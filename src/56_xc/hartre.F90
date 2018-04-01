!{\src2tex{textfont=tt}}
!!****f* ABINIT/hartre
!! NAME
!! hartre
!!
!! FUNCTION
!! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
!! When cplex=1, assume q=(0 0 0), and vhartr will be REAL
!! When cplex=2, q must be taken into account, and vhartr will be COMPLEX
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!         (gsqcut=(boxcut**2)*ecut/(2.d0*(Pi**2))
!!  izero=if 1, unbalanced components of Vhartree(g) are set to zero
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  [qpt(3)=reduced coordinates for a wavevector to be combined with the G vectors (needed if cplex==2).]
!!  rhog(2,nfft)=electron density in G space
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  divgq0= [optional argument] value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq
!!
!! OUTPUT
!!  vhartr(cplex*nfft)=Hartree potential in real space, either REAL or COMPLEX
!!
!! PARENTS
!!      dfpt_rhotov,dfptnl_loop,energy,fock_getghc,m_kxc,nonlinear,nres2vres
!!      odamix,prcref,prcref_PMA,respfn,rhotov,setup_positron,setvtr,tddft
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp,timab,zerosym
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine hartre(cplex,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,rhog,rprimd,vhartr,&
&  divgq0,qpt) ! Optional argument

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_time,     only : timab
 use m_geometry, only : metric
 use m_mpinfo,   only : ptabs_fourdp
 use m_fft,      only : zerosym

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hartre'
 use interfaces_18_timing
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,izero,nfft,paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),rhog(2,nfft)
 real(dp),intent(in),optional :: divgq0
 real(dp),intent(in),optional :: qpt(3)
 real(dp),intent(out) :: vhartr(cplex*nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i2_local,i3,id1,id2,id3
 integer :: ig,ig1min,ig1,ig1max,ig2,ig2min,ig2max,ig3,ig3min,ig3max
 integer :: ii,ii1,ing,n1,n2,n3,qeq0,qeq05,me_fft,nproc_fft
 real(dp),parameter :: tolfix=1.000000001e0_dp
 real(dp) :: cutoff,den,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3,ucvol
 character(len=500) :: message
!arrays
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gmet(3,3),gprimd(3,3),qpt_(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************

 ! Keep track of total time spent in hartre
 call timab(10,1,tsec)

 ! Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message, '(a,i0,a,a)' )&
&   'From the calling routine, cplex=',cplex,ch10,&
&   'but the only value allowed are 1 and 2.'
   MSG_BUG(message)
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)
 nproc_fft = mpi_enreg%nproc_fft; me_fft = mpi_enreg%me_fft

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 ! Initialize a few quantities
 cutoff=gsqcut*tolfix
 if(present(qpt))then
   qpt_=qpt
 else
   qpt_=zero
 end if
 qeq0=0
 if(qpt_(1)**2+qpt_(2)**2+qpt_(3)**2<1.d-15) qeq0=1
 qeq05=0
 if (qeq0==0) then
   if (abs(abs(qpt_(1))-half)<tol12.or.abs(abs(qpt_(2))-half)<tol12.or. &
&   abs(abs(qpt_(3))-half)<tol12) qeq05=1
 end if

 ! If cplex=1 then qpt_ should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message,'(a,3e12.4,a,a)')&
&   'cplex=1 but qpt=',qpt_,ch10,&
&   'qpt should be 0 0 0.'
   MSG_BUG(message)
 end if

 ! If FFT parallelism then qpt should not be 1/2
 if (nproc_fft>1.and.qeq05==1) then
   write(message, '(a,3e12.4,a,a)' )&
&   'FFT parallelism selected but qpt',qpt_,ch10,&
&   'qpt(i) should not be 1/2...'
   MSG_ERROR(message)
 end if

 ! In order to speed the routine, precompute the components of g+q
 ! Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig+qpt_(ii)
   end do
 end do
 ig1max=-1;ig2max=-1;ig3max=-1
 ig1min=n1;ig2min=n2;ig3min=n3

 ABI_ALLOCATE(work1,(2,nfft))
 id1=n1/2+2;id2=n2/2+2;id3=n3/2+2

 ! Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   ! Precompute some products that do not depend on i2 and i1
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2) == me_fft) then
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12

       i2_local = ffti2_local(i2)
       i23=n1*(i2_local-1 +(n2/nproc_fft)*(i3-1))
       ! Do the test that eliminates the Gamma point outside of the inner loop
       ii1=1
       if(i23==0 .and. qeq0==1  .and. ig2==0 .and. ig3==0)then
         ii1=2
         work1(re,1+i23)=zero
         work1(im,1+i23)=zero
         ! If the value of the integration of the Coulomb singularity 4pi\int_BZ 1/q^2 dq is given, use it
         if (PRESENT(divgq0)) then
           work1(re,1+i23)=rhog(re,1+i23)*divgq0*piinv
           work1(im,1+i23)=rhog(im,1+i23)*divgq0*piinv
         end if
       end if

       ! Final inner loop on the first dimension (note the lower limit)
       do i1=ii1,n1
         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
         ii=i1+i23

         if(gs<=cutoff)then
           ! Identify min/max indexes (to cancel unbalanced contributions later)
           ! Count (q+g)-vectors with similar norm
           if ((qeq05==1).and.(izero==1)) then
             ig1=i1-(i1/id1)*n1-1
             ig1max=max(ig1max,ig1); ig1min=min(ig1min,ig1)
             ig2max=max(ig2max,ig2); ig2min=min(ig2min,ig2)
             ig3max=max(ig3max,ig3); ig3min=min(ig3min,ig3)
           end if

           den=piinv/gs
           work1(re,ii)=rhog(re,ii)*den
           work1(im,ii)=rhog(im,ii)*den
         else
           ! gs>cutoff
           work1(re,ii)=zero
           work1(im,ii)=zero
         end if

       end do ! End loop on i1
     end if
   end do ! End loop on i2
 end do ! End loop on i3

 ABI_DEALLOCATE(gq)

 if (izero==1) then
   ! Set contribution of unbalanced components to zero

   if (qeq0==1) then !q=0
     call zerosym(work1,2,n1,n2,n3,comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)

   else if (qeq05==1) then
     !q=1/2; this doesn't work in parallel
     ig1=-1;if (mod(n1,2)==0) ig1=1+n1/2
     ig2=-1;if (mod(n2,2)==0) ig2=1+n2/2
     ig3=-1;if (mod(n3,2)==0) ig3=1+n3/2
     if (abs(abs(qpt_(1))-half)<tol12) then
       if (abs(ig1min)<abs(ig1max)) ig1=abs(ig1max)
       if (abs(ig1min)>abs(ig1max)) ig1=n1-abs(ig1min)
     end if
     if (abs(abs(qpt_(2))-half)<tol12) then
       if (abs(ig2min)<abs(ig2max)) ig2=abs(ig2max)
       if (abs(ig2min)>abs(ig2max)) ig2=n2-abs(ig2min)
     end if
     if (abs(abs(qpt_(3))-half)<tol12) then
       if (abs(ig3min)<abs(ig3max)) ig3=abs(ig3max)
       if (abs(ig3min)>abs(ig3max)) ig3=n3-abs(ig3min)
     end if
     call zerosym(work1,2,n1,n2,n3,ig1=ig1,ig2=ig2,ig3=ig3,&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
   end if
 end if

 ! Fourier Transform Vhartree. Vh in reciprocal space was stored in work1
 call fourdp(cplex,work1,vhartr,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

 ABI_DEALLOCATE(work1)

 call timab(10,2,tsec)

end subroutine hartre
!!***

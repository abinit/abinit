!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spacepar
!! NAME
!! m_spacepar
!!
!! FUNCTION
!!  Relatively Low-level procedures operating on arrays defined on the FFT box (G- or R- space)
!!  Unlike the procedures in m_cgtools, the routines declared in this module can use mpi_type.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2018 ABINIT group (XG, BA, MT, DRH)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_spacepar

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_xmpi

 use m_time,            only : timab
 use defs_abitypes,     only : MPI_type
 use m_geometry,        only : metric
 use m_mpinfo,          only : ptabs_fourdp
 use m_fft,             only : zerosym

 implicit none

 private
!!***

public :: hartre            ! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
public :: meanvalue_g       ! Compute <wf|op|wf> where op is real and diagonal in G-space.
public :: laplacian         ! Compute the laplacian of a function defined in real space
public :: redgr             ! Compute reduced gradients of a real function on the usual unshifted FFT grid.
public :: hartrestr         ! FFT of (rho(G)/pi)*[d(1/G**2)/d(strain) - delta(diagonal strain)*(1/G**2)]
!!***

contains
!!***

!!****f* m_spacepar/hartre
!! NAME
!! hartre
!!
!! FUNCTION
!! Given rho(G), compute Hartree potential (=FFT of rho(G)/pi/(G+q)**2)
!! When cplex=1, assume q=(0 0 0), and vhartr will be REAL
!! When cplex=2, q must be taken into account, and vhartr will be COMPLEX
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

subroutine hartre(cplex,gsqcut,izero,mpi_enreg,nfft,ngfft,paral_kgb,rhog,rprimd,vhartr,&
&  divgq0,qpt) ! Optional argument


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hartre'
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

!!****f* m_spacepar/meanvalue_g
!! NAME
!! meanvalue_g
!!
!! FUNCTION
!!  Compute the mean value of one wavefunction, in reciprocal space,
!!  for an operator that is real, diagonal in G-space: <wf|op|wf>
!!  For the time being, only spin-independent operators are treated.
!!
!! INPUTS
!!  diag(npw)=diagonal operator (real, spin-independent!)
!!  filter= if 1, need to filter on the value of diag, that must be less than huge(0.0d0)*1.d-11
!!          otherwise, should be 0
!!  istwf_k=storage mode of the vectors
!!  npw=number of planewaves of the vector
!!  nspinor=number of spinor components
!!  vect(2,npw*nspinor)=vector
!!  vect1(2,npw*nspinor*use_ndo)=vector1 (=vector in most of the cases)
!!  use_ndo = says if vect=/vect1
!!
!! OUTPUT
!!  ar=mean value
!!
!! PARENTS
!!      dfpt_vtowfk,energy,forstrnps,vtowfk
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine meanvalue_g(ar,diag,filter,istwf_k,mpi_enreg,npw,nspinor,vect,vect1,use_ndo,ar_im)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'meanvalue_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: filter,istwf_k,npw,nspinor,use_ndo
 real(dp),intent(out) :: ar
 real(dp),intent(out),optional :: ar_im
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 real(dp),intent(in) :: diag(npw),vect(2,npw*nspinor)
 real(dp),intent(in) :: vect1(2,npw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: i1,ierr,ipw,jpw,me_g0
 character(len=500) :: message
!arrays

! *************************************************************************
 me_g0 = mpi_enreg%me_g0

 DBG_CHECK(ANY(filter==(/0,1/)),"Wrong filter")
 DBG_CHECK(ANY(nspinor==(/1,2/)),"Wrong nspinor")
 DBG_CHECK(ANY(istwf_k==(/(ipw,ipw=1,9)/)),"Wrong istwf_k")

 if(nspinor==2 .and. istwf_k/=1)then
   write(message,'(a,a,a,i6,a,i6)')&
&   'When istwf_k/=1, nspinor must be 1,',ch10,&
&   'however, nspinor=',nspinor,', and istwf_k=',istwf_k
   MSG_BUG(message)
 end if

 if(use_ndo==1 .and. (istwf_k==2 .and.me_g0==1)) then
   MSG_BUG('use_ndo==1, not tested')
 end if

 ar=zero
 if(present(ar_im)) ar_im=zero

!Normal storage mode
 if(istwf_k==1)then

!  No filter
   if(filter==0)then
!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end do
     end if
     if(use_ndo==1 .and. nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
       end do
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
       end do
     end if

!    !$OMP PARALLEL DO REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end do
!    end if
   else ! will filter

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
     if(nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         if(diag(jpw)<huge(0.0d0)*1.d-11)then
           ar=ar+diag(jpw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
         end if
       end do
     end if
     if(use_ndo==1 .and. nspinor==2)then
!$OMP PARALLEL DO REDUCTION(+:ar_im)
       do ipw=1,npw
         if(diag(ipw)<huge(0.0d0)*1.d-11)then
           ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end if
       end do
!$OMP PARALLEL DO REDUCTION(+:ar_im) PRIVATE(jpw)
       do ipw=1+npw,2*npw
         jpw=ipw-npw
         if(diag(jpw)<huge(0.0d0)*1.d-11)then
           ar_im=ar_im+diag(jpw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
         end if
       end do
     end if


!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1,npw
!    if(diag(ipw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    if(nspinor==2)then
!    !$OMP PARALLEL DO PRIVATE(ipw) REDUCTION(+:ar,ar_im)
!    do ipw=1+npw,2*npw
!    if(diag(ipw-npw)<huge(0.0d0)*1.d-11)then
!    ar=ar+diag(ipw-npw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
!    if(use_ndo==1.and.nspinor==2) ar_im=ar_im+diag(ipw-npw)*(vect1(1,ipw)*vect(2,ipw)-vect1(2,ipw)*vect(1,ipw))
!    end if
!    end do
!    end if ! nspinor==2

   end if ! filter==0

 else if(istwf_k>=2)then

   if(filter==0)then
     i1=1
     if(istwf_k==2 .and. me_g0==1)then ! MPIWF need to know which proc has G=0
       ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
     end do

   else ! filter/=0
     i1=1
     if(istwf_k==2 .and. me_g0==1)then
       if(diag(1)<huge(0.0d0)*1.d-11)then
         ar=half*diag(1)*vect(1,1)*vect1(1,1) ; i1=2
       end if
     end if

!$OMP PARALLEL DO REDUCTION(+:ar)
     do ipw=i1,npw
       if(diag(ipw)<huge(0.0d0)*1.d-11)then
         ar=ar+diag(ipw)*(vect(1,ipw)*vect1(1,ipw)+vect(2,ipw)*vect1(2,ipw))
       end if
     end do
   end if ! filter==0

   ar=two*ar

 end if ! istwf_k

!MPIWF need to make reduction on ar and ai .
 if(mpi_enreg%paral_kgb==1)then
   call xmpi_sum(ar,mpi_enreg%comm_bandspinorfft ,ierr)
   if(present(ar_im))then
     call xmpi_sum(ar_im,mpi_enreg%comm_bandspinorfft,ierr)
   end if
 end if

end subroutine meanvalue_g
!!***

!!****f* m_spacepar/laplacian
!! NAME
!! laplacian
!!
!! FUNCTION
!! compute the laplacian of a function defined in real space
!! the code is written in the way of /3xc/xcden.F90
!!
!! INPUTS
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  mpi_enreg=information about MPI parallelization
!!  nfft=number of points of the fft grid
!!  nfunc=number of functions on the grid for which the laplacian is to be calculated
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  paral_kgb=flag controlling (k,g,bands) parallelization
!!  (optional) rdfuncr(nfft,nfunc)=real(dp) discretized functions in real space
!!  rdfuncg_in TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_in TO BE DESCRIBED SB 090901
!!  (optional) g2cart_in(nfft) = G**2 on the grid
!!
!! OUTPUT
!! (optional) laplacerdfuncr = laplacian in real space of the functions in rdfuncr
!! (optional) rdfuncg = real(dp) discretized functions in fourier space
!! (optional) laplacerdfuncg = real(dp) discretized laplacian of the functions in fourier space
!! (optional) g2cart_out(nfft) = G**2 on the grid
!!  rdfuncg_out TO BE DESCRIBED SB 090901
!!  laplacerdfuncg_out TO BE DESCRIBED SB 090901
!!
!! PARENTS
!!      frskerker1,frskerker2,moddiel_csrb,prcrskerker1,prcrskerker2
!!
!! CHILDREN
!!      fourdp,ptabs_fourdp
!!
!! SOURCE

subroutine laplacian(gprimd,mpi_enreg,nfft,nfunc,ngfft,paral_kgb,rdfuncr,&
&  laplacerdfuncr,rdfuncg_out,laplacerdfuncg_out,g2cart_out,rdfuncg_in,g2cart_in)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'laplacian'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,nfunc,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout),optional :: laplacerdfuncr(nfft,nfunc)
 real(dp),intent(inout),optional,target :: rdfuncr(nfft,nfunc)
 real(dp),intent(in),optional,target :: g2cart_in(nfft) !vz_i
 real(dp),intent(out),optional,target :: g2cart_out(nfft)  !vz_i
 real(dp),intent(out),optional,target :: laplacerdfuncg_out(2,nfft,nfunc)
 real(dp),intent(in),optional,target :: rdfuncg_in(2,nfft,nfunc) !vz_i
 real(dp),intent(out),optional,target :: rdfuncg_out(2,nfft,nfunc)

!Local variables-------------------------------
!scalars
 integer :: count,i1,i2,i3,id1,id2,id3,ifft,ifunc,ig1,ig2,ig3,ii1,n1,n2
 integer :: n3
 real(dp) :: b11,b12,b13,b21,b22,b23,b31,b32,b33
!arrays
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp),ABI_CONTIGUOUS pointer :: g2cart(:),laplacerdfuncg(:,:,:),rdfuncg(:,:,:)

! *************************************************************************

!Keep local copy of fft dimensions
 n1=ngfft(1)
 n2=ngfft(2)
 n3=ngfft(3)

 if(present(laplacerdfuncg_out)) then
   laplacerdfuncg => laplacerdfuncg_out
 else
   ABI_ALLOCATE(laplacerdfuncg,(2,nfft,nfunc))
 end if

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!change the real density rdfuncr on real space on the real density
!rdfuncg in reciprocal space
 if(.not.present(rdfuncg_in)) then
   if(present(rdfuncg_out)) then
     rdfuncg => rdfuncg_out
   else
     ABI_ALLOCATE(rdfuncg,(2,nfft,nfunc))
   end if
   if(present(rdfuncr)) then
     do ifunc=1,nfunc
       call fourdp(1,rdfuncg(:,:,ifunc),rdfuncr(:,ifunc),-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     end do
   end if
 else
   rdfuncg => rdfuncg_in
 end if

!apply the laplacian on laplacerdfuncr
!code from /3xc/xcden.F90
!see also src/5common/hatre.F90 and src/5common/moddiel.F90
!Keep local copy of fft dimensions
!Initialize computation of G^2 in cartesian coordinates
 if(.not.present(g2cart_in)) then
   if(present(g2cart_out)) then
     g2cart => g2cart_out
   else
     ABI_ALLOCATE(g2cart,(nfft))
   end if
   id1=int(n1/2)+2
   id2=int(n2/2)+2
   id3=int(n3/2)+2
   count=0
   do i3=1,n3
     ifft=(i3-1)*n1*(n2/mpi_enreg%nproc_fft)
     ig3=i3-int(i3/id3)*n3-1
     do i2=1,n2
       if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
         ig2=i2-int(i2/id2)*n2-1

         ii1=1
         do i1=ii1,n1
           ig1=i1-int(i1/id1)*n1-1
           ifft=ifft+1

           b11=gprimd(1,1)*real(ig1,dp)
           b21=gprimd(2,1)*real(ig1,dp)
           b31=gprimd(3,1)*real(ig1,dp)
           b12=gprimd(1,2)*real(ig2,dp)
           b22=gprimd(2,2)*real(ig2,dp)
           b32=gprimd(3,2)*real(ig2,dp)
           b13=gprimd(1,3)*real(ig3,dp)
           b23=gprimd(2,3)*real(ig3,dp)
           b33=gprimd(3,3)*real(ig3,dp)

           g2cart(ifft)=( &
&           (b11+b12+b13)**2&
&           +(b21+b22+b23)**2&
&           +(b31+b32+b33)**2&
&           )
           do ifunc=1,nfunc
!            compute the laplacian in Fourier space that is * (i x 2pi x G)**2
             laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
             laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
           end do
         end do
       end if
     end do
   end do
   if(.not.present(g2cart_out))  then
     ABI_DEALLOCATE(g2cart)
   end if
 else
   g2cart => g2cart_in
   do ifunc=1,nfunc
     do ifft=1,nfft
!      compute the laplacian in Fourier space that is * (i x 2pi x G)**2
       laplacerdfuncg(1,ifft,ifunc) = -rdfuncg(1,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
       laplacerdfuncg(2,ifft,ifunc) = -rdfuncg(2,ifft,ifunc)*g2cart(ifft)*two_pi*two_pi
     end do
   end do
 end if

!get the result back into real space
 if(present(laplacerdfuncr)) then
   do ifunc=1,nfunc
     call fourdp(1,laplacerdfuncg(:,:,ifunc),laplacerdfuncr(:,ifunc),1,mpi_enreg,nfft,ngfft,paral_kgb,0)
   end do
 end if

!deallocate pointers
 if((.not.present(rdfuncg_in)).and.(.not.present(rdfuncg_in)))  then
   ABI_DEALLOCATE(rdfuncg)
 end if
 if(.not.present(laplacerdfuncg_out))  then
   ABI_DEALLOCATE(laplacerdfuncg)
 end if

end subroutine laplacian
!!***

!!****f* m_spacepar/redgr
!! NAME
!! redgr
!!
!! FUNCTION
!! Compute reduced gradients of a real function on the usual unshifted
!! fft grid. The gradient directions are the along the primitive
!! reciprocal lattice vectors.
!! The input function is intended to be a single spin component of
!! the valence charge density, the valence + core charge densities
!! or the first-order core charge density for use in frozen wf
!! elastic tensor calculations within the GGA.
!!
!! NOTES
!! Closely linked to xcden, but limited to Q=0, real charge densities,
!! and unshifted grids.
!!
!! INPUTS
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!   see ~abinit/doc/variables/vargs.htm#ngfft
!!  frin(nfft)=real space input function
!!
!! OUTPUT
!!  frredgr(nfft,3)= reduced gradient of input function (same units as frin)
!!
!! PARENTS
!!      dfpt_eltfrxc
!!
!! CHILDREN
!!      fourdp
!!
!! SOURCE

subroutine redgr (frin,frredgr,mpi_enreg,nfft,ngfft,paral_kgb)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'redgr'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: frin(nfft)
 real(dp),intent(out) :: frredgr(nfft,3)

!Local variables-------------------------------
!scalars
 integer :: cplex_tmp,i1,i2,i3,id,idir,ifft,ig,ii,ing,n1,n2,n3
!arrays
 real(dp),allocatable :: gg(:,:),wkcmpx(:,:),work(:),workgr(:,:)

! *************************************************************************

!Only real arrays are treated
 cplex_tmp=1

!Keep local copy of fft dimensions
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!In order to speed the routine, precompute the components of g, including 2pi factor
 ABI_ALLOCATE(gg,(max(n1,n2,n3),3))
 do ii=1,3
   id=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id)*ngfft(ii)-1
     gg(ing,ii)=two_pi*ig
   end do
!  Note that the G <-> -G symmetry must be maintained
   if(mod(ngfft(ii),2)==0)gg(ngfft(ii)/2+1,ii)=zero
 end do

 ABI_ALLOCATE(wkcmpx,(2,nfft))
 ABI_ALLOCATE(work,(nfft))
 ABI_ALLOCATE(workgr,(2,nfft))

!Obtain rho(G) in wkcmpx from input rho(r)
 work(:)=frin(:)

 call fourdp(cplex_tmp,wkcmpx,work,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!Gradient calculation for three reduced components in turn.
!Code duplicated to remove logic from loops.
 do idir=1,3
   if(idir==1) then
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i1,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i1,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   else if(idir==2) then
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i2,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i2,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   else
!$OMP PARALLEL DO PRIVATE(ifft)
     do i3=1,n3
       ifft=(i3-1)*n1*n2
       do i2=1,n2
         do i1=1,n1
           ifft=ifft+1
!          Multiply by i 2pi G(idir)
           workgr(2,ifft)= gg(i3,idir)*wkcmpx(1,ifft)
           workgr(1,ifft)=-gg(i3,idir)*wkcmpx(2,ifft)
         end do
       end do
     end do
   end if !idir

   call fourdp(cplex_tmp,workgr,work,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

!$OMP PARALLEL DO
   do ifft=1,nfft
     frredgr(ifft,idir)=work(ifft)
   end do

 end do !idir

 ABI_DEALLOCATE(gg)
 ABI_DEALLOCATE(wkcmpx)
 ABI_DEALLOCATE(work)
 ABI_DEALLOCATE(workgr)

end subroutine redgr
!!***

!!****f* m_spacepar/hartrestr
!! NAME
!! hartrestr
!!
!! FUNCTION
!! To be called for strain perturbation only
!! Compute the inhomogenous terms generated by the strain derivative of
!! Hartree potential due to the ground state charge rho(G)
!!
!!  FFT of (rho(G)/pi)*[d(1/G**2)/d(strain) - delta(diagonal strain)*(1/G**2)]
!!
!! NOTES
!! *based largely on hartre.f
!! *Modified code to avoid if statements inside loops to skip G=0.
!!  Replaced if statement on G^2>gsqcut to skip G s outside where
!!  rho(G) should be 0.  Effect is negligible but gsqcut should be
!!  used to be strictly consistent with usage elsewhere in code.
!! *The speed-up is provided by doing a few precomputations outside
!!  the inner loop. One variable size array is needed for this (gq).
!!
!! INPUTS
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=number of fft grid points (gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2))
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  vhartr1(nfft)=Inhomogeneous term in strain-perturbation-induced Hartree
!!   potential in real space,
!!
!! PARENTS
!!      dfpt_nselt,dfpt_nstpaw,dfpt_rhotov
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp
!!
!! SOURCE

subroutine hartrestr(gsqcut,idir,ipert,mpi_enreg,natom,nfft,ngfft,&
&  paral_kgb,rhog,rprimd,vhartr1)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hartrestr'
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: idir,ipert,natom,nfft,paral_kgb
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhog(2,nfft),rprimd(3,3)
 real(dp),intent(out) :: vhartr1(nfft)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,id2,id3,ig,ig2,ig3,ii,ii1,ing,istr,ka,kb,n1,n2,n3
 real(dp),parameter :: tolfix=1.000000001_dp
 real(dp) :: cutoff,ddends,den,dgsds,gqg2p3,gqgm12,gqgm13,gqgm23,gs,gs2,gs3
 real(dp) :: term,ucvol
 character(len=500) :: message
!arrays
 integer,save :: idx(12)=(/1,1,2,2,3,3,3,2,3,1,2,1/)
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: dgmetds(3,3),gmet(3,3),gprimd(3,3),gqr(3),rmet(3,3)
 real(dp),allocatable :: gq(:,:),work1(:,:)

! *************************************************************************

 if( .not. (ipert==natom+3 .or. ipert==natom+4))then
   write(message, '(a,i0,a,a)' )&
&   'From the calling routine, ipert=',ipert,ch10,&
&   'so this routine for the strain perturbation should not be called.'
   MSG_BUG(message)
 end if

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Initialize a few quantities
 cutoff=gsqcut*tolfix

 istr=idir + 3*(ipert-natom-3)

 if(istr<1 .or. istr>6)then
   write(message, '(a,i10,a,a,a)' )&
&   'Input dir gives istr=',istr,' not allowed.',ch10,&
&   'Possible values are 1,2,3,4,5,6 only.'
   MSG_BUG(message)
 end if

 ka=idx(2*istr-1);kb=idx(2*istr)
 do ii = 1,3
   dgmetds(:,ii)=-(gprimd(ka,:)*gprimd(kb,ii)+gprimd(kb,:)*gprimd(ka,ii))
 end do
!For historical reasons:
 dgmetds(:,:)=0.5_dp*dgmetds(:,:)

!In order to speed the routine, precompute the components of g+q
!Also check if the booked space was large enough...
 ABI_ALLOCATE(gq,(3,max(n1,n2,n3)))
 do ii=1,3
   id(ii)=ngfft(ii)/2+2
   do ing=1,ngfft(ii)
     ig=ing-(ing/id(ii))*ngfft(ii)-1
     gq(ii,ing)=ig
   end do
 end do

 ABI_ALLOCATE(work1,(2,nfft))
 id2=n2/2+2
 id3=n3/2+2
!Triple loop on each dimension
 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
!  Precompute some products that do not depend on i2 and i1
   gqr(3)=gq(3,i3)
   gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
   gqgm23=gq(3,i3)*gmet(2,3)*2
   gqgm13=gq(3,i3)*gmet(1,3)*2

   do i2=1,n2
     if (fftn2_distrib(i2)==mpi_enreg%me_fft) then
       gqr(2)=gq(2,i2)
       gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
       gqgm12=gq(2,i2)*gmet(1,2)*2
       gqg2p3=gqgm13+gqgm12
       ig2=i2-(i2/id2)*n2-1
!      i23=n1*((i2-1)+n2*(i3-1))
       i23=n1*((ffti2_local(i2)-1)+(n2/mpi_enreg%nproc_fft)*(i3-1))
!      Do the test that eliminates the Gamma point outside
!      of the inner loop
       ii1=1
       if(i23==0  .and. ig2==0 .and. ig3==0)then
         ii1=2
         work1(re,1+i23)=0.0_dp
         work1(im,1+i23)=0.0_dp
       end if

!      Final inner loop on the first dimension
!      (note the lower limit)
       do i1=ii1,n1
         gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
         ii=i1+i23
         if(gs<=cutoff)then
           den=piinv/gs
           gqr(1)=gq(1,i1)
           dgsds=&
&           (gqr(1)*(dgmetds(1,1)*gqr(1)+dgmetds(1,2)*gqr(2)+dgmetds(1,3)*gqr(3))+  &
&           gqr(2)*(dgmetds(2,1)*gqr(1)+dgmetds(2,2)*gqr(2)+dgmetds(2,3)*gqr(3))+  &
&           gqr(3)*(dgmetds(3,1)*gqr(1)+dgmetds(3,2)*gqr(2)+dgmetds(3,3)*gqr(3)) )
           ddends=-piinv*dgsds/gs**2
           if(istr<=3)then
             term=2.0_dp*ddends-den
           else
             term=2.0_dp*ddends
           end if
           work1(re,ii)=rhog(re,ii)*term
           work1(im,ii)=rhog(im,ii)*term
         else
           work1(re,ii)=0.0_dp
           work1(im,ii)=0.0_dp
         end if

       end do ! End loop on i1
     end if
   end do ! End loop on i2
 end do !  End loop on i3

 ABI_DEALLOCATE(gq)

!Fourier Transform Vhartree.
!Vh in reciprocal space was stored in work1
 call fourdp(1,work1,vhartr1,1,mpi_enreg,nfft,ngfft,paral_kgb,0)

 ABI_DEALLOCATE(work1)

end subroutine hartrestr
!!***

end module m_spacepar
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/moddiel
!! NAME
!! moddiel
!!
!! FUNCTION
!! Precondition the residual, using a model dielectric function.
!! When cplex=1, assume q=(0 0 0), and vresid and vrespc will be REAL
!! When cplex=2, q must be taken into account, and vresid and vrespc will be COMPLEX
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, vhartr is REAL, if 2, vhartr is COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  optreal=1 if residual potential is in REAL space, 2 if it is in RECIPROCAL SPACE
!!  optres=0: the array vresid contains a potential residual
!!         1: the array vresid contains a density residual
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  vresid(cplex*nfft,nspden)=residual density/potential in REAL space       (if optreal==1)
!!                            residual density/potential in RECIPROCAL space (if optreal==2)
!!
!! OUTPUT
!!  vrespc(cplex*nfft,nspden)=preconditioned residual of the density/potential
!!                            in REAL space       if optreal==1
!!                            in RECIPROCAL space if optreal==2
!!
!! SIDE EFFECTS
!!
!! NOTES
!! optreal==2 is not compatible with cplex==1
!!
!! PARENTS
!!      dfpt_newvtr,prcref,prcref_PMA
!!
!! CHILDREN
!!      fourdp,metric,ptabs_fourdp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine moddiel(cplex,dielar,mpi_enreg,nfft,ngfft,nspden,optreal,optres,paral_kgb,qphon,rprimd,vresid,vrespc)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

 use m_mpinfo,   only : ptabs_fourdp

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'moddiel'
 use interfaces_41_geometry
 use interfaces_53_ffts
!End of the abilint section

 implicit none

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nspden,optreal,optres,paral_kgb
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dielar(7),qphon(3),rprimd(3,3)
 real(dp),intent(in) :: vresid(cplex*nfft,nspden)
 real(dp),intent(out) :: vrespc(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i23,i3,ifft,ig,ii,ii1,ing,ispden,me_fft,mg,n1,n2,n3,nproc_fft
 integer :: nspden_eff,qeq0
 logical :: magn_precon
 real(dp) :: dielng,diemac,diemac_inv,diemix,diemixmag,diemix_eff,factor,gqg2p3,gqgm12,gqgm13
 real(dp) :: gqgm23,gs,gs2,gs3,l2g2,length2,ucvol
 character(len=500) :: message
!arrays
 integer :: id(3)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)
 real(dp) :: gmet(3,3),gprimd(3,3),potg0(4),rmet(3,3)
 real(dp),allocatable :: gq(:,:),work1(:,:),work2(:)

! *************************************************************************

!Check that cplex has an allowed value
 if(cplex/=1 .and. cplex/=2)then
   write(message,'(a,i0,a,a)')&
&   '  From the calling routine, cplex=',cplex,ch10,&
&   '  but the only value allowed are 1 and 2.'
   MSG_BUG(message)
 end if

 if(cplex==1.and.optreal==2)then
   MSG_BUG('When optreal=2, cplex must be 2.')
 end if

!This is to allow q=0
 qeq0=0
 if(qphon(1)**2+qphon(2)**2+qphon(3)**2<1.d-15) qeq0=1

!If cplex=1 then qphon should be 0 0 0
 if (cplex==1.and. qeq0/=1) then
   write(message,'(a,3e12.4,a)' )' cplex=1 but qphon=',qphon,' qphon should be 0 0 0.'
   MSG_BUG(message)
 end if

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)

!Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 dielng=dielar(2) ; diemac=dielar(3) ; diemix=dielar(4) ; diemixmag=dielar(7)

 magn_precon=(diemixmag>=zero) ! Set to true if magnetization has to be preconditionned
 diemixmag=abs(diemixmag)

!DEBUG
!write(std_out,*)' moddiel : diemac, diemix, diemixmag =',diemac,diemix,diemixmag
!ENDDEBUG

 if(abs(diemac-1.0_dp)<1.0d-6)then

!  Here, simple mixing is required, through macroscopic
!  dielectric constant set to 1.0_dp .
   vrespc(:,1)=diemix*vresid(:,1)
   if (nspden/=1) vrespc(:,2:nspden)=diemixmag*vresid(:,2:nspden)
 else

!  Magnetization is not preconditionned
   if (optres==1.and.nspden>1.and.(.not.magn_precon)) vrespc(:,2:nspden)=diemixmag*vresid(:,2:nspden)

!  Here, model dielectric function (G-diagonal operator)

   length2=(two_pi*dielng)**2
   diemac_inv=1.0_dp/diemac
   ABI_ALLOCATE(work1,(2,nfft))
   if (optreal==1) then
     ABI_ALLOCATE(work2,(cplex*nfft))
   end if

!  In order to speed the routine, precompute the components of g
   mg=maxval(ngfft)
   ABI_ALLOCATE(gq,(3,mg))
   do ii=1,3
     id(ii)=ngfft(ii)/2+2
     do ing=1,ngfft(ii)
       ig=ing-(ing/id(ii))*ngfft(ii)-1
       gq(ii,ing)=ig+qphon(ii)
     end do
   end do

!  Do-loop on spins
!  Note XG 010922 : I doubt the preconditioner is OK for the magnetization
   nspden_eff=nspden;if (optres==1.and.(.not.magn_precon)) nspden_eff=1
   do ispden=1,nspden_eff

     diemix_eff=diemix;if (ispden>1) diemix_eff=diemixmag

!    Do fft from real space (work2) to G space (work1)
     if (optreal==1) then
       work2(:)=vresid(:,ispden)
       call fourdp(cplex,work1,work2,-1,mpi_enreg,nfft,ngfft,paral_kgb,0)
     else
!      work1(:,:)=reshape(vresid(:,ispden),(/2,nfft/))
!      Reshape function does not work with big arrays for some compilers
       do ifft=1,nfft
         work1(1,ifft)=vresid(2*ifft-1,ispden)
         work1(2,ifft)=vresid(2*ifft  ,ispden)
       end do
     end if

!    Store G=0 value
     potg0(ispden)=work1(re,1)

!    Triple loop, for the three dimensions
     do i3=1,n3
!      Precompute some products that do not depend on i2 and i1
       gs3=gq(3,i3)*gq(3,i3)*gmet(3,3)
       gqgm23=gq(3,i3)*gmet(2,3)*2
       gqgm13=gq(3,i3)*gmet(1,3)*2
       do i2=1,n2
         if (fftn2_distrib(i2)==me_fft) then
           gs2=gs3+ gq(2,i2)*(gq(2,i2)*gmet(2,2)+gqgm23)
           gqgm12=gq(2,i2)*gmet(1,2)*2
           gqg2p3=gqgm13+gqgm12
           i23=n1*( ffti2_local(i2)-1+(n2/nproc_fft)*(i3-1))

!          Do the test that eliminates the Gamma point outside
!          of the inner loop
           ii1=1
           if(i2 == 1 .and. i3 == 1 .and. qeq0==1)then
!            if(i23==0 .and. qeq0==1)then: this changes with the number of fft procs...
!            and seems to be wrong.Pls check
             ii1=2
           end if

!          Here, unlike in hartre.f, the G=0 term is not eliminated, albeit
!          not changed.
           do i1=ii1,n1

!            One obtains the square of the norm of q+G (defined by i1,i2,i3)
             gs=gs2+ gq(1,i1)*(gq(1,i1)*gmet(1,1)+gqg2p3)
             ifft=i1+i23

             l2g2=length2*gs
!            The model dielectric function is now computed
             factor = (l2g2+diemac_inv)/(l2g2+1.0_dp) * diemix_eff
             work1(re,ifft)=work1(re,ifft)*factor
             work1(im,ifft)=work1(im,ifft)*factor

           end do
         end if
       end do
     end do

!    Might get rid of the G=0 term
!    if(qeq0==1)then
!    work1(re,1)=0.0_dp
!    work1(im,1)=0.0_dp
!    end if

!    Fourier transform
     if (optreal==1) then
       call fourdp(cplex,work1,work2,1,mpi_enreg,nfft,ngfft,paral_kgb,0)
       vrespc(:,ispden)=work2(:)
     else
!      vrespc(:,ispden)=reshape(work1(:,:),(/nfft*2/))
!      Reshape function does not work with big arrays for some compilers
       do ifft=1,nfft
         vrespc(2*ifft-1,ispden)=work1(1,ifft)
         vrespc(2*ifft  ,ispden)=work1(2,ifft)
       end do
     end if

!    End of loop on spin polarizations
   end do

   ABI_DEALLOCATE(gq)
   ABI_DEALLOCATE(work1)
   if (optreal==1) then
     ABI_DEALLOCATE(work2)
   end if

!  End condition diemac/=1.0
 end if

end subroutine moddiel
!!***

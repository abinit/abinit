!{\src2tex{textfont=tt}}
!!****f* ABINIT/redgr
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
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DRH, DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine redgr (frin,frredgr,mpi_enreg,nfft,ngfft,paral_kgb)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi

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

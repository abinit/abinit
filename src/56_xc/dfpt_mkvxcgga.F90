!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkvxcgga
!! NAME
!! dfpt_mkvxcgga
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2016 ABINIT group (XG, DRH)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=informations about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT
!!  nhat1(cplex*nfft,2*nhat1dim)= -PAW only- 1st-order compensation density
!!  nhat1dim= -PAW only- 1 if nhat1 array is used ; 0 otherwise
!!  nhat1gr(cplex*nfft,2,3*nhat1grdim)= -PAW only- gradients of 1st-order compensation density
!!  nhat1grdim= -PAW only- 1 if nhat1gr array is used ; 0 otherwise
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1tmp(cplex*nfft,2)=array for first-order electron spin-density
!!   in electrons/bohr**3 (second index corresponds to spin-up and spin-down)
!!  usexcnhat= -PAW only- 1 if nhat density has to be taken into account in Vxc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential
!!
!! SIDE EFFECTS
!!
!! NOTES
!!  For the time being, a rather crude coding, to be optimized ...
!!
!! PARENTS
!!      dfpt_mkvxc
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_mkvxcgga(cplex,gprimd,kxc,mpi_enreg,nfft,ngfft,&
&                    nhat1,nhat1dim,nhat1gr,nhat1grdim,nkxc,&
&                    nspden,paral_kgb,qphon,rhor1tmp,usexcnhat,vxc1)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxcgga'
 use interfaces_56_xc, except_this_one => dfpt_mkvxcgga
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,nfft,nhat1dim,nhat1grdim,nkxc,nspden,paral_kgb,usexcnhat
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat1(cplex*nfft,2*nhat1dim)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,2,3*nhat1grdim)
 real(dp),intent(in) :: qphon(3)
 real(dp),intent(in),target :: rhor1tmp(cplex*nfft,2)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ispden,mgga,ngrad,nspdentmp,nspgrad
 logical :: test_nhat
 real(dp) :: coeff_grho_corr,coeff_grho_dn,coeff_grho_up,coeffim_grho_corr
 real(dp) :: coeffim_grho_dn,coeffim_grho_up,gradrho_gradrho1
 real(dp) :: gradrho_gradrho1_dn,gradrho_gradrho1_up,gradrho_gradrho1im
 real(dp) :: gradrho_gradrho1im_dn,gradrho_gradrho1im_up
! character(len=500) :: message
!arrays
 real(dp) :: r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp) :: r1im(3),r1im_dn(3),r1im_up(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:),rhortmp(:,:,:)
 real(dp),allocatable :: vxc1tmp(:,:)
 real(dp),ABI_CONTIGUOUS pointer :: rhor1_ptr(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!metaGGA contributions are not taken into account here
 mgga=0

!Treat all cases in a generic way (to be optimized !!)
 nspdentmp=2

!PAW: substract 1st-order compensation density from 1st-order density
 test_nhat=((nhat1dim==1).and.(usexcnhat==0.or.nhat1grdim==1))
 if (test_nhat) then
   ABI_ALLOCATE(rhor1_ptr,(cplex*nfft,2))
   rhor1_ptr(:,:)=rhor1tmp(:,:)-nhat1(:,:)
 else
   rhor1_ptr => rhor1tmp
 end if

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1_ptr)

!Compute the gradients of the first-order density
!rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density
 ishift=0 ; ngrad=2 ; nspdentmp=2
 ABI_ALLOCATE(rho1now,(cplex*nfft,nspdentmp,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,paral_kgb,qphon,rhor1_ptr,rho1now)

!PAW: add "exact" gradients of compensation density
 if (test_nhat.and.usexcnhat==1) then
   rho1now(:,1:nspdentmp,1)=rho1now(:,1:nspdentmp,1)+nhat1(:,1:nspdentmp)
 end if
 if (nhat1grdim==1) then
   do ii=1,3
     rho1now(:,1:nspdentmp,ii+1)=rho1now(:,1:nspdentmp,ii+1)+nhat1gr(:,1:nspdentmp,ii)
   end do
 end if
 if (test_nhat) then
   ABI_DEALLOCATE(rhor1_ptr)
 end if

!Transfer the ground-state density and its gradient to spin-polarized storage
!rhortmp(:,:,1) contains the GS density, and
!rhortmp(:,:,2:4) contains the gradients of the GS density
 ABI_ALLOCATE(rhortmp,(nfft,nspdentmp,4))
 if(nspden==1)then
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,14+2*ii)*half
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)*half
     end do
   end do
 else
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,15+2*ii)
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)-kxc(ir,15+2*ii)
     end do
   end do
 end if

!Apply the XC kernel
 nspgrad=5
 ABI_ALLOCATE(dnexcdn,(cplex*nfft,nspgrad))

 if(cplex==1)then  ! Treat real case first

   do ir=1,nfft
     r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
     r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
     r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
     r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
     r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
     r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
     gradrho_gradrho1_up=r1_up(1)*r0_up(1)+r1_up(2)*r0_up(2)+r1_up(3)*r0_up(3)
     gradrho_gradrho1_dn=r1_dn(1)*r0_dn(1)+r1_dn(2)*r0_dn(2)+r1_dn(3)*r0_dn(3)
     gradrho_gradrho1   =r1(1)*r0(1)+r1(2)*r0(2)+r1(3)*r0(3)
     dnexcdn(ir,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(ir,1,1)+&
&     kxc(ir,10)*rho1now(ir,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1_up+&
&     kxc(ir,13)*gradrho_gradrho1
     dnexcdn(ir,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(ir,2,1)+&
&     kxc(ir,10)*rho1now(ir,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1_dn+&
&     kxc(ir,14)*gradrho_gradrho1
     coeff_grho_corr=(kxc(ir,13)*rho1now(ir,1,1)+kxc(ir,14)*rho1now(ir,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1
     coeff_grho_up= kxc(ir,5)*rho1now(ir,1,1)+kxc(ir,7)*gradrho_gradrho1_up
     coeff_grho_dn= kxc(ir,6)*rho1now(ir,2,1)+kxc(ir,8)*gradrho_gradrho1_dn
!    Reuse the storage in rho1now
     rho1now(ir,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeff_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(ir,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeff_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
   end do

 else ! if cplex==2

   do ir=1,nfft
     r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
     r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
     r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
     r1_up(:)=rho1now(2*ir-1,1,2:4)   ! grad of spin-up rho1
     r1im_up(:)=rho1now(2*ir,1,2:4)   ! grad of spin-up rho1 , im part
     r1_dn(:)=rho1now(2*ir-1,2,2:4)   ! grad of spin-down rho1
     r1im_dn(:)=rho1now(2*ir,2,2:4)   ! grad of spin-down rho1 , im part
     r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
     r1im(:)=r1im_up(:)+r1im_dn(:)      ! grad of GS rho1, im part
     gradrho_gradrho1_up=r1_up(1)*r0_up(1)+r1_up(2)*r0_up(2)+r1_up(3)*r0_up(3)
     gradrho_gradrho1_dn=r1_dn(1)*r0_dn(1)+r1_dn(2)*r0_dn(2)+r1_dn(3)*r0_dn(3)
     gradrho_gradrho1   =r1(1)*r0(1)+r1(2)*r0(2)+r1(3)*r0(3)
     gradrho_gradrho1im_up=r1im_up(1)*r0_up(1)+r1im_up(2)*r0_up(2)+r1im_up(3)*r0_up(3)
     gradrho_gradrho1im_dn=r1im_dn(1)*r0_dn(1)+r1im_dn(2)*r0_dn(2)+r1im_dn(3)*r0_dn(3)
     gradrho_gradrho1im   =r1im(1)*r0(1)+r1im(2)*r0(2)+r1im(3)*r0(3)
     dnexcdn(2*ir-1,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir-1,1,1)+&
&     kxc(ir,10)*rho1now(2*ir-1,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1_up+&
&     kxc(ir,13)*gradrho_gradrho1
     dnexcdn(2*ir  ,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir,1,1)+&
&     kxc(ir,10)*rho1now(2*ir,2,1)+&
&     kxc(ir,5)*gradrho_gradrho1im_up+&
&     kxc(ir,13)*gradrho_gradrho1im
     dnexcdn(2*ir-1,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir-1,2,1)+&
&     kxc(ir,10)*rho1now(2*ir-1,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1_dn+&
&     kxc(ir,14)*gradrho_gradrho1
     dnexcdn(2*ir  ,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir,2,1)+&
&     kxc(ir,10)*rho1now(2*ir,1,1)+&
&     kxc(ir,6)*gradrho_gradrho1im_dn+&
&     kxc(ir,14)*gradrho_gradrho1im
     coeff_grho_corr=(kxc(ir,13)*rho1now(2*ir-1,1,1)+kxc(ir,14)*rho1now(2*ir-1,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1
     coeffim_grho_corr=(kxc(ir,13)*rho1now(2*ir,1,1)+kxc(ir,14)*rho1now(2*ir,2,1))+&
&     kxc(ir,15)*gradrho_gradrho1im
     coeff_grho_up= kxc(ir,5)*rho1now(2*ir-1,1,1)+kxc(ir,7)*gradrho_gradrho1_up
     coeffim_grho_up= kxc(ir,5)*rho1now(2*ir,1,1)+kxc(ir,7)*gradrho_gradrho1im_up
     coeff_grho_dn= kxc(ir,6)*rho1now(2*ir-1,2,1)+kxc(ir,8)*gradrho_gradrho1_dn
     coeffim_grho_dn= kxc(ir,6)*rho1now(2*ir,2,1)+kxc(ir,8)*gradrho_gradrho1im_dn
!    Reuse the storage in rho1now
     rho1now(2*ir-1,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeff_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(2*ir  ,1,2:4)=r1im_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&     +r1im_dn(:)*kxc(ir,12)               &
&     +r0_up(:)*coeffim_grho_up            &
&     +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
     rho1now(2*ir-1,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeff_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeff_grho_corr
     rho1now(2*ir  ,2,2:4)=r1im_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&     +r1im_up(:)*kxc(ir,12)               &
&     +r0_dn(:)*coeffim_grho_dn            &
&     +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
   end do

 end if

 ABI_ALLOCATE(vxc1tmp,(cplex*nfft,nspdentmp))
 vxc1tmp(:,:)=zero
 call xcpot (cplex,dnexcdn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspdentmp,&
& nspgrad,paral_kgb,qphon,rho1now,vxc1tmp)

!Transfer the data from spin-polarized storage
 do ispden=1,nspden
   do ir=1,cplex*nfft
     vxc1(ir,ispden)=vxc1tmp(ir,ispden)
   end do
 end do

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 ABI_DEALLOCATE(dnexcdn)
 ABI_DEALLOCATE(rhortmp)
 ABI_DEALLOCATE(rho1now)
 ABI_DEALLOCATE(vxc1tmp)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxcgga
!!***

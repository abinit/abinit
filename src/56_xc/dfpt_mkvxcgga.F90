!{\src2tex{textfont=tt}
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
!! Copyright (C) 2001-2017 ABINIT group (XG, DRH)
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
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see below)
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
!!  Content of Kxc array:
!!   ===== if GGA
!!    if nspden==1:
!!       kxc(:,1)= d2Exc/drho2
!!       kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!       kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!       kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!       kxc(:,5)= gradx(rho)
!!       kxc(:,6)= grady(rho)
!!       kxc(:,7)= gradz(rho)
!!    if nspden>=2:
!!       kxc(:,1)= d2Ex/drho_up drho_up
!!       kxc(:,2)= d2Ex/drho_dn drho_dn
!!       kxc(:,3)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,4)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,5)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,6)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,7)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,8)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,9)= d2Ec/drho_up drho_up
!!       kxc(:,10)=d2Ec/drho_up drho_dn
!!       kxc(:,11)=d2Ec/drho_dn drho_dn
!!       kxc(:,12)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,13)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,14)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,15)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,16)=rho_up
!!       kxc(:,17)=rho_dn
!!       kxc(:,18)=gradx(rho_up)
!!       kxc(:,19)=gradx(rho_dn)
!!       kxc(:,20)=grady(rho_up)
!!       kxc(:,21)=grady(rho_dn)
!!       kxc(:,22)=gradz(rho_up)
!!       kxc(:,23)=gradz(rho_dn)
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
&                    nspden,paral_kgb,qphon,rhor1,usexcnhat,vxc1)

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
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: kxc(nfft,nkxc)
 real(dp),intent(in) :: nhat1(cplex*nfft,nspden*nhat1dim)
 real(dp),intent(in) :: nhat1gr(cplex*nfft,nspden,3*nhat1grdim)
 real(dp),intent(in) :: qphon(3)
 real(dp),intent(in),target :: rhor1(cplex*nfft,nspden)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ispden,mgga,ngrad,nspgrad
 logical :: test_nhat
 real(dp) :: coeff_grho,coeff_grho_corr,coeff_grho_dn,coeff_grho_up
 real(dp) :: coeffim_grho,coeffim_grho_corr,coeffim_grho_dn,coeffim_grho_up
 real(dp) :: gradrho_gradrho1,gradrho_gradrho1_dn,gradrho_gradrho1_up
 real(dp) :: gradrho_gradrho1im,gradrho_gradrho1im_dn,gradrho_gradrho1im_up
!arrays
 real(dp) :: r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp) :: r1im(3),r1im_dn(3),r1im_up(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:)
 real(dp),allocatable :: rhortmp(:,:,:)
 real(dp),ABI_CONTIGUOUS pointer :: rhor1_ptr(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

!metaGGA contributions are not taken into account here
 mgga=0

!PAW: substract 1st-order compensation density from 1st-order density
 test_nhat=((nhat1dim==1).and.(usexcnhat==0.or.nhat1grdim==1))
 if (test_nhat) then
   ABI_ALLOCATE(rhor1_ptr,(cplex*nfft,nspden))
   rhor1_ptr(:,:)=rhor1(:,:)-nhat1(:,:)
 else
   rhor1_ptr => rhor1
 end if

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1_ptr)

!Compute the gradients of the first-order density
!rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density
 ishift=0 ; ngrad=2
 ABI_ALLOCATE(rho1now,(cplex*nfft,nspden,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor1_ptr,rho1now)

!PAW: add "exact" gradients of compensation density
 if (test_nhat.and.usexcnhat==1) then
   rho1now(:,1:nspden,1)=rho1now(:,1:nspden,1)+nhat1(:,1:nspden)
 end if
 if (nhat1grdim==1) then
   do ii=1,3
     rho1now(:,1:nspden,ii+1)=rho1now(:,1:nspden,ii+1)+nhat1gr(:,1:nspden,ii)
   end do
 end if
 if (test_nhat) then
   ABI_DEALLOCATE(rhor1_ptr)
 end if

!Transfer the ground-state density and its gradient to spin-polarized storage
!rhortmp(:,:,1) contains the GS density, and
!rhortmp(:,:,2:4) contains the gradients of the GS density
 if(nspden==2)then
   ABI_ALLOCATE(rhortmp,(nfft,nspden,4))
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,15+2*ii)
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)-kxc(ir,15+2*ii)
     end do
   end do
 end if

!Apply the XC kernel
 nspgrad=2; if (nspden==2) nspgrad=5
 ABI_ALLOCATE(dnexcdn,(cplex*nfft,nspgrad))

 if (cplex==1) then  ! Treat real case first
   if (nspden==1) then
     do ir=1,nfft
       r0(:)=kxc(ir,5:7) ; r1(:)=rho1now(ir,1,2:4)
       gradrho_gradrho1=dot_product(r0,r1)
       dnexcdn(ir,1)=kxc(ir,1)*rho1now(ir,1,1) + kxc(ir,3)*gradrho_gradrho1
       coeff_grho=kxc(ir,3)*rho1now(ir,1,1) + kxc(ir,4)*gradrho_gradrho1
  !    Reuse the storage in rho1now
       rho1now(ir,1,2:4)=r1(:)*kxc(ir,2)+r0(:)*coeff_grho
     end do
   else
     do ir=1,nfft
       r0_up(:)=rhortmp(ir,1,2:4)   ! grad of spin-up GS rho
       r0_dn(:)=rhortmp(ir,2,2:4)   ! grad of spin-down GS rho
       r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
       r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
       r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
       r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
       gradrho_gradrho1_up=dot_product(r0_up,r1_up)
       gradrho_gradrho1_dn=dot_product(r0_dn,r1_dn)
       gradrho_gradrho1   =dot_product(r0,r1)
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
  &     +r0(:)*coeff_grho_corr
       rho1now(ir,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
  &     +r1_up(:)*kxc(ir,12)               &
  &     +r0_dn(:)*coeff_grho_dn            &
  &     +r0(:)*coeff_grho_corr
     end do
   end if ! nspden

 else ! if cplex==2
   if (nspden==1) then
     do ir=1,nfft
       r0(:)=kxc(ir,5:7)
       r1(:)  =rho1now(2*ir-1,1,2:4)
       r1im(:)=rho1now(2*ir  ,1,2:4)
       gradrho_gradrho1  =dot_product(r0,r1)
       gradrho_gradrho1im=dot_product(r0,r1im)
       dnexcdn(2*ir-1,1)=kxc(ir,1)*rho1now(2*ir-1,1,1) + kxc(ir,3)*gradrho_gradrho1
       dnexcdn(2*ir  ,1)=kxc(ir,1)*rho1now(2*ir  ,1,1) + kxc(ir,3)*gradrho_gradrho1im
       coeff_grho  =kxc(ir,3)*rho1now(2*ir-1,1,1) + kxc(ir,4)*gradrho_gradrho1
       coeffim_grho=kxc(ir,3)*rho1now(2*ir  ,1,1) + kxc(ir,4)*gradrho_gradrho1im
  !    Reuse the storage in rho1now
       rho1now(2*ir-1,1,2:4)=r1(:)  *kxc(ir,2)+r0(:)*coeff_grho
       rho1now(2*ir  ,1,2:4)=r1im(:)*kxc(ir,2)+r0(:)*coeffim_grho
     end do
   else
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
       gradrho_gradrho1_up  =dot_product(r0_up,r1_up)
       gradrho_gradrho1_dn  =dot_product(r0_dn,r1_dn)
       gradrho_gradrho1     =dot_product(r0,r1)
       gradrho_gradrho1im_up=dot_product(r0_up,r1im_up)
       gradrho_gradrho1im_dn=dot_product(r0_dn,r1im_dn)
       gradrho_gradrho1im   =dot_product(r0,r1im)
       dnexcdn(2*ir-1,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir-1,1,1)+&
&       kxc(ir,10)*rho1now(2*ir-1,2,1)+&
&       kxc(ir,5)*gradrho_gradrho1_up+&
&       kxc(ir,13)*gradrho_gradrho1
       dnexcdn(2*ir  ,1)=(kxc(ir,1)+kxc(ir,9))*rho1now(2*ir,1,1)+&
&       kxc(ir,10)*rho1now(2*ir,2,1)+&
&       kxc(ir,5)*gradrho_gradrho1im_up+&
&       kxc(ir,13)*gradrho_gradrho1im
       dnexcdn(2*ir-1,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir-1,2,1)+&
&       kxc(ir,10)*rho1now(2*ir-1,1,1)+&
&       kxc(ir,6)*gradrho_gradrho1_dn+&
&       kxc(ir,14)*gradrho_gradrho1
       dnexcdn(2*ir  ,2)=(kxc(ir,2)+kxc(ir,11))*rho1now(2*ir,2,1)+&
&       kxc(ir,10)*rho1now(2*ir,1,1)+&
&       kxc(ir,6)*gradrho_gradrho1im_dn+&
&       kxc(ir,14)*gradrho_gradrho1im
       coeff_grho_corr=(kxc(ir,13)*rho1now(2*ir-1,1,1)+kxc(ir,14)*rho1now(2*ir-1,2,1))+&
&       kxc(ir,15)*gradrho_gradrho1
       coeffim_grho_corr=(kxc(ir,13)*rho1now(2*ir,1,1)+kxc(ir,14)*rho1now(2*ir,2,1))+&
&       kxc(ir,15)*gradrho_gradrho1im
       coeff_grho_up= kxc(ir,5)*rho1now(2*ir-1,1,1)+kxc(ir,7)*gradrho_gradrho1_up
       coeffim_grho_up= kxc(ir,5)*rho1now(2*ir,1,1)+kxc(ir,7)*gradrho_gradrho1im_up
       coeff_grho_dn= kxc(ir,6)*rho1now(2*ir-1,2,1)+kxc(ir,8)*gradrho_gradrho1_dn
       coeffim_grho_dn= kxc(ir,6)*rho1now(2*ir,2,1)+kxc(ir,8)*gradrho_gradrho1im_dn
!      Reuse the storage in rho1now
       rho1now(2*ir-1,1,2:4)=r1_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&       +r1_dn(:)*kxc(ir,12)               &
&       +r0_up(:)*coeff_grho_up            &
&       +(r0_up(:)+r0_dn(:))*coeff_grho_corr
       rho1now(2*ir  ,1,2:4)=r1im_up(:)*(kxc(ir,3)+kxc(ir,12))   &
&       +r1im_dn(:)*kxc(ir,12)               &
&       +r0_up(:)*coeffim_grho_up            &
&       +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
       rho1now(2*ir-1,2,2:4)=r1_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&       +r1_up(:)*kxc(ir,12)               &
&       +r0_dn(:)*coeff_grho_dn            &
&       +(r0_up(:)+r0_dn(:))*coeff_grho_corr
       rho1now(2*ir  ,2,2:4)=r1im_dn(:)*(kxc(ir,4)+kxc(ir,12))   &
&       +r1im_up(:)*kxc(ir,12)               &
&       +r0_dn(:)*coeffim_grho_dn            &
&       +(r0_up(:)+r0_dn(:))*coeffim_grho_corr
     end do
   end if ! nspden

 end if

 vxc1(:,:)=zero
 call xcpot(cplex,dnexcdn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden,&
&           nspgrad,paral_kgb,qphon,rho1now,vxc1)

!call filterpot(paral_kgb,cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 ABI_DEALLOCATE(dnexcdn)
 if (nspden==2) then
   ABI_DEALLOCATE(rhortmp)
 end if
 ABI_DEALLOCATE(rho1now)

 DBG_EXIT("COLL")

end subroutine dfpt_mkvxcgga
!!***

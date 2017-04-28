!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkvxcstrgga
!! NAME
!! dfpt_mkvxcstrgga
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! for the strain perturbation in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (DRH, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  dgprimdds(3,3)=strain derivaive of gprimd.
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  istr=index of the strain perturbation (1..6)
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor1tmp(cplex*nfft,2)=array for first-order electron spin-density
!!   in electrons/bohr**3 (second index corresponds to spin-up and spin-down)
!!  str_scale=scaling factor for gradient operator strain-derivative (1. or 2.)
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential
!!
!! NOTES
!!  Closely related to dfpt_mkvxcgga.
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
!!      dfpt_mkvxcstr
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_mkvxcstrgga(cplex,dgprimdds,gprimd,istr,kxc,mpi_enreg,nfft,ngfft,&
& nkxc,nspden,paral_kgb,qphon,rhor1tmp,str_scale,vxc1)

 use defs_basis
 use defs_abitypes
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxcstrgga'
 use interfaces_56_xc
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istr,nfft,nkxc,nspden,paral_kgb
 real(dp) :: str_scale
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: dgprimdds(3,3),gprimd(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: qphon(3),rhor1tmp(cplex*nfft,2)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)
!Local variables-------------------------------
!scalars
 integer :: ii,ir,ishift,ispden,mgga,ngrad,nspgrad
 real(dp) :: coeff_grho,coeff_grho_corr,coeff_grho_dn,coeff_grho_up
 real(dp) :: gradrho_gradrho1,gradrho_gradrho1_dn,gradrho_gradrho1_up
 character(len=500) :: msg
!arrays
 real(dp) :: r0(3),r0_dn(3),r0_up(3),r1(3),r1_dn(3),r1_up(3)
 real(dp),allocatable :: dnexcdn(:,:),rho1now(:,:,:),rhodgnow(:,:,:)
 real(dp),allocatable :: rhordgtmp(:,:),rhortmp(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 if (nspden>2) then
   msg='Not compatible with non-collinear magnetism!'
   MSG_ERROR(msg)
 end if

!metaGGA contributions are not taken into account here
 mgga=0

!if you uncomment the following line, you will have to modify
!the original function call to pass in gmet and gsqcut
!call filterpot(cplex,gmet,gsqcut,nfft,ngfft,2,qphon,rhor1tmp)

!Compute the gradients of the first-order density
!rho1now(:,:,1) contains the first-order density, and
!rho1now(:,:,2:4) contains the gradients of the first-order density
 ishift=0 ; ngrad=2
 ABI_ALLOCATE(rho1now,(cplex*nfft,nspden,ngrad*ngrad))
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,paral_kgb,qphon,rhor1tmp,rho1now)

!Transfer the ground-state density and its gradient to spin-polarized storage
!rhortmp(:,:,1) contains the GS density, and
!rhortmp(:,:,2:4) contains the gradients of the GS density
 ABI_ALLOCATE(rhortmp,(nfft,nspden,4))
 if(nspden==2)then
   do ii=1,4
     do ir=1,nfft
       rhortmp(ir,1,ii)=kxc(ir,15+2*ii)
       rhortmp(ir,2,ii)=kxc(ir,14+2*ii)-kxc(ir,15+2*ii)
     end do
   end do
 end if

!Calculate the 1st-order contribution to grad(n) from the strain derivative
!  acting on the gradient operator acting on the GS charge density,
!Simply use the following formula:
!   (dGprim/ds_alpha_beta)(i,j) = -half.( delta_alpha,i Gprim(beta,j) + delta_beta,i Gprim(alpha,j) )
!To finally get:
!   (nabla)^(alpha,beta)_i[n] = -half ( delta_alpha,i nabla_beta[n] + delta_beta,i nabla_alpha[n] )
 ABI_ALLOCATE(rhodgnow,(cplex*nfft,nspden,ngrad*ngrad))
 rhodgnow(1:nfft,1:nspden,1:4)=zero
 if (nspden==1) then
   if (istr==1) rhodgnow(1:nfft,1,2)=-     kxc(1:nfft,5)
   if (istr==2) rhodgnow(1:nfft,1,3)=-     kxc(1:nfft,6)
   if (istr==3) rhodgnow(1:nfft,1,4)=-     kxc(1:nfft,7)
   if (istr==4) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,7)
   if (istr==4) rhodgnow(1:nfft,1,4)=-half*kxc(1:nfft,6)
   if (istr==5) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,7)
   if (istr==5) rhodgnow(1:nfft,1,4)=-half*kxc(1:nfft,5)
   if (istr==6) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,6)
   if (istr==6) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,5)
 else
   if (istr==1) rhodgnow(1:nfft,1,2)=-     kxc(1:nfft,19)
   if (istr==2) rhodgnow(1:nfft,1,3)=-     kxc(1:nfft,21)
   if (istr==3) rhodgnow(1:nfft,1,4)=-     kxc(1:nfft,23)
   if (istr==4) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,23)
   if (istr==4) rhodgnow(1:nfft,1,4)=-half*kxc(1:nfft,21)
   if (istr==5) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,23)
   if (istr==5) rhodgnow(1:nfft,1,4)=-half*kxc(1:nfft,19)
   if (istr==6) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,21)
   if (istr==6) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,19)
   if (istr==1) rhodgnow(1:nfft,2,2)=-     (kxc(1:nfft,18)-kxc(1:nfft,19))
   if (istr==2) rhodgnow(1:nfft,2,3)=-     (kxc(1:nfft,20)-kxc(1:nfft,21))
   if (istr==3) rhodgnow(1:nfft,2,4)=-     (kxc(1:nfft,22)-kxc(1:nfft,23))
   if (istr==4) rhodgnow(1:nfft,2,3)=-half*(kxc(1:nfft,22)-kxc(1:nfft,23))
   if (istr==4) rhodgnow(1:nfft,2,4)=-half*(kxc(1:nfft,20)-kxc(1:nfft,21))
   if (istr==5) rhodgnow(1:nfft,2,2)=-half*(kxc(1:nfft,22)-kxc(1:nfft,23))
   if (istr==5) rhodgnow(1:nfft,2,4)=-half*(kxc(1:nfft,18)-kxc(1:nfft,19))
   if (istr==6) rhodgnow(1:nfft,2,2)=-half*(kxc(1:nfft,20)-kxc(1:nfft,21))
   if (istr==6) rhodgnow(1:nfft,2,3)=-half*(kxc(1:nfft,18)-kxc(1:nfft,19))
 end if

!Add to the gradients of the first-order density
 do ii=2,4
   do ispden=1,nspden
     rhodgnow(1:nfft,ispden,ii)=str_scale*rhodgnow(1:nfft,ispden,ii)
     rho1now(1:nfft,ispden,ii)=rho1now(1:nfft,ispden,ii)+rhodgnow(1:nfft,ispden,ii)
   end do
 end do

!rho1now(:,:,1) contains the 1st-order density, and rho1now(:,:,2:4) contains the grads of the 1st-order density

!Apply the XC kernel
 nspgrad=2; if (nspden==2) nspgrad=5
 ABI_ALLOCATE(dnexcdn,(cplex*nfft,nspgrad))

!== Non polarized
 if (nspden==1) then
   do ir=1,nfft
     r0(:)=kxc(ir,5:7) ; r1(:)=rho1now(ir,1,2:4)
     gradrho_gradrho1=dot_product(r0,r1)
     dnexcdn(ir,1)=kxc(ir,1)*rho1now(ir,1,1) + kxc(ir,3)*gradrho_gradrho1
     coeff_grho=kxc(ir,3)*rho1now(ir,1,1) + kxc(ir,4)*gradrho_gradrho1
!    Grad strain derivative contribution enters the following term with a
!    factor of two compared to above terms, so add it again.
     r1(:)=r1(:)+rhodgnow(ir,1,2:4)
!    Reuse the storage in rho1now
     rho1now(ir,1,2:4)=r1(:)*kxc(ir,2)+r0(:)*coeff_grho
   end do

!== Spin-polarized
 else ! nspden==2
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

!    grad strain derivative contribution enters the following term with a
!    factor of two compared to above terms, so add it again.
     r1_up(:)=r1_up(:)+rhodgnow(ir,1,2:4)
     r1_dn(:)=r1_dn(:)+rhodgnow(ir,2,2:4)

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
 ABI_DEALLOCATE(rhodgnow)

 vxc1(:,:)=zero
 call xcpot(cplex,dnexcdn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden,&
&           nspgrad,paral_kgb,qphon,rho1now,vxc1)

!if you uncomment the following line, you will have to modify
!the original function call to pass in gmet and gsqcut
!call filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 ABI_DEALLOCATE(dnexcdn)
 if (nspden==2) then
   ABI_DEALLOCATE(rhortmp)
 end if
 ABI_DEALLOCATE(rho1now)

end subroutine dfpt_mkvxcstrgga
!!***

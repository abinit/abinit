!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_mkvxcstr
!! NAME
!!  m_dfpt_mkvxcstr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2019 ABINIT group (DRH,XG)
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

module m_dfpt_mkvxcstr

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_abicore

 use m_time,      only : timab
 use m_symtk,     only : matr3inv
 use m_xctk,      only : xcden, xcpot

 implicit none

 private
!!***

 public :: dfpt_mkvxcstr
!!***

contains
!!***

!!****f* ABINIT/dfpt_mkvxcstr
!! NAME
!! dfpt_mkvxcstr
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to strain: assemble the first-order density change with the
!! frozen-core density change, then use the exchange-correlation kernel.
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhotoxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhat1(cplex*nfft,2nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nkxc=second dimension of the kxc array
!!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  nspden=number of spin-density components
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used, otherwise, nfft
!!  option=if 0, work only with strain-derivative frozen-wavefunction
!!    charge and the XC core-correction,
!!   if 1, treat both density change and XC core correction
!!   if 2, like 0 but multiply gradient strain derivative term by 2.0 for GGA.
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!
!! OUTPUT
!!  vxc1(cplex*nfft,nspden)=change in exchange-correlation potential (including
!!   core-correction, if applicable)
!!
!! PARENTS
!!      dfpt_eltfrxc,dfpt_nselt,dfpt_nstpaw,dfpt_rhotov
!!
!! CHILDREN
!!      dfpt_mkvxcstrgga,matr3inv,timab
!!
!! SOURCE

subroutine dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,&
&                        nkxc,non_magnetic_xc,nspden,n3xccc,option,qphon,&
&                        rhor,rhor1,rprimd,usepaw,usexcnhat,vxc1,xccc3d1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n3xccc,natom,nfft,nkxc,nspden,option
 integer,intent(in) :: usepaw,usexcnhat
 logical,intent(in) :: non_magnetic_xc
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),target,intent(in) :: nhat(nfft,nspden)
 real(dp),target,intent(in) :: nhat1(cplex*nfft,nspden)
 real(dp),intent(in) :: kxc(nfft,nkxc),qphon(3)
 real(dp),target,intent(in) :: rhor(nfft,nspden),rhor1(cplex*nfft,nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xccc3d1(cplex*n3xccc)
 real(dp),intent(out) :: vxc1(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: ii,ir,istr
 real(dp) :: rho1_dn,rho1_up,spin_scale,str_scale
 character(len=500) :: message
!arrays
 real(dp) :: gprimd(3,3),tsec(2)
 real(dp),allocatable :: rhor1tmp(:,:),rhowk1(:,:)
 real(dp),pointer :: rhor_(:,:),rhor1_(:,:)

! *************************************************************************

 call timab(181,1,tsec)

 if(nspden/=1 .and. nspden/=2) then
   message = ' dfpt_mkvxc, Only for nspden==1 and 2.'
   MSG_BUG(message)
 end if

 if (usepaw==1.and.usexcnhat==0) then
   ABI_ALLOCATE(rhor_,(nfft,nspden))
   rhor_(:,:)=rhor(:,:)-nhat(:,:)
 else
   rhor_ => rhor
 end if

 if (usepaw==1.and.usexcnhat==0.and.option==1) then
   ABI_ALLOCATE(rhor1_,(nfft,nspden))
   rhor1_(:,:)=rhor1(:,:)-nhat1(:,:)
 else
   rhor1_ => rhor1
 end if

!Inhomogeneous term for diagonal strain
 ABI_ALLOCATE(rhowk1,(nfft,nspden))
 if(option==0 .or. option==2) then
   if(ipert==natom+3) then
     rhowk1(:,:)=-rhor_(:,:)
   else
     rhowk1(:,:)=zero
   end if
 else if(option==1) then
   if(ipert==natom+3) then
     rhowk1(:,:)=rhor1_(:,:)-rhor_(:,:)
   else
     rhowk1(:,:)=rhor1_(:,:)
   end if
 end if

 if (non_magnetic_xc) then
   if(nspden==2) rhowk1(:,2)=rhowk1(:,1)*half
   if(nspden==4) rhowk1(:,2:4)=zero
 end if

!Treat first LDA
 if(nkxc==1.or.nkxc==3)then

!  Case without non-linear core correction
   if(n3xccc==0)then

!    Non-spin-polarized
     if(nspden==1)then
       do ir=1,nfft
         vxc1(ir,1)=kxc(ir,1)*rhowk1(ir,1)
       end do

!      Spin-polarized
     else
       do ir=1,nfft
         rho1_dn=rhowk1(ir,1)-rhowk1(ir,2)
         vxc1(ir,1)=kxc(ir,1)*rhowk1(ir,2)+kxc(ir,2)*rho1_dn
         vxc1(ir,2)=kxc(ir,2)*rhowk1(ir,2)+kxc(ir,3)*rho1_dn
       end do
     end if ! nspden==1

!    Treat case with non-linear core correction
   else
     if(nspden==1)then
       do ir=1,nfft
         vxc1(ir,1)=kxc(ir,1)*(rhowk1(ir,1)+xccc3d1(ir))
       end do
     else
       do ir=1,nfft
         rho1_dn=rhowk1(ir,1)-rhowk1(ir,2) + xccc3d1(ir)*half
         rho1_up=rhowk1(ir,2)              + xccc3d1(ir)*half
         vxc1(ir,1)=kxc(ir,1)*rho1_up+kxc(ir,2)*rho1_dn
         vxc1(ir,2)=kxc(ir,2)*rho1_up+kxc(ir,3)*rho1_dn
       end do
     end if ! nspden==1

   end if ! n3xccc==0

!  Treat GGA
 else if (nkxc==7.or.nkxc==19) then

!  Generates gprimd and its strain derivative
!  Note that unlike the implicitly symmetric metric tensor strain
!  derivatives, we must explicltly symmetrize the strain derivative
!  here.
   call matr3inv(rprimd,gprimd)
   istr=idir + 3*(ipert-natom-3)
   if(istr<1 .or. istr>6)then
     write(message, '(a,i10,a,a,a)' )&
&     'Input dir gives istr=',istr,' not allowed.',ch10,&
&     'Possible values are 1,2,3,4,5,6 only.'
     MSG_BUG(message)
   end if

!  Rescalling needed for use in dfpt_eltfrxc for elastic tensor (not internal strain tensor).
   str_scale=one;if(option==2) str_scale=two

!  Transfer the data to spin-polarized storage
   ABI_ALLOCATE(rhor1tmp,(cplex*nfft,nspden))
   if(nspden==1)then
     do ir=1,cplex*nfft
       rhor1tmp(ir,1)=rhowk1(ir,1)
     end do
   else
     do ir=1,cplex*nfft
       rho1_dn=rhowk1(ir,1)-rhowk1(ir,2)
       rhor1tmp(ir,1)=rhowk1(ir,2)
       rhor1tmp(ir,2)=rho1_dn
     end do
   end if ! nspden==1
   if(n3xccc/=0)then
     spin_scale=one;if (nspden==2) spin_scale=half
     do ii=1,nspden
       do ir=1,cplex*nfft
         rhor1tmp(ir,ii)=rhor1tmp(ir,ii)+xccc3d1(ir)*spin_scale
       end do
     end do
   end if

   call dfpt_mkvxcstrgga(cplex,gprimd,istr,kxc,mpi_enreg,nfft,ngfft,nkxc,&
&   nspden,qphon,rhor1tmp,str_scale,vxc1)
   ABI_DEALLOCATE(rhor1tmp)

 else
   MSG_BUG('Invalid nkxc!')

 end if ! LDA or GGA

 if (usepaw==1.and.usexcnhat==0) then
   ABI_DEALLOCATE(rhor_)
 end if
 if (usepaw==1.and.usexcnhat==0.and.option==1) then
   ABI_DEALLOCATE(rhor1_)
 end if

 ABI_DEALLOCATE(rhowk1)

 call timab(181,2,tsec)

end subroutine dfpt_mkvxcstr
!!***

!!****f* ABINIT/dfpt_mkvxcstrgga
!! NAME
!! dfpt_mkvxcstrgga
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! for the strain perturbation in case of GGA functionals
!! Use the exchange-correlation kernel.
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  gprimd(3,3)=dimensional primitive translations in reciprocal space (bohr^-1)
!!  istr=index of the strain perturbation (1..6)
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhotoxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
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
!!       kxc(:,1)= d2Exc/drho_up drho_up
!!       kxc(:,2)= d2Exc/drho_up drho_dn
!!       kxc(:,3)= d2Exc/drho_dn drho_dn
!!       kxc(:,4)= 1/|grad(rho_up)| dEx/d|grad(rho_up)|
!!       kxc(:,5)= 1/|grad(rho_dn)| dEx/d|grad(rho_dn)|
!!       kxc(:,6)= 1/|grad(rho_up)| d2Ex/d|grad(rho_up)| drho_up
!!       kxc(:,7)= 1/|grad(rho_dn)| d2Ex/d|grad(rho_dn)| drho_dn
!!       kxc(:,8)= 1/|grad(rho_up)| * d/d|grad(rho_up)| ( 1/|grad(rho_up)| dEx/d|grad(rho_up)| )
!!       kxc(:,9)= 1/|grad(rho_dn)| * d/d|grad(rho_dn)| ( 1/|grad(rho_dn)| dEx/d|grad(rho_dn)| )
!!       kxc(:,10)=1/|grad(rho)| dEc/d|grad(rho)|
!!       kxc(:,11)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_up
!!       kxc(:,12)=1/|grad(rho)| d2Ec/d|grad(rho)| drho_dn
!!       kxc(:,13)=1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dEc/d|grad(rho)| )
!!       kxc(:,14)=gradx(rho_up)
!!       kxc(:,15)=gradx(rho_dn)
!!       kxc(:,16)=grady(rho_up)
!!       kxc(:,17)=grady(rho_dn)
!!       kxc(:,18)=gradz(rho_up)
!!       kxc(:,19)=gradz(rho_dn)
!!
!! PARENTS
!!      dfpt_mkvxcstr
!!
!! CHILDREN
!!      xcden,xcpot
!!
!! SOURCE

subroutine dfpt_mkvxcstrgga(cplex,gprimd,istr,kxc,mpi_enreg,nfft,ngfft,&
& nkxc,nspden,qphon,rhor1tmp,str_scale,vxc1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,istr,nfft,nkxc,nspden
 real(dp) :: str_scale
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: gprimd(3,3),kxc(nfft,nkxc)
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

! *************************************************************************

 DBG_ENTER("COLL")

 if (nkxc/=12*min(nspden,2)-5) then
   msg='Wrong nkxc value for GGA!'
   MSG_BUG(msg)
 end if
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
 call xcden(cplex,gprimd,ishift,mpi_enreg,nfft,ngfft,ngrad,nspden,qphon,rhor1tmp,rho1now)

!Calculate the 1st-order contribution to grad(n) from the strain derivative
!  acting on the gradient operator acting on the GS charge density,
!Simply use the following formula:
!   (dGprim/ds_alpha_beta)(i,j) = -half.( delta_alpha,i Gprim(beta,j) + delta_beta,i Gprim(alpha,j) )
!To finally get:
!   (nabla)^(alpha,beta)_i[n] = -half ( delta_alpha,i nabla_beta[n] + delta_beta,i nabla_alpha[n] )
 ABI_ALLOCATE(rhodgnow,(cplex*nfft,nspden,3))
 rhodgnow(1:nfft,1:nspden,1:3)=zero
 if (nspden==1) then
   if (istr==1) rhodgnow(1:nfft,1,1)=-     kxc(1:nfft,5)
   if (istr==2) rhodgnow(1:nfft,1,2)=-     kxc(1:nfft,6)
   if (istr==3) rhodgnow(1:nfft,1,3)=-     kxc(1:nfft,7)
   if (istr==4) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,7)
   if (istr==4) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,6)
   if (istr==5) rhodgnow(1:nfft,1,1)=-half*kxc(1:nfft,7)
   if (istr==5) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,5)
   if (istr==6) rhodgnow(1:nfft,1,1)=-half*kxc(1:nfft,6)
   if (istr==6) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,5)
 else
   if (istr==1) rhodgnow(1:nfft,1,1)=-     kxc(1:nfft,15)
   if (istr==2) rhodgnow(1:nfft,1,2)=-     kxc(1:nfft,17)
   if (istr==3) rhodgnow(1:nfft,1,3)=-     kxc(1:nfft,19)
   if (istr==4) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,19)
   if (istr==4) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,17)
   if (istr==5) rhodgnow(1:nfft,1,1)=-half*kxc(1:nfft,19)
   if (istr==5) rhodgnow(1:nfft,1,3)=-half*kxc(1:nfft,15)
   if (istr==6) rhodgnow(1:nfft,1,1)=-half*kxc(1:nfft,17)
   if (istr==6) rhodgnow(1:nfft,1,2)=-half*kxc(1:nfft,15)
   if (istr==1) rhodgnow(1:nfft,2,1)=-     (kxc(1:nfft,14)-kxc(1:nfft,15))
   if (istr==2) rhodgnow(1:nfft,2,2)=-     (kxc(1:nfft,16)-kxc(1:nfft,17))
   if (istr==3) rhodgnow(1:nfft,2,3)=-     (kxc(1:nfft,18)-kxc(1:nfft,19))
   if (istr==4) rhodgnow(1:nfft,2,2)=-half*(kxc(1:nfft,18)-kxc(1:nfft,19))
   if (istr==4) rhodgnow(1:nfft,2,3)=-half*(kxc(1:nfft,16)-kxc(1:nfft,17))
   if (istr==5) rhodgnow(1:nfft,2,1)=-half*(kxc(1:nfft,18)-kxc(1:nfft,19))
   if (istr==5) rhodgnow(1:nfft,2,3)=-half*(kxc(1:nfft,14)-kxc(1:nfft,15))
   if (istr==6) rhodgnow(1:nfft,2,1)=-half*(kxc(1:nfft,16)-kxc(1:nfft,17))
   if (istr==6) rhodgnow(1:nfft,2,2)=-half*(kxc(1:nfft,14)-kxc(1:nfft,15))
 end if

!Add to the gradients of the first-order density
 do ii=1,3
   do ispden=1,nspden
     rhodgnow(1:nfft,ispden,ii)=str_scale*rhodgnow(1:nfft,ispden,ii)
     rho1now(1:nfft,ispden,1+ii)=rho1now(1:nfft,ispden,1+ii)+rhodgnow(1:nfft,ispden,ii)
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
     r1(:)=r1(:)+rhodgnow(ir,1,1:3)
!    Reuse the storage in rho1now
     rho1now(ir,1,2:4)=r1(:)*kxc(ir,2)+r0(:)*coeff_grho
   end do

!== Spin-polarized
 else ! nspden==2
   do ir=1,nfft
     do ii=1,3  ! grad of spin-up ans spin_dn GS rho
       r0_up(ii)=kxc(ir,13+2*ii);r0_dn(ii)=kxc(ir,12+2*ii)-kxc(ir,13+2*ii)
     end do
     r0(:)=r0_up(:)+r0_dn(:)      ! grad of GS rho
     r1_up(:)=rho1now(ir,1,2:4)   ! grad of spin-up rho1
     r1_dn(:)=rho1now(ir,2,2:4)   ! grad of spin-down rho1
     r1(:)=r1_up(:)+r1_dn(:)      ! grad of GS rho1
     gradrho_gradrho1_up=dot_product(r0_up,r1_up)
     gradrho_gradrho1_dn=dot_product(r0_dn,r1_dn)
     gradrho_gradrho1   =dot_product(r0,r1)

     dnexcdn(ir,1)=kxc(ir, 1)*rho1now(ir,1,1)     &
&     +kxc(ir, 2)*rho1now(ir,2,1)     &
&     +kxc(ir, 6)*gradrho_gradrho1_up &
&     +kxc(ir,11)*gradrho_gradrho1
     dnexcdn(ir,2)=kxc(ir, 3)*rho1now(ir,2,1)     &
&     +kxc(ir, 2)*rho1now(ir,1,1)     &
&     +kxc(ir, 7)*gradrho_gradrho1_dn &
&     +kxc(ir,12)*gradrho_gradrho1
     coeff_grho_corr=kxc(ir,11)*rho1now(ir,1,1) &
&     +kxc(ir,12)*rho1now(ir,2,1) &
&     +kxc(ir,13)*gradrho_gradrho1
     coeff_grho_up=kxc(ir,6)*rho1now(ir,1,1)+kxc(ir,8)*gradrho_gradrho1_up
     coeff_grho_dn=kxc(ir,7)*rho1now(ir,2,1)+kxc(ir,9)*gradrho_gradrho1_dn

!    Grad strain derivative contribution enters the following term with a
!    factor of two compared to above terms, so add it again.
     r1_up(:)=r1_up(:)+rhodgnow(ir,1,1:3)
     r1_dn(:)=r1_dn(:)+rhodgnow(ir,2,1:3)

!    Reuse the storage in rho1now
     rho1now(ir,1,2:4)=(kxc(ir,4)+kxc(ir,10))*r1_up(:) &
&     +kxc(ir,10)            *r1_dn(:) &
&     +coeff_grho_up         *r0_up(:) &
&     +coeff_grho_corr       *r0(:)
     rho1now(ir,2,2:4)=(kxc(ir,5)+kxc(ir,10))*r1_dn(:) &
&     +kxc(ir,10)            *r1_up(:) &
&     +coeff_grho_dn         *r0_dn(:) &
&     +coeff_grho_corr       *r0(:)
   end do

 end if ! nspden
 ABI_DEALLOCATE(rhodgnow)

 vxc1(:,:)=zero
 call xcpot(cplex,dnexcdn,gprimd,ishift,mgga,mpi_enreg,nfft,ngfft,ngrad,nspden,&
& nspgrad,qphon,rho1now,vxc1)

!if you uncomment the following line, you will have to modify
!the original function call to pass in gmet and gsqcut
!call filterpot(cplex,gmet,gsqcut,nfft,ngfft,nspden,qphon,vxc1)

 ABI_DEALLOCATE(dnexcdn)
 ABI_DEALLOCATE(rho1now)

end subroutine dfpt_mkvxcstrgga
!!***

end module m_dfpt_mkvxcstr
!!***

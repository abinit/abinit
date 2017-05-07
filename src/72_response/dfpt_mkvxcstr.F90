!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_mkvxcstr
!! NAME
!! dfpt_mkvxcstr
!!
!! FUNCTION
!! Compute the first-order change of exchange-correlation potential
!! due to strain: assemble the first-order density change with the
!! frozen-core density change, then use
!! the exchange-correlation kernel.
!!
!! COPYRIGHT
!! Copyright (C) 2001-2017 ABINIT group (DRH,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  cplex= if 1, real space 1-order functions on FFT grid are REAL,
!!     if 2, COMPLEX
!!  idir=direction of the current perturbation
!!  ipert=type of the perturbation
!!  kxc(nfft,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mpi_enreg=information about MPI parallelization
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/input_variables/vargs.htm#ngfft
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhat1(cplex*nfft,2nspden*usepaw)= -PAW only- 1st-order compensation density
!!  nkxc=second dimension of the kxc array
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine dfpt_mkvxcstr(cplex,idir,ipert,kxc,mpi_enreg,natom,nfft,ngfft,nhat,nhat1,&
& nkxc,nspden,n3xccc,option,paral_kgb,qphon,rhor,rhor1,rprimd,usepaw,usexcnhat,vxc1,xccc3d1)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_mkvxcstr'
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_72_response, except_this_one => dfpt_mkvxcstr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,n3xccc,natom,nfft,nkxc,nspden,option
 integer,intent(in) :: paral_kgb,usepaw,usexcnhat
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

!  FTransfer the data to spin-polarized storage
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
&                        nspden,paral_kgb,qphon,rhor1tmp,str_scale,vxc1)
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

!{\src2tex{textfont=tt}}
!!****f* ABINIT/spectral_function
!! NAME
!! spectral_function
!!
!! FUNCTION
!! Print the spectral function computed from Green function in real frequency
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cryst_struc <type(crystal_t)>=crystal structure data
!!  green  <type(green_type)>= green function data 
!!  hu  <type(hu_type)>= datatype of type hu
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  self <type(self_type)>= variables related to self-energy
!!  prtopt= option for printing 
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve
!!
!! CHILDREN
!!      compute_green,copy_green,copy_matlu,dc_self,destroy_green,destroy_self
!!      dyson,hubbard_one,init_green,initialize_self,ldau_self,print_green
!!      rw_self,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine spectral_function(cryst_struc,green,hu,paw_dmft,&
& pawang,pawtab,self_old,prtopt)

 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_crystal, only : crystal_t
 use m_green, only : init_green,green_type,print_green,copy_green,compute_green,destroy_green
 use m_matlu, only : copy_matlu
 use m_paw_dmft, only : paw_dmft_type
 use m_hu, only : hu_type
 use m_self, only : self_type,initialize_self,dc_self,destroy_self,rw_self
 use m_energy, only : energy_type
 use m_pawang, only : pawang_type
 use m_pawtab, only : pawtab_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'spectral_function'
 use interfaces_14_hidewrite
 use interfaces_68_dmft, except_this_one => spectral_function
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
! type(pawang_type), intent(in) :: pawang
 type(crystal_t),intent(in) :: cryst_struc
 type(green_type), intent(in) :: green
 type(hu_type),intent(inout) :: hu(cryst_struc%ntypat)
 !type(MPI_type), intent(inout) :: mpi_enreg
 type(pawang_type), intent(in) :: pawang
 type(pawtab_type),intent(in)  :: pawtab(cryst_struc%ntypat)
 type(self_type), intent(inout) :: self_old
 type(paw_dmft_type), intent(inout)  :: paw_dmft
 integer, intent(in) :: prtopt

!Local variables ------------------------------
 character(len=500) :: message
 type(green_type) :: greenr
 type(green_type) :: weissr
 type(self_type) :: selfr
!scalars
!************************************************************************
!character(len=500) :: message

!   opt_oper_ksloc=3 to be able to compute spectral function.
 call init_green(greenr,paw_dmft,opt_oper_ksloc=3,wtype="real")
 call init_green(weissr,paw_dmft,wtype="real")
 call copy_matlu(green%occup%matlu,greenr%occup%matlu,paw_dmft%natom)
 call initialize_self(selfr,paw_dmft,wtype="real")
!=======================================================================
!== Solve impurity model with green function for real frequency   
!=======================================================================
 write(message,'(2a,i3,13x,a)') ch10,'  ===  Write Spectral function'
 call wrtout(std_out,message,'COLL')
 if(abs(paw_dmft%dmft_solv)==1) then

!  == LDA+U for test
!  -------------------
   call ldau_self(cryst_struc,greenr,paw_dmft,&
&   pawtab,selfr,opt_ldau=1,prtopt=prtopt)
 else if(abs(paw_dmft%dmft_solv)==2) then

!  == Hubbard One
!  -------------------
   call hubbard_one(cryst_struc,greenr,hu,paw_dmft,&
&   pawang,prtopt,self_old%hdc,weissr)

 else if(abs(paw_dmft%dmft_solv)==4) then

!  == Nothing
!  -------------------
   message = "spectral_function: This section of code is disabled!"
   MSG_ERROR(message)
   call copy_green(weissr,greenr,opt_tw=1)

 else if(abs(paw_dmft%dmft_solv)>=5) then

!  == Nothing
!  -------------------
   MSG_ERROR("Stopping before copy_green")
   call copy_green(weissr,greenr,opt_tw=1)

 else if(abs(paw_dmft%dmft_solv)==0) then

!  == Nothing
!  -------------------
!  weiss%occup%has_operks=0 -> only local part is duplicated
   call copy_green(weissr,greenr,opt_tw=2)
 end if

!=======================================================================
!== Integrate green function and printout occupations
!For dmft_solv=-1,0,or 1 , the green function was not computed: it
!cannot be integrated
!=======================================================================
 call dc_self(green%charge_matlu_solver,cryst_struc,hu,selfr,paw_dmft%dmft_dc,prtopt)
 if(abs(paw_dmft%dmft_solv)/=1.and.paw_dmft%dmft_solv/=0) then
   call dyson(greenr,paw_dmft,selfr,weissr,opt_weissself=2) 
 end if
 call compute_green(cryst_struc,greenr,paw_dmft,pawang,1,selfr,opt_self=1)
 call print_green("realw",greenr,4,paw_dmft,pawprtvol=3)
 call rw_self(selfr,paw_dmft,prtopt=2,opt_rw=2)

 if(abs(prtopt)>0) then
 end if
 call destroy_self(selfr)
 call destroy_green(weissr)
 call destroy_green(greenr)

end subroutine spectral_function
!!***

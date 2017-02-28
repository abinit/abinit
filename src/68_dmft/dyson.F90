!{\src2tex{textfont=tt}}
!!****f* ABINIT/dyson
!! NAME
!! dyson
!!
!! FUNCTION
!! Use the Dyson Equation to compute self-energy from green function 
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  istep    =  step of iteration for LDA.
!!  lda_occup
!!  mpi_enreg=informations about MPI parallelization
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  opt_weissgreen= 1: compute weiss from green and self
!!                = 2: compute green from weiss and self
!!                = 4: compute green from weiss and self without
!!                     inversion of weiss  (weiss is previously inverted)
!!
!! OUTPUT
!!  edmft  = energy in DMFT.
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      dmft_solve,spectral_function
!!
!! CHILDREN
!!      add_matlu,copy_green,destroy_green,init_green,inverse_oper,timab,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine dyson(green,paw_dmft,self,weiss,opt_weissself)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_paw_dmft, only: paw_dmft_type
 use m_crystal, only : crystal_t
 use m_green, only : green_type, destroy_green,init_green,copy_green
 use m_oper, only : oper_type,inverse_oper
 use m_matlu, only : matlu_type,add_matlu,print_matlu
 use m_self, only : self_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dyson'
 use interfaces_14_hidewrite
 use interfaces_18_timing
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(green_type),intent(inout)  :: green
 type(paw_dmft_type), intent(inout) :: paw_dmft
 type(self_type),intent(inout)  :: self
 type(green_type),intent(inout)  :: weiss
 integer,intent(in) :: opt_weissself
! type(paw_dmft_type), intent(inout)  :: paw_dmft

!Local variables ------------------------------
!arrays
 real(dp) :: tsec(2)
!scalars
 integer :: ifreq,natom,nspinor,nsppol,weissinv
 type(green_type)  :: greeninv
 character(len=500) :: message
! type
! type(matlu_type), pointer :: matlutemp,matlu1,matlu2
!************************************************************************

 call timab(623,1,tsec)
 DBG_ENTER("COLL")
 natom=green%oper(1)%natom
 nsppol=green%oper(1)%nsppol
 nspinor=green%oper(1)%nspinor
 weissinv=1
 if(opt_weissself==1) then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Use Dyson Equation => weiss '
   call wrtout(std_out,message,'COLL')
 else if(opt_weissself==2) then
   write(message,'(2a,i3,13x,a)') ch10,'  ===  Use Dyson Equation => self'
   call wrtout(std_out,message,'COLL')
 end if

!call xmpi_barrier(spaceComm)
 if(paw_dmft%dmft_solv==2) weissinv=0
 call init_green(greeninv,paw_dmft,opt_oper_ksloc=2,wtype=green%w_type)
 call copy_green(green,greeninv,opt_tw=2)

 do ifreq=1,green%nw
   if(opt_weissself==1) then
     call inverse_oper(greeninv%oper(ifreq),option=1,prtopt=1)
!    warning green is now inversed
     call add_matlu(greeninv%oper(ifreq)%matlu,self%oper(ifreq)%matlu,&
&     weiss%oper(ifreq)%matlu,natom,sign_matlu2=1)
     call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
   else if(opt_weissself==2) then

   ! write(59,*) paw_dmft%omega_lo(ifreq), real(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(green%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(61,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(60,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(62,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
!    call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
!    write(63,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))

!    write(std_out,*) "-----------------------IFREQ",ifreq
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
     call inverse_oper(greeninv%oper(ifreq),option=1,prtopt=1)
!    call print_matlu(greeninv%oper(ifreq)%matlu,paw_dmft%natom,1,opt_diag=-1)
!    If paw_dmft%dmft_solv==2, then inverse of weiss function is
!    computed in hubbard_one.F90
     if(weissinv/=0) then
       call inverse_oper(weiss%oper(ifreq),option=1,prtopt=1)
     end if

!    write(std_out,*) weiss%oper(1)%matlu(ifreq)%mat(1,1,1,1,1),"-",greeninv%oper(ifreq)
     call add_matlu(weiss%oper(ifreq)%matlu,greeninv%oper(ifreq)%matlu,&
&     self%oper(ifreq)%matlu,natom,sign_matlu2=-1)

   ! write(64,*) paw_dmft%omega_lo(ifreq), real(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(greeninv%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(65,*) paw_dmft%omega_lo(ifreq), real(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(weiss%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   ! write(66,*) paw_dmft%omega_lo(ifreq), real(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1)),aimag(self%oper(ifreq)%matlu(1)%mat(1,1,1,1,1))
   else
     message = " BUG in dyson.F90"
     MSG_BUG(message)
   end if
 end do

 call destroy_green(greeninv)


 call timab(623,2,tsec)
 DBG_EXIT("COLL")

end subroutine dyson
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/compute_levels
!! NAME
!! compute_levels
!!
!! FUNCTION
!! Compute levels for ctqmc
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (BAmadon)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  
!! 
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      hubbard_one,qmc_prep_ctqmc
!!
!! CHILDREN
!!      checkdiag_matlu,loc_oper,print_matlu,sym_matlu,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


 subroutine compute_levels(cryst_struc,energy_level,hdc,pawang,paw_dmft,nondiag)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_errors
 use m_profiling_abi

 use m_pawang, only : pawang_type
 use m_crystal, only : crystal_t
 use m_paw_dmft, only : paw_dmft_type
 use m_oper, only : oper_type,loc_oper
 use m_matlu, only : sym_matlu, print_matlu, checkdiag_matlu

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'compute_levels'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(pawang_type), intent(in) :: pawang
 type(oper_type), intent(in) :: hdc
 type(paw_dmft_type), intent(in)  :: paw_dmft
 type(oper_type),intent(inout)  :: energy_level !vz_i
 logical, optional, intent(out) :: nondiag

!Local variables ------------------------------
! scalars
 integer :: iatom,iband,ikpt,im1,isppol,ispinor
 integer :: lpawu,mbandc,natom,nspinor,nsppol,nkpt
 character(len=500) :: message
! arrays
!************************************************************************

 mbandc=paw_dmft%mbandc
 nkpt=paw_dmft%nkpt
 nsppol=paw_dmft%nsppol
 nspinor=paw_dmft%nspinor
 natom=paw_dmft%natom
 if(present(nondiag)) nondiag=.false.

!========================
!Get KS eigenvalues
!========================
 do iband=1,mbandc
   do ikpt=1,nkpt
     do isppol=1,nsppol
!      Take \epsilon_{nks}
!      ========================
       energy_level%ks(isppol,ikpt,iband,iband)=paw_dmft%eigen_lda(isppol,ikpt,iband)
     end do
   end do
 end do


!======================================================================
!Compute atomic levels from projection of \epsilon_{nks} and symetrize
!======================================================================
 call loc_oper(energy_level,paw_dmft,1)
! write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels before sym and only LDA"
! call wrtout(std_out,message,'COLL')
! call print_matlu(energy_level%matlu,natom,1)
 do iatom = 1 , natom
   lpawu=paw_dmft%lpawu(iatom)
   if(lpawu/=-1) then
     do isppol=1,nsppol
       do ispinor=1,nspinor
         do im1=1,2*lpawu+1
           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)= &
&           energy_level%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)&
&           -hdc%matlu(iatom)%mat(im1,im1,isppol,ispinor,ispinor)-paw_dmft%fermie 
         end do
       end do
     end do
!    write(std_out,*) "DC,fermie",hdc%matlu(iatom)%mat(1,1,1,1,1),paw_dmft%fermie
   end if
 end do ! natom
 call sym_matlu(cryst_struc,energy_level%matlu,pawang)
 if(present(nondiag)) call checkdiag_matlu(energy_level%matlu,natom,tol7,nondiag)

 write(message,'(a,2x,a,f13.5)') ch10," == Print Energy levels for Fermi Level=",paw_dmft%fermie
 call wrtout(std_out,message,'COLL')
!call print_oper(energy_level,1,paw_dmft,1)
 call print_matlu(energy_level%matlu,natom,1)



 end subroutine compute_levels
!!***

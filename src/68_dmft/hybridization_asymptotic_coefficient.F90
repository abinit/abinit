!{\src2tex{textfont=tt}}
!!****f* ABINIT/hybridization_asymptotic_coefficient
!! NAME
!! hybridization_asymptotic_coefficient
!!
!! FUNCTION
!! Solve the Impurity problem
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
!!  lda_occup
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!  pawang <type(pawang)>=paw angular mesh and related data
!!  pawtab <type(pawtab)>
!!
!! OUTPUT
!!  paw_dmft =  data for self-consistent LDA+DMFT calculations.
!!
!! NOTES
!!
!! PARENTS
!!      qmc_prep_ctqmc
!!
!! CHILDREN
!!      add_matlu,copy_oper,destroy_oper,init_oper,loc_oper,prod_oper,sym_matlu
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif


#include "abi_common.h"

subroutine hybridization_asymptotic_coefficient(cryst_struc,paw_dmft,pawang,hybri_coeff)


 use defs_basis
 use defs_abitypes
 use m_errors
 use m_profiling_abi
 use m_crystal, only : crystal_t
 use m_oper, only : oper_type,init_oper,upfold_oper,copy_oper,prod_oper,destroy_oper,loc_oper
 use m_matlu, only : matlu_type,init_matlu,add_matlu,destroy_matlu,print_matlu,sym_matlu
 use m_paw_dmft, only: paw_dmft_type
 use m_pawang, only : pawang_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'hybridization_asymptotic_coefficient'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(crystal_t),intent(in) :: cryst_struc
 type(paw_dmft_type), intent(in) :: paw_dmft
 type(pawang_type), intent(in) :: pawang
 type(matlu_type), intent(inout) :: hybri_coeff(paw_dmft%natom)
!Local variables ------------------------------
 type(oper_type)  :: ham_a
 type(oper_type)  :: ham_b
 type(oper_type)  :: ham_squarelocal
 type(oper_type)  :: ham_squareks
 integer :: iband1,iband2,ikpt,isppol
!************************************************************************

! call init_oper(paw_dmft,self_minus_hdc)
 call init_oper(paw_dmft,ham_a)
 call init_oper(paw_dmft,ham_b)
 call init_oper(paw_dmft,ham_squareks)
 call init_oper(paw_dmft,ham_squarelocal)

! Create self_minus_hdc%matlu = Sigma-Hdc in local orbitals
! call add_matlu(self%oper(paw_dmft%dmft_nwlo)%matlu,self%hdc%matlu,&
!&             self_minus_hdc%matlu,cryst_struc%natom,-1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (1)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

!! Create self_minus_hdc%ks = Upfold Sigma-Hdc
! call upfold_oper(self_minus_hdc,paw_dmft,1)
! call loc_oper(self_minus_hdc,paw_dmft,1)

!   write(message,'(a,2x,a)') ch10,        "  == self_minus_hdc (2)"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(self_minus_hdc%matlu,paw_dmft%natom,1,opt_exp=1)

! Create ham_a%ks = H_ks  in KS basis
!----------------------------------------------------
 do iband1=1,paw_dmft%mbandc
   do iband2=1,paw_dmft%mbandc
     do ikpt=1,paw_dmft%nkpt
       do isppol=1,paw_dmft%nsppol
         if(iband1==iband2) then
           ham_a%ks(isppol,ikpt,iband1,iband2) = cmplx(paw_dmft%eigen_lda(isppol,ikpt,iband1),0.d0,kind=dp)
         else 
           ham_a%ks(isppol,ikpt,iband1,iband2) = czero 
         end if
       end do
     end do
   end do
 end do

! Create ham_a%matlu = H_ks in local orbitals
!---------------------------------------------
 call loc_oper(ham_a,paw_dmft,1)

! Symetrise the local quantity (energy levels)
!---------------------------------------------
 call sym_matlu(cryst_struc,ham_a%matlu,pawang)

! Create ham_b%ks : Duplicate both ks and local part of ham_a into ham_b
!-----------------------------------------------------------------------
 call copy_oper(ham_a,ham_b)

! Compute ham_squareks%ks   : Make a product of the two KS version
!------------------------------------------------------------------
 call prod_oper(ham_a,ham_b,ham_squareks,1)

! Compute ham_squareks%matlu 
!---------------------------
 call loc_oper(ham_squareks,paw_dmft,1)

! Symetrise ham_squareks%matlu 
!------------------------------
 call sym_matlu(cryst_struc,ham_squareks%matlu,pawang)

!   write(message,'(a,2x,a)') ch10,        "  == squareks"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squareks%matlu,paw_dmft%natom,1,opt_exp=1)


! Compute ham_squarelocal%matlu  
!-------------------------------
 call prod_oper(ham_a,ham_b,ham_squarelocal,2)

! Compute the product in local orbitals
!--------------------------------------
 call sym_matlu(cryst_struc,ham_squarelocal%matlu,pawang)

!   write(message,'(a,2x,a)') ch10,        "  == squarelocal"
!   call wrtout(std_out,message,'COLL')
!   call print_matlu(ham_squarelocal%matlu,paw_dmft%natom,1,opt_exp=1)

! The difference of ham_squareks and ham_squarelocal
! gives the coefficient that we are looking for ( such that F_ij(iw_n) = -C_ij/(iw_n) ).
!----------------------------------------------------------------------------------------
 call add_matlu(ham_squareks%matlu,ham_squarelocal%matlu,hybri_coeff,cryst_struc%natom,-1)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij before sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

! Symetrise the local quantity 
!------------------------------
 call sym_matlu(cryst_struc,hybri_coeff,pawang)

 !  write(message,'(a,2x,a)') ch10,        "  == Coeff C_ij after sym"
 !  call wrtout(std_out,message,'COLL')
 !  call print_matlu(hybri_coeff,paw_dmft%natom,1,opt_exp=1)

 call destroy_oper(ham_squarelocal)
! call destroy_oper(self_minus_hdc)
 call destroy_oper(ham_a)
 call destroy_oper(ham_b)
 call destroy_oper(ham_squareks)


end subroutine hybridization_asymptotic_coefficient
!!***

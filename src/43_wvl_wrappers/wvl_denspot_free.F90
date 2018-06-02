!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_denspot_free
!! NAME
!!  wvl_denspot_free
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2018 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_denspot_distribution,deallocate_rho_descriptors
!!      denspot_free_history,f_free_ptr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_denspot_free(den)
    
 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
#if defined HAVE_BIGDFT
 use BigDFT_API, only: deallocate_rho_descriptors, &
      & deallocate_denspot_distribution, denspot_free_history
 use dynamic_memory
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_denspot_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvl_denspot_type), intent(inout) :: den

!Local variables-------------------------------
 
! *************************************************************************
 
!DEBUG
!write (std_out,*) ' wvl_denspot_free : enter'
!ENDDEBUG

#if defined HAVE_BIGDFT
 if(associated(den%denspot%rhov)) then
   call f_free_ptr(den%denspot%rhov)
 end if
 if(associated(den%denspot%rho_psi)) then
   call f_free_ptr(den%denspot%rho_psi)
 end if
 if(associated(den%denspot%rho_C)) then
   call f_free_ptr(den%denspot%rho_C)
 end if
 if(associated(den%denspot%V_ext)) then
   call f_free_ptr(den%denspot%V_ext)
 end if
 if(associated(den%denspot%V_XC)) then
   call f_free_ptr(den%denspot%V_XC)
 end if
 if(associated(den%denspot%Vloc_KS)) then
   call f_free_ptr(den%denspot%Vloc_KS)
 end if
 if(associated(den%denspot%f_XC)) then
   call f_free_ptr(den%denspot%f_XC)
 end if
 if(associated(den%denspot%rho_work)) then
   call f_free_ptr(den%denspot%rho_work)
 end if
 if(associated(den%denspot%pot_work)) then
   call f_free_ptr(den%denspot%pot_work)
 end if
 nullify(den%denspot%rhov)
 nullify(den%denspot%rho_psi)
 nullify(den%denspot%rho_C)
 nullify(den%denspot%V_ext)
 nullify(den%denspot%V_XC)
 nullify(den%denspot%Vloc_KS)
 nullify(den%denspot%f_XC)
 nullify(den%denspot%rho_work)
 nullify(den%denspot%pot_work)
 !
 call deallocate_rho_descriptors(den%denspot%rhod)
 call deallocate_denspot_distribution(den%denspot%dpbox)
 call denspot_free_history(den%denspot)
#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) den%symObj
#endif

end subroutine wvl_denspot_free
!!***

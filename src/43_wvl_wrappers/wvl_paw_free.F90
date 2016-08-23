!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_paw_free
!! NAME
!!  wvl_paw_free
!!
!! FUNCTION
!!  Frees memory for WVL+PAW implementation
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2016 ABINIT group (T. Rangel)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  ntypat = number of atom types
!!  wvl= wvl type
!!  wvl_proj= wvl projector type
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_paw_free(wvl)
    
 use defs_basis
 use defs_wvltypes
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_paw_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(wvl_internal_type),intent(inout) :: wvl
  
!Local variables-------------------------------
 
! *************************************************************************
 
#if defined HAVE_DFT_BIGDFT

!PAW objects
 if( associated(wvl%paw%spsi)) then
   ABI_DEALLOCATE(wvl%paw%spsi)
 end if
 if( associated(wvl%paw%indlmn)) then
   ABI_DEALLOCATE(wvl%paw%indlmn)
 end if
 if( associated(wvl%paw%sij)) then
   ABI_DEALLOCATE(wvl%paw%sij)
 end if
 if( associated(wvl%paw%rpaw)) then
   ABI_DEALLOCATE(wvl%paw%rpaw)
 end if

!rholoc
 if( associated(wvl%rholoc%msz )) then
   ABI_DEALLOCATE(wvl%rholoc%msz)
 end if
 if( associated(wvl%rholoc%d )) then
   ABI_DEALLOCATE(wvl%rholoc%d)
 end if
 if( associated(wvl%rholoc%rad)) then
   ABI_DEALLOCATE(wvl%rholoc%rad)
 end if
 if( associated(wvl%rholoc%radius)) then
   ABI_DEALLOCATE(wvl%rholoc%radius)
 end if

#else
 if (.false.) write(std_out,*) wvl%h(1)
#endif

!paw%paw_ij and paw%cprj are allocated and deallocated inside vtorho

end subroutine wvl_paw_free
!!***

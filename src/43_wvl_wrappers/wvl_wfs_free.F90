!{\src2tex{textfont=tt}}
!!****f* ABINIT/wvl_wfs_free
!!
!! NAME
!! wvl_wfs_free
!!
!! FUNCTION
!! Freeing routine.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  wfs <type(wvl_wf_type)>=wavefunctions informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      deallocate_comms,deallocate_lr,deallocate_lzd_except_glr
!!      deallocate_orbs,f_free_ptr
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_wfs_free(wfs)

 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
#if defined HAVE_BIGDFT
 use BigDFT_API, only: deallocate_Lzd_except_Glr, deallocate_lr, &
      & deallocate_orbs, deallocate_comms
 use dynamic_memory
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_wfs_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wvl_wf_type),intent(inout) :: wfs

!Local variables -------------------------

! *********************************************************************

#if defined HAVE_BIGDFT
 call deallocate_Lzd_except_Glr(wfs%ks%lzd)
 call deallocate_lr(wfs%ks%lzd%Glr)
 call deallocate_orbs(wfs%ks%orbs)
 call deallocate_comms(wfs%ks%comms)
 if (associated(wfs%ks%orbs%eval))  then
   ABI_DEALLOCATE(wfs%ks%orbs%eval)
 end if
 ABI_DATATYPE_DEALLOCATE(wfs%ks%confdatarr)

 if (associated(wfs%ks%psi)) then
   call f_free_ptr(wfs%ks%psi)
 end if
 if (associated(wfs%ks%hpsi)) then
   call f_free_ptr(wfs%ks%hpsi)
 end if
 if (associated(wfs%ks%psit)) then
   call f_free_ptr(wfs%ks%psit)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) wfs%ks
#endif

end subroutine wvl_wfs_free
!!***

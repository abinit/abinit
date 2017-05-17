
!!****f* ABINIT/wvl_projectors_free
!!
!! NAME
!! wvl_projectors_free
!!
!! FUNCTION
!! Freeing routine.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DC)
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
!!  proj <type(wvl_projectors_type)>=projectors informations in a wavelet basis.
!!
!! PARENTS
!!      gstate,wvl_wfsinp_reformat
!!
!! CHILDREN
!!      free_dft_psp_projectors
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_projectors_free(proj)

 use defs_basis
 use defs_wvltypes
 use m_profiling_abi
 use m_errors
#if defined HAVE_BIGDFT
  use BigDFT_API, only : free_DFT_PSP_projectors,deallocate_gwf
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_projectors_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(wvl_projectors_type),intent(inout) :: proj

!Local variables -------------------------
#if defined HAVE_BIGDFT
 integer :: ii
#endif

  ! *********************************************************************

#if defined HAVE_BIGDFT

 call free_DFT_PSP_projectors(proj%nlpsp)

 if (allocated(proj%G)) then
   do ii=1,size(proj%G)
!    MT dec 2014: cannot call bigdft deallocation routine
!    because content of proj%G datastructure was created
!    without f_malloc (without memory profiling).
!    call deallocate_gwf(proj%G(ii))
     if (associated(proj%G(ii)%ndoc)) then
       ABI_DEALLOCATE(proj%G(ii)%ndoc)
     end if
     if (associated(proj%G(ii)%nam)) then
       ABI_DEALLOCATE(proj%G(ii)%nam)
     end if
     if (associated(proj%G(ii)%nshell)) then
       ABI_DEALLOCATE(proj%G(ii)%nshell)
     end if
     if (associated(proj%G(ii)%psiat)) then
       ABI_DEALLOCATE(proj%G(ii)%psiat)
     end if
     if (associated(proj%G(ii)%xp)) then
       ABI_DEALLOCATE(proj%G(ii)%xp)
     end if
   end do
   ABI_DATATYPE_DEALLOCATE(proj%G)
 end if

#else
 BIGDFT_NOTENABLED_ERROR()
 if (.false.) write(std_out,*) proj%nlpsp
#endif

end subroutine wvl_projectors_free
!!***

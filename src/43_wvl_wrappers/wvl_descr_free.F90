!!****f* defs_wvltypes/wvl_descr_free
!!
!! NAME
!! wvl_descr_free
!!
!! FUNCTION
!! Free the wvl%atoms% datastructure (deallocate or nullify)
!!
!! INPUTS
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! OUTPUT
!! wvl <type(wvl_internal_type)>=internal variables for wavelets
!!
!! PARENTS
!!      gstate,wvl_memory
!!
!! CHILDREN
!!      deallocate_atoms_data,f_free_ptr
!!
!! SOURCE
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine wvl_descr_free(wvl)

 use m_profiling_abi
 use defs_wvltypes
#if defined HAVE_BIGDFT
 use BigDFT_API, only : deallocate_atoms_data
 use dynamic_memory
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'wvl_descr_free'
!End of the abilint section

  implicit none

!Arguments ------------------------------------
!scalars
  type(wvl_internal_type), intent(inout) :: wvl
!arrays

!Local variables-------------------------------
!scalars

! *********************************************************************

#if defined HAVE_BIGDFT
!These arrays are pointers on memory handled by ABINIT.
 nullify(wvl%atoms%astruct%sym%irrzon)
 nullify(wvl%atoms%astruct%sym%phnons)
 if (associated(wvl%atoms%nlccpar)) then
   call f_free_ptr(wvl%atoms%nlccpar)
 end if
 call deallocate_atoms_data(wvl%atoms)
#endif
 if(allocated(wvl%npspcode_paw_init_guess)) then
   ABI_DEALLOCATE(wvl%npspcode_paw_init_guess)
 end if
end subroutine wvl_descr_free
!!***

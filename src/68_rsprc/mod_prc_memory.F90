!!****m* ABINIT/mod_prc_memory
!! NAME
!! mod_prc_memory
!!
!! FUNCTION
!! This modules defines arrays and data used for the real-space kerker
!! preconditionning of potential residuals.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (PMA).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!  FIXME: this is highly non-kosher. Should be a datastructure which is declared dynamically
!!  MG: I completely agree! We don't use modules to share data and I don't see why we should
!!  break the rule here.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module mod_prc_memory

 use defs_basis
 use m_abicore

 implicit none

private

  real(dp),public,save,allocatable :: rdiemac(:)

  integer,save,public :: cycle=0 ! This is great! A global variable with the same name as a Fortran Statement!
  real(dp),save,public :: energy_min

public :: prc_mem_init
public :: prc_mem_free

! *********************************************************************

 contains
!!***

!!****f* ABINIT/prc_mem_init
!! NAME
!! prc_mem_init
!!
!! FUNCTION
!! This subroutine allocates the module's main component
!!
!! PARENTS
!!      m_prcref
!!
!! CHILDREN
!!
!! SOURCE

subroutine prc_mem_init(nfft)

implicit none

!Arguments -------------------------------
integer, intent(in) :: nfft
!Local variables -------------------------
! *********************************************************************

   if (.not. allocated(rdiemac))  then
     ABI_MALLOC(rdiemac,(nfft))
   end if
   if(nfft.ne.size(rdiemac)) then ! This steps should be done over "istep" instead
     ABI_FREE(rdiemac)
     ABI_MALLOC(rdiemac,(nfft))
     cycle=0
   end if

 end subroutine prc_mem_init
!!***

!!****f* ABINIT/prc_mem_free
!! NAME
!! prc_mem_free
!!
!! FUNCTION
!! This subroutine deallocates the module's main component
!!
!! PARENTS
!!      m_scfcv_core
!!
!! CHILDREN
!!
!! SOURCE

subroutine prc_mem_free()

implicit none

! *********************************************************************

   if (allocated(rdiemac))  then
     ABI_FREE(rdiemac)
   end if

 end subroutine prc_mem_free
!!***

end module mod_prc_memory
!!***

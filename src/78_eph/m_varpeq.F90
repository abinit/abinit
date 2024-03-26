!!****m* ABINIT/m_varpeq
!! NAME
!!  m_varpeq
!!
!! FUNCTION
!!  Description
!!
!! COPYRIGHT
!!  Copyright (C) 2023-2024 ABINIT group (VV)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_varpeq

 use defs_basis
 use m_abicore
 use m_dtset
 use m_crystal
 use m_ebands
 use m_errors
 use m_ifc

 use m_fstrings,            only : sjoin, itoa

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_varpeq/varpeq_t
!! NAME
!!  varpeq_t
!!
!! FUNCTION
!!  Description
!!
!! SOURCE

 type, public :: varpeq_t

  contains

    procedure :: init => varpeq_init

 end type varpeq_t
!!***


contains !=====================================================================
!!***

!!****f* m_varpeq/varpeq_init
!! NAME
!!  varpeq_init
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine varpeq_init(self, dtset)

!Arguments ------------------------------------
!scalars
 class(dataset_type), intent(in) :: dtset
 class(varpeq_t), intent(inout) :: self
!arrays

!Local variables-------------------------------
!scalars

!arrays

! *************************************************************************

end subroutine varpeq_init
!!***


end module m_varpeq
!!***

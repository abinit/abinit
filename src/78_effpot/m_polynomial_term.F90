!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_term
!!
!! NAME
!! m_polynomial_term
!!
!! FUNCTION
!! COPYRIGHT
!! Copyright (C) 2010-2015 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_polynomial_term

 use defs_basis
 use m_errors
 use m_profiling_abi

 implicit none

 public :: polynomial_term_init
 public :: polynomial_term_free
!!***

!!****t* m_polynomial_term/polynomial_term_type
!! NAME
!! polynomial_term_type
!!
!! FUNCTION
!! structure for specific displacements of term
!!
!! SOURCE

 type, public :: polynomial_term_type

   integer :: atindx(2)
!     atindx(2)
!     Indexes of the atoms a and b in the unit cell

   integer :: cell(3,2)
!     cell(3,2)
!     indexes of the cell of the atom a and b

   integer :: ndisp
!     Number of displacement for this terms
!     1 for (X_y-O_y)^3, 2 for (X_y-O_y)(X_x-O_y)^2...

   integer :: power
!      power of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1

   real(dp) :: weight
!     weight of the term


 end type polynomial_term_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_polynomial_term/polynomial_term_init
!!
!! NAME
!! polynomial_term_init
!!
!! FUNCTION
!! Initialize polynomial_term_init for given set of displacements
!!
!! INPUTS
!! atindx(2) = Indexes of the atoms a and b in the unit cell
!! cell(3,2) = Indexes of the cell of the atom a and b
!! ndisp     = Number of displacement for this terms
!! power     = Power of the displacement 2 (X_z-O_z)^2 or 1 for (X_y-O_y)^1
!! weight    = Weight of the term
!!
!! OUTPUT
!! polynomial_term = polynomial_term structure to be initialized
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine polynomial_term_init(atindx,cell,ndisp,polynomial_term,power,weight)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_term_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: ndisp,power
 real(dp),intent(in) :: weight
!arrays
 integer, intent(in) :: atindx(2)
 integer, intent(in) :: cell(3,2)
 type(polynomial_term_type), intent(out) :: polynomial_term
!Local variables-------------------------------
!scalar
 integer :: ii
!arrays
 character(len=500) :: msg

! *************************************************************************

end subroutine polynomial_term_init
!!***


!!****f* m_polynomial_term/polynomial_term_free
!!
!! NAME
!! polynomial_term_free
!!
!! FUNCTION
!! Free polynomial_term
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_term = polynomial_term structure to be free
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine polynomial_term_free(polynomial_term)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'polynomial_term_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(polynomial_term_type), intent(inout) :: polynomial_term
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

 polynomial_term%atindx(:) = zero
 polynomial_term%cell(:,:) = zero
 polynomial_term%ndisp     = zero
 polynomial_term%power     = zero
 polynomial_term%weight    = zero

end subroutine polynomial_term_free
!!***

end module m_polynomial_term
!!***

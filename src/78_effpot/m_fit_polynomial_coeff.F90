!{\src2tex{textfont=tt}}
!!****f* ABINIT/m_polynomial_coeff
!!
!! NAME
!! m_fit_polynomial_coeff
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

module m_fit_polynomial_coeff

 use defs_basis
 use m_errors
 use m_profiling_abi
 use m_polynomial_coeff
 use m_xmpi

 implicit none

 public :: fit_polynomial_coeff_free
 public :: fit_polynomial_coeff_init

!!****t* m_fit_polynomial_coeff/fit_polynomial_coeff_type
!! NAME
!! fit_polynomial_coeff_type
!!
!! FUNCTION
!! structure for a polynomial coefficient
!! contains the value of the coefficient and a 
!! list of terms (displacement) relating to the coefficient
!!
!! SOURCE

 type, public :: fit_polynomial_coeff_type

   character(200) :: name
!     Name of the polynomial_coeff (Sr_y-O1_y)^3) for example

   type(polynomial_coeff_type),dimension(:),allocatable :: coefficients
!     terms(nterm)
!     contains all the displacements for this coefficient      

 end type fit_polynomial_coeff_type
!!***

CONTAINS  !===========================================================================================


!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_init
!!
!! NAME
!! fit_polynomial_coeff_init
!!
!! FUNCTION
!! Initialize scell structure, from unit cell vectors, qpoint chosen, and atoms
!!
!! INPUTS
!!  name     = Name of the polynomial_coeff (Sr_y-O1_y)^3) for example
!!  nterm   = Number of terms (short range interaction) for this polynomial_coeff
!!  coefficient  = Value of the coefficient of this term
!!  termstype(nterm) = Polynomial_term_type contains all the displacements for this coefficient 
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be initialized
!!
!! PARENTS
!!
!!
!! CHILDREN
!!
!!
!! SOURCE

subroutine fit_polynomial_coeff_init(coefficient,name,nterm,polynomial_coeff,terms)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_init'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: nterm
 real(dp),intent(in) :: coefficient
!arrays
 character(200),intent(in) :: name
 type(polynomial_term_type),intent(in) :: terms(nterm)
 type(fit_polynomial_coeff_type), intent(out) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************
 
!First free before initilisation
 call fit_polynomial_coeff_free(polynomial_coeff)

!Initilisation

end subroutine fit_polynomial_coeff_init
!!***

!!****f* m_fit_polynomial_coeff/fit_polynomial_coeff_free
!!
!! NAME
!! fit_polynomial_coeff_free
!!
!! FUNCTION
!! Free polynomial_coeff
!!
!! INPUTS
!!
!! OUTPUT
!! polynomial_coeff = polynomial_coeff structure to be free
!!
!! PARENTS
!!
!!
!! CHILDREN
!!   
!!
!! SOURCE

subroutine fit_polynomial_coeff_free(polynomial_coeff)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fit_polynomial_coeff_free'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
!arrays
 type(fit_polynomial_coeff_type), intent(inout) :: polynomial_coeff
!Local variables-------------------------------
!scalar
!arrays

! *************************************************************************

end subroutine fit_polynomial_coeff_free
!!***

end module m_fit_polynomial_coeff
!!***
 

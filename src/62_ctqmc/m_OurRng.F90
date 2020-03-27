
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!!****m* ABINIT/m_OurRng
!! NAME
!!  m_OurRng
!!
!! FUNCTION
!!  Random number generator module
!!  Should be modify and merge with uniformrandom and zbq
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#include "defs.h"

MODULE m_OurRng
!! Implementation of various RNG with a small footprint

 !use m_numeric_tools,  only : uniformrandom

IMPLICIT NONE

PRIVATE

PUBLIC :: OurRng

CONTAINS
!!***

!!****f* ABINIT/m_OurRng/OurRng
!! NAME
!!  OurRng
!!
!! FUNCTION
!!  Generator given by G. Colin de Verdiere
!!  Efficient on GPU and MIC
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2020 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  xn=seed
!!  rng=random number
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE OurRng(xn,rng)
  ! returns a value between 0. and 1. with a period of 2**31
  ! implements the Marsaglia serie:
  !   xn+1 = (69069 * xn) mod 2^31
!Arguments ------------------------------------
  DOUBLE PRECISION, INTENT(  OUT) :: rng
  INTEGER(8), INTENT(INOUT) :: xn
  !
  INTEGER(8) :: two31  ! 2 ** 31
  INTEGER(8) :: two31m ! 2 ** 31 -1
  INTEGER(8), PARAMETER :: mars   = 69069
  INTEGER(8) :: xn8
  INTRINSIC MOD, REAL, IAND

  two31 = 1
  two31 = two31 * 65536   ! **16
  two31 = two31 * 32768   ! **31
  two31m = two31 - 1

!!$  two31  = z'80000000'
!!$  two31m = z'7FFFFFFF'

  IF (xn == 0) xn = 1
  xn8 = (mars * xn)
  xn8 = IAND(xn8, two31m)
  xn = xn8

  rng = REAL(xn, 8) / REAL(two31m, 8)
  ! guard to avoid pick up one since that sould never happen (otherwise ctqmc
  ! may generate an error and exit the code)
  if ( rng == 1.d0 ) rng = 0.d0
END SUBROUTINE OurRng
!!***

END MODULE m_OurRng
!!***

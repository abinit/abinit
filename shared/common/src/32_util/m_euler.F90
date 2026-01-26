!!****m* ABINIT/m_euler
!! NAME
!!  m_euler
!!
!! FUNCTION
!!  This subroutine calculate the euler angles
!!
!! COPYRIGHT
!! Copyright (C) 2002-2025 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_euler

 use defs_basis
 use m_abicore

 implicit none

 private

 public :: geteuler  ! exponential of a complex matrix


CONTAINS  !===========================================================
 !!***

 !!****f* m_euler/geteuler
 !! NAME
 !! geteuler
 !!
 !! FUNCTION
 !! Compute the Euler angles (alpha, beta) corresponding to the spin quantization axis given in Cartesian coordinates.
 !!
 !! INPUTS
 !! spinaxis(3)=spin quantization axis
 !!
 !! OUTPUT
 !! alpha=Euler angle for rotation around z-axis
 !! beta =Euler angle for rotation around y-axis
 !!
 !! SOURCE

 subroutine geteuler(spinaxis, alpha, beta)

 !Arguments -------------------------------
 !scalars
  real(dp),intent(out) :: alpha, beta
 !arrays
  real(dp),intent(in) :: spinaxis(3)

 !Local variables -------------------------
 !scalars
  real(dp) :: sx, sy, sz, norm, rxy
 !***********************************************************************
  norm = DOT_PRODUCT(spinaxis(:), spinaxis(:))
  if (norm < tol8*tol8) then
    alpha = zero; beta = zero
    return
  end if

  sx = spinaxis(1); sy = spinaxis(2); sz = spinaxis(3)
  rxy = sqrt(sx*sx + sy*sy)
  if (rxy < tol8) then
    alpha = zero
  else
    alpha = atan2(sy, sx)
  end if

  beta  = atan2(rxy, sz)

 end subroutine geteuler
 !!***

END MODULE m_euler
!!***

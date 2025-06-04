!!****m* ABINIT/m_brentq
!! NAME
!! m_brentq
!!
!! FUNCTION
!! This module contains Brent's root-finding method.
!! This was translated from scipy.
!! Originally written by Charles Harris charles.harris@sdl.usu.edu
!!
!! COPYRIGHT
!! Copyright (c) 2001-2002 Enthought, Inc. 2003-2024, SciPy Developers.
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions
!! are met:
!!
!! 1. Redistributions of source code must retain the above copyright
!!    notice, this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above
!!    copyright notice, this list of conditions and the following
!!    disclaimer in the documentation and/or other materials provided
!!    with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived
!!    from this software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
!! "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
!! LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
!! A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
!! OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
!! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
!! LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
!! DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
!! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
!! (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
!! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!!
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_brentq

 use defs_basis
 use m_errors

 implicit none

 private

 public :: brentq

CONTAINS  !========================================================================================
!!***

!!****f* m_brentq/brentq
!! NAME
!! brentq
!!
!! FUNCTION
!!
!! Brent's root finding method (taken from scipy).
!!
!! INPUTS
!! f= subroutine corresponding to the function for which the optimization is to be done, of the form subroutine f(x,fx)
!! xa,xb= interval on which to perform the optimization
!! xtol,rtol=tolerance criterion, the output root x0 will satisfy the criterion
!!           abs(x-x0) < xtol + rtol * x0 with x the true root, default value on scipy are
!!           2e-12 and 4*machine_precision respectively
!! iter= maximum number of iterations, default value on scipy is 100
!!
!! OUTPUT
!! xcur= root
!! ierr= 1 if the algorithm converged, 0 otherwise
!!
!! SOURCE

subroutine brentq(f,xa,xb,xtol,rtol,iter,xcur,ierr)

! Written by Charles Harris charles.harris@sdl.usu.edu

!Arguments ------------------------------------
  interface
    subroutine f(x,fx)
      use defs_basis
      real(dp), intent(in) :: x
      real(dp), intent(out) :: fx
    end subroutine f
  end interface
  real(dp), intent(in) :: xa,xb,xtol,rtol
  integer, intent(in) :: iter
  real(dp), intent(out) :: xcur
  integer, intent(out) :: ierr
!Local variables ------------------------------
  integer :: i
  real(dp) :: dblk,delta,dpre,fblk,fcur,fpre,sbis,scur,spre,stry,xblk,xpre
  character(len=200) :: msg
!************************************************************************

  xpre = xa
  xcur = xb
  xblk = zero
  fblk = zero
  spre = zero
  scur = zero
  ierr = 0

! the tolerance is 2*delta
  call f(xpre,fpre)
  call f(xcur,fcur)

  if (fpre == zero) then
    ierr = 1
    xcur = xpre
    return
  end if

  if (fcur == zero) then
    ierr = 1
    return
  end if

  if (sign(one,fpre) == sign(one,fcur)) then
    msg = 'Sign error in brentq'
    ABI_ERROR(msg)
  end if

  do i=1,iter

    if ((fpre /= zero) .and. (fcur /= zero) .and. (sign(one,fpre) /= sign(one,fcur))) then
      xblk = xpre
      fblk = fpre
      spre = xcur - xpre
      scur = xcur - xpre
    end if

    if (abs(fblk) < abs(fcur)) then
      xpre = xcur
      xcur = xblk
      xblk = xpre
      fpre = fcur
      fcur = fblk
      fblk = fpre
    end if

    delta = (xtol+rtol*abs(xcur)) * half
    sbis = (xblk - xcur) * half
    if ((fcur == zero) .or. (abs(sbis) < delta)) then
      ierr = 1
      return
    end if

    if ((abs(spre) > delta) .and. (abs(fcur) < abs(fpre))) then
      if (xpre == xblk) then
        ! interpolate
        stry = -fcur * (xcur-xpre) / (fcur-fpre)
      else
        ! extrapolate
        dpre = (fpre-fcur) / (xpre-xcur)
        dblk = (fblk-fcur) / (xblk-xcur)
        stry = -fcur * (fblk*dblk-fpre*dpre) / (dblk*dpre*(fblk-fpre))
      end if

      if (two*abs(stry) < min(abs(spre),three*abs(sbis)-delta)) then
        ! good short step
        spre = scur
        scur = stry
      else
        ! bisect
        spre = sbis
        scur = sbis
      end if
    else
      ! bisect
      spre = sbis
      scur = sbis
    end if

    xpre = xcur
    fpre = fcur

    if (abs(scur) > delta) then
      xcur = xcur + scur
    else
      if (sbis > zero) then
        xcur = xcur + delta
      else
        xcur = xcur - delta
      end if
    end if

    call f(xcur,fcur)

  end do ! i

end subroutine brentq
!!***

!----------------------------------------------------------------------

END MODULE m_brentq
!!***


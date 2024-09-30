!!****m* ABINIT/m_brentq
!! NAME
!! m_brentq
!!
!! FUNCTION
!! This module contains Brent's root-finding method (taken from scipy).
!!
!! COPYRIGHT
!! Copyright (C) 2023 
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

 use m_errors
 use defs_basis

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
!! f= subroutine for which the optimization is to be done, of the form subroutine f(x,fx)
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

subroutine brentq(f, xa, xb, xtol, rtol, iter, xcur, ierr)

 implicit none

! Written by Charles Harris charles.harris@sdl.usu.edu 

!Arguments ------------------------------------
!scalars
  real(dp), intent(in)  :: xa, xb, xtol, rtol
  integer, intent(in)   :: iter
  real(dp), intent(out) :: xcur
  external              :: f
  integer, intent(out)  :: ierr

!Local variables ------------------------------
  real(dp) :: xpre, xblk, fpre, fcur, fblk, spre, scur, sbis, delta, stry, dpre, dblk
  integer  :: i
  character(len=200) :: msg

  xpre = xa
  xcur = xb
  xblk = zero
  fpre = zero 
  fcur = zero
  fblk = zero 
  spre = zero 
  scur = zero 
  ierr = 0
 
! the tolerance is 2*delta 
  call f(xpre, fpre)
  call f(xcur, fcur)
 
  if (fpre == zero) then
    xcur = xpre
    return 
  end if

  if (fcur == zero) then
    return
  end if

  if (sign(one,fpre)==sign(one,fcur)) then
    msg='Sign error in brentq'
    ABI_ERROR(msg)
  end if

  do i=1,iter 
    if ((fpre .ne. zero) .and. (fcur .ne. zero) .and. (sign(one,fpre) .ne. sign(one,fcur))) then
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

    delta = (xtol + rtol*abs(xcur))/two
    sbis = (xblk - xcur)/two;
    if ((fcur == zero) .or. (abs(sbis) < delta)) then
      ierr = 1 
      return 
    end if

    if ((abs(spre) > delta) .and. (abs(fcur) < abs(fpre))) then
      if (xpre == xblk) then
        ! interpolate 
        stry = -fcur*(xcur - xpre)/(fcur - fpre)
      else 
        !extrapolate 
        dpre = (fpre - fcur)/(xpre - xcur)
        dblk = (fblk - fcur)/(xblk - xcur)
        stry = -fcur*(fblk*dblk - fpre*dpre) /(dblk*dpre*(fblk - fpre))
      end if

      if (2*abs(stry) < MIN(abs(spre), three*abs(sbis) - delta)) then
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

    call f(xcur, fcur) 

  end do

end subroutine brentq
!!***

!----------------------------------------------------------------------

END MODULE m_brentq
!!***


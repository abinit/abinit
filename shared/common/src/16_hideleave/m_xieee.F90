!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xieee
!! NAME
!!  m_xieee
!!
!! FUNCTION
!!   Debugging tools and helper functions providing access to IEEE exceptions
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!   See F2003 standard and http://www.nag.com/nagware/np/r51_doc/ieee_exceptions.html
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xieee

#ifdef HAVE_FC_IEEE_EXCEPTIONS
 !use, intrinsic :: ieee_exceptions
 use ieee_exceptions
#endif

 implicit none

 private

 public :: xieee_halt_ifexc       ! Halt the code if one of the *usual* IEEE exceptions is raised.
 public :: xieee_signal_ifexc     ! Signal if any IEEE exception is raised.

 integer,private,parameter :: std_out = 6

contains
!!***

!!****f* m_xieee/xieee_halt_ifexc
!! NAME
!!  xieee_halt_ifexc
!!
!! FUNCTION
!!  Halt the code if one of the *usual* IEEE exceptions is raised.
!!
!! INPUTS
!!  halt= If the value is true, the exceptions will cause halting; otherwise, execution will continue after this exception.
!!
!! PARENTS
!!      m_argparse
!!
!! CHILDREN
!!      ieee_set_flag
!!
!! SOURCE

subroutine xieee_halt_ifexc(halt)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: halt
! *************************************************************************

#ifdef HAVE_FC_IEEE_EXCEPTIONS
 ! Possible Flags: ieee_invalid, ieee_overflow, ieee_divide_by_zero, ieee_inexact, and ieee_underflow
 if (ieee_support_halting(ieee_invalid)) then
   call ieee_set_halting_mode(ieee_invalid, halt)
 end if
 if (ieee_support_halting(ieee_overflow)) then
   call ieee_set_halting_mode(ieee_overflow, halt)
 end if
 if (ieee_support_halting(ieee_divide_by_zero)) then
   call ieee_set_halting_mode(ieee_divide_by_zero, halt)
 end if
 !if (ieee_support_halting(ieee_inexact)) then
 !  call ieee_set_halting_mode(ieee_inexact, halt)
 !end if
 !if (ieee_support_halting(ieee_underflow)) then
 !  call ieee_set_halting_mode(ieee_underflow, halt)
 !end if
#else
 write(std_out,*)"Cannot set halting mode to: ",halt
#endif

end subroutine xieee_halt_ifexc
!!***

!----------------------------------------------------------------------

!!****f* m_xieee/xieee_signal_ifexc
!! NAME
!!  xieee_signal_ifexc
!!
!! FUNCTION
!!  Signal if one of the *usual* IEEE exceptions is raised.

!! INPUTS
!!  flag= If the value is true, the exceptions will be signalled
!!
!! PARENTS
!!      m_argparse
!!
!! CHILDREN
!!      ieee_set_flag
!!
!! SOURCE

subroutine xieee_signal_ifexc(flag)

!Arguments ------------------------------------
!scalars
 logical,intent(in) :: flag
! *************************************************************************

#ifdef HAVE_FC_IEEE_EXCEPTIONS
 ! Possible Flags: ieee_invalid, ieee_overflow, ieee_divide_by_zero, ieee_inexact, and ieee_underflow
 call ieee_set_flag(ieee_invalid, flag)
 call ieee_set_flag(ieee_overflow, flag)
 call ieee_set_flag(ieee_divide_by_zero, flag)
 call ieee_set_flag(ieee_inexact, flag)
 call ieee_set_flag(ieee_underflow, flag)
#else
 write(std_out,*)"Cannot set signal flag to: ",flag
#endif

end subroutine xieee_signal_ifexc
!!***

end module m_xieee
!!***

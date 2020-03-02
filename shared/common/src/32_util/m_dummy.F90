!!****f* DUMMY_TESTS/m_dummy_tests
!! NAME
!! m_dummy_tests
!!
!! FUNCTION
!!  Dummy module, to detect unused values
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  used_arg = used argument
!!
!! OUTPUT
!!  dummy_out1 = first output value
!!
!! SIDE EFFECTS
!!  used_arg = used argument
!!
!! PARENTS
!!      dummy_tests
!!
!! CHILDREN
!!
!! SOURCE

 module m_dummy_tests

 implicit none

 !private
 public

 double precision, save :: dummy_value
 public :: test_dummy

 contains
!!***

!!****f* DUMMY_TESTS/test_unused_arg
!! NAME
!! test_unused_arg
!!
!! FUNCTION
!!  Dummy subroutine, to test unused arguments
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  unused_arg = unused argument
!!
!! OUTPUT

!! SIDE EFFECTS
!!  used_arg = used argument
!!
!! PARENTS
!!      dummy_tests
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

 subroutine test_unused_arg(used_arg,unused_arg)

 implicit none

!Local variables-------------------------------
!scalars
 integer, intent(inout) :: used_arg
 integer, intent(in) :: unused_arg

!******************************************************************

 used_arg=used_arg+1

 end subroutine test_unused_arg
!!***

!!****f* DUMMY_TESTS/test_same_actual_arg
!! NAME
!! test_same_actual_arg
!!
!! FUNCTION
!!  Dummy subroutine, to detect when the calling subroutine attributes to one variable results from two different arguments
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  used_arg = used argument
!!
!! OUTPUT
!!  dummy_out1 = first output value
!!  dummy_out2 = second output value
!!
!! SIDE EFFECTS
!!  used_arg = used argument
!!
!! PARENTS
!!      dummy_tests
!!
!! CHILDREN
!!
!! SOURCE

 subroutine test_same_actual_arg(dummy_out1,dummy_out2,used_arg)

 implicit none

!Local variables-------------------------------
!scalars
 integer, intent(in) :: used_arg
 integer, intent(out) :: dummy_out1, dummy_out2

!******************************************************************

 dummy_out1=used_arg+1
 dummy_out2=used_arg+2

 end subroutine test_same_actual_arg
!!***

!!****f* DUMMY_TESTS/test_dummy
!! NAME
!! test_dummy
!!
!! FUNCTION
!!  Dummy subroutine
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2020 ABINIT group (XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  used_arg = used argument
!!
!! OUTPUT
!!  dummy_out1 = first output value
!!
!! PARENTS
!!      dummy_tests
!!
!! CHILDREN
!!
!! SOURCE

 subroutine test_dummy(dummy_out1,used_arg)

 implicit none

!Local variables-------------------------------
!scalars
 integer, intent(in) :: used_arg
 integer, intent(out) :: dummy_out1

!******************************************************************

 dummy_out1=used_arg

 end subroutine test_dummy

 end module m_dummy_tests
!!***

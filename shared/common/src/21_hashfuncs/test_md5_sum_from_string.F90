!{\src2tex{textfont=tt}}
!!****p* ABINIT/21_hashfuncs/tests/test_md5_sum_from_string
!! NAME
!!  test_md5_sum_from_string
!!
!! FUNCTION
!!  Unit test for the test_md5_sum_from_string routine of m_hash_md5.
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2020 ABINIT Group (Yann Pouillon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!!  This program is run by "make check".
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

program test_md5_sum_from_string

use defs_basis, only: ch10,std_out
use m_hash_md5
implicit none

character(len=*), parameter  :: abinit_input_string = "ABINIT" // ch10
character(len=32), parameter :: abinit_md5_ref = &
&                                 "1c57f60933157a3b6e13def40286f0c5"

character(len=32) :: chk_value
logical :: test_ok

#if defined DEBUG_MODE
write(std_out,'(a,a)') "Unitary test: md5_sum_from_string", ch10
write(std_out,'(1x,a," = ",a)') "Reference MD5 sum", abinit_md5_ref
#endif

chk_value = md5_sum_from_string(abinit_input_string)

#if defined DEBUG_MODE
write(std_out,'(1x,a," = ",a)') "Computed  MD5 sum", chk_value
#endif

! Report test result
test_ok = md5_check(abinit_md5_ref, chk_value)
#if defined DEBUG_MODE
write(std_out,*)
if ( test_ok ) then
  write(std_out,'(a,a)') "TEST OK", ch10
else
  write(std_out,'(a,a)') "TEST FAILED", ch10
end if
#else
if ( .not. test_ok ) then
  write(std_out,'(a,a)') "TEST FAILED", ch10
end if
#endif

end program test_md5_sum_from_string
!!***

!!****p* ABINIT/21_hashfuncs/tests/test_md5_sum_from_file
!! NAME
!!  test_md5_sum_from_file
!!
!! FUNCTION
!!  Unit test for the test_md5_sum_from_file routine of m_hash_md5.
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

program test_md5_sum_from_file

use defs_basis, only: ch10,fnlen,std_out
use m_hash_md5
implicit none

character(len=*), parameter  :: abinit_input_string = "ABINIT"
character(len=32), parameter :: abinit_md5_single_line_ref = &
&                                 "1c57f60933157a3b6e13def40286f0c5"
character(len=32), parameter :: abinit_md5_multi_line_ref = &
&                                 "953ac666dcdbf651ee3393d374216499"

character(len=fnlen), parameter :: tmp_file = "test_md5_sum_from_file.tmp"

character(len=32) :: chk_value_single, chk_value_multi
logical :: test_ok

#if defined DEBUG_MODE
write(std_out,'(a,a)') "Unitary test: md5_sum_from_file", ch10
#endif

! File of one line + newline character
#if defined DEBUG_MODE
write(std_out,'(1x,a,1x,a)') "Creating single-line ", trim(tmp_file)
write(std_out,'(1x,a," = ",a)') "Reference MD5 sum", abinit_md5_single_line_ref
#endif

open(unit=1, file=trim(tmp_file), form="formatted", status="replace")
write(1,'(a)') abinit_input_string
close(1)

chk_value_single = md5_sum_from_file(tmp_file)

#if defined DEBUG_MODE
write(std_out,'(1x,a," = ",a)') "Computed  MD5 sum", chk_value_single
#endif

! File of three lines + newline characters
#if defined DEBUG_MODE
write(std_out,*)
write(std_out,'(1x,a,1x,a)') "Creating multi-line ", trim(tmp_file)
write(std_out,'(1x,a," = ",a)') "Reference MD5 sum", abinit_md5_multi_line_ref
#endif

open(unit=1, file=trim(tmp_file), form="formatted", status="replace")
write(1,'(a)') abinit_input_string
write(1,'(a)') abinit_input_string
write(1,'(a)') abinit_input_string
close(1)

chk_value_multi = md5_sum_from_file(tmp_file)

#if defined DEBUG_MODE
write(std_out,'(1x,a," = ",a)') "Computed  MD5 sum", chk_value_multi
#endif

! Report test result
test_ok = ( md5_check(abinit_md5_single_line_ref, chk_value_single) .and. &
& md5_check(abinit_md5_multi_line_ref, chk_value_multi) )
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

end program test_md5_sum_from_file
!!***

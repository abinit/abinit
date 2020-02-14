!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_pair_list
!! NAME
!!  m_pair_list
!!
!! FUNCTION
!!  This module defines an API to build
!!  dictionaries containing string keys and numeric or string values.
!!  It is implemented in C as a simple linked pair list (associative list).
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (TC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! This module provide an implementation of a pair list
!! Possible improvement:
!! - Simplify the usage of get by removing the limit in key and string size
!! - Simplify the usage of get by removing the need for variable for all possible
!!   content when you know what is stored
!!
!! PARENTS
!!   m_yaml_out, m_neat, m_common
!!
!! CHILDREN
!!   m_type_pair_list
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_pair_list

  use iso_c_binding
  use m_type_pair_list
  use m_errors

  use m_fstrings, only : sjoin

  implicit none

  integer,parameter :: TC_EMPTY=-2, TC_NOTFOUND=-1, TC_INT=0, TC_REAL=1, TC_STRING=2

  private
  public :: pair_list_set, pair_list_get, pair_list_free
  public :: pair_list_next, pair_list_look, pair_list_iter, pair_list_restart
  public :: pair_list
  public :: TC_EMPTY, TC_NOTFOUND, TC_INT, TC_REAL, TC_STRING

  type :: pair_list
    type(c_pair_list) :: plc
    contains
      procedure :: set => pair_list_set
      procedure :: set_keys => pair_list_set_keys
      procedure :: set_keys_to_null => pair_list_set_keys_to_null
      procedure :: get => pair_list_get
      procedure :: free => pair_list_free
      procedure :: next => pair_list_next
      procedure :: look => pair_list_look
      procedure :: iter => pair_list_iter
      procedure :: restart => pair_list_restart
      procedure :: length => pair_list_length
  end type pair_list

! -------------------------------------------------------------------------------
! -                                                                             -
! -                        Private C function binding                           -
! -                                                                             -
! -------------------------------------------------------------------------------
  interface

    subroutine pair_list_next_c(pl) bind(C, name="pair_list_next")
      use m_type_pair_list
      type(c_pair_list),intent(in) :: pl
    end subroutine pair_list_next_c

    subroutine pair_list_free_c(pl) bind(C, name="pair_list_free")
      use m_type_pair_list
      type(c_pair_list),intent(inout) :: pl
    end subroutine pair_list_free_c

    subroutine pair_list_seti(pl, key, i, len) bind(C, name="pair_list_seti")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*)
      integer(kind=c_int) :: i, len
    end subroutine pair_list_seti

    subroutine pair_list_setr(pl, key, r, len) bind(C, name="pair_list_setr")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*)
      integer(kind=c_int) :: len
      real(kind=c_double) :: r
    end subroutine pair_list_setr

    subroutine pair_list_sets(pl, key, s, len, len_s) bind(C, name="pair_list_sets")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*), s(*)
      integer(kind=c_int) :: len, len_s
      real(kind=c_double) :: r
    end subroutine pair_list_sets

    subroutine pair_list_get_c(pl, key, type_code, i, r, s, len, len_s) bind(C, name="pair_list_get_")
      use m_type_pair_list
      type(c_pair_list) :: pl
      character(kind=c_char) :: key(*), s(*)
      integer(kind=c_int) :: i, type_code, len, len_s
      real(kind=c_double) :: r
    end subroutine pair_list_get_c

    subroutine pair_list_look_c(pl, key, type_code, i, r, s, len, len_s) bind(C, name="pair_list_look_")
      use m_type_pair_list
      type(c_pair_list) :: pl
      integer(kind=c_int) :: type_code, i, len, len_s
      character(kind=c_char) :: key(len), s(len_s)
      real(kind=c_double) :: r
    end subroutine pair_list_look_c

  end interface

! -------------------------------------------------------------------------------
! -                                                                             -
! -                          Pure Fortran Wrapper                               -
! -                                                                             -
! -------------------------------------------------------------------------------
  contains
!!***

!!****f* m_pair_list/pair_list_length
!! NAME
!! pair_list_length
!!
!! FUNCTION
!!  get the number of pair stored in pl
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function pair_list_length(pl) result(length)
  class(pair_list),intent(in) :: pl
  integer :: length
  length = pl%plc%length
end function pair_list_length
!!***

!!****f* m_pair_list/pair_list_get
!! NAME
!! pair_list_get
!!
!! FUNCTION
!!  Get the value associated with a key, only one of i and r is modified
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!  key <character(kind=c_char,len=*)>=
!!  s <character(kind=c_char,len=*)>=
!!
!! OUTPUT
!!  i <integer(kind=c_int)>=
!!  type_code <integer(kind=c_int)>=
!!      0 if the value was an integer (and so that i is setted)
!!      1 if the value was a real number (and so that r is setted)
!!      2 if the value was a string (and so that s is setted)
!!     -1 if the key was not present (neither i nor r are setted)
!!     -2 if the list is empty (neither i nor r are setted)
!!  r <real(kind=c_double)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_get(pl, key, type_code, i, r, s)
  class(pair_list),intent(in) :: pl
  character(kind=c_char,len=*),intent(in) :: key, s
  integer(kind=c_int),intent(out) :: i, type_code
  real(kind=c_double),intent(out) :: r
  call pair_list_get_c(pl%plc, trim(key), type_code, i, r, s, len_trim(key), len(s))
end subroutine pair_list_get
!!***

!!****f* m_pair_list/pair_list_look
!! NAME
!! pair_list_look
!!
!! FUNCTION
!!  pair_list variables have a cursor wich point onto an arbitrary element
!!  of the list. pair_list_look allow to extract the key-value pair from
!!  that element
!!
!!  If key is shorter than the actual key of the pair, only available space
!!  is used resulting in truncated key
!!  If key is longer than the actual key remaining space is filled with spaces
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!
!! OUTPUT
!!  key <character(kind=c_char,len=*)>=
!!  s <character(kind=c_char,len=*)>=
!!  type_code <integer(kind=c_int)>=
!!      1 if the value was a real number (and so that r is setted)
!!      0 if the value was an integer (and so that i is setted)
!!     -2 if the cursor is null (list is empty or end have been reached)
!!  i <integer(kind=c_int)>=
!!  r <real(kind=c_double)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_look(pl, key, type_code, i, r, s)
  use m_type_pair_list
  class(pair_list),intent(in) :: pl
  character(kind=c_char,len=*),intent(out) :: key, s
  integer(kind=c_int),intent(out) :: type_code, i
  real(kind=c_double),intent(out) :: r
  call pair_list_look_c(pl%plc, key, type_code, i, r, s, len(key), len(s))
end subroutine pair_list_look
!!***

!!****f* m_pair_list/pair_list_next
!! NAME
!! pair_list_next
!!
!! FUNCTION
!!  have the cursor (cf: pair_list_look) moving forward of one element.
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE
    subroutine pair_list_next(pl)
      class(pair_list),intent(in) :: pl
      call pair_list_next_c(pl%plc)
    end subroutine pair_list_next
!!***

!!****f* m_pair_list/pair_list_free
!! NAME
!! pair_list_free
!!
!! FUNCTION
!!  free memory occupied by the list (not the pair_list variable itself !)
!!  and reset the pair_list variable (it can be reused as an empty list)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_free(pl)
  class(pair_list),intent(inout) :: pl
  call pair_list_free_c(pl%plc)
end subroutine pair_list_free
!!***

!!****f* m_pair_list/pair_list_set
!! NAME
!! pair_list_set
!!
!! FUNCTION
!!  set a key-value par into the list. If the key is already present, the
!!  corresponding pair is updated. If not the pair is created.
!!  Only one of i and r should be provided (i is the default if both are
!!  provided). Nothing happen if none of them are provided.
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!  key <character(len=*)>=
!!  i <integer>=optional
!!  r <real(kind=c_double)>=optional
!!  s <character(len=*)>=optional
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_set(pl, key, i, r, s)

 class(pair_list),intent(in) :: pl
 character(len=*),intent(in) :: key
 integer,intent(in),optional :: i
 real(kind=c_double),intent(in),optional :: r
 character(len=*),intent(in),optional :: s

 if (present(i)) then
   call pair_list_seti(pl%plc, trim(key), i, len_trim(key))
 else if (present(r)) then
   call pair_list_setr(pl%plc, trim(key), r, len_trim(key))
 else if (present(s)) then
   call pair_list_sets(pl%plc, trim(key), s, len_trim(key), len_trim(s))
 end if

end subroutine pair_list_set
!!***

!!****f* m_pair_list/pair_list_set_keys
!! NAME
!! pair_list_set_keys
!!
!! FUNCTION
!!  Set the value of a list of comma-separated keys.
!!
!!  Example
!!
!!  d%set_keys("foo, bar", ivals=[1, 2])
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!  keylist <character(len=*)>=
!!  i <integer>=optional
!!  r <real(kind=c_double)>=optional
!!  s <character(len=*)>=optional
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_set_keys(pl, keylist, ivals, rvals) !, svals)

 class(pair_list),intent(in) :: pl
 character(len=*),intent(in) :: keylist
 integer,intent(in),optional :: ivals(:)
 real(kind=c_double),intent(in),optional :: rvals(:)
 !character(len=*),intent(in),optional :: svals(:)

!Local variables-------------------------------
 integer :: i, n, start, stp
 character(len=len(keylist)) :: key
! *************************************************************************

 n = 1
 do i=1,len_trim(keylist)
   if (keylist(i:i) == ",") n = n + 1
 end do

 start = 1
 do i=1,n
   stp = index(keylist(start:), ",")
   if (stp == 0) then
     key = keylist(start:)
   else
     key = keylist(start: start + stp - 2)
     start = start + stp
     ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
   end if
   key = adjustl(key)

   if (present(ivals)) then
     ABI_CHECK(size(ivals) == n, "size(ivals) != n")
     call pair_list_seti(pl%plc, trim(key), ivals(i), len_trim(key))

   else if (present(rvals)) then
     ABI_CHECK(size(rvals) == n, "size(rvals) != n")
     call pair_list_setr(pl%plc, trim(key), rvals(i), len_trim(key))

   !else if (present(svals)) then
   !  TODO: Pass single string with comma-separated tokens.
   !  ABI_CHECK(size(svals) == n, "size(svals) != n")
   !  call pair_list_sets(pl%plc, trim(key), svals(i), len_trim(key), len_trim(svals(i)))
   end if
 end do

end subroutine pair_list_set_keys
!!***

!!****f* m_pair_list/pair_list_set_keys_to_null
!! NAME
!! pair_list_set_keys_to_null
!!
!! FUNCTION
!!  Set the value of a list of comma-separated keys to null
!!
!!  Example:
!!
!!      dict%set_keys_to_null("foo, bar")
!!
!! INPUTS
!!  keylist: List of comma-separated keys
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_set_keys_to_null(pl, keylist)

 class(pair_list),intent(in) :: pl
 character(len=*),intent(in) :: keylist

!Local variables-------------------------------
 integer :: start, stp
! *************************************************************************

 start = 1
 do
   stp = index(keylist(start:), ",")
   if (stp == 0) then
     call pl%set(adjustl(trim(keylist(start:))), s="null")
     exit
   else
     call pl%set(adjustl(trim(keylist(start:start+stp-2))), s="null")
     start = start + stp
     ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
   end if
 end do

end subroutine pair_list_set_keys_to_null
!!***

!!****f* m_pair_list/pair_list_restart
!! NAME
!! pair_list_restart
!!
!! FUNCTION
!!  have the cursor going back to the first element (cf: pair_list_next)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_restart(pl)

  class(pair_list),intent(inout) :: pl
  pl%plc%cursor = pl%plc%first;

end subroutine pair_list_restart
!!***

!!****f* m_pair_list/pair_list_iter
!! NAME
!! pair_list_iter
!!
!! FUNCTION
!!  equivalent to pair_list_look followed by pair_list_next
!!
!! INPUTS
!!  pl <class(pair_list)>=
!!
!! OUTPUT
!!  key <character(len=*)>=
!!  type_code <integer>=
!!  i <integer>=
!!  r <real(kind=c_double)>=
!!  s <character(len=*)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine pair_list_iter(pl, key, type_code, i, r, s)

  class(pair_list),intent(in) :: pl
  character(len=*),intent(out) :: key
  integer,intent(out) :: type_code
  integer,intent(out) :: i
  real(kind=c_double),intent(out) :: r
  character(len=*),intent(out) :: s

  call pair_list_look(pl, key, type_code, i, r, s)
  if(type_code >= 0) call pair_list_next_c(pl%plc)

end subroutine pair_list_iter
!!***

end module m_pair_list
!!***

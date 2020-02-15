!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_yaml
!! NAME
!!  m_yaml
!!
!! FUNCTION
!!  This module defines low-level routines to format data into YAML documents.
!!  Supported data include numeric arrays of one and two dimensions,
!!  strings, numbers, dictionaries from m_pair_list and 1D arrays of dictionaries.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (TC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!   m_stream_string, m_pair_list
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_yaml

 use defs_basis
#ifdef HAVE_FC_IEEE_ARITHMETIC
 use ieee_arithmetic
#endif
 use m_errors
 use m_pair_list
 use m_stream_string

 use m_fstrings, only : sjoin, char_count, itoa, sjoin
 use m_io_tools, only : is_open

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_yaml/yamldoc_t
!! NAME
!! yamldoc_t
!!
!! FUNCTION
!! High-level API to write (simple) Yaml documents.
!!
!! SOURCE

 type,public :: yamldoc_t

   integer :: use_yaml = 1
   ! Temporary flag used to deactivate Yaml output

   integer :: default_keysize = 30
   ! Default key size

   integer :: default_stringsize = 500
   ! Default string size

   integer :: default_width = 0
   ! impose a minimum width of the field name side of the column (padding with spaces)

   integer :: default_multiline_trig = 5
   ! minimum number of elements before switching to multiline representation.

   character(len=20) :: default_ifmt = '(I0)'
   !character(len=20) :: default_ifmt = '(I8)'
   ! Default format for integer

   character(len=20) :: default_rfmt = '(ES16.8)'
   ! Default format for real

   character(len=20) :: default_kfmt = "(A)"
   ! Default format for keys

   character(len=20) :: default_sfmt = "(A)"
   ! Default format for strings

   type(stream_string) :: stream
   ! Stream object used to build yaml string.

 contains

   procedure :: write_and_free => yamldoc_write_and_free
    ! Write Yaml document to unit and free memory.

   procedure :: add_real => yamldoc_add_real
     ! Add a real number field to a document

   procedure :: add_reals => yamldoc_add_reals
     ! Add a list of real number fields to a document

   procedure :: add_paired_real2d => yamldoc_add_paired_real2d
     !  Add a field containing two 2D array of real numbers with the same shape.

   procedure :: add_int => yamldoc_add_int
     ! Add an integer field to a document

   procedure :: add_ints => yamldoc_add_ints
     ! Add an integer field to a document

   procedure :: add_string => yamldoc_add_string
     ! Add a string field to a document

   procedure :: add_real1d => yamldoc_add_real1d
     ! Add a field containing a 1D array of real numbers

   procedure :: add_real2d => yamldoc_add_real2d
     ! Add a field containing a 2D real number array

   procedure :: add_int1d => yamldoc_add_int1d
     ! Add a field containing a 1D integer array

   procedure :: add_int2d => yamldoc_add_int2d
     ! Add a field containing a 2D integer array

   !procedure :: add_tabular => yamldoc_add_tabular
     ! Add a field with a complete table data

   procedure :: open_tabular => yamldoc_open_tabular
     ! Open a field for tabular data

   procedure :: add_tabular_line => yamldoc_add_tabular_line
     ! Add a line of tabular data in an already opened table field

   procedure :: add_dict => yamldoc_add_dict
     ! Add a field containing a dictionary/pair_list

   procedure :: add_dictlist => yamldoc_add_dictlist
     ! Add a field containing a list of dictionaries/array of pair_list

   procedure :: set_keys_to_string => yamldoc_set_keys_to_string
     ! Set all keys to a commong (string) value

 end type yamldoc_t
!!***

 public :: yamldoc_open
  ! Open a yaml document

 public :: yaml_single_dict
  ! Create a full document from a single dictionary

 public :: yaml_iterstart
  ! Set the value of the iteration indices used to build the iteration_state dict in the Yaml documents

 character(len=1),parameter :: eol = char(10)

 ! This is a list of reserved_keywords that shall not be used as keys in Yaml dictionaries.
 character(len=12),parameter :: reserved_keywords(10) = [character(len=12) :: &
   "tol_abs", "tol_rel", "tol_vec", "tol_eq", "ignore", &
   "ceil", "equation", "equations", "callback", "callbacks"]

 ! Global variables used to save the iteration state in Abinit.
 ! Set by yaml_iterstart
 integer,save,protected :: DTSET_IDX = -1
 integer,save,protected :: TIMIMAGE_IDX = -1
 integer,save,protected :: IMAGE_IDX = -1
 integer,save,protected :: ITIME_IDX = -1
 integer,save,protected :: ICYCLE_IDX = -1

contains

!!****f* m_yaml/yaml_iterstart
!! NAME
!! yaml_iterstart
!!
!! FUNCTION
!!  Mark the start of an iteration named by label and numbered by file
!!
!! INPUTS
!!  label=key name
!!  val=value
!!  [newline] = set to false to prevent adding newlines after fields
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_iterstart(label, val, unit, use_yaml, newline)

!Arguments ------------------------------------
 integer,intent(in) :: val, unit, use_yaml
 character(len=*),intent(in) :: label
 logical,intent(in),optional :: newline

!Local variables-------------------------------
 character(len=6) :: tmp_i
 logical :: nl
 type(stream_string) :: stream
! *************************************************************************

 if (unit == dev_null) return

 select case (label)
 case ("dtset")
   DTSET_IDX = val
   TIMIMAGE_IDX = -1
   IMAGE_IDX = -1
   ITIME_IDX = -1
   ICYCLE_IDX = -1
 case ("timimage")
   TIMIMAGE_IDX = val
 case ("image")
   IMAGE_IDX = val
 case ("itime")
   ITIME_IDX = val
 case ("icycle")
   ICYCLE_IDX = val
 case default
   MSG_ERROR(sjoin("Invalid value for label:", label))
 end select

 if (use_yaml == 1) then
   ABI_DEFAULT(nl, newline, .true.)
   write(tmp_i, '(I6)') val
   call stream%push('--- !IterStart'//eol//label//':'//tmp_i//eol//'...')
   if (nl) call stream%push(eol)
   call stream%flush(unit)
 end if

end subroutine yaml_iterstart
!!***

!!****f* m_yaml/yamldoc_open
!! NAME
!! yamldoc_open
!!
!! FUNCTION
!!  Open a yaml document
!!
!! INPUTS
!!  tag: add a tag to the field
!!  comment: string with comment.
!!  [newline]: optional, set to false to prevent adding newlines after fields
!!  [width]: optional, impose a minimum width of the field name side of the column (padding with spaces)
!!  [int_fmt]: Default format for integers.
!!  [real_fmt]: Default format for real.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

type(yamldoc_t) function yamldoc_open(tag, comment, newline, width, int_fmt, real_fmt) result(new)

!Arguments ------------------------------------
 character(len=*),intent(in) :: tag, comment
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width
 character(len=*),optional,intent(in) :: int_fmt, real_fmt

!Local variables-------------------------------
 logical :: nl
 type(pair_list) :: dict
! *************************************************************************

 ABI_DEFAULT(nl, newline, .False.)

 if (present(width)) new%default_width = width
 if (present(int_fmt)) new%default_ifmt = int_fmt
 if (present(real_fmt)) new%default_rfmt = real_fmt

 call new%stream%push('---'//' !'//trim(tag)//ch10)

 if (DTSET_IDX /= -1) then
   ! Write dictionary with iteration state.
   call dict%set('dtset', i=DTSET_IDX)
   if (TIMIMAGE_IDX /= -1) call dict%set("timimage", i=TIMIMAGE_IDX)
   if (IMAGE_IDX /= -1) call dict%set("image", i=IMAGE_IDX)
   if (ITIME_IDX /= -1) call dict%set("itime", i=ITIME_IDX)
   if (ICYCLE_IDX /= -1) call dict%set("icycle", i=ICYCLE_IDX)
   call new%add_dict('iteration_state', dict, int_fmt="(i0)")
   call dict%free()
 end if

 if (comment /= '') then
   call new%stream%push('comment')
   if (new%default_width > 7) call new%stream%push(repeat(' ', new%default_width - 7))
   call new%stream%push(': ')
   call yaml_print_string(new%stream, comment)
   call new%stream%push(eol)
 end if
 if (nl) call new%stream%push(eol)

end function yamldoc_open
!!***

!!****f* m_yaml/yamldoc_add_real
!! NAME
!! yamldoc_add_real
!!
!! FUNCTION
!!  Add a real number field to a document
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [real_fmt] = override the default formatting
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_real(self, label, val, tag, real_fmt, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: label
 real(dp),intent(in) :: val
 character(len=*),intent(in),optional :: tag, real_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 character(len=50) :: tmp_r
 character(len=30) :: rfmt
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call self%stream%push(' ')
 call format_real(val, tmp_r, trim(rfmt))
 call self%stream%push(trim(tmp_r))
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_real
!!***

!!****f* m_yaml/yamldoc_add_reals
!! NAME
!! yamldoc_add_reals
!!
!! FUNCTION
!!  Add a list of real numbers to the document
!!
!! INPUTS
!!  keylist = List of comma-separated keywords
!!  values = List of values
!!  [real_fmt] = override the default formatting
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!  [dict_key]=If present, a dictionary with key `dict_key` is created instead of a list.
!!  [multiline_trig] = optional minimum number of elements before switching to multiline representation
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_reals(self, keylist, values, real_fmt, width, dict_key, multiline_trig)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: keylist
 real(dp),intent(in) :: values(:)
 character(len=*),intent(in),optional :: real_fmt, dict_key
 integer,intent(in),optional :: width, multiline_trig

!Local variables-------------------------------
 integer :: i, n, w, start, stp, vmax
 character(len=30) :: rfmt
 type(pair_list) :: dict
! *************************************************************************

 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)

 n = char_count(keylist, ",") + 1
 ABI_CHECK(size(values) == n, "size of values != len(tokens)")

 start = 1

 if (.not. present(dict_key)) then
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call self%add_real(adjustl(keylist(start:)), values(i), real_fmt=rfmt, width=w)
     else
       call self%add_real(adjustl(keylist(start: start + stp - 2)), values(i), real_fmt=rfmt, width=w)
       start = start + stp
       ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
     end if
   end do

 else

   ! Create and insert dictionary.
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call dict%set(adjustl(keylist(start:)), r=values(i))
     else
       call dict%set(adjustl(keylist(start: start + stp - 2)), r=values(i))
       start = start + stp
       ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
     end if
   end do
   ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)
   call self%add_dict(trim(dict_key), dict, multiline_trig=vmax, real_fmt=rfmt, width=w)
   call dict%free()
 end if

end subroutine yamldoc_add_reals
!!***

!!****f* m_yaml/yamldoc_add_int
!! NAME
!! yamldoc_add_int
!!
!! FUNCTION
!!  Add an integer field to a document
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [int_fmt] = optional  override the default formatting
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_int(self, label, val, tag, int_fmt, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in) :: val
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag, int_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 character(50) :: tmp_i
 character(len=30) :: ifmt
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call self%stream%push(' ')
 write(tmp_i, trim(ifmt)) val
 call self%stream%push(trim(tmp_i))
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_int
!!***

!!****f* m_yaml/yamldoc_add_ints
!! NAME
!! yamldoc_add_ints
!!
!! FUNCTION
!!  Add a list of integer numbers to the document
!!
!! INPUTS
!!  keylist = List of comma-separated keywords
!!  values = List of values
!!  [int_fmt] = override the default formatting
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!  [dict_key]=If present, a dictionary with key `dict_key` is created instead of a list.
!!  [multiline_trig] = optional minimum number of elements before switching to multiline representation
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_ints(self, keylist, values, int_fmt, width, dict_key, multiline_trig)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: keylist
 integer,intent(in) :: values(:)
 character(len=*),intent(in),optional :: int_fmt, dict_key
 integer,intent(in),optional :: width, multiline_trig

!Local variables-------------------------------
 integer :: i, n, w, start, stp, vmax
 character(len=30) :: ifmt
 type(pair_list) :: dict
! *************************************************************************

 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)

 n = char_count(keylist, ",") + 1
 ABI_CHECK(size(values) == n, "size of values != len(tokens)")

 start = 1

 if (.not. present(dict_key)) then
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call self%add_int(adjustl(keylist(start:)), values(i), int_fmt=ifmt, width=w)
     else
       call self%add_int(adjustl(keylist(start: start + stp - 2)), values(i), int_fmt=ifmt, width=w)
       start = start + stp
       ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
     end if
   end do

 else
   ! Create and insert dictionary.
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call dict%set(adjustl(keylist(start:)), i=values(i))
     else
       call dict%set(adjustl(keylist(start: start + stp - 2)), i=values(i))
     end if
   end do
   ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)
   call self%add_dict(trim(dict_key), dict, multiline_trig=vmax, int_fmt=ifmt, width=w)
   call dict%free()
 end if

end subroutine yamldoc_add_ints
!!***

!!****f* m_yaml/yamldoc_add_string
!! NAME
!! yamldoc_add_string
!!
!! FUNCTION
!!  Add a string field to a document
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_string(self, label, val, tag, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: val
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call self%stream%push(' ')
 call yaml_print_string(self%stream, val)
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_string
!!***

!!****f* m_yaml/yamldoc_add_real1d
!! NAME
!! yamldoc_add_real1d
!!
!! FUNCTION
!!  Add a field containing a 1D array of real numbers
!!
!! INPUTS
!!  label = key name
!!  arr(:)
!!  [multiline_trig] = optional minimum number of elements before switching to multiline representation
!!  [tag] = optional, add a tag to the field
!!  [real_fmt] = override the default formatting
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_real1d(self, label, arr, tag, real_fmt, multiline_trig, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in),optional :: multiline_trig
 real(dp),intent(in) :: arr(:)
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag, real_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w, length, vmax
 character(len=30) :: rfmt
 logical :: nl
! *************************************************************************

 length = size(arr)

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call yaml_print_real1d(self%stream, length, arr, trim(rfmt), vmax)
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_real1d
!!***

!!****f* m_yaml/yamldoc_add_int1d
!! NAME
!! yamldoc_add_int1d
!!
!! FUNCTION
!!  Add a field containing a 1D integer array
!!
!! INPUTS
!!  label = key name
!!  arr(:) <integer>=
!!  [multiline_trig] = optional minimum number of elements before switching to multiline representation
!!  [tag] : add a tag to the field
!!  int_fmt <character(len=*)>=optional override the default formatting
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_int1d(self, label, arr, tag, int_fmt, multiline_trig, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in),optional :: multiline_trig
 integer,intent(in) :: arr(:)
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag, int_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 character(len=30) :: ifmt
 integer :: w, length, vmax
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)
 length = size(arr)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call yaml_print_int1d(self%stream, length, arr, trim(ifmt), vmax)
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_int1d
!!***

!!****f* m_yaml/yamldoc_add_dict
!! NAME
!! yamldoc_add_dict
!!
!! FUNCTION
!!  Add a field containing a dictionary/pair_list
!!
!! INPUTS
!!  label <character(len=*)>=
!!  pl <type(pair_list)>=
!!  string_size <integer>=optional maximum storage size for strings found in a pair_list
!!  key_size <integer>=optional maximum storage size for keys of a pair_list
!!  multiline_trig <integer>=optional minimum number of elements before switching to multiline representation
!!  tag <character(len=*)>=optional  add a tag to the field
!!  key_fmt <character(len=*)>=optional  override the default formatting
!!  int_fmt <character(len=*)>=optional  override the default formatting
!!  real_fmt <character(len=*)>=optional  override the default formatting
!!  string_fmt <character(len=*)>=optional  override the default formatting
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!  width <integer>=optional impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  pl <type(pair_list)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_dict(self, label, pl, tag, key_size, string_size, key_fmt, &
                            int_fmt, real_fmt, string_fmt, multiline_trig, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 type(pair_list),intent(inout) :: pl
 character(len=*),intent(in) :: label
 integer,intent(in),optional :: string_size, key_size, multiline_trig
 character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w, vmax, ks, ss
 character(len=30) :: kfmt, ifmt, rfmt, sfmt
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(ks, key_size, self%default_keysize)
 ABI_DEFAULT(ss, string_size, self%default_stringsize)
 ABI_DEFAULT(kfmt, key_fmt, self%default_kfmt)
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)
 ABI_DEFAULT(sfmt, string_fmt, self%default_sfmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 call yaml_print_dict(self%stream, pl, ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_dict
!!***

!!****f* m_yaml/yamldoc_add_real2d
!! NAME
!! yamldoc_add_real2d
!!
!! FUNCTION
!!  Add a field containing a 2D array of real numbers
!!
!! INPUTS
!!  label = key name
!!  arr(:, :) = input array.
!!  [slist(:)]= List of strings (same length as the first dim or second dime of arr, depending on mode).
!!    If present, the string will be included in the the row.
!!  [tag]= add a tag to the field
!!  [real_fmt]= override the default formatting
!!  [multiline_trig]: optional minimum number of elements before switching to multiline representation
!!  [newline]: set to false to prevent adding newlines after fields
!!  [width]: impose a minimum width of the field name side of the column (padding with spaces)
!!  [mode]: "T" to write the transpose of arr i.e columns become rows in output (DEFAULT), "N" for normal order
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_real2d(self, label, arr, slist, tag, real_fmt, multiline_trig, newline, width, mode)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 real(dp),intent(in) :: arr(:, :)
 character(len=*),intent(in) :: label
 character(len=*),optional,intent(in) :: slist(:)
 character(len=*),intent(in),optional :: tag, real_fmt
 integer,intent(in),optional :: multiline_trig
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width
 character(len=1),intent(in),optional :: mode

!Local variables-------------------------------
 integer :: m, n, w, i, vmax
 real(dp) :: line(max(size(arr, dim=1), size(arr, dim=2)))
 character(len=30) :: rfmt
 character(len=1) :: my_mode
 logical :: nl
! *************************************************************************

 m = size(arr, dim=1)
 n = size(arr, dim=2)

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(my_mode, mode, "T")
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 if (my_mode == "T") then
   do i=1,n
     call self%stream%push(eol//'-')
     line(1:m) = arr(:,i)
     if (.not. present(slist)) then
       call yaml_print_real1d(self%stream, m, line, rfmt, vmax)
     else
       call yaml_print_real1d(self%stream, m, line, rfmt, vmax, string=slist(i))
     end if
   end do
 else
   do i=1,m
     call self%stream%push(eol//'-')
     line(1:n) = arr(i,:)
     if (.not. present(slist)) then
       call yaml_print_real1d(self%stream, n, line, rfmt, vmax)
     else
       call yaml_print_real1d(self%stream, n, line, rfmt, vmax, string=slist(i))
     end if
   end do
 end if

 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_real2d
!!***

!!****f* m_yaml/yamldoc_add_paired_real2d
!! NAME
!! yamldoc_add_paired_real2d
!!
!! FUNCTION
!!  Add a field containing two 2D real arrays with the same shape.
!!
!!  Example:
!!    cartesian_forces_and_xred:
!!    - [ [ -0.0000E+00,  -0.0000E+00,  -0.0000E+00, ], [  0.0000E+00,   0.0000E+00,   0.0000E+00, ] ]
!!    - [ [ -0.0000E+00,  -0.0000E+00,  -0.0000E+00, ], [  2.5000E-01,   2.5000E-01,   2.5000E-01, ] ]
!!
!! INPUTS
!!  label = key name
!!  arr1(:,:), arr2(:,:) = input arrays.
!!  [slist(:)]= List of strings (same length as the first dim or second dime of arr, depending on mode).
!!    If present, the string will be included in the the row.
!!  [tag]= add a tag to the field
!!  [real_fmt]= override the default formatting
!!  [multiline_trig]: optional minimum number of elements before switching to multiline representation
!!  [newline]: set to false to prevent adding newlines after fields
!!  [width]: impose a minimum width of the field name side of the column (padding with spaces)
!!  [mode]: "T" to write the transpose of arr i.e columns become rows in output (DEFAULT), "N" for normal order
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_paired_real2d(self, label, arr1, arr2, slist, tag, real_fmt, multiline_trig, newline, width, mode)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 real(dp),intent(in) :: arr1(:, :), arr2(:,:)
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag, real_fmt
 integer,intent(in),optional :: multiline_trig
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width
 character(len=1),intent(in),optional :: mode
 character(len=*),optional,intent(in) :: slist(:)

!Local variables-------------------------------
 integer :: m, n, w, i, vmax
 real(dp) :: line(2 * max(size(arr1, dim=1), size(arr1, dim=2)))
 character(len=30) :: rfmt
 character(len=1) :: my_mode
 logical :: nl
! *************************************************************************

 m = size(arr1, dim=1)
 n = size(arr1, dim=2)

 ABI_CHECK(all(shape(arr1) == shape(arr2)), "arr1 and arr2 must have same shape")

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(my_mode, mode, "T")
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 if (my_mode == "T") then
   if (present(slist)) then
     ABI_CHECK(size(slist) == n, "size(slist) != n")
   end if
   do i=1,n
     call self%stream%push(eol//'- [')
     line(1:m) = arr1(:,i)
     call yaml_print_real1d(self%stream, m, line, rfmt, vmax)
     call self%stream%push(',')
     line(1:m) = arr2(:,i)
     call yaml_print_real1d(self%stream, m, line, rfmt, vmax)
     if (present(slist)) call self%stream%push(', '//trim(slist(i)))
     call self%stream%push(' ]')
   end do
 else
   if (present(slist)) then
     ABI_CHECK(size(slist) == n, "size(slist) != m")
   end if
   do i=1,m
     call self%stream%push(eol//'- [')
     line(1:n) = arr1(i,:)
     call yaml_print_real1d(self%stream, n, line, rfmt, vmax)
     call self%stream%push(',')
     line(1:n) = arr2(i,:)
     call yaml_print_real1d(self%stream, n, line, rfmt, vmax)
     if (present(slist)) call self%stream%push(', '//trim(slist(i)))
     call self%stream%push(']')
   end do
 end if

 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_paired_real2d
!!***

!!****f* m_yaml/yamldoc_add_int2d
!! NAME
!! yamldoc_add_int2d
!!
!! FUNCTION
!!  Add a field containing a 2D integer array
!!
!! INPUTS
!!  label = key name
!!  arr(:, :) <integer>=
!!  [slist(:)]= List of strings (same length as the first dim or second dime of arr, depending on mode).
!!    If present, the string will be included in the the row.
!!  [tag]= add a tag to the field
!!  [int_fmt]: override the default formatting
!!  multiline_trig <integer>=optional minimum number of elements before switching to multiline representation
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!  [mode] = "T" to write the transpose of arr i.e columns become rows in output (DEFAULT), "N" for normal order
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_int2d(self, label, arr, slist, tag, int_fmt, multiline_trig, newline, width, mode)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in) :: arr(:, :)
 character(len=*),intent(in) :: label
 character(len=*),optional,intent(in) :: slist(:)
 character(len=*),intent(in),optional :: tag, int_fmt
 integer,intent(in),optional :: multiline_trig
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width
 character(len=1),intent(in),optional :: mode

!Local variables-------------------------------
 integer :: m, n, w, i, vmax
 integer :: line(max(size(arr, dim=1), size(arr, dim=2)))
 character(len=30) :: ifmt
 character(len=1) :: my_mode
 logical :: nl
! *************************************************************************

 m = size(arr, dim=1)
 n = size(arr, dim=2)

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(my_mode, mode, "T")
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if

 if (my_mode == "T") then
   do i=1,n
     call self%stream%push(eol//'-')
     line(1:m) = arr(:,i)
     if (.not. present(slist)) then
       call yaml_print_int1d(self%stream, m, line, ifmt, vmax)
     else
       call yaml_print_int1d(self%stream, m, line, ifmt, vmax, string=slist(i))
     end if
   end do
 else
   do i=1,m
     call self%stream%push(eol//'-')
     line(1:n) = arr(i,:)
     if (.not. present(slist)) then
       call yaml_print_int1d(self%stream, n, line, ifmt, vmax)
     else
       call yaml_print_int1d(self%stream, n, line, ifmt, vmax, string=slist(i))
     end if
   end do
 end if

 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_int2d
!!***

!!****f* m_yaml/yamldoc_add_dictlist
!! NAME
!! yamldoc_add_dictlist
!!
!! FUNCTION
!!  Add a field containing a list of dictionaries/array of pair_list
!!
!! INPUTS
!!  label = key name
!!  n <integer>=
!!  plarr(n) <type(pair_list)>=
!!  key_size <integer>=optional maximum storage size for keys of a pair_list
!!  string_size <integer>=optional maximum storage size for strings of a pair_list
!!  multiline_trig <integer>=optional minimum number of elements before switching to multiline representation
!!  [tag]= add a tag to the field
!!  key_fmt <character(len=*)>=optional  override the default formatting
!!  int_fmt <character(len=*)>=optional  override the default formatting
!!  real_fmt <character(len=*)>=optional  override the default formatting
!!  string_fmt <character(len=*)>=optional  override the default formatting
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_dictlist(self, label, n, plarr, tag, key_size, string_size, key_fmt, int_fmt, &
                                real_fmt, string_fmt, multiline_trig, newline, width)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in) :: n
 type(pair_list),intent(inout) :: plarr(n)
 character(len=*),intent(in) :: label
 integer,intent(in),optional :: key_size, string_size
 integer,intent(in),optional :: multiline_trig
 character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 character(len=30) :: kfmt, ifmt, rfmt, sfmt
 integer :: vmax, ks, i, ss
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(w, width, self%default_width)
 ABI_DEFAULT(kfmt, key_fmt, self%default_kfmt)
 ABI_DEFAULT(rfmt, real_fmt, self%default_rfmt)
 ABI_DEFAULT(ifmt, int_fmt, self%default_ifmt)
 ABI_DEFAULT(sfmt, string_fmt, self%default_sfmt)
 ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)
 ABI_DEFAULT(ks, key_size, self%default_keysize)
 ABI_DEFAULT(ss, string_size, self%default_keysize)

 if (present(tag)) then
   call yaml_start_field(self%stream, label, width=w, tag=tag)
 else
   call yaml_start_field(self%stream, label, width=w)
 end if
 call self%stream%push(eol)

 do i=1,n
   call self%stream%push('- ')
   call yaml_print_dict(self%stream, plarr(i), ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
   if (nl .or. i/=n) call self%stream%push(eol)
 end do

end subroutine yamldoc_add_dictlist
!!***

!!****f* m_yaml/yamldoc_open_tabular
!! NAME
!! yamldoc_open_tabular
!!
!! FUNCTION
!!  Open a field for tabular data
!!
!! INPUTS
!!  label = key name
!!  tag <character(len=*)>=optional  add a tag to the field
!!  [newline] = set to false to prevent adding newlines after fields
!!  [indent] = optional number of spaces to add to the header
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_open_tabular(self, label, tag, indent, newline)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: indent

!Local variables-------------------------------
 integer :: n
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(n, indent, 4)

 if (n > 4) then
   call self%stream%push(repeat(' ', n-4))
 end if

 if (present(tag)) then
   call yaml_start_field(self%stream, label, tag=tag)
 else
   call yaml_start_field(self%stream, label, tag='Tabular')
 end if
 call self%stream%push(' |'//eol)

end subroutine yamldoc_open_tabular
!!***

!!****f* m_yaml/yamldoc_add_tabular_line
!! NAME
!! yamldoc_add_tabular_line
!!
!! FUNCTION
!!  Add a line of tabular data in an already opened table field
!!
!! INPUTS
!!  line <character(len=*)>=
!!  [newline] = set to false to prevent adding newlines after fields
!!  [indent] = optional number of spaces to add to the header
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_add_tabular_line(self, line, newline, indent)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: line
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: indent

!Local variables-------------------------------
 integer :: n
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(n, indent, 4)

 call self%stream%push(repeat(' ', n)//trim(line))
 if (nl) call self%stream%push(eol)

end subroutine yamldoc_add_tabular_line
!!***

!!****f* m_yaml/yamldoc_add_tabular
!! NAME
!! yamldoc_add_tabular
!!
!! FUNCTION
!!  Add a field with a complete table data
!!
!! INPUTS
!!  label <character(len=*)>=
!!  input <type(stream_string)>=stream containing an already built table
!!  tag <character(len=*)>=optional  add a tag to the field
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!  indent <integer>=optional number of spaces to add to each line
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

!subroutine yamldoc_add_tabular(self, label, input, tag,  newline, indent)
!
!!Arguments ------------------------------------
! class(yamldoc_t),intent(inout) :: self
! character(len=*),intent(in) :: label
! type(stream_string),intent(inout) :: input
! character(len=*),intent(in),optional :: tag
! logical,intent(in),optional :: newline
! integer,intent(in),optional :: indent
!
!!Local variables-------------------------------
! integer :: n
! character(len=100) :: t
! logical :: nl
!! *************************************************************************
!
! ABI_DEFAULT(nl, newline, .true.)
! ABI_DEFAULT(n, indent, 4)
! ABI_DEFAULT(t, tag, 'Tabular')
!
! call yaml_open_tabular(label, tag=t, stream=self%stream, newline=nl)
!
! if (n > 4) call self%stream%push(repeat(' ', n - 4))
!
! call write_indent(input, self%stream, n)
! if (nl) call self%stream%push(eol)
!
!end subroutine yamldoc_add_tabular
!!***

!!****f* m_yaml/yaml_single_dict
!! NAME
!! yaml_single_dict
!!
!! FUNCTION
!!  Create a full document from a single dictionary
!!
!! INPUTS
!!  unit
!!  tag <character(len=*)>=
!!  comment <character(len=*)>=
!!  pl <type(pair_list)>=
!!  key_size <integer>=maximum storage size for the keys of pl
!!  string_size <integer>=maximum storage size for the strings found in pl
!!  tag <character(len=*)>=optional  add a tag to the field
!!  int_fmt <character(len=*)>=optional  override the default formatting
!!  real_fmt <character(len=*)>=optional  override the default formatting
!!  string_fmt <character(len=*)>=optional  override the default formatting
!!  width <integer>=optional impose a minimum width of the field name side of the column (padding with spaces)
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!
!! OUTPUT
!!  pl <type(pair_list)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_single_dict(unit, tag, comment, pl, key_size, string_size, &
                            int_fmt, real_fmt, string_fmt, newline, width)

!Arguments ------------------------------------
 integer,intent(in) :: unit
 type(pair_list),intent(inout) :: pl
 character(len=*),intent(in) :: tag
 character(len=*),intent(in) :: comment
 integer,intent(in) :: key_size, string_size
 character(len=*),intent(in),optional :: int_fmt, real_fmt, string_fmt
 integer,intent(in), optional :: width
 logical,intent(in),optional :: newline

!Local variables-------------------------------
 type(yamldoc_t) :: doc
 character(len=30) :: ifmt, rfmt, sfmt
 character(len=string_size) :: vs, tmp_s
 character(len=key_size) :: key
 integer :: vi, k, type_code, w
 character(len=50) :: tmp_i, tmp_r
 real(dp) :: vr
 logical :: nl
! *************************************************************************

 ABI_DEFAULT(nl, newline, .true.)
 ABI_DEFAULT(rfmt, real_fmt, doc%default_rfmt)
 ABI_DEFAULT(ifmt, int_fmt, doc%default_ifmt)
 ABI_DEFAULT(sfmt, string_fmt, doc%default_sfmt)
 ABI_DEFAULT(w, width, doc%default_width)

 call doc%stream%push('--- !'//tag)

 if (comment /= '') then
   call doc%stream%push(eol)
   call yaml_start_field(doc%stream, 'comment', width=w)
   call yaml_print_string(doc%stream, comment)
 end if
 call doc%stream%push(eol)

 call pl%restart()
 do k=1,pl%length()
   call string_clear(key)
   call string_clear(vs)
   call pl%iter(key, type_code, vi, vr, vs)

   call yaml_start_field(doc%stream, trim(key), width=w)
   call doc%stream%push(' ')
   if (type_code == TC_INT) then
     call string_clear(tmp_i)
     write(tmp_i, ifmt) vi
     call doc%stream%push(trim(tmp_i))
   else if (type_code == TC_REAL) then
     call string_clear(tmp_r)
     call format_real(vr, tmp_r, rfmt)
     call doc%stream%push(trim(tmp_r))
   else if (type_code == TC_STRING) then
     call string_clear(tmp_s)
     write(tmp_s, sfmt) vs
     call yaml_print_string(doc%stream, trim(tmp_s))
   end if
   call doc%stream%push(eol)
 end do

 call doc%write_and_free(unit, newline=nl)

end subroutine yaml_single_dict
!!***

!!****f* m_yaml/yamldoc_write_and_free
!! NAME
!! yamldoc_write_and_free
!!
!! FUNCTION
!!  Close a previously opened document
!!
!! INPUTS
!!  [newline]= set to false to prevent adding newlines after fields. Default: True
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_write_and_free(self, unit, newline)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 integer,intent(in) :: unit
 logical,intent(in),optional :: newline

!Local variables-------------------------------
 logical :: nl
! *************************************************************************

 if (self%stream%length == 0) return
 ABI_DEFAULT(nl, newline, .true.)

 call self%stream%push('...')

 if (self%use_yaml == 1) then
   if (is_open(unit)) then
     write(unit, "(a)")""
     call self%stream%flush(unit, newline=nl)
   else
     call self%stream%free()
   end if
 else
   call self%stream%free()
 end if

end subroutine yamldoc_write_and_free
!!***

!!****f* m_yaml/yamldoc_set_keys_to_string
!! NAME
!! yamldoc_set_keys_to_string
!!
!! FUNCTION
!! Set all keys to a commong (string) value
!!
!! INPUTS
!!  keylist = List of comma-separated keywords
!!  svalue = String Value
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!  [dict_key]=If present, a dictionary with key `dict_key` is created instead of a list.
!!  [multiline_trig] = optional minimum number of elements before switching to multiline representation
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yamldoc_set_keys_to_string(self, keylist, svalue, dict_key, width, multiline_trig)

!Arguments ------------------------------------
 class(yamldoc_t),intent(inout) :: self
 character(len=*),intent(in) :: keylist, svalue
 character(len=*),intent(in),optional :: dict_key
 integer,intent(in),optional :: width, multiline_trig

!Local variables-------------------------------
 integer :: i, n, w, start, stp, vmax
 type(pair_list) :: dict

! *************************************************************************

 ABI_DEFAULT(w, width, self%default_width)

 n = char_count(keylist, ",") + 1
 start = 1

 if (.not. present(dict_key)) then
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call self%add_string(adjustl(keylist(start:)), svalue, width=w)
     else
       call self%add_string(adjustl(keylist(start: start + stp - 2)), svalue, width=w)
       start = start + stp
       ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
     end if
   end do

 else
   ! Create and insert dictionary.
   do i=1,n
     stp = index(keylist(start:), ",")
     if (stp == 0) then
       call dict%set(adjustl(keylist(start:)), s=svalue)
     else
       call dict%set(adjustl(keylist(start: start + stp - 2)), s=svalue)
       start = start + stp
       ABI_CHECK(start < len_trim(keylist), sjoin("Invalid keylist:", keylist))
     end if
   end do
   ABI_DEFAULT(vmax, multiline_trig, self%default_multiline_trig)
   call self%add_dict(trim(dict_key), dict, multiline_trig=vmax, width=w)
   call dict%free()
 end if

end subroutine yamldoc_set_keys_to_string
!!***

! private
subroutine string_clear(string)
  character(len=*),intent(inout) :: string
  string = repeat(' ', len(string))
end subroutine string_clear

subroutine format_real(val, dest, formt)
  real(dp),intent(in) :: val
  character(len=*),intent(out) :: dest
  character(len=*),intent(in) :: formt

#ifdef HAVE_FC_IEEE_ARITHMETIC
  if (ieee_is_nan(val)) then  ! NaN
    write(dest, '(a)') '.nan'
  else if (val == MAGIC_UNDEF) then
#else
  if (val == MAGIC_UNDEF) then
#endif
    write(dest, '(a)') 'undef'
  else
    write(dest, trim(formt)) val
  end if

end subroutine format_real

!subroutine write_indent(input, output, n)
! class(stream_string),intent(inout) :: input, output
! integer,intent(in) :: n
!
! integer :: buffstart, buffstop, length
! character(len=chunk_size) :: buffer
!
! do while (input%length > 0)
!   length = input%length
!   call input%pop_chunk(buffer)
!
!   buffstart = 1
!   buffstop = 1
!   do while (buffstart < min(length, chunk_size))
!     buffstop = index(buffer(buffstart:), eol)
!     if (buffstop > 0) then
!       call output%push(buffer(buffstart:buffstop))
!       call output%push(repeat(' ', n))
!       buffstart = buffstop+1
!     else if (buffstart < min(length, chunk_size)) then
!       call output%push(buffer(buffstart:min(length, chunk_size)))
!       buffstart = chunk_size
!     end if
!   end do
! end do
!end subroutine write_indent

subroutine forbid_reserved_label(label)
 character(len=*),intent(in) :: label
 integer :: i

 do i=1,size(reserved_keywords)
   if (reserved_keywords(i) == label) then
     MSG_ERROR(trim(label)//' is a reserved keyword and cannot be used as a YAML label.')
   end if
 end do
end subroutine forbid_reserved_label

pure function yaml_quote_string(string) result(quoted)
 character(len=*),intent(in) :: string
 character(len=len(string)+2) :: quoted
 logical :: multiline, spec_char, quote
 spec_char = index(string, ':', back=.true.) /= 0
 spec_char = spec_char .or. index(string, '{') /= 0
 spec_char = spec_char .or. index(string, '}') /= 0
 spec_char = spec_char .or. index(string, '[') /= 0
 spec_char = spec_char .or. index(string, ']') /= 0
 spec_char = spec_char .or. index(string, ',') /= 0
 spec_char = spec_char .or. index(string, '&') /= 0
 spec_char = spec_char .or. index(string, '*') /= 0
 spec_char = spec_char .or. index(string, '#') /= 0
 spec_char = spec_char .or. index(string, '?') /= 0
 spec_char = spec_char .or. index(string, '|') /= 0
 spec_char = spec_char .or. index(string, '-') /= 0
 spec_char = spec_char .or. index(string, '<') /= 0
 spec_char = spec_char .or. index(string, '>') /= 0
 spec_char = spec_char .or. index(string, '=') /= 0
 spec_char = spec_char .or. index(string, '!') /= 0
 spec_char = spec_char .or. index(string, '%') /= 0
 spec_char = spec_char .or. index(string, '@') /= 0
 spec_char = spec_char .or. index(string, '`') /= 0

 quote = index(string, "'") /= 0
 multiline = index(string, eol, back=.true.) /= 0

 if (quote) then
   quoted='"'//string//'"'
 else if (multiline .or. spec_char) then
   quoted="'"//string//"'"
 else
   quoted=string
 endif
end function yaml_quote_string

subroutine yaml_start_field(stream, label, tag, width)

 type(stream_string),intent(inout) :: stream
 character(len=*),intent(in) :: label
 integer,optional,intent(in) :: width
 character(len=*),intent(in),optional :: tag

 character(len=len_trim(label)+2) :: quoted

!#ifdef HAVE_DEBUG_MODE
 call forbid_reserved_label(trim(label))
!#endif

 quoted = yaml_quote_string(label)
 if (present(width)) then
   if (width > len_trim(label)) then
     call stream%push(trim(quoted)//repeat(' ', width-len_trim(quoted))//':')
   else
     call stream%push(trim(quoted)//':')
   end if
 else
   call stream%push(trim(quoted)//':')
 end if
 if (present(tag)) call stream%push(' !'//trim(tag))

end subroutine yaml_start_field

subroutine yaml_print_real1d(stream, length, arr, rfmt, vmax, string)

!Arguments ------------------------------------
 type(stream_string),intent(inout) :: stream
 integer,intent(in) :: vmax, length
 real(dp),intent(in) :: arr(length)
 character(len=*),intent(in) :: rfmt
 character(len=*),optional,intent(in) :: string

!Local variables-------------------------------
 integer :: i
 character(len=50) :: tmp_r
! *************************************************************************

 if (length > vmax) then
   call stream%push(' ['//eol//'    ')
 else
   call stream%push(' [')
 end if

 do i=1,length
   call string_clear(tmp_r)
   call format_real(arr(i), tmp_r, rfmt)
   call stream%push(trim(tmp_r))
   if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
     call stream%push(', '//eol//'    ')
   else
     call stream%push(', ')
   end if
 end do

 if (length > vmax) call stream%push(eol)
 if (present(string)) call stream%push(trim(string))
 call stream%push(']')

end subroutine yaml_print_real1d

subroutine yaml_print_int1d(stream, length, arr, ifmt, vmax, string)

!Arguments ------------------------------------
 type(stream_string),intent(inout) :: stream
 integer,intent(in) :: vmax
 integer,intent(in) :: length
 integer,intent(in) :: arr(length)
 character(len=*),intent(in) :: ifmt
 character(len=*),optional,intent(in) :: string

!Local variables-------------------------------
 integer :: i
 character(len=50) :: tmp_i
! *************************************************************************

 if (length > vmax) then
   call stream%push(' ['//eol//'    ')
 else
   call stream%push(' [')
 end if

 do i=1,length
   call string_clear(tmp_i)
   write(tmp_i, ifmt) arr(i)
   call stream%push(trim(tmp_i))
   if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
     call stream%push(', '//eol//'    ')
   else
     call stream%push(', ')
   end if
 end do

 if (length > vmax) call stream%push(eol)
 if (present(string)) call stream%push(trim(string))
 call stream%push(']')

end subroutine yaml_print_int1d

subroutine yaml_print_dict(stream, pl, key_size, s_size, kfmt, ifmt, rfmt, sfmt, vmax)

 type(stream_string),intent(inout) :: stream
 integer,intent(in) :: vmax
 type(pair_list),intent(inout) :: pl
 character(len=*),intent(in) :: ifmt, rfmt, kfmt, sfmt
 integer,intent(in) :: key_size, s_size
 character(len=key_size) :: key
 character(len=key_size+5) :: tmp_key
 character(len=100) :: tmp_r, tmp_i
 character(len=s_size) :: tmp_s

 integer :: i, vi, type_code
 real(dp) :: vr
 character(len=s_size) :: vs

 if (pl%length() > vmax) then
   call stream%push(' {'//eol//'    ')
 else
   call stream%push(' {')
 end if

 call pl%restart()
 do i=1,pl%length()
   call pl%iter(key, type_code, vi, vr, vs)

!#ifdef HAVE_DEBUG_MODE
   call forbid_reserved_label(trim(key))
!#endif

   ! TODO: Should enclose key in double quotation markers only if needed
   !if has_whitespaces(key) then
   !  call string_clear(tmp_key)
   !  write(tmp_key, kfmt) '"'//trim(key)//'"'
   !end if

   write(tmp_key, kfmt) trim(key)
   call stream%push(trim(tmp_key)//': ')

   select case (type_code)
   case (TC_INT)
     call string_clear(tmp_i)
     write(tmp_i, ifmt) vi
     call stream%push(trim(tmp_i))
   case (TC_REAL)
     call string_clear(tmp_r)
     call format_real(vr, tmp_r, rfmt)
     call stream%push(trim(tmp_r))
   case (TC_STRING)
     call string_clear(tmp_s)
     write(tmp_s, sfmt) vs
     call yaml_print_string(stream, trim(tmp_s))
   case default
     MSG_ERROR(sjoin("Invalid type_code:", itoa(type_code)))
   end select

   if (i > 0 .and. mod(i, vmax) == 0 .and. i /= pl%length()) then
     call stream%push(', '//eol//'    ')
   else
     call stream%push(', ')
   end if
 end do

 if (pl%length() > vmax) call stream%push(eol)
 call stream%push('}')

end subroutine yaml_print_dict

subroutine yaml_print_string(stream, string)

 type(stream_string),intent(inout) :: stream
 character(len=*),intent(in) :: string
 character(len=len_trim(string)+2) :: quoted

 quoted = yaml_quote_string(string)
 call stream%push(trim(quoted))

end subroutine yaml_print_string

!pure logical function has_whitespaces(string) result (ans)
!
! character(len=*),intent(in) :: string
! integer :: ii, jj
!
! ans = .False.
! do ii=len_trim(string), 1, -1
!   if (string(ii:ii) == " ") then
!     ans = .True.; exit
!   end if
! end do
!
! if (ans) then
!   do jj=1,ii-1
!     if (string(ii:ii) /= " ") exit
!   end do
!   if (jj == ii) ans = .False.
! end do
!
!end function has_whitespaces

end module m_yaml
!!***

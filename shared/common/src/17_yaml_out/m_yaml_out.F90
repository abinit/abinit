!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_yaml_out
!! NAME
!!  m_yaml_out
!!
!! FUNCTION
!!  This module defines low-level routines to format data into YAML documents.
!!  Supported data include numeric arrays of one and two dimensions,
!!  strings, numbers, dictionaries from m_pair_list and 1D arrays of dictionaries.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (TC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!   m_neat
!!
!! CHILDREN
!!   m_stream_string, m_pair_list
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#define SET_DEFAULT(v, optv, defv) v = defv; if(present(optv)) v = optv
#define ERROR_NO_OUT MSG_ERROR("No output medium has been provided.")

module m_yaml_out

 use defs_basis
 use ieee_arithmetic
 use m_errors
 use m_pair_list
 use m_stream_string

 implicit none

 private
!!***

 public :: yaml_open_doc, yaml_close_doc, yaml_single_dict, yaml_iterstart
 public :: yaml_add_realfield, yaml_add_intfield, yaml_add_stringfield
 public :: yaml_add_real1d, yaml_add_real2d
 public :: yaml_add_dict,  yaml_add_dictlist
 public :: yaml_add_int1d, yaml_add_int2d
 public :: yaml_add_tabular
 public :: yaml_open_tabular, yaml_add_tabular_line

 character(len=1),parameter :: eol=char(10)
 character(len=11),parameter :: default_rfmt='(ES23.15E3)'
 character(len=4),parameter :: default_ifmt='(I8)'
 character(len=13),parameter :: default_kfmt="(A)"
 character(len=13),parameter :: default_sfmt="(A)"
 integer,parameter :: default_keysize=30
 integer,parameter :: default_stringsize=500

 ! This is a list of reserved_keywords that shall not be used as keys in Yaml dictionaries.
 character(len=12),parameter :: reserved_keywords(10) = [character(len=12) :: &
   "tol_abs", "tol_rel", "tol_vec", "tol_eq", "ignore", &
   "ceil", "equation", "equations", "callback", "callbacks"]

contains

!!****f* m_yaml_out/yaml_iterstart
!! NAME
!! yaml_iterstart
!!
!! FUNCTION
!!  Mark the start of an iteration named by label and numbered by file
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label=key name
!!  val=value
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_iterstart(label, val, file_d, string, stream, newline, width)

 integer,intent(in) :: val
 character(len=*),intent(in) :: label
 integer,intent(in), optional :: file_d
 type(stream_string),intent(inout),optional :: stream
 character(len=*),intent(out),optional :: string
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 character(len=6) :: tmp_i
 logical :: nl

 SET_DEFAULT(nl, newline, .true.)
 SET_DEFAULT(w, width, 0)
 write(tmp_i, '(I6)') val

 if(present(stream)) then
   call stream%write('--- !IterStart'//eol//label//':'//tmp_i//eol//'...')
   if(nl) call stream%write(eol)
 else if(present(string)) then
   if(nl) then
     write(string, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'//eol
   else
     write(string, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
   end if
 else if(present(file_d)) then
   if(nl) then
     write(file_d, '(A)') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
   else
     write(file_d, '(A)', advance='no') '--- !IterStart'//eol//label//':'//tmp_i//eol//'...'
   end if
 else
   ERROR_NO_OUT
 end if

end subroutine yaml_iterstart
!!***

!!****f* m_yaml_out/yaml_open_doc
!! NAME
!! yaml_open_doc
!!
!! FUNCTION
!!  Open a yaml document
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  tag: add a tag to the field
!!  comment: string with comment.
!!  [newline]: optional, set to false to prevent adding newlines after fields
!!  [width]: optional, impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_open_doc(tag, comment, file_d, string, stream, newline, width)

 character(len=*),intent(in) :: tag
 character(len=*),intent(in) :: comment
 integer,intent(in), optional :: file_d
 type(stream_string),intent(inout),optional :: stream
 character(len=*),intent(out),optional :: string
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 type(stream_string) :: interm
 logical :: nl

 SET_DEFAULT(nl, newline, .true.)
 SET_DEFAULT(w, width, 0)

 call interm%write('---'//' !'//trim(tag))

 if (comment /= '') then
   call interm%write(eol//'comment')
   if(present(width)) then
     if(width > 7) then
       call interm%write(repeat(' ', width - 7))
     end if
   end if
   call interm%write(': ')
   call yaml_print_string(interm, comment, 70)
 end if
 if(nl) call interm%write(eol)

 if(present(stream)) then
   call interm%transfer(stream)
 else if(present(string)) then
   call interm%to_string(string)
 else if(present(file_d)) then
   call interm%to_file(file_d)
 else
   ERROR_NO_OUT
 end if

end subroutine yaml_open_doc
!!***

!!****f* m_yaml_out/yaml_add_realfield
!! NAME
!! yaml_add_realfield
!!
!! FUNCTION
!!  Add a real number field to a document
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [real_fmt] = override the default formating
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_realfield(label, val, file_d, string, stream, tag, real_fmt, newline, width)

  real(kind=dp),intent(in) :: val
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag, real_fmt
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

!Local variables-------------------------------
  integer :: w
  character(len=50) :: tmp_r
  type(stream_string) :: interm
  character(len=30) :: rfmt
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  rfmt = '                              '
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  call interm%write(' ')
  call format_real(val, tmp_r, trim(rfmt))
  call interm%write(trim(tmp_r))
  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_realfield
!!***

!!****f* m_yaml_out/yaml_add_intfield
!! NAME
!! yaml_add_intfield
!!
!! FUNCTION
!!  Add an integer field to a document
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [int_fmt] = optional  override the default formating
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_intfield(label, val, file_d, string, stream, tag, int_fmt, newline, width)

 integer,intent(in) :: val
 character(len=*),intent(in) :: label
 character(len=*),intent(in),optional :: tag, int_fmt
 integer,intent(in), optional :: file_d
 type(stream_string),intent(inout),optional :: stream
 character(len=*),intent(out),optional :: string
 logical,intent(in),optional :: newline
 integer,intent(in),optional :: width

!Local variables-------------------------------
 integer :: w
 character(50) :: tmp_i
 character(len=30) :: ifmt
 type(stream_string) :: interm
 logical :: nl

 SET_DEFAULT(nl, newline, .true.)
 SET_DEFAULT(w, width, 0)

 ifmt = '                              '
 SET_DEFAULT(ifmt, int_fmt, default_ifmt)
 if(present(tag)) then
   call yaml_start_field(interm, label, width=w, tag=tag)
 else
   call yaml_start_field(interm, label, width=w)
 end if

 call interm%write(' ')
 write(tmp_i, trim(ifmt)) val
 call interm%write(trim(tmp_i))
 if(nl) call interm%write(eol)

 if(present(stream)) then
   call interm%transfer(stream)
 else if(present(string)) then
   call interm%to_string(string)
 else if(present(file_d)) then
   call interm%to_file(file_d)
 else
   ERROR_NO_OUT
 end if

end subroutine yaml_add_intfield
!!***

!!****f* m_yaml_out/yaml_add_stringfield
!! NAME
!! yaml_add_stringfield
!!
!! FUNCTION
!!  Add a string field to a document
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  val = value
!!  [tag] = optional, add a tag to the field
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_stringfield(label, val, file_d, string, stream, tag, newline, width)

  character(len=*),intent(in) :: val
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  integer :: w
  type(stream_string) :: interm
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  call interm%write(' ')
  call yaml_print_string(interm, val, 70)
  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_stringfield
!!***

!!****f* m_yaml_out/yaml_add_real1d
!! NAME
!! yaml_add_real1d
!!
!! FUNCTION
!!  Add a field containing a 1D array of real numbers
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  length = size of arr
!!  arr(length) <real(kind=dp)>=
!!  [multiline_trig] = optional minimun number of elements before switching to multline representation
!!  [tag] = optional, add a tag to the field
!!  [real_fmt] = override the default formating
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_real1d(label, length, arr, file_d, string, stream, tag, real_fmt, multiline_trig, newline, width)

  integer,intent(in) :: length
  integer,intent(in),optional :: multiline_trig
  real(kind=dp),intent(in) :: arr(length)
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag, real_fmt
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  integer :: w
  character(len=30) :: rfmt
  integer :: vmax
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  rfmt = '                              '
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  SET_DEFAULT(vmax, multiline_trig, 5)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  call yaml_print_real1d(interm, length, arr, trim(rfmt), vmax)
  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_real1d
!!***

!!****f* m_yaml_out/yaml_add_int1d
!! NAME
!! yaml_add_int1d
!!
!! FUNCTION
!!  Add a field containing a 1D integer array
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  length = size of arr
!!  arr(length) <integer>=
!!  [multiline_trig] = optional minimun number of elements before switching to multline representation
!!  [tag] : add a tag to the field
!!  int_fmt <character(len=*)>=optional  override the default formating
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_int1d(label, length, arr, file_d, string, stream, tag, int_fmt, multiline_trig, newline, width)

  integer,intent(in) :: length
  integer,intent(in),optional :: multiline_trig
  integer,intent(in) :: arr(length)
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag, int_fmt
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  character(len=30) :: ifmt
  integer :: w
  integer :: vmax
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  ifmt = '                              '
  SET_DEFAULT(ifmt, int_fmt, default_ifmt)
  SET_DEFAULT(vmax, multiline_trig, 5)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  call yaml_print_int1d(interm, length, arr, trim(ifmt), vmax)
  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_int1d
!!***

!!****f* m_yaml_out/yaml_add_dict
!! NAME
!! yaml_add_dict
!!
!! FUNCTION
!!  Add a field containing a dictionary/pair_list
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label <character(len=*)>=
!!  pl <type(pair_list)>=
!!  string_size <integer>=optional maximum storage size for strings found in a pair_list
!!  key_size <integer>=optional maximum storage size for keys of a pair_list
!!  multiline_trig <integer>=optional minimun number of elements before switching to multline representation
!!  tag <character(len=*)>=optional  add a tag to the field
!!  key_fmt <character(len=*)>=optional  override the default formating
!!  int_fmt <character(len=*)>=optional  override the default formating
!!  real_fmt <character(len=*)>=optional  override the default formating
!!  string_fmt <character(len=*)>=optional  override the default formating
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!  width <integer>=optional impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  pl <type(pair_list)>=
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_dict(label, pl, file_d, string, stream, tag, key_size, string_size, key_fmt, &
        int_fmt, real_fmt, string_fmt, multiline_trig, newline, width)

  type(pair_list),intent(inout) :: pl
  character(len=*),intent(in) :: label
  integer,intent(in),optional :: string_size, key_size, multiline_trig
  character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  integer :: w
  character(len=30) :: kfmt, ifmt, rfmt, sfmt
  integer :: vmax, ks, ss
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)
  SET_DEFAULT(ks, key_size, default_keysize)
  SET_DEFAULT(ss, string_size, default_stringsize)

  kfmt = '                              '
  rfmt = '                              '
  ifmt = '                              '
  sfmt = '                              '
  SET_DEFAULT(kfmt, key_fmt, default_kfmt)
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  SET_DEFAULT(ifmt, int_fmt, default_ifmt)
  SET_DEFAULT(sfmt, string_fmt, default_sfmt)
  SET_DEFAULT(vmax, multiline_trig, 5)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  call yaml_print_dict(interm, pl, ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_dict
!!***

!!****f* m_yaml_out/yaml_add_real2d
!! NAME
!! yaml_add_real2d
!!
!! FUNCTION
!!  Add a field containing a 2D real number array
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  m <integer>=
!!  n <integer>=
!!  arr(m, n) <real(kind=dp)>=
!!  [tag]= add a tag to the field
!!  [real_fmt]: override the default formating
!!  [multiline_trig] <integer>=optional minimun number of elements before switching to multline representation
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_real2d(label, m, n, arr, file_d, string, stream, tag, real_fmt, multiline_trig, newline, width)
  integer,intent(in) :: m, n
  real(kind=dp),intent(in) :: arr(m, n)
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag, real_fmt
  integer,intent(in),optional :: multiline_trig
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  integer :: w
  integer :: i, vmax
  real(kind=dp) :: line(n)
  character(len=30) :: rfmt
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  rfmt = '                              '
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  SET_DEFAULT(vmax, multiline_trig, 5)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  do i=1,m
    call interm%write(eol//'-')
    line = arr(i,:)
    call yaml_print_real1d(interm, n, line, rfmt, vmax)
  end do

  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_real2d
!!***

!!****f* m_yaml_out/yaml_add_int2d
!! NAME
!! yaml_add_int2d
!!
!! FUNCTION
!!  Add a field containing a 2D integer array
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  m <integer>=
!!  n <integer>=
!!  arr(m, n) <integer>=
!!  [tag]= add a tag to the field
!!  [int_fmt]: override the default formating
!!  multiline_trig <integer>=optional minimun number of elements before switching to multline representation
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_int2d(interm, label, m, n, arr, file_d, string, stream, tag, int_fmt, multiline_trig, newline, width)

  integer,intent(in) :: m, n
  integer,intent(in) :: arr(m, n)
  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag, int_fmt
  integer,intent(in),optional :: multiline_trig
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  integer :: w
  integer :: i, vmax
  integer :: line(n)
  character(len=30) :: ifmt
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  ifmt = '                              '
  SET_DEFAULT(ifmt, int_fmt, default_ifmt)
  SET_DEFAULT(vmax, multiline_trig, 5)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if

  do i=1,m
    call interm%write(eol//'-')
    line = arr(i,:)
    call yaml_print_int1d(interm, n, line, ifmt, vmax)
  end do

  if(nl) call interm%write(eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_int2d
!!***

!!****f* m_yaml_out/yaml_add_dictlist
!! NAME
!! yaml_add_dictlist
!!
!! FUNCTION
!!  Add a field containing a list of dictionaries/array of pair_list
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  label = key name
!!  n <integer>=
!!  plarr(n) <type(pair_list)>=
!!  key_size <integer>=optional maximum storage size for keys of a pair_list
!!  string_size <integer>=optional maximum storage size for strings of a pair_list
!!  multiline_trig <integer>=optional minimun number of elements before switching to multline representation
!!  [tag]= add a tag to the field
!!  key_fmt <character(len=*)>=optional  override the default formating
!!  int_fmt <character(len=*)>=optional  override the default formating
!!  real_fmt <character(len=*)>=optional  override the default formating
!!  string_fmt <character(len=*)>=optional  override the default formating
!!  [newline] = set to false to prevent adding newlines after fields
!!  [width] = impose a minimum width of the field name side of the column (padding with spaces)
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_dictlist(label, n, plarr, file_d, string, stream, tag, key_size, string_size, key_fmt, int_fmt, &
                             real_fmt, string_fmt, multiline_trig, newline, width)
  integer,intent(in) :: n
  type(pair_list),intent(inout) :: plarr(n)
  character(len=*),intent(in) :: label
  integer,intent(in),optional :: key_size, string_size
  integer,intent(in),optional :: multiline_trig
  character(len=*),intent(in),optional :: tag, key_fmt, int_fmt, real_fmt, string_fmt
  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: width

  type(stream_string) :: interm
  integer :: w
  character(len=30) :: kfmt, ifmt, rfmt, sfmt
  integer :: vmax, ks, i, ss
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(w, width, 0)

  kfmt = '                              '
  rfmt = '                              '
  ifmt = '                              '
  SET_DEFAULT(kfmt, key_fmt, default_kfmt)
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  SET_DEFAULT(ifmt, int_fmt, default_ifmt)
  SET_DEFAULT(sfmt, string_fmt, default_sfmt)
  SET_DEFAULT(vmax, multiline_trig, 5)
  SET_DEFAULT(ks, key_size, default_keysize)
  SET_DEFAULT(ss, string_size, default_keysize)

  if(present(tag)) then
    call yaml_start_field(interm, label, width=w, tag=tag)
  else
    call yaml_start_field(interm, label, width=w)
  end if
  call interm%write(eol)

  do i=1,n
    call interm%write('- ')
    call yaml_print_dict(interm, plarr(i), ks, ss, trim(kfmt), trim(ifmt), trim(rfmt), trim(sfmt), vmax)
    if(nl .or. i/=n) then
      call interm%write(eol)
    end if
  end do


  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if
end subroutine yaml_add_dictlist
!!***

!!****f* m_yaml_out/yaml_open_tabular
!! NAME
!! yaml_open_tabular
!!
!! FUNCTION
!!  Open a field for tabular data
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  label = key name
!!  tag <character(len=*)>=optional  add a tag to the field
!!  [newline] = set to false to prevent adding newlines after fields
!!  [indent] = optional number of spaces to add to the header
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_open_tabular(label, tag, file_d, string, stream, indent, newline)

  character(len=*),intent(in) :: label
  character(len=*),intent(in),optional :: tag
  integer,intent(in),optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: indent

  integer :: n
  type(stream_string) :: interm
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(n, indent, 4)

  if(n > 4) then
    call interm%write(repeat(' ', n-4))
  end if

  if(present(tag)) then
    call yaml_start_field(interm, label, tag=tag)
  else
    call yaml_start_field(interm, label, tag='Tabular')
  end if
  call interm%write(' |'//eol)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_open_tabular
!!***

!!****f* m_yaml_out/yaml_add_tabular_line
!! NAME
!! yaml_add_tabular_line
!!
!! FUNCTION
!!  Add a line of tabular data in an already opened table field
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  line <character(len=*)>=
!!  [newline] = set to false to prevent adding newlines after fields
!!  [indent] = optional number of spaces to add to the header
!!
!! OUTPUT
!!  [file_d] = optional output file descriptor
!!  [string] = optional output string.
!!  [stream] = optional output stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_tabular_line(line, file_d, string, stream, newline, indent)

  character(len=*),intent(in) :: line
  integer,intent(in),optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: indent

  integer :: n
  type(stream_string) :: interm
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(n, indent, 4)

  call interm%write(repeat(' ', n)//trim(line))

  if(nl) then
    call interm%write(eol)
  end if

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if
end subroutine yaml_add_tabular_line
!!***

!!****f* m_yaml_out/yaml_add_tabular
!! NAME
!! yaml_add_tabular
!!
!! FUNCTION
!!  Add a field with a complete table data
!!  One and only one of file_d, stream or string have to be provided as the output destination.
!!
!! INPUTS
!!  label <character(len=*)>=
!!  input <type(stream_string)>=stream containing an already built table
!!  tag <character(len=*)>=optional  add a tag to the field
!!  file_d <integer>=optional output file descriptor
!!  stream <type(stream_string)>=optional output stream
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!  indent <integer>=optional number of spaces to add to each line
!!
!! OUTPUT
!!  input <type(stream_string)>=
!!  stream <type(stream_string)>=optional
!!  string <character(len=*)>=optional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_add_tabular(label, input, tag, file_d, string, stream, newline, indent)

  character(len=*),intent(in) :: label
  type(stream_string),intent(inout) :: input
  character(len=*),intent(in),optional :: tag
  integer,intent(in),optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  integer,intent(in),optional :: indent

  integer :: n
  type(stream_string) :: interm
  character(len=100) :: t
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)
  SET_DEFAULT(n, indent, 4)
  SET_DEFAULT(t, tag, 'Tabular')

  call yaml_open_tabular(label, tag=t, stream=interm, newline=nl)

  if(n > 4) then
    call interm%write(repeat(' ', n - 4))
  end if

  call write_indent(input, interm, n)

  if(nl) then
    call interm%write(eol)
  end if

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_add_tabular
!!***

!!****f* m_yaml_out/yaml_single_dict
!! NAME
!! yaml_single_dict
!!
!! FUNCTION
!!  Create a full document from a single dictionary
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  tag <character(len=*)>=
!!  comment <character(len=*)>=
!!  pl <type(pair_list)>=
!!  key_size <integer>=maximum storage size for the keys of pl
!!  string_size <integer>=maximum storage size for the strings found in pl
!!  tag <character(len=*)>=optional  add a tag to the field
!!  int_fmt <character(len=*)>=optional  override the default formating
!!  real_fmt <character(len=*)>=optional  override the default formating
!!  string_fmt <character(len=*)>=optional  override the default formating
!!  file_d <integer>=optional output file descriptor
!!  width <integer>=optional impose a minimum width of the field name side of the column (padding with spaces)
!!  stream <type(stream_string)>=optional output stream
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!
!! OUTPUT
!!  pl <type(pair_list)>=
!!  stream <type(stream_string)>=optional
!!  string <character(len=*)>=optional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_single_dict(tag, comment, pl, key_size, string_size, file_d, string, stream, &
                           int_fmt, real_fmt, string_fmt, newline, width)

  type(pair_list),intent(inout) :: pl
  character(len=*),intent(in) :: tag
  character(len=*),intent(in) :: comment
  integer,intent(in) :: key_size, string_size
  character(len=*),intent(in),optional :: int_fmt, real_fmt, string_fmt
  integer,intent(in), optional :: file_d, width
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline

  type(stream_string) :: interm
  character(len=30) :: kfmt, ifmt, rfmt, sfmt
  character(len=string_size) :: vs, tmp_s
  character(len=key_size) :: key
  integer :: vi, k, type_code, w
  character(len=50) :: tmp_i, tmp_r
  real(kind=dp) :: vr
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)

  kfmt = '                              '
  rfmt = '                              '
  ifmt = '                              '
  SET_DEFAULT(rfmt, real_fmt, default_rfmt)
  SET_DEFAULT(ifmt, int_fmt, default_ifmt)
  SET_DEFAULT(sfmt, string_fmt, default_sfmt)
  SET_DEFAULT(w, width, 0)

  call interm%write('--- !'//tag)

  if (comment /= '') then
    call interm%write(eol)
    call yaml_start_field(interm, 'comment', width=width)
    call yaml_print_string(interm, comment, 70)
  end if
  call interm%write(eol)

  call pl%restart()
  do k=1,pl%length()
    call string_clear(key)
    call string_clear(vs)
    call pl%iter(key, type_code, vi, vr, vs)

    call yaml_start_field(interm, trim(key), width=width)
    call interm%write(' ')
    if(type_code == TC_INT) then
      call string_clear(tmp_i)
      write(tmp_i, ifmt) vi
      call interm%write(trim(tmp_i))
    else if(type_code == TC_REAL) then
      call string_clear(tmp_r)
      call format_real(vr, tmp_r, rfmt)
      call interm%write(trim(tmp_r))
    else if(type_code == TC_STRING) then
      call string_clear(tmp_s)
      write(tmp_s, sfmt) vs
      call yaml_print_string(interm, trim(tmp_s), 100)
    end if
    call interm%write(eol)
  end do

  call yaml_close_doc(stream=interm, newline=nl)

  if(present(stream)) then
    call interm%transfer(stream)
  else if(present(string)) then
    call interm%to_string(string)
  else if(present(file_d)) then
    call interm%to_file(file_d)
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_single_dict
!!***

!!****f* m_yaml_out/yaml_close_doc
!! NAME
!! yaml_close_doc
!!
!! FUNCTION
!!  Close a previously opened document
!!  One and only one of file_d, stream or string have to be provided as
!!  the output destination.
!!
!! INPUTS
!!  file_d <integer>=optional output file descriptor
!!  stream <type(stream_string)>=optional output stream
!!  newline <logical>=optional  set to false to prevent adding newlines after fields
!!
!! OUTPUT
!!  stream <type(stream_string)>=optional
!!  string <character(len=*)>=optional
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine yaml_close_doc(file_d, string, stream, newline)

  integer,intent(in), optional :: file_d
  type(stream_string),intent(inout),optional :: stream
  character(len=*),intent(out),optional :: string
  logical,intent(in),optional :: newline
  logical :: nl

  SET_DEFAULT(nl, newline, .true.)

  if(present(stream)) then
    call stream%write('...')
    if(nl) call stream%write(eol)
  else if(present(string)) then
    if(nl) then
      write(string, '(A)') '...'//eol
    else
      write(string, '(A)') '...'
    end if
  else if(present(file_d)) then
    if(nl) then
      write(file_d, '(A)') '...'
    else
      write(file_d, '(A)', advance='no') '...'
    end if
  else
    ERROR_NO_OUT
  end if

end subroutine yaml_close_doc
!!***

! private
subroutine string_clear(string)
  character(len=*),intent(inout) :: string
  string = repeat(' ', len(string))
end subroutine string_clear

subroutine format_real(val, dest, formt)
  real(kind=dp),intent(in) :: val
  character(len=*),intent(out) :: dest
  character(len=*),intent(in) :: formt

  if(ieee_is_nan(val)) then  ! NaN
    write(dest, '(a)') '.nan'
  else if (val == MAGIC_UNDEF) then
    write(dest, '(a)') 'undef'
  else
    write(dest, trim(formt)) val
  end if
end subroutine format_real

subroutine write_indent(input, output, n)
  class(stream_string),intent(inout) :: input, output
  integer,intent(in) :: n

  integer :: buffstart, buffstop, length
  character(len=chunk_size) :: buffer

  do while (input%length > 0)
    length = input%length
    call input%get_chunk(buffer)

    buffstart = 1
    buffstop = 1
    do while (buffstart < min(length, chunk_size))
      buffstop = index(buffer(buffstart:), eol)
      if (buffstop > 0) then
        call output%write(buffer(buffstart:buffstop))
        call output%write(repeat(' ', n))
        buffstart = buffstop+1
      else if(buffstart < min(length, chunk_size)) then
        call output%write(buffer(buffstart:min(length, chunk_size)))
        buffstart = chunk_size
      end if
    end do
  end do
end subroutine write_indent

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
  else if(multiline .or. spec_char) then
    quoted="'"//string//"'"
  else
    quoted=string
  endif
end function yaml_quote_string

subroutine yaml_start_field(stream, label, tag, width)
  type(stream_string),intent(inout) :: stream
  character(len=*),intent(in) :: label
  integer,optional :: width
  character(len=*),intent(in),optional :: tag
  character(len=len_trim(label)+2) :: quoted

  call forbid_reserved_label(trim(label))
  quoted = yaml_quote_string(label)
  if(present(width)) then
    if(width > len_trim(label)) then
      call stream%write(trim(quoted)//repeat(' ', width-len_trim(quoted))//':')
    else
      call stream%write(trim(quoted)//':')
    end if
  else
    call stream%write(trim(quoted)//':')
  end if
  if (present(tag)) then
    call stream%write(' !'//trim(tag))
  end if
end subroutine yaml_start_field

subroutine yaml_print_real1d(stream, length, arr, rfmt, vmax)
  type(stream_string),intent(inout) :: stream
  integer,intent(in) :: vmax
  integer,intent(in) :: length
  real(kind=dp),intent(in) :: arr(length)
  character(len=*),intent(in) :: rfmt
  character(len=50) :: tmp_r

  integer :: i

  if (length > vmax) then
    call stream%write(' ['//eol//'    ')
  else
    call stream%write(' [')
  end if

  do i=1,length
    call string_clear(tmp_r)
    call format_real(arr(i), tmp_r, rfmt)
    call stream%write(trim(tmp_r))
    if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
      call stream%write(', '//eol//'    ')
    else
      call stream%write(', ')
    end if
  end do
  if (length > vmax) then
    call stream%write(eol)
  end if
  call stream%write(']')

end subroutine yaml_print_real1d

subroutine yaml_print_int1d(stream, length, arr, ifmt, vmax)

  type(stream_string),intent(inout) :: stream
  integer,intent(in) :: vmax
  integer,intent(in) :: length
  integer,intent(in) :: arr(length)
  character(len=*),intent(in) :: ifmt
  character(len=50) :: tmp_i

  integer :: i

  if (length > vmax) then
    call stream%write(' ['//eol//'    ')
  else
    call stream%write(' [')
  end if

  do i=1,length
    call string_clear(tmp_i)
    write(tmp_i, ifmt) arr(i)
    call stream%write(trim(tmp_i))
    if (i > 0 .and. mod(i, vmax) == 0 .and. i /= length) then
      call stream%write(', '//eol//'    ')
    else
      call stream%write(', ')
    end if
  end do
  if (length > vmax) then
    call stream%write(eol)
  end if
  call stream%write(']')

end subroutine yaml_print_int1d

subroutine yaml_print_dict(stream, pl, key_size, s_size, kfmt, ifmt, rfmt, sfmt, vmax)

  type(stream_string),intent(inout) :: stream
  integer,intent(in) :: vmax
  type(pair_list),intent(inout) :: pl
  character(len=*),intent(in) :: ifmt, rfmt, kfmt, sfmt
  integer,intent(in) :: key_size, s_size
  character(len=key_size) :: key
  character(len=key_size+5) :: tmp_key
  character(len=100) :: tmp_r
  character(len=100) :: tmp_i
  character(len=s_size) :: tmp_s

  integer :: i, vi, type_code
  real(kind=dp) :: vr
  character(len=s_size) :: vs

  if (pl%length() > vmax) then
    call stream%write(' {'//eol//'    ')
  else
    call stream%write(' {')
  end if

  call pl%restart()
  do i=1,pl%length()
    call pl%iter(key, type_code, vi, vr, vs)

    call forbid_reserved_label(trim(key))

    call string_clear(tmp_key)
    write(tmp_key, kfmt) '"'//trim(key)//'"'
    call stream%write(trim(tmp_key)//': ')
    if(type_code == TC_INT) then
      call string_clear(tmp_i)
      write(tmp_i, ifmt) vi
      call stream%write(trim(tmp_i))
    else if(type_code == TC_REAL) then
      call string_clear(tmp_r)
      call format_real(vr, tmp_r, rfmt)
      call stream%write(trim(tmp_r))
    else if(type_code == TC_STRING) then
      call string_clear(tmp_s)
      write(tmp_s, sfmt) vs
      call yaml_print_string(stream, trim(tmp_s), 100)
    end if
    if (i > 0 .and. mod(i, vmax) == 0 .and. i /= pl%length()) then
      call stream%write(', '//eol//'    ')
    else
      call stream%write(', ')
    end if
  end do

  if (pl%length() > vmax) then
    call stream%write(eol)
  end if
  call stream%write('}')

end subroutine yaml_print_dict

subroutine yaml_print_string(stream, string, vmax)

  type(stream_string),intent(inout) :: stream
  integer,intent(in) :: vmax
  character(len=*),intent(in) :: string
  character(len=len_trim(string)+2) :: quoted

  logical :: auto_wrap

  if(len(string) > vmax) then
    auto_wrap = .true.
  else
    auto_wrap = .false.
  end if

  quoted = yaml_quote_string(string)
  call stream%write(trim(quoted))

end subroutine yaml_print_string

end module m_yaml_out
!!***

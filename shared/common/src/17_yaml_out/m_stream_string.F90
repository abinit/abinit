!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_stream_string
!! NAME
!!  m_stream_string
!!
!! FUNCTION
!!  This module define a type representing a variable size
!!  string. It can be used in a file-like way by writing to it or reading it.
!!  Memory is automatically allocated on writing and freed on reading.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (TC, MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Provide tools to manipulate variable size strings in an incremental FIFO way
!! Use `stream%push` to incrementaly fill the string. The required memory space will be allocated
!! automatically when needed.
!! To avoid memory leaks you have to use stream_free on the stream to free the memory space unless
!! you already flushed it using stream%flush, stream%transfer, stream%to_string or stream%to_file.
!! Unlike the last four methods, stream_copy and stream_debug do not modify the source stream
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_stream_string

  use defs_basis
  use m_profiling_abi

  implicit none

  private

  integer,public,parameter :: chunk_size = 248

  type,private :: stream_chunk
    type(stream_chunk), pointer :: next => null()
    character(len=chunk_size) :: chunk = repeat(' ', chunk_size)
  end type stream_chunk

  type,public :: stream_string
    integer :: length = 0
    type(stream_chunk), pointer :: head => null()
    contains
      procedure :: flush => stream_flush_unit
      procedure :: flush_units => stream_flush_units
      procedure :: free => stream_free
      procedure :: copy => stream_copy
      procedure :: push => stream_push
      procedure :: pop_chunk => stream_pop_chunk
      procedure :: to_string => stream_to_string
      procedure :: to_file => stream_to_file
      procedure :: transfer => stream_transfer
      procedure :: debug => stream_debug

  end type stream_string

contains
!!***

subroutine stream_flush_unit(stream, unit, newline)

  class(stream_string),intent(inout) :: stream
  integer,intent(in) :: unit
  logical,optional,intent(in) :: newline

  character(len=stream%length) :: s

  if (unit == dev_null) then
    call stream%free()
    return
  end if

  call stream%to_string(s)
  write(unit, "(a)")trim(s)

  if (present(newline)) then
    if (newline) write(unit, "(a)")""
  end if

  call stream%free()

end subroutine stream_flush_unit


subroutine stream_flush_units(stream, units, newline)

  class(stream_string),intent(inout) :: stream
  integer,intent(in) :: units(:)
  logical,optional,intent(in) :: newline

!Local variables-------------------------------
!scalars
 integer :: ii, cnt
 character(len=stream%length) :: s
!arrays
 integer :: my_units(size(units))

!******************************************************************

 ! Remove duplicated units (if any)
 my_units(1) = units(1); cnt = 1
 do ii=2,size(units)
   if (any(units(ii) == my_units(1:cnt))) cycle
   cnt = cnt + 1
   my_units(cnt) = units(ii)
 end do

 call stream%to_string(s)

 do ii=1,cnt
   if (units(ii) == dev_null) cycle
   write(units(ii), "(a)")trim(s)
   if (present(newline)) then
     if (newline) write(units(ii), "(a)")""
   end if
 end do

 call stream%free()

end subroutine stream_flush_units

!!****f* m_stream_string/stream_free
!! NAME
!! stream_free
!!
!! FUNCTION
!!  free stream. Most of the time this is not needed since
!!  routines to access the content free the stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_free(stream)

  class(stream_string),intent(inout) :: stream
  type(stream_chunk), pointer :: cursor, prev
  cursor => stream%head
  do while (associated(cursor))
    prev => cursor
    cursor => cursor%next
    ABI_FREE_SCALAR(prev)
  end do
  stream%head => NULL()
  stream%length = 0

end subroutine stream_free
!!***

!!****f* m_stream_string/stream_copy
!! NAME
!! stream_copy
!!
!! FUNCTION
!!  copy src content to dest without altering src
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_copy(src, dest)
  class(stream_string),intent(inout) :: src, dest
  type(stream_chunk), pointer :: cursor
  cursor => src%head
  do while (associated(cursor))
    call dest%push(cursor%chunk)
    cursor => cursor%next
  end do
end subroutine stream_copy
!!***

!!****f* m_stream_string/stream_push
!! NAME
!! stream_push
!!
!! FUNCTION
!!  Write string to stream, allocating memory if needed
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_push(stream, string)
  class(stream_string),intent(inout) :: stream
  character(len=*),intent(in) :: string
  integer :: offset, room_left, soffset
  type(stream_chunk), pointer :: cursor

  offset = stream%length

  if (.not.associated(stream%head)) then
    ABI_MALLOC_SCALAR(stream%head)
  end if
  cursor => stream%head

  do while(offset > chunk_size)
    cursor => cursor%next
    offset = offset - chunk_size
  end do

  room_left = chunk_size - offset
  if (room_left < len(string)) then
    cursor%chunk(offset+1:chunk_size) = string(1:room_left)
    soffset = room_left
    do while (soffset < len(string))
      ABI_MALLOC_SCALAR(cursor%next)
      cursor%next%chunk(1:min(chunk_size, len(string)-soffset)) = &
        string(soffset+1:min(soffset+chunk_size,len(string)))
      cursor => cursor%next
      soffset = soffset + chunk_size
    end do
  else
    cursor%chunk(offset+1:offset+len(string)) = string
  end if
  stream%length = stream%length + len(string)

end subroutine stream_push
!!***

!!****f* m_stream_string/stream_pop_chunk
!! NAME
!! stream_pop_chunk
!!
!! FUNCTION
!!  Remove the last chunk of stream an put its content in string
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_pop_chunk(stream, string)
  class(stream_string),intent(inout) :: stream
  character(len=chunk_size),intent(out) :: string
  type(stream_chunk),pointer :: cursor

  string = stream%head%chunk
  if (stream%length > chunk_size) then
    ! copy the next pointer
    cursor => stream%head%next
    ! have next pointing to nothing
    stream%head%next => NULL()
    ! free head
    ABI_FREE_SCALAR(stream%head)
    stream%head => cursor
    stream%length = stream%length - chunk_size
  else
    ABI_FREE_SCALAR(stream%head)
    stream%length = 0
  end if

end subroutine stream_pop_chunk
!!***

!!****f* m_stream_string/stream_to_string
!! NAME
!! stream_to_string
!!
!! FUNCTION
!!  Copy the content of stream to string, freeing stream string must be large enough
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_to_string(stream, string)

  class(stream_string),intent(inout) :: stream
  character(len=*),intent(out) :: string
  character(len=chunk_size) :: stmp
  integer :: offset, length
  offset = 0

  string = repeat(' ', len(string))
  do while (stream%length > 0)
    length = stream%length
    call stream%pop_chunk(stmp)
    string(offset+1:offset+min(length, chunk_size)) = stmp(1:min(length, chunk_size))
    offset = offset + chunk_size
  end do

end subroutine stream_to_string
!!***

!!****f* m_stream_string/stream_to_file
!! NAME
!! stream_to_file
!!
!! FUNCTION
!!  Write the content of stream to the file, freeing stream
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_to_file(stream, file_d)
  class(stream_string),intent(inout) :: stream
  integer,intent(in) :: file_d
  character(len=chunk_size) :: stmp
  integer :: offset, length
  offset = 0

  do while (stream%length > 0)
    length = stream%length
    call stream%pop_chunk(stmp)
    write(file_d, '(A)', advance='no') stmp(1:min(length, chunk_size))
    offset = offset + chunk_size
  end do

end subroutine stream_to_file
!!***

!!****f* m_stream_string/stream_transfer
!! NAME
!! stream_transfer
!!
!! FUNCTION
!!  Copy the content of src to dest, freeing src
!!  If possible does not reallocate memory and just have
!!  dest point to src content
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_transfer(src, dest)
  class(stream_string),intent(inout) :: src, dest
  character(len=chunk_size) :: chunk
  integer :: length
  if(.not.associated(dest%head)) then
    ! if possible just transfer the pointer
    dest%head => src%head
    dest%length = src%length
    src%head => NULL()
  else
    do while (src%length > 0)
      length = src%length
      call src%pop_chunk(chunk)
      if(length > chunk_size) then
        call dest%push(chunk)
      else
        call dest%push(chunk(1:length))
      end if
    end do
  end if

end subroutine stream_transfer
!!***

!!****f* m_stream_string/stream_debug
!! NAME
!! stream_debug
!!
!! FUNCTION
!!  Show the content of the chunks on stdout
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_debug(src)
  class(stream_string),intent(inout) :: src
  type(stream_chunk), pointer :: cursor
  integer :: c
  cursor => src%head
  c = 1
  do while (associated(cursor))
    write(std_out,*) "Chunk no", c
    write(std_out,'(A)') cursor%chunk
    cursor => cursor%next
    c = c + 1
  end do
end subroutine stream_debug
!!***

end module m_stream_string
!!***

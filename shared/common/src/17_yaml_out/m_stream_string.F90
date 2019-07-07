!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_stream_string
!! NAME
!!  m_stream_string
!!
!! FUNCTION
!!  This module define a type representing a variable size
!!  string. It can be used in a file-like way by writing to it
!!  or reading it.
!!  Memory is automatically allocated on writing and freed on reading.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Provide tools to manipulate variable size strings in an incremental FIFO way
!! Use write to incrementaly fill the string. The required memory space will be allocated
!! automatically when needed.
!! To avoid memory leaks you have to use stream_free on the stream to free the memory space unless
!! you already flushed it using stream_transfer, stream_to_string or stream_to_file.
!! Unlike the latter three methods, stream_copy and stream_debug do not modify the source stream
!!
!! PARENTS
!!   m_yaml_out, m_neat
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
  use m_errors
  use m_profiling_abi

  implicit none

  private

  public :: enable_yaml

  integer,public,parameter :: chunk_size = 248

  type,private :: stream_chunk
    type(stream_chunk), pointer :: next => null()
    character(len=chunk_size) :: chunk = repeat(' ', chunk_size)
  end type stream_chunk

  type,public :: stream_string
    integer :: length = 0
    type(stream_chunk), pointer :: head => null()
    contains
      procedure :: dump => wrtout_stream
      procedure :: free => stream_free
      procedure :: copy => stream_copy
      procedure :: write => stream_write
      procedure :: get_chunk => stream_get_chunk
      procedure :: to_string => stream_to_string
      procedure :: to_file => stream_to_file
      procedure :: transfer => stream_transfer
      procedure :: debug => stream_debug
  end type stream_string

  logical,private,save :: enable = .false.

contains
!!***

! Can only be called once, then it lock. This is to prevent unconsistent
! behaviour with activating and disabling on demand (the parser on the python
! side won't like it)
subroutine enable_yaml(yes)
  logical, intent(in) :: yes
  enable = yes
end subroutine enable_yaml

subroutine wrtout_stream(stream, unit, newline)
  class(stream_string),intent(inout) :: stream
  integer,intent(in) :: unit
  logical,optional,intent(in) :: newline

  character(len=stream%length) :: s

  if (enable) then
    call stream%to_string(s)
    if (present(newline)) then
      if (newline) write(unit, "(a)")""
    end if
    write(unit, "(a)")trim(s)
  endif

  call stream%free()
end subroutine wrtout_stream


!!****f* m_stream_string/stream_free
!! NAME
!! stream_free
!!
!! FUNCTION
!!  free stream. Most of the time this is not needed since
!!  routines to access the content free the stream
!!
!! INPUTS
!!  stream <class(stream_string)>=
!!
!! OUTPUT
!!  stream <class(stream_string)>=
!!
!! NOTES
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
!! INPUTS
!!  src <class(stream_string)>=
!!  dest <class(stream_string)>=
!!
!! OUTPUT
!!  src <class(stream_string)>=
!!  dest <class(stream_string)>=
!!
!! NOTES
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
    call stream_write(dest, cursor%chunk)
    cursor => cursor%next
  end do
end subroutine stream_copy
!!***

!!****f* m_stream_string/stream_write
!! NAME
!! stream_write
!!
!! FUNCTION
!!  Write string to stream, allocating memory if needed
!!
!! INPUTS
!!  stream <class(stream_string)>=
!!  string <character(len=*)>=
!!
!! OUTPUT
!!  stream <class(stream_string)>=
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_write(stream, string)
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

end subroutine stream_write
!!***

!!****f* m_stream_string/stream_get_chunk
!! NAME
!! stream_get_chunk
!!
!! FUNCTION
!!  Remove the last chunk of stream an put its content in string
!!
!! INPUTS
!!  stream <class(stream_string)>=
!!
!! OUTPUT
!!  stream <class(stream_string)>=
!!  string <character(len=chunk_size)>=
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine stream_get_chunk(stream, string)
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

end subroutine stream_get_chunk
!!***

!!****f* m_stream_string/stream_to_string
!! NAME
!! stream_to_string
!!
!! FUNCTION
!!  Copy the content of stream to string, freeing stream
!!  string HAVE to be large enough
!!
!! INPUTS
!!  stream <class(stream_string)>=
!!
!! OUTPUT
!!  stream <class(stream_string)>=
!!  string <character(len=*)>=
!!
!! NOTES
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
    call stream_get_chunk(stream, stmp)
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
!! INPUTS
!!  stream <class(stream_string)>=
!!  file_d <integer>=
!!
!! OUTPUT
!!  stream <class(stream_string)>=
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
    call stream_get_chunk(stream, stmp)
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
!! INPUTS
!!  src <class(stream_string)>=
!!  dest <class(stream_string)>=
!!
!! OUTPUT
!!  src <class(stream_string)>=
!!  dest <class(stream_string)>=
!!
!! NOTES
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
      call stream_get_chunk(src, chunk)
      if(length > chunk_size) then
        call stream_write(dest, chunk)
      else
        call stream_write(dest, chunk(1:length))
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
!! INPUTS
!!  src <class(stream_string)>=
!!
!! OUTPUT
!!  src <class(stream_string)>=
!!
!! NOTES
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

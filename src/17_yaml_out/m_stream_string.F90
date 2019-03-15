#if defined HAVE_CONFIG_H
#include "config.h"
#endif

! Provide tools to manipulate variable size strings in an incermental FIFO way
! Use write to incrementaly fill the string. The required memory space will be allocated
! automatically when needed.
! To avoid memory flaws you have to use stream_free on the stream to free the memory space unless
! you already flushed it using stream_transfer, stream_to_string or stream_to_file.
! Unlike the latter three stream_copy and stream_debug do not modify the source stream
module m_stream_string
  
  implicit none

  integer,parameter :: chunk_size=248

  private
  public :: stream_write, stream_get_chunk, stream_free, stream_string
  public :: stream_to_string, stream_copy, stream_transfer, stream_to_file
  public :: stream_debug, stream_chunk

  type stream_chunk
    type(stream_chunk), pointer :: next => NULL()
    character(len=chunk_size) :: chunk = repeat(' ', chunk_size)
  end type stream_chunk

  type stream_string
    integer :: length = 0
    type(stream_chunk), pointer :: head => NULL()
    contains
      procedure :: free => stream_free
      procedure :: copy => stream_copy
      procedure :: write => stream_write
      procedure :: get_chunk => stream_get_chunk
      procedure :: to_string => stream_to_string
      procedure :: to_file => stream_to_file
      procedure :: transfer => stream_transfer
      procedure :: debug => stream_debug
  end type stream_string

  contains

  subroutine stream_free(stream)
    class(stream_string),intent(inout) :: stream
    type(stream_chunk), pointer :: cursor, prev
    cursor => stream%head
    do while (associated(cursor))
      prev => cursor
      cursor => cursor%next
      deallocate(prev)
    end do
  end subroutine stream_free

  subroutine stream_copy(src, dest)
    class(stream_string),intent(inout) :: src, dest
    type(stream_chunk), pointer :: cursor
    cursor => src%head
    do while (associated(cursor))
      call stream_write(dest, cursor%chunk)
      cursor => cursor%next
    end do
  end subroutine stream_copy

  subroutine stream_write(stream, string)
    class(stream_string),intent(inout) :: stream
    character(len=*),intent(in) :: string
    integer :: offset, room_left, soffset
    type(stream_chunk), pointer :: cursor

    offset = stream%length

    if (.not.associated(stream%head)) then
      allocate(stream%head)
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
        allocate(cursor%next)
        cursor%next%chunk(1:min(chunk_size, len(string)-soffset)) = &
&         string(soffset+1:min(soffset+chunk_size,len(string)))
        cursor => cursor%next
        soffset = soffset + chunk_size
      end do
    else
      cursor%chunk(offset+1:offset+len(string)) = string
    end if
    stream%length = stream%length + len(string)
  end subroutine stream_write

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
      deallocate(stream%head)
      stream%head => cursor
      stream%length = stream%length - chunk_size
    else
      deallocate(stream%head)
      stream%length = 0
    end if
  end subroutine stream_get_chunk

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

  subroutine stream_debug(src)
    class(stream_string),intent(inout) :: src
    type(stream_chunk), pointer :: cursor
    integer :: c
    cursor => src%head
    c = 1
    do while (associated(cursor))
      write(*,*) "Chunk no", c
      write(*,'(A)') cursor%chunk
      cursor => cursor%next
      c = c + 1
    end do
  end subroutine stream_debug

end module m_stream_string

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_hash_md5
!! NAME
!!  m_hash_md5
!!
!! FUNCTION
!!  This module provides resources to calculate MD5 checksums.
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2020 ABINIT group (Yann Pouillon)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_hash_md5

  use iso_c_binding
  use m_abicore

  implicit none

  private

  public :: md5_check             ! Checks whether two MD5 sums are identical
  public :: md5_sum_from_file     ! Computes a MD5 sum from a file
  public :: md5_sum_from_string   ! Computes a MD5 sum from a text string

  type :: md5_context_t
    private
    type(c_ptr) :: ptr = C_NULL_PTR
  end type md5_context_t

  interface
    function MD5_Context_New() bind(c, name="MD5_Context_New")
      import
      type(c_ptr) :: MD5_Context_New
    end function MD5_Context_New
  end interface
  interface
    subroutine MD5_Digest_File(fname, retval) bind(c, name="MD5_Digest_File")
      import
      character(kind=c_char) :: fname(*)
      character(kind=c_char) :: retval(33)
    end subroutine MD5_Digest_File
  end interface
  interface
    subroutine MD5_Final(retval, ctx) bind(c, name="MD5_Final")
      import
      character(kind=c_char) :: retval(33)
      type(c_ptr), value :: ctx
    end subroutine MD5_Final
  end interface
  interface
    subroutine MD5_Init(ctx) bind(c, name="MD5_Init")
      import
      type(c_ptr), value :: ctx
    end subroutine MD5_Init
  end interface
  interface
    subroutine MD5_Update(ctx, buffer, bufsize) bind(c, name="MD5_Update")
      import
      type(c_ptr), value :: ctx
      character(kind=c_char) :: buffer(*)
      integer(c_int), value :: bufsize
    end subroutine MD5_Update
  end interface

contains  !===========================================================
!!***

!!****f* m_hash_md5/md5_check
!! NAME
!! md5_check
!!
!! FUNCTION
!!  Checks whether two MD5 sums are identical.
!!
!! INPUTS
!!  sum1 = the first sum to compare
!!  sum2 = the second sum to compare
!!
!! OUTPUT
!!  boolean telling whether the two sums are identical
!!
!! NOTES
!!  Created a function to be able to add more operations than just checking
!!  the equality of the sums.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function md5_check(sum1,sum2)

!Arguments ------------------------------------
  character(len=32),intent(in) :: sum1
  character(len=32),intent(in) :: sum2

!Local variables-------------------------------
  logical :: md5_check

! *********************************************************************

  md5_check = ( sum1 == sum2 )

end function md5_check
!!***

! ---------------------------------------------------------------------

!!****f* m_hash_md5/md5_sum_from_file
!! NAME
!! md5_sum_from_file
!!
!! FUNCTION
!!  Computes a MD5 sum from a file.
!!
!! INPUTS
!!  fname = path to the file
!!
!! OUTPUT
!!  String representing the MD5 sum of the file
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function md5_sum_from_file(fname)

!Arguments ------------------------------------
  character(len=*),intent(in) :: fname

!Local variables-------------------------------
  character(len=32) :: md5_sum_from_file
  character(kind=c_char) :: retval(33)
  character(kind=c_char), allocatable :: path(:)
  integer :: strlen

! *********************************************************************

  ! Translate file name to C
  strlen = len_trim(fname)
  ABI_ALLOCATE(path,(strlen+1))
  call f_to_c_string(fname, path)

  ! Get MD5 sum from C
  call MD5_Digest_file(path, retval)
  call c_to_f_string(retval, md5_sum_from_file)

  ! Clean up the mess
  ABI_DEALLOCATE(path)

end function md5_sum_from_file
!!***

! ---------------------------------------------------------------------

!!****f* m_hash_md5/md5_sum_from_string
!! NAME
!! md5_sum_from_string
!!
!! FUNCTION
!!  Computes a MD5 sum from a string.
!!
!! INPUTS
!!  text = string to process
!!
!! OUTPUT
!!  String representing the MD5 sum of the argument
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function md5_sum_from_string(text)

!Arguments ------------------------------------
  character(len=*), intent(in) :: text

!Local variables-------------------------------
  character(len=32) :: md5_sum_from_string
  type(md5_context_t) :: ctx

! *********************************************************************

  call hash_init(ctx)
  call hash_update(ctx, text, len_trim(text))
  call hash_final(md5_sum_from_string, ctx)

end function md5_sum_from_string
!!***

! ========================= Private functions =========================

!!****f* m_hash_md5/hash_final
!! NAME
!!  hash_final
!!
!! FUNCTION
!!  Builds the final return value of the MD5 checksum.
!!
!! INPUTS
!!  ctx = MD5 context object
!!
!! OUTPUT
!!  retval = string containing the MD5 checksum
!!
!! PARENTS
!!      m_hash_md5
!!
!! CHILDREN
!!
!! SOURCE

subroutine hash_final(retval, ctx)

!Arguments ------------------------------------
  character(len=32), intent(out) :: retval
  type(md5_context_t), intent(inout) :: ctx

!Local variables-------------------------------
  character(kind=c_char) :: c_retval(33)

! *********************************************************************

  call MD5_Final(c_retval, ctx%ptr)
  call c_to_f_string(c_retval, retval)

end subroutine hash_final
!!***

! ---------------------------------------------------------------------

!!****f* m_hash_md5/hash_init
!! NAME
!! hash_init
!!
!! FUNCTION
!!  Checks whether two MD5 sums are identical.
!!
!! INPUTS
!!  ctx = MD5 context object
!!
!! SIDE EFFECTS
!!  ctx is reset to its intial values
!!
!! NOTES
!!  Created a function to be able to add more operations than just checking
!!  the equality of the sums.
!!
!! PARENTS
!!      m_hash_md5
!!
!! CHILDREN
!!
!! SOURCE

subroutine hash_init(ctx)

!Arguments ------------------------------------
  type(md5_context_t), intent(inout) :: ctx

! *********************************************************************

  ctx%ptr = MD5_Context_New()
  call MD5_Init(ctx%ptr)

end subroutine hash_init
!!***

! ---------------------------------------------------------------------

!!****f* m_hash_md5/hash_update
!! NAME
!! hash_update
!!
!! FUNCTION
!!  Updates a MD5 context object.
!!
!! INPUTS
!!  ctx = MD5 context object
!!  buffer = data to process
!!  bufsize = number of bytes to process
!!
!! SIDE EFFECTS
!!  ctx gets updated with the new data
!!
!! PARENTS
!!      m_hash_md5
!!
!! CHILDREN
!!
!! SOURCE

subroutine hash_update(ctx, buffer, bufsize)

!Arguments ------------------------------------
  type(md5_context_t), intent(inout) :: ctx
  character(len=*), intent(in) :: buffer
  integer, intent(in) :: bufsize

!Local variables-------------------------------
  character(kind=c_char), allocatable :: c_buffer(:)
  integer :: strlen

! *********************************************************************

  ! Translate buffer into C
  strlen = len_trim(buffer)
! allocate(c_buffer(strlen+1))
  ABI_ALLOCATE(c_buffer,(strlen+1))
  call f_to_c_string(trim(buffer), c_buffer)

  ! Update C MD5 context
  call MD5_Update(ctx%ptr, c_buffer, bufsize)

  ! Clean up the mess
! deallocate(c_buffer)
  ABI_DEALLOCATE(c_buffer)

end subroutine hash_update
!!***

                    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Helper functions to convert between C and Fortran strings
  ! Based on the routines by Joseph M. Krahn

  subroutine c_to_f_string(c_string, f_string)

    character(kind=c_char,len=1), intent(in)  :: c_string(*)
    character(len=*), intent(out) :: f_string

    integer :: i

    i = 1
    do while(c_string(i) /= C_NULL_CHAR .and. i <= len(f_string))
        f_string(i:i) = c_string(i)
        i = i + 1
    end do
    if (i < len(f_string)) f_string(i:) = ' '

  end subroutine c_to_f_string

  subroutine c_to_f_string_ptr(c_string, f_string)

    type(c_ptr),      intent(in)  :: c_string
    character(len=*), intent(out) :: f_string

    character(len=1, kind=c_char), pointer :: p_chars(:)
    integer :: i

    if (.not. c_associated(c_string)) then
        f_string = ' '
    else
        call c_f_pointer(c_string, p_chars, [huge(0)])
        i = 1
        do while(p_chars(i) /= C_NULL_CHAR .and. i <= len(f_string))
          f_string(i:i) = p_chars(i)
          i = i + 1
        end do
        if (i < len(f_string)) f_string(i:) = ' '
    end if

  end subroutine c_to_f_string_ptr

  subroutine f_to_c_string(f_string, c_string)

    character(len=*), intent(in) :: f_string
    character(kind=c_char,len=1), intent(out) :: c_string(len_trim(f_string)+1)

    integer :: i, strlen

    strlen = len_trim(f_string)

    forall (i=1:strlen)
        c_string(i) = f_string(i:i)
    end forall
    c_string(strlen+1) = C_NULL_CHAR

  end subroutine f_to_c_string

end module m_hash_md5
!!***

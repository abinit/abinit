!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2025 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_clib

 use, intrinsic :: iso_c_binding

 implicit none

 private

 public :: clib_rename        !  Rename a file with a new name using the rename function from C stdlib
 public :: clib_cclock
 public :: clib_etime
 public :: clib_mtrace
 public :: clib_print_mallinfo
 public :: clib_ulimit_stack    ! Set stack size limit to maximum allowed value.
 public :: clib_getpid          ! Get process id
 public :: clib_sleep           ! Sleep for a certain number of seconds.
 public :: clib_setenv          ! Set env variable.
 public :: clib_mkdir_if_needed
 public :: clib_lock_file_by_name
 public :: clib_close_fd

! ===================================================
! ==== Fortran-bindings declared in fsi_posix.c ====
! ===================================================

 interface
    subroutine c_mkdir_if_needed(dirpath, ierr) bind(C, name="c_mkdir_if_needed")
      import
      character(kind=c_char), dimension(*), intent(in) :: dirpath
      integer(c_int),intent(out) :: ierr
    end subroutine c_mkdir_if_needed
 end interface

 interface
   integer(c_int) function c_rename(oldname, newname) bind(C, name='rename')
     import
     character(kind=c_char),intent(in) :: oldname(*)
     character(kind=c_char),intent(in) :: newname(*)
   end function c_rename
 end interface

 interface
    subroutine clib_lock_file_by_name(filename, fd, ierr) bind(C, name="c_lock_file_by_name")
      import
      character(kind=c_char), dimension(*), intent(in) :: filename
      integer(c_int),intent(out) :: fd
      integer(c_int),intent(out) :: ierr
    end subroutine clib_lock_file_by_name
 end interface

 interface
    subroutine clib_close_fd(fd) bind(C, name="c_close_fd")
      import
      integer(c_int),intent(in) :: fd
    end subroutine clib_close_fd
 end interface

 interface
   subroutine clib_cclock(cpu) bind(C, name="cclock")
     import
     real(c_double),intent(out) :: cpu
   end subroutine clib_cclock
 end interface

 interface
   real(c_double) function clib_etime(tt) bind(C, name="etime") result(res)
     import
     real(c_float),intent(out) :: tt(2)
   end function clib_etime
 end interface

 interface
   ! The type of pid_t data is a signed integer type (signed int or we can say int).
   function clib_getpid() bind(C, name='getpid')
     import
     integer(c_int) :: clib_getpid
   end function clib_getpid
 end interface

 interface
   subroutine clib_sleep(seconds) bind(C, name="sleep")
     import
     integer(c_int), value :: seconds  ! This is unsigned int in C
   end subroutine clib_sleep
 end interface


! =================================================
! ==== Fortran-bindings declared in mallinfo.c ====
! =================================================
 interface
   subroutine clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks) bind(C, name="clib_mallinfo")
     import
     integer(c_long),intent(out) :: arena, hblkhd, usmblks, fsmblks, uordblks, fordblks
   end subroutine clib_mallinfo
 end interface

! ==================================================
! ==== Fortran-bindings declared in gnu_tools.c ====
! ==================================================

 interface
   subroutine clib_mtrace(ierr) bind(C, name="clib_mtrace")
     import
     integer(c_int),intent(out) :: ierr
   end subroutine
 end interface

 interface
   subroutine clib_muntrace(ierr) bind(C, name="clib_muntrace")
     import
     integer(c_int),intent(out) :: ierr
   end subroutine
 end interface

 interface
   subroutine clib_mcheck(ierr) bind(C, name="clib_mcheck")
     import
     integer(c_int),intent(out) :: ierr
   end subroutine
 end interface

 interface
   ! Set stack size limit to maximum allowed value. Return soft and hard limit and exit status.
   subroutine clib_ulimit_stack(rlim_cur, rlim_max, ierr) bind(C, name="ulimit_stack")
     import
     integer(c_long),intent(out) :: rlim_cur, rlim_max
     integer(c_int),intent(out) :: ierr
   end subroutine
 end interface

 interface
   integer(C_INT) function setenv(name, value, overwrite) bind(C, name="setenv")
     import
     character(kind=c_char),intent(in) :: name(*), value(*)
     integer(c_int),intent(in) :: overwrite
   end function
 end interface

! ==========================================
! ==== Fortran-bindings for file_lock.c ====
! ==========================================

 !interface
 !  function lock_file(filepath) bind(C)
 !    import
 !    implicit none
 !    character(kind=c_char),intent(in) :: filepath(*)
 !    integer(c_int) :: lock_file
 !  end function lock_file
 !end interface

 !interface
 !  function unlock_fd(fd) bind(C)
 !    import
 !    implicit none
 !    integer(c_int),value,intent(in) :: fd
 !    integer(c_int) unlock_fd
 !  end function unlock_fd
 !end interface

contains
!!***

!!****f* m_clib/clib_print_fmallinfo
!! NAME
!!   clib_print_fmallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine clib_print_mallinfo(unit)

!Arguments ------------------------------------
 integer,intent(in) :: unit

!Local variables-------------------------------
 integer(c_long) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
! *********************************************************************

  call clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks)

  write(unit,*)""
  write(unit,*)"--- !Mallinfo"
  write(unit,*)' Total space in arena: ',arena
  write(unit,*)' Space in holding block headers: ',hblkhd
  write(unit,*)' Space in small blocks in use: ',usmblks
  write(unit,*)' Space in free small blocks: ',fsmblks
  write(unit,*)' Space in ordinary blocks in use: ',uordblks
  write(unit,*)' Space in free ordinary blocks: ',fordblks
  write(unit,*)"..."
  write(unit,*)""

end subroutine clib_print_mallinfo
!!***

!!****f* m_clib/clib_rename
!! NAME
!!  clib_rename
!!
!! FUNCTION
!!  Rename a file with a new name using the rename function from C stdlib
!!
!! SOURCE

integer function clib_rename(old_fname, new_fname) result(ierr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: old_fname, new_fname
! *********************************************************************

 ierr = c_rename(trim(old_fname)//c_null_char, trim(new_fname)//c_null_char)

end function clib_rename
!!***

!!****f* m_clib/clib_mkdir
!! NAME
!!  clib_mkdir
!!
!! FUNCTION
!!  Create a directory if it does not exist. Return 0 on success.
!!
!! SOURCE

subroutine clib_mkdir_if_needed(dirpath, ierr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: dirpath
 integer,intent(out) :: ierr
! *********************************************************************

 call c_mkdir_if_needed(trim(dirpath)//c_null_char, ierr)

end subroutine clib_mkdir_if_needed
!!***

!!****f* m_clib/clib_setenv
!! NAME
!!  clib_setenv
!!
!! FUNCTION
!!   The setenv() function adds the variable name to the environment
!!   with the value value, if name does not already exist.  If name
!!   does exist in the environment, then its value is changed to value
!!   if overwrite is nonzero; if overwrite is zero, then the value of
!!   name is not changed (and setenv() returns a success status).
!!   This function makes copies of the strings pointed to by name and
!!   value (by contrast with putenv(3)).
!!
!! SOURCE

integer function clib_setenv(name, value, overwrite) result(ierr)

!Arguments ------------------------------------
 character(len=*) ,intent(in) :: name, value
 integer(C_INT), intent(in) :: overwrite
! *********************************************************************

 ierr = setenv(trim(name)//C_NULL_CHAR, trim(value)//C_NULL_CHAR, overwrite)

end function clib_setenv
!!***

!integer(c_long) function cache_size()
!  interface
!     function sysconf(name) bind(C, name="sysconf")
!       import :: c_int, c_long
!       integer(c_int), value :: name
!       integer(c_long) :: sysconf
!     end function sysconf
!  end interface
!
!  integer(c_int), parameter :: SC_LEVEL1_DCACHE_SIZE = 190  ! POSIX macro on Linux not portable. should use C directly
!
!  cache_size = sysconf(SC_LEVEL1_DCACHE_SIZE)
!  if (cache_size > 0) then
!     print *, "L1 Data Cache Size:", cache_size, "bytes"
!  else
!     print *, "Could not determine cache size."
!  end if
!end function cache_size

END MODULE m_clib
!!***

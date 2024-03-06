!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2024 ABINIT group (MG)
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
 public :: clib_getpid
 !public :: clib_usleep         ! Suspend calling thread for microseconds of clock time


!FIXME the interfaces below have been commented out since abilint
! JB : because interface must have a name in abilint

! ===================================================
! ==== Fortran-bindings declared in fsi_posix.c ====
! ===================================================
! interface
!   subroutine clib_mkdir(path, ierr)
!     import
!     character(len=*),intent(in) :: path
!     integer(c_int),intent(out) :: ierr
!   end subroutine clib_mkdir
! end interface
!

 interface
   integer(c_int) function c_rename(oldname, newname) bind(C, name='rename')
     import
     character(kind=c_char),intent(in) :: oldname(*)
     character(kind=c_char),intent(in) :: newname(*)
   end function c_rename
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
   ! pid_t getpid().
   ! The type of pid_t data is a signed integer type (signed int or we can say int).
   function clib_getpid() bind(C, name='getpid')
     import
     integer(c_int) :: clib_getpid
   end function clib_getpid
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

  !interface
  !  ! suspend calling thread for microseconds of clock time
  !  ! uses unistd.h for Fortran standard compliant sleep.
  !  ! sleep() is a GNU extension, not standard Fortran
  !  subroutine usleep(us) bind(C)
  !    import
  !    integer(c_int), value :: us
  !  end subroutine usleep
  !end interface

  !interface
  !  ! int usleep(useconds_t useconds)
  !  function clib_usleep(useconds) bind(c, name='usleep')
  !    import
  !    integer(kind=c_int32_t), value :: useconds
  !    integer(kind=c_int)            :: c_usleep
  !  end function clib_usleep
  !end interface

! ==========================================
! ==== Fortran-bindings for file_lock.c ====
! ==========================================

 !interface
 !  function lock_file(path) bind(C)
 !    import
 !    implicit none
 !    character(kind=c_char),intent(in) :: path(*)
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
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

integer function clib_rename(old_fname, new_fname) result(ierr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: old_fname, new_fname

! *********************************************************************

 ierr = c_rename(trim(old_fname)//c_null_char, trim(new_fname)//c_null_char)

end function clib_rename
!!***

END MODULE m_clib
!!***

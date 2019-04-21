!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG)
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

 use iso_c_binding

 implicit none

 private

 public :: clib_rename
 public :: clib_cclock
 public :: clib_etime
 public :: clib_mtrace
 public :: clib_print_mallinfo

!FIXME the interfaces below have been commented out since abilint
! JB : because interface must have a name in abilint

!! ===================================================
!! ==== Fortran-bindings declared in fsi_posix.c ====
!! ===================================================
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
   end function
 end interface

 interface
   subroutine clib_cclock(cpu) bind(C, name="cclock")
     import
     real(c_double),intent(out) :: cpu
   end subroutine
 end interface

 interface
   real(c_double) function clib_etime(tt) bind(C, name="etime") result(res)
     import
     real(c_float),intent(out) :: tt(2)
   end function
 end interface

!! =================================================
!! ==== Fortran-bindings declared in mallinfo.c ====
!! =================================================
 interface
   subroutine clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks) bind(C, name="clib_mallinfo")
     import
     integer(c_long),intent(out) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
   end subroutine clib_mallinfo
 end interface

!! ==================================================
!! ==== Fortran-bindings declared in gnu_tools.c ====
!! ==================================================

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

!!****f* m_clib/fmallinfo
!! NAME
!!   fmallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clib_print_mallinfo(unit)

!Arguments ------------------------------------
 integer,intent(in) :: unit

!Local variables-------------------------------
 integer(c_long) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
! *********************************************************************

  call clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks)

  write(unit, *)""
  write(unit,*)' Total space in arena            : ',arena
  write(unit,*)' Space in holding block headers  : ',hblkhd
  write(unit,*)' Space in small blocks in use    : ',usmblks
  write(unit,*)' Space in free small blocks      : ',fsmblks
  write(unit,*)' Space in ordinary blocks in use : ',uordblks
  write(unit,*)' Space in free ordinary blocks   : ',fordblks
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
!! PARENTS
!!      isfile
!!
!! CHILDREN
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

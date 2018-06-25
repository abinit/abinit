!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_clib
!! NAME
!! m_clib
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (MG)
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

#ifdef HAVE_FC_ISO_C_BINDING
#define USE_MODULE 
#else
#define USE_MODULE use m_iso_c_binding
#endif

 use iso_c_binding

 implicit none

 private 

 integer, parameter :: dp=kind(1.0d0)
 integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  ! Complex should not be used presently

 type,public :: Mallinfo_t
   integer(C_LONG) :: arena 
   integer(C_LONG) :: hblkhd 
   integer(C_LONG) :: usmblks 
   integer(C_LONG) :: fsmblks 
   integer(C_LONG) :: uordblks 
   integer(C_LONG) :: fordblks
 end type Mallinfo_t

 public :: clib_rename

!FIXME the interfaces below have been commented out since abilint
! crashes during the analysis of the file (maybe due to the macro USE_MODULE!)      
! JB : No, it is because interface must have a name in abilint

! ===================================================
! ==== Fortran-bindings declared in intrinsics.c ====
! ===================================================
! interface 
!   subroutine clib_fflush()
!     import 
!     implicit none
!   end subroutine clib_fflush
! end interface
!
! interface 
!   subroutine clib_getenv(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getenv
! end interface
!
!! ===================================================
!! ==== Fortran-bindings declared in fsi_posix.c ====
!! ===================================================
! interface 
!   subroutine clib_mkdir(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_mkdir
! end interface
!
! interface 
!   subroutine clib_chdir(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_chdir
! end interface
!
interface c_rename
    integer(c_int) function rename(oldname,newname) bind(C,name='rename')
      use iso_c_binding
      implicit none
      character(kind=c_char),intent(in) :: oldname(*)
      character(kind=c_char),intent(in) :: newname(*)
    end function rename
 end interface c_rename
!
! interface 
!   subroutine clib_remove(ierr,fname)
!     import 
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_remove
! end interface
!
! interface 
!   subroutine clib_getcwd(ierr,fname)
!     import 
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_getcwd
! end interface
!
! interface 
!   subroutine clib_gethname(ierr,fname)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!     character(len=*),intent(in) :: fname
!   end subroutine clib_gethname
! end interface
!
!! =====================================================
!! ==== Fortran-bindings declared in progress_bar.c ====
!! =====================================================
! interface
!   subroutine clib_progress_bar(actual, max)
!     import
!     implicit none
!     integer(C_INT),intent(in) :: actual
!     integer(C_INT),intent(in) :: max
!   end subroutine clib_progress_bar
! end interface
!
!! =================================================
!! ==== Fortran-bindings declared in mallinfo.c ====
!! =================================================
! interface
!   subroutine clib_mallinfo(arena, hblkhd, usmblks, fsmblks, uordblks, fordblks)
!     import
!     implicit none
!     integer(C_LONG),intent(out) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
!   end subroutine clib_mallinfo
! end interface
!
!
!! ==================================================
!! ==== Fortran-bindings declared in gnu_tools.c ====
!! ==================================================
!
! interface
!   subroutine clib_mtrace(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mtrace
! end interface
!
! interface
!   subroutine clib_muntrace(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_muntrace
! end interface
!
! interface
!   subroutine clib_mcheck(ierr)
!     import
!     implicit none
!     integer(C_INT),intent(out) :: ierr
!   end subroutine clib_mcheck
! end interface

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


CONTAINS  !===========================================================
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
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine fmallinfo(Minfo)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fmallinfo'
!End of the abilint section

 type(Mallinfo_t),intent(out) :: Minfo

!Local variables-------------------------------
 integer(C_LONG) :: arena,hblkhd,usmblks,fsmblks,uordblks,fordblks
! *********************************************************************

  call clib_mallinfo(arena,hblkhd,usmblks,fsmblks,uordblks,fordblks) 

  Minfo%arena    = arena
  Minfo%hblkhd   = hblkhd
  Minfo%usmblks  = usmblks
  Minfo%fsmblks  = fsmblks
  Minfo%uordblks = uordblks
  Minfo%fordblks = fordblks

end subroutine fmallinfo 
!!***

!!****f* m_clib/clib_print_mallinfo
!! NAME
!! clib_print_mallinfo
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clib_print_mallinfo(Minfo,unt)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clib_print_mallinfo'
!End of the abilint section

 integer,intent(in) :: unt
 type(Mallinfo_t),intent(in) :: Minfo
! *********************************************************************

 write(unt,*)' Total space in arena            : ',Minfo%arena
 write(unt,*)' Space in holding block headers  : ',Minfo%hblkhd
 write(unt,*)' Space in small blocks in use    : ',Minfo%usmblks
 write(unt,*)' Space in free small blocks      : ',Minfo%fsmblks
 write(unt,*)' Space in ordinary blocks in use : ',Minfo%uordblks
 write(unt,*)' Space in free ordinary blocks   : ',Minfo%fordblks
 write(unt,*)' End memory statistics '

end subroutine clib_print_mallinfo
!!***

!----------------------------------------------------------------------


!!****f* m_clib/clib_show_fc_alignment
!! NAME
!! clib_show_fc_alignment
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine clib_show_fc_alignment(unt)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clib_show_fc_alignment'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: unt

!Local variables-------------------------------
 integer,parameter :: sz=1024
 integer :: algn
!arrays
 real(dp),allocatable :: arr(:)

! *************************************************************************

 !allocate(arr(sz))
 call clib_alignment_of(arr, 16, algn)
 write(unt,*)" dp_arr: p % 16 = ",algn
 call clib_alignment_of(arr, 64, algn)
 write(unt,*)" dp_arr: p % 64 = ",algn
 !deallocate(arr)

end subroutine clib_show_fc_alignment
!!***

!!****f* m_clib/clib_rename
!! NAME
!!  clib_rename
!!
!! FUNCTION
!!  Rename a file with a new name
!!  It uses the C rename function from stdlib
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

subroutine clib_rename(old_fname, new_fname, ierr)

!Arguments ------------------------------------

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'clib_rename'
!End of the abilint section

 character(len=*),intent(in) :: old_fname, new_fname
 integer, optional, intent(out) :: ierr
 integer :: ier

!Local variables-------------------------------

! *********************************************************************

 ier = rename(trim(old_fname)//c_null_char, trim(new_fname)//c_null_char)
 if ( present(ierr) ) ierr = ier

end subroutine clib_rename
!!***
END MODULE m_clib
!!***


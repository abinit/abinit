!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_copy
!! NAME
!!  m_copy
!!
!! FUNCTION
!!  This module provides a generic interface used to copy pointers:
!!  deep_copy: used to return a deep copy of pointers. The procedure is useful if data types
!!             with several pointers have to be copied.
!!  addr_copy: used to copy the address contained in a pointer
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  * The intent for pointer arguments is not specified since
!!    we have to conform to the F90 specifications. However xval is IN while copy is OUT
!!
!!  * copy is a pointer and is supposed to be *not allocated*.
!!    If the value to be copied points to null(), also the copy will be nullified.
!!
!!  * On the alloc_copy routine:
!!   Since copy is INTENT(OUT), if the associated actual argument is
!!   currently allocated, the actual argument is deallocated on procedure invocation so that the dummy
!!   argument has an allocation status of not currently allocated.
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

MODULE m_copy

 use defs_basis,  only : dp, spc, dpc
 use m_abicore
#if defined HAVE_FC_ISO_C_BINDING
 use iso_c_binding, only : c_ptr,c_loc,c_f_pointer
#endif

 implicit none

 private

 public :: deep_copy     ! Performs deep copy of two pointers
 public :: alloc_copy    ! Allocate an allocable array and copy data. See notes in alloc_copy_int1d
 public :: addr_copy     ! Performs a bitwise copy of a pointer (copy address)

 interface deep_copy
  module procedure deep_copy_int0d
  module procedure deep_copy_int1d
  module procedure deep_copy_int2d
  module procedure deep_copy_int3d
  module procedure deep_copy_int4d
  module procedure deep_copy_rdp0d
  module procedure deep_copy_rdp1d
  module procedure deep_copy_rdp2d
  module procedure deep_copy_rdp3d
  module procedure deep_copy_rdp4d
  module procedure deep_copy_csp0d
  module procedure deep_copy_csp1d
  module procedure deep_copy_csp2d
  module procedure deep_copy_csp3d
  module procedure deep_copy_csp4d
  module procedure deep_copy_cdp0d
  module procedure deep_copy_cdp1d
  module procedure deep_copy_cdp2d
  module procedure deep_copy_cdp3d
  module procedure deep_copy_cdp4d
  module procedure deep_copy_log0d
  module procedure deep_copy_log1d
  module procedure deep_copy_log2d
  module procedure deep_copy_log3d
  module procedure deep_copy_log4d
  !module procedure deep_copy_ch1d  !Does not work on XLF, do not use it for the time being.
 end interface deep_copy

 interface alloc_copy
  module procedure alloc_copy_int1d
  module procedure alloc_copy_int2d
  module procedure alloc_copy_int3d
  module procedure alloc_copy_int4d
  module procedure alloc_copy_rdp1d
  module procedure alloc_copy_rdp2d
  module procedure alloc_copy_rdp3d
  module procedure alloc_copy_rdp4d
  module procedure alloc_copy_rdp5d
  module procedure alloc_copy_rdp6d
  module procedure alloc_copy_csp1d
  module procedure alloc_copy_csp2d
  module procedure alloc_copy_csp3d
  module procedure alloc_copy_csp4d
  module procedure alloc_copy_cdp1d
  module procedure alloc_copy_cdp2d
  module procedure alloc_copy_cdp3d
  module procedure alloc_copy_cdp4d
  module procedure alloc_copy_log1d
  module procedure alloc_copy_log2d
  module procedure alloc_copy_log3d
  module procedure alloc_copy_log4d
 end interface alloc_copy

 interface addr_copy
  module procedure addr_copy_int1d
  module procedure addr_copy_int2d
  module procedure addr_copy_int3d
  module procedure addr_copy_int4d
  module procedure addr_copy_dp1d
  module procedure addr_copy_dp2d
  module procedure addr_copy_dp3d
  module procedure addr_copy_dp4d
  module procedure addr_copy_dp5d
 end interface addr_copy

CONTAINS  !===========================================================
!!***

!!****f* m_copy/deep_copy_int0d
!! NAME
!! deep_copy_int0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_int0d(xval,copy)

!Arguments ------------------------------------
 integer,intent(in) :: xval
 integer,intent(out) :: copy
! *********************************************************************

  copy=xval

end subroutine deep_copy_int0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int1d
!! NAME
!! deep_copy_int1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_int1d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:)
 integer,pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int2d
!! NAME
!! deep_copy_int2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_int2d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:)
 integer,pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int3d
!! NAME
!! deep_copy_int3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_int3d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:,:)
 integer,pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_int4d
!! NAME
!! deep_copy_int4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_int4d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:,:,:)
 integer,pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_int4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp0d
!! NAME
!! deep_copy_rdp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_rdp0d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval
 real(dp),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_rdp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp1d
!! NAME
!! deep_copy_rdp1d

!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_rdp1d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:)
 real(dp),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_rdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp2d
!! NAME
!! deep_copy_rdp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_rdp2d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:)
 real(dp),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_rdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp3d
!! NAME
!! deep_copy_rdp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_rdp3d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:,:)
 real(dp),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_rdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_rdp4d
!! NAME
!! deep_copy_rdp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_rdp4d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:,:,:)
 real(dp),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_rdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp0d
!! NAME
!! deep_copy_csp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_csp0d(xval,copy)

!Arguments ------------------------------------
 complex(spc),intent(in) :: xval
 complex(spc),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_csp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp1d
!! NAME
!! deep_copy_csp1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_csp1d(xval,copy)

!Arguments ------------------------------------
 complex(spc),pointer :: xval(:)
 complex(spc),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_csp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp2d
!! NAME
!! deep_copy_csp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_csp2d(xval,copy)

!Arguments ------------------------------------
 complex(spc),pointer :: xval(:,:)
 complex(spc),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_csp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp3d
!! NAME
!! deep_copy_csp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_csp3d(xval,copy)

!Arguments ------------------------------------
 complex(spc),pointer :: xval(:,:,:)
 complex(spc),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_csp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_csp4d
!! NAME
!! deep_copy_csp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_csp4d(xval,copy)

!Arguments ------------------------------------
 complex(spc),pointer :: xval(:,:,:,:)
 complex(spc),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_csp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp0d
!! NAME
!! deep_copy_cdp0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_cdp0d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: xval
 complex(dpc),intent(out) :: copy
! *********************************************************************
  copy=xval

end subroutine deep_copy_cdp0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp1d
!! NAME
!! deep_copy_cdp1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_cdp1d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),pointer :: xval(:)
 complex(dpc),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_cdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp2d
!! NAME
!! deep_copy_cdp2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_cdp2d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),pointer :: xval(:,:)
 complex(dpc),pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_cdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp3d
!! NAME
!! deep_copy_cdp3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_cdp3d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),pointer :: xval(:,:,:)
 complex(dpc),pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_cdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_cdp4d
!! NAME
!! deep_copy_cdp4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_cdp4d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),pointer :: xval(:,:,:,:)
 complex(dpc),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_MALLOC(copy,(il1:iu1,il2:il2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_cdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_ch1d
!! NAME
!! deep_copy_ch1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! NOTES
!!  This routine segfaults on XLF, disabled for the time being
!!  Should test whether passing slen fixes the problem
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_ch1d(xval,copy,slen)

!Arguments ------------------------------------
 integer,intent(in) :: slen
 character(len=slen),pointer :: xval(:)
 character(len=slen),pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_ch1d
!!***

!!****f* m_copy/deep_copy_log0d
!! NAME
!! deep_copy_log0d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_log0d(xval,copy)

!Arguments ------------------------------------
 logical,intent(in) :: xval
 logical,intent(out) :: copy
! *********************************************************************

  copy=xval

end subroutine deep_copy_log0d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log1d
!! NAME
!! deep_copy_log1d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_log1d(xval,copy)

!Arguments ------------------------------------
 logical,pointer :: xval(:)
 logical,pointer :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 if (associated(xval)) then
  il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
  ABI_MALLOC(copy,(il:iu))
  copy(:)=xval(:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log2d
!! NAME
!! deep_copy_log2d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_log2d(xval,copy)

!Arguments ------------------------------------
 logical,pointer :: xval(:,:)
 logical,pointer :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2))
  copy(:,:)=xval(:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log3d
!! NAME
!! deep_copy_log3d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_log3d(xval,copy)

!Arguments ------------------------------------
 logical,pointer :: xval(:,:,:)
 logical,pointer :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
  copy(:,:,:)=xval(:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/deep_copy_log4d
!! NAME
!! deep_copy_log4d
!!
!! FUNCTION
!!  Performs a deep copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine deep_copy_log4d(xval,copy)

!Arguments ------------------------------------
 logical,pointer :: xval(:,:,:,:)
 logical,pointer :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 if (associated(xval)) then
  il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
  il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
  il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
  il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
  ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
  copy(:,:,:,:)=xval(:,:,:,:)
 else
  nullify(copy)
 end if

end subroutine deep_copy_log4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_int1d
!! NAME
!! alloc_copy_int1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_int1d(xval,copy)

!Arguments ------------------------------------
 integer,intent(in) :: xval(:)
 integer,allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_int1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_int2d
!! NAME
!! alloc_copy_int2d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_int2d(xval,copy)

!Arguments ------------------------------------
 integer,intent(in) :: xval(:,:)
 integer,allocatable,intent(out) :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2))
 copy(:,:)=xval(:,:)

end subroutine alloc_copy_int2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_int3d
!! NAME
!! alloc_copy_int3d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_int3d(xval,copy)

!Arguments ------------------------------------
 integer,intent(in) :: xval(:,:,:)
 integer,allocatable,intent(out) :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
 copy(:,:,:)=xval(:,:,:)

end subroutine alloc_copy_int3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_int4d
!! NAME
!! alloc_copy_int4d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_int4d(xval,copy)

!Arguments ------------------------------------
 integer,intent(in) :: xval(:,:,:,:)
 integer,allocatable,intent(out) :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
 copy(:,:,:,:)=xval(:,:,:,:)

end subroutine alloc_copy_int4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp1d
!! NAME
!! alloc_copy_rdp1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp1d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:)
 real(dp),allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_rdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp2d
!! NAME
!! alloc_copy_rdp2d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp2d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:,:)
 real(dp),allocatable,intent(out) :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2))
 copy(:,:)=xval(:,:)

end subroutine alloc_copy_rdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp3d
!! NAME
!! alloc_copy_rdp3d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp3d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:,:,:)
 real(dp),allocatable,intent(out) :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
 copy(:,:,:)=xval(:,:,:)

end subroutine alloc_copy_rdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp4d
!! NAME
!! alloc_copy_rdp4d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp4d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:,:,:,:)
 real(dp),allocatable,intent(out) :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
 copy(:,:,:,:)=xval(:,:,:,:)

end subroutine alloc_copy_rdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp5d
!! NAME
!! alloc_copy_rdp5d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp5d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:,:,:,:,:)
 real(dp),allocatable,intent(out) :: copy(:,:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4,il5,iu5
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 il5=lbound(xval,DIM=5); iu5=ubound(xval,DIM=5)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4,il5:iu5))
 copy(:,:,:,:,:)=xval(:,:,:,:,:)

end subroutine alloc_copy_rdp5d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_rdp6d
!! NAME
!! alloc_copy_rdp6d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_rdp6d(xval,copy)

!Arguments ------------------------------------
 real(dp),intent(in) :: xval(:,:,:,:,:,:)
 real(dp),allocatable,intent(out) :: copy(:,:,:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4,il5,iu5,il6,iu6
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 il5=lbound(xval,DIM=5); iu5=ubound(xval,DIM=5)
 il6=lbound(xval,DIM=6); iu6=ubound(xval,DIM=6)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4,il5:iu5,il6:iu6))
 copy(:,:,:,:,:,:)=xval(:,:,:,:,:,:)

end subroutine alloc_copy_rdp6d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_csp1d
!! NAME
!! alloc_copy_csp1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_csp1d(xval,copy)

!Arguments ------------------------------------
 complex(spc),intent(in) :: xval(:)
 complex(spc),allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_csp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_csp2d
!! NAME
!! alloc_copy_csp2d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_csp2d(xval,copy)

!Arguments ------------------------------------
 complex(spc),intent(in) :: xval(:,:)
 complex(spc),allocatable,intent(out) :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2))
 copy(:,:)=xval(:,:)

end subroutine alloc_copy_csp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_csp3d
!! NAME
!! alloc_copy_csp3d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_csp3d(xval,copy)

!Arguments ------------------------------------
 complex(spc),intent(in) :: xval(:,:,:)
 complex(spc),allocatable,intent(out) :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
 copy(:,:,:)=xval(:,:,:)

end subroutine alloc_copy_csp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_csp4d
!! NAME
!! alloc_copy_csp4d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_csp4d(xval,copy)

!Arguments ------------------------------------
 complex(spc),intent(in) :: xval(:,:,:,:)
 complex(spc),allocatable,intent(out) :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
 copy(:,:,:,:)=xval(:,:,:,:)

end subroutine alloc_copy_csp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_cdp1d
!! NAME
!! alloc_copy_cdp1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_cdp1d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: xval(:)
 complex(dpc),allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_cdp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_cdp2d
!! NAME
!! alloc_copy_cdp2d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_cdp2d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: xval(:,:)
 complex(dpc),allocatable,intent(out) :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2))
 copy(:,:)=xval(:,:)

end subroutine alloc_copy_cdp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_cdp3d
!! NAME
!! alloc_copy_cdp3d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_cdp3d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: xval(:,:,:)
 complex(dpc),allocatable,intent(out) :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
 copy(:,:,:)=xval(:,:,:)

end subroutine alloc_copy_cdp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_cdp4d
!! NAME
!! alloc_copy_cdp4d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_cdp4d(xval,copy)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: xval(:,:,:,:)
 complex(dpc),allocatable,intent(out) :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 ABI_MALLOC(copy,(il1:iu1,il2:il2,il3:iu3,il4:iu4))
 copy(:,:,:,:)=xval(:,:,:,:)

end subroutine alloc_copy_cdp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_ch1d
!! NAME
!! alloc_copy_ch1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! NOTES
!!  This routine segfaults on XLF, disabled for the time being
!!  Should test whether passing slen fixes the problem
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_ch1d(xval,copy,slen)

!Arguments ------------------------------------
 integer,intent(in) :: slen
 character(len=slen),intent(in) :: xval(:)
 character(len=slen),allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_ch1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_log1d
!! NAME
!! alloc_copy_log1d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_log1d(xval,copy)

!Arguments ------------------------------------
 logical,intent(in) :: xval(:)
 logical,allocatable,intent(out) :: copy(:)

!Local variables-------------------------------
 integer :: il,iu
! *********************************************************************

 il=lbound(xval,DIM=1); iu=ubound(xval,DIM=1)
 ABI_MALLOC(copy,(il:iu))
 copy(:)=xval(:)

end subroutine alloc_copy_log1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_log2d
!! NAME
!! alloc_copy_log2d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_log2d(xval,copy)

!Arguments ------------------------------------
 logical,intent(in) :: xval(:,:)
 logical,allocatable,intent(out) :: copy(:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2))
 copy(:,:)=xval(:,:)

end subroutine alloc_copy_log2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_log3d
!! NAME
!! alloc_copy_log3d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_log3d(xval,copy)

!Arguments ------------------------------------
 logical,intent(in) :: xval(:,:,:)
 logical,allocatable,intent(out) :: copy(:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3))
 copy(:,:,:)=xval(:,:,:)

end subroutine alloc_copy_log3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/alloc_copy_log4d
!! NAME
!! alloc_copy_log4d
!!
!! FUNCTION
!!  Performs a copy of an array.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine alloc_copy_log4d(xval,copy)

!Arguments ------------------------------------
 logical,intent(in) :: xval(:,:,:,:)
 logical,allocatable,intent(out) :: copy(:,:,:,:)

!Local variables-------------------------------
 integer :: il1,iu1,il2,iu2,il3,iu3,il4,iu4
! *********************************************************************

 il1=lbound(xval,DIM=1); iu1=ubound(xval,DIM=1)
 il2=lbound(xval,DIM=2); iu2=ubound(xval,DIM=2)
 il3=lbound(xval,DIM=3); iu3=ubound(xval,DIM=3)
 il4=lbound(xval,DIM=4); iu4=ubound(xval,DIM=4)
 ABI_MALLOC(copy,(il1:iu1,il2:iu2,il3:iu3,il4:iu4))
 copy(:,:,:,:)=xval(:,:,:,:)

end subroutine alloc_copy_log4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_int1d
!! NAME
!! addr_copy_int1d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_int1d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:)
 integer,pointer :: copy(:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(1)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_int1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_int2d
!! NAME
!! addr_copy_int2d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_int2d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:)
 integer,pointer :: copy(:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(2)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_int2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_int3d
!! NAME
!! addr_copy_int3d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_int3d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:,:)
 integer,pointer :: copy(:,:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(3)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_int3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_int4d
!! NAME
!! addr_copy_int4d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_int4d(xval,copy)

!Arguments ------------------------------------
 integer,pointer :: xval(:,:,:,:)
 integer,pointer :: copy(:,:,:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(4)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1,1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0,0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_int4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_dp1d
!! NAME
!! addr_copy_dp1d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_dp1d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:)
 real(dp),pointer :: copy(:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(1)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_dp1d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_dp2d
!! NAME
!! addr_copy_dp2d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_dp2d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:)
 real(dp),pointer :: copy(:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(2)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_dp2d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_dp3d
!! NAME
!! addr_copy_dp3d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_dp3d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:,:)
 real(dp),pointer :: copy(:,:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(3)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_dp3d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_dp4d
!! NAME
!! addr_copy_dp4d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_dp4d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:,:,:)
 real(dp),pointer :: copy(:,:,:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(4)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1,1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0,0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_dp4d
!!***

!----------------------------------------------------------------------

!!****f* m_copy/addr_copy_dp5d
!! NAME
!! addr_copy_dp5d
!!
!! FUNCTION
!!  Performs a bitwise copy of a pointer.
!!
!! PARENTS
!!
!! CHILDREN
!!      c_f_pointer
!!
!! SOURCE

subroutine addr_copy_dp5d(xval,copy)

!Arguments ------------------------------------
 real(dp),pointer :: xval(:,:,:,:,:)
 real(dp),pointer :: copy(:,:,:,:,:)

!Local variables-------------------------------
#if defined HAVE_FC_ISO_C_BINDING
 integer :: shp(5)
 type(C_PTR) :: ham_ptr
#endif

! *********************************************************************

 if (associated(xval)) then
#if defined HAVE_FC_ISO_C_BINDING
   shp=shape(xval)
   if (product(shp)>0) then
     ham_ptr=c_loc(xval(1,1,1,1,1))
     call c_f_pointer(ham_ptr,copy,shp)
   else
     ABI_ALLOCATE(copy,(0,0,0,0,0))
   end if
#else
   copy=transfer(xval,copy)
#endif
 else
  nullify(copy)
 end if

end subroutine addr_copy_dp5d
!!***

END MODULE m_copy
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkrprimd
!!
!! NAME
!! chkrprimd
!!
!! FUNCTION
!! Test if {rprim,acell,rprimd}
!! It means that rprimd can be reconstructed
!! from the rprim and acell
!! Output a message if is not the case
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!  (only writing)
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


subroutine chkrprimd(acell,rprim,rprimd,iout)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkrprimd'
!End of the abilint section

implicit none

!Arguments ------------------------------------
!scalars
integer,intent(in) :: iout
!arrays
real(dp),intent(in) :: rprim(3,3)
real(dp),intent(in) :: rprimd(3,3)
real(dp),intent(in) :: acell(3)

!Local variables-------------------------------
!scalars
integer :: ii,jj
!arrays
real(dp) :: rprimd_test(3,3)
logical :: equal

! ***********************************************************

!###########################################################
!### 1. Compute rprimd from rprim and acell
 do ii=1,3
   do jj=1,3
     rprimd_test(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do


!###########################################################
!### 2. Compare rprimd and rprimd_test

 equal=.TRUE.
 do ii=1,3
   do jj=1,3
     if (abs(rprimd_test(ii,jj)-rprimd(ii,jj))>1.E-12) then
       equal=.FALSE.
     end if
   end do
 end do

 if (equal)then
   write(iout,*) 'chkrprimd: rprimd is consistent'
 else
   write(iout,*) 'chkrprimd: rprimd is NOT consistent ERROR'
 end if

end subroutine chkrprimd
!!***

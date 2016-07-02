!{\src2tex{textfont=tt}}
!!****f* ABINIT/mati3det
!! NAME
!! mati3det
!!
!! FUNCTION
!! COmpute the determinant of a 3x3 matrix of INTEGER elements.
!!
!! COPYRIGHT
!! Copyright (C) 2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! mm = integer matrix 
!!
!! OUTPUT
!! det = determinant of the matrix
!!
!! NOTES
!!
!! PARENTS
!!      get_kpt_fullbz,getspinrot,m_ab7_symmetry,symdet
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine mati3det(mm,det)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mati3det'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 integer,intent(in) :: mm(3,3)
 integer,intent(out) :: det

!Local variables-------------------------------
!scalars

! *************************************************************************
 det=mm(1,1)*(mm(2,2) * mm(3,3) - mm(3,2) * mm(2,3)) &
& + mm(2,1)*(mm(3,2) * mm(1,3) - mm(1,2) * mm(3,3)) &
& + mm(3,1)*(mm(1,2) * mm(2,3) - mm(2,2) * mm(1,3))

end subroutine mati3det
!!***

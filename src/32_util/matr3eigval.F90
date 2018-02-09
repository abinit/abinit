!{\src2tex{textfont=tt}}
!!****f* ABINIT/matr3eigval
!! NAME
!! matr3eigval
!!
!! FUNCTION
!! Find the eigenvalues of a real symmetric 3x3 matrix,
!! entered in full storage mode.
!!
!! COPYRIGHT
!! Copyright (C) 2003-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  matr(3,3)=real symmetric 3x3 matrix
!!
!! OUTPUT
!!  eigval(3)=three eigenvalues
!!
!! PARENTS
!!      chkdilatmx
!!
!! CHILDREN
!!      zhpev
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine matr3eigval(eigval,matr)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matr3eigval'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: matr(3,3)
 real(dp),intent(out) :: eigval(3)

!Local variables-------------------------------
!scalars
 integer :: ier
!arrays
 real(dp) :: eigvec(2,3,3),matrx(2,6),zhpev1(2,2*3-1),zhpev2(3*3-2)

! *************************************************************************

 matrx(1,1)=matr(1,1)
 matrx(1,2)=matr(1,2)
 matrx(1,3)=matr(2,2)
 matrx(1,4)=matr(1,3)
 matrx(1,5)=matr(2,3)
 matrx(1,6)=matr(3,3)
 matrx(2,:)=zero

 call ZHPEV ('V','U',3,matrx,eigval,eigvec,3,zhpev1,zhpev2,ier)
!write(std_out,*)' eigval=',eigval

end subroutine matr3eigval
!!***

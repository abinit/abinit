!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotproduct
!! NAME
!! dotproduct
!!
!! FUNCTION
!! scalar product of two vectors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~ABINIT/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~ABINIT/Infos/contributors .
!!
!! INPUTS
!! v1 and v2: two real(dp) vectors
!!
!! OUTPUT
!! scalar product of the two vectors
!!
!! SIDE EFFECTS
!!
!! WARNINGS
!! vector size is not checked
!!
!! NOTES
!! I've benchmarked this to be speedier than the intrinsic dot_product even on
!! big vectors. The point is that less check is performed.
!!
!! PARENTS
!! cgpr,brent
!!
!! CHILDREN
!!
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


function dotproduct(nv1,nv2,v1,v2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotproduct'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nv1,nv2
 real(dp) :: dotproduct
!arrays
 real(dp),intent(in) :: v1(nv1,nv2),v2(nv1,nv2)

!Local variables-------------------------------
!scalars
 integer :: i,j

! *************************************************************************
 dotproduct=zero
 do j=1,nv2
  do i=1,nv1
   dotproduct=dotproduct+v1(i,j)*v2(i,j)
  end do
 end do
end function dotproduct

!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/dotproductqrc
!! NAME
!! dotproductqrc
!!
!! FUNCTION
!! compute the dot product of two vectors in reciprocal space,
!! one being real and one being complex.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (GMR, VO, LR, RWG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  b1(3),b2(3),b3(3)=the three primitive vectors in reciprocal space
!!  c(3)=the complex vector
!!  r(3)=the real vector
!!
!! OUTPUT
!!  dotproductqrc=The scalar product.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

 function dotproductqrc(r,c,b1,b2,b3)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dotproductqrc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 complex(gwpc) :: dotproductqrc
!arrays
 real(dp),intent(in) :: b1(3),b2(3),b3(3),r(3)
 complex(gwpc),intent(in) :: c(3)

!Local variables-------------------------------
!scalars
 integer :: i

! *************************************************************************
 dotproductqrc=(0.0_gwp,0.0_gwp)
 do i=1,3
   dotproductqrc=dotproductqrc+(r(1)*b1(i)+r(2)*b2(i)+r(3)*b3(i))*&
&                              (c(1)*b1(i)+c(2)*b2(i)+c(3)*b3(i))
 end do

end function dotproductqrc
!!***

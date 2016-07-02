!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont33cso
!! NAME
!! cont33cso
!!
!! FUNCTION
!! Contract symmetric rank 3 tensor gxa1 with symmetric rank 3 tensor
!! gxa2 using metric tensor gmet to produce rank 2 complex tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa1(2,10)=rank 3 complex symmetric tensor
!!  gxa2(2,10)=rank 3 complex symmetric tensor
!!  gmet(3,3)=usual metric tensor (symmetric, real)
!!
!! OUTPUT
!!  rank2c(2,6)=rank 2 complex tensor (pseudo-symmetric storage)
!!
!! NOTES
!! This contraction is used for spin-orbit correction in non-local
!! contribution to stresses.
!!
!! Symmetric gxa1, gxa2 are stored as
!!               111 221 331 321 311 211 222 332 322 333;
!! gmet(3,3) is symmetric but stored fully (9 elements);
!! Output rank2c is not symmetric but since
!!      $rank2c_{gxa1,gxa2}(a,b)=conjg(rank2c_{gxa2,gxa1}(b,a))$
!!       it is stored as 11 22 33 32 31 21.
!!
!! rank2c(1,1), rank2c(2,2), rank3c(3,3) are not needed;
!! They are not calculated.
!!
!!{{\ \begin{equation}
!! rank2c(a,b)=7.5 conjg(gxa1(a,i,j))*r_3(i,j,b) - 1.5 r_{11}(a)*r_{12}(b)
!!\end{equation} }}
!!   where:
!!{{\ \begin{eqnarray}
!! r_3(i,j,b) & = & gxa2(b,l,m) gmet(i,l) gmet(j,m) \nonumber
!! r_{11}(a)  & = & conjg(gxa1(a,l,m)) gmet(l,m)    \nonumber
!! r_{12}(b)  & = & gxa2(b,l,m) gmet(l,m)
!! \end{eqnarray} }}
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine cont33cso(gxa1,gxa2,gmet,rank2c)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont33cso'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa1(2,10),gxa2(2,10)
 real(dp),intent(out) :: rank2c(2,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
!arrays
 real(dp) :: r11(2,3),r12(2,3),r3(2,6,3),r3a(2,6,3)

! *************************************************************************

!Initialize output tensor
 rank2c(:,:)=0.d0

!Compute r3(i,j,b)=gxa2(b,l,m)*gmet(i,l)*gmet(j,m)
!stored as 11 22 33 32 31 21 for (i,j)
!First compute r3a(b,l,j)=gxa2(b,l,m)*gmet(j,m)
!stored as r3a(re/im,(b,l),j)
 r3a(:,1,1)=gxa2(:,1 )*gmet(1,1)+gxa2(:,6 )*gmet(1,2)+gxa2(:,5 )*gmet(1,3)
 r3a(:,2,1)=gxa2(:,2 )*gmet(1,1)+gxa2(:,7 )*gmet(1,2)+gxa2(:,9 )*gmet(1,3)
 r3a(:,3,1)=gxa2(:,3 )*gmet(1,1)+gxa2(:,8 )*gmet(1,2)+gxa2(:,10)*gmet(1,3)
 r3a(:,4,1)=gxa2(:,4 )*gmet(1,1)+gxa2(:,9 )*gmet(1,2)+gxa2(:,8 )*gmet(1,3)
 r3a(:,5,1)=gxa2(:,5 )*gmet(1,1)+gxa2(:,4 )*gmet(1,2)+gxa2(:,3 )*gmet(1,3)
 r3a(:,6,1)=gxa2(:,6 )*gmet(1,1)+gxa2(:,2 )*gmet(1,2)+gxa2(:,4 )*gmet(1,3)
 r3a(:,1,2)=gxa2(:,1 )*gmet(2,1)+gxa2(:,6 )*gmet(2,2)+gxa2(:,5 )*gmet(2,3)
 r3a(:,2,2)=gxa2(:,2 )*gmet(2,1)+gxa2(:,7 )*gmet(2,2)+gxa2(:,9 )*gmet(2,3)
 r3a(:,3,2)=gxa2(:,3 )*gmet(2,1)+gxa2(:,8 )*gmet(2,2)+gxa2(:,10)*gmet(2,3)
 r3a(:,4,2)=gxa2(:,4 )*gmet(2,1)+gxa2(:,9 )*gmet(2,2)+gxa2(:,8 )*gmet(2,3)
 r3a(:,5,2)=gxa2(:,5 )*gmet(2,1)+gxa2(:,4 )*gmet(2,2)+gxa2(:,3 )*gmet(2,3)
 r3a(:,6,2)=gxa2(:,6 )*gmet(2,1)+gxa2(:,2 )*gmet(2,2)+gxa2(:,4 )*gmet(2,3)
 r3a(:,1,3)=gxa2(:,1 )*gmet(3,1)+gxa2(:,6 )*gmet(3,2)+gxa2(:,5 )*gmet(3,3)
 r3a(:,2,3)=gxa2(:,2 )*gmet(3,1)+gxa2(:,7 )*gmet(3,2)+gxa2(:,9 )*gmet(3,3)
 r3a(:,3,3)=gxa2(:,3 )*gmet(3,1)+gxa2(:,8 )*gmet(3,2)+gxa2(:,10)*gmet(3,3)
 r3a(:,4,3)=gxa2(:,4 )*gmet(3,1)+gxa2(:,9 )*gmet(3,2)+gxa2(:,8 )*gmet(3,3)
 r3a(:,5,3)=gxa2(:,5 )*gmet(3,1)+gxa2(:,4 )*gmet(3,2)+gxa2(:,3 )*gmet(3,3)
 r3a(:,6,3)=gxa2(:,6 )*gmet(3,1)+gxa2(:,2 )*gmet(3,2)+gxa2(:,4 )*gmet(3,3)

!Then  compute r3(i,j,b)=r3a(b,l,j)*gmet(i,l)
!stored as r3(re/im,(i,j),b)
 r3(:,1,1)=r3a(:,1,1)*gmet(1,1)+r3a(:,6,1)*gmet(1,2)+r3a(:,5,1)*gmet(1,3)
 r3(:,1,2)=r3a(:,6,1)*gmet(1,1)+r3a(:,2,1)*gmet(1,2)+r3a(:,4,1)*gmet(1,3)
 r3(:,1,3)=r3a(:,5,1)*gmet(1,1)+r3a(:,4,1)*gmet(1,2)+r3a(:,3,1)*gmet(1,3)
 r3(:,2,1)=r3a(:,1,2)*gmet(2,1)+r3a(:,6,2)*gmet(2,2)+r3a(:,5,2)*gmet(2,3)
 r3(:,2,2)=r3a(:,6,2)*gmet(2,1)+r3a(:,2,2)*gmet(2,2)+r3a(:,4,2)*gmet(2,3)
 r3(:,2,3)=r3a(:,5,2)*gmet(2,1)+r3a(:,4,2)*gmet(2,2)+r3a(:,3,2)*gmet(2,3)
 r3(:,3,1)=r3a(:,1,3)*gmet(3,1)+r3a(:,6,3)*gmet(3,2)+r3a(:,5,3)*gmet(3,3)
 r3(:,3,2)=r3a(:,6,3)*gmet(3,1)+r3a(:,2,3)*gmet(3,2)+r3a(:,4,3)*gmet(3,3)
 r3(:,3,3)=r3a(:,5,3)*gmet(3,1)+r3a(:,4,3)*gmet(3,2)+r3a(:,3,3)*gmet(3,3)
 r3(:,4,1)=r3a(:,1,2)*gmet(3,1)+r3a(:,6,2)*gmet(3,2)+r3a(:,5,2)*gmet(3,3)
 r3(:,4,2)=r3a(:,6,2)*gmet(3,1)+r3a(:,2,2)*gmet(3,2)+r3a(:,4,2)*gmet(3,3)
 r3(:,4,3)=r3a(:,5,2)*gmet(3,1)+r3a(:,4,2)*gmet(3,2)+r3a(:,3,2)*gmet(3,3)
 r3(:,5,1)=r3a(:,1,1)*gmet(3,1)+r3a(:,6,1)*gmet(3,2)+r3a(:,5,1)*gmet(3,3)
 r3(:,5,2)=r3a(:,6,1)*gmet(3,1)+r3a(:,2,1)*gmet(3,2)+r3a(:,4,1)*gmet(3,3)
 r3(:,5,3)=r3a(:,5,1)*gmet(3,1)+r3a(:,4,1)*gmet(3,2)+r3a(:,3,1)*gmet(3,3)
 r3(:,6,1)=r3a(:,1,1)*gmet(2,1)+r3a(:,6,1)*gmet(2,2)+r3a(:,5,1)*gmet(2,3)
 r3(:,6,2)=r3a(:,6,1)*gmet(2,1)+r3a(:,2,1)*gmet(2,2)+r3a(:,4,1)*gmet(2,3)
 r3(:,6,3)=r3a(:,5,1)*gmet(2,1)+r3a(:,4,1)*gmet(2,2)+r3a(:,3,1)*gmet(2,3)

!Compute r11(a)=conjg(gxa1(a,l,m))*gmet(l,m)
 r11(:,1)=gxa1(:,1)*gmet(1,1)+gxa1(:,2)*gmet(2,2)+gxa1(:,3 )*gmet(3,3)&
& +2.d0*(gxa1(:,4)*gmet(3,2)+gxa1(:,5)*gmet(3,1)+gxa1(:,6 )*gmet(2,1))
 r11(:,2)=gxa1(:,6)*gmet(1,1)+gxa1(:,7)*gmet(2,2)+gxa1(:,8 )*gmet(3,3)&
& +2.d0*(gxa1(:,9)*gmet(3,2)+gxa1(:,4)*gmet(3,1)+gxa1(:,2 )*gmet(2,1))
 r11(:,3)=gxa1(:,5)*gmet(1,1)+gxa1(:,9)*gmet(2,2)+gxa1(:,10)*gmet(3,3)&
& +2.d0*(gxa1(:,8)*gmet(3,2)+gxa1(:,3)*gmet(3,1)+gxa1(:,4 )*gmet(2,1))
 r11(im,1)=-r11(im,1);r11(im,2)=-r11(im,2);r11(im,3)=-r11(im,3)

!Compute r12(b)=gxa2(b,l,m)*gmet(l,m)
 r12(:,1)=gxa2(:,1)*gmet(1,1)+gxa2(:,2)*gmet(2,2)+gxa2(:,3 )*gmet(3,3)&
& +2.d0*(gxa2(:,4)*gmet(3,2)+gxa2(:,5)*gmet(3,1)+gxa2(:,6 )*gmet(2,1))
 r12(:,2)=gxa2(:,6)*gmet(1,1)+gxa2(:,7)*gmet(2,2)+gxa2(:,8 )*gmet(3,3)&
& +2.d0*(gxa2(:,9)*gmet(3,2)+gxa2(:,4)*gmet(3,1)+gxa2(:,2 )*gmet(2,1))
 r12(:,3)=gxa2(:,5)*gmet(1,1)+gxa2(:,9)*gmet(2,2)+gxa2(:,10)*gmet(3,3)&
& +2.d0*(gxa2(:,8)*gmet(3,2)+gxa2(:,3)*gmet(3,1)+gxa2(:,4 )*gmet(2,1))

!Finally compute rank2c(a,b)=7.5*conjg(gxa1(a,i,j))*r3(i,j,b) - 1.5*r11(a)*r12(b)
!rank2c(re,1)=7.5d0*(gxa1(re,1 )*r3(re,1,1)+gxa1(im,1 )*r3(im,1,1)&
!&                    +gxa1(re,2 )*r3(re,2,1)+gxa1(im,2 )*r3(im,2,1)&
!&                    +gxa1(re,3 )*r3(re,3,1)+gxa1(im,3 )*r3(im,3,1))&
!&             +15.d0*(gxa1(re,4 )*r3(re,4,1)+gxa1(im,4 )*r3(im,4,1)&
!&                    +gxa1(re,5 )*r3(re,5,1)+gxa1(im,5 )*r3(im,5,1)&
!&                    +gxa1(re,6 )*r3(re,6,1)+gxa1(im,6 )*r3(im,6,1))&
!&             -1.5d0*(r11(re,1)*r12(re,1)-r11(im,1)*r12(im,1))
!rank2c(re,2)=7.5d0*(gxa1(re,6 )*r3(re,1,2)+gxa1(im,6 )*r3(im,1,2)&
!&                    +gxa1(re,7 )*r3(re,2,2)+gxa1(im,7 )*r3(im,2,2)&
!&                    +gxa1(re,8 )*r3(re,3,2)+gxa1(im,8 )*r3(im,3,2))&
!&             +15.d0*(gxa1(re,9 )*r3(re,4,2)+gxa1(im,9 )*r3(im,4,2)&
!&                    +gxa1(re,4 )*r3(re,5,2)+gxa1(im,4 )*r3(im,5,2)&
!&                    +gxa1(re,2 )*r3(re,6,2)+gxa1(im,2 )*r3(im,6,2))&
!&             -1.5d0*(r11(re,2)*r12(re,2)-r11(im,2)*r12(im,2))
!rank2c(re,3)=7.5d0*(gxa1(re,5 )*r3(re,1,3)+gxa1(im,5 )*r3(im,1,3)&
!&                    +gxa1(re,9 )*r3(re,2,3)+gxa1(im,9 )*r3(im,2,3)&
!&                    +gxa1(re,10)*r3(re,3,3)+gxa1(im,10)*r3(im,3,3))&
!&             +15.d0*(gxa1(re,8 )*r3(re,4,3)+gxa1(im,8 )*r3(im,4,3)&
!&                    +gxa1(re,3 )*r3(re,5,3)+gxa1(im,3 )*r3(im,5,3)&
!&                    +gxa1(re,4 )*r3(re,6,3)+gxa1(im,4 )*r3(im,6,3))&
!&             -1.5d0*(r11(re,3)*r12(re,3)-r11(im,3)*r12(im,3))
 rank2c(re,4)=7.5d0*(gxa1(re,5 )*r3(re,1,2)+gxa1(im,5 )*r3(im,1,2)&
& +gxa1(re,9 )*r3(re,2,2)+gxa1(im,9 )*r3(im,2,2)&
& +gxa1(re,10)*r3(re,3,2)+gxa1(im,10)*r3(im,3,2))&
& +15.d0*(gxa1(re,8 )*r3(re,4,2)+gxa1(im,8 )*r3(im,4,2)&
& +gxa1(re,3 )*r3(re,5,2)+gxa1(im,3 )*r3(im,5,2)&
& +gxa1(re,4 )*r3(re,6,2)+gxa1(im,4 )*r3(im,6,2))&
& -1.5d0*(r11(re,3)*r12(re,2)-r11(im,3)*r12(im,2))
 rank2c(re,5)=7.5d0*(gxa1(re,5 )*r3(re,1,1)+gxa1(im,5 )*r3(im,1,1)&
& +gxa1(re,9 )*r3(re,2,1)+gxa1(im,9 )*r3(im,2,1)&
& +gxa1(re,10)*r3(re,3,1)+gxa1(im,10)*r3(im,3,1))&
& +15.d0*(gxa1(re,8 )*r3(re,4,1)+gxa1(im,8 )*r3(im,4,1)&
& +gxa1(re,3 )*r3(re,5,1)+gxa1(im,3 )*r3(im,5,1)&
& +gxa1(re,4 )*r3(re,6,1)+gxa1(im,4 )*r3(im,6,1))&
& -1.5d0*(r11(re,3)*r12(re,1)-r11(im,3)*r12(im,1))
 rank2c(re,6)=7.5d0*(gxa1(re,6 )*r3(re,1,1)+gxa1(im,6 )*r3(im,1,1)&
& +gxa1(re,7 )*r3(re,2,1)+gxa1(im,7 )*r3(im,2,1)&
& +gxa1(re,8 )*r3(re,3,1)+gxa1(im,8 )*r3(im,3,1))&
& +15.d0*(gxa1(re,9 )*r3(re,4,1)+gxa1(im,9 )*r3(im,4,1)&
& +gxa1(re,4 )*r3(re,5,1)+gxa1(im,4 )*r3(im,5,1)&
& +gxa1(re,2 )*r3(re,6,1)+gxa1(im,2 )*r3(im,6,1))&
& -1.5d0*(r11(re,2)*r12(re,1)-r11(im,2)*r12(im,1))
!rank2c(im,1)=7.5d0*(gxa1(re,1 )*r3(im,1,1)-gxa1(im,1 )*r3(re,1,1)&
!&                    +gxa1(re,2 )*r3(im,2,1)-gxa1(im,2 )*r3(re,2,1)&
!&                    +gxa1(re,3 )*r3(im,3,1)-gxa1(im,3 )*r3(re,3,1))&
!&             +15.d0*(gxa1(re,4 )*r3(im,4,1)-gxa1(im,4 )*r3(re,4,1)&
!&                    +gxa1(re,5 )*r3(im,5,1)-gxa1(im,5 )*r3(re,5,1)&
!&                    +gxa1(re,6 )*r3(im,6,1)-gxa1(im,6 )*r3(re,6,1))&
!&             -1.5d0*(r11(re,1)*r12(im,1)+r11(im,1)*r12(re,1))
!rank2c(im,2)=7.5d0*(gxa1(re,6 )*r3(im,1,2)-gxa1(im,6 )*r3(re,1,2)&
!&                    +gxa1(re,7 )*r3(im,2,2)-gxa1(im,7 )*r3(re,2,2)&
!&                    +gxa1(re,8 )*r3(im,3,2)-gxa1(im,8 )*r3(re,3,2))&
!&             +15.d0*(gxa1(re,9 )*r3(im,4,2)-gxa1(im,9 )*r3(re,4,2)&
!&                    +gxa1(re,4 )*r3(im,5,2)-gxa1(im,4 )*r3(re,5,2)&
!&                    +gxa1(re,2 )*r3(im,6,2)-gxa1(im,2 )*r3(re,6,2))&
!&             -1.5d0*(r11(re,2)*r12(im,2)+r11(im,2)*r12(re,2))
!rank2c(im,3)=7.5d0*(gxa1(re,5 )*r3(im,1,3)-gxa1(im,5 )*r3(re,1,3)&
!&                    +gxa1(re,9 )*r3(im,2,3)-gxa1(im,9 )*r3(re,2,3)&
!&                    +gxa1(re,10)*r3(im,3,3)-gxa1(im,10)*r3(re,3,3))&
!&             +15.d0*(gxa1(re,8 )*r3(im,4,3)-gxa1(im,8 )*r3(re,4,3)&
!&                    +gxa1(re,3 )*r3(im,5,3)-gxa1(im,3 )*r3(re,5,3)&
!&                    +gxa1(re,4 )*r3(im,6,3)-gxa1(im,4 )*r3(re,6,3))&
!&             -1.5d0*(r11(re,3)*r12(im,3)+r11(im,3)*r12(re,3))
 rank2c(im,4)=7.5d0*(gxa1(re,5 )*r3(im,1,2)-gxa1(im,5 )*r3(re,1,2)&
& +gxa1(re,9 )*r3(im,2,2)-gxa1(im,9 )*r3(re,2,2)&
& +gxa1(re,10)*r3(im,3,2)-gxa1(im,10)*r3(re,3,2))&
& +15.d0*(gxa1(re,8 )*r3(im,4,2)-gxa1(im,8 )*r3(re,4,2)&
& +gxa1(re,3 )*r3(im,5,2)-gxa1(im,3 )*r3(re,5,2)&
& +gxa1(re,4 )*r3(im,6,2)-gxa1(im,4 )*r3(re,6,2))&
& -1.5d0*(r11(re,3)*r12(im,2)+r11(im,3)*r12(re,2))
 rank2c(im,5)=7.5d0*(gxa1(re,5 )*r3(im,1,1)-gxa1(im,5 )*r3(re,1,1)&
& +gxa1(re,9 )*r3(im,2,1)-gxa1(im,9 )*r3(re,2,1)&
& +gxa1(re,10)*r3(im,3,1)-gxa1(im,10)*r3(re,3,1))&
& +15.d0*(gxa1(re,8 )*r3(im,4,1)-gxa1(im,8 )*r3(re,4,1)&
& +gxa1(re,3 )*r3(im,5,1)-gxa1(im,3 )*r3(re,5,1)&
& +gxa1(re,4 )*r3(im,6,1)-gxa1(im,4 )*r3(re,6,1))&
& -1.5d0*(r11(re,3)*r12(im,1)+r11(im,3)*r12(re,1))
 rank2c(im,6)=7.5d0*(gxa1(re,6 )*r3(im,1,1)-gxa1(im,6 )*r3(re,1,1)&
& +gxa1(re,7 )*r3(im,2,1)-gxa1(im,7 )*r3(re,2,1)&
& +gxa1(re,8 )*r3(im,3,1)-gxa1(im,8 )*r3(re,3,1))&
& +15.d0*(gxa1(re,9 )*r3(im,4,1)-gxa1(im,9 )*r3(re,4,1)&
& +gxa1(re,4 )*r3(im,5,1)-gxa1(im,4 )*r3(re,5,1)&
& +gxa1(re,2 )*r3(im,6,1)-gxa1(im,2 )*r3(re,6,1))&
& -1.5d0*(r11(re,2)*r12(im,1)+r11(im,2)*r12(re,1))

end subroutine cont33cso
!!***

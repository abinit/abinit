!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont22cso
!! NAME
!! cont22cso
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa1 with symmetric rank 2 tensor
!! gxa2 using metric tensor gmet to produce rank 2 complex tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa1(2,10)=rank 2 complex symmetric tensor
!!  gxa2(2,10)=rank 2 complex symmetric tensor
!!  gmet(3,3)=usual metric tensor (symmetric, real)
!!
!! OUTPUT
!!  rank2c(2,6)=rank 2 complex tensor (pseudo-symmetric storage)
!!
!! NOTES
!! This contraction is used for spin-orbit correction in non-local
!! contribution to stresses.
!!
!! Symmetric gxa1, gxa2 are stored as 11 22 33 32 31 21;
!! gmet(3,3) is symmetric but stored fully (9 elements);
!! Output rank2c is not symmetric but since
!!      $rank2c_{gxa1,gxa2}(a,b)=conjg(rank2c_{gxa2,gxa1}(b,a))$
!!       it is stored as 11 22 33 32 31 21.
!!
!! rank2c(1,1), rank2c(2,2), rank3c(3,3) are not needed;
!! They are not calculated.
!!
!!{{\ \begin{equation}
!! rank2c(a,b)=3 conjg(gxa1(i,a)) gmet(i,j) gxa2(j,b)
!!\end{equation} }}
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


subroutine cont22cso(gxa1,gxa2,gmet,rank2c)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont22cso'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa1(2,6),gxa2(2,6)
 real(dp),intent(out) :: rank2c(2,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
!arrays
 real(dp) :: r2(2,3,3)

! *************************************************************************

!Initialize output tensor
 rank2c(:,:)=0.d0

!First compute r2(i,j) = gmet(i,k) gxa2(j,k)
 r2(:,1,1)=gmet(1,1)*gxa2(:,1)+gmet(1,2)*gxa2(:,6)+gmet(1,3)*gxa2(:,5)
 r2(:,1,2)=gmet(1,1)*gxa2(:,6)+gmet(1,2)*gxa2(:,2)+gmet(1,3)*gxa2(:,4)
 r2(:,1,3)=gmet(1,1)*gxa2(:,5)+gmet(1,2)*gxa2(:,4)+gmet(1,3)*gxa2(:,3)
 r2(:,2,1)=gmet(2,1)*gxa2(:,1)+gmet(2,2)*gxa2(:,6)+gmet(2,3)*gxa2(:,5)
 r2(:,2,2)=gmet(2,1)*gxa2(:,6)+gmet(2,2)*gxa2(:,2)+gmet(2,3)*gxa2(:,4)
 r2(:,2,3)=gmet(2,1)*gxa2(:,5)+gmet(2,2)*gxa2(:,4)+gmet(2,3)*gxa2(:,3)
 r2(:,3,1)=gmet(3,1)*gxa2(:,1)+gmet(3,2)*gxa2(:,6)+gmet(3,3)*gxa2(:,5)
 r2(:,3,2)=gmet(3,1)*gxa2(:,6)+gmet(3,2)*gxa2(:,2)+gmet(3,3)*gxa2(:,4)
 r2(:,3,3)=gmet(3,1)*gxa2(:,5)+gmet(3,2)*gxa2(:,4)+gmet(3,3)*gxa2(:,3)

!Then compute rank2c(a,b) = 3 conjg(gxa1(a,i)) r2(i,b)
!stored as 11 22 33 32 31 21
!rank2c(re,1) = 3.d0*(gxa1(re,1)*r2(re,1,1)+gxa1(im,1)*r2(im,1,1)&
!&                     +gxa1(re,6)*r2(re,2,1)+gxa1(im,6)*r2(im,2,1)&
!&                     +gxa1(re,5)*r2(re,3,1)+gxa1(im,5)*r2(im,3,1))
!rank2c(re,2) = 3.d0*(gxa1(re,6)*r2(re,1,2)+gxa1(im,6)*r2(im,1,2)&
!&                     +gxa1(re,2)*r2(re,2,2)+gxa1(im,2)*r2(im,2,2)&
!&                     +gxa1(re,4)*r2(re,3,2)+gxa1(im,4)*r2(im,3,2))
!rank2c(re,3) = 3.d0*(gxa1(re,5)*r2(re,1,3)+gxa1(im,5)*r2(im,1,3)&
!&                     +gxa1(re,4)*r2(re,2,3)+gxa1(im,4)*r2(im,2,3)&
!&                     +gxa1(re,3)*r2(re,3,3)+gxa1(im,3)*r2(im,3,3))
 rank2c(re,4) = 3.d0*(gxa1(re,5)*r2(re,1,2)+gxa1(im,5)*r2(im,1,2)&
& +gxa1(re,4)*r2(re,2,2)+gxa1(im,4)*r2(im,2,2)&
& +gxa1(re,3)*r2(re,3,2)+gxa1(im,3)*r2(im,3,2))
 rank2c(re,5) = 3.d0*(gxa1(re,5)*r2(re,1,1)+gxa1(im,5)*r2(im,1,1)&
& +gxa1(re,4)*r2(re,2,1)+gxa1(im,4)*r2(im,2,1)&
& +gxa1(re,3)*r2(re,3,1)+gxa1(im,3)*r2(im,3,1))
 rank2c(re,6) = 3.d0*(gxa1(re,6)*r2(re,1,1)+gxa1(im,6)*r2(im,1,1)&
& +gxa1(re,2)*r2(re,2,1)+gxa1(im,2)*r2(im,2,1)&
& +gxa1(re,4)*r2(re,3,1)+gxa1(im,4)*r2(im,3,1))
!rank2c(im,1) = 3.d0*(gxa1(re,1)*r2(im,1,1)-gxa1(im,1)*r2(re,1,1)&
!&                     +gxa1(re,6)*r2(im,2,1)-gxa1(im,6)*r2(re,2,1)&
!&                     +gxa1(re,5)*r2(im,3,1)-gxa1(im,5)*r2(re,3,1))
!rank2c(im,2) = 3.d0*(gxa1(re,6)*r2(im,1,2)-gxa1(im,6)*r2(re,1,2)&
!&                     +gxa1(re,2)*r2(im,2,2)-gxa1(im,2)*r2(re,2,2)&
!&                     +gxa1(re,4)*r2(im,3,2)-gxa1(im,4)*r2(re,3,2))
!rank2c(im,3) = 3.d0*(gxa1(re,5)*r2(im,1,3)-gxa1(im,5)*r2(re,1,3)&
!&                     +gxa1(re,4)*r2(im,2,3)-gxa1(im,4)*r2(re,2,3)&
!&                     +gxa1(re,3)*r2(im,3,3)-gxa1(im,3)*r2(re,3,3))
 rank2c(im,4) = 3.d0*(gxa1(re,5)*r2(im,1,2)-gxa1(im,5)*r2(re,1,2)&
& +gxa1(re,4)*r2(im,2,2)-gxa1(im,4)*r2(re,2,2)&
& +gxa1(re,3)*r2(im,3,2)-gxa1(im,3)*r2(re,3,2))
 rank2c(im,5) = 3.d0*(gxa1(re,5)*r2(im,1,1)-gxa1(im,5)*r2(re,1,1)&
& +gxa1(re,4)*r2(im,2,1)-gxa1(im,4)*r2(re,2,1)&
& +gxa1(re,3)*r2(im,3,1)-gxa1(im,3)*r2(re,3,1))
 rank2c(im,6) = 3.d0*(gxa1(re,6)*r2(im,1,1)-gxa1(im,6)*r2(re,1,1)&
& +gxa1(re,2)*r2(im,2,1)-gxa1(im,2)*r2(re,2,1)&
& +gxa1(re,4)*r2(im,3,1)-gxa1(im,4)*r2(re,3,1))

end subroutine cont22cso
!!***

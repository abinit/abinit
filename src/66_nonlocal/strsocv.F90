!{\src2tex{textfont=tt}}
!!****f* ABINIT/strsocv
!! NAME
!! strsocv
!!
!! FUNCTION
!! Convert from antisymmetric storage mode 3x3x3 rank3 tensor in reduced
!! coordinates "red" to symmetric storage mode 3x3 rank2 tensor in
!! cartesian coordinates "cart", using metric tensor "gprimd".
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  red(6,3)=3x3x3 tensor in antisymmetric storage mode,
!!           reduced coordinates
!!  gprimd(3,3)=reciprocal space dimensional primitive translations
!!
!! OUTPUT
!!  cart(6)=3x3 tensor in symmetric storage mode,
!!          cartesian coordinates
!!
!! NOTES
!! This routine is used to compute spin-orbit stress tensor.
!!
!! red is antisymmetric in first two indices and stored
!!     as 11 22 33 32 31 21.
!! cart is symmetric and stored as 11 22 33 32 31 21.
!!
!!{{\ \begin{eqnarray}
!! cart(1,1) & = &        & red(i,j,2) G(3,i) G(1,j) + red(i,j,3) G(1,i) G(2,j) \nonumber
!! cart(2,2) & = &        & red(i,j,1) G(2,i) G(3,j) + red(i,j,3) G(1,i) G(2,j) \nonumber
!! cart(3,3) & = &        & red(i,j,1) G(2,i) G(3,j) + red(i,j,2) G(3,i) G(1,j) \nonumber
!! cart(3,2) & = &  0.5 ( & red(i,j,3) G(1,i) G(3,j) + red(i,j,2) G(2,i) G(1,j)) \nonumber
!! cart(3,1) & = &  0.5 ( & red(i,j,3) G(3,i) G(2,j) + red(i,j,1) G(2,i) G(1,j)) \nonumber
!! cart(2,1) & = &  0.5 ( & red(i,j,2) G(3,i) G(2,j) + red(i,j,1) G(1,i) G(3,j))
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


subroutine strsocv(red,gprimd,cart)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strsocv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gprimd(3,3),red(6,3)
 real(dp),intent(out) :: cart(6)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work(3,3,3)

! *************************************************************************

 do ii=1,3
   work(1,1,ii)=0.d0
   work(2,2,ii)=0.d0
   work(3,3,ii)=0.d0
   work(3,2,ii)=red(4,ii) ; work(2,3,ii)=-red(4,ii)
   work(3,1,ii)=red(5,ii) ; work(1,3,ii)=-red(5,ii)
   work(2,1,ii)=red(6,ii) ; work(1,2,ii)=-red(6,ii)
 end do

 cart(:)=0.d0
 do jj=1,3
   do ii=1,3
     cart(1)=cart(1)+work(ii,jj,2)*gprimd(3,ii)*gprimd(1,jj)&
&     +work(ii,jj,3)*gprimd(1,ii)*gprimd(2,jj)
     cart(2)=cart(2)+work(ii,jj,1)*gprimd(2,ii)*gprimd(3,jj)&
&     +work(ii,jj,3)*gprimd(1,ii)*gprimd(2,jj)
     cart(3)=cart(3)+work(ii,jj,1)*gprimd(2,ii)*gprimd(3,jj)&
&     +work(ii,jj,2)*gprimd(3,ii)*gprimd(1,jj)
     cart(4)=cart(4)+work(ii,jj,3)*gprimd(1,ii)*gprimd(3,jj)&
&     +work(ii,jj,2)*gprimd(2,ii)*gprimd(1,jj)
     cart(5)=cart(5)+work(ii,jj,3)*gprimd(3,ii)*gprimd(2,jj)&
&     +work(ii,jj,1)*gprimd(2,ii)*gprimd(1,jj)
     cart(6)=cart(6)+work(ii,jj,2)*gprimd(3,ii)*gprimd(2,jj)&
&     +work(ii,jj,1)*gprimd(1,ii)*gprimd(3,jj)
   end do
 end do
 cart(4)=0.5d0*cart(4)
 cart(5)=0.5d0*cart(5)
 cart(6)=0.5d0*cart(6)

end subroutine strsocv
!!***

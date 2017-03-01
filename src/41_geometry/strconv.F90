!{\src2tex{textfont=tt}}
!!****f* ABINIT/strconv
!! NAME
!! strconv
!!
!! FUNCTION
!! If original gprimd is input, convert from symmetric storage mode
!! 3x3 tensor in reduced coordinates "frac" to symmetric storage mode
!! symmetric tensor in cartesian coordinates "cart".
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  frac(6)=3x3 tensor in symmetric storage mode, reduced coordinates
!!  gprimd(3,3)=reciprocal space dimensional primitive translations (bohr^-1)
!!
!! OUTPUT
!!  cart(6)=symmetric storage mode for symmetric 3x3 tensor in cartesian coords.
!!
!! NOTES
!! $cart(i,j)=G(i,a) G(j,b) frac(a,b)$
!! "Symmetric" storage mode for 3x3 tensor is 6 element array with
!! elements 11, 22, 33, 32, 31, and 21.
!! "cart" may be same array as "frac".
!! If rprimd transpose is input instead of gprimd, then convert tensor
!! in cartesian coordinates to reduced coordinates
!!
!! PARENTS
!!      d2frnl,mkcore,mkcore_paw,mkcore_wvl,nonlop_pl,nonlop_ylm,stresssym
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine strconv(frac,gprimd,cart)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strconv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: frac(6),gprimd(3,3)
 real(dp),intent(inout) :: cart(6) ! alias of frac   !vz_i

!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work1(3,3),work2(3,3)

! *************************************************************************

 work1(1,1)=frac(1)
 work1(2,2)=frac(2)
 work1(3,3)=frac(3)
 work1(3,2)=frac(4) ; work1(2,3)=frac(4)
 work1(3,1)=frac(5) ; work1(1,3)=frac(5)
 work1(2,1)=frac(6) ; work1(1,2)=frac(6)

 do ii=1,3
   work2(:,ii)=zero
   do jj=1,3
     work2(:,ii)=work2(:,ii)+gprimd(ii,jj)*work1(:,jj)
   end do
 end do

 do ii=1,3
   work1(ii,:)=zero
   do jj=1,3
     work1(ii,:)=work1(ii,:)+gprimd(ii,jj)*work2(jj,:)
   end do
 end do

 cart(1)=work1(1,1)
 cart(2)=work1(2,2)
 cart(3)=work1(3,3)
 cart(4)=work1(2,3)
 cart(5)=work1(1,3)
 cart(6)=work1(1,2)

end subroutine strconv
!!***

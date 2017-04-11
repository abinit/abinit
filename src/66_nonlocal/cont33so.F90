!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont33so
!! NAME
!! cont33so
!!
!! FUNCTION
!! Contract symmetric rank 3 tensor gxa1 with symmetric rank 3 tensor
!! gxa2 using metric tensor gmet and antisymmetric tensor amet to
!! produce rank 2 real tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa1(2,10)=rank 3 complex symmetric tensor
!!  gxa2(2,10)=rank 3 complex symmetric tensor
!!  gmet(3,3)=usual metric tensor (symmetric, real)
!!  amet(2,3,3)=antisymmetric complex tensor used for spin-orbit
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (pseudo-symmetric storage)
!!
!! NOTES
!! This contraction is used for spin-orbit correction in non-local
!! contribution to stresses.
!!
!! Symmetric gxa1, gxa2 are stored as
!!               111 221 331 321 311 211 222 332 322 333;
!! gmet(3,3) is symmetric but stored fully (9 elements);
!! amet(3,3) is antisymmetric but stored fully (9 elements);
!! Output rank2 is not symmetric but since
!!      $rank2_{gxa1,gxa2}(a,b)=conjg(rank2_{gxa2,gxa1}(b,a))$
!!       it is stored as 11 22 33 32 31 21.
!! Want 2*Re[contraction].
!!
!!{{\ \begin{equation}
!! rank2(a,b)=2 Re[15 r_{3A}(a,i,j) r_{3G}(b,j,i)-3 r_{3A}(a,b,i) r_1(i)]
!!\end{equation} }}
!!   where:
!!{{\ \begin{eqnarray}
!!   r_1(i)        & = & gxa2(i,j,k) gmet(j,k) \nonumber
!!   r_{3A}(i,j,k) & = & conjg(gxa1(p,i,j)) amet(p,k) \nonumber
!!   r_{3G}(i,j,k) & = & gxa2(p,i,j) gmet(p,k)
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


subroutine cont33so(gxa1,gxa2,gmet,amet,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont33so'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: amet(2,3,3),gmet(3,3),gxa1(2,10),gxa2(2,10)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: ii
!arrays
 real(dp) :: r1(2,3),r3A(2,18),r3G(2,18),s31(6),s33(6)

! *************************************************************************

!Compute r1(i)=gxa2(i,j,k)*gmet(j,k)

 do ii=1,2
   r1(ii,1)=gmet(1,1)*gxa2(ii,1)+gmet(2,2)*gxa2(ii,2)+&
&   gmet(3,3)*gxa2(ii,3)+2.d0*(&
&   gmet(3,2)*gxa2(ii,4)+gmet(3,1)*gxa2(ii,5)+&
&   gmet(2,1)*gxa2(ii,6))
   r1(ii,2)=gmet(1,1)*gxa2(ii,6)+gmet(2,2)*gxa2(ii,7)+&
&   gmet(3,3)*gxa2(ii,8)+2.d0*(&
&   gmet(3,2)*gxa2(ii,9)+gmet(3,1)*gxa2(ii,4)+&
&   gmet(2,1)*gxa2(ii,2))
   r1(ii,3)=gmet(1,1)*gxa2(ii,5)+gmet(2,2)*gxa2(ii,9)+&
&   gmet(3,3)*gxa2(ii,10)+2.d0*(&
&   gmet(3,2)*gxa2(ii,8)+gmet(3,1)*gxa2(ii,3)+&
&   gmet(2,1)*gxa2(ii,4))
 end do

!Compute r3G(i,j,k)=gxa2(p,i,j)*gmet(p,k)
!Write out components for 18 distinct terms, Re and Im
!(r3G is symmetric in first two indices, not in all permutations)
!Store as 111 221 331 321 311 211
!112 222 332 322 312 212
!113 223 333 323 313 213
 do ii=1,2
   r3G(ii, 1)=gmet(1,1)*gxa2(ii,1)+gmet(2,1)*gxa2(ii,6)+gmet(3,1)*gxa2(ii,5)
   r3G(ii, 2)=gmet(1,1)*gxa2(ii,2)+gmet(2,1)*gxa2(ii,7)+gmet(3,1)*gxa2(ii,9)
   r3G(ii, 3)=gmet(1,1)*gxa2(ii,3)+gmet(2,1)*gxa2(ii,8)+gmet(3,1)*gxa2(ii,10)
   r3G(ii, 4)=gmet(1,1)*gxa2(ii,4)+gmet(2,1)*gxa2(ii,9)+gmet(3,1)*gxa2(ii,8)
   r3G(ii, 5)=gmet(1,1)*gxa2(ii,5)+gmet(2,1)*gxa2(ii,4)+gmet(3,1)*gxa2(ii,3)
   r3G(ii, 6)=gmet(1,1)*gxa2(ii,6)+gmet(2,1)*gxa2(ii,2)+gmet(3,1)*gxa2(ii,4)
   r3G(ii, 7)=gmet(2,1)*gxa2(ii,1)+gmet(2,2)*gxa2(ii,6)+gmet(3,2)*gxa2(ii,5)
   r3G(ii, 8)=gmet(2,1)*gxa2(ii,2)+gmet(2,2)*gxa2(ii,7)+gmet(3,2)*gxa2(ii,9)
   r3G(ii, 9)=gmet(2,1)*gxa2(ii,3)+gmet(2,2)*gxa2(ii,8)+gmet(3,2)*gxa2(ii,10)
   r3G(ii,10)=gmet(2,1)*gxa2(ii,4)+gmet(2,2)*gxa2(ii,9)+gmet(3,2)*gxa2(ii,8)
   r3G(ii,11)=gmet(2,1)*gxa2(ii,5)+gmet(2,2)*gxa2(ii,4)+gmet(3,2)*gxa2(ii,3)
   r3G(ii,12)=gmet(2,1)*gxa2(ii,6)+gmet(2,2)*gxa2(ii,2)+gmet(3,2)*gxa2(ii,4)
   r3G(ii,13)=gmet(3,1)*gxa2(ii,1)+gmet(3,2)*gxa2(ii,6)+gmet(3,3)*gxa2(ii,5)
   r3G(ii,14)=gmet(3,1)*gxa2(ii,2)+gmet(3,2)*gxa2(ii,7)+gmet(3,3)*gxa2(ii,9)
   r3G(ii,15)=gmet(3,1)*gxa2(ii,3)+gmet(3,2)*gxa2(ii,8)+gmet(3,3)*gxa2(ii,10)
   r3G(ii,16)=gmet(3,1)*gxa2(ii,4)+gmet(3,2)*gxa2(ii,9)+gmet(3,3)*gxa2(ii,8)
   r3G(ii,17)=gmet(3,1)*gxa2(ii,5)+gmet(3,2)*gxa2(ii,4)+gmet(3,3)*gxa2(ii,3)
   r3G(ii,18)=gmet(3,1)*gxa2(ii,6)+gmet(3,2)*gxa2(ii,2)+gmet(3,3)*gxa2(ii,4)
 end do

!Compute r3A(i,j,k)=conjg(gxa1(p,i,j))*amet(p,k)
!Write out components for 18 distinct terms, Re and Im
!(r3A is symmetric in first two indices, not in all permutations)
!Store as 111 221 331 321 311 211
!112 222 332 322 312 212
!113 223 333 323 313 213
!Note that, since amet is antisymmetric, amet(i,i)=0

 r3A(re, 1)=amet(re,2,1)*gxa1(re,6 )+amet(im,2,1)*gxa1(im,6 )+&
& amet(re,3,1)*gxa1(re,5 )+amet(im,3,1)*gxa1(im,5 )
 r3A(re, 2)=amet(re,2,1)*gxa1(re,7 )+amet(im,2,1)*gxa1(im,7 )+&
& amet(re,3,1)*gxa1(re,9 )+amet(im,3,1)*gxa1(im,9 )
 r3A(re, 3)=amet(re,2,1)*gxa1(re,8 )+amet(im,2,1)*gxa1(im,8 )+&
& amet(re,3,1)*gxa1(re,10)+amet(im,3,1)*gxa1(im,10)
 r3A(re, 4)=amet(re,2,1)*gxa1(re,9 )+amet(im,2,1)*gxa1(im,9 )+&
& amet(re,3,1)*gxa1(re,8 )+amet(im,3,1)*gxa1(im,8 )
 r3A(re, 5)=amet(re,2,1)*gxa1(re,4 )+amet(im,2,1)*gxa1(im,4 )+&
& amet(re,3,1)*gxa1(re,3 )+amet(im,3,1)*gxa1(im,3 )
 r3A(re, 6)=amet(re,2,1)*gxa1(re,2 )+amet(im,2,1)*gxa1(im,2 )+&
& amet(re,3,1)*gxa1(re,4 )+amet(im,3,1)*gxa1(im,4 )
 r3A(re, 7)=amet(re,1,2)*gxa1(re,1 )+amet(im,1,2)*gxa1(im,1 )+&
& amet(re,3,2)*gxa1(re,5 )+amet(im,3,2)*gxa1(im,5 )
 r3A(re, 8)=amet(re,1,2)*gxa1(re,2 )+amet(im,1,2)*gxa1(im,2 )+&
& amet(re,3,2)*gxa1(re,9 )+amet(im,3,2)*gxa1(im,9 )
 r3A(re, 9)=amet(re,1,2)*gxa1(re,3 )+amet(im,1,2)*gxa1(im,3 )+&
& amet(re,3,2)*gxa1(re,10)+amet(im,3,2)*gxa1(im,10)
 r3A(re,10)=amet(re,1,2)*gxa1(re,4 )+amet(im,1,2)*gxa1(im,4 )+&
& amet(re,3,2)*gxa1(re,8 )+amet(im,3,2)*gxa1(im,8 )
 r3A(re,11)=amet(re,1,2)*gxa1(re,5 )+amet(im,1,2)*gxa1(im,5 )+&
& amet(re,3,2)*gxa1(re,3 )+amet(im,3,2)*gxa1(im,3 )
 r3A(re,12)=amet(re,1,2)*gxa1(re,6 )+amet(im,1,2)*gxa1(im,6 )+&
& amet(re,3,2)*gxa1(re,4 )+amet(im,3,2)*gxa1(im,4 )
 r3A(re,13)=amet(re,1,3)*gxa1(re,1 )+amet(im,1,3)*gxa1(im,1 )+&
& amet(re,2,3)*gxa1(re,6 )+amet(im,2,3)*gxa1(im,6 )
 r3A(re,14)=amet(re,1,3)*gxa1(re,2 )+amet(im,1,3)*gxa1(im,2 )+&
& amet(re,2,3)*gxa1(re,7 )+amet(im,2,3)*gxa1(im,7 )
 r3A(re,15)=amet(re,1,3)*gxa1(re,3 )+amet(im,1,3)*gxa1(im,3 )+&
& amet(re,2,3)*gxa1(re,8 )+amet(im,2,3)*gxa1(im,8 )
 r3A(re,16)=amet(re,1,3)*gxa1(re,4 )+amet(im,1,3)*gxa1(im,4 )+&
& amet(re,2,3)*gxa1(re,9 )+amet(im,2,3)*gxa1(im,9 )
 r3A(re,17)=amet(re,1,3)*gxa1(re,5 )+amet(im,1,3)*gxa1(im,5 )+&
& amet(re,2,3)*gxa1(re,4 )+amet(im,2,3)*gxa1(im,4 )
 r3A(re,18)=amet(re,1,3)*gxa1(re,6 )+amet(im,1,3)*gxa1(im,6 )+&
& amet(re,2,3)*gxa1(re,2 )+amet(im,2,3)*gxa1(im,2 )

 r3A(im, 1)=amet(im,2,1)*gxa1(re,6 )-amet(re,2,1)*gxa1(im,6 )+&
& amet(im,3,1)*gxa1(re,5 )-amet(re,3,1)*gxa1(im,5 )
 r3A(im, 2)=amet(im,2,1)*gxa1(re,7 )-amet(re,2,1)*gxa1(im,7 )+&
& amet(im,3,1)*gxa1(re,9 )-amet(re,3,1)*gxa1(im,9 )
 r3A(im, 3)=amet(im,2,1)*gxa1(re,8 )-amet(re,2,1)*gxa1(im,8 )+&
& amet(im,3,1)*gxa1(re,10)-amet(re,3,1)*gxa1(im,10)
 r3A(im, 4)=amet(im,2,1)*gxa1(re,9 )-amet(re,2,1)*gxa1(im,9 )+&
& amet(im,3,1)*gxa1(re,8 )-amet(re,3,1)*gxa1(im,8 )
 r3A(im, 5)=amet(im,2,1)*gxa1(re,4 )-amet(re,2,1)*gxa1(im,4 )+&
& amet(im,3,1)*gxa1(re,3 )-amet(re,3,1)*gxa1(im,3 )
 r3A(im, 6)=amet(im,2,1)*gxa1(re,2 )-amet(re,2,1)*gxa1(im,2 )+&
& amet(im,3,1)*gxa1(re,4 )-amet(re,3,1)*gxa1(im,4 )
 r3A(im, 7)=amet(im,1,2)*gxa1(re,1 )-amet(re,1,2)*gxa1(im,1 )+&
& amet(im,3,2)*gxa1(re,5 )-amet(re,3,2)*gxa1(im,5 )
 r3A(im, 8)=amet(im,1,2)*gxa1(re,2 )-amet(re,1,2)*gxa1(im,2 )+&
& amet(im,3,2)*gxa1(re,9 )-amet(re,3,2)*gxa1(im,9 )
 r3A(im, 9)=amet(im,1,2)*gxa1(re,3 )-amet(re,1,2)*gxa1(im,3 )+&
& amet(im,3,2)*gxa1(re,10)-amet(re,3,2)*gxa1(im,10)
 r3A(im,10)=amet(im,1,2)*gxa1(re,4 )-amet(re,1,2)*gxa1(im,4 )+&
& amet(im,3,2)*gxa1(re,8 )-amet(re,3,2)*gxa1(im,8 )
 r3A(im,11)=amet(im,1,2)*gxa1(re,5 )-amet(re,1,2)*gxa1(im,5 )+&
& amet(im,3,2)*gxa1(re,3 )-amet(re,3,2)*gxa1(im,3 )
 r3A(im,12)=amet(im,1,2)*gxa1(re,6 )-amet(re,1,2)*gxa1(im,6 )+&
& amet(im,3,2)*gxa1(re,4 )-amet(re,3,2)*gxa1(im,4 )
 r3A(im,13)=amet(im,1,3)*gxa1(re,1 )-amet(re,1,3)*gxa1(im,1 )+&
& amet(im,2,3)*gxa1(re,6 )-amet(re,2,3)*gxa1(im,6 )
 r3A(im,14)=amet(im,1,3)*gxa1(re,2 )-amet(re,1,3)*gxa1(im,2 )+&
& amet(im,2,3)*gxa1(re,7 )-amet(re,2,3)*gxa1(im,7 )
 r3A(im,15)=amet(im,1,3)*gxa1(re,3 )-amet(re,1,3)*gxa1(im,3 )+&
& amet(im,2,3)*gxa1(re,8 )-amet(re,2,3)*gxa1(im,8 )
 r3A(im,16)=amet(im,1,3)*gxa1(re,4 )-amet(re,1,3)*gxa1(im,4 )+&
& amet(im,2,3)*gxa1(re,9 )-amet(re,2,3)*gxa1(im,9 )
 r3A(im,17)=amet(im,1,3)*gxa1(re,5 )-amet(re,1,3)*gxa1(im,5 )+&
& amet(im,2,3)*gxa1(re,4 )-amet(re,2,3)*gxa1(im,4 )
 r3A(im,18)=amet(im,1,3)*gxa1(re,6 )-amet(re,1,3)*gxa1(im,6 )+&
& amet(im,2,3)*gxa1(re,2 )-amet(re,2,3)*gxa1(im,2 )

!Compute s33(a,b)=2*Re[r3A(a,i,j)*r3G(b,j,i)]

 s33(1)=2.d0*(r3A(re, 1)*r3G(re, 1)-r3A(im, 1)*r3G(im, 1)+&
& r3A(re,12)*r3G(re,12)-r3A(im,12)*r3G(im,12)+&
& r3A(re,17)*r3G(re,17)-r3A(im,17)*r3G(im,17)+&
& r3A(re,11)*r3G(re,18)-r3A(im,11)*r3G(im,18)+&
& r3A(re,18)*r3G(re,11)-r3A(im,18)*r3G(im,11)+&
& r3A(re, 5)*r3G(re,13)-r3A(im, 5)*r3G(im,13)+&
& r3A(re,13)*r3G(re, 5)-r3A(im,13)*r3G(im, 5)+&
& r3A(re, 6)*r3G(re, 7)-r3A(im, 6)*r3G(im, 7)+&
& r3A(re, 7)*r3G(re, 6)-r3A(im, 7)*r3G(im, 6))
 s33(2)=2.d0*(r3A(re, 6)*r3G(re, 6)-r3A(im, 6)*r3G(im, 6)+&
& r3A(re, 8)*r3G(re, 8)-r3A(im, 8)*r3G(im, 8)+&
& r3A(re,16)*r3G(re,16)-r3A(im,16)*r3G(im,16)+&
& r3A(re,10)*r3G(re,14)-r3A(im,10)*r3G(im,14)+&
& r3A(re,14)*r3G(re,10)-r3A(im,14)*r3G(im,10)+&
& r3A(re, 4)*r3G(re,18)-r3A(im, 4)*r3G(im,18)+&
& r3A(re,18)*r3G(re, 4)-r3A(im,18)*r3G(im, 4)+&
& r3A(re, 2)*r3G(re,12)-r3A(im, 2)*r3G(im,12)+&
& r3A(re,12)*r3G(re, 2)-r3A(im,12)*r3G(im, 2))
 s33(3)=2.d0*(r3A(re, 5)*r3G(re, 5)-r3A(im, 5)*r3G(im, 5)+&
& r3A(re,10)*r3G(re,10)-r3A(im,10)*r3G(im,10)+&
& r3A(re,15)*r3G(re,15)-r3A(im,15)*r3G(im,15)+&
& r3A(re, 9)*r3G(re,16)-r3A(im, 9)*r3G(im,16)+&
& r3A(re,16)*r3G(re, 9)-r3A(im,16)*r3G(im, 9)+&
& r3A(re, 3)*r3G(re,17)-r3A(im, 3)*r3G(im,17)+&
& r3A(re,17)*r3G(re, 3)-r3A(im,17)*r3G(im, 3)+&
& r3A(re, 4)*r3G(re,11)-r3A(im, 4)*r3G(im,11)+&
& r3A(re,11)*r3G(re, 4)-r3A(im,11)*r3G(im, 4))
 s33(4)=2.d0*(r3A(re, 5)*r3G(re, 6)-r3A(im, 5)*r3G(im, 6)+&
& r3A(re,10)*r3G(re, 8)-r3A(im,10)*r3G(im, 8)+&
& r3A(re,15)*r3G(re,16)-r3A(im,15)*r3G(im,16)+&
& r3A(re, 9)*r3G(re,14)-r3A(im, 9)*r3G(im,14)+&
& r3A(re,16)*r3G(re,10)-r3A(im,16)*r3G(im,10)+&
& r3A(re, 3)*r3G(re,18)-r3A(im, 3)*r3G(im,18)+&
& r3A(re,17)*r3G(re, 4)-r3A(im,17)*r3G(im, 4)+&
& r3A(re, 4)*r3G(re,12)-r3A(im, 4)*r3G(im,12)+&
& r3A(re,11)*r3G(re, 2)-r3A(im,11)*r3G(im, 2))
 s33(5)=2.d0*(r3A(re, 5)*r3G(re, 1)-r3A(im, 5)*r3G(im, 1)+&
& r3A(re,10)*r3G(re,12)-r3A(im,10)*r3G(im,12)+&
& r3A(re,15)*r3G(re,17)-r3A(im,15)*r3G(im,17)+&
& r3A(re, 9)*r3G(re,18)-r3A(im, 9)*r3G(im,18)+&
& r3A(re,16)*r3G(re,11)-r3A(im,16)*r3G(im,11)+&
& r3A(re, 3)*r3G(re,13)-r3A(im, 3)*r3G(im,13)+&
& r3A(re,17)*r3G(re, 5)-r3A(im,17)*r3G(im, 5)+&
& r3A(re, 4)*r3G(re, 7)-r3A(im, 4)*r3G(im, 7)+&
& r3A(re,11)*r3G(re, 6)-r3A(im,11)*r3G(im, 6))
 s33(6)=2.d0*(r3A(re, 6)*r3G(re, 1)-r3A(im, 6)*r3G(im, 1)+&
& r3A(re, 8)*r3G(re,12)-r3A(im, 8)*r3G(im,12)+&
& r3A(re,16)*r3G(re,17)-r3A(im,16)*r3G(im,17)+&
& r3A(re,10)*r3G(re,18)-r3A(im,10)*r3G(im,18)+&
& r3A(re,14)*r3G(re,11)-r3A(im,14)*r3G(im,11)+&
& r3A(re, 4)*r3G(re,13)-r3A(im, 4)*r3G(im,13)+&
& r3A(re,18)*r3G(re, 5)-r3A(im,18)*r3G(im, 5)+&
& r3A(re, 2)*r3G(re, 7)-r3A(im, 2)*r3G(im, 7)+&
& r3A(re,12)*r3G(re, 6)-r3A(im,12)*r3G(im, 6))

!Compute s31(a,b)=2*Re[r3A(a,b,i)*r1(i)]

 s31(1)=2.d0*(r1(re,1)*r3A(re, 1)-r1(im,1)*r3A(im, 1)+&
& r1(re,2)*r3A(re, 7)-r1(im,2)*r3A(im, 7)+&
& r1(re,3)*r3A(re,13)-r1(im,3)*r3A(im,13))
 s31(2)=2.d0*(r1(re,1)*r3A(re, 2)-r1(im,1)*r3A(im, 2)+&
& r1(re,2)*r3A(re, 8)-r1(im,2)*r3A(im, 8)+&
& r1(re,3)*r3A(re,14)-r1(im,3)*r3A(im,14))
 s31(3)=2.d0*(r1(re,1)*r3A(re, 3)-r1(im,1)*r3A(im, 3)+&
& r1(re,2)*r3A(re, 9)-r1(im,2)*r3A(im, 9)+&
& r1(re,3)*r3A(re,15)-r1(im,3)*r3A(im,15))
 s31(4)=2.d0*(r1(re,1)*r3A(re, 4)-r1(im,1)*r3A(im, 4)+&
& r1(re,2)*r3A(re,10)-r1(im,2)*r3A(im,10)+&
& r1(re,3)*r3A(re,16)-r1(im,3)*r3A(im,16))
 s31(5)=2.d0*(r1(re,1)*r3A(re, 5)-r1(im,1)*r3A(im, 5)+&
& r1(re,2)*r3A(re,11)-r1(im,2)*r3A(im,11)+&
& r1(re,3)*r3A(re,17)-r1(im,3)*r3A(im,17))
 s31(6)=2.d0*(r1(re,1)*r3A(re, 6)-r1(im,1)*r3A(im, 6)+&
& r1(re,2)*r3A(re,12)-r1(im,2)*r3A(im,12)+&
& r1(re,3)*r3A(re,18)-r1(im,3)*r3A(im,18))

!Finally, compute rank2(a,b)=-15*s33(a,b)+3*s31(a,b)

 rank2(:)=15.d0*s33(:)-3.d0*s31(:)

end subroutine cont33so
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont3
!! NAME
!! cont3
!!
!! FUNCTION
!! Compute several specialized contractions needed for the
!! l=3 part of the stress tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa(2,10)=complex symmetric rank 3 tensor
!!  gmet(3,3)=usual metric tensor, a symmetric matrix stored in
!!            full storage mode (bohr^-2)
!!
!! OUTPUT
!!  rank2(6)=2*Re[contraction] given by
!!   2*Re[(15/2)*r3(a,i,j)*r3(b,j,i)-3*r1(i)*r3(a,b,i)-(3/2)*r1(a)*r1(b)]
!!   where r3(a,i,j)=gmet(j,k) gxa(a,i,k) and r1(a)=gmet(i,j) gxa(i,j,a).
!!  rank2 is stored in the compressed form 11 22 33 32 31 21.
!!
!! NOTES
!! Input gxa is a completely symmetric rank 3 tensor (complex)
!! in compressed storage: 111 221 331 321 311 211 222 332 322 333.
!! The output tensor is completely symmetric rank 2, real, and is given by
!!  $2 Re[{15 \over 2} r3(a,i,j) r3(b,j,i) - 3 r1(i) r3(a,b,i) - {3 \over 2} r1(a) r1(b)]$
!!  where $r3(a,i,j)=gmet(j,k) gxa(a,i,k)$ and $r1(a)=gmet(i,j) gxa(i,j,a)$.
!!  rank2 is stored in the compressed form 11 22 33 32 31 21.
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


subroutine cont3(gxa,gmet,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont3'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa(2,10)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: ii
!arrays
 real(dp) :: r1(2,3),r3(2,18),s13(6),s33(6)

! *************************************************************************

!Compute r1(a) = gmet(i,j) gxa(i,j,a)

!Write out components for 3 distinct terms, Re and Im
 do ii=1,2
   r1(ii,1)=gmet(1,1)*gxa(ii,1)+gmet(2,2)*gxa(ii,2)+&
&   gmet(3,3)*gxa(ii,3)+2.d0*(&
&   gmet(3,2)*gxa(ii,4)+gmet(3,1)*gxa(ii,5)+&
&   gmet(2,1)*gxa(ii,6))
   r1(ii,2)=gmet(1,1)*gxa(ii,6)+gmet(2,2)*gxa(ii,7)+&
&   gmet(3,3)*gxa(ii,8)+2.d0*(&
&   gmet(3,2)*gxa(ii,9)+gmet(3,1)*gxa(ii,4)+&
&   gmet(2,1)*gxa(ii,2))
   r1(ii,3)=gmet(1,1)*gxa(ii,5)+gmet(2,2)*gxa(ii,9)+&
&   gmet(3,3)*gxa(ii,10)+2.d0*(&
&   gmet(3,2)*gxa(ii,8)+gmet(3,1)*gxa(ii,3)+&
&   gmet(2,1)*gxa(ii,4))
 end do

!Compute r3(a,b,k)=gmet(k,n) gxa(a,b,n)

!Write out components for 18 distinct terms, Re and Im
!(symmetric in first two indices, not in all permutations)
!store as 111 221 331 321 311 211
!112 222 332 322 312 212
!113 223 333 323 313 213
 do ii=1,2
   r3(ii, 1)=gmet(1,1)*gxa(ii,1)+gmet(2,1)*gxa(ii,6)+&
&   gmet(3,1)*gxa(ii,5)
   r3(ii, 2)=gmet(1,1)*gxa(ii,2)+gmet(2,1)*gxa(ii,7)+&
&   gmet(3,1)*gxa(ii,9)
   r3(ii, 3)=gmet(1,1)*gxa(ii,3)+gmet(2,1)*gxa(ii,8)+&
&   gmet(3,1)*gxa(ii,10)
   r3(ii, 4)=gmet(1,1)*gxa(ii,4)+gmet(2,1)*gxa(ii,9)+&
&   gmet(3,1)*gxa(ii,8)
   r3(ii, 5)=gmet(1,1)*gxa(ii,5)+gmet(2,1)*gxa(ii,4)+&
&   gmet(3,1)*gxa(ii,3)
   r3(ii, 6)=gmet(1,1)*gxa(ii,6)+gmet(2,1)*gxa(ii,2)+&
&   gmet(3,1)*gxa(ii,4)
   r3(ii, 7)=gmet(2,1)*gxa(ii,1)+gmet(2,2)*gxa(ii,6)+&
&   gmet(3,2)*gxa(ii,5)
   r3(ii, 8)=gmet(2,1)*gxa(ii,2)+gmet(2,2)*gxa(ii,7)+&
&   gmet(3,2)*gxa(ii,9)
   r3(ii, 9)=gmet(2,1)*gxa(ii,3)+gmet(2,2)*gxa(ii,8)+&
&   gmet(3,2)*gxa(ii,10)
   r3(ii,10)=gmet(2,1)*gxa(ii,4)+gmet(2,2)*gxa(ii,9)+&
&   gmet(3,2)*gxa(ii,8)
   r3(ii,11)=gmet(2,1)*gxa(ii,5)+gmet(2,2)*gxa(ii,4)+&
&   gmet(3,2)*gxa(ii,3)
   r3(ii,12)=gmet(2,1)*gxa(ii,6)+gmet(2,2)*gxa(ii,2)+&
&   gmet(3,2)*gxa(ii,4)
   r3(ii,13)=gmet(3,1)*gxa(ii,1)+gmet(3,2)*gxa(ii,6)+&
&   gmet(3,3)*gxa(ii,5)
   r3(ii,14)=gmet(3,1)*gxa(ii,2)+gmet(3,2)*gxa(ii,7)+&
&   gmet(3,3)*gxa(ii,9)
   r3(ii,15)=gmet(3,1)*gxa(ii,3)+gmet(3,2)*gxa(ii,8)+&
&   gmet(3,3)*gxa(ii,10)
   r3(ii,16)=gmet(3,1)*gxa(ii,4)+gmet(3,2)*gxa(ii,9)+&
&   gmet(3,3)*gxa(ii,8)
   r3(ii,17)=gmet(3,1)*gxa(ii,5)+gmet(3,2)*gxa(ii,4)+&
&   gmet(3,3)*gxa(ii,3)
   r3(ii,18)=gmet(3,1)*gxa(ii,6)+gmet(3,2)*gxa(ii,2)+&
&   gmet(3,3)*gxa(ii,4)

 end do

!Now need
!2*Re[(15/2)*r3(a,i,j)*r3(b,j,i)-3*r1(i)*r3(a,b,i)-(3/2)*r1(a)*r1(b)].

!Write out s33(a,b)=2*Re[r3(a,i,j)*r3(b,j,i)]

 s33(1)=2.d0*(r3(re, 1)*r3(re, 1)+r3(im, 1)*r3(im, 1)+&
& r3(re,12)*r3(re,12)+r3(im,12)*r3(im,12)+&
& r3(re,17)*r3(re,17)+r3(im,17)*r3(im,17)+&
& r3(re,11)*r3(re,18)+r3(im,11)*r3(im,18)+&
& r3(re,18)*r3(re,11)+r3(im,18)*r3(im,11)+&
& r3(re, 5)*r3(re,13)+r3(im, 5)*r3(im,13)+&
& r3(re,13)*r3(re, 5)+r3(im,13)*r3(im, 5)+&
& r3(re, 6)*r3(re, 7)+r3(im, 6)*r3(im, 7)+&
& r3(re, 7)*r3(re, 6)+r3(im, 7)*r3(im, 6))

 s33(2)=2.d0*(r3(re, 6)*r3(re, 6)+r3(im, 6)*r3(im, 6)+&
& r3(re, 8)*r3(re, 8)+r3(im, 8)*r3(im, 8)+&
& r3(re,16)*r3(re,16)+r3(im,16)*r3(im,16)+&
& r3(re,10)*r3(re,14)+r3(im,10)*r3(im,14)+&
& r3(re,14)*r3(re,10)+r3(im,14)*r3(im,10)+&
& r3(re, 4)*r3(re,18)+r3(im, 4)*r3(im,18)+&
& r3(re,18)*r3(re, 4)+r3(im,18)*r3(im, 4)+&
& r3(re, 2)*r3(re,12)+r3(im, 2)*r3(im,12)+&
& r3(re,12)*r3(re, 2)+r3(im,12)*r3(im, 2))

 s33(3)=2.d0*(r3(re, 5)*r3(re, 5)+r3(im, 5)*r3(im, 5)+&
& r3(re,10)*r3(re,10)+r3(im,10)*r3(im,10)+&
& r3(re,15)*r3(re,15)+r3(im,15)*r3(im,15)+&
& r3(re, 9)*r3(re,16)+r3(im, 9)*r3(im,16)+&
& r3(re,16)*r3(re, 9)+r3(im,16)*r3(im, 9)+&
& r3(re, 3)*r3(re,17)+r3(im, 3)*r3(im,17)+&
& r3(re,17)*r3(re, 3)+r3(im,17)*r3(im, 3)+&
& r3(re, 4)*r3(re,11)+r3(im, 4)*r3(im,11)+&
& r3(re,11)*r3(re, 4)+r3(im,11)*r3(im, 4))

 s33(4)=2.d0*(r3(re, 5)*r3(re, 6)+r3(im, 5)*r3(im, 6)+&
& r3(re,10)*r3(re, 8)+r3(im,10)*r3(im, 8)+&
& r3(re,15)*r3(re,16)+r3(im,15)*r3(im,16)+&
& r3(re, 9)*r3(re,14)+r3(im, 9)*r3(im,14)+&
& r3(re,16)*r3(re,10)+r3(im,16)*r3(im,10)+&
& r3(re, 3)*r3(re,18)+r3(im, 3)*r3(im,18)+&
& r3(re,17)*r3(re, 4)+r3(im,17)*r3(im, 4)+&
& r3(re, 4)*r3(re,12)+r3(im, 4)*r3(im,12)+&
& r3(re,11)*r3(re, 2)+r3(im,11)*r3(im, 2))

 s33(5)=2.d0*(r3(re, 5)*r3(re, 1)+r3(im, 5)*r3(im, 1)+&
& r3(re,10)*r3(re,12)+r3(im,10)*r3(im,12)+&
& r3(re,15)*r3(re,17)+r3(im,15)*r3(im,17)+&
& r3(re, 9)*r3(re,18)+r3(im, 9)*r3(im,18)+&
& r3(re,16)*r3(re,11)+r3(im,16)*r3(im,11)+&
& r3(re, 3)*r3(re,13)+r3(im, 3)*r3(im,13)+&
& r3(re,17)*r3(re, 5)+r3(im,17)*r3(im, 5)+&
& r3(re, 4)*r3(re, 7)+r3(im, 4)*r3(im, 7)+&
& r3(re,11)*r3(re, 6)+r3(im,11)*r3(im, 6))

 s33(6)=2.d0*(r3(re, 6)*r3(re, 1)+r3(im, 6)*r3(im, 1)+&
& r3(re, 8)*r3(re,12)+r3(im, 8)*r3(im,12)+&
& r3(re,16)*r3(re,17)+r3(im,16)*r3(im,17)+&
& r3(re,10)*r3(re,18)+r3(im,10)*r3(im,18)+&
& r3(re,14)*r3(re,11)+r3(im,14)*r3(im,11)+&
& r3(re, 4)*r3(re,13)+r3(im, 4)*r3(im,13)+&
& r3(re,18)*r3(re, 5)+r3(im,18)*r3(im, 5)+&
& r3(re, 2)*r3(re, 7)+r3(im, 2)*r3(im, 7)+&
& r3(re,12)*r3(re, 6)+r3(im,12)*r3(im, 6))


!Write out s13(a,b)=2*Re[r1(i)*r3(a,b,i)]

 s13(1)=2.d0*(r1(re,1)*r3(re, 1)+r1(im,1)*r3(im, 1)+&
& r1(re,2)*r3(re, 7)+r1(im,2)*r3(im, 7)+&
& r1(re,3)*r3(re,13)+r1(im,3)*r3(im,13))
 s13(2)=2.d0*(r1(re,1)*r3(re, 2)+r1(im,1)*r3(im, 2)+&
& r1(re,2)*r3(re, 8)+r1(im,2)*r3(im, 8)+&
& r1(re,3)*r3(re,14)+r1(im,3)*r3(im,14))
 s13(3)=2.d0*(r1(re,1)*r3(re, 3)+r1(im,1)*r3(im, 3)+&
& r1(re,2)*r3(re, 9)+r1(im,2)*r3(im, 9)+&
& r1(re,3)*r3(re,15)+r1(im,3)*r3(im,15))
 s13(4)=2.d0*(r1(re,1)*r3(re, 4)+r1(im,1)*r3(im, 4)+&
& r1(re,2)*r3(re,10)+r1(im,2)*r3(im,10)+&
& r1(re,3)*r3(re,16)+r1(im,3)*r3(im,16))
 s13(5)=2.d0*(r1(re,1)*r3(re, 5)+r1(im,1)*r3(im, 5)+&
& r1(re,2)*r3(re,11)+r1(im,2)*r3(im,11)+&
& r1(re,3)*r3(re,17)+r1(im,3)*r3(im,17))
 s13(6)=2.d0*(r1(re,1)*r3(re, 6)+r1(im,1)*r3(im, 6)+&
& r1(re,2)*r3(re,12)+r1(im,2)*r3(im,12)+&
& r1(re,3)*r3(re,18)+r1(im,3)*r3(im,18))

!Finally, write out the six terms as final answer
!rank2(a,b)=(15/2)*s33(a,b)-3*s13(a,b)-(3/2)*2*Re[r1(a)*r1(b)]

 rank2(1)=7.5d0*s33(1)-3.d0*s13(1)&
& -3.d0*(r1(re,1)*r1(re,1)+r1(im,1)*r1(im,1))
 rank2(2)=7.5d0*s33(2)-3.d0*s13(2)&
& -3.d0*(r1(re,2)*r1(re,2)+r1(im,2)*r1(im,2))
 rank2(3)=7.5d0*s33(3)-3.d0*s13(3)&
& -3.d0*(r1(re,3)*r1(re,3)+r1(im,3)*r1(im,3))
 rank2(4)=7.5d0*s33(4)-3.d0*s13(4)&
& -3.d0*(r1(re,3)*r1(re,2)+r1(im,3)*r1(im,2))
 rank2(5)=7.5d0*s33(5)-3.d0*s13(5)&
& -3.d0*(r1(re,3)*r1(re,1)+r1(im,3)*r1(im,1))
 rank2(6)=7.5d0*s33(6)-3.d0*s13(6)&
& -3.d0*(r1(re,2)*r1(re,1)+r1(im,2)*r1(im,1))

end subroutine cont3
!!***

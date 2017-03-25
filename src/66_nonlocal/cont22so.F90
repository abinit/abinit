!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont22so
!! NAME
!! cont22so
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa1 with symmetric rank 2 tensor
!! gxa2 using antisymmetric tensor amet to produce rank 2 real tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa1(2,6)=rank 2 complex symmetric tensor
!!  gxa2(2,6)=rank 2 complex symmetric tensor
!!  amet(2,3,3)=antisymmetric complex tensor used for spin-orbit
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (pseudo-symmetric storage)
!!
!! NOTES
!! This contraction is used for spin-orbit correction in non-local
!! contribution to stresses.
!!
!! Symmetric gxa1, gxa2 are stored as 11 22 33 32 31 21;
!! amet(3,3) is antisymmetric but stored fully (9 elements);
!! Output rank2 is not symmetric but since
!!      $rank2_{gxa1,gxa2}(a,b)=conjg(rank2_{gxa2,gxa1}(b,a))$
!!       it is stored as 11 22 33 32 31 21.
!! Want 2*Re[contraction].
!!
!!{{\ \begin{equation}
!! rank2(a,b)=2 Re[conjg(gxa1(i,a)) amet(i,j) gxa2(j,b)]
!!\end{equation} }}
!!
!! Note that, since amet is antisymmetric, amet(i,i)=0
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


subroutine cont22so(gxa1,gxa2,amet,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont22so'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: amet(2,3,3),gxa1(2,6),gxa2(2,6)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1

! *************************************************************************

!Simply write out index summations
!a=1, b=1 in rank2(a,b) --> maps to index 1
 rank2(1)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,5)*gxa2(re,6)+gxa1(im,5)*gxa2(im,6))-&
& amet(2,3,2)*(gxa1(re,5)*gxa2(im,6)-gxa1(im,5)*gxa2(re,6))+&
& amet(1,3,1)*(gxa1(re,5)*gxa2(re,1)+gxa1(im,5)*gxa2(im,1))-&
& amet(2,3,1)*(gxa1(re,5)*gxa2(im,1)-gxa1(im,5)*gxa2(re,1))+&
& amet(1,2,1)*(gxa1(re,6)*gxa2(re,1)+gxa1(im,6)*gxa2(im,1))-&
& amet(2,2,1)*(gxa1(re,6)*gxa2(im,1)-gxa1(im,6)*gxa2(re,1))+&
& amet(1,2,3)*(gxa1(re,6)*gxa2(re,5)+gxa1(im,6)*gxa2(im,5))-&
& amet(2,2,3)*(gxa1(re,6)*gxa2(im,5)-gxa1(im,6)*gxa2(re,5))+&
& amet(1,1,3)*(gxa1(re,1)*gxa2(re,5)+gxa1(im,1)*gxa2(im,5))-&
& amet(2,1,3)*(gxa1(re,1)*gxa2(im,5)-gxa1(im,1)*gxa2(re,5))+&
& amet(1,1,2)*(gxa1(re,1)*gxa2(re,6)+gxa1(im,1)*gxa2(im,6))-&
& amet(2,1,2)*(gxa1(re,1)*gxa2(im,6)-gxa1(im,1)*gxa2(re,6)))

!a=2, b=2 in rank2(a,b) --> maps to index 2
 rank2(2)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,4)*gxa2(re,2)+gxa1(im,4)*gxa2(im,2))-&
& amet(2,3,2)*(gxa1(re,4)*gxa2(im,2)-gxa1(im,4)*gxa2(re,2))+&
& amet(1,3,1)*(gxa1(re,4)*gxa2(re,6)+gxa1(im,4)*gxa2(im,6))-&
& amet(2,3,1)*(gxa1(re,4)*gxa2(im,6)-gxa1(im,4)*gxa2(re,6))+&
& amet(1,2,1)*(gxa1(re,2)*gxa2(re,6)+gxa1(im,2)*gxa2(im,6))-&
& amet(2,2,1)*(gxa1(re,2)*gxa2(im,6)-gxa1(im,2)*gxa2(re,6))+&
& amet(1,2,3)*(gxa1(re,2)*gxa2(re,4)+gxa1(im,2)*gxa2(im,4))-&
& amet(2,2,3)*(gxa1(re,2)*gxa2(im,4)-gxa1(im,2)*gxa2(re,4))+&
& amet(1,1,3)*(gxa1(re,6)*gxa2(re,4)+gxa1(im,6)*gxa2(im,4))-&
& amet(2,1,3)*(gxa1(re,6)*gxa2(im,4)-gxa1(im,6)*gxa2(re,4))+&
& amet(1,1,2)*(gxa1(re,6)*gxa2(re,2)+gxa1(im,6)*gxa2(im,2))-&
& amet(2,1,2)*(gxa1(re,6)*gxa2(im,2)-gxa1(im,6)*gxa2(re,2)))

!a=3, b=3 in rank2(a,b) --> maps to index 3
 rank2(3)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,3)*gxa2(re,4)+gxa1(im,3)*gxa2(im,4))-&
& amet(2,3,2)*(gxa1(re,3)*gxa2(im,4)-gxa1(im,3)*gxa2(re,4))+&
& amet(1,3,1)*(gxa1(re,3)*gxa2(re,5)+gxa1(im,3)*gxa2(im,5))-&
& amet(2,3,1)*(gxa1(re,3)*gxa2(im,5)-gxa1(im,3)*gxa2(re,5))+&
& amet(1,2,1)*(gxa1(re,4)*gxa2(re,5)+gxa1(im,4)*gxa2(im,5))-&
& amet(2,2,1)*(gxa1(re,4)*gxa2(im,5)-gxa1(im,4)*gxa2(re,5))+&
& amet(1,2,3)*(gxa1(re,4)*gxa2(re,3)+gxa1(im,4)*gxa2(im,3))-&
& amet(2,2,3)*(gxa1(re,4)*gxa2(im,3)-gxa1(im,4)*gxa2(re,3))+&
& amet(1,1,3)*(gxa1(re,5)*gxa2(re,3)+gxa1(im,5)*gxa2(im,3))-&
& amet(2,1,3)*(gxa1(re,5)*gxa2(im,3)-gxa1(im,5)*gxa2(re,3))+&
& amet(1,1,2)*(gxa1(re,5)*gxa2(re,4)+gxa1(im,5)*gxa2(im,4))-&
& amet(2,1,2)*(gxa1(re,5)*gxa2(im,4)-gxa1(im,5)*gxa2(re,4)))

!a=3, b=2 in rank2(a,b) --> maps to index 4
 rank2(4)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,3)*gxa2(re,2)+gxa1(im,3)*gxa2(im,2))-&
& amet(2,3,2)*(gxa1(re,3)*gxa2(im,2)-gxa1(im,3)*gxa2(re,2))+&
& amet(1,3,1)*(gxa1(re,3)*gxa2(re,6)+gxa1(im,3)*gxa2(im,6))-&
& amet(2,3,1)*(gxa1(re,3)*gxa2(im,6)-gxa1(im,3)*gxa2(re,6))+&
& amet(1,2,1)*(gxa1(re,4)*gxa2(re,6)+gxa1(im,4)*gxa2(im,6))-&
& amet(2,2,1)*(gxa1(re,4)*gxa2(im,6)-gxa1(im,4)*gxa2(re,6))+&
& amet(1,2,3)*(gxa1(re,4)*gxa2(re,4)+gxa1(im,4)*gxa2(im,4))-&
& amet(2,2,3)*(gxa1(re,4)*gxa2(im,4)-gxa1(im,4)*gxa2(re,4))+&
& amet(1,1,3)*(gxa1(re,5)*gxa2(re,4)+gxa1(im,5)*gxa2(im,4))-&
& amet(2,1,3)*(gxa1(re,5)*gxa2(im,4)-gxa1(im,5)*gxa2(re,4))+&
& amet(1,1,2)*(gxa1(re,5)*gxa2(re,2)+gxa1(im,5)*gxa2(im,2))-&
& amet(2,1,2)*(gxa1(re,5)*gxa2(im,2)-gxa1(im,5)*gxa2(re,2)))

!a=3, b=1 in rank2(a,b) --> maps to index 5
 rank2(5)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,3)*gxa2(re,6)+gxa1(im,3)*gxa2(im,6))-&
& amet(2,3,2)*(gxa1(re,3)*gxa2(im,6)-gxa1(im,3)*gxa2(re,6))+&
& amet(1,3,1)*(gxa1(re,3)*gxa2(re,1)+gxa1(im,3)*gxa2(im,1))-&
& amet(2,3,1)*(gxa1(re,3)*gxa2(im,1)-gxa1(im,3)*gxa2(re,1))+&
& amet(1,2,1)*(gxa1(re,4)*gxa2(re,1)+gxa1(im,4)*gxa2(im,1))-&
& amet(2,2,1)*(gxa1(re,4)*gxa2(im,1)-gxa1(im,4)*gxa2(re,1))+&
& amet(1,2,3)*(gxa1(re,4)*gxa2(re,5)+gxa1(im,4)*gxa2(im,5))-&
& amet(2,2,3)*(gxa1(re,4)*gxa2(im,5)-gxa1(im,4)*gxa2(re,5))+&
& amet(1,1,3)*(gxa1(re,5)*gxa2(re,5)+gxa1(im,5)*gxa2(im,5))-&
& amet(2,1,3)*(gxa1(re,5)*gxa2(im,5)-gxa1(im,5)*gxa2(re,5))+&
& amet(1,1,2)*(gxa1(re,5)*gxa2(re,6)+gxa1(im,5)*gxa2(im,6))-&
& amet(2,1,2)*(gxa1(re,5)*gxa2(im,6)-gxa1(im,5)*gxa2(re,6)))

!a=2, b=1 in rank2(a,b) --> maps to index 6
 rank2(6)=2.0d0*(&
& amet(1,3,2)*(gxa1(re,4)*gxa2(re,6)+gxa1(im,4)*gxa2(im,6))-&
& amet(2,3,2)*(gxa1(re,4)*gxa2(im,6)-gxa1(im,4)*gxa2(re,6))+&
& amet(1,3,1)*(gxa1(re,4)*gxa2(re,1)+gxa1(im,4)*gxa2(im,1))-&
& amet(2,3,1)*(gxa1(re,4)*gxa2(im,1)-gxa1(im,4)*gxa2(re,1))+&
& amet(1,2,1)*(gxa1(re,2)*gxa2(re,1)+gxa1(im,2)*gxa2(im,1))-&
& amet(2,2,1)*(gxa1(re,2)*gxa2(im,1)-gxa1(im,2)*gxa2(re,1))+&
& amet(1,2,3)*(gxa1(re,2)*gxa2(re,5)+gxa1(im,2)*gxa2(im,5))-&
& amet(2,2,3)*(gxa1(re,2)*gxa2(im,5)-gxa1(im,2)*gxa2(re,5))+&
& amet(1,1,3)*(gxa1(re,6)*gxa2(re,5)+gxa1(im,6)*gxa2(im,5))-&
& amet(2,1,3)*(gxa1(re,6)*gxa2(im,5)-gxa1(im,6)*gxa2(re,5))+&
& amet(1,1,2)*(gxa1(re,6)*gxa2(re,6)+gxa1(im,6)*gxa2(im,6))-&
& amet(2,1,2)*(gxa1(re,6)*gxa2(im,6)-gxa1(im,6)*gxa2(re,6)))

end subroutine cont22so
!!***

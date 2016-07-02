!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont22
!! NAME
!! cont22
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa with itself using gmet to
!! produce symmetric rank 2 tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DXA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa(2,6)=rank 2 complex tensor
!!  gmet(3,3)=real symmetric metric tensor (full storage)
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (symmetric storage)
!!
!! NOTES
!! Symmetric gxa is stored as 11 22 33 32 31 21;
!! gmet(3,3) is symmetric but stored fully (9 elements);
!! output rank2 is stored as 11 22 33 32 31 21.
!! Want $2 Re[contraction]$.
!! $rank2(a,b)=2 Re[gxa(i,a)^"*" gmet(i,j) gxa(j,b)]$.
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


subroutine cont22(gxa,gmet,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont22'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa(2,6)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1

! *************************************************************************

!Simply write out index summations
!a=1, b=1 in rank2(a,b) --> maps to index 1
 rank2(1)=2.0d0*(&
& gmet(1,1)*(gxa(re,1)*gxa(re,1)+gxa(im,1)*gxa(im,1))+&
& gmet(2,2)*(gxa(re,6)*gxa(re,6)+gxa(im,6)*gxa(im,6))+&
& gmet(3,3)*(gxa(re,5)*gxa(re,5)+gxa(im,5)*gxa(im,5))+&
& 2.0d0*(&
& gmet(3,2)*(gxa(re,5)*gxa(re,6)+gxa(im,5)*gxa(im,6))+&
& gmet(3,1)*(gxa(re,5)*gxa(re,1)+gxa(im,5)*gxa(im,1))+&
& gmet(2,1)*(gxa(re,6)*gxa(re,1)+gxa(im,6)*gxa(im,1))))

!a=2, b=2 in rank2(a,b) --> maps to index 2
 rank2(2)=2.0d0*(&
& gmet(1,1)*(gxa(re,6)*gxa(re,6)+gxa(im,6)*gxa(im,6))+&
& gmet(2,2)*(gxa(re,2)*gxa(re,2)+gxa(im,2)*gxa(im,2))+&
& gmet(3,3)*(gxa(re,4)*gxa(re,4)+gxa(im,4)*gxa(im,4))+&
& 2.0d0*(&
& gmet(3,2)*(gxa(re,4)*gxa(re,2)+gxa(im,4)*gxa(im,2))+&
& gmet(3,1)*(gxa(re,4)*gxa(re,6)+gxa(im,4)*gxa(im,6))+&
& gmet(2,1)*(gxa(re,2)*gxa(re,6)+gxa(im,2)*gxa(im,6))))

!a=3, b=3 in rank2(a,b) --> maps to index 3
 rank2(3)=2.0d0*(&
& gmet(1,1)*(gxa(re,5)*gxa(re,5)+gxa(im,5)*gxa(im,5))+&
& gmet(2,2)*(gxa(re,4)*gxa(re,4)+gxa(im,4)*gxa(im,4))+&
& gmet(3,3)*(gxa(re,3)*gxa(re,3)+gxa(im,3)*gxa(im,3))+&
& 2.0d0*(&
& gmet(3,2)*(gxa(re,4)*gxa(re,3)+gxa(im,4)*gxa(im,3))+&
& gmet(3,1)*(gxa(re,5)*gxa(re,3)+gxa(im,5)*gxa(im,3))+&
& gmet(2,1)*(gxa(re,5)*gxa(re,4)+gxa(im,5)*gxa(im,4))))

!a=3, b=2 in rank2(a,b) --> maps to index 4
 rank2(4)=2.0d0*(&
& gmet(1,1)*(gxa(re,5)*gxa(re,6)+gxa(im,5)*gxa(im,6))+&
& gmet(2,2)*(gxa(re,4)*gxa(re,2)+gxa(im,4)*gxa(im,2))+&
& gmet(3,3)*(gxa(re,3)*gxa(re,4)+gxa(im,3)*gxa(im,4))+&
& gmet(3,2)*(gxa(re,3)*gxa(re,2)+gxa(im,3)*gxa(im,2))+&
& gmet(3,1)*(gxa(re,3)*gxa(re,6)+gxa(im,3)*gxa(im,6))+&
& gmet(2,1)*(gxa(re,4)*gxa(re,6)+gxa(im,4)*gxa(im,6))+&
& gmet(2,3)*(gxa(re,4)*gxa(re,4)+gxa(im,4)*gxa(im,4))+&
& gmet(1,3)*(gxa(re,5)*gxa(re,4)+gxa(im,5)*gxa(im,4))+&
& gmet(1,2)*(gxa(re,5)*gxa(re,2)+gxa(im,5)*gxa(im,2)))

!a=3, b=1 in rank2(a,b) --> maps to index 5
 rank2(5)=2.0d0*(&
& gmet(1,1)*(gxa(re,5)*gxa(re,1)+gxa(im,5)*gxa(im,1))+&
& gmet(2,2)*(gxa(re,4)*gxa(re,6)+gxa(im,4)*gxa(im,6))+&
& gmet(3,3)*(gxa(re,3)*gxa(re,5)+gxa(im,3)*gxa(im,5))+&
& gmet(3,2)*(gxa(re,3)*gxa(re,6)+gxa(im,3)*gxa(im,6))+&
& gmet(3,1)*(gxa(re,3)*gxa(re,1)+gxa(im,3)*gxa(im,1))+&
& gmet(2,1)*(gxa(re,4)*gxa(re,1)+gxa(im,4)*gxa(im,1))+&
& gmet(2,3)*(gxa(re,4)*gxa(re,5)+gxa(im,4)*gxa(im,5))+&
& gmet(1,3)*(gxa(re,5)*gxa(re,5)+gxa(im,5)*gxa(im,5))+&
& gmet(1,2)*(gxa(re,5)*gxa(re,6)+gxa(im,5)*gxa(im,6)))

!a=2, b=1 in rank2(a,b) --> maps to index 6
 rank2(6)=2.0d0*(&
& gmet(1,1)*(gxa(re,6)*gxa(re,1)+gxa(im,6)*gxa(im,1))+&
& gmet(2,2)*(gxa(re,2)*gxa(re,6)+gxa(im,2)*gxa(im,6))+&
& gmet(3,3)*(gxa(re,4)*gxa(re,5)+gxa(im,4)*gxa(im,5))+&
& gmet(3,2)*(gxa(re,4)*gxa(re,6)+gxa(im,4)*gxa(im,6))+&
& gmet(3,1)*(gxa(re,4)*gxa(re,1)+gxa(im,4)*gxa(im,1))+&
& gmet(2,1)*(gxa(re,2)*gxa(re,1)+gxa(im,2)*gxa(im,1))+&
& gmet(2,3)*(gxa(re,2)*gxa(re,5)+gxa(im,2)*gxa(im,5))+&
& gmet(1,3)*(gxa(re,6)*gxa(re,5)+gxa(im,6)*gxa(im,5))+&
& gmet(1,2)*(gxa(re,6)*gxa(re,6)+gxa(im,6)*gxa(im,6)))

end subroutine cont22
!!***

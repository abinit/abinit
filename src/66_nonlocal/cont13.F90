!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont13
!! NAME
!! cont13
!!
!! FUNCTION
!! Contract rank1 tensor with rank3 symmetric tensor to
!! produce symmetric rank2 tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  rank1(2,3)=rank 1 complex tensor (vector of length 3)
!!  rank3(2,10)=rank 3 complex tensor (symmetric storage)
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (symmetric storage)
!!
!! NOTES
!! Tensors are in "symmetric" storage mode.
!! For rank1 this is 1, 2, 3;
!! for rank2 this is 11, 22, 33, 32, 31, 21;
!! for rank3 this is 111, 221, 331, 321, 311, 211, 222, 332, 322, 333.
!! rank1 and rank3 are complex; rank2 is real.
!! Want $2 Re[contraction]$.
!! $rank2(a,b)=2 Re[rank1(i)^"*" rank3(a,b,i)]$.
!! In typical usage the input rank1 tensor is actually
!! $rank1(i)=gmet(i,j) gxa(j)$
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


subroutine cont13(rank1,rank3,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont13'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: rank1(2,3),rank3(2,10)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1

! *************************************************************************

!Simply write out index summations
!a=1, b=1 in rank2(a,b) --> maps to index 1
 rank2(1)=2.0d0*(&
& (rank1(re,1)*rank3(re,1)+rank1(im,1)*rank3(im,1))+&
& (rank1(re,2)*rank3(re,6)+rank1(im,2)*rank3(im,6))+&
& (rank1(re,3)*rank3(re,5)+rank1(im,3)*rank3(im,5)))

!a=2, b=2 in rank2(a,b) --> maps to index 2
 rank2(2)=2.0d0*(&
& (rank1(re,1)*rank3(re,2)+rank1(im,1)*rank3(im,2))+&
& (rank1(re,2)*rank3(re,7)+rank1(im,2)*rank3(im,7))+&
& (rank1(re,3)*rank3(re,9)+rank1(im,3)*rank3(im,9)))

!a=3, b=3 in rank2(a,b) --> maps to index 3
 rank2(3)=2.0d0*(&
& (rank1(re,1)*rank3(re,3)+rank1(im,1)*rank3(im,3))+&
& (rank1(re,2)*rank3(re,8)+rank1(im,2)*rank3(im,8))+&
& (rank1(re,3)*rank3(re,10)+rank1(im,3)*rank3(im,10)))

!a=3, b=2 in rank2(a,b) --> maps to index 4
 rank2(4)=2.0d0*(&
& (rank1(re,1)*rank3(re,4)+rank1(im,1)*rank3(im,4))+&
& (rank1(re,2)*rank3(re,9)+rank1(im,2)*rank3(im,9))+&
& (rank1(re,3)*rank3(re,8)+rank1(im,3)*rank3(im,8)))

!a=3, b=1 in rank2(a,b) --> maps to index 5
 rank2(5)=2.0d0*(&
& (rank1(re,1)*rank3(re,5)+rank1(im,1)*rank3(im,5))+&
& (rank1(re,2)*rank3(re,4)+rank1(im,2)*rank3(im,4))+&
& (rank1(re,3)*rank3(re,3)+rank1(im,3)*rank3(im,3)))

!a=2, b=1 in rank2(a,b) --> maps to index 6
 rank2(6)=2.0d0*(&
& (rank1(re,1)*rank3(re,6)+rank1(im,1)*rank3(im,6))+&
& (rank1(re,2)*rank3(re,2)+rank1(im,2)*rank3(im,2))+&
& (rank1(re,3)*rank3(re,4)+rank1(im,3)*rank3(im,4)))

end subroutine cont13
!!***

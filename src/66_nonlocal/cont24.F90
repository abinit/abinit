!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont24
!! NAME
!! cont24
!!
!! FUNCTION
!! Contract symmetric rank2 tensor gxa with rank4 symmetric tensor to
!! produce symmetric rank2 tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa(2,6)=rank 2 symmetric complex tensor in order 11 22 33 32 31 21
!!  rank4(2,15)=rank 4 complex tensor (symmetric storage)
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (symmetric storage) 11 22 33 32 31 21.
!!
!! NOTES
!! Tensors are in "symmetric" storage mode.
!! for gxa and rank2 this is 11, 22, 33, 32, 31, 21;
!! for the rank 4 tensor rank4 this is
!! 1111 2211 3311 3211 3111 2111 2221 3321 3221 3331 2222 3322 3222 3332 3333.
!! gxa and rank4 are complex; rank2 is real.
!! Want $2 Re[contraction]$.
!! $rank2(a,b)=2 Re[gxa(i,j)^"*" rank4(a,b,i,j)]$.
!!
!! Note that the input gxa is typically the result of
!! $gxa(i,j)=[{3 \over 2} gmet(i,l) gmet(j,m) - {1 \over 2} gmet(i,j) gmet(l,m)] gxa_old(l,m)$
!! where the subroutine "metcon" already includes weights in the
!! definition of gxa for off-diagonal elements (weight of 2 for
!! symmetry).
!! Components 4, 5, and 6 of gxa have already been multiplied by 2
!! so the expressions below do not carry the 2.
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


subroutine cont24(gxa,rank4,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont24'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gxa(2,6),rank4(2,15)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1

! *************************************************************************

!Simply write out index summations

!a=1, b=1 in rank2(a,b) --> maps to index 1
 rank2(1)=2.0d0*(&
& gxa(re,1)*rank4(re, 1)+gxa(im,1)*rank4(im, 1)+&
& gxa(re,2)*rank4(re, 2)+gxa(im,2)*rank4(im, 2)+&
& gxa(re,3)*rank4(re, 3)+gxa(im,3)*rank4(im, 3)+&
& gxa(re,4)*rank4(re, 4)+gxa(im,4)*rank4(im, 4)+&
& gxa(re,5)*rank4(re, 5)+gxa(im,5)*rank4(im, 5)+&
& gxa(re,6)*rank4(re, 6)+gxa(im,6)*rank4(im, 6))

!a=2, b=2 in rank2(a,b) --> maps to index 2
 rank2(2)=2.0d0*(&
& gxa(re,1)*rank4(re, 2)+gxa(im,1)*rank4(im, 2)+&
& gxa(re,2)*rank4(re,11)+gxa(im,2)*rank4(im,11)+&
& gxa(re,3)*rank4(re,12)+gxa(im,3)*rank4(im,12)+&
& gxa(re,4)*rank4(re,13)+gxa(im,4)*rank4(im,13)+&
& gxa(re,5)*rank4(re, 9)+gxa(im,5)*rank4(im, 9)+&
& gxa(re,6)*rank4(re, 7)+gxa(im,6)*rank4(im, 7))

!a=3, b=3 in rank2(a,b) --> maps to index 3
 rank2(3)=2.0d0*(&
& gxa(re,1)*rank4(re, 3)+gxa(im,1)*rank4(im, 3)+&
& gxa(re,2)*rank4(re,12)+gxa(im,2)*rank4(im,12)+&
& gxa(re,3)*rank4(re,15)+gxa(im,3)*rank4(im,15)+&
& gxa(re,4)*rank4(re,14)+gxa(im,4)*rank4(im,14)+&
& gxa(re,5)*rank4(re,10)+gxa(im,5)*rank4(im,10)+&
& gxa(re,6)*rank4(re, 8)+gxa(im,6)*rank4(im, 8))

!a=3, b=2 in rank2(a,b) --> maps to index 4
 rank2(4)=2.0d0*(&
& gxa(re,1)*rank4(re, 4)+gxa(im,1)*rank4(im, 4)+&
& gxa(re,2)*rank4(re,13)+gxa(im,2)*rank4(im,13)+&
& gxa(re,3)*rank4(re,14)+gxa(im,3)*rank4(im,14)+&
& gxa(re,4)*rank4(re,12)+gxa(im,4)*rank4(im,12)+&
& gxa(re,5)*rank4(re, 8)+gxa(im,5)*rank4(im, 8)+&
& gxa(re,6)*rank4(re, 9)+gxa(im,6)*rank4(im, 9))

!a=3, b=1 in rank2(a,b) --> maps to index 5
 rank2(5)=2.0d0*(&
& gxa(re,1)*rank4(re, 5)+gxa(im,1)*rank4(im, 5)+&
& gxa(re,2)*rank4(re, 9)+gxa(im,2)*rank4(im, 9)+&
& gxa(re,3)*rank4(re,10)+gxa(im,3)*rank4(im,10)+&
& gxa(re,4)*rank4(re, 8)+gxa(im,4)*rank4(im, 8)+&
& gxa(re,5)*rank4(re, 3)+gxa(im,5)*rank4(im, 3)+&
& gxa(re,6)*rank4(re, 4)+gxa(im,6)*rank4(im, 4))

!a=2, b=1 in rank2(a,b) --> maps to index 6
 rank2(6)=2.0d0*(&
& gxa(re,1)*rank4(re, 6)+gxa(im,1)*rank4(im, 6)+&
& gxa(re,2)*rank4(re, 7)+gxa(im,2)*rank4(im, 7)+&
& gxa(re,3)*rank4(re, 8)+gxa(im,3)*rank4(im, 8)+&
& gxa(re,4)*rank4(re, 9)+gxa(im,4)*rank4(im, 9)+&
& gxa(re,5)*rank4(re, 4)+gxa(im,5)*rank4(im, 4)+&
& gxa(re,6)*rank4(re, 2)+gxa(im,6)*rank4(im, 2))

end subroutine cont24
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/cont35
!! NAME
!! cont35
!!
!! FUNCTION
!! Contract symmetric rank3 tensor gxa with rank5 symmetric tensor to
!! produce symmetric rank2 tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  gxa(2,10)=rank 3 symmetric complex tensor in order
!!  rank5(2,21)=rank 5 complex tensor (symmetric storage)
!!
!! OUTPUT
!!  rank2(6)=rank 2 real tensor (symmetric storage) 11 22 33 32 31 21.
!!
!! NOTES
!! Tensors are in "symmetric" storage mode.
!! For rank 3 tensor gxa this is
!!     111 221 331 321 311 211 222 332 322 333;
!! For rank 5 tensor rank5 this is
!!     11111 22111 33111 32111 31111 21111 22211 33211 32211 33311
!!     22221 33221 32221 33321 33331 22222 33222 32222 33322 33332 33333;
!! For rank 2 tensor rank2 this is 11, 22, 33, 32, 31, 21;
!! gxa and rank5 are complex; rank2 is real.
!! Want $2 Re[contraction]$.
!! $rank2(a,b)=2 Re[gxa(i,j,k)^"*" rank5(a,b,i,j,k)]$.
!!
!! Note that the input gxa is typically the result of
!!{{\ \begin{equation}
!! gxa(i,j,k)=[{5 \over 2} gmet(i,l) gmet(j,m) gmet(k,n) - {3 \over 2} gmet(i,j) gmet(l,m) gmet(k,n)] gxa_old(l,m)
!!\end{equation} }}
!! where the subroutine "metcon" already includes weights in the definition
!! of gxa for off-diagonal elements.
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


subroutine cont35(gxa,rank5,rank2)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cont35'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gxa(2,10),rank5(2,21)
 real(dp),intent(out) :: rank2(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1

! *************************************************************************

!Simply write out index summations

!a=1, b=1 in rank2(a,b) --> maps to index 1
 rank2(1)=2.0d0*(&
& gxa(re, 1)*rank5(re, 1)+gxa(im, 1)*rank5(im, 1)+&
& gxa(re, 7)*rank5(re, 7)+gxa(im, 7)*rank5(im, 7)+&
& gxa(re,10)*rank5(re,10)+gxa(im,10)*rank5(im,10)+&
& gxa(re, 2)*rank5(re, 2)+gxa(im, 2)*rank5(im, 2)+&
& gxa(re, 3)*rank5(re, 3)+gxa(im, 3)*rank5(im, 3)+&
& gxa(re, 5)*rank5(re, 5)+gxa(im, 5)*rank5(im, 5)+&
& gxa(re, 6)*rank5(re, 6)+gxa(im, 6)*rank5(im, 6)+&
& gxa(re, 8)*rank5(re, 8)+gxa(im, 8)*rank5(im, 8)+&
& gxa(re, 9)*rank5(re, 9)+gxa(im, 9)*rank5(im, 9)+&
& gxa(re, 4)*rank5(re, 4)+gxa(im, 4)*rank5(im, 4))


!a=2, b=2 in rank2(a,b) --> maps to index 2
 rank2(2)=2.0d0*(&
& gxa(re, 1)*rank5(re, 2)+gxa(im, 1)*rank5(im, 2)+&
& gxa(re, 7)*rank5(re,16)+gxa(im, 7)*rank5(im,16)+&
& gxa(re,10)*rank5(re,19)+gxa(im,10)*rank5(im,19)+&
& gxa(re, 2)*rank5(re,11)+gxa(im, 2)*rank5(im,11)+&
& gxa(re, 3)*rank5(re,12)+gxa(im, 3)*rank5(im,12)+&
& gxa(re, 5)*rank5(re, 9)+gxa(im, 5)*rank5(im, 9)+&
& gxa(re, 6)*rank5(re, 7)+gxa(im, 6)*rank5(im, 7)+&
& gxa(re, 8)*rank5(re,17)+gxa(im, 8)*rank5(im,17)+&
& gxa(re, 9)*rank5(re,18)+gxa(im, 9)*rank5(im,18)+&
& gxa(re, 4)*rank5(re,13)+gxa(im, 4)*rank5(im,13))

!a=3, b=3 in rank2(a,b) --> maps to index 3
 rank2(3)=2.0d0*(&
& gxa(re, 1)*rank5(re, 3)+gxa(im, 1)*rank5(im, 3)+&
& gxa(re, 7)*rank5(re,17)+gxa(im, 7)*rank5(im,17)+&
& gxa(re,10)*rank5(re,21)+gxa(im,10)*rank5(im,21)+&
& gxa(re, 2)*rank5(re,12)+gxa(im, 2)*rank5(im,12)+&
& gxa(re, 3)*rank5(re,15)+gxa(im, 3)*rank5(im,15)+&
& gxa(re, 5)*rank5(re,10)+gxa(im, 5)*rank5(im,10)+&
& gxa(re, 6)*rank5(re, 8)+gxa(im, 6)*rank5(im, 8)+&
& gxa(re, 8)*rank5(re,20)+gxa(im, 8)*rank5(im,20)+&
& gxa(re, 9)*rank5(re,19)+gxa(im, 9)*rank5(im,19)+&
& gxa(re, 4)*rank5(re,14)+gxa(im, 4)*rank5(im,14))

!a=3, b=2 in rank2(a,b) --> maps to index 4
 rank2(4)=2.0d0*(&
& gxa(re, 1)*rank5(re, 4)+gxa(im, 1)*rank5(im, 4)+&
& gxa(re, 7)*rank5(re,18)+gxa(im, 7)*rank5(im,18)+&
& gxa(re,10)*rank5(re,20)+gxa(im,10)*rank5(im,20)+&
& gxa(re, 2)*rank5(re,13)+gxa(im, 2)*rank5(im,13)+&
& gxa(re, 3)*rank5(re,14)+gxa(im, 3)*rank5(im,14)+&
& gxa(re, 5)*rank5(re, 8)+gxa(im, 5)*rank5(im, 8)+&
& gxa(re, 6)*rank5(re, 9)+gxa(im, 6)*rank5(im, 9)+&
& gxa(re, 8)*rank5(re,19)+gxa(im, 8)*rank5(im,19)+&
& gxa(re, 9)*rank5(re,17)+gxa(im, 9)*rank5(im,17)+&
& gxa(re, 4)*rank5(re,12)+gxa(im, 4)*rank5(im,12))

!a=3, b=1 in rank2(a,b) --> maps to index 5
 rank2(5)=2.0d0*(&
& gxa(re, 1)*rank5(re, 5)+gxa(im, 1)*rank5(im, 5)+&
& gxa(re, 7)*rank5(re,13)+gxa(im, 7)*rank5(im,13)+&
& gxa(re,10)*rank5(re,15)+gxa(im,10)*rank5(im,15)+&
& gxa(re, 2)*rank5(re, 9)+gxa(im, 2)*rank5(im, 9)+&
& gxa(re, 3)*rank5(re,10)+gxa(im, 3)*rank5(im,10)+&
& gxa(re, 5)*rank5(re, 3)+gxa(im, 5)*rank5(im, 3)+&
& gxa(re, 6)*rank5(re, 4)+gxa(im, 6)*rank5(im, 4)+&
& gxa(re, 8)*rank5(re,14)+gxa(im, 8)*rank5(im,14)+&
& gxa(re, 9)*rank5(re,12)+gxa(im, 9)*rank5(im,12)+&
& gxa(re, 4)*rank5(re, 8)+gxa(im, 4)*rank5(im, 8))

!a=2, b=1 in rank2(a,b) --> maps to index 6
 rank2(6)=2.0d0*(&
& gxa(re, 1)*rank5(re, 6)+gxa(im, 1)*rank5(im, 6)+&
& gxa(re, 7)*rank5(re,11)+gxa(im, 7)*rank5(im,11)+&
& gxa(re,10)*rank5(re,14)+gxa(im,10)*rank5(im,14)+&
& gxa(re, 2)*rank5(re, 7)+gxa(im, 2)*rank5(im, 7)+&
& gxa(re, 3)*rank5(re, 8)+gxa(im, 3)*rank5(im, 8)+&
& gxa(re, 5)*rank5(re, 4)+gxa(im, 5)*rank5(im, 4)+&
& gxa(re, 6)*rank5(re, 2)+gxa(im, 6)*rank5(im, 2)+&
& gxa(re, 8)*rank5(re,12)+gxa(im, 8)*rank5(im,12)+&
& gxa(re, 9)*rank5(re,13)+gxa(im, 9)*rank5(im,13)+&
& gxa(re, 4)*rank5(re, 9)+gxa(im, 4)*rank5(im, 9))

end subroutine cont35
!!***

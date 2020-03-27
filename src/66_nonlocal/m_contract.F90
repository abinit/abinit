!!****m* ABINIT/m_contract
!! NAME
!!  m_contract
!!
!! FUNCTION
!! Low-level procedeures used in nonlop_pl to contract tensors
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, MT, GZ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_contract

 use defs_basis
 use m_errors

 implicit none

 private
!!***

 public :: cont22cso
 public :: cont22so
 public :: cont24
 public :: cont33cso
 public :: cont33so
 public :: cont35
 public :: cont22
 public :: cont3
 public :: cont13
 public :: metcon
 public :: metcon_so
 public :: metric_so
!!***

contains
!!***

!!****f* m_contract/cont13
!! NAME
!! cont13
!!
!! FUNCTION
!! Contract rank1 tensor with rank3 symmetric tensor to
!! produce symmetric rank2 tensor
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

subroutine cont13(rank1,rank3,rank2)

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


!!****f* m_contract/cont22
!! NAME
!! cont22
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa with itself using gmet to
!! produce symmetric rank 2 tensor.
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

subroutine cont22(gxa,gmet,rank2)

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

!!****f* m_contract/cont22cso
!! NAME
!! cont22cso
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa1 with symmetric rank 2 tensor
!! gxa2 using metric tensor gmet to produce rank 2 complex tensor.
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

subroutine cont22cso(gxa1,gxa2,gmet,rank2c)

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


!!****f* m_contract/cont22so
!! NAME
!! cont22so
!!
!! FUNCTION
!! Contract symmetric rank 2 tensor gxa1 with symmetric rank 2 tensor
!! gxa2 using antisymmetric tensor amet to produce rank 2 real tensor.
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

subroutine cont22so(gxa1,gxa2,amet,rank2)

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


!!****f* m_contract/cont24
!! NAME
!! cont24
!!
!! FUNCTION
!! Contract symmetric rank2 tensor gxa with rank4 symmetric tensor to
!! produce symmetric rank2 tensor.
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

subroutine cont24(gxa,rank4,rank2)

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

!!****f* m_contract/cont3
!! NAME
!! cont3
!!
!! FUNCTION
!! Compute several specialized contractions needed for the
!! l=3 part of the stress tensor.
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

subroutine cont3(gxa,gmet,rank2)

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

!!****f* m_contract/cont33cso
!! NAME
!! cont33cso
!!
!! FUNCTION
!! Contract symmetric rank 3 tensor gxa1 with symmetric rank 3 tensor
!! gxa2 using metric tensor gmet to produce rank 2 complex tensor.
!!
!! INPUTS
!!  gxa1(2,10)=rank 3 complex symmetric tensor
!!  gxa2(2,10)=rank 3 complex symmetric tensor
!!  gmet(3,3)=usual metric tensor (symmetric, real)
!!
!! OUTPUT
!!  rank2c(2,6)=rank 2 complex tensor (pseudo-symmetric storage)
!!
!! NOTES
!! This contraction is used for spin-orbit correction in non-local
!! contribution to stresses.
!!
!! Symmetric gxa1, gxa2 are stored as
!!               111 221 331 321 311 211 222 332 322 333;
!! gmet(3,3) is symmetric but stored fully (9 elements);
!! Output rank2c is not symmetric but since
!!      $rank2c_{gxa1,gxa2}(a,b)=conjg(rank2c_{gxa2,gxa1}(b,a))$
!!       it is stored as 11 22 33 32 31 21.
!!
!! rank2c(1,1), rank2c(2,2), rank3c(3,3) are not needed;
!! They are not calculated.
!!
!!{{\ \begin{equation}
!! rank2c(a,b)=7.5 conjg(gxa1(a,i,j))*r_3(i,j,b) - 1.5 r_{11}(a)*r_{12}(b)
!!\end{equation} }}
!!   where:
!!{{\ \begin{eqnarray}
!! r_3(i,j,b) & = & gxa2(b,l,m) gmet(i,l) gmet(j,m) \nonumber
!! r_{11}(a)  & = & conjg(gxa1(a,l,m)) gmet(l,m)    \nonumber
!! r_{12}(b)  & = & gxa2(b,l,m) gmet(l,m)
!! \end{eqnarray} }}
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine cont33cso(gxa1,gxa2,gmet,rank2c)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gmet(3,3),gxa1(2,10),gxa2(2,10)
 real(dp),intent(out) :: rank2c(2,6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
!arrays
 real(dp) :: r11(2,3),r12(2,3),r3(2,6,3),r3a(2,6,3)

! *************************************************************************

!Initialize output tensor
 rank2c(:,:)=0.d0

!Compute r3(i,j,b)=gxa2(b,l,m)*gmet(i,l)*gmet(j,m)
!stored as 11 22 33 32 31 21 for (i,j)
!First compute r3a(b,l,j)=gxa2(b,l,m)*gmet(j,m)
!stored as r3a(re/im,(b,l),j)
 r3a(:,1,1)=gxa2(:,1 )*gmet(1,1)+gxa2(:,6 )*gmet(1,2)+gxa2(:,5 )*gmet(1,3)
 r3a(:,2,1)=gxa2(:,2 )*gmet(1,1)+gxa2(:,7 )*gmet(1,2)+gxa2(:,9 )*gmet(1,3)
 r3a(:,3,1)=gxa2(:,3 )*gmet(1,1)+gxa2(:,8 )*gmet(1,2)+gxa2(:,10)*gmet(1,3)
 r3a(:,4,1)=gxa2(:,4 )*gmet(1,1)+gxa2(:,9 )*gmet(1,2)+gxa2(:,8 )*gmet(1,3)
 r3a(:,5,1)=gxa2(:,5 )*gmet(1,1)+gxa2(:,4 )*gmet(1,2)+gxa2(:,3 )*gmet(1,3)
 r3a(:,6,1)=gxa2(:,6 )*gmet(1,1)+gxa2(:,2 )*gmet(1,2)+gxa2(:,4 )*gmet(1,3)
 r3a(:,1,2)=gxa2(:,1 )*gmet(2,1)+gxa2(:,6 )*gmet(2,2)+gxa2(:,5 )*gmet(2,3)
 r3a(:,2,2)=gxa2(:,2 )*gmet(2,1)+gxa2(:,7 )*gmet(2,2)+gxa2(:,9 )*gmet(2,3)
 r3a(:,3,2)=gxa2(:,3 )*gmet(2,1)+gxa2(:,8 )*gmet(2,2)+gxa2(:,10)*gmet(2,3)
 r3a(:,4,2)=gxa2(:,4 )*gmet(2,1)+gxa2(:,9 )*gmet(2,2)+gxa2(:,8 )*gmet(2,3)
 r3a(:,5,2)=gxa2(:,5 )*gmet(2,1)+gxa2(:,4 )*gmet(2,2)+gxa2(:,3 )*gmet(2,3)
 r3a(:,6,2)=gxa2(:,6 )*gmet(2,1)+gxa2(:,2 )*gmet(2,2)+gxa2(:,4 )*gmet(2,3)
 r3a(:,1,3)=gxa2(:,1 )*gmet(3,1)+gxa2(:,6 )*gmet(3,2)+gxa2(:,5 )*gmet(3,3)
 r3a(:,2,3)=gxa2(:,2 )*gmet(3,1)+gxa2(:,7 )*gmet(3,2)+gxa2(:,9 )*gmet(3,3)
 r3a(:,3,3)=gxa2(:,3 )*gmet(3,1)+gxa2(:,8 )*gmet(3,2)+gxa2(:,10)*gmet(3,3)
 r3a(:,4,3)=gxa2(:,4 )*gmet(3,1)+gxa2(:,9 )*gmet(3,2)+gxa2(:,8 )*gmet(3,3)
 r3a(:,5,3)=gxa2(:,5 )*gmet(3,1)+gxa2(:,4 )*gmet(3,2)+gxa2(:,3 )*gmet(3,3)
 r3a(:,6,3)=gxa2(:,6 )*gmet(3,1)+gxa2(:,2 )*gmet(3,2)+gxa2(:,4 )*gmet(3,3)

!Then  compute r3(i,j,b)=r3a(b,l,j)*gmet(i,l)
!stored as r3(re/im,(i,j),b)
 r3(:,1,1)=r3a(:,1,1)*gmet(1,1)+r3a(:,6,1)*gmet(1,2)+r3a(:,5,1)*gmet(1,3)
 r3(:,1,2)=r3a(:,6,1)*gmet(1,1)+r3a(:,2,1)*gmet(1,2)+r3a(:,4,1)*gmet(1,3)
 r3(:,1,3)=r3a(:,5,1)*gmet(1,1)+r3a(:,4,1)*gmet(1,2)+r3a(:,3,1)*gmet(1,3)
 r3(:,2,1)=r3a(:,1,2)*gmet(2,1)+r3a(:,6,2)*gmet(2,2)+r3a(:,5,2)*gmet(2,3)
 r3(:,2,2)=r3a(:,6,2)*gmet(2,1)+r3a(:,2,2)*gmet(2,2)+r3a(:,4,2)*gmet(2,3)
 r3(:,2,3)=r3a(:,5,2)*gmet(2,1)+r3a(:,4,2)*gmet(2,2)+r3a(:,3,2)*gmet(2,3)
 r3(:,3,1)=r3a(:,1,3)*gmet(3,1)+r3a(:,6,3)*gmet(3,2)+r3a(:,5,3)*gmet(3,3)
 r3(:,3,2)=r3a(:,6,3)*gmet(3,1)+r3a(:,2,3)*gmet(3,2)+r3a(:,4,3)*gmet(3,3)
 r3(:,3,3)=r3a(:,5,3)*gmet(3,1)+r3a(:,4,3)*gmet(3,2)+r3a(:,3,3)*gmet(3,3)
 r3(:,4,1)=r3a(:,1,2)*gmet(3,1)+r3a(:,6,2)*gmet(3,2)+r3a(:,5,2)*gmet(3,3)
 r3(:,4,2)=r3a(:,6,2)*gmet(3,1)+r3a(:,2,2)*gmet(3,2)+r3a(:,4,2)*gmet(3,3)
 r3(:,4,3)=r3a(:,5,2)*gmet(3,1)+r3a(:,4,2)*gmet(3,2)+r3a(:,3,2)*gmet(3,3)
 r3(:,5,1)=r3a(:,1,1)*gmet(3,1)+r3a(:,6,1)*gmet(3,2)+r3a(:,5,1)*gmet(3,3)
 r3(:,5,2)=r3a(:,6,1)*gmet(3,1)+r3a(:,2,1)*gmet(3,2)+r3a(:,4,1)*gmet(3,3)
 r3(:,5,3)=r3a(:,5,1)*gmet(3,1)+r3a(:,4,1)*gmet(3,2)+r3a(:,3,1)*gmet(3,3)
 r3(:,6,1)=r3a(:,1,1)*gmet(2,1)+r3a(:,6,1)*gmet(2,2)+r3a(:,5,1)*gmet(2,3)
 r3(:,6,2)=r3a(:,6,1)*gmet(2,1)+r3a(:,2,1)*gmet(2,2)+r3a(:,4,1)*gmet(2,3)
 r3(:,6,3)=r3a(:,5,1)*gmet(2,1)+r3a(:,4,1)*gmet(2,2)+r3a(:,3,1)*gmet(2,3)

!Compute r11(a)=conjg(gxa1(a,l,m))*gmet(l,m)
 r11(:,1)=gxa1(:,1)*gmet(1,1)+gxa1(:,2)*gmet(2,2)+gxa1(:,3 )*gmet(3,3)&
& +2.d0*(gxa1(:,4)*gmet(3,2)+gxa1(:,5)*gmet(3,1)+gxa1(:,6 )*gmet(2,1))
 r11(:,2)=gxa1(:,6)*gmet(1,1)+gxa1(:,7)*gmet(2,2)+gxa1(:,8 )*gmet(3,3)&
& +2.d0*(gxa1(:,9)*gmet(3,2)+gxa1(:,4)*gmet(3,1)+gxa1(:,2 )*gmet(2,1))
 r11(:,3)=gxa1(:,5)*gmet(1,1)+gxa1(:,9)*gmet(2,2)+gxa1(:,10)*gmet(3,3)&
& +2.d0*(gxa1(:,8)*gmet(3,2)+gxa1(:,3)*gmet(3,1)+gxa1(:,4 )*gmet(2,1))
 r11(im,1)=-r11(im,1);r11(im,2)=-r11(im,2);r11(im,3)=-r11(im,3)

!Compute r12(b)=gxa2(b,l,m)*gmet(l,m)
 r12(:,1)=gxa2(:,1)*gmet(1,1)+gxa2(:,2)*gmet(2,2)+gxa2(:,3 )*gmet(3,3)&
& +2.d0*(gxa2(:,4)*gmet(3,2)+gxa2(:,5)*gmet(3,1)+gxa2(:,6 )*gmet(2,1))
 r12(:,2)=gxa2(:,6)*gmet(1,1)+gxa2(:,7)*gmet(2,2)+gxa2(:,8 )*gmet(3,3)&
& +2.d0*(gxa2(:,9)*gmet(3,2)+gxa2(:,4)*gmet(3,1)+gxa2(:,2 )*gmet(2,1))
 r12(:,3)=gxa2(:,5)*gmet(1,1)+gxa2(:,9)*gmet(2,2)+gxa2(:,10)*gmet(3,3)&
& +2.d0*(gxa2(:,8)*gmet(3,2)+gxa2(:,3)*gmet(3,1)+gxa2(:,4 )*gmet(2,1))

!Finally compute rank2c(a,b)=7.5*conjg(gxa1(a,i,j))*r3(i,j,b) - 1.5*r11(a)*r12(b)
!rank2c(re,1)=7.5d0*(gxa1(re,1 )*r3(re,1,1)+gxa1(im,1 )*r3(im,1,1)&
!&                    +gxa1(re,2 )*r3(re,2,1)+gxa1(im,2 )*r3(im,2,1)&
!&                    +gxa1(re,3 )*r3(re,3,1)+gxa1(im,3 )*r3(im,3,1))&
!&             +15.d0*(gxa1(re,4 )*r3(re,4,1)+gxa1(im,4 )*r3(im,4,1)&
!&                    +gxa1(re,5 )*r3(re,5,1)+gxa1(im,5 )*r3(im,5,1)&
!&                    +gxa1(re,6 )*r3(re,6,1)+gxa1(im,6 )*r3(im,6,1))&
!&             -1.5d0*(r11(re,1)*r12(re,1)-r11(im,1)*r12(im,1))
!rank2c(re,2)=7.5d0*(gxa1(re,6 )*r3(re,1,2)+gxa1(im,6 )*r3(im,1,2)&
!&                    +gxa1(re,7 )*r3(re,2,2)+gxa1(im,7 )*r3(im,2,2)&
!&                    +gxa1(re,8 )*r3(re,3,2)+gxa1(im,8 )*r3(im,3,2))&
!&             +15.d0*(gxa1(re,9 )*r3(re,4,2)+gxa1(im,9 )*r3(im,4,2)&
!&                    +gxa1(re,4 )*r3(re,5,2)+gxa1(im,4 )*r3(im,5,2)&
!&                    +gxa1(re,2 )*r3(re,6,2)+gxa1(im,2 )*r3(im,6,2))&
!&             -1.5d0*(r11(re,2)*r12(re,2)-r11(im,2)*r12(im,2))
!rank2c(re,3)=7.5d0*(gxa1(re,5 )*r3(re,1,3)+gxa1(im,5 )*r3(im,1,3)&
!&                    +gxa1(re,9 )*r3(re,2,3)+gxa1(im,9 )*r3(im,2,3)&
!&                    +gxa1(re,10)*r3(re,3,3)+gxa1(im,10)*r3(im,3,3))&
!&             +15.d0*(gxa1(re,8 )*r3(re,4,3)+gxa1(im,8 )*r3(im,4,3)&
!&                    +gxa1(re,3 )*r3(re,5,3)+gxa1(im,3 )*r3(im,5,3)&
!&                    +gxa1(re,4 )*r3(re,6,3)+gxa1(im,4 )*r3(im,6,3))&
!&             -1.5d0*(r11(re,3)*r12(re,3)-r11(im,3)*r12(im,3))
 rank2c(re,4)=7.5d0*(gxa1(re,5 )*r3(re,1,2)+gxa1(im,5 )*r3(im,1,2)&
& +gxa1(re,9 )*r3(re,2,2)+gxa1(im,9 )*r3(im,2,2)&
& +gxa1(re,10)*r3(re,3,2)+gxa1(im,10)*r3(im,3,2))&
& +15.d0*(gxa1(re,8 )*r3(re,4,2)+gxa1(im,8 )*r3(im,4,2)&
& +gxa1(re,3 )*r3(re,5,2)+gxa1(im,3 )*r3(im,5,2)&
& +gxa1(re,4 )*r3(re,6,2)+gxa1(im,4 )*r3(im,6,2))&
& -1.5d0*(r11(re,3)*r12(re,2)-r11(im,3)*r12(im,2))
 rank2c(re,5)=7.5d0*(gxa1(re,5 )*r3(re,1,1)+gxa1(im,5 )*r3(im,1,1)&
& +gxa1(re,9 )*r3(re,2,1)+gxa1(im,9 )*r3(im,2,1)&
& +gxa1(re,10)*r3(re,3,1)+gxa1(im,10)*r3(im,3,1))&
& +15.d0*(gxa1(re,8 )*r3(re,4,1)+gxa1(im,8 )*r3(im,4,1)&
& +gxa1(re,3 )*r3(re,5,1)+gxa1(im,3 )*r3(im,5,1)&
& +gxa1(re,4 )*r3(re,6,1)+gxa1(im,4 )*r3(im,6,1))&
& -1.5d0*(r11(re,3)*r12(re,1)-r11(im,3)*r12(im,1))
 rank2c(re,6)=7.5d0*(gxa1(re,6 )*r3(re,1,1)+gxa1(im,6 )*r3(im,1,1)&
& +gxa1(re,7 )*r3(re,2,1)+gxa1(im,7 )*r3(im,2,1)&
& +gxa1(re,8 )*r3(re,3,1)+gxa1(im,8 )*r3(im,3,1))&
& +15.d0*(gxa1(re,9 )*r3(re,4,1)+gxa1(im,9 )*r3(im,4,1)&
& +gxa1(re,4 )*r3(re,5,1)+gxa1(im,4 )*r3(im,5,1)&
& +gxa1(re,2 )*r3(re,6,1)+gxa1(im,2 )*r3(im,6,1))&
& -1.5d0*(r11(re,2)*r12(re,1)-r11(im,2)*r12(im,1))
!rank2c(im,1)=7.5d0*(gxa1(re,1 )*r3(im,1,1)-gxa1(im,1 )*r3(re,1,1)&
!&                    +gxa1(re,2 )*r3(im,2,1)-gxa1(im,2 )*r3(re,2,1)&
!&                    +gxa1(re,3 )*r3(im,3,1)-gxa1(im,3 )*r3(re,3,1))&
!&             +15.d0*(gxa1(re,4 )*r3(im,4,1)-gxa1(im,4 )*r3(re,4,1)&
!&                    +gxa1(re,5 )*r3(im,5,1)-gxa1(im,5 )*r3(re,5,1)&
!&                    +gxa1(re,6 )*r3(im,6,1)-gxa1(im,6 )*r3(re,6,1))&
!&             -1.5d0*(r11(re,1)*r12(im,1)+r11(im,1)*r12(re,1))
!rank2c(im,2)=7.5d0*(gxa1(re,6 )*r3(im,1,2)-gxa1(im,6 )*r3(re,1,2)&
!&                    +gxa1(re,7 )*r3(im,2,2)-gxa1(im,7 )*r3(re,2,2)&
!&                    +gxa1(re,8 )*r3(im,3,2)-gxa1(im,8 )*r3(re,3,2))&
!&             +15.d0*(gxa1(re,9 )*r3(im,4,2)-gxa1(im,9 )*r3(re,4,2)&
!&                    +gxa1(re,4 )*r3(im,5,2)-gxa1(im,4 )*r3(re,5,2)&
!&                    +gxa1(re,2 )*r3(im,6,2)-gxa1(im,2 )*r3(re,6,2))&
!&             -1.5d0*(r11(re,2)*r12(im,2)+r11(im,2)*r12(re,2))
!rank2c(im,3)=7.5d0*(gxa1(re,5 )*r3(im,1,3)-gxa1(im,5 )*r3(re,1,3)&
!&                    +gxa1(re,9 )*r3(im,2,3)-gxa1(im,9 )*r3(re,2,3)&
!&                    +gxa1(re,10)*r3(im,3,3)-gxa1(im,10)*r3(re,3,3))&
!&             +15.d0*(gxa1(re,8 )*r3(im,4,3)-gxa1(im,8 )*r3(re,4,3)&
!&                    +gxa1(re,3 )*r3(im,5,3)-gxa1(im,3 )*r3(re,5,3)&
!&                    +gxa1(re,4 )*r3(im,6,3)-gxa1(im,4 )*r3(re,6,3))&
!&             -1.5d0*(r11(re,3)*r12(im,3)+r11(im,3)*r12(re,3))
 rank2c(im,4)=7.5d0*(gxa1(re,5 )*r3(im,1,2)-gxa1(im,5 )*r3(re,1,2)&
& +gxa1(re,9 )*r3(im,2,2)-gxa1(im,9 )*r3(re,2,2)&
& +gxa1(re,10)*r3(im,3,2)-gxa1(im,10)*r3(re,3,2))&
& +15.d0*(gxa1(re,8 )*r3(im,4,2)-gxa1(im,8 )*r3(re,4,2)&
& +gxa1(re,3 )*r3(im,5,2)-gxa1(im,3 )*r3(re,5,2)&
& +gxa1(re,4 )*r3(im,6,2)-gxa1(im,4 )*r3(re,6,2))&
& -1.5d0*(r11(re,3)*r12(im,2)+r11(im,3)*r12(re,2))
 rank2c(im,5)=7.5d0*(gxa1(re,5 )*r3(im,1,1)-gxa1(im,5 )*r3(re,1,1)&
& +gxa1(re,9 )*r3(im,2,1)-gxa1(im,9 )*r3(re,2,1)&
& +gxa1(re,10)*r3(im,3,1)-gxa1(im,10)*r3(re,3,1))&
& +15.d0*(gxa1(re,8 )*r3(im,4,1)-gxa1(im,8 )*r3(re,4,1)&
& +gxa1(re,3 )*r3(im,5,1)-gxa1(im,3 )*r3(re,5,1)&
& +gxa1(re,4 )*r3(im,6,1)-gxa1(im,4 )*r3(re,6,1))&
& -1.5d0*(r11(re,3)*r12(im,1)+r11(im,3)*r12(re,1))
 rank2c(im,6)=7.5d0*(gxa1(re,6 )*r3(im,1,1)-gxa1(im,6 )*r3(re,1,1)&
& +gxa1(re,7 )*r3(im,2,1)-gxa1(im,7 )*r3(re,2,1)&
& +gxa1(re,8 )*r3(im,3,1)-gxa1(im,8 )*r3(re,3,1))&
& +15.d0*(gxa1(re,9 )*r3(im,4,1)-gxa1(im,9 )*r3(re,4,1)&
& +gxa1(re,4 )*r3(im,5,1)-gxa1(im,4 )*r3(re,5,1)&
& +gxa1(re,2 )*r3(im,6,1)-gxa1(im,2 )*r3(re,6,1))&
& -1.5d0*(r11(re,2)*r12(im,1)+r11(im,2)*r12(re,1))

end subroutine cont33cso
!!***


!!****f* m_contract/cont33so
!! NAME
!! cont33so
!!
!! FUNCTION
!! Contract symmetric rank 3 tensor gxa1 with symmetric rank 3 tensor
!! gxa2 using metric tensor gmet and antisymmetric tensor amet to
!! produce rank 2 real tensor.
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

subroutine cont33so(gxa1,gxa2,gmet,amet,rank2)

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

!!****f* m_contract/cont35
!! NAME
!! cont35
!!
!! FUNCTION
!! Contract symmetric rank3 tensor gxa with rank5 symmetric tensor to
!! produce symmetric rank2 tensor.
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

subroutine cont35(gxa,rank5,rank2)

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

!!****f* m_contract/metcon
!! NAME
!! metcon
!!
!! FUNCTION
!! Carries out specialized metric tensor contractions needed for
!! l=0,1,2,3 nonlocal Kleinman-Bylander pseudopotential operation.
!! Full advantage is taken of the full permutational symmetry of these
!! tensors.
!!
!! INPUTS
!!  rank=0,1,2, or 3 = rank of input tensor aa
!!  gmet(3,3)=metric tensor (array is symmetric but stored as 3x3)
!!  aa(2,(rank+1)*(rank+2)/2)=unique elements of complex input tensor
!!
!! OUTPUT
!!  bb(2,(rank+1)*(rank+2)/2)=unique elements of complex output tensor
!!
!! NOTES
!! All tensors are stored in a compressed storage mode defined below;
!! input and output conform to this scheme.
!! When tensor elements occur repeatedly due to symmetry, the
!! WEIGHT IS INCLUDED in the output tensor element to simplify later
!! contractions with other tensors of the same rank and form, i.e. the
!! next contraction is then simply a dot product over the unique elements.
!!
!! Definitions of the contractions:
!!
!! rank=0:  bb = aa  (simply copy the scalar value, no use of gmet)
!!
!! rank=1:  bb(i)= $gmet(i,l) aa(l)$  (3 elements in, 3 elements out)
!!
!! rank=2:  bb(i,j)= $[{3 \over 2} gmet(i,l) gmet(j,m) - {1 \over 2} gmet(i,j) gmet(l,m)] aa(l,m)$
!!          (6 elements in, 6 elements out)
!!
!! rank=3:  bb(i,j,k)= $[{5 \over 2} g(i,l) g(j,m) g(k,n) - {3 \over 2} g(i,j) g(l,m) g(k,n)] aa(l,m,n)$
!!          (10 elements in, 10 elements out)
!!         In this rank 3 case, the second term is NOT symmetric in all
!!         permutations of i,j,k, but the final tensor b(ijk) may be
!!         symmetrized over all permutations because it will be
!!         contracted with a completely symmetric tensor.
!!
!! The compressed storage scheme is based on storing a symmetric 3x3 matrix as
!!     (1 . .)
!!     (6 2 .)
!!     (5 4 3)
!! which leads to the following mappings for all ranks
!! where the compressed storage index is to the right of the arrow:
!! rank=0 1->1 (only a scalar)
!! rank=1 1->1 2->2 3->3 (standard vector, no compression)
!! rank=2 11->1 22->2 33->3 32->4 31->5 21->6
!!  weights   1     1     1     2     2     2
!! rank=3 111->1 221->2 331->3 321->4 311->5 211->6 222->7 332->8 322->9 333->10
!!  weights    1      3      3      6      3      3      1      3      3       1
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine metcon(rank,gmet,aa,bb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank
!arrays
 real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2),gmet(3,3)
 real(dp),intent(out) :: bb(2,((rank+1)*(rank+2))/2)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: scalar,tmpiii,tmpijk
 character(len=500) :: message
!arrays
 real(dp) :: vector(3)

! *************************************************************************

!This statement function defines the l=3 contraction in
!terms of the free indices of the contracted tensor (re and im)
! coniii(ii,i1,i2,i3)=gmet(i1,1)*gmet(i2,1)*gmet(i3,1)*aa(ii,1)+&
!& gmet(i1,2)*gmet(i2,2)*gmet(i3,2)*aa(ii,7)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,3)*aa(ii,10)
! conijk(ii,i1,i2,i3)=aa(ii,4)*&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,2)*gmet(i2,3)*gmet(i3,1)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,2)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,3)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,3))
! con(ii,i1,i2,i3)=coniii(ii,i1,i2,i3)+conijk(ii,i1,i2,i3)+&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,1)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,1)*gmet(i3,2))*aa(ii,6)+&
!& (gmet(i1,1)*gmet(i2,2)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,1)*gmet(i3,2)+&
!& gmet(i1,2)*gmet(i2,2)*gmet(i3,1))*aa(ii,2)+&
!& (gmet(i1,1)*gmet(i2,3)*gmet(i3,1)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,1)+&
!& gmet(i1,1)*gmet(i2,1)*gmet(i3,3))*aa(ii,5)+&
!& (gmet(i1,1)*gmet(i2,3)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,1)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,1))*aa(ii,3)+&
!& (gmet(i1,2)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,2)*gmet(i2,3)*gmet(i3,2)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,2))*aa(ii,9)+&
!& (gmet(i1,2)*gmet(i2,3)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,2)*gmet(i3,3)+&
!& gmet(i1,3)*gmet(i2,3)*gmet(i3,2))*aa(ii,8)

!DEBUG
!write(std_out,*)' metcon : enter '
!stop
!ENDDEBUG
 if (rank==0) then
!  Simply copy scalar, re and im
   bb(1,1)=aa(1,1)
   bb(2,1)=aa(2,1)

 else if (rank==1) then
!  Apply gmet to input vector, re and im
   do ii=1,2
     do jj=1,3
       bb(ii,jj)=gmet(jj,1)*aa(ii,1)+gmet(jj,2)*aa(ii,2)+gmet(jj,3)*aa(ii,3)
     end do
   end do


 else if (rank==2) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Carry out g(c,d)*aa(c,d) contraction to get scalar
     scalar=gmet(1,1)*aa(ii,1)+gmet(2,2)*aa(ii,2)+&
&     gmet(3,3)*aa(ii,3)+2.0d0*(gmet(2,1)*aa(ii,6)+&
&     gmet(3,1)*aa(ii,5)+gmet(3,2)*aa(ii,4) )
!    Write out components of contraction
!    (1,1)->1
     bb(ii,1)=0.5d0*(3.0d0*(gmet(1,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(1,2)*gmet(1,2)*aa(ii,2)+gmet(1,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(1,2)*gmet(1,1)+gmet(1,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(1,3)*gmet(1,1)+gmet(1,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(1,3)*gmet(1,2)+gmet(1,2)*gmet(1,3))*aa(ii,4) ) &
&     - gmet(1,1)*scalar)
!    (2,2)->2
     bb(ii,2)=0.5d0*(3.0d0*(gmet(2,1)*gmet(2,1)*aa(ii,1)+&
&     gmet(2,2)*gmet(2,2)*aa(ii,2)+gmet(2,3)*gmet(2,3)*aa(ii,3)+&
&     (gmet(2,2)*gmet(2,1)+gmet(2,1)*gmet(2,2))*aa(ii,6)+&
&     (gmet(2,3)*gmet(2,1)+gmet(2,1)*gmet(2,3))*aa(ii,5)+&
&     (gmet(2,3)*gmet(2,2)+gmet(2,2)*gmet(2,3))*aa(ii,4) )&
&     - gmet(2,2)*scalar)
!    (3,3)->3
     bb(ii,3)=0.5d0*(3.0d0*(gmet(3,1)*gmet(3,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(3,2)*aa(ii,2)+gmet(3,3)*gmet(3,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(3,1)+gmet(3,1)*gmet(3,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(3,1)+gmet(3,1)*gmet(3,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(3,2)+gmet(3,2)*gmet(3,3))*aa(ii,4) )&
&     - gmet(3,3)*scalar)
!    (3,2)->4
     bb(ii,4)=0.5d0*(3.0d0*(gmet(3,1)*gmet(2,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(2,2)*aa(ii,2)+gmet(3,3)*gmet(2,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(2,1)+gmet(3,1)*gmet(2,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(2,1)+gmet(3,1)*gmet(2,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(2,2)+gmet(3,2)*gmet(2,3))*aa(ii,4) )&
&     - gmet(3,2)*scalar)
!    (3,1)->5
     bb(ii,5)=0.5d0*(3.0d0*(gmet(3,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(3,2)*gmet(1,2)*aa(ii,2)+gmet(3,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(3,2)*gmet(1,1)+gmet(3,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(3,3)*gmet(1,1)+gmet(3,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(3,3)*gmet(1,2)+gmet(3,2)*gmet(1,3))*aa(ii,4) )&
&     - gmet(3,1)*scalar)
!    (2,1)->6
     bb(ii,6)=0.5d0*(3.0d0*(gmet(2,1)*gmet(1,1)*aa(ii,1)+&
&     gmet(2,2)*gmet(1,2)*aa(ii,2)+gmet(2,3)*gmet(1,3)*aa(ii,3)+&
&     (gmet(2,2)*gmet(1,1)+gmet(2,1)*gmet(1,2))*aa(ii,6)+&
&     (gmet(2,3)*gmet(1,1)+gmet(2,1)*gmet(1,3))*aa(ii,5)+&
&     (gmet(2,3)*gmet(1,2)+gmet(2,2)*gmet(1,3))*aa(ii,4) ) &
&     - gmet(2,1)*scalar)
!    Include appropriate weights for multiplicity
     bb(ii,4)=2.d0*bb(ii,4)
     bb(ii,5)=2.d0*bb(ii,5)
     bb(ii,6)=2.d0*bb(ii,6)
   end do

 else if (rank==3) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Carry out g(l,m)g(j,n)*aa(l,m,n) contraction to get vector(j)
     do jj=1,3
       tmpiii=   gmet(1,1)*gmet(jj,1)*aa(ii,1)+&
&       gmet(2,2)*gmet(jj,2)*aa(ii,7)+&
&       gmet(3,3)*gmet(jj,3)*aa(ii,10)
       tmpijk=  (gmet(1,2)*gmet(jj,3)+&
&       gmet(3,1)*gmet(jj,2)+&
&       gmet(2,3)*gmet(jj,1)+&
&       gmet(3,2)*gmet(jj,1)+&
&       gmet(1,3)*gmet(jj,2)+&
&       gmet(2,1)*gmet(jj,3)) *aa(ii,4)
       vector(jj)=tmpiii + tmpijk +&
&       (gmet(1,2)*gmet(jj,1)+&
&       gmet(2,1)*gmet(jj,1)+&
&       gmet(1,1)*gmet(jj,2)) *aa(ii,6)+&
&       (gmet(1,2)*gmet(jj,2)+&
&       gmet(2,1)*gmet(jj,2)+&
&       gmet(2,2)*gmet(jj,1)) *aa(ii,2)+&
&       (gmet(1,3)*gmet(jj,1)+&
&       gmet(3,1)*gmet(jj,1)+&
&       gmet(1,1)*gmet(jj,3)) *aa(ii,5)+&
&       (gmet(1,3)*gmet(jj,3)+&
&       gmet(3,1)*gmet(jj,3)+&
&       gmet(3,3)*gmet(jj,1)) *aa(ii,3)+&
&       (gmet(2,3)*gmet(jj,2)+&
&       gmet(3,2)*gmet(jj,2)+&
&       gmet(2,2)*gmet(jj,3)) *aa(ii,9)+&
&       (gmet(2,3)*gmet(jj,3)+&
&       gmet(3,2)*gmet(jj,3)+&
&       gmet(3,3)*gmet(jj,2)) *aa(ii,8)
     end do
!    Write out components of contraction
!    (111)->1
     bb(ii,1) =2.5d0*con_met(ii,1,1,1)-1.5d0*(gmet(1,1)*vector(1))
!    (221)->2
     bb(ii,2) =2.5d0*con_met(ii,2,2,1)-0.5d0*(gmet(1,2)*vector(2)+&
&     gmet(1,2)*vector(2)+gmet(2,2)*vector(1))
!    (331)->3
     bb(ii,3) =2.5d0*con_met(ii,3,3,1)-0.5d0*(gmet(1,3)*vector(3)+&
&     gmet(1,3)*vector(3)+gmet(3,3)*vector(1))
!    (321)->4
     bb(ii,4) =2.5d0*con_met(ii,3,2,1)-0.5d0*(gmet(1,3)*vector(2)+&
&     gmet(1,2)*vector(3)+gmet(3,2)*vector(1))
!    (311)->5
     bb(ii,5) =2.5d0*con_met(ii,3,1,1)-0.5d0*(gmet(1,3)*vector(1)+&
&     gmet(1,1)*vector(3)+gmet(3,1)*vector(1))
!    (211)->6
     bb(ii,6) =2.5d0*con_met(ii,2,1,1)-0.5d0*(gmet(1,2)*vector(1)+&
&     gmet(1,1)*vector(2)+gmet(2,1)*vector(1))
!    (222)->7
     bb(ii,7) =2.5d0*con_met(ii,2,2,2)-1.5d0*(gmet(2,2)*vector(2))

!    (332)->8
     bb(ii,8) =2.5d0*con_met(ii,3,3,2)-0.5d0*(gmet(2,3)*vector(3)+&
&     gmet(2,3)*vector(3)+gmet(3,3)*vector(2))
!    (322)->9
     bb(ii,9) =2.5d0*con_met(ii,3,2,2)-0.5d0*(gmet(2,3)*vector(2)+&
&     gmet(2,2)*vector(3)+gmet(3,2)*vector(2))
!    (333)->10
     bb(ii,10)=2.5d0*con_met(ii,3,3,3)-1.5d0*(gmet(3,3)*vector(3))
!    Include appropriate weights for multiplicity
     bb(ii,2)=3.d0*bb(ii,2)
     bb(ii,3)=3.d0*bb(ii,3)
     bb(ii,4)=6.d0*bb(ii,4)
     bb(ii,5)=3.d0*bb(ii,5)
     bb(ii,6)=3.d0*bb(ii,6)
     bb(ii,8)=3.d0*bb(ii,8)
     bb(ii,9)=3.d0*bb(ii,9)
   end do

 else
   write(message, '(a,i0,a,a,a)' )&
&   'Input rank=',rank,' not allowed.',ch10,&
&   'Possible values are 0,1,2,3 only.'
   MSG_BUG(message)
 end if

 contains

   function con_met(ii,i1,i2,i3)

   real(dp) :: con_met
   integer :: ii,i1,i2,i3
   real(dp)::coniii,conijk

   coniii=gmet(i1,1)*gmet(i2,1)*gmet(i3,1)*aa(ii,1)+&
&   gmet(i1,2)*gmet(i2,2)*gmet(i3,2)*aa(ii,7)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,3)*aa(ii,10)
   conijk=aa(ii,4)*&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*gmet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,3)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,3))
   con_met=coniii+conijk+&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,1)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*gmet(i3,2))*aa(ii,6)+&
&   (gmet(i1,1)*gmet(i2,2)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*gmet(i3,2)+&
&   gmet(i1,2)*gmet(i2,2)*gmet(i3,1))*aa(ii,2)+&
&   (gmet(i1,1)*gmet(i2,3)*gmet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*gmet(i3,3))*aa(ii,5)+&
&   (gmet(i1,1)*gmet(i2,3)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,1)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,1))*aa(ii,3)+&
&   (gmet(i1,2)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*gmet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,2))*aa(ii,9)+&
&   (gmet(i1,2)*gmet(i2,3)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,2)*gmet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*gmet(i3,2))*aa(ii,8)

 end function con_met

end subroutine metcon
!!***

!!****f* m_contract/metcon_so
!! NAME
!! metcon_so
!!
!! FUNCTION
!! Carries out specialized metric tensor contractions needed for
!! l=0,1,2,3 nonlocal Kleinman-Bylander pseudopotential operation
!! in the spin-orbit case.
!! Full advantage is taken of the full permutational symmetry of these tensors.
!!
!! INPUTS
!!  rank=0,1,2, or 3 = rank of input tensor aa
!!  gmet(3,3)=metric tensor (array is symmetric but stored as 3x3)
!!  amet(3,3)=real or imaginary part of one spin matrix element of the
!!            "spin metric" tensor
!!  aa(2,(rank+1)*(rank+2)/2)=unique elements of complex input tensor
!!
!! OUTPUT
!!  bb(2,(rank+1)*(rank+2)/2)=unique elements of complex output tensor
!!
!! NOTES
!! All tensors are stored in a compressed storage mode defined below;
!! input and output conform to this scheme.
!! When tensor elements occur repeatedly due to symmetry, the
!! WEIGHT IS INCLUDED in the output tensor element to simplify later
!! contractions with other tensors of the same rank and form, i.e. the
!! next contraction is then simply a dot product over the unique elements.
!!
!! Definitions of the contractions:
!!
!! rank=0:  bb=0
!!
!! rank=1:  bb(i)= $amet(i,l) aa(l)$  (3 elements in, 3 elements out)
!!
!! rank=2:  bb(i,j)= $[3 gmet(i,l) amet(j,m)] aa(l,m)$
!!          (6 elements in, 6 elements out)
!!
!! rank=3:  bb(i,j,k)= $[{15 \over 2} g(i,l) g(j,m) a(k,n) - {3 \over 2} g(i,j) g(l,m) a(k,n)] aa(l,m,n)$
!!          (10 elements in, 10 elements out)
!!          In this rank 3 case, the second term is NOT symmetric in all
!!          permutations of i,j,k, but the final tensor b(ijk) may be
!!          symmetrized over all permutations because it will be
!!          contracted with a completely symmetric tensor.
!!
!! The compressed storage scheme is based on storing
!! a symmetric 3x3 matrix as
!!     (1 . .)
!!     (6 2 .)
!!     (5 4 3)
!! which leads to the following mappings for all ranks
!! where the compressed storage index is to the right of the arrow:
!! rank=0 1->1 (only a scalar)
!! rank=1 1->1 2->2 3->3 (standard vector, no compression)
!! rank=2 11->1 22->2 33->3 32->4 31->5 21->6
!!  weights   1     1     1     2     2     2
!! rank=3 111->1 221->2 331->3 321->4 311->5 211->6 222->7 332->8 322->9 333->10
!!  weights    1      3      3      6      3      3      1      3      3       1
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine metcon_so(rank,gmet,amet,aa,bb)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: rank
!arrays
 real(dp),intent(in) :: aa(2,((rank+1)*(rank+2))/2),amet(3,3),gmet(3,3)
 real(dp),intent(out) :: bb(2,((rank+1)*(rank+2))/2)

!Local variables-------------------------------
!scalars
 integer :: ii,jj
 real(dp) :: tmpiii,tmpijk
 character(len=500) :: message
!arrays
 real(dp) :: vector(3)

! *************************************************************************

 if (rank==0) then
!  Simply copy scalar, re and im
   bb(1,1)=0.d0
   bb(2,1)=0.d0
!
 else if (rank==1) then
!  Apply gmet to input vector, re and im
   do ii=1,2
     do jj=1,3
       bb(ii,jj)=amet(jj,1)*aa(ii,1)+amet(jj,2)*aa(ii,2)+amet(jj,3)*aa(ii,3)
     end do
   end do

 else if (rank==2) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Write out components of contraction
!    (1,1)->1
     bb(ii,1)=3.0d0*(gmet(1,1)*amet(1,1)*aa(ii,1)+&
&     gmet(1,2)*amet(1,2)*aa(ii,2)+gmet(1,3)*amet(1,3)*aa(ii,3)+&
&     (gmet(1,2)*amet(1,1)+gmet(1,1)*amet(1,2))*aa(ii,6)+&
&     (gmet(1,3)*amet(1,1)+gmet(1,1)*amet(1,3))*aa(ii,5)+&
&     (gmet(1,3)*amet(1,2)+gmet(1,2)*amet(1,3))*aa(ii,4))
!    (2,2)->2
     bb(ii,2)=3.0d0*(gmet(2,1)*amet(2,1)*aa(ii,1)+&
&     gmet(2,2)*amet(2,2)*aa(ii,2)+gmet(2,3)*amet(2,3)*aa(ii,3)+&
&     (gmet(2,2)*amet(2,1)+gmet(2,1)*amet(2,2))*aa(ii,6)+&
&     (gmet(2,3)*amet(2,1)+gmet(2,1)*amet(2,3))*aa(ii,5)+&
&     (gmet(2,3)*amet(2,2)+gmet(2,2)*amet(2,3))*aa(ii,4) )
!    (3,3)->3
     bb(ii,3)=3.0d0*(gmet(3,1)*amet(3,1)*aa(ii,1)+&
&     gmet(3,2)*amet(3,2)*aa(ii,2)+gmet(3,3)*amet(3,3)*aa(ii,3)+&
&     (gmet(3,2)*amet(3,1)+gmet(3,1)*amet(3,2))*aa(ii,6)+&
&     (gmet(3,3)*amet(3,1)+gmet(3,1)*amet(3,3))*aa(ii,5)+&
&     (gmet(3,3)*amet(3,2)+gmet(3,2)*amet(3,3))*aa(ii,4) )
!    (3,2)->4
     bb(ii,4)=3.0d0*(gmet(3,1)*amet(2,1)*aa(ii,1)+&
&     gmet(3,2)*amet(2,2)*aa(ii,2)+gmet(3,3)*amet(2,3)*aa(ii,3)+&
&     (gmet(3,2)*amet(2,1)+gmet(3,1)*amet(2,2))*aa(ii,6)+&
&     (gmet(3,3)*amet(2,1)+gmet(3,1)*amet(2,3))*aa(ii,5)+&
&     (gmet(3,3)*amet(2,2)+gmet(3,2)*amet(2,3))*aa(ii,4) )
     bb(ii,4)=bb(ii,4)+3.0d0*(amet(3,1)*gmet(2,1)*aa(ii,1)+&
&     amet(3,2)*gmet(2,2)*aa(ii,2)+amet(3,3)*gmet(2,3)*aa(ii,3)+&
&     (amet(3,2)*gmet(2,1)+amet(3,1)*gmet(2,2))*aa(ii,6)+&
&     (amet(3,3)*gmet(2,1)+amet(3,1)*gmet(2,3))*aa(ii,5)+&
&     (amet(3,3)*gmet(2,2)+amet(3,2)*gmet(2,3))*aa(ii,4) )
!    (3,1)->5
     bb(ii,5)=3.0d0*(gmet(3,1)*amet(1,1)*aa(ii,1)+&
&     gmet(3,2)*amet(1,2)*aa(ii,2)+gmet(3,3)*amet(1,3)*aa(ii,3)+&
&     (gmet(3,2)*amet(1,1)+gmet(3,1)*amet(1,2))*aa(ii,6)+&
&     (gmet(3,3)*amet(1,1)+gmet(3,1)*amet(1,3))*aa(ii,5)+&
&     (gmet(3,3)*amet(1,2)+gmet(3,2)*amet(1,3))*aa(ii,4) )
     bb(ii,5)=bb(ii,5)+3.0d0*(amet(3,1)*gmet(1,1)*aa(ii,1)+&
&     amet(3,2)*gmet(1,2)*aa(ii,2)+amet(3,3)*gmet(1,3)*aa(ii,3)+&
&     (amet(3,2)*gmet(1,1)+amet(3,1)*gmet(1,2))*aa(ii,6)+&
&     (amet(3,3)*gmet(1,1)+amet(3,1)*gmet(1,3))*aa(ii,5)+&
&     (amet(3,3)*gmet(1,2)+amet(3,2)*gmet(1,3))*aa(ii,4) )
!    (2,1)->6
     bb(ii,6)=3.0d0*(gmet(2,1)*amet(1,1)*aa(ii,1)+&
&     gmet(2,2)*amet(1,2)*aa(ii,2)+gmet(2,3)*amet(1,3)*aa(ii,3)+&
&     (gmet(2,2)*amet(1,1)+gmet(2,1)*amet(1,2))*aa(ii,6)+&
&     (gmet(2,3)*amet(1,1)+gmet(2,1)*amet(1,3))*aa(ii,5)+&
&     (gmet(2,3)*amet(1,2)+gmet(2,2)*amet(1,3))*aa(ii,4) )
     bb(ii,6)=bb(ii,6)+3.0d0*(amet(2,1)*gmet(1,1)*aa(ii,1)+&
&     amet(2,2)*gmet(1,2)*aa(ii,2)+amet(2,3)*gmet(1,3)*aa(ii,3)+&
&     (amet(2,2)*gmet(1,1)+amet(2,1)*gmet(1,2))*aa(ii,6)+&
&     (amet(2,3)*gmet(1,1)+amet(2,1)*gmet(1,3))*aa(ii,5)+&
&     (amet(2,3)*gmet(1,2)+amet(2,2)*gmet(1,3))*aa(ii,4) )
!    Include appropriate weights for multiplicity
     bb(ii,4)=bb(ii,4)
     bb(ii,5)=bb(ii,5)
     bb(ii,6)=bb(ii,6)
   end do

 else if (rank==3) then
!  Apply rank2 expression, re and im
   do ii=1,2
!    Carry out g(l,m)g(j,n)*aa(l,m,n) contraction to get vector(j)
     do jj=1,3
       tmpiii=   gmet(1,1)*amet(jj,1)*aa(ii,1)+&
&       gmet(2,2)*amet(jj,2)*aa(ii,7)+&
&       gmet(3,3)*amet(jj,3)*aa(ii,10)
       tmpijk=  (gmet(1,2)*amet(jj,3)+&
&       gmet(3,1)*amet(jj,2)+&
&       gmet(2,3)*amet(jj,1)+&
&       gmet(3,2)*amet(jj,1)+&
&       gmet(1,3)*amet(jj,2)+&
&       gmet(2,1)*amet(jj,3)) *aa(ii,4)
       vector(jj)=tmpiii + tmpijk +&
&       (gmet(1,2)*amet(jj,1)+&
&       gmet(2,1)*amet(jj,1)+&
&       gmet(1,1)*amet(jj,2)) *aa(ii,6)+&
&       (gmet(1,2)*amet(jj,2)+&
&       gmet(2,1)*amet(jj,2)+&
&       gmet(2,2)*amet(jj,1)) *aa(ii,2)+&
&       (gmet(1,3)*amet(jj,1)+&
&       gmet(3,1)*amet(jj,1)+&
&       gmet(1,1)*amet(jj,3)) *aa(ii,5)+&
&       (gmet(1,3)*amet(jj,3)+&
&       gmet(3,1)*amet(jj,3)+&
&       gmet(3,3)*amet(jj,1)) *aa(ii,3)+&
&       (gmet(2,3)*amet(jj,2)+&
&       gmet(3,2)*amet(jj,2)+&
&       gmet(2,2)*amet(jj,3)) *aa(ii,9)+&
&       (gmet(2,3)*amet(jj,3)+&
&       gmet(3,2)*amet(jj,3)+&
&       gmet(3,3)*amet(jj,2)) *aa(ii,8)
     end do
!    Write out components of contraction
!    (111)->1
     bb(ii,1) =7.5d0*con_metso(ii,1,1,1)-1.5d0*(gmet(1,1)*vector(1))
!    (221)->2
     bb(ii,2) =7.5d0*con_metso(ii,2,2,1)-0.5d0*(gmet(1,2)*vector(2)+&
&     gmet(1,2)*vector(2)+gmet(2,2)*vector(1))
!    (331)->3
     bb(ii,3) =7.5d0*con_metso(ii,3,3,1)-0.5d0*(gmet(1,3)*vector(3)+&
&     gmet(1,3)*vector(3)+gmet(3,3)*vector(1))
!    (321)->4
     bb(ii,4) =7.5d0*con_metso(ii,3,2,1)-0.5d0*(gmet(1,3)*vector(2)+&
&     gmet(1,2)*vector(3)+gmet(3,2)*vector(1))
!    (311)->5
     bb(ii,5) =7.5d0*con_metso(ii,3,1,1)-0.5d0*(gmet(1,3)*vector(1)+&
&     gmet(1,1)*vector(3)+gmet(3,1)*vector(1))
!    (211)->6
     bb(ii,6) =7.5d0*con_metso(ii,2,1,1)-0.5d0*(gmet(1,2)*vector(1)+&
&     gmet(1,1)*vector(2)+gmet(2,1)*vector(1))
!    (222)->7
     bb(ii,7) =7.5d0*con_metso(ii,2,2,2)-1.5d0*(gmet(2,2)*vector(2))

!    (332)->8
     bb(ii,8) =7.5d0*con_metso(ii,3,3,2)-0.5d0*(gmet(2,3)*vector(3)+&
&     gmet(2,3)*vector(3)+gmet(3,3)*vector(2))
!    (322)->9
     bb(ii,9) =7.5d0*con_metso(ii,3,2,2)-0.5d0*(gmet(2,3)*vector(2)+&
&     gmet(2,2)*vector(3)+gmet(3,2)*vector(2))
!    (333)->10
     bb(ii,10)=7.5d0*con_metso(ii,3,3,3)-1.5d0*(gmet(3,3)*vector(3))
!    Include appropriate weights for multiplicity
     bb(ii,2)=3.d0*bb(ii,2)
     bb(ii,3)=3.d0*bb(ii,3)
     bb(ii,4)=6.d0*bb(ii,4)
     bb(ii,5)=3.d0*bb(ii,5)
     bb(ii,6)=3.d0*bb(ii,6)
     bb(ii,8)=3.d0*bb(ii,8)
     bb(ii,9)=3.d0*bb(ii,9)
   end do

 else
   write(message, '(a,i0,a,a,a)' )&
&   'Input rank=',rank,' not allowed.',ch10,&
&   'Possible values are 0,1,2,3 only.'
   MSG_BUG(message)
 end if

 contains

!This function defines the l=3 contraction in
!terms of the free indices of the contracted tensor (re and im)

   function cona_metso(ii,i1,i2,i3)

   real(dp) :: cona_metso
   integer,intent(in) :: ii,i1,i2,i3
   real(dp) :: coniii, conijk

   coniii=gmet(i1,1)*gmet(i2,1)*amet(i3,1)*aa(ii,1)+&
&   gmet(i1,2)*gmet(i2,2)*amet(i3,2)*aa(ii,7)+&
&   gmet(i1,3)*gmet(i2,3)*amet(i3,3)*aa(ii,10)
   conijk=aa(ii,4)*&
&   (gmet(i1,1)*gmet(i2,2)*amet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*amet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*amet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*amet(i3,1)+&
&   gmet(i1,1)*gmet(i2,3)*amet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*amet(i3,3))
   cona_metso=coniii+conijk+&
&   (gmet(i1,1)*gmet(i2,2)*amet(i3,1)+&
&   gmet(i1,2)*gmet(i2,1)*amet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*amet(i3,2))*aa(ii,6)+&
&   (gmet(i1,1)*gmet(i2,2)*amet(i3,2)+&
&   gmet(i1,2)*gmet(i2,1)*amet(i3,2)+&
&   gmet(i1,2)*gmet(i2,2)*amet(i3,1))*aa(ii,2)+&
&   (gmet(i1,1)*gmet(i2,3)*amet(i3,1)+&
&   gmet(i1,3)*gmet(i2,1)*amet(i3,1)+&
&   gmet(i1,1)*gmet(i2,1)*amet(i3,3))*aa(ii,5)+&
&   (gmet(i1,1)*gmet(i2,3)*amet(i3,3)+&
&   gmet(i1,3)*gmet(i2,1)*amet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*amet(i3,1))*aa(ii,3)+&
&   (gmet(i1,2)*gmet(i2,2)*amet(i3,3)+&
&   gmet(i1,2)*gmet(i2,3)*amet(i3,2)+&
&   gmet(i1,3)*gmet(i2,2)*amet(i3,2))*aa(ii,9)+&
&   (gmet(i1,2)*gmet(i2,3)*amet(i3,3)+&
&   gmet(i1,3)*gmet(i2,2)*amet(i3,3)+&
&   gmet(i1,3)*gmet(i2,3)*amet(i3,2))*aa(ii,8)
 end function cona_metso


   function con_metso(ii,i1,i2,i3)

   real(dp) :: con_metso
   integer,intent(in) :: ii,i1,i2,i3

   con_metso=(cona_metso(ii,i3,i1,i2)+cona_metso(ii,i2,i3,i1)+cona_metso(ii,i1,i2,i3))/3.d0

 end function con_metso

end subroutine metcon_so
!!***

!!****f* m_contract/metric_so
!! NAME
!! metric_so
!!
!! FUNCTION
!! Computes Pauli matrices and antisymmetric tensor needed for
!! spin-orbit.
!!
!! INPUTS
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space
!!              (bohr**-1)
!!
!! OUTPUT
!!  amet(2,3,3,2,2)=the antisymmetric tensor A(Re/Im,y,y'',s,s'')
!!  pauli(2,2,2,3)=Pauli matrixes
!!
!! NOTES
!! (the preprocessing based on CPP does not allow isolated quotes,
!!  so that they have been doubled in the following latex equations)
!!{{\ \begin{eqnarray}
!! $2 Pauli(s,s'',1) & = & 0 & 1 \nonumber
!!                   &   & 1 & 0 \nonumber
!!  2 Pauli(s,s'',2) & = & 0 & -i \nonumber
!!                   &   & i &  0 \nonumber
!!  2 Pauli(s,s'',3) & = & 1 &  0 \nonumber
!!                   &   & 0 & -1
!! \end{eqnarray} }}
!!{{\ \begin{eqnarray}
!!    $Amet(y,y'',s,s'') & = & -i Pauli(s,s'',n) E(n,m,m'') Gprimd(m,y) Gprimd(m'',y'')
!! \end{eqnarray} }}
!!
!! E(n,m,m''):   full antisymetric tensor
!! s,s'':        spin indices (1..2)
!! y,y'',m,m'',n: metric indices (1..3)
!! a,b:         strain indices (1..3)
!! Amet and Pauli are complex
!!
!! PARENTS
!!      nonlop_pl
!!
!! CHILDREN
!!
!! SOURCE

subroutine metric_so(amet,gprimd,pauli)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: amet(2,3,3,2,2),pauli(2,2,2,3)

!Local variables-------------------------------
!scalars
 integer :: iy1,iy2,m1,m2,n
!arrays
 real(dp) :: buffer1(3,3,2,2) !,buffer2(3,3,3,3,2,2)

! **********************************************************************

!Fill in Pauli matrices and make them spin matrices:
!Pauli(Re/Im,up,down,coord):
 pauli(:,:,:,:)=0.d0
 pauli(1,1,2,1)= 1.d0;pauli(1,2,1,1)= 1.d0
 pauli(2,1,2,2)=-1.d0;pauli(2,2,1,2)= 1.d0
 pauli(1,1,1,3)= 1.d0;pauli(1,2,2,3)=-1.d0
 pauli(:,:,:,:)= 0.5d0*pauli(:,:,:,:)

!Construct the antisymmetric tensor:
 amet(:,:,:,:,:)=0.d0
 do iy2=1,3
   do iy1=1,3
     do n=1,3
       m1=mod(n ,3)+1    !  n,m1,m2 is an even permutation
       m2=mod(m1,3)+1
       amet(1:2,iy1,iy2,1:2,1:2) = amet(:,iy1,iy2,:,:) &
&       + pauli(:,:,:,n) &
&       *(gprimd(m1,iy1)*gprimd(m2,iy2) &
&       -gprimd(m2,iy1)*gprimd(m1,iy2))
     end do
   end do
 end do
!amet= -i amet
 buffer1(:,:,:,:)=-amet(1,:,:,:,:)
 amet(1,:,:,:,:) = amet(2,:,:,:,:)
 amet(2,:,:,:,:) =buffer1(:,:,:,:)

!DEBUG
!!Eventually construct the gradients of Amet wrt strains:
!! DAmet(y,y'',a,b,s,s'')= d[Amet(y,y'',s,s'')]/d[strain(a,b)]
!!                        -i Pauli(s,s'',n)*
!!                          ( E(n,a,m'')*Gprimd(b,y )*Gprimd(m'',y'')
!!                           +E(n,m,a )*Gprimd(b,y'')*Gprimd(m ,y ) )
!damet(:,:,:,:,:,:,:)=0.d0
!do ib=1,3
!do ia=1,3
!m1=mod(ia,3)+1    !  ia,m1,m2 is an even permutation
!m2=mod(m1,3)+1
!do iy2=1,3
!do iy1=1,3
!damet(:,iy1,iy2,ia,ib,:,:) = damet(:,iy1,iy2,ia,ib,:,:) &
!&        + (pauli(:,:,:,m2)*gprimd(m1,iy2) &
!-pauli(:,:,:,m1)*gprimd(m2,iy2))*gprimd(ib,iy1) &
!&        + (pauli(:,:,:,m1)*gprimd(m2,iy1) &
!-pauli(:,:,:,m2)*gprimd(m1,iy1))*gprimd(ib,iy2)
!end do
!end do
!end do
!end do
!! damet= i damet
!buffer2(:,:,:,:,:,:)= damet(1,:,:,:,:,:,:)
!damet(1,:,:,:,:,:,:)= -damet(2,:,:,:,:,:,:)
!damet(2,:,:,:,:,:,:)=buffer2(:,:,:,:,:,:)
!! Symetrize damet(:,:,:,a,b,:,:)
!damet(:,:,:,1,2,:,:)=0.5d0*(damet(:,:,:,1,2,:,:)+damet(:,:,:,2,1,:,:))
!damet(:,:,:,1,3,:,:)=0.5d0*(damet(:,:,:,1,3,:,:)+damet(:,:,:,3,1,:,:))
!damet(:,:,:,2,3,:,:)=0.5d0*(damet(:,:,:,2,3,:,:)+damet(:,:,:,3,2,:,:))
!damet(:,:,:,2,1,:,:)=damet(:,:,:,1,2,:,:)
!damet(:,:,:,3,1,:,:)=damet(:,:,:,1,3,:,:)
!damet(:,:,:,3,2,:,:)=damet(:,:,:,2,3,:,:)
!ENDDEBUG

end subroutine metric_so
!!***

end module m_contract
!!***

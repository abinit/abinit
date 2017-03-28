!{\src2tex{textfont=tt}}
!!****f* ABINIT/symsgtetra
!! NAME
!! symsgtetra
!!
!! FUNCTION
!! Yields all the TETRAGONAL symmetry operations starting from the space group symbol.
!! according to the International Tables of Crystallography, 1983.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! msym = default number of symmetries
!! nsym = the number of symmetry operations
!! shubnikov= magnetic type of the space group to be generated
!! spgorig = the origin choice (1 or 2) for the axes system
!! spgroup = the numeric symbol of the space groups
!! spgroupma= number of the magnetic space group
!!
!! OUTPUT
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym) = 3D matrix containg symmetry operations
!! tnons(3,nsym) = 2D matrix containing translations associated
!!
!! PARENTS
!!      gensymspgr
!!
!! CHILDREN
!!      bldgrp,spgdata
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symsgtetra'
 use interfaces_41_geometry, except_this_one => symsgtetra
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym,shubnikov,spgaxor,spgorig,spgroup
 integer,intent(in) :: spgroupma
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables ------------------------------
!integer :: isym
!scalars
 integer :: nogen,sporder
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
!arrays
 integer :: gen4m(3,3),genmmm(3,3),genmmp(3,3),genmpm(3,3),genmpp(3,3)
 integer :: genpmm(3,3),genpmp(3,3),genppm(3,3)

! *************************************************************************
!DEBUG
!write(std_out,*) ' symsgtetra: ',spgroup,shubnikov,spgroupma
!ENDDEBUG
 nogen=0

 tnons(:,:)=zero
!The identity operation belongs to all space groups
 symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1

!The next operation belongs to most TETRAGONAL space groups, except 81,82,111,121.
 symrel(:,:,2)=0 ; symrel(1,2,2)=-1 ; symrel(2,1,2)=1 ; symrel(3,3,2)=1

!Predefine some generators
 genmpp(:,:)=0 ; genmpp(1,1)=-1 ; genmpp(2,2)= 1 ; genmpp(3,3)= 1
 genpmp(:,:)=0 ; genpmp(1,1)= 1 ; genpmp(2,2)=-1 ; genpmp(3,3)= 1
 genppm(:,:)=0 ; genppm(1,1)= 1 ; genppm(2,2)= 1 ; genppm(3,3)=-1
 genpmm(:,:)=0 ; genpmm(1,1)= 1 ; genpmm(2,2)=-1 ; genpmm(3,3)=-1
 genmpm(:,:)=0 ; genmpm(1,1)=-1 ; genmpm(2,2)= 1 ; genmpm(3,3)=-1
 genmmp(:,:)=0 ; genmmp(1,1)=-1 ; genmmp(2,2)=-1 ; genmmp(3,3)= 1
 genmmm(:,:)=0 ; genmmm(1,1)=-1 ; genmmm(2,2)=-1 ; genmmm(3,3)=-1
 gen4m(:,:)=0  ; gen4m(1,2)=-1 ; gen4m(2,1)=1 ; gen4m(3,3)=-1

!Default non-magnetic behaviour
 symafm(1:nsym)=1


!assigns the generators to each space group
 select case (spgroup)
!  TETRAGONAL space groups
 case (75,79,83,87)        !P4, I4, P4/m, I4/m
   nogen=2
 case (76,80)                !P41, I41
   tnons(:,2)=(/0.d0,0.d0,0.25d0/)
   nogen=2
 case (77,84)                !P42, P42/m
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   nogen=2
 case (78)                !P43
   tnons(:,2)=(/0.d0,0.d0,0.75d0/)
   nogen=2
 case (81,82)                !PB4, IB4
   symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1
   symrel(:,:,2)=0 ; symrel(1,1,2)=-1; symrel(2,2,2)=-1; symrel(3,3,2)=1
   symrel(:,:,3)=0 ; symrel(1,2,3)=1 ; symrel(2,1,3)=-1 ; symrel(3,3,3)=-1
   symrel(:,:,4)=0 ; symrel(1,2,4)=-1 ; symrel(2,1,4)=1 ; symrel(3,3,4)=-1
   nogen=0
   if (shubnikov==3) then
     symafm(3:4)=-1
   end if
 case (85,86,88)                !P4/n
   symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1 ; symrel(:,:,5) = -1*symrel(:,:,1)
   symrel(:,:,2)=0 ; symrel(1,1,2)=-1; symrel(2,2,2)=-1; symrel(3,3,2)=1 ; symrel(:,:,6) = -1*symrel(:,:,2)
   symrel(:,:,3)=0 ; symrel(1,2,3)=-1 ; symrel(2,1,3)=1 ; symrel(3,3,3)=1 ; symrel(:,:,7) = -1*symrel(:,:,3)
   symrel(:,:,4)=0 ; symrel(1,2,4)=1 ; symrel(2,1,4)=-1 ; symrel(3,3,4)=1 ; symrel(:,:,8) = -1*symrel(:,:,4)
   nogen=0
   if(spgorig==1) then
     select case (spgroup)
     case (85)                !P4/n
       tnons(:,3)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     case (86)                !P42/n
       tnons(:,3)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (88)                 !I41/a
       tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,3)=(/0.0d0,0.5d0,0.25d0/)
       tnons(:,4)=(/0.5d0,0.0d0,0.75d0/)
       tnons(:,5)=(/0.0d0,0.5d0,0.25d0/)
       tnons(:,6)=(/0.5d0,0.0d0,0.75d0/)
       tnons(:,7)=(/0.5d0,0.5d0,0.5d0/)
     end select
   else
     select case (spgroup)
     case (85)          !P4/n
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ;  tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
       tnons(:,3)=(/0.5d0,0.0d0,0.d0/) ;  tnons(:,7)=(/0.5d0,0.0d0,0.d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.d0/) ;  tnons(:,8)=(/0.0d0,0.5d0,0.d0/)
     case (86)          !P42/n
       tnons(:,2)=(/0.5d0,0.5d0,0.0d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.0d0/)
       tnons(:,3)=(/0.0d0,0.5d0,0.5d0/) ; tnons(:,7)=(/0.0d0,0.5d0,0.5d0/)
       tnons(:,4)=(/0.5d0,0.0d0,0.5d0/) ; tnons(:,8)=(/0.5d0,0.0d0,0.5d0/)
     case (88)          !I41/a
       tnons(:,2)=(/0.5d0,0.0d0,0.5d0/) ;  tnons(:,6)=(/0.5d0,0.0d0,0.5d0/)
       tnons(:,3)=(/0.75d0,0.25d0,0.25d0/)
       tnons(:,7)=(/0.25d0,0.75d0,0.75d0/)
       tnons(:,4)=(/0.75d0,0.75d0,0.75d0/)
       tnons(:,8)=(/0.25d0,0.25d0,0.25d0/)
     end select
   end if
   if (shubnikov==3) then
     select case (spgroupma)
     case(61,69,83)
       symafm(3)=-1;symafm(4)=-1;symafm(7)=-1; symafm(8)=-1
     case(62,70,84)
       symafm(5)=-1;symafm(6)=-1;symafm(7)=-1; symafm(8)=-1
     case(63,71,85)
       symafm(3)=-1;symafm(5)=-1;symafm(4)=-1; symafm(6)=-1
     end select
   end if
 case (89,97)                !P422, I422
   symrel(:,:,3) = genmpm(:,:)
   nogen=3
 case (90)                !P4212
   tnons(:,2)=(/0.5d0,0.5d0,0.0d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (91)                !P4122
   tnons(:,2)=(/0.d0,0.d0,0.25d0/)
   symrel(:,:,3) = genmpm(:,:)
   nogen=3
 case (92)                !P41212
   tnons(:,2)=(/0.5d0,0.5d0,0.25d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.25d0/)
   nogen=3
 case (93)                !P4222
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   nogen=3
 case (94)                !P42212
   tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (95)                !P4322
   tnons(:,2)=(/0.d0,0.d0,0.75d0/)
   symrel(:,:,3) = genmpm(:,:)
   nogen=3
 case (96)                !P43212
   tnons(:,2)=(/0.5d0,0.5d0,0.75d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.75d0/)
   nogen=3
 case (98)                !I4122
   tnons(:,2)=(/0.d0,0.5d0,0.25d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.0d0,0.75d0/)
   nogen=3
 case (99,107,123,139)        !P4mm, I4mm, P4/mmm, I4/mmm
   symrel(:,:,3) = genpmp(:,:)
   nogen=3
 case (100)                !P4bm
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (101,132)                !P42cm, P42/mcm
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 case (102)                !P42nm
   tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (103,124)                !P4cc, P4/mcc
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 case (104)                !P4nc
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (105,131)                !P42mc, P42/mmc
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   nogen=3
 case (106,135)                !P42bc, P42/mbc
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (108)                !I4cm
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 case (109)                !I41md
   tnons(:,2)=(/0.d0,0.5d0,0.25d0/)
   symrel(:,:,3) = genpmp(:,:)
   symrel(:,:,4) = genmmp(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
   nogen=4
 case (110)                !I41cd
   tnons(:,2)=(/0.d0,0.5d0,0.25d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genmmp(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
   nogen=4
 case (111,121)                !PB42m, IB42m
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genmpm(:,:)
   nogen=3
 case (112)                !PB42c
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 case (113)                !PB421m
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (114)                !PB421c
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (115,119)                !PB4m2,IB4m2
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genpmp(:,:)
   nogen=3
 case (116,120)                !PB4c2, IB4c2
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 case (117)                !PB4b2
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (118)                !PB4n2
   symrel(:,:,2) = gen4m(:,:)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (122)                !IB42d
   symrel(:,:,2)=0 ; symrel(1,2,2)=-1 ; symrel(2,1,2)=1 ; symrel(3,3,2)=-1
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.d0,0.75d0/)
   nogen=3
 case (125,126,129,130,133,134,137,138,141,142)
   symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1
   symrel(:,:,2)=0 ; symrel(1,1,2)=-1; symrel(2,2,2)=-1; symrel(3,3,2)=1
   symrel(:,:,3)=0 ; symrel(1,2,3)=-1 ; symrel(2,1,3)=1 ; symrel(3,3,3)=1
   symrel(:,:,4)=0 ; symrel(1,2,4)=1 ; symrel(2,1,4)=-1 ; symrel(3,3,4)=1
   symrel(:,:,5) = genmpm(:,:)
   symrel(:,:,6) = genmmm(:,:)
   if (spgorig==1) then
     select case (spgroup)
     case (125)                !P4/nbm
       tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     case (126)         !P4/nnc
       tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (129)                !P4/nmm
       tnons(:,3)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     case (130)                !P4/ncc
       tnons(:,3)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     case (133)         !P42/nbc
       tnons(:,3)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,5)=(/0.d0,0.d0,0.5d0/)
       tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (134)         !P42/nnm
       tnons(:,3)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (137)         !P42/nmc
       tnons(:,3)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (138)         !P42/ncm
       tnons(:,3)=(/0.5d0,0.5d0,0.5d0/) ; tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,5)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     case (141)         !I41/amd
       tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,3)=(/0.d0,0.5d0,0.25d0/)
       tnons(:,4)=(/0.5d0,0.d0,0.75d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.75d0/)
       tnons(:,6)=(/0.d0,0.5d0,0.25d0/)
     case (142)         !I41/acd
       tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
       tnons(:,3)=(/0.d0,0.5d0,0.25d0/)
       tnons(:,4)=(/0.5d0,0.d0,0.75d0/)
       tnons(:,5)=(/0.5d0,0.d0,0.25d0/)
       tnons(:,6)=(/0.d0,0.5d0,0.25d0/)
     end select
   else
     select case (spgroup)
     case (125)         !P4/nbm
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.d0/)
     case (126)                !P4/nnc
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.5d0/)
     case (129)                !P4/nmm
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.d0/) ;  tnons(:,5)=(/0.0d0,0.5d0,0.d0/)
     case (130)         !P4/ncc
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ;  tnons(:,3)=(/0.5d0,0.d0,0.d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.d0/) ; tnons(:,5)=(/0.0d0,0.5d0,0.5d0/)
     case (133)        !P42/nbc
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.5d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.d0/)
     case (134)                !P42/nnm
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.5d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.5d0/)
     case (137)         !P42/nmc
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.5d0/) ; tnons(:,5)=(/0.d0,0.5d0,0.d0/)
     case (138)         !P42/ncm
       tnons(:,2)=(/0.5d0,0.5d0,0.d0/) ; tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
       tnons(:,4)=(/0.0d0,0.5d0,0.5d0/) ; tnons(:,5)=(/0.d0,0.5d0,0.5d0/)
     case (141)         !I41/amd
       tnons(:,2)=(/0.5d0,0.d0,0.5d0/) ; tnons(:,3)=(/0.25d0,0.75d0,0.25d0/)
       tnons(:,4)=(/0.25d0,0.25d0,0.75d0/) ;  tnons(:,5)=(/0.5d0,0.d0,0.5d0/)
     case (142)         !I41/acd
       tnons(:,2)=(/0.5d0,0.d0,0.5d0/) ; tnons(:,3)=(/0.25d0,0.75d0,0.25d0/)
       tnons(:,4)=(/0.25d0,0.25d0,0.75d0/) ; tnons(:,5)=(/0.5d0,0.d0,0.d0/)
     end select
   end if
   nogen=6
 case (127)                !P4/mbm
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   nogen=3
 case (128)                !P4/mnc
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (136)                !P42/mnm
   tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   nogen=3
 case (140)                !I4/mcm
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   nogen=3
 end select

 if(shubnikov==3)then
   select case(spgroupma)
   case(3,9,15,21,27,31,45,47,53,55,77,79,&
&     89,97,105,113,121,129,137,145,153,159,166,174,182,190,198,206,&
&     214,222,230,236,242,248,254,262,270,278,285,293,301,309,317,&
&     323,330,336,343,346,355,358,391,392,403,404,&
&     439,442,451,454,487,490,499,500,535,&
&     538,545,546)
     symafm(2)=-1
   case(90,98,106,114,122,130,138,146,154,160,167,175,183,191,199,&
&     207,215,223,231,237,243,249,255,263,271,279,287,295,303,311,319,&
&     325,331,337,345,347,357,359,389,393,401,405,441,443,453,455,&
&     489,491,497,501,537,539,543,547)
     symafm(3)=-1
   case(91,99,107,115,123,131,139,147,155,161,165,173,181,189,&
&     197,205,213,221,229,235,241,247,253,261,269,277,286,294,302,310,&
&     318,324,329,335,342,344,354,356,390,394,402,406,&
&     438,440,450,452,486,488,498,502,534,536,544,548)
     symafm(2:3)=-1
   case(365,377,413,425,461,473,509,521)
     symafm(2)=-1
     symafm(5:6)=-1
   case(366,378,414,426,462,474,510,522,554,564)
     symafm(3:5)=-1
   case(367,379,415,427,463,475,511,523,555,565)
     symafm(3:4)=-1
   case(368,380,416,428,464,476,512,524,556,566)
     symafm(3:4)=-1
     symafm(6)=-1
   case(369,381,417,429,465,477,513,525)
     symafm(2)=-1
     symafm(5)=-1
   case(370,382,418,430,466,478,514,526,558,568)
     symafm(3:6)=-1
   case(371,383,419,431,467,479,515,527,559,569)
     symafm(6)=-1
   case(553,563)
     symafm(5:6)=-1
   case(557,567)
     symafm(5)=-1
   end select
 end if

!DEBUG
!write(std_out,*)' symsgtetra : symrel(:,:,6)=',symrel(:,:,6)
!ENDDEBUG

 if (nogen>1) then
   call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
 end if

 call spgdata(brvsb,intsb,intsbl,ptintsb,&
& ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

!DEBUG
!write(std_out,*)' symsgtetra : exit'
!do isym=1,nsym
!write(std_out,'(i3,2x,9i3,3es13.3,i3)') isym,symrel(:,:,isym),tnons(:,isym),symafm(isym)
!end do
!ENDDEBUG

end subroutine symsgtetra
!!***

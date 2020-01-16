!!****m* ABINIT/m_symsg
!! NAME
!!  m_symsg
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (RC,XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_symsg

 use defs_basis
 use m_errors
 use m_abicore

 use m_time,     only : timab
 use m_spgdata,  only : spgdata

 implicit none

 private
!!***

 public :: symsgcube
 public :: symsghexa
 public :: symsgmono
 public :: symsgortho
 public :: symsgtetra
!!***

contains
!!***

!!****f* m_symsg/symsgcube
!! NAME
!! symsgcube
!!
!! FUNCTION
!! Generate all the symmetry operations starting from the space group symbol
!! for the cubic groups (according to the International Tables of Crystallography, 1983)
!!
!! INPUTS
!! msym = default number of symmetries
!! shubnikov= magnetic type of the space group to be generated
!! spgaxor = the possible orientation of the axes system
!! spgorig = the origin choice (1 or 2) for the axes system
!! spgroup = the numeric symbol of the space groups
!! spgroupma = number of the magnetic space group
!!
!! OUTPUT
!! nsym = the number of symmetry operations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym) = 3D matrix containg symmetry operations
!! tnons(3,nsym) = 2D matrix containing translations associated
!!
!! PARENTS
!!      gensymspgr
!!
!! CHILDREN
!!      bldgrp,spgdata,timab
!!
!! SOURCE

subroutine symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma,symafm,symrel,tnons)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,shubnikov,spgaxor,spgorig,spgroup,spgroupma
 integer,intent(inout) :: nsym !vz_i
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(out) :: tnons(3,msym)

!Local variables -----------------------------
!scalars
 integer :: ii,nogen,sporder
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
!arrays
 integer :: gen1(3,3),gen2(3,3),gen3(3,3),gen4(3,3),gen5(3,3),gen6(3,3)
 integer :: gen7(3,3),gen8(3,3),gen9(3,3),genmmm(3,3),genmmp(3,3),genmpm(3,3)
 integer :: genmpp(3,3),genpmm(3,3),genpmp(3,3),genppm(3,3),genrot(3,3)
 integer :: genswm(3,3),genswp(3,3)
 real(dp) :: tsec(2)

!*************************************************************************

!DEBUG
!write(std_out,*) ' symsgcube : enter with spgroup ',spgroup,' and ',' origin choice ',spgorig
!write(std_out,*) ' msym,nsym=',msym,nsym
!ENDDEBUG


!The identity operation belongs to all space groups
 symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1

!Initialize the associated translations matrix to 0
 do ii=1,msym
   tnons(:,ii)= 0.0d0
 end do
 nogen=0

!Predefine some generators
 genmpp(:,:)=0 ; genmpp(1,1)=-1 ; genmpp(2,2)= 1 ; genmpp(3,3)= 1
 genpmp(:,:)=0 ; genpmp(1,1)= 1 ; genpmp(2,2)=-1 ; genpmp(3,3)= 1
 genppm(:,:)=0 ; genppm(1,1)= 1 ; genppm(2,2)= 1 ; genppm(3,3)=-1
 genpmm(:,:)=0 ; genpmm(1,1)= 1 ; genpmm(2,2)=-1 ; genpmm(3,3)=-1
 genmpm(:,:)=0 ; genmpm(1,1)=-1 ; genmpm(2,2)= 1 ; genmpm(3,3)=-1
 genmmp(:,:)=0 ; genmmp(1,1)=-1 ; genmmp(2,2)=-1 ; genmmp(3,3)= 1
 genmmm(:,:)=0 ; genmmm(1,1)=-1 ; genmmm(2,2)=-1 ; genmmm(3,3)=-1
 genrot(:,:)=0 ; genrot(1,3)=1 ; genrot(3,2)=1 ; genrot(2,1)=1  !reshape((/0,0,1,1,0,0,0,1,0/),(/3,3/),(/0,0/),(/2,1/) )
 genswm(:,:)=0 ; genswm(2,1)=1 ; genswm(1,2)=1 ; genswm(3,3)=-1 !reshape((/0,1,0,1,0,0,0,0,-1/),(/3,3/),(/0,0/),(/2,1/) )
 genswp(:,:)=0 ; genswp(2,1)=1 ; genswp(1,2)=1 ; genswp(3,3)=1  !reshape((/0,1,0,1,0,0,0,0,1/),(/3,3/),(/0,0/),(/2,1/) )

!Because of problems with the IBM compiler, that does not like reshape
!operations, define 9 basic matrices
 gen1(:,:)=0 ; gen1(1,1)=1
 gen2(:,:)=0 ; gen2(1,2)=1
 gen3(:,:)=0 ; gen3(1,3)=1
 gen4(:,:)=0 ; gen4(2,1)=1
 gen5(:,:)=0 ; gen5(2,2)=1
 gen6(:,:)=0 ; gen6(2,3)=1
 gen7(:,:)=0 ; gen7(3,1)=1
 gen8(:,:)=0 ; gen8(3,2)=1
 gen9(:,:)=0 ; gen9(3,3)=1

!Default non-magnetic behaviour
 symafm(1:msym)=1

!*************************************************************************

!Treat CUBIC groups

 if(195<=spgroup .and. spgroup<=230)then

   select case(spgroup)
   case (195,196,197)        !P23, F23, I23
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genpmm(:,:)
     nogen=5
   case (200,202,204)                !PmB3, FmB3, ImB3
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     nogen=4
   case (198,199)                !P213, I213
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genpmm(:,:)
     tnons(:,5)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,6) = gen3(:,:)-gen4(:,:)-gen8(:,:)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,7) = -gen3(:,:)-gen4(:,:)+gen8(:,:)
     tnons(:,7)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,8) = -gen3(:,:)+gen4(:,:)-gen8(:,:)
     tnons(:,8)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,10) = -gen2(:,:)+gen6(:,:)-gen7(:,:)
     tnons(:,10)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,11) =  gen2(:,:)-gen6(:,:)-gen7(:,:)
     tnons(:,11)=(/0.5d0,0.5d0,0.0d0/)
     symrel(:,:,12) = -gen2(:,:)-gen6(:,:)+gen7(:,:)
     tnons(:,12)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,9) =  gen2(:,:)+gen6(:,:)+gen7(:,:)
   case (201)                !Pn-3
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       symrel(:,:,4) = genmpm(:,:)
       symrel(:,:,5) = genmmm(:,:)
       tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
       nogen=5
       if(shubnikov==3)symafm(5)=-1
       call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
       symrel(:,:,5) = genpmm(:,:)
       tnons(:,5)=(/0.d0,0.5d0,0.5d0/)
       symrel(:,:,6) =  gen3(:,:)-gen4(:,:)-gen8(:,:)
       tnons(:,6)=(/0.d0,0.5d0,0.5d0/)
       symrel(:,:,7) = -gen3(:,:)-gen4(:,:)+gen8(:,:)
       tnons(:,7)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,8) = -gen3(:,:)+gen4(:,:)-gen8(:,:)
       tnons(:,8)=(/0.5d0,0.d0,0.5d0/)
       symrel(:,:,10) = -gen2(:,:)+gen6(:,:)-gen7(:,:)
       tnons(:,10)=(/0.5d0,0.d0,0.5d0/)
       symrel(:,:,11) =  gen2(:,:)-gen6(:,:)-gen7(:,:)
       tnons(:,11)=(/0.d0,0.5d0,0.5d0/)
       symrel(:,:,12) = -gen2(:,:)-gen6(:,:)+gen7(:,:)
       tnons(:,12)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,9) =  gen2(:,:)+gen6(:,:)+gen7(:,:)
       do ii=1,12
         symrel(:,:,ii+12)= - symrel(:,:,ii)
         tnons(:,ii+12)=tnons(:,ii)
         if(shubnikov==3)symafm(ii+12)=-1
       end do
     end if
     nogen=0
   case (205,206)                !Pa-3,Ia-3
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genpmm(:,:)
     tnons(:,5)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,6) =  gen3(:,:)-gen4(:,:)-gen8(:,:)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,7) = -gen3(:,:)-gen4(:,:)+gen8(:,:)
     tnons(:,7)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,8) = -gen3(:,:)+gen4(:,:)-gen8(:,:)
     tnons(:,8)=(/0.0d0,0.5d0,0.5d0/)
     symrel(:,:,9) =  gen2(:,:)+gen6(:,:)+gen7(:,:)
     symrel(:,:,10) = -gen2(:,:)+gen6(:,:)-gen7(:,:)
     tnons(:,10)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,11) =  gen2(:,:)-gen6(:,:)-gen7(:,:)
     tnons(:,11)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,12) = -gen2(:,:)-gen6(:,:)+gen7(:,:)
     tnons(:,12)=(/0.5d0,0.d0,0.5d0/)
     nogen=0
   case (203)                !FdB3
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       symrel(:,:,4) = genmpm(:,:)
       nogen=4
       nsym=12
       call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
       do ii=1,12
         symrel(:,:,ii+12) = - symrel(:,:,ii)
         tnons(:,ii+12)=(/0.25d0,0.25d0,0.25d0/)
         if(shubnikov==3)symafm(ii+12)=-1
       end do
       nsym=24
       nogen=0
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.25d0,0.25d0,0.d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.25d0,0.d0,0.25d0/)
       symrel(:,:,5) = genpmm(:,:)
       tnons(:,5)=(/0.d0,0.25d0,0.25d0/)
       symrel(:,:,6) =  gen3(:,:)-gen4(:,:)-gen8(:,:)
       tnons(:,6)=(/0.d0,0.25d0,0.25d0/)
       symrel(:,:,7) = -gen3(:,:)-gen4(:,:)+gen8(:,:)
       tnons(:,7)=(/0.25d0,0.25d0,0.d0/)
       symrel(:,:,8) = -gen3(:,:)+gen4(:,:)-gen8(:,:)
       tnons(:,8)=(/0.25d0,0.d0,0.25d0/)
       symrel(:,:,10) = -gen2(:,:)+gen6(:,:)-gen7(:,:)
       tnons(:,10)=(/0.25d0,0.d0,0.25d0/)
       symrel(:,:,11) =  gen2(:,:)-gen6(:,:)-gen7(:,:)
       tnons(:,11)=(/0.d0,0.25d0,0.25d0/)
       symrel(:,:,12) = -gen2(:,:)-gen6(:,:)+gen7(:,:)
       tnons(:,12)=(/0.25d0,0.25d0,0.d0/)
       symrel(:,:,9) =  gen2(:,:)+gen6(:,:)+gen7(:,:)
       do ii=1,12
         symrel(:,:,ii+12) = - symrel(:,:,ii)
         tnons(:,ii+12)=tnons(:,ii)
         if(shubnikov==3)symafm(ii+12)=-1
       end do
     end if
     nsym=24
     nogen=0
   case (207,209,211)        !P432, F432, I432, PmB3m, FmB3m, ImB3m
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswm(:,:)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (208)                !P4232
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (210)                !F4132
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.75d0,0.25d0,0.75d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (212)                !P4332
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.25d0,0.75d0,0.75d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (213,214)                !P4132, I4132
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.75d0,0.25d0,0.25d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (215,216,217)        !PB43m, FB43m, IB43m
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswp(:,:)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (218,219)                !PB43n, FB343c
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswp(:,:)
     tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (220)                !IB43d
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genswp(:,:)
     tnons(:,5)=(/0.25d0,0.25d0,0.25d0/)
     if(shubnikov==3)symafm(5)=-1
     nogen=5
   case (221,225,229)        !Pm3m,Fm3m,Im3m
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswm(:,:)
     if(shubnikov==3)then
       if(spgroupma==94 .or. spgroupma==95 .or. spgroupma==118 .or. &
&       spgroupma==119 .or. spgroupma==142 .or. spgroupma==143   )symafm(5)=-1
     end if
     nogen=5
   case (222)                !PnB3n
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.d0,0.d0,0.5d0/)
       if(shubnikov==3)then
         if(spgroupma==100 .or. spgroupma==101)symafm(5)=-1
       end if
       nogen=5
       nsym=24
       call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
       do ii=1,24
         symrel(:,:,ii+24)=-symrel(:,:,ii)
         tnons(:,ii+24)=tnons(:,ii)
         if(shubnikov==3)then
           if(spgroupma==100 .or. spgroupma==102)symafm(ii+24)=-symafm(ii)
           if(spgroupma==101)symafm(ii+24)=symafm(ii)
         end if
       end do
       nogen=0 ; nsym=48
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       symrel(:,:,4) = genmpm(:,:)
       symrel(:,:,5) = genswm(:,:)
       if(shubnikov==3)then
         if(spgroupma==100 .or. spgroupma==101)symafm(5)=-1
       end if
       nogen=5
       nsym=24
       call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
       do ii=1,24
         symrel(:,:,ii+24)=-symrel(:,:,ii)
         tnons(:,ii+24)=(/0.5d0,0.5d0,0.5d0/)
         if(shubnikov==3)then
           if(spgroupma==100 .or. spgroupma==102)symafm(ii+24)=-symafm(ii)
           if(spgroupma==101)symafm(ii+24)=symafm(ii)
         end if
       end do
       nogen=0 ; nsym=48
     end if
   case (223,226)           ! PmB3n, FmB3c
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     symrel(:,:,4) = genmpm(:,:)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
     if(shubnikov==3)then
       if(spgroupma==106 .or. spgroupma==124 .or. &
&       spgroupma==107 .or. spgroupma==125     )symafm(5)=-1
     end if
     nogen=5
   case (224)                !PnB3m
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       symrel(:,:,4) = genmpm(:,:)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
       symrel(:,:,6) = genmmm(:,:)
       tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.5d0,0.5d0,0.0d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.5d0,0.0d0,0.5d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,6) = genmmm(:,:)
     end if
     if(shubnikov==3)then
       if(spgroupma==112 .or. spgroupma==113)symafm(5)=-1
       if(spgroupma==112 .or. spgroupma==114)symafm(6)=-1
     end if
     nogen=6
   case (227)                !FdB3m
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.75d0,0.25d0,0.75d0/)
       symrel(:,:,6) = genmmm(:,:)
       tnons(:,6)=(/0.25d0,0.25d0,0.25d0/)
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.75d0,0.25d0,0.5d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.25d0,0.5d0,0.75d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.75d0,0.25d0,0.5d0/)
       symrel(:,:,6) = genmmm(:,:)
     end if
     if(shubnikov==3)then
       if(spgroupma==130 .or. spgroupma==131)symafm(5)=-1
       if(spgroupma==130 .or. spgroupma==132)symafm(6)=-1
     end if
     nogen=6
   case (228)                !FdB3c
     if (spgorig==1) then
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.75d0,0.25d0,0.75d0/)
       symrel(:,:,6) = genmmm(:,:)
       tnons(:,6)=(/0.75d0,0.75d0,0.75d0/)
     else
       symrel(:,:,2) = genrot(:,:)
       symrel(:,:,3) = genmmp(:,:)
       tnons(:,3)=(/0.25d0,0.75d0,0.5d0/)
       symrel(:,:,4) = genmpm(:,:)
       tnons(:,4)=(/0.75d0,0.5d0,0.25d0/)
       symrel(:,:,5) = genswm(:,:)
       tnons(:,5)=(/0.75d0,0.25d0,0.d0/)
       symrel(:,:,6) = genmmm(:,:)
     end if
     if(shubnikov==3)then
       if(spgroupma==136 .or. spgroupma==137)symafm(5)=-1
       if(spgroupma==136 .or. spgroupma==138)symafm(6)=-1
     end if
     nogen=6
   case (230)                !IaB3d
     symrel(:,:,2) = genrot(:,:)
     symrel(:,:,3) = genmmp(:,:)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     symrel(:,:,4) = genmpm(:,:)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     symrel(:,:,5) = genswm(:,:)
     tnons(:,5)=(/0.75d0,0.25d0,0.25d0/)
     if(shubnikov==3)then
       if(spgroupma==147 .or. spgroupma==148)symafm(5)=-1
     end if
     nogen=5
   end select

!  End CUBIC space groups
 end if

!***************************************************************************

 call timab(47,1,tsec)

 if (nogen>0) then
   call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
 end if

 call timab(47,2,tsec)

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

end subroutine symsgcube
!!***

!!****f* m_symsg/symsghexa
!! NAME
!! symsghexa
!!
!! FUNCTION
!! Yields all the TRIGONAL & HEXAGONAL symmetry operations starting from the space group symbol.
!! according to the International Tables of Crystallography, 1983.
!!
!! INPUTS
!! msym = default number of symmetries
!! nsym = the number of symmetry operations
!! shubnikov= magnetic type of the space group to be generated
!! spgaxor = ossible orientation of the axes system
!! spgorig = the origin choice (1 or 2) for the axes system
!! spgroup = the numeric symbol of the space groups
!! spgroupma= number of the magnetic space group
!!
!! OUTPUT
!! brvltt = bravais lattice type, here, only for rhombohedral groups
!!  with hexagonal axes
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

subroutine symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma,symafm,symrel,tnons)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma
 integer,intent(inout) :: brvltt !vz_i
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables -----------------------------
!scalars
 integer :: ii,nogen,sporder
 real(dp),parameter :: fivesixth=5.0d0/6.0d0,twothird=2.0d0/3.0d0
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
!arrays
 integer :: genm(3,3),genmmp(3,3),genswm(3,3),genswmmm(3,3),genswmmp(3,3)
 integer :: genswp(3,3)

!*************************************************************************

!DEBUG
!write(std_out,*) 'symsghexa',spgroup,shubnikov,spgroupma
!ENDDEBUG

!The identity operation belongs to all space groups
 symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1

!Predefine some generators
 genswm(:,:)=0 ; genswm(2,1)=1 ; genswm(1,2)=1 ; genswm(3,3)=-1 !reshape((/0,1,0,1,0,0,0,0,-1/),(/3,3/),(/0,0/),(/2,1/) )
 genswmmm(:,:)=0 ; genswmmm(2,1)=-1 ; genswmmm(1,2)=-1 ; genswmmm(3,3)=-1
!reshape((/0,-1,0,-1,0,0,0,0,-1/),(/3,3/),(/0,0/),(/2,1/) )
 genswmmp(:,:)=0 ; genswmmp(2,1)=-1 ; genswmmp(1,2)=-1 ; genswmmp(3,3)=1
!reshape((/0,-1,0,-1,0,0,0,0,1/),(/3,3/),(/0,0/),(/2,1/) )
 genswp(:,:)=0 ; genswp(2,1)=1 ; genswp(1,2)=1 ; genswp(3,3)=1
!reshape((/0,1,0,1,0,0,0,0,1/),(/3,3/),(/0,0/),(/2,1/) )
 genmmp(:,:)=0 ; genmmp(1,1)=-1 ; genmmp(2,2)=-1 ; genmmp(3,3)=1
!reshape((/-1,0,0,0,-1,0,0,0,1/),(/3,3/),(/0,0/),(/2,1/) )

!Initialize the associated translations matrix to 0
 do ii=1,nsym
   tnons(:,ii)= 0.0d0
 end do
 nogen=0

!Default non-magnetic behaviour
 symafm(1:nsym)=1

!*************************************************************************

!Treat TRIGONAL case
 if(143<=spgroup .and. spgroup<=167)then

!  The hexagonal axis choice (orientation) is first treated
   if (spgaxor == 1) then

!    This matrix is common to ALL trigonal spatial groups in this orientation
!    (Note : this is the 3- symmetry operation)
     symrel(:,:,2)=0 ; symrel(1,1,2)=-1 ; symrel(1,2,2)=1 ; symrel(2,1,2)=-1 ; symrel(3,3,2)=1
!    reshape((/-1,1,0,-1,0,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )

!    Assigns the generators to each space group
     select case (spgroup)
     case (143,146,147,148)        !P3, R3, PB3, RB3
       symrel(:,:,3)=0 ; symrel(1,2,3)=-1 ; symrel(2,1,3)=1 ; symrel(2,2,3)=-1 ; symrel(3,3,3)=1
!        reshape((/0,-1,0,1,-1,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )
       nogen=0
     case (144)                !P31
       tnons(:,2)=(/0.d0,0.d0,twothird/)
       symrel(:,:,3)=0 ; symrel(1,2,3)=-1 ; symrel(2,1,3)=1 ; symrel(2,2,3)=-1 ; symrel(3,3,3)=1
!        reshape((/0,-1,0,1,-1,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )
       tnons(:,3)=(/0.d0,0.d0,third/)
       nogen=0
     case (145)                !P32
       tnons(:,2)=(/0.d0,0.d0,third/)
       symrel(:,:,3)=0 ; symrel(1,2,3)=-1 ; symrel(2,1,3)=1 ; symrel(2,2,3)=-1 ; symrel(3,3,3)=1
!        reshape((/0,-1,0,1,-1,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )
       tnons(:,3)=(/0.d0,0.d0,twothird/)
       nogen=0
     case (149)                !P312
       symrel(:,:,3) = genswmmm(:,:)
       nogen=3
     case (150,155)                !P321, R32
       symrel(:,:,3) = genswm(:,:)
       nogen=3
     case (151)                !P3112
       tnons(:,2)=(/0.d0,0.d0,twothird/)
       symrel(:,:,3) = genswmmm(:,:)
       nogen=3
     case (152)                !P3121
       tnons(:,2)=(/0.d0,0.d0,twothird/)
       symrel(:,:,3) = genswm(:,:)
       nogen=3
     case (153)                !P3212
       tnons(:,2)=(/0.d0,0.d0,third/)
       symrel(:,:,3) = genswmmm(:,:)
       nogen=3
     case (154)                !P3221
       tnons(:,2)=(/0.d0,0.d0,third/)
       symrel(:,:,3) = genswm(:,:)
       nogen=3
     case (156,160,164,166)        !P3m1, R3m, PB3m1, RB3m
       symrel(:,:,3) = genswmmp(:,:)
       nogen=3
     case (157,162)                !P31m, PB31m
       symrel(:,:,3) = genswp(:,:)
       nogen=3
     case (158,161,165,167)        !P3c1, R3c, PB3c1, RB3c
       symrel(:,:,3) = genswmmp(:,:)
       tnons(:,3)=(/0.d0,0.d0,0.5d0/)
       nogen=3
     case (159,163)                !P31c, PB31c
       symrel(:,:,3) = genswp(:,:)
       tnons(:,3)=(/0.d0,0.d0,0.5d0/)
       nogen=3
     end select

     select case (spgroup)
     case (146,148,155,160,166,167)
       brvltt=7
     end select

!    Quite simple, because the generator of even order is always the third one.
     if(shubnikov==3)then
       select case(spgroupma)
       case (23,27,31,35,39,43,47,51,55,59,63,67,71,76,77,82,83,88,89,94,95,&
&         100,101,106,107)
         symafm(3)=-1
       end select
     end if

   else if (spgaxor == 2) then
!    The rhombohedral axis choice (orientation) is now treated
     write(std_out,*)'rhombohedral axes'
!    Assignment of common three-fold rotation
     symrel(:,:,2)=0 ; symrel(1,3,2)=1 ; symrel(3,2,2)=1 ; symrel(2,1,2)=1
!    reshape((/0,0,1,1,0,0,0,1,0/),(/3,3/),(/0,0/),(/2,1/) )
     symrel(:,:,3)=0 ; symrel(3,1,3)=1 ; symrel(2,3,3)=1 ; symrel(1,2,3)=1
!    reshape((/0,1,0,0,0,1,1,0,0/), (/3,3/), (/0,0/), (/2,1/) )

     select case (spgroup)
     case (146,148)       !R3
     case (155,166)       !R32, RB3m
       symrel(:,:,4) = genswmmm(:,:)
       nogen=4
     case (160)           !R3m
       symrel(:,:,4) = genswp(:,:)
       nogen=4
     case (161,167)       !R3c, RB3c
       symrel(:,:,4) = genswp(:,:)
       tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
       nogen=4
     end select

     if(shubnikov==3)then
       select case(spgroupma)
       case (47,67,71,99,101,106,107)
         symafm(4)=-1
       end select
     end if

!    End selection of axis orientation
   end if

!  End trigonal groups
 end if

!*************************************************************************

!Treat HEXAGONAL case
 if(168<=spgroup .and. spgroup<=194)then

!  This matrix (6) is common to most hexagonal spatial groups, except 174,187,188,189,190
   symrel(:,:,2)=0 ; symrel(1,1,2)=1 ; symrel(1,2,2)=-1 ; symrel(2,1,2)=1 ; symrel(3,3,2)=1
!  reshape((/1,-1,0,1,0,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )
!  This one (6 bar) is present in the other cases
   genm(:,:)=0 ; genm(1,2)=-1 ; genm(2,1)=1 ; genm(2,2)=-1 ; genm(3,3)=-1
!  reshape((/0,-1,0,1,-1,0,0,0,-1/), (/3,3/), (/0,0/), (/2,1/) )
   select case(spgroup)
   case (168,175)        !P6
     nogen=2
   case (169)                !P61
     tnons(:,2)=(/0.d0,0.d0,sixth/)
     nogen=2
   case (170)                !P65
     tnons(:,2)=(/0.d0,0.d0,fivesixth/)
     nogen=2
   case (171)                !P62
     tnons(:,2)=(/0.d0,0.d0,third/)
     nogen=2
   case (172)                !P64
     tnons(:,2)=(/0.d0,0.d0,twothird/)
     nogen=2
   case (173,176)                !P63, P63/m
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     nogen=2
   case (174)                !PB6
     symrel(:,:,2) = genm(:,:)
     nogen=2
   case (177)                !P622
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (178)                !P6122
     tnons(:,2)=(/0.d0,0.d0,sixth/)
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (179)                !P6522
     tnons(:,2)=(/0.d0,0.d0,fivesixth/)
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (180)                !P6222
     tnons(:,2)=(/0.d0,0.d0,third/)
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (181)                !P6422
     tnons(:,2)=(/0.d0,0.d0,twothird/)
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (182)                !P6322
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     symrel(:,:,3) =  genswm(:,:)
     nogen=3
   case (183,191)                !P6mm, P6/mmm
     symrel(:,:,3) = genswp(:,:)
     nogen=3
   case (184,192)                !P6cc, P6/mcc
     symrel(:,:,3) = genswp(:,:)
     tnons(:,3)=(/0.d0,0.d0,0.5d0/)
     nogen=3
   case (185,193)                !P63cm, P63/mcm
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     symrel(:,:,3) = genswp(:,:)
     nogen=3
   case (186,194)                !P63mc, P63/mmc
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     symrel(:,:,3) = genswp(:,:)
     tnons(:,3)=(/0.d0,0.d0,0.5d0/)
     nogen=3
   case (187)                !PB6m2
     symrel(:,:,2)=0 ; symrel(1,2,2)=-1 ; symrel(2,1,2)=1 ; symrel(2,2,2)=-1 ; symrel(3,3,2)=-1
!      reshape((/0,-1,0,1,-1,0,0,0,-1/), (/3,3/), (/0,0/), (/2,1/) )
     symrel(:,:,3)=0 ; symrel(1,1,3)=-1 ; symrel(1,2,3)=1 ; symrel(2,2,3)=1 ; symrel(3,3,3)=1
!      reshape((/-1,1,0,0,1,0,0,0,1/), (/3,3/), (/0,0/), (/2,1/) )
     nogen=3
     if (shubnikov==3) then
       if (spgroupma==211) symafm(2:3)=-1
       if (spgroupma==212) symafm(2)=-1
       if (spgroupma==213) symafm(3)=-1
     end if
   case (188)                !PB6c2
     symrel(:,:,2) = genm(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     symrel(:,:,3) = genswmmm(:,:)
     nogen=3
   case (189)                !PB62m
     symrel(:,:,2) = genm(:,:)
     symrel(:,:,3) = genswp(:,:)
     nogen=3
   case (190)                !PB62c
     symrel(:,:,2) = genm(:,:)
     symrel(:,:,3) = genswp(:,:)
     tnons(:,3)=(/0.d0,0.d0,0.5d0/)
     nogen=3
   end select

   if(shubnikov==3)then
     select case(spgroupma)
!      spgroup from 168 to 176 are OK, 177 to 194 are not done
     case (111,115,119,123,127,131,135,139,141,145,147,152,158,164,170,&
&       176,182,187,193,199,205,217,224,230,237,239,247,249,257,259,267,269)
       symafm(2)=-1
     case(153,159,165,171,177,183,189,195,&
&       201,207,219,225,231,240,241,250,251,260,261,270,271)
       symafm(3)=-1
     case(151,157,163,169,175,181,188,194,200,206,218,223,229,236,238,246,248,256,258,266,268)
       symafm(2:3)=-1
     end select
   end if

   call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

!  End HEXAGONAL groups
 end if

!***************************************************************************

!DEBUG
!write(std_out,*) 'symsghexa : out with nogen = ',nogen
!ENDDEBUG


 if (nogen>0) then
   call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
 end if

!DEBUG
!write(std_out,*)'symrel:'
!write(std_out,*) symrel(:,:,1:nsym)
!ENDDEBUG

end subroutine symsghexa
!!***

!!****f* m_symsg/symsgmono
!! NAME
!! symsgmono
!!
!! FUNCTION
!! Yields all the MONOCLINIC symmetry operations starting from the space group symbol.
!! according to the International Tables of Crystallography, 1983.
!! It solves also the problem of the axes orientation
!! according to the spgaxor
!!
!! INPUTS
!! msym = default number of symmetries
!! nsym = the number of symmetry operations
!! shubnikov= magnetic type of the space group to be generated
!! spgaxor = the orientation choice of the unit cell
!! spgorig = possible origin of the axes system
!! spgroup = the numeric symbol of the space groups
!! spgroupma= number of the magnetic space group
!!
!! OUTPUT
!! brvltt = bravais lattice type, here, only for rhombohedral groups
!!  with hexagonal axes (1=P; 2=I; 3=F; 4=C; 5=A; 6=B; 7=R)
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym) = 3D matrix containg symmetry operations
!! tnons(3,nsym) = 2D matrix containing translations associated
!!
!! PARENTS
!!      gensymspgr
!!
!! CHILDREN
!!      spgdata
!!
!! SOURCE

subroutine symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma,symafm,symrel,tnons)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma
 integer,intent(inout) :: brvltt !vz_i
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables -----------------------------
!scalars
 integer :: sporder
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
!arrays
 integer :: genmmm(3,3),genmmp(3,3),genmpm(3,3),genmpp(3,3),genpmm(3,3)
 integer :: genpmp(3,3),genppm(3,3)

! *************************************************************************
!the identity operation belonging to all space groups
 symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1

!Predefine some generators
 genmpp(:,:)=0 ; genmpp(1,1)=-1 ; genmpp(2,2)= 1 ; genmpp(3,3)= 1
 genpmp(:,:)=0 ; genpmp(1,1)= 1 ; genpmp(2,2)=-1 ; genpmp(3,3)= 1
 genppm(:,:)=0 ; genppm(1,1)= 1 ; genppm(2,2)= 1 ; genppm(3,3)=-1
 genpmm(:,:)=0 ; genpmm(1,1)= 1 ; genpmm(2,2)=-1 ; genpmm(3,3)=-1
 genmpm(:,:)=0 ; genmpm(1,1)=-1 ; genmpm(2,2)= 1 ; genmpm(3,3)=-1
 genmmp(:,:)=0 ; genmmp(1,1)=-1 ; genmmp(2,2)=-1 ; genmmp(3,3)= 1
 genmmm(:,:)=0 ; genmmm(1,1)=-1 ; genmmm(2,2)=-1 ; genmmm(3,3)=-1

!Default non-magnetic behaviour
 symafm(1:nsym)=1

!assigns the generators to each space group
 select case (spgroup)
 case (3)                 ! P2
   select case (spgaxor)
   case (1)                ! 3:b, P2_b = P2
     symrel(:,:,2) = genmpm(:,:)
   case (2)                ! 3:a, P2_a = P2
     symrel(:,:,2) = genpmm(:,:)
   case (3)                ! 3:c, P2_c = P2
     symrel(:,:,2) = genmmp(:,:)
   end select
   if(shubnikov==3)symafm(2)=-1
 case (4)                ! P21
   select case (spgaxor)
   case (1)                ! 3:b, P21_b = P21
     symrel(:,:,2) = genmpm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
   case (2)                ! 3:a, P21_a = P21
     symrel(:,:,2) = genpmm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   case (3)                ! 3:c, P21_c = P21
     symrel(:,:,2) = genmmp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   end select
   if(shubnikov==3)symafm(2)=-1
 case (5)                ! C2
   select case (spgaxor)
   case (1)                ! 5:b1, C2  = C2
     symrel(:,:,2) = genmpm(:,:)
     brvltt=4
   case (2)                ! 5:a1, B2_a = C2
     symrel(:,:,2) = genpmm(:,:)
     brvltt=6
   case (3)                ! 5:a2, C2_a = C2
     symrel(:,:,2) = genpmm(:,:)
     brvltt=4
   case (4)                ! 5:a3, I2_a = C2
     symrel(:,:,2) = genpmm(:,:)
     brvltt=2
   case (5)                ! 5:b2, A2_b = C2
     symrel(:,:,2) = genmpm(:,:)
     brvltt=5
   case (6)                ! 5:b3, I2_b = C2
     symrel(:,:,2) = genmpm(:,:)
     brvltt=2
   case (7)                ! 5:c1, A2_c = C2
     symrel(:,:,2) = genmmp(:,:)
     brvltt=5
   case (8)                ! 5:c2, B2_c = C2
     symrel(:,:,2) = genmmp(:,:)
     brvltt=6
   case (9)                ! 5:c3, I2_c = C2
     symrel(:,:,2) = genmmp(:,:)
     brvltt=2
   end select
   if(shubnikov==3)symafm(2)=-1
 case (6)                ! Pm
   select case (spgaxor)
   case (1)                ! Pm_b = Pm
     symrel(:,:,2) = genpmp(:,:)
   case (2)                ! Pm_a = Pm
     symrel(:,:,2) = genmpp(:,:)
   case (3)                ! Pm_c = Pm
     symrel(:,:,2) = genppm(:,:)
   end select
   if(shubnikov==3)symafm(2)=-1
 case (7)                ! Pc
   select case (spgaxor)
   case (1)                ! 7:b1, Pc_b = Pc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   case (2)                ! 7:a1, Pb_a = Pc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
   case (3)                ! 7:a2, Pn_a = Pc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
   case (4)                ! 7:a3, Pc_a = Pc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   case (5)                ! 7:b2, Pn_b = Pc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   case (6)                ! 7:b3, Pa_b = Pc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   case (7)                ! 7:c1, Pa_c = Pc
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   case (8)                ! 7:c2, Pn_c = Pc
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
   case (9)                ! 7:c3, Pb_c = Pb
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
   end select
   if(shubnikov==3)symafm(2)=-1
 case (8)                ! Cm
   select case (spgaxor)
   case (1)                ! 8:b1, Cm = Cm
     symrel(:,:,2) = genpmp(:,:)
     brvltt=4
   case (2)                ! 8:a1, Bm_a = Cm
     symrel(:,:,2) = genmpp(:,:)
     brvltt=6
   case (3)                ! 8:a2, Cm_a = Cm
     symrel(:,:,2) = genmpp(:,:)
     brvltt=4
   case (4)                ! 8:a3, Im_a = Cm
     symrel(:,:,2) = genmpp(:,:)
     brvltt=2
   case (5)                ! 8:b2, Am_b = Cm
     symrel(:,:,2) = genpmp(:,:)
     brvltt=5
   case (6)                ! 8:b3, Im_b = Cm
     symrel(:,:,2) = genpmp(:,:)
     brvltt=2
   case (7)                ! 8:c1, Am_c = Cm
     symrel(:,:,2) = genppm(:,:)
     brvltt=5
   case (8)                ! 8:c2, Bm_c = Bm
     symrel(:,:,2) = genppm(:,:)
     brvltt=6
   case (9)                ! 8:c3, Im_c = Cm
     symrel(:,:,2) = genppm(:,:)
     brvltt=2
   end select
   if(shubnikov==3)symafm(2)=-1
 case (9)                ! Cc
   select case (spgaxor)
   case (1)                ! 9:b1, Cc_b = Cc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     brvltt=4
   case (2)                ! 9:a1, Bb_a = Cc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     brvltt=6
   case (3)                ! 9:a2, Cn_a = Cc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     brvltt=4
   case (4)                ! 9:a3, Ic_a = Cc
     symrel(:,:,2) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     brvltt=2
   case (5)                ! 9:b2, An_b = Cc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
     brvltt=5
   case (6)                ! 9:b3, Ia_b = Cc
     symrel(:,:,2) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     brvltt=2
   case (7)                ! 9:c1, Aa_c = Cc
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     brvltt=5
   case (8)                ! 9:c2, B(b+c)_c = Cc
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     brvltt=6
   case (9)                ! 9:c3, Ib_c = Cc
     symrel(:,:,2) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     brvltt=2
   end select
   if(shubnikov==3)symafm(2)=-1
 case (10)                ! P2/m
   select case (spgaxor)
   case (1)                ! 10:b, P2/m = P2/m
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
   case (2)                ! 10:a, P2/m_a = P2/m
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
   case (3)                ! 10:c, P2/m_c = P2/m
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==44)symafm(4)=1
     if(spgroupma==45)symafm(2)=1
     if(spgroupma==46)symafm(3)=1
   end if
 case (11)                ! P21/m
   select case (spgaxor)
   case (1)                ! 11:b, P21/m = P21/m
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
   case (2)                ! 11:a, P21/m_a = P21/m
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
   case (3)                ! 11:c, P21/m_c = P21/m
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.d0,0.5d0/)
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==52)symafm(4)=1
     if(spgroupma==53)symafm(2)=1
     if(spgroupma==54)symafm(3)=1
   end if
 case (12)                ! C2/m
   select case (spgaxor)
   case (1)                ! 12:b1, C2/m = C2/m
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     brvltt=4
   case (2)                ! 12:a1, B2/m_a = C2/m
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     brvltt=6
   case (3)                ! 12:a2, C2/m_a = C2/m
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     brvltt=4
   case (4)                ! 12:a3, I2/m_a = C2/m
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     brvltt=2
   case (5)                ! 12:b2, A2/m_b = C2/m
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     brvltt=5
   case (6)                ! 12:b3, I2/m_b = C2/m
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     brvltt=2
   case (7)                ! 12:c1, A2/m_c = C2/m
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     brvltt=5
   case (8)                ! 12:c2, B2/m_c = B2/m
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     brvltt=6
   case (9)                ! 12:c3, I2/m_c = C2/m
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     brvltt=2
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==60)symafm(4)=1
     if(spgroupma==61)symafm(2)=1
     if(spgroupma==62)symafm(3)=1
   end if
 case (13)                ! P2/c
   select case (spgaxor)
   case (1)                ! 13:b1, P2/c = P2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.d0,0.5d0/)
   case (2)                ! 13:a1, P2/b_a = P2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
   case (3)                ! 13:a2, P2/n_a = P2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
   case (4)                ! 13:a3, P2/c_a = P2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.d0,0.5d0/)
   case (5)                ! 13:b2, P2/n_b = P2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
   case (6)                ! 13:b3, P2/a_b = P2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
   case (7)                ! 13:c1, P2/a_c = P2/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
   case (8)                ! 13:c2, P2/n_c = P2/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
   case (9)                ! 13:c3, P2/b_c = P2/b
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==67)symafm(4)=1
     if(spgroupma==68)symafm(2)=1
     if(spgroupma==69)symafm(3)=1
   end if
 case (14)              ! P21/c
   select case (spgaxor)
   case (1)             ! 14:b1, P21/c_b = P21/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
   case (2)             ! 14:a1, P21/a_b = P21/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
   case (3)                ! 14:a2, P21/n_a = P21/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
   case (4)                ! 14:a3, P21/c_a = P21/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
   case (5)                ! 14:b2, P21/n_b = P21/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
   case (6)                ! 14:b3, P21/a_b = P21/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
   case (7)                ! 14:c1, P21/a_c = P21/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
   case (8)                ! 14:c2, P21/n_c = P21/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
   case (9)                ! 14/c3, P21/b_c = P21/b
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==77)symafm(4)=1
     if(spgroupma==78)symafm(2)=1
     if(spgroupma==79)symafm(3)=1
   end if
 case (15)                ! C2/c
   select case (spgaxor)
   case (1)                ! 15:b1, C2/c_b = C2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.d0,0.5d0/)
     brvltt = 4
   case (2)                ! 15:a1, B2/b_a = C2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
     brvltt = 6
   case (3)                ! 15:a2, C2/n_a = C2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     brvltt = 4
   case (4)                ! 15:a3, I2/c_a = C2/c
     symrel(:,:,2) = genpmm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genmpp(:,:)
     tnons(:,2)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.d0,0.5d0/)
     brvltt = 2
   case (5)                ! 15:b2, A2/n_b = C2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
     brvltt = 5
   case (6)                ! 15:b3, I2/a_b = C2/c
     symrel(:,:,2) = genmpm(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genpmp(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
     brvltt = 2
   case (7)                ! 15:c1, A2/a_c = C2/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
     brvltt = 5
   case (8)                ! 15:c2, B21/b_c = C2/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     brvltt = 6
   case (9)                ! 15:c3, I2/b_c = C2/c
     symrel(:,:,2) = genmmp(:,:)
     symrel(:,:,3) = genmmm(:,:)
     symrel(:,:,4) = genppm(:,:)
     tnons(:,2)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
     brvltt = 2
   end select
   if(shubnikov==3)then
     symafm(2:4)=-1 ! Default
     if(spgroupma==87)symafm(4)=1
     if(spgroupma==88)symafm(2)=1
     if(spgroupma==89)symafm(3)=1
   end if
 end select

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

end subroutine symsgmono
!!***

!!****f* m_symsg/symsgortho
!! NAME
!! symsgortho
!!
!! FUNCTION
!! Yields all the ORTHORHOMBIC symmetry operations starting from the space group symbol.
!! It deals only with the orthorhombic groups
!! taken in the standard orientation
!! according to the International Tables of Crystallography, 1983.
!!
!! INPUTS
!! msym = default number of symmetries
!! nsym = the number of symmetry operations
!! shubnikov= magnetic type of the space group to be generated
!! spgorig = the origin choice (1 or 2) for the axes system
!! spgaxor = the possible orientation of the axes system
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

subroutine symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym,shubnikov,spgaxor,spgorig,spgroup
 integer,intent(in) :: spgroupma
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables ------------------------------
! nogen = number of generators selected
!scalars
 integer :: nogen,sporder
 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
!arrays
 integer :: genmmm(3,3),genmmp(3,3),genmpm(3,3),genmpp(3,3),genpmm(3,3)
 integer :: genpmp(3,3),genppm(3,3)

! *************************************************************************

!DEBUG
!write(std_out,*)'symsgortho ( orthorhombic groups) : enter with space group ',spgroup
!ENDDEBUG

!The orientation of the space group:
!first we will permute the input coordinates of the atoms, xred
!then we will make the calculation in the "normal" space group
!then the coordinates are reoriented to match the initial orientation
!and finally the symrel is reoriented to correspond to the new orientation
!further all the calculations are performed into the space group
!with the user-defined orientation

!The identity operation belongs to all space groups
 symrel(:,:,1)=0 ; symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1

 nogen=4

!Predefine some generators
 genmpp(:,:)=0 ; genmpp(1,1)=-1 ; genmpp(2,2)= 1 ; genmpp(3,3)= 1
 genpmp(:,:)=0 ; genpmp(1,1)= 1 ; genpmp(2,2)=-1 ; genpmp(3,3)= 1
 genppm(:,:)=0 ; genppm(1,1)= 1 ; genppm(2,2)= 1 ; genppm(3,3)=-1
 genpmm(:,:)=0 ; genpmm(1,1)= 1 ; genpmm(2,2)=-1 ; genpmm(3,3)=-1
 genmpm(:,:)=0 ; genmpm(1,1)=-1 ; genmpm(2,2)= 1 ; genmpm(3,3)=-1
 genmmp(:,:)=0 ; genmmp(1,1)=-1 ; genmmp(2,2)=-1 ; genmmp(3,3)= 1
 genmmm(:,:)=0 ; genmmm(1,1)=-1 ; genmmm(2,2)=-1 ; genmmm(3,3)=-1


!For all the groups in this routine symrel(:,:,2) is the same
 symrel(:,:,2)=genmmp(:,:)

!Default non-magnetic behaviour
 symafm(1:nsym)=1

!DEBUG
!write(std_out,*) 'symsgortho:',spgroup,shubnikov,spgroupma
!ENDDEBUG

!assigns the generators to each space group
 select case (spgroup)
!  ORTHORHOMBIC space groups
 case (16,21,22,23)        !P222, C222, F222, I222
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   if(shubnikov==3)symafm(3:4)=-1
 case (17,20)                !P2221, C2221
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmm(:,:)
   symrel(:,:,4) = genmpm(:,:)
   tnons(:,4)=(/0.d0,0.d0,0.5d0/)
   if(shubnikov==3)then
     symafm(4)=-1
     if(spgroupma==9) symafm(3)=-1
     if(spgroupma==10)symafm(2)=-1
     if(spgroupma==33)symafm(3:4)=-1
     if(spgroupma==34)symafm(2)=-1
   end if
 case (18)                !P21212
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
   if(shubnikov==3)then
     symafm(3)=-1
     if(spgroupma==18)symafm(4)=-1
     if(spgroupma==19)symafm(2)=-1
   end if
 case (19,24)                !P212121, I212121
   tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
   if(shubnikov==3)symafm(3:4)=-1
 case (25,35,38,42,44)        !Pmm2, Cmm2, Amm2, Fmm2, Imm2
   symrel(:,:,3) = genpmp(:,:)
   symrel(:,:,4) = genmpp(:,:)
 case (26,36)                !Pmc21, Cmc21
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genmpp(:,:)
 case (27,37)                !Pcc2, Ccc2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.d0,0.d0,0.5d0/)
 case (28,40,46)        !Pma2, Ama2, Ima2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.d0,0.d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.5d0,0.d0,0.d0/)
 case (29)                 !Pca21
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.d0,0.d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
 case (30)                !Pnc2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
 case (31)                !Pmn21
   tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,4) = genmpp(:,:)
 case (32,41,45)        !Pba2, Aba2, Iba2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
 case (33)                !Pna21
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
 case (34)                !Pnn2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
 case (39)                !Abm2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.d0,0.5d0,0.d0/)
 case (43)                !Fdd2
   symrel(:,:,3) = genpmp(:,:)
   tnons(:,3)=(/0.25d0,0.25d0,0.25d0/)
   symrel(:,:,4) = genmpp(:,:)
   tnons(:,4)=(/0.25d0,0.25d0,0.25d0/)
 case (47,65,69,71)        !Pmmm, Cmmm, Fmmm, Immm
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
 case (48)                !Pnnn
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,5) = genmmm(:,:)
   symrel(:,:,6) = genppm(:,:)
   symrel(:,:,7) = genpmp(:,:)
   symrel(:,:,8) = genmpp(:,:)
   if (spgorig==1) then
     tnons(:,5)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,7)=(/0.5d0,0.5d0,0.5d0/)
     tnons(:,8)=(/0.5d0,0.5d0,0.5d0/)
   else
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,7)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,8)=(/0.d0,0.5d0,0.5d0/)
   end if
   if(shubnikov==3)then
     if(spgroupma==259)then
       symafm(2:3)=-1 ; symafm(5)=-1 ; symafm(8)=-1
     else if(spgroupma==260)then
       symafm(3:4)=-1 ; symafm(7:8)=-1
     else if(spgroupma==261)then
       symafm(5:8)=-1
     end if
   end if
   nogen=0
 case (49,66)                !Pccm, Cccm
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
 case (50)                !Pban
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,5) = genmmm(:,:)
   symrel(:,:,6) = genppm(:,:)
   symrel(:,:,7) = genpmp(:,:)
   symrel(:,:,8) = genmpp(:,:)
   if (spgorig==1) then
     tnons(:,5)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,7)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,8)=(/0.5d0,0.5d0,0.d0/)
   else
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,3)=(/0.5d0,0.d0,0.d0/)
     tnons(:,7)=(/0.5d0,0.d0,0.d0/)
     tnons(:,4)=(/0.d0,0.5d0,0.d0/)
     tnons(:,8)=(/0.d0,0.5d0,0.d0/)
   end if
   if(shubnikov==3)then
     if(spgroupma==279)then
       symafm(2:3)=-1 ; symafm(5)=-1 ; symafm(8)=-1
     else if(spgroupma==280)then
       symafm(3:4)=-1 ; symafm(5:6)=-1
     else if(spgroupma==281)then
       symafm(3:4)=-1 ; symafm(7:8)=-1
     else if(spgroupma==282)then
       symafm(2:8:2)=-1
     else if(spgroupma==283)then
       symafm(5:8)=-1
     end if
   end if
   nogen=0
 case (51)                !Pmma
   tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.d0,0.d0/)
 case (52)                !Pnna
   tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.d0,0.5d0,0.5d0/)
 case (53)                 !Pmna
   tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
 case (54)                !Pcca
   tnons(:,2)=(/0.5d0,0.d0,0.d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
 case (55,72)                !Pbam, Ibam
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
 case (56)                !Pccn
   tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
 case (57)                !Pbcm
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.d0,0.5d0,0.d0/)
 case (58)                !Pnnm
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
 case (59)                !Pmmn
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,5) = genmmm(:,:)
   symrel(:,:,6) = genppm(:,:)
   symrel(:,:,7) = genpmp(:,:)
   symrel(:,:,8) = genmpp(:,:)
   if (spgorig==1) then
     tnons(:,3)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,5)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
   else
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,3)=(/0.d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.d0/)
     tnons(:,6)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,7)=(/0.d0,0.5d0,0.d0/)
     tnons(:,8)=(/0.5d0,0.d0,0.d0/)
   end if
   if(shubnikov==3)then
     if(spgroupma==407)then
       symafm(2:3)=-1 ; symafm(5)=-1 ; symafm(8)=-1
     else if(spgroupma==408)then
       symafm(3:4)=-1 ; symafm(5:6)=-1
     else if(spgroupma==409)then
       symafm(3:4)=-1 ; symafm(7:8)=-1
     else if(spgroupma==410)then
       symafm(2:8:2)=-1
     else if(spgroupma==411)then
       symafm(5:8)=-1
     end if
   end if
   nogen=0
 case (60)                !Pbcn
   tnons(:,2)=(/0.5d0,0.5d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
 case (61,73)                !Pbca, Ibca
   tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
 case (62)                !Pnma
   tnons(:,2)=(/0.5d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.d0/)
   symrel(:,:,4) = genpmm(:,:)
   tnons(:,4)=(/0.5d0,0.5d0,0.5d0/)
 case (63)                !Cmcm
   tnons(:,2)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
   if(shubnikov==3)then
     if(spgroupma==459 .or. spgroupma==463)symafm(2:3)=-1
     if(spgroupma==460 .or. spgroupma==464)symafm(4)=-1
     if(spgroupma==460 .or. spgroupma==464)symafm(2)=-1
     if(spgroupma==461 .or. spgroupma==462)symafm(3:4)=-1
   end if
 case (64)                !Cmca
   tnons(:,2)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.5d0/)
   symrel(:,:,4) = genpmm(:,:)
 case (67)                !Cmma
   tnons(:,2)=(/0.d0,0.5d0,0.d0/)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.d0/)
 case (68)                !Ccca
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,5) = genmmm(:,:)
   symrel(:,:,6) = genppm(:,:)
   symrel(:,:,7) = genpmp(:,:)
   symrel(:,:,8) = genmpp(:,:)
   if (spgorig==1) then
     tnons(:,2)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,4)=(/0.5d0,0.5d0,0.d0/)
     tnons(:,5)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,6)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,7)=(/0.d0,0.5d0,0.5d0/)
     tnons(:,8)=(/0.5d0,0.d0,0.5d0/)
   else
     tnons(:,2)=(/0.5d0,0.d0,0.d0/)
     tnons(:,3)=(/0.d0,0.d0,0.5d0/)
     tnons(:,4)=(/0.5d0,0.d0,0.5d0/)
     tnons(:,6)=(/0.5d0,0.d0,0.d0/)
     tnons(:,7)=(/0.d0,0.d0,0.5d0/)
     tnons(:,8)=(/0.5d0,0.d0,0.5d0/)
   end if
   if(shubnikov==3)then
     if(spgroupma==513)then
       symafm(2:3)=-1 ; symafm(5)=-1 ; symafm(8)=-1
     else if(spgroupma==514)then
       symafm(3:4)=-1 ; symafm(5:6)=-1
     else if(spgroupma==515)then
       symafm(3:4)=-1 ; symafm(7:8)=-1
     else if(spgroupma==516)then
       symafm(2:8:2)=-1
     else if(spgroupma==517)then
       symafm(5:8)=-1
     end if
   end if
   nogen=0
 case (70)                !Fddd
   symrel(:,:,3) = genmpm(:,:)
   symrel(:,:,4) = genpmm(:,:)
   symrel(:,:,5) = genmmm(:,:)
   symrel(:,:,6) = genppm(:,:)
   symrel(:,:,7) = genpmp(:,:)
   symrel(:,:,8) = genmpp(:,:)
   if (spgorig==1) then
     tnons(:,5)=(/0.25d0,0.25d0,0.25d0/)
     tnons(:,6)=(/0.25d0,0.25d0,0.25d0/)
     tnons(:,7)=(/0.25d0,0.25d0,0.25d0/)
     tnons(:,8)=(/0.25d0,0.25d0,0.25d0/)
   else
     tnons(:,2)=(/0.75d0,0.75d0,0.d0/)
     tnons(:,3)=(/0.75d0,0.d0,0.75d0/)
     tnons(:,4)=(/0.d0,0.75d0,0.75d0/)
!      JWZ DEBUGGING BEGIN
!      the following lines were present in 5.7 as of Dec 10 2008 but
!      gave wrong results in a spgroup 70 spgorig 2 case (Na2SO4)
!      tnons(:,5)=(/0.25d0,0.25d0,0.d0/) ! original code
!      tnons(:,6)=(/0.25d0,0.d0,0.25d0/) ! original code
!      tnons(:,7)=(/0.d0,0.25d0,0.25d0/) ! original code
!      here are the corrected values of tnons for this case
     tnons(:,5)=(/0.0d0,0.0d0,0.0d0/)
     tnons(:,6)=(/0.25d0,0.25d0,0.0d0/)
     tnons(:,7)=(/0.25d0,0.0d0,0.25d0/)
     tnons(:,8)=(/0.0d0,0.25d0,0.25d0/)
!      JWZ DEBUGGING END
   end if
   if(shubnikov==3)then
     if(spgroupma==529)then
       symafm(2:3)=-1 ; symafm(5)=-1 ; symafm(8)=-1
     else if(spgroupma==530)then
       symafm(3:4)=-1 ; symafm(7:8)=-1
     else if(spgroupma==531)then
       symafm(5:8)=-1
     end if
   end if
   nogen=0
 case (74)                !Imma
   tnons(:,2)=(/0.d0,0.5d0,0.d0/)
   symrel(:,:,3) = genmpm(:,:)
   tnons(:,3)=(/0.d0,0.5d0,0.d0/)
   symrel(:,:,4) = genpmm(:,:)
 end select

 if (shubnikov==3) then
   select case (spgroupma)
   case (59,68,80,89,101,113,125,137,146,158,167,174,182,189,197,&
&     205,213,221,226,231,237,243,270,292,296,308,312,324,328,340,344,&
&     358,370,380,384,420,424,444,448,460,464,472,476)
     symafm(2)=-1
     symafm(4)=-1
   case (69,90,102,114,126,147,175,190,198,206,214,244)
     symafm(2:3)=-1
   case (60,70,81,91,103,115,127,138,148,159,168,176,183,191,199,&
&     207,215,222,227,232,238,245,252,268,269,293,294,309,310,&
&     325,326,341,342,356,357,368,369,381,382,396,397,421,422,436,445,&
&     446,461,462,473,474,484,485,494,495,504,505,524,536,&
&     542,543,551,557,558)
     symafm(3:4)=-1
   case (251,267,291,295,307,311,323,327,339,343,355,367,379,383,&
&     395,398,419,423,435,443,447,459,463,&
&     471,475,483,486,493,496,503,506,523,535,541,544,550,556,559)
     symafm(2:3)=-1

   end select
 end if

 if (nogen>1)then
   call bldgrp(msym,nogen,nsym,symafm,symrel,tnons)
 end if

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

!DEBUG
!write(std_out,*)'symsgortho : end of symmetry assignement'
!ENDDEBUG

end subroutine symsgortho
!!***

!!****f* m_symsg/symsgtetra
!! NAME
!! symsgtetra
!!
!! FUNCTION
!! Yields all the TETRAGONAL symmetry operations starting from the space group symbol.
!! according to the International Tables of Crystallography, 1983.
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

subroutine symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,spgroupma,symafm,symrel,tnons)

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

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

!DEBUG
!write(std_out,*)' symsgtetra : exit'
!do isym=1,nsym
!write(std_out,'(i3,2x,9i3,3es13.3,i3)') isym,symrel(:,:,isym),tnons(:,isym),symafm(isym)
!end do
!ENDDEBUG

end subroutine symsgtetra
!!***

!!****f* m_symsg/bldgrp
!! NAME
!! bldgrp
!!
!! FUNCTION
!! Yields all the symmetry operations starting from the generators.
!! Applies all the generators onto themselves, and obtains all the other operations.
!! Iterates until it reaches nsym.
!!
!! INPUTS
!! msym = default number of symmetry operations
!! nsym = number of symmetry operations
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym) = 3D matrix containg symmetry operations
!! tnons(3,msym) = 2D matrix containing translations of the symmery operations
!!
!! OUTPUT
!!
!! symafm(msym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,msym) = 3D matrix containg symmetry operations
!! tnons(3,msym) = 2D matrix containing translations of the symmery operations
!!
!! SIDE EFFECTS
!! nogen = number of generators, number of operations to be applied onto themselves
!!
!! PARENTS
!!      symsgcube,symsghexa,symsgortho,symsgtetra
!!
!! CHILDREN
!!
!! SOURCE

subroutine bldgrp(msym,nogen,nsym,symafm,symrel,tnons)

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,nsym
 integer,intent(inout) :: nogen
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym)
 real(dp),intent(inout) :: tnons(3,msym)

!Local variables ------------------------------
!matrintoper(3,3) & matrinttransl(3) are intermediate arrays of the new
!      symmetry operations obtained, in order to check their uniqueness.
!flagop,flagtr = flags used during the checking of the similarity between
!      the obtained operation and the already existent ones
!ii,ijk,ijkl,jjj,kk = counters in the cycles
!scalars
 integer :: flagop,flagtr,ii,ijk,ijkl,jjj,kk,matrintsymafm,nogen_new
 real(dp) :: nastyzero
 character(len=500) :: message
!arrays
 integer :: bcksymafm(2*msym),bcksymrel(3,3,2*msym),matrintoper(3,3)
 real(dp) :: bcktnons(3,2*msym),matrinttransl(3)

! *************************************************************************

 nastyzero=0.1

!DEBUG
! write(std_out,*)' bldgrp : enter, builds the space group symmetry '
! write(std_out,*)' bldgrp : number of generators : ',nogen
! write(std_out,*)' bldgrp : nsym,msym=',nsym,msym
!ENDDEBUG

 if (nogen<1) then
   write(message, '(a,i4,a,a,a,a,a)' )&
&   'The number of generators nogen is ',nogen,&
&   'and it should be greater than one',ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: Contact ABINIT group '
   MSG_ERROR(message)
 end if

!Transfer the generators to bcksymrel
 do ii=1,nogen
   bcksymrel(:,:,ii)=symrel(:,:,ii)
   bcktnons(:,ii)=tnons(:,ii)
   bcksymafm(ii)=symafm(ii)
 end do

!DEBUG
 write(std_out,*)' Describe the different generators (index,symrel,tnons,symafm)'
 do ii=1,nogen
   write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
 end do
!ENDDEBUG

!Simply iterate until the group is complete
 do ijkl=1,nsym

!  DEBUG
!  write(std_out,*)' bldgrp : in loop, ijkl,nogen=',ijkl,nogen
!  ENDDEBUG

   nogen_new=nogen

   do jjj=2,nogen
     do kk=2,nogen

!      Computing block of the new symmetry operation according to:
!      !   $ { R1 | v1 }{ R2 | v2 } = { R1.R2 | v1+R1.v2 } $
       matrintoper(:,:) = matmul(bcksymrel(:,:,jjj),bcksymrel(:,:,kk))
       matrinttransl(:) = bcktnons(:,jjj)+matmul(bcksymrel(:,:,jjj),bcktnons(:,kk))
       matrintsymafm    = bcksymafm(jjj)*bcksymafm(kk)

!      Rescaling translation between 0 and 1
       do ii=1,3
         if (matrinttransl(ii)>=0.9) then
           do while (matrinttransl(ii)>=0.9)
             matrinttransl(ii)=matrinttransl(ii)-1.0
           end do
         end if
         if (matrinttransl(ii)<0.0) then
           do while (matrinttransl(ii)<0.0)
             matrinttransl(ii)=matrinttransl(ii)+1.0
           end do
         end if
         if ( abs(matrinttransl(ii))<nastyzero) matrinttransl(ii)=0.0
         if ( abs(matrinttransl(ii)-1.0)<nastyzero) matrinttransl(ii)=0.0
       end do

!      Cheking block to validate the new symmetry operation
       do ijk=1,nogen_new

         flagop=0 ; flagtr=0

!        Check for rotation similarity
         if(sum((matrintoper-bcksymrel(:,:,ijk))**2)==0)flagop=1

!        Check for translation similarity
         if(maxval((matrinttransl-bcktnons(:,ijk))**2)<nastyzero**2)flagtr=1

         if(flagop+flagtr==2)exit

       end do

!      Add the new determined symmetry if it is unique
       if (flagtr+flagop<2) then
         nogen_new=nogen_new+1
!        DEBUG
!         write(std_out,*)' added one more symmetry : nogen_new=',nogen_new
!        ENDDEBUG
         bcksymrel(:,:,nogen_new)=matrintoper(:,:)
         bcktnons(:,nogen_new)=matrinttransl(:)
         bcksymafm(nogen_new)=matrintsymafm
       end if

     end do
   end do

   nogen=nogen_new

   if(nogen==nsym)exit

 end do

!Transfer of the calculated symmetry to the routine output
 if (nogen==nsym) then
   symrel(:,:,1:nsym)=bcksymrel(:,:,1:nsym)
   tnons(:,1:nsym)=bcktnons(:,1:nsym)
   symafm(1:nsym)=bcksymafm(1:nsym)
 else
!  Problem with the generation of the symmetry operations
   write(message, '(a,i7,a,a,i7)' )&
&   'The symmetries obtained are  ',nogen,ch10,&
&   'and they should be ',nsym
   MSG_BUG(message)
 end if

!DEBUG
!write(std_out,*)' bldgrp : exit with  ',nogen,' operation symmetries'
!ENDDEBUG

end subroutine bldgrp
!!***

end module m_symsg
!!***

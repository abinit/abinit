!{\src2tex{textfont=tt}}
!!****f* ABINIT/symsgcube
!! NAME
!! symsgcube
!!
!! FUNCTION
!! Generate all the symmetry operations starting from the space group symbol
!! for the cubic groups
!! (according to the International Tables of Crystallography, 1983)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!! TODO
!!
!! PARENTS
!!      gensymspgr
!!
!! CHILDREN
!!      bldgrp,spgdata,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symsgcube'
 use interfaces_18_timing
 use interfaces_41_geometry, except_this_one => symsgcube
!End of the abilint section

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

 call spgdata(brvsb,intsb,intsbl,ptintsb,&
& ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

end subroutine symsgcube

!!***

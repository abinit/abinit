!{\src2tex{textfont=tt}}
!!****f* ABINIT/symsgortho
!! NAME
!! symsgortho
!!
!! FUNCTION
!! Yields all the ORTHORHOMBIC symmetry operations starting from the space group symbol.
!! It deals only with the orthorhombic groups
!! taken in the standard orientation
!! according to the International Tables of Crystallography, 1983.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2018 ABINIT group (RC,XG)
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symsgortho'
 use interfaces_41_geometry, except_this_one => symsgortho
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

 call spgdata(brvsb,intsb,intsbl,ptintsb,&
& ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

!DEBUG
!write(std_out,*)'symsgortho : end of symmetry assignement'
!ENDDEBUG

end subroutine symsgortho
!!***

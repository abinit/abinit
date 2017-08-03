!{\src2tex{textfont=tt}}
!!****f* ABINIT/symsghexa
!! NAME
!! symsghexa
!!
!! FUNCTION
!! Yields all the TRIGONAL & HEXAGONAL symmetry operations starting from the space group symbol.
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symsghexa'
 use interfaces_41_geometry, except_this_one => symsghexa
!End of the abilint section

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

   call spgdata(brvsb,intsb,intsbl,ptintsb,&
&   ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

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

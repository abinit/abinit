!{\src2tex{textfont=tt}}
!!****f* ABINIT/symsgmono
!! NAME
!! symsgmono
!!
!! FUNCTION
!! Yields all the MONOCLINIC symmetry operations starting from the space group symbol.
!! according to the International Tables of Crystallography, 1983.
!! It solves also the problem of the axes orientation
!! according to the spgaxor
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)


 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symsgmono'
 use interfaces_41_geometry, except_this_one => symsgmono
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

 call spgdata(brvsb,intsb,intsbl,ptintsb,&
& ptschsb,schsb,spgaxor,spgroup,sporder,spgorig)

end subroutine symsgmono
!!***

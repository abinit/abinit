!{\src2tex{textfont=tt}}
!!****f* ABINIT/gensymspgr
!! NAME
!! gensymspgr
!!
!! FUNCTION
!! Give all the symmetry operations starting from the space group symbol.
!! Suppose we are working in a conventional cell
!! If brvltt 0 or positive, the pure translations of the
!! Bravais lattice are included as generator of the group.
!! In brvltt=-1, no pure translation is present, and the
!! cell should be changed from conventional to primitive,
!! outside of this routine.
!! Treat also Shubnikov type III space groups.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (RC,XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! brvltt = input variable giving Bravais lattice
!! msym = default number of symmetry operations
!! shubnikov= magnetic type of the space group to be generated
!! spgaxor = orientation of the cell axis (might be needed)
!! spgorig = second choice of origin for certain groups
!! (might be needed if nsym==0)
!! spgroup = number of space group
!! spgroupma= number of the magnetic space group
!!
!! OUTPUT
!! nsym = number of symmetry operations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space in terms
!!                  of primitive translations
!! tnons(3,nsym)=nonsymmorphic translations for symmetry operations
!!
!! SIDE EFFECTS
!! brvltt = input variable giving Bravais lattice
!!
!! NOTES
!!
!! PARENTS
!!      ingeo
!!
!! CHILDREN
!!      chkgrp,print_symmetries,symsgcube,symsghexa,symsgmono,symsgortho
!!      symsgtetra
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine gensymspgr(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,&
& spgroup,spgroupma,symafm,symrel,tnons)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_crystal, only : print_symmetries

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'gensymspgr'
 use interfaces_41_geometry, except_this_one => gensymspgr
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: msym,shubnikov,spgaxor,spgorig,spgroup,spgroupma
 integer,intent(inout) :: brvltt
 integer,intent(out) :: nsym
!arrays
 integer,intent(inout) :: symafm(msym),symrel(3,3,msym) !vz_i
 real(dp),intent(inout) :: tnons(3,msym) !vz_i

!Local variables ------------------------------
! intsym,inttn = intermediate real to swap the columns of the symmetry matrix
! bckbrvltt = backup to brvltt to compare the assigned and the input values
!real(dp) :: tsec(2)
!scalars
 integer :: bckbrvltt,ii,inversion,jj,kk,ierr
 integer :: intsym
 real(dp) :: inttn
 character(len=500) :: message

! *************************************************************************

!List of the input parameters
!DEBUG
!write(std_out,*)' gensymspgr : enter with:'
!write(std_out,*)' spgroup = ',spgroup
!write(std_out,*)' spgaxor = ',spgaxor
!write(std_out,*)' spgorig = ',spgorig
!write(std_out,*)' brvltt  = ',brvltt
!ENDDEBUG

!Assume that the value of spgroupma is consistent with the one of spgroup
!(this has been checked earlier)

!Tests for consistency first the space group number and then the orientation
!Checks the space group number
 if (.not.(spgroup>0 .and. spgroup<231) ) then
   write(message, '(a,i0,a,a,a,a)' )&
&   'spgroup must be between 1 to 230, but is ',spgroup,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify spgroup in the input file.'
   MSG_ERROR(message)
 end if

!Checks the orientation
 if (.not.(spgaxor>0 .and. spgaxor<10)) then
   write(message, '(a,i12,a,a,a,a)' )&
&   'spgaxor must be from 1 to 9, but is',spgaxor,ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify spgaxor in the input file.'
   MSG_ERROR(message)
 end if

!Checks the consistency between the origin and space group
 if (spgorig==1 .or. spgorig==2) then
 else
   write(message, '(a,i4,a,a,a,a,a,a)' )&
&   'spgorig is',spgorig,ch10,&
&   'while it should be 0 or 1',ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify spgorig in the input file.'
   MSG_ERROR(message)
 end if

 if (spgorig>1) then
   select case (spgroup)
   case (48,50,59,68,70,85,86,88,125,126,129,130,133,134,137,138,&
&     141,142,201,203,222,224,227,228)
   case default
     write(message, '(a,a,a,a,a)' )&
&     'spgroup does not accept several origin choices',ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action: modify spgorig in the input file.'
     MSG_ERROR(message)
   end select
 end if

!Checks for consistency between the orientation and space group
 if (spgaxor>1) then
   select case (spgroup)
   case (3:74,146,148,155,160,161,166,167)
   case default
     write(message, '(a,a,a,a,a)' )&
&     'spgroup does not accept several orientations',ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action: modify spgaxor or spgroup in the input file.'
     MSG_ERROR(message)
   end select
 end if

 if (brvltt<-1 .or. brvltt>7)then
   write(message, '(a,i4,a,a,a,a,a,a)' )&
&   'The input brvltt was ',brvltt,ch10,&
&   'and it should be an integer from -1 to 7',ch10,&
&   'This is not allowed.  ',ch10,&
&   'Action: modify brvltt  in the input file.'
   MSG_ERROR(message)
 end if

!Assign nsym for each group according first to the order of the group
!Note that this value might be modified later:
!first because of the product with the inversion,
!second because of the centering operations
 select case (spgroup)
 case (1,2)
   nsym=1
 case (3:9)
   nsym=2
 case (143:148)
   nsym=3
 case (10:42,44:47,49,51:58,60:67,69,71:84,87)
   nsym=4
 case (149:176)
   nsym=6
 case (48,50,59,68,70,85,86,88:121,123,124,127,128,131,132,135,136,139,140)
   nsym=8
 case (177:200,202,204:206)
   nsym=12
 case (43,122,125,126,129,130,133,134,137,138)
   nsym=16
 case (201,207:219,221,223,225,226,229,230)
   nsym=24
 case (141,142)
   nsym=32
 case (203,220,222,224)
   nsym=48
 case (227,228)
   nsym=192
 end select

!DEBUG
!write(std_out,*)'gensymspgr :  assigns nsym = ',nsym
!ENDDEBUG

!Makes a backup to the brvltt for further comparison with the assigned value
 bckbrvltt=brvltt
!Default brvltt
 brvltt=1

!call timab(47,1,tsec)

!Assigns the first part of the symmetry operations:
!Rotation axis and mirror planes with or without translations,
!and sometimes also inversion operations.
 select case (spgroup)
 case (1:2)
   symrel(:,:,1)=0
   symrel(1,1,1)=1 ; symrel(2,2,1)=1 ; symrel(3,3,1)=1
   tnons(:,1)=zero ; symafm(1)=1
 case (3:15)
   call symsgmono(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)
 case (16:74)
   call symsgortho(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)
 case (75:142)
   call symsgtetra(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)
 case (143:194)
   call symsghexa(brvltt,msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)
 case (195:230)
   call symsgcube(msym,nsym,shubnikov,spgaxor,spgorig,spgroup,&
&   spgroupma,symafm,symrel,tnons)
 end select

!call timab(47,2,tsec)

!Assign the inversion center (if necessary).
!Note that for monoclinic space groups, the inversion was already
!assigned in symsgmono.f. Some other inversions have also been assigned in the
!corresponding system routine
 inversion=0
 select case (spgroup)
 case (2,47,49,51:58,60:67,69,71:74,83,84,87,123,124,127,128,131,132,&
&   135,136,139,140,147,148,162:167,175,176,191:194,200,202,204:206,&
&   221,223,225,226,229,230)
   inversion=1
!    Treat the magnetic part
   if(shubnikov==3)then
     select case (spgroup)
     case(2) ! Triclinic
       inversion=-1
     case(47,49,51:58,60:67,69,71:74) ! Orthorhombic
       select case (spgroupma)
       case(251,253,259,261,267,268,271,279,280,283,291,292,293,297,307,&
&         308,309,313,323,324,325,329,339,340,341,345,355,356,359,367,368,371,&
&         379,380,381,385,395,396,399,407,408,411,419,420,421,425,435,437,443,&
&         444,445,449,459,460,461,465,471,472,473,477,483,484,487,493,494,497,&
&         503,504,507,513,514,517,523,525,529,531,535,537,541,542,545,550,552,556,557,560)
         inversion=-1
       end select
     case(83,84,87,123,124,127,128,131,132,135,136,139,140) ! Tetragonal
       select case (spgroupma)
       case(46,47,54,55,62,63,70,71,78,79,84,85,341,344,346,347,353,356,358,359,365,&
&         368,370,371,377,380,382,383,389,392,394,395,401,404,406,407,&
&         413,416,418,419,425,428,430,431,437,440,442,443,449,452,454,455,&
&         461,464,466,467,473,476,478,479,485,488,490,491,497,500,502,503,&
&         509,512,514,515,521,524,526,527,533,536,538,539,543,546,548,549)
         inversion=-1
       end select
     case(147,148,162:167,175,176,191:194) ! Hexagonal or rhombohedral
       select case (spgroupma)
       case(15,19,75,76,81,82,87,88,93,94,99,100,105,106,139,140,145,146,&
&         235,236,237,241,245,246,247,251,255,256,257,261,265,266,267,271)
         inversion=-1
       end select
     case(200,202,204:206,221,223,225,226,229,230) ! Cubic
       select case (spgroupma)
       case(16,20,24,28,32,35,39,94,96,100,102,106,108,112,114,118,120,124,126,&
&         130,132,136,138,142,144,147,149)
         inversion=-1
       end select
     end select
   end if
 end select

!DEBUG
!write(std_out,*)' gensymspgr : before inversion'
!write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!do ii=1,nsym
!write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!end do
!ENDDEBUG

 if(inversion/=0)then
   do ii=1,nsym        ! visit all the symmetries assigned before
     do jj=1,3        ! visit the 3x3 matrix corresponding to the symmetry i
       tnons(jj,nsym+ii)=-tnons(jj,ii)
       do kk=1,3
         symrel(jj,kk,nsym+ii)=-symrel(jj,kk,ii)
       end do
     end do
     symafm(nsym+ii)=inversion*symafm(ii)
   end do
   nsym=nsym*2
 end if

!DEBUG
!write(std_out,*)' gensymspgr : after inversion'
!write(std_out,*)' Describe the different symmetry operations (index,symrel,tnons,symafm)'
!do ii=1,nsym
!write(std_out,'(i3,2x,9i3,3es12.2,i3)')ii,symrel(:,:,ii),tnons(:,ii),symafm(ii)
!end do
!ENDDEBUG

!Assign the Bravais lattice to each space group to which it has not yet
!been assigned
 select case (spgroup)
 case (38:41)
   brvltt=5                ! A
 case (20,21,35:37,63:68)
   brvltt=4                ! C
 case (22,42,43,69,70,196,202,203,209,210,216,219,225:228)
   brvltt=3                ! F
 case (23,24,44:46,71:74,79,80,82,87,88,97,98,107:110,119:122,&
&   139:142,197,199,204,206,211,214,217,220,229,230)
   brvltt=2                ! I
 case (146,148,155,160,161,166,167)
   if (spgaxor==1) then
     brvltt=7
   end if
 end select

 if (bckbrvltt/=0 .and. bckbrvltt/=-1) then
   if (bckbrvltt/=brvltt) then
     write(message, '(a,i8,a,a,a,i8,a,a)' )&
&     'The assigned brvltt ',brvltt,' is not equal',ch10,&
&     'to the input value ',bckbrvltt,ch10,&
&     'Assume experienced user. Execution will continue.'
     MSG_WARNING(message)
   end if
 end if

!if(bckbrvltt>=0)then
!Complete the set of primitive symmetries by translations
!associated with brvltt.
 select case (brvltt)
!  Bravais lattice type : A ! translation associated: b/2+c/2
 case (5)
   do ii=1,nsym
     tnons(1,nsym+ii)=tnons(1,ii)
     tnons(2,nsym+ii)=tnons(2,ii)+0.5
     tnons(3,nsym+ii)=tnons(3,ii)+0.5
     symrel(:,:,nsym+ii)=symrel(:,:,ii)
     symafm(nsym+ii)=symafm(ii)
   end do
   nsym=nsym*2

!    Bravais lattice type : B ! translation associated: a/2+c/2
 case (6)
   do ii=1,nsym
     tnons(1,nsym+ii)=tnons(1,ii)+0.5
     tnons(2,nsym+ii)=tnons(2,ii)
     tnons(3,nsym+ii)=tnons(3,ii)+0.5
     symrel(:,:,nsym+ii)=symrel(:,:,ii)
     symafm(nsym+ii)=symafm(ii)
   end do
   nsym=nsym*2

!    Bravais lattice type : C ! translation associated: a/2+b/2
 case (4)
   do ii=1,nsym
     tnons(1,nsym+ii)=tnons(1,ii)+0.5
     tnons(2,nsym+ii)=tnons(2,ii)+0.5
     tnons(3,nsym+ii)=tnons(3,ii)
     symrel(:,:,nsym+ii)=symrel(:,:,ii)
     symafm(nsym+ii)=symafm(ii)
   end do
   nsym=nsym*2

!    Bravais lattice type : F ! translations associated: a/2+b/2,b/2+c/2,b/2+c/2
 case (3)
!    For space groups containing d elements, all the symmetry operations
!    have already been obtained
   if(spgroup/=43 .and. spgroup/=203 .and. spgroup/=227 .and. &
&   spgroup/=228)then
     do ii=1,nsym
!        First translation: a/2+b/2
       tnons(1,nsym+ii)=tnons(1,ii)+0.5
       tnons(2,nsym+ii)=tnons(2,ii)+0.5
       tnons(3,nsym+ii)=tnons(3,ii)
       symrel(:,:,nsym+ii)=symrel(:,:,ii)
       symafm(nsym+ii)=symafm(ii)
     end do
!      Second translation: b/2+c/2
     do ii=1,nsym
       tnons(1,2*nsym+ii)=tnons(1,ii)
       tnons(2,2*nsym+ii)=tnons(2,ii)+0.5
       tnons(3,2*nsym+ii)=tnons(3,ii)+0.5
       symrel(:,:,2*nsym+ii)=symrel(:,:,ii)
       symafm(2*nsym+ii)=symafm(ii)
     end do
!      Third translation: a/2+c/2
     do ii=1,nsym
       tnons(1,3*nsym+ii)=tnons(1,ii)+0.5
       tnons(2,3*nsym+ii)=tnons(2,ii)
       tnons(3,3*nsym+ii)=tnons(3,ii)+0.5
       symrel(:,:,3*nsym+ii)=symrel(:,:,ii)
       symafm(3*nsym+ii)=symafm(ii)
     end do
     nsym=nsym*4
   end if

!    Bravais lattice type: I ! translation associated: a/2+b/2+c/2
 case (2)
!    For space groups containing d elements, all the symmetry operations
!    have already been obtained
   if(spgroup/=122 .and. spgroup/=141 .and. spgroup/=142 .and. &
&   spgroup/=220 )then
     do ii=1,nsym        ! visit all the symmetries assigned before
       tnons(:,nsym+ii)=tnons(:,ii)+0.5
       symrel(:,:,nsym+ii)=symrel(:,:,ii)
       symafm(nsym+ii)=symafm(ii)
     end do
     nsym=nsym*2
   end if

!    Bravais lattice type: R
!    translations for hexagonal axes ONLY: (2/3,1/3,1/3) & (1/3,2/3,2/3)
!    first translation (2/3,1/3,1/3)
 case (7)
   do ii=1,nsym
     tnons(1,nsym+ii)=tnons(1,ii)+two_thirds
     tnons(2,nsym+ii)=tnons(2,ii)+third
     tnons(3,nsym+ii)=tnons(3,ii)+third
     symrel(:,:,nsym+ii)=symrel(:,:,ii)
     symafm(nsym+ii)=symafm(ii)
   end do
!    Second translation (1/3,2/3,2/3)
   do ii=1,nsym
     tnons(1,2*nsym+ii)=tnons(1,ii)+third
     tnons(2,2*nsym+ii)=tnons(2,ii)+two_thirds
     tnons(3,2*nsym+ii)=tnons(3,ii)+two_thirds
     symrel(:,:,2*nsym+ii)=symrel(:,:,ii)
     symafm(2*nsym+ii)=symafm(ii)
   end do
   nsym=nsym*3

 end select

!end if

!Translate tnons in the ]-0.5,0.5] interval
 tnons(:,1:nsym)=tnons(:,1:nsym)-nint(tnons(:,1:nsym)-1.0d-8)

!Orientations for the orthorhombic space groups
!WARNING : XG 000620 : I am not sure that this coding is correct !!
 if (spgroup>15 .and. spgroup <75) then
   select case (spgaxor)
   case (1)             ! abc
     write(std_out,*)' the choosen orientation corresponds to: abc; the proper one'
   case (2)             ! cab
     do ii=1,nsym
       intsym=symrel(1,1,ii)
       symrel(1,1,ii)=symrel(2,2,ii)
       symrel(2,2,ii)=intsym
       inttn=tnons(1,ii)
       tnons(1,ii)=tnons(2,ii)
       tnons(2,ii)=inttn
     end do
     do ii=1,nsym
       intsym=symrel(1,1,ii)
       symrel(1,1,ii)=symrel(3,3,ii)
       symrel(3,3,ii)=intsym
       inttn=tnons(1,ii)
       tnons(1,ii)=tnons(3,ii)
       tnons(3,ii)=inttn
     end do
     write(std_out,*)' the choosen orientation corresponds to:  cab'
   case (3)             ! bca
     do ii=1,nsym
       intsym=symrel(1,1,ii)
       symrel(1,1,ii)=symrel(2,2,ii)
       symrel(2,2,ii)=intsym
       inttn=tnons(1,ii)
       tnons(1,ii)=tnons(2,ii)
       tnons(2,ii)=inttn
     end do
     do ii=1,nsym
       intsym=symrel(2,2,ii)
       symrel(2,2,ii)=symrel(3,3,ii)
       symrel(3,3,ii)=intsym
       inttn=tnons(2,ii)
       tnons(2,ii)=tnons(3,ii)
       tnons(3,ii)=inttn
     end do
     write(std_out,*)' the choosen orientation corresponds to:  bca'
   case (4)             ! acb
     do ii=1,nsym
       intsym=symrel(2,2,ii)
       symrel(2,2,ii)=symrel(3,3,ii)
       symrel(3,3,ii)=intsym
       inttn=tnons(1,ii)
       tnons(2,ii)=tnons(3,ii)
       tnons(3,ii)=inttn
     end do
     write(std_out,*)' the choosen orientation corresponds to:  acb'
   case (5)             ! bac
     do ii=1,nsym
       intsym=symrel(1,1,ii)
       symrel(1,1,ii)=symrel(2,2,ii)
       symrel(2,2,ii)=intsym
       inttn=tnons(1,ii)
       tnons(1,ii)=tnons(2,ii)
       tnons(2,ii)=inttn
     end do
     write(std_out,*)' the choosen orientation corresponds to:  bac'
   case (6)             ! cba
     do ii=1,nsym
       intsym=symrel(1,1,ii)
       symrel(1,1,ii)=symrel(3,3,ii)
       symrel(3,3,ii)=intsym
       inttn=tnons(1,ii)
       tnons(1,ii)=tnons(3,ii)
       tnons(3,ii)=inttn
     end do
     write(std_out,*)' the choosen orientation corresponds to:  cba'
   end select
 end if

!DEBUG
!write(std_out,*)' gensymspgr  : out of the Bravais lattice, nsym is',nsym
!ENDDEBUG

 call chkgrp(nsym,symafm,symrel,ierr)
 if (ierr/=0) then
   call print_symmetries(nsym,symrel,tnons,symafm)
 end if

 ABI_CHECK(ierr==0,"Error in group closure")

end subroutine gensymspgr
!!***

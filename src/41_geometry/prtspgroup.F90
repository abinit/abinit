!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtspgroup
!! NAME
!! prtspgroup
!!
!! FUNCTION
!! Print the space group (first, the dataset)
!!
!! COPYRIGHT
!! Copyright (C) 2002-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  bravais(11)=characteristics of Bravais lattice (see symlatt.f)
!!  genafm(3)=generator of magnetic translations, in case of
!!            Shubnikov type IV magnetic groups (if zero, the group is
!!            not a type IV magnetic group)
!!  iout=unit number of output file
!!  jdtset= actual number of the dataset to be read
!!  ptgroupma=magnetic point group, in case of
!!            Shubnikov type III magnetic groups (if zero, the group is
!!            not a type III magnetic group)
!!  spgroup=space group number
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      memory_eval
!!
!! CHILDREN
!!      ptgmadata,spgdata,wrtout,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtspgroup(bravais,genafm,iout,jdtset,ptgroupma,spgroup)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtspgroup'
 use interfaces_14_hidewrite
 use interfaces_41_geometry, except_this_one => prtspgroup
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iout,jdtset,ptgroupma,spgroup
!arrays
 integer,intent(in) :: bravais(11)
 real(dp),intent(inout) :: genafm(3)

!Local variables -------------------------------
!scalars
 integer :: center,iholohedry,ii,shubnikov,spgaxor,spgorig,sporder,sumgen
 character(len=1) :: brvsb
 character(len=10) :: ptgrpmasb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl
 character(len=500) :: message
 character(len=80) :: bravais_name
!arrays
 integer :: genafmint(3)
 real(dp) :: genafmconv(3),rprimdconv(3,3)

!*************************************************************************

!DEBUG
!write(std_out,*)' prtspgroup : enter '
!write(std_out,*)' ptgroupma=',ptgroupma
!write(std_out,*)' genafm(:)=',genafm(:)
!ENDDEBUG

 center=bravais(2)
 iholohedry=bravais(1)

!Determine the magnetic type
 shubnikov=1
 if(ptgroupma/=0)shubnikov=3
 if(sum(abs(genafm(:)))>tol6)then
   shubnikov=4
!  Produce genafm in conventional axes,
   rprimdconv(:,1)=bravais(3:5)
   rprimdconv(:,2)=bravais(6:8)
   rprimdconv(:,3)=bravais(9:11)
   if(center/=0)rprimdconv(:,:)=rprimdconv(:,:)*half
   call xred2xcart(1,rprimdconv,genafmconv,genafm)
!  Gives the associated translation, with components in the
!  interval ]-0.5,0.5] .
   genafmconv(:)=genafmconv(:)-nint(genafmconv(:)-tol6)
   do ii=1,3
     genafmint(ii)=-1
     if(abs(genafmconv(ii)-zero)<tol6)genafmint(ii)=0
     if(abs(genafmconv(ii)-half)<tol6)genafmint(ii)=1
   end do
   if(minval(genafmint(:))==-1)then
     write(message, '(3a,3es12.2,a)' )&
&     'The magnetic translation generator,',ch10,&
&     'genafmconv(:)=',genafmconv(:),&
&     'could not be identified.'
     MSG_BUG(message)
   end if
 end if

!Determine whether the space group can be printed
 if(iholohedry<=0)then
   if(jdtset/=0)then
     write(message,'(a,a,i5,a)')ch10,&
&     ' DATASET',jdtset,' : the unit cell is not primitive'
   else
     write(message,'(a,a)')ch10,&
&     ' Symmetries : the unit cell is not primitive'
   end if
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
 else if(spgroup==0)then
   if(jdtset/=0)then
     write(message,'(a,a,i5,a)')ch10,&
&     ' DATASET',jdtset,' : the space group has not been recognized'
   else
     write(message,'(a,a)')ch10,&
&     ' Symmetries : the space group has not been recognized'
   end if
   call wrtout(std_out,message,'COLL')
   call wrtout(iout,message,'COLL')
 else

!  ------------------------------------------------------------------
!  The space group can be printed

!  Determine the Bravais lattice

   bravais_name=' (the Bravais lattice could not be identified)'

   if(iholohedry==7)then ! Cubic

     if(center==0) then
       if(shubnikov/=4)bravais_name='cP (primitive cubic)'
       if(shubnikov==4)bravais_name='cP_I (primitive cubic, inner magnetic, #33)'
     else if(center==-1) then
       if(shubnikov/=4)bravais_name='cI (body-center cubic)' ! Only non-magnetic is possible
     else if(center==-3) then
       if(shubnikov/=4)bravais_name='cF (face-center cubic)'
       if(shubnikov==4)bravais_name='cF_s (face-center cubic, simple cubic magnetic, #35)'
     end if

   else if(iholohedry==4)then ! Tetragonal

     if(center==0) then
       if(shubnikov/=4)bravais_name='tP (primitive tetrag.)'
       if(shubnikov==4)then
         sumgen=sum(genafmint(:))
         if(sumgen==1)bravais_name='tP_c (primitive tetrag., c-magnetic, #23)'
         if(sumgen==2)bravais_name='tP_C (primitive tetrag., C-magnetic, #24)'
         if(sumgen==3)bravais_name='tP_I (primitive tetrag., centered magnetic, #25)'
       end if
     else if(center==-1)then
       if(shubnikov/=4)bravais_name='tI (body-center tetrag.)'
       if(shubnikov==4)bravais_name='tI_c (body-center tetrag., simple tetragonal magnetic, #27)'
     end if

   else if(iholohedry==3)then ! Orthorhombic

     if(center==0) then
       if(shubnikov/=4)bravais_name='oP (primitive ortho.)'
       if(shubnikov==4)then
         sumgen=sum(genafmint(:))
         if(sumgen==1)then
           if(genafmint(1)==1)bravais_name='oP_a (primitive ortho., a-magnetic, #11)'
           if(genafmint(2)==1)bravais_name='oP_b (primitive ortho., b-magnetic, #11)'
           if(genafmint(3)==1)bravais_name='oP_c (primitive ortho., c-magnetic, #11)'
         else if(sumgen==2)then
           if(genafmint(1)==0)bravais_name='oP_A (primitive ortho., A-magnetic, #12)'
           if(genafmint(2)==0)bravais_name='oP_B (primitive ortho., B-magnetic, #12)'
           if(genafmint(3)==0)bravais_name='oP_C (primitive ortho., C-magnetic, #12)'
         else if(sumgen==3)then
           bravais_name='oP_I (primitive ortho., centered magnetic, #13)'
         end if
       end if
     else if(center==-1)then
       if(shubnikov/=4)bravais_name='oI (body-center ortho.)'
       if(shubnikov==4)bravais_name='oI_c (body-center ortho., simple ortho. magn., #21)'
     else if(center==1 .or. center==2 .or. center==3)then
       if(shubnikov/=4) bravais_name='oC (1-face-center ortho.)'
       if(shubnikov==4)then
         sumgen=sum(genafmint(:))
         if(sumgen==1)then
!          One should work more to distinguish these magnetic groups
           bravais_name='oC_(a,b,c) (1-face-cent. ortho., 1-magn., #15 or 16)'
         else if(sumgen==2)then
           bravais_name='oC_A (1-face-centered ortho., 1-face-magnetic, #17)'
         else if(sumgen==3)then
           bravais_name='oC_c (C-face-centered ortho., c-magnetic, #15)'
         end if
       end if
     else if(center==-3)then
       if(shubnikov/=4)bravais_name='oF (face-center ortho.)'
       if(shubnikov==4)bravais_name='oF_s (face-center ortho., simple ortho. magnetic, #19)'
     end if

   else if(iholohedry==6)then ! Hexagonal

     if(shubnikov/=4)bravais_name='hP (primitive hexag.)'
     if(shubnikov==4)bravais_name='hP_c (primitive hexag., c-magnetic, #29)'

   else if(iholohedry==5)then ! Rhombohedral

     if(shubnikov/=4)bravais_name='hR (rhombohedral)'
     if(shubnikov==4)bravais_name='hR_I (rhombohedral, centered magnetic, #31)'

   else if(iholohedry==2)then ! Monoclinic

     if(center==0)then
       if(shubnikov/=4)bravais_name='mP (primitive monocl.)'
       if(shubnikov==4)then
         sumgen=sum(genafmint(:))
         if(sumgen==1)then
           if(genafmint(1)==1)bravais_name='mP_a (primitive monocl., a-magnetic, #5)'
           if(genafmint(2)==1)bravais_name='mP_b (primitive monocl., b-magnetic, #4)'
           if(genafmint(3)==1)bravais_name='mP_c (primitive monocl., c-magnetic, #5)'
         else if(sumgen==2)then
           if(genafmint(1)==0)bravais_name='mP_A (primitive monocl., A-magnetic, #6)'
           if(genafmint(2)==0)bravais_name='mP_B (primitive monocl., B-magnetic, #6)'
           if(genafmint(3)==0)bravais_name='mP_C (primitive monocl., C-magnetic, #6)'
         end if
       end if
     else if(center==3)then
       if(shubnikov/=4)bravais_name='mC (1-face-center monocl.)'
       if(shubnikov==4)then
         if(genafmint(3)==1)bravais_name='mC_c (C-face-center monocl., c-magnetic, #8)'
         if(genafmint(3)/=1)bravais_name='mC_a (C-face-center monocl., a-magnetic, #9)'
       end if
     else if(center==-3)then
       if(shubnikov/=4)bravais_name='(reduction of face-center)'
     end if

   else if(iholohedry==1)then ! Triclinic

     if(shubnikov/=4)bravais_name='aP (primitive triclinic)'
     if(shubnikov==4)bravais_name='aP_s (primitive triclinic, simple magnetic, #2)'

   end if

!  Determine the symbol of the Fedorov space group
   spgaxor=1 ; spgorig=1
   call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&
&   schsb,spgaxor,spgroup,sporder,spgorig)

!  Prepare print of the dataset, symmetry point group, Bravais lattice
   if(shubnikov==1)then

     if(jdtset/=0)then
       write(message,'(a,a,i5,a,a,a,a,i3,a,a,a)' )ch10,&
&       ' DATASET',jdtset,' : space group ',trim(brvsb),trim(intsb),' (#',spgroup,')',&
&       '; Bravais ',trim(bravais_name)
     else
       write(message,'(a,a,a,a,a,i3,a,a,a)' )ch10,&
&       ' Symmetries : space group ',trim(brvsb),trim(intsb),' (#',spgroup,')',&
&       '; Bravais ',trim(bravais_name)
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

   else if(shubnikov==3)then

     if(jdtset/=0)then
       write(message,'(a,a,i5,a)' )ch10,&
&       ' DATASET',jdtset,' : magnetic group, Shubnikov type III '
     else
       write(message,'(2a)' )ch10,&
&       ' Magnetic group, Shubnikov type III '
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

     write(message,'(a,a,a,a,i3,a,a,a)' )&
&     ' Fedorov space group ',trim(brvsb),trim(intsb),' (#',spgroup,')',&
&     '; Bravais ',trim(bravais_name)
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

     call ptgmadata(ptgroupma,ptgrpmasb)

     write(message,'(3a,i3,a)' )&
&     ' Magnetic point group ',trim(ptgrpmasb),' (#',ptgroupma,')'
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

   else if(shubnikov==4)then

     if(jdtset/=0)then
       write(message,'(a,a,i5,a)' )ch10,&
&       ' DATASET',jdtset,' : magnetic group, Shubnikov type IV '
     else
       write(message,'(2a)' )ch10,&
&       ' Magnetic group, Shubnikov type IV '
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

     write(message,'(a,a,a,a,i3,a)' )&
&     ' Fedorov space group ',trim(brvsb),trim(intsb),' (#',spgroup,')'
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

     write(message,'(2a)' )&
&     ' Magnetic Bravais lattice ',trim(bravais_name)
     call wrtout(std_out,message,'COLL')
     call wrtout(iout,message,'COLL')

   end if
 end if

end subroutine prtspgroup
!!***

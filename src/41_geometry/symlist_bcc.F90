!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_bcc
!! NAME
!! symlist_bcc
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! BCC case.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2018 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! additional_info=information that is needed beyond n_axes, in order
!!  to discriminate between specific space groups
!! nsym=actual number of symmetries
!! n_axes(31)=array containing the number of all the possible symmetry operations
!!
!! OUTPUT
!! spgroup=space group number ; returns 0 if not found
!!
!! NOTES
!!
!! The list of symmetry operations is for the conventional cell
!!
!! TODO
!! For the time being there are several groups where uncertainties still exist
!! This will be solved in the very next ABINIT version
!!
!! PARENTS
!!      symspgr
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symlist_bcc(additional_info,nsym,n_axes,spgroup)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symlist_bcc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: additional_info,nsym
 integer,intent(out) :: spgroup
!arrays
 integer,intent(in) :: n_axes(31)

!Local variables-------------------------------
!character(len=500) :: message
!scalars
 integer :: ii
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(std_out,*) ' symlist_bcc : enter '
!write(std_out,*) ' nsym = ', nsym
!write(std_out,'(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(std_out,'(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(std_out,'(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 select case(nsym)

 case(4)
!    XG021015, from LSi.
!    The coding of this case is different because of a compiler problem on
!    SUN. Seems that an internal limit is reached. Do not know why.
   do ii=1,3
     if(ii==1)then
       n_axest(:)=0 ; n_axest(7)=1 ; n_axest(8)=1 ; n_axest(9)=1 ; n_axest(20)=1
       spgroup=5
     else if(ii==2)then
       n_axest(:)=0 ; n_axest(7)=1 ; n_axest(8)=1 ; n_axest(15)=1 ; n_axest(18)=1
       spgroup=8
     else if(ii==3)then
       n_axest(:)=0 ; n_axest(7)=1 ; n_axest(8)=1 ; n_axest(16)=1 ; n_axest(18)=1
       spgroup=9
     end if
     if(sum((n_axes-n_axest)**2)==0)exit
     spgroup=0
   end do

 case(8)

   n_axest=(/0,0,0,0,2,0,1,1,1,0,  0,0,0,0,0,2,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=15
   n_axest=(/0,0,0,0,0,0,1,1,3,0,  0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) then
     if(additional_info==1) spgroup=23
     if(additional_info==2) spgroup=24
   end if
   n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,2,0,0,2,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=44
   n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,0,4,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=45
   n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,1,2,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=46
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,2,0,0,0,0,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=79
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=80

   n_axest=(/0,4,0,0,0,0,1,1,2,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=82

 case(16)

   n_axest=(/0,0,0,0,2,0,1,1,3,0,  0,0,0,0,3,0,0,3,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=71
   n_axest=(/0,0,0,0,2,0,1,1,3,0,  0,0,0,0,1,4,0,1,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=72
   n_axest=(/0,0,0,0,2,0,1,1,3,0,  0,0,0,0,0,6,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=73
   n_axest=(/0,0,0,0,2,0,1,1,3,0,  0,0,0,0,2,2,0,2,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=74
   n_axest=(/0,0,0,0,0,0,1,1,6,0,  0,2,0,0,0,0,0,0,0,0,  4,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=97
   n_axest=(/0,0,0,0,0,0,1,1,6,0,  0,0,0,0,0,0,0,0,0,0,  4,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=98
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,2,0,0,8,0,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=107
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,2,0,0,4,4,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=108
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,4,0,4,0,0,0,  0,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=109
   n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,4,4,0,0,0,  0,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=110

   n_axest=(/0,4,0,0,2,0,1,1,2,0,  0,2,0,0,2,0,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=87
   n_axest=(/0,4,0,0,2,0,1,1,2,0,  0,0,0,0,0,2,0,0,0,0,  0,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=88
   n_axest=(/0,4,0,0,0,0,1,1,2,0,  0,0,0,0,4,0,0,0,0,0,  4,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=119
   n_axest=(/0,4,0,0,0,0,1,1,2,0,  0,0,0,0,0,4,0,0,0,0,  4,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=120
   n_axest=(/0,4,0,0,0,0,1,1,6,0,  0,0,0,0,4,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=121
   n_axest=(/0,4,0,0,0,0,1,1,6,0,  0,0,0,0,0,0,4,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=122

 case(24)

   n_axest=(/0,0,0,0,0,0,1,1,3,16, 0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) then
     if(additional_info==1) spgroup=197
     if(additional_info==2) spgroup=199
   end if

 case(32)

   n_axest=(/0,4,0,0,2,0,1,1,6,0,  0,2,0,0,10,0,0,0,0,0, 4,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=139
   n_axest=(/0,4,0,0,2,0,1,1,6,0,  0,2,0,0,6,4,0,0,0,0,  4,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=140
   n_axest=(/0,4,0,0,2,0,1,1,6,0,  0,0,0,0,4,2,4,0,0,0,  4,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=141
   n_axest=(/0,4,0,0,2,0,1,1,6,0,  0,0,0,0,0,6,4,0,0,0,  4,0,0,4,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=142

 case(48)

   n_axest=(/0,0,16,0,2,0,1,1,3,16,  0,0,0,0,3,0,0,3,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=204
   n_axest=(/0,0,16,0,2,0,1,1,3,16,  0,0,0,0,0,0,0,6,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=206
   n_axest=(/0,0,0,0,0,0,1,1,12,16,  0,6,0,0,0,0,0,0,0,6,  0,0,0,0,6,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=211
   n_axest=(/0,0,0,0,0,0,1,1,12,16,  0,0,0,0,0,0,0,0,0,6,  0,0,0,12,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=214
   n_axest=(/0,12,0,0,0,0,1,1,3,16,  0,0,0,0,6,0,0,6,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=217
   n_axest=(/0,12,0,0,0,0,1,1,3,16,  0,0,0,0,0,0,12,0,0,3, 0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=220

 case(96)

!    Note that the identification of the mirror planes is still ambiguous for cI
   n_axest=(/0,12,16,0,2,0,1,1,12,16,  0,6,0,0,9,0,0,9,0,6,  0,0,0,0,6,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=229
   n_axest=(/0,12,16,0,2,0,1,1,12,16,  0,0,0,0,0,0,12,6,0,6, 0,0,0,12,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=230

 end select

end subroutine symlist_bcc
!!***

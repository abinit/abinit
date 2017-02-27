!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_fcc
!! NAME
!! symlist_fcc
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! FCC case
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
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


subroutine symlist_fcc(nsym,n_axes,spgroup)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symlist_fcc'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: spgroup
!arrays
 integer,intent(in) :: n_axes(31)

!Local variables-------------------------------
!character(len=500) :: message
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(std_out,*) ' symlist_fcc : enter '
!write(std_out,*) ' nsym = ', nsym
!write(std_out,'(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(std_out,'(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(std_out,'(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 select case(nsym)

 case(16)

   n_axest=(/0,0,0,0,0,0,3,1,6,0,  0,0,0,0,0,0,0,0,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=22
   n_axest=(/0,0,0,0,0,0,3,1,2,0,  0,0,0,0,2,4,0,2,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=42
   n_axest=(/0,0,0,0,0,0,3,1,2,0,  0,0,0,0,0,0,8,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=43


 case(32)

   n_axest=(/0,0,0,0,4,0,3,1,6,0,  0,0,0,0,3,6,0,3,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=69
   n_axest=(/0,0,0,0,4,0,3,1,6,0,  0,0,0,0,0,0,12,0,0,6, 0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=70

 case(48)

   n_axest=(/0,0,0,0,0,0,3,1,6,32, 0,0,0,0,0,0,0,0,0,6, 0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=196

 case(96)

   n_axest=(/0,0,32,0,4,0,3,1,6,32,  0,0,0,0,3,0,0,9,0,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=202
   n_axest=(/0,0,32,0,4,0,3,1,6,32,  0,0,0,0,0,0,12,0,0,6,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=203
   n_axest=(/0,0,0,0,0,0,3,1,18,32,  0,12,0,0,0,0,0,0,0,18, 0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=209
   n_axest=(/0,0,0,0,0,0,3,1,18,32,  0,0,0,0,0,0,0,0,0,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=210
   n_axest=(/0,24,0,0,0,0,3,1,6,32,  0,0,0,0,24,0,0,0,0,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=216
   n_axest=(/0,24,0,0,0,0,3,1,6,32,  0,0,0,0,0,0,0,24,0,6,   0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=219


 case(192)

   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,12,0,0,27,0,0,9,0,18, 0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=225
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,12,0,0,3,0,0,33,0,18,  0,0,0,0,12,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=226
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,24,0,12,0,0,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=227
   n_axest=(/0,24,32,0,4,0,3,1,18,32,  0,0,0,0,0,0,12,24,0,18,  0,0,0,24,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=228

 end select

end subroutine symlist_fcc
!!***

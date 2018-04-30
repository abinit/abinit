!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_prim
!! NAME
!! symlist_prim
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
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


subroutine symlist_prim(additional_info,nsym,n_axes,spgroup)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symlist_prim'
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
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(std_out,*) ' symlist_prim : enter '
!write(std_out,*) ' nsym = ', nsym
!write(std_out,'(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(std_out,'(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(std_out,'(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 select case(nsym)

 case(1)

   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=1

 case(2)

   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=2
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=3
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=4
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,1,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=6
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,1,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=7
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=7

 case(3)

   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=143
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,2,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=144
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,0,2,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=145

 case(4)

   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,1,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=10
   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,1,0,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=11
   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,0,1,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=13
   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=13
   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,0,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=14
   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,0,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=14

   n_axest=(/0,0,0,0,0,0,0,1,3,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=16
   n_axest=(/0,0,0,0,0,0,0,1,2,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=17
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=18
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=19
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)

   if(sum((n_axes-n_axest)**2)==0) spgroup=25
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,1,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=26
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,2,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=29
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) then
     if(additional_info==1) spgroup=27
     if(additional_info==2) spgroup=32
   end if
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,1,1,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=28
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,1,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=30
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,1,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=31
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=33
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,2,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=34

   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=75
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,2,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=76
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=77
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,2,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=78


   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=81

 case(6)

   n_axest=(/0,0,2,0,1,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=147

   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=149
   n_axest=(/0,0,0,3,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=150
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  3,2,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=151
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  3,0,2,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=153

   n_axest=(/0,0,0,3,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,2,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=152
   n_axest=(/0,0,0,3,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,  0,0,2,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=154

   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=156
   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,3,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=157
   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=158
   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,3,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=159

   n_axest=(/0,0,0,0,0,0,0,1,1,2,  0,0,0,2,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=168
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  0,2,0,0,0,0,2,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=169
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  0,0,2,0,0,0,0,0,0,0,2/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=170

   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  0,0,2,0,0,0,0,2,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=171
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  0,2,0,0,0,0,0,0,0,2,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=172

   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=173
   n_axest=(/2,0,0,0,0,0,0,1,0,2,  0,0,0,0,1,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=174

 case(8)

   n_axest=(/0,0,0,0,1,0,0,1,3,0,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=47
   n_axest=(/0,0,0,0,1,0,0,1,3,0,  0,0,0,0,0,0,0,3,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=48
   n_axest=(/0,0,0,0,1,0,0,1,3,0,  0,0,0,0,1,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=49
   n_axest=(/0,0,0,0,1,0,0,1,3,0,  0,0,0,0,0,2,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=50

   n_axest=(/0,0,0,0,1,0,0,1,2,0,  0,0,0,0,2,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=51
   n_axest=(/0,0,0,0,1,0,0,1,2,0,  0,0,0,0,0,1,0,2,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=52
   n_axest=(/0,0,0,0,1,0,0,1,2,0,  0,0,0,0,1,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=53
   n_axest=(/0,0,0,0,1,0,0,1,2,0,  0,0,0,0,0,3,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=54

   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,1,2,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0)  then
     if(additional_info==1) spgroup=55
     if(additional_info==2) spgroup=57
   end if
   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,0,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0)  then
     if(additional_info==1) spgroup=56
     if(additional_info==2) spgroup=60
   end if
   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,1,0,0,2,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=58
   n_axest=(/0,0,0,0,1,0,0,1,1,0,  0,0,0,0,2,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=59

   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,0,3,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=61
   n_axest=(/0,0,0,0,1,0,0,1,0,0,  0,0,0,0,1,1,0,1,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=62

   n_axest=(/0,0,0,0,0,0,0,1,3,0,  0,2,0,0,0,0,0,0,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=89
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,0,0,0,0,0,2,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=90

   n_axest=(/0,0,0,0,0,0,0,1,2,0,  0,0,0,0,0,0,0,0,0,1,  2,0,0,2,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=91
   n_axest=(/0,0,0,0,0,0,0,1,2,0,  0,0,0,0,0,0,0,0,0,1,  2,0,0,0,0,2,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=95

   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,3,  2,0,0,2,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=92
   n_axest=(/0,0,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,3,  2,0,0,0,0,2,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=96

   n_axest=(/0,0,0,0,0,0,0,1,3,0,  0,0,0,0,0,0,0,0,0,0,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=93
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,2,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=94

   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,4,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=99
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,2,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=100
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,2,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=101
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,0,2,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=102
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,0,0,2,0,2,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=103
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,2,0,0,0,0,0,2,2,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=104
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,0,0,2,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=105
   n_axest=(/0,0,0,0,0,0,0,1,1,0,  0,0,0,0,0,2,0,0,2,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=106

   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,1,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=83
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,1,0,0,0,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=84
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=85
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,0,0,0,1,0,0,  0,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=86

   n_axest=(/0,2,0,0,0,0,0,1,3,0,  0,0,0,0,2,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=111
   n_axest=(/0,2,0,0,0,0,0,1,3,0,  0,0,0,0,0,0,0,0,2,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=112
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=113
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,2,2,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=114
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,2,0,0,0,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=115
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,2,0,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=116
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,0,2,0,0,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=117
   n_axest=(/0,2,0,0,0,0,0,1,1,0,  0,0,0,0,0,0,0,2,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=118

 case(12)

   n_axest=(/0,0,2,0,1,0,0,1,0,2,  0,0,0,0,0,0,3,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=162
   n_axest=(/0,0,2,0,1,0,0,1,0,2,  0,0,0,0,0,0,0,3,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=163
   n_axest=(/0,0,2,3,1,0,0,1,0,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=164
   n_axest=(/0,0,2,3,1,0,0,1,0,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=165

   n_axest=(/2,0,2,0,1,0,0,1,1,2,  0,0,0,2,1,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=175
   n_axest=(/2,0,2,0,1,0,0,1,0,2,  0,0,0,0,1,0,0,0,0,1,  0,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=176

   n_axest=(/0,0,0,3,0,0,0,1,1,2,  0,0,0,2,0,0,0,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=177

   n_axest=(/0,0,0,3,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  3,2,0,0,0,0,2,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=178
   n_axest=(/0,0,0,3,0,0,0,1,0,0,  0,0,0,0,0,0,0,0,0,1,  3,0,2,0,0,0,0,0,0,0,2/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=179

   n_axest=(/0,0,0,3,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  3,0,2,0,0,0,0,2,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=180
   n_axest=(/0,0,0,3,0,0,0,1,1,0,  0,0,0,0,0,0,0,0,0,0,  3,2,0,0,0,0,0,0,0,2,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=181

   n_axest=(/0,0,0,3,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,1,  3,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=182

   n_axest=(/0,0,0,0,0,0,0,1,1,2,  0,0,0,2,3,0,3,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=183
   n_axest=(/0,0,0,0,0,0,0,1,1,2,  0,0,0,2,0,3,0,3,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=184
   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,3,3,0,0,1,  0,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=185
   n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,3,0,0,3,0,1,  0,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=186

   n_axest=(/2,0,0,0,0,0,0,1,0,2,  0,0,0,0,4,0,0,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=187
   n_axest=(/2,0,0,0,0,0,0,1,0,2,  0,0,0,0,1,3,0,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=188
   n_axest=(/2,0,0,3,0,0,0,1,0,2,  0,0,0,0,1,0,3,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=189
   n_axest=(/2,0,0,3,0,0,0,1,0,2,  0,0,0,0,1,0,0,3,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=190

   n_axest=(/0,0,0,0,0,0,0,1,3,8,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=195
   n_axest=(/0,0,0,0,0,0,0,1,0,8,  0,0,0,0,0,0,0,0,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=198

 case(16)

   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,2,0,0,5,0,0,0,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=123
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,2,0,0,1,0,2,0,2,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=124
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,2,0,0,2,2,0,1,0,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=125
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,2,0,0,0,0,0,3,2,0,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=126

   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,3,2,0,0,0,2,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=127
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,1,0,0,2,2,2,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=128
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,4,0,0,1,0,2,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=129
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,2,0,0,0,0,2,1,2,2,  2,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=130

   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,0,0,0,3,0,0,0,2,0,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=131
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,0,0,0,3,0,2,0,0,0,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=132
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,0,0,0,0,2,0,1,2,0,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=133
   n_axest=(/0,2,0,0,1,0,0,1,3,0,  0,0,0,0,2,0,0,3,0,0,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=134

   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,1,2,0,0,2,2,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=135
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,3,0,0,2,0,2,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=136
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,2,0,0,1,2,2,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=137
   n_axest=(/0,2,0,0,1,0,0,1,1,0,  0,0,0,0,2,0,2,1,0,2,  2,0,0,0,2,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=138

 case(24)

   n_axest=(/2,0,2,3,1,0,0,1,1,2,  0,0,0,2,4,0,3,0,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=191
   n_axest=(/2,0,2,3,1,0,0,1,1,2,  0,0,0,2,1,3,0,3,0,0,  3,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=192
   n_axest=(/2,0,2,3,1,0,0,1,0,2,  0,0,0,0,1,3,3,0,0,1,  3,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=193
   n_axest=(/2,0,2,3,1,0,0,1,0,2,  0,0,0,0,4,0,0,3,0,1,  3,0,0,0,0,0,0,0,2,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=194


   n_axest=(/0,0,8,0,1,0,0,1,3,8,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=200
   n_axest=(/0,0,8,0,1,0,0,1,3,8,  0,0,0,0,0,0,0,3,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=201
   n_axest=(/0,0,8,0,1,0,0,1,0,8,  0,0,0,0,0,0,0,3,0,3,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=205
   n_axest=(/0,0,0,0,0,0,0,1,3,8,  0,6,0,0,0,0,0,0,0,0,  6,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=207
   n_axest=(/0,0,0,0,0,0,0,1,3,8,  0,0,0,0,0,0,0,0,0,0,  6,0,0,0,6,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=208
   n_axest=(/0,0,0,0,0,0,0,1,0,8,  0,0,0,0,0,0,0,0,0,3,  6,0,0,0,0,6,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=212
   n_axest=(/0,0,0,0,0,0,0,1,0,8,  0,0,0,0,0,0,0,0,0,3,  6,0,0,6,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=213

   n_axest=(/0,6,0,0,0,0,0,1,3,8,  0,0,0,0,6,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=215
   n_axest=(/0,6,0,0,0,0,0,1,3,8,  0,0,0,0,0,0,0,6,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=218

 case(48)

   n_axest=(/0,6,8,0,1,0,0,1,3,8,  0,6,0,0,9,0,0,0,0,0,  6,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=221
   n_axest=(/0,6,8,0,1,0,0,1,3,8,  0,6,0,0,0,0,0,9,0,0,  6,0,0,0,0,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=222
   n_axest=(/0,6,8,0,1,0,0,1,3,8,  0,0,0,0,3,0,0,6,0,0,  6,0,0,0,6,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=223
   n_axest=(/0,6,8,0,1,0,0,1,3,8,  0,0,0,0,6,0,0,3,0,0,  6,0,0,0,6,0,0,0,0,0,0/)
   if(sum((n_axes-n_axest)**2)==0) spgroup=224

 end select

end subroutine symlist_prim
!!***

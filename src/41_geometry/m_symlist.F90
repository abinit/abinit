!!****m* ABINIT/m_symlist
!! NAME
!!  m_symlist
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group (RC)
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

module m_symlist

 use defs_basis

 implicit none

 private
!!***

 public :: symlist_prim
 public :: symlist_bcc
 public :: symlist_fcc
 public :: symlist_others
!!***

contains
!!***

!!****f* m_symlist/symlist_prim
!! NAME
!! symlist_prim
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
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

subroutine symlist_prim(additional_info,nsym,n_axes,spgroup)

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

!!****f* m_symlist/symlist_bcc
!! NAME
!! symlist_bcc
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! BCC case.
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

subroutine symlist_bcc(additional_info,nsym,n_axes,spgroup)

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

!!****f* m_symlist/symlist_fcc
!! NAME
!! symlist_fcc
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! FCC case
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

subroutine symlist_fcc(nsym,n_axes,spgroup)

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

!!****f* m_symlist/symlist_others
!! NAME
!! symlist_others
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! Non primitive, non BCC, non FCC case : rhombohedral or face centered
!!
!! INPUTS
!! brvltt=Bravais lattice type
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

subroutine symlist_others(brvltt,nsym,n_axes,spgroup)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: brvltt,nsym
 integer,intent(out) :: spgroup
!arrays
 integer,intent(in) :: n_axes(31)

!Local variables-------------------------------
!character(len=500) :: message
!arrays
 integer :: n_axest(31)

!**************************************************************************

!DEBUG
!write(std_out,*) ' symlist_others : enter '
!write(std_out,*) ' nsym = ', nsym
!write(std_out,*) ' brvltt = ',brvltt
!write(std_out,'(a,10i3)' ) ' n_axes(1:10) =',n_axes(1:10)
!write(std_out,'(a,10i3)' ) ' n_axes(11:20)=',n_axes(11:20)
!write(std_out,'(a,11i3)' ) ' n_axes(21:31)=',n_axes(21:31)
!ENDDEBUG

 spgroup=0

 if(brvltt==4 .or. brvltt==5 .or. brvltt==6)then       !  A, B, C face centered

   select case(nsym)

   case(4)

     n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,0,0,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=5
     n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,1,1,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=8
     n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,0,1,0,1,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=9

   case(8)

     n_axest=(/0,0,0,0,0,0,1,1,0,0,  0,0,0,0,1,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=36

     n_axest=(/0,0,0,0,2,0,1,1,1,0,  0,0,0,0,1,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=12
     n_axest=(/0,0,0,0,2,0,1,1,1,0,  0,0,0,0,0,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=15
     n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,2,1,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=38
     n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,1,3,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=39
     n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,1,1,0,2,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=40
     n_axest=(/0,0,0,0,0,0,1,1,1,0,  0,0,0,0,0,3,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=41

     n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,0,0,0,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=20
     n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,2,2,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=35
     n_axest=(/0,0,0,0,0,0,1,1,2,0,  0,0,0,0,0,2,0,2,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=37

     n_axest=(/0,0,0,0,0,0,1,1,4,0,  0,0,0,0,0,0,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=21

   case(16)

     n_axest=(/0,0,0,0,2,0,1,1,2,0,  0,0,0,0,2,2,0,2,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=63
     n_axest=(/0,0,0,0,2,0,1,1,2,0,  0,0,0,0,1,4,0,1,0,4,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=64

     n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,3,2,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=65
     n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,2,4,0,0,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=67
     n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,0,4,0,2,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=68
     n_axest=(/0,0,0,0,2,0,1,1,4,0,  0,0,0,0,1,2,0,3,0,2,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=66

   end select

 else if(brvltt==7)then                ! Rhombohedral lattice

   select case(nsym)

   case(3)

     n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=146

   case(6)

     n_axest=(/0,0,2,0,1,0,0,1,0,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=148
     n_axest=(/0,0,0,0,0,0,0,1,3,2,  0,0,0,0,0,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=155
     n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=160
     n_axest=(/0,0,0,0,0,0,0,1,0,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=161

   case(12)

     n_axest=(/0,0,2,0,1,0,0,1,3,2,  0,0,0,0,3,0,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=166
     n_axest=(/0,0,2,0,1,0,0,1,3,2,  0,0,0,0,0,3,0,0,0,0,  0,0,0,0,0,0,0,0,0,0,0/)
     if(sum((n_axes-n_axest)**2)==0) spgroup=167

   end select

 end if

end subroutine symlist_others
!!***

end module m_symlist
!!***

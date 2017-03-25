!{\src2tex{textfont=tt}}
!!****f* ABINIT/symlist_others
!! NAME
!! symlist_others
!!
!! FUNCTION
!! Determine the space group from the number and type of symmetry operations
!! Non primitive, non BCC, non FCC case : rhombohedral or face centered
!!
!! COPYRIGHT
!! Copyright (C) 2000-2017 ABINIT group (RC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symlist_others(brvltt,nsym,n_axes,spgroup)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symlist_others'
!End of the abilint section

 implicit none

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

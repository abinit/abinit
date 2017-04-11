!{\src2tex{textfont=tt}}
!!****f* ABINIT/proc_distrb_cycle
!! NAME 
!!  proc_distrb_cycle
!!
!! FUNCTION
!!  test a condition to cycle
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2017 ABINIT group (FJ)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
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

function proc_distrb_cycle(distrb,ikpt,iband1,iband2,isppol,me) 

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'proc_distrb_cycle'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ikpt,iband1,iband2,isppol,me
 integer,allocatable,intent(in) :: distrb(:,:,:)
 logical :: proc_distrb_cycle

! *************************************************************************
 proc_distrb_cycle=.false.
 if (allocated(distrb)) then
   if (isppol==-1) then
     proc_distrb_cycle=(minval(abs(distrb(ikpt,iband1:iband2,:)-me))/=0)
   else
     proc_distrb_cycle=(minval(abs(distrb(ikpt,iband1:iband2,isppol)-me))/=0)
   end if
 end if

end function proc_distrb_cycle
!!***

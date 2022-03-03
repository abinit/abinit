!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfptlw_pert
!! NAME
!!  m_dfptlw_pert
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! COPYRIGHT
!!  Copyright (C) 2022 ABINIT group (FIXME: add author)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_dfptlw_pert
    
 use defs_basis
 use m_profiling_abi
 use m_errors

 implicit none

 private

! *************************************************************************

contains 
!!***

!!****f* ABINIT/m_dfptlw_pert/dummy_routine
!! NAME
!!  dummy_routine
!!
!! FUNCTION
!!  FIXME: add description.
!!
!! INPUTS
!!  argin(sizein)=description
!!
!! OUTPUT
!!  argout(sizeout)=description
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine dummy_routine(argin,argout,option,sizein,sizeout)
    
 use defs_basis

 implicit none

!Arguments ------------------------------------
 integer , intent(in)  :: option,sizein,sizeout
 integer , intent(in)  :: argin(sizein)
 integer , intent(out) :: argout(sizeout)
 real(dp), intent(out) ::                        

!Local variables-------------------------------
 integer ::                                      
 real(dp) ::                                     
!character(len=500) :: msg                   
 
! *************************************************************************

 DBG_ENTER("COLL")
 
! if (option/=1 .and. option/=2 ) then
!  write(msg,'(3a,i0)')&
!&  'The argument option should be 1 or 2,',ch10,&
!&  'however, option=',option
!  MSG_BUG(msg)
! end if
!
! if (sizein<1) then
!  write(msg,'(3a,i0)')&
!&  '  The argument sizein should be a positive number,',ch10,&
!&  '  however, sizein=',sizein
!  MSG_ERROR(msg)
! end if

 DBG_EXIT("COLL")

end subroutine dummy_routine
!!***

end module m_dfptlw_pert
!!***

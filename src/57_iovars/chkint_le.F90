!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkint_le
!! NAME
!! chkint_le
!!
!! FUNCTION
!! Checks the value of an input integer variable, 
!! expected to be lower than some value, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkint_le,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint_le.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! minmax_value=see the description of minmax_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECT
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number) or list_values(list_number)
!! must be between -99 and 999 to be printed correctly.
!!
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      chkint_prt
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkint_le(advice_change_cond,cond_number,cond_string,cond_values,&
&  ierr,input_name,input_value,&
&  minmax_value,unit)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkint_le'
 use interfaces_57_iovars, except_this_one => chkint_le
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value
 integer,intent(in) :: minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: list_number,minmax_flag,ok
 integer, allocatable :: list_values(:)

!******************************************************************

!Checks the allowed values
 ok=0
 minmax_flag=-1
 if(input_value<=minmax_value)ok=1
!write(std_out,*)' chkint_le : input_value,minmax_value=',input_value,minmax_value

 list_number=1
 ABI_ALLOCATE(list_values,(1))
 list_values=minmax_value

!If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
&   ierr,input_name,input_value,&
&   list_number,list_values,minmax_flag,minmax_value,unit)
 end if

 ABI_DEALLOCATE(list_values)

! reset all cond_strings
 cond_string(:)='#####'

end subroutine chkint_le
!!***

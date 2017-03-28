!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkint
!! NAME
!! chkint
!!
!! FUNCTION
!! Checks the value of an input integer variable, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkint,
!! and these are mentioned in the error message.
!!
!! See the examples in the NOTES
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! advice_change_cond= if 1, and if an error is detected, will
!!  advice to change the value of the conditions.
!! cond_number= number of conditions checked before calling chkint.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! list_number=number of allowed values (maximum 40).
!! list_values=list of allowed values
!! minmax_flag=if 0, only values in the list are allowed
!!              if 1, admit values larger or equal to minmax_value
!!              if -1, admit values smaller or equal to minmax_value
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
!! in order to ask only for a minimal value, set list_number
!! as well as minmax_flag to 1, and put the minimal value in both
!! list_values and minmax_value.
!!
!! Examples :
!!  List of values - ionmov must be equal to 0, 1, 3, 8, or 9
!!   call chkint(0,0,cond_string,cond_values,ierr,&
!!  & 'ionmov',ionmov,5,(/0,1,3,8,9/),0,0,iout)
!!
!!  Larger or equal to a given value - nberry >= limit
!!   call chkint(0,0,cond_string,cond_values,ierr,&
!!  & 'nberry',nberry,1,(/limit/),1,limit,iout)
!!
!!  Smaller or equal to a given value - nberry <= limit
!!   call chkint(0,0,cond_string,cond_values,ierr,&
!!  & 'nberry',nberry,1,(/limit/),-1,limit,iout)
!!
!!  Conditional cases (examples to be provided - see chkinp.f for the
!!  time being)
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


subroutine chkint(advice_change_cond,cond_number,cond_string,cond_values,&
&  ierr,input_name,input_value,&
&  list_number,list_values,minmax_flag,minmax_value,unit)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkint'
 use interfaces_57_iovars, except_this_one => chkint
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,input_value,list_number
 integer,intent(in) :: minmax_flag,minmax_value,unit
 integer,intent(inout) :: ierr
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4),list_values(list_number)
 character(len=*),intent(inout) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: ilist,ok

!******************************************************************

!Checks the allowed values
 ok=0
 if(list_number>0)then
   do ilist=1,list_number
     if(input_value == list_values(ilist))ok=1
   end do
 end if
 if(minmax_flag==1 .and. input_value>=minmax_value)ok=1
 if(minmax_flag==-1 .and. input_value<=minmax_value)ok=1

!If there is something wrong, compose the message, and print it
 if(ok==0)then
   call chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
&   ierr,input_name,input_value,&
&   list_number,list_values,minmax_flag,minmax_value,unit)
 end if
 
! reset all cond_strings
 cond_string(:)='#####'


end subroutine chkint
!!***

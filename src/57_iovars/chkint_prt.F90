!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkint_prt
!! NAME
!! chkint_prt
!!
!! FUNCTION
!! During the checking of the value of a variable,
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkval,
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
!!              if 2, values in the list are not allowed
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
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,&
!!  & 'ionmov',ionmov,5,(/0,1,3,8,9/),0,0,iout)
!!
!!  Larger or equal to a given value - nberry >= limit
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,&
!!  & 'nberry',nberry,1,(/limit/),1,limit,iout)
!!
!!  Smaller or equal to a given value - nberry <= limit
!!   call chkint_prt(0,0,cond_string,cond_values,ierr,&
!!  & 'nberry',nberry,1,(/limit/),-1,limit,iout)
!!
!!  Conditional cases (examples to be provided - see chkinp.f for the
!!  time being)
!!
!! PARENTS
!!      chkint,chkint_eq,chkint_ge,chkint_le,chkint_ne
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkint_prt(advice_change_cond,cond_number,cond_string,cond_values,&
&  ierr,input_name,input_value,&
&  list_number,list_values,minmax_flag,minmax_value,unit)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkint_prt'
 use interfaces_14_hidewrite
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
 character(len=*),intent(in) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: icond
 character(len=500) :: message

!******************************************************************

 if(cond_number<0 .or. cond_number>4)then
   write(message,'(a,i0,a)' )&
&   'The value of cond_number is ',cond_number,' but it should be positive and < 5.'
   MSG_BUG(message)
 end if

 if(list_number<0 .or. list_number>40)then
   write(message,'(a,i0,a)' )&
&   'The value of list_number is',list_number,' but it should be between 0 and 40.'
   MSG_BUG(messagE)
 end if

!Compose the message, and print it
 ierr=1
 write(message, '(2a)' ) ch10,' chkint_prt: ERROR -'
 if(cond_number/=0)then
   do icond=1,cond_number
!    The following format restricts cond_values(icond) to be between -99 and 999
     write(message, '(5a,i0,a)' ) trim(message),ch10,&
&     ' Context: the value of the variable ',trim(cond_string(icond)),' is ',cond_values(icond),'.'
   end do
 end if
 write(message, '(5a,i0,a)' ) trim(message),ch10,&
& '  The value of the input variable ',trim(input_name),' is ',input_value,', while it must be'
 if(minmax_flag==2)then
   write(message, '(3a,20(i0,1x))' ) trim(message),ch10,&
   '  different from one of the following:',list_values(1:list_number)
 else if(list_number>1 .or. &
&   minmax_flag==0 .or. list_values(1)/=minmax_value )then
!  The following format restricts list_values to be between -99 and 999
   if(list_number/=1)then
     write(message, '(3a,40(i0,1x))' ) trim(message),ch10,&
     '  equal to one of the following: ',list_values(1:list_number)
   else
     write(message, '(3a,40(i0,1x))' ) trim(message),ch10,&
     '  equal to ',list_values(1)
   end if
   if(minmax_flag==1)then
!    The following format restricts minmax_value to be between -99 and 999
     write(message, '(3a,i0,a)' ) trim(message),ch10,&
&     '  or it must be larger or equal to ',minmax_value,'.'
   else if(minmax_flag==-1)then
     write(message, '(3a,i0,a)' ) trim(message),ch10,&
&     '  or it must be smaller or equal to ',minmax_value,'.'
   end if
 else if(minmax_flag==1)then
!  The following format restricts minmax_value to be between -99 and 999
   write(message, '(3a,i0,a)' ) trim(message),ch10,&
&   '  larger or equal to ',minmax_value,'.'
 else if(minmax_flag==-1)then
!  The following format restricts minmax_value to be between -99 and 999
   write(message, '(3a,i0,a)' ) trim(message),ch10,&
&   '  smaller or equal to ',minmax_value,'.'
 end if
 if(cond_number==0 .or. advice_change_cond==0)then
   write(message, '(5a)' ) trim(message),ch10,&
&   '  Action: you should change the input variable ',trim(input_name),'.'
 else if(cond_number==1)then
   write(message, '(7a)' ) trim(message),ch10,&
&   '  Action: you should change the input variables ',trim(input_name),' or ',trim(cond_string(1)),'.'
 else if(cond_number==2)then
   write(message, '(11a)' ) trim(message),ch10,&
&   '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
&   '   ',trim(cond_string(1)),' or ',trim(cond_string(2)),'.'
 else if(cond_number==3)then
   write(message, '(13a)' ) trim(message),ch10,&
&   '  Action: you should change one of the input variables ',trim(input_name),',',ch10,&
&   '   ',trim(cond_string(1)),', ',trim(cond_string(2)),' or ',trim(cond_string(3)),'.'
 end if
 call wrtout(unit   ,message,'COLL')
 call wrtout(std_out,message,'COLL')

end subroutine chkint_prt
!!***

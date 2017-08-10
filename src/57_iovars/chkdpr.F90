!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkdpr
!! NAME
!! chkdpr
!!
!! FUNCTION
!! Checks the value of an input real(dp) variable, and
!! write a sophisticated error message when it is erroneous.
!! A few conditions might have been checked before calling chkdpr,
!! and these are mentioned in the error message.
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
!! cond_number= number of conditions checked before calling chkdpr.
!! cond_string(cond_number)= name of the variables associated to the conditions.
!! cond_values(cond_number)= value of the variables associated to the conditions.
!! input_name=name of the input variable to be checked
!! input_value=value of the input variable to be checked
!! minimal_flag=if 0, the reference_value must be matched within 1.0d-10
!!              if 1, admit values larger or equal to reference_value
!!              if -1, admit values smaller or equal to reference_value
!! reference_value=see the description of minimal_flag
!! unit=unit number for clean output file
!!
!! OUTPUT
!!  (only side effect)
!!
!! SIDE EFFECTS
!! ierr= switch it to 1 if an error was detected. No action otherwise.
!!
!! NOTES
!! cond_values(cond_number)
!! must be between -99 and 999 to be printed correctly.
!! for the time being, at most 3 conditions are allowed.
!!
!! PARENTS
!!      chkinp
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkdpr(advice_change_cond,cond_number,cond_string,cond_values,&
&  ierr,input_name,input_value,minimal_flag,reference_value,unit)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkdpr'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: advice_change_cond,cond_number,minimal_flag,unit
 integer,intent(inout) :: ierr
 real(dp),intent(in) :: input_value,reference_value
 character(len=*),intent(in) :: input_name
!arrays
 integer,intent(in) :: cond_values(4)
 character(len=*),intent(in) :: cond_string(4)

!Local variables-------------------------------
!scalars
 integer :: icond,ok
 character(len=500) :: message

!******************************************************************

 if(cond_number<0 .or. cond_number>4)then
   write(message,'(a,i6,a)' )&
&   'The value of cond_number is',cond_number,&
&   'but it should be positive and < 5.'
   MSG_BUG(message)
 end if

!Checks the allowed values
 ok=0
 if(minimal_flag==1 .and. input_value>=reference_value-tol10)              ok=1
 if(minimal_flag==-1 .and. input_value<=reference_value+tol10)             ok=1
 if(minimal_flag==0 .and. abs(input_value-reference_value)<=tol10) ok=1

!If there is something wrong, compose the message, and print it
 if(ok==0)then
   ierr=1
   write(message, '(a,a)' ) ch10,' chkdpr: ERROR -'
   if(cond_number/=0)then
     do icond=1,cond_number
!      The following format restricts cond_values(icond) to be between -99 and 999
       write(message, '(2a,a,a,a,i4,a)' ) trim(message),ch10,&
&       '  Context : the value of the variable ',&
&       trim(cond_string(icond)),' is',cond_values(icond),'.'
     end do
   end if
   write(message, '(2a,a,a,a,es20.12,a)' ) trim(message),ch10,&
&   '  The value of the input variable ',trim(input_name),&
&   ' is',input_value,','
   if(minimal_flag==0)then
     write(message, '(2a,a,es20.12,a)' ) trim(message),ch10,&
     '  while it must be equal to ',reference_value,'.'
   else if(minimal_flag==1)then
     write(message, '(2a,a,es20.12,a)' ) trim(message),ch10,&
&     '  while it must be larger or equal to',reference_value,'.'
   else if(minimal_flag==-1)then
     write(message, '(2a,a,es20.12,a)' ) trim(message),ch10,&
&     '  while it must be smaller or equal to',reference_value,'.'
   end if

   if(cond_number==0 .or. advice_change_cond==0)then
     write(message, '(2a,a,a,a)' ) trim(message),ch10,&
&     '  Action : you should change the input variable ',trim(input_name),'.'
   else if(cond_number==1)then
     write(message, '(2a,a,a,a,a,a)' ) trim(message),ch10,&
&     '  Action : you should change the input variables ',trim(input_name),&
&     ' or ',trim(cond_string(1)),'.'
   else if(cond_number==2)then
     write(message, '(2a,a,a,a,a,a,a,a,a,a)' ) trim(message),ch10,&
&     '  Action : you should change one of the input variables ',&
&     trim(input_name),',',ch10,&
&     '   ',trim(cond_string(1)),' or ',trim(cond_string(2)),'.'
   else if(cond_number==3)then
     write(message, '(2a,a,a,a,a,a,a,a,a,a,a,a)' ) trim(message),ch10,&
&     '  Action : you should change one of the input variables ',&
&     trim(input_name),',',ch10,&
&     '   ',trim(cond_string(1)),', ',trim(cond_string(2)),&
&     ' or ',trim(cond_string(3)),'.'
   end if

   call wrtout(unit,message,'COLL')
   !call wrtout(std_out,  message,'COLL')
   MSG_WARNING(message)
 end if

end subroutine chkdpr
!!***

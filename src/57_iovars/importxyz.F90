!{\src2tex{textfont=tt}}
!!****f* ABINIT/importxyz
!! NAME
!! importxyz
!!
!! FUNCTION
!! Examine the input string, to see whether data from xyz
!! file(s) has to be incorporated.
!! For each such xyz file, translate the relevant
!! information into intermediate input variables compatible
!! with the usual ABINIT formatting, then append it
!! to the input string.
!!
!! COPYRIGHT
!! Copyright (C) 2002-2016 ABINIT group (MJV).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  string_raw*(strln)=raw string of character from input file (with original case)
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of character in string
!!  string_upper*(strln)=string of character
!!   the string (with upper case) from the input file, to which the xyz data are appended to it
!!
!! PARENTS
!!      m_ab7_invars_f90,parsefile
!!
!! CHILDREN
!!      append_xyz,incomprs,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine importxyz (lenstr,string_raw,string_upper,strln)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'importxyz'
 use interfaces_14_hidewrite
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => importxyz
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=*),intent(in) :: string_raw
 character(len=*),intent(inout) :: string_upper

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: dtset_len,ixyz,ii,index_already_done,index_xyz_fname
 integer :: index_xyz_fname_end,index_xyz_token,kk
 character(len=2) :: dtset_char
 character(len=500) :: message
 character(len=fnlen) :: xyz_fname

!************************************************************************

 index_already_done=1
 ixyz=0

 do    ! Infinite do-loop, to identify the presence of the xyzFILE token

   index_xyz_token=index(string_upper(index_already_done:lenstr),"XYZFILE")
   if(index_xyz_token==0)exit

   ixyz=ixyz+1
   if(ixyz==1)then
     write(message,'(80a)')('=',ii=1,80)
     call wrtout(ab_out,message,'COLL')
   end if

!  The xyzFILE token has been identified
   index_xyz_token=index_already_done+index_xyz_token-1

!  Find the related dataset tag, and length
   dtset_char=string_upper(index_xyz_token+7:index_xyz_token+8)
   if(dtset_char(1:1)==blank)dtset_char(2:2)=blank
   dtset_len=len_trim(dtset_char)

!  Find the name of the xyz file
   index_xyz_fname=index_xyz_token+8+dtset_len
   index_xyz_fname_end=index(string_upper(index_xyz_fname:lenstr),blank)

   if(index_xyz_fname_end ==0 )then
     write(message, '(5a,i4,2a)' )&
&     'Could not find the name of the xyz file.',ch10,&
&     'index_xyz_fname_end should be non-zero, while it is :',ch10,&
&     'index_xyz_fname_end=',index_xyz_fname_end,ch10,&
&     'Action: check the filename that was provided after the XYZFILE input variable keyword.'
     MSG_ERROR(message)
   end if

   index_xyz_fname_end=index_xyz_fname_end+index_xyz_fname-1

   index_already_done=index_xyz_fname_end

   xyz_fname=repeat(blank,fnlen)                  ! Initialize xyz_fname to a blank line
   xyz_fname=string_raw(index_xyz_fname:index_xyz_fname_end-1)

   write(message, '(3a)') ch10,&
&   ' importxyz : Identified token XYZFILE, referring to file ',trim(xyz_fname)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Append the data from the xyz file to the string, and update the length of the string
   call append_xyz(dtset_char,lenstr,string_upper,xyz_fname,strln)

!  erase the file name from string_upper
   string_upper(index_xyz_fname:index_xyz_fname_end-1) = blank

 end do

 if (index_already_done > 1) then 
   xyz_fname=repeat(blank,fnlen) ! Initialize xyz_fname to a blank line
   call append_xyz("-1",lenstr,string_upper,xyz_fname,strln)
 end if

 if(ixyz/=0)then
   call incomprs(string_upper,lenstr)
!  A blank is needed at the beginning of the string
   do kk=lenstr,1,-1
     string_upper(kk+1:kk+1)=string_upper(kk:kk)
   end do
   string_upper(1:1)=blank
   lenstr=lenstr+1
   write(message,'(a,80a,a)')ch10,('=',ii=1,80),ch10
   call wrtout(ab_out,message,'COLL')
 end if

end subroutine importxyz
!!***

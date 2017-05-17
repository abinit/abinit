!{\src2tex{textfont=tt}}
!!****f* ABINIT/incomprs
!! NAME
!! incomprs
!!
!! FUNCTION
!! Compresses input character string into the following form:
!! (1) Replaces tabs and all other characters lexically less than
!! SP (blank) with SP (blank), where lexically less than refers to
!! the ASCII collating sequence (SP is hex 20, dec 32).
!! The use of llt is needed e.g. on the IBM 9000 because it does not
!! handle tab characters sensibly in its AIX fortran.
!! Also replace occurences of '=' by a SP.
!! (2) Removes all repeated blanks, ignoring trailing blanks
!! after first (returns nontrailing final length in arg 'length').
!! (3) Makes first character in string NONBLANK.  This is done
!! to prevent double blanks from occurring when compressed string
!! is concatenated with other compressed strings.
!! (4) Makes last character (string(length:length)) a blank.
!! If input string is entirely blank or tabs, simply returns with length=0.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (see side effects)
!!
!! OUTPUT
!!  length=nonblank, nontab length of string as defined above
!!
!! SIDE EFFECT
!!  string=at input:  character string
!!         at output: repeated blanks and tabs have been removed and
!!                    remaining tabs have been replaced by blanks
!!
!! PARENTS
!!      importxyz,instrng
!!
!! CHILDREN
!!      inreplsp
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine incomprs(string,length)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'incomprs'
 use interfaces_42_parser, except_this_one => incomprs
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(out) :: length
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
 character(len=1) :: blank=' '
!scalars
 integer :: bb,f1,ii,jj,kk,l1,lbef,lcut,lold,stringlen
 character(len=500) :: message

! *************************************************************************

!
!String length determined by calling program declaration of "string"
 stringlen=len(string)
 length=stringlen
!
!Only proceed if string has nonzero length
 if (length>0) then
!  Find last nonblank character (i.e. nonblank and nontab length)
   length=len_trim(string)
   if (length==0) then
!    Line is all blanks or tabs so do not proceed
!    write(std_out,*)' incomprs: blank line encountered'
   else

!    Replace all characters lexically less than SP, and '=', by SP (blank)
     call inreplsp(string(1:length))

!    Continue with parsing
!    l1 is set to last nonblank, nontab character position
     l1=length
     do ii=1,l1
       if (string(ii:ii)/=blank) exit
     end do

!    f1 is set to first nonblank, nontab character position
     f1=ii
!    lbef is number of characters in string starting at
!    first nonblank, nontab and going to last
     lbef=l1-f1+1

!    Process characters one at a time from right to left:
     bb=0
     lcut=lbef
     do ii=1,lbef
       jj=lbef+f1-ii
!      set bb=position of next blank coming in from right
       if (string(jj:jj)==blank) then
         if (bb==0) then
           bb=jj
         end if
       else
         if (bb/=0) then
!          if several blanks in a row were found, cut from string
           if (jj<bb-1) then
!            lold becomes string length before cutting blanks
             lold=lcut
!            lcut will be new string length
             lcut=lcut-(bb-1-jj)
!            redefine string with repeated blanks gone
             do kk=1,f1+lcut-1-jj
               string(jj+kk:jj+kk)=string(kk+bb-1:kk+bb-1)
             end do
           end if
           bb=0
         end if
       end if
     end do
!    
!    Remove initial blanks in string if any
     if (f1>1) then
       string(1:lcut)=string(f1:f1+lcut-1)
     end if
!    
!    Add blank on end unless string had no extra space
     if (lcut==stringlen) then
       write(message,'(a,i7,a,a,a,a,a,a,a,a)')&
&       'For input file, with data forming a string of',stringlen,' characters,',ch10,&
&       'no double blanks or tabs were found.',ch10,&
&       'This is unusual for an input file (or any file),',ch10,&
&       'and may cause parsing trouble.  Is this a binary file?',ch10
       MSG_WARNING(message)
     else
       length=lcut+1
       string(length:length)=blank
     end if
   end if
 end if

end subroutine incomprs
!!***

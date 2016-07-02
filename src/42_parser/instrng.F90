!{\src2tex{textfont=tt}}
!!****f* ABINIT/instrng
!! NAME
!! instrng
!!
!! FUNCTION
!! Read the input file, and product a string of character,
!! with all data, to be analyzed in later routines. The length
!! of this string is lenstr. This number is checked to be smaller
!! than the dimension of the string of character, namely strln .
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  filnam=name of the input file, to be read
!!  option= if 0, simple storing of the character string,
!!             no special treatment for ABINIT (comment delimiters, checks, include ...)
!!          if 1, suppresses text after an ABINIT comment delimiter (! or #),
!!             checks that a minus sign is followed by a number ...
!!                check for INCLUDE statement:
!!                if present, add string from included file
!!  strln=maximal number of character of string, as declared in the calling routine
!!
!! OUTPUT
!!  lenstr=actual number of character in string
!!  string*(strln)=string of character
!!
!! PARENTS
!!      anaddb,importcml,localorb_S,lwf,parsefile
!!
!! CHILDREN
!!      incomprs,instrng,wrtout
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


recursive subroutine instrng(filnam,lenstr,option,strln,string)

 use defs_basis
 use m_profiling_abi
 use m_errors

 use m_io_tools,  only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'instrng'
 use interfaces_14_hidewrite
 use interfaces_42_parser, except_this_one => instrng
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: option,strln
 integer,intent(out) :: lenstr
 character(len=*),intent(in) :: filnam
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer,save :: include_level=-1
 integer :: ii,ii1,ii2,ij,iline,ios,iost,lenc,lenstr_inc,mline,nline1,input_unit
 logical :: include_found,ex
 character(len=1) :: string1
 character(len=3) :: string3
 character(len=500) :: filnam_inc,msg
 character(len=fnlen+20) :: line
 character(len=strlen),pointer :: string_inc

!************************************************************************

 DBG_ENTER("COLL")

!%%%%%%%%%%%%%%%%%%%%%%%%
!read in string from file
!%%%%%%%%%%%%%%%%%%%%%%%%

!The file can be included in another (prevent too many include levels)
 include_level=include_level+1
 if (include_level>2) then
   write(msg, '(3a)' ) &
&   'At least 4 levels of included files are present in input file !',ch10,&
&   'This is not allowed. Action: change your input file.'
   MSG_ERROR(msg)
 end if

!Open data file and read one line at a time, compressing data
!and concatenating into single string:
 if (open_file(filnam,msg,newunit=input_unit,form="formatted",status="old",action="read") /= 0) then
   MSG_ERROR(msg)
 end if
 rewind (unit=input_unit)

!Initialize string to blanks
 string=blank
 lenstr=1

!Set maximum number lines to be read to some large number
 mline=50000
 do iline=1,mline

!  Keeps reading lines until end of input file
   read (unit=input_unit,fmt= '(a)' ,iostat=ios) line(1:fnlen+20)
!  Hello ! This is a commentary. Please, do not remove me.
!  In fact, this commentary protect tests_v4 t47 for miscopying
!  the input file into the output string. It _is_ strange.
!  The number of lines in the commentary is also resulting from
!  a long tuning..

!  DEBUG
!  write(std_out,*)' instrng, iline=',iline,' ios=',ios,' echo :',trim(line(1:fnlen+20))
!  ENDDEBUG

!  Exit the reading loop when arrived at the end
   if(ios/=0)then
     backspace(input_unit)
     read (unit=input_unit,fmt= '(a1)' ,iostat=ios) string1
     if(ios/=0)exit
     backspace(input_unit)
     read (unit=input_unit,fmt= '(a3)' ,iostat=ios) string3
     if(string3=='end')exit
     write(msg, '(3a,i0,11a)' ) &
&     'It is observed in the input file: ',TRIM(filnam),', line number ',iline,',',ch10,&
&     'that there is a non-zero IO signal.',ch10,&
&     'This is normal when the file is completely read.',ch10,&
&     'However, it seems that the error appears while your file has not been completely read.',ch10,&
&     'Action: correct your file. If your file seems correct, then,',ch10,&
&     'add the keyword ''end'' at the very beginning of the last line of your input file.'
     MSG_ERROR(msg)
   end if

!  Find length of input line ignoring delimiter characters (# or !)
!  and any characters beyond it (allows for comments beyond # or !)
   ii1=index(line(1:fnlen+20),'#')
   ii2=index(line(1:fnlen+20),'!')
   if ( (ii1==0 .and. ii2==0) .or. option==0 ) then
!    delimiter character was not found on line so use full line
     ii=fnlen+20
   else if(ii1==0)then
!    ii will represent length of line up to but not including !
     ii=ii2-1
   else if(ii2==0)then
!    ii will represent length of line up to but not including #
     ii=ii1-1
   else
     ii=min(ii1,ii2)-1
   end if

!  Checks that nothing is left beyond fnlen
   if(ii>fnlen)then
     do ij=fnlen+1,ii
       if(line(ij:ij)/=' ')then
         write(msg,'(3a,i0,3a,i0,3a)' ) &
&         'It is observed in the input file: ',TRIM(filnam),' line number ',iline,',',ch10,&
&         'that more than ',fnlen,' columns are used.',ch10,&
&         'This is not allowed. Change this line of your input file.'
         MSG_ERROR(msg)
       end if
     end do
   end if

   if (ii>0) then
!    Check for the occurence of a minus sign followed by a blank
     ij=index(line(1:ii),'- ')
     if (ij>0 .and. option==1) then
       write(msg, '(3a,i0,11a)' ) &
&       'It is observed in the input file:, ',TRIM(filnam),' line number ',iline,',',ch10,&
&       'the occurence of a minus sign followed',ch10,&
&       'by a blank. This is forbidden.',ch10,&
&       'If the minus sign is meaningful, do not leave a blank',ch10,&
&       'between it and the number to which it applies.',ch10,&
&       'Otherwise, remove it.'
       MSG_ERROR(msg)
     end if
!    Check for the occurence of a tab
     ij=index(line(1:ii),char(9))
     if (ij>0 .and. option==1 ) then
       write(msg, '(3a,i0,3a)' ) &
&       'The occurence of a tab, in the input file: ',TRIM(filnam),' line number ',iline,',',ch10,&
&       'is observed. This sign is confusing, and has been forbidden.'
       MSG_ERROR(msg)
     end if

!    Check for the occurence of a include statement
     include_found=.false.
     if (option==1) then
!      Look for include statement
       ii1=index(line(1:ii),"include");ii2=index(line(1:ii),"INCLUDE")
       include_found=(ii1>0.or.ii2>0)
       if (include_found) then
         ij=max(ii1,ii2);ii1=0;ii2=0
!        Look for quotes (ascii 34)
         ii1=index(line(ij+7:ii),char(34))
         if (ii1>1) ii2=index(line(ij+7+ii1:ii),char(34))
!        Look for quotes (ascii 39)
         if (ii1==0.and.ii2==0) then
           ii1=index(line(ij+7:ii),char(39))
           if (ii1>1) ii2=index(line(ij+7+ii1:ii),char(39))
         end if
!        Check if quotes are correctly set
         ex=(ii1<=1.or.ii2<=1)
         if (.not.ex) then
           msg=line(ij+7:ij+5+ii1)
           call incomprs(msg(1:ii1-1),lenc)
           ex=(len(trim(msg))/=0)
         end if
         if (ex) then
           write(msg, '(6a)' ) &
&           'A "include" statement has been found in input file: ',TRIM(filnam),ch10,&
&           'but there must be a problem with the quotes.',ch10,&
&           'Action: change your input file.'
           MSG_ERROR(msg)
         end if
!        Store included file name
         filnam_inc=line(ij+7+ii1:ij+5+ii1+ii2)
!        Extract include statement from line
         lenc=ii1+ii2+7
         msg(1:ii-lenc)=line(1:ij-1)//line(ij+lenc:ii)
         ii=ii-lenc;line(1:ii)=msg(1:ii)
       end if
     end if

!    Compress: remove repeated blanks, make all ASCII characters
!    less than a blank (and '=') to become a blank.
     call incomprs(line(1:ii),lenc)

   else
!    ii=0 means line starts with #, is entirely a comment line
     lenc=0;include_found=.false.
   end if

!  Check resulting total string length
   if (lenstr+lenc>strln) then
     write(msg, '(8a)' ) &
&     'The size of your input file: ',TRIM(filnam),' is such that the internal',ch10,&
&     'character string that should contain it is too small.',ch10,&
&     'Action: decrease the size of your input file,',ch10,&
&     'or contact the ABINIT group.'
     MSG_ERROR(msg)
   end if

   if (lenc>0) then
!    Concatenate new compressed characters
!    with previous part of compressed string (unless all blank)
     string(lenstr+1:lenstr+lenc)=line(1:lenc)
   end if
!  Keep track of total string length
   lenstr=lenstr+lenc

!  Eventually (recursively) read included file
   if (include_found) then
!    Check file existence
     inquire(file=filnam_inc ,iostat=iost,exist=ex)
     if ((.not.ex).or.(iost/=0)) then
       write(msg, '(5a)' ) &
&       'Input file: ',TRIM(filnam),' reading: the included file ',trim(filnam_inc),' cannot be found !'
       MSG_ERROR(msg)
     end if
!    Read included file (warning: recursive call !)
     ABI_ALLOCATE(string_inc,)
     call instrng(trim(filnam_inc),lenstr_inc,option,strln-lenstr,string_inc)
!    Check resulting total string length
     if (lenstr+lenstr_inc>strln) then
       write(msg, '(6a)' ) &
&       'The size of your input file: ',TRIM(filnam),' (including included files) is such that',ch10,&
&       'the internal character string that should contain it is too small !',ch10,&
&       'Action : decrease the size of your input file.'
       MSG_ERROR(msg)
     end if
!    Concatenate total string
     string(lenstr+1:lenstr+lenstr_inc)=string_inc(1:lenstr_inc)
     lenstr=lenstr+lenstr_inc
     ABI_DEALLOCATE(string_inc)
   end if

!  If mline is reached, something is wrong
   if (iline>=mline) then
     write(msg, '(a,i0,2a,i0,4a)' ) &
&     'The number of lines already read from input file=',iline,ch10,&
&     'is equal or greater than maximum allowed mline=',mline,ch10,&
&     'Action: you could decrease the length of the input file, or',ch10,&
&     'contact the ABINIT group.'
     MSG_ERROR(msg)
   end if

!  End loop on iline. Note that there is an "exit" instruction in the loop
 end do

 nline1=iline-1
 close (unit=input_unit)

 write(msg,'(a,i0,3a)')'-instrng: ',nline1,' lines of input have been read from file ',trim(filnam),ch10
 call wrtout(std_out,msg,'COLL')

 include_level=include_level-1

 DBG_EXIT("COLL")

end subroutine instrng
!!***

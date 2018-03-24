!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_parser
!! NAME
!! m_parser
!!
!! FUNCTION
!! This module contains routines and functions used to
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (XG, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_parser

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_atomdata
 use m_xmpi

 use m_io_tools,  only : open_file
 use m_fstrings,  only : sjoin, itoa, inupper

 implicit none

 private

 public :: parsefile
 public :: inread
 public :: instrng
 public :: incomprs
 public :: intagm
 public :: importxyz

CONTAINS  !===========================================================
!!***

!!****f* m_parser/parsefile
!! NAME
!! parsefile
!!
!! FUNCTION
!!  Glue function, to read the given file, put it into a string,
!!  change everything to uppercase, remove carriage returns and
!!  non significant blank characters. May also read a XYZ input
!!  file if specified. Finally read ndtset input variable.
!!
!! INPUTS
!!  filnamin= the file to read
!!  comm=MPI communicator
!!
!! OUTPUT
!!  lenstr= the length of the resulting string.
!!  ndtset= the number of declared datasets.
!!  string= contains on output the content of the file, ready for parsing.
!!
!! PARENTS
!!      abinit,m_ab7_invars_f90,ujdet
!!
!! CHILDREN
!!      importxyz,instrng,intagm,inupper,xmpi_bcast
!!
!! SOURCE

subroutine parsefile(filnamin,lenstr,ndtset,string,comm)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'parsefile'
 use interfaces_42_parser
 use interfaces_57_iovars, except_this_one => parsefile
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 character(len=*),intent(in) :: filnamin
 integer,intent(in) :: comm
 integer,intent(out) :: ndtset,lenstr
 character(len=strlen),intent(out) :: string

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: option,marr,tread,lenstr_noxyz,ierr
 character(len=strlen) :: string_raw
 character(len=500) :: message
!arrays
 integer :: intarr(1)
 real(dp) :: dprarr(1)

! *************************************************************************

 ! Read the input file, and store the information in a long string of characters
 ! Note: this is done only by me=0, and then string and other output vars are BCASTED

 if (xmpi_comm_rank(comm) == master) then
   !strlen from defs_basis module
   option=1
   call instrng (filnamin,lenstr,option,strlen,string)

   ! Copy original file, without change of case
   string_raw=string

   ! To make case-insensitive, map characters of string to upper case:
   call inupper(string(1:lenstr))

   ! Might import data from xyz file(s) into string
   ! Need string_raw to deal properly with xyz filenames
   lenstr_noxyz = lenstr
   call importxyz(lenstr,string_raw,string,strlen)

   !6) Take ndtset from the input string
   ndtset=0; marr=1
   call intagm(dprarr,intarr,0,marr,1,string(1:lenstr),"ndtset",tread,'INT')
   if (tread==1) ndtset=intarr(1)
   ! Check that ndtset is not negative
   if (ndtset<0 .or. ndtset>9999) then
     write(message, '(a,i0,a,a,a,a)' )&
&     'Input ndtset must be non-negative and < 10000, but was ',ndtset,ch10,&
&     'This is not allowed.  ',ch10,&
&     'Action : modify ndtset in the input file.'
     MSG_ERROR(message)
   end if
 end if ! master

 if (xmpi_comm_size(comm) > 1) then
   ! Broadcast data.
   call xmpi_bcast(lenstr,master,comm,ierr)
   call xmpi_bcast(ndtset,master,comm,ierr)
   call xmpi_bcast(string,master,comm,ierr)
 end if

end subroutine parsefile
!!***

!!****f* m_parser/inread
!! NAME
!! inread
!!
!! FUNCTION
!! Carry out internal read from input character string, starting
!! at first character in string, reading ndig digits (including possible
!! sign, decimal, and exponent) by computing the appropriate format and
!! performing a formatted read (list-directed read would be perfect for
!! this application but is inconsistent with internal read according to
!! Fortran90 standard).
!! In case of a real number, this routine
!! is also able to read SQRT(number): return the square root of the number.
!!
!! INPUTS
!!  string=character string.
!!  ndig=length of field to be read (including signs, decimals, and exponents).
!!  typevarphys=variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR','LEN','ENE'=>real(dp) (no special treatment)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in token
!!
!! OUTPUT
!!  outi or outr (integer or real respectively)
!!  errcod, =0 for success, 1,2 for ini, inr failure resp.
!!
!! PARENTS
!!      adini,inarray
!!
!! CHILDREN
!!
!! SOURCE

subroutine inread(string,ndig,typevarphys,outi,outr,errcod)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inread'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ndig
 integer,intent(out) :: errcod,outi
 real(dp),intent(out) :: outr
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: typevarphys

!Local variables-------------------------------
!scalars
 integer :: done,idig,index_slash,sign
 real(dp) :: den,num
 logical :: logi
 character(len=500) :: msg

! *************************************************************************

!write(std_out,*)'inread : enter '
!write(std_out,*)'string(1:ndig)=',string(1:ndig)
!write(std_out,*)'typevarphys=',typevarphys

 if (typevarphys=='INT') then

!  integer input section
   read (unit=string(1:ndig),fmt=*,iostat=errcod) outi
   if(errcod/=0)then
!    integer reading error
     write(std_out,'(/,a,/,a,i0,a)' ) &
&     ' inread : ERROR -',&
&     '  Attempted to read ndig=',ndig,' integer digits,'
     write(std_out,'(a,a,a)' ) '   from string(1:ndig)= ',string(1:ndig),&
&     ', to initialize an integer variable'
     errcod=1
   end if

 else if (typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI') then

!  real(dp) input section

!  Special treatment of SQRT(xxx) or -SQRT(xxx) chains of characters, where xxx can be a fraction
   done=0
   if (ndig>5) then
     if(string(1:5)=='SQRT(' .and. string(ndig:ndig)==')')then
       done=1 ; sign=1
     else if(string(1:6)=='-SQRT(' .and. string(ndig:ndig)==')')then
       done=1 ; sign=2
     end if
     if(done==1)then
       index_slash=index(string(5+sign:ndig-1),'/')
       if(index_slash==0)then
         read (unit=string(5+sign:ndig-1),fmt=*,iostat=errcod) outr
       else if(index_slash/=0)then
         read (unit=string(5+sign:5+sign+index_slash-2),fmt=*,iostat=errcod) num
         if(errcod==0)then
           read (unit=string(5+sign+index_slash:ndig-1),fmt=*,iostat=errcod) den
           if(errcod==0)then
             if(abs(den)<tol12)then
               errcod=1
             else
               outr=num/den
             end if
           end if
         end if
       end if
       if(outr<-tol12)then
         errcod=1
       else
         outr=sqrt(outr)
         if(sign==2)outr=-outr
       end if
     end if
   end if

!  Special treatment of fractions
   if(done==0)then
     index_slash=index(string(1:ndig),'/')
     if(index_slash/=0)then
       done=1
       read (unit=string(1:index_slash-1),fmt=*,iostat=errcod) num
       if(errcod==0)then
         read (unit=string(index_slash+1:ndig),fmt=*,iostat=errcod) den
         if(errcod==0)then
           if(abs(den)<tol12)then
             errcod=1
           else
             outr=num/den
           end if
         end if
       end if
     end if
   end if

!  Normal treatment of floats
   if(done==0)then ! Normal treatment of float numbers
     read (unit=string(1:ndig),fmt=*,iostat=errcod) outr
   end if

!  Treatment of errors
   if(errcod/=0)then
!    real(dp) data reading error
     write(std_out,'(/,a,/,a,i0,a)' ) &
&     'inread : ERROR -',&
&     'Attempted to read ndig=',ndig,' floating point digits,'
     write(std_out,'(a,a,a)' ) '   from string(1:ndig) ',string(1:ndig),&
&     ', to initialize a floating variable.'
     errcod=2
   end if

 else if (typevarphys=='LOG') then

   read (unit=string(1:ndig),fmt=*,iostat=errcod) logi
   if(errcod/=0)then
!    integer reading error
     write(std_out,'(/,a,/,a,i0,a)' ) &
&     'inread : ERROR -',&
&     'Attempted to read ndig=',ndig,' integer digits,'
     write(std_out,'(a,a,a)' ) '   from string(1:ndig)= ',string(1:ndig),', to initialize a logical variable.'
     errcod=3
   end if
   if(logi)outi=1
   if(.not.logi)outi=0

 else
   write(msg,'(4a)' ) &
&   'Argument typevarphys must be INT,DPR,LEN,ENE,BFI or LOG ',ch10,&
&   'but input value was: ',trim(typevarphys)
   MSG_ERROR(msg)
 end if

 if(errcod /= 0)then
   do idig=1,ndig
     if( string(idig:idig) == 'O' )then
       write(std_out,'(/,a,/,a,a,a)' ) &
&       'inread : WARNING -',&
&       'Note that this string contains the letter O. ',ch10,&
&       'It is likely that this letter should be replaced by the number 0.'
       exit
     end if
   end do
 end if

end subroutine inread
!!***

!!****f* m_parser/instrng
!! NAME
!! instrng
!!
!! FUNCTION
!! Read the input file, and product a string of character,
!! with all data, to be analyzed in later routines. The length
!! of this string is lenstr. This number is checked to be smaller
!! than the dimension of the string of character, namely strln .
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

recursive subroutine instrng(filnam,lenstr,option,strln,string)

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

!!****f* m_parser/inreplsp
!! NAME
!! inreplsp
!!
!! FUNCTION
!! Replace all occurrences of characters lexically less than SP (blank)
!! by SP in the input string, returning modified string of same length.
!! Also replace a '=' by a SP.
!!
!! INPUTS
!!  string=character string to be modified
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  string=same character string with ASCII (decimal) 0-31 replaced by 32.
!!
!! PARENTS
!!      incomprs
!!
!! CHILDREN
!!
!! SOURCE

subroutine inreplsp(string)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inreplsp'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ilenth,length

! *************************************************************************

!Get length of string
 length=len(string)

!Proceed only if string has nonzero length
 if (length>0) then

!  Do replacement by going through input
!  character string one character at a time
   do ilenth=1,length
     if (llt(string(ilenth:ilenth),' ')) then
       string(ilenth:ilenth)=' '
     end if
     if(string(ilenth:ilenth)=='=')then
       string(ilenth:ilenth)=' '
     end if
   end do

 end if

end subroutine inreplsp
!!***

!!****f* m_parser/incomprs
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

subroutine incomprs(string,length)

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

!!****f* m_parser/intagm
!! NAME
!! intagm
!!
!! FUNCTION
!! Search input 'string' for specific 'token'. Search depends on
!! input dataset through 'jdtset'. Then, return the information
!! mentioned after 'token'.
!! See the "notes" section
!!
!! INPUTS
!!  jdtset=see the notes section
!!  marr=dimension of the intarr and dprarr arrays, as declared in the calling subroutine.
!!  narr=actual size of array to be read in.
!!  string=character string containing 'tags' and data.
!!  token=character string for 'tag'.
!!  typevarphys= variable type (might indicate the physical meaning of for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (expect a "length", identify bohr, au or angstrom,
!!       and return in au -atomic units=bohr- )
!!   'ENE'=>real(dp) (expect a "energy", identify Ha, hartree, eV, Ry, Rydberg)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in key_value
!!
!! OUTPUT
!!  intarr(1:narr), dprarr(1:narr)
!!   integer or real(dp) arrays, respectively (see typevarphys),
!!   into which data is read if typevarphys/='KEY'. Use these arrays even for scalars.
!!  tread is an integer : tread = 0 => no data was read
!!                        tread = 1 => data was read
!!  ds_input is an optional integer flag:
!!           ds_input = 0 => value was found which is not specific to jdtset
!!           ds_input > 0 => value was found which is specific to jdtset
!!   one could add more information, eg whether a ? or a : was used, etc...
!!   [key_value]=Stores the value of key if typevarphys=="KEY"
!!               String of len fnlen
!!
!! NOTES
!!
!! If jdtset==0:
!!
!!  Search compressed 'string' for blank//'token'//blank and
!!  read input data beside 'token', to be read into appropriate variable.
!!  For this routine to find a given token, the token has to be preceded
!!  and followed by blanks--i.e. the first token should not start out as
!!  the first character in the input file.  This is checked in the calling
!!  subroutine 'input'. Calls inread which performs internal read from
!!  specified string.  Also calls upper which maps characters to all upper case.
!!  Also checks whether there is an occurence of blank//'token'//digit,
!!  in which case the input file might be erroneous, so stops.
!!
!! If jdtset is a positive number:
!!
!!  (1) First search for modified string, blank//'token'//jdtset//blank
!!
!!  (2a) if the occurence of (1) is not found,
!!       look for other modified strings,
!!       blank//'token'//'?'//unities//blank
!!       or
!!       blank//'token'//dozens//'?'//blank
!!       (issue an error message if more than one occurs)
!!       where jdtset=dozens*10+unities (decimal decomposition of jdtset)
!!       if one of them exists, just take the value
!!       Note that unities is a one-digit number, while dozens might be bigger than 9.
!!
!!  (2b-2c) search for a series, with the following tokens :
!!       (issue an error message if more than one occurs, or
!!       goto (3) if none exist)
!!
!!      blank//'token'//':'//blank
!!      if it exists, then a series might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'+'//blank,
!!      or the multiplicative factor, blank//'token'//'*'//blank
!!
!!      blank//'token'//'?'//':'//blank
!!      if it exists, then a series for the inner loop
!!      might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'?'//'+'//blank,
!!      or the multiplicative factor, blank//'token'//'?'//'*'//blank
!!
!!      blank//'token'//':'//'?'//blank
!!      if it exists, then a series for the outer loop
!!      might have been defined in the input file
!!      must thus find either the increment, blank//'token'//'+'//'?'//blank,
!!      or the multiplicative factor, blank//'token'//'*'//'?'//blank
!!
!!  (3) if neither (1) nor (2) are found, search for the 'normal'
!!       string, blank//'token'//blank
!!
!!
!! PARENTS
!!      ingeo,ingeobld,inkpts,inqpt,invacuum,invars0,invars1,invars2
!!      m_ab7_invars_f90,m_anaddb_dataset,m_band2eps_dataset,m_ingeo_img
!!      m_multibinit_dataset,macroin,mpi_setup,parsefile,ujdet
!!
!! CHILDREN
!!      appdig,inarray,inupper,wrtout
!!
!! SOURCE

subroutine intagm(dprarr,intarr,jdtset,marr,narr,string,token,tread,typevarphys,ds_input,key_value)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'intagm'
 use interfaces_14_hidewrite
 use interfaces_32_util
 use interfaces_42_parser, except_this_one => intagm
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: jdtset,marr,narr
 integer,intent(out) :: tread
 integer,intent(out),optional :: ds_input
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: token
 character(len=*),intent(in) :: typevarphys
 character(len=fnlen),optional,intent(out) :: key_value
!arrays
 integer,intent(inout) :: intarr(marr) !vz_i
 real(dp),intent(inout) :: dprarr(marr) !vz_i

!Local variables-------------------------------
 character(len=1), parameter :: blank=' '
!scalars
 integer :: b1,cs1len,cslen,dozens,ier,itoken,itoken1,itoken2,itoken2_1colon
 integer :: itoken2_1plus,itoken2_1times,itoken2_2colon,itoken2_2plus
 integer :: itoken2_2times,itoken2_colon,itoken2_plus,itoken2_times
 integer :: itoken_1colon,itoken_1plus,itoken_1times,itoken_2colon,itoken_2plus
 integer :: itoken_2times,itoken_colon,itoken_plus,itoken_times,number,opttoken
 integer :: sum_token,toklen,trial_cslen,trial_jdtset,unities
 integer :: ds_input_
 character(len=4) :: appen
 character(len=3) :: typevar
 character(len=500) :: message
 character(len=fnlen) :: cs,cs1,cs1colon,cs1plus,cs1times,cs2colon,cs2plus
 character(len=fnlen) :: cs2times,cscolon,csplus,cstimes,trial_cs
!arrays
 integer,allocatable :: int1(:),int2(:)
 real(dp),allocatable :: dpr1(:),dpr2(:)

! *************************************************************************

 ABI_CHECK(marr >= narr, sjoin("marr", itoa(marr)," < narr ", itoa(narr)))

 ds_input_ = -1

 dozens=jdtset/10
 unities=jdtset-10*dozens

 if(jdtset<0)then
   write(message,'(a,i0,a)')' jdtset=',jdtset,', while it should be non-negative.'
   MSG_ERROR(message)
 end if

 if(jdtset>9999)then
   write(message,'(a,i0,a)')' jdtset=',jdtset,', while it must be lower than 10000.'
   MSG_ERROR(message)
 end if

!Default values : nothing has been read
 itoken=0
 opttoken=0
!Initialise flags in case of opttoken >= 2 later.
 itoken_times=0
 itoken_plus=0
 itoken_colon=0
 cslen=1

 if (narr/=0) then

   toklen=len_trim(token)

!  --------------------------------------------------------------------------
!  (1) try to find the token with dataset number appended
   if(jdtset>0)then

     call appdig(jdtset,'',appen)
     cs=blank//token(1:toklen)//trim(appen)//blank
     if(jdtset<10) then
       cslen=toklen+3
     else if(jdtset<100) then
       cslen=toklen+4
     else if(jdtset<1000) then
       cslen=toklen+5
     else if(jdtset<10000)then
       cslen=toklen+6
     end if
!    Map token to all upper case (make case-insensitive):
     call inupper(cs)
!    Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))
!    Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if(itoken/=itoken2)then
       write(message, '(7a)' )&
&       'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
&       'This is confusing, so it has been forbidden.',ch10,&
&       'Action: remove one of the two occurences.'
       MSG_ERROR(message)
     end if

     if(itoken/=0) then
       opttoken=1
       ds_input_=jdtset
     end if
   end if

!  --------------------------------------------------------------------------
!  (2a) try to find the token appended with a string that contains the metacharacter "?".

   if(jdtset>0 .and. opttoken==0)then

!    Use the metacharacter for the dozens, and save in cs and itoken
     write(appen,'(i1)')unities
     cs=blank//token(1:toklen)//'?'//trim(appen)//blank
     cslen=toklen+4
!    Map token to all upper case (make case-insensitive):
     call inupper(cs)
!    Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))
!    Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if(itoken/=itoken2)then
       write(message, '(7a)' )&
&       'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
&       'This is confusing, so it has been forbidden.',ch10,&
&       'Action: remove one of the two occurences.'
       MSG_ERROR(message)
     end if
     if(itoken/=0) then
       opttoken=1
       ds_input_=jdtset
     end if

!    Use the metacharacter for the unities, and save in cs1 and itoken1
     write(appen,'(i1)')dozens
     cs1=blank//token(1:toklen)//trim(appen)//'?'//blank
!    Map token to all upper case (make case-insensitive):
     call inupper(cs1)
!    Absolute index of blank//token//blank in string:
     itoken1=index(string,cs1(1:cslen))
!    Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs1(1:cslen), BACK=.true. )
     if(itoken1/=itoken2)then
       write(message, '(7a)' )&
&       'There are two occurences of the keyword "',cs1(1:cslen),'" in the input file.',ch10,&
&       'This is confusing, so it has been forbidden.',ch10,&
&       'Action: remove one of the two occurences.'
       MSG_ERROR(message)
     end if

     if(itoken/=0 .and. itoken1/=0)then
       write(message, '(9a)' )&
&       'The keywords "',cs(1:cslen),'" and "',cs1(1:cslen),'"',ch10,&
&       'cannot be used together in the input file.',ch10,&
&       'Action: remove one of the two keywords.'
       MSG_ERROR(message)
     end if

     if(itoken1/=0)then
       opttoken=1
       itoken=itoken1
       cs=cs1
       ds_input_=jdtset
     end if

   end if

!  --------------------------------------------------------------------------
!  (2b) try to find the tokens defining a series
   if(opttoken==0)then

     cs=token(1:toklen)

     cslen=toklen+3
     cs1len=toklen+4

     cscolon=blank//token(1:toklen)//':'//blank
     csplus=blank//token(1:toklen)//'+'//blank
     cstimes=blank//token(1:toklen)//'*'//blank

     cs1colon=blank//token(1:toklen)//'?'//':'//blank
     cs1plus=blank//token(1:toklen)//'?'//'+'//blank
     cs1times=blank//token(1:toklen)//'?'//'*'//blank

     cs2colon=blank//token(1:toklen)//':'//'?'//blank
     cs2plus=blank//token(1:toklen)//'+'//'?'//blank
     cs2times=blank//token(1:toklen)//'*'//'?'//blank

!    Map token to all upper case (make case-insensitive):
     call inupper(cscolon)
     call inupper(csplus)
     call inupper(cstimes)
     call inupper(cs1colon)
     call inupper(cs1plus)
     call inupper(cs1times)
     call inupper(cs2colon)
     call inupper(cs2plus)
     call inupper(cs2times)

!    Absolute index of tokens in string:
     itoken_colon=index(string,cscolon(1:cslen))
     itoken_plus=index(string,csplus(1:cslen))
     itoken_times=index(string,cstimes(1:cslen))
     itoken_1colon=index(string,cs1colon(1:cs1len))
     itoken_1plus=index(string,cs1plus(1:cs1len))
     itoken_1times=index(string,cs1times(1:cs1len))
     itoken_2colon=index(string,cs2colon(1:cs1len))
     itoken_2plus=index(string,cs2plus(1:cs1len))
     itoken_2times=index(string,cs2times(1:cs1len))

!    Look for another occurence of the same tokens in string
     itoken2_colon=index(string,cscolon(1:cslen), BACK=.true. )
     itoken2_plus=index(string,csplus(1:cslen), BACK=.true. )
     itoken2_times=index(string,cstimes(1:cslen), BACK=.true. )
     itoken2_1colon=index(string,cs1colon(1:cs1len), BACK=.true. )
     itoken2_1plus=index(string,cs1plus(1:cs1len), BACK=.true. )
     itoken2_1times=index(string,cs1times(1:cs1len), BACK=.true. )
     itoken2_2colon=index(string,cs2colon(1:cs1len), BACK=.true. )
     itoken2_2plus=index(string,cs2plus(1:cs1len), BACK=.true. )
     itoken2_2times=index(string,cs2times(1:cs1len), BACK=.true. )

     if(jdtset==0)then

!      If the multi-dataset mode is not used, no token should have been found
       if(itoken_colon+itoken_plus+itoken_times+&
&       itoken_2colon+itoken_2plus+itoken_2times > 0 ) then
         write(message,'(a,a,a,a,a,a,a,a,a,a,a,a, a)' )&
&         'Although the multi-dataset mode is not activated,',ch10,&
&         'the keyword "',trim(cs),'" has been found',ch10,&
&         'appended with  + * or :  .',ch10,&
&         'This is not allowed.',ch10,&
&         'Action: remove the appended keyword, or',ch10,&
&         'use the multi-dataset mode (ndtset/=0).'
         MSG_ERROR(message)
       end if
       if(itoken_1colon+itoken_1plus+itoken_1times > 0 ) then
         write(message, '(a,a,a,a,a,a,a,a,a,a,a,a,a)' )&
&         'Although the multi-dataset mode is not activated,',ch10,&
&         'the keyword "',trim(cs),'" has been found',ch10,&
&         'appended with ? , then + * or :  .',ch10,&
&         'This is not allowed.',ch10,&
&         'Action: remove the appended keyword, or',ch10,&
&         'use the multi-dataset mode (ndtset/=0).'
         MSG_ERROR(message)
       end if

     else

!      If the multi-dataset mode is used, exactly zero or two token must be found
       sum_token=0
       if(itoken_colon/=0)sum_token=sum_token+1
       if(itoken_plus /=0)sum_token=sum_token+1
       if(itoken_times/=0)sum_token=sum_token+1
       if(itoken_1colon/=0)sum_token=sum_token+1
       if(itoken_1plus /=0)sum_token=sum_token+1
       if(itoken_1times/=0)sum_token=sum_token+1
       if(itoken_2colon/=0)sum_token=sum_token+1
       if(itoken_2plus /=0)sum_token=sum_token+1
       if(itoken_2times/=0)sum_token=sum_token+1

       if(sum_token/=0 .and. sum_token/=2) then
         write(message, '(a,a,a,a,a,i3,a,a,a,a,a,a,a)' )&
&         'The keyword "',trim(cs),'" has been found to take part',ch10,&
&         'to series definition in the multi-dataset mode',sum_token,' times.',ch10,&
&         'This is not allowed, since it should be used once with ":",',ch10,&
&         'and once with "+" or "*".',ch10,&
&         'Action: change the number of occurences of this keyword.'
         MSG_ERROR(message)
       end if

!      If the multi-dataset mode is used, make sure that
!      no twice the same combined keyword happens
       ier=0
       if(itoken_colon/=itoken2_colon)then
         ier=1 ; cs=cscolon
       end if
       if(itoken_plus/=itoken2_plus)then
         ier=1 ; cs=csplus
       end if
       if(itoken_times/=itoken2_times)then
         ier=1 ; cs=cstimes
       end if
       if(itoken_1colon/=itoken2_1colon)then
         ier=1 ; cs=cs1colon
       end if
       if(itoken_1plus/=itoken2_1plus)then
         ier=1 ; cs=cs1plus
       end if
       if(itoken_1times/=itoken2_1times)then
         ier=1 ; cs=cs1times
       end if
       if(itoken_2colon/=itoken2_2colon)then
         ier=1 ; cs=cs2colon
       end if
       if(itoken_2plus/=itoken2_2plus)then
         ier=1 ; cs=cs2plus
       end if
       if(itoken_2times/=itoken2_2times)then
         ier=1 ; cs=cs2times
       end if
       if(ier==1)then
         write(message, '(a,a,a,a,a,a,a)' )&
&         'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
&         'This is confusing, so it has been forbidden.',ch10,&
&         'Action: remove one of the two occurences.'
         MSG_ERROR(message)
       end if

!      Select the series according to the presence of a colon flag
       if(itoken_colon>0)then
         opttoken=2
         ds_input_=jdtset
       else if(itoken_1colon>0)then
         opttoken=3
         cscolon=cs1colon ; csplus=cs1plus ; cstimes=cs1times
         itoken_colon=itoken_1colon
         itoken_plus=itoken_1plus ; itoken_times=itoken_1times
         cslen=cs1len
         ds_input_=jdtset
       else if(itoken_2colon>0)then
         opttoken=4
         cscolon=cs2colon ; csplus=cs2plus ; cstimes=cs2times
         itoken_colon=itoken_2colon
         itoken_plus=itoken_2plus ; itoken_times=itoken_2times
         cslen=cs1len
         ds_input_=jdtset
       end if

!      Make sure that the proper combination of : + and * is found .
       if(itoken_colon > 0 .and. (itoken_plus==0 .and. itoken_times==0) )then
         write(message, '(13a)' )&
&         'The keyword "',cscolon(1:cslen),'" initiate a series,',ch10,&
&         'but there is no occurence of "',csplus(1:cslen),'" or "',cstimes(1:cslen),'".',ch10,&
&         'Action: either suppress the series, or make the increment',ch10,&
&         'or the factor available.'
         MSG_ERROR(message)
       end if
       if(itoken_plus/=0 .and. itoken_times/=0)then
         write(message, '(a,a, a,a,a,a,a)' )&
&         'The combined occurence of keywords "',csplus(1:cslen),'" and "',cstimes(1:cslen),'" is not allowed.',ch10,&
&         'Action: suppress one of them in your input file.'
         MSG_ERROR(message)
       end if
       if(itoken_colon==0 .and. (itoken_plus/=0 .or. itoken_times/=0) ) then
         cs=csplus
         if(itoken_times/=0)cs=cstimes
         write(message, '(a,a,a,a,a,a,a,a,a,a,a)' )&
&         'The keyword "',cscolon(1:cslen),'" does not appear in the input file.',ch10,&
&         'However, the keyword "',cs(1:cslen),'" appears.',ch10,&
&         'This is forbidden.',ch10,&
&         'Action: make the first appear, or suppress the second.'
         MSG_ERROR(message)
       end if

!      At this stage, either
!      - itoken_colon vanish as well as itoken_plus and itoken_times
!      - itoken_colon does not vanish,
!      as well as one of itoken_plus or itoken_times

!      End the condition of multi-dataset mode
     end if

!    End the check on existence of a series
   end if

!  --------------------------------------------------------------------------
!  (3) if not found, try to find the token with non-modified string
   if(opttoken==0)then

     cs=blank//token(1:toklen)//blank
     cslen=toklen+2

!    Map token to all upper case (make case-insensitive):
     call inupper(cs)

!    Absolute index of blank//token//blank in string:
     itoken=index(string,cs(1:cslen))

!    Look for another occurence of the same token in string, if so, leaves:
     itoken2=index(string,cs(1:cslen), BACK=.true. )
     if(itoken/=itoken2)then
       write(message, '(a,a,a,a,a,a,a)' )&
&       'There are two occurences of the keyword "',cs(1:cslen),'" in the input file.',ch10,&
&       'This is confusing, so it has been forbidden.',ch10,&
&       'Action: remove one of the two occurences.'
       MSG_ERROR(message)
     end if

     if(itoken/=0) then
       opttoken=1
       ds_input_=0
     end if

   end if

!  --------------------------------------------------------------------------
!  If jdtset==0, means that the multi-dataset mode is not used, so
!  checks whether the input file contains a multi-dataset keyword,
!  and if this occurs, stop. Check also the forbidden occurence of
!  use of 0 as a multi-dataset index.
!  Note that the occurence of series initiators has already been checked.

   do trial_jdtset=0,9
     if(jdtset==0 .or. trial_jdtset==0)then
       write(appen,'(i1)')trial_jdtset
       trial_cs=blank//token(1:toklen)//trim(appen)
       trial_cslen=toklen+2
!      Map token to all upper case (make case-insensitive):
       call inupper(trial_cs)
!      Look for an occurence of this token in string, if so, leaves:
       itoken2=index(string,trial_cs(1:trial_cslen))
!      If itoken2/=0
       if(itoken2/=0)then
         if(trial_jdtset==0)then
           write(message, '(a,a,a,a,a,a,a)' )&
&           'There is an occurence of the keyword "',trim(token),'" appended with 0 in the input file.',ch10,&
&           'This is forbidden.',ch10,&
&           'Action: remove this occurence.'
           call wrtout(std_out,message,'COLL')
         else
           write(message, '(a,a,a,a,a,i1,a,a,a,a,a)' )&
&           'In the input file, there is an occurence of the ',ch10,&
&           'keyword "',trim(token),'", appended with the digit "',trial_jdtset,'".',ch10,&
&           'This is forbidden when ndtset==0 .',ch10,&
&           'Action: remove this occurence, or change ndtset.'
           call wrtout(std_out,message,'COLL')
         end if
         MSG_ERROR(message)
       end if
     end if
   end do

 end if

!===========================================================================
!At this stage, the location of the keyword string is known, as well
!as its length. So, can read the data.
!Usual reading if opttoken==1 (need itoken).
!If opttoken>=2, the characteristics of a series must be read
!(need itoken_colon and either itoken_plus or itoken_times)

 tread = 0
 typevar='INT'
 if(typevarphys=='LOG')typevar='INT'
 if(typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI')typevar='DPR'
 if(typevarphys=='KEY')then
   if(opttoken>=2)then
     write(message, '(9a)' )&
&     'For the keyword "',cs(1:cslen),'", of KEY type,',ch10,&
&     'a series has been defined in the input file.',ch10,&
&     'This is forbidden.',ch10,&
&     'Action: check your input file.'
     MSG_ERROR(message)
   end if
   if(narr>=2)then
     write(message, '(9a)' )&
&     'For the keyword "',cs(1:cslen),'", of KEY type,',ch10,&
&     'the number of data requested is larger than 1.',ch10,&
&     'This is forbidden.',ch10,&
&     'Action: check your input file.'
     MSG_ERROR(message)
   end if
   typevar='KEY'
!  write(std_out,*)' intagm : will read cs=',trim(cs)
!  stop
 end if

!There is something to be read if opttoken>=1
 if(opttoken==1)then

!  DEBUG
!  write(std_out,*)' intagm : opttoken==1 , token has been found, will read '
!  ENDDEBUG

!  Absolute location in string of blank which follows token:
   b1=itoken+cslen-1

!  Read the array (or eventual scalar) that follows the blank
!  In case of typevarphys='KEY', the chain of character will be returned in cs.
   call inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)

   if(typevarphys=='KEY')then
     if (.not. PRESENT(key_value)) then
       MSG_ERROR("typevarphys == KEY requires the optional argument key_value")
     end if
     !token=trim(cs)
     !write(std_out,*)' intagm : after inarray, token=',trim(token)
     key_value = TRIM(cs)
   end if

!  if this point is reached then data has been read in successfully
   tread = 1

 else if(opttoken>=2)then

!  write(std_out,*)' intagm : opttoken>=2 , token has been found, will read '
   ABI_ALLOCATE(dpr1,(narr))
   ABI_ALLOCATE(dpr2,(narr))
   ABI_ALLOCATE(int1,(narr))
   ABI_ALLOCATE(int2,(narr))

!  Absolute location in string of blank which follows token//':':
   b1=itoken_colon+cslen-1
   call inarray(b1,cscolon,dpr1,int1,narr,narr,string,typevarphys)

!  Initialise number even if the if series treat all cases.
   number=1
!  Define the number of the term in the series
   if(opttoken==2)number=jdtset-1
   if(opttoken==3)number=unities-1
   if(opttoken==4)number=dozens-1

!  Distinguish additive and multiplicative series
   if(itoken_plus/=0)then

     b1=itoken_plus+cslen-1
     call inarray(b1,csplus,dpr2,int2,narr,narr,string,typevarphys)

     if(typevar=='INT')then
       intarr(1:narr)=int1(:)+int2(:)*number
     else if(typevar=='DPR')then
       dprarr(1:narr)=dpr1(:)+dpr2(:)*number
     end if

   else if(itoken_times/=0)then

     b1=itoken_times+cslen-1
     call inarray(b1,cstimes,dpr2,int2,narr,narr,string,typevarphys)
     if(typevar=='INT')then
       intarr(1:narr)=int1(:)*int2(:)**number
     else if(typevar=='DPR')then
       dprarr(1:narr)=dpr1(:)*dpr2(:)**number
     end if

   end if

   tread = 1

   ABI_DEALLOCATE(dpr1)
   ABI_DEALLOCATE(dpr2)
   ABI_DEALLOCATE(int1)
   ABI_DEALLOCATE(int2)
 end if

 if(present(ds_input)) then
   ds_input = ds_input_
 end if

!DEBUG
!write(std_out,*) ' intagm : exit value tread=',tread
!write(std_out,*) ' intarr =',intarr(1:narr)
!write(std_out,*) ' dprarr =',dprarr(1:narr)
!stop
!ENDDEBUG

end subroutine intagm
!!***

!!****f* m_parser/inarray
!! NAME
!! inarray
!!
!! FUNCTION
!! Read the array of narr numbers located immediately after
!! a specified blank in a string of character.
!! Might read instead one word, after the specified blank.
!! Takes care of multipliers.
!!
!! INPUTS
!!  cs=character token
!!  marr=dimension of the intarr and dprarr arrays, as declared in the
!!   calling subroutine.
!!  narr=actual size of array to be read in  (if typevarphys='KEY', only narr=1 is allowed)
!!  string=character string containing the data.
!!  typevarphys=variable type (might indicate the physical meaning of
!!   for dimensionality purposes)
!!   'INT'=>integer
!!   'DPR'=>real(dp) (no special treatment)
!!   'LEN'=>real(dp) (expect a "length", identify bohr, au or angstrom,
!!       and return in au -atomic units=bohr- )
!!   'ENE'=>real(dp) (expect a "energy", identify Ha, hartree, eV, Ry, Rydberg)
!!   'BFI'=>real(dp) (expect a "magnetic field", identify T, Tesla)
!!   'LOG'=>integer, but read logical variable T,F,.true., or .false.
!!   'KEY'=>character, returned in token cs
!!
!! OUTPUT
!!  intarr(1:narr), dprarr(1:narr)
!!   integer or real(dp) arrays, respectively,
!!   into which data is read. Use these arrays even for scalars.
!!  errcod: if /=0, then something went wrong in subroutine "inread"
!!
!! SIDE EFFECT
!!   b1=absolute location in string of blank which follows the token
!!             (will be modified in the execution)
!!   cs=at input  : character token
!!      at output : chain of character replacing the token (only if typevarphys='KEY')
!!
!! PARENTS
!!      intagm
!!
!! CHILDREN
!!      inread,wrtout
!!
!! SOURCE

subroutine inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inarray'
 use interfaces_14_hidewrite
 use interfaces_42_parser, except_this_one => inarray
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: marr,narr
 integer,intent(inout) :: b1
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: typevarphys
 character(len=fnlen),intent(inout) :: cs
!arrays
 integer,intent(inout) :: intarr(marr) !vz_i
 real(dp),intent(out) :: dprarr(marr)

!Local variables-------------------------------
 character(len=1), parameter :: blank=' '
!scalars
 integer :: asciichar,b2,errcod,ii,integ,istar,nrep,strln
 real(dp) :: factor,real8
 character(len=3) :: typevar
 character(len=500) :: message

! *************************************************************************

!DEBUG
!write(std_out,'(2a)' )' inarray : token=',trim(cs)
!write(std_out,'(a,i4)' )' inarray : narr=',narr
!write(std_out,'(2a)' )' inarray : typevarphys=',typevarphys
!ENDDEBUG

 ii=0
 typevar='INT'
 if(typevarphys=='LOG')typevar='INT'
 if(typevarphys=='DPR' .or. typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI')typevar='DPR'
 strln=len_trim(string)

 do while (ii<narr)

!  Relative location of next blank after data
   if(b1>=strln)exit   ! b1 is the last character of the string
   b2=index(string(b1+1:),blank)
!  If no second blank is found put the second blank just beyond strln
   if(b2==0) b2=strln-b1+1

   if(typevarphys=='KEY')then
     cs=string(b1+1:b1+b2-1)
     errcod=0
     exit
   end if

!  nrep tells how many times to repeat input in array:
   nrep=1

!  Check for *, meaning repeated input (as in list-directed input):
   istar=index(string(b1+1:b1+b2-1),'*')
   if (istar/=0) then
     if (istar==1) then ! Simply fills the array with the data, repeated as many times as needed
       nrep=narr-ii
       errcod=0
     else
       call inread(string(b1+1:b1+istar-1),istar-1,'INT',nrep,real8,errcod)
     end if
     if (errcod/=0) exit
!    Shift starting position of input field:
     b1=b1+istar
     b2=b2-istar
   end if

!  Read data internally by calling inread at entry ini:
   call inread(string(b1+1:b1+b2-1),b2-1,typevarphys,integ,real8,errcod)
   if (errcod/=0) exit

!  Allow for list-directed input with repeat number nrep:
   if(typevar=='INT')then
     intarr(1+ii:min(nrep+ii,narr))=integ
   else if(typevar=='DPR')then
     dprarr(1+ii:min(nrep+ii,narr))=real8
   else
     MSG_BUG('Disallowed typevar='//typevar)
   end if
   ii=min(ii+nrep,narr)

!  Find new absolute location of next element of array:
   b1=b1+b2

!  End do while (ii<narr). Note "exit" instructions within loop.
 end do

!if (ii>narr) then
!write(message, '(a,a,a,a,a,a,a,a,a,a,i4,a,i4,a,a,a,a,a,a,a,a)' ) ch10,&
!' inarray : ERROR -',ch10,&
!&  '  Too many data are provided in the input file for',ch10,&
!&  '  the keyword "',cs,'" :',ch10,&
!&  '  attempted to read',ii,' elements for array length',narr,ch10,&
!&  '  This might be due to an erroneous value for the size ',ch10,&
!&  '  of this array, in the input file.',ch10,&
!&  '  Action : check the data provided for this keyword,',ch10,&
!&  '  as well as its declared dimension. They do not match.'
!call wrtout(std_out,message,'COLL')
!end if

 if(errcod/=0)then

   write(message, '(8a,i0,a)' ) ch10,&
&   ' inarray : ',ch10,&
&   '  An error occurred reading data for keyword "',trim(cs),'",',ch10,&
&   '  looking for ',narr,' array elements.'
   call wrtout(std_out,message,do_flush=.true.)

   write(message,'(8a)')&
&   'There is a problem with the input file : maybe  ',ch10,&
&   'a disagreement between the declared dimension of the array,',ch10,&
&   'and the number of data actually provided. ',ch10,&
&   'Action: correct your input file, and especially the keywork', trim(cs)
   MSG_ERROR(message)
 end if

!In case of 'LEN', 'ENE', or 'BFI', try to identify the unit
 if(typevarphys=='LEN' .or. typevarphys=='ENE' .or. typevarphys=='BFI')then
   do

!    Relative location of next blank after data
     if(b1>=strln)exit   ! b1 is the last character of the string
     b2=index(string(b1+1:),blank)
!    If no second blank is found put the second blank just beyond strln
     if(b2==0) b2=strln-b1+1

!    DEBUG
!    write(std_out,*)' inarray : string(b1+1:)=',string(b1+1:)
!    write(std_out,*)' inarray : b2=',b2
!    write(std_out,*)' typevarphys==',typevarphys
!    ENDDEBUG

!    Identify the presence of a non-digit character
     asciichar=iachar(string(b1+1:b1+1))
     if(asciichar<48 .or. asciichar>57)then
       factor=one
       if(typevarphys=='LEN' .and. b2>=7)then
         if(string(b1+1:b1+6)=='ANGSTR')then
           factor=one/Bohr_Ang
         end if
       else if(typevarphys=='ENE' .and. b2>=3)then
         if(string(b1+1:b1+3)=='RY ')then
           factor=half
         else if(string(b1+1:b1+3)=='EV ')then
           factor=one/Ha_eV
         end if
       else if(typevarphys=='ENE' .and. b2>=2)then
         if(string(b1+1:b1+2)=='K ')then
           factor=kb_HaK
         end if
       else if(typevarphys=='BFI' .and. b2>=2)then
         if(string(b1+1:b1+2)=='T ' .or. string(b1+1:b1+2)=='TE')then
           factor=BField_Tesla
         end if
       end if
       dprarr(1:narr)=dprarr(1:narr)*factor
       exit
     else
!      A digit has been observed, go to the next sequence
       b1=b2
       cycle
     end if

   end do
 end if

!DEBUG
!write(std_out,*)' inarray : exit '
!write(std_out,*)' dprarr(1:narr)==',dprarr(1:narr)
!ENDDEBUG

end subroutine inarray
!!***

!!****f* m_parser/importxyz
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

subroutine importxyz(lenstr,string_raw,string_upper,strln)

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

!!****f* m_parser/append_xyz
!! NAME
!! append_xyz
!!
!! FUNCTION
!! Translate the data from a xyz file (xyz_fname),
!! and add it at the end of the usual ABINIT input data string (string),
!! taking into account the dtset (dtset_char)
!!
!! INPUTS
!!  dtset_char*2=possible dtset label
!!  xyz_fname = name of the xyz file
!!  strln=maximal number of characters of string, as declared in the calling routine
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  lenstr=actual number of characters in string
!!  string*(strln)=string of characters  (upper case) to which the xyz data are appended
!!
!! PARENTS
!!      importxyz
!!
!! CHILDREN
!!      atomdata_from_symbol,wrtout
!!
!! SOURCE

subroutine append_xyz(dtset_char,lenstr,string,xyz_fname,strln)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'append_xyz'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: strln
 integer,intent(inout) :: lenstr
 character(len=2),intent(in) :: dtset_char
 character(len=fnlen),intent(in) :: xyz_fname
 character(len=strln),intent(inout) :: string

!Local variables-------------------------------
 character :: blank=' '
!scalars
 integer :: unitxyz, iatom, natom, mu
 integer :: lenstr_new
 integer :: lenstr_old
 integer :: ntypat
 real(dp) :: znucl
 character(len=5) :: string5
 character(len=20) :: string20
 character(len=500) :: message
 type(atomdata_t) :: atom
!arrays
 real(dp),allocatable :: xangst(:,:)
 integer, save :: atomspecies(200) = 0
 character(len=500), save :: znuclstring = ""
 character(len=2),allocatable :: elementtype(:)

!************************************************************************

 lenstr_new=lenstr

 if (dtset_char == "-1") then
!  write znucl
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+7+len_trim(znuclstring)+1
   string(lenstr_old+1:lenstr_new)=" ZNUCL"//blank//trim(znuclstring)//blank

!  write ntypat
   ntypat = sum(atomspecies)
   write(string20,'(i10)') ntypat
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+8+len_trim(string20)+1
   string(lenstr_old+1:lenstr_new)=" NTYPAT"//blank//trim(string20)//blank

   return
 end if

!open file with xyz data
 if (open_file(xyz_fname, message, newunit=unitxyz, status="unknown") /= 0) then
   MSG_ERROR(message)
 end if
 write(message, '(3a)')' importxyz : Opened file ',trim(xyz_fname),'; content stored in string_xyz'
 call wrtout(std_out,message,'COLL')

!check number of atoms is correct
 read(unitxyz,*) natom

 write(string5,'(i5)')natom
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1+5
 string(lenstr_old+1:lenstr_new)=" _NATOM"//trim(dtset_char)//blank//string5

 ABI_ALLOCATE(xangst,(3,natom))
 ABI_ALLOCATE(elementtype,(natom))

!read dummy line
 read(unitxyz,*)

!read atomic types and positions
 do iatom = 1, natom
   read(unitxyz,*) elementtype(iatom), xangst(:,iatom)
!  extract znucl for each atom type
   call atomdata_from_symbol(atom,elementtype(iatom))
   znucl = atom%znucl
   if (znucl > 200) then
     write (message,'(5a)')&
&     'Error: found element beyond Z=200 ', ch10,&
&     'Solution: increase size of atomspecies in append_xyz', ch10
     MSG_ERROR(message)
   end if
!  found a new atom type
   if (atomspecies(int(znucl)) == 0) then
     write(string20,'(f10.2)') znucl
     znuclstring = trim(znuclstring) // " " // trim(string20) // " "
   end if
   atomspecies(int(znucl)) = 1
 end do
 close (unitxyz)


!Write the element types
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+7+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _TYPAX"//trim(dtset_char)//blank
 do iatom=1,natom
   lenstr_old=lenstr_new
   lenstr_new=lenstr_new+3
   string(lenstr_old+1:lenstr_new)=elementtype(iatom)//blank
 end do
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+3
 string(lenstr_old+1:lenstr_new)="XX " ! end card for TYPAX

!Write the coordinates
 lenstr_old=lenstr_new
 lenstr_new=lenstr_new+8+len_trim(dtset_char)+1
 string(lenstr_old+1:lenstr_new)=" _XANGST"//trim(dtset_char)//blank

 do iatom=1,natom
   do mu=1,3
     write(string20,'(f20.12)')xangst(mu,iatom)
     lenstr_old=lenstr_new
     lenstr_new=lenstr_new+20
     string(lenstr_old+1:lenstr_new)=string20
   end do
 end do

 ABI_DEALLOCATE(elementtype)
 ABI_DEALLOCATE(xangst)

!Check the length of the string
 if(lenstr_new>strln)then
   write(message,'(3a)')&
&   'The maximal size of the input variable string has been exceeded.',ch10,&
&   'The use of a xyz file is more character-consuming than the usual input file. Sorry.'
   MSG_BUG(message)
 end if

!Update the length of the string
 lenstr=lenstr_new

end subroutine append_xyz
!!***

end module m_parser
!!***

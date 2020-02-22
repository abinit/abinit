!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_fstrings
!! NAME
!!  m_fstrings
!!
!! FUNCTION
!!  This module contains basic tools to operate on Fortran strings.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MG, XG, MT, DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
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

MODULE m_fstrings

 use defs_basis, only : dp, std_out, ch10

 implicit none

 private

 public :: is_letter       ! Returns .TRUE. if ch is a letter and .FALSE. otherwise
 public :: is_digit        ! Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise
 public :: upper           ! Convert lower case letters to UPPER CASE
 public :: toupper         ! Convert lower case letters to UPPER CASE (function version)
 public :: lower           ! Convert UPPER CASE letters to lower case
 public :: tolower         ! Convert UPPER CASE letters to lower case  (function version)
 public :: removesp        ! Removes spaces, tabs, and control characters in string str
 public :: replace_ch0     ! Replace final '\0' with whitespaces
 public :: lstrip          ! Remove leading spaces from string
 public :: replace         ! Replace chars in string.
 public :: ljust           ! Return a left-justified string of length width.
 public :: lpad            ! Pad a string adding repeat characters fillchar on the left side.
 public :: round_brackets  ! Return a new string enclosed in parentheses if not already present.
 public :: quote           ! Return a new string enclosed by quotation marks.
 public :: rmquotes        ! Remove quotation marks from a string. Return new string
 public :: write_num       ! Writes a number to a string using format fmt
 public :: trimzero        ! Deletes nonsignificant trailing zeroes from a number string.
 public :: writeq          ! Writes a string of the form <name> = value to unit
 public :: strcat          ! Concatenate strings (function version)
 public :: sjoin           ! Joins strings with a space separator.
 public :: yesno           ! Convert boolean to "yes", "no"
 public :: itoa            ! Convert an integer into a string
 public :: ftoa            ! Convert a float into a string
 public :: ktoa            ! Convert a k-point into a string.
 public :: ltoa            ! Convert a list into a string.
 public :: atoi            ! Convert a string into a integer
 public :: atof            ! Convert a string into a floating-point number.
 public :: basename        ! Returns the final component of a pathname.
 public :: firstchar       ! Returns .TRUE. is the first character in a string belongs to a gives set.
 public :: startswith      ! Returns .TRUE. is the string starts with the specified prefix.
 public :: endswith        ! Returns .True if the string ends with the specified suffix.
 public :: indent          ! Indent text
 public :: prep_dash       ! Prepend `-` to each line in a string.
 public :: int2char4       ! Convert a positive integer number (zero included) to a character(len=*)
                           ! with trailing zeros if the number is <=9999
 public :: int2char10      ! Convert a positive integer number (zero included) to a character(len=10)
                           ! with trailing blanks
 public :: char_count      ! Count the occurrences of a character in a string.
 public :: next_token      ! Tokenize a string made of whitespace-separated tokens.
 public :: inupper         ! Maps all characters in string to uppercase except for tokens between quotation marks.

 !TODO method to center a string

 interface write_num
   module procedure write_rdp_0D
   module procedure write_int_0D
 end interface write_num

 interface writeq
   module procedure writeq_rdp_0D
   module procedure writeq_int_0D
 end interface writeq

 interface is_digit
   module procedure is_digit_0D
 end interface is_digit

 interface firstchar
   module procedure firstchar_0d
   module procedure firstchar_1d
 end interface firstchar

 interface sjoin
   module procedure sjoin_2
   module procedure sjoin_3
   module procedure sjoin_4
   module procedure sjoin_5
   module procedure sjoin_6
   module procedure sjoin_7
 end interface sjoin

 interface strcat
   module procedure strcat_2
   module procedure strcat_3
   module procedure strcat_4
   module procedure strcat_5
 end interface strcat

 interface ltoa
   module procedure ltoa_int
   module procedure ltoa_dp
 end interface ltoa

 character(len=1),parameter :: BLANK=' '
 character(len=1),parameter :: NCHAR = char(10)
 character(len=1),parameter :: DIR_SEPARATOR = '/'

 integer,parameter :: ASCII_A=ICHAR('A')
 integer,parameter :: ASCII_Z=ICHAR('Z')
 integer,parameter :: ASCII_aa=ICHAR('a')
 integer,parameter :: ASCII_zz=ICHAR('z')
 integer,parameter :: SHIFT=ASCII_aa-ASCII_A ! Capital letters have smaller Dec value in the ASCII table.
 integer,parameter :: ASCII_0=ICHAR('0')
 integer,parameter :: ASCII_9=ICHAR('9')

 integer,parameter :: MAX_SLEN = 500


CONTAINS  !===========================================================
!!***

!!****f* m_fstrings/is_letter
!! NAME
!!  is_letter
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a letter and .FALSE. otherwise.
!!
!! SOURCE

pure function is_letter(ch) result(ans)

 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case (ASCII_A:ASCII_Z,ASCII_aa:ASCII_zz)
   ans=.TRUE.
 case DEFAULT
   ans=.FALSE.
 end select

end function is_letter
!!***

!!****f* m_fstrings/is_digit_0D
!! NAME
!!  is_digit_0D
!!
!! FUNCTION
!!  Returns .TRUE. if ch is a digit (0,1,...,9) and .FALSE. otherwise.
!!
!! SOURCE

pure function is_digit_0D(ch) result(ans)

!Arguments ------------------------------------
 character(len=1),intent(in) :: ch
 logical :: ans
! *********************************************************************

 select case (ICHAR(ch))
 case(ASCII_0:ASCII_9)
   ans=.TRUE.
 case default
   ans=.FALSE.
 end select

end function is_digit_0D
!!***

!!****f* m_fstrings/upper
!! NAME
!!  upper
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE.
!!
!! SOURCE

pure subroutine upper(str)

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str)
   iasc=IACHAR(str(ic:ic))
   if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) str(ic:ic)=ACHAR(iasc-SHIFT)
 end do

end subroutine upper
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/toupper
!! NAME
!!  toupper
!!
!! FUNCTION
!!  Convert lower case letters to UPPER CASE (function version).
!!
!! SOURCE

pure function toupper(str_in) result(str_out)

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str_in)
   iasc=IACHAR(str_in(ic:ic))
   if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) then
     str_out(ic:ic)=ACHAR(iasc-SHIFT)
   else
     str_out(ic:ic)=str_in(ic:ic)
   end if
 end do

end function toupper
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/lower
!! NAME
!!  lower
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case.
!!
!! PARENTS
!!      fftprof,ioprof,lapackprof
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine lower(str)

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str)
   iasc=IACHAR(str(ic:ic))
   if (iasc>=ASCII_A.and.iasc<=ASCII_Z) str(ic:ic)=ACHAR(iasc+SHIFT)
 end do

end subroutine lower
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/tolower
!! NAME
!!  tolower
!!
!! FUNCTION
!!  Convert UPPER CASE letters to lower case (function version).
!!
!! SOURCE

pure function tolower(str_in) result(str_out)

 character(len=*),intent(in) :: str_in
 character(len=LEN_TRIM(str_in)) :: str_out

!Local variables-------------------------------
 integer :: ic,iasc
! *********************************************************************

 do ic=1,LEN_TRIM(str_in)
   iasc=IACHAR(str_in(ic:ic))
   if (iasc>=ASCII_A.and.iasc<=ASCII_Z) then
     str_out(ic:ic)=ACHAR(iasc+SHIFT)
   else
     str_out(ic:ic)=str_in(ic:ic)
   end if
 end do

end function tolower
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/removesp
!! NAME
!!  removesp
!!
!! FUNCTION
!!  Removes spaces, tabs, and control characters in string str.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE

subroutine removesp(str)

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,k,lenstr,ich
 character(len=1):: ch
 character(len=LEN_TRIM(str)):: outstr
! *********************************************************************

 str=ADJUSTL(str) ; lenstr=LEN_TRIM(str)

 outstr=BLANK ; k=0
 do i=1,lenstr
   ch=str(i:i)
   ich=IACHAR(ch)
   select case(ich)
   case(0:32)  ! space, tab, or control character
     CYCLE
   case(33:)
     k=k+1
     outstr(k:k)=ch
   end select
 end do

 str=ADJUSTL(outstr)

end subroutine removesp
!!***

!!****m* m_fstrings/replace_ch0
!! NAME
!!  replace_ch0
!!
!! FUNCTION
!!  Little tool to change all final '\0' (end of string in C) characters to ' ' (space).
!!
!! SIDE EFFECTS
!!  * string = the string to convert. It is done in-place.
!!
!! SOURCE

elemental subroutine replace_ch0(string)

  character(len=*), intent(inout) :: string

  integer :: i, l

  i = index(string, char(0))
  if (i > 0) then
     l = len(string)
     string(i:l) = repeat(" ", l - i + 1)
  end if

end subroutine replace_ch0
!!***

!!****m* m_fstrings/replace
!! NAME
!!  replace
!!
!! FUNCTION
!!  Replace `text` with `rep` in string `s`. Return new string.
!!
!! NOTES:
!!  The length of the output string is increased by 500 but this could not be enough
!!  if len_trim(text) > len_trim(re) and there are several occurrences of `text` in s.
!!
!! SOURCE

function replace(s, text, rep)  result(outs)

 character(len=*),intent(in) :: s,text,rep
 character(len(s)+500) :: outs     ! provide outs with extra 500 char len

!Local variables-------------------------------
 integer :: i, j, nt, nr, last
! *********************************************************************

 outs = s; nt = len_trim(text); nr = len_trim(rep); last = 1
 do
   i = index(outs(last:), text(1:nt)); if (i == 0) exit
   j = last + i - 1; last = j + nr
   if (j - 1 < 1) then
     outs = rep(:nr) // outs(j+nt:)
   else
     outs = outs(:j-1) // rep(:nr) // outs(j+nt:)
   end if
 end do

end function replace
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/lstrip
!! NAME
!!  lstrip
!!
!! FUNCTION
!!  Removes leading spaces from the input string.
!!
!! SOURCE

pure function lstrip(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=len(istr)) :: ostr

!Local variables-------------------------------
 integer :: ii,jj,lg
! *********************************************************************

 lg=LEN(istr)
 do ii=1,lg
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 ostr = " "
 do jj=1,lg-ii+1
   ostr(jj:jj) = istr(ii:ii)
   ii=ii+1
 end do

end function lstrip
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ljust
!! NAME
!!  ljust
!!
!! FUNCTION
!!  Return S left-justified in a string of length width. Padding is
!!  done using the specified fill character (default is a space).
!!
!! SOURCE

pure function ljust(istr, width, fillchar) result(ostr)

 character(len=*),intent(in) :: istr
 integer,intent(in) :: width
 character(len=width) :: ostr
 character(len=1),optional,intent(in) :: fillchar

!Local variables-------------------------------
 integer :: ii
! *********************************************************************

 ostr = ADJUSTL(istr)

 if (PRESENT(fillchar)) then
   do ii=LEN_TRIM(ostr)+1,width
     ostr(ii:ii) = fillchar
   end do
 end if

end function ljust
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/lpad
!! NAME
!!  lpad
!!
!! FUNCTION
!!  Pad a string adding repeat characters fillchar on the left side.
!!  Padding is done using the specified fill character (default is a blanck character).
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function lpad(istr, repeat, fillchar) result(ostr)

 character(len=*),intent(in) :: istr
 integer,intent(in) :: repeat
 character(len=LEN_TRIM(istr) + repeat) :: ostr
 character(len=1),optional,intent(in) :: fillchar

!Local variables-------------------------------
 integer :: ii
 character(len=1) :: ch
! *********************************************************************

 ostr(repeat+1:) = TRIM(istr)

 ch = " "; if (PRESENT(fillchar)) ch = fillchar
 do ii=1,repeat
   ostr(ii:ii) = ch
 end do

end function lpad
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/round_brackets
!! NAME
!!  round_brackets
!!
!! FUNCTION
!!  Return a new string enclosed in parentheses if not already present.
!!
!! SOURCE

pure function round_brackets(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)+2) :: ostr

!Local variables-------------------------------
 integer :: ii
 character(len=1) :: qq
 character(len=LEN(istr)+2) :: tmp

! *********************************************************************

 do ii=1,LEN(istr)
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 qq = istr(ii:ii)

 if (qq == "(") then
   ! Don't add quotation marks if they already present.
   tmp = istr
   ii = LEN_TRIM(tmp)
   ! If the string is not closed, fix it.
   if (tmp(ii:ii) /= ")") tmp(ii+1:ii+1) = ")"
   ostr = TRIM(tmp)

 else
   qq = '('
   ostr(1:1) = qq
   ostr(2:) = TRIM(istr)
   ii = LEN_TRIM(ostr)+1
   ostr(ii:ii) = ")"
 end if

end function round_brackets
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/quote
!! NAME
!!  quote
!!
!! FUNCTION
!!  Return a new string enclosed by quotation marks.
!!
!! SOURCE

pure function quote(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)+2) :: ostr

!Local variables-------------------------------
 integer :: ii
 character(len=1) :: qq
 character(len=LEN(istr)+2) :: tmp

! *********************************************************************

 do ii=1,LEN(istr)
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 qq = istr(ii:ii)

 if (qq == "'" .or. qq == '"') then
   ! Don't add quotation marks if they already present.
   tmp = istr
   ii = LEN_TRIM(tmp)
   ! If the string is not closed, fix it.
   if (tmp(ii:ii) /= qq) tmp(ii+1:ii+1) = qq
   ostr = TRIM(tmp)

 else
   qq = '"'
   ostr(1:1) = qq
   ostr(2:) = TRIM(istr)
   ii = LEN_TRIM(ostr)+1
   ostr(ii:ii) = qq
 end if

end function quote
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/rmquotes
!! NAME
!!  rmquotes
!!
!! FUNCTION
!!  Remove quotation marks from a string. Return new string
!!
!! SOURCE

pure function rmquotes(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=len(istr)) :: ostr

!Local variables-------------------------------
 integer :: ii,cnt

! *********************************************************************

 ostr = ""; cnt = 0
 do ii=1,len_trim(istr)
   if (any(istr(ii:ii) == ["'", '"'])) cycle
   cnt = cnt + 1
   ostr(cnt:cnt) = istr(ii:ii)
 end do

end function rmquotes
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/write_rdp_0d
!! NAME
!!  write_rdp_0d
!!
!! FUNCTION
!!  Writes a number to a string using format fmt.
!!
!! PARENTS
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE

subroutine write_rdp_0d(rnum,str,fmt)

!Arguments ------------------------------------
 real(dp),intent(in) :: rnum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt)rnum
 str=ADJUSTL(str)

end subroutine write_rdp_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/write_int_0d
!! NAME
!!  write_int_0d
!!
!! FUNCTION
!!  Writes a number to a string using format fmt.
!!
!! PARENTS
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE

subroutine write_int_0D(inum,str,fmt)

!Arguments ------------------------------------
 integer,intent(in) :: inum
 character(len=*),intent(in) :: fmt
 character(len=*),intent(out) :: str

!Local variables-------------------------------
 character(len=LEN(fmt)+2) :: formt
! *********************************************************************

 formt='('//TRIM(fmt)//')'
 write(str,formt) inum
 str=ADJUSTL(str)

end subroutine write_int_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/trimzero
!! NAME
!!  trimzero
!!
!! FUNCTION
!! Deletes nonsignificant trailing zeroes from number string str. If number
!! string ends in a decimal point, one trailing zero is added.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!      m_fstrings
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE
! NOT sure it will work

subroutine trimzero(str)

 character(len=*),intent(inout) :: str

!Local variables-------------------------------
 integer :: i,ipos,lstr
 character :: ch
 character(len=10) :: sexp
! *********************************************************************

 ipos=SCAN(str,'eE')
 if (ipos>0) then
  sexp=str(ipos:)
  str=str(1:ipos-1)
 end if
 lstr=LEN_TRIM(str)
 do i=lstr,1,-1
  ch=str(i:i)
  if (ch=='0') CYCLE
  if (ch=='.') then
   str=str(1:i)//'0'
   if (ipos>0) str=TRIM(str)//TRIM(sexp)
   EXIT
  end if
  str=str(1:i)
  EXIT
 end do

 if (ipos>0) str=TRIM(str)//TRIM(sexp)

end subroutine trimzero
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/writeq_rdp_0D
!! NAME
!!  writeq_rdp_0D
!!
!! FUNCTION
!!  Writes a string of the form <name> = value to unit.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE
subroutine writeq_rdp_0D(unit,namestr,value,fmt)

 real(dp),intent(in) :: value
 integer,intent(in) :: unit
 character(len=*),intent(in) :: fmt
 character(len=*),intent(in) :: namestr

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(value,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_rdp_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/writeq_int_0D
!! NAME
!!  writeq_int_0D
!!
!! FUNCTION
!!  Writes a string of the form <name> = value to unit.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      trimzero,write_num
!!
!! SOURCE

subroutine writeq_int_0D(unit,namestr,ivalue,fmt)

 integer,intent(in) :: ivalue
 integer,intent(in) :: unit
 character(len=*),intent(in) :: namestr
 character(len=*),intent(in) :: fmt

!Local variables-------------------------------
 character(len=32) :: tempstr
! *********************************************************************

 call write_num(ivalue,tempstr,fmt)
 call trimzero(tempstr)
 write(unit,*)TRIM(namestr)//' = '//TRIM(tempstr)

end subroutine writeq_int_0D
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_2
!! NAME
!! sjoin_2
!!
!! FUNCTION
!!  Joins two strings with a space separator except if first string is empty.
!!

pure function sjoin_2(str1,str2) result(ostr)

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+1) :: ostr

! *********************************************************************

 if (len_trim(str1) > 0) then
   ostr=TRIM(str1)//" "//TRIM(str2)
 else
   ostr=TRIM(str2)
 end if

end function sjoin_2
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_3
!! NAME
!! sjoin_3
!!
!! FUNCTION
!!  Joins three strings with a space separator.
!!

pure function sjoin_3(str1,str2,str3) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+2) :: ostr

! *********************************************************************

 ostr = sjoin_2(sjoin_2(str1, str2), str3)

end function sjoin_3
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_4
!! NAME
!! sjoin_4
!!
!! FUNCTION
!!  Joins four strings with a space separator.
!!

pure function sjoin_4(str1,str2,str3,str4) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+len_trim(str4)+3) :: ostr

! *********************************************************************

 ostr = sjoin_2(str1, sjoin_3(str2, str3, str4))

end function sjoin_4
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_5
!! NAME
!! sjoin_5
!!
!! FUNCTION
!!  Joins five strings with a space separator.
!!

pure function sjoin_5(str1,str2,str3,str4,str5) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4,str5
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+len_trim(str4)+len_trim(str5)+4) :: ostr

! *********************************************************************

 ostr = sjoin_2(str1, sjoin_4(str2, str3, str4, str5))

end function sjoin_5
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_6
!! NAME
!! sjoin_6
!!
!! FUNCTION
!!  Joins six strings with a space separator.
!!

pure function sjoin_6(str1,str2,str3,str4,str5,str6) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4,str5,str6
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+len_trim(str4)+len_trim(str5)+len_trim(str6)+5) :: ostr

! *********************************************************************

 ostr = sjoin_2(str1, sjoin_5(str2, str3, str4, str5, str6))

end function sjoin_6
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/sjoin_7
!! NAME
!! sjoin_7
!!
!! FUNCTION
!!  Joins six strings with a space separator.
!!

pure function sjoin_7(str1,str2,str3,str4,str5,str6,str7) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4,str5,str6,str7
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+len_trim(str4)+len_trim(str5)+len_trim(str6)+len(str7)+6) &
&  :: ostr

! *********************************************************************

 ostr = sjoin_2(str1, sjoin_6(str2, str3, str4, str5, str6, str7))

end function sjoin_7
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_2
!! NAME
!! strcat_2
!!
!! FUNCTION
!!  Returns two concatenated strings.
!!

pure function strcat_2(str1,str2) result(ostr)

 character(len=*),intent(in) :: str1,str2
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)) :: ostr

! *********************************************************************

 ostr=TRIM(str1)//TRIM(str2)

end function strcat_2
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_3
!! NAME
!! strcat_3
!!
!! FUNCTION
!!  Concatenate 3 strings
!!

pure function strcat_3(str1, str2, str3) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)) :: ostr

! *********************************************************************

 ostr = TRIM(str1)//TRIM(str2)//TRIM(str3)

end function strcat_3
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_4
!! NAME
!! strcat_3
!!
!! FUNCTION
!!  Concatenate 4 strings
!!

pure function strcat_4(str1, str2, str3, str4) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+LEN_TRIM(str4)) :: ostr

! *********************************************************************

 ostr = TRIM(str1)//TRIM(str2)//TRIM(str3)//TRIM(str4)

end function strcat_4
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/strcat_5
!! NAME
!! strcat_5
!!
!! FUNCTION
!!  Concatenate 5 strings
!!

pure function strcat_5(str1, str2, str3, str4, str5) result(ostr)

 character(len=*),intent(in) :: str1,str2,str3,str4,str5
 character(len=LEN_TRIM(str1)+LEN_TRIM(str2)+LEN_TRIM(str3)+LEN_TRIM(str4)+LEN_TRIM(str5)) :: ostr

! *********************************************************************

 ostr = TRIM(str1)//TRIM(str2)//TRIM(str3)//TRIM(str4)//trim(str5)

end function strcat_5
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/yesno
!! NAME
!! yesno
!!
!! FUNCTION
!!  Convert boolean into "yes" or "no"
!!

character(len=3) pure function yesno(bool)

!Arguments ------------------------------------
 logical,intent(in) :: bool

! *********************************************************************

 if (bool) then
   yesno = "yes"
 else
   yesno = "no"
 end if

end function yesno
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/atoi
!! NAME
!! atoi
!!
!! FUNCTION
!!  Convert a string into a integer
!!

integer function atoi(string)

!Arguments ------------------------------------
 character(len=*),intent(in) :: string

! *********************************************************************

 read(string,*,err=10)atoi
 return
 10 write(std_out,*)"Error while trying to convert string to integer. string: ",trim(string)

end function atoi
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/atof
!! NAME
!! atof
!!
!! FUNCTION
!!  Convert a string into a floating-point number
!!

real(dp) function atof(string)

!Arguments ------------------------------------
 character(len=*),intent(in) :: string

! *********************************************************************

 read(string,*,err=10)atof
 return
 10 write(std_out,*)"Error while trying to convert string to floating-point. string: ",trim(string)

end function atof
!!***

!!****f* m_fstrings/itoa
!! NAME
!! itoa
!!
!! FUNCTION
!!  Convert an integer into a string
!!
!! INPUTS
!!  value=The integer
!!
!! PARENTS
!!
!! CHILDREN
!!

pure function itoa(value)

 integer,intent(in) :: value
 character(len=22) :: itoa

! *********************************************************************

 ! len=22 is large enough to contain integer*8
 write(itoa,"(i0)")value
 itoa = ADJUSTL(itoa)

end function itoa
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ftoa
!! NAME
!! ftoa
!!
!! FUNCTION
!!  Convert an float into a string using format fmt  (es16.6 if fmt is not given).
!!
!! PARENTS
!!
!! CHILDREN
!!

pure function ftoa(value,fmt)

 real(dp),intent(in) :: value
 character(len=*),optional,intent(in) :: fmt
 character(len=MAX_SLEN) :: ftoa

! *********************************************************************

 if (present(fmt)) then
   write(ftoa,round_brackets(fmt))value
 else
   write(ftoa,"(es16.6)")value
 end if
 ftoa = ADJUSTL(ftoa)

end function ftoa
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ktoa
!! NAME
!! ktoa
!!
!! FUNCTION
!!  Convert an k-point into a string using format fmt  (es.16.6 if fmt is not given).
!!
!! PARENTS
!!
!! CHILDREN
!!

pure function ktoa(kpt,fmt)

 real(dp),intent(in) :: kpt(3)
 character(len=*),optional,intent(in) :: fmt
 character(len=MAX_SLEN) :: ktoa

! *********************************************************************

 if (present(fmt)) then
   write(ktoa,fmt)kpt
 else
   write(ktoa,"(a,3(es11.4,a))")"[",kpt(1),", ",kpt(2),", ",kpt(3),"]"
 end if
 ktoa = ADJUSTL(ktoa)

end function ktoa
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ltoa_int
!! NAME
!! ltoa_int
!!
!! FUNCTION
!!  Convert a list of integers into a string.
!!
!! PARENTS
!!
!! CHILDREN

pure function ltoa_int(list) result(str)

 integer,intent(in) :: list(:)
 character(len=MAX_SLEN) :: str

!Local variables-------------------------------
 integer :: ii,base,sz
 character(len=MAX_SLEN) :: temp

! *********************************************************************

 sz = size(list)

 if (any(sz == [0, 1])) then
   if (sz == 0) str = "[]"
   if (sz == 1) write(str, "(a,i0,a)")"[",list(1),"]"
   return
 end if

 str = ""; base = 1
 do ii=1,sz

   ! Write to temp string and copy it to str if we have enough chars.
   ! Return if MAX_SLEN is too short.
   if (ii == 1) then
     write(temp, "(a,i0,a)")"[",list(1),", "
   else if (ii == sz) then
     write(temp, "(i0,a)")list(ii),"]"
   else
     write(temp, "(i0,a)")list(ii),", "
   end if

   if (base + len_trim(temp) - 1 <= MAX_SLEN) then
     str(base:) = trim(temp)
     base = len_trim(str) + 1
   else
     return
   end if
 end do

end function ltoa_int
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/ltoa_dp
!! NAME
!! ltoa_dp
!!
!! FUNCTION
!!  Convert a list of double precision numbers into a string.
!!  fmt specifies the format to be used ("es13.4" by default)
!!
!! PARENTS
!!
!! CHILDREN

pure function ltoa_dp(list, fmt) result(str)

 real(dp),intent(in) :: list(:)
 character(len=*),optional,intent(in) :: fmt
 character(len=MAX_SLEN) :: str

!Local variables-------------------------------
 integer :: ii,base,sz
 character(len=MAX_SLEN) :: temp,myfmt,fa

! *********************************************************************

 myfmt = "es13.4"; if (present(fmt)) myfmt = fmt
 sz = size(list)

 if (any(sz == [0, 1])) then
   if (sz == 0) str = "[]"
   if (sz == 1) write(str, sjoin("(a,",myfmt,",a)")) "[",list(1),"]"
   return
 end if

 str = ""; base = 1; fa = sjoin("(",myfmt,",a)")
 do ii=1,sz

   ! Write to temp string and copy it to str if we have enough chars.
   ! Return if MAX_SLEN is too short.
   if (ii == 1) then
     write(temp, sjoin("(a,",myfmt,",a)")) "[",list(1),","
   else if (ii == sz) then
     write(temp, fa)list(ii),"]"
   else
     write(temp, fa) list(ii),","
   end if

   if (base + len_trim(temp) - 1 <= MAX_SLEN) then
     str(base:) = trim(temp)
     base = len_trim(str) + 1
   else
     return
   end if
 end do

end function ltoa_dp
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/basename
!! NAME
!! basename
!!
!! FUNCTION
!!  Returns the final component of a pathname.
!!
!! INPUTS
!!  string=The input string
!!
!! NOTES
!!  * If the input string in not a valid path to a file (i.e not in the form foo/name)
!!    a blank strink is returned
!!  * We do a backward search becase we want to optimize the algorithm for Fortran strings.
!!
!! SOURCE

pure function basename(string)

 character(len=*),intent(in) :: string
 character(len=LEN_TRIM(string)) :: basename

!Local variables-------------------------------
 integer :: ic,nch_trim,nch
!************************************************************************

 nch     =LEN     (string)
 nch_trim=LEN_TRIM(string)

 ic = INDEX (TRIM(string), DIR_SEPARATOR, back=.TRUE.)
 !write(*,*)'DEBUG ',TRIM(string),ic

 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the separator.
  basename = string(ic+1:nch_trim)
  return
 else if (ic==0 .or. ic == nch_trim+1) then ! no separator in string or zero length string,
  basename = TRIM(string)                   ! return trimmed string.
  return
 else              ! (ic == nch_trim) separator is the last char.
  basename= BLANK  ! This is not a valid path to a file, return blank.
  return
 end if

end function basename
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/firstchar_0d
!! NAME
!! firstchar_0d
!!
!! FUNCTION
!!   Return True if string starts with the specified character
!!
!! INPUTS
!!  string=The string whose first character has to be cheched
!!  ch=Character
!!  [csens]=.TRUE. if comparison is done regardless of case. Defaults to .FALSE.
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

pure function firstchar_0d(string,ch,csens) result(ans)

 logical :: ans
 logical,optional,intent(in) :: csens
 character(len=*),intent(in) :: string
 character(len=1),intent(in) :: ch

!Local variables-------------------------------
 logical :: my_csens
!************************************************************************

 my_csens=.FALSE.; if (PRESENT(csens)) my_csens = csens

 if (.not.my_csens) then
   ans = ( string(1:1) == ch)
 else
   ans = ( toupper(string(1:1)) == toupper(ch))
 end if

end function firstchar_0d
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/firstchar_1d
!! NAME
!! firstchar_1d
!!
!! FUNCTION
!!  Returns .TRUE. is the first character of the string belongs to a given list.
!!
!! INPUTS
!!  string=The string whose first character has to be cheched
!!  char_list=The list of characters.
!!  [csens]=.TRUE. if comparison is done regardless of case. Defaults to .FALSE.
!!
!! PARENTS
!!
!! CHILDREN
!!
!!
!! SOURCE

pure function firstchar_1d(string,char_list,csens) result(ans)

 logical :: ans
 logical,optional,intent(in) :: csens
 character(len=*),intent(in) :: string
 character(len=1),intent(in) :: char_list(:)

!Local variables-------------------------------
 integer :: ii
 logical :: my_csens
 character(len=1) :: first_ch
!************************************************************************

 my_csens=.FALSE.; if (PRESENT(csens)) my_csens = csens

 first_ch = string(1:1)

 ans=.FALSE.

 if (.not.my_csens) then
   do ii=1,SIZE(char_list)
     ans = ( first_ch == char_list(ii) ); if (ans) EXIT
   end do
 else
   do ii=1,SIZE(char_list)
     ans = ( toupper(first_ch) == toupper(char_list(ii)) ); if (ans) EXIT
   end do
 end if

end function firstchar_1d
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/startswith
!! NAME
!! startswith
!!
!! FUNCTION
!!  Returns .TRUE. is the string starts with the specified prefix.
!!
!! SOURCE

pure function startswith(string, prefix) result(ans)

 logical :: ans
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: prefix

!Local variables-------------------------------
 integer :: ii,lenstr,lenpre
!************************************************************************

 ans = .False.
 lenstr = len_trim(string); lenpre = len_trim(prefix)
 if (lenpre > lenstr) return

 do ii=1,lenpre
   if (prefix(ii:ii) /= string(ii:ii)) return
 end do
 ans = .True.

end function startswith
!!***

!----------------------------------------------------------------------

!!****f* m_fstring/endswith
!! NAME
!! endswith
!!
!! FUNCTION
!!  Returns .TRUE. is the string ends with the specified suffix
!!
!! SOURCE

pure function endswith(string, suffix) result(ans)

 logical :: ans
 character(len=*),intent(in) :: string
 character(len=*),intent(in) :: suffix

!Local variables-------------------------------
 integer :: ii,p,lenstr,lensuf
!************************************************************************

 ans = .False.
 lenstr = len_trim(string); lensuf = len_trim(suffix)
 if (lensuf > lenstr) return

 do ii=1,lensuf
   p = lenstr - lensuf + ii
   if (suffix(ii:ii) /= string(p:p)) return
 end do
 ans = .True.

end function endswith
!!***

!!****f* m_fstrings/indent
!! NAME
!!  indent
!!
!! FUNCTION
!!  Indent text
!!
!! INPUTS
!!   istr=Input string
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function indent(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=len(istr)*4+4) :: ostr

!Local variables-------------------------------
 integer,parameter :: n=4 ! ostr is large enough to allocate all the possible indentations.
 integer :: ii,jj,kk
 character(len=1) :: ch

! *********************************************************************

 ostr = " "
 jj = n
 do ii=1,LEN_TRIM(istr)
   ch = istr(ii:ii)
   jj = jj + 1
   if (ch == NCHAR) then
      ostr(jj:jj) = NCHAR
      do kk=jj+1,jj+n
        ostr(kk:kk) = " "
      end do
      jj = jj+n
   else
     ostr(jj:jj) = ch
   end if
 end do
 !ostr(jj+1:) = "H"

end function indent
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/prep_dash
!! NAME
!!  prep_dash
!!
!! FUNCTION
!!  Prepend `-` to each line in a string.
!!
!! INPUTS
!!   istr=Input string
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function prep_dash(istr) result(ostr)

 character(len=*),intent(in) :: istr
 character(len=2*len(istr)) :: ostr

!Local variables-------------------------------
 integer :: ii,jj
 character(len=1) :: ch

! *********************************************************************
 ostr = ""
 jj = 1; ostr(jj:jj) = "-"
 !jj = 0

 do ii=1,LEN_TRIM(istr)
   ch = istr(ii:ii)
   jj = jj + 1
   if (ch == ch10) then
      ostr(jj:jj) = ch10
      ostr(jj+1:jj+1) = "-"
      jj = jj+1
   else
     ostr(jj:jj) = ch
   end if
 end do
 !ostr(jj+1:) = "H"

end function prep_dash
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/int2char4
!! NAME
!! int2char4
!!
!! FUNCTION
!! Convert an integer number to ("2") a character(len=*)
!! with trailing zeros if the number is <=9999.
!! Exemple : 123 will be mapped to "0123" ; 12345 will be mapped to "12345"
!! Makes sure that the integer fits the string length
!! (ex.: between 0 and 99999 if the string is a character(len=5)).
!!
!! INPUTS
!!  iint=integer to be converted
!!
!! OUTPUT
!!  string=character string ('####...' if error)
!!
!! PARENTS
!!      aim,anaddb,dtfil_init1,gaus_dos,get_all_gkq,iofn1,m_atprj,m_green
!!      m_io_redirect,m_phonon_supercell,m_self,mrgscr,optic,pawmkaewf
!!      prtfatbands,read_wfrspa,scfcv,tetrahedron
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine int2char4(iint,string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=*),intent(out) :: string

!Local variables-------------------------------
 integer :: lenstr

! *************************************************************************

 lenstr=min(len(string),25)
 if(iint<0 .or. iint>10._dp**(lenstr-1))then
   string=repeat('#',lenstr)
   return
 end if
 if(iint<10)then
   write(string,'("000",i1)')iint
 else if(iint<100)then
   write(string,'("00",i2)')iint
 else if(iint<1000)then
   write(string,'("0",i3)')iint
 else if(iint<10000)then
   write(string,'(i4)')iint
 else if(iint<1.0d5)then
   write(string,'(i5)')iint
 else if(iint<1.0d6)then
   write(string,'(i6)')iint
 else if(iint<1.0d7)then
   write(string,'(i7)')iint
 else if(iint<1.0d8)then
   write(string,'(i8)')iint
 else if(iint<1.0d9)then
   write(string,'(i9)')iint
 else if(iint<1.0d9)then
   write(string,'(i10)')iint
 else
   string=repeat('#',lenstr)
 end if

end subroutine int2char4
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/int2char10
!! NAME
!! int2char10
!!
!! FUNCTION
!! Convert a positive integer number (zero included) to a character(len=10),
!! with blanks to COMPLETE the string.
!! Exemple : 1234 will be mapped to "1234      "
!! Makes sure that the integer is between 0 and 9 999 999 999
!! Should be enough for integer*4
!!
!! INPUTS
!!  iint=integer to be converted
!!
!! OUTPUT
!!  string=character string ('##########' if error)
!!
!! PARENTS
!!      handle_ncerr,m_esymm,m_dfti,m_dyson_solver,m_qparticles,m_wfd
!!      prt_cif,wffile
!!
!! CHILDREN
!!
!! SOURCE

pure subroutine int2char10(iint,string)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iint
 character(len=10),intent(out) :: string

! *************************************************************************

!Note the use of floating numbers instead of large integers, for portability
 if(iint<0 .or. iint>=1.d10)then
   string='####'
   return
 end if
 if(iint<10)then
   write(string,'(i1,9x)')iint
 else if(iint<100)then
   write(string,'(i2,8x)')iint
 else if(iint<1.0d3)then
   write(string,'(i3,7x)')iint
 else if(iint<1.0d4)then
   write(string,'(i4,6x)')iint
 else if(iint<1.0d5)then
   write(string,'(i5,5x)')iint
 else if(iint<1.0d6)then
   write(string,'(i6,4x)')iint
 else if(iint<1.0d7)then
   write(string,'(i7,3x)')iint
 else if(iint<1.0d8)then
   write(string,'(i8,2x)')iint
 else if(iint<1.0d9)then
   write(string,'(i9,1x)')iint
 else
   write(string,'(i10)')iint
 end if

end subroutine int2char10
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/char_count
!! NAME
!! chcount
!!
!! FUNCTION
!!   Count the occurrences of a character in a string.
!!
!! PARENTS
!!
!! SOURCE

integer pure function char_count(string, char)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string
 character(len=1),intent(in) :: char
 integer :: i

! *************************************************************************

 char_count = 0
 do i=1,len(string)
   if (string(i:i) == char) char_count = char_count + 1
 end do

end function char_count
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/next_token
!! NAME
!! next_token
!!
!! FUNCTION
!!  Assume a string with whitespace-separated tokens.
!!  Find the next token starting from `start`, return it in `ostr` and update `start`
!!  so that one can call the function inside a loop.
!!  Return exit status.
!!
!! PARENTS
!!
!! SOURCE

integer function next_token(string, start, ostr) result(ierr)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: string
 character(len=*),intent(out) :: ostr
 integer,intent(inout) :: start

!Local variables-------------------------------
 integer :: ii,beg

! *************************************************************************
 !print *, "string:", trim(string(start:))

 ierr = 1; beg = 0
 ! Find first non-empty char.
 do ii=start,len_trim(string)
   if (string(ii:ii) /= " ") then
     beg = ii; exit
   end if
 end do
 if (beg == 0) return

 ! Find end of token.
 start = 0
 do ii=beg,len_trim(string)
   if (string(ii:ii) == " ") then
     start = ii; exit
   end if
 end do
 ! Handle end of string.
 if (start == 0) start = len_trim(string) + 1

 ierr = 0
 ostr = string(beg:start-1)

end function next_token
!!***

!----------------------------------------------------------------------

!!****f* m_fstrings/inupper
!! NAME
!! inupper
!!
!! FUNCTION
!! Maps all characters in string to uppercase except for tokens between quotation marks.
!! Uses fortran90 character string manipulation but should work
!! independent of EBCDIC or ASCII assumptions--only relies on
!! 'index' intrinsic character string matching function.
!! Makes sure that the string 'lolett' remains defined as the lower
!! case 26-character alphabet string and 'uplett' remains upper case.
!!
!! INPUTS
!!  string= character string with arbitrary case
!!
!! OUTPUT
!!  string= same character string mapped to upper case
!!
!! SIDE EFFECTS
!!  string= (input) character string with arbitrary case
!!          (output) same character string mapped to upper case
!!
!! PARENTS
!!      anaddb,band2eps,chkvars,intagm,invars1,m_ab7_invars_f90
!!      m_anaddb_dataset,m_exit,multibinit,parsefile
!!
!! CHILDREN
!!
!! SOURCE

subroutine inupper(string)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ii,indx,inquotes
 logical,save :: first=.true.
 character(len=1) :: cc
 !character(len=500) :: message
 character(len=26), parameter :: uplett='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 character(len=26), parameter :: lolett='abcdefghijklmnopqrstuvwxyz'

! *************************************************************************
!
!On first entry make sure lower case letters stayed
!lower case and upper case letters stayed upper case
 if (first) then
   do ii=1,26
     ! Look for occurrence of each upper case character
     ! anywhere in string of all lower case letters
     indx=index(lolett,uplett(ii:ii))
     ! If found then print error message and quit
     if (indx>0) then
       write(std_out, '(a,a,a,a,a,a,a,a,a)' )&
        'Upper case string = ',uplett,ch10,&
        'Lower case string = ',lolett,ch10,&
        'Upper case character ',uplett(ii:ii),'found in supposedly lower case string.'
       stop
     end if
   end do
   first=.false.
 end if

 inquotes = 0
 do ii=1,len_trim(string)
   !  Pick off single character of string (one byte):
   cc=string(ii:ii)

   ! Ignore tokens between quotation marks.
   if (cc == "'" .or. cc == '"') inquotes = inquotes + 1
   if (inquotes == 1) cycle
   if (inquotes == 2) then
     inquotes = 0; cycle
   end if
   ! determine whether a lowercase letter:
   indx=index(lolett,cc)
   ! Map to uppercase:
   if (indx>0) string(ii:ii)=uplett(indx:indx)
 end do

end subroutine inupper
!!***

end module m_fstrings
!!***

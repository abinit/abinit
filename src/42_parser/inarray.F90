!{\src2tex{textfont=tt}}
!!****f* ABINIT/inarray
!! NAME
!! inarray
!!
!! FUNCTION
!! Read the array of narr numbers located immediately after
!! a specified blank in a string of character.
!! Might read instead one word, after the specified blank.
!! Takes care of multipliers.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (DCA, XG).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inarray(b1,cs,dprarr,intarr,marr,narr,string,typevarphys)

 use defs_basis
 use m_errors
 use m_profiling_abi

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

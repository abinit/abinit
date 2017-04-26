!{\src2tex{textfont=tt}}
!!****f* ABINIT/inread
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
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine inread(string,ndig,typevarphys,outi,outr,errcod)

 use defs_basis
 use m_profiling_abi
 use m_errors

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

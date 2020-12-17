!!****m* ABINIT/m_abicore
!! NAME
!!  m_abicore
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group
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

module m_abicore

 !use defs_basis
 !use m_build_info
 use m_profiling_abi
 use m_specialmsg !,  only : herald, specialmsg_setcount, specialmsg_getcount, specialmsg_mpisum, wrtout
 !use m_errors

 implicit none

 public

contains  !=====================================================
!!***

! TODO: Replace with F90 i0

!!****f* ABINIT/appdig
!! NAME
!! appdig
!!
!! FUNCTION
!! Using input string "string" and integer "integ", make a string
!! named 'strinn' by concatenating digits of "integ" with characters
!! of "string"; return final string in "strinn".
!! Can also treat initial empty string, then simply returns the integer in the form of a string
!!
!! INPUTS
!! integ=nonnegative integer whose digits will be appended to string
!! string=string to which digits will be appended
!!
!! OUTPUT
!! strinn=string//nn
!!
!! PARENTS
!!      m_berryphase_new,m_d2frnl,m_dfpt_looppert,m_dfpt_lw,m_dfpt_nstwf
!!      m_dfpt_scfcv,m_dfptnl_loop,m_dtfil,m_elpolariz,m_gkk,m_ifc
!!      m_io_redirect,m_outvar_o_z,m_parser,m_pead_nl_loop
!!
!! CHILDREN
!!
!! SOURCE

subroutine appdig(integ,string,strinn)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: integ
 character(len=*),intent(in) :: string
 character(len=*),intent(out) :: strinn

!Local variables-------------------------------
!scalars
 integer :: i,length,ndig
 character(len=2) :: ncha
 character(len=8) :: form
 !character(len=500) :: msg

! *************************************************************************

 ! Check that integer is nonnegative
 !if (integ<0) then
 !  write(msg,'(a,i0,a)') &
 !  'Input integer =',integ,' must not be <0. Argument integ was input as negative.'
 !  MSG_BUG(msg)
 !end if

 ! Fill output string initially with blanks to end of dimensioned length
 length=len(strinn)
 do i=1,length
   strinn(i:i)=' '
 end do

!Find nonwhitespace length of string
 length=len_trim(string)
!Copy input character string into first part of output string
 if(length>0)then
   strinn(1:length)=string(1:length)
 end if

!Find how many digits "integ" has
 ndig=int(log10(real(integ)+0.50))+1

!Create a format for exact number of digits using internal write
 write(unit=ncha,fmt='(i2)') ndig
 form='(i'//ncha//')'
!Do internal write to get digits of integer into character string,
!placing digits into appropriate end of string.
 write(unit=strinn(1+length:length+ndig),fmt=form) integ
!(Note that present version writes "1" or "2" for single digit,
!not "01" or "02".  Latter may be preferable.  Can be amended.)
!

end subroutine appdig
!!***

end module m_abicore
!!***

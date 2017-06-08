!{\src2tex{textfont=tt}}
!!****f* ABINIT/inupper
!! NAME
!! inupper
!!
!! FUNCTION
!! Maps all characters in string to uppercase except for tokes between quotation marks.
!! Uses fortran90 character string manipulation but should work
!! independent of EBCDIC or ASCII assumptions--only relies on
!! 'index' intrinsic character string matching function.
!! Makes sure that the string 'lolett' remains defined as the lower
!! case 26-character alphabet string and 'uplett' remains upper case.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR).
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine inupper(string)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'inupper'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=*),intent(inout) :: string

!Local variables-------------------------------
!scalars
 integer :: ii,indx,inquotes
 logical,save :: first=.true.
 character(len=1) :: cc
 character(len=500) :: message
 character(len=26), parameter :: uplett='ABCDEFGHIJKLMNOPQRSTUVWXYZ'
 character(len=26), parameter :: lolett='abcdefghijklmnopqrstuvwxyz'

! *************************************************************************
!
!On first entry make sure lower case letters stayed
!lower case and upper case letters stayed upper case
 if (first) then
   do ii=1,26
!    Look for occurrence of each upper case character
!    anywhere in string of all lower case letters
     indx=index(lolett,uplett(ii:ii))
!    If found then print error message and quit
     if (indx>0) then
       write(message, '(a,a,a,a,a,a,a,a,a)' )&
&       'Upper case string = ',uplett,ch10,&
&       'Lower case string = ',lolett,ch10,&
&       'Upper case character ',uplett(ii:ii),'found in supposedly lower case string.'
       MSG_BUG(message)
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

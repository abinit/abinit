!{\src2tex{textfont=tt}}
!!****f* ABINIT/fappnd
!!
!! NAME
!! fappnd
!!
!! FUNCTION
!! Create the modified root name to be used for output of density, potential,
!! and geometry files. See the description of the iapp input variable.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam= generic output root name
!! iapp=indicates the eventual suffix to be appended to the generic output root
!!      (the suffixe depends on the presence of the suff (optional) argument.
!!        if 0 : no suffix to be appended (called directly from gstate)
!!        if positive : append "_SUF//iapp" (called from move or brdmin)
!!        if -1 : append "_SUF0" (called from brdmin)
!!        if -2, -3, -4, -5: append "_SUFA", ... ,"_SUFD", (called from move)
!!      SUF can be TIM (default) or IMG
!! suff= --optional argument--indicates the suffixe to be appended:
!!       SUF=TIM (default) or SUF=IMG or ...
!! OUTPUT
!! filapp= filename with appended string
!!
!! PARENTS
!!      dtfil_init_time
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fappnd(filapp,filnam,iapp,&
&                 suff) ! optional argument

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fappnd'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iapp
 character(len=fnlen),intent(in) :: filnam
 character(len=fnlen),intent(out) :: filapp
 character(len=3),optional,intent(in) :: suff

!Local variables-------------------------------
!scalars
 integer :: ndig
 character(len=3) :: suffixe
 character(len=8) :: nchar
 character(len=500) :: msg

! *************************************************************************

 if(iapp==0)then
   filapp=trim(filnam)
 else
   suffixe="TIM"
   if (present(suff)) suffixe=trim(suff(1:3))
   if(iapp>0)then
!    Create character string for filename. Determine the number of digits in iapp.
     ndig=int(log10(dble(iapp)+0.5_dp))+1
!    Make integer format field of exact size (internal write)
!    for assumed nchar string of 8 characters
     write(nchar, '(i8)' ) iapp
     if (ndig>8) then
       write(msg,'(5a,i0,2a,i0,2a)')&
&       'Requested file name extension has more than the allowed 8 digits.',ch10,&
&       'Action : resubmit the job with smaller value for ntime.',ch10,&
&       'Value computed here was ndig=',ndig,ch10,&
&       'iapp=',iapp,' filnam=',TRIM(filnam)
       MSG_ERROR(msg)
     end if
!    Concatenate into character string, picking off exact number of digits
!    The potential or density label will be appended in ioarr
     filapp=trim(filnam)//'_'//suffixe(1:3)//nchar(9-ndig:8)
   else if(iapp==-1)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'0'
   else if(iapp==-2)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'A'
   else if(iapp==-3)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'B'
   else if(iapp==-4)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'C'
   else if(iapp==-5)then
     filapp=trim(filnam)//'_'//suffixe(1:3)//'D'
   end if
 end if

end subroutine fappnd
!!***

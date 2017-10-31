!{\src2tex{textfont=tt}}
!!****f* ABINIT/isfile
!! NAME
!! isfile
!!
!! FUNCTION
!! Inquire Status of FILE
!! Checks that for status =
!! 'old': file already exists
!! 'new': file does not exist; if file exists,
!! filnam is modified to filnam.A or filnam.B,....
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, JJ)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! filnam=character string to specify filename
!! status='old' or 'new'
!!
!! OUTPUT
!! stops processing if old file does not exist; changes name
!! and returns new name in redefined filnam if new file already exists.
!!
!! PARENTS
!!      anaddb,iofn1,m_effective_potential,m_polynomial_coeff,m_vcoul
!!      multibinit,ujdet
!!
!! CHILDREN
!!      clib_rename,int2char4
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine isfile(filnam,status)

 use defs_basis
 use m_errors

 use m_clib, only : clib_rename
 use m_fstrings, only : int2char4

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'isfile'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=3),intent(in) :: status
 character(len=fnlen),intent(inout) :: filnam

!Local variables-------------------------------
!scalars
 logical :: ex,found
 integer :: ii,ios, ioserr
 character(len=500) :: message
 character(len=fnlen) :: filnam_tmp
 character(len=fnlen) :: trialnam

! *************************************************************************

 filnam_tmp=filnam

 if (status=='old') then !  Check that old file exists
   inquire(file=filnam,iostat=ios,exist=ex)

   if (ios/=0) then
     write(message,'(4a,i0,2a)')&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE statement returns error code',ios,ch10,&
&     'Action: identify which problem appears with this file.'
     MSG_ERROR(message)
   else if (.not.ex) then
     write(message, '(5a)' )&
&     'Checks for existence of file  ',trim(filnam),ch10,&
&     'but INQUIRE finds file does not exist.',&
&     'Action: check file name and re-run.'
     MSG_ERROR(message)
   end if

 else if (status=='new') then

   ! Check that new output file does NOT exist
   ioserr = 0
   trialnam = filnam
   ii = 0
   inquire(file=trim(trialnam),iostat=ios,exist=ex)
   if ( ios /= 0 ) then
     write(message,'(4a)') 'Something is wrong with permissions for ', &
&     'reading/writing on this filesystem.',ch10,&
&     'Action : Check permissions.'
     MSG_ERROR(message)
   end if

   if ( ex .eqv. .true. ) then
     write(message,'(4a)')'Output file ',trim(trialnam),ch10,' already exists.'
     MSG_COMMENT(message)
     found=.false.

     ii=1
     do while ( (found .eqv. .false.) .and. (ii < 10000) )
       call int2char4(ii,message)
       trialnam=trim(trim(filnam_tmp)//message)
       inquire(file=trim(trialnam),iostat=ios,exist=ex)
       if ( (ex .eqv. .false.) .and. (ios == 0)) then
         found  = .true.
       end if
       if ( ios /= 0 )  ioserr=ioserr+1
       if ( ioserr > 10 ) then
!        There is a problem => stop
         write(message, '(2a,i0,2a)' )&
&         'Check for permissions of reading/writing files on the filesystem', &
&         '10 INQUIRE statements returned an error code like ',ios,ch10,&
&         'Action: Check permissions'
         MSG_ERROR(message)
       end if
       ii=ii+1
     end do
     if ( found .eqv. .true. ) then
       write(message,'(4a)') 'Renaming old ',trim(filnam),' to ',trim(trialnam)
       MSG_COMMENT(message)
       call clib_rename(filnam,trialnam,ioserr)
       if ( ioserr /= 0 ) then
         write(message,'(4a)') 'Failed to rename file ', trim(filnam),' to ',trim(trialnam)
         MSG_ERROR(message)
       end if
     else
       write(message,'(3a)')&
&       'Have used all names of the form filenameXXXX, X in [0-9]',ch10,&
&       'Action: clean up your directory and start over.'
       MSG_ERROR(message)
     end if
   end if

   ! if ii > 0 we iterated so rename abi_out to abi_outXXXX
   ! and just write to abi_out
 else ! status not recognized
   write(message,'(3a)')'  Input status= ',status,' not recognized.'
   MSG_BUG(message)
 end if

end subroutine isfile
!!***

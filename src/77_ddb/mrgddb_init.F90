!{\src2tex{textfont=tt}}
!!****f* ABINIT/mrgddb_init
!!
!! NAME
!! mrgddb_init
!!
!! FUNCTION
!! Initialize the code : write heading and make the first i/os
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! dscrpt=character string that describe the derivative database
!! filnam(3)=character strings giving file names
!! nddb=(=1 => will initialize the ddb, using an input GS file)
!!  (>1 => will merge the whole set of ddbs listed)
!!
!! OUTPUT
!!  (None)
!!
!! NOTES
!! 1. To be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Output Derivative Database
!!    if nddb==1,
!!     (2) Formatted input file for the Corning ground-state code
!!    if nddb>1,
!!     (2 ... nddb+1) Derivative Databases to be added
!!
!! PARENTS
!!      mrgddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mrgddb_init(dscrpt,filnam,mddb,nddb)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgddb_init'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: mddb
 integer,intent(out) :: nddb
 character(len=*),intent(out) :: dscrpt
!arrays
 character(len=*),intent(out) :: filnam(mddb+1)

!Local variables -------------------------
!scalars
 integer :: iddb
 character(len=500) :: message

! *********************************************************************

!Read the name of the output ddb
 write(std_out,*)' Give name for output derivative database : '
 read(std_in, '(a)' ) filnam(1)
 write(std_out,'(a,a)' )' ',trim(filnam(1))

!Read the description of the derivative database
 write(std_out,*)' Give short description of the derivative database :'
 read(std_in, '(a)' )dscrpt
 write(std_out,'(a,a)' )' ',trim(dscrpt)

!Read the number of input ddbs, and check its value
 write(std_out,*)' Give number of input ddbs, or 1 if input GS file'
 read(std_in,*)nddb
 write(std_out,*)nddb
 if(nddb<=0.or.nddb>mddb)then
   write(message, '(a,a,i0,a,i0,a,a,a)' )&
&   'nddb should be positive, >1 , and lower',&
&   'than mddb= ',mddb,' while the input nddb is ',nddb,'.',ch10,&
&   'Action: correct the input nddb.'
   MSG_ERROR(message)
 end if

!Read the file names
 if(nddb==1)then
   write(std_out,*)' Give name for ABINIT input file : '
   read(std_in, '(a)' ) filnam(2)
   write(std_out,'(a,a)' )' ',trim(filnam(2))
 else
   do iddb=1,nddb
     write(std_out,*)' Give name for derivative database number',iddb,' : '
     read(std_in, '(a)' ) filnam(iddb+1)
     write(std_out,'(a,a)' )' ',trim(filnam(iddb+1))
   end do
 end if

end subroutine mrgddb_init
!!***

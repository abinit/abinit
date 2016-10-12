!{\src2tex{textfont=tt}}
!!****f* ABINIT/init10
!!
!! NAME
!! init10
!!
!! FUNCTION
!! Initialize the code multibinit: write heading and make the first i/os
!!
!! COPYRIGHT
!! Copyright (C) 1999-2015 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! character(len=fnlen) filnam(4)=character strings giving file names
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Formatted input file
!!     (2) Formatted output file
!!     (3) Input Derivative Database (DDB file)
!!     (4-12) Input Derivative Database (XML format)
!! 
!! PARENTS
!!      anaddb
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine init10(filnam,comm)

 use defs_basis
 use m_xmpi
 use m_errors
 use m_ab7_invars

 use m_fstrings,     only : int2char4
 use m_io_tools,     only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init10'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm
!arrays
 character(len=*),intent(out) :: filnam(15)

!Local variables--------------------------
!scalars
 integer,parameter :: master=0
 integer :: me,nproc,ierr
 integer :: ii,io
!arrays
! *********************************************************************
 
!Determine who I am in comm
 me = xmpi_comm_rank(comm)
 nproc = xmpi_comm_size(comm)

 filnam(:) = ""

 if (me==master)then
!Read the file names
   write(std_out,*)' Give name for      formatted input file : '
   read(std_in, '(a)',IOSTAT=io) filnam(1)
   write(std_out,'(a,a)' )'-   ',trim(filnam(1))
   write(std_out,*)' Give name for     formatted output file : '
   read(std_in, '(a)',IOSTAT=io) filnam(2)
   write(std_out,'(a,a)' )'-   ',trim(filnam(2))
   write(std_out,*)' Give name for input derivative database of reference structure',&
&                  ' (DDB or XML file): '
   read(std_in, '(a)',IOSTAT=io) filnam(3)
   write(std_out,'(a,a)' )'-   ',trim(filnam(3))
   ii = 4
   do while (io>=0 .and. ii<16)
     write(std_out,*)' Give name for input derivative database (DDB or XML file): '
     read(std_in, '(a)',IOSTAT=io) filnam(ii)
     write(std_out,'(a,a)' )'-   ',trim(filnam(ii))
     if(trim(filnam(ii))==" ") exit
     ii = ii + 1
   end do
 end if

!Communicate filenames to all processors
 call xmpi_bcast (filnam, master, comm, ierr)

end subroutine init10
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/init10
!!
!! NAME
!! init10
!!
!! FUNCTION
!! Initialize the code epigene: write heading and make the first i/os
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


subroutine init10(filnam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'init10'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!scalars
 integer :: ii,io
!arrays
 character(len=*),intent(out) :: filnam(15)

! *********************************************************************
 
 filnam(:) = ""

!Read the file names
 write(std_out,*)' Give name for      formatted input file : '
 read(std_in, '(a)',IOSTAT=io) filnam(1)
 write(std_out,'(a,a)' )'-   ',trim(filnam(1))
 write(std_out,*)' Give name for     formatted output file : '
 read(std_in, '(a)',IOSTAT=io) filnam(2)
 write(std_out,'(a,a)' )'-   ',trim(filnam(2))
 write(std_out,*)' Give name for input derivative database of reference structure',&
& ' (DDB or XML file): '
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
end subroutine init10
!!***

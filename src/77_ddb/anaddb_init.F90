!{\src2tex{textfont=tt}}
!!****f* ABINIT/anaddb_init
!!
!! NAME
!! anaddb_init
!!
!! FUNCTION
!! Initialize the code ppddb9: write heading and make the first i/os
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!
!! OUTPUT
!! character(len=fnlen) filnam(7)=character strings giving file names
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Formatted input file
!!     (2) Formatted output file
!!     (3) Input Derivative Database
!!     (4) Output Molecular Dynamics
!!     (5) Input electron-phonon matrix elements
!!     (6) Root name for electron-phonon file names
!!     (7) Name of file containing the 3 ddk filenames and the GS wf file name
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


subroutine anaddb_init(filnam)

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'anaddb_init'
!End of the abilint section

 implicit none

!Arguments -------------------------------
!arrays
 character(len=*),intent(out) :: filnam(7)

! *********************************************************************

!Read the file names
 write(std_out,*)' Give name for      formatted input file : '
 read(std_in, '(a)' ) filnam(1)
 write(std_out,'(a,a)' )'-   ',trim(filnam(1))
 write(std_out,*)' Give name for     formatted output file : '
 read(std_in, '(a)' ) filnam(2)
 write(std_out,'(a,a)' )'-   ',trim(filnam(2))
 write(std_out,*)' Give name for input derivative database : '
 read(std_in, '(a)' ) filnam(3)
 write(std_out,'(a,a)' )'-   ',trim(filnam(3))
 write(std_out,*)' Give name for output molecular dynamics : '
 read(std_in, '(a)' ) filnam(4)
 write(std_out,'(a,a)' )'-   ',trim(filnam(4))
 write(std_out,*)' Give name for input elphon matrix elements  (GKK file) : '
 read(std_in, '(a)' ) filnam(5)
 write(std_out,'(a,a)' )'-   ',trim(filnam(5))
 write(std_out,*)' Give root name for elphon output files : '
 read(std_in, '(a)' ) filnam(6)
 write(std_out,'(a,a)' )'-   ',trim(filnam(6))
 write(std_out,*)' Give name for file containing ddk filenames for elphon/transport : '
 read(std_in, '(a)' ) filnam(7)
 write(std_out,'(a,a)' )'-   ',trim(filnam(7))

end subroutine anaddb_init
!!***

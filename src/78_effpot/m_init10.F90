!!****m* ABINIT/m_init10
!! NAME
!!  m_init10
!!
!! FUNCTION
!!   It should be "contained" in multibinit but abilint does not accept "contains" in programs.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group ()
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
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

module m_init10

 use defs_basis
 use m_errors
 use m_abicore
 use m_fstrings,     only : int2char4, rmquotes, sjoin, strcat, basename

 implicit none

 private
!!***

 public :: init10
!!***

contains
!!***


!!****f* ABINIT/init10
!!
!! NAME
!! init10
!!
!! FUNCTION
!! Initialize the code multibinit: write heading and make the first i/os
!!
!! INPUTS
!!  character(len=fnlen) input_path: gives the name of the input file.
  !! If not given, the files file are then used with a deprecation message.
!! OUTPUT
!! character(len=fnlen) filnam(6)=character strings giving file names
!!
!! NOTES
!! 1. Should be executed by one processor only.
!! 2. File names refer to following files, in order:
!!     (1) Formatted input file
!!     (2) Formatted output file
!!     (3) Input for reference structure and harmonic part (DDB file or XML)
!!     (4) Input for XML with polynomial coefficients (DDB file)
!!     (5) Input for HIST Training-Set file (netcdf .nc format)
!!     (6) Input for HIST Test-Set file (netcdf .nc format)
!!     (7-18) Input Derivative Database (XML format)
!!
!! PARENTS
!!      multibinit
!!
!! CHILDREN
!!      xmpi_bcast
!!
!! SOURCE

subroutine init10(input_path, filnam,comm)

 use defs_basis
 use m_xmpi
 use m_errors

 use m_fstrings,     only : int2char4
 use m_io_tools,     only : open_file

 implicit none

!Arguments -------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=fnlen), intent(in) :: input_path
!arrays
 character(len=fnlen),intent(out) :: filnam(18)

!Local variables--------------------------
!scalars
 integer,parameter :: master=0
 integer :: me,nproc,ierr
 integer :: ii,io, i1, i2

 character(len=fnlen) :: fname
! *********************************************************************

!Determine who I am in comm
 me = xmpi_comm_rank(comm)
 nproc = xmpi_comm_size(comm)

 filnam(:) = ""

 if (me==master)then
    print *, "input_path:", input_path
    print *, "input_path:", trim(input_path)
    print *, "input_path:", len_trim(input_path)
    if (len_trim(input_path) == 0) then
       ! Legacy Files file mode.
       write(std_out, "(2a)")" DeprecationWarning: ",ch10
       write(std_out, "(a)") "     The files file has been deprecated in Abinit9 and will be removed in Abinit10."
       write(std_out, "(a)")"     Use the syntax `multibinit t01.abi` where t01.abi is an input."
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
       write(std_out,*)' Give name for input coefficients from fitted polynomial',&
            &                  ' (XML file or enter no): '
       read(std_in, '(a)',IOSTAT=io) filnam(4)
       write(std_out,'(a,a)' )'-   ',trim(filnam(4))
       write(std_out,*)' Give name for training-set file',&
            &                  ' (netcdf file or enter no): '
       read(std_in, '(a)',IOSTAT=io) filnam(5)
       write(std_out,'(a,a)' )'-   ',trim(filnam(5)) 
       write(std_out,*)' Give name for test-set file',&
            &                  ' (netcdf file or enter no): '
       read(std_in, '(a)',IOSTAT=io) filnam(6)
       write(std_out,'(a,a)' )'-   ',trim(filnam(6)) 
       ii = 7
       do while (io>=0 .and. ii<19)
          write(std_out,*)' Give name for input derivative database (DDB or XML file): '
          read(std_in, '(a)',IOSTAT=io) filnam(ii)
          write(std_out,'(a,a)' )'-   ',trim(filnam(ii))
          if(trim(filnam(ii))==" ") exit
          ii = ii + 1
       end do
    else
       filnam(1)=input_path
       filnam(2) = trim(input_path)//".abo"
       filnam(3) = "i"
       filnam(4) = "o"
       filnam(5) = "t"

       fname = basename(input_path)
       i1 = index(fname, ".")
       !if (i1 /= 0) then
       if (i1 > 1) then
          ! file ext is present --> use prefix to initialize filnam
          !TODO:  read the input file and set the filnam array
          i2 = index(input_path, ".", back=.True.)
          filnam(2) = input_path(:i2) // "abo"
          ! file ext is present --> use prefix to initialize filnam
          !TODO:  read the input file and set the filnam array
       
       end if


    end if
 end if

!Communicate filenames to all processors
 call xmpi_bcast (filnam, master, comm, ierr)

end subroutine init10
!!***

end module m_init10
!!***

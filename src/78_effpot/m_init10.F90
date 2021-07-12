!!****m* ABINIT/m_init10
!! NAME
!!  m_init10
!!
!! FUNCTION
!!   It should be "contained" in multibinit but abilint does not accept "contains" in programs.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2021 ABINIT group ()
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
 use m_multibinit_dataset, only: multibinit_dtset_type

 implicit none

 private
!!***

 public :: init10
 public :: postfix_fnames
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
!! character(len=fnlen) filnam(18)=character strings giving file names
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

       fname = basename(input_path)
       i1 = index(fname, ".")
       if (i1 > 1) then
          i2 = index(input_path, ".", back=.True.)
          filnam(2) = input_path(:i2) // "abo"
          ! The rest filnam(3:) are set after reading from input file.
          ! as in postfix_filenames subroutine
       end if
    end if
 end if


!Communicate filenames to all processors
 call xmpi_bcast (filnam, master, comm, ierr)

end subroutine init10
!!***

!===============================================================
! Change the extension of a filename. e.g. run.abo -> run_hist.nc
!> @ fname: the old fname
!> @ new_ext: the new extension. e.g. .nc or _hist.nc
!===============================================================
function change_extension(fname, new_ext) result (ret)
    character(len=fnlen), intent(in) :: fname, new_ext
    character(len=fnlen) :: tmp_fname
    character(len=fnlen) :: ret
    integer :: i1, i2
    tmp_fname=basename(fname)
    i1 = index(tmp_fname, ".")
    if (i1 > 1) then
       i2 = index(fname, ".", back=.True.)
       ret = fname(:i2) // trim(new_ext)
    end if
end function change_extension

!===============================================================
! Sync the filename info from files file and from input file
!> @ input_path: input file
!> @ filnam : the filenames from files file
!> @ params: data read from inputfile
!===============================================================
subroutine postfix_fnames(input_path, filnam, params)
  character(len=fnlen), intent(inout) :: input_path, filnam(18)
  type(multibinit_dtset_type), intent(inout) :: params
  ! if using files file, set the params%*_fname according to filnam
  integer :: i
  if (len_trim(input_path)==0) then
     if (len_trim(params%spin_pot_fname)==0) params%spin_pot_fname=filnam(3)
     if (len_trim(params%latt_pot_fname)==0) params%latt_pot_fname=filnam(3)
     if (len_trim(params%slc_pot_fname)==0) params%slc_pot_fname=filnam(3)
     ! TODO lattice model
     if (len_trim(params%latt_inp_ddb_fname)==0) params%latt_inp_ddb_fname=filnam(3)
     if (len_trim(params%latt_inp_coeff_fname)==0) params%latt_inp_coeff_fname=filnam(4)
     if (len_trim(params%latt_training_set_fname)==0) params%latt_training_set_fname=filnam(5)
     if (len_trim(params%latt_test_set_fname)==0) params%latt_test_set_fname=filnam(6)
     do i=1, 12
        if (len_trim(params%latt_ddb_fnames(i))==0) params%latt_ddb_fnames(i)=filnam(6+i)
     end do
  else
     filnam(2)=params%outdata_prefix
     filnam(3)=params%latt_inp_ddb_fname
     filnam(4)=params%latt_inp_coeff_fname
     filnam(5)=params%latt_training_set_fname
     filnam(6)=params%latt_test_set_fname
     do i=1, 12
        filnam(i+6)=params%latt_ddb_fnames(i)
     end do
  end if

end subroutine postfix_fnames



end module m_init10
!!***

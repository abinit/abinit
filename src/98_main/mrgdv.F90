!{\src2tex{textfont=tt}}
!!****p* ABINIT/mrgdv
!! NAME
!! mrgdv
!!
!! FUNCTION
!! This program merges DFPT potentials for different q-vectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2018 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! DVDB file format:
!!   version (integer)
!!   number of potentials (integer)
!!   for each potential:
!!     Abinit header with info on the perturbation and the FFT mesh
!!     potential on the FFT mesh
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,dvdb_free,dvdb_init
!!      dvdb_list_perts,dvdb_merge_files,dvdb_print,dvdb_test_ftinterp
!!      dvdb_test_v1complete,dvdb_test_v1rsym,get_command_argument,herald
!!      prompt,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgdv

 use defs_basis
 use defs_abitypes
 use m_xmpi
 use m_errors
 use m_build_info
 use m_abicore
 use m_dvdb

 use m_specialmsg,      only : specialmsg_getcount, herald
 use m_fstrings,        only : sjoin, itoa, ltoa
 use m_numeric_tools,   only : vdiff_eval, vdiff_print
 use m_io_tools,        only : file_exists, prompt

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mrgdv'
!End of the abilint section

 implicit none

!Local variables-------------------------------
!scalars
 integer :: ii,nargs,nfiles,comm,prtvol,my_rank,ierr
 character(len=24) :: codename
 character(len=500) :: command,arg, msg
 character(len=fnlen) :: db_path,dump_file
 type(dvdb_t) :: db
!arrays
 integer :: ngqpt(3)
 character(len=fnlen),allocatable :: v1files(:)

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world
 my_rank = xmpi_comm_rank(comm)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 ! write greating,read the file names, etc.
 codename='MRGDV'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

 ABI_CHECK(xmpi_comm_size(comm) == 1, "Not programmed for parallel execution")
 prtvol = 0
 !prtvol = 10

 nargs = command_argument_count()

 if (nargs == 0) then
   ! We are reading from stdin
   call prompt("Enter name of output file:", db_path)
   call prompt("Enter total number of DFPT POT files:", nfiles)
   ABI_MALLOC(v1files, (nfiles))
   do ii=1,nfiles
     call prompt(sjoin("Enter name of POT file",itoa(ii),":"), v1files(ii))
   end do
   call dvdb_merge_files(nfiles, v1files, db_path, prtvol)
   ABI_FREE(v1files)

 else
   ! Command line options.
   do ii=1,command_argument_count()
     call get_command_argument(ii, arg)
     if (arg == "-v" .or. arg == "--version") then
       write(std_out,"(a)") trim(abinit_version); goto 100

     else if (arg == "-h" .or. arg == "--help") then
       ! Document the options.
       write(std_out,*)"-v, --version              Show version number and exit."
       write(std_out,*)"merge out_DVDB POT1 POT2   Merge list of POT files, produce out_DVDB file."
       write(std_out,*)"info out_DVDB              Print information on DVDB file"
       write(std_out,*)"-h, --help                 Show this help and exit."
       write(std_out,*)" "
       write(std_out,*)"Options for developers:"
       write(std_out,*)"test_v1complete [file]     Test symmetrization of DFPT potentials."
       write(std_out,*)"                           Assume DVDB with all 3*natom perturbations for each q (prep_gkk)."
       write(std_out,*)"test_v1rsym                Test symmetries of DFPT potentials in real space."
       write(std_out,*)"test_ftinterp [n1,n2,n3]   Test Fourier interpolation of DFPT potentials."
       goto 100
     end if
   end do

   call get_command_argument(1, command)

   select case (command)
   case ("merge")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, db_path)
     if (file_exists(db_path)) then
       MSG_ERROR(sjoin("Cannot overwrite existing file:", db_path))
     end if

     nfiles = nargs - 2
     ABI_MALLOC(v1files, (nfiles))
     do ii=1,nfiles
       call get_command_argument(ii+2, v1files(ii))
     end do

     ! Merge POT files.
     call dvdb_merge_files(nfiles, v1files, db_path, prtvol)
     ABI_FREE(v1files)

   case ("info")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, db_path)

     call dvdb_init(db, db_path, comm)
     call dvdb_print(db, prtvol=prtvol)
     call dvdb_list_perts(db, [-1,-1,-1])
     call dvdb_free(db)

   case ("test_v1comp", "test_v1complete")
     call wrtout(std_out," Testing symmetries (assuming overcomplete DVDB, pass extra argument to dump v1(r)) to file")
     call get_command_argument(2, db_path)
     dump_file = ""; if (nargs > 2) call get_command_argument(3, dump_file)
     call dvdb_test_v1complete(db_path, dump_file, comm)

   case ("test_v1rsym")
     call wrtout(std_out," Testing symmetries of V1(r) in real space.")
     call get_command_argument(2, db_path)
     call dvdb_test_v1rsym(db_path, comm)

   case ("test_ftinterp")
     call get_command_argument(2, db_path)
     ngqpt = [2,2,2]
     if (nargs > 2) then
       call get_command_argument(3, arg)
       read(arg, *, iostat=ierr, iomsg=msg)ngqpt
       ABI_CHECK(ierr == 0, msg)
     end if

     write(std_out,"(a)")sjoin("Testing Fourier interpolation of V1(r) with ngqpt:", ltoa(ngqpt))
     call dvdb_test_ftinterp(db_path, ngqpt, comm)

   case default
     MSG_ERROR(sjoin("Unknown command:", command))
   end select

 end if

 call wrtout(std_out," Done",'COLL')

 call abinit_doctor("__mrgdv")

 100 call xmpi_end()

 end program mrgdv
!!***

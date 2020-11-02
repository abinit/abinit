!!****p* ABINIT/mrgdv
!! NAME
!! mrgdv
!!
!! FUNCTION
!! This program merges DFPT potentials for different q-vectors and perturbations.
!!
!! COPYRIGHT
!! Copyright (C) 2004-2020 ABINIT group (MG)
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
!!      abi_io_redirect,abimem_init,abinit_doctor,dvdb%free,dvdb%list_perts
!!      dvdb%open_read,dvdb%print,dvdb%qdownsample,dvdb_merge_files
!!      dvdb_test_ftinterp,dvdb_test_v1complete,dvdb_test_v1rsym
!!      get_command_argument,herald,ngfft_seq,prompt,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgdv

 use defs_basis
 use m_xmpi
 use m_errors
 use m_build_info
 use m_abicore
 use m_dvdb

 use m_specialmsg,      only : specialmsg_getcount, herald
 use m_fstrings,        only : sjoin, itoa, ltoa
 use m_numeric_tools,   only : vdiff_eval, vdiff_print
 use m_io_tools,        only : file_exists, prompt
 use m_argparse,        only : get_arg, get_arg_list
 use m_fftcore,         only : ngfft_seq

 implicit none

!Local variables-------------------------------
!scalars
 integer :: ii, nargs, nfiles, comm, prtvol, my_rank, lenr, dvdb_add_lr, rspace_cell, symv1scf, npert_miss, abimem_level
 real(dp) :: dvdb_qdamp, abimem_limit_mb
 character(len=24) :: codename
 character(len=500) :: command, arg, msg
 character(len=fnlen) :: dvdb_filepath, dump_file, ddb_filepath
 type(dvdb_t) :: dvdb
!arrays
 integer :: ngqpt(3), coarse_ngqpt(3), ngfftf(18)
 character(len=fnlen),allocatable :: v1files(:)

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world
 my_rank = xmpi_comm_rank(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
 ABI_CHECK(get_arg("abimem-level", abimem_level, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("abimem-limit-mb", abimem_limit_mb, msg, default=20.0_dp) == 0, msg)
#ifdef HAVE_MEM_PROFILING
 call abimem_init(abimem_level, limit_mb=abimem_limit_mb)
#endif

 ! write greating,read the file names, etc.
 codename='MRGDV'//repeat(' ',18)
 call herald(codename, abinit_version, std_out)

 ABI_CHECK(xmpi_comm_size(comm) == 1, "Not programmed for parallel execution")
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)

 nargs = command_argument_count()

 if (nargs == 0) then
   ! We are reading from stdin
   ! Prepend prompt with `-` to bypass bug in intel18-19
   call prompt("Enter name of output file:", dvdb_filepath)
   call prompt("Enter total number of DFPT POT files:", nfiles)
   ABI_MALLOC(v1files, (nfiles))
   do ii=1,nfiles
     call prompt(sjoin("Enter name of POT file", itoa(ii), ":"), v1files(ii))
   end do
   call dvdb_merge_files(nfiles, v1files, dvdb_filepath, prtvol)
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
       write(std_out,*)"=== Options for developers ==="
       write(std_out,*)" "
       write(std_out,*)"test_v1complete FILE [--symv1scf 1] [--potfile foo.nc]"
       write(std_out,*)"                           Test symmetrization of DFPT potentials with symv1scf option"
       write(std_out,*)"                           Assume DVDB with all 3*natom perturbations for each q (use prep_gkk)."
       write(std_out,*)"test_v1rsym [--symv1scf]   Test symmetries of DFPT potentials in real space."
       write(std_out,*)"test_ftinterp in_DVDB --ngqpt 4 4 4 [--ddb-path] [--dvdb-add-lr 0] [--qdamp -1]"
       write(std_out,*)"                                    [--symv1scf] [--coarse-ngqpt 2 2 2]"
       write(std_out,*)"                           Test Fourier interpolation of DFPT potentials."
       write(std_out,*)"downsample in_DVDB out_DVDB [n1, n2, n3] Produce new DVDB with q-subsmesh"
       !write(std_out,*)"convert in_old_DVDB out_DVDB.nc  Convert old DVDB format to new DVDB in netcdf format"
       goto 100
     end if
   end do

   call get_command_argument(1, command)

   select case (command)
   case ("merge")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, dvdb_filepath)
     if (file_exists(dvdb_filepath)) then
       MSG_ERROR(sjoin("Cannot overwrite existing file:", dvdb_filepath))
     end if

     nfiles = nargs - 2
     ABI_MALLOC(v1files, (nfiles))
     do ii=1,nfiles
       call get_command_argument(ii+2, v1files(ii))
     end do

     ! Merge POT files.
     call dvdb_merge_files(nfiles, v1files, dvdb_filepath, prtvol)
     ABI_FREE(v1files)

   case ("info")
     ! Get name of output database and list of v1 files.
     ABI_CHECK(nargs > 1, "Additional arguments are missing")
     call get_command_argument(2, dvdb_filepath)

     dvdb = dvdb_new(dvdb_filepath, comm)
     if (prtvol > 0) call dvdb%print(prtvol=prtvol)
     call dvdb%list_perts([-1, -1, -1], npert_miss)
     call dvdb%free()

   case ("test_v1comp", "test_v1complete")
     call wrtout(std_out," Testing symmetries (assuming overcomplete DVDB, pass extra argument to dump v1(r)) to file")
     call get_command_argument(2, dvdb_filepath)
     ABI_CHECK(get_arg("symv1scf", symv1scf, msg, default=0) == 0, msg)
     ABI_CHECK(get_arg("potfile", dump_file, msg, default="") == 0, msg)
     call dvdb_test_v1complete(dvdb_filepath, symv1scf, dump_file, comm)

   case ("test_v1rsym")
     call wrtout(std_out," Testing symmetries of V1(r) in real space.")
     call get_command_argument(2, dvdb_filepath)
     ABI_CHECK(get_arg("symv1scf", symv1scf, msg, default=0) == 0, msg)
     call dvdb_test_v1rsym(dvdb_filepath, symv1scf, comm)

   case ("test_ftinterp")
     call get_command_argument(2, dvdb_filepath)
     ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, default=2, want_len=3) == 0, msg)
     ABI_CHECK(get_arg("ddb-path", ddb_filepath, msg, default="") == 0, msg)
     ABI_CHECK(get_arg("rspace_cell", rspace_cell, msg, default=0) == 0, msg)
     ABI_CHECK(get_arg("symv1scf", symv1scf, msg, default=0) == 0, msg)
     ABI_CHECK(get_arg("dvdb-add-lr", dvdb_add_lr, msg, default=1) == 0, msg)
     ABI_CHECK(get_arg("qdamp", dvdb_qdamp, msg, default=0.1_dp) == 0, msg)
     ABI_CHECK(get_arg_list("coarse-ngqpt", coarse_ngqpt, lenr, msg, default=0, want_len=3) == 0, msg)
     call dvdb_test_ftinterp(dvdb_filepath, rspace_cell, symv1scf, ngqpt, dvdb_add_lr, dvdb_qdamp, &
                             ddb_filepath, prtvol, coarse_ngqpt, comm)

   case ("downsample")
     call get_command_argument(2, dvdb_filepath)
     call get_command_argument(3, dump_file)
     ABI_CHECK(get_arg_list("ngqpt", ngqpt, lenr, msg, want_len=3) == 0, msg)
     write(std_out,"(a)")sjoin(" Downsampling q-mesh with ngqpt:", ltoa(ngqpt))
     write(std_out,"(a)")trim(dvdb_filepath), " --> ", trim(dump_file)

     dvdb = dvdb_new(dvdb_filepath, xmpi_comm_self)
     call ngfft_seq(ngfftf, dvdb%ngfft3_v1(:, 1))
     call dvdb%open_read(ngfftf, xmpi_comm_self)
     if (prtvol > 0) call dvdb%print(prtvol=prtvol)
     call dvdb%list_perts([-1,-1,-1], npert_miss, unit=std_out)
     call dvdb%qdownsample(dump_file, ngqpt, comm)
     call dvdb%free()

   !case ("convert")
   !  call get_command_argument(2, dvdb_filepath)
   !  call get_command_argument(3, dump_file)
   !  call dvdb_convert_fort2nc(dvdb_filepath, dump_file, comm)

   case default
     MSG_ERROR(sjoin("Unknown command:", command))
   end select

 end if

 call wrtout(std_out," Done")
 call abinit_doctor("__mrgdv")

 100 call xmpi_end()

 end program mrgdv
!!***

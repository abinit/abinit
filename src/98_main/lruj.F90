!!****p* ABINIT/lruj
!! NAME
!! lruj
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2022 ABINIT group ()
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


program lruj

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_build_info
 use m_errors
 use m_argparse
 use m_crystal
 use netcdf
 use m_nctk

 use m_fstrings,    only : itoa, sjoin, ltoa
 !use defs_abitypes, only : MPI_type
 use m_specialmsg,  only : specialmsg_getcount, herald
 !use m_io_tools,    only : open_file
 use m_sort,        only : sort_dp
 !use m_parser,      only : intagm, parsefile
 use m_common,      only : crystal_from_file
 use m_mpinfo,      only : destroy_mpi_enreg, initmpi_seq
 !use m_dtfil,       only : isfile
 use m_paw_uj,      only : pawuj_ini,pawuj_free,pawuj_det, macro_uj_type

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter     :: master=0 !ndtpawuj=4,nwfchr=6,
 integer               :: ii,jdtset,lenstr,marr,ndtset,tread,comm, ncid, nnat, prtvol
 integer               :: nat,nproc,my_rank,nspden,macro_uj,pawujat,pawprtvol, nargs, nfiles, ndtpawuj
 real(dp)              :: dfact
 logical               :: iam_master
 character(len=24)     :: codename
 !character(len=strlen) :: string
 !character(len=500)    :: message
 !type(MPI_type)        :: mpi_enreg
 real(dp)              :: ures
 type(crystal_t) :: cryst
 character(len=500) :: command, arg, msg
!arrays
 character(len=fnlen),allocatable :: file_paths(:)
 integer,allocatable   :: iperm(:), pawujat_file(:)
 real(dp),allocatable  :: dpdumar(:),dprarr(:), uj_perts(:), luocc(:,:), luocc_nnat(:,:)

! *********************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI (not used but necessary)
 call xmpi_init()
 comm = xmpi_world
 nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm); iam_master = (my_rank == master)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 ! Default for sequential use
 !call initmpi_seq(mpi_enreg)

 if (my_rank /= master) goto 100

 ! Command line options.
 ! Syntax:

 !   lruj  FILE1 FILE2 FILE3 ... [-d 5]

 ! i.e. list of netcdf files come first, followed by options
 nargs = command_argument_count()
 ABI_MALLOC(file_paths, (nargs))
 nfiles = 0
 do ii=1,nargs
   call get_command_argument(ii, arg)
   if (arg(1:1) == "-") exit
   nfiles = nfiles + 1
   file_paths(nfiles) = trim(arg)
 end do

 if (nfiles == 0) then
   write(std_out, *) "Empty file list!"
   goto 100
 end if

 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     call lruj_show_help()
     goto 100

   else
     ! NOP
     !nfiles = nfiles + 1
     !file_paths(nfiles) = trim(arg)
   end if
 end do

 ! Get other options from the CLI. e.g. --prtvol 0 -d 3.0
 ! Should be documented in lruj_show_help
 ABI_CHECK(get_arg("prtvol", prtvol, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("d", dfact, msg, default=5.0_dp) == 0, msg)

 ! Print header
 codename='LRUJ'//repeat(' ',18)
 call herald(codename, abinit_version, std_out)

 ! Read uj_pert from each file and sort file_paths accordingly.
 ABI_MALLOC(uj_perts, (nfiles))

 do ii=1,nfiles
   NCF_CHECK(nctk_open_read(ncid, file_paths(ii), xmpi_comm_self))

   NCF_CHECK(nf90_get_var(ncid, vid("uj_pert"), uj_perts(ii)))

   ! Make sure ndtpawuj is always 4.
   NCF_CHECK(nctk_get_dim(ncid, "ndtpawuj", ndtpawuj))
   ABI_CHECK_IEQ(ndtpawuj, 4, "Wrong ndtpawuj")

   NCF_CHECK(nctk_get_dim(ncid, "nnat", nnat))
   NCF_CHECK(nf90_close(ncid))
 end do

 ! Comment this section if you don't need to sort files by uj_pert
 ABI_MALLOC(iperm, (nfiles))
 iperm = [(ii, ii=1,nfiles)]
 call sort_dp(nfiles, uj_perts, iperm, tol12)
 file_paths(1:nfiles) = file_paths(iperm(:))
 ABI_FREE(iperm)

 ! Now read data and perform basic consistency check
 ABI_MALLOC(luocc, (ndtpawuj, nfiles))
 ABI_MALLOC(luocc_nnat, (ndtpawuj, nnat))
 ABI_MALLOC(pawujat_file, (nfiles))

 do ii=1,nfiles
   NCF_CHECK(nctk_open_read(ncid, file_paths(ii), xmpi_comm_self))

   NCF_CHECK(nf90_get_var(ncid, vid("pawujat"), pawujat_file(ii)))
   NCF_CHECK(nf90_get_var(ncid, vid("luocc"), luocc_nnat))
   luocc(:,ii) = luocc_nnat(:, pawujat)
   write(std_out, *) "luocc for file:", ii, "with uj_pert", uj_perts(ii), "(Ha)", ch10, luocc(:,ii)
   NCF_CHECK(nf90_close(ncid))
 end do

 ! Consistency check.
 ! Very basic. Additional checks might be added.
 if (any(pawujat_file /= pawujat_file(1))) then
   ABI_ERROR(sjoin("Found differerent values of pawujat in files:", ltoa(pawujat_file)))
 end if

 ABI_FREE(luocc_nnat)
 ABI_FREE(pawujat_file)

 !call cryst%print()
 !call cryst%free()

 ABI_SFREE(luocc)
 ABI_SFREE(uj_perts)
 ABI_SFREE(file_paths)

 !call destroy_mpi_enreg(mpi_enreg)

 ! Write information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor("__lruj")

 ! Close files
 if (iam_master) close(ab_out)

 100 call xmpi_end()

contains

!!****f* abitk/lruj_show_help
!! NAME
!! lruj_show_help
!!
!! FUNCTION
!!  Show command line help
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine lruj_show_help()

  write(std_out,"(a)")" --version              Show version number and exit."
  write(std_out,"(a)")" -h, --help             Show this help and exit."
  write(std_out,"(a)")" -v                     Increase verbosity level"

  !write(std_out,"(2a)")ch10,"=== HEADER ==="
  !write(std_out,"(a)")"hdr_print FILE [--prtvol 0]  Print ABINIT header."

  !write(std_out,"(2a)")ch10,"=== KPOINTS ==="
  !write(std_out,"(a)")"ibz FILE --ngkpt 2 2 2 or --kptrlatt [--kptopt 1] [--shiftk 0.5 0.5 0.5] [--chksymbreak 1]"

end subroutine lruj_show_help
!!***

 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
end function vid

end program lruj
!!***

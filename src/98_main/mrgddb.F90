!!****p* ABINIT/mrgddb
!! NAME
!! mrgddb
!!
!! FUNCTION
!! This code merges the derivative databases.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2025 ABINIT group (DCA, XG, GMR, SP, GA)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! NOTES
!! The heading of the constituted database is read,
!! then the heading of the temporary database to be added is read,
!! the code check their compatibility, and create a new
!! database that mixes the old and the temporary ones.
!! This process can be iterated.
!! The whole database will be stored in
!! central memory. One could introduce a third mode in which
!! only the temporary DDB is in central memory, while the
!! input DDB is read twice: first to make a table of blocks,
!! counting the final number of blocks, and second to merge
!! the two DDBs. This would save memory.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program mrgddb

 use defs_basis
 use m_abicore
 use m_errors
 use m_xmpi
 use m_ddb_hdr
 use m_ddb

 use m_build_info,   only : abinit_version
 use m_specialmsg,   only : herald
 use m_time ,        only : asctime, timein
 use m_io_tools,     only : file_exists
 use m_fstrings,     only : sjoin

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: ddbun=2
 integer :: chkopt,ios, mddb
 integer :: iddb,ii,nddb,nfiles_cli,nargs,comm,my_rank
 real(dp) :: tcpu,tcpui,twall,twalli
 logical :: cannot_overwrite=.True.
 character(len=24) :: codename
 character(len=fnlen) :: dscrpt, outname
!arrays
 real(dp) :: tsec(2)
 character(len=fnlen),allocatable :: filnam(:),copy_filnam(:)
 character(len=500) :: msg,arg

!******************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 call timein(tcpui,twalli)

 codename='MRGDDB'//repeat(' ',18)
 call herald(codename,abinit_version,std_out)

 ABI_CHECK(xmpi_comm_size(comm)==1, "mrgddb not programmed for parallel execution")

 nargs = command_argument_count()

 mddb = 5000 ! maximum number of databases (initial guess)

 chkopt = 1; nfiles_cli = 0
 do ii=1,nargs
   call get_command_argument(ii, arg)
   if (arg == "-v" .or. arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "--nostrict") then
     ! Disable consistency checks
     chkopt = 0

   else if (arg == "-f") then
     cannot_overwrite = .False.

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     write(std_out,*)"Usage:"
     write(std_out,*)"    mrgddb                           Interactive prompt."
     write(std_out,*)"    mrgddb < run.files               Read arguments from run.files."
     write(std_out,*)"    mrgddb out_DDB in1_DDB in2_DDB   Merge list of input DDB files, produce new out_DDB file."
     write(std_out,*)"    mrgddb out_DDB in*_DDB           Same as above but use shell wildcards instead of file list."
     write(std_out,*)"    mrgddb out_DDB in_DDB.nc         Convert DDB from NetCDF format to plain text format."
     write(std_out,*)"    mrgddb out_DDB.nc in_DDB         Convert DDB from plain text to NetCDF format."
     write(std_out,*)" "
     write(std_out,*)"Available options:"
     write(std_out,*)"    -v, --version      Show version number and exit."
     write(std_out,*)"    -f                 Overwrite output DDB if file already exists."
     write(std_out,*)"    --nostrict         Disable consistency checks"
     write(std_out,*)"    -h, --help         Show this help and exit."
     goto 100

   else
     if (ii == 1) then
       ! First arg is the name of the output file.
       outname = trim(arg); cycle
     end if

     ! Save filenames passed via command-line.
     if (.not. allocated(filnam)) then
       ABI_MALLOC(filnam, (mddb+1))
     end if
     nfiles_cli = nfiles_cli + 1
     if (nfiles_cli > mddb + 1) then
       ! Extend filnam
       ABI_MALLOC(copy_filnam, (mddb+1))
       copy_filnam = filnam
       iddb = mddb + 1
       mddb = 2 * mddb
       ABI_FREE(filnam)
       ABI_MALLOC(filnam, (mddb+1))
       filnam(:iddb) = copy_filnam(:iddb)
       ABI_FREE(copy_filnam)
     end if
     filnam(nfiles_cli) = trim(arg)
   end if
 end do ! nargs

 if (nfiles_cli == 0) then
   ! Read names of files, operating mode (also check its value),
   ! and short description of new database.

   ! Read the name of the output ddb
   write(std_out,*)' Give name for output derivative database : '
   read(std_in, '(a)' ) outname
   write(std_out,'(2a)' )' ',trim(outname)

   ! Read the description of the derivative database
   write(std_out,*)' Give short description of the derivative database :'
   read(std_in, '(a)' )dscrpt
   write(std_out,'(2a)' )' ',trim(dscrpt)

   ! Read the number of input ddbs, and check its value
   ! MG NOTE: In the documentation of mrgddb_init I found:
   !
   ! nddb = (=1 => will initialize the ddb, using an input GS file)
   !        (>1 => will merge the whole set of ddbs listed)
   !    if nddb==1,
   !     (2) Formatted input file for the Corning ground-state code
   !
   ! but the case nddb=1 with input file is not supported anymore!

   write(std_out,*)' Give number of input ddbs'
   read(std_in,*)nddb
   write(std_out,*)nddb
   ABI_MALLOC(filnam, (nddb))

   ! Read the file names
   do iddb=1,nddb
     !Added to catch error message if the number of input ddbs is greater than the
     !actual number of ddb files entered by the user.
     read(std_in, '(a)',IOSTAT =ios ) filnam(iddb)
     if (ios < 0) then
       write(msg, '(a,i0,4a)' )&
       'The number of input ddb files: ',nddb,' exceeds the number ',&
       'of ddb file names.', ch10, &
       'Action: change the number of ddb files in the mrgddb input file.'
       ABI_ERROR(msg)
     else
       write(std_out,*)' Give name for derivative database number',iddb,' : '
       write(std_out,'(2a)' )' ',trim(filnam(iddb))
     end if
   end do

 else
   if (cannot_overwrite .and. file_exists(outname)) then
     ABI_ERROR(sjoin("Cannot overwrite existing file:", outname))
   end if
   nddb = nfiles_cli
   dscrpt = sjoin("Generated by mrgddb on:", asctime())
 end if

 ! Call the main merging routine
 call merge_ddb(nddb, filnam, outname, dscrpt, chkopt)

 ABI_FREE(filnam)

!**********************************************************************

 call timein(tcpu,twall)

 tsec(1)=tcpu-tcpui
 tsec(2)=twall-twalli

 write(std_out, '(3a,f13.1,a,f13.1)' ) '-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 call wrtout(std_out,'+mrgddb : the run completed successfully ','COLL', do_flush=.True.)

 call abinit_doctor("__mrgddb")

 100 call xmpi_end()

 end program mrgddb
!!***

!!****p* ABINIT/ioprof
!! NAME
!! ioprof
!!
!! FUNCTION
!! Tool for profiling and and testing the IO routines used in abinit
!!
!! COPYRIGHT
!! Copyright (C) 2004-2025 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program ioprof

 use defs_basis
 use m_errors
 use m_xmpi
 use m_wfk
 use m_abicore
 use m_hdr
 use netcdf
 use m_nctk

 use m_build_info,     only : abinit_version
 use m_specialmsg,     only : specialmsg_getcount, herald
 use m_fstrings,       only : lower, sjoin, itoa
 use m_io_tools,       only : delete_file, file_exists, iomode_from_fname, get_unit
 use m_argparse,       only : get_arg !, get_arg_list, get_start_step_num

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0, MAX_NFILES=50
 integer :: comm,my_rank, rank, nprocs,iomode,formeig,ierr, test_ierr, fform, ncid, ncerr
 integer :: ii,io,check_iomode,method,feg,ount,nband2read,nargs, abimem_level, nkibz, nband
 logical,parameter :: verbose=.FALSE.
 logical :: independent
 real(dp) :: abimem_limit_mb
 character(len=24) :: codename
 character(len=500) :: msg,command,arg
 character(len=fnlen) :: new_filename, wfk_source, wfk_dest, wfk_path
 type(hdr_type) :: hdr
 type(wfk_t) :: wfk
!arrays
 integer,parameter :: formeigs(2) = [0,1]
 integer,parameter :: io_modes(1) = [IO_MODE_FORTRAN]
 !integer,parameter :: io_modes(1) = [IO_MODE_MPI]
 !integer,parameter :: io_modes(1) = [IO_MODE_ETSF]
 !integer,parameter :: io_modes(2) = [IO_MODE_FORTRAN, IO_MODE_MPI]
 !integer,parameter :: io_modes(3) = [IO_MODE_FORTRAN, IO_MODE_MPI, IO_MODE_ETSF]
 integer :: start3(3), count3(3)
 real(dp) :: cwtimes(2)
 real(dp),allocatable :: waves(:,:), read_waves(:,:)
 type(kvars_t),allocatable :: Kvars(:)
 character(len=fnlen) :: hdr_fnames(MAX_NFILES) = ABI_NOFILE
! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
 ABI_CHECK(get_arg("abimem-level", abimem_level, msg, default=0) == 0, msg)
 ABI_CHECK(get_arg("abimem-limit-mb", abimem_limit_mb, msg, default=20.0_dp) == 0, msg)
#ifdef HAVE_MEM_PROFILING
 call abimem_init(abimem_level, limit_mb=abimem_limit_mb)
#endif

 codename='IOPROF'//REPEAT(' ',17)
 call herald(codename, abinit_version, std_out)

 ount = dev_null; if (verbose) ount = std_out

 nargs = command_argument_count()

 ! Command line options.
 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "-v" .or. arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     write(std_out,*)"-v, --version                 Show version number and exit."
     write(std_out,*)"info WFK_FILE                 Print info on WFK file"
     write(std_out,*)"nc2fort IN_WFK.nc OUT_WFK     Convert IN_WFK ncfile to Fortran OUT_WFK file"
     write(std_out,*)"nc_tests                      Perform basis tests of netcdf API with MPI-IO"
     write(std_out,*)"-h, --help                    Show this help and exit."
     goto 100
   end if
 end do

 call get_command_argument(1, command)

 select case (command)

 case ("info")
   ! Print info
   ABI_CHECK(nargs == 2, "Usage: info wfk_file")
   call get_command_argument(2, wfk_path)
   if (my_rank == master) then
     formeig = 0
     call wfk%open_read(wfk_path, formeig, iomode_from_fname(wfk_path), get_unit(), xmpi_comm_self)
     call wfk%print()
     call wfk%close()
   end if

 case ("nc2fort")
   ! Converter: netcdf --> Fortran
   if (my_rank /= master) goto 100
   ABI_CHECK(nargs == 3, "Usage: nc2form IN_WFK.nc OUT_WFK")
   call get_command_argument(2, wfk_source)
   call get_command_argument(3, wfk_dest)

   if (file_exists(wfk_dest)) then
     ABI_ERROR(sjoin("Cannot overwrite:", wfk_dest, ". Remove file and rerun"))
   end if
   call wfk_nc2fort(wfk_source, wfk_dest)

 case ("bench")
   ! Profile
   ABI_CHECK(nargs > 1, "Usage: bench WFK_FILE")
   call get_command_argument(2, wfk_path)
   !if (my_rank /= master) goto 100
   ! To create a new file.
   !call wfk_create_wfkfile(new_filename, hdr, iomode, formeig, Kvars, cwtimes, comm)

   formeig = 0; nband2read = 1500
   call wfk_prof(wfk_path, formeig, nband2read, comm)

 case ("nc_tests")

   call nctk_test_mpiio(print_warning=.True.)
   if (nctk_has_mpiio) then
     call wrtout(std_out, " netcdf library supports MPI-IO: [OK] ")
   else
     call wrtout(std_out, " netcdf library does not support MPI-IO: [FAILED]")
     call xmpi_abort()
   end if

   new_filename = "__IOPROF__.nc"
   NCF_CHECK(nctk_open_create(ncid, new_filename, comm))

   ! Define dimensions.
   nkibz = nprocs; nband = 101

   ncerr = nctk_def_dims(ncid, [ &
      nctkdim_t("nkibz", nkibz), &
      nctkdim_t("nband", nband) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Define arrays.
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("waves", "dp", "two, nband, nkibz"), &
     nctkarr_t("compressed_waves", "dp", "two, nband, nkibz"), &
     nctkarr_t("collective_waves", "dp", "two, nband, nkibz"), &
     nctkarr_t("collective_compressed_waves", "dp", "two, nband, nkibz") &
   ])

   ! Compress waves to reduce size on disk.
   NCF_CHECK(nf90_def_var_deflate(ncid, vid("compressed_waves"), shuffle=1, deflate=1, deflate_level=5))

   ! Begin writing.
   NCF_CHECK(nctk_set_datamode(ncid))

   ABI_CALLOC(read_waves, (2, nband))
   ABI_CALLOC(waves, (2, nband))
   waves(1, :) = [(one * ii, ii=1,nband)]
   waves(2, :) = [(ten * ii, ii=1,nband)]
   waves = waves * (my_rank + 1)
   start3 = [1,1,my_rank+1]; count3 = [2,nband,1]

   call wrtout(std_out, "Writing data in indipendent mode...")
   NCF_CHECK(nctk_set_collective(ncid, vid("waves"), independent=.True.))
   do rank=0, nprocs-1
     if (my_rank == rank) then
       call wrtout(std_out, sjoin("MPI rank:", itoa(my_rank)))
       NCF_CHECK(nf90_put_var(ncid, vid("waves"), waves, start=start3, count=count3))
     end if
     call xmpi_barrier(comm)
   end do

   ! It seems indipendent mode cannot be used with deflated variables.
   !NCF_CHECK(nctk_set_collective(ncid, vid("compressed_waves"), independent=.True.))
   call wrtout(std_out, "Writing compressed data in indipendent mode...")
   NCF_CHECK(nf90_put_var(ncid, vid("compressed_waves"), waves, start=start3, count=count3))

   ! Write variables in collective mode.
   call wrtout(std_out, "Writing data in collective mode...")
   independent = .False.
   NCF_CHECK(nctk_set_collective(ncid, vid("collective_waves"), independent=independent))
   NCF_CHECK(nf90_put_var(ncid, vid("collective_waves"), waves, start=start3, count=count3))
   NCF_CHECK(nctk_set_collective(ncid, vid("collective_compressed_waves"), independent=independent))
   NCF_CHECK(nf90_put_var(ncid, vid("collective_compressed_waves"), waves, start=start3, count=count3))
   NCF_CHECK(nf90_close(ncid))
   call xmpi_barrier(comm)

   ! Now read the data and performs consistency check.
   NCF_CHECK(nctk_open_read(ncid, new_filename, comm))
   test_ierr = 0

   call wrtout(std_out, "Reading data in indipendent mode...")
   NCF_CHECK(nctk_set_collective(ncid, vid("waves"), independent=.True.))
   do rank=0, nprocs-1
     if (my_rank == rank) then
       call wrtout(std_out, sjoin("MPI rank:", itoa(my_rank)))
       NCF_CHECK(nf90_get_var(ncid, vid("waves"), read_waves, start=start3, count=count3))
       call check_waves("Reading waves", test_ierr)
       NCF_CHECK(nf90_get_var(ncid, vid("compressed_waves"), read_waves, start=start3, count=count3))
       call check_waves("Reading compressed waves", test_ierr)
     end if
     call xmpi_barrier(comm)
   end do

   call wrtout(std_out, "Reading compressed data in indipendent mode...")
   NCF_CHECK(nf90_get_var(ncid, vid("compressed_waves"), read_waves, start=start3, count=count3))
   call check_waves("Reading compressed waves", test_ierr)

   call wrtout(std_out, "Reading data in collective mode...")
   NCF_CHECK(nctk_set_collective(ncid, vid("collective_waves"), independent=independent))
   NCF_CHECK(nf90_get_var(ncid, vid("collective_waves"), read_waves, start=start3, count=count3))
   call check_waves("Reading collective_waves", test_ierr)

   call wrtout(std_out, "Reading compressed data in collective mode...")
   NCF_CHECK(nctk_set_collective(ncid, vid("collective_compressed_waves"), independent=independent))
   NCF_CHECK(nf90_get_var(ncid, vid("collective_compressed_waves"), read_waves, start=start3, count=count3))
   call check_waves("Reading collective_compressed waves", test_ierr)

   NCF_CHECK(nf90_close(ncid))
   ABI_FREE(waves)
   ABI_FREE(read_waves)
   if (my_rank == master) call delete_file(new_filename, ierr)

   ABI_CHECK_IEQ(test_ierr, 0, "test_ierr /= 0")

 case ("wfk_tests")
   ! Unit tests for WFK file.
   if (my_rank /= master) goto 100

   do ii=1,count(hdr_fnames/=ABI_NOFILE)

     ! Read the header from an external netcdf file
     ! This trick is needed because the initialization of
     ! the header is *IMPOSSIBLE* if we don't have a full ABINIT input file!
     call hdr%from_fname(hdr_fnames(ii), fform, comm)
     ABI_CHECK(fform /= 0, "fform==0")

     call hdr%echo(fform,4,unit=std_out)

     do feg=1,size(formeigs)
       formeig = formeigs(feg)
       do io=1,size(io_modes)
         iomode = io_modes(io)
         new_filename = "NEW_WFK"
         if (iomode == IO_MODE_ETSF) new_filename = "NEW_WFK.nc"
         if (file_exists(new_filename)) call delete_file(new_filename,ierr)

         ! TODO
         if (formeig==1 .and. iomode==IO_MODE_ETSF) then
           ABI_WARNING("iomode==1 with ETSF_IO not coded yet")
           CYCLE
         end if

         ABI_MALLOC(Kvars, (hdr%nkpt))

         if (my_rank == master) then
           call wrtout(ount, "Calling wfk_create_wfkfile")
           call wfk_create_wfkfile(new_filename,hdr,iomode,formeig,Kvars,cwtimes,xmpi_comm_self)
           !call wfk_create_wfkfile(new_filename,hdr,IO_MODE_FORTRAN,formeig,Kvars,cwtimes,xmpi_comm_self)
           call wrtout(ount, "Done wfk_create_wfkfile")
         end if

         if (nprocs > 1) then
           ABI_ERROR("Not coded")
           !call kvars_mpicast(Kvars,master,comm)
         end if

         method = 0
         call wrtout(std_out,sjoin("Checking file:", new_filename, ", with iomode=",itoa(iomode)))
         call wfk_check_wfkfile(new_filename,hdr,iomode,method,formeig,Kvars,cwtimes,comm,ierr)

         write(msg,'(2(a,i0),2(a,f8.2))')" Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
         call wrtout(ount, msg)

         ABI_CHECK(ierr == 0, sjoin("wfk_check_wfkfile returned ierr ",itoa(ierr)))

         ! If not a netcdf file, try to read the file with the other mode that is compatible with it.
         check_iomode = -100
         if (iomode == IO_MODE_FORTRAN) check_iomode = IO_MODE_MPI
         if (iomode == IO_MODE_MPI)     check_iomode = IO_MODE_FORTRAN

         if (check_iomode /= -100) then
           call wrtout(std_out,sjoin("Trying to read file:",new_filename,", with check_iomode=",itoa(check_iomode)))
           call wfk_check_wfkfile(new_filename,hdr,check_iomode,method,formeig,Kvars,cwtimes,comm,ierr)

           write(msg,'(2(a,i0),2(a,f8.2))')" Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
           call wrtout(ount, msg)
           ABI_CHECK(ierr == 0, sjoin("wfk_check_wfkfile returned ierr: ", itoa(ierr)))
         end if

         ABI_FREE(Kvars)
       end do ! iomode
     end do ! formeig

     call hdr%free()
   end do

 case default
   ABI_ERROR(sjoin("Unknown command:", command))
 end select

 !call wrtout(std_out, ch10//" Analysis completed.)

 call abinit_doctor("__ioprof")
 100 call xmpi_end()

contains
 subroutine check_waves(header, ierr)
   character(len=*),intent(in) :: header
   integer,intent(inout) :: ierr
   if (all(waves == read_waves)) then
     write(std_out, "(1x,2a,i0,a)")trim(header), ", MPI rank: ", my_rank, " [OK] "
   else
     write(std_out, "(1x,2a,i0,a)")trim(header), ", MPI rank: ", my_rank, " [FAILED] "
     ierr = ierr + 1
   end if
 end subroutine check_waves

 integer function vid(var_name)
   character(len=*),intent(in) :: var_name
   vid = nctk_idname(ncid, var_name)
 end function vid

 end program ioprof
!!***

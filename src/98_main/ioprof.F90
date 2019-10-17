!{\src2tex{textfont=tt}}
!!****p* ABINIT/ioprof
!! NAME
!! ioprof
!!
!! FUNCTION
!! Tool for frofiling and and testing the IO routines used in abinit
!!
!! COPYRIGHT
!! Copyright (C) 2004-2019 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  (main program)
!!
!! NOTES
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,delete_file
!!      get_command_argument,hdr_echo,hdr_free,hdr_read_from_fname,herald
!!      wfk_check_wfkfile,wfk_close,wfk_create_wfkfile,wfk_nc2fort
!!      wfk_open_read,wfk_print,wfk_prof,wrtout,xmpi_init
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program ioprof

 use defs_basis
 use m_build_info
 use m_errors
 use m_xmpi
 use m_wfk
 use m_abicore
 use m_hdr
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_specialmsg,     only : specialmsg_getcount, herald
 use m_fstrings,       only : lower, sjoin, itoa
 use m_io_tools,       only : delete_file, file_exists, iomode_from_fname, get_unit

 implicit none

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0,MAX_NFILES=50
 integer :: comm,my_rank,nprocs,iomode,formeig,ierr,fform
 integer :: ii,io,check_iomode,method,feg,ount,nband2read,nargs
 logical :: verbose=.FALSE.
 character(len=24) :: codename
 character(len=500) :: msg,command,arg
 character(len=fnlen) :: new_fname,wfk_source,wfk_dest,wfk_path
 type(hdr_type) :: hdr
 type(wfk_t) :: wfk
!arrays
 integer,parameter :: formeigs(2) = (/0,1/)
 integer,parameter :: io_modes(1) = (/IO_MODE_FORTRAN/)
 !integer,parameter :: io_modes(1) = (/IO_MODE_MPI/)
 !integer,parameter :: io_modes(1) = (/IO_MODE_ETSF/)
 !integer,parameter :: io_modes(2) = (/IO_MODE_FORTRAN, IO_MODE_MPI/)
 !integer,parameter :: io_modes(3) = (/IO_MODE_FORTRAN, IO_MODE_MPI, IO_MODE_ETSF/)
 real(dp) :: cwtimes(2)
 type(kvars_t),allocatable :: Kvars(:)
 character(len=fnlen) :: hdr_fnames(MAX_NFILES) = ABI_NOFILE

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 ! Initialize MPI
 call xmpi_init()
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

!Initialize memory profiling if it is activated
!if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
!note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 codename='ABIWFK'//REPEAT(' ',17)
 call herald(codename,abinit_version,std_out)

 !verbose = .TRUE.
 ount = dev_null; if (verbose) ount = std_out

 if (nprocs>1) then
   MSG_ERROR("not programmed for parallel execution, Run it in sequential")
 end if

 nargs = command_argument_count()

 ! Command line options.
 do ii=1,command_argument_count()
   call get_command_argument(ii, arg)
   if (arg == "-v" .or. arg == "--version") then
     write(std_out,"(a)") trim(abinit_version); goto 100

   else if (arg == "-h" .or. arg == "--help") then
     ! Document the options.
     write(std_out,*)"-v, --version                 Show version number and exit."
     write(std_out,*)"info wfk_file                 Print info on WFK file"
     write(std_out,*)"nc2fort ncfile fortran_file   Convert WFK ncfile to Fortran WFK file"
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
     call wfk_open_read(wfk, wfk_path, formeig, iomode_from_fname(wfk_path), get_unit(), xmpi_comm_self)
     call wfk%print()
     call wfk%close()
   end if

 case ("nc2fort")
   ! Converter: netcdf --> Fortran
   ABI_CHECK(nargs == 3, "Usage: nc2form wfk_source wfk_dest")
   call get_command_argument(2, wfk_source)
   call get_command_argument(3, wfk_dest)

   if (file_exists(wfk_dest)) then
     MSG_ERROR(sjoin("Cannot overwrite:", wfk_dest, ". Remove file and rerun"))
   end if

   if (my_rank == master) call wfk_nc2fort(wfk_source, wfk_dest)

 case ("prof")
   ! Profile
   ABI_CHECK(nargs > 1, "Usage: bench wfk_path")
   call get_command_argument(2, wfk_path)
   formeig = 0
   nband2read = 1500
   call wfk_prof(wfk_path, formeig, nband2read, comm)
   ! TO CREATE new file
   !call wfk_create_wfkfile(new_fname,hdr,iomode,formeig,Kvars,cwtimes,xmpi_comm_self)

 case ("unittests")
   ! Unit tests.
   do ii=1,count(hdr_fnames/=ABI_NOFILE)

     ! Read the header from an external netcdf file
     ! This trick is needed because the initialization of
     ! the header is *IMPOSSIBLE* if we don't have a full ABINIT input file!
     !
     call hdr_read_from_fname(hdr,hdr_fnames(ii),fform,comm)
     ABI_CHECK(fform/=0,"fform==0")

     call hdr%echo(fform,4,unit=std_out)

     do feg=1,size(formeigs)
       formeig = formeigs(feg)
       do io=1,size(io_modes)
         iomode = io_modes(io)

         new_fname = "NEW_WFK"
         if (iomode==IO_MODE_ETSF) new_fname = "NEW_WFK.nc"

         if (file_exists(new_fname)) then
           call delete_file(new_fname,ierr)
         end if

         ! TODO
         if (formeig==1 .and. iomode==IO_MODE_ETSF) then
           MSG_WARNING("iomode==1 with ETSF_IO not coded yet")
           CYCLE
         end if

         ABI_DT_MALLOC(Kvars, (hdr%nkpt))

         if (my_rank == master) then
           call wrtout(ount,"Calling wfk_create_wfkfile","COLL")
           call wfk_create_wfkfile(new_fname,hdr,iomode,formeig,Kvars,cwtimes,xmpi_comm_self)
           !call wfk_create_wfkfile(new_fname,hdr,IO_MODE_FORTRAN,formeig,Kvars,cwtimes,xmpi_comm_self)
           call wrtout(ount,"Done wfk_create_wfkfile","COLL")
         end if

         if (nprocs > 1) then
           MSG_ERROR("Not coded")
           !call kvars_mpicast(Kvars,master,comm)
         end if

         method = 0

         call wrtout(std_out,sjoin("Checking file:",new_fname,", with iomode=",itoa(iomode)))
         call wfk_check_wfkfile(new_fname,hdr,iomode,method,formeig,Kvars,cwtimes,comm,ierr)

         write(msg,'(2(a,i2),2(a,f8.2))')&
         " Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
         call wrtout(ount,msg,"COLL")

         if (ierr/=0) then
           MSG_ERROR(sjoin("wfk_check_wfkfile returned ierr ",itoa(ierr)))
         end if

         ! If not netcdf file, try to read the file with the other mode that is compatible with it.
         check_iomode = -100
         if (iomode == IO_MODE_FORTRAN) check_iomode = IO_MODE_MPI
         if (iomode == IO_MODE_MPI)     check_iomode = IO_MODE_FORTRAN

         if (check_iomode /= -100) then
           call wrtout(std_out,sjoin("Trying to read file:",new_fname,", with check_iomode=",itoa(check_iomode)))

           call wfk_check_wfkfile(new_fname,hdr,check_iomode,method,formeig,Kvars,cwtimes,comm,ierr)

           write(msg,'(2(a,i2),2(a,f8.2))')&
           " Read with iomode: ",iomode,", nproc: ",nprocs,", cpu: ",cwtimes(1),", wall:",cwtimes(2)
           call wrtout(ount,msg,"COLL")

           if (ierr/=0) then
             MSG_ERROR(sjoin("wfk_check_wfkfile returned ierr ",itoa(ierr)))
           end if
         end if

         ABI_DT_FREE(Kvars)
       end do ! iomode
     end do ! formeig

     call hdr%free()
   end do

 case default
   MSG_ERROR(sjoin("Unknown command:", command))
 end select

 !call wrtout(std_out,ch10//" Analysis completed.","COLL")

!Writes information on file about the memory before ending mpi module, if memory profiling is enabled
 !ABI_ALLOCATE(nband, (2))
 call abinit_doctor("__ioprof")

 100 call xmpi_end()

 end program ioprof
!!***

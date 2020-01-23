!!****p* ABINIT/multibinit
!! NAME
!! multibinit
!!
!! FUNCTION
!! Main routine MULTIBINIT.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2020 ABINIT group (AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (main routine)
!!
!! OUTPUT
!!  (main routine)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abihist_bcast,abihist_free
!!      abimem_init,abinit_doctor,compute_anharmonics
!!      effective_potential_file_getdimsystem,effective_potential_file_gettype
!!      effective_potential_file_maphisttoref,effective_potential_file_read
!!      effective_potential_file_readmdfile,effective_potential_free
!!      effective_potential_setconfinement,effective_potential_writenetcdf
!!      effective_potential_writexml,fit_polynomial_coeff_fit
!!      fit_polynomial_printsystemfiles,flush_unit,herald,init10,instrng
!!      inupper,invars10,isfile,mover_effpot,multibinit_dtset_free
!!      outvars_multibinit,timein,wrtout,xmpi_bcast,xmpi_init,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program multibinit

  use defs_basis
  use defs_abitypes
  use m_build_info
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors
  use m_argparse
  use m_effective_potential
  use m_fit_polynomial_coeff
  use m_multibinit_dataset
  use m_effective_potential_file
  use m_abihist

  use m_specialmsg, only : specialmsg_getcount, herald
  use m_io_tools,   only : flush_unit, open_file
  use m_time,       only : asctime, timein
  use m_parser,     only : instrng
  use m_dtset,      only : chkvars
  use m_dtfil,      only : isfile

  !use m_generate_training_set, only : generate_training_set
  !use m_compute_anharmonics, only : compute_anharmonics
  use m_init10,              only : init10
  use m_multibinit_unittest, only: mb_test_main
  use m_multibinit_driver
  use m_multibinit_main2, only: multibinit_main2

  implicit none

  !Arguments -----------------------------------

  !Local variables-------------------------------
  ! Set array dimensions
  real(dp) :: tcpu,tcpui,twall,twalli
  real(dp) :: tsec(2)
  character(len=24) :: codename,start_datetime
  character(len=fnlen) :: filnam(18),tmpfilename
  character(len=500) :: message 
  type(args_t) :: args
  integer :: ii
  integer :: master, my_rank, comm, nproc, ierr
  logical :: iam_master

  !TEST_AM
  ! integer :: natom_sp
  ! real(dp),allocatable :: dynmat(:,:,:,:,:)
  !TEST_AM
  !******************************************************************

  !Change communicator for I/O (mandatory!)
  call abi_io_redirect(new_io_comm=xmpi_world)

  !Initialize MPI
  call xmpi_init()

  master = 0
  comm = xmpi_world
  nproc = xmpi_comm_size(comm)
  my_rank = xmpi_comm_rank(comm)
  iam_master = (my_rank == master)

  ! Parse command line arguments.
  args = args_parser(); if (args%exit /= 0) goto 100

 ! nargs = command_argument_count()
 ! do iarg=1,nargs
 !    call get_command_argument(number=iarg, value=arg)
 !    if (arg == "-v" .or. arg == "--version") then
 !       write(std_out,"(a)") trim(abinit_version); goto 100
 !       goto 100
 !    else if (arg == "--unittest") then
 !       unittest=.True.
 !       call mb_test_main()
 !       goto 100
 !    else if(arg== "-F03") then
 !       use_f03=.True.
 !    endif
 ! end do

 ! Initialize memory profiling if activated at configure time.
 ! if a full report is desired, set the argument of abimem_init to "2" instead of "0" via the command line.
 ! note that the file can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level, limit_mb=args%abimem_limit_mb)
#endif

  !Initialisation of the timing
  call timein(tcpui,twalli)

  if (iam_master) then
     codename='MULTIBINIT'//repeat(' ',14)
     call herald(codename,abinit_version,std_out)
  end if

  start_datetime = asctime()
  !Print the number of cpus in log file
  write(message,'(a,i5,a)') '-  nproc =',nproc,ch10
  call wrtout(std_out,message,'COLL')

  !Initialise the code : write heading, and read names of files.
  call init10(filnam,comm)

  !******************************************************************

  call timein(tcpu,twall)

  write(message, '(a,f11.3,a,f11.3,a)' )'-begin at tcpu',tcpu-tcpui,'  and twall',twall-twalli,' sec'
  call wrtout(std_out,message,'COLL')

  ! Open output files and ab_out (might change its name if needed)
  ! MJV 1/2010 : now output file is open, but filnam(2) continues unmodified
  ! so the other output files are overwritten instead of accumulating.
  if (iam_master) then
     tmpfilename = filnam(2)
     call isfile(tmpfilename,'new')
     if (open_file(tmpfilename,message,unit=ab_out,form="formatted",status="new",&
          &   action="write") /= 0) then
        MSG_ERROR(message)
     end if
     !  Call open_file(unit=ab_out,file=tmpfilename,form='formatted',status='new')
     rewind (unit=ab_out)
     call herald(codename,abinit_version,ab_out)
     !  Print the number of cpus in output
     write(message,'(a,i5,a)') '-  nproc =',nproc
     call wrtout(ab_out,message,'COLL')
  else
     ab_out = dev_null
  end if

  write(message, '(a,(80a),a)' ) ch10,&
       & ('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  !***************************************************************************************
  !***************************************************************************************
  if(args%multibinit_F03_mode==1) then
     ! Use the F03 mode, which has only spin and a simple harmonic lattice now
     ! After everything is migrated, it will becomes default and multibinit_main will be deprecated.
     call multibinit_main2(filnam)
  else
     call multibinit_main(filnam, args%dry_run)
  end if
  ! Final message
  !****************************************************************************************

  write(message,'(a,a,a,(80a))') ch10,('=',ii=1,80),ch10
  call wrtout(ab_out,message,'COLL')
  call wrtout(std_out,message,'COLL')

  call timein(tcpu,twall)
  tsec(1)=tcpu-tcpui
  tsec(2)=twall-twalli

  write(message, '(a,i4,a,f13.1,a,f13.1)' )' Proc.',my_rank,' individual time (sec): cpu=',&
       &   tsec(1),'  wall=',tsec(2)
  call wrtout(std_out,message,"COLL")

  if (iam_master) then
     write(ab_out, '(a,a,a,i4,a,f13.1,a,f13.1)' )'-',ch10,&
          &     '- Proc.',my_rank,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
  end if

   call xmpi_sum(tsec,comm,ierr)

   write(message, '(a,(80a),a,a,a,f11.3,a,f11.3,a,a,a,a)' ) ch10,&
&   ('=',ii=1,80),ch10,ch10,&
&   '+Total cpu time',tsec(1),&
&   '  and wall time',tsec(2),' sec',ch10,ch10,&
&   ' multibinit : the run completed succesfully.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   if (iam_master) then
   ! Write YAML document with the final summary.
   ! we use this doc to test whether the calculation is completed.
     write(std_out,"(a)")"--- !FinalSummary"
     write(std_out,"(a)")"program: multibinit"
     write(std_out,"(2a)")"version: ",trim(abinit_version)
     write(std_out,"(2a)")"start_datetime: ",start_datetime
     write(std_out,"(2a)")"end_datetime: ",asctime()
     write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
     write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
     write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
     write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   !write(std_out,"(a,i0)")"num_warnings: ",nwarning
   !write(std_out,"(a,i0)")"num_comments: ",ncomment
     write(std_out,"(a)")"..."
     call flush_unit(std_out)
   end if

!Write information on file about the memory before ending mpi module, if memory profiling is enabled
   call abinit_doctor("__multibinit")

   call flush_unit(ab_out)
   call flush_unit(std_out)

   if (iam_master) close(ab_out)

   100 call xmpi_end()

   end program multibinit
!!***

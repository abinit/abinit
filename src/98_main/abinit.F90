!!****p* ABINIT/abinit
!! NAME
!! abinit
!!
!! FUNCTION
!! Main routine for conducting Density-Functional Theory calculations or Many-Body Perturbation Theory calculations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, MKV, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! NOTES
!! The new user is strongly adviced to read the
!! latest version of the file ~abinit/doc/users/new_user_guide.html
!! before trying to modify or even use the code.
!! Even experienced users of the code should also be careful in coding,
!! please read the latest version of the file ~abinit/doc/developers/rules_coding
!!
!! The present main routine drives the following operations :
!!
!! 1) Eventually initialize MPI
!! 2) Initialize overall timing of run
!! 3) Print greeting for interactive user and
!!    Read names of files (input, output, rootinput, rootoutput, roottemporaries),
!!    create the name of the status file, initialize the status subroutine.
!! 4) Open output file and print herald at top of output and log files
!! 5) Read the input file, and store the information in a long string of characters
!! 6) Take ndtset from the input string, then allocate
!!    the arrays whose dimensions depends only on ndtset
!! 7) Continue to analyze the input string, and allocate the remaining arrays.
!!    Also modulate the timing according to timopt.
!! 8) Finish to read the "file" file completely,
!!    and also initialize pspheads (the pseudopotential header information)
!! 9) Provide defaults for the variables that have not yet been initialized.
!! 10) Perform some global initialization, depending on the value of
!! pseudopotentials, parallelism variables, or macro input variables
!! 11) Call the main input routine.
!! 12) Echo input data to output file and log file
!! 13) Perform additional checks on input data
!!  At this stage, all the information from the "files" file and "input" file
!!  have been read and checked.
!! 14) Print more information, and activate GPU
!! ___________________________________________
!! 15) Perform main calculation  (call driver)
!! -------------------------------------------
!!
!! 16) Give final echo of coordinates, etc.
!! 17) Timing analysis
!! 18) Bibliographical recommendations
!! 19) Delete the status file, and, for build-in tests, analyse the correctness of results
!! 20) Write the final timing, close the output file, and write a final line to the log file
!! 21) Eventual cleaning of MPI run
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
!!      ab7_invars_free,ab7_invars_get_abinit_vars,ab7_invars_load
!!      ab7_invars_set_flags,abi_io_redirect,abimem_init,abinit_doctor
!!      bigdft_init_errors,bigdft_init_timing_categories,chkinp,chkvars
!!      clnmpi_atom,clnmpi_grid,clnmpi_img,clnmpi_pert,date_and_time
!!      delete_file,destroy_mpi_enreg,destroy_results_out,driver,dump_config
!!      dump_cpp_options,dump_optim,f_lib_finalize,f_lib_initialize,flush_unit
!!      gather_results_out,herald,init_results_out,iofn1,libpaw_spmsg_getcount
!!      memory_eval,mpi_allreduce,mpi_setup,nctk_test_mpiio,out_acknowl,outvars
!!      outxml_finalise,outxml_open,parsefile,paw_setup_free,print_kinds
!!      setdevice_cuda,specialmsg_getcount,status,testfi,timab,timana
!!      time_set_papiopt,timein,unsetdevice_cuda,wrtout,xmpi_init
!!      xmpi_show_info,xomp_show_info,xpapi_init,xpapi_show_info,xpapi_shutdown
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program abinit

 use defs_basis
 use m_build_info
 use m_cppopts_dumper
 use m_optim_dumper
 use m_abicore
 use m_dtset
 use m_results_out
 use m_xmpi
 use m_xomp
 use m_xpapi
 use m_errors
 use m_argparse
 use m_nctk
#if defined HAVE_MPI2
 use mpi
#endif

 use defs_datatypes, only : pspheader_type
 use defs_abitypes, only : MPI_type
 use m_parser,      only : ab_dimensions
 use m_time ,       only : asctime, sec2str, timein, time_set_papiopt, timab
 use m_fstrings,    only : sjoin, strcat, itoa, yesno, ljust
 use m_io_tools,    only : flush_unit, delete_file
 use m_specialmsg,  only : specialmsg_getcount, herald
 use m_exit,        only : get_timelimit_string
 use m_atomdata,    only : znucl2symbol
 use m_libpaw_tools,only : libpaw_spmsg_getcount
 use m_mpinfo,      only : destroy_mpi_enreg, clnmpi_img, clnmpi_grid, clnmpi_atom, clnmpi_pert
 use m_memeval,     only : memory_eval
 use m_chkinp,      only : chkinp
 use m_dtfil,       only : iofn1
 use m_outxml,      only : outxml_open, outxml_finalise
 use m_out_acknowl, only : out_acknowl
 use m_timana,      only : timana
 use m_builtin_tests, only : testfi
 use m_mpi_setup,     only : mpi_setup
 use m_outvars,       only : outvars
 use m_driver,       only : driver

#ifdef HAVE_GPU_CUDA
 use m_gpu_toolbox
#endif

#if defined HAVE_BIGDFT
 use BigDFT_API,    only : bigdft_init_errors,bigdft_init_timing_categories
#endif

 use m_common, only : get_dtsets_pspheads

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

!Arguments -----------------------------------
!Local variables-------------------------------
!
!===============================================================================
!  abinit_version designate overall code version
!  mpw=maximum number of planewaves in basis sphere
!  unit numbers (ab_in,ab_out,std_out,tmp_unit) have been defined in defs_basis.f .
!  The array filnam is used for the name of input and output files,
!  and roots for generic input, output or temporary files.
!  Pseudopotential file names are set in iofn2, and are contained in pspheads.
!  The name filstat will be needed beyond gstate to check
!  the appearance of the "exit" flag, to make a hasty exit, as well as
!  in order to output the status of the computation.
!==============================================================================
! Declarations
! Define "level of the routine", for debugging purposes
 integer,parameter :: level=1
 integer :: choice,dmatpuflag,ierr,ii,iounit,ios
 integer :: lenstr,me,print_mem_report
 integer :: mu,natom,ncomment,ncomment_paw,ndtset
 integer :: ndtset_alloc,nexit,nexit_paw,nfft,nkpt,npsp
 integer :: nsppol,nwarning,nwarning_paw,prtvol,timopt,use_gpu_cuda
 integer,allocatable :: nband(:),npwtot(:)
 real(dp) :: etotal, tcpui, twalli
 real(dp) :: strten(6),tsec(2)
 real(dp),allocatable :: fred(:,:),xred(:,:)
 character(len=24) :: codename
 character(len=24) :: start_datetime
 character(len=5000) :: msg
 character(len=strlen) :: string
 character(len=fnlen) :: filstat
 character(len=fnlen) :: filnam(5)
 type(args_t) :: args
 type(dataset_type),allocatable  :: dtsets(:)
 type(MPI_type),allocatable :: mpi_enregs(:)
 type(pspheader_type),allocatable :: pspheads(:)
 type(results_out_type),allocatable,target :: results_out(:)
 type(results_out_type),pointer :: results_out_all(:)
 type(ab_dimensions) :: mx
 logical :: test_img,test_exit,use_results_all,xml_output=.false.
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=13) :: warn_fmt
#ifdef HAVE_GPU_CUDA
 integer :: gpu_devices(5)
#endif

!******************************************************************

!0) Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)
 !call xlf_set_sighandler()

!------------------------------------------------------------------------------

!1) Eventually initialize MPI. Pay attention: me and comm may be initialzed again in finddistrproc
 call xmpi_init()
 me = xmpi_comm_rank(xmpi_world)

 ! Parse command line arguments.
 args = args_parser(); if (args%exit /= 0) goto 100

 ! Initialize memory profiling if activated at configure time.
 ! if a full report is desired, set the argument of abimem_init to "2" instead of "0" via the command line.
 ! note that the file can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(args%abimem_level, limit_mb=args%abimem_limit_mb)
#endif

!------------------------------------------------------------------------------

 ! 2) Initialize overall timing of run:
 call xpapi_init()
 call xpapi_show_info(unit=std_out,mode_paral="COLL")

 start_datetime = asctime()
 call timein(tcpui,twalli)
 call timab(1,0,tsec)

 ! Start to accumulate time for the entire run. The end of accumulation is in timana.f
 call timab(1,1,tsec)

!------------------------------------------------------------------------------

!3) Print greeting for interactive user,
!read names of files (input, output, rootinput, rootoutput, roottemporaries),
!create the name of the status file, initialize the status subroutine.

 call timab(41,3,tsec)
 call iofn1(args%input_path, filnam, filstat, xmpi_world)

!------------------------------------------------------------------------------

!4) Open output file and print herald at top of output and log files

 if (me==0) then
#ifdef FC_NAG
   open(unit=ab_out,file=filnam(2),form='formatted',status='new', action="write", recl=ABI_RECL, iomsg=msg, iostat=ios)
#else
   open(unit=ab_out,file=filnam(2),form='formatted',status='new', action="write", iomsg=msg, iostat=ios)
#endif
   ABI_CHECK(ios == 0, msg)
   rewind (unit=ab_out)
   codename='ABINIT'//repeat(' ',18)
   call herald(codename,abinit_version,ab_out)
   call herald(codename,abinit_version,std_out)
   call dump_config(std_out)
   call dump_optim(std_out)
   call dump_cpp_options(std_out)
   ! Write names of files
   write(msg, '(a,a,a,a,a,a,a,a,a,a,a,a)' )&
    '- input  file    -> ',trim(filnam(1)),ch10,&
    '- output file    -> ',trim(filnam(2)),ch10,&
    '- root for input  files -> ',trim(filnam(3)),ch10,&
    '- root for output files -> ',trim(filnam(4)),ch10
   call wrtout([std_out, ab_out], msg)
 end if

 call timab(44,1,tsec)

 ! Test if the netcdf library supports MPI-IO
 call nctk_test_mpiio()

 call get_dtsets_pspheads(args%input_path, filnam(1), ndtset, lenstr, string, &
                          timopt, dtsets, pspheads, mx, dmatpuflag, xmpi_world)

 ndtset_alloc = size(dtsets) - 1
 npsp = size(pspheads)

#if defined HAVE_BIGDFT
 call f_lib_initialize()
 call bigdft_init_errors()
 call bigdft_init_timing_categories()
#endif

 ABI_MALLOC(mpi_enregs, (0:max(1,ndtset)))
 call mpi_setup(dtsets,filnam,lenstr,mpi_enregs,ndtset,ndtset_alloc,string)

 call memory_eval(dtsets,ab_out,mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads)

!------------------------------------------------------------------------------

!12) Echo input data to output file and log file

 ! For evolving variables, and results
 ABI_MALLOC(results_out, (0:ndtset_alloc))

 ! Initialize results_out datastructure
 call init_results_out(dtsets,1,1,mpi_enregs, mx%natom, mx%mband_upper, mx%nkpt,npsp,&
  mx%nsppol, mx%ntypat, results_out)

 ! Gather contributions to results_out from images of the cell, if needed
 test_img = (mx%nimage/=1.and.maxval(dtsets(:)%npimage)>1)
 use_results_all=.false.
 if (test_img) then
   use_results_all=(me==0)
   if (use_results_all) then
     ABI_MALLOC(results_out_all, (0:ndtset_alloc))
   end if

   call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all, allgather=.false.,master=0)

 else
   results_out_all => results_out
 end if

 if (me == 0) then
   ! Echo input to output file on unit ab_out, and to log file on unit 06 :
   choice=1
   do ii=1,2
     if(ii==1)iounit=ab_out
     if(ii==2)iounit=std_out

     call outvars(choice,dmatpuflag,dtsets, filnam(4), iounit, mx, ndtset,ndtset_alloc,npsp,results_out_all,timopt)
   end do

   if (dtsets(1)%prtxml == 1) then
     call outxml_open(trim(filnam(4)))
     call date_and_time(strdat,strtime,strzone,values)
     xml_output = .true.
   else
     xml_output = .false.
   end if

 end if ! me==0

 ! Clean memory
 if (test_img.and.me==0) then
   call destroy_results_out(results_out_all)
   ABI_FREE(results_out_all)
 end if

!This synchronization is not strictly needed, but without it,
!there are problems with Tv1#93 in parallel, PGI compiler, on Intel/PC
 call abi_io_redirect(new_io_comm=xmpi_world)

 call timab(44,2,tsec)

!------------------------------------------------------------------------------

!13) Perform additional checks on input data
 call timab(45,3,tsec)
 call chkinp(dtsets, ab_out, mpi_enregs, ndtset, ndtset_alloc, npsp, pspheads, xmpi_world)

 ! Check whether the string only contains valid keywords
 call chkvars(string)

!At this stage, all the information from the "files" file and "input" file have been read and checked.

!------------------------------------------------------------------------------

!14) Print more information, and activate GPU

 if (me == 0) then
   call print_kinds(std_out)     ! Printout of kinds and precisions.
   call xomp_show_info(std_out)  ! Info on the openMP environment.
   call xmpi_show_info(std_out)  ! Info on the MPI environment.
 end if

!Eventually activate GPU
 use_gpu_cuda=0
#if defined HAVE_GPU_CUDA
 gpu_devices(:)=-1
 do ii=1,ndtset_alloc
   if (dtsets(ii)%use_gpu_cuda==1) then
     use_gpu_cuda=1
     gpu_devices(:)=dtsets(ii)%gpu_devices(:)
   end if
 end do
 call setdevice_cuda(gpu_devices,use_gpu_cuda)
#endif

!------------------------------------------------------------------------------

!15) Perform main calculation
 call timab(45,2,tsec)

 test_exit=.false.
 prtvol=dtsets(1)%prtvol
 if (prtvol == -level .or. prtvol == -2 .or. args%dry_run /= 0) then
   write(msg,'(a,a,i0,a)')ch10,' abinit : before driver, prtvol=',prtvol,', debugging mode => will skip driver '
   call wrtout([std_out, ab_out], msg)
   test_exit=.true.
 end if

 if(.not.test_exit)then
   call driver(abinit_version,tcpui,dtsets,filnam,filstat, mpi_enregs,ndtset,ndtset_alloc,npsp,pspheads,results_out)
 end if

!------------------------------------------------------------------------------

 ! 16) Give final echo of coordinates, etc.
 call timab(46,1,tsec)

 write(msg,'(a,a,a,62a,80a)') ch10,'== END DATASET(S) ',('=',mu=1,62),ch10,('=',mu=1,80)
 call wrtout([std_out, ab_out], msg)

 ! Gather contributions to results_out from images of the cell, if needed
 if (test_img) then
   if (use_results_all)  then
     ABI_MALLOC(results_out_all,(0:ndtset_alloc))
   end if

   call gather_results_out(dtsets,mpi_enregs,results_out,results_out_all,use_results_all,allgather=.false.,master=0)
 end if

 if(me==0) then
   if(test_exit)then
     write(msg,'(a,a,i0,a)')ch10,' abinit : before driver, prtvol=',prtvol,', debugging mode => will skip outvars '
     call wrtout([std_out, ab_out], msg)
   else
     ! Echo input to output file on unit ab_out, and to log file on unit std_out.
     choice=2
     do ii=1,2
       if(ii==1)iounit=ab_out
       if(ii==2)iounit=std_out
       write(iounit,*)' '
       call outvars (choice,dmatpuflag,dtsets, filnam(4), iounit,mx,ndtset,ndtset_alloc,npsp,results_out_all,timopt)
       if(ii==2)write(std_out,*)' '
     end do
   end if
 end if ! me==0

 ! Clean memory
 if (test_img.and.me==0) then
   call destroy_results_out(results_out_all)
   ABI_FREE(results_out_all)
 else
   nullify(results_out_all)
 end if

 ! In prevision of the next two calls, some variables need to be transfered.
 ! They concern the case ndtset<2, and nimage=1 so take first value.
 natom=dtsets(1)%natom ; nkpt=dtsets(1)%nkpt ; nsppol=dtsets(1)%nsppol
 nfft=dtsets(1)%nfft

 ABI_MALLOC(nband,(nkpt*nsppol))
 ABI_MALLOC(npwtot,(nkpt))
 ABI_MALLOC(fred,(3,natom))
 ABI_MALLOC(xred,(3,natom))

 etotal=results_out(1)%etotal(1)
 fred(:,:)  =results_out(1)%fred(:,1:natom,1)
 nband(:)   =dtsets(1)%nband(1:nkpt*nsppol)
 npwtot(:)  =results_out(1)%npwtot(1:nkpt,1)
 strten(:)  =results_out(1)%strten(:,1)
 xred(:,:)  =results_out(1)%xred(:,1:natom,1)

 call timab(46,2,tsec)

!------------------------------------------------------------------------------

 ! 17) Timing analysis
 if(timopt/=0)then
   call timana (mpi_enregs(1), natom, nband, ndtset, nfft, nkpt, npwtot, nsppol, timopt)
 else
#if defined HAVE_MPI
   if(me==0)then ! This is for the automatic tests
     write(ab_out,'(5a)')ch10,ch10,'- Timing analysis has been suppressed with timopt=0',ch10,ch10
   end if
#endif
 end if

!------------------------------------------------------------------------------

 ! 18) Bibliographical recommendations
 if(me==0) then
   if(test_exit)then
     write(msg,'(a,a,i0,a)')ch10,' abinit : before driver, prtvol=',prtvol,', debugging mode => will skip acknowledgments'
     call wrtout([std_out, ab_out], msg)
   else
     do ii=1,2
       if(ii==1)iounit=ab_out
       if(ii==2)iounit=std_out
       call out_acknowl(dtsets,iounit,ndtset_alloc,npsp,pspheads)
     end do
   end if
 end if ! me==0

!------------------------------------------------------------------------------

 ! 19) Delete the status file, and, for build-in tests, analyse the correctness of results.
 if (ndtset == 0) then
   call testfi(dtsets(1)%builtintest,etotal,filstat,fred,natom,strten,xred)
 end if

 ! One should have here the explicit deallocation of all arrays
 call destroy_results_out(results_out)

 ABI_FREE(fred)
 ABI_FREE(nband)
 ABI_FREE(npwtot)
 ABI_FREE(results_out)
 ABI_FREE(xred)

 ! 20) Write the final timing, close the output file, and write a final line to the log file
 call timein(tsec(1),tsec(2))
 tsec(1)=tsec(1)-tcpui
 tsec(2)=tsec(2)-twalli

 ! Get number of comments/warnings
 call specialmsg_getcount(ncomment,nwarning,nexit)
 call libpaw_spmsg_getcount(ncomment_paw,nwarning_paw,nexit_paw)
 ncomment=ncomment+ncomment_paw;nwarning=nwarning+nwarning_paw;nexit=nexit+nexit_paw
 warn_fmt='(a,i6,a,i6,a)'
 if (nwarning<10000.and.ncomment<10000) warn_fmt='(a,i5,a,i5,a)'
 if (nwarning<1000 .and.ncomment<1000 ) warn_fmt='(a,i4,a,i4,a)'

#if defined HAVE_MPI
 write(std_out,'(a,i4,a,f13.1,a,f13.1)')' Proc.',mpi_enregs(1)%me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 if(me==0)then
   write(ab_out,'(3a,i4,a,f13.1,a,f13.1)')'-',ch10,'- Proc.',me,' individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
 end if
 call xmpi_sum(tsec, xmpi_world, ierr)
#else
 write(ab_out, '(a,a,a,f13.1,a,f13.1)' )'-',ch10,'- Proc.   0 individual time (sec): cpu=',tsec(1),'  wall=',tsec(2)
#endif

 write(msg,'(a,80a,a,a,a)' ) ch10,('=',mu=1,80),ch10,ch10,' Calculation completed.'
 call wrtout(ab_out, msg)
 write(msg,fmt=warn_fmt) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
 if (nexit/=0) write(msg,'(3a)') trim(msg),ch10,' Note : exit requested by the user.'
 call wrtout(ab_out, msg)

 if (me==0) then
   write(ab_out, '(a,f13.1,a,f13.1)' )'+Overall time at end (sec) : cpu=',tsec(1),'  wall=',tsec(2)
   write(msg, '(a,a)' ) ch10,' Calculation completed.'
   call wrtout(std_out, msg)
   write(msg,fmt=warn_fmt) '.Delivered',nwarning,' WARNINGs and',ncomment,' COMMENTs to log file.'
   if (nexit/=0) write(msg,'(3a)') trim(msg),ch10,' Note : exit requested by the user.'
   call wrtout(std_out, msg)
 end if

 if (me==0) then
   ! Write YAML document with the final summary.
   ! We use this doc to test whether the calculation is completed.
   write(std_out,"(a)")
   write(std_out,"(a)")"--- !FinalSummary"
   write(std_out,"(a)")"program: abinit"
   write(std_out,"(2a)")"version: ",trim(abinit_version)
   write(std_out,"(2a)")"start_datetime: ",start_datetime
   write(std_out,"(2a)")"end_datetime: ",asctime()
   write(std_out,"(a,f13.1)")"overall_cpu_time: ",tsec(1)
   write(std_out,"(a,f13.1)")"overall_wall_time: ",tsec(2)
   write(std_out,"(2a)")"exit_requested_by_user: ",yesno(nexit /= 0)
   write(std_out,"(2a)")"timelimit: ",trim(get_timelimit_string())
   write(std_out,"(a)")"pseudos: "
   do ii=1,npsp
     write(std_out,"(4a)")"    ",ljust(znucl2symbol(pspheads(ii)%znuclpsp), 4),": ",trim(pspheads(ii)%md5_checksum)
   end do
   write(std_out,"(a,i0)")"usepaw: ",dtsets(1)%usepaw
   write(std_out,"(a,i0)")"mpi_procs: ",xmpi_comm_size(xmpi_world)
   write(std_out,"(a,i0)")"omp_threads: ",xomp_get_num_threads(open_parallel=.True.)
   write(std_out,"(a,i0)")"num_warnings: ",nwarning
   write(std_out,"(a,i0)")"num_comments: ",ncomment
   write(std_out,"(a)")"..."
   call flush_unit(std_out)
 end if

 if (me==0) then
   if (xml_output) call outxml_finalise(tsec, values)
#ifndef HAVE_MEM_PROFILING
   close(unit=ab_out)
#endif
 end if

 ! 21) Eventual cleaning of MPI (and/or GPU) run
 call clnmpi_img(mpi_enregs(0))
 do ii=1,ndtset_alloc
   if(mpi_enregs(ii)%me<0) cycle
   call clnmpi_img(mpi_enregs(ii))
   call clnmpi_grid(mpi_enregs(ii))
   call clnmpi_atom(mpi_enregs(ii))
   call clnmpi_pert(mpi_enregs(ii))
 end do
 do ii=0,max(1,ndtset)
   call destroy_mpi_enreg(mpi_enregs(ii))
 end do
 ABI_FREE(mpi_enregs)

 ! If memory profiling is activated, check if bigdft plugin is used or not
 print_mem_report = 1
 do ii=1,ndtset_alloc
   if ((dtsets(ii)%usewvl == 1) .or. (dtsets(ii)%icoulomb > 0)) then
     print_mem_report = 0
     exit
   end if
 end do

 ! Here we deallocate dtsets. Do not access dtsets after this line!
 do ii=0,size(dtsets)-1,1
   call dtsets(ii)%free()
 end do
 ABI_FREE(dtsets)
 ABI_FREE(pspheads)

#if defined HAVE_BIGDFT
 call f_lib_finalize()
#endif

#if defined HAVE_GPU_CUDA
 call unsetdevice_cuda(use_gpu_cuda)
#endif

 call xpapi_shutdown()

 ! Writes information on file about the memory before ending mpi module, if memory profiling is enabled
 call abinit_doctor(filnam(4), print_mem_report=print_mem_report)

 call flush_unit(std_out)
 call flush_unit(ab_out)

 if (me == 0) close(unit=ab_out)

 100 call xmpi_end()

 end program abinit
!!***

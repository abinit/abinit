!{\src2tex{textfont=tt}}
!!****p* ABINIT/fftprof
!! NAME
!! fftprof
!!
!! FUNCTION
!!  Utility for profiling the FFT libraries supported by ABINIT.
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
!! OUTPUT
!!  Timing analysis of the different FFT libraries and algorithms.
!!
!! NOTES
!!  Point-group symmetries are not taken into account in getng during the generation
!!  of the FFT mesh. Therefore the FFT mesh might differ from the one
!!  found by abinit for the same cutoff and Bravais lattice (actually it might be smaller).
!!
!! Description of the input file (Fortran NAMELIST):
!!
!!   &CONTROL
!!     tasks = string specifying the tasks to perform i.e. the routines that should be tested or profiled.
!!             allowed values:
!!                 fourdp --> Test FFT transforms of density and potentials on the full box.
!!                 fourwf --> Test FFT transforms of wavefunctions using the zero-padded algorithm.
!!                 gw_fft --> Test the FFT transforms used in the GW code.
!!                 all    --> Test all FFT routines (DEFAULT)
!!     fftalgs = list of fftalg values (used to select the FFT libraries to use, see abinit doc for more info)
!!     ncalls = integer defining the number of calls for each tests. The final Wall time and CPU time
!!              are computed by averaging the final results over ncalls executions.
!!     max_nthreads = Maximum number of threads (DEFAULT 1, meaningful only if the code
!!                    uses threaded external libraries or OpenMP parallelization)
!!     ndat   = integer specifying how many FFT transforms should be executed for each call to the FFT routine
!!              (same meaning as the ndat input variable passed to fourwf)
!!     necut  = Used if tasks = "bench". Specifies the number of cutoff energies to profile (see also ecut_arth)
!!     ecut_arth = Used if tasks = "bench". Used to generate an arithmetic progression of cutoff energies
!!                 whose starting value is ecut_arth(1) and whose step is ecut_arth(2)
!!
!!  &SYSTEM
!!     ecut = cutoff energy for wavefunctions (real, Hartree units)
!!     rprimd = Direct lattice vectors in Bohr. (3,3) matrix in Fortran column-major order
!!     kpoint = real(3) vector specifying the reduced coordinates of the k-point of the wavefunction (used if tasks = "fourwf").
!!               The value of the k-point defines the storage scheme (istwfk) of the u(G) coefficients and therefore
!!               the FFT algorithm used to perform the transform u(G) <--> u(R) in fourwf.
!!     osc_ecut  = cutoff energy (Hartree) for the oscillator matrix elements computed in the GW code
!!                 Corresponds to the input variables (ecuteps, ecutsigx) used in the main code.
!!     nsym     =Number of symmetry operations (DEFAULT 1)
!!     symrel(3,3,nsym) = Symmetry operation in real space used to select the FFT mesh in the routine getng (default: Identity matrix)
!!
!! PARENTS
!!
!! CHILDREN
!!      abi_io_redirect,abimem_init,abinit_doctor,destroy_mpi_enreg
!!      fft_test_free,fft_test_init,fft_test_nullify,fft_test_print
!!      fft_use_lib_threads,fftprof_free,fftprof_ncalls_per_test,fftprof_print
!!      fftw3_init_threads,flush_unit,get_kg,getng,herald,initmpi_seq,lower
!!      metric,prof_fourdp,prof_fourwf,prof_rhotwg,time_fftbox,time_fftu
!!      time_fourdp,time_fourwf,time_rhotwg,wrtout,xmpi_bcast,xmpi_init
!!      xomp_show_info
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

program fftprof

 use defs_basis
 use defs_abitypes
 use m_build_info
 use m_xmpi
 use m_xomp
 use m_errors
 use m_FFT_prof
 use m_abicore
 use m_dfti

 use m_fstrings,   only : lower
 use m_specialmsg, only : specialmsg_getcount, herald
 use m_io_tools,   only : flush_unit
 use m_geometry,   only : metric
 use m_fftcore,    only : get_cache_kb, get_kg, fftalg_isavailable, fftalg_has_mpi, getng
 use m_fft,        only : fft_use_lib_threads, fftbox_utests, fftu_utests, fftbox_mpi_utests, fftu_mpi_utests
 use m_fftw3,      only : fftw3_init_threads
 use m_mpinfo,     only : destroy_mpi_enreg, initmpi_seq

 implicit none

!Arguments -----------------------------------
!Local variables-------------------------------
!scalars
 integer,parameter :: MAX_NFFTALGS=50,MAX_NSYM=48
 integer,parameter :: paral_kgb0=0,me_fft0=0,nproc_fft1=1,master=0
 integer :: ii,fftcache,it,cplex,ntests,option_fourwf,osc_npw
 integer :: map2sphere,use_padfft,isign,nthreads,comm,nprocs,my_rank
 integer :: iset,iall,inplace,nsets,avail,ith,idx,ut_nfft,ut_mgfft
 integer :: nfftalgs,alg,fftalg,fftalga,fftalgc,nfailed,ierr,paral_kgb
 character(len=24) :: codename
 character(len=500) :: header,msg
 real(dp),parameter :: boxcutmin2=two
 real(dp) :: ucvol
 logical :: test_fourdp,test_fourwf,test_gw,do_seq_bench,do_seq_utests
 logical :: iam_master,do_mpi_utests !,do_mpi_bench,
 type(MPI_type) :: MPI_enreg
!arrays
 integer :: ut_ngfft(18)
 real(dp),parameter :: gamma_point(3)=(/zero,zero,zero/)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3)
 type(FFT_test_t),allocatable :: Ftest(:)
 type(FFT_prof_t),allocatable :: Ftprof(:)
 integer,allocatable :: osc_gvec(:,:)
 integer,allocatable :: fft_setups(:,:),fourwf_params(:,:)
! ==========  INPUT FILE ==============
 integer :: ncalls=10,max_nthreads=1,ndat=1,necut=0,nsym=1
 character(len=500) :: tasks="all"
 integer :: fftalgs(MAX_NFFTALGS) = 0
 integer :: symrel(3,3,MAX_NSYM) = 0
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: ecut=30,osc_ecut=3
 real(dp) :: ecut_arth(2)=zero
 real(dp) :: rprimd(3,3)
 real(dp) :: kpoint(3) = (/0.1,0.2,0.3/)
 logical :: use_lib_threads = .FALSE.
 namelist /CONTROL/ tasks, ncalls, max_nthreads, ndat, fftalgs, necut, ecut_arth, use_lib_threads
 namelist /SYSTEM/ ecut, rprimd, kpoint, osc_ecut, nsym,symrel

! *************************************************************************

 ! Change communicator for I/O (mandatory!)
 call abi_io_redirect(new_io_comm=xmpi_world)

 call xmpi_init()
 comm = xmpi_world
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 iam_master = (my_rank == master)

 ! Initialize memory profiling if it is activated
 ! if a full abimem.mocc report is desired, set the argument of abimem_init to "2" instead of "0"
 ! note that abimem.mocc files can easily be multiple GB in size so don't use this option normally
#ifdef HAVE_MEM_PROFILING
 call abimem_init(0)
#endif

 if (iam_master) then
   ! Print header and read input file.
   codename='FFTPROF'//REPEAT(' ',17)
   call herald(codename,abinit_version,std_out)

   write(std_out,'(a)')" Tool for profiling and testing the FFT libraries used in ABINIT."
   write(std_out,'(a)')" Allowed options are: "
   write(std_out,'(a)')"   fourdp --> Test FFT transforms of density and potentials on the full box."
   write(std_out,'(a)')"   fourwf --> Test FFT transforms of wavefunctions using the zero-pad algorithm."
   write(std_out,'(a)')"   gw_fft --> Test the FFT transforms used in the GW code."
   write(std_out,'(a)')"   all    --> Test all FFT routines."
   !write(std_out,'(a)')"   seq-utests --> Units tests for sequential FFTs."
   !write(std_out,'(a)')"   mpi-utests --> Units tests for MPI FFTs."
   !write(std_out,'(a)')"   seq-bench  --> Benchmark mode for sequential FFTs."
   !write(std_out,'(a)')"   mpi-bench  --> Benchmark mode for MPI FFTs."
   write(std_out,'(a)')" "

   !Read input file
   read(std_in, NML=CONTROL)
   read(std_in, NML=SYSTEM)
   !write(std_out, NML=CONTROL)
   !write(std_out, NML=SYSTEM)
   call lower(tasks)
 end if

 ! Broadcast variables.
 if (nprocs > 1) then
   ! CONTROL
   call xmpi_bcast(tasks,master,comm,ierr)
   call xmpi_bcast(ncalls,master,comm,ierr)
   call xmpi_bcast(max_nthreads,master,comm,ierr)
   call xmpi_bcast(ndat,master,comm,ierr)
   call xmpi_bcast(fftalgs,master,comm,ierr)
   call xmpi_bcast(necut,master,comm,ierr)
   ! SYSTEM
   call xmpi_bcast(ecut_arth,master,comm,ierr)
   call xmpi_bcast(use_lib_threads,master,comm,ierr)
   call xmpi_bcast(ecut,master,comm,ierr)
   call xmpi_bcast(rprimd,master,comm,ierr)
   call xmpi_bcast(kpoint,master,comm,ierr)
   call xmpi_bcast(osc_ecut,master,comm,ierr)
   call xmpi_bcast(nsym,master,comm,ierr)
   call xmpi_bcast(symrel,master,comm,ierr)
 end if

 ! replace "+" with white spaces
 do ii=1,LEN_TRIM(tasks)
   if (tasks(ii:ii) == "+") tasks(ii:ii) = " "
 end do
 !%call str_replace_chars(tasks,"+"," ")
 tasks = " "//TRIM(tasks)//""

 iall=INDEX (tasks,"all")
 test_fourdp = (iall>0 .or. INDEX(tasks," fourdp")>0 )
 test_fourwf = (iall>0 .or. INDEX(tasks," fourwf")>0 )
 test_gw     = (iall>0 .or. INDEX(tasks," gw_fft")>0 )

 do_seq_utests = INDEX(tasks," utests")>0
 do_seq_bench  = INDEX(tasks," bench")>0

 do_mpi_utests = INDEX(tasks," mpi-utests")>0
 !do_mpi_bench  = INDEX(tasks," mpi-bench")>0

#ifdef HAVE_FFT_FFTW3_THREADS
 call fftw3_init_threads()
#endif

#ifndef HAVE_OPENMP
 ABI_CHECK(max_nthreads<=1,"nthreads>1 but OMP support is not enabled!")
#endif
 if (max_nthreads>1.and. iam_master) then
   call xomp_show_info(std_out)
 end if

 call fft_use_lib_threads(use_lib_threads)
 !write(std_out,*)"use_lib_threads: ",use_lib_threads

 if (do_mpi_utests) then
   ! Execute unit tests for MPI FFTs and terminate execution.
   call wrtout(std_out,"=== MPI FFT Unit Tests ===","COLL")

   ! List the FFT libraries that will be tested.
   ! Goedecker FFTs are always available, other libs are optional.
   nfftalgs = COUNT(fftalgs/=0)
   ABI_CHECK(nfftalgs > 0,"fftalgs must be specified")

   nfailed = 0; nthreads = 0
   do ii=1,nfftalgs
     fftalg = fftalgs(ii)
     do paral_kgb=1,1
       write(msg,"(5(a,i0))")&
&       "MPI fftu_utests with fftalg = ",fftalg,", paral_kgb = ",paral_kgb," ndat = ",ndat,", nthreads = ",nthreads
       call wrtout(std_out,msg,"COLL")
       nfailed = nfailed + fftu_mpi_utests(fftalg,ecut,rprimd,ndat,nthreads,comm,paral_kgb)
     end do
   end do
   ! TEMPORARY
   !goto 1

   nfailed = 0; nthreads = 0
   do ii=1,nfftalgs
     fftalg = fftalgs(ii)
     do cplex=1,2
       write(msg,"(4(a,i0))")"MPI fftbox_utests with fftalg = ",fftalg,", cplex = ",cplex," ndat = ",ndat,", nthreads = ",nthreads
       call wrtout(std_out,msg,"COLL")
       nfailed = nfailed + fftbox_mpi_utests(fftalg=fftalg,cplex=cplex,ndat=ndat,nthreads=nthreads,comm_fft=comm)
     end do
   end do

   write(msg,'(a,i0)')"Total number of failed tests = ",nfailed
   call wrtout(std_out,msg,"COLL")
   goto 100 ! Jump to xmpi_end
 end if

 call initmpi_seq(MPI_enreg)

 call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)

 ! symmetries (if given) are used for defining the FFT mesh.
 if (nsym==1 .and. ALL(symrel==0)) then
   symrel(:,:,1) = RESHAPE((/1,0,0,0,1,0,0,0,1/),(/3,3/))
 end if

 ! Set the number of calls for test.
 ncalls=ABS(ncalls); if (ncalls==0) ncalls=10
 call fftprof_ncalls_per_test(ncalls)

 ! List the FFT libraries that will be tested.
 ! Goedecker FFTs are always available, other libs are optional.
 ! write(std_out,*)"ndat = ",ndat
 nfftalgs = COUNT(fftalgs/=0)
 ABI_CHECK(nfftalgs > 0,"fftalgs must be specified")

 ntests = max_nthreads * nfftalgs

 ! First dimension contains (fftalg,fftcache,ndat,nthreads,available).
 ABI_MALLOC(fft_setups,(5,ntests))

 ! Default Goedecker library.
 idx=0
 do alg=1,nfftalgs
   fftalg = fftalgs(alg)
   fftalga=fftalg/100; fftalgc=mod(fftalg,10)
   avail = 1
   if (.not.fftalg_isavailable(fftalg)) avail = 0
   do ith=1,max_nthreads
     !fftcache is machine-dependent.
     fftcache = get_cache_kb()
     idx = idx + 1
     fft_setups(:,idx) = (/fftalg,fftcache,ndat,ith,avail/)
   end do
 end do

 ! Compute FFT box.
 ABI_DT_MALLOC(Ftest,(ntests))
 ABI_DT_MALLOC(Ftprof,(ntests))

 do it=1,ntests ! Needed for crappy compilers that do not support => null() in type declarations.
   call fft_test_nullify(Ftest(it))
 end do

 do it=1,ntests
   call fft_test_init(Ftest(it),fft_setups(:,it),kpoint,ecut,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)
    ! Ftest%results is allocated using nfftot and the consistency btw libs is tested assuming an equal number of FFT-points.
   if ( ANY(Ftest(it)%ngfft(1:3) /= Ftest(1)%ngfft(1:3)) ) then
     MSG_ERROR("Consistency check assumes equal FFT meshes. Cannot continue")
   end if
   if (it == 1) then
     call fft_test_print(Ftest(it)) !,header)
   else if (fft_setups(1,it) /= fft_setups(1,it-1)) then
     call fft_test_print(Ftest(it)) !,header)
   end if
 end do
 !
 ! =======================
 ! ==== fourdp timing ====
 ! =======================
 if (test_fourdp) then
   do isign=-1,1,2
     do cplex=1,2
       do it=1,ntests
         call time_fourdp(Ftest(it),isign,cplex,header,Ftprof(it))
       end do
       call fftprof_print(Ftprof,header)
       call fftprof_free(Ftprof)
     end do
   end do
 end if
 !
 ! =======================
 ! ==== fourwf timing ====
 ! =======================
 if (test_fourwf) then
   ! possible combinations of (option, cplex) supported in fourwf.
   ! (cplex=2 only allowed for option=2, and istwf_k=1)
   nsets=4; if (Ftest(1)%istwf_k==1) nsets=5
   ABI_MALLOC(fourwf_params,(2,nsets))
   fourwf_params(:,1) = (/0,0/)
   fourwf_params(:,2) = (/1,1/)
   fourwf_params(:,3) = (/2,1/)
   fourwf_params(:,4) = (/3,0/)
   if (nsets==5) fourwf_params(:,5) = (/2,2/)

   do iset=1,nsets
     option_fourwf = fourwf_params(1,iset)
     cplex         = fourwf_params(2,iset)
     do it=1,ntests
       call time_fourwf(Ftest(it),cplex,option_fourwf,header,Ftprof(it))
     end do
     call fftprof_print(Ftprof,header)
     call fftprof_free(Ftprof)
   end do
   ABI_FREE(fourwf_params)
 end if
 !
 ! ==========================
 ! ==== Test GW routines ====
 ! ==========================
 ! These routines are used in the GW part, FFTW3 is expected to
 ! be more efficient as the conversion complex(:) <--> real(2,:) is not needed.
 if (test_gw) then
   !
   ! ==== fourdp timing with complex arrays ====
   do isign=-1,1,2
     do inplace=0,1
       do it=1,ntests
         call time_fftbox(Ftest(it),isign,inplace,header,Ftprof(it))
       end do
       call fftprof_print(Ftprof,header)
       call fftprof_free(Ftprof)
     end do
   end do
   !
   ! ==== zero padded with complex arrays ====
   do isign=-1,1,2
     do it=1,ntests
       call time_fftu(Ftest(it),isign,header,Ftprof(it))
     end do
     call fftprof_print(Ftprof,header)
     call fftprof_free(Ftprof)
   end do
   !
   ! ==== rho_tw_g timing ====
   ABI_CHECK(osc_ecut > zero,"osc_ecut <= zero!")

   call get_kg(gamma_point,1,osc_ecut,gmet,osc_npw,osc_gvec)
   ! TODO should reorder by shells to be consistent with the GW part!
   ! Moreover I guess this ordering is more efficient when we have
   ! to map the box to the G-sphere!
   map2sphere=1; !map2sphere=0

   do use_padfft=0,1
     do it=1,ntests
       call time_rhotwg(Ftest(it),map2sphere,use_padfft,osc_npw,osc_gvec,header,Ftprof(it))
     end do
     call fftprof_print(Ftprof,header)
     call fftprof_free(Ftprof)
   end do

   ABI_FREE(osc_gvec)
 end if ! test_gw

 if (do_seq_utests) then
   call wrtout(std_out,"=== FFT Unit Tests ===","COLL")

   nfailed = 0
   do idx=1,ntests
     ! fft_setups(:,idx) = (/fftalg,fftcache,ndat,ith,avail/)
     fftalg   = fft_setups(1,idx)
     fftcache = fft_setups(2,idx)
     ndat     = fft_setups(3,idx)
     nthreads = fft_setups(4,idx)
     ! Skip the test if library is not available.
     if (fft_setups(5,idx) == 0) CYCLE

     write(msg,"(3(a,i0))")"fftbox_utests with fftalg = ",fftalg,", ndat = ",ndat,", nthreads = ",nthreads
     call wrtout(std_out,msg,"COLL")

     nfailed = nfailed + fftbox_utests(fftalg,ndat,nthreads)

     ! Warning: This routine is not thread safe hence we have to initialize ngfft it here.
     ut_ngfft = -1
     ut_ngfft(7) = fftalg
     ut_ngfft(8) = fftcache

     call getng(boxcutmin2,ecut,gmet,k0,me_fft0,ut_mgfft,ut_nfft,ut_ngfft,nproc_fft1,nsym,&
&     paral_kgb0,symrel,unit=dev_null)

     write(msg,"(3(a,i0))")"fftu_utests with fftalg = ",fftalg,", ndat = ",ndat,", nthreads = ",nthreads
     call wrtout(std_out,msg,"COLL")

     nfailed = nfailed + fftu_utests(ecut,ut_ngfft,rprimd,ndat,nthreads)
   end do

   write(msg,'(a,i0)')"Total number of failed tests = ",nfailed
   call wrtout(std_out,msg,"COLL")
 end if

! Benchmarks for the sequential version.
 if (do_seq_bench) then
   call wrtout(std_out,"Entering benchmark mode","COLL")
   write(std_out,*)"ecut_arth",ecut_arth,", necut ",necut

   if (INDEX(tasks,"bench_fourdp")>0) then
     isign=1; cplex=1
     call prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)

     isign=1; cplex=2
     call prof_fourdp(fft_setups,isign,cplex,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)
   end if

   if (INDEX(tasks,"bench_fourwf")>0) then
     cplex=2; option_fourwf=0
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)

     cplex=1; option_fourwf=1
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)

     cplex=1; option_fourwf=2
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)

     cplex=2; option_fourwf=2
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)

     cplex=0; option_fourwf=3
     call prof_fourwf(fft_setups,cplex,option_fourwf,kpoint,necut,ecut_arth,boxcutmin2,rprimd,nsym,symrel,MPI_enreg)
   end if

   if (INDEX(tasks,"bench_rhotwg")>0) then
     map2sphere=1; use_padfft=1
     call prof_rhotwg(fft_setups,map2sphere,use_padfft,necut,ecut_arth,osc_ecut,boxcutmin2,&
&     rprimd,nsym,symrel,gmet,MPI_enreg)

     map2sphere=1; use_padfft=0
     call prof_rhotwg(fft_setups,map2sphere,use_padfft,necut,ecut_arth,osc_ecut,boxcutmin2,&
&     rprimd,nsym,symrel,gmet,MPI_enreg)
   end if
 end if
 !
 !===============================
 !=== End of run, free memory ===
 !===============================
 call wrtout(std_out,ch10//" Analysis completed.","COLL")

 ABI_FREE(fft_setups)

 call fft_test_free(Ftest)
 ABI_DT_FREE(Ftest)
 call fftprof_free(Ftprof)
 ABI_DT_FREE(Ftprof)
 call destroy_mpi_enreg(MPI_enreg)

 call flush_unit(std_out)

 call abinit_doctor("__fftprof")

 100 call xmpi_end()

 end program fftprof
!!***

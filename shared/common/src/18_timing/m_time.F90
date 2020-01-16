!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_time
!! NAME
!! m_time
!!
!! FUNCTION
!! This module contains accumulators for the timer.
!! and functions to get cpu and wall time.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2019 ABINIT group (MG, XG, MT, TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_time

 use defs_basis
 use m_abicore
 use m_errors
 use iso_c_binding
 use m_xmpi
#if defined HAVE_MPI2
 use mpi
#endif
 use m_clib

 use m_xpapi,    only: xpapi_flops
 use m_fstrings, only: char_count, sjoin

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private

 public :: asctime       ! Build a 24-character string of the following form: 'Sun Jun 20 23:21:05 1993'.
 public :: sec2str       ! Convert time data in seconds to string
 public :: str2sec       ! Convert a string with time (Slurm form) in seconds
 public :: abi_wtime     ! Returns wall clock time in seconds since some arbitrary start.
 public :: abi_cpu_time  ! Returns cpu time in seconds since some arbitrary start.
 public :: cwtime        ! Returns cpu, wall clock time and gflops
 public :: cwtime_report ! Stop timers, write message, reinit counters.

 ! FIXME: Deprecated Should be replaced by cwtime
 public :: timein
 public :: time_accu
 public :: timab
 public :: time_set_papiopt
 public :: time_get_papiopt

!!***

 ! papiopt is a flag which indicates if there is or not an analysis of speed execution is made.
 ! By defaut the analysis is not done
 integer,private,save :: papiopt=0

 !==================
 ! Counter variables
 !==================

 ! TIMER_SIZE determines the maximum number of "timing slots" available
 integer,public,parameter :: TIMER_SIZE=1999

 ! timeopt is a flag which indicates the suppression or not of the timing.
 integer,private,save :: timopt=1

 ! Number of times that the routine has been called
 integer,private,save :: ncount(TIMER_SIZE)=0

 ! Accumulating cpu time (1) and wall to wall time (2) for each "timing slots"
 real(dp),private,save  :: acctim(2,TIMER_SIZE)=zero,tzero(2,TIMER_SIZE)=zero

 ! Accumulating number of floating point operation and cpu time (1) and wall time (2) for each "performance slot"
 real(dp),private,save :: papi_accflops(TIMER_SIZE)=zero, papi_acctim(2,TIMER_SIZE)=zero

 ! Reference value for number of floating point operation and time (cpu and wall) for each performance slot
 real(dp),private,save :: papi_flops(TIMER_SIZE)=zero , papi_tzero(2,TIMER_SIZE)=zero

 ! Elapsed time and elapsed number of floating point operation since a reference
#ifdef HAVE_PAPI
 real(dp),private,save :: papi_tottim(2,TIMER_SIZE)=zero, papi_totflops(TIMER_SIZE)=zero
#endif

CONTAINS
!!***

!!****f* m_time/asctime
!! NAME
!!  asctime
!!
!! FUNCTION
!!   Build a 24-character string of the following form: 'Sun Jun 20 23:21:05 1993'.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function asctime()

!Arguments ------------------------------------
 character(len=24) :: asctime

!Local variables-------------------------------
 integer :: day,dd,ja,jy,jm,jdn,mm,year
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=3),parameter :: day_names(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!Get year, month and day
 call date_and_time(strdat,strtime,strzone,values)

 year=values(1)
 mm=values(2)
 dd=values(3)

!Get day of the week
 if (mm > 2) then
   jy=year
   jm=mm+1
 else
   jy=year-1
   jm=mm+13
 end if

 jdn=int(365.25d0*jy)+int(30.6001d0*jm)+dd+1720995
 ja=int(0.01d0*jy)
 jdn=jdn+2-ja+int(quarter*ja)
 day=mod(jdn,7)+1

 ! Build a 24-character string of the following form: 'Sun Jun 20 23:21:05 1993'.
 write(asctime, '(a,1x,a,1x,i0.2,1x,2(i0.2,a),i0.2,1x,i4)')&
   day_names(day),month_names(mm),dd,values(5),":",values(6),":",values(7),year

end function asctime
!!***

!----------------------------------------------------------------------

!!****f* m_time/sec2str
!! NAME
!!  sec2str
!!
!! FUNCTION
!!  Convert time data in seconds to string
!!
!! INPUTS
!!   time_s=Time in seconds
!!
!! OUTPUT
!!   string with time displayed in the form: [days-][hours:][minutes:]seconds
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

pure function sec2str(time_s) result(str)

!Arguments ------------------------------------
!scalars
 real(dp),intent(in) :: time_s
 character(len=500) :: str

!Local variables-------------------------------
 integer :: days,hours,minutes,seconds

! *************************************************************************

 days    = time_s / 86400
 hours   = MOD(time_s,86400._dp) / 3600
 minutes = MOD(time_s,3600._dp) / 60
 seconds = MOD(time_s,60._dp)

 if (days > 0) then
   write(str,'(i0,3(a,i0.2),a)')days,"-",hours,":",minutes,":",seconds, " [days]"
 else if (hours > 0) then
   write(str,'(i0.2,2(a,i0.2),a)')hours,":",minutes,":",seconds, " [hours]"
 else if (minutes > 0) then
   write(str,'(i0.2,a,i0.2,a)')minutes,":",seconds, " [minutes]"
 else
   write(str,'(f5.2,a)')time_s," [s]"
 end if

end function sec2str
!!***

!----------------------------------------------------------------------

!!****f* m_time/str2sec
!! NAME
!!  str2sec
!!
!! FUNCTION
!!  Convert a string to time data in seconds. Return negative value if not valid string
!!  Accepts a string in one the following (SLURM) forms:
!!
!!     # "days-hours",
!!     # "days-hours:minutes",
!!     # "days-hours:minutes:seconds".
!!     # "minutes",
!!     # "minutes:seconds",
!!     # "hours:minutes:seconds",
!!
!! PARENTS
!!
!! SOURCE

real(dp) pure function str2sec(str) result(time)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: str

!Local variables-------------------------------
 integer :: days,hours,minutes,seconds,dash,i,j

! *************************************************************************

 days = 0; hours = 0; minutes = 0; seconds = 0
 dash = index(str, "-")
 if (dash /= 0) read(str(:dash-1),*,err=1) days

 select case (char_count(str, ":"))
 case (0)
   if (dash /= 0) then
     read(str(dash+1:),*,err=1)hours
   else
     read(str(dash+1:),*,err=1)minutes
   end if

 case (1)
   i = index(str, ":")
   if (dash /= 0) then
     read(str(dash+1:i-1),*,err=1)hours
     read(str(i+1:),*,err=1)minutes
   else
     read(str(:i-1),*,err=1)minutes
     read(str(i+1:),*,err=1)seconds
   end if

 case(2)
   i = index(str, ":")
   read(str(dash+1:i-1),*,err=1)hours
   j = index(str(i+1:), ":") + i
   read(str(i+1:j-1),*,err=1)minutes
   read(str(j+1:),*,err=1)seconds

 case default
   time = -one; return
 end select

 time = 24 * 3600 * days + hours * 3600 + minutes * 60 + seconds
 return

1 time = -one

end function str2sec
!!***

!----------------------------------------------------------------------

!!****f* m_time/abi_cpu_time
!! NAME
!!  abi_cpu_time
!!
!! FUNCTION
!!  Timing routine. Returns cpu time in seconds since some arbitrary start.
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  cpu_time= cpu time in seconds
!!
!! NOTES
!!  For CPU time, contains machine-dependent code (choice will be selected by c preprocessor).
!!  Note that all supported machines are listed explicitly below; there
!!  is no "else" which covers "other".  The C preprocessor will place
!!  a spurious line of code (see below) into the fortran source unless
!!  preprocessed with -Dflag where flag refers to one of the supported machines.
!!
!!  WARNING: the following list is no more accurate (YP 20060530)
!!
!!  Presently supported flags: "ibm", "hp", "P6", "dec_alpha", "sgi", "vpp", "sun", "mac", "nec", "sr8k"
!!  Previously supported flags:  "ultrix". Might still work !
!!
!!  Calls machine-dependent "mclock" for "ibm" .
!!  Calls ANSI C subroutine "cclock" for "hp" and "sgi".
!!  Calls machine-dependent "etime" for "P6", "mac", "dec_alpha", "sun", "nec" .
!!  Calls machine-dependent "clock" for "vpp"
!!  Calls machine-dependent "xclock" for "sr8k"
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function abi_cpu_time() result(cpu)

!Arguments ------------------------------------
 real(dp) :: cpu

!Local variables-------------------------------
#ifdef HAVE_FC_CPUTIME
 real :: cpu_sp
#elif defined FC_IBM
 integer :: mclock
#elif defined HAVE_OS_MACOSX
 real :: tmp(2) !real array only needed by etime
 real(dp) :: etime
#else
 integer :: count_now,count_max,count_rate
#endif

! *************************************************************************

!Machine-dependent timers
#ifdef HAVE_CCLOCK
 call clib_cclock(cpu)

#elif defined HAVE_FC_CPUTIME
!This is the F95 standard subroutine.
 call cpu_time(cpu_sp)
 cpu = cpu_sp

#elif defined FC_IBM
 cpu = mclock()*0.01d0

#elif defined HAVE_OS_MACOSX
 cpu = clib_etime(tmp)

#else
!This is the Fortran90 standard subroutine, might not always be sufficiently accurate
 call system_clock(count_now,count_rate,count_max)
 cpu=dble(count_now)/dble(count_rate)
#endif

end function abi_cpu_time
!!***

!----------------------------------------------------------------------

!!****f* m_time/abi_wtime
!! NAME
!!  abi_wtime
!!
!! FUNCTION
!!  Return wall clock time in seconds since some arbitrary start.
!!  Call the F90 intrinsic date_and_time .
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  wall= wall clock time in seconds
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function abi_wtime() result(wall)

!Arguments ------------------------------------
!scalars
 real(dp) :: wall

!Local variables-------------------------------
!scalars
#ifndef HAVE_MPI
 integer,parameter :: nday(24)=(/31,28,31,30,31,30,31,31,30,31,30,31,&
&                                31,28,31,30,31,30,31,31,30,31,30,31/)
 integer,save :: month_init,month_now,start=1,year_init
 integer :: months
 character(len=8)   :: date
 character(len=10)  :: time
 character(len=5)   :: zone
 character(len=500) :: msg
!arrays
 integer :: values(8)
#endif

! *************************************************************************

#ifndef HAVE_MPI

!The following section of code is standard F90, but it is useful only if the intrinsics
!date_and_time is accurate at the 0.01 sec level, which is not the case for a P6 with the pghpf compiler ...
!Year and month initialisation
 if(start==1)then
   start=0
   call date_and_time(date,time,zone,values)
   year_init=values(1)
   month_init=values(2)
 end if

!Uses intrinsic F90 subroutine Date_and_time for
!wall clock (not correct when a change of year happen)
 call date_and_time(date,time,zone,values)

!Compute first the number of seconds from the beginning of the month
 wall=(values(3)*24.0d0+values(5))*3600.0d0+values(6)*60.0d0+values(7)+values(8)*0.001d0

!If the month has changed, compute the number of seconds
!to be added. This fails if the program ran one year !!
 month_now=values(2)
 if(month_now/=month_init)then
   if(year_init+1==values(1))then
     month_now=month_now+12
   end if
   if(month_now<=month_init)then
     msg = 'Problem with month and year numbers.'
     MSG_BUG(msg)
   end if
   do months=month_init,month_now-1
     wall=wall+86400.0d0*nday(months)
   end do
 end if

!Now take into account bissextile years (I think 2000 is bissextile, but I am not sure ...)
 if(mod(year_init,4)==0 .and. month_init<=2 .and. month_now>2)   wall=wall+3600.0d0
 if(mod(values(1),4)==0 .and. month_init<=14 .and. month_now>14) wall=wall+3600.0d0

#else
!Use the timer provided by MPI1.
 wall = MPI_WTIME()
#endif

end function abi_wtime
!!***

!----------------------------------------------------------------------

!!****f* m_time/cwtime
!! NAME
!!  cwtime
!!
!! FUNCTION
!!  Timing routine. Returns cpu and wall clock time in seconds.
!!
!! INPUTS
!!  start_or_stop=
!!    "start" to start the timers
!!    "stop" to stop the timers and return the final cpu_time and wall_time
!!  [msg]: Optional message printed to std_out
!!  [comm]: MPI communicator. If values averaged inside comm are wanted. Only for "stop"
!!
!! OUTPUT
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!  gflops = Gigaflops
!!
!! NOTES
!!  Example:
!!  ! Init cpu and wall
!!  call cwtime(cpu,wall,gflops,"start")
!!
!!  do_stuff()
!!
!!  ! stop the counters, return cpu- and wall-time spent in do_stuff()
!!  call cwtime(cpu,wall,gflops,"stop")
!!
!! PARENTS
!!      calc_sigc_me,calc_sigx_me,cchi0,cchi0q0,eph,exc_build_block
!!      exc_build_ham,lapackprof,m_abilasi,m_bse_io,m_dvdb,m_epjdos
!!      m_exc_itdiago,m_fft,m_fft_prof,m_fstab,m_gkk,m_gruneisen,m_ifc,m_ioarr
!!      m_iowf,m_phgamma,m_phonons,m_phpi,m_shirley,m_sigmaph,m_skw,m_wfd,m_wfk
!!      partial_dos_fractions
!!
!! CHILDREN
!!      xpapi_flops
!!
!! SOURCE

subroutine cwtime(cpu, wall, gflops, start_or_stop, msg, comm)

!Arguments ------------------------------------
!scalars
 real(dp),intent(inout) :: cpu,wall
 real(dp),intent(out) :: gflops
 character(len=*),intent(in) :: start_or_stop
 character(len=*),intent(in),optional :: msg
 integer,intent(in),optional :: comm

!Local variables-------------------------------
#ifndef HAVE_PAPI
 logical,parameter :: use_papi=.FALSE.
#else
 logical,parameter :: use_papi=.TRUE.
#endif
 integer :: ierr
 integer(C_INT)  :: check
 integer(C_LONG_LONG) :: flops
 real(C_FLOAT) :: real_time,proc_time,mflops
 real(dp) :: vals(3)

! *************************************************************************

 if (present(msg)) call wrtout(std_out, msg)

 SELECT CASE (start_or_stop)
 CASE ("start")
 if (use_papi) then
   call xpapi_flops(real_time,proc_time,flops,mflops,check)
   cpu = proc_time; wall = real_time; gflops = mflops / 1000
 else
   cpu = abi_cpu_time(); wall = abi_wtime(); gflops = -one
 end if

 CASE ("stop")
 if (use_papi) then
   call xpapi_flops(real_time,proc_time,flops,mflops,check)
   cpu = proc_time - cpu; wall = real_time - wall; gflops = mflops / 1000
 else
   cpu = abi_cpu_time() - cpu; wall = abi_wtime() - wall; gflops = -one
 end if
 if (present(comm)) then
   vals = [cpu, wall, gflops]
   call xmpi_sum(vals, comm, ierr)
   vals = vals / xmpi_comm_size(comm)
   cpu = vals(1); wall = vals(2); gflops = vals(3)
 end if

 CASE DEFAULT
   MSG_ERROR("Wrong option for start_or_stop: "//TRIM(start_or_stop))
 END SELECT

end subroutine cwtime
!!***

!----------------------------------------------------------------------

!!****f* m_time/cwtime_report
!! NAME
!!  cwtime_report
!!
!! FUNCTION
!! Stop timers, write message, reinit counters.
!!
!! INPUT
!!  [pre_str], [end_str]: String to print before and after the timing section
!!  [comm]: MPI communicator. If values averaged inside comm is wanted. Only for "stop"
!!
!! SIDE EFFECTS
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!  gflops = Gigaflops
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cwtime_report(tag, cpu, wall, gflops, pre_str, end_str, comm)

!Arguments ------------------------------------
!scalars
 real(dp),intent(inout) :: cpu,wall
 real(dp),intent(out) :: gflops
 integer,intent(in),optional :: comm
 character(len=*),intent(in) :: tag
 character(len=*),optional,intent(in) :: pre_str, end_str

!Local variables-------------------------------
!scalars
 character(len=500) :: avg_type

! *************************************************************************

 if (present(comm)) then
   call cwtime(cpu, wall, gflops, "stop", comm=comm)
   avg_type = "(MPI average)"
 else
   call cwtime(cpu, wall, gflops, "stop")
   avg_type = ""
 end if
 if (present(pre_str)) call wrtout(std_out, pre_str)
 call wrtout(std_out, sjoin(tag, "completed. cpu:", sec2str(cpu), ", wall:", sec2str(wall), avg_type), &
     do_flush=.True.)
 if (present(end_str)) call wrtout(std_out, end_str)
 call cwtime(cpu, wall, gflops, "start")

end subroutine cwtime_report
!!***

!!****f* m_time/timein
!! NAME
!!  timein
!!
!! FUNCTION
!!  Timing routine. Returns cpu and wall clock time in seconds since some arbitrary start.
!!  For wall clock time, call the F90 intrinsic date_and_time.
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!
!! NOTES
!!  For CPU time, contains machine-dependent code (choice will be selected
!!  by C preprocessor, see abi_cpu_time).
!!
!! TODO
!!  Should be replaced by cwtime
!!
!! PARENTS
!!      abinit,aim,aim_follow,anaddb,bsepostproc,conducti,cpdrv,cut3d,drvaim
!!      elphon,first_rec,m_exit,mrgddb,mrgscr,multibinit,optic,rsurf,surf,timab
!!
!! CHILDREN
!!
!! SOURCE

subroutine timein(cpu,wall)

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: cpu,wall
! *************************************************************************

 ! CPU time
 cpu = abi_cpu_time()
 ! Wall time
 wall = abi_wtime()

end subroutine timein
!!***

!!****f* m_time/time_accu
!! NAME
!!  time_accu
!!
!! FUNCTION
!!  Return the number of times the counter has been called
!!  and corresponding data for given index
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed);
!!
!! OUTPUT
!!  tottim(2)=accumulated time for accumulator nn
!!  totftimes(2)=accumulated time for accumulator nn evaluated by papi
!!  totffops =accumulated number of flops for accumulator nn evaluated by papi
!!  return_ncount gives the number of times that the accumulator has been incremented
!!
!! PARENTS
!!      timana
!!
!! CHILDREN
!!
!! SOURCE

subroutine time_accu(nn,return_ncount,tottim,totflops,totftimes)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn
 integer,intent(out) :: return_ncount
 real(dp),intent(out) :: totflops
!arrays
 real(dp),intent(out) :: totftimes(2),tottim(2)

!Local variables-------------------------------
!scalars
 character(len=500) :: message

! *************************************************************************

!Check that nn lies in sensible bounds
 if (nn<0.or.nn>TIMER_SIZE) then
   write(message,'(a,i6,a,i8,a)')' dim TIMER_SIZE=',TIMER_SIZE,' but input nn=',nn,'.'
   MSG_BUG(message)
 end if

!return accumulated time for nn
 tottim(1)=acctim(1,nn)
 tottim(2)=acctim(2,nn)

!return accumulated number flops for nn
 totflops = papi_accflops(nn)

!return accumulated time for nn evaluated by papi
 totftimes(1) = papi_acctim(1,nn)
 totftimes(2) = papi_acctim(2,nn)
 return_ncount=ncount(nn)

end subroutine time_accu
!!***

!!****f* m_time/time_set_papiopt
!! NAME
!!  time_set_papiopt
!!
!! FUNCTION
!!  Set the value of papiopt
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine time_set_papiopt(opt)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: opt

! *************************************************************************

 papiopt = opt

end subroutine time_set_papiopt
!!***

!----------------------------------------------------------------------

!!****f* m_time/time_get_papiopt
!! NAME
!!  time_get_papiopt
!!
!! FUNCTION
!!  Return the value of papiopt
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function time_get_papiopt()

!Arguments ------------------------------------
!scalars
 integer :: time_get_papiopt

! *************************************************************************

 time_get_papiopt = papiopt

end function time_get_papiopt
!!***

!!****f* m_time/timab
!! NAME
!!  timab
!!
!! FUNCTION
!!  Timing subroutine. Calls machine-dependent "timein" which returns elapsed cpu and wall clock times in sec.
!!  Depending on value of "option" routine will:
!!
!!  (0) zero all accumulators
!!  (1) start with new incremental time slice for accumulator n using explicit call to timein (or PAPI)
!!  (2) stop time slice; add time to accumulator n also increase by one the counter for this accumulator
!!  (3) start with new incremental time slice for accumulator n
!!        using stored values for cpu, wall, and PAPI infos ( ! do not use for stop )
!!  (4) report time slice for accumlator n (not full time accumlated)
!!  (5) option to suppress timing (nn should be 0) or reenable it (nn /=0)
!!
!!  If, on first entry, subroutine is not being initialized, it
!!  will automatically initialize as well as rezero accumulator n.
!!  However, initialization SHOULD be done explicitly by the user
!!  so that it can be done near the top of his/her main routine.
!!
!! INPUTS
!!  nn=index of accumulator (distinguish what is being timed); NOT used if option=0
!!  option=see comment above
!!
!! OUTPUT
!!  on option=4 :
!!    tottim(2,nn)=accumulated time for accumulator nn; otherwise
!!     tottim is a dummy variable.
!!    option gives the number of times that the
!!     accumulator has been incremented
!!
!! PARENTS
!!      abinit,afterscfloop,atm2fft,bethe_salpeter,calc_sigc_me,calc_sigx_me
!!      calcdenmagsph,cchi0,cgq_builder,cgwf,chebfi,cohsex_me,corrmetalwf1,d2frnl
!!      density_rec,dfpt_cgwf,dfpt_dyfro,dfpt_dyxc1,dfpt_eltfrhar,dfpt_eltfrkin
!!      dfpt_eltfrloc,dfpt_eltfrxc,dfpt_ewald,dfpt_looppert,dfpt_mkrho
!!      dfpt_mkvxc,dfpt_mkvxc_noncoll,dfpt_mkvxcstr,dfpt_newvtr,dfpt_nstdy
!!      dfpt_nstpaw,dfpt_nstwf,dfpt_rhofermi,dfpt_rhotov,dfpt_scfcv,dfpt_vtorho
!!      dfpt_vtowfk,dfpt_wfkfermi,dfptnl_loop,dielmt,dieltcel,m_dmft
!!      dotprodm_v,dotprodm_vn,driver,dyson,eig2stern,eig2tot,elt_ewald
!!      eltxccore,energy,entropyrec,etotfor,exc_build_block,exc_build_ham
!!      fermisolverec,first_rec,fock2ACE,fock_getghc,forces,forstr,forstrnps
!!      fourdp,fourwf,fxphas,getgh1c,getghc,getgsc,getngrec,gran_potrec
!!      green_kernel,gstate,gstateimg,gwls_ComputeCorrelationEnergy
!!      gwls_DielectricArray,gwls_QR_factorization,gwls_lineqsolver
!!      gwls_model_polarisability,gwls_polarisability,gwls_sternheimer,hartre
!!      impurity_solve,initberry,initorbmag,initwf,inkpts,invars2,inwffil
!!      listkk,lobpcgwf,m_ab7_invars_f90,m_ab7_mixing,m_cgtools,m_dyson_solver
!!      m_fftcore,m_fftw3,m_fock,m_green,m_haydock,m_hexc,m_invovl,m_iowf
!!      m_lobpcg,m_lobpcg2,m_lobpcgwf,m_paral_pert,m_sg2002,m_wfutils,m_xg
!!      m_xgScalapack,mag_penalty,mkcore,mkcore_paw,mkcore_wvl,mkffnl
!!      mklocl_realspace,mklocl_recipspace,mkresi,mkrho,newkpt,newocc,newrho
!!      newvtr,nhatgrid,nlenergyrec,nonlinear,nonlop,odamix,opernla_ylm
!!      optics_paw,optics_paw_core,optics_vloc,outkss,outscfcv,pareigocc
!!      partial_dos_fractions_paw,pawdenpot,pawdfptenergy,pawinit,pawmknhat
!!      pawmknhat_psipsi,pawmkrho,pawpolev,prep_bandfft_tabs,prep_calc_ucrpa
!!      prep_fourwf,prep_getghc,prep_nonlop,pspatm,pspheads_comm,pspini
!!      pw_orthon,rayleigh_ritz,recursion,recursion_nl,respfn,rhotov,rhotoxc
!!      rwwf,scfcv,screening,setsym,setvtr,sigma,sqnormm_v,status,stress,strhar
!!      suscep_stat,susk,suskmm,symrhg,symsgcube,tddft,timana,vn_nl_rec,vtorho
!!      vtorhorec,vtorhotf,vtowfk,wf_mixing,wfconv,wfk_analyze,wfsinp
!!      wvl_nhatgrid,xcden,xcpot
!!
!! CHILDREN
!!      papif_flops,papif_perror,timein
!!
!! SOURCE
!!

subroutine timab(nn,option,tottim)

#ifdef HAVE_PAPI
#include "f90papi.h"
#endif

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nn,option
!arrays
 real(dp),intent(out) :: tottim(2)

!Local variables-------------------------------
!scalars
 real(dp),save :: cpu,wall
 character(len=500) :: message
#ifdef HAVE_PAPI
 integer(C_INT) :: check
 integer(C_LONG_LONG),save :: flops1
 real(C_FLOAT),save :: real_time,proc_time
 real(C_FLOAT) :: mflops1
 character(len=PAPI_MAX_STR_LEN) :: papi_errstr
#endif
! *************************************************************************

 if (option==5) timopt=nn

!If timopt was set to zero by a call with option=5, suppress
!all action of this routine (might as well return at this point !)
 if(timopt/=0 .and. option/=5)then
   ! Check that nn lies in sensible bounds
   if (nn<1.or.nn>TIMER_SIZE) then
     write(message,'(a,i0,a,i0)')'  TIMER_SIZE = ',TIMER_SIZE,' but input nn = ',nn
     MSG_BUG(message)
   end if

#ifdef HAVE_PAPI
   ! for all active options for time if papi analysis has been selected.
   if (option/=3.and.time_get_papiopt()==1) then
     call PAPIf_flops(real_time, proc_time, flops1, mflops1, check)
     if (check /= PAPI_OK) then
       call papif_perror(check,papi_errstr,check)
       write(std_out,*) 'Problem to initialize papi high level inteface'
       write(std_out,*) 'Error code', papi_errstr
     end if
     if (flops1 < 0) then
       MSG_WARNING("Number of floating point instruction Overflow")
       papi_flops(:)=-1
     end if
   end if
#endif

   select case (option)
   case (0)
     ! Zero out all accumulators of time and init timers
     acctim(:,:)      = 0.0d0
     tzero(:,:)       = 0.0d0
     ncount(:)        = 0
     papi_flops(:)    = 0
     papi_acctim(:,:) = 0.
     papi_accflops(:) = 0.
     papi_tzero(:,:)  = 0.

   case (1)
     ! Initialize timab for nn
     call timein(cpu,wall)
     tzero(1,nn)=cpu
     tzero(2,nn)=wall
#ifdef HAVE_PAPI
     papi_flops(nn)   = flops1       ! Initialize megaflops for nn
     papi_tzero(1,nn) = proc_time
     papi_tzero(2,nn) = real_time
#endif

   case (2)
     ! Accumulate time for nn (also keep the values of cpu, wall, proc_time, real_time, flops1)
     call timein(cpu,wall)
     acctim(1,nn)=acctim(1,nn)+cpu -tzero(1,nn)
     acctim(2,nn)=acctim(2,nn)+wall-tzero(2,nn)
     ncount(nn)=ncount(nn)+1
#ifdef HAVE_PAPI
     ! accumulate time and flops for nn Difference between 2 calls to Papif_flops
     papi_acctim(1,nn)=papi_acctim(1,nn)+ proc_time - papi_tzero(1,nn)
     papi_acctim(2,nn)=papi_acctim(2,nn)+ real_time - papi_tzero(2,nn)
     papi_accflops(nn)=papi_accflops(nn)+ flops1- papi_flops(nn)
#endif

   case (3)
     ! Use previously obtained values to initialize timab for nn
     tzero(1,nn)=cpu
     tzero(2,nn)=wall
#ifdef HAVE_PAPI
     papi_flops(nn)=flops1
     papi_tzero(1,nn) = proc_time
     papi_tzero(2,nn) = real_time
#endif

   case (4)
     ! Return elapsed time for nn (do not accumulate)
     call timein(cpu,wall)
     tottim(1)=cpu-tzero(1,nn)
     tottim(2)=wall-tzero(2,nn)
#ifdef HAVE_PAPI
     ! return elapsed floating point operationfor nn (do not accumulate)
     papi_tottim(1,nn)= proc_time - papi_tzero(1,nn)
     papi_tottim(2,nn)= real_time - papi_tzero(2,nn)
     papi_totflops(nn)= flops1 - papi_flops(nn)
#endif

   case default
     write(message,'(a,i10,a)')'  Input option not valid, =',option,'.'
     MSG_BUG(message)
   end select
 end if

end subroutine timab
!!***

END MODULE m_time

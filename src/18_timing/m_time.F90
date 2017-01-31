!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_time
!! NAME
!! m_time
!!
!! FUNCTION
!!
!! COPYRIGHT
!! Copyright (C) 2009-2016 ABINIT group (MG, XG, MT)
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
 use m_profiling_abi
 use m_errors
 use iso_c_binding
#if defined HAVE_MPI2
 use mpi
#endif

 use m_xpapi,    only: xpapi_flops
 use m_fstrings, only: char_count 

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
!!***

!----------------------------------------------------------------------

CONTAINS  !===========================================================
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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'asctime'
!End of the abilint section

 implicit none

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sec2str'
!End of the abilint section

 implicit none

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

 if (days>0) then
   write(str,'(i0,3(a,i0.2))')days,"-",hours,":",minutes,":",seconds
 else if (hours>0) then
   write(str,'(i0.2,2(a,i0.2))')hours,":",minutes,":",seconds
 else if (minutes>0) then
   write(str,'(i0.2,a,i0.2)')minutes,":",seconds
 else
   write(str,'(i0.2,a)')seconds," [s]"
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

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'str2sec'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_cpu_time'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 real(dp) :: cpu

!Local variables-------------------------------
#ifdef HAVE_FC_CPUTIME
 real :: cpu_sp
#elif defined FC_IBM
 integer :: mclock
#elif defined FC_SUN
 real :: tmp(2)
 real :: etime
#elif defined FC_COMPAQ || defined HAVE_OS_MACOSX 
 real :: tmp(2) !real array only needed by etime
 real(dp) :: etime
#else
 integer :: count_now,count_max,count_rate
#endif

! *************************************************************************

!Machine-dependent timers
#ifdef HAVE_CCLOCK
 call cclock(cpu)

#elif defined HAVE_FC_CPUTIME
!This is the F95 standard subroutine.
 call cpu_time(cpu_sp)
 cpu = cpu_sp

#elif defined FC_IBM
 cpu = mclock()*0.01d0

#elif defined HAVE_OS_MACOSX || defined FC_COMPAQ || defined FC_SUN
 cpu = etime(tmp)

#elif defined FC_FUJITSU
 call clock(cpu,0,2)

#elif defined FC_HITACHI
 call xclock(cpu,5)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abi_wtime'
!End of the abilint section

 implicit none

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
!!      exc_build_ham,lapackprof,m_abilasi,m_bse_io,m_epjdos,m_exc_itdiago
!!      m_fft,m_fft_prof,m_fstab,m_gkk,m_ifc,m_ioarr,m_iowf,m_phgamma,m_phpi
!!      m_shirley,m_sigmaph,m_wfd,m_wfk,partial_dos_fractions
!!
!! CHILDREN
!!      xpapi_flops
!!
!! SOURCE

subroutine cwtime(cpu,wall,gflops,start_or_stop)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cwtime'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(inout) :: cpu,wall
 real(dp),intent(out) :: gflops
 character(len=*),intent(in) :: start_or_stop

!Local variables-------------------------------
#ifndef HAVE_PAPI
 logical,parameter :: use_papi=.FALSE.
#else
 logical,parameter :: use_papi=.TRUE.
#endif
 integer(C_INT)  :: check 
 integer(C_LONG_LONG) :: flops
 real(C_FLOAT) :: real_time,proc_time,mflops

! *************************************************************************

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

 CASE DEFAULT
   MSG_ERROR("Wrong option for start_or_stop: "//TRIM(start_or_stop))
 END SELECT

end subroutine cwtime
!!***

END MODULE m_time

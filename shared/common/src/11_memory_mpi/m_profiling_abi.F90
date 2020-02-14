!!****m* ABINIT/m_profiling_abi
!! NAME
!! m_profiling_abi
!!
!! FUNCTION
!!  This module is used for tracing memory allocations/deallocations
!!  when we compile the code with --enable-memory-profiling="yes" that,
!!  in turn, defines the CPP macro HAVE_MEM_PROFILE in abi_common.h
!!  The main entry point is abimem_init. abimem_record is interfaced via CPP macros
!!  defined in abi_common
!!
!! COPYRIGHT
!! Copyright (C) 2010-2020 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! This module is not thread-safe. However, this should not represent
!! a significant limitation since memory-tracing is only enabled in debug mode.

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_profiling_abi

 use defs_basis
 !use m_clib
#ifdef HAVE_MPI2
 use mpi
#endif

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

#define _ABORT(msg) call abimem_abort(msg, __FILE__, __LINE__)

 public :: abimem_get_info
 public :: abimem_init              ! Initialize memory profiling.
 public :: abimem_set_snapshot_time ! Set time interface for snapshots.
 public :: abimem_shutdown          ! Final cleanup.
 public :: abimem_report            ! Print allocation status.
 public :: abimem_record            ! Central routine to be used for allocation/deallocation.
                                    ! Interfaced via CPP macros defined in abi_common.h

 integer,private,parameter :: slen = 500
 character(fnlen),parameter :: NONE_STRING = "__NONE_STRING__"

!!****t* m_profiling_abi/abimem_t
!! NAME
!!  minfo_t
!!
!! FUNCTION
!!  Internal datastructure storing information on the memory allocated at run-time
!!
!! SOURCE

 type :: abimem_t

   integer :: level = huge(1)
   ! Integer selecting the operation mode
   ! The initial value is set to huge so that main executables that don't call abimem_init will
   ! produce an Error when the first alocation is performed and abimem_record is invoked.

   integer(i8b) :: memory = 0
   ! Total memory allocated so far in bytes.

   integer(i8b) :: peak = 0
   ! Memory peak in bytes.

   integer :: peak_fileline = -1
   ! Line in peak_file

   integer(i8b) :: num_alloc = 0
   ! Total numer of allocations performed so far.

   integer(i8b) :: num_free = 0
   ! Total numer of deallocations performed so far.

   integer :: logunt = 99
   ! Unit number of logfile (hardcoded)

   integer :: my_rank = 0
   ! Rank of this processor.

   logical :: iwrite = .False.
   ! True if this MPI rank should write to file.

   !real(dp),private,save :: start_time
   ! Origin of time in seconds.

   real(dp) :: last_snapshot = -one
   ! time of the last snapshot in seconds.

   real(dp) :: dt_snapshot = -one
   ! time between two consecutive snapshots in seconds.

   real(dp) :: limit_mb = 20_dp
   ! Optional memory limit in Mb. used when level == 3

   character(len=slen) :: peak_vname = "_vname"
   ! Name of the last variable for which the memory peak occurred.

   character(len=slen) :: peak_file = NONE_STRING
   ! Name of the file in which peak occurred.

   ! Selective memory tracing
   character(fnlen) :: select_file = NONE_STRING

   character(len=fnlen) :: logfile
   ! File used for logging allocations/deallocations.

 end type abimem_t
!!***

 type(abimem_t),private,save :: minfo
 ! Internal datastructure storing memory profiling data.

contains

!!****f* m_profiling_abi/abimem_init
!! NAME
!! abimem_init
!!
!! FUNCTION
!!  Initialize memory profiling module.
!!
!! INPUT
!!  level = Integer selecting the operation mode:
!!       0 -> no file abimem.mocc is created, only memory allocation counters running
!!       1 -> light version. Only memory peaks are written.
!!       2 -> file abimem.mocc is created with full information inside.
!!       3 -> Write info only if allocation/deallocation is larger that limit_mb
!!    NOTE: By default, only master node writes, use negative values to make all MPI procs write info to disk.
!!  [delta_time]=Interval in second for snapshots. Will write report to std_out evety delta_time seconds.
!!  [filename] = If present, activate memory logging only inside filename (basename).
!!  [limit_mb]= Set memory limit in Mb if level == 3. Print allocation/deallocation only above this limit.

subroutine abimem_init(level, delta_time, filename, limit_mb)

!Arguments ------------------------------------
 integer, intent(in) :: level
 real(dp),optional,intent(in) :: delta_time
 real(dp),optional,intent(in) :: limit_mb
 character(len=*),optional,intent(in) :: filename

!Local variables-------------------------------
 integer :: ierr
 logical :: file_exists
 character(len=500) :: msg
! *************************************************************************

 minfo%level = level
 !start_time = abimem_wtime()

 ! Optionally, selects functions or files to be profiled.
 if (present(filename)) then
   if (len_trim(filename) > 0) minfo%select_file = filename
 end if

 ! Optionally, set max limit in Mb used if level == 2
 if (present(limit_mb)) minfo%limit_mb = limit_mb

 ! Build name of file used for logging.
 minfo%my_rank = 0
#if defined HAVE_MPI
 call MPI_COMM_RANK(MPI_COMM_WORLD, minfo%my_rank, ierr)
#endif
 write(minfo%logfile,"(a,i0,a)")"abimem_rank",minfo%my_rank,".mocc"

 ! Clean the file if it already exists.
 inquire(file=minfo%logfile, exist=file_exists)
 if (file_exists) then
   open(unit=minfo%logunt, file=minfo%logfile, status="old", iostat=ierr)
   if (ierr==0) close(unit=minfo%logunt, status="delete", iostat=ierr)
 end if

 ! Activate snapshots.
 if (present(delta_time)) then
   minfo%dt_snapshot = delta_time
   minfo%last_snapshot = zero
   if (delta_time < 1.0e-6) then
     _ABORT("delta_time is too small")
   end if
 end if

 select case (abs(minfo%level))
 case (0)
   ! No action required

 case (1, 2, 3)
   minfo%iwrite = .False.
   if (minfo%my_rank == 0) minfo%iwrite = .True.
   if (minfo%level < 0) minfo%iwrite = .False.
   if (minfo%iwrite) then
     open(unit=minfo%logunt, file=minfo%logfile, status='unknown', action='write', iostat=ierr)
     if (ierr /= 0) then
       _ABORT("Opening abimem file")
     end if
     if (minfo%level == 1) call write_header("# Write memory allocations larger than previous peak")
     if (minfo%level == 2) call write_header("# To be used for inspecting a variable which is not deallocated")
     if (minfo%level == 3) then
       write(msg, "(a,f9.1,a)")"# Write memory allocations/deallocations larger than ", minfo%limit_mb, "(MB)"
       call write_header(msg)
     end if
  end if

 case default
   write(msg, "(a,i0,2a)") &
    "Invalid value for abimem_level:", minfo%level, ch10, &
    "Make sure you are calling abimem_init and abinit_doctor in main!"
   _ABORT(msg)
 end select

 contains

 subroutine write_header(info)
   character(len=*),intent(in) :: info
   write(minfo%logunt, "(a, i0)")"# memocc file generated by Abinit compiled with HAVE_MEM_PROFILE."
   write(minfo%logunt, "(a)") trim(info)
   write(minfo%logunt, "(2(a,i0),a)")"# {level: ", level, ", rank: ", minfo%my_rank, "}"
   write(minfo%logunt,'(a,t60,a)')&
       '# Variable name', 'Action Address Size[b] File Line Total Memory [bits]'
 end subroutine write_header

end subroutine abimem_init
!!***

!!****f* m_profiling_abi/abimem_set_snapshot_time
!! NAME
!! abimem_set_snapshot_time
!!
!! FUNCTION
!!  Set time interface for snapshots.
!!
!! INPUT
!!  delta_time=Interval in second for snapshots.
!!
!! INPUT

subroutine abimem_set_snapshot_time(delta_time)

!Arguments ------------------------------------
 real(dp),optional,intent(in) :: delta_time
! *************************************************************************

 minfo%dt_snapshot = delta_time
 minfo%last_snapshot = zero
 if (delta_time < 1.0e-6) then
   _ABORT("delta_time is too small")
 end if

end subroutine abimem_set_snapshot_time
!!***

!!****f* m_profiling_abi/abimem_shutdown
!! NAME
!! abimem_shutdown
!!
!! FUNCTION
!! Perform final cleanup of the module and close files.
!!
!! INPUT

subroutine abimem_shutdown()

!Local variables-------------------------------
 integer :: unt_found
 logical :: isopen
! *************************************************************************

 minfo%level = 0

 ! Close the file if it's connected
 inquire(file=minfo%logfile, number=unt_found, opened=isopen)
 if (isopen .and. unt_found == minfo%logunt) close(minfo%logunt)

end subroutine abimem_shutdown
!!***

!!****f* m_profiling_abi/abimem_report
!! NAME
!! abimem_report
!!
!! FUNCTION
!!  Print info about memory usage to unit `unt`.
!!  Add mallinfo values if `with_mallinfo` (default: False)
!!
!! INPUT
!!

subroutine abimem_report(unt, with_mallinfo)

!Arguments ------------------------------------
 integer,intent(in) :: unt
 logical,optional,intent(in) :: with_mallinfo

!Local variables-------------------------------
 integer,save :: icall = 0
 integer(i8b),save :: prev_memory
 real(dp) :: diff_mb

! *************************************************************************

 if (minfo%level == huge(one)) return
 icall = icall + 1
 write(unt,"(a)")"------------------------- MEMORY CONSUMPTION REPORT -----------------------------"
 write(unt,"(3(a,i0))")" Malloc: ",minfo%num_alloc,", Free: ", minfo%num_free, ", M-F: ", minfo%num_alloc - minfo%num_free
 write(unt,"(a,f8.1,a)")" Memory allocated so far: ", minfo%memory * b2Mb, " (Mb)"
 write(unt,"(a,f8.1,5a,i0)")" Peak: ", minfo%peak * b2Mb," (MB) for variable: ", trim(minfo%peak_vname), &
   "at:", trim(abimem_basename(minfo%peak_file)),":",minfo%peak_fileline
 diff_mb = zero; if (icall > 1) diff_mb = (minfo%memory - prev_memory) * b2Mb
 write(unt,"(a,f8.1,a)")" Memory allocated wrt previous call: ", diff_mb, " (Mb)"
 prev_memory = minfo%memory

 if (present(with_mallinfo)) then
   !if (with_mallinfo) call clib_print_mallinfo(unit=unt)
 end if

end subroutine abimem_report
!!***

!!****f* m_profiling_abi/abimem_get_info
!! NAME
!! abimem_get_info
!!
!! FUNCTION
!!  Function that returns the number of allocations and deallocations that have
!!  been performed in Fortran and the memory currently used
!!
!! OUTPUT
!!  nalloc: number of allocations that have been done
!!  nfree:  number of deallocations that have been done
!!  allocmemory:  total memory used

subroutine abimem_get_info(nalloc, nfree, allocmemory)

!Arguments ------------------------------------
 integer(i8b),intent(out) :: nalloc, nfree, allocmemory
! *************************************************************************

 nalloc = minfo%num_alloc; nfree = minfo%num_free; allocmemory = minfo%memory

end subroutine abimem_get_info
!!***

!!****f* m_profiling_abi/abimem_record
!! NAME
!! abimem_record
!!
!! FUNCTION
!!  Control the memory occupation by calculating the overall size of the allocated arrays
!!  At the end of the calculation a short report is printed on the screen,
!!  some information can be also written on disk following the needs
!!
!! PARENTS
!!
!! CHILDREN
!!      date_and_time,mpi_abort
!!
!! SOURCE

subroutine abimem_record(istat, vname, addr, act, isize, file, line)

!Arguments ------------------------------------
 integer,intent(in) :: istat,line
 integer(i8b), intent(in) :: isize,addr
 character(len=*), intent(in) :: vname,act,file

!Local variables-------------------------------
 !integer :: ierr
 !real(dp) :: now
 logical :: do_log, new_peak
 character(len=500) :: msg
! *************************************************************************

 ! Handle possible allocate/deallocate failures.
 if (istat /= 0) then
   if (isize >= 0) then
     write(msg,*)" Problem of allocation of variable: ",trim(vname),', error code= ',istat
     _ABORT(msg)
   else if (isize<0) then
     write(msg,*)" Problem of deallocation of variable ",trim(vname),', error code= ',istat
     _ABORT(msg)
   end if
 end if

 ! Increase total counter
 minfo%memory = minfo%memory + isize
 new_peak = .False.
 if (isize > minfo%peak) then
   ! New peak
   new_peak = .True.
   minfo%peak = isize
   minfo%peak_vname = vname
   minfo%peak_file = file
   minfo%peak_fileline = line
 end if

 if (isize > 0) then
   minfo%num_alloc = minfo%num_alloc + 1
 else if (isize < 0) then
   minfo%num_free = minfo%num_free + 1
 else
   ! This is the correct check but tests fail!
   !if (act == "A") then
   !  minfo%num_alloc = minfo%num_alloc + 1
   !else if (act == "D") then
   !  minfo%num_free = minfo%num_free + 1
   !else
   !  _ABORT("Wrong action: "//trim(act))
   !end if
 end if

 ! Selective memory tracing
 do_log = .True.
 if (minfo%select_file /= NONE_STRING) do_log = (minfo%select_file == file)
 do_log = do_log .and. minfo%iwrite

 ! Snapshot (write to std_out)
 !if (do_log .and. minfo%last_snapshot >= zero) then
 !  now = abimem_wtime()
 !  if (now - minfo%last_snapshot >= minfo%dt_snapshot) then
 !    call abimem_report(std_out)
 !    minfo%last_snapshot = now
 !  end if
 !end if

 ! IMPORTANT:
 ! Remember to change the pyton code in ~abinit/tests/pymods/memprof.py to account for changes in the format
 if (do_log) then
   select case (minfo%level)
   case (0)
     ! No action required

   case (1)
     ! Write only if we have a new peak
     if (new_peak) then
       write(minfo%logunt,'(a,t60,a,1x,2(i0,1x),a,1x,2(i0,1x))') &
         trim(vname), trim(act), addr, isize, trim(abimem_basename(file)), line, minfo%memory
     end if

   case (2)
     ! To be used for inspecting a variable which is not deallocated
     write(minfo%logunt,'(a,t60,a,1x,2(i0,1x),a,1x,2(i0,1x))') &
       trim(vname), trim(act), addr, isize, trim(abimem_basename(file)), line, minfo%memory

   case (3)
     ! Write memory allocations larger than limit_mb
     if (abs(isize * b2Mb) > minfo%limit_mb) then
       write(minfo%logunt,'(a,t60,a,1x,2(i0,1x),a,1x,2(i0,1x))') &
         trim(vname), trim(act), addr, isize, trim(abimem_basename(file)), line, minfo%memory
     end if

   case default
     _ABORT("Invalid abimem_level")
   end select
 end if

end subroutine abimem_record
!!***


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Private routine providing services already implemented in other higher level modules.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!****f* m_abimem/abimem_abort
!! NAME
!!  abimem_abort
!!
!! FUNCTION
!!  Stop the code if an error occurs.
!!
!! INPUT
!!  msg=Error message
!!  file=File name
!!  line=Line number
!!
!! PARENTS
!!
!! CHILDREN
!!      date_and_time,mpi_abort
!!
!! SOURCE

subroutine abimem_abort(msg, file, line)

!Arguments ------------------------------------
 integer,intent(in) :: line
 character(len=*),intent(in) :: msg,file

!Local variables-------------------------------
 integer :: ierr

!Local variables-------------------------------
 integer :: unt_found
 logical :: isopen
! *************************************************************************

 write(std_out,*)trim(msg),", file: ", trim(file), ", line: ", line

 ! Close logfile if it's connected to flush io buffers and avoid file corruption
 inquire(file=minfo%logfile, number=unt_found, opened=isopen)
 if (isopen .and. (unt_found == minfo%logunt)) close(unit=minfo%logunt)

 ierr = 0
#ifdef HAVE_MPI
 call MPI_ABORT(MPI_COMM_WORLD, MPI_ERR_UNKNOWN, ierr)
#endif
 stop

end subroutine abimem_abort
!!***

!----------------------------------------------------------------------

!!****f* m_profiling_abi/abimem_basename
!! NAME
!! abimem_basename
!!
!! FUNCTION
!!  Returns the final component of a pathname.
!!
!! INPUTS
!!  string=The input string
!!
!! NOTES
!!  * If the input string in not a valid path to a file (i.e not in the form foo/name)
!!    a blank strink is returned
!!  * We do a backward search becase we want to optimize the algorithm for Fortran strings.
!!
!! SOURCE

pure function abimem_basename(string) result(basename)

 character(len=*),intent(in) :: string
 character(len=LEN_TRIM(string)) :: basename

!Local variables-------------------------------
 integer :: ic,nch_trim,nch
 character(len=1),parameter :: DIR_SEPARATOR = '/'
 character(len=1),parameter :: BLANK=' '
!************************************************************************

 nch     =LEN     (string)
 nch_trim=LEN_TRIM(string)

 ic = INDEX (TRIM(string), DIR_SEPARATOR, back=.TRUE.)
 !write(*,*)'DEBUG ',TRIM(string),ic

 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the separator.
   basename = string(ic+1:nch_trim)
   return
 else if (ic==0 .or. ic == nch_trim+1) then ! no separator in string or zero length string,
   basename = TRIM(string)                   ! return trimmed string.
   return
 else              ! (ic == nch_trim) separator is the last char.
   basename= BLANK  ! This is not a valid path to a file, return blank.
   return
 end if

end function abimem_basename
!!***

!----------------------------------------------------------------------

!!****f* m_time/abimem_wtime
!! NAME
!!  abimem_wtime
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

function abimem_wtime() result(wall)

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
     _ABORT('Problem with month and year numbers.')
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

end function abimem_wtime
!!***

end module m_profiling_abi
!!***

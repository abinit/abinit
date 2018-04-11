!!****m* ABINIT/m_profiling_abi
!! NAME
!! m_profiling_abi
!!
!! FUNCTION
!!  This module is used for tracing memory allocations/deallocations
!!  when we compile the code with --enable-memory-profiling="yes" that, in turns, defines HAVE_MEM_PROFILE.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  *  This module is not thread-safe. However, this should not represent
!!     a significant limitation since memory-tracing is only enabled in debug mode.

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_profiling_abi

 use defs_basis
#ifdef HAVE_MPI2
 use mpi
#endif

 implicit none

 private
!!***

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

#ifdef HAVE_MEM_PROFILING

#define _ABORT(msg) call abimem_abort(msg, __FILE__, ABI_FUNC, __LINE__)

 public :: abimem_get_info
 public :: abimem_init
 public :: abimem_shutdown
 public :: abimem_set_opts
 public :: abimem_report
 public :: abimem_record               ! Central routine to be used for allocation/deallocation
 !public :: abimem_enable
 !public :: abimem_disable
 !public :: abimem_reset

 integer,private,parameter :: slen = 500

!!****t* m_profiling_abi/abimem_t
!! NAME
!!  minfo_t
!!
!! FUNCTION
!!  Store information on the memory allocated at run-time
!!
!! SOURCE

 !Memory profiling
 type :: abimem_t
   integer(kind=8) :: memory = int(0, kind=8)
   integer(kind=8) :: peak = int(0, kind=8)
   character(len=slen) :: func = "_func"
   character(len=slen) :: vname = "_vname"
 end type abimem_t
!!***

 ! PRIVATE STUFF
 ! Selective memory tracing
 character(fnlen),parameter :: NONE_STRING = "__NONE_STRING__"
 character(fnlen),save :: select_file = NONE_STRING
 character(fnlen),save :: select_func = NONE_STRING

 ! Save values for memocc_abi.
 real(dp),save :: abimem_limit = zero
 logical,save :: abimem_isinit = .False.
 integer,parameter :: logunt = 99
 character(len=fnlen),save :: abimem_file
 type(abimem_t),save :: memloc_abi, memtot_abi
 integer,save :: num_alloc = 0
 integer,save :: num_free = 0
 integer,save :: my_rank = 0
 !Debug option for memocc_abi, set in the input file
 !logical,parameter :: abimem_debug=.True.
 !logical,save :: abimem_ilog = .False.
 integer,save :: abimem_level = 0

 !real(dp),private,save :: start_time
! Origin of time in seconds.
 real(dp),private,save :: last_snapshot = -one
! time of the last snapshot in seconds.
 real(dp),private,save :: dt_snapshot = -one
! time between two consecutive snapshots in seconds.

contains

!!****f* m_profiling_abi/abimem_init
!! NAME
!! abimem_init
!!
!! FUNCTION
!!
!! INPUT
!!  level = Integer selecting the operation mode:
!!       0 no file abimem.mocc is created, only memory allocation counters running
!!       1 file abimem.mocc is created in a light version (only current information is written)
!!       2 file abimem.mocc is created with full information inside (default state if not specified)
!!       The status can only be downgraded. A stop signal is produced if status is increased
!!  [deltat]
!!  [filename] = If present, activate memory logging only inside filename.
!!  [funcname] = If present, activate memory logging only inside funcname.

 subroutine abimem_init(level, deltat, filename, funcname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_init'
!End of the abilint section

  implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_init'
!End of the abilint section

!Arguments ------------------------------------
 integer, intent(in) :: level
 real(dp),optional,intent(in) :: deltat
 character(len=*),optional,intent(in) :: filename, funcname

!Local variables-------------------------------
 integer :: ierr
 logical :: file_exists
! *************************************************************************

 !if (level > abimem_level) stop 'abimem_level can be only downgraded'
 abimem_level = level
 abimem_isinit = .True.
 !start_time = abimem_wtime()

 ! Build name of file used for logging.
 my_rank = 0
#if defined HAVE_MPI
 call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
#endif
 write(abimem_file,"(a,i0,a)")"abimem_rank",my_rank,".mocc"

 ! Optionally, selects functions or files to be profiled.
 if (present(filename)) then
   if (len_trim(filename) > 0) select_file = filename
 end if
 if (present(funcname)) then
    if (len_trim(funcname) > 0) select_func = funcname
 endif

 ! Clean the file if it already exists.
 ! The file should be deleted
 inquire(file=abimem_file, exist=file_exists)
 if (file_exists) then
   open(unit=logunt, file=abimem_file, status="old", iostat=ierr)
   if (ierr==0) close(unit=logunt, status="delete", iostat=ierr)
 end if

 ! Activate snapshots
 if (present(deltat)) then
   dt_snapshot = deltat
   last_snapshot = zero
   if (deltat < 1.0e-6) then
     _ABORT("deltat is too small")
   end if
 end if

 select case (level)
 case (0)
   ! No action required

 case (1)
   open(unit=logunt, file=abimem_file, status='unknown', action='write', iostat=ierr)
   if (ierr /= 0) then
     _ABORT("Opening abimem file")
   end if
   write(logunt,'(a,t60,a,t90,4(1x,a12))')&
      '(Data in KB) Routine','Array name    ','Array size','Total Memory'

 case (2)
   open(unit=logunt, file=abimem_file, status='unknown', action='write', iostat=ierr)
   if (ierr /= 0) then
     _ABORT("Opening abimem file")
   end if
   write(logunt,'(a,t60,a,t90,4(1x,a12))')&
      '(Data in KB) Routine','Array name    ','Array size','Total Memory'

 case default
   _ABORT("invalid abimem_level")
 end select

end subroutine abimem_init
!!***

!!****f* m_profiling_abi/abimem_shutdown
!! NAME
!! abimem_shutdowns
!!
!! FUNCTION
!! Perform final cleanup of the module and close files.
!!
!! INPUT

subroutine abimem_shutdown()


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_shutdown'
!End of the abilint section

  implicit none

!Local variables-------------------------------
 integer :: unt_found
 logical :: isopen
! *************************************************************************

 abimem_level = 0; abimem_isinit = .False.

 ! Close the file if it's connected
 inquire(file=abimem_file, number=unt_found, opened=isopen)
 if (isopen .and. (unt_found==logunt)) close(logunt)

end subroutine abimem_shutdown
!!***

!!****f* m_profiling_abi/abimem_set_opts
!! NAME
!! abimem_set_opts
!!
!! FUNCTION
!!
!! INPUT
!!  level = Integer selecting the operation mode:
!!       0 no file abimem.mocc is created, only memory allocation counters running
!!       1 file abimem.mocc is created in a light version (only current information is written)
!!       2 file abimem.mocc is created with full information inside (default state if not specified)
!!       The status can only be downgraded. A stop signal is produced if status is increased
!!  [limit]= Give a memory limit above which the code will stop properly. The unit is bytes.
!!  [filename] = If present, activate memory logging only inside filename.
!!  [funcname] = If present, activate memory logging only inside funcname.
!!  [deltat]


subroutine abimem_set_opts(level, limit, deltat, filename, funcname)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_set_opts'
!End of the abilint section

  implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_setlogopts'
!End of the abilint section

!Arguments ------------------------------------
 integer, intent(in) :: level
 real(dp),optional,intent(in) :: limit
 real(dp),optional,intent(in) :: deltat
 character(len=*),optional,intent(in) :: filename, funcname

! *************************************************************************

 abimem_level = level

 ! Optionally, set max limit on total memory allocated
 if (present(limit)) abimem_limit = limit

 ! Optionally, selects functions or files to be profiles.
 if (present(filename)) then
   if (len_trim(filename) > 0) select_file = filename
 end if
 if (present(funcname)) then
    if (len_trim(funcname) > 0) select_func = funcname
 endif

 ! Activate snapshots
 if (present(deltat)) then
   dt_snapshot = deltat
   last_snapshot = zero
   if (deltat < 1.0e-6) then
     _ABORT("deltat is too small")
   end if
 end if

end subroutine abimem_set_opts
!!***

subroutine abimem_report(unit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_report'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_report'
!End of the abilint section

!Arguments ------------------------------------
 integer,optional,intent(in) :: unit

!Local variables-------------------------------
 integer :: unt
! *************************************************************************

 unt = std_out; if (present(unit)) unt = unit

#if 0
 if (trim(func)=='stop' .and. my_rank == 0) then
   if (abimem_level > 0) then
     if (abimem_level == 1) rewind(logunt)
     write(logunt,'(a,t60,a,t90,4(1x,i12))')&
       trim(memloc_abi%func),trim(memloc_abi%vname),&
       memloc_abi%memory/int(1024,kind=8),memloc_abi%peak/int(1024,kind=8),&
       memtot_abi%memory/int(1024,kind=8),&
       (memtot_abi%peak+memloc_abi%peak-memloc_abi%memory)/int(1024,kind=8)
     close(unit=logunt)
   end if

   write(unt,'(1x,a)')'-------------------------MEMORY CONSUMPTION REPORT-----------------------------'
   write(unt,'(1x,2(i0,a,1x),i0)')&
        num_alloc,' allocations and',num_free,' deallocations, remaining memory(B):',memtot_abi%memory
   write(unt,'(1x,a,i0,a)') 'memory occupation peak: ',memtot_abi%peak/int(1048576,kind=8),' MB'
   write(unt,'(4(1x,a))') 'for the variable ',trim(memtot_abi%vname),'in the function',trim(memtot_abi%func)
   !here we can add a func which open the abimem.mocc file in case of some
   !memory allocation problem, and which eliminates it for a successful run
   f (abimem_level == 1 .and. num_alloc == num_free .and. memtot_abi%memory==int(0,kind=8)) then
   !remove file should be put here
   open(unit=logunt,file=abimem_file,status='unknown',action='write')
   write(unit=logunt,fmt='()',advance='no')
   close(unit=logunt)
   else
     call abimem_check(num_alloc,num_free)
   end if

 else if (trim(func)/='stop') then
   write(unt,*) "memocc_abi: ",array," ",func
   write(unt,"(a,i0,a)") "Error[",my_rank,"]: Use memocc_abi and the word 'count' only with the word 'stop'."
   _ABORT("Exit requested by user")
 end if
#endif

end subroutine abimem_report
!!***

!!****f* m_profiling_abi/abimem_get_info
!! NAME
!! abimem_get_info
!!
!! FUNCTION
!!  Function that returns the number of allocations and deallocations that have
!!  been done and the memory currently used
!!
!! OUTPUT
!!  nalloc       number of allocations that have been done
!!  ndealloc     number of deallocations that have been done
!!  allocmemory  total memory used

subroutine abimem_get_info(nalloc, ndealloc, allocmemory)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_get_info'
!End of the abilint section

 implicit none

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_get_info'
!End of the abilint section

!Arguments ------------------------------------
 integer(kind=8), intent(out) :: allocmemory
 integer, intent(out) :: nalloc,ndealloc
! *************************************************************************

 nalloc = num_alloc; ndealloc = num_free; allocmemory = memtot_abi%memory

end subroutine abimem_get_info
!!***

!!****f* m_profiling_abi/abimem_record
!! NAME
!! abimem_record
!!
!! FUNCTION
!!  Control the memory occupation by calculating the overall size of the allocated arrays
!!
!!  At the end of the calculation a short report is printed on the screen,
!!  some information can be also written on disk following the needs
!!
!!  The file abimem.mocc is not deleted if the final total memory is not equal to zero.
!!  abimem_debug (parameter)
!!    == .true.  verbose format (useful with tests/scripts/abimem.py)
!!               then display a line per allocation or deallocation
!!               a routine at the end parses the file
!!    == .false. compact format
!!
!! PARENTS
!!      abinit,aim,anaddb,band2eps,conducti,cut3d,dummy_tests,fftprof
!!      fold2Bloch,ioprof,lapackprof,macroave,mrgddb,mrgdv,mrggkk,mrgscr
!!      multibinit,optic,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!      date_and_time,mpi_abort
!!
!! SOURCE

subroutine abimem_record(istat, vname, addr, act, isize, file, func, line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_record'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: istat,line
 integer(kind=8), intent(in) :: isize,addr
 character(len=*), intent(in) :: vname,act,file,func

!Local variables-------------------------------
 integer :: ierr
 real(dp) :: now
 logical :: do_log
 character(len=500) :: msg
! *************************************************************************

 ! Handle allocate/deallocate failures
 !if (istat /= 0) then
 !  write(msg,('(5a,i0,/,2a,/,a,i0,a)'))&
 !   'Procedure: ',trim(funcname),"@",trim(filename),":",line,&
 !   'problem of allocation of variable: ',trim(arr_name),&
 !   'Error code = ',istat,' Aborting now...'
 !  call abimem_abort(istat, msg, filename, funcname, line)
 ! end if

 !control of the allocation/deallocation status
 if (istat /= 0) then
   if (isize >= 0) then
     write(msg,*)&
       trim(func),': problem of allocation of variable ',trim(vname),', error code= ',istat
     _ABORT(msg)
   else if (isize<0) then
     write(msg,*)&
       trim(func),': problem of deallocation of variable ',trim(vname),', error code= ',istat
     _ABORT(msg)
   end if
 end if

 ! total counter, for all the processes
 memtot_abi%memory = memtot_abi%memory + isize
 if (memtot_abi%memory > memtot_abi%peak) then
   memtot_abi%peak = memtot_abi%memory
   memtot_abi%func = func
   memtot_abi%vname = vname
 end if

 if (isize > 0) then
   num_alloc = num_alloc + 1
 else if (isize < 0) then
   num_free = num_free + 1
 end if

 ! This is the correct check but tests fail!
 !if (act == "A") then
 !  num_alloc = num_alloc + 1
 !else if (act == "D") then
 !  num_free = num_free + 1
 !else
 !  _ABORT("Wrong action: "//trim(act))
 !end if

 ! Check on memory limit.
 if (abimem_limit /= zero .and. memtot_abi%memory > int(real(abimem_limit,kind=8)*1073741824.d0,kind=8)) then
   ! memory limit is in GB
   write(msg,'(a,f7.3,2(a,i0),a,2(a,i0))')&
     'Memory limit of ',abimem_limit,' GB reached for rank ',my_rank,', total memory is ',memtot_abi%memory,' B.',&
     'this happened for variable '//trim(memtot_abi%vname)//' in func '//trim(memtot_abi%func)
   _ABORT(msg)
 end if

 ! Selective memory tracing
 do_log = .True.
 if (select_file /= NONE_STRING) do_log = (select_file == file)
 if (select_func /= NONE_STRING) do_log = do_log .and. (select_func == func)
 !do_log = (do_log .and. my_rank == 0)

 ! Snapshot
 if (do_log .and. last_snapshot >= zero) then
   now = abimem_wtime()
   if ((now - last_snapshot) >= dt_snapshot) then
      last_snapshot = now
   else
      do_log = .False.
   end if
 end if

 if (do_log) then
   select case (abimem_level)
   case (0)
     ! No action required

   case (1)
     ! Compact format
     if (trim(memloc_abi%func) /= func) then
       if (memloc_abi%memory /= int(0,kind=8)) then
         rewind(logunt)
         write(logunt,'(a,t60,a,t90,4(1x,i12))')&
           trim(memloc_abi%func),trim(memloc_abi%vname),&
           memloc_abi%memory/int(1024,kind=8),memloc_abi%peak/int(1024,kind=8),&
           memtot_abi%memory/int(1024,kind=8),&
           (memtot_abi%memory+memloc_abi%peak-memloc_abi%memory)/int(1024,kind=8)
       end if
       memloc_abi%func = func
       memloc_abi%vname = vname
       memloc_abi%memory = isize
       memloc_abi%peak = isize

     else
       memloc_abi%memory=memloc_abi%memory+isize
       if (memloc_abi%memory > memloc_abi%peak) then
         memloc_abi%peak = memloc_abi%memory
         memloc_abi%vname = vname
       end if
     end if

   case (2)
     !to be used for inspecting a variable which is not deallocated
     write(logunt,'(a,t60,a,1x,2(i0,1x),2(a,1x),2(i0,1x))')&
       trim(vname), trim(act), addr, isize, trim(abimem_basename(file)), trim(func), line, memtot_abi%memory

   case default
     _ABORT("invalid abimem_level")
   end select
 end if

end subroutine abimem_record
!!***

!!****f* m_profiling_abi/abimem_check
!! NAME
!! abimem_check
!!
!! FUNCTION
!!   Check the abimem.mocc file (verbose format)
!!
!! PARENTS
!!      m_profiling_abi
!!
!! CHILDREN
!!      date_and_time,mpi_abort
!!
!! SOURCE

subroutine abimem_check(nalloc, ndealloc)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_check'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nalloc,ndealloc
! *************************************************************************

 if (abimem_level==2 .and. nalloc /= ndealloc) then
   write(std_out,"(3a)")"Use the python script 'abimem.py' in ~abinit/tests/scripts to check ",trim(abimem_file)," file"
   write(std_out,"(a)")"Note that one can use the command line option `abinit --abimem-level=2"
 end if

end subroutine abimem_check
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
!!  func=Function name.
!!  line=Line number
!!
!! PARENTS
!!
!! CHILDREN
!!      date_and_time,mpi_abort
!!
!! SOURCE

subroutine abimem_abort(msg, file, func, line)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_abort'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(in) :: line
 character(len=*),intent(in) :: msg,file,func

!Local variables-------------------------------
 integer :: ierr

!Local variables-------------------------------
 integer :: unt_found
 logical :: isopen
! *************************************************************************

 write(std_out,*)msg,file,func,line

 ! Close abimem_file if it's connected to flush io buffers and avoid file corruption
 inquire(file=abimem_file, number=unt_found, opened=isopen)
 if (isopen .and. (unt_found == logunt)) close(unit=logunt)

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_basename'
!End of the abilint section

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


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'abimem_wtime'
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
     _ABORT(msg)
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

#endif
! HAVE_MEM_PROFILING
end module m_profiling_abi
!!***

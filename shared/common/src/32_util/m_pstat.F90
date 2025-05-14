!!****m* ABINIT/m_pstat
!! NAME
!! m_pstat
!!
!! FUNCTION
!! Interface to the /proc/{pid}/status file available on Linux.
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_pstat

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors
 use m_yaml

 use m_fstrings, only : find_and_select, basename
 use m_clib,     only : clib_getpid

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_pstat/pstat_t
!! NAME
!! pstat_t
!!
!! FUNCTION
!! This object stores the most important quantites reported in the /proc/{pid}/status file
!! in particular the virtual memory VmRSS. See https://docs.kernel.org/filesystems/proc.html
!!
!! NB: This file is only available on Linux hence one should always check the value of pstat%ok
!! before using quantities such as vmrss_mb.
!!
!! SOURCE

 type, private :: pstat_t

  logical :: ok = .False.
  ! False if stat file is not available

  integer :: pid = -1
  ! Process identifier.

  integer :: threads = -1
  ! number of threads.

  integer :: fdsize = -1
   ! Number of file descriptor slots currently allocated.

  real(dp) :: vmrss_mb = -one
  ! Actual physical RAM used. It contains the three following parts (VmRSS = RssAnon + RssFile + RssShmem)

  real(dp) :: vmpeak_mb = -one
  ! Peak virtual memory size

  real(dp) :: vmstk_mb = -one
  ! Size of stack segments

  character(len=fnlen) :: filepath = ""
  ! Path of status file

  character(len=500) :: iomsg = ""
  ! Error message returned when parsing filepath

 contains
   procedure :: from_pid => pstat_from_pid     ! Init object from process identifier (main entry point).
   procedure :: from_file => pstat_from_file   ! Init object from file (useful for debugging).
   procedure :: print => pstat_print           ! Print object.
   procedure :: min_mem_mb_per_proc => pstat_min_mem_mb_per_proc
 end type pstat_t
!!***

 type(pstat_t),save,public :: pstat_proc

!----------------------------------------------------------------------

contains
!!***

!!****f* m_pstat/pstat_from_pid
!! NAME
!!  pstat_from_pid
!!
!! FUNCTION
!!   Init object from process identifier (main entry point for client code).
!!
!! SOURCE

subroutine pstat_from_pid(pstat)

!Arguments ------------------------------------
 class(pstat_t),intent(out) :: pstat

!Local variables-------------------------------
 integer(c_int) :: pid
 character(len=500) :: spid
! *************************************************************************

 pid = clib_getpid()
 write(spid, "(i0)") pid
 spid = adjustl(spid)
 call pstat%from_file("/proc/"//trim(spid)//"/status")

end subroutine pstat_from_pid
!!***

!!****f* m_pstat/pstat_from_file
!! NAME
!!  pstat_from_file
!!
!! FUNCTION
!! Init object from file (useful for debugging).
!!
!! SOURCE

subroutine pstat_from_file(pstat, filepath)

!Arguments ------------------------------------
 class(pstat_t),intent(inout) :: pstat
 character(len=*),intent(in) :: filepath

!Local variables-------------------------------
 integer :: unit, ierr
 character(len=500) :: line
 integer :: istart, istop, iostat
! *************************************************************************

 pstat%ok = .False.
 pstat%filepath = filepath

 open(newunit=unit, file=trim(pstat%filepath), action="read", status="old", iostat=ierr, iomsg=pstat%iomsg)
 if (ierr /= 0) then
   close(unit); return
 end if

 do
   read(unit, "(a)", iostat=ierr, end=10, iomsg=pstat%iomsg) line
   if (ierr > 0) then ! EOF
     close(unit); return
   end if

   ! Parse useful integers
   if (index(line, "Pid:") == 1) call get_int(line, pstat%pid)
   if (index(line, "Threads:") == 1) call get_int(line, pstat%threads)
   if (index(line, "FDSize:") == 1) call get_int(line, pstat%fdsize)

   ! Parse memory entries
   if (index(line, "VmRSS:") == 1) call get_mem_mb(line, pstat%vmrss_mb)
   if (index(line, "VmPeak:") == 1) call get_mem_mb(line, pstat%vmpeak_mb)
   if (index(line, "VmStk:") == 1) call get_mem_mb(line, pstat%vmstk_mb)
 end do

10 close(unit)
  pstat%ok = .True.
  pstat%iomsg = ""

contains

subroutine get_mem_mb(str, mem_mb)

 character(len=*),intent(in) :: str
 real(dp),intent(out) :: mem_mb

 ! Generic mem entry has format `VmRSS: 2492 kB`
 real(dp) :: mem_fact
 istart = index(str, ":") + 1
 istop = find_and_select(str, &
                        ["kB", "mB"], &
                        [one/1024._dp, one], mem_fact, pstat%iomsg) !default=one,
 ABI_CHECK(istop /= -1, pstat%iomsg)
 read(str(istart+1:istop-1), fmt=*, iostat=iostat, iomsg=pstat%iomsg) mem_mb
 ABI_CHECK(iostat == 0, pstat%iomsg)
 mem_mb = mem_mb * mem_fact

end subroutine get_mem_mb

subroutine get_int(str, out_ival)

 character(len=*),intent(in) :: str
 integer,intent(out) :: out_ival
 istart = index(str, ":") + 1
 read(str(istart+1:), fmt=*, iostat=iostat, iomsg=pstat%iomsg) out_ival
 !ABI_CHECK(iostat == 0, pstat%iomsg)

end subroutine get_int

end subroutine pstat_from_file
!!***

!!****f* m_pstat/pstat_print
!! NAME
!!  pstat_print
!!
!! FUNCTION
!!  Print object in Yaml format to std_out
!!
!! SOURCE

subroutine pstat_print(pstat, comm, file, line)

 class(pstat_t),intent(inout) :: pstat
 character(len=*),optional,intent(in) :: file
 integer,optional,intent(in) :: line, comm

!Local variables-------------------------------
 integer :: units(1), ierr
 integer :: f90line = 0
 real(dp) :: min_mpicomm_vmrss_mb, max_mpicomm_vmrss_mb
 character(len=500) :: f90name='Subroutine Unknown'
 type(yamldoc_t) :: ydoc
! *************************************************************************

 if (pstat%pid == -1) return
 units(1) = std_out
 if (std_out < 1) return

 call pstat%from_file(pstat%filepath)

 if (present(line)) f90line = line
 if (present(file)) f90name = basename(file)

 min_mpicomm_vmrss_mb = pstat%vmrss_mb
 max_mpicomm_vmrss_mb = pstat%vmrss_mb
 if (present(comm)) then
   call xmpi_min(pstat%vmrss_mb, min_mpicomm_vmrss_mb, comm, ierr)
   call xmpi_max(pstat%vmrss_mb, max_mpicomm_vmrss_mb, comm, ierr)
 end if

#ifndef FC_NVHPC
 ydoc = yamldoc_open("PstatData")
 call ydoc%add_int("pid", pstat%pid)
 call ydoc%add_string("file", f90name)
 call ydoc%add_int("line", f90line)
 call ydoc%add_real("vmrss_mb", pstat%vmrss_mb)
 call ydoc%add_real("min_mpicomm_vmrss_mb", min_mpicomm_vmrss_mb)
 call ydoc%add_real("max_mpicomm_vmrss_mb", max_mpicomm_vmrss_mb)
 call ydoc%add_real("vmpeak_mb", pstat%vmpeak_mb)
 call ydoc%add_real("vmstk_mb", pstat%vmstk_mb)
 if (len_trim(pstat%iomsg) > 0) call ydoc%add_string("iomsg", trim(pstat%iomsg))
 call ydoc%write_units_and_free(units)
#else
 ! Yet another wild NVHPC bug (only on eos_nvhpc_23.9_elpa)
 write(std_out, "(a)")"--- !PstatData"
 write(std_out, *)"vmrss_mb: ", pstat%vmrss_mb
 write(std_out, *)"min_mpicomm_vmrss: ", min_mpicomm_vmrss_mb
 write(std_out, *)"max_mpicomm_vmrss: ", max_mpicomm_vmrss_mb
 write(std_out, "(a)")"..."
#endif

end subroutine pstat_print
!!***

!!****f* m_pstat/pstat_min_mem_mb_per_proc
!! NAME
!!  pstat_min_mem_mb_per_proc
!!
!! FUNCTION
!!  This function estimates the available memory (in MB) per process, within the MPI communicator comm
!!  based on process statistics (pstat).
!!
!! SOURCE

real(dp) function pstat_min_mem_mb_per_proc(pstat, comm) result(min_mem_mb)

!Arguments ------------------------------------
 class(pstat_t),intent(inout) :: pstat
 integer,intent(in) :: comm

!Local variables-------------------------------
 integer :: ierr
 logical :: all_ok
! *************************************************************************

 call pstat%from_file(pstat%filepath)

 all_ok = pstat%ok
 call xmpi_land(all_ok, comm)

 if (.not. all_ok) then
   ! Handle case in which pstat is not available or something went wrong when reading.
   min_mem_mb = mem_per_cpu_mb * half
   return
 end if

 ! Compute min inside comm
 call xmpi_min(pstat%vmrss_mb, min_mem_mb, comm, ierr)

 min_mem_mb = (mem_per_cpu_mb - min_mem_mb)

 ! Fallback for too small values
 if (min_mem_mb <= tol1) min_mem_mb = mem_per_cpu_mb * half

end function pstat_min_mem_mb_per_proc
!!***

end module m_pstat
!!***

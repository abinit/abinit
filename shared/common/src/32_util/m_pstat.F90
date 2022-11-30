!!****m* ABINIT/m_pstat
!! NAME
!! m_pstat
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2017-2022 ABINIT group (MG)
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

 use defs_basis
 use m_abicore
 use m_errors
 use iso_c_binding
 use m_yaml

 use m_fstrings, only : find_and_select
 use m_clib, only : clib_getpid

 implicit none

 private

!----------------------------------------------------------------------

!!****t* m_pstat/pstat_t
!! NAME
!! pstat_t
!!
!! FUNCTION
!! This object stores the most important quantites reported in the /proc/{pid}/status file
!! in particular the virtual memory VmRSS that can be used at runtime to define block sizes
!! See https://docs.kernel.org/filesystems/proc.html
!!
!! NB: This file is only available on Linux hence one should always check the value of pstat%ok
!! before using quantities such as vmrss_mb.
!!
!! SOURCE

 type, public :: pstat_t

  logical :: ok = .False.

  integer :: pid = -1

  integer :: threads = -1
  ! number of threads

  integer :: fdsize = -1
   ! number of file descriptor slots currently allocated

  real(dp) :: vmrss_mb = -one
  ! size of memory portions. It contains the three following parts (VmRSS = RssAnon + RssFile + RssShmem)

  real(dp) :: vmpeak_mb = -one
  ! peak virtual memory size

  real(dp) :: vmstk_mb = -one
  ! size of stack segments

 contains
   procedure :: from_pid => pstat_from_pid     ! Init object from process identifier (main entry point).
   procedure :: from_file => pstat_from_file   ! Init object from file (useful for debugging).
   procedure :: print => pstat_print           ! Print object.
 end type pstat_t

!----------------------------------------------------------------------

contains
!!***

!!****f* m_gwr/pstat_from_pid
!! NAME
!!  pstat_from_pid
!!
!! FUNCTION
!!   Init object from process identifier (main entry point).
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






 call pstat%from_file('/proc/'//trim(spid)//'/status')

end subroutine pstat_from_pid
!!***

!!****f* m_gwr/pstat_from_file
!! NAME
!!  pstat_from_file
!!
!! FUNCTION
!! Init object from file (useful for debugging).
!!
!! SOURCE

subroutine pstat_from_file(pstat, filepath)

!Arguments ------------------------------------
 class(pstat_t),intent(out) :: pstat
 character(len=*),intent(in) :: filepath

!Local variables-------------------------------
 integer :: unit, ierr
 character(len=500) :: line, spid, iomsg
 integer :: istart, istop, iostat
! *************************************************************************

 !inquire(file="/proc/"//trim(spid)//'/status', exist=exist)
 !if (.not. exist) then
 !  ierr = 1; return
 !end if

 pstat%ok = .False.
 open(newunit=unit, file=filepath, action="read", status='old', iostat=ierr)
 if (ierr /= 0) return

 do
   read(unit, "(a)", iostat=ierr, end=10) line
   if (ierr > 0) then
     close(unit)
     return
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

10  close(unit)
  pstat%ok = .True.

contains

subroutine get_mem_mb(str, mem_mb)

  character(len=*),intent(in) :: str
  real(dp),intent(out) :: mem_mb

  real(dp) :: mem_fact

  ! Generic mem entry has format:
  !VmRSS: 2492 kB
  istart = index(str, ":") + 1
  istop = find_and_select(str, &
                         ["kB", "mB"], &
                         [one/1024._dp, one], mem_fact, iomsg) !default=one,
  ABI_CHECK(istop /= -1, iomsg)
  read(str(istart+1:istop-1), fmt=*, iostat=iostat, iomsg=iomsg) mem_mb
  ABI_CHECK(iostat == 0, iomsg)
  mem_mb = mem_mb * mem_fact
end subroutine get_mem_mb

subroutine get_int(str, out_ival)
  character(len=*),intent(in) :: str
  integer,intent(out) :: out_ival
  istart = index(str, ":") + 1
  read(str(istart+1:), fmt=*, iostat=iostat, iomsg=iomsg) out_ival
  ABI_CHECK(iostat == 0, iomsg)
end subroutine get_int

end subroutine pstat_from_file
!!***

!!****f* m_gwr/pstat_print
!! NAME
!!  pstat_print
!!
!! FUNCTION
!!  Print object in Yaml format.
!!
!! SOURCE

subroutine pstat_print(pstat, units, header)

 class(pstat_t),intent(in) :: pstat
 integer,intent(in) :: units(:)
 character(len=*),optional,intent(in) :: header

!Local variables-------------------------------
 character(len=500) :: header__
 type(yamldoc_t) :: ydoc
! *************************************************************************

 header__ = "unknown"; if (present(header)) header__ = header
 ydoc = yamldoc_open(header__) !, width=11, real_fmt='(3f8.3)')

 call ydoc%add_int("pid", pstat%pid)
 call ydoc%add_int("threads", pstat%threads)
 call ydoc%add_int("fdsize", pstat%fdsize)
 call ydoc%add_real("vmrss_mb", pstat%vmrss_mb)
 call ydoc%add_real("vmpeak_mb", pstat%vmpeak_mb)
 call ydoc%add_real("vmstk_mb", pstat%vmstk_mb)

 call ydoc%write_units_and_free(units)

end subroutine pstat_print
!!***

!!****f* m_gwr/pstat_gather
!! NAME
!!  pstat_gather
!!
!! FUNCTION
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine pstat_gather(pstat, vmrss_mb, comm)

 class(pstat_t),intent(out) :: pstat
 real(dp),intent(out) :: vmrss_mb
 integer,intent(in) :: comm
 !integer,intent(out) :: ierr

 integer :: ierr, int_list(5)
 real(dp) :: real_list(3)

 call pstat%from_pid()
 real_list = [pstat%vmrss_mb, pstat%vmpeak_mb, pstat%vmstk_mb]
 !call xmpi_max_ip(real_list, comm, ierr)

end subroutine pstat_gather
!!***

end module m_pstat
!!***

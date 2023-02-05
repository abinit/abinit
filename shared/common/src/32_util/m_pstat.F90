!!****m* ABINIT/m_pstat
!! NAME
!! m_pstat
!!
!! FUNCTION
!! Interface to the /proc/{pid}/status file available on Linux.
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
!! Usage:
!!
!!   call ps%from_pid()
!!   call ps%print([std_out])
!!   ! reload data and print it
!!   call ps%update()
!!   call ps%print([std_out])
!!
!! SOURCE

 type, public :: pstat_t

  logical :: ok = .False.
  ! False if stat file is not available

  integer :: pid = -1
  ! Process identifier.

  integer :: threads = -1
  ! number of threads.

  integer :: fdsize = -1
   ! Number of file descriptor slots currently allocated

  real(dp) :: vmrss_mb = -one
  ! Size of memory portions.
  ! It contains the three following parts (VmRSS = RssAnon + RssFile + RssShmem)

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
 end type pstat_t
!!***

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

subroutine pstat_from_pid(ps)

!Arguments ------------------------------------
 class(pstat_t),intent(out) :: ps

!Local variables-------------------------------
 integer(c_int) :: pid
 character(len=500) :: spid
! *************************************************************************

 pid = clib_getpid()
 write(spid, "(i0)") pid
 spid = adjustl(spid)
 call ps%from_file('/proc/'//trim(spid)//'/status')

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

subroutine pstat_from_file(ps, filepath)

!Arguments ------------------------------------
 class(pstat_t),intent(inout) :: ps
 character(len=*),intent(in) :: filepath

!Local variables-------------------------------
 integer :: unit, ierr
 character(len=500) :: line
 integer :: istart, istop, iostat
! *************************************************************************

 ps%ok = .False.
 ps%filepath = filepath

 open(newunit=unit, file=trim(ps%filepath), action="read", status="old", iostat=ierr, iomsg=ps%iomsg)
 if (ierr /= 0) then
   close(unit); return
 end if

 do
   read(unit, "(a)", iostat=ierr, end=10, iomsg=ps%iomsg) line
   if (ierr > 0) then ! EOF
     close(unit); return
   end if

   ! Parse useful integers
   if (index(line, "Pid:") == 1) call get_int(line, ps%pid)
   if (index(line, "Threads:") == 1) call get_int(line, ps%threads)
   if (index(line, "FDSize:") == 1) call get_int(line, ps%fdsize)

   ! Parse memory entries
   if (index(line, "VmRSS:") == 1) call get_mem_mb(line, ps%vmrss_mb)
   if (index(line, "VmPeak:") == 1) call get_mem_mb(line, ps%vmpeak_mb)
   if (index(line, "VmStk:") == 1) call get_mem_mb(line, ps%vmstk_mb)
 end do

10 close(unit)
  ps%ok = .True.

contains

subroutine get_mem_mb(str, mem_mb)
 character(len=*),intent(in) :: str
 real(dp),intent(out) :: mem_mb

 ! Generic mem entry has format:
 !VmRSS: 2492 kB
 real(dp) :: mem_fact
 istart = index(str, ":") + 1
 istop = find_and_select(str, &
                        ["kB", "mB"], &
                        [one/1024._dp, one], mem_fact, ps%iomsg) !default=one,
 ABI_CHECK(istop /= -1, ps%iomsg)
 read(str(istart+1:istop-1), fmt=*, iostat=iostat, iomsg=ps%iomsg) mem_mb
 ABI_CHECK(iostat == 0, ps%iomsg)
 mem_mb = mem_mb * mem_fact
end subroutine get_mem_mb

subroutine get_int(str, out_ival)
 character(len=*),intent(in) :: str
 integer,intent(out) :: out_ival
 istart = index(str, ":") + 1
 read(str(istart+1:), fmt=*, iostat=iostat, iomsg=ps%iomsg) out_ival
 ABI_CHECK(iostat == 0, ps%iomsg)
end subroutine get_int

end subroutine pstat_from_file
!!***

!!****f* m_pstat/pstat_print
!! NAME
!!  pstat_print
!!
!! FUNCTION
!!  Print object in Yaml format.
!!
!! SOURCE

subroutine pstat_print(ps, units, header, reload)

 class(pstat_t),intent(inout) :: ps
 integer,intent(in) :: units(:)
 character(len=*),optional,intent(in) :: header
 logical,optional,intent(in) :: reload

!Local variables-------------------------------
 character(len=500) :: header__
 type(yamldoc_t) :: ydoc
! *************************************************************************

 if (present(reload)) then
   if (reload) call ps%from_file(ps%filepath)
 end if

 header__ = "unknown"; if (present(header)) header__ = header
 ydoc = yamldoc_open(header__) !, width=11, real_fmt='(3f8.3)')

 call ydoc%add_int("pid", ps%pid)
 call ydoc%add_int("threads", ps%threads) !, comment="")
 call ydoc%add_int("fdsize", ps%fdsize) !, comment="")
 call ydoc%add_real("vmrss_mb", ps%vmrss_mb) !, comment="")
 call ydoc%add_real("vmpeak_mb", ps%vmpeak_mb) !, comment="")
 call ydoc%add_real("vmstk_mb", ps%vmstk_mb) !, comment="")
 if (len_trim(ps%iomsg) > 0) call ydoc%add_string("iomsg", trim(ps%iomsg))

 call ydoc%write_units_and_free(units)

end subroutine pstat_print
!!***

!!****f* m_pstat/pstat_gather
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

#if 0

subroutine pstat_gather(ps, vmrss_mb, comm)

 class(pstat_t),intent(inout) :: ps
 real(dp),intent(out) :: vmrss_mb
 integer,intent(in) :: comm

 !integer :: ierr, int_list(5)
 real(dp) :: real_list(3)

 call ps%from_pid()
 real_list = [ps%vmrss_mb, ps%vmpeak_mb, ps%vmstk_mb]
 !call xmpi_max_ip(real_list, comm, ierr)

 ps%vmrss_mb = real_list(1)
 ps%vmpeak_mb = real_list(2)
 ps%vmstk_mb = real_list(3)

end subroutine pstat_gather
!!***
#endif

end module m_pstat
!!***

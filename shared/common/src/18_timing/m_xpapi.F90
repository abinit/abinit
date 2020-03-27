!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xpapi
!! NAME
!! m_xpapi
!!
!! FUNCTION
!!  Thin wrapper for the PAPI library.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (MG,DC)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_xpapi

 use defs_basis
 use m_errors
 use iso_c_binding

 implicit none

 private

#ifdef HAVE_PAPI
#include "f90papi.h"
#endif

 public :: xpapi_init
 public :: xpapi_show_info
 public :: xpapi_flops
 public :: xpapi_shutdown
!!***

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_xpapi/xpapi_init
!! NAME
!!  xpapi_init
!!
!! FUNCTION
!!  initialize the PAPI library. It must be called before any low level PAPI functions can be used.
!!  If your application is making use of threads PAPI_thread_init (3) also be called prior to making
!!  any calls to the library other than PAPI_library_init().
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      papif_is_initialized,papif_perror
!!
!! SOURCE

subroutine xpapi_init()

#ifdef HAVE_PAPI
!Local variables-------------------------------
 character(len=PAPI_MAX_STR_LEN) :: papi_errstr
 integer(C_INT) :: check
 real(C_FLOAT) :: real_time,proc_time,mflops
 integer(C_LONG_LONG) :: flpops

! *************************************************************************

 check = PAPI_VER_CURRENT
 call PAPIf_library_init(check)

 if ( check /= PAPI_VER_CURRENT .and. check >0 ) then
   MSG_WARNING(" PAPI library version mismatch!")
 end if
 !ABI_CHECK(check>=0," PAPI Initialization error!")
 !XPAPI_CHECK(check," PAPI Initialization error!")

!#ifdef HAVE_OPENMP
 !call PAPIf_thread_init(pthread_self, check)
 !XPAPI_CHECK(check>=0," PAPIf_thread_init")
!#endif

! First pass. Set up counters to monitor PAPI_FP_OPS and PAPI_TOT_CYC events and start the counters
! Subsequent calls will read the counters and return total real time,
! total process time, total floting point instructions or operations
! since the start of the mesurement and the Mflop/s rate since latests call to PAPI_flops
 call PAPIf_flops(real_time, proc_time, flpops, mflops, check)
 XPAPI_CHECK(check,"Problem in PAPIf_flops")

#endif

end subroutine xpapi_init
!!***

!----------------------------------------------------------------------

!!****f* m_xpapi/xpapi_show_info
!! NAME
!!  xpapi_show_info
!!
!! FUNCTION
!!
!! INPUTS
!!  [unit]=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!         Default = std_out
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine: DEFAULT = "COLL"
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      papif_is_initialized,papif_perror
!!
!! SOURCE

subroutine xpapi_show_info(unit,mode_paral)

!Arguments-------------------------
 integer,optional,intent(in) :: unit
 character(len=*),optional,intent(in) :: mode_paral

!Local variables-------------------
#ifdef HAVE_PAPI
 integer :: unt
 integer(C_INT) :: num_hwcntrs,ncpu,nnodes,totalcpus,vendor,model
 real(C_FLOAT) :: revision,mhz
 character(len=PAPI_MAX_STR_LEN) :: vendor_string,model_string
 character(len=500) :: msg,my_mode
#endif

! *************************************************************************

#ifdef HAVE_PAPI
 unt = std_out;    if (PRESENT(unit)) unt = unit
 my_mode = "COLL"; if (PRESENT(mode_paral)) my_mode = mode_paral

 call PAPIf_num_counters(num_hwcntrs)
 if (num_hwcntrs  < 0) then
   MSG_WARNING(" The installation does not support PAPI")
 end if

 if (num_hwcntrs == 0) then
   msg = " The installation supports PAPI, but this machine does not provide hardware counters."
   MSG_WARNING(msg)
 end if

 call PAPIF_get_hardware_info (ncpu,nnodes,totalcpus,vendor,vendor_string,model,model_string,revision,mhz)

 write(msg,"(a,i0)")" PAPI Version ",PAPI_VER_CURRENT
 call wrtout(unt,msg,my_mode)
 write(msg,"(a,i0)")" Number of hardware counters: ",num_hwcntrs
 call wrtout(unt,msg,my_mode)
 write(msg,"(a,i0)")" Number of CPUs in an SMP Node: ",ncpu
 call wrtout(unt,msg,my_mode)
 write(msg,"(a,i0)")" Number of nodes in the entire system: ",nnodes
 call wrtout(unt,msg,my_mode)
 write(msg,"(a,i0)")" Total number of CPUs in the entire system: ",totalcpus
 call wrtout(unt,msg,my_mode)
 !write(msg,"(a,i0)")" Vendor id number of CPU: ",vendor
 !call wrtout(unt,msg,my_mode)
 !write(msg,"(a)")   " Vendor id string of CPU: "//TRIM(vendor_string)
 !call wrtout(unt,msg,my_mode)
 !write(msg,"(a,i0)")" Model number of CPU: ",model
 !call wrtout(unt,msg,my_mode)
 write(msg,"(a)")   " Model string of CPU: "//TRIM(model_string)
 !call wrtout(unt,msg,my_mode)
 !write(msg,"(a,f5.1)")" Revision number of CPU: ",revision
 !write(msg,"(a,f8.2")" Cycle time of this CPU; *may* be an estimate generated at init time with a quick timing routine",mhz

#else
 ABI_UNUSED(mode_paral)
 ABI_UNUSED(unit)
#endif

end subroutine xpapi_show_info
!!***

!----------------------------------------------------------------------

!!****f* m_xpapi/xpapi_flops
!! NAME
!!  xpapi_flops
!!
!! FUNCTION
!!   PAPI High Level: Simplified call to get Mflops/s, real and processor time
!!
!! OUTPUT
!!  real_time -- total realtime since the first PAPI_flops() call
!!  proc_time -- total process time since the first PAPI_flops() call
!!  flops -- total floating point instructions or operations since the first call
!!  mflops -- Mflop/s achieved since the previous call
!!  check = exit status
!!
!! PARENTS
!!      m_time
!!
!! CHILDREN
!!      papif_is_initialized,papif_perror
!!
!! SOURCE

subroutine xpapi_flops(real_time,proc_time,flops,mflops,check)

!Arguments-------------------------
 integer(C_INT),intent(out) :: check
 integer(C_LONG_LONG),intent(out) :: flops
 real(C_FLOAT),intent(out) :: real_time,proc_time,mflops

! *************************************************************************

#ifdef HAVE_PAPI
 call PAPIf_flops(real_time, proc_time, flops, mflops, check)
#endif

end subroutine xpapi_flops
!!***

!----------------------------------------------------------------------

!!****f* m_xpapi/xpapi_shutdown
!! NAME
!!  xpapi_shutdown
!!
!! FUNCTION
!!   exit function used by the PAPI Library to free resources and shut down when certain error conditions arise.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      papif_is_initialized,papif_perror
!!
!! SOURCE

subroutine xpapi_shutdown()

! *************************************************************************

#ifdef HAVE_PAPI
 call PAPIf_shutdown()
#endif

end subroutine xpapi_shutdown
!!***

!----------------------------------------------------------------------

!!****f* m_xpapi/xpapi_handle_error
!! NAME
!!  xpapi_handle_error
!!
!! FUNCTION
!!
!! PARENTS
!!
!! CHILDREN
!!      papif_is_initialized,papif_perror
!!
!! SOURCE

subroutine xpapi_handle_error(check,err_msg,file,line)

!Arguments ------------------------------------
!scalars
 integer(C_INT),intent(in) :: check
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: err_msg
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line
 character(len=10) :: lnum
 character(len=500) :: f90name
 character(len=500) :: my_msg
#ifdef HAVE_PAPI
 integer(C_INT) :: ierr,ans
 character(len=PAPI_MAX_STR_LEN) :: papi_errstr
#endif

! *************************************************************************

 if (PRESENT(line)) then
   f90line=line
 else
   f90line=0
 end if
 write(lnum,'(i0)')f90line

 if (PRESENT(file)) then
   !f90name = basename(file)
   f90name = file
 else
   f90name='Subroutine Unknown'
 end if

 my_msg=TRIM(f90name)//":"//TRIM(lnum)//":"

#ifdef HAVE_PAPI
 if (check /= PAPI_OK) then
   write(std_out,*) " Error in papi library at: "//TRIM(my_msg)
   write(std_out,*) " User message: "//TRIM(err_msg)
   call PAPIF_is_initialized(ans)
   if (ans == PAPI_LOW_LEVEL_INITED) then
     call papif_perror(check,papi_errstr,ierr)
     write(std_out,*) 'Error code: ',TRIM(papi_errstr)
   else
     write(std_out,*) "Papi library is not initialized!"
   end if
 end if
 MSG_ERROR("Fatal error")
#else
 ABI_UNUSED(err_msg)
 if (.FALSE.) write(std_out,*) check
#endif

end subroutine xpapi_handle_error
!!***

END MODULE m_xpapi
!!***

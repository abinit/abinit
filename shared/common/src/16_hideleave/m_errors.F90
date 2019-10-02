!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_errors
!! NAME
!!  m_errors
!!
!! FUNCTION
!!  This module contains low-level procedures to check assertions and handle errors.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2019 ABINIT group (MG,YP,NCJ,MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_errors

 use defs_basis
 use m_profiling_abi
 use m_xmpi
 use m_specialmsg, only : wrtout
#ifdef HAVE_NETCDF
 use netcdf
#endif
#ifdef HAVE_ETSF_IO
 use etsf_io_low_level
 use etsf_io
#endif
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef FC_NAG
 use f90_unix_proc
#endif
#ifdef FC_INTEL
 use ifcore
#endif

 use m_io_tools,        only : flush_unit, lock_and_write, file_exists, num_opened_units, show_units, open_file
 use m_fstrings,        only : toupper, basename, indent, lstrip, atoi, strcat, itoa
 use m_build_info,      only : dump_config, abinit_version
 use m_cppopts_dumper,  only : dump_cpp_options
 use m_optim_dumper,    only : dump_optim

 implicit none

#if defined HAVE_MPI1
include 'mpif.h'
#endif

#ifdef FC_IBM
include "fexcp.h"
#endif

 private
!!***

!Public procedures
 public :: assert_eq        ! Report and die gracefully if integers not all equal (used for size checking).
 public :: assert           ! Report and die if any logical is false (used for argument range checking).
 public :: sentinel         ! Announce the entering or the exiting from a procedure.
 public :: die              ! Stop execution in case of unexpected events.
 public :: msg_hndl         ! Basic Error handlers.
 public :: netcdf_check     ! Stop execution after a NetCDF I/O error
 public :: check_mpi_ierr   ! Error handler for MPI routines.
 public :: set_backtrace_onerr ! Activate show_backtrace call in msg_hndl. 0 to disable it.
 !public :: show_backtrace   ! Shows a backtrace at an arbitrary place in user code. (Gfortran/Ifort extension)
 public :: unused_var       ! Helper function used to silence compiler warnings due to unused variables.
#if defined HAVE_ETSF_IO
 public :: abietsf_msg_hndl ! Error handler for ETSF-IO routines.
 public :: abietsf_warn     ! Write warnings reported by ETSF-IO routines.
#endif
 public :: bigdft_lib_error
 public :: xlf_set_sighandler
 public :: abinit_doctor         ! Perform checks on memory leaks and leaking file descriptors
                                 ! at the end of the run.
 public :: abi_abort             ! Abort the code
 public :: abi_cabort            ! C-interoperable version.

 ! This flag activate the output of the backtrace in msg_hndl
 ! Unfortunately, gcc4.9 seems to crash inside this routine
 ! hence, for the time being, this optional feature has been disabled
 integer, save, private :: m_errors_show_backtrace = 0

 interface assert_eq
   module procedure assert_eq2
   module procedure assert_eq3
   module procedure assert_eq4
   module procedure assert_eqn
 end interface assert_eq

 interface assert
   module procedure assert1
   module procedure assert2
   module procedure assert3
   module procedure assert4
   module procedure assert_v
 end interface assert

 interface unused_var
   module procedure unused_int
   module procedure unused_real_dp
   module procedure unused_real_sp
   module procedure unused_cplx_dpc
   module procedure unused_cplx_spc
   module procedure unused_logical
   module procedure unused_ch
 end interface unused_var

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq2
!! NAME
!!  assert_eq2
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! SOURCE

function assert_eq2(l1,l2,message,file,line)

!Arguments ------------------------------------
 integer,intent(in) :: l1,l2
 integer,optional,intent(in) :: line
 integer :: assert_eq2
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'

! *************************************************************************

 if (l1==l2) then
  assert_eq2=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,line)
 end if

end function assert_eq2
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq3
!! NAME
!!  assert_eq3
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! SOURCE

function assert_eq3(l1,l2,l3,message,file,line)

!Arguments ------------------------------------
 integer,intent(in) :: l1,l2,l3
 integer,optional,intent(in) :: line
 integer :: assert_eq3
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (l1==l2.and.l2==l3) then
  assert_eq3=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,line)
 end if

end function assert_eq3
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eq4
!! NAME
!!  assert_eq4
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! INPUTS
!!  l1,l2,.. Integers to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! SOURCE

function assert_eq4(l1,l2,l3,l4,message,file,line)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: l1,l2,l3,l4
 integer,optional,intent(in) :: line
 integer :: assert_eq4
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (l1==l2.and.l2==l3.and.l3==l4) then
  assert_eq4=l1
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,line)
 end if

end function assert_eq4
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_eqn
!! NAME
!!  assert_eqn
!!
!! FUNCTION
!!  Report and die gracefully if integers not all equal (used for size checking).
!!
!! SOURCE

function assert_eqn(nn,message,file,line)

!Arguments ------------------------------------
!scalars
 integer,optional,intent(in) :: line
 integer :: assert_eqn
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
!arrays
 integer,intent(in) :: nn(:)

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (ALL(nn(2:)==nn(1))) then
  assert_eqn=nn(1)
 else
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,line)
 end if

end function assert_eqn
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert1
!! NAME
!!  assert1
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if
!!  any logical is false (used for arg range checking).
!!
!! INPUTS
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additiona information
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine assert1(l1,message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.l1) then
   if (PRESENT(line)) f90line=line
   if (PRESENT(file)) f90name= basename(file)
   call msg_hndl(message,'ERROR','PERS',f90name,f90line)
 end if

end subroutine assert1
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert2
!! NAME
!!  assert2
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if
!   any logical is false (used for arg range checking).
!!
!! INPUTS
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine assert2(l1,l2,message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,f90line)
 end if

end subroutine assert2
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert3
!! NAME
!!  assert3
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if
!!  any logical is false (used for arg range checking).
!!
!! INPUTS
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine assert3(l1,l2,l3,message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,f90line)
 end if

end subroutine assert3
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert4
!! NAME
!!  assert4
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if
!!  any logical is false (used for arg range checking).
!!
!! INPUTS
!!  l1,l2,.. logical values to be checked (array version is also provided)
!!  message(len=*)=tag with additional information
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine assert4(l1,l2,l3,l4,message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: l1,l2,l3,l4

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.(l1.and.l2.and.l3.and.l4)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,f90line)
 end if

end subroutine assert4
!!***

!----------------------------------------------------------------------

!!****f* m_errors/assert_v
!! NAME
!!  assert_v
!!
!! FUNCTION
!!  Routines for argument checking and error handling. Report and die if
!!  any logical is false (used for arg range checking).
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine assert_v(n,message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file
 logical,intent(in) :: n(:)

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
! *************************************************************************

 if (.not.ALL(n)) then
  if (PRESENT(line)) f90line=line
  if (PRESENT(file)) f90name= basename(file)
  call msg_hndl(message,'ERROR','PERS',f90name,f90line)
 end if

end subroutine assert_v
!!***

!----------------------------------------------------------------------

!!****f* m_errors/netcdf_check
!! NAME
!!  netcdf_check
!!
!! FUNCTION
!!  Error handler for Netcdf calls.
!!
!! INPUTS
!!  ncerr=Status error returned by the Netcdf library.
!!  msg=User-defined string with info on the action that was performed
!!  file= name of the file.
!!  line= line number.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine netcdf_check(ncerr,msg,file,line)

!Arguments ------------------------------------
 integer,intent(in) :: ncerr
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: file
 integer,optional,intent(in) :: line

!Local variables-------------------------------
 integer :: f90line
 character(len=500) :: f90name
 character(len=1024) :: nc_msg
 character(len=2048) :: my_msg

! *************************************************************************

#ifdef HAVE_NETCDF
 if (ncerr /= NF90_NOERR) then
   if (PRESENT(line)) then
     f90line=line
   else
     f90line=0
   end if
   if (PRESENT(file)) then
     f90name = basename(file)
   else
     f90name='Subroutine Unknown'
   end if
   !
   ! Append Netcdf string to user-defined message.
   write(nc_msg,'(a,2x,a)')' - NetCDF library returned:',TRIM(nf90_strerror(ncerr))
   !write(std_out,*)TRIM(nf90_strerror(ncerr))
   my_msg = TRIM(msg) // TRIM(nc_msg)

   call msg_hndl(my_msg,"ERROR","PERS",f90name,f90line)
 end if
#endif

end subroutine netcdf_check
!!***

!----------------------------------------------------------------------

!!****f* m_errors/sentinel
!! NAME
!!  sentinel
!!
!! FUNCTION
!!  Announce the entering and the exiting from a function. Useful for poor-man debugging.
!!
!! INPUTS
!!  level=1 when entering, 2 for exit.
!!  mode_paral= ['COLL'|'PERS'|'COLL_SILENT|PERS_SILENT']
!!   'COLL' and 'PERS' refer to the output mode used in wrtout to report the message.
!!   'COLL_SILENT' and 'PERS_SILENT' can be used if the procedure is called several times inside a loop.
!!   In this case sentinel will report only the first entry and the first exit using either 'COLL' or 'PERS' mode.
!!  file=File name
!!  func=Name of the procedure to be tested (passed through ABI_FUNC macro)
!!  [line]=Line number. Defaults to 0.
!!
!! NOTES
!!  This routine is usually interfaced with the macros defined in abi_common.h
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine sentinel(level,mode_paral,file,func,line)

!Arguments ------------------------------------
 integer,intent(in) :: level
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: mode_paral
 character(len=*),optional,intent(in) :: func
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer,save :: level_save=0
 integer :: ii
 integer :: f90line
 character(len=500),save :: func_save
 character(len=4) :: my_mode
 character(len=10) :: lnum
 character(len=500) :: my_func, my_file
 character(len=500) :: msg

! *********************************************************************

 ! initialize the variable
 my_func = 'Function Unknown'; if (PRESENT(func)) my_func = basename(func)
 my_file = "File Unknown"; if (PRESENT(file)) my_file = basename(file)

 level_save=level; func_save=my_func

 f90line=0; if (PRESENT(line)) f90line=line

 if (toupper(mode_paral)=='COLL_SILENT'.or.toupper(mode_paral)=='PERS_SILENT') then
    ! * Silent mode, check if we are inside a loop.
    if (level==level_save .and. my_func==func_save) RETURN
    ii = index( toupper(mode_paral), '_SILENT')
    my_mode=toupper(mode_paral(1:ii-1))
 else ! * Normal mode.
    my_mode=mode_paral
 end if

 if (my_mode/='COLL'.or.my_mode/='PERS') my_mode='COLL'

 write(lnum,"(i0)")f90line
 my_func= TRIM(my_func)//"@"//TRIM(my_file)//":"//TRIM(lnum)

 if (level==1) then
    msg = ' '//TRIM(my_func)//' >>>>> ENTER'//ch10
 else if (level==2) then
    msg = ' '//TRIM(my_func)//' >>>>> EXIT '//ch10
 else
    call die('Wrong level', &
&   __FILE__,&
&   __LINE__)
 end if

 call wrtout(std_out,msg,my_mode)
 call flush_unit(std_out)

end subroutine sentinel
!!***

!----------------------------------------------------------------------

!!****f* m_errors/die
!! NAME
!!  die
!!
!! FUNCTION
!!  Stop smoothly the execution in case of unexpected events reporting the
!!  line number and the file name where the error occurred as well as the
!!  MPI rank of the processor. This routine is usually interfaced through
!!  some macro defined in abi_common.h
!!
!! INPUTS
!!  message=String containing additional information on the nature of the problem
!!  line=Line number of the file where problem occurred
!!  f90name=Name of the f90 file containing the caller
!!
!! PARENTS
!!      m_errors,m_xc_vdw
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine die(message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer :: rank
 integer :: f90line=0
 character(len=10) :: lnum,strank
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=500) :: msg

! *********************************************************************

 if (PRESENT(line)) f90line=line
 write(lnum,"(i0)")f90line

 ! === Determine my rank inside MPI_COMM_WORLD ===
 rank = xmpi_comm_rank(xmpi_world)
 write(strank,"(i0)")rank

 if (PRESENT(file)) f90name= basename(file)
 msg=TRIM(f90name)//':'//TRIM(lnum)//' P'//TRIM(strank)

 write(msg,'(a,2x,2a,2x,a)')ch10,TRIM(msg),ch10,TRIM(message)

 call wrtout(std_out,msg,'PERS')
 !if is_connected(ab_out)) call wrtout(ab_out,msg,'PERS')
 call abi_abort('PERS')

end subroutine die
!!***

!----------------------------------------------------------------------

!!****f* m_errors/msg_hndl
!! NAME
!!  msg_hndl
!!
!! FUNCTION
!!  Basic error handler for abinit. This routine is usually interfaced through some macro defined in abi_common.h
!!
!! INPUTS
!!  message=string containing additional information on the nature of the problem
!!  level=string defining the type of problem. Possible values are
!!   COMMENT
!!   WARNING
!!   ERROR
!!   STOP
!!   BUG
!!  mode_paral=Either "COLL" or "PERS".
!!  [line] = line number of the file where problem occurred
!!  [file] = name of the f90 file containing the caller
!!  [NODUMP]= if present dump config before stopping
!!  [NOSTOP]= if present don't stop even in the case of an error or a bug
!!  [unit]= Unit number (defaults to std_out)
!!
!! OUTPUT
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine msg_hndl(message,level,mode_paral,file,line,NODUMP,NOSTOP,unit)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line, unit
 logical,optional,intent(in) :: NODUMP,NOSTOP
 character(len=*),intent(in) :: level,message
 character(len=*),optional,intent(in) :: file
 character(len=*),intent(in) :: mode_paral

!Local variables-------------------------------
 integer :: f90line,ierr,unit_
 character(len=10) :: lnum
 character(len=500) :: f90name
 character(len=LEN(message)) :: my_msg
 character(len=MAX(4*LEN(message),2000)) :: sbuf ! Increase size and keep fingers crossed!

! *********************************************************************
 unit_ = std_out; if (present(unit)) unit_ = unit

 if (PRESENT(line)) then
   f90line=line
 else
   f90line=0
 end if
 ! TODO: fldiff.py should ignore f90line when comparing files (we don't want to
 ! update ref files if a new line is added to F90 source file!
 if (unit_ == ab_out) f90line = 0
 write(lnum,"(i0)")f90line

 if (PRESENT(file)) then
   f90name = basename(file)
 else
   f90name='Subroutine Unknown'
 end if

 my_msg = lstrip(message)

 select case (toupper(level))

 case ('COMMENT','WARNING')

   write(sbuf,'(8a,i0,7a)')ch10,&
     "--- !",TRIM(level),ch10,&
     "src_file: ",TRIM(f90name),ch10,&
     "src_line: ",f90line,ch10,&
     "message: |",ch10,TRIM(indent(my_msg)),ch10,&
     "...",ch10
   call wrtout(unit_, sbuf, mode_paral)

 case ('STOP')

   write(sbuf,'(8a)')ch10,&
     "--- !",TRIM(level),ch10,&
     "message: |",ch10,TRIM(indent(my_msg)),ch10
   call wrtout(unit_, sbuf, mode_paral)
   if (.not.present(NOSTOP)) then
     call abi_abort(mode_paral,print_config=.FALSE.)
   end if

 ! ERROR' or 'BUG'
 case default

   if ((.not.present(NOSTOP)).and.(.not.present(NODUMP))) then
     !call print_kinds()
     !call xmpi_show_info()
     !call dump_config(std_out)
     ! Dump the backtrace if the compiler supports it.
     if (m_errors_show_backtrace == 1) call show_backtrace()
   end if

   write(sbuf,'(8a,i0,2a,i0,7a)')ch10,&
     "--- !",TRIM(level),ch10,&
     "src_file: ",TRIM(f90name),ch10,&
     "src_line: ",f90line,ch10,&
     "mpi_rank: ",xmpi_comm_rank(xmpi_world),ch10,&
     "message: |",ch10,TRIM(indent(my_msg)),ch10,&
     "...",ch10
   call wrtout(unit_, sbuf, mode_paral)

   if (.not.present(NOSTOP)) then
     ! The first MPI proc that gets here, writes the ABI_MPIABORTFILE with the message!
     ! The file is written only if nprocs > 1. Do not change this behaviour!
     if (.not. file_exists(ABI_MPIABORTFILE) .and. xmpi_comm_size(xmpi_world) > 1) then
        call lock_and_write(ABI_MPIABORTFILE, sbuf, ierr)
     end if
     ! And now we die!
     call abi_abort(mode_paral,print_config=.FALSE.)
   end if

 end select

end subroutine msg_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_errors/set_backtrace_onerr
!! NAME
!! set_backtrace_onerr
!!
!! FUNCTION
!!  1 to activate show_backtrace call in msg_hndl. 0 to disable it
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine set_backtrace_onerr(iflag)

!Arguments ------------------------------------
 integer,intent(in) :: iflag
! *********************************************************************

  m_errors_show_backtrace = iflag

end subroutine set_backtrace_onerr
!!***

!----------------------------------------------------------------------

!!****f* m_errors/show_backtrace
!! NAME
!! show_backtrace
!!
!! FUNCTION
!!  shows a backtrace at an arbitrary place in user code.
!!  Program execution continues normally afterwards.
!!  The backtrace information is printed to the unit corresponding to ERROR_UNIT in ISO_FORTRAN_ENV.
!!  This is a (Gfortran extension| Ifort Extension)
!!
!! PARENTS
!!      m_errors
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine show_backtrace()


#if defined FC_GNU && defined HAVE_FC_BACKTRACE
  call backtrace()  ! Gfortran extension

#elif defined FC_INTEL
  call TRACEBACKQQ(USER_EXIT_CODE=-1)  ! Ifort extension
#endif

end subroutine show_backtrace
!!***

!----------------------------------------------------------------------

!!****f* m_errors/check_mpi_ierr
!! NAME
!!  check_mpi_ierr
!!
!! FUNCTION
!!  Basic error handler for MPI calls. This routine is usually interfaced through some macro defined in abi_common.h
!!
!! INPUTS
!!  ierr=Exit status reported by an MPI call.
!!  line=line number of the file where problem occurred
!!  file=name of the f90 file containing the caller
!!
!! OUTPUT
!!  Write error message thep stop execution.
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine check_mpi_ierr(ierr,msg,file,line)

!Arguments ------------------------------------
 integer,intent(in) :: ierr
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 integer,parameter :: mpi_msg_len=1000
 integer :: f90line,ilen,ierr2
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=mpi_msg_len) :: mpi_msg_error
 character(len=mpi_msg_len+500) :: my_msg
! *********************************************************************

#ifdef HAVE_MPI
 if (ierr==MPI_SUCCESS) RETURN
 call MPI_ERROR_STRING(ierr, mpi_msg_error, ilen, ierr2)
#else
 ilen=0; ierr2=0
 mpi_msg_error = " Check_mpi_ierr should not be called in non-MPI mode!"
 if (ierr==0) RETURN
#endif

 if (ilen>mpi_msg_len) write(std_out,*)" Warning_ MPI message has been truncated!"
 if (ierr2/=0) write(std_out,*)" Warning: MPI_ERROR_STRING returned ierr2= ",ierr2

 f90line=0; if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name = basename(file)
 my_msg = TRIM(msg)//ch10//TRIM(mpi_msg_error)

 call msg_hndl(my_msg,"ERROR","PERS",file=f90name,line=f90line)

end subroutine check_mpi_ierr
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_int
!! NAME
!!  unused_int
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar integer value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! SOURCE

elemental subroutine unused_int(var)

!Arguments ------------------------------------
 integer,intent(in) :: var

!Local variables-------------------------------
 integer :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_int
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_real_dp
!! NAME
!!  unused_real_dp
!!
!! FUNCTION
!!  Helper function used to silence warning messages due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar real value.
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      signal
!!
!! SOURCE

elemental subroutine unused_real_dp(var)

!Arguments ------------------------------------
 real(dp),intent(in) :: var

!Local variables-------------------------------
 real(dp) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_real_dp
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_real_sp
!! NAME
!!  unused_real_sp
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro. Target: one-dimensional real(dp) vector.
!!
!! SOURCE

elemental subroutine unused_real_sp(var)

!Arguments ------------------------------------
 real(sp),intent(in) :: var

!Local variables-------------------------------
 real(sp) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_real_sp
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_spc
!! NAME
!!  unused_cplx_spc
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar complex value
!!
!! OUTPUT
!!  None
!!
!! SOURCE

elemental subroutine unused_cplx_spc(var)

!Arguments ------------------------------------
 complex(spc),intent(in) :: var

!Local variables-------------------------------
 complex(spc) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_cplx_spc
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_cplx_dpc
!! NAME
!!  unused_cplx_dpc
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar complex value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      signal
!!
!! SOURCE

elemental subroutine unused_cplx_dpc(var)

!Arguments ------------------------------------
 complex(dpc),intent(in) :: var

!Local variables-------------------------------
 complex(dpc) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_cplx_dpc
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_logical
!! NAME
!!  unused_logical
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar logical value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      signal
!!
!! SOURCE

elemental subroutine unused_logical(var)

!Arguments ------------------------------------
 logical,intent(in) :: var

!Local variables-------------------------------
 logical :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_logical
!!***

!----------------------------------------------------------------------

!!****f* m_errors/unused_ch
!! NAME
!!  unused_ch
!!
!! FUNCTION
!!  Helper function used to silence compiler warnings due to unused variables.
!!  Interfaced via the ABI_UNUSED macro.
!!
!! INPUTS
!!  var=Scalar character value
!!
!! OUTPUT
!!  None
!!
!! PARENTS
!!
!! CHILDREN
!!      signal
!!
!! SOURCE

elemental subroutine unused_ch(var)

!Arguments ------------------------------------
 character(len=*),intent(in) :: var

!Local variables-------------------------------
 character(len=LEN(var)) :: dummy
! *********************************************************************

 dummy = var

end subroutine unused_ch
!!***

!----------------------------------------------------------------------

!!****f* m_errors/abietsf_msg_hndl
!! NAME
!!  abietsf_msg_hndl
!!
!! FUNCTION
!!  Wrapper to interface the abinint error handlers with the error handling routines used in etsf-io.
!!  It is usually interfaced via the macro ETSF_* defined in abi_common.h
!!
!! INPUTS
!!  lstat=Logical flag returned by etsf-io routines.
!!  Error_data<ETSF_io_low_error>=Structure storing the error returned by etsf-io calls.
!!  [line]=line number of the file where the problem occurred
!!  [file]=name of the f90 file containing the caller
!!  mode_paral=Either "COLL" or "PERS".
!!
!! OUTPUT
!!  Only writing, then the code is stopped.
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

#if defined HAVE_ETSF_IO

subroutine abietsf_msg_hndl(lstat,Error_data,mode_paral,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),optional,intent(in) :: file
 character(len=*),intent(in) :: mode_paral
 logical,intent(in) :: lstat
 type(ETSF_io_low_error),intent(in) :: Error_data

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=etsf_io_low_error_len) :: errmess
! *********************************************************************

 if (lstat) RETURN

 if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name = file
 call etsf_io_low_error_to_str(errmess,Error_data)

 call msg_hndl(errmess,"ERROR",mode_paral,f90name,f90line)

end subroutine abietsf_msg_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_errors/abietsf_warn
!! NAME
!!  abietsf_warn
!!
!! FUNCTION
!!  Wrapper to write warning messages, only used for ETSF_IO routines
!!  It is usually interfaced via the macro ETSF_WARN defined in abi_common.h
!!
!! INPUTS
!!  lstat=status error.
!!  Error_data<ETSF_io_low_error>=Structure storing the error returned by etsf-io calls.
!!  [line]=line number of the file where the problem occurred
!!  [file]=name of the f90 file containing the caller
!!  mode_paral=Either "COLL" or "PERS".
!!
!! OUTPUT
!!  Only writing.
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE


subroutine abietsf_warn(lstat,Error_data,mode_paral,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 logical,intent(in) :: lstat
 character(len=*),optional,intent(in) :: file
 character(len=*),intent(in) :: mode_paral
 type(ETSF_io_low_error),intent(in) :: Error_data

!Local variables-------------------------------
 integer :: f90line=0
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=etsf_io_low_error_len) :: errmess
! *********************************************************************

 if (lstat) RETURN

 if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name = file
 call etsf_io_low_error_to_str(errmess,Error_data)

 call msg_hndl(errmess,"WARNING",mode_paral,f90name,f90line)

end subroutine abietsf_warn
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_errors/bigdft_lib_error
!! NAME
!!  bigdft_lib_error
!!
!! FUNCTION
!!  Stop the code if bigdft library has not been enabled.
!!  Interfaced with the CPP macro BIGDFT_NOTENABLED_ERROR
!!
!! INPUTS
!!  line=line number of the file where problem occurred
!!  file=name of the f90 file containing the caller
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine bigdft_lib_error(file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),optional,intent(in) :: file

!Local variables-------------------------------
 character(len=500) :: message

! *********************************************************************

  write(message,'(4a)') ch10,&
&  ' BigDFT support has not been enabled.', ch10, &
&  ' Action, used the flag --enable-bigdft when configuring.'

 if (PRESENT(file) .and. PRESENT(line)) then
   call msg_hndl(message,"ERROR","PERS",file=file,line=line)
 else
   call msg_hndl(message,"ERROR", "PERS")
 end if

end subroutine bigdft_lib_error
!!***

!----------------------------------------------------------------------

!!****f* m_errors/xlf_set_sighandler
!! NAME
!!  xlf_set_sighandler
!!
!! FUNCTION
!!   Set the signal handler for IBM XLF
!!
!! NOTES
!!   See http://publib.boulder.ibm.com/infocenter/compbgpl/v9v111/index.jsp?topic=/com.ibm.xlf111.bg.doc/xlfopg/fptrap.htm
!!   The XL Fortran exception handlers and related routines are:
!!   xl__ieee
!!   Produces a traceback and an explanation of the signal and continues execution by supplying the default IEEE result
!!   for the failed computation. This handler allows the program to produce the same results as if exception detection was not turned on.
!!   xl__trce
!!   Produces a traceback and stops the program.
!!   xl__trcedump
!!   Produces a traceback and a core file and stops the program.
!!   xl__sigdump
!!   Provides a traceback that starts from the point at which it is called and provides information about the signal.
!!   You can only call it from inside a user-written signal handler.
!!   It does not stop the program. To successfully continue, the signal handler must perform some cleanup after calling this subprogram.
!!   xl__trbk
!!   Provides a traceback that starts from the point at which it is called.
!!   You call it as a subroutine from your code, rather than specifying it with the -qsigtrap option. It requires no parameters. It does not stop the program.
!!
!! PARENTS
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine xlf_set_sighandler()

! *************************************************************************

#ifdef FC_IBM
 call SIGNAL(SIGTRAP, xl__trcedump)
 call SIGNAL(SIGFPE, xl__trcedump)
#endif

end subroutine xlf_set_sighandler
!!***

!----------------------------------------------------------------------

!!****f* m_errors/abinit_doctor
!! NAME
!!  abinit_doctor
!!
!! FUNCTION
!! Perform checks on memory leaks and leaking file descriptors at the end of the run.
!!
!! INPUTS
!!  prefix=Prefix for output file  (usually "__nameofprogram" e.g. __cut3d)
!!  [print_mem_report]=0 to disable the test on memory leaks (used in Abinit if bigdft is activated).
!!    Default: 1, i.e. memory check is always activated.
!!
!! PARENTS
!!      abinit,anaddb,conducti,cut3d,fftprof,fold2Bloch,ioprof,lapackprof
!!      macroave,mrgddb,mrgdv,mrggkk,mrgscr,multibinit,optic,ujdet
!!      vdw_kernelgen
!!
!! CHILDREN
!!      abimem_get_info,abimem_shutdown,show_units,wrtout
!!
!! SOURCE

subroutine abinit_doctor(prefix, print_mem_report)

!Arguments ------------------------------------
 integer,optional,intent(in) :: print_mem_report
 character(len=*),intent(in) :: prefix

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer :: do_mem_report, my_rank
 character(len=5000) :: msg
#ifdef HAVE_MEM_PROFILING
 integer :: ii,ierr,unt
 integer(i8b) :: memtot, nalloc, nfree
 character(len=fnlen) :: path
 character(len=5000) :: errmsg
#endif

! *************************************************************************

 do_mem_report = 1; if (present(print_mem_report)) do_mem_report = print_mem_report
 my_rank = xmpi_comm_rank(xmpi_world)

#ifdef HAVE_MEM_PROFILING
 errmsg = ""; ierr = 0

 ! Test on memory leaks.
 call abimem_get_info(nalloc, nfree, memtot)
 call abimem_shutdown()

 if (do_mem_report == 1) then
   if (nalloc == nfree .and. memtot == 0) then
     write(msg,'(3a,i0,a,i0,3a,i0)') &
       '- MEMORY CONSUMPTION REPORT:',ch10, &
       '-   There were ',nalloc,' allocations and ',nfree,' deallocations',ch10, &
       '-   Remaining memory at the end of the calculation is ',memtot
   else
     ! This msg will make the test fail if the memory leak occurs on master (no dash in the first column)
     write(msg,'(2a,2(a,i0),3a,f12.4,1x,11a)') &
       'MEMORY CONSUMPTION REPORT:',ch10, &
       '   There were ',nalloc,' allocations and ',nfree,' deallocations',ch10, &
       '   Remaining memory at the end of the calculation: ',memtot * b2Mb, " (Mb)", ch10, &
       '   As a help for debugging, you might set call abimem_init(2) in the main program,', ch10, &
       '   or use the command line option `abinit --abimem-level 2`', ch10, &
       '   then use tests/Scripts/abimem.py to analyse the file abimem_rank[num].mocc that has been created,',ch10, &
       '   e.g. from tests/Scripts issue the command: ./abimem.py leaks ../<dir>/<subdir>/abimem_rank0.mocc .',ch10, &
       '   Note that abimem files can easily be multiple GB in size so do not use this option normally!'
     ! And this will make the code call mpi_abort if the leak occurs on my_rank != master
     ierr = ierr + 1
     errmsg = strcat(errmsg, ch10, msg)
   end if

 else
   write(msg,'(3a)') &
     '- MEMORY CONSUMPTION REPORT:',ch10, &
     '- Memory profiling is activated but not yet usable when bigdft is used'
 end if
 if (my_rank == master) call wrtout(ab_out, msg)

 ! Test whether all logical units have been closed.
 ! If you wonder why I'm doing this, remember that there's a per-user
 ! limit on the maximum number of open file descriptors. Hence descriptors
 ! represent a precious resource and we should close them as soon as possible.
 ii = num_opened_units(ignore=[std_err, std_in, std_out, ab_out])
 if (ii > 0) then
   path = strcat(prefix, "_lunits_rank", itoa(my_rank), ".flun")
   if (open_file(path, msg, newunit=unt) /= 0) then
     MSG_ERROR(msg)
   end if
   call show_units(unt)
   close(unt)
   write(msg, "(a,i0,2a)")"Leaking ",ii," Fortran logical units. See: ",trim(path)
   errmsg = strcat(errmsg, ch10, msg)
   ierr = ierr + 1
 end if

 if (my_rank == master) call wrtout(ab_out, msg)
 if (ierr /= 0) then
   MSG_ERROR(errmsg)
 end if

#else
 ABI_UNUSED(prefix)
#endif

 ! Check for pending requests.
 if (xmpi_count_requests /= 0) then
   write(msg, "(a,i0,a)")"Leaking ", xmpi_count_requests, " MPI requests at the end of the run"
   MSG_WARNING(msg)
   !MSG_ERROR(msg)
#ifdef HAVE_MEM_PROFILING
   MSG_ERROR(msg)
#endif
 end if

end subroutine abinit_doctor
!!***

!!****f* m_errors/abi_abort
!! NAME
!!  abi_abort
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible parallelization.
!!
!!  Note the this routine is private and should never be called explicitly.
!!  Please, use the macros:
!!    MSG_ERROR, MSG_BUG
!!  defined in abi_common.h to abort the execution.
!!  XG : this is not true, in very rare cases, ABINIT has to exit without giving an error (e.g. for non-zero prtkpt )
!!
!! INPUTS
!!  exit_status=(optional, default=1 or -1, see below) the return code of the routine
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same message to be
!!     written once only or
!!   'PERS' if the procs are calling the routine with different mesgs
!!     each to be written, or if one proc is calling the routine
!!  print_config=(optional, default=true)
!!       if true print out several information before leaving
!!
!! OUTPUT
!!  (only writing, then stop)
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      m_errors,testkgrid,vtorho
!!
!! CHILDREN
!!      dump_config,print_kinds,wrtout,xmpi_abort,xmpi_show_info
!!
!! SOURCE

subroutine abi_abort(mode_paral,exit_status,print_config)

!Arguments ------------------------------------
 character(len=4),intent(in) :: mode_paral
 integer,intent(in),optional :: exit_status
 logical,intent(in),optional :: print_config

!Local variables-------------------------------
 logical :: print_config_

! **********************************************************************

 call wrtout(std_out,ch10//' abi_abort: decision taken to exit ...','PERS')

! Caveat: Do not use MPI collective calls!
 if (mode_paral == "COLL") then
   call wrtout(std_out,"Why are you using COLL? Are you sure that ALL the processors are calling abi_abort?")
 end if

!Dump configuration before exiting
 print_config_=.False.; if (present(print_config)) print_config_=print_config
 if (print_config_) then
   call print_kinds()
   call xmpi_show_info()
   call dump_config(std_out)
 end if

 if (present(exit_status)) then
   call xmpi_abort(exit_status=exit_status)
 else
   call xmpi_abort()
 end if

end subroutine abi_abort
!!***

!!****f* m_errors/abi_cabort
!! NAME
!!  abi_cabort
!!
!! FUNCTION
!!  C-interoperable version of abi_abort

subroutine abi_cabort() bind(C, name='abi_cabort')

  call abi_abort("COLL", exit_status=1, print_config=.False.)

end subroutine abi_cabort
!!***

END MODULE m_errors
!!***

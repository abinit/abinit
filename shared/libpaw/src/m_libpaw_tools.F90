!!****m* ABINIT/m_libpaw_tools
!! NAME
!!  m_libpaw_tools
!!
!! FUNCTION
!!  Several libPAW tools: message printing, error handling, string handling...
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (MT, MG, ...)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Parts of this module come directly from hide_write & hide_leave src files delivered with ABINIT.
!!
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_libpaw_tools

 USE_DEFS
 USE_MPI_WRAPPERS

#if defined HAVE_YAML
  use yaml_output
#endif
#ifdef LIBPAW_HAVE_NETCDF
  use netcdf
#endif

 implicit none

 private

!PUBLIC FUNCTIONS - MESSAGE HANDLING
 public  :: libpaw_wrtout         ! Parallel output of messages
 public  :: libpaw_msg_hndl       ! Basic error handler
 public  :: libpaw_flush          ! Wrapper for the standard flush routine
 public  :: libpaw_spmsg_getcount ! Get number of special messages (WARNING/COMMENT) already printed
 public  :: libpaw_spmsg_mpisum   ! Reduce number of special messages (WARNING/COMMENT) over MPI comm
 public  :: libpaw_write_comm_set ! Set the MPI communicator used for parallel write
 public  :: libpaw_log_flag_set   ! Set the flag controlling the filling of the LOG file
 public  :: libpaw_netcdf_check   ! Stop execution after a NetCDF I/O error

!PUBLIC FUNCTIONS - STRING HANDLING
 public  :: libpaw_basename       ! String, base name extraction from path
 public  :: libpaw_to_upper       ! String conversion to uppercase
 public  :: libpaw_lstrip         ! String right blanks removal
 public  :: libpaw_indent         ! String indentation

!PUBLIC FUNCTIONS - IO TOOLS
 public :: libpaw_get_free_unit  ! Get a free I/O unit

!PRIVATE FUNCTIONS
 private :: libpaw_wrtout_myproc  ! Sequential output of messages
 private :: libpaw_write_lines    ! OS-compatible string output
 private :: libpaw_leave          ! Clean exit of F90 routines
 private :: libpaw_die            ! Clean exit
 private :: libpaw_lock_and_write ! Write a string to a file with locking mechanism

!PRIVATE VARIABLES
 integer,save :: LIBPAW_WRITE_COMM=xmpi_world ! Communicator used for the parallel write
 integer,save :: LIBPAW_COMMENT_COUNT=0           ! Number of COMMENTs printed in log file
 integer,save :: LIBPAW_WARNING_COUNT=0           ! Number of WARNINGs printed in log file
 integer,save :: LIBPAW_EXIT_FLAG=0               ! Flag set to 1 if an exit is requested
 logical,save :: LIBPAW_HAS_LOG_FILE=.TRUE.       ! Flag: True if std output exists

!PRIVATE PARAMETERS
 integer,parameter :: LIBPAW_NULL_UNIT=-1     ! Fake null unit
 character(len=25),parameter :: LIBPAW_MPIABORTFILE="__LIBPAW_MPIABORTFILE__"
#if defined HAVE_OS_WINDOWS
 character(len=3),parameter :: LIBPAW_NULL_FILE="NUL"
#else
 character(len=9),parameter :: LIBPAW_NULL_FILE="/dev/null"
#endif

!!***

CONTAINS !===========================================================

!!****f* m_libpaw_tools/libpaw_wrtout
!! NAME
!! libpaw_wrtout
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRTOUT routine delivered with ABINIT.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_wrtout(unit,msg,mode_paral)

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral

!Local variables ------------------------------
 integer :: comm,me,nproc
 integer,save :: master=0
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if ((unit==std_out).and.(.not.LIBPAW_HAS_LOG_FILE)) RETURN
 if (unit==LIBPAW_NULL_UNIT) RETURN

 my_mode_paral = "COLL"; if (PRESENT(mode_paral)) my_mode_paral = mode_paral

!Communicator used for the parallel write
 comm=LIBPAW_WRITE_COMM
 nproc = xmpi_comm_size(comm)
 me    = xmpi_comm_rank(comm)

 if ((my_mode_paral=='COLL').or.(nproc==1)) then
   if (me==master) then
     call libpaw_wrtout_myproc(unit,msg)
   end if
 else if (my_mode_paral=='PERS') then
   call libpaw_write_lines(unit,msg)
 else if (my_mode_paral=='INIT') then
   master=unit
 else
   write(string,'(7a)')ch10,&
&   'libpaw_wrtout: ERROR -',ch10,&
&   '  Unknown write mode: ',my_mode_paral,ch10,&
&   '  Continuing anyway ...'
   write(unit,'(A)') trim(string)
 end if

end subroutine libpaw_wrtout
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_wrtout_myproc
!! NAME
!!  libpaw_wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc.
!!
!! INPUTS
!!  unit=unit number for writing
!!  msg=(character(len=*)) message to be written
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRTOUT_MYPROC routine delivered with ABINIT.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_wrtout_myproc(unit,msg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg

!Local variables ------------------------------
!scalars
 logical :: print_std_err
!arrays

!******************************************************************

 print_std_err=(unit==std_out.and.(index(trim(msg),'BUG')/=0.or.index(trim(msg),'ERROR')/=0))

!Print message
 call libpaw_write_lines(unit,msg)
 if (print_std_err) then
   call libpaw_write_lines(std_err,msg)
 end if

!Append "Contact Abinit group" to BUG messages
 if (index(trim(msg),'BUG')/=0) then
   write(unit,'(a)') '  Action: contact libPAW developers.'
   if (print_std_err) write(std_err, '(a)' ) '  Action: contact libPAW developers.'
   write(unit,*); if (print_std_err) write(std_err,*)
 end if

!Count the number of warnings and comments. Only take into
!account unit std_out, in order not to duplicate these numbers.
 if (index(trim(msg),'WARNING')/=0 .and. unit==std_out) LIBPAW_WARNING_COUNT=LIBPAW_WARNING_COUNT+1
 if (index(trim(msg),'COMMENT')/=0 .and. unit==std_out) LIBPAW_COMMENT_COUNT=LIBPAW_COMMENT_COUNT+1
 if (index(trim(msg),'Exit'   )/=0) LIBPAW_EXIT_FLAG=1

end subroutine libpaw_wrtout_myproc
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_write_lines
!! NAME
!!  libpaw_write_lines
!!
!! FUNCTION
!!  This routine receives a string, split the message in lines according to the
!!  ch10 character and output the text to the specified unit.
!!  Allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! INPUTS
!!  unit=unit number for writing
!!  msg=(character(len=*)) message to be written
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the WRITE_LINES routine delivered with ABINIT.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_write_lines(unit,msg)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg

!Local variables ------------------------------
!scalars
 integer :: msg_size,ii,jj,rtnpos
#if defined HAVE_YAML
 character(len = len_trim(msg)) :: msg_out
#endif

!******************************************************************

 msg_size=len_trim(msg)

#if defined HAVE_YAML
 if (msg_size>0 .and. unit==std_out) then
    ! Change any carriage return into space.
    do ii = 1, msg_size
       if (msg(ii:ii) /= char(10)) then
          msg_out(ii:ii) = msg(ii:ii)
       else
          msg_out(ii:ii) = " "
       end if
    end do
    call yaml_comment(msg_out)
 end if
 return
#endif

 if (msg_size==0) then
   write(unit,*) ; return
 end if

!Here, split the message, according to the char(10) characters (carriage return).
!This technique is portable accross different OS.
 rtnpos=index(msg,ch10)
 if (rtnpos==0) then
   write(unit,"(a)") msg(1:msg_size) ; return
 end if

 ii=1; jj=rtnpos
 do
   if (ii==jj) then
     write(unit,*)
   else
     write(unit,'(a)') msg(ii:jj-1)
   end if
   ii=jj+1 ; if (ii>msg_size) exit
   jj=index(msg(ii:msg_size),ch10)
   if (jj==0) then
     jj=msg_size+1
   else
     jj=jj+ii-1
   end if
 end do

 if (msg(msg_size:msg_size)==ch10) write(unit,*)

end subroutine libpaw_write_lines
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_msg_hndl
!! NAME
!!  libpaw_msg_hndl
!!
!! FUNCTION
!!  Basic error handler.
!!
!! INPUTS
!!  msg=string containing additional information on the nature of the problem
!!  level=string defining the type of problem. Possible values are:
!!   COMMENT, WARNING, ERROR,BUG
!!  mode_paral=Either "COLL" or "PERS".
!!  [line]=line number of the file where problem occurred (optional)
!!  [file]=name of the f90 file containing the caller (optional)
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the MSG_HNDL routine delivered with ABINIT.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_msg_hndl(msg,level,mode_paral,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: level,msg,mode_paral
 character(len=*),optional,intent(in) :: file

!Local variables ------------------------------
 logical :: file_exists
 character(len=500) :: f90name='Unknown'
 character(len=LEN(msg)) :: my_msg
 character(len=MAX(4*LEN(msg),2000)) :: sbuf

! *********************************************************************

 my_msg=libpaw_lstrip(msg)

 write(sbuf,'(3a)') ch10,"--- !",TRIM(level)
 if (PRESENT(file)) then
   f90name=libpaw_basename(file)
   write(sbuf,'(4a)') trim(sbuf),ch10,"src_file: ",TRIM(f90name)
 end if
 if (PRESENT(line)) then
   write(sbuf,'(3a,i0)') trim(sbuf),ch10,"src_line: ",line
 end if
 write(sbuf,'(8a)') trim(sbuf),ch10,&
&  "message: |",ch10,trim(libpaw_indent(my_msg)),ch10,&
&  "...",ch10

 select case (libpaw_to_upper(level))
 case ('COMMENT','WARNING')
   call libpaw_wrtout(std_out,sbuf,mode_paral)
 case ('ERROR','BUG')
   call libpaw_wrtout(std_out,sbuf,mode_paral)
   inquire(file=LIBPAW_MPIABORTFILE,exist=file_exists)
   if ((.not.file_exists).and.xmpi_comm_size(xmpi_world)>1) then
     call libpaw_lock_and_write(LIBPAW_MPIABORTFILE,sbuf)
   end if
   call libpaw_leave(mode_paral)
 case default
   write(sbuf,'(4a)') ch10,' libpaw_msg_hndl: BUG**2 - ',ch10,' Wrong value for level!'
   call libpaw_die(sbuf)
 end select

end subroutine libpaw_msg_hndl
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_spmsg_getcount
!! NAME
!!  libpaw_spmsg_getcount
!!
!! FUNCTION
!!  Get the values of the counters of special messages (WARNING, COMMENT)
!!
!! INPUTS
!!  ncomment= number of COMMENTs in log file
!!  nwarning= number of WARNINGs in log file
!!  nexit=    1 if exit requested
!!
!! OUTPUT
!!  (only counters updated)
!!
!! NOTES
!!  This routine comes directly from the SPECIALMSG_GETCOUNT routine delivered with ABINIT.
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_spmsg_getcount(ncomment,nwarning,nexit)

!Arguments ------------------------------------
 integer,intent(out) :: ncomment,nexit,nwarning

!Local variables ------------------------------

! **********************************************************************

 ncomment=LIBPAW_COMMENT_COUNT
 nwarning=LIBPAW_WARNING_COUNT
 nexit   =LIBPAW_EXIT_FLAG

end subroutine libpaw_spmsg_getcount
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_spmsg_mpisum
!! NAME
!!  libpaw_spmsg_mpisum
!!
!! FUNCTION
!!  Reduce the counters of special messages (WARNING, COMMENTS, EXIT) over a MPI communicator
!!
!! INPUTS
!!  mpicomm= MPI communicator
!!
!! OUTPUT
!!  (only counters updated)
!!
!! NOTES
!!  This routine comes directly from the SPECIALMSG_MPISUM routine delivered with ABINIT.
!!
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_spmsg_mpisum(mpicomm)

!Arguments ------------------------------------
 integer,intent(in) :: mpicomm

!Local variables ------------------------------
 integer :: ierr
 integer :: buf(3)

! **********************************************************************

  buf(1)=LIBPAW_COMMENT_COUNT;buf(2)=LIBPAW_WARNING_COUNT;buf(3)=LIBPAW_EXIT_FLAG

  call xmpi_sum(buf,mpicomm,ierr)

  LIBPAW_COMMENT_COUNT=buf(1)
  LIBPAW_WARNING_COUNT=buf(2)
  LIBPAW_EXIT_FLAG=buf(3) ; if (LIBPAW_EXIT_FLAG/=0) LIBPAW_EXIT_FLAG=1

end subroutine libpaw_spmsg_mpisum
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_write_comm_set
!! NAME
!!  libpaw_write_comm_set
!!
!! FUNCTION
!!  Set the MPI communicator used for parallel write
!!
!! INPUTS
!!  new_write_comm= new value for the parallel write MPI communicator
!!
!! OUTPUT
!!
!! PARENTS
!!      m_driver,m_io_redirect,m_memeval,m_mpi_setup,m_mpinfo
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_write_comm_set(new_write_comm)

!Arguments ------------------------------------
 integer,intent(in) :: new_write_comm

!Local variables ------------------------------

! **********************************************************************

 LIBPAW_WRITE_COMM=new_write_comm

end subroutine libpaw_write_comm_set
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_log_flag_set
!! NAME
!!  libpaw_log_flag_set
!!
!! FUNCTION
!!  Set the flag controlling the filling of the LOG file
!!
!! INPUTS
!!  log_flag= new value for LOG file flag
!!            True: the log file is filled; False: no the log file
!!
!! OUTPUT
!!
!! PARENTS
!!      m_argparse,m_dtfil
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_log_flag_set(log_flag)

!Arguments ------------------------------------
 logical,intent(in) :: log_flag

!Local variables ------------------------------

! **********************************************************************

 LIBPAW_HAS_LOG_FILE=log_flag

end subroutine libpaw_log_flag_set
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tool/libpaw_netcdf_check
!! NAME
!!  libpaw_netcdf_check
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
!!  This routine is usually interfaced with the macros defined in libpaw.h
!!
!! PARENTS
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_netcdf_check(ncerr,msg,file,line)

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

#ifdef LIBPAW_HAVE_NETCDF
 if (ncerr /= NF90_NOERR) then
   if (PRESENT(line)) then
     f90line=line
   else
     f90line=0
   end if
   if (PRESENT(file)) then
     f90name = libpaw_basename(file)
   else
     f90name='Subroutine Unknown'
   end if
   !
   ! Append Netcdf string to user-defined message.
   write(nc_msg,'(a,3x,a)')' - NetCDF library returned:',TRIM(nf90_strerror(ncerr))
   my_msg = TRIM(msg) // TRIM(nc_msg)

   call libpaw_msg_hndl(my_msg,"ERROR","PERS",f90name,f90line)
 end if
#else
 call libpaw_die("LIBPAW_HAVE_NETCDF is not defined!")
#endif

end subroutine libpaw_netcdf_check
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_leave
!! NAME
!!  libpaw_leave
!!
!! FUNCTION
!!  Routine for clean exit of f90 code, taking into account possible parallelization.
!!
!! INPUTS
!!  mode_paral=
!!   'COLL' if all procs are calling the routine with the same msg to be written once only
!!   'PERS' if the procs are calling the routine with different msgs each to be written,
!!          or if one proc is calling the routine
!!  [exit_status]=(optional, default=1 or -1, see below) the return code of the routine
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!!  This routine comes directly from the LEAVE_NEW routine delivered with ABINIT.
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_leave(mode_paral,exit_status)

!Arguments ------------------------------------
 integer,intent(in),optional :: exit_status
 character(len=4),intent(in) :: mode_paral

!Local variables ------------------------------

! **********************************************************************

 call libpaw_wrtout(std_out,ch10//' leave_new : decision taken to exit ...','PERS')

!Caveat: Do not use MPI collective calls!
 if (mode_paral=="COLL") then
   call libpaw_wrtout(std_out,"Why COLL? Are you sure that ALL the processors are calling leave_new?")
 end if

 if (present(exit_status)) then
   call xmpi_abort(exit_status=exit_status)
 else
   call xmpi_abort()
 end if

end subroutine libpaw_leave
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_die
!! NAME
!!  libpaw_die
!!
!! FUNCTION
!!  Stop smoothly the execution in case of unexpected events reporting the
!!  line number and the file name where the error occurred as well as the
!!  MPI rank of the processor.
!!
!! INPUTS
!!  msg=String containing additional information on the nature of the problem
!!  [file]=Name of the f90 file containing the caller
!!  [line]=Line number of the file where problem occurred
!!
!! NOTES
!!  This routine comes directly from the DIE routine delivered with ABINIT.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_die(message,file,line)

!Arguments ------------------------------------
 integer,optional,intent(in) :: line
 character(len=*),intent(in) :: message
 character(len=*),optional,intent(in) :: file

!Local variables ------------------------------
 integer :: rank
 integer :: f90line=0
 character(len=10) :: lnum,strank
 character(len=500) :: f90name='Subroutine Unknown'
 character(len=500) :: msg

! *********************************************************************

 if (PRESENT(line)) f90line=line
 if (PRESENT(file)) f90name= libpaw_basename(file)

 rank=xmpi_comm_rank(xmpi_world) !Determine my rank inside world communicator

 write(lnum,"(i0)") f90line
 write(strank,"(i0)") rank
 msg=TRIM(f90name)//':'//TRIM(lnum)//' P'//TRIM(strank)
 write(msg,'(a,2x,2a,2x,a)') ch10,TRIM(msg),ch10,TRIM(message)

 call libpaw_wrtout(std_out,msg,'PERS')
 call libpaw_leave('PERS')

end subroutine libpaw_die
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_lock_and_write
!! NAME
!!  libpaw_lock_and_write
!!
!! FUNCTION
!!  Writes a string to filename with locking mechanism.
!!
!! INPUTS
!!  filename= Name of the file.
!!  string= Input string.
!!
!! PARENTS
!!      m_libpaw_tools
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_lock_and_write(filename,string)

!Arguments ------------------------------------
 character(len=*),intent(in) :: filename,string

!Local variables-------------------------------
 integer :: lock_unit,file_unit
 character(len=len(filename)+5) :: lock

! *********************************************************************

 !Try to acquire the lock.
 lock=trim(filename)//".lock"
 lock_unit=libpaw_get_free_unit()
 open(unit=lock_unit,file=trim(lock),status='new',err=99)

 file_unit=libpaw_get_free_unit()
 open(unit=file_unit,file=trim(filename),form="formatted")
 call libpaw_write_lines(file_unit,string)
 close(lock_unit,status="delete")
 close(file_unit)
 return

99 continue

end subroutine libpaw_lock_and_write
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_get_free_unit
!! NAME
!!  libpaw_get_free_unit
!!
!! FUNCTION
!!  Obtain a free logical Fortran unit.
!!
!! OUTPUT
!!  The unit number (free unit)
!!  Raises:
!!   -1 if no logical unit is free (!)
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function libpaw_get_free_unit()

!Local variables-------------------------------
 integer,parameter :: MIN_UNIT_NUMBER=10
#ifdef FC_NAG
  integer,parameter :: MAX_UNIT_NUMBER=64    ! There's a serious problem in Nag6.0. In principle
                                             ! Maximum unit number: 2147483647
#else
 integer,parameter :: MAX_UNIT_NUMBER=1024
#endif
 integer :: iunt
 logical :: isopen

! *********************************************************************

 do iunt=MAX_UNIT_NUMBER,MIN_UNIT_NUMBER,-1
   inquire(unit=iunt,opened=isopen)
   if (.not.isopen) then
      libpaw_get_free_unit=iunt; return
   end if
 end do
 libpaw_get_free_unit=-1

end function libpaw_get_free_unit
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_flush
!! NAME
!!  libpaw_flush
!!
!! FUNCTION
!!  Wrapper for the standard flush routine
!!  Available only if the compiler implements this intrinsic procedure.
!!
!! INPUTS
!!  unit=Fortran logical Unit number
!!
!! NOTES
!!  This routine comes directly from the FLUSH_UNIT routine delivered with ABINIT.
!!
!! PARENTS
!!      m_pawrhoij
!!
!! CHILDREN
!!      flush,flush_
!!
!! SOURCE

subroutine libpaw_flush(unit)

!Arguments ------------------------------------
 integer,intent(in) :: unit

!Local variables ------------------------------
 integer, parameter :: dev_null=-1
 logical :: isopen

!************************************************************************

 if (unit==dev_null) return

!FLUSH on unconnected unit is illegal: F95 std., 9.3.5.
 inquire(unit=unit,opened=isopen)

#if defined HAVE_FC_FLUSH
 if (isopen) then
   call flush(unit)
 endif
#elif defined HAVE_FC_FLUSH_
 if (isopen) then
   call flush_(unit)
  end if
#endif

end subroutine libpaw_flush
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_basename
!! NAME
!! libpaw_basename
!!
!! FUNCTION
!!  Returns the final component of a pathname (function version).
!!
!! INPUTS
!!  string=The input string
!!
!! NOTES
!!  This routine comes directly from the BASENAME routine delivered with ABINIT.
!!  If the input string in not a valid path to a file, a blank strink is returned
!!
!! SOURCE

pure function libpaw_basename(istr) result(ostr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)) :: ostr

!Local variables ------------------------------
 integer :: ic,nch_trim,nch
 character(len=1),parameter :: BLANK=' '
 character(len=1),parameter :: DIR_SEPARATOR = '/'

!************************************************************************

 nch     =LEN     (istr)
 nch_trim=LEN_TRIM(istr)

 ic = INDEX (TRIM(istr), DIR_SEPARATOR, back=.TRUE.)
 if (ic >= 1 .and. ic <= nch_trim-1) then ! there is stuff after the separator.
   ostr = istr(ic+1:nch_trim)
 else if (ic==0 .or. ic == nch_trim+1) then ! no separator in string or zero length string,
   ostr = TRIM(istr)     ! return trimmed string.
 else                    ! (ic == nch_trim) separator is the last char.
   ostr = BLANK ! This is not a valid path to a file, return blank.
 end if
 return

end function libpaw_basename
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/to_upper
!! NAME
!!  libpaw_to_upper
!!
!! FUNCTION
!!  Convert a string to UPPER CASE (function version).
!!
!! INPUTS
!!   istr=Input string
!!
!! NOTES
!!  This routine comes directly from the TOUPPER routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_to_upper(istr) result(ostr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=LEN_TRIM(istr)) :: ostr

!Local variables ------------------------------
 integer,parameter :: ASCII_aa=ICHAR('a')
 integer,parameter :: ASCII_zz=ICHAR('z')
 integer,parameter :: SHIFT=ICHAR('a')-ICHAR('A')
 integer :: ic,iasc

! *********************************************************************

 do ic=1,LEN_TRIM(istr)
   iasc=IACHAR(istr(ic:ic))
   if (iasc>=ASCII_aa.and.iasc<=ASCII_zz) then
     ostr(ic:ic)=ACHAR(iasc-SHIFT)
   else
     ostr(ic:ic)=istr(ic:ic)
   end if
 end do

end function libpaw_to_upper
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_lstrip
!! NAME
!!  libpaw_lstrip
!!
!! FUNCTION
!!  Removes leading spaces from the input string.
!!
!! NOTES
!!  This routine comes directly from the LSTRIP routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_lstrip(istr) result(ostr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=len(istr)) :: ostr

!Local variables ------------------------------
 integer :: ii,jj,lg
 character(len=1),parameter :: BLANK=' '

! *********************************************************************

 lg=LEN(istr)
 do ii=1,lg
   if (istr(ii:ii)/=BLANK) EXIT
 end do

 ostr = " "
 do jj=1,lg-ii+1
   ostr(jj:jj) = istr(ii:ii)
   ii=ii+1
 end do

end function libpaw_lstrip
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_tools/libpaw_indent
!! NAME
!!  libpaw_indent
!!
!! FUNCTION
!!  Indent text (function version).
!!
!! INPUTS
!!   istr=Input string
!!
!! NOTES
!!  This routine comes directly from the INDENT routine delivered with ABINIT.
!!
!! SOURCE

pure function libpaw_indent(istr) result(ostr)

!Arguments ------------------------------------
 character(len=*),intent(in) :: istr
 character(len=len(istr)*4+4) :: ostr

!Local variables-------------------------------
 character(len=1),parameter :: NCHAR = char(10)
 integer,parameter :: n=4
 integer :: ii,jj,kk
 character(len=1) :: ch

! *********************************************************************

 ostr=" "
 jj=n
 do ii=1,LEN_TRIM(istr)
   ch=istr(ii:ii)
   jj=jj+1
   if (ch==NCHAR) then
      ostr(jj:jj)=NCHAR
      do kk=jj+1,jj+n
        ostr(kk:kk)=" "
      end do
      jj=jj+n
   else
     ostr(jj:jj)=ch
   end if
 end do

end function libpaw_indent
!!***

!----------------------------------------------------------------------

end module m_libpaw_tools
!!***

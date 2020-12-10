!!****m* ABINIT/m_exit
!! NAME
!! m_exit
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG, DCA, XG, GMR)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_exit

 use defs_basis
 use m_xmpi
 use m_abicore
 use m_errors

 use m_time,      only : abi_wtime, sec2str, timein
 use m_fstrings,  only : inupper
 use m_io_tools,  only : open_file

 implicit none

 private
!!***

 public :: exit_init                  ! Initialize the global variables of the module
 public :: get_start_time             ! Return the origin of time (begin of execution ) in seconds
 public :: have_timelimit_in          ! True if the run must completed before timelimit.
 public :: disable_timelimit          ! Disable time limit handlers.
 public :: enable_timelimit_in        ! Eable time limit handler in the given function.
 public :: get_timelimit              ! Return the time limit in seconds
 public :: get_timelimit_string       ! Return the time limit in string form.
 public :: exit_check                 ! Test if we should try to stop the code gracefully and to create a restart point.

 real(dp),private,save :: WALL0
! Origin of time in seconds.

 real(dp),private,save :: WTIME_LIMIT=-one
! Wall time limit in seconds. Negative value if not set.

 character(len=fnlen),private,save :: TIMELIMIT_INABIFUNC = "__None__"
! Name of the abinit function in which the time limit is handled.
! Note that there's **only one caller** in charge of the check.
! In structure relaxations, for example, only mover decides whether the itime loop must be exited
! and similar time limit handlers located in the children (e.g. scfcv) are automatically deactivated.
! This approach facilitates the treatment of nested parallelism e.g. image parallelism and reduce the number of MPI synchronizations.
! The drawback is that we loose the possibility of controlling the exit at a fine-grained level.
! If, for example, mover underestimates the time needed to complete the itime step, scfcv wont' be
! able to return if the time limit is approaching.

!----------------------------------------------------------------------

CONTAINS
!!***

!!****f* m_exit/exit_init
!! NAME
!! exit_init
!!
!! FUNCTION
!!  Initialize the global variables of the modules.
!!  This is a collective function that should be called by all the nodes in COMM_WORLD
!!
!! INPUTS
!!
!! PARENTS
!!      m_argparse
!!
!! CHILDREN
!!      inupper,timein,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine exit_init(time_limit)

!Arguments ------------------------------------
 real(dp),intent(in) :: time_limit

! *************************************************************************

 WTIME_LIMIT = time_limit
 WALL0 = abi_wtime()

end subroutine exit_init
!!***

!----------------------------------------------------------------------

!!****f* m_exit/disable_timelimit
!! NAME
!! disable_timelimit
!!
!! FUNCTION
!!  Disable the time limit handler. This function should be called by a driver
!!  routine that is not able to handle time limit events and wants to prevent
!!  its children from installing their handlers.
!!
!! PARENTS
!!      m_dfpt_looppert
!!
!! CHILDREN
!!      inupper,timein,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine disable_timelimit()

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 WTIME_LIMIT = -one

 if (TIMELIMIT_INABIFUNC /= "__None__") then
   msg = "Timelimit is already activated in function: "//trim(TIMELIMIT_INABIFUNC)
   ABI_WARNING(msg)
   !ABI_ERROR(msg)
 end if

end subroutine disable_timelimit
!!***

!----------------------------------------------------------------------

!!****f* m_exit/have_timelimit_in
!! NAME
!! have_timelimit_in
!!
!! FUNCTION
!!  Return .True. if timelimit is enabled in this caller
!!
!! PARENTS
!!
!! SOURCE

logical pure function have_timelimit_in(abifunc) result(ans)

!Arguments -----------------------------------
 character(len=*),intent(in) :: abifunc

! *************************************************************************

 ans = WTIME_LIMIT > zero .and. abifunc == TIMELIMIT_INABIFUNC

end function have_timelimit_in
!!***

!----------------------------------------------------------------------

!!****f* m_exit/enable_timelimit_in
!! NAME
!! enable_timelimit_in
!!
!! FUNCTION
!!  Eable time limit handler in the given function if not already done in one of the callers.
!!  Return the name of procedure that is handling the time limit.
!!  Example:
!!
!!    ! enable time limit handler if not done in callers.
!!    if (enable_timelimit_in(FUNC_NAME) == FUNC_NAME) then
!!      write(std_out,*)"Enabling timelimit check in function: ",trim(FUNC_NAME)," with timelimit: ",trim(sec2str(get_timelimit()))
!!    end if
!!
!! PARENTS
!!
!! SOURCE

function enable_timelimit_in(abifunc) result(prev_func)

!Arguments -----------------------------------
 character(len=*),intent(in) :: abifunc
 character(len=fnlen) :: prev_func

! *************************************************************************

 if (WTIME_LIMIT > zero .and. TIMELIMIT_INABIFUNC == "__None__") TIMELIMIT_INABIFUNC = abifunc
 prev_func = TIMELIMIT_INABIFUNC

end function enable_timelimit_in
!!***

!----------------------------------------------------------------------

!!****f* m_exit/get_timelimit
!! NAME
!! get_timelimit
!!
!! FUNCTION
!!  Return the time limit in seconds
!!
!! PARENTS
!!
!! SOURCE

real(dp) pure function get_timelimit()

 get_timelimit = WTIME_LIMIT

end function get_timelimit
!!***

!----------------------------------------------------------------------

!!****f* m_exit/get_timelimit_string
!! NAME
!! get_timelimit_string
!!
!! FUNCTION
!!  Return the time limit in string form.
!!
!! PARENTS
!!
!! SOURCE

pure function get_timelimit_string() result(string)

!Local variables-------------------------------
!scalars
 real(dp) :: timelimit
 character(len=500) :: string

! *************************************************************************

 ! Handle negative values
 timelimit = get_timelimit()
 if (timelimit > zero) then
   string = sec2str(timelimit)
 else
   string = "0"
 end if

end function get_timelimit_string
!!***

!!****f* m_exit/get_start_time
!! NAME
!! get_start_time
!!
!! FUNCTION
!!  Return the origin of execution time in seconds
!!
!! PARENTS
!!
!! SOURCE

real(dp) pure function get_start_time()

 get_start_time = WALL0

end function get_start_time
!!***

!!****f* m_exit/exit_check
!! NAME
!! exit_check
!!
!! FUNCTION
!! This routine checks whether the CPU time limit is exceeded or not.
!! If openexit is non-zero, it also checks the "filename" file
!! for the "exit" character string in its first line and returns the location
!! of the string on the line (0 if not found).  Maps both strings to upper case
!! before attempting to match them. Also checks for the existence
!! of the "abinit.exit" file in the directory where the job was started.
!! Finally, checks whether the CPU time limit was not exceeded.
!! If one of these conditions occurs, will induce graceful exit of iterations.
!!
!! INPUTS
!!  cpus = CPU time limit
!!  filename = character string giving name of file to be opened
!!  iout = unit number to print output to
!!  openexit = if 1, open the "filename" and "abinit.exit" files
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  iexit = index of "exit" on first line of file (0 if not found),
!!      or -1 if the exit was ordered through the existence of the "exit" file
!!      or -2 if the exit was ordered through the CPU time limit.
!!
!! PARENTS
!!      m_common,m_dfpt_looppert,m_driver,m_gstate,m_respfn_driver
!!
!! CHILDREN
!!      inupper,timein,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine exit_check(cpus,filename,iexit,iout,comm,openexit)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: cpus
 character(len=*),intent(in) :: filename
 integer,intent(in) :: openexit,iout
 integer,intent(out) :: iexit

!Local variables-------------------------------
!scalars
 integer,parameter :: master=0
 integer,save :: iexit_save=0
 integer :: ierr,temp_unit,ierrmpi
 logical :: ex
 real(dp),save :: tcpu_last=zero
 character(len=500) :: message
 character(len=fnlen) :: line
 character(len=4), parameter :: string='EXIT'
!arrays
 real(dp) :: tsec(2)

! *************************************************************************

 if (iexit_save==0) then
   ! ABINIT will pass again in this routine even after exit call has been detected

   if (xmpi_comm_rank(comm)==master) then
     ! Master tests and broadcast the result to others
     iexit=0

     ! Is it worth to test the cpu time ?
     tsec = zero
     if (abs(cpus)>1.0d-5 .or. openexit==1) then
       call timein(tsec(1),tsec(2))
     end if

     ! A first way of exiting: the cpu time limit
     if (abs(cpus)>1.0d-5) then
       if(cpus<tsec(1))iexit=-2
     end if

     ! Test the content of files only when sufficient time (2 sec) has elapsed from last time it was tested.
     if (openexit==1 .and. iexit==0 .and. tsec(1)-tcpu_last>two ) then
       ! TODO Remove this approach. Use abinit.exit!
       tcpu_last=tsec(1)
       ! Open file and read first line as character string
       if (open_file(filename,message,newunit=temp_unit,form='formatted',status='old') /= 0) then
         ABI_ERROR(message)
       end if
       rewind (unit=temp_unit)
       read (unit=temp_unit,fmt='(a)',iostat=ierr) line
       if(ierr/=0)then
         write(message, '(a,a,a,i5,a,a)' )&
&         'Problem when reading file=',TRIM(filename),'iostat =',ierr,ch10,&
&         'Action: check whether this file is OK.'
         ABI_ERROR(message)
       end if
       ! Make a local copy of matching string of length equal to nonblank length of input string
       ! Map to upper case
       call inupper(line)
       iexit=index(line,string)
       close (unit=temp_unit)

       ! This is another way of exiting : the previous one does not work
       ! on some machines, may be because they keep a copy of the initial input file.
       if(iexit==0)then
         inquire(file='abinit.exit',exist=ex)
         if(ex)iexit=-1
       end if

     end if
   end if

   call xmpi_bcast(iexit,master,comm,ierrmpi)

 else
   ! In case the exit mechanism has already been activated
   iexit=iexit_save
 end if

 if (iexit/=0) then
   if (iexit>0) write(message, '(a,a,a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit has been requested from file ',trim(filename),'.',ch10
   if (iexit==-1) write(message, '(a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit has been requested from file "abinit.exit".',ch10
   if (iexit==-2) write(message, '(a,a,a,a,a)' ) ch10,&
&   ' chkexi: WARNING -',ch10,&
&   '  Exit due to cpu time limit exceeded.',ch10
   if (iout/=std_out) then
     call wrtout(iout,message,'COLL')
   end if
   call wrtout(std_out,  message,'COLL')
 end if

 iexit_save=iexit

end subroutine exit_check
!!***

END MODULE m_exit
!!***

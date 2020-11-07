!!****m* ABINIT/m_specialmsg
!! NAME
!!  m_specialmsg
!!
!! FUNCTION
!!  This module contains tools to deal with special messages counters.
!!  Special messages= WARNING, COMMENT, EXIT
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MT)
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

module m_specialmsg

 use defs_basis
 use m_build_info
 use m_xmpi

 use m_io_tools,   only : flush_unit, write_lines, is_open

 implicit none

 private
!!***

 public :: herald        ! Prints message to unit iout giving info about current
                         ! code, version of code, platform, and starting date.

!Number of WARNINGS/COMMENTS printed in log file
 integer,save :: COMMENT_COUNT = 0
 integer,save :: WARNING_COUNT = 0
 integer,save :: EXIT_FLAG = 0

!Public procedures
 public :: specialmsg_setcount ! Update number of special messages (WARNING/COMMENT) present in log file
 public :: specialmsg_getcount ! Get number of special messages (WARNING/COMMENT) present in log file
 public :: specialmsg_mpisum   ! Reduce number of special messages (WARNING/COMMENT) over MPI comm

 public :: wrtout

 interface wrtout
   module procedure wrtout_unit
   module procedure wrtout_units
 end interface wrtout

CONTAINS  !===========================================================
!!***

!----------------------------------------------------------------------

!!****f* m_specialmsg/specialmsg_setcount
!! NAME
!!  specialmsg_setcount
!!
!! FUNCTION
!!  Update the counters of special messages (WARNING, COMMENTS, EXIT) printed in log file
!!
!! INPUTS
!!  [n_add_comment]= (optional) number of comments to add to the counter
!!  [n_add_exit]   = (optional) number of exit messages to add to the counter
!!  [n_add_warning]= (optional) number of warnings to add to the counter
!!
!! OUTPUT
!!  (only counters updated)
!!
!! PARENTS
!!      m_specialmsg
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine specialmsg_setcount(n_add_comment,n_add_warning,n_add_exit)

!Arguments ------------------------------------
 integer,optional,intent(in) :: n_add_comment,n_add_warning,n_add_exit

! *********************************************************************

 if (PRESENT(n_add_comment)) COMMENT_COUNT=COMMENT_COUNT+n_add_comment
 if (PRESENT(n_add_warning)) WARNING_COUNT=WARNING_COUNT+n_add_warning
 if (PRESENT(n_add_exit)) then
   EXIT_FLAG=EXIT_FLAG+n_add_exit
   if (EXIT_FLAG>1) EXIT_FLAG=1
 end if

end subroutine specialmsg_setcount
!!***

!----------------------------------------------------------------------

!!****f* m_specialmsg/specialmsg_getcount
!! NAME
!!  specialmsg_getcount
!!
!! FUNCTION
!!  Get the values of the counters of special messages (WARNING, COMMENT)
!!
!! INPUTS
!!
!! OUTPUT
!!  ncomment= number of COMMENTs in log file
!!  nwarning= number of WARNINGs in log file
!!  nexit= 1 if exit requested
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine specialmsg_getcount(ncomment,nwarning,nexit)

!Arguments ------------------------------------
 integer,intent(out) :: ncomment,nexit,nwarning

! *********************************************************************

 ncomment=COMMENT_COUNT
 nwarning=WARNING_COUNT
 nexit   =EXIT_FLAG

end subroutine specialmsg_getcount
!!***

!----------------------------------------------------------------------

!!****f* m_specialmsg/specialmsg_mpisum
!! NAME
!!  specialmsg_mpisum
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
!! PARENTS
!!      m_gstateimg
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine specialmsg_mpisum(mpicomm)

!Arguments ------------------------------------
 integer,intent(in) :: mpicomm

!Local variables-------------------------------
 integer :: ierr
 integer :: buf(3)

! *********************************************************************

  buf(1)=COMMENT_COUNT;buf(2)=WARNING_COUNT;buf(3)=EXIT_FLAG

  call xmpi_sum(buf,mpicomm,ierr)

  COMMENT_COUNT=buf(1)
  WARNING_COUNT=buf(2)
  EXIT_FLAG=buf(3) ; if (EXIT_FLAG/=0) EXIT_FLAG=1

end subroutine specialmsg_mpisum
!!***

!----------------------------------------------------------------------

!!****f* m_specialmsg/herald
!! NAME
!!  herald
!!
!! FUNCTION
!!  Prints out a message to unit iout giving info about current
!!  code, version of code, platform, and starting date.
!!
!! INPUTS
!!  code_name= code name
!!  code_version= code version
!!  iout=unit number for output
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      abinit,aim,anaddb,cut3d,fftprof,ioprof,lapackprof,mrgddb,mrgdv,mrggkk
!!      mrgscr,multibinit,optic,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine herald(code_name,code_version,iout)

!Arguments ------------------------------------
 integer,intent(in) :: iout
 character(len=*),intent(in) :: code_name
 character(len=*),intent(in) :: code_version

!Local variables-------------------------------
 integer :: day,dd,ja,jy,jm,jdn,mm,mm_rel,year,year_rel
 integer :: values(8)
 character(len=5) :: strzone
 character(len=8) :: strdat
 character(len=10) :: strtime
 character(len=500) :: msg
 character(len=3),parameter :: day_names(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
                                                 'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!RELEASE TIME FROM ABIRULES
 year_rel=2020
 mm_rel=10
!END OF RELEASE TIME

!The technique used hereafter is the only one that we have found to obtain
!perfect transferability across platforms and OS.
 write(iout, '(/,a,a,a,a,a)' ) '.Version ',trim(code_version),' of ',trim(code_name),' '
#if defined HAVE_MPI
 write(iout, '(a,a,a,/)' ) '.(MPI version, prepared for a ',build_target,' computer) '
#else
 write(iout, '(a,a,a,/)' ) '.(sequential version, prepared for a ',build_target,' computer) '
#endif

!GNU GPL license
 write(iout, '(a,/,a,a,a,/,a,/,a,/,a,/)' ) &
 '.Copyright (C) 1998-2020 ABINIT group . ',&
 ' ',trim(code_name),' comes with ABSOLUTELY NO WARRANTY.',&
 ' It is free software, and you are welcome to redistribute it',&
 ' under certain conditions (GNU General Public License,',&
 ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'

 if(trim(code_name)=='OPTIC')then
   write(iout, '(a,a,a,/,a,/,a,/,a,/,a,/,a,/,a,/)' ) &
   ' ',trim(code_name),' has originally been developed by',&
   ' Sangeeta Sharma and incorporated in ABINIT with the help of M. Verstraete.',&
   ' Please refer to : ',&
   ' S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B 67, 165332 (2003), and',&
   ' S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 (2004).',&
   '- URLs and DOI at https://docs.abinit.org/theory/bibliography/#sharma2003',&
   '- and https://docs.abinit.org/theory/bibliography/#sharma2004'
 end if

 write(iout, '(a,/,a,/,a,/,a,/,a)' ) &
 ' ABINIT is a project of the Universite Catholique de Louvain,',&
 ' Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt .',&
 ' Please read https://docs.abinit.org/theory/acknowledgments for suggested',&
 ' acknowledgments of the ABINIT effort.',&
 ' For more information, see https://www.abinit.org .'

!Get year, month and day
 call date_and_time(strdat,strtime,strzone,values)
 year=values(1)
 mm=values(2)
 dd=values(3)

!Get day of the week
 if (mm.gt.2) then
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

!Print date in nice format (* new format *)
 write(iout, '(/,a,a,1x,i2,1x,a,1x,i4,a,/,a,i2.2,a,i2.2,a)' ) &
 '.Starting date : ',day_names(day),dd,month_names(mm),year,'.','- ( at ',values(5),'h',values(6),' )'
 write(iout,*)' '

!Impose a maximal life cycle of 3 years
 if(year>year_rel+3 .or. (year==year_rel+3 .and. mm>mm_rel) ) then
   write(msg, '(5a,i4,5a)' )&
   '- The starting date is more than 3 years after the initial release',ch10,&
   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
   '- This version of ABINIT is not supported anymore.',ch10,&
   '- Action: please, switch to a more recent version of ABINIT.'
   call wrtout(iout,msg,'COLL')

!  Gives a warning beyond 2 years
 else if(year>year_rel+2 .or. (year==year_rel+2 .and. mm>mm_rel) ) then
   write(msg, '(5a,i4,6a)' )&
   '- The starting date is more than 2 years after the initial release',ch10,&
   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
   '- Note that the use beyond 3 years after the release will not be supported.',ch10,&
   '- Action: please, switch to a more recent version of ABINIT.',ch10
   call wrtout(iout,msg,'COLL')
 end if

end subroutine herald
!!***

!!****f* m_specialmsg/wrtout_unit
!! NAME
!!  wrtout_unit
!!
!! FUNCTION
!!  Organizes the sequential or parallel version of the write intrinsic
!!  Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  unit=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!  [newlines]: Number of new lines added after message. Default 0
!!  [pre_newlines]: Number of new lines added vefore message. Default 0
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_specialmsg
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine wrtout_unit(unit, msg, mode_paral, do_flush, newlines, pre_newlines)

!Arguments ------------------------------------
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: do_flush
 integer,optional,intent(in) :: newlines, pre_newlines

!Local variables-------------------------------
 integer,save :: master=0
 integer :: comm, me, nproc, my_newlines, ii,  my_pre_newlines
 logical :: my_flush
 character(len=len(msg)+50) :: string
 character(len=500) :: my_mode_paral

!******************************************************************

 if (unit == std_out .and. .not. do_write_log) return
 if (unit == dev_null) return
 !if (.not. is_open(unit)) return

 my_mode_paral = "COLL"; if (present(mode_paral)) my_mode_paral = mode_paral
 my_flush = .false.; if (present(do_flush)) my_flush = do_flush
 my_newlines = 0; if (present(newlines)) my_newlines = newlines
 my_pre_newlines = 0; if (present(pre_newlines)) my_pre_newlines = pre_newlines

 ! Communicator is xmpi_world by default, except for the parallelization over images
 if (abinit_comm_output /= -1) then
   comm = abinit_comm_output
 else
   comm = xmpi_world
 end if

 ! Determine who I am in COMM_WORLD
 me = xmpi_comm_rank(comm); nproc = xmpi_comm_size(comm)

 if (my_mode_paral == 'COLL' .or. nproc == 1) then
   if (me == master) then
     if (my_pre_newlines /= 0) then
       do ii=1,my_pre_newlines; write(unit, "(a)")""; end do
     end if
     call wrtout_myproc(unit, msg, do_flush=my_flush)
     if (my_newlines /= 0) then
       do ii=1,my_newlines; write(unit, "(a)")""; end do
     end if
   end if

 else if (my_mode_paral == 'PERS') then
   if (my_pre_newlines /= 0) then
     do ii=1,my_pre_newlines; write(unit, "(a)")""; end do
   end if
   call write_lines(unit,msg)
   if (my_newlines /= 0) then
     do ii=1,my_newlines; write(unit, "(a)")""; end do
   end if
   ! Flush unit
   if (my_flush) call flush_unit(unit)

 else if (my_mode_paral == 'INIT') then
   master = unit

 else
   write(string,'(7a)')ch10,&
   'wrtout_unit: ERROR -',ch10,&
   '  Unknown write mode: ',trim(my_mode_paral),ch10,&
   '  Continuing anyway ...'
   write(unit, '(A)' ) trim(string)
 end if

end subroutine wrtout_unit
!!***

!!****f* m_specialmsg/wrtout_units
!! NAME
!!  wrtout_units
!!
!! FUNCTION
!!  Write string to multiple units. Wraps wrtout_unit
!!
!! INPUTS
!!  msg=(character(len=*)) message to be written
!!  units=unit number for writing. The named constant dev_null defined in defs_basis can be used to avoid any printing.
!!  [mode_paral]= --optional argument--
!!   'COLL' if all procs are calling the routine with the same message to be written once only. Default.
!!   'PERS' if the procs are calling the routine with different messages each to be written,
!!          or if one proc is calling the routine
!!   "INIT" to change the rank of the master node that prints the message if "COLL" is used.
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!  [newlines]: Number of new lines added after message. Default 0
!!  [pre_newlines]: Number of new lines added vefore message. Default 0
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine wrtout_units(units, msg, mode_paral, do_flush, newlines, pre_newlines)

!Arguments ------------------------------------
 integer,intent(in) :: units(:)
 character(len=*),intent(in) :: msg
 character(len=*),optional,intent(in) :: mode_paral
 logical,optional,intent(in) :: do_flush
 integer,optional,intent(in) :: newlines, pre_newlines

!Local variables-------------------------------
!scalars
 integer :: ii, cnt, my_newlines, my_pre_newlines
 logical :: my_flush
 character(len=500) :: my_mode_paral
!arrays
 integer :: my_units(size(units))

!******************************************************************

 my_mode_paral = "COLL"; if (present(mode_paral)) my_mode_paral = mode_paral
 my_flush = .false.; if (present(do_flush)) my_flush = do_flush
 my_newlines = 0; if (present(newlines)) my_newlines = newlines
 my_pre_newlines = 0; if (present(pre_newlines)) my_pre_newlines = pre_newlines

 ! Remove duplicated units (if any)
 my_units(1) = units(1); cnt = 1
 do ii=2,size(units)
   if (any(units(ii) == my_units(1:cnt))) cycle
   cnt = cnt + 1
   my_units(cnt) = units(ii)
 end do

 do ii=1,cnt
   call wrtout_unit(my_units(ii), msg, mode_paral=my_mode_paral, &
                    do_flush=my_flush, newlines=my_newlines, pre_newlines=my_pre_newlines)
 end do

end subroutine wrtout_units
!!***

!----------------------------------------------------------------------

!!****f* m_specialmsg/wrtout_myproc
!! NAME
!!  wrtout_myproc
!!
!! FUNCTION
!!  Do the output for one proc. For parallel or sequential output use wrtout()
!!  instead. Also allows to treat correctly the write operations for Unix (+DOS) and MacOS.
!!
!! INPUTS
!!  unit=unit number for writing
!!  msg=(character(len=*)) message to be written
!!  [do_flush]=True to flush the unit. Defaults to .False.
!!
!! OUTPUT
!!  (only writing)
!!
!! PARENTS
!!      m_specialmsg
!!
!! CHILDREN
!!      flush_unit,specialmsg_setcount,write_lines
!!
!! SOURCE

subroutine wrtout_myproc(unit, msg, do_flush) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: unit
 character(len=*),intent(in) :: msg
 logical,optional,intent(in) :: do_flush

!Local variables-------------------------------
!scalars
 logical :: print_std_err

!******************************************************************

 print_std_err = (unit == std_out .and. std_out /= std_err .and. &
   (index(trim(msg), 'BUG') /= 0 .or. index(trim(msg), 'ERROR') /= 0))

 ! Print message
 call write_lines(unit, msg)
 if (print_std_err) call write_lines(std_err, msg)

 ! Append "Contact Abinit group" to BUG messages
 if (index(trim(msg), 'BUG') /= 0 )then
   write(unit, '(a)' ) '  Action: contact ABINIT group (please attach the output of `abinit -b`)'
   write(unit,*)
   if (print_std_err) then
     write(std_err, '(a)' ) '  Action: contact ABINIT group (please attach the output of `abinit -b`)'
     write(std_err,*)
   end if
 end if

 ! Count the number of warnings and comments. Only take into
 ! account unit std_out, in order not to duplicate these numbers.
 if (index(trim(msg), 'WARNING') /= 0 .and. unit==std_out) call specialmsg_setcount(n_add_warning=1)
 if (index(trim(msg), 'COMMENT') /= 0 .and. unit==std_out) call specialmsg_setcount(n_add_comment=1)
 if (index(trim(msg), 'Exit') /= 0 ) call specialmsg_setcount(n_add_exit=1)

 ! Flush unit
 if (present(do_flush)) then
   if (do_flush) call flush_unit(unit)
 end if
#ifdef DEBUG_MODE
 call flush_unit(unit)
 if (print_std_err) call flush_unit(std_err)
#endif

end subroutine wrtout_myproc
!!***

end module m_specialmsg
!!***

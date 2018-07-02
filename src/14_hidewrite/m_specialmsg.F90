!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_specialmsg
!! NAME
!!  m_specialmsg
!!
!! FUNCTION
!!  This module contains tools to deal with special messages counters.
!!  Special messages= WARNING, COMMENT, EXIT
!!
!! COPYRIGHT
!! Copyright (C) 2008-2018 ABINIT group (MT)
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

MODULE m_specialmsg

 use defs_basis
 use m_build_info
 use m_xmpi

 implicit none

 private
!!***

 public :: herald        ! Prints message to unit iout giving info about current
                         ! code, version of code, platform, and starting date.

!Number of WARNINGS/COMMENTS printed in log file
 integer,save :: COMMENT_COUNT=0
 integer,save :: WARNING_COUNT=0
 integer,save :: EXIT_FLAG=0

!Public procedures
 public :: specialmsg_setcount ! Update number of special messages (WARNING/COMMENT) present in log file
 public :: specialmsg_getcount ! Get number of special messages (WARNING/COMMENT) present in log file
 public :: specialmsg_mpisum   ! Reduce number of special messages (WARNING/COMMENT) over MPI comm

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
!!      wrtout
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine specialmsg_setcount(n_add_comment,n_add_warning,n_add_exit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'specialmsg_setcount'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,optional,intent(in) :: n_add_comment,n_add_warning,n_add_exit

!Local variables-------------------------------

! *********************************************************************

 if (PRESENT(n_add_comment)) COMMENT_COUNT=COMMENT_COUNT+n_add_comment
 if (PRESENT(n_add_warning)) WARNING_COUNT=WARNING_COUNT+n_add_warning
 if (PRESENT(n_add_exit   )) then
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
!!  nexit=    1 if exit requested
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine specialmsg_getcount(ncomment,nwarning,nexit)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'specialmsg_getcount'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 integer,intent(out) :: ncomment,nexit,nwarning

!Local variables-------------------------------

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
!!      gstateimg
!!
!! CHILDREN
!!      xmpi_sum
!!
!! SOURCE

subroutine specialmsg_mpisum(mpicomm)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'specialmsg_mpisum'
!End of the abilint section

 implicit none

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
!!      abinit,aim,anaddb,bsepostproc,cut3d,fftprof,ioprof,lapackprof,mrgddb
!!      mrgdv,mrggkk,mrgscr,multibinit,optic,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!      date_and_time,wrtout
!!
!! SOURCE

subroutine herald(code_name,code_version,iout)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'herald'
 use interfaces_14_hidewrite
!End of the abilint section

 implicit none

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
 character(len=500) :: message
 character(len=3),parameter :: day_names(7)=(/'Mon','Tue','Wed','Thu','Fri','Sat','Sun'/)
 character(len=3),parameter :: month_names(12)=(/'Jan','Feb','Mar','Apr','May','Jun',&
&                                                'Jul','Aug','Sep','Oct','Nov','Dec'/)

! *************************************************************************

!RELEASE TIME FROM ABIRULES
 year_rel=2018
 mm_rel=6
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
& '.Copyright (C) 1998-2018 ABINIT group . ',&
& ' ',trim(code_name),' comes with ABSOLUTELY NO WARRANTY.',&
& ' It is free software, and you are welcome to redistribute it',&
& ' under certain conditions (GNU General Public License,',&
& ' see ~abinit/COPYING or http://www.gnu.org/copyleft/gpl.txt).'

 if(trim(code_name)=='OPTIC')then
   write(iout, '(a,a,a,/,a,/,a,/,a,/,a,/)' ) &
&   ' ',trim(code_name),' has originally been developed by',&
&   ' Sangeeta Sharma and incorporated in ABINIT with the help of M. Verstraete.',&
&   ' Please refer to : ',&
&   ' S. Sharma, J. K. Dewhurst and C. Ambrosch-Draxl, Phys. Rev. B 67, 165332 (2003), and',&
&   ' S. Sharma and C. Ambrosch-Draxl, Physica Scripta T 109 (2004).'
 end if

 write(iout, '(a,/,a,/,a,/,a,/,a)' ) &
& ' ABINIT is a project of the Universite Catholique de Louvain,',&
& ' Corning Inc. and other collaborators, see ~abinit/doc/developers/contributors.txt .',&
& ' Please read ~abinit/doc/users/acknowledgments.html for suggested',&
& ' acknowledgments of the ABINIT effort.',&
& ' For more information, see https://www.abinit.org .'

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
& '.Starting date : ',day_names(day),dd,month_names(mm),year,'.','- ( at ',values(5),'h',values(6),' )'
 write(iout,*)' '

!Impose a maximal life cycle of 3 years
 if(year>year_rel+3 .or. (year==year_rel+3 .and. mm>mm_rel) ) then
   write(message, '(5a,i4,5a)' )&
&   '- The starting date is more than 3 years after the initial release',ch10,&
&   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
&   '- This version of ABINIT is not supported anymore.',ch10,&
&   '- Action: please, switch to a more recent version of ABINIT.'
   call wrtout(iout,message,'COLL')

!  Gives a warning beyond 2 years
 else if(year>year_rel+2 .or. (year==year_rel+2 .and. mm>mm_rel) ) then
   write(message, '(5a,i4,6a)' )&
&   '- The starting date is more than 2 years after the initial release',ch10,&
&   '- of this version of ABINIT, namely ',month_names(mm_rel),' ',year_rel,'.',ch10,&
&   '- Note that the use beyond 3 years after the release will not be supported.',ch10,&
&   '- Action: please, switch to a more recent version of ABINIT.',ch10
   call wrtout(iout,message,'COLL')
 end if

end subroutine herald
!!***

END MODULE m_specialmsg
!!***

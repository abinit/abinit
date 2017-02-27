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
!! Copyright (C) 2008-2017 ABINIT group (MT)
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

 use m_xmpi, only : xmpi_sum

 implicit none
 
 private
!!***

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

END MODULE m_specialmsg
!!***

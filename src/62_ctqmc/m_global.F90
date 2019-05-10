
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_global
!! NAME
!!  m_global
!!
!! FUNCTION
!!  Manage error and warnings for the ctqmc
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

#include "defs.h"

MODULE m_Global

#if defined HAVE_CONFIG_H
! we are in abinit
USE defs_basis
USE m_profiling_abi
USE m_errors
USE m_xmpi
#endif
#ifdef HAVE_MPI2
USE mpi
#endif

IMPLICIT NONE

PUBLIC

PUBLIC :: ERROR
PUBLIC :: WARN
PUBLIC :: WARNALL

CONTAINS
!!***

!!****f* ABINIT/m_global/ERROR
!! NAME
!!  ERROR
!!
!! FUNCTION
!!  error dectected => leave
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  message=error message to display
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE ERROR(message)

!Arguments ------------------------------------
#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  CHARACTER(LEN=*), INTENT(IN) :: message
!Local variables ------------------------------
  CHARACTER(LEN=500)            :: messend
#ifdef HAVE_MPI
  INTEGER                       :: ierr
  INTEGER                       :: rank
  CALL MPI_Comm_rank(MY_WORLD, rank, ierr)
  WRITE(messend,'(A,i5,A,A)') "ERROR in QMC rank ", rank, " : ",TRIM(message)
  myERROR(TRIM(messend))
  CALL MPI_Finalize(ierr) ! IF in abinit, does nothing since killed in _myERROR_
#else
  WRITE(messend,'(A,A)') "ERROR in QMC : ", TRIM(message)
  myERROR(TRIM(messend))
#endif
  !CALL FLUSH(0)
  STOP
END SUBROUTINE ERROR
!!***

!!****f* ABINIT/m_global/WARN
!! NAME
!!  WARN
!!
!! FUNCTION
!!  on cpu wants to tell something
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  message=warning message
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE WARN(message)

!Arguments ------------------------------------
#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  CHARACTER(LEN=*), INTENT(IN) :: message
!Local variables ------------------------------
  CHARACTER(LEN=500)           :: messend
#ifdef HAVE_MPI
  INTEGER                      :: ierr
  INTEGER                      :: rank
  CALL MPI_Comm_rank(MY_WORLD, rank, ierr)
  WRITE(messend,'(A,I6,A)') "WARNING in QMC rank ", rank, " : ", TRIM(message)
#else
  WRITE(messend,'(A,A)') "WARNING in QMC : ", TRIM(message)
#endif
  myWARN(TRIM(messend))
  !CALL FLUSH(0)

END SUBROUTINE WARN
!!***

!!****f* ABINIT/m_global/WARNALL
!! NAME
!!  WARNALL
!!
!! FUNCTION
!!  collective warning function
!!
!! COPYRIGHT
!!  Copyright (C) 2013-2019 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  message=message to display
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! NOTES
!!
!! PARENTS
!!  Will be filled automatically by the parent script
!!
!! CHILDREN
!!  Will be filled automatically by the parent script
!!
!! SOURCE

SUBROUTINE WARNALL(message)

!Arguments ------------------------------------
#ifdef HAVE_MPI1
include 'mpif.h'
#endif
  CHARACTER(LEN=*), INTENT(IN) :: message
!Local variables ------------------------------------
  CHARACTER(LEN=500)           :: messend
#ifdef HAVE_MPI
  INTEGER                      :: ierr
  INTEGER                      :: rank
  CALL MPI_Comm_rank(MY_WORLD, rank, ierr)
  IF ( rank .EQ. 0) THEN
#endif
    WRITE(messend,'(A,A)') "WARNING in QMC : ", TRIM(message)
    myWARNALL(TRIM(messend))
#ifdef HAVE_MPI
  END IF
#endif
  !CALL FLUSH(0)

END SUBROUTINE WARNALL
!!***

END MODULE m_global
!!***

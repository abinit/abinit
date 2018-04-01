!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_time
!! NAME
!! defs_time
!!
!! FUNCTION
!! This module contains accumulators for the timer.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2018 ABINIT group (TD)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! TODO
!! Include the name of all routines: better modularity
!!
!! PARENTS
!!    timab,time_accu
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module defs_time

 use defs_basis
 use iso_c_binding

 implicit none

 private

!mtim determines the maximum number of "timing slots" available
 integer,public,parameter :: mtim=1999

! timeopt is a flag which indicates the suppression or not of the timing.
 integer,public,save :: timopt=1

! Number of times that the routine has been called
 integer,public,save :: ncount(mtim)=0

! Accumulating cpu time (1) and wall to wall time (2) for each "timing slots"
 real(dp),public,save  :: acctim(2,mtim)=zero,tzero(2,mtim)=zero

! Accumulating number of floating point operation and cpu time (1) and wall time (2) for each "performance slot"
 real(dp),public,save :: papi_accflops(mtim)=zero, papi_acctim(2,mtim)=zero

! Reference value for number of floating point operation and time (cpu and wall) for each performance slot
 real(dp),public,save :: papi_flops(mtim)=zero , papi_tzero(2,mtim)=zero

! Elapsed time and elapsed number of floating point operation since a reference
 real(dp),public,save :: papi_tottim(2,mtim)=zero, papi_totflops(mtim)=zero

 public :: time_set_papiopt
 public :: time_get_papiopt

! papiopt is a flag which indicates if there is or not an analysis of speed execution is made.
! By defaut the analysis is not done
 integer,private,save :: papiopt=0

CONTAINS  !=========================================================================================================================
!!***

!----------------------------------------------------------------------

!!****f* defs_time/time_set_papiopt
!! NAME
!!  time_set_papiopt
!!
!! FUNCTION
!!  Set the value of papiopt
!!
!! PARENTS
!!      abinit
!!
!! CHILDREN
!!
!! SOURCE

subroutine time_set_papiopt(opt)

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_set_papiopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: opt

! *************************************************************************

 papiopt = opt

end subroutine time_set_papiopt
!!***

!----------------------------------------------------------------------

!!****f* defs_time/time_get_papiopt
!! NAME
!!  time_get_papiopt
!!
!! FUNCTION
!!  Return the value of papiopt
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

function time_get_papiopt()

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'time_get_papiopt'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer :: time_get_papiopt

! *************************************************************************

 time_get_papiopt = papiopt

end function time_get_papiopt
!!***

!----------------------------------------------------------------------

end module defs_time
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/timein
!! NAME
!!  timein
!!
!! FUNCTION
!!  Timing routine. Returns cpu and wall clock time in seconds since some arbitrary start.
!!  For wall clock time, call the F90 intrinsic date_and_time .
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, LSI, MM, MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!  For the initials of contributors, see
!!  ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  (no inputs)
!!
!! OUTPUT
!!  cpu= cpu time in seconds
!!  wall= wall clock time in seconds
!!
!! NOTES
!!  For CPU time, contains machine-dependent code (choice will be selected
!!  by C preprocessor, see abi_cpu_time).
!!
!! TODO
!!  Should be replaced by cwtime
!!
!! PARENTS
!!      abinit,aim,aim_follow,anaddb,bsepostproc,conducti,cpdrv,cut3d,drvaim
!!      elphon,first_rec,m_exit,mrgddb,mrgscr,multibinit,optic,rsurf,surf,timab
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine timein(cpu,wall)

 use defs_basis
 use m_time

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'timein'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 real(dp),intent(out) :: cpu,wall
! *************************************************************************

!CPU time
 cpu = abi_cpu_time()

!Wall time
 wall = abi_wtime()

end subroutine timein
!!***

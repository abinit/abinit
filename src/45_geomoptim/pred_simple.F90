!{\src2tex{textfont=tt}}
!!****f* ABINIT/pred_simple
!! NAME
!! pred_simple
!!
!! FUNCTION
!! Ionmov predictors (4 & 6) Internal to SCFV
!!
!! IONMOV 4:
!! Conjugate gradient algorithm for simultaneous optimization
!! of potential and ionic degrees of freedom. It can be used with
!! iscf=2 and iscf=5 or 6
!!
!! IONMOV 5:
!! Simple relaxation of ionic positions according to (converged)
!! forces. Equivalent to ionmov=1 with zero masses, albeit the
!! relaxation coefficient is not vis, but iprcfc.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, JCC, SE)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors,
!! see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! ab_mover <type(abimover)> : Datatype with all the information
!!                                needed by the preditor
!! zDEBUG : if true print some debugging information
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!! hist <type(abihist)> : History of positions,forces
!!                               acell, rprimd, stresses
!! NOTES
!!
!! PARENTS
!!      mover
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pred_simple(ab_mover,hist,iexit)

 use defs_basis
 use m_profiling_abi
 use m_abimover
 use m_abihist

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pred_simple'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 type(abimover),intent(in) :: ab_mover
 type(abihist),intent(inout) :: hist
 integer,intent(in) :: iexit

!Local variables-------------------------------
!scalars
 integer  :: kk,jj

!***************************************************************************
!Beginning of executable session
!***************************************************************************

 if(iexit/=0)then
   return
 end if

!All the operations are internal to scfcv.F90

!XRED, FCART and VEL
 do kk=1,ab_mover%natom
   do jj=1,3
     hist%histXF(jj,kk,1,hist%ihist+1)=hist%histXF(jj,kk,1,hist%ihist)
     hist%histXF(jj,kk,2,hist%ihist+1)=hist%histXF(jj,kk,2,hist%ihist)
     hist%histV (jj,kk,  hist%ihist+1)=hist%histV (jj,kk,  hist%ihist)
   end do ! jj=1,3
 end do ! kk=1,ab_mover%natom

!ACELL
!For the 3 dimensions of space
 do jj=1,3
   hist%histA(jj,hist%ihist+1)=hist%histA(jj,hist%ihist)
 end do

!RPRIMD
!For the 3 dimensions of space
 do kk=1,3
!  For the 3 dimensions of space
   do jj=1,3
     hist%histR(jj,kk,hist%ihist+1)=hist%histR(jj,kk,hist%ihist)
   end do
 end do

 hist%ihist=hist%ihist+1

end subroutine pred_simple
!!***

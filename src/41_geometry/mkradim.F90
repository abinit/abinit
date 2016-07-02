!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkradim
!! NAME
!! mkradim
!!
!! FUNCTION
!!  Not so trivial subroutine to make dimensionless real space
!!  primitive translations rprim(3,3) from dimentional rprimd(3).
!!  also make acell(3).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! OUTPUT
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! PARENTS
!!      gstate,ingeo,m_dvdb,m_pimd,m_use_ga,mover,pred_steepdesc,predict_pimd
!!      wvl_memory,xfpack_vin2x
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkradim(acell,rprim,rprimd)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkradim'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(out) :: acell(3),rprim(3,3)
 real(dp),intent(in) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii

! *************************************************************************

!Use a representation based on normalised rprim vectors
 do ii=1,3
   acell(ii)=sqrt(rprimd(1,ii)**2+rprimd(2,ii)**2+rprimd(3,ii)**2)
   rprim(:,ii)=rprimd(:,ii)/acell(ii)
 end do

end subroutine mkradim
!!***

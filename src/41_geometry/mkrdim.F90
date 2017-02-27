!{\src2tex{textfont=tt}}
!!****f* ABINIT/mkrdim
!! NAME
!! mkrdim
!!
!! FUNCTION
!!  Trivial subroutine to make dimensional real space
!!  primitive translations from length scales acell(3)
!!  and dimensionless translations rprim(3,3).
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  acell(3)=unit cell length scales (bohr)
!!  rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!              where: rprimd(i,j)=rprim(i,j)*acell(j)
!!
!! PARENTS
!!      bethe_salpeter,dfpt_looppert,dfpt_symph,driver,finddistrproc
!!      get_npert_rbz,harmonic_thermo,ingeo,invars1,invars2m,m_ddb,m_ifc
!!      m_results_img,m_use_ga,memory_eval,mpi_setup,outvar_o_z,pred_bfgs
!!      pred_delocint,pred_diisrelax,pred_isothermal,pred_lbfgs,pred_steepdesc
!!      pred_verlet,predict_pimd,randomcellpos,screening,setup1,setup_bse
!!      setup_screening,setup_sigma,sigma,thmeig,wvl_setboxgeometry
!!      xfpack_x2vin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine mkrdim(acell,rprim,rprimd)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'mkrdim'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 real(dp),intent(out) :: rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,jj

! *************************************************************************

 do ii=1,3
   do jj=1,3
     rprimd(ii,jj)=rprim(ii,jj)*acell(jj)
   end do
 end do

end subroutine mkrdim
!!***

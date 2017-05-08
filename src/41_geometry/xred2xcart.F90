!{\src2tex{textfont=tt}}
!!****f* ABINIT/xred2xcart
!! NAME
!! xred2xcart
!!
!! FUNCTION
!! Convert from dimensionless reduced coordinates xred(3,natom)
!! to cartesian coordinates xcart(3,natom) in bohr by using
!! xcart(mu,ia)=rprimd(mu,1)*xred(1,ia)
!!             +rprimd(mu,2)*xred(2,ia)
!!             +rprimd(mu,3)*xred(3,ia)
!! Note that the reverse operation is done by xcart2xred.F90
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! OUTPUT
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      afterscfloop,berryphase,berryphase_new,bonds_lgth_angles,constrf,cut3d
!!      denfgr,driver,evdw_wannier,forstr,ingeo,ionion_realspace,ionion_surface
!!      m_abihist,m_crystal,m_ddb,m_effective_potential,m_fit_polynomial_coeff
!!      m_mep,m_pred_lotf,m_results_img,make_efg_el,make_efg_ion,mkcore_paw
!!      mkcore_wvl,mkgrid_fft,mklocl,mklocl_realspace,mlwfovlp_projpaw
!!      mover_effpot,out1dm,outqmc,outvar_o_z,outxml,pimd_langevin_npt
!!      pimd_langevin_nvt,pimd_nosehoover_npt,pimd_nosehoover_nvt,prec_simple
!!      pred_delocint,pred_diisrelax,pred_hmc,pred_isokinetic,pred_isothermal
!!      pred_langevin,pred_moldyn,pred_nose,pred_srkna14,pred_steepdesc
!!      pred_velverlet,pred_verlet,prtimg,prtspgroup,prtxfase,randomcellpos
!!      rhotov,setvtr,spin_current,symspgr,thmeig,vso_realspace_local,vtorho
!!      wrt_moldyn_netcdf,wvl_denspot_set,wvl_initro,wvl_memory,wvl_nhatgrid
!!      wvl_projectors_set,wvl_rwwf,wvl_setboxgeometry,wvl_wfs_set
!!      wvl_wfsinp_reformat,wvl_wfsinp_scratch,xfh_recover_deloc
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine xred2xcart(natom,rprimd,xcart,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xred2xcart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: rprimd(3,3),xred(3,natom)
 real(dp),intent(out) :: xcart(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
!arrays

! *************************************************************************

 do iatom=1,natom
   do mu=1,3
     xcart(mu,iatom)=rprimd(mu,1)*xred(1,iatom)+rprimd(mu,2)*xred(2,iatom)+rprimd(mu,3)*xred(3,iatom)
   end do
 end do

end subroutine xred2xcart
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/xcart2xred
!! NAME
!! xcart2xred
!!
!! FUNCTION
!! Convert from cartesian coordinates xcart(3,natom) in bohr to
!! dimensionless reduced coordinates xred(3,natom) by using
!! xred(mu,ia)=gprimd(1,mu)*xcart(1,ia)
!!            +gprimd(2,mu)*xcart(2,ia)
!!            +gprimd(3,mu)*xcart(3,ia)
!! where gprimd is the inverse of rprimd
!! Note that the reverse operation is deon by xred2xcart
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in unit cell
!!  rprimd(3,3)=dimensional real space primitive translations (bohr)
!!  xcart(3,natom)=cartesian coordinates of atoms (bohr)
!!
!! OUTPUT
!!  xred(3,natom)=dimensionless reduced coordinates of atoms
!!
!! PARENTS
!!      driver,evdw_wannier,ingeo,m_cut3d,m_dens,m_effective_potential
!!      m_effective_potential_file,m_mep,m_paw_pwaves_lmn,m_pred_lotf
!!      mkcore_paw,mkcore_wvl,mover_effpot,pawmkaewf,pimd_langevin_npt
!!      pimd_langevin_nvt,pimd_nosehoover_npt,pimd_nosehoover_nvt,prcref
!!      prcref_PMA,pred_delocint,pred_diisrelax,pred_isokinetic,pred_isothermal
!!      pred_langevin,pred_moldyn,pred_nose,pred_srkna14,pred_steepdesc
!!      pred_velverlet,pred_verlet,relaxpol,wrt_moldyn_netcdf
!!      wvl_setboxgeometry
!!
!! CHILDREN
!!      matr3inv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine xcart2xred(natom,rprimd,xcart,xred)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xcart2xred'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom
!arrays
 real(dp),intent(in) :: rprimd(3,3),xcart(3,natom)
 real(dp),intent(out) :: xred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
!arrays
 real(dp) :: gprimd(3,3)

! *************************************************************************

 call matr3inv(rprimd,gprimd)
 do iatom=1,natom
   do mu=1,3
     xred(mu,iatom)= gprimd(1,mu)*xcart(1,iatom)+gprimd(2,mu)*xcart(2,iatom)+&
&     gprimd(3,mu)*xcart(3,iatom)
   end do
 end do

end subroutine xcart2xred
!!***

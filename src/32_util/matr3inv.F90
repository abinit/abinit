!{\src2tex{textfont=tt}}
!!****f* ABINIT/matr3inv
!! NAME
!! matr3inv
!!
!! FUNCTION
!! Invert and transpose general 3x3 matrix of real*8 elements.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!! aa = 3x3 matrix to be inverted
!!
!! OUTPUT
!! ait = inverse of aa input matrix
!!
!! NOTES
!! Returned array is TRANSPOSE of inverse, as needed to get g from r.
!!
!! PARENTS
!!      berryphase,chkdilatmx,conducti_nc,ddb_hybrid,dfpt_mkvxc,dfpt_mkvxcstr
!!      dfpt_symph,electrooptic,ep_el_weights,ep_fs_weights,ep_ph_weights
!!      get_kpt_fullbz,getghcnd,getkgrid,getspinrot,gstate,harmonic_thermo
!!      invars2,inwffil,m_cut3d,m_ddb,m_ddk,m_double_grid,m_dynmat
!!      m_effective_potential,m_esymm,m_ewald,m_fock,m_fstab,m_ifc
!!      m_phonon_supercell,m_psps,m_strain,make_efg_el,make_efg_ion,metric
!!      optic,outwant,pimd_langevin_npt,prtxf,relaxpol,smpbz,stresssym,symbrav
!!      symlatt,symmetrize_rprimd,symrelrot,symrhg,tddft,testkgrid,thmeig
!!      uderiv,xcart2xred,xfpack_x2vin
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine matr3inv(aa,ait)

 use defs_basis
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'matr3inv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,det,t1,t2,t3
 character(len=500) :: message

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det  = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3 

!Make sure matrix is not singular
 if (abs(det)>tol16) then
   dd=one/det
 else
   write(message, '(2a,2x,9es16.8,a,a,es16.8,a)' )&
&   ' Attempting to invert real(8) 3x3 array',ch10,&
&   aa(:,:),ch10,&
&   '   ==> determinant=',det,' is zero.'
   MSG_BUG(message)
 end if

 ait(1,1) = t1 * dd
 ait(2,1) = t2 * dd
 ait(3,1) = t3 * dd
 ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
 ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
 ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
 ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
 ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
 ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd

end subroutine matr3inv
!!***

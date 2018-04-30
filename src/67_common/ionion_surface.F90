!{\src2tex{textfont=tt}}
!!****f* ABINIT/ionion_surface
!!
!! NAME
!! ionion_surface
!!
!! FUNCTION
!! Compute the ion/ion interaction energies and forces in real space
!! case. Use ewald() instead if computations are done in reciprocal
!! space since it also includes the correction for the shift done in
!! potentials calculations and includes replica interactions.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!  rmet(3,3)=metric tensor in real space (bohr^2)
!!  xred(3,natom)=relative coords of atoms in unit cell (dimensionless)
!!  zion(ntypat)=charge on each type of atom (real number)
!!
!! OUTPUT
!!  eew=final ion/ion energy in hartrees
!!  grewtn(3,natom)=grads of ion/ion wrt xred(3,natom), hartrees.
!!
!! PARENTS
!!      setvtr
!!
!! CHILDREN
!!      ionicenergyandforces,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine ionion_surface(dtset, eew, grewtn, me, nproc, rprimd, wvl, wvl_den, xred)

 use defs_basis
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
#if defined HAVE_BIGDFT
 use BigDFT_API, only: IonicEnergyandForces
#endif

 use m_geometry,    only : xred2xcart

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'ionion_surface'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: me, nproc
 real(dp),intent(out) :: eew
 type(dataset_type),intent(in) :: dtset
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
!arrays
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: xred(3,dtset%natom)
 real(dp),intent(out) :: grewtn(3,dtset%natom)

!Local variables-------------------------------
!scalars
 integer :: dispersion, iatom, igeo
 real(dp) :: psoffset
!arrays
 real(dp),allocatable :: xcart(:,:)
 real(dp),pointer :: grew_cart(:,:),fdisp(:,:)
#if defined HAVE_BIGDFT
 real(dp) :: edisp
 real(dp) :: ewaldstr(6)
#endif

! *************************************************************************

!Store xcart for each atom
 ABI_ALLOCATE(xcart,(3, dtset%natom))
 call xred2xcart(dtset%natom, rprimd, xcart, xred)

 nullify(fdisp)
 nullify(grew_cart)
 dispersion = 0
 psoffset = 0._dp
#if defined HAVE_BIGDFT
 call IonicEnergyandForces(me, nproc, wvl_den%denspot%dpbox,&
& wvl%atoms, dtset%efield, xcart, &
& eew, grew_cart, dispersion, edisp, fdisp,&
& ewaldstr,wvl%Glr%d%n1,wvl%Glr%d%n2,wvl%Glr%d%n3,&
& wvl_den%denspot%V_ext, wvl_den%denspot%pkernel,psoffset)

 if (associated(fdisp)) then
   ABI_DEALLOCATE(fdisp)
 end if
#endif

 ABI_DEALLOCATE(xcart)

!Transform cartesian gradients to reduced gradients.
 do iatom = 1, dtset%natom, 1
   do igeo = 1, 3, 1
     grewtn(igeo, iatom) = -rprimd(1, igeo) * grew_cart(1, iatom) - &
&     rprimd(2, igeo) * grew_cart(2, iatom) - &
&     rprimd(3, igeo) * grew_cart(3, iatom)
   end do
 end do
 if (associated(grew_cart)) then
   ABI_DEALLOCATE(grew_cart)
 end if

#if !defined HAVE_BIGDFT
 if (.false.) write(std_out,*) me,nproc,wvl%h(1),wvl_den%symObj
#endif

end subroutine ionion_surface
!!***

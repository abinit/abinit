!{\src2tex{textfont=tt}}
!!****f* ABINIT/pspcor
!! NAME
!! pspcor
!!
!! FUNCTION
!! Compute ecore pseudoion-pseudoion correction energy from epsatm for
!! different types of atoms in unit cell.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  ntypat=number of types of atoms
!!  typat(natom)=integer label of 'typat' for each atom in cell
!!  epsatm(ntypat)=pseudoatom energy for each type of atom
!!  zion(ntypat)=valence charge on each type of atom in cell
!!
!! OUTPUT
!!  ecore=resulting psion-psion energy in Hartrees
!!
!! PARENTS
!!      pspini
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine pspcor(ecore,epsatm,natom,ntypat,typat,zion)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'pspcor'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,ntypat
 real(dp),intent(out) :: ecore
!arrays
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: epsatm(ntypat),zion(ntypat)

!Local variables-------------------------------
!scalars
 integer :: ia
 real(dp) :: charge,esum

! *************************************************************************

 charge = 0.d0
 esum = 0.d0
 do ia=1,natom
!  compute pseudocharge:
   charge=charge+zion(typat(ia))
!  add pseudocore energies together:
   esum = esum + epsatm(typat(ia))
 end do

 ecore=charge*esum

end subroutine pspcor
!!***

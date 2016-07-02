!{\src2tex{textfont=tt}}
!!****f* ABINIT/fred2fcart
!!
!! NAME
!! fred2fcart
!!
!! FUNCTION
!! Convert reduced forces into cartesian forces
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR, FJ, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  natom=Number of atoms in the unitary cell
!!  jellslab=
!!  gprimd(3,3)=
!!
!! OUTPUT
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!
!! SIDE EFFECTS
!!
!! NOTES
!!    Unlike fred, fcart has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! PARENTS
!!      forces,m_mep,mover
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fred2fcart(favg,fcart,fred,gprimd,jellslab,natom)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fred2fcart'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,jellslab
!arrays
 real(dp),intent(out) :: fcart(3,natom)
 real(dp),intent(in) :: fred(3,natom)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(out) :: favg(3)

!Local variables-------------------------------
!scalars
 integer :: iatom,mu
!arrays

! *************************************************************************

!Note conversion to cartesian coordinates (bohr) AND
!negation to make a force out of a gradient
 favg(:)=zero
 do iatom=1,natom
   do mu=1,3
     fcart(mu,iatom)= - (gprimd(mu,1)*fred(1,iatom)+&
&     gprimd(mu,2)*fred(2,iatom)+&
&     gprimd(mu,3)*fred(3,iatom))
     favg(mu)=favg(mu)+fcart(mu,iatom)
   end do
 end do

!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
 favg(:)=favg(:)/dble(natom)
 if(jellslab/=0) favg(3)=zero
 do iatom=1,natom
   fcart(:,iatom)=fcart(:,iatom)-favg(:)
 end do

end subroutine fred2fcart
!!***

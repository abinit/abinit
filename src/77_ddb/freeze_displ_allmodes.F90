!!****f* ABINIT/freeze_displ_allmodes
!!
!! NAME
!! freeze_displ_allmodes
!!
!! FUNCTION
!!  From a given set of phonon modes, generate and output supercells and
!!  displaced configurations of atoms.
!!  Typically useful to follow soft modes and see distorsions of crystal structures
!!
!! COPYRIGHT
!! Copyright (C) 2010-2017 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public Licence, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! amu(ntypat) = mass of the atoms (atomic mass unit)
!! displ(2,3*natom,3*natom) = phonon mode displacements (complex)
!! freeze_displ = amplitude of the displacement to freeze into the supercell
!! natom = number of atoms in the unit cell
!! ntypat = number of atom types
!! phfrq(3*natom) = phonon frequencies
!! qphnrm = norm of phonon q vector (should be 1 or 0)
!! qphon = phonon wavevector
!! rprimd(3,3) = dimensionfull primitive translations in real space
!! typat(natom) = integer label of each type of atom (1,2,...)
!! xcart(3,natom) = cartesian coords of atoms in unit cell (bohr)
!!
!! OUTPUT
!! for the moment only prints to file, but could also return pointer to supercell object, with
!! rprimd and atomic positions, for further use
!!
!! NOTES
!! freeze_displ could be determined automatically from a temperature and the phonon frequency,
!! as the average displacement of the mode with a Bose distribution.
!!
!! PARENTS
!!      m_phonons
!!
!! CHILDREN
!!      destroy_supercell,freeze_displ_supercell,init_supercell_for_qpt
!!      prt_supercell_for_qpt
!!
!! SOURCE
!!

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine freeze_displ_allmodes(displ, freeze_displ, natom, outfile_radix, phfreq,  &
&         qphon, rprimd, typat, xcart, znucl)


 use defs_basis
 use m_profiling_abi
 use m_supercell

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'freeze_displ_allmodes'
!End of the abilint section

 implicit none

! arguments
! scalar
 integer,intent(in) :: natom
 character(len=*),intent(in) :: outfile_radix
 real(dp), intent(in) :: freeze_displ

!arrays
 integer,intent(in) :: typat(natom)

 real(dp),intent(in) :: displ(2,3*natom,3*natom)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(in) :: phfreq(3*natom)
 real(dp),intent(in) :: qphon(3)
 real(dp),intent(in) :: xcart(3,natom)
 real(dp),intent(in) :: znucl(:) 

! local vars
 integer :: jmode
 type(supercell_type) :: scell

! *************************************************************************

!determine supercell needed to freeze phonon
 call init_supercell_for_qpt(natom, qphon, rprimd, typat, xcart, znucl, scell)

 do jmode = 1, 3*natom
! reset positions
   scell%xcart = scell%xcart_ref

!  displace atoms according to phonon jmode
   call freeze_displ_supercell(displ(:,:,jmode), freeze_displ, scell)

!  print out everything for this wavevector and mode
   call prt_supercell_for_qpt (phfreq(jmode), jmode, outfile_radix, scell)
 end do

 call destroy_supercell (scell)

end subroutine freeze_displ_allmodes
!!***

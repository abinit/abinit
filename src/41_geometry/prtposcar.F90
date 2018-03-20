!{\src2tex{textfont=tt}}
!!****f* ABINIT/prtposcar
!! NAME
!!  prtposcar
!!
!! FUNCTION
!!  output VASP style POSCAR and FORCES files for use with frozen phonon codes, like
!!  PHON from Dario Alfe' or frophon
!!  IMPORTANT: the order of atoms is fixed such that typat is re-grouped. First typat=1
!!  then typat=2, etc...
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  fcart = forces on atoms in cartesian coordinates
!!  natom = number of atoms
!!  ntypat = number of types of atoms
!!  rprimd = lattice vectors for the primitive cell
!!  typat = type for each of the natom atoms
!!  ucvol = unit cell volume 
!!  xred = reduced positions of the atoms
!!  znucl = nuclear charge of each atomic type
!!
!! OUTPUTS
!!   Only files written
!!
!! PARENTS
!!      afterscfloop
!!
!! CHILDREN
!!      atomdata_from_znucl
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prtposcar(fcart, fnameradix, natom, ntypat, rprimd, typat, ucvol, xred, znucl)

 use defs_basis
 use m_profiling_abi
 use m_atomdata
 use m_errors

 use m_io_tools,     only : open_file

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prtposcar'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: natom, ntypat
 real(dp), intent(in) :: ucvol
!arrays
 integer, intent(in) :: typat(natom)
 real(dp), intent(in) :: fcart(3,natom)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,natom)
 real(dp), intent(in) :: znucl(ntypat)
 character(len=fnlen), intent(in) :: fnameradix

!Local variables-------------------------------
!scalars
 integer :: iatom, itypat, iout
 type(atomdata_t) :: atom
! arrays
 integer :: natoms_this_type(ntypat)
 character(len=2) :: symbol
 character(len=7) :: natoms_this_type_str
 character(len=100) :: chem_formula, natoms_all_types
 character(len=500) :: msg
 character(len=fnlen) :: fname

!************************************************************************

!output POSCAR file for positions, atom types etc
 fname = trim(fnameradix)//"_POSCAR"
 if (open_file(fname,msg,newunit=iout) /= 0) then
   MSG_ERROR(msg)
 end if

 natoms_this_type = 0
 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) == itypat) natoms_this_type(itypat) = natoms_this_type(itypat)+1
   end do
 end do

 chem_formula = ""
 do itypat=1, ntypat
   call atomdata_from_znucl(atom,znucl(itypat))
   symbol = atom%symbol
   if (natoms_this_type(itypat) < 10) then
     write (natoms_this_type_str, '(I1)') natoms_this_type(itypat)
   else if (natoms_this_type(itypat) < 100) then
     write (natoms_this_type_str, '(I2)') natoms_this_type(itypat)
   else if (natoms_this_type(itypat) < 1000) then
     write (natoms_this_type_str, '(I3)') natoms_this_type(itypat)
   end if
   chem_formula = trim(chem_formula) // symbol // trim(natoms_this_type_str)
 end do
 write (iout,'(2a)') "ABINIT generated POSCAR file. Title string - should be chemical formula... ",&
& trim(chem_formula)

 write (iout,'(E24.14)') -ucvol*Bohr_Ang*Bohr_Ang*Bohr_Ang
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,1) ! (angstr? bohr?)
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,2) 
 write (iout,'(3E24.14,1x)') Bohr_Ang*rprimd(:,3) 
 natoms_all_types = "   "
 do itypat=1, ntypat
   write (natoms_this_type_str, '(I7)') natoms_this_type(itypat)
   natoms_all_types = trim(natoms_all_types) // "   " // trim(natoms_this_type_str)
 end do
 write (iout,'(a)') trim(natoms_all_types)
 write (iout,'(a)') "Direct"
 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) /= itypat) cycle
     write (iout,'(3(E24.14,1x))') xred(:,iatom)
   end do
 end do
 close (iout)


!output FORCES file for forces in same order as positions above
 fname = trim(fnameradix)//"_FORCES"
 if (open_file(fname,msg,newunit=iout) /= 0 ) then
   MSG_ERROR(msg)
 end if 
!ndisplacements
!iatom_displaced displacement_red_coord(3)
!forces_cart_ev_Angstr(3)
!...
!<repeat for other displaced atoms>
 write (iout,'(I7)') 1
 write (iout,'(a)') '1 0 0 0        ! TO BE FILLED IN '
 do itypat=1, ntypat
   do iatom=1,natom
     if (typat(iatom) /= itypat) cycle
     write (iout,'(3(E24.14,1x))') Ha_eV/Bohr_Ang*fcart(:,iatom)
   end do
 end do
 close (iout)

end subroutine prtposcar
!!***

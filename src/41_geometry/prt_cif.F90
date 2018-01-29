!{\src2tex{textfont=tt}}
!!****f* ABINIT/prt_cif
!! NAME
!! prt_cif
!!
!! FUNCTION
!!   print out CIF format file, from abinit internal data
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      outscfcv
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine prt_cif(brvltt, ciffname, natom, nsym, ntypat, rprimd, &
&   spgaxor, spgroup, spgorig, symrel, tnon, typat, xred, znucl)

 use defs_basis
 use m_profiling_abi
 use m_errors
 use m_atomdata

 use m_io_tools,       only : open_file
 use m_fstrings,       only : int2char10

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'prt_cif'
 use interfaces_41_geometry, except_this_one => prt_cif
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom, ntypat, nsym
 integer, intent(in) :: brvltt, spgaxor, spgroup, spgorig
!arrays
 integer, intent(in) :: typat(natom)
 integer, intent(in) :: symrel(3,3,nsym)
 character(len=*), intent(in) :: ciffname
 real(dp), intent(in) :: tnon(3,nsym)
 real(dp), intent(in) :: rprimd(3,3)
 real(dp), intent(in) :: xred(3,natom)
 real(dp), intent(in) :: znucl(ntypat)

!Local variables -------------------------------
!scalars
 integer :: unitcif
 integer :: iatom, isym
 integer :: sporder
 integer :: itypat, nat_this_type
 real(dp) :: ucvol
 type(atomdata_t) :: atom

!arrays
 character(len=80) :: tmpstring

 character(len=1) :: brvsb
 character(len=15) :: intsb,ptintsb,ptschsb,schsb
 character(len=35) :: intsbl

 character(len=10) :: str_nat_type
 character(len=100) :: chemformula
 character(len=500) :: msg

 real(dp) :: angle(3)
 real(dp) :: gprimd(3,3)
 real(dp) :: rmet(3,3), gmet(3,3)

!*************************************************************************

!open file in append mode xlf and other compilers refuse append mode
 if (open_file(ciffname,msg,newunit=unitcif) /=0) then
   MSG_WARNING(msg)
   return
 end if

!print title for dataset
 write (unitcif,'(a)') 'data_set'

!print cell parameters a,b,c, angles, volume
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 angle(1)=acos(rmet(2,3)/sqrt(rmet(2,2)*rmet(3,3)))/two_pi*360.0_dp
 angle(2)=acos(rmet(1,3)/sqrt(rmet(1,1)*rmet(3,3)))/two_pi*360.0_dp
 angle(3)=acos(rmet(1,2)/sqrt(rmet(1,1)*rmet(2,2)))/two_pi*360.0_dp

 write (unitcif,'(a,E20.10)') '_cell_length_a                     ', sqrt(rmet(1,1))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_length_b                     ', sqrt(rmet(2,2))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_length_c                     ', sqrt(rmet(3,3))*Bohr_Ang
 write (unitcif,'(a,E20.10)') '_cell_angle_alpha                  ', angle(1)
 write (unitcif,'(a,E20.10)') '_cell_angle_beta                   ', angle(2)
 write (unitcif,'(a,E20.10)') '_cell_angle_gamma                  ', angle(3)
 write (unitcif,'(a,E20.10)') '_cell_volume                       ', ucvol*(Bohr_Ang)**3

!print reduced positions
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a,E20.10)') '  _atom_site_label                   '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_x                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_y                 '
 write (unitcif,'(a,E20.10)') '  _atom_site_fract_z                 ' 
 do iatom = 1, natom
   call atomdata_from_znucl(atom,znucl(typat(iatom)))
   write (unitcif,'(2a,3E20.10)') '  ', atom%symbol, xred(:,iatom)
 end do

!
!other specs in CIF dictionary which may be useful:
!GEOM_BOND GEOM_ANGLE GEOM_TORSION
!

!print chemical composition in simplest form
 chemformula = "'"
 do itypat = 1, ntypat
   nat_this_type = 0
   do iatom = 1, natom
     if (typat(iatom) == itypat) nat_this_type = nat_this_type+1
   end do
   call atomdata_from_znucl(atom,znucl(itypat))
   call int2char10(nat_this_type, str_nat_type)
   chemformula = trim(chemformula) // atom%symbol // trim(str_nat_type) // "  " 
 end do
 chemformula = trim(chemformula) // "'"
 write (unitcif,'(2a)') '_chemical_formula_analytical              ', chemformula

!FIXME: check that brvltt is correctly used here - is it equal to bravais(1) in the invars routines?
 if     (brvltt==1)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             triclinic'
 else if(brvltt==2)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             monoclinic'
 else if(brvltt==3)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             orthorhombic'
 else if(brvltt==4)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             tetragonal'
 else if(brvltt==5)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             rhombohedral'
 else if(brvltt==6)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             hexagonal'
 else if(brvltt==7)then 
   write (unitcif,'(a)') '_symmetry_cell_setting             cubic'
 end if

 call spgdata(brvsb,intsb,intsbl,ptintsb,ptschsb,&
& schsb,spgaxor,spgroup,sporder,spgorig)

!print symmetry operations
 write (unitcif,'(a,I6)') "_symmetry_Int_Tables_number          ", spgroup
 write (unitcif,'(5a)') "_symmetry_space_group_name_H-M        '", brvsb, " ", trim(intsb), "'"
 write (unitcif,'(a)') ''
 write (unitcif,'(a)') 'loop_'
 write (unitcif,'(a)') '  _symmetry_equiv_pos_as_xyz           '
 do isym = 1, nsym
   call  symrel2string(symrel(:,:,isym), tnon(:,isym), tmpstring)
   write (unitcif,'(2a)') '  ', trim(tmpstring)
 end do

 close (unitcif)

end subroutine prt_cif
!!***

!!****f* ABINIT/symrel2string
!! NAME
!! symrel2string
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 2010-2018 ABINIT group (MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! NOTES
!!
!! PARENTS
!!      prt_cif
!!
!! CHILDREN
!!
!! SOURCE

subroutine symrel2string(symrel1, tnon, string)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symrel2string'
!End of the abilint section

 implicit none
 
!Arguments ------------------------------------
!scalars
!arrays
 integer, intent(in) :: symrel1(3,3)
 real(dp), intent(in) :: tnon(3)
 character(len=80), intent(out) :: string

!Local variables -------------------------
!scalars
 integer :: i1,i2
 character(len=1) :: xyz(3)

! *********************************************************************

 xyz(1) = 'x'
 xyz(2) = 'y'
 xyz(3) = 'z'

 string = ''
 do i1=1,3
   if (abs(tnon(i1)) > tol10) then
!    find fraction 1/n for tnon, otherwise do not know what to print
     if (abs(one-two*tnon(i1)) < tol10) string = trim(string)//'1/2'
     if (abs(one+two*tnon(i1)) < tol10) string = trim(string)//'-1/2'

     if (abs(one-three*tnon(i1)) < tol10) string = trim(string)//'1/3'
     if (abs(one+three*tnon(i1)) < tol10) string = trim(string)//'-1/3'
     if (abs(two-three*tnon(i1)) < tol10) string = trim(string)//'2/3'
     if (abs(two+three*tnon(i1)) < tol10) string = trim(string)//'-2/3'

     if (abs(one-six*tnon(i1)) < tol10) string = trim(string)//'1/6'
     if (abs(one+six*tnon(i1)) < tol10) string = trim(string)//'-1/6'
     if (abs(five-six*tnon(i1)) < tol10) string = trim(string)//'5/6'
     if (abs(five+six*tnon(i1)) < tol10) string = trim(string)//'-5/6'
   end if
   do i2=1,3
!    FIXME: check if this is correct ordering for symrel(i1,i2) looks ok
     if (symrel1(i1,i2) == 1)  string = trim(string)//'+'//xyz(i2)
     if (symrel1(i1,i2) == -1) string = trim(string)//'-'//xyz(i2)
   end do
   if (i1 /= 3) string = trim(string)//','
 end do

end subroutine symrel2string
!!***

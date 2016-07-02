!{\src2tex{textfont=tt}}
!!****f* ABINIT/fixsym
!!
!! NAME
!! fixsym
!!
!! FUNCTION
!! Using input indsym which tells which atoms are related by symmetry,
!! check that iatfix consistently fixes (freezes) all atoms which are
!! related by symmetry, i.e. that iatfix does not break symmetry.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! iatfix(3,natom)=integer array with 1 in every position for which
!!  the atom is to be kept fixed
!!  NOTE that this is not the input data structure for iatfix but it is
!!  the internal data structure used through most of the subroutines
!! indsym(4,nsym,natom)=indirect indexing array for symmetrically related
!!  atoms; 4th element is label of symmetrically related atom
!! natom=number of atoms
!! nsym=number of symmetries (should be > 1 when this is called)
!!
!! OUTPUT
!!  (only checking)
!!
!! NOTE
!!  Stops execution with an error message if iatfix breaks symmetry.
!!
!! PARENTS
!!      gstate
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine fixsym(iatfix,indsym,natom,nsym)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'fixsym'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: iatfix(3,natom),indsym(4,nsym,natom)

!Local variables-------------------------------
!scalars
 integer :: iatom,isym,jatom
 character(len=500) :: message

! *************************************************************************

 if (nsym > 1) then
   do iatom=1,natom
     do isym=1,nsym
!      jatom is the label of a symmetrically related atom
       jatom=indsym(4,isym,iatom)
!      Thus the atoms jatom and iatom must be fixed along the same directions
       if ( iatfix(1,jatom) /=  iatfix(1,iatom) .or. &
&       iatfix(2,jatom) /=  iatfix(2,iatom) .or. &
&       iatfix(3,jatom) /=  iatfix(3,iatom)       ) then
         write(message, '(a,i6,a,a,i6,a,a,a,a,a,a,a)' )&
&         '  Atom number ',jatom,' is symmetrically',&
&         ' equivalent to atom number ',iatom,',',ch10,&
&         '  but according to iatfix, iatfixx, iatfixy and iatfixz, they',ch10,&
&         '  are not fixed along the same directions, which is forbidden.',ch10,&
&         '  Action : modify either the symmetry or iatfix(x,y,z) and resubmit.'
         MSG_ERROR(message)
       end if
     end do
   end do
 end if

end subroutine fixsym
!!***

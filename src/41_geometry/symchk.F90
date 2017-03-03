!{\src2tex{textfont=tt}}
!!****f* ABINIT/symchk
!! NAME
!! symchk
!!
!! FUNCTION
!! Symmetry checker for atomic coordinates.
!! Checks for translated atomic coordinate tratom(3) to agree
!! with some coordinate xred(3,iatom) where atomic types agree too.
!! All coordinates are "reduced", i.e. given in terms of primitive
!! reciprocal translations.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! natom=number of atoms in unit cell
!! tratom(3)=reduced coordinates for a single atom which presumably
!!   result from the application of a symmetry operation to an atomic
!!   coordinate
!! trtypat=type of atom (integer) translated to tratom
!! typat(natom)=types of all atoms in unit cell (integer)
!! xred(3,natom)=reduced coordinates for all atoms in unit cell
!!
!! OUTPUT
!! difmin(3)=minimum difference between apparently equivalent atoms
!!   (give value separately for each coordinate)--note that value
!!   may be NEGATIVE so take abs later if needed
!! eatom=atom label of atom which is SAME as tratom to within a primitive
!!   cell translation ("equivalent atom")
!! transl(3)=primitive cell translation to make iatom same as tratom (integers)
!!
!! PARENTS
!!      symatm
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symchk(difmin,eatom,natom,tratom,transl,trtypat,typat,xred)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symchk'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,trtypat
 integer,intent(out) :: eatom
!arrays
 integer,intent(in) :: typat(natom)
 integer,intent(out) :: transl(3)
 real(dp),intent(in) :: tratom(3),xred(3,natom)
 real(dp),intent(out) :: difmin(3)

!Local variables-------------------------------
!scalars
 integer :: iatom,jatom,trans1,trans2,trans3
 real(dp) :: test,test1,test2,test3,testmn

! *************************************************************************

!Start testmn out at large value
 testmn=1000000.d0

!Loop through atoms--
!when types agree, check for agreement after primitive translation
 jatom=1
 do iatom=1,natom
   if (trtypat/=typat(iatom)) cycle

!  Check all three components
   test1=tratom(1)-xred(1,iatom)
   test2=tratom(2)-xred(2,iatom)
   test3=tratom(3)-xred(3,iatom)
!  Find nearest integer part of difference
   trans1=nint(test1)
   trans2=nint(test2)
   trans3=nint(test3)
!  Check whether, after translation, they agree
   test1=test1-dble(trans1)
   test2=test2-dble(trans2)
   test3=test3-dble(trans3)

   test=abs(test1)+abs(test2)+abs(test3)
   if (test<tol10) then
!    Note that abs() is not taken here
     difmin(1)=test1
     difmin(2)=test2
     difmin(3)=test3
     jatom=iatom
     transl(1)=trans1
     transl(2)=trans2
     transl(3)=trans3
!    Break out of loop when agreement is within tolerance
     exit
   else
!    Keep track of smallest difference if greater than tol10
     if (test<testmn) then
       testmn=test
!      Note that abs() is not taken here
       difmin(1)=test1
       difmin(2)=test2
       difmin(3)=test3
       jatom=iatom
       transl(1)=trans1
       transl(2)=trans2
       transl(3)=trans3
     end if
   end if

!  End loop over iatom. Note a "cycle" and an "exit" inside the loop
 end do

 eatom=jatom

end subroutine symchk
!!***

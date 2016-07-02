!{\src2tex{textfont=tt}}
!!****f* ABINIT/symdet
!! NAME
!! symdet
!!
!! FUNCTION
!! Compute determinant of each input symmetry matrix sym(3,3,i)
!! and check that the determinant is always +/- 1.  Integer arithmetic.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2016 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym=number of symmetry operations
!! sym(3,3,nsym)=integer symmetry array
!!
!! OUTPUT
!! determinant(nsym)=determinant of each symmetry operation
!!
!! PARENTS
!!      remove_inversion,setsym,symptgroup,symspgr
!!
!! CHILDREN
!!      mati3det
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine symdet(determinant,nsym,sym)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'symdet'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 integer,intent(out) :: determinant(nsym)

!Local variables-------------------------------
!scalars
 integer :: det,isym
 character(len=500) :: message

! *************************************************************************

 do isym=1,nsym
   call mati3det(sym(:,:,isym),det)
   determinant(isym)=det
   if (abs(det)/=1) then
     write(message,'(a,i5,a,i10,a,a,a,a,a)')&
&     'Abs(determinant) for symmetry number',isym,' is',det,' .',ch10,&
&     'For a legitimate symmetry, abs(determinant) must be 1.',ch10,&
&     'Action: check your symmetry operations (symrel) in input file.'
     MSG_ERROR(message)
   end if
 end do

end subroutine symdet
!!***

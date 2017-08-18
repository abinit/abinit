!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkgrp
!! NAME chkgrp
!! chkgrp
!!
!! FUNCTION
!! Checks that a set of input symmetries constitutes a group.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! nsym = number of symmetry operations
!! symafm = (anti)ferromagnetic part of symmetry operations
!! symrel = 3D matrix containg symmetry operations
!!
!! OUTPUT
!!  ierr=Status error.
!!
!! TODO
!! SHOULD ALSO CHECK THE tnons !
!!
!! PARENTS
!!      chkinp,gensymspgr,m_bz_mesh,m_esymm,m_sigmaph,setsym
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine chkgrp(nsym,symafm,symrel,ierr)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkgrp'
 use interfaces_32_util
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)

!Local variables-------------------------------
!scalars
 integer :: isym,jsym,ksym,symafmchk,testeq=1
 logical :: found_inv
 character(len=500) :: msg
!arrays
 integer :: chk(3,3)

! *************************************************************************

!DEBUG
!write(std_out,*)' chkgrp : enter'
!write(std_out,*)'     isym         symrel            symafm '
!do isym=1,nsym
!write(std_out,'(i3,a,9i3,a,i3)' )isym,'   ',symrel(:,:,isym),'   ',symafm(isym)
!end do
!ENDDEBUG

 ierr = 0

!1) Identity must be the first symmetry.
 if (ANY(symrel(:,:,1) /= identity_3d .or. symafm(1)/=1 )) then
   MSG_WARNING("First operation must be the identity operator")
   ierr = ierr+1
 end if
!
!2) The inverse of each element must belong to the group.
 do isym=1,nsym
   call mati3inv(symrel(:,:,isym),chk)
   chk = TRANSPOSE(chk)
   found_inv = .FALSE.
   do jsym=1,nsym
     if ( ALL(symrel(:,:,jsym) == chk) .and. (symafm(jsym)*symafm(isym) == 1 )) then
       found_inv = .TRUE.; EXIT
     end if
   end do

   if (.not.found_inv) then
     write(msg,'(a,i0,2a)')&
&     "Cannot find the inverse of symmetry operation ",isym,ch10,&
&     "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr+1
   end if

 end do
!
!Closure relation under composition.
 do isym=1,nsym
   do jsym=1,nsym
!
!    Compute the product of the two symmetries
     chk = MATMUL(symrel(:,:,jsym), symrel(:,:,isym))
     symafmchk=symafm(jsym)*symafm(isym)
!
!    Check that product array is one of the original symmetries.
     do ksym=1,nsym
       testeq=1
       if ( ANY(chk/=symrel(:,:,ksym) )) testeq=0
#if 0
!      FIXME this check make v4/t26 and v4/t27 fails.
!      The rotational part is in the group but with different magnetic part!
       if (symafmchk/=symafm(ksym))testeq=0
#endif
       if (testeq==1) exit ! The test is positive
     end do
!
     if(testeq==0) then ! The test is negative
       write(msg, '(a,2i3,a,7a)' )&
&       'product of symmetries',isym,jsym,' is not in group.',ch10,&
&       'This indicates that the input symmetry elements',ch10,&
&       'do not possess closure under group composition.',ch10,&
&       'Action: check symrel, symafm and fix them.'
       MSG_WARNING(msg)
       ierr = ierr+1
     end if

   end do ! jsym
 end do ! isym

end subroutine chkgrp
!!***


!!****f* ABINIT/sg_multable
!! NAME
!! sg_multable
!!
!! FUNCTION
!! Checks that a set of input symmetries constitutes a group.
!!
!! INPUTS
!! nsym=number of symmetry operations
!! symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!! symrel(3,3,nsym)=symmetry operations in real space.
!! tnons(3,nsym)=Fractional translations.
!!
!! OUTPUT
!!  ierr=Status error. A non-zero value signals a failure.
!!  [multable(4,nsym,nsym)]= Optional output.
!!    multable(1,sym1,sym2) gives the index of the symmetry product S1 * S2 in the symrel array. 0 if not found.
!!    multable(2:4,sym1,sym2)= the lattice vector that has to added to the fractional translation
!!      of the operation of index multable(1,sym1,sym2) to obtain the fractional traslation of the product S1 * S2.
!!  [toinv(4,nsym)]= Optional output.
!!    toinv(1,sym1)=Gives the index of the inverse of the symmetry operation.
!!     S1 * S1^{-1} = {E, L} with E identity and L a lattice vector
!!    toinv(2:4,sym1)=The lattice vector L
!!      Note that toinv can be easily obtained from multable but sometimes we do not need the full table.
!!
!! TODO
!!  This improved version should replace chkgrp.
!!
!! PARENTS
!!      m_crystal,m_shirley
!!
!! CHILDREN
!!
!! SOURCE

subroutine sg_multable(nsym,symafm,symrel,tnons,tnons_tol,ierr,multable,toinv)

 use defs_basis
 use m_errors
 use m_profiling_abi

 use m_numeric_tools,  only : isinteger

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sg_multable'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(out) :: ierr
 real(dp),intent(in) :: tnons_tol
!arrays
 integer,intent(in) :: symafm(nsym),symrel(3,3,nsym)
 integer,optional,intent(out) :: multable(4,nsym,nsym)
 integer,optional,intent(out) :: toinv(4,nsym)
 real(dp),intent(in) :: tnons(3,nsym)

!Local variables-------------------------------
!scalars
 integer :: sym1,sym2,sym3,prd_symafm
 logical :: found_inv,iseq
 character(len=500) :: msg
!arrays
 integer :: prd_symrel(3,3)
 real(dp) :: prd_tnons(3)

! *************************************************************************

 ierr = 0

!1) Identity must be the first symmetry. Do not check tnons, cell might not be primitive.
 if (ANY(symrel(:,:,1) /= identity_3d .or. symafm(1)/=1 )) then
   MSG_WARNING("First operation must be the identity operator")
   ierr = ierr+1
 end if
!
!2) The inverse of each element must belong to the group.
 do sym1=1,nsym
   found_inv = .FALSE.
   do sym2=1,nsym
     prd_symrel = MATMUL(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_tnons = tnons(:,sym1) + MATMUL(symrel(:,:,sym1),tnons(:,sym2))
     prd_symafm = symafm(sym1)*symafm(sym2)
     if ( ALL(prd_symrel == identity_3d) .and. isinteger(prd_tnons,tnons_tol) .and. prd_symafm == 1 ) then
       found_inv = .TRUE.
       if (PRESENT(toinv)) then
         toinv(1,sym1) = sym2
         toinv(2:4,sym1) = NINT(prd_tnons)
       end if
       EXIT
     end if
   end do

   if (.not.found_inv) then
     write(msg,'(a,i0,2a)')&
&     "Cannot find the inverse of symmetry operation ",sym1,ch10,&
&     "Input symmetries do not form a group "
     MSG_WARNING(msg)
     ierr = ierr+1
   end if
 end do
!
!Check closure relation under composition and construct multiplication table.
 do sym1=1,nsym
   do sym2=1,nsym
!
!    Compute the product of the two symmetries. Convention {A,a} {B,b} = {AB, a+Ab}
     prd_symrel = MATMUL(symrel(:,:,sym1), symrel(:,:,sym2))
     prd_symafm = symafm(sym1)*symafm(sym2)
     prd_tnons = tnons(:,sym1) + MATMUL(symrel(:,:,sym1),tnons(:,sym2))
!
     iseq=.FALSE.
     do sym3=1,nsym ! Check that product array is one of the original symmetries.
       iseq = ( ALL(prd_symrel==symrel(:,:,sym3) )           .and. &
&       isinteger(prd_tnons-tnons(:,sym3),tnons_tol) .and. &
&       prd_symafm==symafm(sym3) )  ! Here v4/t26 and v4/t27 will fail.
!      The rotational part is in the group but with different magnetic part!

       if (iseq) then ! The test is positive
         if (PRESENT(multable)) then
           multable(1,sym1,sym2) = sym3
           multable(2:4,sym1,sym2) = NINT(prd_tnons-tnons(:,sym3))
         end if
         EXIT
       end if
     end do
!
     if (.not.iseq) then ! The test is negative
       write(msg, '(a,2(i0,1x),a,7a)' )&
&       'product of symmetries:',sym1,sym2,' is not in group.',ch10,&
&       'This indicates that the input symmetry elements',ch10,&
&       'do not possess closure under group composition.',ch10,&
&       'Action: check symrel, symafm and fix them.'
       MSG_WARNING(msg)
       ierr = ierr+1
       if (PRESENT(multable)) then
         multable(1,sym1,sym2) = 0
         multable(2:4,sym1,sym2) = HUGE(0)
       end if
     end if

   end do ! sym2
 end do ! sym1

end subroutine sg_multable
!!***

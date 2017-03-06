!{\src2tex{textfont=tt}}
!!****f* ABINIT/chkorthsy
!! NAME chkorthsy
!! chkorthsy
!!
!! FUNCTION
!! Check the orthogonality of the symmetry operations
!! (lengths and absolute values of scalar products should be preserved)
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive transl. for reciprocal space (bohr**-1)
!! rmet=Real space metric.
!! nsym=actual number of symmetries
!! rprimd(3,3)=dimensional primitive translations for real space (bohr)
!! symrel(3,3,1:nsym)=symmetry operations in real space in terms of primitive translations
!!
!! SIDE EFFECTS
!! iexit= if 0 at input, will do the check, and stop if there is a problem, return 0 if no problem
!!        if 1 at input, will always input, return 0 if no problem, -1 if there is a problem,
!!                       also, suppresses printing of problem
!!
!! PARENTS
!!      chkinp,ingeo,symmetrize_rprimd
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine chkorthsy(gprimd,iexit,nsym,rmet,rprimd,symrel)

 use defs_basis
 use m_errors
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chkorthsy'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
 integer,intent(inout) :: iexit
!arrays
 integer,intent(in) :: symrel(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3),rmet(3,3),rprimd(3,3)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,jj
 real(dp),parameter :: tol=2.0d-8
 real(dp) :: residual,rmet2
 character(len=500) :: message
!arrays
 real(dp) :: prods(3,3),rmet_sym(3,3),rprimd_sym(3,3)

! *************************************************************************

!  DEBUG
!  write(std_out,*)' chkorthsy : enter '
!  write(std_out,*)rprimd(1:3,1)
!  write(std_out,*)rprimd(1:3,2)
!  write(std_out,*)rprimd(1:3,3)
!  ENDDEBUG

 rmet2=zero
 do ii=1,3
   do jj=1,3
     rmet2=rmet2+rmet(ii,jj)*2
   end do
 end do

!Loop over all symmetry operations
 do isym=1,nsym

!  Compute symmetric of primitive vectors under point symmetry operations
   do ii=1,3
     rprimd_sym(:,ii)=symrel(1,ii,isym)*rprimd(:,1)+&
&     symrel(2,ii,isym)*rprimd(:,2)+&
&     symrel(3,ii,isym)*rprimd(:,3)
   end do

!  DEBUG
   !write(std_out,*)' chkorthsy : isym=',isym
!  write(std_out,*)rprimd(1:3,1)
!  write(std_out,*)rprimd(1:3,2)
!  write(std_out,*)rprimd(1:3,3)
!  write(std_out,*)rprimd_sym(1:3,1)
!  write(std_out,*)rprimd_sym(1:3,2)
!  write(std_out,*)rprimd_sym(1:3,3)
!  ENDDEBUG

!  If the new lattice is the same as the original one,
!  the lengths and angles are preserved
   do ii=1,3
     rmet_sym(ii,:)=rprimd_sym(1,ii)*rprimd_sym(1,:)+&
&     rprimd_sym(2,ii)*rprimd_sym(2,:)+&
&     rprimd_sym(3,ii)*rprimd_sym(3,:)
   end do

!  DEBUG
!  write(std_out,*)' rmet :'
!  write(std_out,*)rmet(1:3,1)
!  write(std_out,*)rmet(1:3,2)
!  write(std_out,*)rmet(1:3,3)
!  write(std_out,*)rmet_sym(1:3,1)
!  write(std_out,*)rmet_sym(1:3,2)
!  write(std_out,*)rmet_sym(1:3,3)
!  ENDDEBUG

   residual=zero
   do ii=1,3
     do jj=1,3
       residual=residual+(rmet_sym(ii,jj)-rmet(ii,jj))**2
     end do
   end do

   if(sqrt(residual)>tol*sqrt(rmet2))then
     if(iexit==0)then
       write(message, '(a,i5,a,a,a,a,a,es12.4,a,a,a,a,a,a,a)' )&
&       'The symmetry operation number',isym,' does not preserve',ch10,&
&       'vector lengths and angles.',ch10,&
&       'The value of the residual is',residual,'.',ch10,&
&       'Action : modify rprim, acell and/or symrel so that',ch10,&
&       'vector lengths and angles are preserved.',ch10,&
&       'Beware, the tolerance on symmetry operations is very small.'
       MSG_ERROR(message)
     else
       iexit=-1
     end if
   end if

!  Also, the scalar product of rprimd_sym and gprimd must give integer numbers
   do ii=1,3
     prods(ii,:)=rprimd_sym(1,ii)*gprimd(1,:)+ &
&     rprimd_sym(2,ii)*gprimd(2,:)+ &
&     rprimd_sym(3,ii)*gprimd(3,:)
   end do

   do ii=1,3
     do jj=1,3
       residual=prods(ii,jj)-anint(prods(ii,jj))
       if(abs(residual)>tol)then
         if(iexit==0)then
           write(message, '(a,i0,a,a,a,a,a,a,a)' )&
&           'The symmetry operation number',isym,' generates',ch10,&
&           'a different lattice.',ch10,&
&           'Action : modify rprim, acell and/or symrel so that',ch10,&
&           'the lattice is preserved.'
           MSG_ERROR(message)
         else
           iexit=-1
         end if
       end if
     end do
   end do

   if(iexit==-1) exit

 end do ! isym

 if(iexit==1)iexit=0

end subroutine chkorthsy
!!***

!{\src2tex{textfont=tt}}
!!****f* ABINIT/sygrad
!!
!! NAME
!! sygrad
!!
!! FUNCTION
!! Symmetrize derivatives of energy with respect to coordinates.
!! Unsymmetrized gradients are input as dedt; symmetrized grads are then
!! placed in fred.
!! If nsym=1 simply copy dedt into fred (only symmetry is identity).
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
!!  dedt(3,natom)=unsymmetrized gradients wrt dimensionless tn (hartree)
!!  nsym=number of symmetry operators in group
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!    reciprocal space primitive translations--see comments below
!!  indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!   label which gets rotated into given atom by given symmetry
!!   (first three elements are related primitive translation--
!!   see symatm where this is computed)
!!
!! OUTPUT
!! fred(3,3,natom)=symmetrized gradients wrt reduced coordinates (hartree)
!!
!! NOTES
!! symmetrization of gradients with respect to reduced
!! coordinates tn is conducted according to the expression
!! $[d(e)/d(t(n,a))]_{symmetrized} = (1/Nsym)*Sum(S)*symrec(n,m,S)*
!!              [d(e)/d(t(m,b))]_{unsymmetrized}$
!! where $t(m,b)= (symrel^{-1})(m,n)*(t(n,a)-tnons(n))$ and tnons
!! is a possible nonsymmorphic translation.  The label "b" here
!! refers to the atom which gets rotated into "a" under symmetry "S".
!! symrel is the symmetry matrix in real space, which is the inverse
!! transpose of symrec.  symrec is the symmetry matrix in reciprocal
!! space.  $sym_{cartesian} = R * symrel * R^{-1} = G * symrec * G^{-1}$
!! where the columns of R and G are the dimensional primitive translations
!! in real and reciprocal space respectively.
!! Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine sygrad(fred,natom,dedt,nsym,symrec,indsym)

 use m_profiling_abi

 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sygrad'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym)
 real(dp),intent(in) :: dedt(3,natom)
 real(dp),intent(out) :: fred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ind,isym,mu
 real(dp),parameter :: tol=1.0d-30
 real(dp) :: summ

! *************************************************************************
!
 if (nsym==1) then
!  only symmetry is identity so simply copy
   fred(:,:)=dedt(:,:)
 else
!  actually conduct symmetrization
   do ia=1,natom
     do mu=1,3
       summ=0._dp
       do isym=1,nsym
         ind=indsym(4,isym,ia)
         summ=summ+dble(symrec(mu,1,isym))*dedt(1,ind)+&
&         dble(symrec(mu,2,isym))*dedt(2,ind)+&
&         dble(symrec(mu,3,isym))*dedt(3,ind)
       end do
       fred(mu,ia)=summ/dble(nsym)
       if(abs(fred(mu,ia))<tol)fred(mu,ia)=0.0_dp
     end do
   end do
 end if

end subroutine sygrad
!!***

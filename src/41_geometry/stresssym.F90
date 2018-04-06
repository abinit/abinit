!{\src2tex{textfont=tt}}
!!****f* ABINIT/stresssym
!! NAME
!! stresssym
!!
!! FUNCTION
!! For given order of point group, symmetrizes the stress tensor,
!! in symmetrized storage mode and cartesian coordinates, using input
!! 3x3 symmetry operators in reduced coordinates.
!! symmetrized tensor replaces input tensor.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2018 ABINIT group (DCA, XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!! nsym=order of group.
!! sym(3,3,nsym)=symmetry operators (usually symrec=expressed in terms
!!               of action on reciprocal lattice primitive translations);
!!               integers.
!!
!! OUTPUT
!! stress(6)=stress tensor, in cartesian coordinates, in symmetric storage mode
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      dfpt_nselt,dfpt_nstpaw,forstrnps,littlegroup_pert,pawgrnl,stress
!!
!! CHILDREN
!!      matr3inv,strconv
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


subroutine stresssym(gprimd,nsym,stress,sym)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'stresssym'
 use interfaces_32_util
 use interfaces_41_geometry, except_this_one => stresssym
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nsym
!arrays
 integer,intent(in) :: sym(3,3,nsym)
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(inout) :: stress(6)

!Local variables-------------------------------
!scalars
 integer :: ii,isym,mu,nu
 real(dp) :: summ,tmp
!arrays
 real(dp) :: rprimd(3,3),rprimdt(3,3),strfrac(6),tensor(3,3),tt(3,3)

!*************************************************************************

!Obtain matrix of real space dimensional primitive translations
!(inverse tranpose of gprimd), and its transpose
 call matr3inv(gprimd,rprimd)
 rprimdt=transpose(rprimd)

!Compute stress tensor in reduced coordinates
 call strconv(stress,rprimdt,strfrac)

!Switch to full storage mode
 tensor(1,1)=strfrac(1)
 tensor(2,2)=strfrac(2)
 tensor(3,3)=strfrac(3)
 tensor(3,2)=strfrac(4)
 tensor(3,1)=strfrac(5)
 tensor(2,1)=strfrac(6)
 tensor(2,3)=tensor(3,2)
 tensor(1,3)=tensor(3,1)
 tensor(1,2)=tensor(2,1)

 do nu=1,3
   do mu=1,3
     tt(mu,nu)=tensor(mu,nu)/dble(nsym)
     tensor(mu,nu)=0.0_dp
   end do
 end do

!loop over all symmetry operations:
 do isym=1,nsym
   do mu=1,3
     do nu=1,3
       summ=0._dp
       do ii=1,3
         tmp=tt(ii,1)*sym(nu,1,isym)+tt(ii,2)*sym(nu,2,isym)+&
&         tt(ii,3)*sym(nu,3,isym)
         summ=summ+sym(mu,ii,isym)*tmp
       end do
       tensor(mu,nu)=tensor(mu,nu)+summ
     end do
   end do
 end do

!Switch back to symmetric storage mode
 strfrac(1)=tensor(1,1)
 strfrac(2)=tensor(2,2)
 strfrac(3)=tensor(3,3)
 strfrac(4)=tensor(3,2)
 strfrac(5)=tensor(3,1)
 strfrac(6)=tensor(2,1)

!Convert back stress tensor (symmetrized) in cartesian coordinates
 call strconv(strfrac,gprimd,stress)

end subroutine stresssym
!!***

!!****f* ABINIT/strconv
!! NAME
!! strconv
!!
!! FUNCTION
!! If original gprimd is input, convert from symmetric storage mode
!! 3x3 tensor in reduced coordinates "frac" to symmetric storage mode
!! symmetric tensor in cartesian coordinates "cart".
!!
!! INPUTS
!!  frac(6)=3x3 tensor in symmetric storage mode, reduced coordinates
!!  gprimd(3,3)=reciprocal space dimensional primitive translations (bohr^-1)
!!
!! OUTPUT
!!  cart(6)=symmetric storage mode for symmetric 3x3 tensor in cartesian coords.
!!
!! NOTES
!! $cart(i,j)=G(i,a) G(j,b) frac(a,b)$
!! "Symmetric" storage mode for 3x3 tensor is 6 element array with
!! elements 11, 22, 33, 32, 31, and 21.
!! "cart" may be same array as "frac".
!! If rprimd transpose is input instead of gprimd, then convert tensor
!! in cartesian coordinates to reduced coordinates
!!
!! PARENTS
!!      ctocprj,d2frnl,mkcore,mkcore_paw,mkcore_wvl,nonlop_pl,nonlop_ylm
!!      stresssym
!!
!! CHILDREN
!!
!! SOURCE

subroutine strconv(frac,gprimd,cart)

 use defs_basis
 use m_profiling_abi

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'strconv'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: frac(6),gprimd(3,3)
 real(dp),intent(inout) :: cart(6) ! alias of frac   !vz_i

!Local variables-------------------------------
!scalars
 integer :: ii,jj
!arrays
 real(dp) :: work1(3,3),work2(3,3)

! *************************************************************************

 work1(1,1)=frac(1)
 work1(2,2)=frac(2)
 work1(3,3)=frac(3)
 work1(3,2)=frac(4) ; work1(2,3)=frac(4)
 work1(3,1)=frac(5) ; work1(1,3)=frac(5)
 work1(2,1)=frac(6) ; work1(1,2)=frac(6)

 do ii=1,3
   work2(:,ii)=zero
   do jj=1,3
     work2(:,ii)=work2(:,ii)+gprimd(ii,jj)*work1(:,jj)
   end do
 end do

 do ii=1,3
   work1(ii,:)=zero
   do jj=1,3
     work1(ii,:)=work1(ii,:)+gprimd(ii,jj)*work2(jj,:)
   end do
 end do

 cart(1)=work1(1,1)
 cart(2)=work1(2,2)
 cart(3)=work1(3,3)
 cart(4)=work1(2,3)
 cart(5)=work1(1,3)
 cart(6)=work1(1,2)

end subroutine strconv
!!***

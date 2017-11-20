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
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR)
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

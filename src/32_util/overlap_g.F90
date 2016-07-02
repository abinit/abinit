!{\src2tex{textfont=tt}}
!!****f* ABINIT/overlap_g
!! NAME
!! overlap_g
!!
!! FUNCTION
!! Compute the scalar product between WF at two different k-points
!! < u_{n,k1} | u_{n,k2} >
!!
!! COPYRIGHT
!! Copyright (C) 1999-2016 ABINIT group ()
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! mpw = maximum dimensioned size of npw
!! npw_k1 = number of plane waves at k1
!! npw_k2 = number of plane waves at k2
!! nspinor = 1 for scalar, 2 for spinor wavefunctions
!! pwind_k = array required to compute the scalar product (see initberry.f)
!! vect1 = wavefunction at k1: | u_{n,k1} >
!! vect2 = wavefunction at k1: | u_{n,k2} >
!!
!! OUTPUT
!! doti = imaginary part of the scalarproduct
!! dotr = real part of the scalarproduct
!!
!! NOTES
!! In case a G-vector of the basis sphere of plane waves at k1
!! does not belong to the basis sphere of plane waves at k2,
!! pwind = 0. Therefore, the dimensions of vect1 &
!! vect2 are (1:2,0:mpw) and the element (1:2,0) MUST be
!! set to zero.
!!
!! The current implementation if not compatible with TR-symmetry (i.e. istwfk/=1) !
!!
!! PARENTS
!!      dfptff_bec,dfptff_die,dfptff_ebp,dfptff_edie,dfptff_gbefd
!!      dfptff_gradberry,qmatrix,smatrix
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine overlap_g(doti,dotr,mpw,npw_k1,npw_k2,nspinor,pwind_k,vect1,vect2)

 use defs_basis
 use m_profiling_abi
 use m_errors

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'overlap_g'
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mpw,npw_k1,npw_k2,nspinor
 real(dp),intent(out) :: doti,dotr
!arrays
 integer,intent(in) :: pwind_k(mpw)
 real(dp),intent(in) :: vect1(1:2,0:mpw*nspinor),vect2(1:2,0:mpw*nspinor)

!Local variables-------------------------------
!scalars
 integer :: ipw,ispinor,jpw,spnshft1,spnshft2
 character(len=500) :: message

! *************************************************************************

!Check if vect1(:,0) = 0 and vect2(:,0) = 0
 if ((abs(vect1(1,0)) > tol12).or.(abs(vect1(2,0)) > tol12).or. &
& (abs(vect2(1,0)) > tol12).or.(abs(vect2(2,0)) > tol12)) then
   message = ' vect1(:,0) and/or vect2(:,0) are not equal to zero'
   MSG_BUG(message)
 end if

!Compute the scalar product
 dotr = zero; doti = zero
 do ispinor = 1, nspinor
   spnshft1 = (ispinor-1)*npw_k1
   spnshft2 = (ispinor-1)*npw_k2
!$OMP PARALLEL DO PRIVATE(jpw) REDUCTION(+:doti,dotr) 
   do ipw = 1, npw_k1
     jpw = pwind_k(ipw)
     dotr = dotr + vect1(1,spnshft1+ipw)*vect2(1,spnshft2+jpw) + vect1(2,spnshft1+ipw)*vect2(2,spnshft2+jpw)
     doti = doti + vect1(1,spnshft1+ipw)*vect2(2,spnshft2+jpw) - vect1(2,spnshft1+ipw)*vect2(1,spnshft2+jpw)
   end do
 end do

end subroutine overlap_g
!!***

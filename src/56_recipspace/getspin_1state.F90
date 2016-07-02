!{\src2tex{textfont=tt}}
!!****f* ABINIT/getspin_1state
!! NAME
!! getspin_1state
!!
!! FUNCTION
!!  sandwich a single wave function on the Pauli matrices
!!
!! COPYRIGHT
!! Copyright (C) 2014-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  npw_k = number of plane waves
!!  cgcband = coefficients of spinorial wave function
!!
!! OUTPUT
!!  spin = 3-vector of spin components for this state
!!  cgcmat_ = outer spin product of spinorial wf with itself
!!
!! PARENTS
!!      m_cut3d,partial_dos_fractions
!!
!! CHILDREN
!!      cg_outer_spin_product
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine  getspin_1state(cgcband, npw_k, spin, cgcmat_)

 use m_errors
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'getspin_1state'
 use interfaces_56_recipspace, except_this_one => getspin_1state
!End of the abilint section

 implicit none

!input var
 integer, intent(in) :: npw_k
 real(dp), intent(in) :: cgcband(2,2*npw_k)
 complex(dpc), intent(out),optional :: cgcmat_(2,2)
 real(dp), intent(out) :: spin(3)


!local var
 complex(dpc) :: pauli_0(2,2) 
 complex(dpc) :: pauli_x(2,2) 
 complex(dpc) :: pauli_y(2,2) 
 complex(dpc) :: pauli_z(2,2)
 complex(dpc) :: cspin(0:3)
 complex(dpc) :: cgcmat(2,2)
! code
! ***********************************************************************

 pauli_0 = reshape((/cone,czero,czero,cone/), (/2,2/))
 pauli_x = reshape((/czero,cone,cone,czero/), (/2,2/))
 pauli_y = reshape((/czero,j_dpc,-j_dpc,czero/), (/2,2/))
 pauli_z = reshape((/cone,czero,czero,-cone/), (/2,2/))
! print *,  "pauli_0 :"
! print *, pauli_0(1,1), pauli_0(1,2)
! print *, pauli_0(2,1), pauli_0(2,2)
! print *,  "pauli_x :"
! print *, pauli_x(1,1), pauli_x(1,2)
! print *, pauli_x(2,1), pauli_x(2,2)
! print *,  "pauli_y :"
! print *, pauli_y(1,1), pauli_y(1,2)
! print *, pauli_y(2,1), pauli_y(2,2)
! print *,  "pauli_z :"
! print *, pauli_z(1,1), pauli_z(1,2)
! print *, pauli_z(2,1), pauli_z(2,2)
 
 call cg_outer_spin_product(cgcband, npw_k, cgcmat)

 cspin = czero 
! spin(*)  = sum_{si sj pi} cgcband(si,pi)^* pauli_*(si,sj) cgcband(sj,pi)
 cspin(0) = cgcmat(1,1)*pauli_0(1,1) + cgcmat(2,1)*pauli_0(2,1) &
& + cgcmat(1,2)*pauli_0(1,2) + cgcmat(2,2)*pauli_0(2,2)
 cspin(1) = cgcmat(1,1)*pauli_x(1,1) + cgcmat(2,1)*pauli_x(2,1) &
& + cgcmat(1,2)*pauli_x(1,2) + cgcmat(2,2)*pauli_x(2,2)
 cspin(2) = cgcmat(1,1)*pauli_y(1,1) + cgcmat(2,1)*pauli_y(2,1) &
& + cgcmat(1,2)*pauli_y(1,2) + cgcmat(2,2)*pauli_y(2,2)
 cspin(3) = cgcmat(1,1)*pauli_z(1,1) + cgcmat(2,1)*pauli_z(2,1) &
& + cgcmat(1,2)*pauli_z(1,2) + cgcmat(2,2)*pauli_z(2,2)

!print *, 'real(spin) ', real(cspin)
!print *, 'aimag(spin) ', aimag(cspin)

 spin = real(cspin(1:3))

 if (present(cgcmat_)) then
   cgcmat_ = cgcmat
 end if

end subroutine getspin_1state
!!***


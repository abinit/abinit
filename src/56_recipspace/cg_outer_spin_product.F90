!{\src2tex{textfont=tt}}
!!****f* ABINIT/cg_outer_spin_product
!! NAME
!! cg_outer_spin_product
!!
!! FUNCTION
!!  trace a single wave function, keeping spin matrix of degrees of freedom
!!
!! COPYRIGHT
!! Copyright (C) 2014-2016 ABINIT group (MVer)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cgcband(2,nspinor*npw_k) = wave function
!!  npw_k = number of plane waves
!!
!! OUTPUT
!!  cgcmat = 2x2 matrix of spin components (dpcomplex)
!!
!! PARENTS
!!      getspin_1state
!!
!! CHILDREN
!!      zgemm
!!
!! NOTES
!!  this routine is forced for nspinor = 2!!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine cg_outer_spin_product(cgcband, npw_k, cgcmat)

 use m_errors
 use defs_basis

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cg_outer_spin_product'
!End of the abilint section

 implicit none

!input var
 integer, intent(in) :: npw_k
 real(dp), intent(in) :: cgcband(2,2*npw_k)
 complex(dpc), intent(out) :: cgcmat(2,2)

! ***********************************************************************

! cgcmat = cgcband * cgcband^T*
 cgcmat = czero
 call zgemm('n','c',2,2,npw_k,cone,cgcband,2,cgcband,2,czero,cgcmat,2)

! debug
 write(std_out,'(a)') " cgcmat :"
 write(std_out,'(a, 2(2E20.10, 2x))') "  ", cgcmat(1,1), cgcmat(1,2) 
 write(std_out,'(a, 2(2E20.10, 2x))') "  ", cgcmat(2,1), cgcmat(2,2)
!enddebug

end subroutine cg_outer_spin_product
!!***

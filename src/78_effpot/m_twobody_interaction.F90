!!****m* ABINIT/m_twobody_interaction
!! NAME
!! m_twobody_interaction
!!
!! FUNCTION
!! This module contains an twobody any interger order interaction
!!
!! Datatypes:
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_twobody_interaction
  use defs_basis
  use m_errors
  use m_abicore
  use m_spmat_ndcoo, only: ndcoo_mat_t
  use m_dynamic_array, only: int2d_array_type, real_array_type, int_array_type
  implicit none

  private
!!***

  !public :: build_index
  public :: get_twobody_dEdx
  public :: get_twobody_delta_E

contains

  subroutine get_twobody_dEdx(coeff, x, f, E)
    type(ndcoo_mat_t), intent(in) :: coeff
    real(dp), intent(in) :: x(:)
    real(dp), intent(inout) :: f(:), E
    integer :: inz, i, j, oi, oj
    real(dp) :: val
    do inz =1, coeff%nnz
       i=coeff%ind%data(1,inz)
       j=coeff%ind%data(2,inz)
       oi=coeff%ind%data(3,inz)
       oj=coeff%ind%data(4,inz)
       val=coeff%val%data(inz)
       f(i)=f(i)+ oi* val * x(i)**(oi-1)* x(j)**oj
       f(j)=f(j)+ oj* val * x(i)**oi* x(j)**(oj-1)
       E=E+val*x(i)**oi*x(j)**oj
    end do
  end subroutine get_twobody_dEdx

  subroutine get_twobody_delta_E(coeff, x, ix, dE)
    type(ndcoo_mat_t), intent(in) :: coeff
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: ix
    real(dp), intent(inout) :: dE
    integer :: inz, i, j, oi, oj
    ABI_UNUSED_A(inz)
    ABI_UNUSED_A(i)
    ABI_UNUSED_A(j)
    ABI_UNUSED_A(oi)
    ABI_UNUSED_A(oj)
    ABI_UNUSED_A(coeff)
    ABI_UNUSED_A(x)
    ABI_UNUSED_A(ix)
    ABI_UNUSED_A(dE)
  end subroutine get_twobody_delta_E

end module m_twobody_interaction




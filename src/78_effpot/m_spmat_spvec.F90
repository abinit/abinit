!!****m* ABINIT/m_spmat_spvec
!! NAME
!! m_spmat_spvec
!!
!! FUNCTION
!! This module contains the a sparse vector.
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

module m_spmat_spvec
  use defs_basis  
  use m_xmpi
  use m_errors
  use m_abicore
  use m_dynamic_array, only: real_array_type, int_array_type
  implicit none
  private
  !!***

  type, public :: sp_real_vec
     integer :: size=0
     integer :: nnz=0
     type(int_array_type) :: ids
     type(real_array_type) :: vals
   contains
     procedure :: initialize
     procedure :: push
     procedure :: finalize
     procedure :: dot => vec_spvec_dot
     procedure :: plus_Ax => spvec_plus_Ax
  end type sp_real_vec

  public :: vec_spvec_dot
  public :: spvec_plus_Ax

contains

  subroutine initialize(self, size)
    class(sp_real_vec), intent(inout) :: self
    integer, intent(in) :: size
    self%size=size
  end subroutine initialize


  subroutine push(self, id, val)
    class(sp_real_vec), intent(inout) :: self
    integer, intent(in) :: id
    real(dp), intent(in) :: val
    call self%ids%push(id)
    call self%vals%push(val)
    self%nnz=self%nnz+1
  end subroutine push

  subroutine finalize(self)
    class(sp_real_vec), intent(inout) :: self
    call self%ids%finalize()
    call self%vals%finalize()
    self%size=0
    self%nnz=0
  end subroutine finalize

  function vec_spvec_dot(spvec, vec) result (ret)
    class(sp_real_vec), intent(in) :: spvec
    real(dp), intent(in) :: vec(spvec%size)
    real(dp) :: ret
    integer :: i
    ret=0.0_dp
    do i=1, spvec%nnz
       ret=ret+vec(spvec%ids%data(i)) * spvec%vals%data(i)
    end do
  end function vec_spvec_dot

  subroutine spvec_plus_Ax(spvec, a, y)
    class(sp_real_vec), intent(in) :: spvec
    real(dp), intent(inout) :: y(spvec%size)
    real(dp), intent(in) :: a
    integer :: i, j
    do i =1, spvec%nnz
       j=spvec%ids%data(i)
       y(j)=y(j) + a * spvec%vals%data(j)
    end do
  end subroutine spvec_plus_Ax



end module m_spmat_spvec

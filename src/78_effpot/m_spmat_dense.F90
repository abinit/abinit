!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spmat_dense
!! NAME
!! m_spmat_dense
!!
!! FUNCTION
!! This module contains the dense matrix as sparse matrix type.
!!
!! Datatypes:
!!  dense_mat_t: dense matrix pretending to be sparse, mostly for test. Sometimes when the matrix is not so sparse, it is then more efficient to use this.
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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


module m_spmat_dense
  use m_errors
  use defs_basis
  use m_spmat_base, only: base_mat2d_t
  implicit none
!!***
  private
  type, extends(base_mat2d_t), public :: dense_mat_t
     real(dp), allocatable :: mat(:,:)
   contains
     procedure :: initialize => dense_mat_t_initialize
     procedure :: finalize => dense_mat_t_finalize
     procedure :: mv => dense_mat_t_mv
  end type dense_mat_t

contains
  subroutine dense_mat_t_initialize(self,  mshape)
    class(dense_mat_t), intent(inout) :: self
    integer, intent(in)::mshape(:)
    if(size(mshape)/=2) MSG_ERROR("mshape should be size 2 for dense_mat_t")
    self%nrow=mshape(1)
    self%ncol=mshape(2)
    ABI_ALLOCATE(self%mat, (self%nrow, self%ncol))
    self%mat(:,:)=0.0d0
  end subroutine dense_mat_t_initialize

  subroutine dense_mat_t_finalize(self)
    class(dense_mat_t), intent(inout) :: self
    if (allocated(self%mat))  ABI_DEALLOCATE(self%mat)
    if (allocated(self%mshape)) ABI_DEALLOCATE(self%mat)
    self%ncol=0
    self%nrow=0
    self%ndim=0
  end subroutine dense_mat_t_finalize

  subroutine dense_mat_insert(self, irow, icol, val)
    class(dense_mat_t), intent(inout) :: self
    integer, intent(inout) :: irow, icol
    real(dp), intent(in) :: val
    self%mat(irow, icol)=val
  end subroutine dense_mat_insert

  ! dense matrix-vector multiplication, using blas DGEMV
  subroutine dense_mat_t_mv(self, x, b)
    class(dense_mat_t), intent(in) :: self
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
    call dgemv("N", self%nrow, self%ncol, 1.0d0,self%mat , 2,  x, 1, 0.0d0,  b, 1)
  end subroutine dense_mat_t_mv

end module m_spmat_dense

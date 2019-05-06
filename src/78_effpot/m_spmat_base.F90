!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spmat_base
!! NAME
!! m_spmat_base
!!
!! FUNCTION
!! This module contains the base type for sparse matrix. 
!!
!! Datatypes:
!!  base_mat_t: base sparse matrix.
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
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


module m_spmat_base
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  implicit none
  private
!!***
  type, public :: base_mat_t
     integer :: ndim
     integer, allocatable:: mshape(:)
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: add_entry
  end type base_mat_t

  type, public, extends(base_mat_t) ::  base_mat2d_t
     integer :: nrow, ncol
   contains
     procedure :: initialize => base_mat2d_t_initialize
     procedure :: mv => base_mat2d_t_mv
  end type base_mat2d_t

contains

  subroutine initialize(self, mshape)
    class(base_mat_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    ABI_ALLOCATE(self%mshape, (size(mshape)))
    self%mshape=mshape
    self%ndim=size(mshape)
  end subroutine initialize

  subroutine finalize(self)
    class(base_mat_t), intent(inout) :: self
    self%ndim=0
    if (allocated(self%mshape)) then
       ABI_DEALLOCATE(self%mshape)
    endif
  end subroutine finalize

  subroutine add_entry(self, ind, val)
    class(base_mat_t), intent(inout) :: self
    integer, intent(in) :: ind(self%ndim)
    real(dp), intent(in) :: val
  end subroutine add_entry

  subroutine base_mat2d_t_initialize(self,mshape)
    class(base_mat2d_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    call self%base_mat_t%initialize(mshape)
    self%nrow=mshape(1)
    self%ncol=mshape(2)
  end subroutine base_mat2d_t_initialize

  subroutine base_mat2d_t_mv(self, x, b)
    class(base_mat2d_t), intent(in) :: self
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
  end subroutine base_mat2d_t_mv

end module m_spmat_base

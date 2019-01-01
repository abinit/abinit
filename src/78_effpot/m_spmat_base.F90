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
!! Copyright (C) 2001-2018 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

module m_spmat_base
  use defs_basis
  implicit none
  private
  type, public ::  base_mat_t
     integer :: nrow, ncol
   contains
     procedure :: mv => base_mat_t_mv
  end type base_mat_t

contains
  subroutine base_mat_t_mv(self, x, b)
    class(base_mat_t), intent(in) :: self
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
  end subroutine base_mat_t_mv

end module m_spmat_base

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

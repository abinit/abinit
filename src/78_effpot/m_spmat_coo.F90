!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_spmat_coo
!! NAME
!! m_spmat_coo
!!
!! FUNCTION
!! This module contains the a COO (coordinate) format of sparse matrix.
!! The efficiency of mat vec multiplication is fine but not as good as CSR
!! Datatypes:
!!  COO_mat_t: COO matrix
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

module m_spmat_COO
  use defs_basis  
  use m_xmpi
  use m_errors
  use m_abicore
  use m_spmat_base
  use m_spmat_ndcoo, only: ndcoo_mat_t
  implicit none
  private
!!***

  !!----------- COO ------------------------
  ! COO sparse matrix.
  ! i, j, val are the row index, col index and value of each entry.
  ! nnz: number of non-zeros.
  type, public, extends(ndcoo_mat_t) :: COO_mat_t
   contains
     procedure :: initialize
     procedure :: mv
  end type COO_mat_t


contains

  subroutine initialize(self, mshape)
    class(coo_mat_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    if (size(mshape)/=2) stop 1
    call self%ndcoo_mat_t%initialize(mshape)
  end subroutine initialize

  ! COO sparse matrix-vector multiplication. naive implementation.
  subroutine mv(self, x, b)
    class(COO_mat_t), intent(in) :: self
    real(dp), intent(in) :: x(self%mshape(1))
    real(dp), intent(out) :: b(self%mshape(2))
    integer:: ind, ind_i, ind_j
    b(:)=0.0D0
    do ind = 1, self%nnz, 1
       ind_i=self%ind%data(1, ind)
       ind_j=self%ind%data(2, ind)
       b(ind_i)=b(ind_i)+self%val%data(ind)*x(ind_j)
    end do
  end subroutine  mv

  subroutine COO_mat_t_mv_mpi(self, x ,b)
    class(coo_mat_t), intent(in) :: self
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: b(:)
    !real(dp):: my_b(self%nrow)
    integer :: ierr, irow, icol
    call xmpi_bcast(x, 0, xmpi_world, ierr)
    b(:)=0.0_dp

    ! TODO implement.
    MSG_ERROR("mpi COO mv Not implemented yet")
    ! TODO : use gather instead of reduce.
    !call mpi_reduce(my_b, b, self%nrow, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call xmpi_sum_master(b, 0, xmpi_world, ierr )
  end subroutine COO_mat_t_mv_mpi

end module m_spmat_COO

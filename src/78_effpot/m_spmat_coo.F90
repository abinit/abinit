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
!! Copyright (C) 2001-2018 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#include "abi_common.h"
module m_spmat_COO
  use defs_basis  
  use m_xmpi
  use m_spmat_base
  implicit none
!!***

  !!----------- COO ------------------------
  ! COO sparse matrix.
  ! i, j, val are the row index, col index and value of each entry.
  ! nnz: number of non-zeros.
  type, public, extends(base_mat_t) :: COO_mat_t
     integer :: nnz
     integer, allocatable :: i(:), j(:)
     real(dp), allocatable:: val(:)
   contains
     procedure :: initialize => coo_mat_t_initialize
     procedure :: finalize => coo_mat_t_finalize
     procedure :: mv => coo_mat_t_mv
  end type COO_mat_t

contains

  ! COO matrix
  subroutine COO_mat_t_initialize(self, nrow, ncol, nnz, i, j, val)
    class(COO_mat_t), intent(inout) :: self
    integer, intent(in) :: nnz, nrow, ncol
    integer , intent(in), optional :: i(:), j(:)
    real(dp), intent(in), optional:: val(:)
    integer :: err
    allocate(self%i(nnz), stat=err)
    allocate(self%j(nnz), stat=err)
    allocate(self%val(nnz), stat=err)
    self%nrow=nrow
    self%ncol=ncol
    self%nnz=nnz
    if(present(i)) then
       self%i(:)=i(:)
    end if
 
    if(present(j)) then
       self%j(:)=j(:)
    end if
    if(present(val)) then
       self%val(:)=val(:)
    endif
  end subroutine COO_mat_t_initialize

  subroutine COO_mat_t_finalize(self)
    class(COO_mat_t), intent(inout) :: self
    self%nrow=0
    self%ncol=0
    self%nnz=0
    if(allocated(self%i)) then
       ABI_DEALLOCATE(self%i)
    end if
    if(allocated(self%j)) then
       ABI_DEALLOCATE(self%j)
    end if
    if(allocated(self%val)) then
       ABI_DEALLOCATE(self%val)
    end if

  end subroutine COO_mat_t_finalize

  ! COO sparse matrix-vector multiplication. naive implementation.
  subroutine COO_mat_t_mv(self, x, b)
    class(COO_mat_t), intent(in) :: self
    real(dp), intent(in):: x(self%ncol)
    real(dp), intent(out):: b(self%nrow)
    integer:: ind, ind_i, ind_j
    b(:)=0.0D0
    do ind = 1, self%nnz, 1
       ind_i=self%i(ind)
       ind_j=self%j(ind)
       b(ind_i)=b(ind_i)+self%val(ind)*x(ind_j)
    end do
  end subroutine  COO_mat_t_mv

  subroutine COO_mat_t_mv_mpi(self, x ,b)
    class(coo_mat_t), intent(in) :: self
    real(dp), intent(inout) :: x(:)
    real(dp), intent(out) :: b(:)
    !real(dp):: my_b(self%nrow)
    integer :: ierr, irow, icol
    !call mpi_bcast(x, self%ncol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call xmpi_bcast(x, 0, xmpi_world, ierr)
    !my_b(:)=0.0_dp
    b(:)=0.0_dp

    ! TODO implement
    print *, "mpi COO mv Not implemented yet"
    ! TODO : use gather instead of reduce.
    !call mpi_reduce(my_b, b, self%nrow, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

    call xmpi_sum_master(b, 0, xmpi_world, ierr )
  end subroutine COO_mat_t_mv_mpi


end module m_spmat_COO

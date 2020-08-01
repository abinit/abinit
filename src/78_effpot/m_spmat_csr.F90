!!****m* ABINIT/m_spmat_csr
!! NAME
!! m_spmat_csr
!!
!! FUNCTION
!! This module contains the a CSR (compressed row) format of sparse matrix.
!! Efficient for mat vec multiplication.
!! MPI matvec is also implemented.
!! Datatypes:
!!  CSR_mat_t: CSR matrix
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

module m_spmat_csr
  use defs_basis  
  use m_abicore
  use m_xmpi
  use m_errors
  use m_spmat_base, only: base_mat2d_t
  use m_mpi_scheduler
  implicit none
!!***
  private
  !!----------- CSR ------------------------
  ! CSR sparse matrix
  ! nnz: number of non-zeros.
  ! icol : column index of entries .size:nnz
  ! row_shift: row_shift(irow) to row_shift(irow+1)-1  are the index
  !            of column and values for entries. size : nrow+1
  ! val: values of non-zero entries. size(nnz)
  !type, public, extends(base_mat_t) :: CSR_mat_t
  type, public, extends(base_mat2d_t) :: CSR_mat_t
     integer :: nnz
     integer, allocatable :: icol(:), row_shift(:)
     real(dp), allocatable:: val(:)
     type(mpi_scheduler_t) :: mps
   contains
     procedure :: initialize
     procedure :: set
     procedure :: finalize=> csr_mat_t_finalize
     procedure :: mv => csr_mat_t_mv ! mat vec multiplication serial version
     procedure :: sync ! sync data to all mpi ranks
     procedure :: mv_mpi => csr_mat_t_mv_mpi ! mpi version of mv
     procedure :: mv_select_row =>csr_mat_t_mv_select_row ! mv of selected rows
     procedure :: mv_one_row =>csr_mat_t_mv_one_row ! mv of one rows
  end type CSR_mat_t

contains

  !-----------------------------------------------------------------------
  !> @brief initialization
  !> @param [in] mshape: shape of matrix, should be size 2.
  !-----------------------------------------------------------------------
  subroutine initialize(self, mshape)
    class(csr_mat_t), intent(inout) :: self
    integer, intent(in) :: mshape(:)
    if (size(mshape)/=2) stop "mshape should be size 2"
    call self%base_mat_t%initialize(mshape)
    self%nrow=mshape(1)
    self%ncol=mshape(2)
  end subroutine initialize

  !-----------------------------------------------------------------------
  !> @brief set the full csr matrix
  !> @param [in] nnz: number of entries
  !> @param [in] icol:
  !> @param [in] row_shift: 
  !> @param [in] val: values:
  !-----------------------------------------------------------------------
  subroutine set(self,  nnz, icol, row_shift, val)
    class(CSR_mat_t), intent(inout) :: self
    integer, intent(in) :: nnz 
    ! i: col number of each entry
    ! j: first 0, n1, n1+n2, ...
    ! val(irow, irow+1) are the values of entries in row irow.
    integer , intent(in), optional :: icol(:), row_shift(:)
    real(dp), intent(in), optional :: val(:)
    integer :: iproc

    iproc=xmpi_comm_rank(xmpi_world)
    if (iproc/=0) then
       MSG_ERROR("This function (CSR_MAT%set) should be only used on root node")
    end if

    self%nnz=nnz
    if(.not. allocated(self%icol)) then
       ABI_ALLOCATE(self%icol, (self%nnz))
    endif

    if(.not. allocated(self%row_shift)) then
       ABI_ALLOCATE(self%row_shift, (self%nrow+1))
    endif

    if(.not. allocated(self%val)) then
       ABI_ALLOCATE(self%val, (self%nnz))
    endif

    if (present(icol)) then
       self%icol(:)=icol(:)
    end if
    if (present(row_shift)) then
       self%row_shift(:)=row_shift(:)
    end if
    if (present(icol)) then
       self%val(:)=val(:)
    end if

  end subroutine set

  !-----------------------------------------------------------------------
  !> @brief sync the matrix from master node to all nodes
  !> @param [in] master: the id of master node
  !> @param [in] comm: the communicator
  !> @param [in] nblock: the minimal block for assigning the tasks to nodes.
  !-----------------------------------------------------------------------
  subroutine sync(self, master, comm, nblock)
    class (csr_mat_t), intent(inout) :: self
    integer , intent(in) :: master, comm, nblock
    integer :: ierr, iproc
    iproc=xmpi_comm_rank(xmpi_world)

    call xmpi_barrier(xmpi_world)
    call xmpi_bcast(self%ndim, master, comm, ierr)
    call xmpi_bcast(self%ncol, master, comm, ierr)
    call xmpi_bcast(self%nrow, master, comm, ierr)
    call xmpi_bcast(self%nnz, master, comm, ierr)

    call self%mps%initialize(self%nrow/nblock, master, comm, nblock)
    if (.not. self%mps%irank==master) then
       if(.not. allocated(self%mshape)) then
          ABI_ALLOCATE(self%mshape, (self%ndim))
       endif

       if(.not. allocated(self%icol)) then
          ABI_ALLOCATE(self%icol, (self%nnz))
       endif

       if(.not. allocated(self%row_shift)) then
          ABI_ALLOCATE(self%row_shift, (self%nrow+1))
       endif
       if(.not. allocated(self%val)) then
          ABI_ALLOCATE(self%val, (self%nnz))
       endif
    end if

    ! TODO: no need to send all the data to all the
    ! only the corresponding row data.
    call xmpi_bcast(self%mshape, master, comm, ierr)
    call xmpi_bcast(self%icol, master, comm, ierr)
    call xmpi_bcast(self%row_shift, master, comm, ierr)
    call xmpi_bcast(self%val, master, comm, ierr)
  end subroutine sync


  !-----------------------------------------------------------------------
  !> @brief Finalize
  !-----------------------------------------------------------------------
  subroutine CSR_mat_t_finalize(self)
    class(CSR_mat_t), intent(inout) :: self
    self%ncol=0
    self%nrow=0
    self%nnz=0
    self%ndim=0
    call self%mps%finalize()
    if(allocated(self%icol)) then
       ABI_DEALLOCATE(self%icol)
    endif
    if(allocated(self%row_shift)) then
       ABI_DEALLOCATE(self%row_shift)
    endif
    if(allocated(self%val)) then
       ABI_DEALLOCATE(self%val)
    endif
    if(allocated(self%mshape)) then
       ABI_DEALLOCATE(self%mshape)
    endif
  end subroutine CSR_mat_t_finalize


  !-----------------------------------------------------------------------
  !> @brief Matrix vector multiplication
  !> @param [in] x : M x = b
  !> @param [out] b: M x = b
  !-----------------------------------------------------------------------
  subroutine CSR_mat_t_mv(self, x, b)
    class(CSR_mat_t), intent(in):: self
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
    integer::irow, i1, i2, i
    b(:)=0.0d0
    !$OMP PARALLEL DO private(i, i1, i2)
    do irow=1, self%nrow
        i1=self%row_shift(irow)
        i2=self%row_shift(irow+1)-1
        do i=i1, i2
            b(irow)=b(irow)+ self%val(i)*x(self%icol(i))
        end do
    enddo
    !$OMP END PARALLEL DO
  end subroutine CSR_mat_t_mv

  !-----------------------------------------------------------------------
  !> @brief Matrix vector multiplication (mpi version)
  !> @param [in] x : M x = b
  !> @param [out] b: M x = b
  !-----------------------------------------------------------------------
  subroutine CSR_mat_t_mv_mpi(self, x, b, bcastx, syncb)
    class(CSR_mat_t), intent(in) :: self
    real(dp), intent(inout) :: x(self%ncol)
    logical, intent(in) :: bcastx, syncb
    real(dp), intent(out) :: b(self%nrow)
    !real(dp):: my_b(self%nrow)
    integer :: ierr, irow,  i1, i2, i
    if (bcastx) then
        call xmpi_bcast(x, 0, xmpi_world, ierr)
    end if
    !if (.not. iam_master) b(:)=0.0_dp
    !my_b(:)=0.0_dp
    b(:)=0.0_dp
    do irow= self%mps%istart, self%mps%iend
       i1=self%row_shift(irow)
       i2=self%row_shift(irow+1)-1
       do i=i1, i2
          b(irow)=b(irow)+ self%val(i)*x(self%icol(i))
       end do
    enddo
    ! TODO : use gather instead of reduce?
    !call mpi_reduce(my_b, b, self%nrow, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    if (syncb) then
       call xmpi_sum_master(b, 0, xmpi_world, ierr )
    endif
  end subroutine CSR_mat_t_mv_mpi

  !-----------------------------------------------------------------------
  !> @brief multiple a row, indexed by  of M by x:   y_i=  Mij x_j
  !> @param [in] nrow: number of rows.
  !> @param [in] ind_row: indices i
  !> @param [in]  x: x
  !> @param [out]  y: y
  !-----------------------------------------------------------------------
  subroutine CSR_mat_t_mv_one_row(self,  j, x, y)
    class(CSR_mat_t), intent(in)::self 
    integer, intent(in) ::j
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: y
    integer :: i, irow, i1, i2 
    y=0.0_dp
    i1=self%row_shift(j)
    i2=self%row_shift(j+1)-1
    do i=i1, i2
       y= y+self%val(i)*x(self%icol(i))
    end do
  end subroutine CSR_mat_t_mv_one_row


  !-----------------------------------------------------------------------
  !> @brief multiple a submatrix (rows indexed by i) of M by x:   y_i=\sum M_ij x_j
  !> @param [in] nrow: number of rows.
  !> @param [in] ind_row: indices i
  !> @param [in]  x: x
  !> @param [out]  y: y
  !-----------------------------------------------------------------------
  subroutine CSR_mat_t_mv_select_row(self, nrow, id_row, x, y)
    class(CSR_mat_t), intent(in)::self 
    integer, intent(in) :: nrow,  id_row(nrow)
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: y(nrow)
    integer :: i, irow, i1, i2, j
    y(:)=0.0_dp
    do j=1, nrow
       irow=id_row(j)
       i1=self%row_shift(irow)
       i2=self%row_shift(irow+1)-1
       do i=i1, i2
          y(j)= y(j)+self%val(i)*x(self%icol(i))
       end do
    end do
  end subroutine CSR_mat_t_mv_select_row

!  subroutine print(self)
!    class(CSR_mat_t), intent(in) :: self
!    print *, "icol:", self%icol
!    print *, "row_shift:", self%row_shift
!    print *, "val:", self%val
!  end subroutine print

end module m_spmat_csr

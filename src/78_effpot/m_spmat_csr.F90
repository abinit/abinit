#include "abi_common.h"
module m_spmat_csr
  use defs_basis  
  use m_xmpi
  use m_spmat_base
  use m_mpi_scheduler
  implicit none
  private
  !!----------- CSR ------------------------
  ! CSR sparse matrix
  ! nnz: number of non-zeros.
  ! icol : column index of entries .size:nnz
  ! row_shift: row_shift(irow) to row_shift(irow+1)-1  are the index
  !            of column and values for entries. size : nrow+1
  ! val: values of non-zero entries. size(nnz)
  !type, public, extends(base_mat_t) :: CSR_mat_t
     type, public, extends(base_mat_t) :: CSR_mat_t
     integer :: nnz
     integer, allocatable :: icol(:), row_shift(:)
     real(dp), allocatable:: val(:)
     type(mpi_scheduler_t) :: mpi_scheduler
   contains
     procedure :: initialize => csr_mat_t_initialize
     procedure :: finalize=> csr_mat_t_finalize
     procedure :: mv => csr_mat_t_mv
     procedure :: sync => csr_mat_t_sync
     procedure :: mv_mpi => csr_mat_t_mv_mpi
     procedure :: mv_select_row =>csr_mat_t_mv_select_row
  end type CSR_mat_t
contains

  ! COO matrix
  subroutine CSR_mat_t_initialize(self, nrow, ncol, nnz, icol, row_shift, val)

    class(CSR_mat_t), intent(inout) :: self
    integer, intent(in) :: nnz, nrow, ncol
    ! i: col number of each entry
    ! j: first 0, n1, n1+n2, ...
    ! val(irow, irow+1) are the values of entries in row irow.
    integer , intent(in), optional :: icol(:), row_shift(:)
    real(dp), intent(in), optional :: val(:)
    integer :: ierr, iproc

    !call MPI_COMM_RANK(MPI_COMM_WORLD, iproc,ierr)
    iproc=xmpi_comm_rank(xmpi_world)
    if (iproc/=0) then
       print *, "This function should be only used on root node"
    end if

    self%nrow=nrow
    self%ncol=ncol
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

  end subroutine CSR_mat_t_initialize

  subroutine CSR_mat_t_sync(self)
    class (csr_mat_t), intent(inout) :: self
    integer :: ierr, iproc
    !call mpi_comm_rank(MPI_COMM_WORLD, iproc, ierr)
    iproc=xmpi_comm_rank(xmpi_world)

    call xmpi_barrier(xmpi_world)
    !call mpi_barrier(MPI_COMM_WORLD, ierr)
    !call mpi_bcast(self%ncol, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call mpi_bcast(self%nrow, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call mpi_bcast(self%nnz, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    call xmpi_bcast(self%ncol, 0, xmpi_world, ierr)
    call xmpi_bcast(self%nrow, 0, xmpi_world, ierr)
    call xmpi_bcast(self%nnz, 0, xmpi_world, ierr)
    call self%mpi_scheduler%initialize(self%nrow, xmpi_world)
    if (.not. self%mpi_scheduler%irank==0) then
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
    !call mpi_bcast(self%icol, self%nnz, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call mpi_bcast(self%row_shift, self%nrow+1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    !call mpi_bcast(self%val, self%nnz, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    call xmpi_bcast(self%icol, 0, xmpi_world, ierr)
    call xmpi_bcast(self%row_shift, 0, xmpi_world, ierr)
    call xmpi_bcast(self%val, 0, xmpi_world, ierr)

  end subroutine CSR_mat_t_sync


  subroutine CSR_mat_t_finalize(self)
    class(CSR_mat_t), intent(inout) :: self
    self%ncol=0
    self%nrow=0
    self%nnz=0
    if(allocated(self%icol)) then
       ABI_DEALLOCATE(self%icol)
    endif
    if(allocated(self%row_shift)) then
       ABI_DEALLOCATE(self%row_shift)
    endif
    if(allocated(self%val)) then
       ABI_DEALLOCATE(self%val)
    endif
  end subroutine CSR_mat_t_finalize



  subroutine CSR_mat_t_mv(self, x, b)

    class(CSR_mat_t), intent(in):: self
    real(dp), intent(in) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
    integer::irow, i1, i2, i
    b(:)=0.0d0
    !!$OMP PARALLEL DO private(i, i1, i2)
    do irow=1, self%nrow
        i1=self%row_shift(irow)
        i2=self%row_shift(irow+1)-1
        do i=i1, i2
            b(irow)=b(irow)+ self%val(i)*x(self%icol(i))
        end do
    enddo
    !!$OMP END PARALLEL DO
  end subroutine CSR_mat_t_mv


  subroutine CSR_mat_t_mv_mpi(self, x, b)
    class(CSR_mat_t), intent(in) :: self
    real(dp), intent(inout) :: x(self%ncol)
    real(dp), intent(out) :: b(self%nrow)
    !real(dp):: my_b(self%nrow)
    integer :: ierr, irow,  i1, i2, i
    !call mpi_bcast(x, self%ncol, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call xmpi_bcast(x, 0, xmpi_world, ierr)
    !my_b(:)=0.0_dp
    b(:)=0.0_dp

    do irow= self%mpi_scheduler%get_istart(), self%mpi_scheduler%get_iend()
       i1=self%row_shift(irow)
       i2=self%row_shift(irow+1)-1
       do i=i1, i2
          b(irow)=b(irow)+ self%val(i)*x(self%icol(i))
       end do
    enddo
    ! TODO : use gather instead of reduce.
    !call mpi_reduce(my_b, b, self%nrow, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
    call xmpi_sum_master(b, 0, xmpi_world, ierr )
  end subroutine CSR_mat_t_mv_mpi

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

end module m_spmat_csr

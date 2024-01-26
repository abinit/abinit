program have_mpi_buggy_interfaces

  use mpi

  implicit none

  integer :: comm,ierr,ival,isum
  real*8 :: xval,xsum

  call mpi_init(ierr)
  call mpi_allreduce(ival,isum,1,MPI_INTEGER,MPI_SUM,comm,ierr)
  call mpi_allreduce(xval,xsum,1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call mpi_allgather(xval,1,MPI_DOUBLE_PRECISION,xsum,1,MPI_DOUBLE_PRECISION,comm,ierr)
  call mpi_finalize(ierr)

end program have_mpi_buggy_interfaces

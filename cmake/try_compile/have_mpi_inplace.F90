program have_mpi_inplace

  use mpi

  implicit none

  integer, parameter :: ii = MPI_INTEGER16
  integer :: comm,ierr,counts(3),displs(3)
  real*8 :: xval(5)
 
  call mpi_init(ierr)
  call mpi_allreduce([MPI_IN_PLACE],xval,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call mpi_allgather([MPI_IN_PLACE],1,MPI_DOUBLE_PRECISION,xval,5,MPI_DOUBLE_PRECISION,comm,ierr)
  call mpi_allgatherv([MPI_IN_PLACE],1,MPI_DOUBLE_PRECISION,xval,counts,displs,MPI_DOUBLE_PRECISION,comm,ierr)
  call mpi_finalize(ierr)

end program have_mpi_inplace

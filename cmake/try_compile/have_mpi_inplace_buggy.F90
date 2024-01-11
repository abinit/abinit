program have_mpi_inplace_buggy

  use mpi

  implicit none

  integer :: comm,ierr,counts(3),displs(3),in_place(1)
  real*8 :: xval(5)
 
  call mpi_init(ierr)
  in_place(1)=MPI_IN_PLACE
  call mpi_allreduce(in_place,xval,5,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
  call mpi_allgather(in_place,1,MPI_DOUBLE_PRECISION,xval,5,MPI_DOUBLE_PRECISION,comm,ierr)
  call mpi_allgatherv(in_place,1,MPI_DOUBLE_PRECISION,xval,counts,displs,MPI_DOUBLE_PRECISION,comm,ierr)
  call mpi_finalize(ierr)

end program have_mpi_inplace_buggy

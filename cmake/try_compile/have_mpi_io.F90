program have_mpi_io

  use mpi

  implicit none

  integer :: ierr, fh, rank, buf(10)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)

  call MPI_File_open(MPI_COMM_WORLD, 'test.out', MPI_MODE_CREATE + MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)
  if (rank==0) then
    call MPI_File_write(fh, buf, 10, MPI_INTEGER, MPI_STATUS_IGNORE, ierr)
  end if
  call MPI_File_close(fh, ierr)

  call MPI_Finalize(ierr)

end program have_mpi_io

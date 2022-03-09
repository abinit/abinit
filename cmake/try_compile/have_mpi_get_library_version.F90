program have_fc_mpi_integer16

  use mpi

  implicit none

  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: info
  integer :: ilen, ierr
  call mpi_init(ierr)
  call mpi_get_library_version(info, ilen, ierr)
  call mpi_finalize(ierr)

end program have_fc_mpi_integer16

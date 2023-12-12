program have_netcdf_fortran_mpi

  use mpi
  use netcdf

  implicit none

  integer :: process_Rank, size_Of_Cluster, ierror
  integer :: ierr, ncid
  size_Of_Cluster = 1

  call MPI_INIT(ierror)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size_Of_Cluster, ierror)
  call MPI_COMM_RANK(MPI_COMM_WORLD, process_Rank, ierror)

  ierr = nf90_create("conftest.nc", ior(NF90_NETCDF4, NF90_MPIPOSIX), &
    ncid, comm=MPI_COMM_WORLD, info=MPI_INFO_NULL)

  ! DEBUG
  !write(*,*) "JMB netcdf-f nfcreate -> ierr= ", ierr

  call MPI_FINALIZE(ierror)

  if(ierr /= 0) stop 1;

end program have_netcdf_fortran_mpi

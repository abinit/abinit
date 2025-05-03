program have_mpi3

  use mpi
  use iso_c_binding

  implicit none

  integer :: win, disp_unit, ierr, comm
  type(c_ptr) :: baseptr
  integer(kind=MPI_ADDRESS_KIND) :: my_size, address

  call MPI_Init(ierr)

  my_size = 0; disp_unit = 1
  call MPI_WIN_ALLOCATE_SHARED(my_size, disp_unit, MPI_INFO_NULL, comm, address, win, ierr)
  call MPI_WIN_SHARED_QUERY(win, 0, my_size, disp_unit, baseptr, ierr)

  call MPI_Finalize(ierr)

end program have_mpi3

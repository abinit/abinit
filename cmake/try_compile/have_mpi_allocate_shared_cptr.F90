program have_mpi_allocate_shared_cptr

  use mpi
  use iso_c_binding

  implicit none

  integer :: win, disp_unit, ierr, comm
  type(c_ptr) :: baseptr
  integer(kind=MPI_ADDRESS_KIND) :: my_size

  call MPI_Init(ierr)

  my_size = 0; disp_unit = 1
  ! Here we pass a c_ptr although the "standard" API defined in mpi module expects integer(MPI_ADDRESS_KIND).
  ! For real world applications we need an API that supports C_PTR so that we can convert to a Fortran pointer.
  ! This problem was fixed in mpi_f08 but Abinit is not ready for that.
  ! See also https://github.com/pmodels/mpich/issues/2659
  call MPI_WIN_ALLOCATE_SHARED(my_size, disp_unit, MPI_INFO_NULL, comm, baseptr, win, ierr)
  call MPI_WIN_SHARED_QUERY(win, 0, my_size, disp_unit, baseptr, ierr)

  call MPI_Finalize(ierr)

end program have_mpi_allocate_shared_cptr

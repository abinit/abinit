program write_irreg_array_ok

  use mpi

  implicit none

  integer, parameter :: nb = 2

  integer :: fh, ierr, nboct_int, nboct_dp
  integer :: bufsize, myrank, nbproc, ii, jj, kk
  integer :: my_status(mpi_status_size)
  integer(kind=mpi_offset_kind) :: offset
  integer :: blockl(2+nb), blockd(2+nb)
  integer :: filetype, blockt(2+nb)
  double precision :: buf(nb)
  double precision, allocatable :: bufftot(:)
  character(len=500) :: errstrg

  ! Initializations
  call mpi_init(ierr)
  call printerr(ierr, "mpi_init")
  call mpi_comm_rank(mpi_comm_world, myrank, ierr)
  call printerr(ierr, "mpi_comm_rank")
  call mpi_comm_size(mpi_comm_world, nbproc, ierr)
  call printerr(ierr, "mpi_comm_size")
  call mpi_type_size(mpi_integer, nboct_int, ierr)
  call printerr(ierr, "mpi_type_size 1")
  call mpi_type_size(mpi_double_precision, nboct_dp, ierr)
  call printerr(ierr, "mpi_type_size 2")

  ! Build a derived datatype (as a mask for the file access)
  blockl(1:nb+2) = (/1, 1, 1, 1/)
  blockt(1:nb+2) = &
&   (/mpi_lb, mpi_double_precision, mpi_double_precision, mpi_ub/)
  if ( myrank == 0 ) blockd(1:nb+2) = (/0, 0       , 2*nboct_dp, 4*nboct_dp/)
  if ( myrank == 1 ) blockd(1:nb+2) = (/0, nboct_dp, 3*nboct_dp, 4*nboct_dp/)
  call mpi_type_struct(nb+2, blockl, blockd, blockt, filetype, ierr)
  call printerr(ierr, "mpi_type_struct")
  call mpi_type_commit(filetype, ierr)
  call printerr(ierr, "mpi_type_commit")

  ! Store 4 numbers in writing buffer
  if ( myrank == 0 ) then
    buf(1) = 0.
    buf(2) = 2.
  end if
  if ( myrank == 1 ) then
    buf(1) = 1.
    buf(2) = 3.
  endif

  ! Write the 4 numbers using a "view"
  offset  = 0
  bufsize = 2
  call mpi_file_open(mpi_comm_world, "conftest.mpi", &
&   mpi_mode_wronly+mpi_mode_create, mpi_info_null, fh, ierr)
  call printerr(ierr, "mpi_file_open 1")
  call mpi_file_set_view(fh, offset, mpi_double_precision, filetype, &
&   'native', mpi_info_null, ierr)
  call printerr(ierr, "mpi_file_set_view")
  call mpi_file_write_all(fh, buf, bufsize, mpi_double_precision, &
&   my_status, ierr)
  call printerr(ierr, "mpi_file_write_all")
  call mpi_file_close(fh, ierr)
  call printerr(ierr, "mpi_file_close 1")
  call mpi_barrier(mpi_comm_world, ierr)
  call printerr(ierr, "mpi_barrier")

  ! Read the 4 numbers (proc 0 only)
  call mpi_file_open(mpi_comm_world, "conftest.mpi", mpi_mode_rdonly, &
&   mpi_info_null, fh, ierr)
  call printerr(ierr, "mpi_file_open 2")
  if ( myrank == 0 ) then
    offset = 0
    allocate(bufftot(nbproc*bufsize))
    bufftot(:) = 0.d0
    call mpi_file_read_at(fh, offset, bufftot, nbproc*bufsize, &
&     mpi_double_precision, my_status, ierr)
    call printerr(ierr, "mpi_file_read_at")
    kk = 0
    do ii = 1, nbproc
      do jj = 1, bufsize
        kk = kk+1
        print *, 'i = ', ii, 'j =', jj, 'buff =', bufftot(kk)
      end do
    end do
    deallocate(bufftot)
  endif
  call mpi_file_close(fh, ierr)
  call printerr(ierr, "mpi_file_close 2")

  ! End of main program
  call mpi_finalize(ierr)
  call printerr(ierr, "mpi_finalize")

contains

  subroutine printerr(ierror, routstrg)

    integer, intent(in) :: ierror
    character(len=*), intent(in) :: routstrg
    integer :: ierr, ilen

    if ( ierror == 0 ) return
    call mpi_error_string(ierror, errstrg, ilen, ierr)
    print *, "error in ", routstrg
    print *, "   ", trim(errstrg)
    stop

  end subroutine printerr

end program write_irreg_array_ok

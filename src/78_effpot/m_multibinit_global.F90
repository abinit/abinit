  ! global
  module m_multibinit_global
    use defs_basis
    use m_random_xoroshiro128plus, only: rng_t
    use m_xmpi
    implicit none
    integer, parameter :: master =0
    logical :: iam_master =.False.
    integer :: my_rank, comm, nproc, ierr
  contains
    subroutine init_multibinit_global()
      comm = xmpi_world
      nproc = xmpi_comm_size(comm)
      my_rank = xmpi_comm_rank(comm)
      iam_master = (my_rank == master)
    end subroutine init_multibinit_global

    subroutine test_bcast()
       real(dp) :: x
       if (iam_master) call random_number(x) 
       call xmpi_bcast(x, master, comm, ierr)
       print *, "myrank :", my_rank, x
    end subroutine test_bcast
  end module m_multibinit_global

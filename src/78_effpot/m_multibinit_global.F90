  ! global
  module m_multibinit_global
    use m_random_xoroshiro128plus, only: rng_t
    use m_xmpi
    implicit none
    integer, parameter :: master =0
    logical, save :: iam_master =.False.
    integer, save :: my_rank, comm, nproc, ierr
  contains
    subroutine init_multibinit_global()
      comm = xmpi_world
      nproc = xmpi_comm_size(comm)
      my_rank = xmpi_comm_rank(comm)
      iam_master = (my_rank == master)
    end subroutine init_multibinit_global
  end module m_multibinit_global

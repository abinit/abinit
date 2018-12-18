#include "abi_common.h"
module m_mpi_scheduler
  !use mpi
  use defs_basis
  use m_xmpi
  implicit none
  private
  type, public :: mpi_scheduler_t
     integer :: nproc, ntasks, irank, comm
     integer,  allocatable :: istart(:), iend(:), ntask_proc(:)
   contains
     procedure :: initialize => mpi_scheduler_t_initialize
     procedure :: finalize => mpi_scheduler_t_finalize
     procedure :: get_iproc => mpi_scheduler_t_get_iproc
     procedure :: get_istart
     procedure :: get_iend
     procedure :: get_ntask_proc
     procedure :: gatherv_dp1d
     procedure :: gatherv_dp2d
  end type mpi_scheduler_t

contains
  subroutine mpi_scheduler_t_initialize(self, ntasks, comm)
    class(mpi_scheduler_t), intent(inout) :: self
    integer, intent(in) :: ntasks
    integer, intent(in) :: comm
    integer :: i,nmore, n, ierr

    !call MPI_COMM_SIZE(comm, self%nproc, ierr)
    self%nproc = xmpi_comm_size(comm)
    !call MPI_COMM_RANK(comm, self%iproc, ierr)
    self%irank=xmpi_comm_rank(comm)
    self%comm=comm
    self%ntasks = ntasks
    !call mpi_bcast(self%ntasks, 1, MPI_INTEGER, 0, comm, ierr)
    call xmpi_bcast(self%ntasks, 0, comm, ierr )
    if (.not. allocated(self%istart)) then
       ABI_ALLOCATE(self%istart, (self%nproc))
    end if

    if (.not. allocated(self%iend)) then
       ABI_ALLOCATE(self%iend, (self%nproc) )
    end if

    if (.not. allocated(self%ntask_proc)) then
       ABI_ALLOCATE(self%ntask_proc, (self%nproc))
    endif

    ! number of procs which has one more
    nmore=mod(self%ntasks, self%nproc)

    if(nmore==0) then
       n=(self%ntasks-nmore)/self%nproc
       do i = 1, self%nproc
          self%istart(i)= 1+(i-1)*n
          self%iend(i)=i*n
          self%ntask_proc(i)=n
       end do
    else
       n=(self%ntasks-nmore)/self%nproc
       do i = 1, nmore
          self%ntask_proc(i)=n+1
          self%istart(i)= 1+(i-1) *(n+1)
          self%iend(i)=self%istart(i)+self%ntask_proc(i)-1
       end do

       do i = nmore+1, self%nproc
          self%ntask_proc(i)=n
          self%istart(i)= self%iend(i-1)+1
          self%iend(i)=self%istart(i)+self%ntask_proc(i)-1
       end do
    end if

  end subroutine mpi_scheduler_t_initialize



  ! find the proc id of which the index of task belong to
  function mpi_scheduler_t_get_iproc(self, i) result(iproc)
    class(mpi_scheduler_t), intent(in) :: self
    integer, intent(in) :: i
    integer :: iproc
    iproc=i/(self%ntasks/self%nproc)
  end function mpi_scheduler_t_get_iproc

  function get_istart(self, rank) result(i)
    class(mpi_scheduler_t), intent(in) :: self
    integer, optional, intent(in):: rank
    integer :: i, r
    if (present(rank)) then
       r=rank
    else
       r=self%irank
    end if
    i=self%istart(r+1)
  end function get_istart

  function get_iend(self, rank) result(i)
    class(mpi_scheduler_t), intent(in) :: self
    integer, optional, intent(in):: rank
    integer :: i, r
    if (present(rank)) then
       r=rank
    else
       r=self%irank
    end if
    i=self%iend(r+1)
  end function get_iend

  
  function get_ntask_proc(self, rank) result(i)
    class(mpi_scheduler_t), intent(in) :: self
    integer, optional, intent(in):: rank
    integer :: i, r
    if (present(rank)) then
       r=rank
    else
       r=self%irank
    end if
    i=self%ntask_proc(r+1)
  end function get_ntask_proc

  subroutine gatherv_dp1d(self, data)
    class(mpi_scheduler_t), intent(inout) :: self
    real(dp), intent(inout) :: data(self%ntasks)
    integer :: ierr
    call xmpi_gatherv(data(self%get_istart(): self%get_iend()), &
         & self%get_ntask_proc(), &
         & data ,&
         & self%ntask_proc, &
         & self%istart-1, &
         & 0, self%comm, ierr )
  end subroutine gatherv_dp1d

  subroutine gatherv_dp2d(self, data, nrow)
    class(mpi_scheduler_t), intent(inout) :: self
    integer, intent(in) :: nrow
    real(dp), intent(inout) :: data(nrow,self%ntasks)
    integer :: ierr
    call xmpi_gatherv(data(:,self%get_istart():self%get_iend()), &
         & self%get_ntask_proc()*nrow, &
         & data, &
         & self%ntask_proc*nrow, &
         & (self%istart-1)*nrow, &
         & 0, self%comm, ierr)

  end subroutine gatherv_dp2d


  subroutine mpi_scheduler_t_finalize(self)
    class(mpi_scheduler_t), intent(inout) :: self
    if (allocated(self%istart)) then
       ABI_DEALLOCATE(self%istart)
    endif
    if (allocated(self%iend)) then
       ABI_DEALLOCATE(self%iend)
    endif
    if (allocated(self%ntask_proc)) then
       ABI_DEALLOCATE(self%ntask_proc)
    endif
  end subroutine mpi_scheduler_t_finalize


end module m_mpi_scheduler

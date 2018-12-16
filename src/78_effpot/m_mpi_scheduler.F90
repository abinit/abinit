#include "abi_common.h"
module m_mpi_scheduler
  !use mpi
  use m_xmpi
  implicit none

  type, public :: mpi_scheduler_t
     integer :: nproc, ntasks, iproc
     integer,  allocatable :: istart(:), iend(:), ntask_proc(:)
   contains
     procedure :: initialize => mpi_scheduler_t_initialize
     procedure :: finalize => mpi_scheduler_t_finalize
     procedure :: get_iproc => mpi_scheduler_t_get_iproc
  end type mpi_scheduler_t

contains
  subroutine mpi_scheduler_t_initialize(self, ntasks, comm)
    class(mpi_scheduler_t), intent(inout) :: self
    integer, intent(inout) :: ntasks
    integer, intent(in) :: comm
    integer :: i,nmore, n, ierr

    !call MPI_COMM_SIZE(comm, self%nproc, ierr)
    self%nproc = xmpi_comm_size(comm)
    !call MPI_COMM_RANK(comm, self%iproc, ierr)
    self%iproc=xmpi_comm_rank(comm)
    self%ntasks = ntasks
    !call mpi_bcast(self%ntasks, 1, MPI_INTEGER, 0, comm, ierr)
    call xmpi_bcast(self%ntasks, 0, comm, ierr )
    if (.not. allocated(self%istart)) then
       ABI_ALLOCATE(self%istart, (0:self%nproc-1))
    end if

    if (.not. allocated(self%iend)) then
       ABI_ALLOCATE(self%iend, (0:self%nproc-1) )
    end if

    if (.not. allocated(self%ntask_proc)) then
       ABI_ALLOCATE(self%ntask_proc, (0:self%nproc-1))
    endif

    ! number of procs which has one more 
    nmore=mod(self%ntasks, self%nproc)

    if(nmore==0) then
       n=(self%ntasks-nmore)/self%nproc
       do i = 0, self%nproc-1
          self%istart(i)= 1+i*n
          self%iend(i)=(i+1)*n
          self%ntask_proc(i)=n
       end do
    else
       n=(self%ntasks-nmore)/self%nproc 
       do i = 0, nmore-1
          self%ntask_proc(i)=n+1
          self%istart(i)= 1+i*(n+1)
          self%iend(i)=self%istart(i)+self%ntask_proc(i)-1
       end do

       do i = nmore, self%nproc-1
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


  subroutine mpi_scheduler_t_finalize(self)
    class(mpi_scheduler_t), intent(inout) :: self
    if (allocated(self%istart)) then
       ABI_DEALLOCATE(self%istart)
    endif
    if (allocated(self%iend)) then
       ABI_DEALLOCATE(self%iend)
    endif
    if (allocated(self%istart)) then
       ABI_DEALLOCATE(self%iend)
    endif
  end subroutine mpi_scheduler_t_finalize


end module m_mpi_scheduler

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mpi_scheduler
!! NAME
!! m_mpi_scheduler
!!
!! FUNCTION
!! This module contains the mpi scheduler for spin dynamics
!! It provide the function to assign site to mpi nodes, and methods for scattering (TODO),
!! and gathering data from nodes.
!!
!! Datatypes:
!!
!! * mpi_scheduler_t
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2019 ABINIT group (hexu)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"


module m_mpi_scheduler
  use defs_basis
  use m_xmpi
  use m_errors
  use m_abicore
 
  implicit none
!!***

  private
  type, public :: mpi_scheduler_t
     integer :: nproc, ntasks, irank, comm, istart, iend, ntask, nblock
     ! ntasks:  total number of tasks
     ! istart: first task id in this proc
     ! iend: last task id in this proc
     ! ntask: number of tasks in this proc
     ! nblock: sometimes it is useful to group task into blocks(eg. each spin has 3 components and its better to put them together.)
     integer,  allocatable :: istart_list(:), iend_list(:), ntask_list(:)
     ! istart_list: istart for all nodes
     ! iend_list: iend for all nodes
     ! ntask_list: ntast for all nodes
   contains
     procedure :: initialize => mpi_scheduler_t_initialize
     procedure :: finalize => mpi_scheduler_t_finalize
     procedure :: get_iproc => mpi_scheduler_t_get_iproc
     procedure :: get_istart
     procedure :: get_iend
     procedure :: get_ntask
     procedure :: gatherv_dp1d ! helper function to gather 1d real(dp) array from nodes.
     procedure :: gatherv_dp2d ! helper function to gather 2d real(dp) array from nodes.
  end type mpi_scheduler_t

contains
  subroutine mpi_scheduler_t_initialize(self, ntasks, comm, nblock)
    ! assign ntasks to ranks in mpi comm.
    ! ntask: number of tasks
    ! nblock: number of subtask per task. TODO: should improve the naming.
    class(mpi_scheduler_t), intent(inout) :: self
    integer, intent(in) :: ntasks
    integer, intent(in) :: comm
    integer, optional, intent(in):: nblock
    integer :: i,nmore, n, ierr

    if( present(nblock)) then
       self%nblock=nblock
    else
       self%nblock=1
    end if
    call xmpi_bcast(self%nblock, 0, comm, ierr)

    !call MPI_COMM_SIZE(comm, self%nproc, ierr)
    self%nproc = xmpi_comm_size(comm)
    !call MPI_COMM_RANK(comm, self%iproc, ierr)
    self%irank=xmpi_comm_rank(comm)
    self%comm=comm
    self%ntasks = ntasks
    call xmpi_bcast(self%ntasks, 0, comm, ierr )
    if (.not. allocated(self%istart_list)) then
       ABI_ALLOCATE(self%istart_list, (self%nproc))
    end if

    if (.not. allocated(self%iend_list)) then
       ABI_ALLOCATE(self%iend_list, (self%nproc) )
    end if

    if (.not. allocated(self%ntask_list)) then
       ABI_ALLOCATE(self%ntask_list, (self%nproc))
    endif

    ! number of procs which has one more
    nmore=mod(self%ntasks, self%nproc)

    if(nmore==0) then
       n=(self%ntasks-nmore)/self%nproc
       do i = 1, self%nproc
          self%istart_list(i)= 1+(i-1)*n
          self%iend_list(i)=i*n
          self%ntask_list(i)=n
       end do
    else
       n=(self%ntasks-nmore)/self%nproc
       do i = 1, nmore
          self%ntask_list(i)=n+1
          self%istart_list(i)= 1+(i-1) *(n+1)
          self%iend_list(i)=self%istart_list(i)+self%ntask_list(i)-1
       end do

       do i = nmore+1, self%nproc
          self%ntask_list(i)=n
          self%istart_list(i)= self%iend_list(i-1)+1
          self%iend_list(i)=self%istart_list(i)+self%ntask_list(i)-1
       end do
    end if

    do i=1, self%nproc
       self%ntask_list(i)=self%ntask_list(i) * self%nblock
       self%istart_list(i) = (self%istart_list(i)-1) *self%nblock + 1
       self%iend_list(i) = self%iend_list(i)*self%nblock
    end do

    self%istart=self%istart_list(self%irank  + 1)
    self%iend=self%iend_list(self%irank  + 1)
    self%ntask=self%ntask_list(self%irank  + 1)

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
    i=self%istart_list(r+1)
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
    i=self%iend_list(r+1)
  end function get_iend

  
  function get_ntask(self, rank) result(i)
    class(mpi_scheduler_t), intent(in) :: self
    integer, optional, intent(in):: rank
    integer :: i, r
    if (present(rank)) then
       r=rank
    else
       r=self%irank
    end if
    i=self%ntask_list(r+1)
  end function get_ntask

  subroutine gatherv_dp1d(self, data)
    class(mpi_scheduler_t), intent(inout) :: self
    real(dp), intent(inout) :: data(self%ntasks)
    integer :: ierr
    call xmpi_gatherv(data(self%istart: self%iend), &
         & self%ntask, &
         & data ,&
         & self%ntask_list, &
         & self%istart_list-1, &
         & 0, self%comm, ierr )
  end subroutine gatherv_dp1d

  subroutine gatherv_dp2d(self, data, nrow)
    class(mpi_scheduler_t), intent(inout) :: self
    integer, intent(in) :: nrow
    real(dp), intent(inout) :: data(nrow,self%ntasks)
    integer :: ierr
    call xmpi_gatherv(data(:,self%istart:self%iend), &
         & self%ntask*nrow, &
         & data, &
         & self%ntask_list*nrow, &
         & (self%istart_list-1)*nrow, &
         & 0, self%comm, ierr)
  end subroutine gatherv_dp2d


  subroutine mpi_scheduler_t_finalize(self)
    class(mpi_scheduler_t), intent(inout) :: self
    if (allocated(self%istart_list)) then
       ABI_DEALLOCATE(self%istart_list)
    endif
    if (allocated(self%iend_list)) then
       ABI_DEALLOCATE(self%iend_list)
    endif
    if (allocated(self%ntask_list)) then
       ABI_DEALLOCATE(self%ntask_list)
    endif
  end subroutine mpi_scheduler_t_finalize


end module m_mpi_scheduler

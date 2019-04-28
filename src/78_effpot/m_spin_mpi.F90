  !{\src2tex{textfont=tt}}
  !!****m* ABINIT/m_spin_mpi
  !! NAME
  !! m_spin_mpi
  !!
  !! FUNCTION
  !! This module contains the module to assign blocks to mpi procs.
  !! and utilities for mpi.
  !!
  !! Datatypes:
  !!  spin_mpi_t
  !!
  !! Subroutines:
  !!
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

module m_spin_mpi

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  implicit none
  private

  !!***

  type, public :: spin_mpi_t
     integer :: nproc, nspins
     integer, allocatable :: istart(:), iend(:), nspin_proc(:)
   contains
     procedure :: initialize => spin_mpi_t_initialize
     procedure :: get_iproc => spin_mpi_t_get_iproc
     procedure :: finalize => spin_mpi_t_finalize
  end type spin_mpi_t

contains
  subroutine spin_mpi_t_initialize(self, nspins, comm)
    class(spin_mpi_t), intent(inout) :: self
    integer, intent(in) :: nspins, comm
    integer :: i,nlast, n

    self%nspins = nspins
    self%nproc = xmpi_comm_size(comm)
    ABI_ALLOCATE(self%istart, (self%nproc))
    ABI_ALLOCATE(self%iend, (self%nproc))
    ABI_ALLOCATE(self%nspin_proc, (self%nproc))


    n=(self%nspins-nlast)/nspins
    do i = 1, self%nproc
       self%istart(i)= 1+(i-1)*n
       self%iend(i)=i*n
       self%nspin_proc(i)=n
    end do
    self%iend(self%nproc)=nspins
    self%nspin_proc(self%nproc)=n
  end subroutine spin_mpi_t_initialize

  ! find the proc id of which the index of spin belong to
  function spin_mpi_t_get_iproc(self, i) result(iproc)
    class(spin_mpi_t), intent(in) :: self
    integer, intent(in) :: i
    integer :: iproc
    iproc=i/(self%nspins/self%nproc)
  end function spin_mpi_t_get_iproc

  subroutine sync_S(self, Slocal, S)
    class(spin_mpi_t), intent(in) :: self
    real(dp), intent(inout) :: Slocal(:,:), S(:,:)
    ABI_UNUSED((/self%nproc/))
    ABI_UNUSED((/Slocal(1,1), S(1,1)/))
    ! TODO 
  end subroutine sync_S

  subroutine spin_mpi_t_finalize(self)
    class(spin_mpi_t), intent(inout) :: self
    ABI_SFREE(self%istart)
    ABI_SFREE(self%iend)
    ABI_SFREE(self%nspin_proc)
  end subroutine spin_mpi_t_finalize

end module m_spin_mpi

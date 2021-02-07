!!****m* ABINIT/m_spin_supercell
!! NAME
!! m_spin_supercell
!!
!! FUNCTION
!! This module define the spin_supercell_t type, which contains the supercell for spin.
!!
!! Datatypes:
!!  spin_supercell_t
!!
!! Subroutines:
!! 
!!
!! COPYRIGHT
!! Copyright (C) 2001-2020 ABINIT group (hexu)
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

module m_spin_supercell
  use defs_basis
  use m_abicore
  use m_errors
  use m_mpi_scheduler, only: init_mpi_info
  implicit none
  private
!!***

  type, public :: spin_supercell_t
     integer :: nspin
     real(dp), allocatable :: ms(:), pos(:, :), gyro_ratio(:), gilbert_damping(:)
     integer, allocatable :: iatoms(:), ispin_prim(:), rvec(:,:)
   contains
     procedure :: initialize
     procedure :: set
     procedure :: finalize
  end type spin_supercell_t
contains
  subroutine initialize(self, nspin)
    class(spin_supercell_t),intent(inout) :: self
    integer, intent(in) :: nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 



    if (iam_master) self%nspin=nspin
    call xmpi_bcast(self%nspin, master, comm, ierr)

    ABI_MALLOC(self%iatoms, (self%nspin))
    ABI_MALLOC(self%pos, (3,self%nspin) )
    ABI_MALLOC(self%ms, (self%nspin))
    ABI_MALLOC(self%ispin_prim, (self%nspin))
    ABI_MALLOC(self%rvec, (3, self%nspin))
    ABI_MALLOC(self%gyro_ratio, (self%nspin))
    ABI_MALLOC(self%gilbert_damping, (self%nspin) )
  end subroutine initialize

  subroutine set(self, ispin_prim, rvec, iatoms, pos, ms, gyro_ratio, damping)
    class(spin_supercell_t),intent(inout) :: self
    integer, intent(in) :: ispin_prim(:), rvec(:,:)
    integer, intent(in) :: iatoms(:)
    real(dp), intent(in)::pos(:,:), ms(:), gyro_ratio(:), damping(:)
    !Local variables-------------------------------
    integer :: nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


    nspin=self%nspin
    if (iam_master) self%iatoms(:)=iatoms(:)
    call xmpi_bcast(self%iatoms, master, comm, ierr)

    if (iam_master) self%pos(:,:)=pos(:,:)
    call xmpi_bcast(self%pos, master, comm, ierr)

    if (iam_master) self%ms(:)=ms(:)
    call xmpi_bcast(self%ms, master, comm, ierr)

    if (iam_master) self%ispin_prim(:)=ispin_prim(:)
    call xmpi_bcast(self%ispin_prim, master, comm, ierr)

    if (iam_master) self%rvec(:,:)=rvec(:,:)
    call xmpi_bcast(self%rvec, master, comm, ierr)

    ! Defautl gyro_ratio
    if(iam_master)  self%gyro_ratio(:)=gyro_ratio
    call xmpi_bcast(self%gyro_ratio, master, comm, ierr)

    if(iam_master) self%gilbert_damping(:)=damping
    call xmpi_bcast(self%gilbert_damping, master, comm, ierr)
  end subroutine set


  subroutine finalize(self)
    class(spin_supercell_t), intent(inout) :: self
    self%nspin=0
    if (allocated(self%ms))  then
       ABI_FREE(self%ms)
    endif

    if (allocated(self%pos))  then
       ABI_FREE(self%pos)
    endif

    if (allocated(self%iatoms)) then
       ABI_FREE(self%iatoms)
    endif
    if (allocated(self%ispin_prim)) then
       ABI_FREE(self%ispin_prim)
    endif
    if (allocated(self%rvec)) then
       ABI_FREE(self%rvec)
    endif

    if (allocated(self%gyro_ratio)) then
       ABI_FREE(self%gyro_ratio)
    endif

    if (allocated(self%gilbert_damping))  then
       ABI_FREE(self%gilbert_damping)
    endif

  end subroutine finalize

end module m_spin_supercell

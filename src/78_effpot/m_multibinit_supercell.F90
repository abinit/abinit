!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_multibinit_supercell
!! NAME
!! m_multibinit_supercell
!!
!! FUNCTION
!! This module define the supercell_t type, which contains all the information of the supercell.
!!
!! Datatypes:
!!  mb_supercell_t
!!
!! Subroutines:
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

module m_multibinit_supercell
  use defs_basis
  use m_abicore
  use m_errors
  use m_mpi_scheduler, only: init_mpi_info
  use m_xmpi
  implicit none
  private
!!***

  type ,public :: mb_supercell_t
     integer :: sc_matrix(3,3)
     integer :: ncell
     logical :: has_lattice=.False.
     logical :: has_spin=.False.
     logical :: has_lwf=.False.
     logical :: has_electron=.False.
     ! lattice related
     integer :: natom
     integer :: ntypat
     real(dp):: rprimd(3,3)
     integer, allocatable :: typat(:)
     real(dp), allocatable :: xred(:, :)
     ! spin related
     integer :: nspin
     real(dp), allocatable :: ms(:), spin_positions(:, :), gyro_ratio(:), gilbert_damping(:)
     integer, allocatable ::  ispin_prim(:), rvec(:,:)

     ! lwf related
     integer :: nlwf

   contains
     procedure:: initialize
     procedure :: finalize
  end type mb_supercell_t
contains


  subroutine initialize(self, natom, nspin, nlwf)
    class(mb_supercell_t),intent(inout) :: self
    integer, intent(in), optional :: nspin, natom, nlwf
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 



    if (present(nspin)) then
       if (iam_master) self%nspin=nspin
       call xmpi_bcast(self%nspin, master, comm, ierr)
       ABI_ALLOCATE(self%spin_positions, (3,self%nspin) )
       ABI_ALLOCATE(self%ms, (self%nspin))
       ABI_ALLOCATE(self%ispin_prim, (self%nspin))
       ABI_ALLOCATE(self%rvec, (3, self%nspin))
       ABI_ALLOCATE(self%gyro_ratio, (self%nspin))
       ABI_ALLOCATE(self%gilbert_damping, (self%nspin) )
    end if
  end subroutine initialize

  subroutine set(self, ispin_prim, rvec,  pos, ms, gyro_ratio, damping)
    class(mb_supercell_t),intent(inout) :: self
    integer, intent(in) :: ispin_prim(:), rvec(:,:)
    real(dp), intent(in)::pos(:,:), ms(:), gyro_ratio(:), damping(:)
    !Local variables-------------------------------
    integer :: nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


    nspin=self%nspin

    if (iam_master) self%spin_positions(:,:)=pos(:,:)
    call xmpi_bcast(self%spin_positions, master, comm, ierr)

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
    class(mb_supercell_t), intent(inout) :: self
    self%nspin=0
    if (allocated(self%ms))  then
       ABI_DEALLOCATE(self%ms)
    endif

    if (allocated(self%spin_positions))  then
       ABI_DEALLOCATE(self%spin_positions)
    endif

    if (allocated(self%ispin_prim)) then
       ABI_DEALLOCATE(self%ispin_prim)
    endif
    if (allocated(self%rvec)) then
       ABI_DEALLOCATE(self%rvec)
    endif

    if (allocated(self%gyro_ratio)) then
       ABI_DEALLOCATE(self%gyro_ratio)
    endif

    if (allocated(self%gilbert_damping))  then
       ABI_DEALLOCATE(self%gilbert_damping)
    endif

  end subroutine finalize



end module m_multibinit_supercell

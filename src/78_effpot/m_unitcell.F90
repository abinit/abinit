!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_unitcell
!! NAME
!! m_unitcell
!!
!! FUNCTION
!! This module define the m_unitcell type, which defines an atomic structure, with spin, lwf, electron
!!
!! Datatypes:
!!  unitcell_t
!!
!! Subroutines:
!!
!! !! COPYRIGHT !! Copyright (C) 2001-2019 ABINIT group (hexu) !! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"


module m_unitcell
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_mpi_scheduler, only: init_mpi_info
  use m_multibinit_supercell, only: mb_supercell_t
  use m_supercell_maker , only: supercell_maker_t
  use m_multibinit_dataset, only : multibinit_dtset_type
  implicit none

!!***
  private
  ! TODO: use crystal_t
  type, public :: unitcell_lattice_t
  !type, public, extends(crystal_t) :: unitcell_lattice_t
   contains
     procedure :: initialize => latt_initialize
     procedure :: finalize => latt_finalize
     procedure :: fill_supercell=>latt_fill_supercell
  end type unitcell_lattice_t

  type, public :: unitcell_spin_t
     integer :: nspin
     real(dp) :: rprimd(3,3)
     real(dp), allocatable :: ms(:)
     real(dp), allocatable :: gyro_ratio(:)
     real(dp), allocatable :: damping_factor(:)
     real(dp), allocatable :: spin_positions(:,:)
   contains
     procedure :: initialize => spin_initialize
     procedure :: finalize => spin_finalize
     procedure :: fill_supercell=>spin_fill_supercell
  end type unitcell_spin_t

  type, public :: unitcell_lwf_t
     integer :: nlwf
   contains
     procedure :: initialize => lwf_initialize
     procedure :: finalize => lwf_finalize
     procedure :: fill_supercell => lwf_fill_supercell
  end type unitcell_lwf_t

  type ,public :: unitcell_t
     logical :: has_lattice=.False.
     logical :: has_spin=.False.
     logical :: has_lwf=.False.
     logical :: has_electron=.False.
     type(unitcell_lattice_t) :: lattice
     type(unitcell_spin_t) :: spin
     type(unitcell_lwf_t) :: lwf
     ! shared 
   contains
     procedure:: initialize
     procedure :: finalize
     procedure :: set_lattice
     procedure :: set_spin
     procedure :: fill_supercell
  end type unitcell_t

contains

  subroutine initialize(self)
    class(unitcell_t), intent(inout) :: self
  end subroutine initialize

  !TODO: Implement
  subroutine set_lattice(self)
    class(unitcell_t), intent(inout) :: self
    self%has_lattice=.True.
    call self%lattice%initialize()
  end subroutine set_lattice


  subroutine set_spin(self,nspin, ms, spin_positions, gyro_ratio, damping_factor)
    class(unitcell_t) , intent(inout):: self
    integer, intent(in) :: nspin
    real(dp), intent(in) :: ms(nspin), spin_positions(3, nspin), gyro_ratio(nspin), damping_factor(nspin)
    self%has_spin=.True.
    call self%spin%initialize(nspin, ms, spin_positions, gyro_ratio, damping_factor)
  end subroutine set_spin

  subroutine set_lwf(self)
    class(unitcell_t) , intent(inout):: self
    self%has_lwf=.True.
    call self%lwf%initialize()
  end subroutine set_lwf

  subroutine read_from_file(self, params, fnames)
    class(unitcell_t), intent(inout) :: self
    type(multibinit_dtset_type), intent(in) :: params
    character(len=fnlen), intent(in) :: fnames(17)
  end subroutine read_from_file


  subroutine finalize(self)
    class(unitcell_t), intent(inout) :: self
    call self%lattice%finalize()
    call self%spin%finalize()
    call self%lwf%finalize()
  end subroutine finalize

  subroutine fill_supercell(self, sc_maker, supercell)
    class(unitcell_t), intent(inout) :: self
    type(supercell_maker_t), intent(inout) :: sc_maker
    type(mb_supercell_t), intent(inout) :: supercell
    call supercell%initialize()
    if (self%has_lattice) then
       call self%lattice%fill_supercell(sc_maker, supercell)
    endif

    if (self%has_spin) then
       call self%spin%fill_supercell(sc_maker, supercell)
    endif

    if (self%has_lwf) then
       call self%lwf%fill_supercell(sc_maker, supercell)
    endif
  end subroutine fill_supercell

  !=========================================================================
  subroutine latt_initialize(self)
    class(unitcell_lattice_t) :: self
  end subroutine latt_initialize

  subroutine latt_finalize(self)
    class(unitcell_lattice_t) :: self
  end subroutine latt_finalize

  subroutine latt_fill_supercell(self, sc_maker, supercell)
    class(unitcell_lattice_t), intent(inout):: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mb_supercell_t), intent(inout):: supercell
  end subroutine latt_fill_supercell


  !========================= SPIN =================================


  Subroutine spin_initialize(self, nspin, ms, spin_positions, gyro_ratio, damping_factor)
    class(unitcell_spin_t) , intent(inout):: self
    integer, intent(in) :: nspin
    real(dp), intent(in) :: ms(nspin), spin_positions(3, nspin), gyro_ratio(nspin), damping_factor(nspin)
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


    self%nspin=nspin
    call xmpi_bcast(self%nspin, master, comm, ierr)
    ABI_ALLOCATE(self%spin_positions, (3, self%nspin))
    ABI_ALLOCATE(self%ms, (self%nspin))
    ABI_ALLOCATE(self%gyro_ratio, (self%nspin))
    ABI_ALLOCATE(self%damping_factor, (self%nspin))
    if (iam_master) then
       self%ms(:) = ms(:)
       self%spin_positions(:,:)=spin_positions(:,:)
       self%gyro_ratio(:)=gyro_ratio(:)
       self%damping_factor(:)=damping_factor(:)
    endif
    call xmpi_bcast(self%spin_positions, master, comm, ierr)
    call xmpi_bcast(self%ms, master, comm, ierr)
    call xmpi_bcast(self%gyro_ratio, master, comm, ierr)
    call xmpi_bcast(self%damping_factor, master, comm, ierr)
  end subroutine spin_initialize

  subroutine spin_finalize(self)
    class(unitcell_spin_t) :: self
    self%nspin=0
    if (allocated(self%ms)) then
       ABI_DEALLOCATE(self%ms)
    end if
    if (allocated(self%spin_positions)) then
       ABI_DEALLOCATE(self%spin_positions)
    end if
    if (allocated(self%gyro_ratio)) then
       ABI_DEALLOCATE(self%gyro_ratio)
    end if
    if (allocated(self%damping_factor)) then
       ABI_DEALLOCATE(self%damping_factor)
    end if
  end subroutine spin_finalize

  subroutine spin_fill_supercell(self, sc_maker, supercell)
    class(unitcell_spin_t),intent(inout) :: self
    type(supercell_maker_t), intent(inout):: sc_maker
    type(mb_supercell_t), intent(inout) :: supercell
    integer :: i, nspin
    integer :: master, my_rank, comm, nproc, ierr
    logical :: iam_master
    call init_mpi_info(master, iam_master, my_rank, comm, nproc) 


    nspin=sc_maker%ncells*self%nspin
    call supercell%initialize(nspin=nspin)
    if(iam_master) then
       call sc_maker%repeat(self%ms, supercell%ms)
       call sc_maker%repeat(self%gyro_ratio, supercell%gyro_ratio)
       call sc_maker%repeat(self%damping_factor, supercell%gilbert_damping)
       call sc_maker%repeat([(i ,i=1, self%nspin)], supercell%ispin_prim)
       supercell%rprimd(:,:)=sc_maker%sc_cell(self%rprimd)
       call sc_maker%trans_xcart(self%rprimd, self%spin_positions, supercell%spin_positions)
       call sc_maker%rvec_for_each(self%nspin, supercell%rvec)
    end if
    call xmpi_bcast(supercell%nspin, master, comm, ierr)
    call xmpi_bcast(supercell%ms, master, comm, ierr)
    call xmpi_bcast(supercell%gyro_ratio, master, comm, ierr)
    call xmpi_bcast(supercell%gilbert_damping, master, comm, ierr)
    call xmpi_bcast(supercell%ispin_prim, master, comm, ierr)
    call xmpi_bcast(supercell%spin_positions, master, comm, ierr)
    call xmpi_bcast(supercell%rvec, master, comm, ierr)
    call xmpi_bcast(supercell%rprimd, master, comm, ierr)
  end subroutine spin_fill_supercell




  !========================= SPIN =================================
  Subroutine lwf_initialize(self)
    class(unitcell_lwf_t) :: self
  end subroutine lwf_initialize

  subroutine lwf_finalize(self)
    class(unitcell_lwf_t) :: self
  end subroutine lwf_finalize

  subroutine lwf_fill_supercell(self, sc_maker,supercell)
    class(unitcell_lwf_t) :: self
    type(supercell_maker_t):: sc_maker
    type(mb_supercell_t) :: supercell
  end subroutine lwf_fill_supercell


end module m_unitcell

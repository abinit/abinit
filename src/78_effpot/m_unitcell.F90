!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_atoms
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
  use m_multibinit_supercell, only: mb_supercell_t
  use m_supercell_maker , only: supercell_maker_t
  use m_multibinit_dataset, only : multibinit_dtset_type

  implicit none

  private
  ! TODO: use crystal_t
  type, public :: unitcell_lattice_t
   contains
     procedure :: initialize => latt_initialize
     procedure :: finalize => latt_finalize
     procedure :: fill_supercell=>latt_fill_supercell
  end type unitcell_lattice_t

  type, public :: unitcell_spin_t
     integer :: nspin
     real(dp), allocatable :: ms(:)
     real(dp), allocatable :: gyroratios(:)
     real(dp), allocatable :: damping_factors(:)
     real(dp), allocatable :: positions(:,:)
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
     integer :: nlwf
     type(unitcell_lattice_t) :: lattice
     type(unitcell_spin_t) :: spin
     type(unitcell_lwf_t) :: lwf
   contains
     procedure:: initialize
     procedure :: finalize
     procedure :: set_lattice
     procedure :: fill_supercell
  end type unitcell_t

contains

  subroutine initialize(self)
    class(unitcell_t), intent(inout) :: self
  end subroutine initialize

  !TODO: Implement
  subroutine set_lattice(self)
    class(unitcell_t), intent(inout) :: self
  end subroutine set_lattice

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
    integer :: sc_natom
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
  Subroutine spin_initialize(self, nspin, ms, positions, gyroratios, damping_factors)
    class(unitcell_spin_t) , intent(inout):: self
    integer, intent(in) :: nspin
    real(dp), intent(in) :: ms(nspin), positions(3, nspin), gyroratios(nspin), damping_factors(nspin)

    self%nspin=nspin
    ABI_ALLOCATE(self%positions, (3, nspin))
    ABI_ALLOCATE(self%ms, (nspin))
    ABI_ALLOCATE(self%gyroratios, (nspin))
    ABI_ALLOCATE(self%damping_factors, (nspin))

    self%ms(:) = ms(:)
    self%positions(:,:)=positions(:,:)
    self%gyroratios(:)=gyroratios(:)
    self%damping_factors(:)=damping_factors(:)
  end subroutine spin_initialize

  subroutine spin_finalize(self)
    class(unitcell_spin_t) :: self
    self%nspin=0
    if (allocated(self%ms)) then
       ABI_DEALLOCATE(self%ms)
    end if
    if (allocated(self%positions)) then
       ABI_DEALLOCATE(self%positions)
    end if
    if (allocated(self%gyroratios)) then
       ABI_DEALLOCATE(self%gyroratios)
    end if
    if (allocated(self%damping_factors)) then
       ABI_DEALLOCATE(self%damping_factors)
    end if
  end subroutine spin_finalize

  subroutine spin_fill_supercell(self, sc_maker, supercell)
    class(unitcell_spin_t) :: self
    type(supercell_maker_t):: sc_maker
    type(mb_supercell_t) :: supercell
    integer:: sc_nspin
    sc_nspin=sc_maker%ncells*self%nspin

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

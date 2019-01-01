!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_effpot_api
!! NAME
!! m_effpot_api
!!
!! FUNCTION
!! This module contains the manager type, which is a thin layer above ALL 
!! TODO: the structure of this is yet to be discussed
!!
!!
!! Datatypes:
!!
!! * base_manager_t: 
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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

module m_effpot_manager
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_supercell, only: supercell_type
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_effpot_api, only: effpot_list_t
  use m_lattice_effpot, only : lattice_effpot_t
  use m_spin_terms, only : spin_terms_t
  use m_lattice_mover, only : lattice_mover_t
  use m_spin_mover, only : spin_mover_t
  use m_spin_lattice_coupling_effpot, only : spin_lattice_coupling_effpot_t
  implicit none
  private

!!***

  type, public :: base_manager_t
   contains
     procedure :: run => base_manager_run
  end type base_manager_t

  type, public, extends(base_manager_t) :: global_manager_t
     type (lattice_effpot_t) :: lattice_effpot
     type (spin_terms_t) :: spin_effpot
     type (spin_lattice_coupling_effpot_t) :: slc_effpot
     type (effpot_list_t) :: effpots
     type (lattice_mover_t) :: lattice_mover
     type (spin_mover_t) :: spin_mover
     type (multibinit_dtset_type):: params
     type (supercell_type) :: supercell
   contains
     procedure :: initialize => global_manager_initialize
     procedure :: finalize => global_manager_finalize
     procedure :: run_one_step => global_manager_run_one_step
     procedure :: run => global_manager_run
  end type global_manager_t

contains

  subroutine base_manager_run(self)
    class (base_manager_t), intent(inout) :: self
  end subroutine base_manager_run

  subroutine global_manager_initialize(self, params, fnames)
    class (global_manager_t), intent(inout) :: self
    type (multibinit_dtset_type), intent(inout) :: params
    character(*), intent(in) :: fnames(:)
    self%params = params
    call self%lattice_effpot%initialize(params, fnames)
    !call self%spin_effpot%initialize(params, fnames)
    call self%slc_effpot%initialize(params, fnames)

    call self%lattice_effpot%make_supercell(self%supercell)
    call self%spin_effpot%make_supercell(self%supercell)
    call self%slc_effpot%make_supercell(self%supercell)

    call self%effpots%initialize(natoms=self%lattice_effpot%natoms, nspins=self%spin_effpot%nspins)
    call self%effpots%append(self%lattice_effpot)
    call self%effpots%append(self%spin_effpot)
    call self%effpots%append(self%slc_effpot)

    call self%spin_mover%initialize(params, nspins=self%spin_effpot%nspins)
    call self%lattice_mover%initialize(params, fnames )

    call self%lattice_mover%set_initial_state()
    call self%spin_mover%set_initial_state()



  end subroutine global_manager_initialize

  subroutine global_manager_finalize(self)

    class (global_manager_t), intent(inout) :: self
    call self%lattice_effpot%finalize()
    call self%spin_effpot%finalize()
    call self%slc_effpot%finalize()

    call self%lattice_mover%finalize()
    call self%spin_mover%finalize()
  end subroutine global_manager_finalize

  subroutine global_manager_run_one_step(self)
    class (global_manager_t), intent(inout) :: self
    ! TODO: steps with only spin or lattice
    ! one step of lattice
    call self%lattice_mover%run_one_step(self%effpots)
    ! one step of spin
    call self%spin_mover%run_one_step(self%effpots)
    ! update changes to spin-lattice_coupling
    !call self%slc_effpot%set_distortion(displacement=self%current_displacement, strain=self%current_strain)
    !call self%slc_effpot%set_spin(spin=self%current_spin)
  end subroutine global_manager_run_one_step

  subroutine global_manager_run(self)
    class (global_manager_t), intent(inout) :: self
  end subroutine global_manager_run

end module m_effpot_manager

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abstract_potential
!! NAME
!! m_abstract_potential
!!
!! FUNCTION
!! This module contains the base type for all effective potentials. 
!!
!!
!! Datatypes:
!!
!! * abstract_potential_t: defines the base api of effective potentials.
!! * potential_list_t: list of abstract_potential_t, which is essentially a list of pointer to abstract_potential_t
!!    itself is also a effpot type, and its energy, 1st derivative to energy are the sum of all items.
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

module m_abstract_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_supercell, only: mb_supercell_t

  implicit none
!!***
  private
  type ,public :: abstract_potential_t
     ! This is the abstract class of effective potential.
     ! It do the following things:
     ! calculate 0th, 1st,... derivative of energy related, e.g. Force for lattice model
     ! labels for variables
     logical :: has_displacement=.False.
     logical :: has_strain=.False.
     logical :: has_spin=.False.
     logical :: has_lwf = .False.
     logical :: is_null=.True.   ! if is_null, this term does not exist.
     type(mb_supercell_t) ,pointer :: supercell => null()
     !real(dp), allocatable :: ms(:)
     character (len=200) :: label="Abstract Potential"
   contains
     procedure :: set_supercell   ! set_supercell
     procedure :: finalize        ! finalize
     procedure :: set_params      ! parameters from input file
     procedure :: calculate       ! get energy and 1st derivative from input state
     procedure :: get_delta_E     ! calculate energy diffence if one component is changed for Monte carlo algorithm
  end type abstract_potential_t

contains

  subroutine set_supercell(self, supercell)
    class(abstract_potential_t), intent(inout) :: self
    type(mb_supercell_t), target, intent(inout) :: supercell
    self%supercell => supercell
    MSG_ERROR("Every potential should override this set_supercell method to avoid mistake.")
  end subroutine set_supercell

  subroutine finalize(self)
    class(abstract_potential_t), intent(inout) :: self
    self%is_null=.True.
    nullify(self%supercell)
    self%label="Destroyed potential"
  end subroutine finalize

  subroutine set_params(self, params)
    class(abstract_potential_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    MSG_ERROR("Every potential should override this method set_params to avoid mistakes")
  end subroutine set_params

  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy)
    ! This function calculate the energy and its first derivative
    ! the inputs and outputs are optional so that each effpot can adapt to its
    ! own.
    ! In principle, the 1st derivatives are only calculated if asked to (present). However, they can be computed if it is simply convinient to do.
    class(abstract_potential_t), intent(inout) :: self  ! the effpot may save the states.

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    ! if present in input
    ! calculate if required
    MSG_ERROR("calculate not implemented for this effpot.")
  end subroutine calculate

  subroutine get_delta_E(self, S, ispin, Snew, deltaE)
    ! for spin monte carlo
    ! calculate energy difference if one spin is moved.
    class(abstract_potential_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), intent(inout) :: S(:,:),  Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(inout) :: deltaE
    MSG_ERROR("get_delta_E currenlty only implemented for spin potential.")
  end subroutine get_delta_E

end module m_abstract_potential

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

module m_abstract_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_hashtable_strval, only: hash_table_t

  implicit none
!!***
  private


  type ,public :: abstract_potential_t
     ! This is the abstract class of effective potential.
     ! It do the following things:
     ! calculate 0th, 1st,... derivative of energy related, e.g. Force for lattice model
     ! labels for variables. If a new degree of freedom is to be added, a label
     ! should be added here.
     logical :: has_displacement=.False.
     logical :: has_strain=.False.
     logical :: has_spin=.False.
     logical :: has_lwf = .False.
     logical :: is_null=.True.   ! if is_null, this term does not exist.
     ! it is important to set it to True before initialization.
     ! Because it is used as the tag for deallocatign memory. 
     type(mbsupercell_t) ,pointer :: supercell => null()
     ! every supercell potential has a pointer to the supercell,
     ! which could be used for like reference structure.
     character (len=200) :: label="Abstract Potential"  !
     ! the label is used for printing information.
   contains
     procedure :: set_supercell   ! set_supercell
     procedure :: finalize        ! finalize
     procedure :: set_params      ! parameters from input file
     procedure :: calculate       ! get energy and 1st derivative from input state
     procedure :: get_delta_E     ! calculate energy diffence if one component is changed for Monte carlo algorithm
  end type abstract_potential_t

contains

  !----------------------------------------------------------------------
  !> @brief link a supercell to the potential
  !>
  !> @param[in]  supercell: supercell object
  !----------------------------------------------------------------------
  subroutine set_supercell(self, supercell)
    class(abstract_potential_t), intent(inout) :: self
    type(mbsupercell_t), target, intent(inout) :: supercell
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(supercell)
    MSG_ERROR("Every potential should override this set_supercell method to avoid mistake.")
  end subroutine set_supercell

  !----------------------------------------------------------------------
  !> @brief finalize
  !>
  !----------------------------------------------------------------------
  subroutine finalize(self)
    class(abstract_potential_t), intent(inout) :: self
    self%is_null=.True.
    nullify(self%supercell)
    self%label="Destroyed potential"
  end subroutine finalize

  !----------------------------------------------------------------------
  !> @brief set_params: set the parameters from input file parameters
  !>
  !> @param[in]  params: multibinit_dtset_type: from input file
  !----------------------------------------------------------------------
  subroutine set_params(self, params)
    class(abstract_potential_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)
    MSG_ERROR("Every potential should override this method set_params to avoid mistakes")
  end subroutine set_params

  ! hexu comment : which one is better, more general variables,
  !  or one function for each type of var?
  !subroutine set_variables(self, displacements, strain, spin)
  !  class(lattice_api_t), intent(inout) :: self
  !  real(dp), optional, intent(in) :: displacements(:,:), strain(:,:), spin(:,:)
  !end subroutine set_variables

  !subroutine get_1st_deriv(self, force, stress, bfield)
  !  class(lattice_api_t), intent(inout) :: self
  !  real(dp), optional, intent(out) :: force(:,:), stress(:,:), bfield(:,:)
  !end subroutine get_1st_deriv

  ! subroutine set_distortion(self, displacement, strain)
  !   class(abstract_potential_t), intent(inout) :: self
  !   real(dp), optional, intent(in) :: displacement(:,:), strain(:,:)
  !   MSG_ERROR("set_distortion not implemented.")
  ! end subroutine set_distortion

  ! subroutine set_spin(self, spin)
  !   class(abstract_potential_t), intent(inout) :: self
  !   real(dp), optional, intent(in) :: spin
  !   MSG_ERROR("set_spin not implemented.")
  ! end subroutine set_spin

  !----------------------------------------------------------------------
  !> @brief calculate energy and derivatives with given state.
  !> This function calculate the energy and its first derivative
  !> the inputs and outputs are optional so that each effpot can adapt to its
  !> own.
  !> In principle, the 1st derivatives are only calculated if
  !> asked to (present).
  !> However, they can be computed if it is simply convinient to do.
  !> @param[in]  displacement (optional)
  !> @param[in]  strain (optional)
  !> @param[in]  spin (optional)
  !> @param[in]  lwf (optional)
  !> @param[out] force (optional)
  !> @param[out] stress (optional)
  !> @param[out] bfield (optional)
  !> @param[out] energy_table (optional)
  !----------------------------------------------------------------------
  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, &
       & energy, energy_table)
      class(abstract_potential_t), intent(inout) :: self  ! the effpot may save the states.

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    type(hash_table_t),optional, intent(inout) :: energy_table
    ! if present in input
    ! calculate if required
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(force)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(lwf_force)
    ABI_UNUSED_A(energy)
    ABI_UNUSED_A(energy_table)
    MSG_ERROR("calculate not implemented for this effpot.")
  end subroutine calculate

  !----------------------------------------------------------------------
  !> @brief get_delta_E: calculate the energy difference when a given spin
  !> is changed. This is to be used for spin Monte Carlo. Currently the
  !> only supported is the spin model. 
  !>
  !> @param[in]  S: spin of full structure. array of (3, nspin)
  !> @param[in]  ispin: the index of spin changed. integer
  !> @param[in]  snew: the new value of the changed spin. 
  !> @param[out] deltaE: the energy difference
  !----------------------------------------------------------------------
  subroutine get_delta_E(self, S, ispin, Snew, deltaE)
    ! for spin monte carlo
    ! calculate energy difference if one spin is moved.
    class(abstract_potential_t), intent(inout) :: self  ! the effpot may save the states.
    real(dp), intent(inout) :: S(:,:),  Snew(:)
    integer, intent(in) :: ispin
    real(dp), intent(inout) :: deltaE
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(S)
    ABI_UNUSED_A(ispin)
    ABI_UNUSED_A(Snew)
    ABI_UNUSED_A(deltaE)
    MSG_ERROR("get_delta_E not implemented for this effpot.")
  end subroutine get_delta_E

!   subroutine get_energy(self, energy)
!     class(abstract_potential_t), intent(inout) :: self
!     real(dp) , intent(inout) :: energy
!   end subroutine get_energy


!   subroutine get_force(self, force)
!     class(abstract_potential_t), intent(inout) :: self
!     real(dp), intent(out) :: force(:,:)
!   end subroutine get_force

!   subroutine get_stress(self, stress)
!     class(abstract_potential_t), intent(inout) :: self
!     real(dp), intent(out) :: stress(:,:)
!   end subroutine get_stress

!   subroutine get_effective_Bfield(self, spin,bfield)
!     class(abstract_potential_t), intent(in) :: self
!     real(dp), intent(in) :: spin(:,:)
!     real(dp), intent(inout) :: bfield(:,:)
!   end subroutine get_effective_Bfield

end module m_abstract_potential

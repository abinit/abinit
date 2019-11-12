!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abstract_mover
!! NAME
!! m_abstract_mover
!!
!! FUNCTION
!! This module contains the base type for all mover types.
!!
!!
!! Datatypes:
!!
!! * abstract_mover_t: defines the base api of movers.
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

module m_abstract_mover
  use defs_basis
  use m_abicore
  use m_errors

  use m_random_xoroshiro128plus, only: rng_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_abstract_potential, only: abstract_potential_t
  use m_multibinit_cell, only: mbcell_t, mbsupercell_t
  use m_hashtable, only: hash_table_t
  implicit none
!!***

  private

  !-------------------------------------------------------------------!
  !>@brief Abstract_mover_t:
  !>  abstract mover. All the movers should be derived from this.
  !> A mover defines the evolution of a structure.
  !> it will call the potentials to calculate the energy derivative and E
  !> and use that to move the states.
  !
  !-------------------------------------------------------------------!
  type ,public :: abstract_mover_t
     ! This is the abstract class of mover
     ! It do the following things:
     ! calculate d(var)/dt and integrate new var. 
     ! call functions to calculate observables.
     ! interact with hist file.

     ! a pointer to the supercell structure
     type(mbsupercell_t), pointer:: supercell=>null() 
     type(hash_table_t), pointer :: etable=>null()
     ! a label for each mover. For printing out information
     character (len=200) :: label="Abstract Mover"
     ! time step.
     real(dp) :: dt= 0.0
     ! total time
     real(dp) :: total_time =0.0
     ! temperature (NOTE: It is not the real temperature of
     ! the structure, but a preset value, even when it is not
     ! a constant temperature mover!!!)
     real(dp) :: temperature =0.0
     ! time for themalization
     real(dp) :: thermal_time =0.0
     ! A pointer to the random number generator
     ! It is initialized outside the mover (by the manager).
     type(rng_t), pointer :: rng=>null()

   contains
     !procedure:: initialize       ! perhaps each effpot type should have own
     !procedure :: finalize
     procedure :: set_rng          ! set the pointer to an random number generator
     procedure :: set_params       ! set the parameters from input
     procedure :: set_initial_state ! initial state
     procedure:: run_one_step      ! defines how the system evolves.
     !procedure:: run_time
     !procedure:: run_temperature
     procedure :: reset            ! reset the mover
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
     procedure :: set_energy_table
  end type abstract_mover_t

contains

  !-------------------------------------------------------------------!
  ! set_rng:
  ! set the random number generator. The rng is already initialize
  ! outside. 
  !-------------------------------------------------------------------!
  subroutine set_rng(self, rng)
    class(abstract_mover_t), intent(inout) :: self
    type(rng_t), target, intent(in) :: rng
    self%rng => rng
  end subroutine set_rng

  !-------------------------------------------------------------------!
  ! set_parms:
  !  set the variables using the input params
  !-------------------------------------------------------------------!
  subroutine set_params(self, params)
    ! set parameters from input file. (something else, like temperature for MvT calculation?)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(params)
    MSG_ERROR("set_params not implemented for this mover")
  end subroutine set_params

  !----------------------------------------------------------------------
  !> @brief set the pointer to supercell
  !> @param[in]  supercell
  !-------------------------------------------------------------------!
  subroutine set_supercell(self, supercell)
    class(abstract_mover_t), intent(inout) :: self
    type(mbsupercell_t), target :: supercell
    self%supercell=>supercell
  end subroutine set_supercell

  !----------------------------------------------------------------------
  !> @brief set initial state.
  !>
  !> @param[in]  mode: a integer to define the kind of initial state.
  !----------------------------------------------------------------------
  subroutine set_initial_state(self, mode, restart_hist_fname)
    ! set initial positions, spin, etc
    class(abstract_mover_t), intent(inout) :: self
    integer, optional, intent(in) :: mode
    character(len=fnlen), optional, intent(in) :: restart_hist_fname

    MSG_ERROR("set_initial_state not implemented for this mover")
    ABI_UNUSED_A(self)
    ABI_UNUSED(mode)
    ABI_UNUSED(restart_hist_fname)
  end subroutine set_initial_state

  !-------------------------------------------------------------------!
  !> @brief: Run_one_step
  !>   effpot: the potential (which do the calculation of E and dE/dvar)
  !> param[in]: effpot
  !> param[in]: (optional) displacement
  !> param[in]: (optional) strain
  !> param[in]: (optional) spin
  !> param[in]: (optional) lwf
  ! NOTE: No need to pass the variable already saved in the mover.
  !     e.g. For spin mover, do NOT pass the spin to it.
  !    The other variables are only required if there is coupling with
  !    the mover variable.
  !-------------------------------------------------------------------!
  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf)
    ! run one step. (For MC also?)
    class(abstract_mover_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    class(abstract_potential_t), intent(inout) :: effpot
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    MSG_ERROR("run_one_step not implemented for this mover")
  end subroutine run_one_step

  !-------------------------------------------------------------------!
  ! Reset:
  ! Reset the mover.
  ! It is used, when multiple sets of move are needed, e.g. for various
  ! temperature.
  !-------------------------------------------------------------------!
  subroutine reset(self)
    ! reset the state of mover (e.g. counter->0)
    ! so it can be reused.
    class(abstract_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
    MSG_ERROR("reset not implemented for this mover")
  end subroutine reset

  !-------------------------------------------------------------------!
  !> @brief: calc_observables: calculate observables
  !-------------------------------------------------------------------!
  subroutine calc_observables(self)
    ! call functions to calculate observables.
    class(abstract_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
    MSG_ERROR("calc_observables not implemented for this mover")
  end subroutine calc_observables

  !-------------------------------------------------------------------!
  !> @brief Write_hist:
  !   write to hist file.
  !-------------------------------------------------------------------!
  subroutine write_hist(self)
    ! write to hist file
    class(abstract_mover_t), intent(inout) :: self
    ABI_UNUSED_A(self)
    MSG_ERROR("write_hist not implemented for this mover")
  end subroutine write_hist

  !-------------------------------------------------------------------!
  !> @breif: Get_state: the current state
  !-------------------------------------------------------------------!
  subroutine get_state(self, displacement, strain, spin, lwf, ihist)
    ! get the state of the ihist(th) step. ihist can be 0 (current), -1 (last), ... -maxhist..
    !Note that the params are optional so it will be returned only if asked for.
    class(abstract_mover_t), intent(in):: self
    real(dp), optional, intent(inout) :: displacement, strain, spin, lwf
    integer, optional, intent(in):: ihist
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(ihist)
    MSG_ERROR("get_state not implemented for this mover")
  end subroutine get_state



end module m_abstract_mover

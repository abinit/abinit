!!****m* ABINIT/m_lwf_dummy_mover
!! NAME
!! m_lwf_dummy_mover
!!
!! FUNCTION
!! This module contains the lwf Markov chain Monte Carlo functions for lwf mover .
!!
!!
!! Datatypes:
!!
!! * lwf_dummy_mover_t : This mover does not move lwf, instead, it allow the lwf
!!  to be moved by the potential. It is used to monitor the move of lwf when the lattice mover
!!  is actually used.
!!
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
    module m_lwf_dummy_mover
    use defs_basis
    use m_abicore
    use m_errors
    use m_abstract_potential, only: abstract_potential_t
    use m_random_xoroshiro128plus, only: rng_t
    use m_hashtable_strval, only: hash_table_t
    use m_multibinit_dataset, only: multibinit_dtset_type
    use m_lwf_mover, only: lwf_mover_t
    use m_multibinit_cell, only: mbcell_t, mbsupercell_t
    implicit none
!!***
    private
    !----------------------------------------------------------------------
    !> @brief An dummy mover for lwf
    !----------------------------------------------------------------------
    type,public, extends(lwf_mover_t) :: lwf_dummy_mover_t
     contains
       procedure  :: run_one_step
    end type lwf_dummy_mover_t
  contains


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
  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf,  energy_table)
    class(lwf_dummy_mover_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    type(hash_table_t), optional, intent(inout) :: energy_table
    class(abstract_potential_t), intent(inout) :: effpot
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(energy_table)
  end subroutine run_one_step
end module m_lwf_dummy_mover

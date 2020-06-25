!!****m* ABINIT/m_lwf_mover
!! NAME
!! m_lwf_mover
!!
!! FUNCTION
!! This module contains the lwf mover, which controls how the lattice wannier function move.
!!
!!
!! Datatypes:
!!
!! * lwf_mover_t
!!
!! Subroutines:
!!
!! * lwf_mover_t_initialize
!! * lwf_mover_t_run_one_step
!! * lwf_mover_t_run_time
!! * TODO: update this when F2003 documentation format decided.
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
module m_lwf_mover
  use defs_basis
  use m_errors
  use m_abicore
  use m_xmpi
  use m_mpi_scheduler, only: mpi_scheduler_t
  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_random_xoroshiro128plus, only: set_seed, rand_normal_array, rng_t
  use m_abstract_potential, only: abstract_potential_t
  use m_abstract_mover, only: abstract_mover_t
  use m_hashtable_strval, only: hash_table_t
  
  implicit none
  private
  !!***

  type, public, extends(abstract_mover_t) :: lwf_mover_t
   contains
     procedure :: run_one_step
  end type lwf_mover_t

contains
  subroutine run_one_step(self, effpot, displacement, strain, spin, lwf, energy_table)
    ! run one step. (For MC also?)
    class(lwf_mover_t), intent(inout) :: self
    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    class(abstract_potential_t), intent(inout) :: effpot
    type(hash_table_t),optional, intent(inout) :: energy_table
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(effpot)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(energy_table)
    MSG_ERROR("run_one_step not implemented for this mover")
  end subroutine run_one_step


end module m_lwf_mover


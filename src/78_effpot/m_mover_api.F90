!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_mover_api
!! NAME
!! m_mover_api
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

module m_mover_api
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type
  use m_effpot_api, only: effpot_t

  implicit none
!!***

  private
  type ,public :: abstract_mover_t
     ! This is the abstract class of mover
     ! It do the following things:
     ! calculate d(var)/dt and integrate new var. 
     ! call functions to calculate observables.
     ! interact with hist file.
   contains
     !procedure:: initialize       ! perhaps each effpot type should have own 
     !procedure :: finalize
     procedure :: set_params
     procedure :: set_initial_state ! initial state
     procedure:: run_one_step 
     procedure :: reset            ! reset the mover
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
  end type abstract_mover_t

contains

  subroutine set_params(self, params)
    ! set parameters from input file. (something else, like temperature for MvT calculation?)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
    MSG_ERROR("set_params not implemented for this mover")
  end subroutine set_params

  subroutine set_initial_state(self)
    ! set initial positions, spin, etc
    class(abstract_mover_t), intent(inout) :: self
    MSG_ERROR("set_initial_state not implemented for this mover")
  end subroutine set_initial_state

  subroutine run_one_step(self, effpot)
    ! run one step. (For MC also?)
    class(abstract_mover_t), intent(inout) :: self
    ! should we use array of effective potentials so that there can be multiple of them?
    class(effpot_t), intent(inout) :: effpot
    MSG_ERROR("run_one_step not implemented for this mover")
  end subroutine run_one_step

  subroutine reset(self)
    ! reset the state of mover (e.g. counter->0)
    ! so it can be reused.
    class(abstract_mover_t), intent(inout) :: self
    MSG_ERROR("reset not implemented for this mover")
  end subroutine reset

  subroutine calc_observables(self)
    ! call functions to calculate observables.
    class(abstract_mover_t), intent(inout) :: self

    MSG_ERROR("calc_observables not implemented for this mover")
  end subroutine calc_observables

  subroutine write_hist(self)
    ! write to hist file
    class(abstract_mover_t), intent(inout) :: self
    MSG_ERROR("write_hist not implemented for this mover")
  end subroutine write_hist

  subroutine get_state(self, position, strain, spin, ihist)
    ! get the state of the ihist(th) step. ihist can be 0 (current), -1 (last), ... -maxhist..
    !Note that the params are optional so it will be returned only if asked for.
    class(abstract_mover_t), intent(in):: self
    real(dp), optional, intent(inout) :: position, strain, spin
    integer, optional, intent(in):: ihist
    MSG_ERROR("get_state not implemented for this mover")
  end subroutine get_state

end module m_mover_api

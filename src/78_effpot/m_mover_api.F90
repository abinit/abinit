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
     procedure :: run_one_step
     procedure :: run_nstep
     procedure :: reset            ! reset the mover
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
  end type abstract_mover_t

contains

  subroutine initialize(self, params, fnames, effpot)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params  ! read from input file
    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).
    ! hexu comment  It's better to get rid of it and put everything into input file!?)
    class(effpot_t) :: effpot ! may use information from effective potential
  end subroutine initialize

  subroutine finalize(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine finalize

  subroutine set_params(self, params)
    ! set parameters from input file. (something else, like temperature for MvT calculation?)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
  end subroutine set_params

  subroutine set_initial_state(self)
    ! set initial positions, spin, etc
    class(abstract_mover_t), intent(inout) :: self
  end subroutine set_initial_state

  subroutine set_force(self, force, add)
    ! set force. Should not be here but in lattice mover only.
    class(abstract_mover_t), intent(inout) :: self
    real(dp), intent(in) :: force(:,:)
    logical, optional :: add ! add or reset force.
  end subroutine set_force

  subroutine run_one_step(self)
    ! run one step. (For MC also?)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine run_one_step

  subroutine run_nstep(self, nstep)
    class(abstract_mover_t), intent(inout) :: self
    integer, intent(in) :: nstep
  end subroutine run_nstep

  subroutine reset(self)
    ! reset the state of mover (e.g. counter->0)
    ! so it can be reused.
    class(abstract_mover_t), intent(inout) :: self
  end subroutine reset

  subroutine calc_observables(self)
    ! call functions to calculate observables.
    class(abstract_mover_t), intent(inout) :: self
  end subroutine calc_observables

  subroutine write_hist(self)
    ! write to hist file
    class(abstract_mover_t), intent(inout) :: self
  end subroutine write_hist

  subroutine get_state(self, position, strain, spin, ihist)
    ! get the state of the ihist(th) step. ihist can be 0 (current), -1 (last), ... -maxhist..
    class(abstract_mover_t), intent(in):: self
    real(dp), optional, intent(inout) :: position, strain, spin
    integer, optional, intent(in):: ihist
  end subroutine get_state

end module m_mover_api

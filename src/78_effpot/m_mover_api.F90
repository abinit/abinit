#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_mover_api
  use defs_basis
  use m_abicore
  use m_errors

  use m_multibinit_dataset, only: multibinit_dtset_type

  implicit none
  private
  type ,public :: abstract_mover_t
     ! This is the abstract class of mover
     ! It do the following things:
     ! calculate d(var)/dt and integrate new var. 
     ! call functions to calculate observables.
     ! interact with hist file.
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: set_params
     procedure :: set_initial_state ! initial state
     procedure :: run_one_step
     procedure :: run_nstep
     procedure :: calc_observables ! call functions to calculate observables
     procedure :: write_hist       ! write hist file
  end type abstract_mover_t

contains

  subroutine initialize(self, params, fnames)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params  ! read from input file
    character(len=*), intent(in) :: fnames(:)  !  files file (xml, DDB, etc).
    ! hexu comment  It's better to get rid of it and put everything into input file!?)
  end subroutine initialize

  subroutine finalize(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine finalize

  subroutine set_params(self, params)
    class(abstract_mover_t), intent(inout) :: self
    type(multibinit_dtset_type) :: params
  end subroutine set_params

  subroutine set_initial_state(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine set_initial_state

  subroutine run_one_step(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine run_one_step

  subroutine run_nstep(self, nstep)
    class(abstract_mover_t), intent(inout) :: self
    integer, intent(in) :: nstep
  end subroutine run_nstep

  subroutine calc_observables(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine calc_observables

  subroutine write_hist(self)
    class(abstract_mover_t), intent(inout) :: self
  end subroutine write_hist

 



end module m_mover_api

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_multibinit_main2
  use defs_basis
  use defs_abitypes
  use m_build_info
  use m_xmpi
  use m_xomp
  use m_abicore
  use m_errors
 
  use m_effective_potential
  use m_fit_polynomial_coeff
  use m_multibinit_dataset
  use m_multibinit_global
  use m_effective_potential_file
  use m_spin_model, only: spin_model_t
  use m_abihist

  use m_multibinit_manager, only: mb_manager_t

  !use m_generate_training_set, only : generate_training_set
  use m_compute_anharmonics, only : compute_anharmonics
  use m_init10,              only : init10
  use m_fstrings,   only : replace, inupper

  implicit none
contains
subroutine multibinit_main2(filnam)
  character(len=fnlen), intent(inout) :: filnam(17)
  type(mb_manager_t) :: manager
  call manager%run_all(filnam)
end subroutine multibinit_main2

end module m_multibinit_main2

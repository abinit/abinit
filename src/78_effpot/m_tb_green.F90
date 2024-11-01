#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_green_J
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_multibinit_dataset, only: multibinit_dtset_type
  implicit none
  private

  type, public :: green_J_t

   contains
     procedure :: initialize
     procedure :: finalize
  end type green_J_t

contains

  subroutine initialize(self)
    class (green_J_t) :: self
  end subroutine initialize

  subroutine finalize(self)
    class (green_J_t) :: self
  end subroutine finalize


end module m_green_J

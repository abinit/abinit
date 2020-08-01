!!****m* ABINIT/m_lwf_observables
!! NAME
!! m_lwf_observables
!!
!! FUNCTION
!! This module contains the subroutines for utilities of calculating the observables
!!
!!
!! Datatypes:
!! lwf_observable_t: store data to calculate observables
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

module m_lwf_observables

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_observables, only: Cv_t

  implicit none

  private
  !!***

  type, public ::lwf_observables_t
     type(Cv_t) :: Cv
   contains
     procedure :: initialize
     procedure :: finalize
  end type lwf_observables_t
contains
  subroutine initialize(self)
    class(lwf_observables_t) :: self
    ABI_UNUSED_A(self)
  end subroutine initialize

  subroutine finalize(self)
    class(lwf_observables_t) :: self
    ABI_UNUSED_A(self)
  end subroutine finalize

end module m_lwf_observables

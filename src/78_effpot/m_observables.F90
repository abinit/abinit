!!****m* ABINIT/m_observables
!! NAME
!! m_observables
!!
!! FUNCTION
!! This module contains the subroutines for utilities of calculating the observables
!!
!!
!! Datatypes:
!! observable_t: store data to calculate observables
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

module m_observables

  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  implicit none

  private
  !!***
  type, public :: Cv_t
     real(dp):: Cv
  end type Cv_t

end module m_observables

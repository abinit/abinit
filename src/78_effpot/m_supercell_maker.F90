!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_supercell_maker
!! NAME
!! m_supercell_maker
!!
!! FUNCTION
!! This module define the supercell_maker file, which provide functions to help build
!! potentials in supercell.
!!
!! Datatypes:
!!  supercell_maker_t
!!
!! Subroutines:
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

module m_supercell_maker
  use defs_basis
  use m_abicore
  use m_errors
  use m_supercell
  implicit none

  type ,public :: supercell_maker_t
     integer :: sc_matrix(3,3)
     integer :: ncell
   contains
     procedure:: initialize       ! perhaps each effpot type should have own 
     procedure :: finalize
  end type supercell_maker_t
contains

  subroutine initialize(self, sc_matrix)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: sc_matrix(3,3)
    self%sc_matrix(:,:)=sc_matrix(:,:)
  end subroutine initialize


  subroutine finalize(self)
    class(supercell_maker_t), intent(inout) :: self
  end subroutine finalize

end module m_supercell_maker

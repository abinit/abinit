!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_supercell
!! NAME
!! m_supercell
!!
!! FUNCTION
!! This module define the supercell_t type, which contains all the information of the supercell.
!!
!! Datatypes:
!!  supercell_t
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

module m_supercell
  use defs_basis
  use m_abicore
  use m_errors
  use m_supercell
  implicit none

  type ,public :: supercell_t
     integer :: sc_matrix(3,3)
     integer :: ncell
     logical :: has_lattice=.False.
     logical :: has_spin=.False.
     logical :: has_lwf=.False.
     logical :: has_electron=.False.
     ! lattice related
     integer :: natom
     integer :: ntypat
     real(dp):: rprimd(3,3)
     integer, allocatable :: typat(:)
     real(dp), allocatable :: xred(:, :)
     ! spin related
     integer :: nspin
     real(dp), allocatable :: ms(:)
     ! lwf related
     integer :: nlwf

   contains
     procedure:: initialize
     procedure :: finalize
  end type supercell_t
contains

  subroutine initialize(self)
    class(supercell_t), intent(inout) :: self
  end subroutine initialize


  subroutine finalize(self)
    class(supercell_t), intent(inout) :: self
  end subroutine finalize

end module m_supercell

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_atoms
!! NAME
!! m_unitcell
!!
!! FUNCTION
!! This module define the m_unitcell type, which defines an atomic structure, with spin, lwf, electron
!!
!! Datatypes:
!!  unitcell_t
!!
!! Subroutines:
!! 
!! !! COPYRIGHT !! Copyright (C) 2001-2019 ABINIT group (hexu) !! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE


#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"

module m_unitcell
  use defs_basis
  use m_abicore
  use m_errors
  use m_supercell, only: supercell
  use m_supercell_maker , only: supercell_maker

  implicit none

  type ,public :: unitcell_t
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
  end type unitcell_t
contains

  subroutine initialize(self)
    class(unitcell_t), intent(inout) :: self
  end subroutine initialize


  subroutine finalize(self)
    class(unitcell_t), intent(inout) :: self
  end subroutine finalize

  subroutine fill_supercell(self, supercell_maker, supercell)
    class(unitcell_t), intent(inout) :: self
    type(supercell_maker_t), intent(inout) :: supercell_maker
  end subroutine fill_supercell

end module m_unitcell

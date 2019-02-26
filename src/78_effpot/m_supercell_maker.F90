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
  use m_mathfuncs , only: mat33det
  implicit none

  type ,public :: supercell_maker_t
     integer :: sc_matrix(3,3)
     integer :: ncell
   contains
     generic:: initialize => initialize_2d, initialize_1d      ! perhaps each effpot type should have own 
     procedure:: initialize_1d
     procedure:: initialize_2d
     procedure :: finalize
  end type supercell_maker_t
contains

  subroutine initialize_2d(self, sc_matrix)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: sc_matrix(3, 3)
    self%sc_matrix(:,:)=sc_matrix
    self%ncell=mat33det(self%sc_matrix)
  end subroutine initialize_2d

  subroutine initialize_1d(self, ncell)
    class(supercell_maker_t), intent(inout) :: self
    integer, intent(in) :: ncell(3)
    integer :: i
    do i =1, 3
       self%sc_matrix(i,i) = ncell(i)
    end do
    self%ncell=mat33det(self%sc_matrix)
  end subroutine initialize_1d

  subroutine finalize(self)
    class(supercell_maker_t), intent(inout) :: self
    self%sc_matrix=0
    self%ncell=0
  end subroutine finalize

end module m_supercell_maker

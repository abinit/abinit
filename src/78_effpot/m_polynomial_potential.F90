!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_polynomial_potential
!! NAME
!! m_polynomial_potential
!!
!! FUNCTION
!! This module contains the polynomial type potential
!!
!!
!! Datatypes:
!!
!! polynomial_potential_t
!!
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
!!
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (hexu)
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

module m_polynomail_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_abstract_potential, only : abstract_potential_t
  use m_supercell_maker, only: supercell_maker_t

  implicit none
!!***

  integer, parameter :: displacement = 0, strain=1, spin=2, lwf=3, electron=4

  private
  type ,public, extends(abstract_potential_t) :: polynomial_potential_t
     integer, allocatable :: nature(:)
     integer :: order=-1
     ! abstract_matrix :: coeff ! TODO : add this after abstract matrix is made.
   contains
     procedure :: initialize
     procedure :: finalize
  end type polynomial_potential_t

contains

  subroutine initialize(self, nature, order)
    class(polynomial_potential_t), intent(inout) :: self  ! the effpot may save the states.
    integer, intent(in) :: order
    integer, intent(in) :: nature(order)
    ABI_ALLOCATE(self%nature, (order))
    self%nature(:)=nature(:)
    self%order=order
  end subroutine initialize

  subroutine finalize(self)
    class(polynomial_potential_t), intent(inout) :: self  ! the effpot may save the states.
    if (allocated(self%nature)) then
       ABI_DEALLOCATE(self%nature)
    end if
    self%order=-1
  end subroutine finalize


end module m_polynomail_potential

!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_primitive_bilinear_potential
!! NAME
!! m_primitive_bilinear_potential
!!
!! FUNCTION
!! This module contains the base type for all primitive_bilinear_potential.
!! It has the format of Rij:val, i.e. For each supercell vector R, it can be represented
!! by a matrix with index i, j.
!!
!!
!! Datatypes:
!!
!! * primitive_bilinear_potential_t
!! Subroutines:
!! TODO: add this when F2003 doc style is determined.
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

module m_primitive_bilinear_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi
  use m_primitive_potential, only: primitive_potential_t
  use m_supercell_maker, only: supercell_maker_t
  use m_polynomial_potential, only: polynomial_potential_t
  !use m_spmat_base

  implicit none
!!***
  private

  type, public, extends(primitive_potential_t) ::  primitive_bilinear_potential_t
     integer :: nR
     integer :: order=2
     integer ::nature(2)
     integer, allocatable :: Rlist(:, :)      ! 3(xyz), indR
     real(dp), allocatable :: mat(:, :, :) ! index: i, j, indR
   contains
     procedure :: initialize
     procedure :: finalize
  end type primitive_bilinear_potential_t
contains

  subroutine initialize(self, nR, ni, nj, nature)
    class(primitive_bilinear_potential_t) :: self
    integer, intent(in) :: nR, ni, nj, nature(2)
    self%nR=nR
    self%nature(:)=nature(:)
    ABI_ALLOCATE(self%Rlist, (3, self%nR))
    ABI_ALLOCATE(self%mat, (ni, nj, nR))
  end subroutine initialize

  subroutine finalize(self)
    class(primitive_bilinear_potential_t) :: self
    self%nR=0
    if ( allocated(self%Rlist)) then
       ABI_DEALLOCATE(self%Rlist)
    end if
    if (allocated(self%mat)) then
       ABI_DEALLOCATE(self%mat)
    end if
  end subroutine finalize

  subroutine fill_supercell(self, sc_maker, scpot)
    class(primitive_bilinear_potential_t) :: self
    type(supercell_maker_t),intent(inout) :: sc_maker
    type(polynomial_potential_t) :: scpot
  end subroutine fill_supercell

end module m_primitive_bilinear_potential

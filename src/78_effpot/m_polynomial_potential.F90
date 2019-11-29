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

module m_polynomial_potential
  use defs_basis
  use m_abicore
  use m_errors
  use m_xmpi

  use m_mathfuncs, only: find_int
  use m_abstract_potential, only : abstract_potential_t
  use m_supercell_maker, only: supercell_maker_t
  use m_spmat_ndcoo, only: NDCOO_mat_t

  implicit none
!!***

  !integer, parameter :: displacement = 0, strain=1, spin=2, lwf=3, electron=4
  enum, bind(c)
     enumerator :: null_nature=0, displacement=1, strain=2, spin=3, lwf=4, electron=5
  endenum

  private
  type ,public, extends(abstract_potential_t) :: polynomial_potential_t
     integer(kind(null_nature)), allocatable :: nature(:)
     integer :: order=0
     ! abstract_matrix :: coeff ! TODO : add this after abstract matrix is made.
     ! type(NDCOO_mat_t) :: coeff
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: calculate
  end type polynomial_potential_t

contains

  subroutine initialize(self, nature, order, mshape)
    class(polynomial_potential_t), intent(inout) :: self  ! the effpot may save the states.
    integer, intent(in) :: order
    integer, intent(in) :: nature(order), mshape(order)
    ABI_ALLOCATE(self%nature, (order))
    self%nature(:)=nature(:)
    self%order=order
    if (find_int(nature, displacement)/=0) self%has_displacement=.True.
    if (find_int(nature, strain)/=0) self%has_strain=.True.
    if (find_int(nature, spin)/=0) self%has_spin=.True.
    if (find_int(nature, lwf)/=0) self%has_lwf=.True.
    !call self%coeff%initialize(mshape=mshape)
    ABI_UNUSED_A(mshape)
    self%label="PolynomialPotential"
  end subroutine initialize

  subroutine finalize(self)
    class(polynomial_potential_t), intent(inout) :: self  ! the effpot may save the states.
    if (allocated(self%nature)) then
       ABI_DEALLOCATE(self%nature)
    end if
    !call self%coeff%finalize()
    self%order=0
    call self%abstract_potential_t%finalize()
  end subroutine finalize

  subroutine calculate(self, displacement, strain, spin, lwf, force, stress, bfield, lwf_force, energy)
    class(polynomial_potential_t), intent(inout) :: self  ! the effpot may save the states.

    real(dp), optional, intent(inout) :: displacement(:,:), strain(:,:), spin(:,:), lwf(:)
    real(dp), optional, intent(inout) :: force(:,:), stress(:,:), bfield(:,:), lwf_force(:), energy
    ABI_UNUSED_A(self)
    ABI_UNUSED_A(displacement)
    ABI_UNUSED_A(strain)
    ABI_UNUSED_A(spin)
    ABI_UNUSED_A(lwf)
    ABI_UNUSED_A(force)
    ABI_UNUSED_A(stress)
    ABI_UNUSED_A(bfield)
    ABI_UNUSED_A(lwf_force)
    ABI_UNUSED_A(energy)
    MSG_ERROR("calculate not implemented for this effpot.")
  end subroutine calculate


end module m_polynomial_potential

!!****m* ABINIT/m_lattice_lwf_map
!! NAME
!! m_lattice_lwf_map
!!
!! FUNCTION
!! This module contains the functions to map between Lattice and lwf, including the amplitude, and the force.
!!
!!
!! Subroutines:
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

module m_lattice_lwf_map
  use, intrinsic :: iso_c_binding
  !use m_dynamic_array, only: int_array_type, real_array_type, int2d_array_type
  use defs_basis
  use m_abicore
  use m_errors
  use m_nctk
  !use m_spmat_spvec, only: sp_real_vec
  use m_spmat_coo, only: COO_mat_t
!#if defined HAVE_NETCDF
  use netcdf
!#endif
  implicit none
  private
  !!***

  type, public :: lwf_latt_coeff_t
     type(COO_mat_t), public :: coeffs
   contains
     procedure :: initialize
     procedure :: finalize
     procedure :: lattice_to_lwf_projection
     procedure :: lwf_force_to_lattice
     procedure :: lwf_amp_to_displacements
  end type lwf_latt_coeff_t

contains

  subroutine initialize(self, nlwf, natom)
    class(lwf_latt_coeff_t), intent(inout) :: self
    integer, intent(in) :: nlwf, natom
    call self%coeffs%initialize([natom, nlwf])
  end subroutine initialize

  subroutine finalize(self)
    class(lwf_latt_coeff_t), intent(inout) :: self
    call self%coeffs%finalize()
  end subroutine finalize


!===============================================================
! Project lattice distortion amplitudes to lwf amplitudes
!> @
!===============================================================
  subroutine lattice_to_lwf_projection(self, latt_amp, lwf_amp)
    class(lwf_latt_coeff_t), intent(in) :: self
    real(dp), intent(in) :: latt_amp(:,:)
    real(dp), intent(inout) :: lwf_amp(:)
    ! integer :: ilwf, nlwf
    ! nlwf=size(coeffs)
    ! do ilwf=1, nlwf
    !    lwf_amp(ilwf) = coeffs(ilwf)%dot(latt_amp)
    ! end do
    call self%coeffs%mv(latt_amp, lwf_amp)
  end subroutine lattice_to_lwf_projection

  !===============================================================
  ! Map force on lwf to force on lattice
  !> @
  !===============================================================
  subroutine lwf_force_to_lattice(self, lwf_force, latt_force)
    class(lwf_latt_coeff_t), intent(in) :: self
    real(dp), intent(in) :: lwf_force(:)
    real(dp), intent(out) :: latt_force(:,:)
    !nlwf=size(coeffs)
    !do ilwf=1, nlwf
    !   call coeffs(ilwf)%plus_Ax(lwf_force(ilwf), latt_force)
    !end do
    call self%coeffs%mv_left(lwf_force, latt_force)
  end subroutine lwf_force_to_lattice

  !===============================================================
  ! Map lwf amplitudes to lattice amplitudes
  !> @
  !===============================================================
  subroutine lwf_amp_to_displacements(self, lwf, displacement)
    class(lwf_latt_coeff_t), intent(in) :: self
    real(dp), intent(in) :: lwf(:)
    real(dp), intent(inout) :: displacement(:,:)
    ! integer :: ilwf, nlwf
    ! nlwf=size(coeffs)
    ! do ilwf=1, nlwf
    !    call coeffs(ilwf)%plus_Ax(lwf(ilwf), displacement)
    ! end do
    call self%coeffs%mv_left(lwf, displacement)
  end subroutine lwf_amp_to_displacements

end module m_lattice_lwf_map


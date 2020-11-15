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
  use iso_c_binding
  !use m_dynamic_array, only: int_array_type, real_array_type, int2d_array_type
  use defs_basis
  use m_abicore
  use m_errors
  use m_nctk
  use m_spmat_spvec, only: sp_real_vec
!#if defined HAVE_NETCDF
  use netcdf
!#endif
  implicit none
  private
  !!***

  type, public :: lwf_latt_coeff_t
     !! for each lwf, map to the 
     type(sp_real_vec), allocatable :: coeffs(:)
  end type lwf_latt_coeff_t
contains


!===============================================================
! Project lattice distortion amplitudes to lwf amplitudes
!> @
!===============================================================
  subroutine lattice_to_lwf_projection(coeffs, latt_amp, lwf_amp)
    type(sp_real_vec) , intent(in) :: coeffs(:)
    real(dp), intent(in) :: latt_amp(:,:)
    real(dp), intent(inout) :: lwf_amp(:)
    integer :: ilwf, nlwf
    nlwf=size(coeffs)
    do ilwf=1, nlwf
       lwf_amp(ilwf) = coeffs(ilwf)%dot(latt_amp)
    end do
  end subroutine lattice_to_lwf_projection

  !===============================================================
  ! Map force on lwf to force on lattice
  !> @
  !===============================================================
  subroutine lwf_force_to_lattice(lwf_force, coeffs, latt_force)
    real(dp), intent(in) :: lwf_force(:)
    type(sp_real_vec) , intent(in) :: coeffs(:)
    real(dp), intent(out) :: latt_force(:,:)
    integer :: ilwf, nlwf
    nlwf=size(coeffs)
    do ilwf=1, nlwf
       call coeffs(ilwf)%plus_Ax(lwf_force(ilwf), latt_force)
    end do
  end subroutine lwf_force_to_lattice


end module m_lattice_lwf_map


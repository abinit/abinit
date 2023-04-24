! ***************************************************************************************************
!  Copyright (C) 2020-2023 Green-X library
!  This file is distributed under the terms of the APACHE2 License.
!
! ***************************************************************************************************
!> \brief This module contains auxiliary procedures and data structures for the main minimax routines
! ***************************************************************************************************
module minimax_utils
#if defined HAVE_CONFIG_H
#include "config.h"
#endif
#include "abi_common.h"
#include "gx_common.h"
  use defs_basis,     only: dp
  use m_errors
  !use kinds, only: dp
  implicit none

  private

  type, public :: er_aw_aux
     ! Sorted array of the energy ranges
     real(kind=dp), dimension(:), allocatable :: energy_range
     ! Matrices with coefficients and weights per energy region
     real(kind=dp), dimension(:, :), allocatable :: aw_erange_matrix
   contains
     procedure :: get_coeff_weight => coeffs_and_weights
  end type er_aw_aux

  !> Transformation types
  integer, parameter, public :: cosine_tw = 1
  integer, parameter, public :: cosine_wt = 2
  integer, parameter, public :: sine_tw = 3

contains

  !> \brief Find first element in unsorted array that is strictly greater than a given value
  !>        This algorithm is O(n), difficult to do better with unsorted arrays
  !! @param[in] lenght - lenght of sorted array
  !! @param[in] einter - sorted array of the energy intervals
  !! @param[in] eval - the energy value
  function find_erange(length, einter, eval) result(idx)
    integer, intent(in)                          :: length
    real(kind=dp), dimension(length), intent(in) :: einter
    real(kind=dp), intent(in)                    :: eval
    integer                                      :: idx

    ! Auxiliary variables
    integer                                      :: jdx
    real(kind=dp)                                :: tmp_min_max

    ! Begin work
    tmp_min_max = huge(0.0_dp)
    idx = length + 1

    do jdx = 1, length
       if (eval < einter(jdx) .and. einter(jdx) < tmp_min_max) then
          idx = jdx
          tmp_min_max = einter(jdx)
       end if
    end do

  end function find_erange

  !> \brief Selects the energy region and scales weights and coefficients
  !! @param[in] grid_size - the grid size
  !! @param[in] bup - length of the energy region array
  !! @param[in] e_range - the selected energy range
  !! @param[inout] e_ratio - an heuristic correction factor
  !! @param[inout] ac_we - vector containing coefficients and weights
  subroutine coeffs_and_weights(this, grid_size, bup, e_range, ac_we, e_ratio)
    class(er_aw_aux), intent(in)               :: this
    integer, intent(in)                        :: grid_size
    integer, intent(in)                        :: bup
    real(kind=dp), intent(in)                  :: e_range
    real(kind=dp), dimension(:), intent(inout) :: ac_we
    real(kind=dp), intent(inout)               :: e_ratio

    ! Internal variables
    integer                                    :: ien

    ! Select energy region
    ien = find_erange(bup, this%energy_range, e_range)

    ! Scale grids for large sizes when erange falls in the first energy range
    if (ien == 1 .and. grid_size > 20) then
       e_ratio = this%energy_range(1) / e_range
       if (e_ratio > 1.5_dp) then
          e_ratio = e_ratio / 1.5_dp
       endif
    end if

    ac_we(:) = this%aw_erange_matrix(:, ien)

  end subroutine coeffs_and_weights

end module minimax_utils

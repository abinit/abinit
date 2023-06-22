!!****m* ABINIT/m_xg_kokkos
!! NAME
!!  m_xg_kokkos
!!
!! FUNCTION
!! This provides iso_c_binding wrappers to kokkos implementation of some of the
!! computational routines located in m_xg
!!
!! COPYRIGHT
!!  Copyright (C) 2016-2022 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xg_kokkos

  use defs_basis
  use, intrinsic :: iso_c_binding

#if defined HAVE_YAKL
 use gator_mod
#endif

  implicit none

  interface

    ! ========================================================================
    ! ========================================================================
    subroutine computeBatchedDotProduct_scalar(x_ptr, y_ptr, res_ptr, nx, ny, ldim) &
      & bind(c, name='computeBatchedDotProduct_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: x_ptr
      type(c_ptr)            , value             :: y_ptr
      type(c_ptr)            , value             :: res_ptr
      integer(kind=c_int32_t), value, intent(in) :: nx
      integer(kind=c_int32_t), value, intent(in) :: ny
      integer(kind=c_int32_t), value, intent(in) :: ldim
    end subroutine computeBatchedDotProduct_scalar

    subroutine computeBatchedDotProduct_cplx(x_ptr, y_ptr, res_ptr, nx, ny, ldim) &
      & bind(c, name='computeBatchedDotProduct_cplx_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: x_ptr
      type(c_ptr)            , value             :: y_ptr
      type(c_ptr)            , value             :: res_ptr
      integer(kind=c_int32_t), value, intent(in) :: nx
      integer(kind=c_int32_t), value, intent(in) :: ny
      integer(kind=c_int32_t), value, intent(in) :: ldim
    end subroutine computeBatchedDotProduct_cplx

    subroutine computeBatchedDotProduct_cplx_scalar(x_ptr, y_ptr, res_ptr, nx, ny, ldim) &
      & bind(c, name='computeBatchedDotProduct_cplx_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: x_ptr
      type(c_ptr)            , value             :: y_ptr
      type(c_ptr)            , value             :: res_ptr
      integer(kind=c_int32_t), value, intent(in) :: nx
      integer(kind=c_int32_t), value, intent(in) :: ny
      integer(kind=c_int32_t), value, intent(in) :: ldim
    end subroutine computeBatchedDotProduct_cplx_scalar

    ! ========================================================================
    ! ========================================================================
    subroutine computeMax_scalar(x_ptr, size, res) &
      & bind(c, name='computeMax_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      real(c_double),                 intent(inout) :: res
    end subroutine computeMax_scalar

    subroutine computeMax_complex(x_ptr, size, res) &
      & bind(c, name='computeMax_complex_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      real(c_double),                 intent(inout) :: res
    end subroutine computeMax_complex

    ! ========================================================================
    ! ========================================================================
    subroutine computeMin_scalar(x_ptr, size, res) &
      & bind(c, name='computeMin_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      real(c_double),                 intent(inout) :: res
    end subroutine computeMin_scalar

    subroutine computeMin_complex(x_ptr, size, res) &
      & bind(c, name='computeMin_complex_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      real(c_double),                 intent(inout) :: res
    end subroutine computeMin_complex

    ! ========================================================================
    ! ========================================================================
    subroutine computeMaxloc_scalar(x_ptr, size, res) &
      & bind(c, name='computeMaxloc_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      integer(c_int32_t),             intent(inout) :: res
    end subroutine computeMaxloc_scalar

    subroutine computeMaxloc_complex(x_ptr, size, res) &
      & bind(c, name='computeMaxloc_complex_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      integer(c_int32_t),             intent(inout) :: res
    end subroutine computeMaxloc_complex

    subroutine computeMaxloc_scalar_2d(x_ptr, nx, ny, res_ptr) &
      & bind(c, name='computeMaxloc_scalar_2d_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr    ! double
      integer(kind=c_int32_t), value, intent(in)    :: nx
      integer(kind=c_int32_t), value, intent(in)    :: ny
      type(c_ptr)            , value                :: res_ptr  ! int32_t
    end subroutine computeMaxloc_scalar_2d

    subroutine computeMaxloc_complex_2d(x_ptr, nx, ny, res_ptr) &
      & bind(c, name='computeMaxloc_complex_2d_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr    ! cplx_t
      integer(kind=c_int32_t), value, intent(in)    :: nx
      integer(kind=c_int32_t), value, intent(in)    :: ny
      type(c_ptr)            , value                :: res_ptr  ! int32_t
    end subroutine computeMaxloc_complex_2d

    ! ========================================================================
    ! ========================================================================
    subroutine computeMinloc_scalar(x_ptr, size, res) &
      & bind(c, name='computeMinloc_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      integer(c_int32_t),             intent(inout) :: res
    end subroutine computeMinloc_scalar

    subroutine computeMinloc_complex(x_ptr, size, res) &
      & bind(c, name='computeMinloc_complex_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      integer(c_int32_t),             intent(inout) :: res
    end subroutine computeMinloc_complex

    subroutine computeMinloc_scalar_2d(x_ptr, nx, ny, res_ptr) &
      & bind(c, name='computeMinloc_scalar_2d_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr    ! double
      integer(kind=c_int32_t), value, intent(in)    :: nx
      integer(kind=c_int32_t), value, intent(in)    :: ny
      type(c_ptr)            , value                :: res_ptr  ! int32_t
    end subroutine computeMinloc_scalar_2d

    subroutine computeMinloc_complex_2d(x_ptr, nx, ny, res_ptr) &
      & bind(c, name='computeMinloc_complex_2d_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr    ! cplx_t
      integer(kind=c_int32_t), value, intent(in)    :: nx
      integer(kind=c_int32_t), value, intent(in)    :: ny
      type(c_ptr)            , value                :: res_ptr  ! int32_t
    end subroutine computeMinloc_complex_2d

    ! ========================================================================
    ! ========================================================================
    subroutine computeColwiseDivision_scalar(x_ptr, y_ptr, size, res_ptr) &
      & bind(c, name='computeColwiseDivision_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      type(c_ptr)            , value                :: y_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      type(c_ptr)            , value                :: res_ptr
    end subroutine computeColwiseDivision_scalar

    subroutine computeColwiseDivision_complex(x_ptr, y_ptr, size, res_ptr) &
      & bind(c, name='computeColwiseDivision_complex_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value                :: x_ptr
      type(c_ptr)            , value                :: y_ptr
      integer(kind=c_int32_t), value, intent(in)    :: size
      type(c_ptr)            , value                :: res_ptr
    end subroutine computeColwiseDivision_complex

    ! ========================================================================
    ! ========================================================================
    subroutine compute_colwiseCymax_scalar(A_ptr, da_ptr, B_ptr, W_ptr, rows, cols, ldim) &
      & bind(c, name='compute_colwiseCymax_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: A_ptr
      type(c_ptr)            , value             :: da_ptr
      type(c_ptr)            , value             :: B_ptr
      type(c_ptr)            , value             :: W_ptr
      integer(kind=c_int32_t), value, intent(in) :: rows, cols, ldim
    end subroutine compute_colwiseCymax_scalar

    subroutine compute_colwiseCymax_cplx(A_ptr, da_ptr, B_ptr, W_ptr, rows, cols, ldim) &
      & bind(c, name='compute_colwiseCymax_cplx_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: A_ptr
      type(c_ptr)            , value             :: da_ptr
      type(c_ptr)            , value             :: B_ptr
      type(c_ptr)            , value             :: W_ptr
      integer(kind=c_int32_t), value, intent(in) :: rows, cols, ldim
    end subroutine compute_colwiseCymax_cplx

    ! ========================================================================
    ! ========================================================================
    subroutine compute_colwiseMul_scalar_scalar(data_ptr, vec_ptr, shift, rows, cols, ldim, vec_size) &
      & bind(c, name='compute_colwiseMul_scalar_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: data_ptr
      type(c_ptr)            , value             :: vec_ptr
      integer(kind=c_int32_t), value, intent(in) :: shift, rows, cols, ldim, vec_size
    end subroutine compute_colwiseMul_scalar_scalar

    subroutine compute_colwiseMul_cplx_scalar(data_ptr, vec_ptr, shift, rows, cols, ldim, vec_size) &
      & bind(c, name='compute_colwiseMul_cplx_scalar_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: data_ptr
      type(c_ptr)            , value             :: vec_ptr
      integer(kind=c_int32_t), value, intent(in) :: shift, rows, cols, ldim, vec_size
    end subroutine compute_colwiseMul_cplx_scalar

    subroutine compute_colwiseMul_cplx_cplx(data_ptr, vec_ptr, shift, rows, cols, ldim, vec_size) &
      & bind(c, name='compute_colwiseMul_cplx_cplx_kokkos_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr)            , value             :: data_ptr
      type(c_ptr)            , value             :: vec_ptr
      integer(kind=c_int32_t), value, intent(in) :: shift, rows, cols, ldim, vec_size
    end subroutine compute_colwiseMul_cplx_cplx

  end interface

contains
  !!***

end module m_xg_kokkos
!!***

!!****m* ABINIT/m_gpu_toolbox
!! NAME
!!  m_gpu_toolbox
!!
!! FUNCTION
!!  Fake module to dupe the build system and allow it to include cuda files
!!   in the chain of dependencies.
!!
!! COPYRIGHT
!!  Copyright (C) 2000-2022 ABINIT group (MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gpu_toolbox

  use m_initcuda
  use m_gpu_detect

! MG: I had to comment this import to avoid the following error on buda2_gnu_8.5_cuda
!
!    type(c_ptr),intent(inout) :: blockvectorbx_gpu, blockvectorx_gpu,sqgram_gpu
!            1
!    Error: Type name 'c_ptr' at (1) is ambiguous
!    abi_gpu_linalg.f90:374:47:
!
! I believe this is due to a misconfiguration issue in the Fortran compilers used by the bot.

#ifdef HAVE_FC_ISO_C_BINDING
 use, intrinsic :: iso_c_binding, only : C_INT32_T,C_SIZE_T
#endif

  implicit none

  !Interfaces for C bindings --- To be completed
#ifdef HAVE_FC_ISO_C_BINDING
#if defined HAVE_GPU_CUDA

  interface

    !  integer(C_INT) function cuda_func() bind(C)
    !    use iso_c_binding, only : C_INT,C_PTR
    !    type(C_PTR) :: ptr
    !  end function cuda_func

    subroutine gpu_device_synchronize() bind(c, name='gpu_device_synchronize_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
    end subroutine gpu_device_synchronize

    subroutine gpu_get_device(deviceId) bind(c, name='gpu_get_device_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      integer(kind=C_INT32_T), intent(inout) :: deviceId
    end subroutine gpu_get_device

    subroutine gpu_data_prefetch_async_f(dev_ptr, count, deviceId) bind(c, name='gpu_data_prefetch_async_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),             value :: dev_ptr
      integer(kind=C_SIZE_T),  value :: count
      integer(kind=C_INT32_T), value :: deviceId
    end subroutine gpu_data_prefetch_async_f

  end interface

  public
  ! -1 is a special value for device id, it represent actually the host (CPU)
  ! see https://docs.nvidia.com/cuda/cuda-runtime-api/group__CUDART__TYPES.html
  integer(C_INT32_T), public, parameter :: CPU_DEVICE_ID = -1

#endif
#endif

contains
  !!***

#ifdef HAVE_FC_ISO_C_BINDING
#if defined HAVE_GPU_CUDA

  ! prefetch data (memory managed pointer) to device
  ! device can be a GPU (deviceId >=0)
  ! device can be   CPU deviceId = -1)
  subroutine gpu_data_prefetch_async(dev_ptr, count, deviceId)
    use, intrinsic :: iso_c_binding
    implicit none
    type(c_ptr),             value          :: dev_ptr
    integer(kind=C_SIZE_T),  value          :: count
    integer(kind=C_INT32_T), value,optional :: deviceId

    integer(kind=C_INT32_T)                 :: currentDevId

    ! if a device id is provided, use it
    ! if not, just probe driver to get current device id
    if (present(deviceId)) then

      if ( deviceId >= CPU_DEVICE_ID ) then
        call gpu_data_prefetch_async_f(dev_ptr, count, deviceId)
      end if

    else

      call gpu_get_device(currentDevId)
      call gpu_data_prefetch_async_f(dev_ptr, count, currentDevId)

    end if

  end subroutine gpu_data_prefetch_async

#endif
#endif

end module m_gpu_toolbox
!!***

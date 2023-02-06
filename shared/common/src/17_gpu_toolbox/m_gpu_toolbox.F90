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

  ! mirroring cuda enum cudaMemoryAdvise usually defined in
  ! /usr/local/cuda/targets/x86_64-linux/include/driver_types.h
  !
  ! to be used as 3rd arg of gpu_memory_advise_f
  !
  ! I didn't find a clean way of using an existing enum defined in C
  ! without redefining it in fortran
  ! I didn't found a way to do it through iso_c_binding, strange...
  ! It means this enum will have to be updated if ever the C defined enum
  ! changes.
  enum, bind(c)
    ! Data will mostly be read and only occassionally be written to
    enumerator :: CUDA_MEM_ADVISE_SET_READ_MOSTLY          = 1

    ! Undo the effect of ::cudaMemAdviseSetReadMostly
    enumerator :: CUDA_MEM_ADVISE_UNSET_READ_MOSTLY        = 2

    ! Set the preferred location for the data as the specified device
    enumerator :: CUDA_MEM_ADVISE_SET_PREFERRED_LOCATION   = 3

    ! Clear the preferred location for the data
    enumerator :: CUDA_MEM_ADVISE_UNSET_PREFERRED_LOCATION = 4

    ! Data will be accessed by the specified device, so prevent page faults as much as possible
    enumerator :: CUDA_MEM_ADVISE_SET_ACCESSED_BY          = 5

    ! Let the Unified Memory subsystem decide on the page faulting policy for the specified device
    enumerator :: CUDA_MEM_ADVISE_UNSET_ACCESSED_BY        = 6
  end enum

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

    subroutine gpu_memory_advise_f(dev_ptr, count, advice, deviceId) bind(c, name='gpu_memory_advise_cpp')
      use, intrinsic :: iso_c_binding
      implicit none
      type(c_ptr),                       value :: dev_ptr
      integer(kind=C_SIZE_T),            value :: count
      integer(kind=C_INT),               value :: advice
      integer(kind=C_INT32_T),           value :: deviceId
    end subroutine gpu_memory_advise_f

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

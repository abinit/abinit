!!****m* ABINIT/m_gputk
!! NAME
!!  m_gputk
!!
!! FUNCTION
!!  Low-level procedures for GPUs.
!!
!! COPYRIGHT
!!  Copyright (C) 2012-2025 ABINIT group (LNguyen,FDahm,MT)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gputk

 use, intrinsic :: iso_c_binding
 USE_MPI
 use defs_basis
 use m_errors
 use m_abicore
 !use m_xomp
!#if defined HAVE_GPU
! use m_gpu_toolbox
!#endif
!
#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 implicit none

 private
!!***

 !----------------------------------------------------------------------

 public :: gpu_set_to_zero
 public :: gpu_set_to_zero_sp
 public :: gpu_set_to_zero_complex
 public :: gpu_set_to_zero_complex_sp
 public :: gpu_copy
 public :: gpu_copy_sp
 public :: gpu_copy_complex
 public :: gpu_copy_complex_sp

 !----------------------------------------------------------------------

#ifdef HAVE_GPU
  interface
    subroutine check_gpu_mem(str) bind(c, name="check_gpu_mem_")
      use, intrinsic :: iso_c_binding
      character (KIND=c_char), intent(in)  :: str(*)
    end subroutine check_gpu_mem

    subroutine alloc_on_gpu(gpu_ptr,size_in_bytes) bind(c, name="alloc_on_gpu_")
      use, intrinsic :: iso_c_binding
      type(c_ptr),                    intent(inout)  :: gpu_ptr
      integer(kind=c_size_t),         intent(in)     :: size_in_bytes
    end subroutine alloc_on_gpu

    subroutine dealloc_on_gpu(gpu_ptr) bind(c, name="dealloc_on_gpu_")
      use, intrinsic :: iso_c_binding
      type(c_ptr),                    intent(inout)  :: gpu_ptr
    end subroutine dealloc_on_gpu

    subroutine copy_gpu_to_gpu(dest_gpu_ptr, src_gpu_ptr, size_in_bytes) bind(c, name="copy_gpu_to_gpu_cpp_")
      use, intrinsic :: iso_c_binding
      type(c_ptr)                                   :: dest_gpu_ptr
      type(c_ptr)                                   :: src_gpu_ptr
      integer(kind=c_size_t),        intent(in)    :: size_in_bytes
    end subroutine copy_gpu_to_gpu

    subroutine gpu_memset(gpu_ptr, val, size_in_bytes) bind(c, name="gpu_memset_cpp_")
      use, intrinsic :: iso_c_binding
      type(c_ptr),                    intent(in) :: gpu_ptr
      integer(kind=c_int32_t),        intent(in)    :: val
      integer(kind=c_size_t),         intent(in)    :: size_in_bytes
    end subroutine gpu_memset

    ! logical(kind=c_bool) function gpu_allocated(gpu_ptr) bind(c, name="gpu_allocated_")
    !   use, intrinsic :: iso_c_binding
    !   type(c_ptr),                    intent(in) :: gpu_ptr
    ! end function gpu_allocated

    subroutine gpu_allocated_impl(gpu_ptr, is_allocated) bind(c, name="gpu_allocated_impl_")
      use, intrinsic :: iso_c_binding
      type(c_ptr),                    intent(in)  :: gpu_ptr
      logical(kind=c_bool),           intent(out) :: is_allocated
    end subroutine gpu_allocated_impl

    subroutine gpu_managed_ptr_status(gpu_ptr, str) bind(c, name="gpu_managed_ptr_status_")
      use, intrinsic :: iso_c_binding
      type(c_ptr),                    intent(in)  :: gpu_ptr
      character (KIND=c_char),        intent(in)  :: str(*)
    end subroutine gpu_managed_ptr_status

  end interface

#else
 !dummy routines replace gpu helper routines
 public :: gpu_device_synchronize
 public :: check_gpu_mem
 public :: alloc_on_gpu
 public :: copy_from_gpu
 public :: copy_on_gpu
 public :: dealloc_on_gpu
 public :: gpu_allocated_impl
 public :: gpu_managed_ptr_status
#endif

 public :: copy_gpu_to_gpu
 public :: gpu_memset
 public :: gpu_allocated


CONTAINS  !===========================================================
!!***

!!
!! this is just a wrapper arround gpu_allocated_cuda, because (strangely)
!! I can't manage to bind a function (not a subroutine) through iso_c_binding
!!
function gpu_allocated(gpu_ptr) result(is_allocated)

  !Arguments ------------------------------------
  type(c_ptr),                    intent(in) :: gpu_ptr
  logical(kind=c_bool)                       :: is_allocated

  call gpu_allocated_impl(gpu_ptr, is_allocated)

end function gpu_allocated

!----------------------------------------------------------------------

#ifndef HAVE_GPU

!!****f* m_gputk/gpu_device_synchronize
!! NAME
!!  gpu_device_synchronize
!!
!! FUNCTION
!!  Wait for any running operation, compute and memory transfer, to complete on GPU.
!!
!! INPUTS
!!  None
!!
!! OUTPUT
!!  None
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine gpu_device_synchronize()
  use, intrinsic :: iso_c_binding
  implicit none
end subroutine gpu_device_synchronize
!!***


!!****f* m_gputk/check_gpu_mem
!! NAME
!!  check_gpu_mem
!!
!! FUNCTION
!!  Print information about amount of free memory on GPU and total amount of memory on GPU (current device).
!!
!! INPUTS
!!  str is a string message (character array).
!!
!! OUTPUT
!!  None
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine check_gpu_mem(str)

  !Arguments ------------------------------------
  character (KIND=c_char), intent(in), target  :: str(*)
  !Local variables ------------------------------
  type(c_ptr)                                  :: dummy

  if(.false.) dummy=c_loc(str)

end subroutine check_gpu_mem
!!***

!!****f* m_gputk/alloc_on_gpu
!! NAME
!!  alloc_on_gpu
!!
!! FUNCTION
!!  Allocate size byte in gpu memory and returns in gpu_ptr this location
!!
!! INPUTS
!!  size= size in byte to allocate
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location that has been allocated
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine alloc_on_gpu(gpu_ptr,size)

!Arguments ------------------------------------
 type(c_ptr),                intent(inout) :: gpu_ptr
 integer(kind=c_size_t),     intent(in)    :: size ! size in bytes to allocate

 ABI_UNUSED(gpu_ptr)
 ABI_UNUSED(size)

end subroutine alloc_on_gpu
!!***

!!****f* m_gputk/copy_from_gpu
!! NAME
!!  copy_from_gpu
!!
!! FUNCTION
!!  copy size byte from gpu memory (pointed by gpu_ptr) to cpu memory (pointed by cpu_ptr)
!!
!! INPUTS
!!  size_in_bytes = size in bytes to allocate
!!  gpu_ptr = C_PTR on gpu memory location that has been allocated
!!
!! OUTPUT
!!  dtab = fortran tab which will contains data
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine copy_from_gpu(dtab,gpu_ptr,size_in_bytes)

!Arguments ------------------------------------
 real(dp),dimension(*)               :: dtab
 type(c_ptr)                         :: gpu_ptr
 integer(kind=c_size_t), intent(in)  :: size_in_bytes ! size in byte (to be transfered)

!Local variables ------------------------------
 type(c_ptr)                         :: cpu_ptr

 if(.false.) write(std_out,*) dtab(1)
 ABI_UNUSED(cpu_ptr)
 ABI_UNUSED(gpu_ptr)
 ABI_UNUSED(size_in_bytes)

end subroutine copy_from_gpu
!!***

!!****f* m_gputk/copy_on_gpu
!! NAME
!!  copy_on_gpu
!!
!! FUNCTION
!!  copy size byte from cpu (pointed by cpu_ptr) to gpu memory (pointed by gpu_ptr)
!!
!! INPUTS
!!  size_in_bytes = size in bytes to allocate
!!  dtab = fortran tab to copy
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine copy_on_gpu(dtab,gpu_ptr,size_in_bytes)

  !Arguments ------------------------------------
  real(dp),dimension(*)               :: dtab
  type(c_ptr)                         :: gpu_ptr
  integer(kind=c_size_t), intent(in)  :: size_in_bytes ! size in byte (to be transfered)

  !Local variables ------------------------------
  type(c_ptr)                         :: cpu_ptr

  if(.false.) write(std_out,*) dtab(1)
  ABI_UNUSED(cpu_ptr)
  ABI_UNUSED(gpu_ptr)
  ABI_UNUSED(size_in_bytes)

end subroutine copy_on_gpu
!!***

!!****f* m_gputk/copy_gpu_to_gpu
!! NAME
!!  copy_gpu_to_gpu
!!
!! FUNCTION
!!  copy size byte from gpu (src) to gpu (dest)
!!
!! INPUTS
!!  size_in_bytes = size in bytes to copy
!!  src_gpu_ptr = C_PTR on gpu memory
!!
!! OUTPUT
!!  dest_gpu_ptr = C_PTR on gpu memory
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine copy_gpu_to_gpu(cpu_ptr,gpu_ptr,size_in_bytes)

  !Arguments ------------------------------------
  type(c_ptr)                         :: cpu_ptr
  type(c_ptr)                         :: gpu_ptr
  integer(kind=c_size_t), intent(in)  :: size_in_bytes ! size in byte (to be transfered)

  ABI_UNUSED(cpu_ptr)
  ABI_UNUSED(gpu_ptr)
  ABI_UNUSED(size_in_bytes)

end subroutine copy_gpu_to_gpu
!!***

!!****f* m_gputk/dealloc_on_gpu
!! NAME
!!  dealloc_on_gpu
!!
!! FUNCTION
!!  free memory location pointed by gpu_ptr
!!
!! INPUTS
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location that has been allocated
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine dealloc_on_gpu(gpu_ptr)

  !Arguments ------------------------------------
  type(c_ptr) :: gpu_ptr

  ABI_UNUSED(gpu_ptr)

end subroutine dealloc_on_gpu
!!***

!!****f* m_gputk/gpu_memset
!! NAME
!!  gpu_memset
!!
!! FUNCTION
!!  Initializes or sets device memory to a value.
!!
!! INPUTS
!!  gpu_ptr= C_PTR on gpu memory location
!!  val= value used to initialized each bytes
!!  size= number of bytes to initialize
!!
!! OUTPUT
!!  gpu_ptr= C_PTR on gpu memory location
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine gpu_memset(gpu_ptr, val, array_size)

  !Arguments ------------------------------------
  type(c_ptr)                         :: gpu_ptr
  integer(kind=c_int32_t), intent(in) :: val
  integer(kind=c_size_t),  intent(in) :: array_size

  ABI_UNUSED(gpu_ptr)
  ABI_UNUSED(val)
  ABI_UNUSED(array_size)

end subroutine gpu_memset
!!***

!!****f* m_gputk/gpu_allocated_impl
!! NAME
!!  gpu_allocated_impl
!!
!! FUNCTION
!!  Check if pointer points to allocated gpu device memory.
!!
!! INPUTS
!!  gpu_ptr= C_PTR on gpu memory location
!!
!! OUTPUT
!!  is_allocate= logical(c_bool) : true (if allocated), false (if not allocated)
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine gpu_allocated_impl(gpu_ptr, is_allocated)

  !Arguments ------------------------------------
  type(c_ptr)                       :: gpu_ptr
  logical(kind=c_bool), intent(out) :: is_allocated

  ABI_UNUSED(gpu_ptr)

  is_allocated = .false.

end subroutine gpu_allocated_impl
!!***

!!****f* m_gputk/gpu_managed_ptr_status
!! NAME
!!  gpu_managed_ptr_status_impl
!!
!! FUNCTION
!!  Print information about a managed pointer (host or device address when accessible).
!!
!! INPUTS
!!  gpu_ptr= C_PTR on gpu memory location
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/dev_spec.cu
!!
!! SOURCE

subroutine gpu_managed_ptr_status(gpu_ptr, str)

  !Arguments ------------------------------------
  type(c_ptr)                                  :: gpu_ptr
  character (KIND=c_char), intent(in), target  :: str(*)
  !Local variables ------------------------------
  type(c_ptr)                                  :: dummy

  ABI_UNUSED(gpu_ptr)
  if(.false.) dummy=c_loc(str)

end subroutine gpu_managed_ptr_status
!!***
#endif

!------------------------------------------------------------------------------
!!****f* m_gputk/gpu_set_to_zero
!! NAME
!!  gpu_set_to_zero
!!
!! FUNCTION
!!  Set array content to zero
!!
!! INPUTS
!!  size = size of array
!!
!! OUTPUT
!!  array  = array to be set to zero
!!
!! SOURCE

subroutine gpu_set_to_zero(array, sizea)
 integer(c_size_t),intent(in)  :: sizea
 real(dp),target,intent(out) :: array(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(array)
 call gpu_memset(c_loc(array), 0, sizea*dp)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:array)
 do i=1,sizea
   array(i)=zero
 end do
#endif

#endif

end subroutine gpu_set_to_zero
!!***

!------------------------------------------------------------------------------
!!****f* m_gputk/gpu_set_to_zero_sp
!! NAME
!!  gpu_set_to_zero_sp
!!
!! FUNCTION
!!  Set array content to zero
!!
!! INPUTS
!!  size = size of array
!!
!! OUTPUT
!!  array  = array to be set to zero
!!
!! SOURCE

subroutine gpu_set_to_zero_sp(array, sizea)
 integer(c_size_t),intent(in)  :: sizea
 real(sp),target,intent(out) :: array(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(array)
 call gpu_memset(c_loc(array), 0, sizea*sp)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:array)
 do i=1,sizea
   array(i)=zero_sp
 end do
#endif

#endif

end subroutine gpu_set_to_zero_sp

!------------------------------------------------------------------------------
!!****f* m_gputk/gpu_set_to_zero_complex
!! NAME
!!  gpu_set_to_zero_complex
!!
!! FUNCTION
!!  Set array content to zero
!!
!! INPUTS
!!  size = size of array
!!
!! OUTPUT
!!  array = array to be set to zero
!!
!! SOURCE

subroutine gpu_set_to_zero_complex(array, sizea)
 integer(c_size_t),intent(in)  :: sizea
 complex(dp),target,intent(out) :: array(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(array)
 call gpu_memset(c_loc(array), 0, sizea*dp*2)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:array)
 do i=1,sizea
   array(i)=czero
 end do
#endif

#endif

end subroutine gpu_set_to_zero_complex
!!***

!!****f* m_gputk/gpu_set_to_zero_complex_sp
!! NAME
!!  gpu_set_to_zero_complex_sp
!!
!! FUNCTION
!!  Set array content to zero
!!
!! INPUTS
!!  size = size of array
!!
!! OUTPUT
!!  array = array to be set to zero
!!
!! SOURCE

subroutine gpu_set_to_zero_complex_sp(array, sizea)
 integer(c_size_t),intent(in)  :: sizea
 complex(sp),target,intent(out) :: array(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(array)
 call gpu_memset(c_loc(array), 0, sizea*sp*2)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:array)
 do i=1,sizea
   array(i)=czero_sp
 end do
#endif

#endif

end subroutine gpu_set_to_zero_complex_sp
!!***

!------------------------------------------------------------------------------

!!****f* m_gputk/gpu_copy
!! NAME
!!  gpu_copy
!!
!! FUNCTION
!!  Copy array content on GPU to another
!!
!! INPUTS
!!  src  = array to be copied
!!  size = size of src and dest
!!
!! OUTPUT
!!  dest = array to be set
!!
!! SOURCE

subroutine gpu_copy(dest, src, sizea)
 integer(c_size_t),intent(in)  :: sizea
 real(dp),target,intent(in)  :: src(sizea)
 real(dp),target,intent(out) :: dest(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(dest,src)
 call copy_gpu_to_gpu(c_loc(dest), c_loc(src), sizea*dp)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:src,dest)
 do i=1,sizea
   dest(i)=src(i)
 end do
#endif

#else
 ! Make testfarm happy
 ABI_UNUSED((/src,dest/))
#endif

end subroutine gpu_copy
!!***

!!****f* m_gputk/gpu_copy_sp
!! NAME
!!  gpu_copy_sp
!!
!! FUNCTION
!!  Copy array content on GPU to another (single precision version)
!!
!! INPUTS
!!  src  = array to be copied
!!  size = size of src and dest
!!
!! OUTPUT
!!  dest = array to be set
!!
!! SOURCE

subroutine gpu_copy_sp(dest, src, sizea)
 integer(c_size_t),intent(in)  :: sizea
 real(sp),target,intent(in)  :: src(sizea)
 real(sp),target,intent(out) :: dest(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(dest,src)
 call copy_gpu_to_gpu(c_loc(dest), c_loc(src), sizea*sp)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:src,dest)
 do i=1,sizea
   dest(i)=src(i)
 end do
#endif

#else
 ! Make testfarm happy
 ABI_UNUSED((/src,dest/))
#endif

end subroutine gpu_copy_sp
!!***

!!****f* m_gputk/gpu_copy_complex
!! NAME
!!  gpu_copy_complex
!!
!! FUNCTION
!!  Copy array content on GPU to another
!!
!! INPUTS
!!  src  = array to be copied
!!  size = size of src and dest
!!
!! OUTPUT
!!  dest = array to be set
!!
!! SOURCE

subroutine gpu_copy_complex(dest, src, sizea)
 integer(c_size_t),intent(in)  :: sizea
 complex(dp),target,intent(in)  :: src(sizea)
 complex(dp),target,intent(out) :: dest(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(dest,src)
 call copy_gpu_to_gpu(c_loc(dest), c_loc(src), sizea*dp*2)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:src,dest)
 do i=1,sizea
   dest(i)=src(i)
 end do
#endif

#else
 ! Make testfarm happy
 ABI_UNUSED((/src,dest/))
#endif

end subroutine gpu_copy_complex
!!***

!!****f* m_gputk/gpu_copy_complex_sp
!! NAME
!!  gpu_copy_complex_sp
!!
!! FUNCTION
!!  Copy array content on GPU to another (single precision version)
!!
!! INPUTS
!!  src  = array to be copied
!!  size = size of src and dest
!!
!! OUTPUT
!!  dest = array to be set
!!
!! SOURCE

subroutine gpu_copy_complex_sp(dest, src, sizea)
 integer(c_size_t),intent(in)  :: sizea
 complex(sp),target,intent(in)  :: src(sizea)
 complex(sp),target,intent(out) :: dest(sizea)
! *********************************************************************

#if defined HAVE_OPENMP_OFFLOAD
 integer(c_size_t)  :: i

#if defined HAVE_GPU_CUDA
 !$OMP TARGET DATA USE_DEVICE_ADDR(dest,src)
 call copy_gpu_to_gpu(c_loc(dest), c_loc(src), sizea*sp*2)
 !$OMP END TARGET DATA
#elif defined HAVE_GPU_HIP
 !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO PRIVATE(i) MAP(to:src,dest)
 do i=1,sizea
   dest(i)=src(i)
 end do
#endif

#else
 ! Make testfarm happy
 ABI_UNUSED((/src,dest/))
#endif

end subroutine gpu_copy_complex_sp
!!***

end module m_gputk
!!***

!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xaxpy
!! NAME
!!  abi_xaxpy
!!
!! FUNCTION
!!  abi_xaxpy is the generic function for BLAS-1 AXPY (y = alpha * x + y)
!!  with optional GPU dispatch.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2026 ABINIT group
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_daxpy
!! NAME
!! abi_daxpy
!!
!! FUNCTION
!!  Real AXPY: y = alpha * x + y, with real(dp) alpha and real(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_daxpy(n, alpha, x, incx, y, incy, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 real(dp), intent(in) :: alpha
 real(dp), target, intent(in) :: x(*)
 integer, intent(in) :: incx
 real(dp), target, intent(inout) :: y(*)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zaxpy(n, dcmplx(alpha, 0.0_dp), x, incx, y, incy, x_cplx=1, gpu_option=gpu_option_)

end subroutine abi_daxpy
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_daxpy_2d
!! NAME
!! abi_daxpy_2d
!!
!! FUNCTION
!!  Real AXPY: y = alpha * x + y, with real(dp) alpha and real(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_daxpy_2d(n, alpha, x, incx, y, incy, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 real(dp), intent(in) :: alpha
 real(dp), target, intent(in) :: x(:,:)
 integer, intent(in) :: incx
 real(dp), target, intent(inout) :: y(:,:)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zaxpy(n, dcmplx(alpha, 0.0_dp), x, incx, y, incy, x_cplx=1, gpu_option=gpu_option_)

end subroutine abi_daxpy_2d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zaxpy_2d
!! NAME
!! abi_d2zaxpy_2d
!!
!! FUNCTION
!!  Real AXPY: y = alpha * x + y, with real(dp) alpha and real(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zaxpy_2d(n, alpha, x, incx, y, incy, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 real(dp), target, intent(in) :: x(:,:)
 integer, intent(in) :: incx
 real(dp), target, intent(inout) :: y(:,:)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zaxpy(n, alpha, x, incx, y, incy, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zaxpy_2d
!!***


!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zaxpy_5d
!! NAME
!! abi_d2zaxpy_5d
!!
!! FUNCTION
!!  Real AXPY: y = alpha * x + y, with real(dp) alpha and real(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zaxpy_5d(n, alpha, x, incx, y, incy, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 real(dp), target, intent(in) :: x(:,:,:,:,:)
 integer, intent(in) :: incx
 real(dp), target, intent(inout) :: y(:,:,:,:,:)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zaxpy(n, alpha, x, incx, y, incy, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zaxpy_5d
!!***


!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zaxpy
!! NAME
!! abi_d2zaxpy
!!
!! FUNCTION
!!  Real AXPY: y = alpha * x + y, with real(dp) alpha and real(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zaxpy(n, alpha, x, incx, y, incy, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 real(dp), target, intent(in) :: x(*)
 integer, intent(in) :: incx
 real(dp), target, intent(inout) :: y(*)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(x)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(y)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
   if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
     call abi_gpu_xaxpy_cptr(cplx_, n, alpha, c_loc(x), incx, c_loc(y), incy)
     !$OMP END TARGET DATA
#endif
   else
     call abi_gpu_xaxpy_cptr(cplx_, n, alpha, c_loc(x), incx, c_loc(y), incy)
   end if
 else
   if(cplx_==1) then
     call daxpy(n, alpha, x, incx, y, incy)
   else
     call zaxpy(n, alpha, x, incx, y, incy)
   end if
 end if

end subroutine abi_d2zaxpy
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zaxpy
!! NAME
!! abi_zaxpy
!!
!! FUNCTION
!!  Complex AXPY: y = alpha * x + y, with complex(dp) alpha and complex(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zaxpy(n, alpha, x, incx, y, incy, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(in) :: x(*)
 integer, intent(in) :: incx
 complex(dp), target, intent(inout) :: y(*)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(x)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(y)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
   if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(x,y)
     call abi_gpu_xaxpy_cptr(2, n, alpha, c_loc(x), incx, c_loc(y), incy)
     !$OMP END TARGET DATA
#endif
   else
     call abi_gpu_xaxpy_cptr(2, n, alpha, c_loc(x), incx, c_loc(y), incy)
   end if
 else
   call zaxpy(n, alpha, x, incx, y, incy)
 end if

end subroutine abi_zaxpy
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zaxpy_2d
!! NAME
!! abi_zaxpy_2d
!!
!! FUNCTION
!!  Complex AXPY: y = alpha * x + y, with complex(dp) alpha and complex(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zaxpy_2d(n, alpha, x, incx, y, incy, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(in) :: x(:,:)
 integer, intent(in) :: incx
 complex(dp), target, intent(inout) :: y(:,:)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zaxpy(n, alpha, x, incx, y, incy, gpu_option=gpu_option_)

end subroutine abi_zaxpy_2d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zaxpy_3d
!! NAME
!! abi_zaxpy_3d
!!
!! FUNCTION
!!  Complex AXPY: y = alpha * x + y, with complex(dp) alpha and complex(dp) arrays.
!!  On GPU, calls abi_gpu_xaxpy_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zaxpy_3d(n, alpha, x, incx, y, incy, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(in) :: x(:,:,:)
 integer, intent(in) :: incx
 complex(dp), target, intent(inout) :: y(:,:,:)
 integer, intent(in) :: incy
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zaxpy(n, alpha, x, incx, y, incy, gpu_option=gpu_option_)

end subroutine abi_zaxpy_3d
!!***

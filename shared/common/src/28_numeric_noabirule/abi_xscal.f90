!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xscal
!! NAME
!!  abi_xscal
!!
!! FUNCTION
!!  abi_xscal is the generic function for BLAS-1 SCAL (x = alpha * x)
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

!!****f* m_abi_linalg/abi_d2zscal
!! NAME
!! abi_d2zscal
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zscal(n, alpha, x, incx, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(*)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx,gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(x)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
   if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(x)
     call abi_gpu_xscal_cptr(cplx_, n, alpha, c_loc(x), incx)
     !$OMP END TARGET DATA
#endif
   else
     call abi_gpu_xscal_cptr(cplx_, n, alpha, c_loc(x), incx)
   end if
 else
   if(cplx_==1) then
     call dscal(n, alpha, x, incx)
   else
     call zscal(n, alpha, x, incx)
   end if
 end if

end subroutine abi_d2zscal
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zscal_3d
!! NAME
!! abi_d2zscal_3d
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zscal_3d(n, alpha, x, incx, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(:,:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx,gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zscal(n, alpha, x, incx, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zscal_3d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zscal_4d
!! NAME
!! abi_d2zscal_4d
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zscal_4d(n, alpha, x, incx, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(:,:,:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx,gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zscal(n, alpha, x, incx, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zscal_4d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zscal_5d
!! NAME
!! abi_d2zscal_5d
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zscal_5d(n, alpha, x, incx, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(:,:,:,:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx,gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zscal(n, alpha, x, incx, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zscal_5d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zscal_7d
!! NAME
!! abi_d2zscal_7d
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zscal_7d(n, alpha, x, incx, x_cplx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dpc), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(:,:,:,:,:,:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx,gpu_option

!Local variables-------------------------------
 integer :: cplx_,gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option
 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx

 call abi_d2zscal(n, alpha, x, incx, x_cplx=cplx_, gpu_option=gpu_option_)

end subroutine abi_d2zscal_7d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dscal
!! NAME
!! abi_dscal
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_dscal(n, alpha, x, incx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 real(dp), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(*)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zscal(n, dcmplx(alpha, 0.0_dp), x, incx, x_cplx=1, gpu_option=gpu_option_)

end subroutine abi_dscal
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_dscal_2d
!! NAME
!! abi_dscal_2d
!!
!! FUNCTION
!!  Real SCAL: x = alpha * x, with real(dp) alpha and real(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=1.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_dscal_2d(n, alpha, x, incx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 real(dp), intent(in) :: alpha
 real(dp), target, intent(inout) :: x(:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zscal(n, dcmplx(alpha, 0.0_dp), x, incx, x_cplx=1, gpu_option=gpu_option_)

end subroutine abi_dscal_2d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zscal
!! NAME
!! abi_zscal
!!
!! FUNCTION
!!  Complex SCAL: x = alpha * x, with complex(dp) alpha and complex(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zscal(n, alpha, x, incx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(inout) :: x(*)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(x)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
   if(gpu_option_==ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(x) IF(gpu_option_==ABI_GPU_OPENMP)
     call abi_gpu_xscal_cptr(2, n, alpha, c_loc(x), incx)
     !$OMP END TARGET DATA
#endif
   else
     call abi_gpu_xscal_cptr(2, n, alpha, c_loc(x), incx)
   end if
 else
   call zscal(n, alpha, x, incx)
 end if

end subroutine abi_zscal
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zscal_2d
!! NAME
!! abi_zscal
!!
!! FUNCTION
!!  Complex SCAL: x = alpha * x, with complex(dp) alpha and complex(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zscal_2d(n, alpha, x, incx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(inout) :: x(:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zscal(n, alpha, x, incx, gpu_option=gpu_option_)

end subroutine abi_zscal_2d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zscal_3d
!! NAME
!! abi_zscal_3d
!!
!! FUNCTION
!!  Complex SCAL: x = alpha * x, with complex(dp) alpha and complex(dp) array.
!!  On GPU, calls abi_gpu_xscal_cptr with cplx=2.
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zscal_3d(n, alpha, x, incx, gpu_option)

!Arguments ------------------------------------
 integer, intent(in) :: n
 complex(dp), intent(in) :: alpha
 complex(dp), target, intent(inout) :: x(:,:,:)
 integer, intent(in) :: incx
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zscal(n, alpha, x, incx, gpu_option=gpu_option_)

end subroutine abi_zscal_3d
!!***

!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xpotrf
!! NAME
!!  abi_xpotrf
!!
!! FUNCTION
!!  abi_xpotrf is the generic function for computing the
!!  Cholesky factorization of a real symmetric (or hermitian)
!!    positive definite matrix A.
!!    The factorization has the form
!!      A = U**T * U,  if UPLO = 'U', or
!!      A = L  * L**T,  if UPLO = 'L',
!!    where U is an upper triangular matrix and L is lower triangular.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2026 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_linalg/abi_dpotrf
!! NAME
!! abi_dpotrf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_dpotrf(uplo,n,a,lda,info,gpu_option)

 !Arguments ------------------------------------
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 integer, intent(out) :: info
 real(dp), target, intent(inout) :: a(*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

 !Local variables ------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 if(gpu_option_/=ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_ADDR(a) IF(gpu_option_==ABI_GPU_OPENMP)
#endif
   call abi_gpu_xpotrf_cptr(1, uplo, n, c_loc(a), lda, info)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP END TARGET DATA
#endif
 else
#ifdef HAVE_LINALG_PLASMA
   if (ABI_LINALG_PLASMA_ISON) then
     call PLASMA_dpotrf(uplo_plasma(uplo),n,a,lda,info)
     return
   end if
#endif
   call dpotrf(uplo,n,a,lda,info)
 end if

end subroutine abi_dpotrf
!!***

!!****f* m_abi_linalg/abi_zpotrf_2d
!! NAME
!! abi_zpotrf_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zpotrf_2d(uplo,n,a,lda,info,gpu_option)

 !Arguments ------------------------------------
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: lda,n
 integer, intent(out) :: info
 complex(dp), target, intent(inout) :: a(lda,*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

! *********************************************************************

 call abi_zpotrf(uplo,n,a(1,1),lda,info,gpu_option=gpu_option)

end subroutine abi_zpotrf_2d
!!***

!!****f* m_abi_linalg/abi_d2zpotrf_3d
!! NAME
!! abi_d2zpotrf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zpotrf_3d(uplo,n,a,lda,info,x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 integer, intent(out) :: info
 integer, intent(in), optional :: x_cplx
 integer, intent(in), optional :: gpu_option
 real(dp),target, intent(inout) :: a(:,:,:)

 !Local Variables -----------------------------
 integer  :: cplx_, gpu_option_

! *********************************************************************

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zpotrf(uplo,n,a,lda,info,x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zpotrf_3d
!!***

!!****f* m_abi_linalg/abi_d2zpotrf
!! NAME
!! abi_d2zpotrf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_d2zpotrf(uplo,n,a,lda,info,x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: n,lda
 integer, intent(out) :: info
 integer, intent(in), optional :: x_cplx
 integer, intent(in), optional :: gpu_option
 real(dp),target, intent(inout) :: a(lda,*)  ! FIXME should be x_cplx * lda

 !Local Variables -----------------------------
 integer  :: cplx_, gpu_option_

! *********************************************************************

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 if(gpu_option_/=ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_ADDR(a) IF(gpu_option_==ABI_GPU_OPENMP)
#endif
   call abi_gpu_xpotrf_cptr(cplx_, uplo, n, c_loc(a), lda, info)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP END TARGET DATA
#endif
 else
#ifdef HAVE_LINALG_PLASMA
   if (ABI_LINALG_PLASMA_ISON) then
     if(cplx_ == 2) then
        info = PLASMA_zpotrf_c(uplo_plasma(uplo),n,c_loc(a),lda)
     else
        info = PLASMA_dpotrf_c(uplo_plasma(uplo),n,c_loc(a),lda)
     end if
     return
   end if
#endif
   if(cplx_ == 2) then
      call zpotrf(uplo,n,a,lda,info)
   else
      call dpotrf(uplo,n,a,lda,info)
   end if
 end if

end subroutine abi_d2zpotrf
!!***

!!****f* m_abi_linalg/abi_zpotrf
!! NAME
!! abi_zpotrf
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine abi_zpotrf(uplo,n,a,lda,info,gpu_option)

 !Arguments ------------------------------------
 character(len=1), intent(in) :: uplo
 integer, intent(in) :: lda,n
 integer, intent(out) :: info
 complex(dp), target, intent(inout) :: a(*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

 !Local variables ------------------------------
 integer :: gpu_option_

! *********************************************************************

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 if(gpu_option_/=ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_ADDR(a) IF(gpu_option_==ABI_GPU_OPENMP)
#endif
   call abi_gpu_xpotrf_cptr(2, uplo, n, c_loc(a), lda, info)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP END TARGET DATA
#endif
 else
#ifdef HAVE_LINALG_PLASMA
   if (ABI_LINALG_PLASMA_ISON) then
     call PLASMA_zpotrf(uplo_plasma(uplo),n,a,lda,info)
     return
   end if
#endif
   call zpotrf(uplo,n,a,lda,info)
 end if

end subroutine abi_zpotrf
!!***

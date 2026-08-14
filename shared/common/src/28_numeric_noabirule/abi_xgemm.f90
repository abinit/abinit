!{\src2tex{textfont=tt}}
!!****f* m_abi_linalg/abi_xgemm
!! NAME
!!  abi_xgemm
!!
!! FUNCTION
!!  abi_xgemm is the generic function that solves:
!!
!!    C := alpha*op( A )*op( B ) + beta*C,
!!
!!  where  op( X ) is one of
!!
!!    op( X ) = X   or   op( X ) = X**T,
!!
!! alpha and beta are scalars, and A, B and C are matrices, with op( A )
!! an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
!!
!! COPYRIGHT
!!  Copyright (C) 2001-2026 ABINIT group (LNguyen,FDahm (CS))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

!!****f* m_abi_gpu_linalg/abi_xgemm_gpu_cptr
!! NAME
!!  abi_gpu_xgemm_gpu_cptr
!!
!! FUNCTION
!!  Compute a scalar-matrix-matrix product and return a scalar-matrix product on GPU
!!  c = alpha * op(a) * op(b) + beta * c
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  transa= from of op(a) to be used in the matrix multiplication
!!  transb= from of op(b) to be used in the matrix multiplication
!!  m     = number of rows of the matrix op(a) and of the matrix c
!!  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
!!  k     = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
!!  alpha = alpha scalar coefficient for matrix op(a)
!!  a = pointer to gpu memory location of  matrix a
!!  lda   = first dimension of a
!!  b = pointer to gpu memory location of  matrix b
!!  ldb   = first dimension of b
!!  beta  = beta scalar coefficient for matrix c
!!  c = pointer to gpu memory location of  matrix c
!!  ldc   = first dimension of c
!!
!! OUTPUT
!!  c     = c matrix
!!
!! SOURCE

subroutine abi_xgemm_gpu_cptr(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dp),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 type(c_ptr),intent(in) :: a,b
 type(c_ptr),intent(in) :: c
! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU
  call gpu_xgemm(cplx,transa,transb,m,n,k,alpha,&
      a,lda,b,ldb,beta,c,ldc)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    ! potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,lda,ldb,ldc,m,n,k/))
  ABI_UNUSED((/alpha,beta/))
  ABI_UNUSED((/transa,transb/))
  ABI_UNUSED((/a,b,c/))
#endif

end subroutine abi_xgemm_gpu_cptr
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zgemm_2dd
!! NAME
!! abi_zgemm_2dd
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

 subroutine abi_zgemm_2dd(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
     gpu_option)

 !Arguments------------------------------------
 character(len=1),intent(in) :: TRANSA
 character(len=1),intent(in) :: TRANSB
 integer,intent(in) :: K
 integer,intent(in) :: LDA
 integer,intent(in) :: LDB
 integer,intent(in) :: LDC
 integer,intent(in) :: M
 integer,intent(in) :: N
 complex(dp),intent(in) :: ALPHA
 complex(dp),intent(in) :: BETA
 complex(dp),target,intent(in) :: A(*)
 complex(dp),target,intent(in) :: B(*)
 complex(dp),target,intent(inout) :: C(*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

 !Local variables-------------------------------
 integer  :: gpu_option_,info
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(a)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(b)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(c)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
   if(gpu_option == ABI_GPU_LEGACY .or. gpu_option == ABI_GPU_KOKKOS) then
     call abi_xgemm_gpu_cptr(2,transa,transb,m,n,k,alpha,&
         c_loc(a),lda,&
         c_loc(b),ldb,&
         beta,&
         c_loc(c),ldc)
   else if(gpu_option == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
     !$OMP TARGET DATA USE_DEVICE_ADDR(a,b,c)
     call abi_xgemm_gpu_cptr(2,transa,transb,m,n,k,alpha,&
       c_loc(a),lda,&
       c_loc(b),ldb,&
       beta,&
       c_loc(c),ldc)
     !$OMP END TARGET DATA
#endif
   else
     ABI_BUG(sjoin("Unhandled GPU mode:", itoa(gpu_option)))
   end if
 else if (ABI_LINALG_PLASMA_ISON) then
   info = -1
#ifdef HAVE_LINALG_PLASMA
   !write(std_out,*)"Will call PLASMA_zgemm_c"
   info = PLASMA_zgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,ALPHA,&
&    c_loc(A),LDA,c_loc(B),LDB,BETA,c_loc(C),LDC)
#endif
   ABI_CHECK(info==0,"PLASMA_zgemm_c returned info !=0")
 else
   if (use_zgemm3m(m,n,k)) then
     call _ZGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   else
     call zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XGEMM,2,tsec)
#endif

end subroutine abi_zgemm_2dd
!!***

!----------------------------------------------------------------------
 subroutine abi_zgemm_2d(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
     gpu_option)

 !Arguments------------------------------------
 character(len=1),intent(in) :: TRANSA
 character(len=1),intent(in) :: TRANSB
 integer,intent(in) :: K
 integer,intent(in) :: LDA
 integer,intent(in) :: LDB
 integer,intent(in) :: LDC
 integer,intent(in) :: M
 integer,intent(in) :: N
 complex(dp),intent(in) :: ALPHA
 complex(dp),intent(in) :: BETA
 complex(dp),target,intent(in) :: A(lda,*)
 complex(dp),target,intent(in) :: B(ldb,*)
 complex(dp),target,intent(inout) :: C(ldc,*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

 !Local variables-------------------------------
 integer  :: gpu_option_
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zgemm_2dd(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    gpu_option=gpu_option_)

end subroutine abi_zgemm_2d
!!***

!----------------------------------------------------------------------
 subroutine abi_zgemm_3d(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
     gpu_option)

 !Arguments------------------------------------
 character(len=1),intent(in) :: TRANSA
 character(len=1),intent(in) :: TRANSB
 integer,intent(in) :: K
 integer,intent(in) :: LDA
 integer,intent(in) :: LDB
 integer,intent(in) :: LDC
 integer,intent(in) :: M
 integer,intent(in) :: N
 complex(dp),intent(in) :: ALPHA
 complex(dp),intent(in) :: BETA
 complex(dp),target,intent(in) :: A(:,:)
 complex(dp),target,intent(in) :: B(:,:,:)
 complex(dp),target,intent(inout) :: C(:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

 !Local variables-------------------------------
 integer  :: gpu_option_
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_zgemm_2dd(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    gpu_option=gpu_option_)

end subroutine abi_zgemm_3d
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm
!! NAME
!! abi_d2zgemm
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target,intent(in)    :: a(*)    ! FIXME should be lda * x_cplx
 real(dp),target,intent(in)    :: b(*)
 real(dp),target,intent(inout) :: c(*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_,info
#ifdef DEV_LINALG_TIMING
 real(dp) :: tsec(2)
 call timab(TIMAB_XGEMM,1,tsec)
#endif

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

#if defined(DEBUG_VERBOSE) && defined(HAVE_OPENMP_OFFLOAD)
 if ( gpu_option_ == ABI_GPU_OPENMP ) then
   ABI_CHECK(xomp_target_is_present(c_loc(a)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(b)), "Array isn't mapped on GPU")
   ABI_CHECK(xomp_target_is_present(c_loc(c)), "Array isn't mapped on GPU")
 end if
#endif

 if(gpu_option_/=ABI_GPU_DISABLED) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_ADDR(a,b,c) IF(gpu_option_==ABI_GPU_OPENMP)
#endif
   call abi_xgemm_gpu_cptr(cplx_,transa,transb,m,n,k,alpha,&
       c_loc(a),lda,&
       c_loc(b),ldb,&
       beta,&
       c_loc(c),ldc)
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP END TARGET DATA
#endif
 else if (ABI_LINALG_PLASMA_ISON) then
   info = -1
#ifdef HAVE_LINALG_PLASMA
   !write(std_out,*) "Will call plasma_[zd]gemm"
   if (cplx_ == 2) then
      info = PLASMA_zgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,&
&      ALPHA,A,LDA,B,LDB,BETA,C,LDC)
   else
      info = PLASMA_dgemm_c(trans_plasma(TRANSA),trans_plasma(TRANSB),M,N,K,&
&      real(ALPHA,dp),A,LDA,B,LDB,real(BETA,dp),C,LDC)
   end if
#endif
   ABI_CHECK(info==0,"PLASMA_[z,d]gemm_c returned info !=0")
 else
   if (cplx_ == 2) then
     if (use_zgemm3m(m,n,k)) then
       call _ZGEMM3M(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     else
       call zgemm(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC)
     end if
   else
     call dgemm(TRANSA,TRANSB,M,N,K,real(ALPHA,dp),A,LDA,B,LDB,real(BETA,dp),C,LDC)
   end if
 end if

#ifdef DEV_LINALG_TIMING
 call timab(TIMAB_XGEMM,2,tsec)
#endif

end subroutine abi_d2zgemm
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_zgemm_2r
!! NAME
!! abi_zgemm_2r
!!
!! FUNCTION
!!
!! INPUTS
!!
!!
!! SOURCE

 subroutine abi_zgemm_2r(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
&                        gpu_option)

 !Arguments------------------------------------
 character(len=1),intent(in) :: TRANSA
 character(len=1),intent(in) :: TRANSB
 integer,intent(in) :: K
 integer,intent(in) :: LDA
 integer,intent(in) :: LDB
 integer,intent(in) :: LDC
 integer,intent(in) :: M
 integer,intent(in) :: N
 complex(dp),intent(in) :: ALPHA
 complex(dp),intent(in) :: BETA
 real(dp),target,intent(in) :: A(*)
 real(dp),target,intent(in) :: B(*)
 real(dp),target,intent(inout) :: C(*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: gpu_option

!Local variables-------------------------------
 integer  :: gpu_option_

 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=2,gpu_option=gpu_option_)

end subroutine abi_zgemm_2r
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_222
!! NAME
!! abi_d2zgemm_333
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_222(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:)
 real(dp),target, intent(in)    :: b(:,:)
 real(dp),target, intent(inout) :: c(:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_222
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_233
!! NAME
!! abi_d2zgemm_333
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_233(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:)
 real(dp),target, intent(in)    :: b(:,:,:)
 real(dp),target, intent(inout) :: c(:,:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_233
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_333
!! NAME
!! abi_d2zgemm_333
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_333(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:,:)
 real(dp),target, intent(in)    :: b(:,:,:)
 real(dp),target, intent(inout) :: c(:,:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_333
!!***


!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_313
!! NAME
!! abi_d2zgemm_313
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_313(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:,:)    ! FIXME should be lda * x_cplx
 real(dp),target, intent(in)    :: b(:)
 real(dp),target, intent(inout) :: c(:,:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_313
!!***


!!****f* m_abi_linalg/abi_d2zgemm_331
!! NAME
!! abi_d2zgemm_331
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_331(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:,:)    ! FIXME should be lda * x_cplx
 real(dp),target, intent(in)    :: b(:,:,:)
 real(dp),target, intent(inout) :: c(:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_331
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_331
!! NAME
!! abi_d2zgemm_331
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_334(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(:,:,:)    ! FIXME should be lda * x_cplx
 real(dp),target, intent(in)    :: b(:,:,:)
 real(dp),target, intent(inout) :: c(:,:,:,:)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_334
!!***

!----------------------------------------------------------------------

!!****f* m_abi_linalg/abi_d2zgemm_d
!! NAME
!! abi_d2zgemm_2d
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE
!!
subroutine abi_d2zgemm_2d(TRANSA,TRANSB,M,N,K,ALPHA,A,LDA,B,LDB,BETA,C,LDC,&
                       x_cplx,gpu_option)

!Arguments ------------------------------------
 character(len=1), intent(in) :: transa
 character(len=1), intent(in) :: transb
 integer, intent(in) :: lda
 integer, intent(in) :: ldb
 integer, intent(in) :: ldc
 integer, intent(in) :: m
 integer, intent(in) :: n
 integer, intent(in) :: k
 complex(dp),intent(in) :: alpha
 complex(dp),intent(in) :: beta
 real(dp),target, intent(in)    :: a(lda,*)    ! FIXME should be lda * x_cplx
 real(dp),target, intent(in)    :: b(ldb,*)
 real(dp),target, intent(inout) :: c(ldc,*)
 !Optionals -----------------------------------
 integer, intent(in), optional :: x_cplx, gpu_option

!Local variables-------------------------------
 integer  :: cplx_,gpu_option_

 cplx_=1 ; if(PRESENT(x_cplx)) cplx_ = x_cplx
 gpu_option_=ABI_GPU_DISABLED ; if(PRESENT(gpu_option)) gpu_option_ = gpu_option

 call abi_d2zgemm(transa,transb,m,n,k,alpha,&
 &    a,lda,&
 &    b,ldb,&
 &    beta,&
 &    c,ldc,&
 &    x_cplx=cplx_,gpu_option=gpu_option_)

end subroutine abi_d2zgemm_2d
!!***


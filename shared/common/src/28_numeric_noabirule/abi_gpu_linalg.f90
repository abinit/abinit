!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_abi_gpu_linalg
!! NAME
!!  m_abi_gpu_linalg
!!
!! FUNCTION
!!  Interfaces of GPU subroutines wrapper
!!
!! COPYRIGHT
!!  Copyright (C) 2011-2024 ABINIT group (FDahm ))
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~ABINIT/Infos/copyright
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

!!***

#ifndef HAVE_GPU

!!****f* m_abi_gpu_linalg/gpu_device_synchronize
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


!!****f* m_abi_gpu_linalg/check_gpu_mem
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


!!****f* m_abi_gpu_linalg/alloc_on_gpu
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

!!****f* m_abi_gpu_linalg/copy_from_gpu
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

!!****f* m_abi_gpu_linalg/copy_on_gpu
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

!!****f* m_abi_gpu_linalg/copy_gpu_to_gpu
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
!! PARENTS
!!      lobpcgwf,m_abi_gpu_linalg
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

!!****f* m_abi_gpu_linalg/dealloc_on_gpu
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

!!****f* m_abi_gpu_linalg/gpu_memset
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
!! PARENTS
!!      lobpcgwf
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

!!****f* m_abi_gpu_linalg/gpu_allocated_impl
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
!! PARENTS
!!      lobpcgwf
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

!!****f* m_abi_gpu_linalg/gpu_managed_ptr_status
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
!! PARENTS
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

!!****f* m_abi_gpu_linalg/gpu_linalg_init
!! NAME
!!  gpu_linalg_init
!!
!! FUNCTION
!!  initialisation of linalg environnement on GPU
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/gpu_linalg.cu
!!
!! SOURCE

subroutine gpu_linalg_init()


end subroutine gpu_linalg_init
!!***

!!****f* m_abi_gpu_linalg/gpu_linalg_shutdown
!! NAME
!!  gpu_linalg_shutdown
!!
!! FUNCTION
!!  close linalg environnement on GPU
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/gpu_linalg.cu
!!
!! SOURCE
subroutine gpu_linalg_shutdown()

end subroutine gpu_linalg_shutdown
!!***

!!****f* m_abi_gpu_linalg/gpu_xgemm
!! NAME
!!  gpu_xgemm
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
!!  a_gpu = pointer to gpu memory location of  matrix a
!!  lda   = first dimension of a
!!  b_gpu = pointer to gpu memory location of  matrix b
!!  ldb   = first dimension of b
!!  beta  = beta scalar coefficient for matrix c
!!  c_gpu = pointer to gpu memory location of  matrix c
!!  ldc   = first dimension of c
!!
!! OUTPUT
!!  c     = c matrix
!!
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/gpu_linalg.cu
!!
!! SOURCE

subroutine gpu_xgemm(cplx,transa,transb,m,n,k,alpha,a_gpu,lda,b_gpu,ldb,beta,c_gpu,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 type(c_ptr),intent(in) :: a_gpu,b_gpu
 type(c_ptr),intent(inout) :: c_gpu
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=a_gpu;cptr=b_gpu;cptr=c_gpu
   write(std_out,*) transa,transb,cplx,lda,ldb,ldc,m,n,k,alpha,beta
 end if

end subroutine gpu_xgemm
!!***

!!****f* m_abi_gpu_linalg/gpu_xtrsm
!! NAME
!!  gpu_xtrsm
!!
!! FUNCTION
!! Solves a matrix equation (one matrix operand is triangular) on GPU.
!! The xtrsm routines solve one of the following matrix equations
!! op(a)*x = alpha*b
!! or
!! x*op(a) = alpha*b,
!!
!! INPUTS
!! cplx= 1 if real 2 if complex
!! side= Specifies whether op(a) appears on the left or right of x for
!!      the operation to be performed as follows:
!!      L or l op(a)*x = alpha*b
!!      R or r x*op(a) = alpha*b
!! uplo= Specifies whether the matrix a is an upper or lower triangular
!!      matrix as follows:
!!      U or u Matrix a is an upper triangular matrix.
!!      L or l Matrix a is a lower triangular matrix
!! transa= Specifies the form of op(a) to be used in the matrix
!!      multiplication as follows:
!!      N or n op(a) = a
!!      T or t op(a) = a'
!!      C or c op(a) = conjg(a')
!! diag= Specifies whether or not a is unit triangular as follows:
!!      U or u Matrix a is assumed to be unit triangular.
!!      N or n Matrix a is not assumed to be unit triangular.
!! m= Specifies the number of rows of b. The value of m must be at least zero
!! n= Specifies the number of columns of b. The value of n must be at least zero
!! alpha= Specifies the scalar alpha. When alpha is zero, then a is not referenced and b
!!      need not be set before entry.
!!  a_gpu = pointer to gpu memory location of  array a, DIMENSION (lda, k), where k is m when side = 'L' or 'l' and is n
!!      when side = 'R' or 'r'.
!! lda= Specifies the first dimension of a as declared in the calling
!!     (sub)program. When side = 'L' or 'l', then lda must be at least max(1,
!!      m), when side = 'R' or 'r', then lda must be at least max(1, n).
!!  b_gpu = pointer to gpu memory location of  b Array, DIMENSION (ldb,n). Before entry, the leading m-by-n part of the array
!!     b must contain the right-hand side matrix b.
!! ldb= Specifies the first dimension of b as declared in the calling
!!     (sub)program. The value of ldb must be at least max(1, m).
!!
!! OUTPUT
!!  b_gpu
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU is not enabled
!!   the correct one is in 17_gpu_toolbox/gpu_linalg.cu
!!
!! SOURCE
subroutine gpu_xtrsm(cplx,side,uplo,transa,diag,m,n,alpha,a_gpu,lda,b_gpu,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 type(c_ptr),intent(in) :: a_gpu
 type(c_ptr),intent(inout) :: b_gpu
!Local variables ------------------------------
 type(c_ptr) :: cptr

! *********************************************************************

 if (.false.) then
   cptr=a_gpu;cptr=b_gpu
   write(std_out,*) side,uplo,transa,diag,cplx,lda,ldb,m,n,alpha
 end if

end subroutine gpu_xtrsm
!!***

!!****f* m_abi_gpu_linalg/gpu_xaxpy
!! NAME
!!  gpu_xaxpy
!!
!! FUNCTION
!!  Compute a BLAS-1 AXPY operation on GPU
!!  y = alpha * x + y
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  size  = vector size
!!  alpha = scalar complex value
!!  x_gpu = pointer to gpu memory location of array x
!! incrx  = stride between consecutive elements of x
!!  y_gpu = pointer to gpu memory location of array y
!! incry  = stride between consecutive elements of y
!!
!! SOURCE
subroutine gpu_xaxpy(cplx, size, alpha, x_gpu, incrx, y_gpu, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  type(c_ptr),  intent(in)    :: x_gpu
  integer,      intent(in)    :: incrx
  type(c_ptr),  intent(inout) :: y_gpu
  integer,      intent(in)    :: incry

  ABI_UNUSED((/cplx,size,incrx,incry/))
  ABI_UNUSED(alpha)
  ABI_UNUSED_A(x_gpu)
  ABI_UNUSED_A(y_gpu)
end subroutine gpu_xaxpy
!!***

!!****f* m_abi_gpu_linalg/gpu_xcopy
!! NAME
!!  gpu_xcopy
!!
!! FUNCTION
!!  Compute a BLAS-1 COPY operation on GPU
!!  y = x (copy x into y)
!!
!! INPUTS
!!  cplx   = 1 if real 2 if complex
!!  size   = input vector size
!!  x_gpu  = pointer to gpu memory location of array x
!!  incrx  = stride between consecutive elements of x
!!  y_gpu  = pointer to gpu memory location of array y
!!  incry  = stride between consecutive elements of y
!!
!! SOURCE
subroutine gpu_xcopy(cplx, size, x_gpu, incrx, y_gpu, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  type(c_ptr),  intent(in)    :: x_gpu
  integer,      intent(in)    :: incrx
  type(c_ptr),  intent(inout) :: y_gpu
  integer,      intent(in)    :: incry

  ABI_UNUSED((/cplx,size,incrx,incry/))
  ABI_UNUSED_A(x_gpu)
  ABI_UNUSED_A(y_gpu)
end subroutine gpu_xcopy
!!***

!!****f* m_abi_gpu_linalg/gpu_xscal
!! NAME
!!  gpu_xscal
!!
!! FUNCTION
!!  Compute a BLAS-1 SCAL operation on GPU
!!  x = alpha * x
!!
!! INPUTS
!!  cplx   = 1 if real 2 if complex
!!  size   = vector size
!!  alpha  = scalar complex value
!!  x_gpu  = pointer to gpu memory location of array x
!!  incrx  = stride between consecutive elements of x
!!
!! SOURCE
subroutine gpu_xscal(cplx, size, alpha, x_gpu, incrx)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  type(c_ptr),  intent(in)    :: x_gpu
  integer,      intent(in)    :: incrx

  ABI_UNUSED((/cplx,size,incrx/))
  ABI_UNUSED(alpha)
  ABI_UNUSED_A(x_gpu)
end subroutine gpu_xscal
!!***

!!****f* m_abi_gpu_linalg/gpu_xsygvd
!! NAME
!!  gpu_xsygvd
!!
!! FUNCTION
!!  Compute a LAPACK SYGVD operation on GPU
!!  compute eigen values/vectors of a real generalized
!!  symmetric-definite eigenproblem
!!
!! See cusolver documentation
!! https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-sygvd
!!
!! See also LAPACK doc in reference implementation:
!! https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dsygvd.f
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  itype = integer, type of problem
!!  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
!!  uplo  = character, 'u' or 'l'
!!  A_nrows = matrix size
!!  A_ptr = pointer to gpu memory location of matrix A
!!  lda   = leading dimension of matrix A
!!  B_ptr = pointer to gpu memory location of matrix B
!!  ldb   = leading dimension of matrix B
!!  W_ptr = pointer to gpu memory location of matrix W (output eigen values)
!!  work_ptr =
!!  lwork =
!!  devInfo  =
!!
!! SOURCE
subroutine gpu_xsygvd(cplx, itype, jobz, uplo, A_nrows, &
  &                   A_ptr, lda, &
  &                   B_ptr, ldb, &
  &                   W_ptr,      &
  &                   work_ptr, lwork, &
  &                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  type(c_ptr),     intent(in   ) :: A_ptr
  integer,         intent(in   ) :: lda
  type(c_ptr),     intent(in   ) :: B_ptr
  integer,         intent(in   ) :: ldb
  type(c_ptr),     intent(inout) :: W_ptr
  type(c_ptr),     intent(inout) :: work_ptr
  integer,         intent(in   ) :: lwork
  integer,         intent(inout) :: devInfo

  ABI_UNUSED((/cplx,itype,A_nrows,lda,ldb,lwork,devInfo/))
  ABI_UNUSED((/jobz,uplo/))
  ABI_UNUSED((/A_ptr,B_ptr,W_ptr,work_ptr/))
end subroutine gpu_xsygvd
!!***

!!****f* m_abi_gpu_linalg/gpu_xsygvd_bufferSize
!! NAME
!!  gpu_xsygvd_bufferSize
!!
!! FUNCTION
!!  Compute required size for auxiliary work buffer used internally by cusolver
!!  in cusolverDnDsygvd / cusolverDnZhegvd
!!
!! See cusolver documentation
!! https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-sygvd
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  itype = integer, type of problem
!!  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
!!  uplo  = character, 'u' or 'l'
!!  A_nrows = matrix size
!!  A_ptr = pointer to gpu memory location of matrix A
!!  lda   = leading dimension of matrix A
!!  B_ptr = pointer to gpu memory location of matrix B
!!  ldb   = leading dimension of matrix B
!!  W_ptr = pointer to gpu memory location of matrix W (output eigen values)
!!  work_ptr =
!!  lwork =
!!  devInfo  =
!!
!! SOURCE
subroutine gpu_xsygvd_buffersize(cplx, itype, jobz, uplo, A_nrows, &
  &                              A_ptr, lda, &
  &                              B_ptr, ldb, &
  &                              W_ptr,      &
  &                              lwork)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  type(c_ptr),     intent(in   ) :: A_ptr
  integer,         intent(in   ) :: lda
  type(c_ptr),     intent(in   ) :: B_ptr
  integer,         intent(in   ) :: ldb
  type(c_ptr),     intent(inout) :: W_ptr
  integer,         intent(in   ) :: lwork

  ABI_UNUSED((/cplx,itype,A_nrows,lda,ldb,lwork/))
  ABI_UNUSED((/jobz,uplo/))
  ABI_UNUSED((/A_ptr,B_ptr,W_ptr/))
end subroutine gpu_xsygvd_bufferSize
!!***


#endif


!------------------------------------------------------------------------------
!                         abi_gpu_xgemm
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xgemm
!! NAME
!!  abi_gpu_xgemm
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

subroutine abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
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
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,lda,ldb,ldc,m,n,k/))
  ABI_UNUSED((/alpha,beta/))
  ABI_UNUSED((/transa,transb/))
  ABI_UNUSED((/a,b,c/))
#endif

end subroutine abi_gpu_xgemm_cptr
!!***

subroutine abi_gpu_xgemm_d(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),   intent(in),target :: a(*),b(*)
 real(dp),   intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_d
!!***

subroutine abi_gpu_xgemm_z(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 complex(dpc),intent(in),target :: a(*),b(*)
 complex(dpc),intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_z
!!***

subroutine abi_gpu_xgemm_2d(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),   intent(in),target :: a(lda,*),b(ldb,*)
 real(dp),   intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_2d
!!***

subroutine abi_gpu_xgemm_2z(cplx,transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 complex(dpc),intent(in),target :: a(lda,*),b(ldb,*)
 complex(dpc),intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xgemm_strided
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xgemm_strided
!! NAME
!!  abi_gpu_xgemm_strided
!!
!! FUNCTION
!!  Compute a batched scalar-matrix-matrix product and return a scalar-matrix product on GPU.
!!  Meant to be used on non-contiguous matrixes with data is uniformly split in the same number of batches in each matrix.
!!  c = alpha * op(a) * op(b) + beta * c
!!
!! INPUTS
!!  cplx       = 1 if real 2 if complex
!!  transa     = from of op(a) to be used in the matrix multiplication
!!  transb     = from of op(b) to be used in the matrix multiplication
!!  m          = number of rows of the matrix op(a) and of the matrix c
!!  n          = number of rows of the matrix op(b) and the number of columns of the matrix c
!!  k          = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
!!  alpha      = alpha scalar coefficient for matrix op(a)
!!  a          = pointer to gpu memory location of  matrix a
!!  lda        = first dimension of a
!!  strideC    = stride between each batch in matrix a
!!  b          = pointer to gpu memory location of  matrix b
!!  ldb        = first dimension of b
!!  strideC    = stride between each batch in matrix b
!!  beta       = beta scalar coefficient for matrix c
!!  c          = pointer to gpu memory location of  matrix c
!!  ldc        = first dimension of c
!!  strideC    = stride between each batch in matrix c
!!  batchCount = number of batches in any matrix
!!
!! OUTPUT
!!  c     = c matrix
!!
!! SOURCE

subroutine abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 integer,intent(in) :: strideA,strideB,strideC,batchCount
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 type(c_ptr),intent(in) :: a,b
 type(c_ptr),intent(in) :: c

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_xgemm_strided_batched(cplx,transa,transb,m,n,k,alpha,&
      a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,lda,ldb,ldc,m,n,k/))
  ABI_UNUSED((/strideA,strideB,strideC,batchCount/))
  ABI_UNUSED((/alpha,beta/))
  ABI_UNUSED((/transa,transb/))
  ABI_UNUSED((/a,b,c/))
#endif

end subroutine abi_gpu_xgemm_strided_cptr
!!***

subroutine abi_gpu_xgemm_strided_d(cplx,transa,transb,m,n,k,alpha,a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 integer,intent(in) :: strideA,strideB,strideC,batchCount
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),   intent(in),target :: a(*),b(*)
 real(dp),   intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,strideA,&
        c_loc(b),ldb,strideB,&
        beta,&
        c_loc(c),ldc,strideC,batchCount)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,strideA,&
        c_loc(b),ldb,strideB,&
        beta,&
        c_loc(c),ldc,strideC,batchCount)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_strided_d
!!***

subroutine abi_gpu_xgemm_strided_z(cplx,transa,transb,m,n,k,alpha,a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 integer,intent(in) :: strideA,strideB,strideC,batchCount
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 complex(dpc),intent(in),target :: a(*),b(*)
 complex(dpc),intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,strideA,&
        c_loc(b),ldb,strideB,&
        beta,&
        c_loc(c),ldc,strideC,batchCount)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,strideA,&
      c_loc(b),ldb,strideB,&
      beta,&
      c_loc(c),ldc,strideC,batchCount)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_strided_z
!!***

subroutine abi_gpu_xgemm_strided_2d(cplx,transa,transb,m,n,k,alpha,a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 integer,intent(in) :: strideA,strideB,strideC,batchCount
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 real(dp),   intent(in),target :: a(lda,*),b(ldb,*)
 real(dp),   intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,strideA,&
        c_loc(b),ldb,strideB,&
        beta,&
        c_loc(c),ldc,strideC,batchCount)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,strideA,&
      c_loc(b),ldb,strideB,&
      beta,&
      c_loc(c),ldc,strideC,batchCount)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_strided_2d
!!***

subroutine abi_gpu_xgemm_strided_2z(cplx,transa,transb,m,n,k,alpha,a,lda,strideA,b,ldb,strideB,beta,c,ldc,strideC,batchCount)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n,k
 integer,intent(in) :: strideA,strideB,strideC,batchCount
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: transa,transb
 complex(dpc),intent(in),target :: a(lda,*),b(ldb,*)
 complex(dpc),intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
        c_loc(a),lda,strideA,&
        c_loc(b),ldb,strideB,&
        beta,&
        c_loc(c),ldc,strideC,batchCount)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xgemm_strided_cptr(cplx,transa,transb,m,n,k,alpha,&
      c_loc(a),lda,strideA,&
      c_loc(b),ldb,strideB,&
      beta,&
      c_loc(c),ldc,strideC,batchCount)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xgemm_strided_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xsymm
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xsymm
!! NAME
!!  abi_gpu_xsymm
!!
!! FUNCTION
!!  Compute a symmetric scalar-matrix-matrix product and return a scalar-matrix product on GPU
!!  c = alpha * op(a) * op(b) + beta * c    , if side == L
!!  c = alpha * op(b) * op(a) + beta * c    , if side == R
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  side  = Specifies whether op(a) appears on the left or right of x for
!!          the operation to be performed as follows:
!!              L or l op(a)*x = alpha*b
!!              R or r x*op(a) = alpha*b
!!  uplo  = Specifies whether the matrix a is an upper or lower triangular
!!          matrix as follows:
!!              U or u Matrix a is an upper triangular matrix.
!!              L or l Matrix a is a lower triangular matrix
!!  m     = number of rows of the matrix op(a) and of the matrix c
!!  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
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

subroutine abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 type(c_ptr),intent(in) :: a,b
 type(c_ptr),intent(in) :: c

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_xsymm(cplx,side,uplo,m,n,alpha,&
      a,lda,b,ldb,beta,c,ldc)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,lda,ldb,ldc,m,n/))
  ABI_UNUSED((/alpha,beta/))
  ABI_UNUSED((/side,uplo/))
  ABI_UNUSED((/a,b,c/))
#endif

end subroutine abi_gpu_xsymm_cptr
!!***

subroutine abi_gpu_xsymm_d(cplx,side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 real(dp),   intent(in),target :: a(*),b(*)
 real(dp),   intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xsymm_d
!!***

subroutine abi_gpu_xsymm_z(cplx,side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 complex(dpc),intent(in),target :: a(*),b(*)
 complex(dpc),intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xsymm_z
!!***

subroutine abi_gpu_xsymm_2d(cplx,side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 real(dp),   intent(in),target :: a(lda,*),b(ldb,*)
 real(dp),   intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xsymm_2d
!!***

subroutine abi_gpu_xsymm_2z(cplx,side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: cplx,lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 complex(dpc),intent(in),target :: a(lda,*),b(ldb,*)
 complex(dpc),intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_xsymm_cptr(cplx,side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xsymm_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_zhemm
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_zhemm
!! NAME
!!  abi_gpu_xhemm
!!
!! FUNCTION
!!  Compute a Hermitian scalar-matrix-matrix product and return a scalar-matrix product on GPU
!!  c = alpha * op(a) * op(b) + beta * c    , if side == L
!!  c = alpha * op(b) * op(a) + beta * c    , if side == R
!!
!! INPUTS
!!  side  = Specifies whether op(a) appears on the left or right of x for
!!          the operation to be performed as follows:
!!              L or l op(a)*x = alpha*b
!!              R or r x*op(a) = alpha*b
!!  uplo  = Specifies whether the matrix a is an upper or lower triangular
!!          matrix as follows:
!!              U or u Matrix a is an upper triangular matrix.
!!              L or l Matrix a is a lower triangular matrix
!!  m     = number of rows of the matrix op(a) and of the matrix c
!!  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
!!  alpha = alpha scalar coefficient for matrix op(a)
!!  a_gpu = pointer to gpu memory location of  matrix a
!!  lda   = first dimension of a
!!  b_gpu = pointer to gpu memory location of  matrix b
!!  ldb   = first dimension of b
!!  beta  = beta scalar coefficient for matrix c
!!  c_gpu = pointer to gpu memory location of  matrix c
!!  ldc   = first dimension of c
!!
!! OUTPUT
!!  c     = c matrix
!!
!! SOURCE

subroutine abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 type(c_ptr),intent(in) :: a,b
 type(c_ptr),intent(in) :: c

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_zhemm(side,uplo,m,n,alpha,&
      a,lda,b,ldb,beta,c,ldc)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/lda,ldb,ldc,m,n/))
  ABI_UNUSED((/alpha,beta/))
  ABI_UNUSED((/side,uplo/))
  ABI_UNUSED((/a,b,c/))
#endif

end subroutine abi_gpu_zhemm_cptr
!!***

subroutine abi_gpu_zhemm_d(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 real(dp),intent(in),target :: a(*),b(*)
 real(dp),intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_zhemm_d
!!***

subroutine abi_gpu_zhemm_z(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 complex(dpc),intent(in),target :: a(*),b(*)
 complex(dpc),intent(inout),target :: c(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_zhemm_z
!!***

subroutine abi_gpu_zhemm_2d(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 real(dp),intent(in),target :: a(lda,*),b(ldb,*)
 real(dp),intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_zhemm_2d
!!***

subroutine abi_gpu_zhemm_2z(side,uplo,m,n,alpha,a,lda,b,ldb,beta,c,ldc)

!Arguments ------------------------------------
 integer,intent(in) :: lda,ldb,ldc,m,n
 complex(dpc),intent(in) :: alpha,beta
 character(len=1),intent(in) :: side,uplo
 complex(dpc),intent(in),target :: a(lda,*),b(ldb,*)
 complex(dpc),intent(inout),target :: c(ldc,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
        c_loc(a),lda,&
        c_loc(b),ldb,&
        beta,&
        c_loc(c),ldc)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b,c)
    call abi_gpu_zhemm_cptr(side,uplo,m,n,alpha,&
      c_loc(a),lda,&
      c_loc(b),ldb,&
      beta,&
      c_loc(c),ldc)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_zhemm_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xscal
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xscal
!! NAME
!!  abi_gpu_xscal
!!
!! FUNCTION
!!  Compute a BLAS-1 SCAL operation on GPU
!!  x = alpha * x
!!
!! INPUTS
!!  cplx   = 1 if real 2 if complex
!!  size   = vector size
!!  alpha  = scalar complex value
!!  x_gpu  = pointer to gpu memory location of array x
!!  incrx  = stride between consecutive elements of x
!!
!! SOURCE
subroutine abi_gpu_xscal_cptr(cplx, size, alpha, x, incrx)

 ! !Arguments ------------------------------------
 integer,      intent(in)    :: cplx
 integer,      intent(in)    :: size
 complex(dpc), intent(in)    :: alpha
 type(c_ptr),  intent(in)    :: x
 integer,      intent(in)    :: incrx

! *************************************************************************

 if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
   ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
 end if

#ifdef HAVE_GPU

 call gpu_xscal(cplx, size, alpha, x, incrx)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,incrx,size/))
  ABI_UNUSED((/alpha/))
  ABI_UNUSED((/x/))
#endif

end subroutine abi_gpu_xscal_cptr
!!***

subroutine abi_gpu_xscal_d(cplx, size, alpha, x, incrx)

 ! !Arguments ------------------------------------
 integer,      intent(in)    :: cplx
 integer,      intent(in)    :: size
 complex(dpc), intent(in)    :: alpha
 real(dp),     intent(inout), target :: x(*)
 integer,      intent(in)    :: incrx

! *************************************************************************

 if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
   ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
 end if

 if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
 else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_PTR(x)
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
   !$OMP END TARGET DATA
#endif
 else
   ABI_BUG("Unhandled GPU mode !")
 end if

end subroutine abi_gpu_xscal_d
!!***

subroutine abi_gpu_xscal_z(cplx, size, alpha, x, incrx)

 ! !Arguments ------------------------------------
 integer,      intent(in)    :: cplx
 integer,      intent(in)    :: size
 complex(dpc), intent(in)    :: alpha
 complex(dpc), intent(inout), target :: x(*)
 integer,      intent(in)    :: incrx

! *************************************************************************

 if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
   ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
 end if

 if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
 else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_PTR(x)
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
   !$OMP END TARGET DATA
#endif
 else
   ABI_BUG("Unhandled GPU mode !")
 end if

end subroutine abi_gpu_xscal_z
!!***

subroutine abi_gpu_xscal_2d(cplx, size, alpha, x, incrx)

 ! !Arguments ------------------------------------
 integer,      intent(in)    :: cplx
 integer,      intent(in)    :: size
 complex(dpc), intent(in)    :: alpha
 real(dp),     intent(inout), target :: x(size,*)
 integer,      intent(in)    :: incrx

! *************************************************************************

 if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
   ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
 end if

 if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
 else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_PTR(x)
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
   !$OMP END TARGET DATA
#endif
 else
   ABI_BUG("Unhandled GPU mode !")
 end if

end subroutine abi_gpu_xscal_2d
!!***

subroutine abi_gpu_xscal_2z(cplx, size, alpha, x, incrx)

 ! !Arguments ------------------------------------
 integer,      intent(in)    :: cplx
 integer,      intent(in)    :: size
 complex(dpc), intent(in)    :: alpha
 complex(dpc), intent(inout), target :: x(size,*)
 integer,      intent(in)    :: incrx

! *************************************************************************

 if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
   ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
 end if

 if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
 else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
   !$OMP TARGET DATA USE_DEVICE_PTR(x)
   call abi_gpu_xscal_cptr(cplx, size, alpha, c_loc(x), incrx)
   !$OMP END TARGET DATA
#endif
 else
   ABI_BUG("Unhandled GPU mode !")
 end if

end subroutine abi_gpu_xscal_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xaxpy
!------------------------------------------------------------------------------

!!****f* m_abi_gpu_linalg/abi_gpu_xaxpy
!! NAME
!!  abi_gpu_xaxpy
!!
!! FUNCTION
!!  Compute a BLAS-1 AXPY operation on GPU
!!  y = alpha * x + y
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  size  = vector size
!!  alpha = scalar complex value
!!  x_gpu = pointer to gpu memory location of array x
!! incrx  = stride between consecutive elements of x
!!  y_gpu = pointer to gpu memory location of array y
!! incry  = stride between consecutive elements of y
!!
!! SOURCE
subroutine abi_gpu_xaxpy_cptr(cplx, size, alpha, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  type(c_ptr),  intent(in)    :: x
  integer,      intent(in)    :: incrx
  type(c_ptr),  intent(in)    :: y
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_xaxpy(cplx, size, alpha, x, incrx, y, incry)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,incrx,incry,size/))
  ABI_UNUSED((/alpha/))
  ABI_UNUSED((/x,y/))
#endif

end subroutine abi_gpu_xaxpy_cptr
!!***

subroutine abi_gpu_xaxpy_d(cplx, size, alpha, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  real(dp),     intent(in), target    :: x(*)
  integer,      intent(in)    :: incrx
  real(dp),     intent(inout), target :: y(*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xaxpy_d
!!***

subroutine abi_gpu_xaxpy_z(cplx, size, alpha, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  complex(dpc), intent(in), target    :: x(*)
  integer,      intent(in)    :: incrx
  complex(dpc), intent(inout), target :: y(*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xaxpy_z
!!***

subroutine abi_gpu_xaxpy_2d(cplx, size, alpha, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  real(dp),     intent(in), target    :: x(size,*)
  integer,      intent(in)    :: incrx
  real(dp),     intent(inout), target :: y(size,*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xaxpy_2d
!!***

subroutine abi_gpu_xaxpy_2z(cplx, size, alpha, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in)    :: alpha
  complex(dpc), intent(in), target    :: x(size,*)
  integer,      intent(in)    :: incrx
  complex(dpc), intent(inout), target :: y(size,*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xaxpy_cptr(cplx, size, alpha, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xaxpy_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xcopy
!------------------------------------------------------------------------------

!!****f* m_abi_gpu_linalg/abi_gpu_xcopy
!! NAME
!!  abi_gpu_xcopy
!!
!! FUNCTION
!!  Compute a BLAS-1 COPY operation on GPU
!!  y = x (copy x into y)
!!
!! INPUTS
!!  cplx   = 1 if real 2 if complex
!!  size   = input vector size
!!  x      = pointer to gpu memory location of array x
!!  incrx  = stride between consecutive elements of x
!!  y      = pointer to gpu memory location of array y
!!  incry  = stride between consecutive elements of y
!!
!! SOURCE
subroutine abi_gpu_xcopy_cptr(cplx, size, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  type(c_ptr),  intent(in)    :: x
  integer,      intent(in)    :: incrx
  type(c_ptr),  intent(in)    :: y
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_xcopy(cplx, size, x, incrx, y, incry)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,incrx,incry,size/))
  ABI_UNUSED((/x,y/))
#endif

end subroutine abi_gpu_xcopy_cptr
!!***

subroutine abi_gpu_xcopy_d(cplx, size, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  real(dp),     intent(in),target    :: x(*)
  integer,      intent(in)    :: incrx
  real(dp),     intent(inout),target :: y(*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xcopy_d
!!***

subroutine abi_gpu_xcopy_z(cplx, size, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in),target    :: x(*)
  integer,      intent(in)    :: incrx
  complex(dpc), intent(inout),target :: y(*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xcopy_z
!!***

subroutine abi_gpu_xcopy_2d(cplx, size, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  real(dp),     intent(in),target    :: x(size,*)
  integer,      intent(in)    :: incrx
  real(dp),     intent(inout),target :: y(size,*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xcopy_2d
!!***

subroutine abi_gpu_xcopy_2z(cplx, size, x, incrx, y, incry)

  ! !Arguments ------------------------------------
  integer,      intent(in)    :: cplx
  integer,      intent(in)    :: size
  complex(dpc), intent(in),target    :: x(size,*)
  integer,      intent(in)    :: incrx
  complex(dpc), intent(inout),target :: y(size,*)
  integer,      intent(in)    :: incry

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(x,y)
    call abi_gpu_xcopy_cptr(cplx, size, c_loc(x), incrx, c_loc(y), incry)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xcopy_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xtrsm
!------------------------------------------------------------------------------

!!****f* m_abi_gpu_linalg/abi_gpu_xtrsm
!! NAME
!!  abi_gpu_xtrsm
!!
!! FUNCTION
!! Solves a matrix equation (one matrix operand is triangular) on GPU.
!! The xtrsm routines solve one of the following matrix equations
!! op(a)*x = alpha*b
!! or
!! x*op(a) = alpha*b,
!!
!! INPUTS
!! cplx= 1 if real 2 if complex
!! side= Specifies whether op(a) appears on the left or right of x for
!!      the operation to be performed as follows:
!!      L or l op(a)*x = alpha*b
!!      R or r x*op(a) = alpha*b
!! uplo= Specifies whether the matrix a is an upper or lower triangular
!!      matrix as follows:
!!      U or u Matrix a is an upper triangular matrix.
!!      L or l Matrix a is a lower triangular matrix
!! transa= Specifies the form of op(a) to be used in the matrix
!!      multiplication as follows:
!!      N or n op(a) = a
!!      T or t op(a) = a'
!!      C or c op(a) = conjg(a')
!! diag= Specifies whether or not a is unit triangular as follows:
!!      U or u Matrix a is assumed to be unit triangular.
!!      N or n Matrix a is not assumed to be unit triangular.
!! m= Specifies the number of rows of b. The value of m must be at least zero
!! n= Specifies the number of columns of b. The value of n must be at least zero
!! alpha= Specifies the scalar alpha. When alpha is zero, then a is not referenced and b
!!      need not be set before entry.
!!  a = pointer to gpu memory location of  array a, DIMENSION (lda, k), where k is m when side = 'L' or 'l' and is n
!!      when side = 'R' or 'r'.
!! lda= Specifies the first dimension of a as declared in the calling
!!     (sub)program. When side = 'L' or 'l', then lda must be at least max(1,
!!      m), when side = 'R' or 'r', then lda must be at least max(1, n).
!!  b = pointer to gpu memory location of  b Array, DIMENSION (ldb,n). Before entry, the leading m-by-n part of the array
!!     b must contain the right-hand side matrix b.
!! ldb= Specifies the first dimension of b as declared in the calling
!!     (sub)program. The value of ldb must be at least max(1, m).
!!
!! OUTPUT
!!  b
!!
!! SIDE EFFECTS
!!   WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
!!   the correct one is in 17_toolbox/gpu_linalg.cu
!!
!! SOURCE
subroutine abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 type(c_ptr),intent(in) :: a
 type(c_ptr),intent(in) :: b

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  call gpu_xtrsm(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,lda,ldb,m,n/))
  ABI_UNUSED((/alpha/))
  ABI_UNUSED((/side,uplo,transa,diag/))
  ABI_UNUSED((/a,b/))
#endif

end subroutine abi_gpu_xtrsm_cptr
!!***

subroutine abi_gpu_xtrsm_d(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 real(dp),   intent(in),target :: a(*)
 real(dp),   intent(inout),target :: b(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b)
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xtrsm_d
!!***

subroutine abi_gpu_xtrsm_z(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 complex(dpc),intent(in),target :: a(*)
 complex(dpc),intent(inout),target :: b(*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b)
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xtrsm_z
!!***

subroutine abi_gpu_xtrsm_2d(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 real(dp),   intent(in),target :: a(lda,*)
 real(dp),   intent(inout),target :: b(ldb,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b)
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xtrsm_2d
!!***

subroutine abi_gpu_xtrsm_2z(cplx,side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)

! !Arguments ------------------------------------
 integer, intent(in) :: cplx,lda,ldb,m,n
 complex(dpc), intent(in) :: alpha
 character(len=1), intent(in) :: side,uplo,transa,diag
 complex(dpc),intent(in),target :: a(lda,*)
 complex(dpc),intent(inout),target :: b(ldb,*)

! *********************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(a,b)
    call abi_gpu_xtrsm_cptr(cplx,side,uplo,transa,&
      diag,m,n,alpha,c_loc(a),lda,c_loc(b),ldb)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xtrsm_2z
!!***


!!****f* m_abi_gpu_linalg/abi_gpu_work_resizeI
!!
!! NAME
!! abi_gpu_work_resizeI
subroutine abi_gpu_work_resizeI(array,array_managed,current_dim,asked_dim)

  integer, allocatable, intent(inout) :: array(:)
  integer(kind=c_int32_t), ABI_CONTIGUOUS pointer, intent(inout) :: array_managed(:)
  integer, intent(inout)  :: current_dim
  integer, intent(in   )  :: asked_dim

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then

#ifdef HAVE_YAKL
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array_managed) ) then
        ABI_FREE_MANAGED(array_managed)
      end if
      ABI_MALLOC_MANAGED(array_managed,(/asked_dim/))
    end if
#endif

  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then

#ifdef HAVE_OPENMP_OFFLOAD
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        !$OMP TARGET EXIT DATA MAP(release:array)
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
      !$OMP TARGET ENTER DATA MAP(alloc:array)
    end if
#endif

  else
    ABI_BUG("Unhandled GPU mode !")
  end if

#ifndef HAVE_GPU
  ! Unused if GPU code disabled
  ABI_UNUSED(array)
  ABI_UNUSED(array_managed)
  ABI_UNUSED((/current_dim,asked_dim/))
#endif

end subroutine abi_gpu_work_resizeI
!!***

!!****f* m_abi_gpu_linalg/abi_gpu_work_resizeR
!!
!! NAME
!! abi_gpu_work_resizeR

subroutine abi_gpu_work_resizeR(array,array_managed,current_dim,asked_dim)

  double precision, allocatable, intent(inout) :: array(:)
  real(kind=c_double), ABI_CONTIGUOUS pointer, intent(inout) :: array_managed(:)
  integer, intent(inout) :: current_dim
  integer, intent(in   ) :: asked_dim

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then

#ifdef HAVE_YAKL
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array_managed) ) then
        ABI_FREE_MANAGED(array_managed)
      end if
      ABI_MALLOC_MANAGED(array_managed,(/asked_dim/))
    end if
#endif

  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then

#ifdef HAVE_OPENMP_OFFLOAD
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        !$OMP TARGET EXIT DATA MAP(release:array)
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
      !$OMP TARGET ENTER DATA MAP(alloc:array)
    end if
#endif

  else
    ABI_BUG("Unhandled GPU mode !")
  end if

#ifndef HAVE_GPU
  ! Unused if GPU code disabled
  ABI_UNUSED(array)
  ABI_UNUSED(array_managed)
  ABI_UNUSED((/current_dim,asked_dim/))
#endif

end subroutine abi_gpu_work_resizeR
!!***

!!****f* m_abi_gpu_linalg/abi_gpu_work_resizeC
!!
!! NAME
!! abi_gpu_work_resizeC

subroutine abi_gpu_work_resizeC(array,array_managed,current_dim,asked_dim)

  complex(kind=8), allocatable, intent(inout) :: array(:)
  complex(kind=c_double_complex), ABI_CONTIGUOUS pointer, intent(inout) :: array_managed(:)
  integer, intent(inout)  :: current_dim
  integer, intent(in   )  :: asked_dim

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then

#ifdef HAVE_YAKL
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( associated(array_managed) ) then
        ABI_FREE_MANAGED(array_managed)
      end if
      ABI_MALLOC_MANAGED(array_managed,(/asked_dim/))
    end if
#endif

  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then

#ifdef HAVE_OPENMP_OFFLOAD
    if ( current_dim < asked_dim  ) then
      current_dim = asked_dim
      if ( allocated(array) ) then
        !$OMP TARGET EXIT DATA MAP(release:array)
        ABI_FREE(array)
      end if
      ABI_MALLOC(array,(asked_dim))
      !$OMP TARGET ENTER DATA MAP(alloc:array)
    end if
#endif

  else
    ABI_BUG("Unhandled GPU mode !")
  end if

#ifndef HAVE_GPU
  ! Unused if GPU code disabled
  ABI_UNUSED(array)
  ABI_UNUSED(array_managed)
  ABI_UNUSED((/current_dim,asked_dim/))
#endif

end subroutine abi_gpu_work_resizeC
!!***

!!****f* m_abi_gpu_linalg/abi_gpu_work_resizeCptr
!!
!! NAME
!! abi_gpu_work_resizeCptr

subroutine abi_gpu_work_resizeCptr(array,current_dim,asked_dim)

  type(c_ptr), intent(inout) :: array
  integer(c_size_t), intent(inout)  :: current_dim
  integer(c_size_t), intent(in   )  :: asked_dim

! *************************************************************************

  if ( current_dim < asked_dim  ) then
    if(current_dim == 0) then
      call dealloc_on_gpu(array)
    end if
    current_dim = asked_dim
    call alloc_on_gpu(array, asked_dim)
  end if

#ifndef HAVE_GPU
  ! Unused if GPU code disabled
  ABI_UNUSED(array)
  ABI_UNUSED((/current_dim,asked_dim/))
#endif

end subroutine abi_gpu_work_resizeCptr
!!***

subroutine abi_gpu_work_finalize()

#ifdef HAVE_GPU
  !FIXME Assuming managed here ?
  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then

#ifdef HAVE_YAKL
    if ( associated(i_work_managed) ) then
      ABI_FREE_MANAGED(i_work_managed)
    end if
    if ( associated(r_work_managed) ) then
      ABI_FREE_MANAGED(r_work_managed)
    end if
    if ( associated(c_work_managed) ) then
      ABI_FREE_MANAGED(c_work_managed)
    end if
#endif

  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then

#ifdef HAVE_OPENMP_OFFLOAD
    if ( allocated(i_work) ) then
      !$OMP TARGET EXIT DATA MAP(release:i_work)
      ABI_FREE(i_work)
    end if
    if ( allocated(r_work) ) then
      !$OMP TARGET EXIT DATA MAP(release:r_work)
      ABI_FREE(r_work)
    end if
    if ( allocated(c_work) ) then
      !$OMP TARGET EXIT DATA MAP(release:c_work)
      ABI_FREE(c_work)
    end if
#endif

    if(gpu_work_len > 0) then
      call dealloc_on_gpu(gpu_work)
    end if

  end if

  i_work_len = 0
  r_work_len = 0
  c_work_len = 0
  gpu_work_len = 0

#endif

end subroutine abi_gpu_work_finalize
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xhegvd
!------------------------------------------------------------------------------

!!****f* m_abi_gpu_linalg/abi_gpu_xhegvd
!! NAME
!!  abi_gpu_xhegvd
!!
!! FUNCTION
!!  Compute a LAPACK SYGVD operation on GPU
!!  compute eigen values/vectors of a real generalized
!!  symmetric-definite eigenproblem
!!
!! See cusolver documentation
!! https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-hegvd
!!
!! See also LAPACK doc in reference implementation:
!! https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dhegvd.f
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  itype = integer, type of problem
!!  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
!!  uplo  = character, 'u' or 'l'
!!  A_nrows = matrix size
!!  A = pointer to gpu memory location of matrix A
!!  lda   = leading dimension of matrix A
!!  B = pointer to gpu memory location of matrix B
!!  ldb   = leading dimension of matrix B
!!  W = pointer to gpu memory location of matrix W (output eigen values)
!!  devInfo  =
!!
!! SOURCE
subroutine abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, A_nrows, &
                   A, lda, &
                   B, ldb, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  type(c_ptr),     intent(in   ) :: A
  integer,         intent(in   ) :: lda
  type(c_ptr),     intent(in   ) :: B
  integer,         intent(in   ) :: ldb
  type(c_ptr),     intent(in   ) :: W
  integer,         intent(inout) :: devInfo

  ! Local variables ------------------------------
  integer     :: bufferSize
  type(c_ptr) :: gpu_ptr

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  ! probe needed bufferSize
  call gpu_xsygvd_buffersize(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 A, lda, &
                 B, ldb, &
                 W, &
                 bufferSize)

  select case(cplx)

  case (1)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = c_loc(r_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(r_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(1,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  case (2)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = c_loc(c_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(c_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(2,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  end select

  ! and compute (finally)
  call gpu_xsygvd(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 A, lda, &
                 B, ldb, &
                 W, &
                 gpu_ptr, bufferSize, devInfo)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,itype,A_nrows,lda,ldb,devInfo,bufferSize/))
  ABI_UNUSED((/jobz,uplo/))
  ABI_UNUSED((/A,B,W,gpu_ptr/))
#endif

end subroutine abi_gpu_xhegvd_cptr
!!***

subroutine abi_gpu_xhegvd_d(cplx, itype, jobz, uplo, A_nrows, &
                   A, lda, &
                   B, ldb, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  real(dp),        intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  real(dp),        intent(in   ),target :: B(*)
  integer,         intent(in   ) :: ldb
  real(dp),        intent(inout),target :: W(*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W)
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xhegvd_d
!!***

subroutine abi_gpu_xhegvd_z(cplx, itype, jobz, uplo, A_nrows, &
                   A, lda, &
                   B, ldb, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  complex(dpc),    intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  complex(dpc),    intent(in   ),target :: B(*)
  integer,         intent(in   ) :: ldb
  real(dp),        intent(inout),target :: W(*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W)
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xhegvd_z
!!***

subroutine abi_gpu_xhegvd_2d(cplx, itype, jobz, uplo, A_nrows, &
                   A, lda, &
                   B, ldb, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda,ldb
  real(dp),        intent(in   ),target :: A(lda,*)
  real(dp),        intent(in   ),target :: B(ldb,*)
  real(dp),        intent(inout),target :: W(A_nrows,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W)
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xhegvd_2d
!!***

subroutine abi_gpu_xhegvd_2z(cplx, itype, jobz, uplo, A_nrows, &
                   A, lda, &
                   B, ldb, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  integer,         intent(in   ) :: itype
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda,ldb
  complex(dpc),    intent(in   ),target :: A(lda,*)
  complex(dpc),    intent(in   ),target :: B(ldb,*)
  real(dp),        intent(inout),target :: W(A_nrows,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,B,W)
    call abi_gpu_xhegvd_cptr(cplx, itype, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(B), ldb, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xhegvd_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xheevd
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xheevd
!! NAME
!!  abi_gpu_xheevd
!!
!! FUNCTION
!!  Compute a LAPACK SYEVD operation on GPU
!!  compute eigen values/vectors of a real
!!  symmetric-definite eigenproblem
!!
!! See cusolver documentation
!! https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-heevd
!!
!! See also LAPACK doc in reference implementation:
!! https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dheevd.f
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
!!  uplo  = character, 'u' or 'l'
!!  A_nrows = matrix size
!!  A = pointer to gpu memory location of matrix A
!!  lda   = leading dimension of matrix A
!!  W = pointer to gpu memory location of matrix W (output eigen values)
!!  devInfo  =
!!
!! SOURCE
subroutine abi_gpu_xheevd_cptr(cplx, jobz, uplo, A_nrows, &
                   A, lda, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  type(c_ptr),     intent(in   ) :: A
  integer,         intent(in   ) :: lda
  type(c_ptr),     intent(in) :: W
  integer,         intent(inout) :: devInfo

  ! Local variables ------------------------------
  integer     :: bufferSize
  type(c_ptr) :: gpu_ptr

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  ! probe needed bufferSize
  call gpu_xsyevd_buffersize(cplx, jobz, uplo, &
                 A_nrows, &
                 A, lda, &
                 W, &
                 bufferSize)

  select case(cplx)

  case (1)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = c_loc(r_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(r_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(1,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  case (2)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = c_loc(c_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(c_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(2,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  end select

  ! and compute (finally)
  call gpu_xsyevd(cplx, jobz, uplo, &
                 A_nrows, &
                 A, lda, &
                 W, &
                 gpu_ptr, bufferSize, devInfo)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,A_nrows,lda,devInfo,bufferSize/))
  ABI_UNUSED((/jobz,uplo/))
  ABI_UNUSED((/A,W,gpu_ptr/))
#endif

end subroutine abi_gpu_xheevd_cptr
!!***

subroutine abi_gpu_xheevd_d(cplx, jobz, uplo, A_nrows, &
                   A, lda, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  real(dp),        intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  real(dp),        intent(inout),target :: W(*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,W)
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xheevd_d
!!***

subroutine abi_gpu_xheevd_z(cplx, jobz, uplo, A_nrows, &
                   A, lda, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  complex(dpc),    intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  real(dp),        intent(inout),target :: W(*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,W)
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xheevd_z
!!***

subroutine abi_gpu_xheevd_2d(cplx, jobz, uplo, A_nrows, &
                   A, lda, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda
  real(dp),        intent(in   ),target :: A(lda,*)
  real(dp),        intent(inout),target :: W(A_nrows,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,W)
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xheevd_2d
!!***

subroutine abi_gpu_xheevd_2z(cplx, jobz, uplo, A_nrows, &
                   A, lda, &
                   W,      &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: jobz
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda
  complex(dpc),    intent(in   ),target :: A(lda,*)
  real(dp),        intent(inout),target :: W(A_nrows,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A,W)
    call abi_gpu_xheevd_cptr(cplx, jobz, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 c_loc(W), &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xheevd_2z
!!***

!------------------------------------------------------------------------------
!                         abi_gpu_xpotrf
!------------------------------------------------------------------------------
!!****f* m_abi_gpu_linalg/abi_gpu_xpotrf
!! NAME
!!  abi_gpu_xpotrf
!!
!! FUNCTION
!!  Compute a LAPACK SYGVD operation on GPU
!!  compute eigen values/vectors of a real
!!  symmetric-definite eigenproblem
!!
!! See cusolver documentation
!! https://docs.nvidia.com/cuda/cusolver/index.html#cuSolverDN-lt-t-gt-potrf
!!
!! See also LAPACK doc in reference implementation:
!! https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dpotrf.f
!!
!! INPUTS
!!  cplx  = 1 if real 2 if complex
!!  uplo  = character, 'u' or 'l'
!!  A_nrows = matrix size
!!  A = pointer to gpu memory location of matrix A
!!  lda   = leading dimension of matrix A
!!  devInfo  =
!!
!! SOURCE
subroutine abi_gpu_xpotrf_cptr(cplx, uplo, A_nrows, &
                   A, lda, &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  type(c_ptr),     intent(in   ) :: A
  integer,         intent(in   ) :: lda
  integer,         intent(inout) :: devInfo

  ! Local variables ------------------------------
  integer     :: bufferSize
  type(c_ptr) :: gpu_ptr

! *************************************************************************

  if (abi_linalg_gpu_mode == ABI_GPU_DISABLED) then
    ABI_BUG("You requested to run on CPU to a GPU wrapper :/")
  end if

#ifdef HAVE_GPU

  ! probe needed bufferSize
  call gpu_xpotrf_buffersize(cplx, uplo, &
                 A_nrows, &
                 A, lda, &
                 bufferSize)

  select case(cplx)

  case (1)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = c_loc(r_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(r_work,r_work_managed,r_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(r_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(1,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  case (2)
    ! resize work array if needed and retrieve work pointer to use
    if(abi_linalg_gpu_mode == ABI_GPU_LEGACY &
          .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = c_loc(c_work_managed)
    else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_GET_MAPPED_PTR
      call abi_gpu_work_resize(c_work,c_work_managed,c_work_len,bufferSize)
      gpu_ptr = xomp_get_mapped_ptr(c_loc(c_work))
#else
      call abi_gpu_work_resizeCptr(gpu_work,gpu_work_len,INT(1,c_size_t)*bufferSize*dp)
      gpu_ptr = gpu_work
#endif
    end if

  end select

  ! and compute (finally)
  call gpu_xpotrf(cplx, uplo, &
                 A_nrows, &
                 A, lda, &
                 gpu_ptr, bufferSize, devInfo)

  if (abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
    ! CUDA/HIP linalg calls are run asynchronously and OpenMP is unaware of them.
    ! Therefore, we issue a stream sync here to avoid
    !potential mistakes in calling context.
    call gpu_linalg_stream_synchronize()
  end if

#else
  ! Unused if GPU code disabled
  ABI_UNUSED((/cplx,A_nrows,lda,devInfo,bufferSize/))
  ABI_UNUSED((/uplo/))
  ABI_UNUSED((/A,gpu_ptr/))
#endif

end subroutine abi_gpu_xpotrf_cptr
!!***

subroutine abi_gpu_xpotrf_d(cplx, uplo, A_nrows, &
                   A, lda, &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  real(dp),        intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A)
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xpotrf_d
!!***

subroutine abi_gpu_xpotrf_z(cplx, uplo, A_nrows, &
                   A, lda, &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows
  complex(dpc),    intent(in   ),target :: A(*)
  integer,         intent(in   ) :: lda
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A)
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xpotrf_z
!!***

subroutine abi_gpu_xpotrf_2d(cplx, uplo, A_nrows, &
                   A, lda, &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda
  real(dp),        intent(in   ),target :: A(lda,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A)
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xpotrf_2d
!!***

subroutine abi_gpu_xpotrf_2z(cplx, uplo, A_nrows, &
                   A, lda, &
                   devInfo)

  ! Arguments ------------------------------------
  integer,         intent(in   ) :: cplx
  character(len=1),intent(in   ) :: uplo
  integer,         intent(in   ) :: A_nrows,lda
  complex(dpc),    intent(in   ),target :: A(lda,*)
  integer,         intent(inout) :: devInfo

! *************************************************************************

  if(abi_linalg_gpu_mode == ABI_GPU_LEGACY .or. abi_linalg_gpu_mode == ABI_GPU_KOKKOS) then
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
  else if(abi_linalg_gpu_mode == ABI_GPU_OPENMP) then
#ifdef HAVE_OPENMP_OFFLOAD
    !$OMP TARGET DATA USE_DEVICE_PTR(A)
    call abi_gpu_xpotrf_cptr(cplx, uplo, &
                 A_nrows, &
                 c_loc(A), lda, &
                 devInfo)
    !$OMP END TARGET DATA
#endif
  else
    ABI_BUG("Unhandled GPU mode !")
  end if

end subroutine abi_gpu_xpotrf_2z
!!***

!!****f* m_abi_linalg/gpu_xorthonormalize
!! NAME
!!  gpu_xorthonormalize
!!
!! FUNCTION
!!  This routine computes the overlap of two complex wavefunctions (for a given number of bands)
!!  and orthonormalizes it using gpu:
!!      - Computes the products of two rectangular matrices
!!         containing the wavefunctions psi and S.psi (where S is the
!!         overlap (with the PAW terms if necessary)).
!!      - Does a Cholesky decomposition of this overlap
!!      - rotates the initial matrix blockvectorx by the triangular matrix to
!!         have an orthonormal set of wavefunctions
!!
!! INPUTS
!!  blockvectorbx = matrix of dimension (blocksize,vectsize) as a GPU ptr
!!                  (e.g. block of overlap*wavefunction)
!!  blocksize     = dimension of matrices (e.g number of bands)
!!  spaceComm     = communicator used for  MPI parallelization
!!  vectsize      = dimension of matrices (e.g number of G vector)
!!
!! OUTPUT
!!  sqgram        = Choleski decomposition of transpose(blockvector)*blockvectorx as a GPU ptr
!!
!! SIDE EFFECTS
!!  blockvectorx  = on input, matrix of dimension (vectsize,blocksize) as a GPU ptr
!!                  (e.g block of wavefunction)
!!  blockvectorx  = on output, orthonormalized wavefunction. as a GPU ptr
!!
!!
!! SOURCE

subroutine gpu_xorthonormalize(blockvectorx_gpu,blockvectorbx_gpu,blocksize,spaceComm,&
&                              sqgram_gpu,vectsize,&
&                              x_cplx,timopt,tim_xortho) ! optional arguments

  use, intrinsic :: iso_c_binding
  implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: blocksize,spaceComm,vectsize,x_cplx
 integer, intent(in), optional :: timopt,tim_xortho
!arrays
 type(c_ptr),intent(inout) :: blockvectorbx_gpu, blockvectorx_gpu, sqgram_gpu
!Local variables-------------------------------
#if defined HAVE_GPU
 integer :: ierr,info
 real(dp),    dimension(:,:),allocatable, target :: d_sqgram
 complex(dpc),dimension(:,:),allocatable, target :: z_sqgram
 character :: tr
 real(dp) :: tsec(2)
 integer(c_size_t) :: size
#else
 type(c_ptr) :: cptr_a
#endif
 character(len=500) :: message

! *********************************************************************

#if defined HAVE_GPU
 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xortho,1,tsec)
   end if
 end if

 if ( x_cplx == 1 ) then
   tr='t'
   ABI_MALLOC(d_sqgram,(blocksize,blocksize))
 else
   tr='c'
   ABI_MALLOC(z_sqgram,(blocksize,blocksize))
 end if

 call gpu_xgemm(x_cplx,tr,'n',blocksize,blocksize,vectsize, &
& cone,blockvectorx_gpu,vectsize,blockvectorbx_gpu,vectsize,czero,sqgram_gpu,blocksize)
 call gpu_device_synchronize()
   size=x_cplx*dp*blocksize*blocksize

 if ( x_cplx == 1 ) then
   call copy_from_gpu(d_sqgram, sqgram_gpu, INT(x_cplx, c_size_t)*dp*blocksize*blocksize)
   call xmpi_sum(d_sqgram,spaceComm,ierr)
   call abi_xpotrf('u',blocksize,d_sqgram,blocksize,info)
   call copy_on_gpu(d_sqgram, sqgram_gpu, INT(x_cplx, c_size_t)*dp*blocksize*blocksize)
 else
   call copy_from_gpu(z_sqgram, sqgram_gpu, size)
   call xmpi_sum(z_sqgram,spaceComm,ierr)
   call abi_xpotrf('u',blocksize,z_sqgram,blocksize,info)
   call copy_on_gpu(z_sqgram, sqgram_gpu, size)
 end if

 if (info /= 0 ) then
   write(message,'(a,i3)') '  xpotrf, info=',info
   ABI_WARNING(message)
 end if

 call gpu_xtrsm(x_cplx,'r','u','n','n',vectsize,blocksize,cone,sqgram_gpu,blocksize,&
&               blockvectorx_gpu,vectsize)

 if(x_cplx==1) then
   ABI_FREE(d_sqgram)
 else
   ABI_FREE(z_sqgram)
 end if
 if (present(tim_xortho).and.present(timopt)) then
   if(abs(timopt)==3) then
     call timab(tim_xortho,2,tsec)
   end if
 end if
 return

#else
 message='  This routine is not allowed when running on GPU is disabled !'
 ABI_BUG(message)
 if (.false.) then
   write(std_out,*) blocksize,vectsize,spaceComm,x_cplx
   if (present(timopt))  write(std_out,*) timopt
   if (present(tim_xortho))  write(std_out,*) tim_xortho
   cptr_a=blockvectorbx_gpu;cptr_a=blockvectorx_gpu;cptr_a=sqgram_gpu
 end if
#endif

end subroutine gpu_xorthonormalize
!!***

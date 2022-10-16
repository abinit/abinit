/* gpu_linalg.cu */

/*
 * Copyright (C) 2008-2022 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain cublas and magma encapsulation routines,
 * that will be callable from fortran routines
 *
 */

#include <gpu_linalg.h>

cublasHandle_t cublas_handle;
cusolverDnHandle_t cusolverDn_handle;

// utility function for compatiblity between cublas v1/v2 API
cublasOperation_t select_cublas_op(char *c)
{
  cublasOperation_t op;

  if (*c == 'n' or *c == 'N')
    op = CUBLAS_OP_N;
  else if (*c == 't' or *c == 'T')
    op = CUBLAS_OP_T;
  else if (*c == 'c' or *c == 'C')
    op = CUBLAS_OP_C;
  else
    printf("CUBLAS API error, character can't be converted to a valid cublas operation !\n");

  return op;
}

// utility function for compatiblity between cublas v1/v2 API
cublasSideMode_t select_cublas_side(char *c)
{
  cublasSideMode_t mode;

  if (*c == 'l' or *c == 'L')
    mode = CUBLAS_SIDE_LEFT;
  else if (*c == 'r' or *c == 'R')
    mode = CUBLAS_SIDE_RIGHT;
  else
    printf("CUBLAS API error, character can't be converted to a valid cublas side mode !\n");

  return mode;
}

// utility function for compatiblity between cublas v1/v2 API
cublasFillMode_t select_cublas_fill_mode(char *c)
{
  cublasFillMode_t mode;

  if (*c == 'u' or *c == 'U')
    mode = CUBLAS_FILL_MODE_UPPER;
  else if (*c == 'l' or *c == 'L')
    mode = CUBLAS_FILL_MODE_LOWER;
  else
    printf("CUBLAS API error, character can't be converted to a valid cublas fill mode !\n");

  return mode;
}

// utility function for compatiblity between cublas v1/v2 API
cublasDiagType_t select_cublas_diag_type(char *c)
{
  cublasDiagType_t diag;

  if (*c == 'n' or *c == 'N')
    diag = CUBLAS_DIAG_NON_UNIT;
  else if (*c == 'u' or *c == 'U')
    diag = CUBLAS_DIAG_UNIT;
  else
    printf("CUBLAS API error, character can't be converted to a valid cublas diag type !\n");

  return diag;
}


/*=========================================================================*/
// NAME
//  gpu_linalg_init
//
// FUNCTION
//  Initialisation of linalg environnement on GPU
//
//  WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
//             the correct one is in xx_gpu_toolbox/gpu_linalg.cu
/*=========================================================================*/

extern "C" void gpu_linalg_init_()
{
  CUDA_API_CHECK( cublasCreate(&cublas_handle) );
  CUDA_API_CHECK( cusolverDnCreate(&cusolverDn_handle) );
}

/*=========================================================================*/
// NAME
//  gpu_linalg_shutdown
//
// FUNCTION
//  Close linalg environnement on GPU
//
//  WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
//             the correct one is in xx_gpu_toolbox/gpu_linalg.cu
/*=========================================================================*/

extern "C" void gpu_linalg_shutdown_()
{
  CUDA_API_CHECK( cublasDestroy(cublas_handle) );
  CUDA_API_CHECK( cusolverDnDestroy(cusolverDn_handle) );
}

/*=========================================================================*/
// NAME
//  gpu_xgemm
//
// FUNCTION
//  Compute a scalar-matrix-matrix product and return a scalar-matrix product on GPU
//  c = alpha * op(a) * op(b) + beta * c
//
// INPUTS
//  cplx  = 1 if real 2 if complex
//  transa= from of op(a) to be used in the matrix multiplication
//  transb= from of op(b) to be used in the matrix multiplication
//  m     = number of rows of the matrix op(a) and of the matrix c
//  n     = number of rows of the matrix op(b) and the number of columns of the matrix c
//  k     = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
//  alpha = alpha scalar coefficient for matrix op(a)
//  a_gpu = pointer to gpu memory location of  matrix a
//  lda   = first dimension of a
//  b_gpu = pointer to gpu memory location of  matrix b
//  ldb   = first dimension of b
//  beta  = beta scalar coefficient for matrix c
//  c_gpu = pointer to gpu memory location of  matrix c
//  ldc   = first dimension of c
//
// OUTPUT
//  c_gpu     = c_gpu matrix
//
// WARNING! : this routine is a dummy one when HAVE_GPU_CUDA is not enabled
//            the correct one is in xx_gpu_toolbox/gpu_linalg.cu
/*=========================================================================*/

extern "C" void gpu_xgemm_(int* cplx, char *transA, char *transB, int *N, int *M, int *K,
                           cuDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           cuDoubleComplex *beta, void** C_ptr, int *ldc)
{

  cublasOperation_t opA = select_cublas_op(transA);
  cublasOperation_t opB = select_cublas_op(transB);

  (*cplx==1)?
    cublasDgemm(cublas_handle, opA, opB, *N, *M, *K, &((*alpha).x),
                (double *)(*A_ptr), *lda,
                (double *)(*B_ptr), *ldb,
                &((*beta).x),
                (double *)(*C_ptr), *ldc) :
    cublasZgemm(cublas_handle, opA, opB, *N, *M, *K, alpha,
                (cuDoubleComplex *)(*A_ptr), *lda,
                (cuDoubleComplex *)(*B_ptr), *ldb,
                beta,
                (cuDoubleComplex *)(*C_ptr), *ldc);

} // gpu_xgemm

/*=========================================================================*/
// NAME
//  gpu_xtrsm
//
// FUNCTION
// Solves a matrix equation (one matrix operand is triangular) on GPU.
// The xtrsm routines solve one of the following matrix equations
// op(a)*x = alpha*b
// or
// x*op(a) = alpha*b,
//
// INPUTS
// cplx= 1 if real 2 if complex
// side= Specifies whether op(a) appears on the left or right of x for
//      the operation to be performed as follows:
//      L or l op(a)*x = alpha*b
//      R or r x*op(a) = alpha*b
// uplo= Specifies whether the matrix a is an upper or lower triangular
//      matrix as follows:
//      U or u Matrix a is an upper triangular matrix.
//      L or l Matrix a is a lower triangular matrix
// transa= Specifies the form of op(a) to be used in the matrix
//      multiplication as follows:
//      N or n op(a) = a
//      T or t op(a) = a'
//      C or c op(a) = conjg(a')
// diag= Specifies whether or not a is unit triangular as follows:
//      U or u Matrix a is assumed to be unit triangular.
//      N or n Matrix a is not assumed to be unit triangular.
// m= Specifies the number of rows of b. The value of m must be at least zero
// n= Specifies the number of columns of b. The value of n must be at least zero
// alpha= Specifies the scalar alpha. When alpha is zero, then a is not referenced and b
//      need not be set before entry.
//  a_gpu = pointer to gpu memory location of  array a, DIMENSION (lda, k), where k is m when side = 'L' or 'l' and is n
//      when side = 'R' or 'r'.
// lda= Specifies the first dimension of a as declared in the calling
//     (sub)program. When side = 'L' or 'l', then lda must be at least max(1,
//      m), when side = 'R' or 'r', then lda must be at least max(1, n).
//  b_gpu = pointer to gpu memory location of  b Array, DIMENSION (ldb,n). Before entry, the leading m-by-n part of the array
//     b must contain the right-hand side matrix b.
// ldb= Specifies the first dimension of b as declared in the calling
//     (sub)program. The value of ldb must be at least max(1, m).
//
// OUTPUT
//  b_gpu
/*=========================================================================*/

extern "C" void gpu_xtrsm_(int* cplx, char *side, char *uplo, char *transA, char *diag,
                           int *N, int *M, cuDoubleComplex *alpha,
                           void **A_ptr, int *ldA,
                           void** B_ptr, int *ldB)
{

  cublasSideMode_t sideMode = select_cublas_side(side);
  cublasFillMode_t fillMode = select_cublas_fill_mode(uplo);
  cublasDiagType_t diagType = select_cublas_diag_type(diag);
  cublasOperation_t opA     = select_cublas_op(transA);

  (*cplx==1) ?
    cublasDtrsm(cublas_handle, sideMode, fillMode, opA, diagType,
                *N, *M, &((*alpha).x),
                (double *)(*A_ptr), *ldA,
                (double *)(*B_ptr), *ldB):
    cublasZtrsm(cublas_handle, sideMode, fillMode, opA, diagType,
                *N, *M, alpha,
                (cuDoubleComplex *)(*A_ptr), *ldA,
                (cuDoubleComplex *)(*B_ptr), *ldB);
} // gpu_xtrsm_

/*=========================================================================*/
// NAME
//  gpu_xaxpy
//
// FUNCTION
//  Compute blas-1 AXPY on GPU
//  y = alpha * x + y
//
// INPUTS
//  cplx  = 1 if real 2 if complex
//  alpha = scalar complex value (when doing computation with real value, only the real part is used)
//  x =     input vector
//  incrx = increment between two consecutive values of x
//  y =     in/out vector
//  incry = increment between two consecutive values of y
//
// OUTPUT
//  y
/*=========================================================================*/

extern "C" void gpu_xaxpy_(int* cplx, int *N,
                           cuDoubleComplex *alpha,
                           void **X_ptr, int *incrx, void** Y_ptr, int *incry)
{

  CUDA_API_CHECK( (*cplx==1) ?
                  cublasDaxpy(cublas_handle,*N,
                              &((*alpha).x),
                              (double *)(*X_ptr), *incrx,
                              (double *)(*Y_ptr), *incry) :
                  cublasZaxpy(cublas_handle, *N,
                              alpha,
                              (cuDoubleComplex *)(*X_ptr), *incrx,
                              (cuDoubleComplex *)(*Y_ptr), *incry) );
} // gpu_xaxpy_

/*=========================================================================*/
// NAME
//  gpu_xscal
//
// FUNCTION
//  Compute blas-1 SCAL on GPU
//  x = alpha * x
//
// INPUTS
//  cplx  = 1 if real 2 if complex
//  alpha = scalar complex value (when doing computation with real value, only the real part is used)
//  x =     input/output vector
//  incrx = increment between two consecutive values of x
//
// OUTPUT
//  x
/*=========================================================================*/

extern "C" void gpu_xscal_(int* cplx, int *N,
                           cuDoubleComplex *alpha,
                           void **X_ptr, int *incrx)
{

  CUDA_API_CHECK( (*cplx==1) ?
                  cublasDscal(cublas_handle,*N,
                              &((*alpha).x),
                              (double *)(*X_ptr), *incrx) :
                  cublasZscal(cublas_handle, *N,
                              alpha,
                              (cuDoubleComplex *)(*X_ptr), *incrx) );
} // gpu_xscal_

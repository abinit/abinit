/* gpu_linalg.cu */

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain hipblas and magma encapsulation routines,
 * that will be callable from fortran routines
 *
 */

#include <gpu_linalg.h>
#include <rocblas/rocblas.h>
#include <rocfft/rocfft.h>

// SYEVG/SYGVD are too slow when using HIP (on CRAY)
// We choose instead SYEVJ/SYGVJ (Jacobi) instead
//  with a ridiculous number of iterations
// This has to be tested forward...
hipblasHandle_t hipblas_handle;
hipsolverDnHandle_t hipsolverDn_handle;
static hipStream_t stream_compute;

// Parameters to HIP SYEVJ and SYGVJ (Jacobi-based eigensolver)
hipsolverSyevjInfo_t syevj_params = nullptr;
const int    MAX_SWEEPS  = 1;
const double TOLERANCE   = 1.e-5;

//! utility function for compatiblity between hipblas v1/v2 API
hipblasOperation_t select_hipblas_op(char *c)
{
  hipblasOperation_t op;

  if (*c == 'n' or *c == 'N')
    op = HIPBLAS_OP_N;
  else if (*c == 't' or *c == 'T')
    op = HIPBLAS_OP_T;
  else if (*c == 'c' or *c == 'C')
    op = HIPBLAS_OP_C;
  else
    printf("HIPBLAS API error, character can't be converted to a valid hipblas operation !\n");

  return op;
}

//! utility function for compatiblity between hipblas v1/v2 API
hipblasSideMode_t select_hipblas_side(char *c)
{
  hipblasSideMode_t mode;

  if (*c == 'l' or *c == 'L')
    mode = HIPBLAS_SIDE_LEFT;
  else if (*c == 'r' or *c == 'R')
    mode = HIPBLAS_SIDE_RIGHT;
  else
    printf("HIPBLAS API error, character can't be converted to a valid hipblas side mode !\n");

  return mode;
}

//! utility function for compatiblity between hipblas v1/v2 API
hipblasFillMode_t select_hipblas_fill_mode(const char *c)
{
  hipblasFillMode_t mode;

  if (*c == 'u' or *c == 'U')
    mode = HIPBLAS_FILL_MODE_UPPER;
  else if (*c == 'l' or *c == 'L')
    mode = HIPBLAS_FILL_MODE_LOWER;
  else
    printf("HIPBLAS API error, character can't be converted to a valid hipblas fill mode !\n");

  return mode;
}

//! utility function for compatiblity between hipblas v1/v2 API
hipblasDiagType_t select_hipblas_diag_type(char *c)
{
  hipblasDiagType_t diag;

  if (*c == 'n' or *c == 'N')
    diag = HIPBLAS_DIAG_NON_UNIT;
  else if (*c == 'u' or *c == 'U')
    diag = HIPBLAS_DIAG_UNIT;
  else
    printf("HIPBLAS API error, character can't be converted to a valid hipblas diag type !\n");

  return diag;
}

//! utility function to select eigen type
hipsolverEigType_t select_eigen_type(const int eigType_int)
{
  hipsolverEigType_t eigType = HIPSOLVER_EIG_TYPE_1;

  if (eigType_int < 1 or eigType_int > 3)
    printf("HIPSOLVER API error, input can't be converted to a valid hipsolver eigen type !\n");
  else if(eigType_int == 1)
    eigType = HIPSOLVER_EIG_TYPE_1;
  else if(eigType_int == 2)
    eigType = HIPSOLVER_EIG_TYPE_2;
  else if(eigType_int == 3)
    eigType = HIPSOLVER_EIG_TYPE_3;

  return eigType;
}

//! utility function to select eigen mode
hipsolverEigMode_t select_eigen_mode(const char* c)
{

  hipsolverEigMode_t eigMode = HIPSOLVER_EIG_MODE_NOVECTOR;

  if (*c =='n' or *c=='N')
    eigMode = HIPSOLVER_EIG_MODE_NOVECTOR;
  else if (*c =='v' or *c=='V')
    eigMode = HIPSOLVER_EIG_MODE_VECTOR;
  else
    printf("HIPSOLVER API error, character can't be converted to a valid eigen mode !\n");

  return eigMode;
}

hipsolverFillMode_t select_hipsolver_fill_mode(const char *c)
{
  hipsolverFillMode_t mode;

  if (*c == 'u' or *c == 'U')
    mode = HIPSOLVER_FILL_MODE_UPPER;
  else if (*c == 'l' or *c == 'L')
    mode = HIPSOLVER_FILL_MODE_LOWER;
  else
    printf("HIPBLAS API error, character can't be converted to a valid hipsolver fill mode !\n");

  return mode;
}



/*=========================================================================*/
// NAME
//  gpu_linalg_init
//
// FUNCTION
//  Initialisation of linalg environnement on GPU
//
/*=========================================================================*/

extern "C" void gpu_linalg_init_()
{
  HIP_API_CHECK( hipblasCreate(&hipblas_handle) );
  HIP_API_CHECK( hipsolverDnCreate(&hipsolverDn_handle) );
  HIP_API_CHECK( hipStreamCreate(&stream_compute) );
  HIP_API_CHECK( hipsolverDnSetStream(hipsolverDn_handle,stream_compute) );
  HIP_API_CHECK( hipblasSetStream(hipblas_handle,stream_compute) );

  HIP_API_CHECK(hipsolverCreateSyevjInfo(&syevj_params));
  HIP_API_CHECK(hipsolverXsyevjSetTolerance(syevj_params, TOLERANCE));
  HIP_API_CHECK(hipsolverXsyevjSetMaxSweeps(syevj_params, MAX_SWEEPS));
  //fprintf(stdout, "Initializing hipBLAS, this may take 10 seconds...\n");
  //fflush(stdout);
  rocblas_initialize();
  rocfft_setup();
}

/*=========================================================================*/
// NAME
//  gpu_linalg_shutdown
//
// FUNCTION
//  Close linalg environnement on GPU
//
/*=========================================================================*/

extern "C" void gpu_linalg_shutdown_()
{
  //FIXME
  //HIP_API_CHECK( hipblasDestroy(hipblas_handle) );
  hipblasDestroy(hipblas_handle);
  HIP_API_CHECK( hipsolverDestroySyevjInfo(syevj_params) );
  HIP_API_CHECK( hipsolverDnDestroy(hipsolverDn_handle) );
  HIP_API_CHECK( hipStreamDestroy(stream_compute) );
}

/*=========================================================================*/
// NAME
//  gpu_linalg_stream_synchronize
//
// FUNCTION
//  Wait for any linalg operations still running on stream
//
/*=========================================================================*/

extern "C" void gpu_linalg_stream_synchronize_()
{
  HIP_API_CHECK( hipStreamSynchronize(stream_compute) );
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
/*=========================================================================*/

extern "C" void gpu_xgemm_(int* cplx, char *transA, char *transB, int *N, int *M, int *K,
                           hipDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           hipDoubleComplex *beta, void** C_ptr, int *ldc)
{

  hipblasOperation_t opA = select_hipblas_op(transA);
  hipblasOperation_t opB = select_hipblas_op(transB);

  (*cplx==1)?
    hipblasDgemm(hipblas_handle, opA, opB, *N, *M, *K, &((*alpha).x),
                (double *)(*A_ptr), *lda,
                (double *)(*B_ptr), *ldb,
                &((*beta).x),
                (double *)(*C_ptr), *ldc) :
    hipblasZgemm(hipblas_handle, opA, opB, *N, *M, *K, (hipblasDoubleComplex *) alpha,
                (hipblasDoubleComplex *)(*A_ptr), *lda,
                (hipblasDoubleComplex *)(*B_ptr), *ldb,
                (hipblasDoubleComplex *) beta,
                (hipblasDoubleComplex *)(*C_ptr), *ldc);

} // gpu_xgemm

/*=========================================================================*/
// NAME
//  gpu_xgemm_strided_batched
//
// FUNCTION
//  Compute a batched scalar-matrix-matrix product and return a scalar-matrix product on GPU.
//  Meant to be used on non-contiguous matrixes with data is uniformly split in the same number of batches in each matrix.
//  c = alpha * op(a) * op(b) + beta * c
//
// INPUTS
//  cplx       = 1 if real 2 if complex
//  transa     = from of op(a) to be used in the matrix multiplication
//  transb     = from of op(b) to be used in the matrix multiplication
//  m          = number of rows of the matrix op(a) and of the matrix c
//  n          = number of rows of the matrix op(b) and the number of columns of the matrix c
//  k          = number of columns of the matrix op(a) and the number of rows of the matrix op(b)
//  alpha      = alpha scalar coefficient for matrix op(a)
//  a_gpu      = pointer to gpu memory location of  matrix a
//  lda        = first dimension of a
//  strideC    = stride between each batch in matrix a
//  b_gpu      = pointer to gpu memory location of  matrix b
//  ldb        = first dimension of b
//  strideC    = stride between each batch in matrix b
//  beta       = beta scalar coefficient for matrix c
//  c_gpu      = pointer to gpu memory location of  matrix c
//  ldc        = first dimension of c
//  strideC    = stride between each batch in matrix c
//  batchCount = number of batches in any matrix
//
// OUTPUT
//  c_gpu     = c_gpu matrix
//
/*=========================================================================*/

extern "C" void gpu_xgemm_strided_batched_(int* cplx, char *transA, char *transB,
                           int *N, int *M, int *K,
                           hipDoubleComplex *alpha,
                           void **A_ptr, int *lda, int *strideA,
                           void** B_ptr, int *ldb, int *strideB,
                           hipDoubleComplex *beta,
                           void** C_ptr, int *ldc, int *strideC, int *batchCount)
{

  hipblasOperation_t opA = select_hipblas_op(transA);
  hipblasOperation_t opB = select_hipblas_op(transB);

  (*cplx==1)?
    HIP_API_CHECK( hipblasDgemmStridedBatched(hipblas_handle, opA, opB,
                *N, *M, *K, &((*alpha).x),
                (double *)(*A_ptr), *lda, *strideA,
                (double *)(*B_ptr), *ldb, *strideB,
                &((*beta).x),
                (double *)(*C_ptr), *ldc, *strideC, *batchCount)) :
    HIP_API_CHECK( hipblasZgemmStridedBatched(hipblas_handle, opA, opB,
                *N, *M, *K, (hipblasDoubleComplex *) alpha,
                (hipblasDoubleComplex *)(*A_ptr), *lda, *strideA,
                (hipblasDoubleComplex *)(*B_ptr), *ldb, *strideB,
                (hipblasDoubleComplex *) beta,
                (hipblasDoubleComplex *)(*C_ptr), *ldc, *strideC, *batchCount));
} // gpu_xgemm_strided_batched

/*=========================================================================*/
// NAME
//  gpu_xsymm
//
// FUNCTION
//  Compute a symmetric scalar-matrix-matrix product and return a scalar-matrix product on GPU
//  c = alpha * op(a) * op(b) + beta * c    , if side == L
//  c = alpha * op(b) * op(a) + beta * c    , if side == R
//
// INPUTS
//  cplx  = 1 if real 2 if complex
//  side  = Specifies whether op(a) appears on the left or right of x for
//          the operation to be performed as follows:
//              L or l op(a)*x = alpha*b
//              R or r x*op(a) = alpha*b
//  uplo  = Specifies whether the matrix a is an upper or lower triangular
//          matrix as follows:
//              U or u Matrix a is an upper triangular matrix.
//              L or l Matrix a is a lower triangular matrix
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
/*=========================================================================*/

extern "C" void gpu_xsymm_(int* cplx, char *side, char *uplo, int *N, int *M,
                           hipDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           hipDoubleComplex *beta, void** C_ptr, int *ldc)
{

  hipblasSideMode_t sideMode = select_hipblas_side(side);
  hipblasFillMode_t fillMode = select_hipblas_fill_mode(uplo);

  (*cplx==1)?
    HIP_API_CHECK( hipblasDsymm(hipblas_handle, sideMode, fillMode, *N, *M, &((*alpha).x),
                (double *)(*A_ptr), *lda,
                (double *)(*B_ptr), *ldb,
                &((*beta).x),
                (double *)(*C_ptr), *ldc)) :
    HIP_API_CHECK( hipblasZsymm(hipblas_handle, sideMode, fillMode, *N, *M, (hipblasDoubleComplex *) alpha,
                (hipblasDoubleComplex *)(*A_ptr), *lda,
                (hipblasDoubleComplex *)(*B_ptr), *ldb,
                (hipblasDoubleComplex *) beta,
                (hipblasDoubleComplex *)(*C_ptr), *ldc));
} // gpu_xsymm_omp

/*=========================================================================*/
// NAME
//  gpu_xhemm
//
// FUNCTION
//  Compute a Hermitian scalar-matrix-matrix product and return a scalar-matrix product on GPU
//  c = alpha * op(a) * op(b) + beta * c    , if side == L
//  c = alpha * op(b) * op(a) + beta * c    , if side == R
//
// INPUTS
//  cplx  = 1 if real 2 if complex
//  side  = Specifies whether op(a) appears on the left or right of x for
//          the operation to be performed as follows:
//              L or l op(a)*x = alpha*b
//              R or r x*op(a) = alpha*b
//  uplo  = Specifies whether the matrix a is an upper or lower triangular
//          matrix as follows:
//              U or u Matrix a is an upper triangular matrix.
//              L or l Matrix a is a lower triangular matrix
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
/*=========================================================================*/

extern "C" void gpu_zhemm_(char *side, char *uplo, int *N, int *M,
                           hipDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           hipDoubleComplex *beta, void** C_ptr, int *ldc)
{

  hipblasSideMode_t sideMode = select_hipblas_side(side);
  hipblasFillMode_t fillMode = select_hipblas_fill_mode(uplo);

  HIP_API_CHECK( hipblasZhemm(hipblas_handle, sideMode, fillMode, *N, *M, (hipblasDoubleComplex *) alpha,
              (hipblasDoubleComplex *)(*A_ptr), *lda,
              (hipblasDoubleComplex *)(*B_ptr), *ldb,
              (hipblasDoubleComplex *) beta,
              (hipblasDoubleComplex *)(*C_ptr), *ldc));
} // gpu_zhemm_omp

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
                           int *N, int *M, hipDoubleComplex *alpha,
                           void **A_ptr, int *ldA,
                           void** B_ptr, int *ldB)
{

  hipblasSideMode_t sideMode = select_hipblas_side(side);
  hipblasFillMode_t fillMode = select_hipblas_fill_mode(uplo);
  hipblasDiagType_t diagType = select_hipblas_diag_type(diag);
  hipblasOperation_t opA     = select_hipblas_op(transA);

  (*cplx==1) ?
    hipblasDtrsm(hipblas_handle, sideMode, fillMode, opA, diagType,
                *N, *M, &((*alpha).x),
                (double *)(*A_ptr), *ldA,
                (double *)(*B_ptr), *ldB):
    hipblasZtrsm(hipblas_handle, sideMode, fillMode, opA, diagType,
                *N, *M, (hipblasDoubleComplex *) alpha,
                (hipblasDoubleComplex *)(*A_ptr), *ldA,
                (hipblasDoubleComplex *)(*B_ptr), *ldB);
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
                           hipDoubleComplex *alpha,
                           void **X_ptr, int *incrx, void** Y_ptr, int *incry)
{
  HIP_API_CHECK( (*cplx==1) ?
                  hipblasDaxpy(hipblas_handle,*N,
                              &((*alpha).x),
                              (double *)(*X_ptr), *incrx,
                              (double *)(*Y_ptr), *incry) :
                  hipblasZaxpy(hipblas_handle, *N,
                              (hipblasDoubleComplex *) alpha,
                              (hipblasDoubleComplex *)(*X_ptr), *incrx,
                              (hipblasDoubleComplex *)(*Y_ptr), *incry) );
} // gpu_xaxpy_

/*=========================================================================*/
// NAME
//  gpu_xcopy
//
// FUNCTION
//  Compute blas-1 COPY on GPU
//  y = x (copy x into y)
//
// INPUTS / OUTPUTS
//  cplx  = 1 if real 2 if complex
//  N     = input vector size
//  x     = input vector
//  incrx = increment between two consecutive values of x
//  y     = output vector
//  incry = increment between two consecutive values of y
//
/*=========================================================================*/

extern "C" void gpu_xcopy_(int* cplx, int *N,
                           void **X_ptr, int *incrx,
                           void **Y_ptr, int *incry)
{

  HIP_API_CHECK( (*cplx==1) ?
                  hipblasDcopy(hipblas_handle,*N,
                              (const double *)(*X_ptr), *incrx,
                              (      double *)(*Y_ptr), *incry) :
                  hipblasZcopy(hipblas_handle, *N,
                              (const hipblasDoubleComplex *)(*X_ptr), *incrx,
                              (      hipblasDoubleComplex *)(*Y_ptr), *incry) );
} // gpu_xcopy_

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
//  N     = vector size
//  alpha = scalar complex value (when doing computation with real value, only the real part is used)
//  x     = input/output vector
//  incrx = increment between two consecutive values of x
//
// OUTPUT
//  x
/*=========================================================================*/

extern "C" void gpu_xscal_(int* cplx, int *N,
                           hipDoubleComplex *alpha,
                           void **X_ptr, int *incrx)
{
  HIP_API_CHECK( (*cplx==1) ?
                  hipblasDscal(hipblas_handle,*N,
                              &((*alpha).x),
                              (double *)(*X_ptr), *incrx) :
                  hipblasZscal(hipblas_handle, *N,
                              (hipblasDoubleComplex *) alpha,
                              (hipblasDoubleComplex *)(*X_ptr), *incrx) );
} // gpu_xscal_

/*=========================================================================*/
// NAME
//  gpu_xpotrf
//
// FUNCTION
//  Compute a LAPACK SYGVD operation on GPU
//  computes the Cholesky factorization of a Hermitian positive-definite matrix
//
// See LAPACK doc in reference implementation:
// https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dpotrf.f
//
// INPUTS
//  uplo  = character, 'u' or 'l'
//  n     = matrix size
//  A_ptr = pointer to gpu memory location of matrix A
//  lda   = leading dimension of matrix A
//  work_ptr =
//  lwork =
//  devInfo  =
//
/*=========================================================================*/

extern "C" void gpu_xpotrf_(const int* cplx,
                            const char* uplo,
                            const int* n,
                            void **A_ptr,
                            const int *lda,
                            void** work_ptr,
                            int* lwork,
                            int* info)
{

  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  int* devInfo;
  HIP_API_CHECK( hipMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDpotrf(hipsolverDn_handle,
                                       fillMode, *n,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZpotrf(hipsolverDn_handle,
                                       fillMode, *n,
                                       (hipDoubleComplex*)(*A_ptr), *lda,
                                       (hipDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo) );
    }

  HIP_API_CHECK( hipMemcpy(info, devInfo, 1*sizeof(int), hipMemcpyDefault) );

  HIP_API_CHECK( hipFree(devInfo) );

} // gpu_xpotrf_

extern "C" void gpu_xpotrf_buffersize_(const int* cplx,
                                       const char* uplo,
                                       const int* n,
                                       void **A_ptr,
                                       const int *lda,
                                       int* lwork)
{

  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDpotrf_bufferSize(hipsolverDn_handle,
                                                  fillMode, *n,
                                                  (double*)(*A_ptr), *lda,
                                                  lwork) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZpotrf_bufferSize(hipsolverDn_handle,
                                                  fillMode, *n,
                                                  (hipDoubleComplex*)(*A_ptr), *lda,
                                                  lwork) );
    }

} // gpu_xpotrf_bufferSize_

/*=========================================================================*/
// NAME
//  gpu_xsygvd
//
// FUNCTION
//  Compute a LAPACK SYGVD operation on GPU
//  compute eigen values/vectors of a real generalized
//  symmetric-definite eigenproblem
//
// See LAPACK doc in reference implementation:
// https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dsygvd.f
//
// INPUTS
//  itype = integer, type of problem
//  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
//  uplo  = character, 'u' or 'l'
//  A_nrows = matrix size
//  A_ptr = pointer to gpu memory location of matrix A
//  lda   = leading dimension of matrix A
//  B_ptr = pointer to gpu memory location of matrix B
//  ldb   = leading dimension of matrix B
//  W_ptr = pointer to gpu memory location of matrix W (output eigen values)
//  work_ptr =
//  lwork =
//  devInfo  =
//
/*=========================================================================*/

extern "C" void gpu_xsygvd_(const int* cplx,
                            const int* itype,
                            const char* jobz,
                            const char* uplo,
                            const int* A_nrows,
                            void **A_ptr,
                            const int *lda,
                            void** B_ptr,
                            const int* ldb,
                            void** W_ptr,
                            void** work_ptr,
                            int* lwork,
                            int* info)
{

  const hipsolverEigType_t itype_cu = select_eigen_type(*itype);
  const hipsolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  int* devInfo;
  HIP_API_CHECK( hipMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDsygvd(hipsolverDn_handle, itype_cu,
                                       jobz_cu, fillMode, *A_nrows,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*B_ptr), *ldb,
                                       (double*)(*W_ptr),
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZhegvj(hipsolverDn_handle, itype_cu,
                                       jobz_cu, fillMode, *A_nrows,
                                       (hipDoubleComplex*)(*A_ptr), *lda,
                                       (hipDoubleComplex*)(*B_ptr), *ldb,
                                       (double*)(*W_ptr),
                                       (hipDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo, syevj_params) );
    }

  HIP_API_CHECK( hipMemcpy(info, devInfo, 1*sizeof(int), hipMemcpyDefault) );

  HIP_API_CHECK( hipFree(devInfo) );

} // gpu_xsygvd_


extern "C" void gpu_xsygvd_buffersize_(const int* cplx,
                                       const int* itype,
                                       const char* jobz,
                                       const char* uplo,
                                       const int* A_nrows,
                                       void **A_ptr,
                                       const int *lda,
                                       void** B_ptr,
                                       const int* ldb,
                                       void** W_ptr,
                                       int* lwork)
{

  const hipsolverEigType_t itype_cu = select_eigen_type(*itype);
  const hipsolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDsygvd_bufferSize(hipsolverDn_handle, itype_cu,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (double*)(*A_ptr), *lda,
                                                  (double*)(*B_ptr), *ldb,
                                                  (double*)(*W_ptr), lwork) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZhegvj_bufferSize(hipsolverDn_handle, itype_cu,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (hipDoubleComplex*)(*A_ptr), *lda,
                                                  (hipDoubleComplex*)(*B_ptr), *ldb,
                                                  (double*)(*W_ptr), lwork, syevj_params) );
    }

} // gpu_xsygvd_bufferSize_


/*=========================================================================*/
// NAME
//  gpu_xsyevd
//
// FUNCTION
//  Compute a LAPACK SYGVD operation on GPU
//  compute eigen values/vectors of a real generalized
//  symmetric-definite eigenproblem
//
// See LAPACK doc in reference implementation:
// https://github.com/Reference-LAPACK/lapack/blob/master/SRC/dsyevd.f
//
// INPUTS
//  itype = integer, type of problem
//  jobz  = character, 'n'(eigenvalues only) or 'v' (eigenvalues + eigenvectors)
//  uplo  = character, 'u' or 'l'
//  A_nrows = matrix size
//  A_ptr = pointer to gpu memory location of matrix A
//  lda   = leading dimension of matrix A
//  B_ptr = pointer to gpu memory location of matrix B
//  ldb   = leading dimension of matrix B
//  W_ptr = pointer to gpu memory location of matrix W (output eigen values)
//  work_ptr =
//  lwork =
//  devInfo  =
//
/*=========================================================================*/

extern "C" void gpu_xsyevd_(const int* cplx,
                            const char* jobz,
                            const char* uplo,
                            const int* A_nrows,
                            void **A_ptr,
                            const int *lda,
                            void** W_ptr,
                            void** work_ptr,
                            int* lwork,
                            int* info)
{

  const hipsolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  int* devInfo;
  HIP_API_CHECK( hipMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDsyevd(hipsolverDn_handle,
                                       jobz_cu, fillMode, *A_nrows,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*W_ptr),
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZheevd(hipsolverDn_handle,
                                       jobz_cu, fillMode, *A_nrows,
                                       (hipDoubleComplex*)(*A_ptr), *lda,
                                       (double*)(*W_ptr),
                                       (hipDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo) );
    }

  HIP_API_CHECK( hipMemcpy(info, devInfo, 1*sizeof(int), hipMemcpyDefault) );

  HIP_API_CHECK( hipFree(devInfo) );

} // gpu_xsyevd_

extern "C" void gpu_xsyevd_buffersize_(const int* cplx,
                                       const char* jobz,
                                       const char* uplo,
                                       const int* A_nrows,
                                       void **A_ptr,
                                       const int *lda,
                                       void** W_ptr,
                                       int* lwork)
{

  const hipsolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const hipsolverFillMode_t  fillMode = select_hipsolver_fill_mode(uplo);

  if (*cplx == 1)
    {
      HIP_API_CHECK( hipsolverDnDsyevd_bufferSize(hipsolverDn_handle,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (double*)(*A_ptr), *lda,
                                                  (double*)(*W_ptr), lwork) );
    }
  else
    {
      HIP_API_CHECK( hipsolverDnZheevd_bufferSize(hipsolverDn_handle,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (hipDoubleComplex*)(*A_ptr), *lda,
                                                  (double*)(*W_ptr), lwork) );
    }

} // gpu_xsyevd_bufferSize_

/* gpu_linalg.cu */

/*
 * Copyright (C) 2008-2024 ABINIT Group (MMancini,FDahm)
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
static cudaStream_t stream_compute;

//! utility function for compatiblity between cublas v1/v2 API
cublasOperation_t select_cublas_op(const char *c)
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

//! utility function for compatiblity between cublas v1/v2 API
cublasSideMode_t select_cublas_side(const char *c)
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

//! utility function for compatiblity between cublas v1/v2 API
cublasFillMode_t select_cublas_fill_mode(const char *c)
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

//! utility function for compatiblity between cublas v1/v2 API
cublasDiagType_t select_cublas_diag_type(const char *c)
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

//! utility function to select eigen type
cusolverEigType_t select_eigen_type(const int eigType_int)
{
  cusolverEigType_t eigType = CUSOLVER_EIG_TYPE_1;

  if (eigType_int < 1 or eigType_int > 3)
    printf("CUSOLVER API error, input can't be converted to a valid cusolver eigen type !\n");
  else
    eigType = static_cast<cusolverEigType_t>(eigType_int);

  return eigType;
}

//! utility function to select eigen mode
cusolverEigMode_t select_eigen_mode(const char* c)
{

  cusolverEigMode_t eigMode = CUSOLVER_EIG_MODE_NOVECTOR;

  if (*c =='n' or *c=='N')
    eigMode = CUSOLVER_EIG_MODE_NOVECTOR;
  else if (*c =='v' or *c=='V')
    eigMode = CUSOLVER_EIG_MODE_VECTOR;
  else
    printf("CUSOLVER API error, character can't be converted to a valid eigen mode !\n");

  return eigMode;
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
  CUDA_API_CHECK( cublasCreate(&cublas_handle) );
  CUDA_API_CHECK( cusolverDnCreate(&cusolverDn_handle) );
  CUDA_API_CHECK( cudaStreamCreate(&stream_compute) );
  CUDA_API_CHECK( cusolverDnSetStream(cusolverDn_handle,stream_compute) );
  CUDA_API_CHECK( cublasSetStream(cublas_handle,stream_compute) );
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
  CUDA_API_CHECK( cublasDestroy(cublas_handle) );
  CUDA_API_CHECK( cusolverDnDestroy(cusolverDn_handle) );
  CUDA_API_CHECK( cudaStreamDestroy(stream_compute) );
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
  CUDA_API_CHECK( cudaStreamSynchronize(stream_compute) );
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
                           cuDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           cuDoubleComplex *beta, void** C_ptr, int *ldc)
{

  cublasOperation_t opA = select_cublas_op(transA);
  cublasOperation_t opB = select_cublas_op(transB);

  (*cplx==1)?
    CUDA_API_CHECK( cublasDgemm(cublas_handle, opA, opB, *N, *M, *K, &((*alpha).x),
                (double *)(*A_ptr), *lda,
                (double *)(*B_ptr), *ldb,
                &((*beta).x),
                (double *)(*C_ptr), *ldc)) :
    CUDA_API_CHECK( cublasZgemm(cublas_handle, opA, opB, *N, *M, *K, alpha,
                (cuDoubleComplex *)(*A_ptr), *lda,
                (cuDoubleComplex *)(*B_ptr), *ldb,
                beta,
                (cuDoubleComplex *)(*C_ptr), *ldc));

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
                           cuDoubleComplex *alpha,
                           void **A_ptr, int *lda, int *strideA,
                           void** B_ptr, int *ldb, int *strideB,
                           cuDoubleComplex *beta,
                           void** C_ptr, int *ldc, int *strideC, int *batchCount)
{

  cublasOperation_t opA = select_cublas_op(transA);
  cublasOperation_t opB = select_cublas_op(transB);

  (*cplx==1)?
    CUDA_API_CHECK( cublasDgemmStridedBatched(cublas_handle, opA, opB,
                *N, *M, *K, &((*alpha).x),
                (double *)(*A_ptr), *lda, *strideA,
                (double *)(*B_ptr), *ldb, *strideB,
                &((*beta).x),
                (double *)(*C_ptr), *ldc, *strideC, *batchCount)) :
    CUDA_API_CHECK( cublasZgemmStridedBatched(cublas_handle, opA, opB,
                *N, *M, *K, alpha,
                (cuDoubleComplex *)(*A_ptr), *lda, *strideA,
                (cuDoubleComplex *)(*B_ptr), *ldb, *strideB,
                beta,
                (cuDoubleComplex *)(*C_ptr), *ldc, *strideC, *batchCount));
} // gpu_xgemm_strided_batched


/*=========================================================================*/
// NAME
//  gpu_xginv_strided_batched
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

extern "C" void gpu_xgetrf_buffersize_(int* cplx,
                           int *N,
                           void **A_ptr, int *lda,
                           size_t *bufferSize_host, size_t *bufferSize_device)
{
  cudaDataType type;
  if(*cplx==1) type=CUDA_R_64F;
  if(*cplx==2) type=CUDA_C_64F;
  CUDA_API_CHECK( cusolverDnXgetrf_bufferSize(cusolverDn_handle,
    NULL,
    *N,
    *N,
    type,
    (void*)A_ptr,
    *lda,
    type,
    bufferSize_device,
    bufferSize_host));
}

extern "C" void gpu_xginv_(int* cplx,
                           int *N,
                           void **A_ptr, int *lda, void **W_ptr,
                           size_t *bufferSize_host,
                           void **work_ptr, size_t *bufferSize_device,
                           int *info)
{
  cudaDataType type;
  if(*cplx==1) type=CUDA_R_64F;
  if(*cplx==2) type=CUDA_C_64F;
  int *devInfo;
  void *bufferOnHost;
  //Complex *matrix = (Complex*) malloc(*N * (*N) * sizeof(Complex));
  //create_identity_matrix(*N, matrix);
  int64_t *ipiv;
  //cusolverDnSetAdvOptions(params, CUSOLVERDN_GETRF, CUSOLVER_ALG_0);
  CUDA_API_CHECK( cudaMalloc( (void**) &ipiv,     *N*sizeof(int64_t)) );
  CUDA_API_CHECK( cudaMalloc( (void**) &devInfo, 1*sizeof(int)) );
  //CUDA_API_CHECK( cudaMemcpy(W_ptr, matrix, *N*(*N)*sizeof(cuDoubleComplex), cudaMemcpyDefault) );

  bufferOnHost = malloc(*bufferSize_host);
  CUDA_API_CHECK( cusolverDnXgetrf(cusolverDn_handle,
    NULL,
    *N,
    *N,
    type,
    (void*)A_ptr,
    *lda,
    ipiv,
    type,
    work_ptr,
    *bufferSize_device,
    bufferOnHost,
    *bufferSize_host,
    devInfo ));
  CUDA_API_CHECK( cudaMemcpy(info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );
  if (*info < 0) {
    fprintf(stderr, "Xgetrf:  The %d-th argument of ZGETRF had an illegal value.", *info);
    fflush(stderr);
    abi_cabort();
  } else if (*info > 0) {
     fprintf(stderr, "Xgetrf: The matrix that has been passed in argument is probably either singular or nearly singular.\n\
     U(i,i) in the P*L*U factorization is exactly zero for i = %d \n\
     The factorization has been completed but the factor U is exactly singular.\n\
     Division by zero will occur if it is used to solve a system of equations.", *info);
     fflush(stderr);
     abi_cabort();
  } else { fprintf(stderr, "Xgetrf: OK");}
  free(bufferOnHost);
  cusolverDnXgetrs(
    cusolverDn_handle,
    NULL,
    CUBLAS_OP_N,
    *N,
    *N,
    type,
    (void*)A_ptr,
    *lda,
    ipiv,
    type,
    W_ptr,
    *lda,
    devInfo );
  CUDA_API_CHECK( cudaMemcpy(info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );
  if (*info < 0) {
    fprintf(stderr, "Xgetrs: The %d-th argument of ZGETRS had an illegal value.", *info);
    fflush(stderr);
    abi_cabort();
  } else if (*info > 0) {
     fprintf(stderr, "Xgetrs: The matrix that has been passed in argument is probably either singular or nearly singular.\n\
     U(i,i) in the P*L*U factorization is exactly zero for i = %d \n\
     The factorization has been completed but the factor U is exactly singular.\n\
     Division by zero will occur if it is used to solve a system of equations.", *info);
     fflush(stderr);
     abi_cabort();
  }
  CUDA_API_CHECK( cudaMemcpy(A_ptr, W_ptr, *N*(*N)*sizeof(cuDoubleComplex), cudaMemcpyDefault) );
  CUDA_API_CHECK( cudaFree(ipiv) );
  CUDA_API_CHECK( cudaFree(devInfo) );
} // gpu_xginv

typedef struct {
    double real;
    double imag;
} Complex;

void create_identity_matrix(int N, Complex *matrix) {
    // Remplir la matrice avec des valeurs complexes
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                matrix[i * N + j].real = 1.0;  // Partie réelle de 1
                matrix[i * N + j].imag = 0.0;  // Partie imaginaire de 0
            } else {
                matrix[i * N + j].real = 0.0;  // Partie réelle de 0
                matrix[i * N + j].imag = 0.0;  // Partie imaginaire de 0
            }
        }
    }
}

extern "C" void gpu_xginv_old_(int* cplx,
                           int *N,
                           void **A_ptr, int *lda, void **W_ptr)
{
  int info;
  int *devInfo;
  size_t workspaceInBytesOnDevice,workspaceInBytesOnHost;
  void *bufferOnDevice, *bufferOnHost;
  int64_t *ipiv;
  //cusolverDnSetAdvOptions(params, CUSOLVERDN_GETRF, CUSOLVER_ALG_0);
  CUDA_API_CHECK( cudaMalloc( (void**) &ipiv,     *N*sizeof(int64_t)) );
  CUDA_API_CHECK( cudaMalloc( (void**) &devInfo, 1*sizeof(int)) );

  CUDA_API_CHECK( cusolverDnXgetrf_bufferSize(cusolverDn_handle,
    NULL,
    *N,
    *N,
    CUDA_C_64F,
    (void*)A_ptr,
    *lda,
    CUDA_C_64F,
    &workspaceInBytesOnDevice,
    &workspaceInBytesOnHost));
  bufferOnHost = malloc(workspaceInBytesOnHost);
  CUDA_API_CHECK( cudaMalloc( (void**) &bufferOnDevice, workspaceInBytesOnDevice) );
  CUDA_API_CHECK( cusolverDnXgetrf(cusolverDn_handle,
    NULL,
    *N,
    *N,
    CUDA_C_64F,
    (void*)A_ptr,
    *lda,
    ipiv,
    CUDA_C_64F,
    bufferOnDevice,
    workspaceInBytesOnDevice,
    bufferOnHost,
    workspaceInBytesOnHost,
    devInfo ));
  CUDA_API_CHECK( cudaMemcpy(&info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );
  if (info < 0) {
    fprintf(stderr, " The %d-th argument of ZGETRF had an illegal value.", info);
    fflush(stderr);
    abi_cabort();
  } else if (info > 0) {
     fprintf(stderr, "The matrix that has been passed in argument is probably either singular or nearly singular.\n\
     U(i,i) in the P*L*U factorization is exactly zero for i = %d \n\
     The factorization has been completed but the factor U is exactly singular.\n\
     Division by zero will occur if it is used to solve a system of equations.", info);
     fflush(stderr);
     abi_cabort();
  }
  CUDA_API_CHECK( cudaFree(bufferOnDevice) );
  free(bufferOnHost);
  cusolverDnXgetrs(
    cusolverDn_handle,
    NULL,
    CUBLAS_OP_N,
    *N,
    *N,
    CUDA_C_64F,
    (void*)A_ptr,
    *lda,
    ipiv,
    CUDA_C_64F,
    W_ptr,
    *lda,
    devInfo );
  CUDA_API_CHECK( cudaMemcpy(&info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );
  if (info < 0) {
    fprintf(stderr, " The %d-th argument of ZGETRF had an illegal value.", info);
    fflush(stderr);
    abi_cabort();
  } else if (info > 0) {
     fprintf(stderr, "The matrix that has been passed in argument is probably either singular or nearly singular.\n\
     U(i,i) in the P*L*U factorization is exactly zero for i = %d \n\
     The factorization has been completed but the factor U is exactly singular.\n\
     Division by zero will occur if it is used to solve a system of equations.", info);
     fflush(stderr);
     abi_cabort();
  }
  CUDA_API_CHECK( cudaMemcpy(A_ptr, W_ptr, *N*(*N)*sizeof(cuDoubleComplex), cudaMemcpyDefault) );
  CUDA_API_CHECK( cudaFree(ipiv) );
  CUDA_API_CHECK( cudaFree(devInfo) );
	/*
  int *infoArray, *pivotArray, *info;
	info = (int*) malloc(*batchCount*sizeof(int));
  CUDA_API_CHECK( cudaMalloc( (void**) &infoArray,  *batchCount*sizeof(int)) );
  CUDA_API_CHECK( cudaMalloc( (void**) &pivotArray, *batchCount*sizeof(int)) );
  cuDoubleComplex **A_array = (cuDoubleComplex**) malloc(*batchCount*sizeof(cuDoubleComplex*));
  cuDoubleComplex **C_array = (cuDoubleComplex**) malloc(*batchCount*sizeof(cuDoubleComplex*));
	int i;
	for(i=0; i<*batchCount; ++i) {
		fprintf(stderr, "5\n"); fflush(stderr);
    A_array[i]=(cuDoubleComplex *)(A_ptr+i*(*strideA)) ;
		CUDA_API_CHECK( cudaMalloc( (void**) &(C_array[i]), *strideA*sizeof(cuDoubleComplex)) );
		fprintf(stderr, "56\n"); fflush(stderr);
		toto((void**)&(A_ptr),"toto");
		toto((void**)&(A_array[i]),"totoA");
		toto((void**)&(C_array[i]),"totoC");
	}
	fprintf(stderr, "6\n"); fflush(stderr);

  if(*cplx==1) {
  *  CUDA_API_CHECK( cublasDgetrfBatched(cublas_handle,
                *N,
                A_array, *lda,
                pivotArray, infoArray, *batchCount));
    CUDA_API_CHECK( cublasDgetriBatched(cublas_handle,
                *N,
                (double *)(*A_ptr), *lda,
                pivotArray,
                (double *)(*A_ptr), *lda,
                infoArray, *batchCount));
  *
  } else {
    CUDA_API_CHECK( cublasZgetrfBatched(cublas_handle,
                *N,
                A_array, *lda,
                pivotArray, infoArray, *batchCount));
		CUDA_API_CHECK( cudaMemcpy(info, infoArray, *batchCount*sizeof(int), cudaMemcpyDefault) );
		for(i=0; i<*batchCount; ++i) {
			if(info[i]!=0){
				fprintf(stderr, "ERROR: zgetrf on batch %d / %d : %d \n",i,*batchCount,info[i]);
				fflush(stderr);
				abi_cabort();
			} else { fprintf(stderr, "OK: zgetrf on batch %d / %d : %d \n",i,*batchCount,info[i]);}
		}
    CUDA_API_CHECK( cublasZgetriBatched(cublas_handle,
                *N,
                A_array, *lda,
                pivotArray,
                C_array, *lda,
                infoArray, *batchCount));
		CUDA_API_CHECK( cudaMemcpy(info, infoArray, *batchCount*sizeof(int), cudaMemcpyDefault) );
		for(i=0; i<*batchCount; ++i) {
			if(info[i]!=0){
				fprintf(stderr, "ERROR: zgetri on batch %d / %d : %d \n",i,*batchCount,info[i]);
				fflush(stderr);
				abi_cabort();
			} else { fprintf(stderr, "OK: zgetri on batch %d / %d : %d \n",i,*batchCount,info[i]);}
		}
  }


  CUDA_API_CHECK( cudaFree(infoArray) );
  CUDA_API_CHECK( cudaFree(pivotArray) );
	free(info);
	for(i=0; i<*batchCount; ++i) {
		CUDA_API_CHECK( cudaMemcpy(A_array[i], C_array[i], *strideA*sizeof(cuDoubleComplex), cudaMemcpyDefault) );
		CUDA_API_CHECK( cudaFree(C_array[i]) );
	}
	free(A_array);
	free(C_array);
	*/
} // gpu_xginv_strided_batched

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
                           cuDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           cuDoubleComplex *beta, void** C_ptr, int *ldc)
{

  cublasSideMode_t sideMode = select_cublas_side(side);
  cublasFillMode_t fillMode = select_cublas_fill_mode(uplo);

  (*cplx==1)?
    CUDA_API_CHECK( cublasDsymm(cublas_handle, sideMode, fillMode, *N, *M, &((*alpha).x),
                (double *)(*A_ptr), *lda,
                (double *)(*B_ptr), *ldb,
                &((*beta).x),
                (double *)(*C_ptr), *ldc)) :
    CUDA_API_CHECK( cublasZsymm(cublas_handle, sideMode, fillMode, *N, *M, alpha,
                (cuDoubleComplex *)(*A_ptr), *lda,
                (cuDoubleComplex *)(*B_ptr), *ldb,
                beta,
                (cuDoubleComplex *)(*C_ptr), *ldc));
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
                           cuDoubleComplex *alpha,
                           void **A_ptr, int *lda, void** B_ptr, int *ldb,
                           cuDoubleComplex *beta, void** C_ptr, int *ldc)
{

  cublasSideMode_t sideMode = select_cublas_side(side);
  cublasFillMode_t fillMode = select_cublas_fill_mode(uplo);

  CUDA_API_CHECK( cublasZhemm(cublas_handle, sideMode, fillMode, *N, *M, alpha,
              (cuDoubleComplex *)(*A_ptr), *lda,
              (cuDoubleComplex *)(*B_ptr), *ldb,
              beta,
              (cuDoubleComplex *)(*C_ptr), *ldc));
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

  CUDA_API_CHECK( (*cplx==1) ?
                  cublasDcopy(cublas_handle,*N,
                              (const double *)(*X_ptr), *incrx,
                              (      double *)(*Y_ptr), *incry) :
                  cublasZcopy(cublas_handle, *N,
                              (const cuDoubleComplex *)(*X_ptr), *incrx,
                              (      cuDoubleComplex *)(*Y_ptr), *incry) );
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

  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  int* devInfo;
  CUDA_API_CHECK( cudaMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDpotrf(cusolverDn_handle,
                                       fillMode, *n,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZpotrf(cusolverDn_handle,
                                       fillMode, *n,
                                       (cuDoubleComplex*)(*A_ptr), *lda,
                                       (cuDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo) );
    }

  CUDA_API_CHECK( cudaMemcpy(info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );

  CUDA_API_CHECK( cudaFree(devInfo) );

} // gpu_xpotrf_

extern "C" void gpu_xpotrf_buffersize_(const int* cplx,
                                       const char* uplo,
                                       const int* n,
                                       void **A_ptr,
                                       const int *lda,
                                       int* lwork)
{

  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDpotrf_bufferSize(cusolverDn_handle,
                                                  fillMode, *n,
                                                  (double*)(*A_ptr), *lda,
                                                  lwork) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZpotrf_bufferSize(cusolverDn_handle,
                                                  fillMode, *n,
                                                  (cuDoubleComplex*)(*A_ptr), *lda,
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

  const cusolverEigType_t itype_cu = select_eigen_type(*itype);
  const cusolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  int* devInfo;
  CUDA_API_CHECK( cudaMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDsygvd(cusolverDn_handle, itype_cu,
                                       jobz_cu, fillMode, *A_nrows,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*B_ptr), *ldb,
                                       (double*)(*W_ptr),
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZhegvd(cusolverDn_handle, itype_cu,
                                       jobz_cu, fillMode, *A_nrows,
                                       (cuDoubleComplex*)(*A_ptr), *lda,
                                       (cuDoubleComplex*)(*B_ptr), *ldb,
                                       (double*)(*W_ptr),
                                       (cuDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo) );
    }

  CUDA_API_CHECK( cudaMemcpy(info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );

  CUDA_API_CHECK( cudaFree(devInfo) );

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

  const cusolverEigType_t itype_cu = select_eigen_type(*itype);
  const cusolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDsygvd_bufferSize(cusolverDn_handle, itype_cu,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (const double*)(*A_ptr), *lda,
                                                  (const double*)(*B_ptr), *ldb,
                                                  (const double*)(*W_ptr), lwork) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZhegvd_bufferSize(cusolverDn_handle, itype_cu,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (const cuDoubleComplex*)(*A_ptr), *lda,
                                                  (const cuDoubleComplex*)(*B_ptr), *ldb,
                                                  (const double*)(*W_ptr), lwork) );
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

  const cusolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  int* devInfo;
  CUDA_API_CHECK( cudaMalloc( (void**) &devInfo, 1*sizeof(int)) );

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDsyevd(cusolverDn_handle,
                                       jobz_cu, fillMode, *A_nrows,
                                       (double*)(*A_ptr), *lda,
                                       (double*)(*W_ptr),
                                       (double*)(*work_ptr),
                                       *lwork, devInfo) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZheevd(cusolverDn_handle,
                                       jobz_cu, fillMode, *A_nrows,
                                       (cuDoubleComplex*)(*A_ptr), *lda,
                                       (double*)(*W_ptr),
                                       (cuDoubleComplex*)(*work_ptr),
                                       *lwork, devInfo) );
    }

  CUDA_API_CHECK( cudaMemcpy(info, devInfo, 1*sizeof(int), cudaMemcpyDefault) );

  CUDA_API_CHECK( cudaFree(devInfo) );

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

  const cusolverEigMode_t jobz_cu  = select_eigen_mode(jobz);
  const cublasFillMode_t  fillMode = select_cublas_fill_mode(uplo);

  if (*cplx == 1)
    {
      CUDA_API_CHECK( cusolverDnDsyevd_bufferSize(cusolverDn_handle,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (const double*)(*A_ptr), *lda,
                                                  (const double*)(*W_ptr), lwork) );
    }
  else
    {
      CUDA_API_CHECK( cusolverDnZheevd_bufferSize(cusolverDn_handle,
                                                  jobz_cu, fillMode, *A_nrows,
                                                  (const cuDoubleComplex*)(*A_ptr), *lda,
                                                  (const double*)(*W_ptr), lwork) );
    }

} // gpu_xsyevd_bufferSize_

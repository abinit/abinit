/* gpu_linalg.cu */

/*
 * Copyright (C) 2008-2019 ABINIT Group (MMancini,FDahm)
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain cublas and magma encapsulation routines,
 * that will be callable from fortran routines
 *
 */

#include <cublas.h>


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

extern "C" void gpu_linalg_init_(){
  cublasInit();
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

extern "C" void gpu_linalg_shutdown_(){
  cublasShutdown();
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

extern "C" void gpu_xgemm_(int* cplx,char *transA,char *transB,int *N,int *M,int *K,cuDoubleComplex *alpha,
			   void **A_ptr,int *lda,void** B_ptr,int *ldb,cuDoubleComplex *beta,void** C_ptr,int *ldc){
  (*cplx==1)?
    cublasDgemm(*transA,*transB,*N,*M,*K,(*alpha).x,(double *)(*A_ptr),*lda,(double *)(*B_ptr),*ldb,(*beta).x,(double *)(*C_ptr),*ldc):
    cublasZgemm(*transA,*transB,*N,*M,*K,*alpha,(cuDoubleComplex *)(*A_ptr),*lda,
		(cuDoubleComplex *)(*B_ptr),*ldb,*beta,(cuDoubleComplex *)(*C_ptr),*ldc);
}

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

extern "C" void gpu_xtrsm_(int* cplx,char *side,char *uplo,char *transA,char *diag,int *N,int *M,cuDoubleComplex *alpha,
			   void **A_ptr,int *ldA,void** B_ptr,int *ldB){
  (*cplx==1)?
    cublasDtrsm(*side,*uplo,*transA,*diag,*N,*M,(*alpha).x,(double *)(*A_ptr),*ldA,(double *)(*B_ptr),*ldB):
    cublasZtrsm(*side,*uplo,*transA,*diag,*N,*M,*alpha,(cuDoubleComplex *)(*A_ptr),*ldA,(cuDoubleComplex *)(*B_ptr),*ldB);
}

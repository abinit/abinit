// NAME
// gpu_gemm_nonlop
//
// FUNCTION
// Helper function: compute non-local operation on GPU using BLAS operators
// This is a direct equivalent of routines found in 66_nonlocal/m_gemm_nonlop.F90
//
// INPUTS
//
// PARENTS
//      m_gemm_nonlop
//
// CHILDREN
//      TODO
//
// SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>

#include "gpu_gemm_nonlop.h"

#include "gpu_gemm_nonlop_kernel.h"


#include "cuda_api_error_check.h"

extern "C" void cuda_gemm_nonlop_gpu_allocate(gemm_nonlop_gpu_t *data, int32_t npw, int32_t nprojs, int32_t istwf_k)
{
  CHECK_CUDA_ERROR( cudaMalloc( &(data->projs  ),   2*npw*nprojs*sizeof(double) ) );
  if (istwf_k > 1) {
    CHECK_CUDA_ERROR( cudaMalloc( &(data->projs_r),     npw*nprojs*sizeof(double) ) );
    CHECK_CUDA_ERROR( cudaMalloc( &(data->projs_i),     npw*nprojs*sizeof(double) ) );
  }
}

extern "C" void cuda_gemm_nonlop_gpu_deallocate(gemm_nonlop_gpu_t *data)
{
  CHECK_CUDA_ERROR( cudaFree( data->projs   ) );

  // don't worry if pointer are actually null (this might append when istwf_k <= 1)
  // cudaFree won't do anything

  CHECK_CUDA_ERROR( cudaFree( data->projs_r ) );
  CHECK_CUDA_ERROR( cudaFree( data->projs_i ) );
}

extern "C" void cuda_extract_real_part(double* data_out, const double* data_in, int size)
{

  {
    constexpr int offset = 0;
    dim3 blockSize {256,1,1};
    dim3 gridSize {(size+blockSize.x-1)/blockSize.x,1,1};
    kernel_extract_real_imag<offset><<<gridSize,blockSize>>>(data_out, data_in, size);
    GET_LAST_CUDA_ERROR("cuda_extract_real_part");
  }

} // cuda_extract_real_part

extern "C" void cuda_extract_imag_part(double* data_out, const double* data_in, int size)
{

  {
    constexpr int offset = 1;
    dim3 blockSize {256,1,1};
    dim3 gridSize {(size+blockSize.x-1)/blockSize.x,1,1};
    kernel_extract_real_imag<offset><<<gridSize,blockSize>>>(data_out, data_in, size);
    GET_LAST_CUDA_ERROR("cuda_extract_imag_part");
  }

} // cuda_extract_imag_part

extern "C" void cuda_insert_real_part(double* data_out, const double* data_in, int size)
{

  {
    constexpr int offset = 0;
    dim3 blockSize {256,1,1};
    dim3 gridSize {(size+blockSize.x-1)/blockSize.x,1,1};
    kernel_insert_real_imag<offset><<<gridSize,blockSize>>>(data_out, data_in, size);
    GET_LAST_CUDA_ERROR("cuda_insert_real_part");
  }

} // cuda_insert_real_part

extern "C" void cuda_insert_imag_part(double* data_out, const double* data_in, int size)
{

  {
    constexpr int offset = 1;
    dim3 blockSize {256,1,1};
    dim3 gridSize {(size+blockSize.x-1)/blockSize.x,1,1};
    kernel_insert_real_imag<offset><<<gridSize,blockSize>>>(data_out, data_in, size);
    GET_LAST_CUDA_ERROR("cuda_insert_imag_part");
  }

} // cuda_insert_imag_part

extern "C" void cuda_fix_realvec(double** data, int32_t npw_in, int32_t ndat_nspinor, int32_t option)
{

  {
    //constexpr int offset = 1;
    dim3 blockSize {256,1,1};
    dim3 gridSize {(ndat_nspinor+blockSize.x-1)/blockSize.x,1,1};
    kernel_fix_realvec<<<gridSize,blockSize>>>(*data, npw_in, ndat_nspinor, option);
    GET_LAST_CUDA_ERROR("cuda_fix_realvec");
  }

} // cuda_fix_realvec

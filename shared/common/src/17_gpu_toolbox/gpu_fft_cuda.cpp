/* gpu_fft_cuda.cpp */

/*
 * Copyright (C) 2008-2025 ABINIT Group
 * this file is distributed under the terms of the
 * gnu general public license, see ~abinit/COPYING
 * or http://www.gnu.org/copyleft/gpl.txt .
 * for the initials of contributors, see ~abinit/doc/developers/contributors.txt.
 *
 * The main goal of this file is to contain cublas and cufft encapsulation routines,
 * that will be callable from fortran routines
 *
 */

#include<assert.h>
#include<gpu_fft.h>
#include "stdio.h"
#include "abi_common.h"

cufftHandle plan_fft[2];
static cudaStream_t stream_compute[2];

//! utility function to select FFT type
static cufftType select_cufft_type(const int fftType_int)
{
  switch(fftType_int){
    case(CUFFT_R2C): return CUFFT_R2C;
    case(CUFFT_C2R): return CUFFT_C2R;
    case(CUFFT_C2C): return CUFFT_C2C;
    case(CUFFT_D2Z): return CUFFT_D2Z;
    case(CUFFT_Z2D): return CUFFT_Z2D;
    case(CUFFT_Z2Z): return CUFFT_Z2Z;
    default:
      fprintf(stderr, "Provided wrong enum value for FFT type: %d\n",fftType_int);
      fflush(stderr);
      abi_cabort();
      return CUFFT_R2C;
  }
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_plan_many_cpp
 *
 * FUNCTION
 *  Initialize a FFT plan with custom dimension, strided and batch size.
 *
 * INPUTS
 *   rank      Dimensionality of the transform (1, 2, or 3).
 *   n         Array of size rank, describing the size of each dimension,
 *             n[0] being the size of the outermost and n[rank-1] innermost
 *             (contiguous) dimension of a transform.
 *   inembed   Pointer of size rank that indicates the storage dimensions
 *             of the input data in memory.
 *             If set to NULL all other advanced data layout parameters are ignored.
 *   istride   Indicates the distance between two successive input elements
 *             in the least significant (i.e., innermost) dimension
 *   idist     Indicates the distance between the first element of two
 *             consecutive signals in a batch of the input data
 *   onembed   Pointer of size rank that indicates the storage dimensions of
 *             the output data in memory.
 *             If set to NULL all other advanced data layout parameters are ignored.
 *   ostride   Indicates the distance between two successive output elements in
 *             the output array in the least significant (i.e., innermost) dimension
 *   odist     Indicates the distance between the first element of two
 *             consecutive signals in a batch of the output data
 *   type      The transform data type
 *             (e.g., FFT_R2C for single precision real to complex)
 *   batch     Batch size for this transform
 */
/*=========================================================================*/

extern "C"
void gpu_fft_plan_many_cpp(int *fft_plan_id, int *rank, int **n, int **inembed,
                           int *istride, int *idist, int **onembed, int *ostride,
                           int *odist, int *fft_type, int *batch){

  assert(CUFFT_Z2Z==0x69 && "cuFFT_Type enum value mismatch !(CUDA update?)");
  assert(CUFFT_FORWARD==-1 && "cuFFT direction enum value mismatch (CUDA update?)");
  assert(CUFFT_INVERSE== 1 && "cuFFT direction enum value mismatch (CUDA update?)");

  cufftType type = select_cufft_type(*fft_type);
  CUDA_API_CHECK(cufftPlanMany(
        &plan_fft[*fft_plan_id],
        *rank,
        *n,
        *inembed,
        *istride,
        *idist,
        *onembed,
        *ostride,
        *odist,
        type,
        *batch));
  CUDA_API_CHECK( cudaStreamCreate(&stream_compute[*fft_plan_id]) );
  CUDA_API_CHECK( cufftSetStream(plan_fft[*fft_plan_id],stream_compute[*fft_plan_id]) );
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_stream_synchronize_cpp
 *
 * FUNCTION
 *  Wait for any FFT operations still running on stream
 */
/*=========================================================================*/

extern "C"
void gpu_fft_stream_synchronize_cpp(int *fft_plan_id)
{
  CUDA_API_CHECK( cudaStreamSynchronize(stream_compute[*fft_plan_id]) );
}


/*=========================================================================*/
// NAME
//  gpu_fft_plan_destroy_cpp
//
// FUNCTION
//  Destroy FFT plan
//
/*=========================================================================*/

extern "C"
void gpu_fft_plan_destroy_cpp(int *fft_plan_id){
  CUDA_API_CHECK(cufftDestroy(plan_fft[*fft_plan_id]));
  CUDA_API_CHECK(cudaStreamDestroy(stream_compute[*fft_plan_id]) );
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_exec_z2z_cpp
 *
 * FUNCTION
 *  Run a Fast Fourier Transform on double-complex input and output
 *
 * INPUTS
 *   idata       Pointer to the complex input data (in GPU memory) to transform
 *   odata       Pointer to the complex output data (in GPU memory)
 *   direction   The transform direction: FFT_FORWARD or FFT_INVERSE
 *
 * OUTPUT
 *   odata       Contains the complex Fourier coefficients
 */
/*=========================================================================*/

extern "C"
void gpu_fft_exec_z2z_cpp(int *fft_plan_id, void **idata, void **odata, int *direction){

  CUDA_API_CHECK(cufftExecZ2Z(plan_fft[*fft_plan_id], (cufftDoubleComplex*) (*idata),
                 (cufftDoubleComplex*) (*odata), *direction));
}


/*=========================================================================*/
/* NAME
 *  gpu_fft_exec_c2c_cpp
 *
 * FUNCTION
 *  Run a Fast Fourier Transform on float complex input and output
 *
 * INPUTS
 *   idata       Pointer to the complex input data (in GPU memory) to transform
 *   odata       Pointer to the complex output data (in GPU memory)
 *   direction   The transform direction: FFT_FORWARD or FFT_INVERSE
 *
 * OUTPUT
 *   odata       Contains the complex Fourier coefficients
 */
/*=========================================================================*/

extern "C"
void gpu_fft_exec_c2c_cpp(int *fft_plan_id, void **idata, void **odata, int *direction){
  CUDA_API_CHECK(cufftExecC2C(plan_fft[*fft_plan_id], (cufftComplex*) *idata,
                 (cufftComplex*) *odata, *direction));
}


#if defined HAVE_GPU_CUDA
/**
 * @struct gpu_context_t
 * @brief Opaque GPU execution context.
 *
 * This structure groups all GPU resources required for FFT-based
 * operations:
 *   - cuFFT plan for batched 3D complex-to-complex transforms,
 *   - CUDA stream used for all operations,
 *   - cuBLAS handle associated with the same stream.
 *
 * The structure is opaque to Fortran and accessed only through
 * a `void *` handle.
 */

typedef struct {
  cufftHandle     fft_plan;
  cudaStream_t    stream;
  cublasHandle_t  cublas_handle;
} gpu_context_t;

/**
 * @brief Convert Fortran FFT sign to cuFFT direction.
 *
 * @param[in]  isign     FFT sign convention (Fortran-style).
 *                       +1 : inverse FFT
 *                       -1 : forward FFT
 * @param[out] direction Corresponding cuFFT direction constant.
 *
 * @abort On invalid value of isign.
 */

static void
_get_direction(int isign, int *direction){

  switch (isign) {
  case 1:
    *direction = CUFFT_INVERSE;
    break;
  case -1:
    *direction = CUFFT_FORWARD;
    break;
  default:
    printf("Invalid isign: %d\n", isign);
    abi_cabort();
  }
}

/**
 * @brief Determine cuFFT type and buffer size from precision.
 *
 * @param[in]  dist_batch Total number of complex elements.
 * @param[in]  kind       Fortran kind (4 = single, 8 = double).
 * @param[out] type       cuFFT transform type (C2C or Z2Z).
 * @param[out] nbytes     Required buffer size in bytes.
 *
 * @abort On unsupported kind.
 */

static void
_get_type_nbytes(int dist_batch, int kind, cufftType *type, size_t *nbytes){

  switch (kind) {
  case 4:
    *type = CUFFT_C2C;
    *nbytes = dist_batch * sizeof(cufftComplex);
    break;
  case 8:
    *type = CUFFT_Z2Z;
    *nbytes = dist_batch * sizeof(cufftDoubleComplex);
    break;
  default:
    printf("Invalid kind: %d\n", kind);
    abi_cabort();
  }
}

/**
 * @brief Initialize a GPU FFT execution context.
 *
 * This routine:
 *   - allocates a new GPU context,
 *   - creates a CUDA stream,
 *   - builds a batched 3D cuFFT plan,
 *   - creates a cuBLAS handle bound to the same stream.
 *
 * Array dimensions are passed from Fortran and converted to C row-major order.
 *
 * @param[out] void_ctx Opaque GPU context pointer (Fortran c_ptr).
 * @param[in]  f_dims   FFT dimensions (Fortran order, size 3).
 * @param[in]  f_embed  Embedded dimensions (Fortran order, size 3).
 * @param[in]  batch    Number of FFT batches.
 * @param[in]  kind     Fortran kind (4 = single precision, 8 = double).
 *
 * @note Ownership of the context is transferred to the caller.
 *       The context must be released with gpu_ctx_free().
 */

extern "C" void
gpu_ctx_init_cpp(void **void_ctx, int *f_dims, int *f_embed, int batch, int kind) {

  const int RANK3 = 3, stride1 = 1;
  size_t nbytes;
  int c_dims[RANK3], c_embed[RANK3];
  // Fortran to C
  c_dims[0] = f_dims[2]; c_dims[1] = f_dims[1]; c_dims[2] = f_dims[0];
  c_embed[0] = f_embed[2]; c_embed[1] = f_embed[1]; c_embed[2] = f_embed[0];
  int dist = f_embed[0] * f_embed[1] * f_embed[2];

  //printf("in C gpu_fftbox_plan_init\n");

  //cufftResult = cufftPlanMany(cufftHandle *plan, int rank, int *c_dims,
  //                            int *inembed, int istride, int idist,
  //                            int *onembed, int ostride, int odist,
  //                            cufftType type, int batch);

  cufftType type;
  _get_type_nbytes(dist * batch, kind, &type, &nbytes);

  gpu_context_t *ctx = (gpu_context_t *) calloc(1, sizeof(*ctx));

  CUDA_API_CHECK(cufftPlanMany(&ctx->fft_plan, RANK3, c_dims, c_embed, stride1, dist, c_embed, stride1, dist, type, batch));
  //printf("Creating new GPU plan_p: %d @ %p\n", plan, *plan_p);

  /* Associate plan with stream */
  CUDA_API_CHECK(cudaStreamCreate(&ctx->stream));
  CUDA_API_CHECK(cufftSetStream(ctx->fft_plan, ctx->stream));

  /* cuBLAS */
  CUDA_API_CHECK(cublasCreate(&ctx->cublas_handle));
  CUDA_API_CHECK(cublasSetStream(ctx->cublas_handle, ctx->stream));

  // Return void pointers to Fortran
  *void_ctx = (void *) ctx;
}

/**
 * @brief Synchronize the GPU stream associated with a context.
 *
 * Blocks until all previously issued operations on the stream
 * have completed.
 *
 * @param[in] void_ctx Opaque GPU context pointer.
 */

extern "C" void
gpu_ctx_synch_cpp(void *void_ctx) {

  gpu_context_t *ctx = (gpu_context_t *) void_ctx;
  CUDA_API_CHECK(cudaStreamSynchronize(ctx->stream));
}

/**
 * @brief Destroy a GPU context and release all associated resources.
 *
 * This routine:
 *   - destroys the cuFFT plan,
 *   - destroys the cuBLAS handle,
 *   - destroys the CUDA stream,
 *   - frees the context structure,
 *   - sets the caller pointer to NULL.
 *
 * Safe to call multiple times.
 *
 * @param[in,out] void_ctx Pointer to opaque GPU context handle.
 */

extern "C" void
gpu_ctx_free_cpp(void **void_ctx)
{
  //printf("In gpu_ctx_free_cpp\n");
  if (!void_ctx || !*void_ctx) return;

  gpu_context_t *ctx = (gpu_context_t *)(*void_ctx);

  cufftDestroy(ctx->fft_plan);
  cublasDestroy(ctx->cublas_handle);
  cudaStreamDestroy(ctx->stream);

  free(ctx);
  *void_ctx = NULL;
}

/**
 * @brief In-place complex-to-complex FFT on the GPU.
 *
 * Executes a batched 3D FFT in-place using the context's cuFFT plan.
 * Optional scaling is applied for forward transforms using cuBLAS.
 *
 * @param[in]     void_ctx Opaque GPU context pointer.
 * @param[in]     nfft     Number of grid points per FFT.
 * @param[in]     batch   Number of FFT batches.
 * @param[in]     isign   FFT sign (+1 inverse, -1 forward).
 * @param[in]     iscale  Apply scaling if non-zero.
 * @param[in]     kind    Precision kind (4 = single, 8 = double).
 * @param[in,out] d_ff    Device pointer to input/output data.
 */

extern "C" void
gpu_fftbox_c2c_ip_cpp(void *void_ctx, int nfft, int batch, int isign, int iscale, int kind, void **d_ff) {

  //printf("in gpu_fftbox_c2c_ip_cpp");
  cufftType type;
  int direction;
  size_t nbytes;
  _get_direction(isign, &direction);
  _get_type_nbytes(0, kind, &type, &nbytes);

  gpu_context_t *ctx = (gpu_context_t *) void_ctx;

  // Transform the signal in place.
  if (type == CUFFT_C2C) {
    CUDA_API_CHECK(cufftExecC2C(ctx->fft_plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CUDA_API_CHECK(cublasCsscal(ctx->cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_ff, 1));
    }
  }

  if (type == CUFFT_Z2Z) {
    CUDA_API_CHECK(cufftExecZ2Z(ctx->fft_plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CUDA_API_CHECK(cublasZdscal(ctx->cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_ff, 1));
    }
  }

  //CUDA_API_CHECK(cudaStreamSynchronize(ctx->stream));
}

/**
 * @brief Out-of-place complex-to-complex FFT on the GPU.
 *
 * Executes a batched 3D FFT using distinct input and output buffers.
 * Optional scaling is applied for forward transforms using cuBLAS.
 *
 * @param[in]  void_ctx Opaque GPU context pointer.
 * @param[in]  nfft     Number of grid points per FFT.
 * @param[in]  batch   Number of FFT batches.
 * @param[in]  isign   FFT sign (+1 inverse, -1 forward).
 * @param[in]  iscale  Apply scaling if non-zero.
 * @param[in]  kind    Precision kind (4 = single, 8 = double).
 * @param[in]  d_ff    Device pointer to input data.
 * @param[out] d_gg    Device pointer to output data.
 */

extern "C" void
gpu_fftbox_c2c_op_cpp(void *void_ctx, int nfft, int batch, int isign, int iscale, int kind,
                  void **d_ff, void **d_gg) {

  //printf("in gpu_fftbox_c2c_op_cpp");
  cufftType type;
  int direction;
  size_t nbytes;
  _get_direction(isign, &direction);
  _get_type_nbytes(0, kind, &type, &nbytes);

  gpu_context_t *ctx = (gpu_context_t *) void_ctx;

  // Transform the signal out of place.
  if (type == CUFFT_C2C) {
    CUDA_API_CHECK(cufftExecC2C(ctx->fft_plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_gg, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CUDA_API_CHECK(cublasCsscal(ctx->cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_gg, 1));
    }
  }
  if (type == CUFFT_Z2Z) {
    CUDA_API_CHECK(cufftExecZ2Z(ctx->fft_plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_gg, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CUDA_API_CHECK(cublasZdscal(ctx->cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_gg, 1));
    }
  }

  //CUDA_API_CHECK(cudaStreamSynchronize(ctx->stream));
}

#endif


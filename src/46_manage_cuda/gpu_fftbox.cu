/**
 * @file gpu_fftbox.cu
 * @brief CUDA GPU context and FFT execution utilities callable from Fortran.
 *
 * This file provides a minimal GPU execution context encapsulating:
 *   - a CUDA stream,
 *   - a cuFFT plan,
 *   - a cuBLAS handle bound to the same stream.
 *
 * The API is designed to be called from Fortran via `ISO_C_BINDING`,
 * using opaque `type(c_ptr)` handles on the Fortran side.
 *
 * All GPU resources are owned by the context and must be explicitly
 * released via gpu_ctx_free().
 *
 * This code is compiled only when HAVE_GPU_CUDA is enabled.
 */

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"

#if defined HAVE_GPU_CUDA

#include "cuda_header.h"
#include <gpu_linalg.h>
#include "cuda_api_error_check.h"


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

  CHECK_CUDA_ERROR(cufftPlanMany(&ctx->fft_plan, RANK3, c_dims, c_embed, stride1, dist, c_embed, stride1, dist, type, batch));
  //printf("Creating new GPU plan_p: %d @ %p\n", plan, *plan_p);

  /* Associate plan with stream */
  CHECK_CUDA_ERROR(cudaStreamCreate(&ctx->stream));
  CHECK_CUDA_ERROR(cufftSetStream(ctx->fft_plan, ctx->stream));

  /* cuBLAS */
  CHECK_CUDA_ERROR(cublasCreate(&ctx->cublas_handle));
  CHECK_CUDA_ERROR(cublasSetStream(ctx->cublas_handle, ctx->stream));

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
  CHECK_CUDA_ERROR(cudaStreamSynchronize(ctx->stream));
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
    CHECK_CUDA_ERROR(cufftExecC2C(ctx->fft_plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CHECK_CUDA_ERROR(cublasCsscal(ctx->cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_ff, 1));
    }
  }

  if (type == CUFFT_Z2Z) {
    CHECK_CUDA_ERROR(cufftExecZ2Z(ctx->fft_plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CHECK_CUDA_ERROR(cublasZdscal(ctx->cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_ff, 1));
    }
  }

  //CHECK_CUDA_ERROR(cudaStreamSynchronize(ctx->stream));
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
    CHECK_CUDA_ERROR(cufftExecC2C(ctx->fft_plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_gg, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CHECK_CUDA_ERROR(cublasCsscal(ctx->cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_gg, 1));
    }
  }
  if (type == CUFFT_Z2Z) {
    CHECK_CUDA_ERROR(cufftExecZ2Z(ctx->fft_plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_gg, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CHECK_CUDA_ERROR(cublasZdscal(ctx->cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_gg, 1));
    }
  }

  //CHECK_CUDA_ERROR(cudaStreamSynchronize(ctx->stream));
}

#endif

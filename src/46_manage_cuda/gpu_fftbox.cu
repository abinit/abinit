
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"

#if defined HAVE_GPU_CUDA

#include "cuda_header.h"
#include <gpu_linalg.h>
#include "cuda_api_error_check.h"


typedef struct {
  cufftHandle     fft_plan;
  cudaStream_t    stream;
  cublasHandle_t  cublas_handle;
} gpu_context_t;


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

extern "C" void
gpu_ctx_init(void **void_ctx, int *f_dims, int *f_embed, int batch, int kind) {

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


extern "C" void
gpu_ctx_free(void **void_ctx)
{
  //printf("In gpu_ctx_free\n");
  if (!void_ctx || !*void_ctx) return;

  gpu_context_t *ctx = (gpu_context_t *)(*void_ctx);

  cufftDestroy(ctx->fft_plan);
  cublasDestroy(ctx->cublas_handle);
  cudaStreamDestroy(ctx->stream);

  free(ctx);
  *void_ctx = NULL;
}


extern "C" void
gpu_fftbox_c2c_ip(void *void_cxt, int nfft, int batch, int isign, int iscale, int kind, void **d_ff) {

  //printf("in gpu_fftbox_c2c_ip");
  cufftType type;
  int direction;
  size_t nbytes;
  _get_direction(isign, &direction);
  _get_type_nbytes(0, kind, &type, &nbytes);

  gpu_context_t *ctx = (gpu_context_t *) void_cxt;

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

  CHECK_CUDA_ERROR(cudaStreamSynchronize(ctx->stream));
}


extern "C" void
gpu_fftbox_c2c_op(void *void_ctx, int nfft, int batch, int isign, int iscale, int kind,
                  void **d_ff, void **d_gg) {

  //printf("in gpu_fftbox_c2c_op");
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

  CHECK_CUDA_ERROR(cudaStreamSynchronize(ctx->stream));
}

#endif

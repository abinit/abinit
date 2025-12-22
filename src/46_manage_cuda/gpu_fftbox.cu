
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"

#if defined HAVE_GPU_CUDA

#include "cuda_header.h"
#include <gpu_linalg.h>
#include "cuda_api_error_check.h"

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
gpu_fft_plan_init(void **plan_pp, void **stream_pp, int *f_dims, int *f_embed, int batch, int kind) {

  const int RANK3 = 3, stride1 = 1;
  size_t nbytes;
  int c_dims[RANK3], c_embed[RANK3];
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

  // Allocate the handle
  cufftHandle *plan_p = (cufftHandle *) malloc(sizeof(cufftHandle));

  CHECK_CUDA_ERROR(cufftPlanMany(plan_p, RANK3, c_dims, c_embed, stride1, dist, c_embed, stride1, dist, type, batch));
  //printf("Creating new GPU plan_p: %d @ %p\n", plan, *plan_p);
  //printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);

  /* Associate plan with stream */
  // see gpu_fft.cu
  cudaStream_t *fft_stream = (cudaStream_t *) malloc(sizeof(cudaStream_t));
  CHECK_CUDA_ERROR(cudaStreamCreate(fft_stream));
  CHECK_CUDA_ERROR(cufftSetStream(*plan_p, *fft_stream));

  // Return void pointers to Fortran
  *plan_pp = (void *) plan_p;
  *stream_pp = (void *) fft_stream;
}


extern "C" void
gpu_fft_plan_free(void *void_ptr)
{
  cufftHandle *plan = (cufftHandle *) void_ptr;
  //printf("In gpu_fft_plan_free. About to free GPU plan: %d @ %p\n", *plan, plan);

  if (plan) {
    CHECK_CUDA_ERROR(cufftDestroy(*plan));
    free(plan);
  }
}


extern "C" void
gpu_stream_free(void *void_ptr)
{
  cudaStream_t *stream =  (cudaStream_t *) void_ptr;
  //printf("In gpu_stream_free. About to free GPU stream: %d @ %p\n", *stream, stream);

  if (stream) {
    CHECK_CUDA_ERROR(cudaStreamDestroy(*stream));
    free(stream);
  }
}


extern "C" void
gpu_fftbox_c2c_ip(void **plan_pp, void *stream, int nfft, int batch, int isign, int iscale, int kind, void **d_ff) {

  //printf("in gpu_fftbox_c2c_ip");
  cufftType type;
  int direction;
  size_t nbytes;
  _get_direction(isign, &direction);
  _get_type_nbytes(0, kind, &type, &nbytes);

  cufftHandle plan = *(cufftHandle *) (*plan_pp);

  // Transform the signal in place.
  if (type == CUFFT_C2C) {
    CHECK_CUDA_ERROR(cufftExecC2C(plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CHECK_CUDA_ERROR(cublasCsscal(cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_ff, 1));
    }
  }

  if (type == CUFFT_Z2Z) {
    CHECK_CUDA_ERROR(cufftExecZ2Z(plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CHECK_CUDA_ERROR(cublasZdscal(cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_ff, 1));
    }
  }

  cudaStream_t *fft_stream = (cudaStream_t *) stream;
  CHECK_CUDA_ERROR(cudaStreamSynchronize(*fft_stream));
}


extern "C" void
gpu_fftbox_c2c_op(void **plan_pp, void *stream, int nfft, int batch, int isign, int iscale, int kind,
                  void **d_ff, void **d_gg) {

  //printf("in gpu_fftbox_c2c_op");

  cufftType type;
  int direction;
  size_t nbytes;
  _get_direction(isign, &direction);
  _get_type_nbytes(0, kind, &type, &nbytes);

  cufftHandle plan = *(cufftHandle *) (*plan_pp);

  // Transform the signal out of place.
  if (type == CUFFT_C2C) {
     CHECK_CUDA_ERROR(cufftExecC2C(plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_gg, direction));
     if (direction == CUFFT_FORWARD and iscale != 0){
         float alpha_sp = 1.0f / nfft;
         CHECK_CUDA_ERROR(cublasCsscal(cublas_handle, nfft*batch, &alpha_sp, (cuComplex *) *d_gg, 1));
     }
  }
  if (type == CUFFT_Z2Z) {
     CHECK_CUDA_ERROR(cufftExecZ2Z(plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_gg, direction));
     if (direction == CUFFT_FORWARD and iscale != 0){
         double alpha_dp = 1.0 / nfft;
         CHECK_CUDA_ERROR(cublasZdscal(cublas_handle, nfft*batch, &alpha_dp, (cuDoubleComplex *) *d_gg, 1));
     }
  }

  cudaStream_t *fft_stream = (cudaStream_t *) stream;
  CHECK_CUDA_ERROR(cudaStreamSynchronize(*fft_stream));
}

#endif

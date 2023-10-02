//#if defined HAVE_CONFIG_H
#include "config.h"
//#endif

#include "abi_common.h"
#include "stdio.h"
#if defined HAVE_GPU_CUDA
#include "cuda_header.h"
#include <gpu_linalg.h>
#include "cuda_api_error_check.h"
#endif

static void
_cuget_type_nbytes_dir_cc(int dist_ndat, int kind, int isign,
                          cufftType *type, size_t *nbytes, int *direction){

  switch (kind) {
  case 4:
    *type = CUFFT_C2C;
    *nbytes = dist_ndat * sizeof(cufftComplex);
    break;
  case 8:
    *type = CUFFT_Z2Z;
    *nbytes = dist_ndat * sizeof(cufftDoubleComplex);
    break;
  default:
    printf("Invalid kind: %d\n", kind);
    abi_cabort();
  }

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

extern "C" void
gpu_planpp_free(void **plan_pp) {
  if (*plan_pp == NULL) return;
  //cufftHandle plan = * ((cufftHandle *) (*plan_pp));
  //printf("About to free GPU plan: %d @ %p\n", plan, &plan);
  //printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);
  //CHECK_CUDA_ERROR(cufftDestroy(plan));
  *plan_pp = NULL;
}


extern "C" void
devpp_free(void **dev_pp) {
  //printf("About to free devptr %p\n", *dev_pp);
  if (*dev_pp == NULL) return;
  CHECK_CUDA_ERROR(cudaFree(*dev_pp));
  *dev_pp = NULL;
}

extern "C" void
xgpu_fftbox_c2c_ip(int *f_dims, int *f_embed, int ndat, int isign, int kind, int iscale,
                   void **h_ff, void **plan_pp, void **d_ff) {

  const int RANK = 3, stride = 1;
  size_t nbytes;
  int c_dims[RANK], c_embed[RANK], direction;
  c_dims[0] = f_dims[2]; c_dims[1] = f_dims[1]; c_dims[2] = f_dims[0];
  c_embed[0] = f_embed[2]; c_embed[1] = f_embed[1]; c_embed[2] = f_embed[0];
  int dist = f_embed[0] * f_embed[1] * f_embed[2];
  int nfft = c_dims[0] * c_dims[1] * c_dims[2];

#if defined HAVE_GPU_CUDA
  cufftType type;
  cufftHandle plan;
  _cuget_type_nbytes_dir_cc(dist * ndat, kind, isign, &type, &nbytes, &direction);

  if (*d_ff == NULL) {
    //printf("Calling cudaMalloc\n");
    //CHECK_CUDA_ERROR(cudaMalloc((void**) &d_ff, nbytes));
    CHECK_CUDA_ERROR(cudaMalloc(d_ff, nbytes));}
  else {
    //printf("Reusing devptr %p\n", *d_ff);
  }

  CHECK_CUDA_ERROR(cudaMemcpy(*d_ff, *h_ff, nbytes, cudaMemcpyHostToDevice));

  /* Create a 3D FFT plan.
  cufftResult = cufftPlanMany(cufftHandle *plan, int rank, int *c_dims,
                              int *inembed, int istride, int idist,
                              int *onembed, int ostride, int odist,
                              cufftType type, int batch);
  */

  //printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);

  if (*plan_pp == NULL) {
    CHECK_CUDA_ERROR(cufftPlanMany(&plan, RANK, c_dims, c_embed, stride, dist, c_embed, stride, dist, type, ndat));
    *plan_pp = &plan;
    //printf("Creating new GPU plan: %d @ %p\n", plan, &plan);
  }
  else {
    plan = * ((cufftHandle *) (*plan_pp));
    //printf("Reusing GPU plan: %d @ %p\n", plan, &plan);
  }

  /* Transform the signal in place. */
  if (type == CUFFT_C2C) {
    CHECK_CUDA_ERROR(cufftExecC2C(plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        float alpha_sp = 1.0f / nfft;
        CHECK_CUDA_ERROR(cublasCsscal(cublas_handle, dist*ndat, &alpha_sp, (cuComplex *) *d_ff, 1));
    }
  }
  if (type == CUFFT_Z2Z) {
    CHECK_CUDA_ERROR(cufftExecZ2Z(plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_ff, direction));
    if (direction == CUFFT_FORWARD and iscale != 0){
        double alpha_dp = 1.0 / nfft;
        CHECK_CUDA_ERROR(cublasZdscal(cublas_handle, dist*ndat, &alpha_dp, (cuDoubleComplex *) *d_ff, 1));
    }
  }

  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(*h_ff, *d_ff, nbytes, cudaMemcpyDeviceToHost));

  //CHECK_CUDA_ERROR(cudaFree(*d_ff)); *d_ff = NULL;
  CHECK_CUDA_ERROR(cufftDestroy(plan)); *plan_pp = NULL;
#endif
}


extern "C" void
xgpu_fftbox_c2c_op(int *f_dims, int *f_embed, int ndat, int isign, int kind, int iscale,
                   void **h_ff, void **h_gg, void **plan_pp, void **d_ff, void **d_gg) {

  const int RANK = 3, stride = 1;
  size_t nbytes;
  int c_dims[RANK], c_embed[RANK], direction;
  c_dims[0] = f_dims[2]; c_dims[1] = f_dims[1]; c_dims[2] = f_dims[0];
  c_embed[0] = f_embed[2]; c_embed[1] = f_embed[1]; c_embed[2] = f_embed[0];
  int dist = f_embed[0] * f_embed[1] * f_embed[2];
  int nfft = c_dims[0] * c_dims[1] * c_dims[2];

#if defined HAVE_GPU_CUDA
  cufftType type;
  cufftHandle plan;
  _cuget_type_nbytes_dir_cc(dist * ndat, kind, isign, &type, &nbytes, &direction);

  if (*d_ff == NULL) {
    //printf("Calling cudaMalloc for d_ff\n");
    CHECK_CUDA_ERROR(cudaMalloc(d_ff, nbytes));}
  else {
    //printf("Reusing *d_ff devptr %p\n", *d_ff);
  }
  CHECK_CUDA_ERROR(cudaMemcpy(*d_ff, *h_ff, nbytes, cudaMemcpyHostToDevice));

  if (*d_gg == NULL) {
    //printf("Calling cudaMalloc for d_gg\n");
    CHECK_CUDA_ERROR(cudaMalloc(d_gg, nbytes));}
  else {
    //printf("Reusing *d_gg devptr %p\n", *d_gg);
  }
  CHECK_CUDA_ERROR(cudaMemcpy(*d_gg, *h_gg, nbytes, cudaMemcpyHostToDevice));

  /* Create a 3D FFT plan. */
  //printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);

  if (*plan_pp == NULL) {
    CHECK_CUDA_ERROR(cufftPlanMany(&plan, RANK, c_dims, c_embed, stride, dist, c_embed, stride, dist, type, ndat));
    *plan_pp = &plan;
    //printf("Creating new GPU plan: %d @ %p\n", plan, &plan);
  }
  else {
    plan = * ((cufftHandle *) (*plan_pp));
    //printf("Reusing GPU plan: %d @ %p\n", plan, &plan);
  }

  /* Transform the signal out of place. */
  if (type == CUFFT_C2C) {
     CHECK_CUDA_ERROR(cufftExecC2C(plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_gg, direction));
     if (direction == CUFFT_FORWARD and iscale != 0){
         float alpha_sp = 1.0f / nfft;
         CHECK_CUDA_ERROR(cublasCsscal(cublas_handle, dist*ndat, &alpha_sp, (cuComplex *) *d_gg, 1));
     }
  }
  if (type == CUFFT_Z2Z) {
     CHECK_CUDA_ERROR(cufftExecZ2Z(plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_gg, direction));
     if (direction == CUFFT_FORWARD and iscale != 0){
         double alpha_dp = 1.0 / nfft;
         CHECK_CUDA_ERROR(cublasZdscal(cublas_handle, dist*ndat, &alpha_dp, (cuDoubleComplex *) *d_gg, 1));
     }
  }

  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(*h_gg, *d_gg, nbytes, cudaMemcpyDeviceToHost));

  //CHECK_CUDA_ERROR(cudaFree(*d_ff)); *d_ff = NULL;
  //CHECK_CUDA_ERROR(cudaFree(*d_gg)); *d_gg = NULL;
  CHECK_CUDA_ERROR(cufftDestroy(plan)); *plan_pp = NULL;
#endif
}

//#if defined HAVE_CONFIG_H
#include "config.h"
//#endif

#include "abi_common.h"
#include "stdio.h"
#if defined HAVE_GPU_CUDA
//#include "abi_gpu_header.h"
#include "cuda_header.h"
#include "cuda_api_error_check.h"
#endif


extern "C" void gpu_planpp_free(void **plan_pp) {
  if (*plan_pp == NULL) return;
  cufftHandle plan = * ((cufftHandle *) (*plan_pp));
  //cufftHandle plan = * ((cufftHandle *) (plan_pp));
  printf("About to free GPU plan: %d @ %p\n", plan, &plan);
  printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);
  //CHECK_CUDA_ERROR(cufftDestroy(plan));
  *plan_pp = NULL;
}


extern "C" void devpp_free(void **dev_pp) {
  printf("About to free devptr %p\n", *dev_pp);
  if (*dev_pp == NULL) return;
  CHECK_CUDA_ERROR(cudaFree(*dev_pp));
  *dev_pp = NULL;
}

// For scaling see:
// https://github.com/zchee/cuda-sample/blob/master/7_CUDALibraries/simpleCUFFT_callback/simpleCUFFT_callback.cu

extern "C" void xgpu_fftbox_c2c_ip(int *f_dims, int *f_embed, int ndat, int isign, int kind,
                                   void **h_ff, void **plan_pp, void **d_ff) {
                                   //cufftComplex *h_ff, cufftHandle *plan_pp, cufftComplex *d_ff) {

  const int RANK = 3, stride = 1;
  size_t nbytes;
  int c_dims[RANK], c_embed[RANK], dist, direction;
  c_dims[0] = f_dims[2]; c_dims[1] = f_dims[1]; c_dims[2] = f_dims[0];
  c_embed[0] = f_embed[2]; c_embed[1] = f_embed[1]; c_embed[2] = f_embed[0];
  dist = f_embed[0] * f_embed[1] * f_embed[2];

#if defined HAVE_GPU_CUDA
  cufftType type;
  cufftHandle plan;
  switch (kind) {
  case 4:
    type = CUFFT_C2C;
    nbytes = dist * ndat * sizeof(cufftComplex);
    break;
  case 8:
    type = CUFFT_Z2Z;
    nbytes = dist * ndat * sizeof(cufftDoubleComplex);
    break;
  default:
    printf("Invalid kind: %d\n", kind);
    abi_cabort();
  }

  switch (isign) {
  case 1:
    direction = CUFFT_INVERSE;
    break;
  case -1:
    direction = CUFFT_FORWARD;
    break;
  default:
    printf("Invalid isign: %d\n", isign);
    abi_cabort();
  }

  if (*d_ff == NULL) {
    printf("Calling cudaMalloc\n");
    //CHECK_CUDA_ERROR(cudaMalloc((void**) &d_ff, nbytes));
    CHECK_CUDA_ERROR(cudaMalloc(d_ff, nbytes));}
  else {
    printf("Reusing devptr %p\n", *d_ff);
  }

  CHECK_CUDA_ERROR(cudaMemcpy(*d_ff, *h_ff, nbytes, cudaMemcpyHostToDevice));

  /* Create a 3D FFT plan.
  cufftResult = cufftPlanMany(cufftHandle *plan, int rank, int *c_dims,
                              int *inembed, int istride, int idist,
                              int *onembed, int ostride, int odist,
                              cufftType type, int batch);
  */

  printf("plan_pp: %p, *plan_pp: %p\n", plan_pp, *plan_pp);

  if (*plan_pp == NULL) {
    CHECK_CUDA_ERROR(cufftPlanMany(&plan, RANK, c_dims, c_embed, stride, dist, c_embed, stride, dist, type, ndat));
    *plan_pp = &plan;
    printf("Creating new GPU plan: %d @ %p\n", plan, &plan);
  }
  else {
    plan = * ((cufftHandle *) (*plan_pp));
    printf("Reusing GPU plan: %d @ %p\n", plan, &plan);
  }

  /* Transform the signal in place. */
  if (type == CUFFT_C2C) {
		CHECK_CUDA_ERROR(cufftExecC2C(plan, (cufftComplex *) *d_ff, (cufftComplex *) *d_ff, direction));
	}
  if (type == CUFFT_Z2Z) {
		CHECK_CUDA_ERROR(cufftExecZ2Z(plan, (cufftDoubleComplex *) *d_ff, (cufftDoubleComplex *) *d_ff, direction));
	}

  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(*h_ff, *d_ff, nbytes, cudaMemcpyDeviceToHost));

  //CHECK_CUDA_ERROR(cudaFree(*d_ff));
  //*d_ff = NULL;
  //CHECK_CUDA_ERROR(cufftDestroy(plan));
  //*plan_pp = NULL;
#endif
}

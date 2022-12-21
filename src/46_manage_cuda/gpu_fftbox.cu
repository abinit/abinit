#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"
#include "stdio.h"
#if defined HAVE_GPU_CUDA
//#include "abi_gpu_header.h"
#include "cuda_header.h" 
#include "cuda_api_error_check.h"
#endif


//extern "C" int  _isign2cuda(int isign){
//  switch (isign){
//  case 1:
//    return cuFFTINVERSE;
//  case -1:
//    return cuFFTFORWARD;
//  default
//    printf("Invalid isign: %d\n", *isign);
//    abi_cabort();
//}


extern "C" void xgpu_fftbox_c2c_ip(int *f_dims, int *f_embed, int ndat, int isign, int kind, 
                                   void *h_ff, void *plan_ptr, void *d_ff) {
                                   //cufftComplex *h_ff, cufftHandle *plan_ptr, cufftComplex *d_ff) {

  const int RANK = 3, stride = 1;
  size_t nbytes;
  int c_dims[RANK], c_embed[RANK], dist, type;
  c_dims[0] = f_dims[2]; c_dims[1] = f_dims[1]; c_dims[2] = f_dims[0];
  c_embed[0] = f_embed[2]; c_embed[1] = f_embed[1]; c_embed[2] = f_embed[0];
  dist = f_embed[0] * f_embed[1] * f_embed[2];

#if defined HAVE_GPU_CUDA
  nbytes = dist * ndat * sizeof(cufftComplex);
  //switch (kind){
  //  case 4:
  //    type = CUFFT_C2C;
  //    break;
  //  case 8:
  //    type = CUFFT_Z2Z;
  //    break;
  //  default
  //    printf("Invalid kind: %d\n", kind);
  //    abi_cabort();
  //}

  type = CUFFT_C2C;
  //type = CUFFT_Z2Z;

  //cufftComplex *d_ff; cufftHandle plan_ptr;
  //
  if (d_ff == NULL){ 
    CHECK_CUDA_ERROR(cudaMalloc((void**) &d_ff, nbytes));
  }

  CHECK_CUDA_ERROR(cudaMemcpy(d_ff, h_ff, nbytes, cudaMemcpyHostToDevice));

  /* Create a 3D FFT plan. 
  cufftResult = cufftPlanMany(cufftHandle *plan, int rank, int *c_dims, 
                              int *inembed, int istride, int idist, 
                              int *onembed, int ostride, int odist, 
                              cufftType type, int batch);
  */

  if (plan_ptr == NULL){
      CHECK_CUDA_ERROR(cufftPlanMany(plan_ptr, RANK, c_dims, c_embed, stride, dist, c_embed, stride, dist, type, ndat));
  }

  /* Transform the signal in place. */
  CHECK_CUDA_ERROR(cufftExecC2C(*plan_ptr, d_ff, d_ff, CUFFT_FORWARD));
  CHECK_CUDA_ERROR(cudaDeviceSynchronize());
  CHECK_CUDA_ERROR(cudaMemcpy(h_ff, d_ff, nbytes, cudaMemcpyDeviceToHost));

  //CHECK_CUDA_ERROR(cudaFree(d_ff));
  //CHECK_CUDA_ERROR(cufftDestroy(*plan_ptr));
#endif
}

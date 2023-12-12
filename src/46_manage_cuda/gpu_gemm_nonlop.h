#ifndef ABI_MANAGE_CUDA_GEMM_NONLOP_H
#define ABI_MANAGE_CUDA_GEMM_NONLOP_H

#include <stdint.h>

#ifdef  __cplusplus
extern "C" {
#endif

  //! type structure used for cuda/c/fortran interoperability (iso_c_bindings)
  //! this struct must be directly compatible (in the iso_c_bindings meaning) with
  //! fortran derived type gemm_nonlop_gpu_type
  typedef struct gemm_nonlop_gpu_t
  {

    int32_t npw;
    int32_t nprojs;

    double* projs;
    double* projs_r;
    double* projs_i;

  } gemm_nonlop_gpu_t;

#ifdef  __cplusplus
}
#endif

#endif /* ABI_MANAGE_CUDA_GEMM_NONLOP_H */

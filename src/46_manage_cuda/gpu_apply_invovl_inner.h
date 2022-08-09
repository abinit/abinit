#ifndef ABI_MANAGE_CUDA_APPLY_INVOVL_INNER_H
#define ABI_MANAGE_CUDA_APPLY_INVOVL_INNER_H

#include <stdint.h>

#ifdef  __cplusplus
extern "C" {
#endif

  //! type structure used for cuda/c/fortran interoperability (iso_c_bindings)
  //! this struct must be directly compatible (in the iso_c_bindings meaning) with
  //! fortran derived type invovl_kpt_gpu_type
  typedef struct invovl_kpt_gpu_t
  {

    int32_t nprojs;

    double* gram_projs;
    int32_t gram_projs_dim[3];

    double* inv_sij;
    int32_t inv_sij_dim[4];

    double* inv_s_approx;
    int32_t inv_s_approx_dim[4];

  } invovl_kpt_gpu_t;

#ifdef  __cplusplus
}
#endif

#endif /* ABI_MANAGE_CUDA_APPLY_INVOVL_INNER_H */

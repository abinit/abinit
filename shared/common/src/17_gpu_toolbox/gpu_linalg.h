#ifndef ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H
#define ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H

#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_api_error_check.h>

extern cublasHandle_t cublas_handle;
extern cusolverDnHandle_t cusolverDn_handle;

#endif // ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H

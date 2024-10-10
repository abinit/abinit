#ifndef ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H
#define ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#ifdef HAVE_GPU_CUDA
#include <cublas_v2.h>
#include <cusolverDn.h>
#include <cuda_api_error_check.h>

extern cublasHandle_t cublas_handle;
extern cusolverDnHandle_t cusolverDn_handle;
#endif

#ifdef HAVE_GPU_HIP
#include <hipblas/hipblas.h>
#include <hipsolver/hipsolver.h>
#include <hip/hip_complex.h>
#include <hip_api_error_check.h>

extern hipblasHandle_t hip_handle;
extern hipsolverDnHandle_t hipsolverDn_handle;
#endif

#endif // ABINIT_SHARED_COMMON_SRC_17_GPU_LINALG_H

/**
 * HIP API error checking utilities.
 *
 * This header is loosely translated from its CUDA counterpart,
 * which itself is slightly adapted from cuda samples
 * https://github.com/Nvidia/cuda-samples
 *
 * Copyright (c) 2022, NVIDIA CORPORATION. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  * Neither the name of NVIDIA CORPORATION nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS ``AS IS'' AND ANY
 * EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
 * OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
#ifndef HIP_API_ERROR_CHECK_H_
#define HIP_API_ERROR_CHECK_H_

#pragma once

#include <hip/hip_runtime.h>
#include <hip/driver_types.h>
#include <hipblas/hipblas.h>
#include <hipfft/hipfft.h>
#include <hipsolver/hipsolver.h>
//#include <hipsparse.h>
#include <string>
#include <stdio.h> // for fflush
#include <cassert>

#ifndef assertm
#define assertm(exp, msg) assert(((void)msg, exp))
#endif

#ifdef __cplusplus
extern "C" {
#endif

  // this is defined in shared/common/src/16_hideleave/m_errors.F90
  void abi_cabort();

#ifdef __cplusplus
}
#endif

// when this header is used inside abinit, call abi_cabort when an error happens
// otherwise, do a regular exit
// we assume here that if config.h (from abinit is included first, then ABINIT_VERSION will be defined)
#ifdef ABINIT_VERSION
#define ABORT() abi_cabort()
#else
#define ABORT() {       \
    exit(EXIT_FAILURE); \
  }
#endif


/*
 * If symbol ALWAYS_SYNC_GPU is defined in your build system
 * then we will always call hipDeviceSynchronize, before checking
 * last error from a hip kernel call.
 *
 * It may help in case of debugging to define ALWAYS_SYNC_GPU, but for a
 * regular run, it should remains to 0.
 *
 * Note: you can also enforce syncing CPU/GPU by set environment variable HIP_LAUNCH_BLOCKING to 1
 *
 * See also: https://www.olcf.ornl.gov/wp-content/uploads/2021/06/hip_training_series_hip_debugging.pdf
 */
#ifdef ALWAYS_SYNC_GPU
#define FORCE_SYNC_GPU 1
#else
#define FORCE_SYNC_GPU 0
#endif


// HIP Runtime error messages
static const char *_hipGetErrorEnum(hipError_t error) {
  return hipGetErrorName(error);
}

// cuBLAS API errors
static const char *_hipGetErrorEnum(hipblasStatus_t error) {
  return hipblasStatusToString(error);
}

#ifdef HIP_DRIVER_API
// HIP Driver API errors
static const char *_hipGetErrorEnum(CUresult error) {
  static char unknown[] = "<unknown>";
  const char *ret = NULL;
  cuGetErrorName(error, &ret);
  return ret ? ret : unknown;
}
#endif

#ifdef HIPBLAS_API_H_
// cuBLAS API errors
static const char *_hipGetErrorEnum(cublasStatus_t error) {
  switch (error) {
    case HIPBLAS_STATUS_SUCCESS:
      return "HIPBLAS_STATUS_SUCCESS";

    case HIPBLAS_STATUS_NOT_INITIALIZED:
      return "HIPBLAS_STATUS_NOT_INITIALIZED";

    case HIPBLAS_STATUS_ALLOC_FAILED:
      return "HIPBLAS_STATUS_ALLOC_FAILED";

    case HIPBLAS_STATUS_INVALID_VALUE:
      return "HIPBLAS_STATUS_INVALID_VALUE";

    case HIPBLAS_STATUS_ARCH_MISMATCH:
      return "HIPBLAS_STATUS_ARCH_MISMATCH";

    case HIPBLAS_STATUS_MAPPING_ERROR:
      return "HIPBLAS_STATUS_MAPPING_ERROR";

    case HIPBLAS_STATUS_EXEHIPTION_FAILED:
      return "HIPBLAS_STATUS_EXEHIPTION_FAILED";

    case HIPBLAS_STATUS_INTERNAL_ERROR:
      return "HIPBLAS_STATUS_INTERNAL_ERROR";

    case HIPBLAS_STATUS_NOT_SUPPORTED:
      return "HIPBLAS_STATUS_NOT_SUPPORTED";

    case HIPBLAS_STATUS_LICENSE_ERROR:
      return "HIPBLAS_STATUS_LICENSE_ERROR";
  }

  return "<unknown>";
}
#endif

#ifdef HIPFFT_H_
// hipFFT API errors
static const char *_hipGetErrorEnum(hipfftResult error) {
  switch (error) {
    case HIPFFT_SUCCESS:
      return "HIPFFT_SUCCESS";

    case HIPFFT_INVALID_PLAN:
      return "HIPFFT_INVALID_PLAN";

    case HIPFFT_ALLOC_FAILED:
      return "HIPFFT_ALLOC_FAILED";

    case HIPFFT_INVALID_TYPE:
      return "HIPFFT_INVALID_TYPE";

    case HIPFFT_INVALID_VALUE:
      return "HIPFFT_INVALID_VALUE";

    case HIPFFT_INTERNAL_ERROR:
      return "HIPFFT_INTERNAL_ERROR";

    case HIPFFT_EXEC_FAILED:
      return "HIPFFT_EXEC_FAILED";

    case HIPFFT_SETUP_FAILED:
      return "HIPFFT_SETUP_FAILED";

    case HIPFFT_INVALID_SIZE:
      return "HIPFFT_INVALID_SIZE";

    case HIPFFT_UNALIGNED_DATA:
      return "HIPFFT_UNALIGNED_DATA";

    case HIPFFT_INCOMPLETE_PARAMETER_LIST:
      return "HIPFFT_INCOMPLETE_PARAMETER_LIST";

    case HIPFFT_INVALID_DEVICE:
      return "HIPFFT_INVALID_DEVICE";

    case HIPFFT_PARSE_ERROR:
      return "HIPFFT_PARSE_ERROR";

    case HIPFFT_NO_WORKSPACE:
      return "HIPFFT_NO_WORKSPACE";

    case HIPFFT_NOT_IMPLEMENTED:
      return "HIPFFT_NOT_IMPLEMENTED";

    case HIPFFT_NOT_SUPPORTED:
      return "HIPFFT_NOT_SUPPORTED";
  }

  return "<unknown>";
}
#endif

#ifdef _HIPSPARSE_H_
// cuSPARSE API errors
static const char *_hipGetErrorEnum(hipsparseStatus_t error) {
  switch (error) {
    case HIPSPARSE_STATUS_SUCCESS:
      return "HIPSPARSE_STATUS_SUCCESS";

    case HIPSPARSE_STATUS_NOT_INITIALIZED:
      return "HIPSPARSE_STATUS_NOT_INITIALIZED";

    case HIPSPARSE_STATUS_ALLOC_FAILED:
      return "HIPSPARSE_STATUS_ALLOC_FAILED";

    case HIPSPARSE_STATUS_INVALID_VALUE:
      return "HIPSPARSE_STATUS_INVALID_VALUE";

    case HIPSPARSE_STATUS_ARCH_MISMATCH:
      return "HIPSPARSE_STATUS_ARCH_MISMATCH";

    case HIPSPARSE_STATUS_MAPPING_ERROR:
      return "HIPSPARSE_STATUS_MAPPING_ERROR";

    case HIPSPARSE_STATUS_EXECUTION_FAILED:
      return "HIPSPARSE_STATUS_EXECUTION_FAILED";

    case HIPSPARSE_STATUS_INTERNAL_ERROR:
      return "HIPSPARSE_STATUS_INTERNAL_ERROR";

    case HIPSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
      return "HIPSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";

    case HIPSPARSE_STATUS_ZERO_PIVOT:
      return "HIPSPARSE_STATUS_ZERO_PIVOT";

    case HIPSPARSE_STATUS_NOT_SUPPORTED:
      return "HIPSPARSE_STATUS_NOT_SUPPORTED";

    case HIPSPARSE_STATUS_INSUFFICIENT_RESOURCES:
      return "HIPSPARSE_STATUS_INSUFFICIENT_RESOURCES";
  }

  return "<unknown>";
}
#endif

#ifdef HIPSOLVER_H
// cuSOLVER API errors
static const char *_hipGetErrorEnum(hipsolverStatus_t error) {
  switch (error) {
    case HIPSOLVER_STATUS_SUCCESS:
      return "HIPSOLVER_STATUS_SUCCESS";
    case HIPSOLVER_STATUS_NOT_INITIALIZED:
      return "HIPSOLVER_STATUS_NOT_INITIALIZED";
    case HIPSOLVER_STATUS_ALLOC_FAILED:
      return "HIPSOLVER_STATUS_ALLOC_FAILED";
    case HIPSOLVER_STATUS_INVALID_VALUE:
      return "HIPSOLVER_STATUS_INVALID_VALUE";
    case HIPSOLVER_STATUS_ARCH_MISMATCH:
      return "HIPSOLVER_STATUS_ARCH_MISMATCH";
    case HIPSOLVER_STATUS_MAPPING_ERROR:
      return "HIPSOLVER_STATUS_MAPPING_ERROR";
    case HIPSOLVER_STATUS_EXECUTION_FAILED:
      return "HIPSOLVER_STATUS_EXECUTION_FAILED";
    case HIPSOLVER_STATUS_INTERNAL_ERROR:
      return "HIPSOLVER_STATUS_INTERNAL_ERROR";
    case HIPSOLVER_STATUS_HANDLE_IS_NULLPTR:
      return "HIPSOLVER_STATUS_HANDLE_IS_NULLPTR";
    case HIPSOLVER_STATUS_INVALID_ENUM:
      return "HIPSOLVER_STATUS_INVALID_ENUM";
    case HIPSOLVER_STATUS_UNKNOWN:
      return "HIPSOLVER_STATUS_UNKNOWN";
    case HIPSOLVER_STATUS_NOT_SUPPORTED:
      return "HIPSOLVER_STATUS_NOT_SUPPORTED ";
  }

  return "<unknown>";
}
#endif


#ifdef HIPRAND_H_
// hipRAND API errors
static const char *_hipGetErrorEnum(hiprandStatus_t error) {
  switch (error) {
    case HIPRAND_STATUS_SUCCESS:
      return "HIPRAND_STATUS_SUCCESS";

    case HIPRAND_STATUS_VERSION_MISMATCH:
      return "HIPRAND_STATUS_VERSION_MISMATCH";

    case HIPRAND_STATUS_NOT_INITIALIZED:
      return "HIPRAND_STATUS_NOT_INITIALIZED";

    case HIPRAND_STATUS_ALLOCATION_FAILED:
      return "HIPRAND_STATUS_ALLOCATION_FAILED";

    case HIPRAND_STATUS_TYPE_ERROR:
      return "HIPRAND_STATUS_TYPE_ERROR";

    case HIPRAND_STATUS_OUT_OF_RANGE:
      return "HIPRAND_STATUS_OUT_OF_RANGE";

    case HIPRAND_STATUS_LENGTH_NOT_MULTIPLE:
      return "HIPRAND_STATUS_LENGTH_NOT_MULTIPLE";

    case HIPRAND_STATUS_DOUBLE_PRECISION_REQUIRED:
      return "HIPRAND_STATUS_DOUBLE_PRECISION_REQUIRED";

    case HIPRAND_STATUS_LAUNCH_FAILURE:
      return "HIPRAND_STATUS_LAUNCH_FAILURE";

    case HIPRAND_STATUS_PREEXISTING_FAILURE:
      return "HIPRAND_STATUS_PREEXISTING_FAILURE";

    case HIPRAND_STATUS_INITIALIZATION_FAILED:
      return "HIPRAND_STATUS_INITIALIZATION_FAILED";

    case HIPRAND_STATUS_ARCH_MISMATCH:
      return "HIPRAND_STATUS_ARCH_MISMATCH";

    case HIPRAND_STATUS_NOT_IMPLEMENTED:
      return "HIPRAND_STATUS_NOT_IMPLEMENTED";

    case HIPRAND_STATUS_INTERNAL_ERROR:
      return "HIPRAND_STATUS_INTERNAL_ERROR";
  }

  return "<unknown>";
}
#endif

template <typename T>
void check_hip_error(T result, char const *const func, const char *const file,
           int const line) {
  if (result) {
    fprintf(stderr, "HIP error at :\n\t%s:%d \n\treturn code = %d (%s)\n\t\"%s\" \n\n",
        file, line,
        static_cast<unsigned int>(result), _hipGetErrorEnum(result), func);
    fflush(stderr);
    ABORT();
  }
}

// This will output the proper HIP error strings in the event
// that a HIP host call returns an error
#define CHECK_HIP_ERROR(value) check_hip_error((value), #value, __FILE__, __LINE__)
#define HIP_API_CHECK(value) CHECK_HIP_ERROR(value)


// This will output the proper error string when calling hipGetLastError
#define getLastHipError(msg) __getLastHipError(msg, __FILE__, __LINE__)
#define GET_LAST_HIP_ERROR(msg) getLastHipError(msg)

inline void __getLastHipError(const char *errorMessage, const char *file,
                               const int line) {
  hipError_t err = hipGetLastError();

  if (hipSuccess != err) {
    fprintf(stderr,
            "%s(%i) : getLastHipError() HIP error :"
            " %s : (%d) %s.\n",
            file, line, errorMessage, static_cast<int>(err),
            hipGetErrorString(err));

    // Make sure we call HIP Device Reset before exiting
    err = hipDeviceReset();

    ABORT();
  }
}

// This will only print the proper error string when calling hipGetLastError
// but not exit program incase error detected.
#define printLastHipError(msg) __printLastHipError(msg, __FILE__, __LINE__)
#define PRINT_LAST_HIP_ERROR(msg) printLastHipError(msg)

inline void __printLastHipError(const char *errorMessage, const char *file,
                                 const int line) {
  hipError_t err = hipGetLastError();

  if (hipSuccess != err) {
    fprintf(stderr,
            "%s(%i) : getLastHipError() HIP error :"
            " %s : (%d) %s.\n",
            file, line, errorMessage, static_cast<int>(err),
            hipGetErrorString(err));
  }
}


/**
 * enum used below; can be used as the second argument of macro
 * HIP_KERNEL_CHECK
 */
enum device_sync_t {
  DEVICE_NO_SYNC = 0,
  DEVICE_SYNC = 1
};


/**
 * a simple macro helper:
 * GET_KERNEL_CHECK_MACRO always picks the 3rd arg
 *
 * see https://stackoverflow.com/questions/11761703/overloading-macro-on-number-of-arguments
 */
#define GET_KERNEL_CHECK_MACRO(_1,_2,NAME,...) NAME

/**
 * another simple macro helper :
 * - if HIP_KERNEL_CHECK is called with only 1 argument, then HIP_KERNEL_CHECK1 is chosen
 * - if HIP_KERNEL_CHECK is called with      2 arguments, then HIP_KERNEL_CHECK2 is chosen
 *
 *
 * this is the macro we want to call
 */
#define HIP_KERNEL_CHECK(...) GET_KERNEL_CHECK_MACRO(__VA_ARGS__, HIP_KERNEL_CHECK2, HIP_KERNEL_CHECK1)(__VA_ARGS__)

/**
 * Preprocessor macro helping to retrieve the exact code
 * location where the error was emitted.
 *
 * Default behavior, don't synchronize device
 */
#define HIP_KERNEL_CHECK1(msg) hip_kernel_check((msg), __FILE__, __LINE__, DEVICE_NO_SYNC)

/**
 * Same as above, but let the user chose if we want to synchronize device.
 */
#define HIP_KERNEL_CHECK2(msg,sync) hip_kernel_check((msg), __FILE__, __LINE__, sync)

/**
 * Check last HIP kernel call status.
 * If it was not successfull then print error message.
 *
 * \param[in] errstr error message to print
 * \param[in] file source filename where error occured
 * \param[in] line line number where error occured
 * \param[in] sync integer, 0 means no device synchronization
 */
static void hip_kernel_check(const char* errstr,
                              const char* file,
                              const int   line,
                              const int   sync)
{

  auto status = hipGetLastError();

  if (sync or FORCE_SYNC_GPU) {
    //fprintf(stderr, "syncing device\n");
    hipDeviceSynchronize();
  }

  if (status != hipSuccess) {
    fprintf(stderr,
            "%s(%i) : getLastHipError() HIP error :"
            " %s : (%d) %s.\n",
            file, line, errstr, static_cast<int>(status),
            hipGetErrorString(status));

    //hipDeviceReset();
    //exit(EXIT_FAILURE);
  }

} // hip_kernel_check

inline const char* _ConvertSMVer2ArchName(int major, int minor) {
  // Defines for GPU Architecture types (using the SM version to determine
  // the GPU Arch name)
  typedef struct {
    int SM;  // 0xMm (hexidecimal notation), M = SM Major version,
    // and m = SM minor version
    const char* name;
  } sSMtoArchName;

  sSMtoArchName nGpuArchNameSM[] = {
      {0x30, "Kepler"},
      {0x32, "Kepler"},
      {0x35, "Kepler"},
      {0x37, "Kepler"},
      {0x50, "Maxwell"},
      {0x52, "Maxwell"},
      {0x53, "Maxwell"},
      {0x60, "Pascal"},
      {0x61, "Pascal"},
      {0x62, "Pascal"},
      {0x70, "Volta"},
      {0x72, "Xavier"},
      {0x75, "Turing"},
      {0x80, "Ampere"},
      {0x86, "Ampere"},
      {0x89, "AdaLovelace"},
      {0x90, "Hopper"},
      {-1, "Graphics Device"}};

  int index = 0;

  while (nGpuArchNameSM[index].SM != -1) {
    if (nGpuArchNameSM[index].SM == ((major << 4) + minor)) {
      return nGpuArchNameSM[index].name;
    }

    index++;
  }

  // If we don't find the values, we default use the previous one
  // to run properly
  printf(
      "MapSMtoArchName for SM %d.%d is undefined."
      "  Default to use %s\n",
      major, minor, nGpuArchNameSM[index - 1].name);
  return nGpuArchNameSM[index - 1].name;
}


#ifdef __HIP_RUNTIME_H__

// General check for HIP GPU SM Capabilities
inline bool checkHipCapabilities(int major_version, int minor_version) {
  int dev;
  int major = 0, minor = 0;

  CHECK_HIP_ERROR(hipGetDevice(&dev));
  CHECK_HIP_ERROR(hipDeviceGetAttribute(&major, hipDevAttrComputeCapabilityMajor, dev));
  CHECK_HIP_ERROR(hipDeviceGetAttribute(&minor, hipDevAttrComputeCapabilityMinor, dev));

  if ((major > major_version) ||
      (major == major_version &&
       minor >= minor_version)) {
    printf("  Device %d: <%16s >, Compute SM %d.%d detected\n", dev,
           _ConvertSMVer2ArchName(major, minor), major, minor);
    return true;
  } else {
    printf(
        "  No GPU device was found that can support "
        "HIP compute capability %d.%d.\n",
        major_version, minor_version);
    return false;
  }
}
#endif /* __HIP_RUNTIME_H__ */


#endif /* HIP_API_ERROR_CHECK_H_ */

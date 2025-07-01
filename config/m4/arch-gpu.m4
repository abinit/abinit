# -*- Autoconf -*-
#
# Copyright (C) 2009-2025 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for GPU
#


# _ABI_GPU_CHECK_CUDA()
# ------------------------
#
# Check whether the Cuda library is working.
#
AC_DEFUN([_ABI_GPU_CHECK_CUDA],[
  # Init
  abi_gpu_cuda_serial="no"
  abi_gpu_cuda_old="unknown"

  # Display variables
  AC_MSG_NOTICE([Cuda incs: ${abi_gpu_cuda_incs}])
  AC_MSG_NOTICE([Cuda libs: ${abi_gpu_cuda_libs}])

  # Prepare environment
  ABI_ENV_BACKUP
  CPPFLAGS="${CPPFLAGS} ${abi_gpu_cuda_incs}"
  LDFLAGS="${CC_LDFLAGS} ${CC_LDFLAGS_GPU}"
  abi_saved_LIBS="${LIBS}"
  LIBS="${abi_gpu_cuda_libs} ${LIBS}"
  AC_LANG_PUSH([C])

  # Check usability of headers
  AC_CHECK_HEADERS([cuda_runtime_api.h cufft.h cublas.h cusolver_common.h])

  # Look for libraries and routines
  AC_MSG_CHECKING([whether Cuda programs can be compiled])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#if defined HAVE_CUDA_RUNTIME_API_H
#include "cuda_runtime_api.h"
#endif
    ]],
    [[
      cudaError_t err;
      int *count;
      err = cudaGetDeviceCount(count);
    ]])], [abi_gpu_cuda_serial="yes"], [])
  AC_MSG_RESULT([${abi_gpu_cuda_serial}])

  # Do we have an old version of Cuda?
  AC_MSG_CHECKING([whether we have Cuda < 4])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#if defined HAVE_CUDA_RUNTIME_API_H
#include "cuda_runtime_api.h"
#endif
    ]],
    [[
      cudaDeviceReset();
    ]])], [abi_gpu_cuda_old="no"], [abi_gpu_cuda_old="yes"])
  AC_MSG_RESULT([${abi_gpu_cuda_old}])

  # Check ISO C Binding (Fortran)
  if test "${fc_has_iso_c_binding}" != "yes"; then
    AC_MSG_WARN([your Fortran compiler does not provide any ISO C binding module])
  fi

  # call this macro here to make sure variable `abi_gpu_cuda_version_10`
  #
  # check if we are using CUDA runtime version at least 10
  # version 10 of CUDA is required for NVTX (header only)
  #
  AC_MSG_CHECKING([whether we have Cuda >= 10])
  AC_LANG_PUSH(C)
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [[
      #include <cuda_runtime_api.h>
      ]],
      [[
#if CUDART_VERSION < 10000
#error
#endif
      ]])], [abi_gpu_cuda_version_10="yes"], [abi_gpu_cuda_version_10="no"])
  AC_MSG_RESULT([${abi_gpu_cuda_version_10}])

  # call this macro here to make sure variable `abi_gpu_cuda_version_129`
  #
  # check if we are using CUDA runtime version at least 12.9
  # From version 12.9 of CUDA, NVTX library is named differently
  #
  AC_MSG_CHECKING([whether we have Cuda >= 12.9])
  AC_LANG_PUSH(C)
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM(
      [[
      #include <cuda_runtime_api.h>
      ]],
      [[
#if CUDART_VERSION < 12090
#error
#endif
      ]])], [abi_gpu_cuda_version_129="yes"], [abi_gpu_cuda_version_129="no"])
  AC_MSG_RESULT([${abi_gpu_cuda_version_129}])

  if test "${abi_gpu_markers_enable}" = "yes"; then
    if test "${abi_gpu_cuda_version_10}" = "yes"; then
      AC_MSG_CHECKING([for NVTX v3 GPU markers])
      nvtx_libname="nvToolsExt"
      if test "${abi_gpu_cuda_version_129}" = "yes"; then
        nvtx_libname="nvtx3interop"
      fi
      if test -e "${abi_gpu_cuda_libdir}/lib${nvtx_libname}.${abi_so_ext}"; then
        # always add link flags to nvtx if available
        if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
          abi_gpu_cuda_libs="-l${nvtx_libname} ${abi_gpu_cuda_libs}"
        else
          abi_gpu_cuda_libs="${abi_gpu_cuda_libs} -l${nvtx_libname}"
        fi
        abi_gpu_nvtx_v3="yes"
        abi_result="${abi_result} nvtx_v3"
      else
        AC_MSG_ERROR([Cuda NVTX: ${abi_gpu_cuda_libdir}/lib${nvtx_libname}.${abi_so_ext} not found])
      fi
      AC_MSG_RESULT([${abi_result}])
    else
      AC_MSG_ERROR([Cuda NVTX was requested but is not available for CUDA < v10])
    fi
  fi

  # Restore build environment
  AC_LANG_POP([C])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE
]) # _ABI_GPU_CHECK_CUDA



                    ########################################



# _ABI_GPU_INIT_CUDA()
# --------------------
#
# Looks for an implementation of Cuda, using the provided prefix.
#
AC_DEFUN([_ABI_GPU_INIT_CUDA],[
  # Init
  abi_gpu_cuda_has_cc="no"
  abi_gpu_cuda_has_common="no"
  abi_gpu_cuda_has_fft="no"
  abi_gpu_cuda_has_cusolver="no"
  abi_gpu_cuda_has_incs="no"
  abi_gpu_cuda_has_libs="no"
  abi_gpu_cuda_has_linalg="no"
  abi_gpu_cuda_has_runtime="no"
  abi_gpu_cuda_libdir=""
  abi_gpu_cuda_incs="${GPU_CPPFLAGS}"
  abi_gpu_cuda_libs="${GPU_LIBS}"
  abi_gpu_cuda_root="${abi_gpu_prefix}"
  abi_gpu_cuda_math_root=""
  abi_gpu_cuda_version_10="unknown"
  abi_gpu_nvtx_v3="unknown"

  test -z "${abi_gpu_cuda_incs}" && abi_gpu_cuda_incs="${CUDA_CPPFLAGS}"
  test -z "${abi_gpu_cuda_libs}" && abi_gpu_cuda_libs="${CUDA_LIBS}"

  # Make use of the CUDA_ROOT, CUDA_HOME, CUDA_PREFIX, CUDA_PATH and NVHPC_CUDA_HOME
  # environment variables to find CUDA from user environment
  if test "${abi_gpu_cuda_root}" = ""; then
    abi_gpu_cuda_root="${CUDA_ROOT}"
  fi
  if test "${abi_gpu_cuda_root}" = ""; then
    abi_gpu_cuda_root="${CUDA_HOME}"
  fi
  if test "${abi_gpu_cuda_root}" = ""; then
    abi_gpu_cuda_root="${CUDA_PREFIX}"
  fi
  if test "${abi_gpu_cuda_root}" = ""; then
    abi_gpu_cuda_root="${NVHPC_CUDA_HOME}"
  fi

  # If CUDA is provided by NVHPC SDK, math libraries such as cuBLAS, cuFFT or cuSOLVER
  # are located in a separate directory
  tmp_math_root=`echo ${abi_gpu_cuda_root} | sed 's/cuda/math_libs/'`
  if test -s "${tmp_math_root}"; then
    abi_gpu_cuda_math_root="${tmp_math_root}"
  else
    abi_gpu_cuda_math_root="${abi_gpu_cuda_root}"
  fi

  # Check whether to look for generic files
  if test "${abi_gpu_cuda_root}" = ""; then

    # nVidia C compiler
    if test "${NVCC}" = ""; then
      AC_CHECK_PROGS(NVCC,[nvcc])
    fi
    AC_MSG_CHECKING([for the nVidia C compiler])
    if test "${NVCC}" = ""; then
      AC_MSG_RESULT([none found])
    else
      abi_gpu_cuda_has_cc="yes"
      AC_MSG_RESULT([${NVCC}])
    fi

  else

    # nVidia C compiler
    AC_MSG_CHECKING([for the nVidia C compiler])
    if test "${NVCC}" = ""; then
      if test -x "${abi_gpu_cuda_root}/bin/nvcc"; then
        abi_gpu_cuda_has_cc="yes"
        NVCC="${abi_gpu_cuda_root}/bin/nvcc"
      fi
    fi
    if test "${NVCC}" = ""; then
      AC_MSG_RESULT([none found])
    else
      AC_MSG_RESULT([${NVCC}])
    fi

  fi

  # Populate NVCC C flags with interesting values
  abi_gpu_nvcc_flags_base="-O3 -Xptxas=-v --use_fast_math --compiler-options -O3,-fPIC -DHAVE_CONFIG_H --forward-unknown-opts"

  # Check whether to look for generic files
  if test "${abi_gpu_cuda_root}" != ""; then

    # Headers
    AC_MSG_CHECKING([for Cuda headers])
    abi_result=""
    if test -s "${abi_gpu_cuda_root}/include/cuda_runtime_api.h"; then
      if test "${GPU_CPPFLAGS}" = "" -a "${CUDA_CPPFLAGS}" = ""; then
        if test "${abi_gpu_cuda_root}" == "${abi_gpu_cuda_math_root}"; then
          abi_gpu_cuda_incs="-I${abi_gpu_cuda_root}/include"
        else
          abi_gpu_cuda_incs="-I${abi_gpu_cuda_root}/include -I${abi_gpu_cuda_math_root}/include"
        fi
      fi
      abi_gpu_cuda_has_incs="yes"
      abi_result="${abi_result} run-time"
    fi
    if test -s "${abi_gpu_cuda_math_root}/include/cufft.h"; then
      abi_result="${abi_result} fft"
    fi
    if test -s "${abi_gpu_cuda_math_root}/include/cublas.h"; then
      abi_result="${abi_result} blas"
    fi
    if test -s "${abi_gpu_cuda_math_root}/include/cusolver_common.h"; then
      abi_result="${abi_result} cusolver"
    fi
    if test -s "${abi_gpu_cuda_root}/SDK/C/common/inc/cutil.h"; then
      if test "${GPU_CPPFLAGS}" = "" -a "${CUDA_CPPFLAGS}" = ""; then
        abi_gpu_cuda_incs="-I${abi_gpu_cuda_root}/SDK/C/common/inc ${abi_gpu_cuda_incs}"
      fi
      abi_result="${abi_result} sdk"
    fi
    if test "${abi_result}" = ""; then
      abi_result="none"
    fi
    AC_MSG_RESULT([${abi_result}])

    # Libraries
    AC_MSG_CHECKING([for Cuda libraries])
    abi_result=""
    if test "${abi_cpu_64bits}" = "yes"; then
      abi_gpu_cuda_libdir="${abi_gpu_cuda_root}/lib64"
      abi_gpu_cuda_math_libdir="${abi_gpu_cuda_math_root}/lib64"
    else
      abi_gpu_cuda_libdir="${abi_gpu_cuda_root}/lib"
      abi_gpu_cuda_math_libdir="${abi_gpu_cuda_math_root}/lib"
    fi
    for tmp_cuda_dir in \
        /usr/lib \
        /usr/local/lib \
        /usr/local/lib${abi_cpu_bits} \
        /usr/lib/nvidia-current \
        /usr/local/lib/nvidia-current \
        /usr/local/lib${abi_cpu_bits}/nvidia-current; do
      if test "${abi_gpu_cuda_has_common}" = "no"; then
        if test -e "${tmp_cuda_dir}/libcuda.${abi_so_ext}"; then
          abi_gpu_cuda_has_libs="yes"
          abi_gpu_cuda_has_common="yes"
          abi_result="${abi_result} common"
        fi
      fi
    done
    if test -e "${abi_gpu_cuda_libdir}/libcudart.${abi_so_ext}"; then
      if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
        abi_gpu_cuda_libs="-lcudart"
      fi
      abi_gpu_cuda_has_libs="yes"
      abi_gpu_cuda_has_runtime="yes"
      abi_result="${abi_result} run-time"
    fi
    if test "${abi_gpu_cuda_has_libs}" = "yes"; then
      if test -e "${abi_gpu_cuda_math_libdir}/libcufft.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcufft ${abi_gpu_cuda_libs}"
        fi
        abi_gpu_cuda_has_fft="yes"
        abi_result="${abi_result} fft"
      fi
      if test -e "${abi_gpu_cuda_math_libdir}/libcublas.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcublas ${abi_gpu_cuda_libs}"
        fi
        abi_gpu_cuda_has_linalg="yes"
        abi_result="${abi_result} blas"
      fi
      if test -e "${abi_gpu_cuda_math_libdir}/libcublasLt.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcublasLt ${abi_gpu_cuda_libs}"
        fi
      fi
      if test -e "${abi_gpu_cuda_math_libdir}/libcusolver.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcusolver ${abi_gpu_cuda_libs}"
        fi
        abi_gpu_cuda_has_cusolver="yes"
        abi_result="${abi_result} cusolver"
      fi
      if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
        if test "${abi_gpu_cuda_root}" == "${abi_gpu_cuda_math_root}"; then
          abi_gpu_cuda_libs="-L${abi_gpu_cuda_libdir} ${abi_gpu_cuda_libs}"
        else
          abi_gpu_cuda_libs="-L${abi_gpu_cuda_libdir} -L${abi_gpu_cuda_math_libdir} ${abi_gpu_cuda_libs}"
        fi
      fi
    fi

    if test -s "${abi_gpu_cuda_root}/SDK/C/lib/libcutil.a"; then
      if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
        abi_gpu_cuda_libs="-L${abi_gpu_cuda_root}/SDK/C/lib -lcutil ${abi_gpu_cuda_libs}"
      fi
      abi_result="${abi_result} sdk"
    fi
    if test "${abi_result}" = ""; then
      abi_result="none"
    fi
    AC_MSG_RESULT([${abi_result}])
    if test "${GPU_LIBS}" = "" -a "${CUDA_LIBS}" = ""; then
      abi_gpu_cuda_libs="${abi_gpu_cuda_libs} -lcuda"
    fi
    if test "${abi_gpu_cuda_has_common}" = "no"; then
      AC_MSG_WARN([could not find libcuda.${abi_so_ext}])
    fi

    # add standart libc++ link flags
    abi_gpu_cuda_libs="${abi_gpu_cuda_libs} -lstdc++"

    # C and C++ link flags
    AC_MSG_CHECKING([for Cuda link flags])
    if test "${CC_LDFLAGS_GPU}" = ""; then
      if test "${abi_cpu_64bits}" = "yes"; then
        CC_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_cuda_root}/lib64"
      else
        CC_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_cuda_root}/lib"
      fi
    fi
    if test "${CXX_LDFLAGS_GPU}" = ""; then
      if test "${abi_cpu_64bits}" = "yes"; then
        CXX_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_cuda_root}/lib64"
      else
        CXX_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_cuda_root}/lib"
      fi
    fi
    AC_MSG_RESULT([${CC_LDFLAGS_GPU}])

  fi # abi_gpu_cuda_root

  AC_MSG_NOTICE([Cuda incs: ${abi_gpu_cuda_incs}])
  AC_MSG_NOTICE([Cuda libs: ${abi_gpu_cuda_libs}])

  AC_SUBST(abi_gpu_nvtx_v3)
  AC_SUBST(NVCC_CFLAGS_ARCH)
  AC_SUBST(abi_gpu_nvcc_flags_base)

]) # _ABI_GPU_INIT_CUDA



                    ########################################


# _ABI_GPU_CHECK_HIP()
# ------------------------
#
# Check whether AMD ROCm/HIP libraries are working.
#
AC_DEFUN([_ABI_GPU_CHECK_HIP],[
  # Init
  abi_gpu_hip_serial="no"
  abi_gpu_hip_old="unknown"

  # Display variables
  AC_MSG_NOTICE([ROCm/HIP incs: ${abi_gpu_hip_incs}])
  AC_MSG_NOTICE([ROCm/HIP libs: ${abi_gpu_hip_libs}])

  # Prepare environment
  ABI_ENV_BACKUP
  CPPFLAGS="${CPPFLAGS} ${abi_gpu_hip_incs}"
  LDFLAGS="${CC_LDFLAGS} ${CC_LDFLAGS_GPU}"
  abi_saved_LIBS="${LIBS}"
  LIBS="${abi_gpu_hip_libs} ${LIBS}"
  AC_LANG_PUSH([C])

  # Check usability of headers
  AC_CHECK_HEADERS([hip/hip_runtime_api.h hipfft/hipfft.h hipblas/hipblas.h hipsolver/hipsolver.h])

  # Look for libraries and routines
  AC_MSG_CHECKING([whether HIP programs can be compiled])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <hip/hip_runtime_api.h>
    ]],
    [[
      hipError_t err;
      int *count;
      err = hipGetDeviceCount(count);
    ]])], [abi_gpu_hip_serial="yes"], [])
  AC_MSG_RESULT([${abi_gpu_hip_serial}])

  # Check ISO C Binding (Fortran)
  if test "${fc_has_iso_c_binding}" != "yes"; then
    AC_MSG_WARN([your Fortran compiler does not provide any ISO C binding module])
  fi

  if test "${abi_gpu_markers_enable}" = "yes"; then
    if test -e "${abi_gpu_hip_libdir}/libroctx64.${abi_so_ext}"; then
      # always add link flags to roctx if available
      if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
        abi_gpu_hip_libs="-lroctx64 ${abi_gpu_hip_libs}"
      else
        abi_gpu_hip_libs="${abi_gpu_hip_libs} -lroctx64"
      fi
      abi_gpu_roctx="yes"
      abi_result="${abi_result} roctx"
    else
      AC_MSG_ERROR([AMD ROCtx: ${abi_gpu_hip_libdir}/libroctx64.${abi_so_ext} not found])
    fi
  fi

  # Restore build environment
  AC_LANG_POP([C])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE
]) # _ABI_GPU_CHECK_HIP



                    ########################################



# _ABI_GPU_INIT_HIP()
# --------------------
#
# Looks for an installation of ROCm/HIP, using the provided prefix.
#
AC_DEFUN([_ABI_GPU_INIT_HIP],[
  # Init
  abi_gpu_hip_has_cc="no"
  abi_gpu_hip_has_common="no"
  abi_gpu_hip_has_fft="no"
  abi_gpu_hip_has_cusolver="no"
  abi_gpu_hip_has_incs="no"
  abi_gpu_hip_has_libs="no"
  abi_gpu_hip_has_linalg="no"
  abi_gpu_hip_has_runtime="no"
  abi_gpu_hip_libdir=""
  abi_gpu_hip_incs="${GPU_CPPFLAGS}"
  abi_gpu_hip_libs="${GPU_LIBS}"
  abi_gpu_hip_root="${abi_gpu_prefix}"
  abi_gpu_hip_version_10="unknown"
  abi_gpu_nvtx_v3="unknown"

  test -z "${abi_gpu_hip_incs}" && abi_gpu_hip_incs="${ROCM_CPPFLAGS}"
  test -z "${abi_gpu_hip_libs}" && abi_gpu_hip_libs="${ROCM_LIBS}"

  # Make use of the ROCM_ROOT, ROCM_PATH and ROCM_PREFIX environment variables
  if test "${abi_gpu_hip_root}" = ""; then
    abi_gpu_hip_root="${ROCM_ROOT}"
  fi
  if test "${abi_gpu_hip_root}" = ""; then
    abi_gpu_hip_root="${ROCM_PATH}"
  fi
  if test "${abi_gpu_hip_root}" = ""; then
    abi_gpu_hip_root="${ROCM_PREFIX}"
  fi

  # Check whether to look for generic files
  if test "${abi_gpu_hip_root}" = ""; then

   # AMD HIP hipcc is not currently used in the code
   AC_MSG_NOTICE("Skipping AMD HIPCC detection (not used)")
   # # AMD HIP C compiler
   # if test "${NVCC}" = ""; then
   #   AC_CHECK_PROGS(NVCC,[hipcc])
   # fi
   # AC_MSG_CHECKING([for the AMD C compiler])
   # if test "${NVCC}" = ""; then
   #   AC_MSG_RESULT([none found])
   # else
   #   abi_gpu_hip_has_cc="yes"
   #   AC_MSG_RESULT([${NVCC}])
   # fi

  else

   # AMD HIP hipcc is not currently used in the code
   AC_MSG_NOTICE("Skipping AMD HIPCC detection (not used)")
   # # AMD HIP C compiler
   # AC_MSG_CHECKING([for the AMD HIP C compiler])
   # if test "${NVCC}" = ""; then
   #   if test -x "${abi_gpu_hip_root}/bin/hipcc"; then
   #     abi_gpu_hip_has_cc="yes"
   #     NVCC="${abi_gpu_hip_root}/bin/hipcc"
   #   fi
   # fi
   # if test "${NVCC}" = ""; then
   #   AC_MSG_RESULT([none found])
   # else
   #   AC_MSG_RESULT([${NVCC}])
   # fi

    # Headers
    AC_MSG_CHECKING([for HIP headers])
    abi_result=""
    if test -s "${abi_gpu_hip_root}/include/hip/hip_runtime_api.h"; then
      if test "${GPU_CPPFLAGS}" = "" -a "${ROCM_CPPFLAGS}" = ""; then
        # Use HIP with binding to AMD ROCM by default.
        # Binding to CUDA, with __HIP_PLATFORM_NVIDIA__ isn't tested
        abi_gpu_hip_incs="-D__HIP_PLATFORM_AMD__ -I${abi_gpu_hip_root}/include"
      fi
      abi_gpu_hip_has_incs="yes"
      abi_result="${abi_result} run-time"
    fi
    if test -s "${abi_gpu_hip_root}/include/hipfft/hipfft.h"; then
      abi_result="${abi_result} fft"
    fi
    if test -s "${abi_gpu_hip_root}/include/hipblas/hipblas.h"; then
      abi_result="${abi_result} blas"
    fi
    if test -s "${abi_gpu_hip_root}/include/hipsolver/hipsolver.h"; then
      abi_result="${abi_result} hipsolver"
    fi
    if test "${abi_result}" = ""; then
      abi_result="none"
    fi
    AC_MSG_RESULT([${abi_result}])

    # Libraries
    AC_MSG_CHECKING([for ROCm/HIP libraries])
    abi_result=""
    abi_gpu_hip_libdir="${abi_gpu_hip_root}/lib"
    if test -e "${abi_gpu_hip_libdir}/libamdhip64.${abi_so_ext}"; then
      if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
        abi_gpu_hip_libs="-lamdhip64"
      fi
      abi_gpu_hip_has_libs="yes"
      abi_gpu_hip_has_runtime="yes"
      abi_result="${abi_result} run-time"
    fi
    if test "${abi_gpu_hip_has_libs}" = "yes"; then
      if test -e "${abi_gpu_hip_libdir}/libhipfft.${abi_so_ext}" -a "${abi_gpu_hip_libdir}/librocfft.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
          abi_gpu_hip_libs="-lrocfft -lhipfft ${abi_gpu_hip_libs}"
        fi
        abi_gpu_hip_has_fft="yes"
        abi_result="${abi_result} fft"
      fi
      if test -e "${abi_gpu_hip_libdir}/libhipblas.${abi_so_ext}" -a "${abi_gpu_hip_libdir}/librocblas.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
          abi_gpu_hip_libs="-lrocblas -lhipblas ${abi_gpu_hip_libs}"
        fi
        abi_gpu_hip_has_linalg="yes"
        abi_result="${abi_result} blas"
      fi
      if test -e "${abi_gpu_hip_libdir}/libhipsolver.${abi_so_ext}" -a "${abi_gpu_hip_libdir}/librocsolver.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
          abi_gpu_hip_libs="-lrocsolver -lhipsolver ${abi_gpu_hip_libs}"
        fi
        abi_gpu_hip_has_hipsolver="yes"
        abi_result="${abi_result} hipsolver"
      fi
      if test "${GPU_LIBS}" = "" -a "${ROCM_LIBS}" = ""; then
        abi_gpu_hip_libs="-L${abi_gpu_hip_libdir} ${abi_gpu_hip_libs}"
      fi
    fi

    if test "${abi_result}" = ""; then
      abi_result="none"
    fi
    AC_MSG_RESULT([${abi_result}])

    # add standart libc++ link flags
    abi_gpu_hip_libs="${abi_gpu_hip_libs} -lstdc++"

    # C and C++ link flags
    AC_MSG_CHECKING([for Cuda link flags])
    if test "${CC_LDFLAGS_GPU}" = ""; then
      CC_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_hip_root}/lib"
    fi
    if test "${CXX_LDFLAGS_GPU}" = ""; then
      CXX_LDFLAGS_GPU="-Wl,-rpath=${abi_gpu_hip_root}/lib"
    fi
    AC_MSG_RESULT([${CC_LDFLAGS_GPU}])

  fi # abi_gpu_hip_root

  AC_MSG_NOTICE([ROCm/HIP incs: ${abi_gpu_hip_incs}])
  AC_MSG_NOTICE([ROCm/HIP libs: ${abi_gpu_hip_libs}])

  AC_SUBST(abi_gpu_roctx)
  AC_SUBST(abi_gpu_hip_libdir)

]) # _ABI_GPU_INIT_HIP



                    ########################################



# ABI_GPU_INIT()
# --------------
#config/m4/arch-gpu.m4
# Looks for an implementation of GPU, using the provided prefix.
# Note 1: this is a convenience feature, purely for comfort.
# Note 2: it should be run as early as possible.
#
AC_DEFUN([ABI_GPU_INIT],[
  # Delegate most of the initialization to Steredeg
  SD_GPU_INIT([optional warn])

  # Init
  abi_gpu_complete="unknown"
  abi_gpu_enable="${sd_gpu_enable}"
  abi_gpu_markers_enable="${sd_gpu_markers_enable}"
  abi_gpu_has_cc="no"
  abi_gpu_has_fft="no"
  abi_gpu_has_incs="no"
  abi_gpu_has_libs="no"
  abi_gpu_has_linalg="no"
  abi_gpu_usable="no"
  abi_gpu_fcflags=""
  abi_gpu_ldflags=""
  abi_gpu_prefix="${sd_gpu_prefix}"
  abi_gpu_flavor="${sd_gpu_flavor}"
  abi_gpu_incs="${GPU_CPPFLAGS}"
  abi_gpu_libs="${GPU_LIBS}"
  abi_gpu_arch="${GPU_ARCH}"

  if test "${abi_gpu_enable}" = "yes"; then
    if test "${abi_gpu_arch}" = "" -a "${NVCC_CFLAGS_ARCH}" = ""; then
      AC_MSG_ERROR([GPU support is enabled but no GPU target architecture has been provided.
                    Please set GPU_ARCH variable accordingly to your target GPU with either
                    NVIDIA compute capability or AMDGPU target, for example:
                    GPU_ARCH=80       , for NVIDIA A100
                    GPU_ARCH=gfx90a   , for AMD Instinct MI250])
    fi
  fi

  if test "${abi_gpu_enable}" = "yes" -o "${abi_gpu_enable}" = "auto"; then

    # Banner
    AC_MSG_NOTICE([Initializing GPU support])
    AC_MSG_CHECKING([which kind of GPU we want])
    AC_MSG_RESULT([${abi_gpu_flavor}])

    # Look for prerequisites
    case "${abi_gpu_flavor}" in

      cuda*)
        _ABI_GPU_INIT_CUDA
        abi_gpu_has_cc="${abi_gpu_cuda_has_cc}"
        abi_gpu_has_fft="${abi_gpu_cuda_has_fft}"
        abi_gpu_has_cublas="${abi_gpu_cuda_has_cublas}"
        abi_gpu_has_incs="${abi_gpu_cuda_has_incs}"
        abi_gpu_has_libs="${abi_gpu_cuda_has_libs}"
        abi_gpu_has_linalg="${abi_gpu_cuda_has_linalg}"
        if test "${abi_gpu_has_cc}" = "yes" -a \
                "${abi_gpu_has_incs}" = "yes" -a \
                "${abi_gpu_has_libs}" = "yes"; then
          abi_gpu_complete="yes"
        else
          abi_gpu_complete="no"
        fi
        ;;

      hip*)
        _ABI_GPU_INIT_HIP
        abi_gpu_has_cc="${abi_gpu_hip_has_cc}"
        abi_gpu_has_fft="${abi_gpu_hip_has_fft}"
        abi_gpu_has_cublas="${abi_gpu_hip_has_cublas}"
        abi_gpu_has_incs="${abi_gpu_hip_has_incs}"
        abi_gpu_has_libs="${abi_gpu_hip_has_libs}"
        abi_gpu_has_linalg="${abi_gpu_hip_has_linalg}"
        if test "${abi_gpu_has_cc}" = "yes" -a \
                "${abi_gpu_has_incs}" = "yes" -a \
                "${abi_gpu_has_libs}" = "yes"; then
          abi_gpu_complete="yes"
        else
          abi_gpu_complete="no"
        fi
        ;;

    esac

  else

    AC_MSG_NOTICE([GPU support disabled from command-line])
    abi_gpu_flavor="none"

  fi # abi_gpu_enable

  # Enable substitution
  AC_SUBST(abi_gpu_fcflags)
  AC_SUBST(abi_gpu_ldflags)
  AC_SUBST(abi_gpu_flavor)
  AC_SUBST(abi_gpu_incs)
  AC_SUBST(abi_gpu_libs)
  AC_SUBST(abi_gpu_arch)
]) # ABI_GPU_INIT



                    ########################################



# ABI_GPU_DETECT()
# ----------------
#
# Sets all variables needed to handle the GPU libraries.
#
AC_DEFUN([ABI_GPU_DETECT],[
  # Initial setup
  abi_gpu_serial="no"
  abi_gpu_precision=`echo "${abi_gpu_flavor}" | cut -d- -f2`
  test "${abi_gpu_precision}" = "" && abi_gpu_precision="single"

  # Display user requests
  AC_MSG_CHECKING([whether to activate GPU support])
  AC_MSG_RESULT([${abi_gpu_enable}])

  # Look for GPU libraries
  if test "${abi_gpu_enable}" = "yes"; then

    # Check whether we have a working gpu environment
    AC_MSG_CHECKING([for the requested GPU support])
    AC_MSG_RESULT([${abi_gpu_flavor}])

    case "${abi_gpu_flavor}" in

      cuda*)
        _ABI_GPU_CHECK_CUDA
        abi_gpu_serial="${abi_gpu_cuda_serial}"
        if test "${abi_gpu_serial}" = "yes"; then
          AC_DEFINE([HAVE_GPU_CUDA],1,[Define to 1 if you have the Cuda library.])
          if test "${abi_gpu_cuda_old}" = "yes"; then
            AC_DEFINE([HAVE_GPU_CUDA3],1,[Define to 1 if you have a Cuda version < 4.])
          fi
          if test "${abi_gpu_cuda_version_10}" = "yes"; then
            AC_DEFINE([HAVE_GPU_CUDA10],1,[Define to 1 if you have a Cuda version >= 10 (for nvtx v3).])
          fi
          case "${abi_gpu_precision}" in
            single)
              AC_DEFINE(HAVE_GPU_CUDA_SP,1,[Define to 1 if you want to perform single-precision Cuda calculations.])
              ;;
            double)
              AC_DEFINE(HAVE_GPU_CUDA_DP,1,[Define to 1 if you want to perform double-precision Cuda calculations.])
              ;;
          esac
          abi_gpu_fcflags="${abi_gpu_cuda_fcflags}"
          abi_gpu_ldflags="${abi_gpu_cuda_ldflags}"
          abi_gpu_incs="${abi_gpu_cuda_incs}"
          abi_gpu_libs="${abi_gpu_cuda_libs}"
        fi
        ;;

      hip*)
        _ABI_GPU_CHECK_HIP
        abi_gpu_serial="${abi_gpu_hip_serial}"
        if test "${abi_gpu_serial}" = "yes"; then
          AC_DEFINE([HAVE_GPU_HIP],1,[Define to 1 if you have the HIP library.])
        fi
        abi_gpu_fcflags="${abi_gpu_hip_fcflags}"
        abi_gpu_ldflags="${abi_gpu_hip_ldflags}"
        abi_gpu_incs="${abi_gpu_hip_incs}"
        abi_gpu_libs="${abi_gpu_hip_libs}"
        ;;
    esac

    if test "${abi_gpu_serial}" = "no"; then
      AC_MSG_ERROR([GPU support is broken])
    fi

  fi

  # Transmit serial status to the source code
  if test "${abi_gpu_serial}" = "yes"; then
    AC_DEFINE([HAVE_GPU],1,[Define to 1 if you have a GPU library.])
    AC_DEFINE([HAVE_GPU_SERIAL],1,[Define to 1 if you have a serial GPU library.])
    abi_gpu_flavor="${abi_gpu_flavor}"
  fi

  # Transmit the possible use of NVTX/ROCTX
  if test "${abi_gpu_nvtx_v3}" = "yes" -o "${abi_gpu_roctx}" = "yes"; then
    AC_DEFINE([HAVE_GPU_MARKERS],1,[Define to 1 if you have library for GPU range markers.])
    if test "${abi_gpu_nvtx_v3}" = "yes"; then
        AC_DEFINE([HAVE_GPU_MARKERS_NVTX],1,[Define to 1 if you have CUDA library for GPU range markers.])
    fi
    if test "${abi_gpu_roctx}" = "yes"; then
        AC_DEFINE([HAVE_GPU_MARKERS_ROCTX],1,[Define to 1 if you have ROCTX library for GPU range markers.])
    fi
  fi

  # Output final flavor
  if test "${abi_gpu_enable}" = "yes"; then
    AC_MSG_CHECKING([for the actual GPU support])
    AC_MSG_RESULT([${abi_gpu_flavor}])
  fi

  # FIXME: Update GPU libraries and includes (in fcflags)
  sd_gpu_libs="${sd_gpu_libs} ${abi_gpu_libs}"
  sd_gpu_fcflags="${sd_gpu_fcflags} ${abi_gpu_incs}"

  # Inform Automake
  AM_CONDITIONAL(DO_BUILD_17_GPU_TOOLBOX,[test "${abi_gpu_flavor}" != "none" -o "${abi_gpu_markers_enable}" = "yes"])
  AM_CONDITIONAL(DO_BUILD_46_MANAGE_CUDA,[test "${abi_gpu_flavor}" = "cuda-double" -o "${abi_gpu_flavor}" = "cuda-single"])
  AM_CONDITIONAL(DO_BUILD_NVTX,[test "${abi_gpu_nvtx_v3}" = "yes"])

]) # ABI_GPU_DETECT

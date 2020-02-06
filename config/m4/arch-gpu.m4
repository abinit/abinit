# -*- Autoconf -*-
#
# Copyright (C) 2009-2020 ABINIT Group (Yann Pouillon)
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
  abi_gpu_cuda_mpi="no"
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
  AC_CHECK_HEADERS([cuda_runtime_api.h cufft.h cublas.h])

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
  abi_gpu_cuda_has_incs="no"
  abi_gpu_cuda_has_libs="no"
  abi_gpu_cuda_has_linalg="no"
  abi_gpu_cuda_has_runtime="no"
  abi_gpu_cuda_libdir=""
  abi_gpu_cuda_incs="${GPU_CPPFLAGS}"
  abi_gpu_cuda_libs="${GPU_LIBS}"
  abi_gpu_cuda_root="${abi_gpu_prefix}"

  # Make use of the CUDA_ROOT environment variable
  if test "${abi_gpu_cuda_root}" = ""; then
    abi_gpu_cuda_root="${CUDA_ROOT}"
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
    if test -x "${abi_gpu_cuda_root}/bin/nvcc"; then
      abi_gpu_cuda_has_cc="yes"
      NVCC="${abi_gpu_cuda_root}/bin/nvcc"
    fi
    if test "${NVCC}" = ""; then
      AC_MSG_RESULT([none found])
    else
      AC_MSG_RESULT([${NVCC}])
    fi

    # Headers
    AC_MSG_CHECKING([for Cuda headers])
    abi_result=""
    if test -s "${abi_gpu_cuda_root}/include/cuda_runtime_api.h"; then
      if test "${GPU_CPPFLAGS}" = ""; then
        abi_gpu_cuda_incs="-I${abi_gpu_cuda_root}/include"
      fi
      abi_gpu_cuda_has_incs="yes"
      abi_result="${abi_result} run-time"
    fi
    if test -s "${abi_gpu_cuda_root}/include/cufft.h"; then
      abi_result="${abi_result} fft"
    fi
    if test -s "${abi_gpu_cuda_root}/include/cublas.h"; then
      abi_result="${abi_result} blas"
    fi
    if test -s "${abi_gpu_cuda_root}/SDK/C/common/inc/cutil.h"; then
      if test "${GPU_CPPFLAGS}" = ""; then
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
    else
      abi_gpu_cuda_libdir="${abi_gpu_cuda_root}/lib"
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
      if test "${GPU_LIBS}" = ""; then
        abi_gpu_cuda_libs="-lcudart"
      fi
      abi_gpu_cuda_has_libs="yes"
      abi_gpu_cuda_has_runtime="yes"
      abi_result="${abi_result} run-time"
    fi
    if test "${abi_gpu_cuda_has_libs}" = "yes"; then
      if test -e "${abi_gpu_cuda_libdir}/libcufft.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcufft ${abi_gpu_cuda_libs}"
        fi
        abi_gpu_cuda_has_fft="yes"
        abi_result="${abi_result} fft"
      fi
      if test -e "${abi_gpu_cuda_libdir}/libcublas.${abi_so_ext}"; then
        if test "${GPU_LIBS}" = ""; then
          abi_gpu_cuda_libs="-lcublas ${abi_gpu_cuda_libs}"
        fi
        abi_gpu_cuda_has_linalg="yes"
        abi_result="${abi_result} blas"
      fi
      if test "${GPU_LIBS}" = ""; then
        abi_gpu_cuda_libs="-L${abi_gpu_cuda_libdir} ${abi_gpu_cuda_libs}"
      fi
    fi
    if test -s "${abi_gpu_cuda_root}/SDK/C/lib/libcutil.a"; then
      if test "${GPU_LIBS}" = ""; then
        abi_gpu_cuda_libs="-L${abi_gpu_cuda_root}/SDK/C/lib -lcutil ${abi_gpu_cuda_libs}"
      fi
      abi_result="${abi_result} sdk"
    fi
    if test "${abi_result}" = ""; then
      abi_result="none"
    fi
    AC_MSG_RESULT([${abi_result}])
    if test "${GPU_LIBS}" = ""; then
      abi_gpu_cuda_libs="${abi_gpu_cuda_libs} -lcuda"
    fi
    if test "${abi_gpu_cuda_has_common}" = "no"; then
      AC_MSG_WARN([could not find libcuda.${abi_so_ext}])
    fi

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
]) # _ABI_GPU_INIT_CUDA



                    ########################################



# ABI_GPU_INIT()
# --------------
#
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
  abi_gpu_has_cc="no"
  abi_gpu_has_fft="no"
  abi_gpu_has_incs="no"
  abi_gpu_has_libs="no"
  abi_gpu_has_linalg="no"
  abi_gpu_usable="no"
  abi_gpu_fcflags=""
  abi_gpu_ldflags=""
  abi_gpu_flavor="${sd_gpu_flavor}"
  abi_gpu_incs="${GPU_CPPFLAGS}"
  abi_gpu_libs="${GPU_LIBS}"

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
  abi_gpu_mpi="no"
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
        abi_gpu_mpi="${abi_gpu_cuda_mpi}"
        if test "${abi_gpu_serial}" = "yes"; then
          AC_DEFINE([HAVE_GPU_CUDA],1,[Define to 1 if you have the Cuda library.])
          if test "${abi_gpu_cuda_old}" = "yes"; then
            AC_DEFINE([HAVE_GPU_CUDA3],1,[Define to 1 if you have a Cuda version < 4.])
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

  # Transmit MPI status to the source code
  if test "${abi_gpu_mpi}" = "yes"; then
    AC_DEFINE([HAVE_GPU_MPI],1,[Define to 1 if you have a MPI-aware GPU library.])
  fi

  # Output final flavor
  if test "${abi_gpu_enable}" = "yes"; then
    AC_MSG_CHECKING([for the actual GPU support])
    AC_MSG_RESULT([${abi_gpu_flavor}])
  fi

  # FIXME: Update GPU libraries
  sd_gpu_libs="${sd_gpu_libs} ${abi_gpu_libs}"

  # Inform Automake
  AM_CONDITIONAL(DO_BUILD_17_GPU_TOOLBOX,[test "${abi_gpu_flavor}" != "none"])
  AM_CONDITIONAL(DO_BUILD_52_MANAGE_CUDA,[test "${abi_gpu_flavor}" = "cuda-double" -o "${abi_gpu_flavor}" = "cuda-single"])
]) # ABI_GPU_DETECT

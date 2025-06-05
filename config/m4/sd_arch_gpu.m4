## Copyright (C) 2019-2025 ABINIT group (Yann Pouillon)

#
# GPU support for Steredeg
#


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_GPU_INIT], [
  # Init
  sd_gpu_cppflags=""
  sd_gpu_cflags=""
  sd_gpu_fcflags=""
  sd_gpu_ldflags=""
  sd_gpu_libs=""
  sd_gpu_enable=""
  sd_gpu_markers_enable=""
  sd_gpu_init="unknown"
  sd_gpu_ok="unknown"
  sd_gpu_prefix=""
  sd_cuda_prefix=""
  sd_cuda_enable=""
  sd_cuda_init=""
  sd_rocm_prefix=""
  sd_rocm_enable=""
  sd_rocm_init=""

  # Set adjustable parameters
  sd_gpu_options="$1"
  sd_gpu_libs_def="$2"
  sd_gpu_cppflags_def="$3"
  sd_gpu_cflags_def="$4"
  sd_gpu_cxxflags_def="$5"
  sd_gpu_fcflags_def="$6"
  sd_gpu_ldflags_def="$7"
  sd_gpu_enable_def=""
  sd_gpu_flavor_def=""
  sd_gpu_markers_enable_def=""
  sd_gpu_policy=""
  sd_gpu_status=""

  # Process options
  for kwd in ${sd_gpu_options}; do
    case "${kwd}" in
      auto|no|yes)
        sd_gpu_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_gpu_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_gpu_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg GPU option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_gpu_enable_def}" && sd_gpu_enable_def="no"
  test -z "${sd_gpu_flavor_def}" && sd_gpu_flavor_def="cuda-double"
  test -z "${sd_gpu_markers_enable_def}" && sd_gpu_markers_enable_def="no"
  test -z "${sd_gpu_policy}" && sd_gpu_policy="warn"
  test -z "${sd_cuda_status}" && sd_gpu_status="${sd_cuda_status}"
  test -z "${sd_rocm_status}" && sd_gpu_status="${sd_rocm_status}"
  test -z "${sd_gpu_status}" && sd_gpu_status="optional"
  # FIXME: improve the setting mechanism
  #test -z "${sd_gpu_libs_def}" && sd_gpu_libs_def="-lopencl"

  # Declare main configure option
  AC_ARG_WITH([gpu],
    [AS_HELP_STRING(
      [--with-gpu],
      [ Install prefix of a GPU SDK . You may use --with-gpu without argument to force GPU detection, in which case detection failures will result in errors, and --without-gpu to disable GPU support completely.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_gpu_enable="${withval}"
        sd_gpu_init="yon"
      else
        sd_gpu_prefix="${withval}"
        sd_gpu_enable="yes"
        sd_gpu_init="dir"
      fi],
    [ sd_gpu_enable="${sd_gpu_enable_def}"; sd_gpu_init="def"])

  AC_ARG_WITH([cuda],
    [AS_HELP_STRING(
      [--with-cuda],
      [Install prefix of NVIDIA CUDA SDK. You may use --with-cuda without argument to force GPU detection, in which case detection failures will result in errors.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_cuda_enable="${withval}"
        sd_cuda_init="yon"
        sd_gpu_flavor="cuda-double"
      else
        sd_cuda_prefix="${withval}"
        sd_cuda_enable="yes"
        sd_cuda_init="dir"
        sd_gpu_flavor_def="cuda-double"
        sd_gpu_flavor="${sd_gpu_flavor_def}"
      fi],
    [ sd_cuda_enable="${sd_gpu_enable_def}"; sd_cuda_init="def"])

  AC_ARG_WITH([rocm],
    [AS_HELP_STRING(
      [--with-rocm],
      [Install prefix of AMD ROCM SDK. You may use --with-rocm without argument to force GPU detection, in which case detection failures will result in errors.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_rocm_enable="${withval}"
        sd_rocm_init="yon"
        sd_gpu_flavor="hip-double"
      else
        sd_rocm_prefix="${withval}"
        sd_rocm_enable="yes"
        sd_rocm_init="dir"
        sd_gpu_flavor_def="hip-double"
        sd_gpu_flavor="${sd_gpu_flavor_def}"
      fi],
    [ sd_rocm_enable="${sd_gpu_enable_def}"; sd_rocm_init="def"])

  if test "${sd_cuda_enable}" != "${sd_gpu_enable_def}"; then
    sd_gpu_prefix="${sd_cuda_prefix}"
    sd_gpu_enable="${sd_cuda_enable}"
    sd_gpu_init="${sd_cuda_init}"
  fi

  if test "${sd_rocm_enable}" != "${sd_gpu_enable_def}"; then
    sd_gpu_prefix="${sd_rocm_prefix}"
    sd_gpu_enable="${sd_rocm_enable}"
    sd_gpu_init="${sd_rocm_init}"
  fi

  # Declare main configure option
  AC_ARG_WITH([gpu_markers],
    [AS_HELP_STRING(
      [--with-gpu-markers],
      [Enable GPU markers such as NVTX for NVIDIA CUDA or ROCTX for AMD ROCm.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_gpu_markers_enable="${withval}"
      else
        sd_gpu_markers_enable="no"
      fi],
    [ sd_gpu_markers_enable="${sd_gpu_markers_enable_def}"; sd_gpu_markers_init="def"])

  # Declare flavor option
  sd_gpu_flavors_supported="cuda-double cuda-single hip-double"
  AC_ARG_WITH([gpu-flavor],
    [AS_HELP_STRING(
      [--with-gpu-flavor],
      [GPU flavor to select. Use it without argument to get a list of supported flavors.])],
    [ case "${withval}" in
        no|yes)
          AC_MSG_NOTICE([])
          AC_MSG_NOTICE([available GPU flavors:])
          AC_MSG_NOTICE([])
          for tmp_gpu_flavor in ${sd_gpu_flavors_supported}; do
            AC_MSG_NOTICE([  * ${tmp_gpu_flavor}])
          done
          AC_MSG_NOTICE([])
          AC_MSG_ERROR([please select a valid GPU flavor])
          ;;
        *)
          tmp_gpu_flavor_ok="no"
          for tmp_gpu_flavor in ${sd_gpu_flavors_supported}; do
            if test "${withval}" = "${tmp_gpu_flavor}"; then
              tmp_gpu_flavor_ok="yes"
              break
            fi
          done
          if test "${tmp_gpu_flavor_ok}" = "yes"; then
            sd_gpu_flavor="${withval}"
            sd_gpu_flavor_init="sel"
          else
            AC_MSG_ERROR([invalid GPU flavor: '${withval}'])
          fi
          unset tmp_gpu_flavor
          unset tmp_gpu_flavor_ok
          ;;
      esac],
    [ sd_gpu_flavor="${sd_gpu_flavor_def}"
      sd_gpu_flavor_init="def"])

  # Declare environment variables
  AC_ARG_VAR([GPU_CPPFLAGS], [C preprocessing flags for GPU, deprecated.])
  AC_ARG_VAR([GPU_CFLAGS], [C flags for GPU, deprecated.])
  AC_ARG_VAR([GPU_CXXFLAGS], [C++ flags for GPU, deprecated.])
  AC_ARG_VAR([GPU_FCFLAGS], [Fortran flags for GPU, deprecated.])
  AC_ARG_VAR([GPU_FFLAGS], [Fortran flags for GPU, deprecated (better use GPU_FCFLAGS).])
  AC_ARG_VAR([GPU_LDFLAGS], [Linker flags for GPU, deprecated.])
  AC_ARG_VAR([GPU_LIBS], [Library flags for GPU, deprecated.])

  # Detect use of environment variables
  if test "${sd_gpu_enable}" = "yes" -o "${sd_gpu_enable}" = "auto"; then
    if test "${sd_cuda_enable}" = "yes"; then
      tmp_gpu_vars="${CUDA_CPPFLAGS}${CUDA_CFLAGS}${CUDA_CXXFLAGS}${CUDA_FFLAGS}${CUDA_FCFLAGS}${CUDA_LDFLAGS}${CUDA_LIBS}"
      if test "${sd_gpu_init}" = "def" -a ! -z "${tmp_gpu_vars}"; then
        sd_gpu_enable="yes"
        sd_gpu_init="env"
      fi
    elif test "${sd_rocm_enable}" = "yes"; then
      tmp_gpu_vars="${ROCM_CPPFLAGS}${ROCM_CFLAGS}${ROCM_CXXFLAGS}${ROCM_FFLAGS}${ROCM_FCFLAGS}${ROCM_LDFLAGS}${ROCM_LIBS}"
      if test "${sd_gpu_init}" = "def" -a ! -z "${tmp_gpu_vars}"; then
        sd_gpu_enable="yes"
        sd_gpu_init="env"
      fi
    else
      tmp_gpu_vars="${GPU_CPPFLAGS}${GPU_CFLAGS}${GPU_CXXFLAGS}${GPU_FFLAGS}${GPU_FCFLAGS}${GPU_LDFLAGS}${GPU_LIBS}"
      if test "${sd_gpu_init}" = "def" -a ! -z "${tmp_gpu_vars}"; then
        sd_gpu_enable="yes"
        sd_gpu_init="env"
      fi
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_GPU_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_gpu_enable}" = "yes" -o "${sd_gpu_enable}" = "auto"; then

    # Set GPU-specific flags
    case "${sd_gpu_init}" in

      def|yon)
        sd_gpu_cppflags="${sd_gpu_cppflags_def}"
        sd_gpu_cflags="${sd_gpu_cflags_def}"
        sd_gpu_cxxflags="${sd_gpu_cxxflags_def}"
        sd_gpu_fcflags="${sd_gpu_fcflags_def}"
        sd_gpu_ldflags="${sd_gpu_ldflags_def}"
        sd_gpu_libs="${sd_gpu_libs_def}"
        ;;

      dir)
        sd_gpu_cppflags="-I${sd_gpu_prefix}/include"
        sd_gpu_cflags="${sd_gpu_cflags_def}"
        sd_gpu_cxxflags="${sd_gpu_cxxflags_def}"
        sd_gpu_fcflags="${sd_gpu_fcflags_def} -I${sd_gpu_prefix}/include"
        sd_gpu_ldflags="${sd_gpu_ldflags_def}"
        sd_gpu_libs="-L${sd_gpu_prefix}/lib ${sd_gpu_libs_def}"
        ;;

      env)
        sd_gpu_cppflags="${sd_gpu_cppflags_def}"
        sd_gpu_cflags="${sd_gpu_cflags_def}"
        sd_gpu_cxxflags="${sd_gpu_cxxflags_def}"
        sd_gpu_fcflags="${sd_gpu_fcflags_def}"
        sd_gpu_ldflags="${sd_gpu_ldflags_def}"
        sd_gpu_libs="${sd_gpu_libs_def}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for GPU: '${sd_gpu_init}'])
        ;;

    esac

    if test "${sd_cuda_enable}" = "yes"; then

      test ! -z "${CUDA_CPPFLAGS}" && sd_gpu_cppflags="${CUDA_CPPFLAGS}"
      test ! -z "${CUDA_CFLAGS}" && sd_gpu_cflags="${CUDA_CFLAGS}"
      if test "${sd_gpu_enable_cxx}" = "yes"; then
        test ! -z "${CUDA_CXXFLAGS}" && sd_gpu_cxxflags="${CUDA_CXXFLAGS}"
      fi
      if test "${sd_gpu_enable_fc}" = "yes"; then
        test ! -z "${CUDA_FFLAGS}" && sd_gpu_fcflags="${CUDA_FFLAGS}"
        test ! -z "${CUDA_FCFLAGS}" && sd_gpu_fcflags="${CUDA_FCFLAGS}"
      fi
      test ! -z "${CUDA_LDFLAGS}" && sd_gpu_ldflags="${CUDA_LDFLAGS}"
      test ! -z "${CUDA_LIBS}" && sd_gpu_libs="${CUDA_LIBS}"

    elif test "${sd_rocm_enable}" = "yes"; then

      test ! -z "${ROCM_CPPFLAGS}" && sd_gpu_cppflags="${ROCM_CPPFLAGS}"
      test ! -z "${ROCM_CFLAGS}" && sd_gpu_cflags="${ROCM_CFLAGS}"
      if test "${sd_gpu_enable_cxx}" = "yes"; then
        test ! -z "${ROCM_CXXFLAGS}" && sd_gpu_cxxflags="${ROCM_CXXFLAGS}"
      fi
      if test "${sd_gpu_enable_fc}" = "yes"; then
        test ! -z "${ROCM_FFLAGS}" && sd_gpu_fcflags="${ROCM_FFLAGS}"
        test ! -z "${ROCM_FCFLAGS}" && sd_gpu_fcflags="${ROCM_FCFLAGS}"
      fi
      test ! -z "${ROCM_LDFLAGS}" && sd_gpu_ldflags="${ROCM_LDFLAGS}"
      test ! -z "${ROCM_LIBS}" && sd_gpu_libs="${ROCM_LIBS}"

    else

      test ! -z "${GPU_CPPFLAGS}" && sd_gpu_cppflags="${GPU_CPPFLAGS}"
      test ! -z "${GPU_CFLAGS}" && sd_gpu_cflags="${GPU_CFLAGS}"
      if test "${sd_gpu_enable_cxx}" = "yes"; then
        test ! -z "${GPU_CXXFLAGS}" && sd_gpu_cxxflags="${GPU_CXXFLAGS}"
      fi
      if test "${sd_gpu_enable_fc}" = "yes"; then
        test ! -z "${GPU_FFLAGS}" && sd_gpu_fcflags="${GPU_FFLAGS}"
        test ! -z "${GPU_FCFLAGS}" && sd_gpu_fcflags="${GPU_FCFLAGS}"
      fi
      test ! -z "${GPU_LDFLAGS}" && sd_gpu_ldflags="${GPU_LDFLAGS}"
      test ! -z "${GPU_LIBS}" && sd_gpu_libs="${GPU_LIBS}"

    fi


    # Toggle use GPU unified memory feature, specific to NVHPC with OpenMP offload, on NVIDIA GPUs
    # This flag is supported for OpenMP since NVHPC v24.3
    # On AMD GPU, this feature seem to be controled using env variables.
    if test "${sd_gpu_nvidia_unified_memory_enable}" == "yes"; then
      if test "${abi_fc_vendor}" != "nvhpc" -o "${abi_openmp_offload_enable}" = "yes"; then
        AC_MSG_ERROR([Unified memory setting is only supported with NVHPC SDK and OpenMP offload enabled !])
      fi
      nvhpc_version=`echo ${abi_fc_version} | sed 's/\.//;s/-//'`
      if test ${nvhpc_version} -lt 2430; then
        AC_MSG_ERROR([Unified memory setting is only supported since NVHPC SDK version 24.3. Your NVHPC is too old.])
      fi
      gpu_unified_flag="-gpu=mem:unified"
      # Use older flag for NVHPC 24.3, deprecated in newer versions
      if test ${nvhpc_version} -eq 2430; then
        gpu_unified_flag="-gpu=unified"
      fi
      sd_gpu_cflags="${sd_gpu_cflags} {gpu_unified_flag}"
      sd_gpu_cxxflags="${sd_gpu_cxxflags} {gpu_unified_flag}"
      sd_gpu_fcflags="${sd_gpu_fcflags} {gpu_unified_flag}"
      sd_gpu_ldflags="${sd_gpu_ldflags} {gpu_unified_flag}"
    fi
  fi

  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_GPU_DUMP_CONFIG
  fi

  # Export configuration
  AC_SUBST(sd_gpu_options)
  AC_SUBST(sd_gpu_enable_def)
  AC_SUBST(sd_gpu_policy)
  AC_SUBST(sd_gpu_prefix)
  AC_SUBST(sd_gpu_status)
  AC_SUBST(sd_gpu_enable)
  AC_SUBST(sd_gpu_markers_enable)
  AC_SUBST(sd_gpu_init)
  AC_SUBST(sd_gpu_ok)
  AC_SUBST(sd_gpu_cppflags)
  AC_SUBST(sd_gpu_cflags)
  AC_SUBST(sd_gpu_cxxflags)
  AC_SUBST(sd_gpu_fcflags)
  AC_SUBST(sd_gpu_ldflags)
  AC_SUBST(sd_gpu_libs)

  # Clean-up
  unset tmp_gpu_vars
]) # SD_GPU_INIT


AC_DEFUN([SD_GPU_DETECT], [
  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
  _SD_GPU_DUMP_CONFIG
  fi

  # Check compilers or APIs
  if test "${sd_gpu_enable}" = "yes" -o "${sd_gpu_enable}" = "auto"; then

    # FIXME: search for GPU implementations
    AC_MSG_WARN([GPU detection not implemented!])

    # Validate implementation status
    tmp_gpu_ok="no"

    # Take decision according to policy
    if test "${sd_gpu_ok}" = "yes"; then
      sd_gpu_enable="yes"
    else
      if test "${sd_gpu_enable}" = "yes"; then
        case "${sd_gpu_policy}" in
          fail)
            AC_MSG_FAILURE([GPU support does not work])
            ;;
          skip)
            sd_gpu_enable="no"
            ;;
          warn)
            sd_gpu_enable="no"
            AC_MSG_WARN([GPU support does not work and has been disabled])
	    ;;
        esac
      else
        sd_gpu_enable="no"
      fi
    fi

  fi

  # Make GPU status available to the source code
  if test "${sd_gpu_enable}" = "yes" -a "${sd_gpu_ok}" = "yes"; then
    AC_DEFINE([HAVE_GPU], 1,
      [Define to 1 if you have a working GPU installation.])
  fi
]) # SD_GPU_DETECT


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Private macros
#


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_GPU_CHECK_CONFIG], [
  # Main trigger must be yes, no, or auto
  tmp_gpu_invalid="no"
  if test "${sd_gpu_enable}" != "auto" -a \
          "${sd_gpu_enable}" != "no" -a \
          "${sd_gpu_enable}" != "yes"; then
    case "${sd_gpu_policy}" in
      fail)
        AC_MSG_ERROR([invalid trigger value: sd_gpu_enable = '${sd_gpu_enable}'])
        ;;
      skip)
        tmp_gpu_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid trigger value: sd_gpu_enable = '${sd_gpu_enable}'])
        tmp_gpu_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_gpu_invalid}" = "yes"; then
    if test "${sd_gpu_status}" = "required"; then
      sd_gpu_enable="yes"
    else
      sd_gpu_enable="no"
    fi
    tmp_gpu_invalid="no"
    AC_MSG_NOTICE([setting sd_gpu_enable to '${sd_gpu_enable}'])
  fi

  # Check consistency between trigger value and package status
  tmp_gpu_invalid="no"
  if test "${sd_gpu_status}" = "implicit" -o \
          "${sd_gpu_status}" = "required"; then
    if test "${sd_gpu_enable}" = "no"; then
      case "${sd_gpu_policy}" in
        fail)
          AC_MSG_ERROR([The GPU package is required and cannot be disabled.])
          ;;
        skip)
          tmp_gpu_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The GPU package is required and cannot be disabled.])
          tmp_gpu_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_gpu_enable}" = "auto"; then
      AC_MSG_NOTICE([setting GPU trigger to yes])
      sd_gpu_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_gpu_invalid}" = "yes"; then
    case "${sd_gpu_status}" in
      implicit|required)
        sd_gpu_enable="yes"
        ;;
      optional)
        sd_gpu_enable="no"
        ;;
    esac
    tmp_gpu_invalid="no"
    AC_MSG_NOTICE([setting sd_gpu_enable to '${sd_gpu_enable}'])
  fi

  # When using environment variables, triggers must be set to yes
  if test "${sd_gpu_init}" = "env" -a "${sd_gpu_enable}" = "no"; then
    if test "${sd_gpu_policy}" != "skip"; then
      AC_MSG_WARN([GPU environment variables will be ignored])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_gpu_status}" = "implicit"; then
    sd_gpu_enable="yes"
    if test "${sd_gpu_cppflags}" != ""; then
      sd_gpu_cppflags=""
      AC_MSG_NOTICE([resetting GPU CPP flags (implicit package)])
    fi
    if test "${sd_gpu_cflags}" != ""; then
      sd_gpu_cflags=""
      AC_MSG_NOTICE([resetting GPU C flags (implicit package)])
    fi
    if test "${sd_gpu_cxxflags}" != ""; then
      sd_gpu_cxxflags=""
      AC_MSG_NOTICE([resetting GPU C++ flags (implicit package)])
    fi
    if test "${sd_gpu_fcflags}" != ""; then
      sd_gpu_fcflags=""
      AC_MSG_NOTICE([resetting GPU Fortran flags (implicit package)])
    fi
    if test "${sd_gpu_ldflags}" != ""; then
      sd_gpu_ldflags=""
      AC_MSG_NOTICE([resetting GPU linker flags (implicit package)])
    fi
    if test "${sd_gpu_libs}" != ""; then
      sd_gpu_libs=""
      AC_MSG_NOTICE([resetting GPU library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_gpu_enable}" = "no"; then
    sd_gpu_cppflags=""
    sd_gpu_cflags=""
    sd_gpu_cxxflags=""
    sd_gpu_fcflags=""
    sd_gpu_ldflags=""
    sd_gpu_libs=""
    sd_gpu_ok="no"
  fi

  # Clean-up
  unset tmp_gpu_invalid
]) # _SD_GPU_CHECK_CONFIG


AC_DEFUN([_SD_GPU_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable GPU])
  AC_MSG_RESULT([${sd_gpu_enable}])
  if test "${sd_gpu_enable}" != "no"; then
    AC_MSG_CHECKING([how GPU parameters have been set])
    AC_MSG_RESULT([${sd_gpu_init}])
    AC_MSG_CHECKING([for GPU C preprocessing flags])
    AC_MSG_RESULT([${sd_gpu_cppflags}])
    AC_MSG_CHECKING([for GPU C flags])
    AC_MSG_RESULT([${sd_gpu_cflags}])
    AC_MSG_CHECKING([for GPU C++ flags])
    AC_MSG_RESULT([${sd_gpu_cxxflags}])
    AC_MSG_CHECKING([for GPU Fortran flags])
    AC_MSG_RESULT([${sd_gpu_fcflags}])
    AC_MSG_CHECKING([for GPU linker flags])
    AC_MSG_RESULT([${sd_gpu_ldflags}])
    AC_MSG_CHECKING([for GPU library flags])
    AC_MSG_RESULT([${sd_gpu_libs}])
  fi
]) # _SD_GPU_DUMP_CONFIG

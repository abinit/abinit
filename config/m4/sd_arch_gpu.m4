## Copyright (C) 2019 Yann Pouillon

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
  sd_gpu_init="unknown"
  sd_gpu_ok="unknown"

  # Set adjustable parameters
  sd_gpu_options="$1"
  sd_gpu_libs_def="$2"
  sd_gpu_cppflags_def="$3"
  sd_gpu_cflags_def="$4"
  sd_gpu_cxxflags_def="$5"
  sd_gpu_fcflags_def="$6"
  sd_gpu_ldflags_def="$7"
  sd_gpu_enable_def=""
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
  test -z "${sd_gpu_policy}" && sd_gpu_policy="warn"
  test -z "${sd_gpu_status}" && sd_gpu_status="optional"
  test -z "${sd_gpu_libs_def}" && sd_gpu_libs_def="-lopencl"

  # Declare main configure option
  AC_ARG_WITH([gpu],
    [AS_HELP_STRING(
      [--with-gpu],
      [Install prefix of a GPU SDK (e.g. /usr/local). You may use --with-gpu without argument to force GPU detection, in which case detection failures will result in errors, and --without-gpu to disable GPU support completely.])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_gpu_enable="${withval}"
        sd_gpu_init="yon"
      else
        sd_gpu_enable="yes"
        sd_gpu_init="dir"
      fi],
    [ sd_gpu_enable="${sd_gpu_enable_def}"; sd_gpu_init="def"])

  # Declare flavor option
  sd_gpu_flavors_supported="cuda-double cuda-single"
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
            test "${withval}" = "${tmp_gpu_flavor}" && tmp_gpu_flavor_ok="yes"
            break
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
    [ sd_gpu_flavor=`echo "${sd_gpu_flavors_supported}" | cut -d' ' -f1`
      sd_gpu_flavor_init="def"])

  # Declare environment variables
  AC_ARG_VAR([GPU_CPPFLAGS], [C preprocessing flags for GPU.])
  AC_ARG_VAR([GPU_CFLAGS], [C flags for GPU.])
  AC_ARG_VAR([GPU_CXXFLAGS], [C++ flags for GPU.])
  AC_ARG_VAR([GPU_FCFLAGS], [Fortran flags for GPU.])
  AC_ARG_VAR([GPU_LDFLAGS], [Linker flags for GPU.])
  AC_ARG_VAR([GPU_LIBS], [Library flags for GPU.])

  # Detect use of environment variables
  if test "${sd_gpu_enable}" = "yes" -o "${sd_gpu_enable}" = "auto"; then
    tmp_gpu_vars="${GPU_CPPFLAGS}${GPU_CFLAGS${GPU_CXXFLAGS}}${GPU_FCFLAGS}${GPU_LDFLAGS}${GPU_LIBS}"
    if test "${sd_gpu_init}" = "def" -a ! -z "${tmp_gpu_vars}"; then
      sd_gpu_enable="yes"
      sd_gpu_init="env"
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
        sd_gpu_cppflags="-I${with_gpu}/include"
        sd_gpu_cflags="${sd_gpu_cflags_def}"
        sd_gpu_cxxflags="${sd_gpu_cxxflags_def}"
        sd_gpu_fcflags="${sd_gpu_fcflags_def} -I${with_gpu}/include"
        sd_gpu_ldflags="${sd_gpu_ldflags_def}"
        sd_gpu_libs="-L${with_gpu}/lib ${sd_gpu_libs_def}"
        ;;

      env)
        sd_gpu_cppflags="${sd_gpu_cppflags_def}"
        sd_gpu_cflags="${sd_gpu_cflags_def}"
        sd_gpu_cxxflags="${sd_gpu_cxxflags_def}"
        sd_gpu_fcflags="${sd_gpu_fcflags_def}"
        sd_gpu_ldflags="${sd_gpu_ldflags_def}"
        sd_gpu_libs="${sd_gpu_libs_def}"
        test ! -z "${GPU_CPPFLAGS}" && sd_gpu_cppflags="${GPU_CPPFLAGS}"
        test ! -z "${GPU_CFLAGS}" && sd_gpu_cflags="${GPU_CFLAGS}"
        if test "${sd_gpu_enable_cxx}" = "yes"; then
          test ! -z "${GPU_CXXFLAGS}" && sd_gpu_cxxflags="${GPU_CXXFLAGS}"
        fi
        if test "${sd_gpu_enable_fc}" = "yes"; then
          test ! -z "${GPU_FCFLAGS}" && sd_gpu_fcflags="${GPU_FCFLAGS}"
        fi
        test ! -z "${GPU_LDFLAGS}" && sd_gpu_ldflags="${GPU_LDFLAGS}"
        test ! -z "${GPU_LIBS}" && sd_gpu_libs="${GPU_LIBS}"
        ;;

      *)
        AC_MSG_ERROR([invalid init type for GPU: '${sd_gpu_init}'])
        ;;

    esac

  fi

  # Display configuration
  if test "${STEREDEG_CONFIG_BYPASS_CHECKS}" != "yes"; then
    _SD_GPU_DUMP_CONFIG
  fi

  # Export configuration
  AC_SUBST(sd_gpu_options)
  AC_SUBST(sd_gpu_enable_def)
  AC_SUBST(sd_gpu_policy)
  AC_SUBST(sd_gpu_status)
  AC_SUBST(sd_gpu_enable)
  AC_SUBST(sd_gpu_init)
  AC_SUBST(sd_gpu_ok)
  AC_SUBST(sd_gpu_cppflags)
  AC_SUBST(sd_gpu_cflags)
  AC_SUBST(sd_gpu_cxxflags)
  AC_SUBST(sd_gpu_fcflags)
  AC_SUBST(sd_gpu_ldflags)
  AC_SUBST(sd_gpu_libs)
  AC_SUBST(with_gpu)

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

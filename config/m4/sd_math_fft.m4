## Copyright (C) 2019 Yann Pouillon <devops@materialsevolution.es>

#
# Multi-flavor Fast Fourier Transform support
#


                    # ------------------------------------ #


AC_DEFUN([SD_FFT_INIT], [
  # Init
  sd_fft_cflags=""
  sd_fft_cppflags=""
  sd_fft_fcflags=""
  sd_fft_ldflags=""
  sd_fft_libs=""
  sd_fft_flavor=""
  sd_fft_init="unknown"
  sd_fft_ok="unknown"

  # Set adjustable parameters
  sd_fft_flavors="$1"
  sd_fft_policy="$3"
  sd_fft_status="$2"
  sd_fft_flavor_def=`echo "${sd_fft_flavors}" | ${AWK} '{print [$]1}'`

  # Make sure default values are correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFT_CHECK_DEFAULTS
  fi

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([fft-flavor],
    [AS_HELP_STRING([--with-fft-flavor],
      [FFT flavor to select among fftw3, and pfft])],
    [ tmp_fft_flavor_ok=`echo "${sd_fft_flavors}" | grep "${withval}"`
      if test "${tmp_fft_flavor_ok}" = ""; then
        AC_MSG_ERROR([invalid FFT flavor: '${withval}'])
      fi
      sd_fft_init="cmd"
      unset tmp_fft_flavor_ok],
    [ sd_fft_flavor="${sd_fft_flavor_def}"; sd_fft_init="def"])

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFT_CHECK_CONFIG
  fi

  # Display configuration
  _SD_FFT_DUMP_CONFIG

  # Init packages for permitted flavors
  for tmp_fft_flavor in ${sd_fft_flavors}; do
    tmp_fft_enable="no"
    test "${tmp_fft_flavor}" = "${sd_fft_flavor}" && tmp_fft_enable="yes"
    case "${tmp_fft_flavor}" in
      fftw3)
        SD_FFTW3_INIT([${tmp_fft_enable}], [], [], [-lfftw3], [optional], [skip])
        ;;
      pfft)
        SD_PFFT_INIT([${tmp_fft_enable}], [], [], [-lpfft], [optional], [skip])
        ;;
      *)
        AC_MSG_ERROR([unsupported FFT flavor: '${tmp_fft_flavor}'])
        ;;
    esac
  done

  # Export configuration
  AC_SUBST(sd_fft_cppflags)
  AC_SUBST(sd_fft_cflags)
  AC_SUBST(sd_fft_fcflags)
  AC_SUBST(sd_fft_ldflags)
  AC_SUBST(sd_fft_libs)
  AC_SUBST(sd_fft_flavor)
  AC_SUBST(sd_fft_init)
  AC_SUBST(sd_fft_ok)
  AC_SUBST(with_fft_flavor)
]) # SD_FFT_INIT


AC_DEFUN([SD_FFT_SELECT_FLAVOR], [
  AC_MSG_CHECKING([for the FFT flavor to use])
  AC_MSG_RESULT([${sd_fft_flavor}])
  case "${sd_fft_flavor}" in
    fftw3)
      SD_FFTW3_DETECT
      if test "${sd_fftw3_ok}" = "yes"; then
        sd_fft_cppflags="${sd_fftw3_cppflags}"
        sd_fft_cflags="${sd_fftw3_cflags}"
        sd_fft_fcflags="${sd_fftw3_fcflags}"
        sd_fft_ldflags="${sd_fftw3_ldflags}"
        sd_fft_libss="${sd_fftw3_libs}"
      else
        AC_MSG_ERROR([FFTW3 is not available
                    Please adjust configure options to point to a working FFTW3
                    installation or change the requested FFT flavor through the
                    --with-fft-flavor option.])
      fi
      ;;
    pfft)
      SD_PFFT_DETECT
      if test "${sd_pfft_ok}" = "yes"; then
        sd_fft_cppflags="${sd_pfft_cppflags}"
        sd_fft_cflags="${sd_pfft_cflags}"
        sd_fft_fcflags="${sd_pfft_fcflags}"
        sd_fft_ldflags="${sd_pfft_ldflags}"
        sd_fft_libss="${sd_pfft_libs}"
      else
        AC_MSG_ERROR([PFFT is not available
                    Please adjust configure options to point to a working PFFT
                    installation or change the requested FFT flavor through the
                    --with-fft-flavor option.])
      fi
      ;;
    *)
      AC_MSG_ERROR([unsupported FFT flavor: '${sd_fft_flavor}'])
      ;;
  esac
  sd_fft_ok="yes"
]) # SD_FFT_SELECT_FLAVOR


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_FFT_CHECK_CONFIG], [
  # Check consistency between trigger value and package status
  tmp_fft_invalid="no"
  tmp_fft_flavor_chk=`echo "${sd_fft_flavors}" | grep "${sd_fft_flavor}"`
  if test "${sd_fft_status}" = "implicit" -o \
          "${sd_fft_status}" = "required"; then
    if test "${tmp_fft_flavor_chk}" = ""; then
      case "${sd_fft_policy}" in
        fail)
          AC_MSG_ERROR([invalid FFT flavor: sd_fft_flavor = '${sd_fft_flavor}'])
          ;;
        skip)
          tmp_fft_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The FFT package is required and cannot be disabled])
          tmp_fft_invalid="yes"
          ;;
      esac
    fi
  fi

  # Fix wrong flavor value
  if test "${tmp_fft_invalid}" = "yes"; then
    sd_fft_flavor="${sd_fft_flavor_def}"
    AC_MSG_NOTICE([setting sd_fft_flavor to '${sd_fft_flavor}'])
  fi

  # Implicit status overrides everything
  if test "${sd_fft_status}" = "implicit"; then
    if test "${sd_fft_cppflags}" != ""; then
      sd_fft_cppflags=""
      AC_MSG_NOTICE([resetting FFT C preprocessing flags (implicit package)])
    fi
    if test "${sd_fft_cflags}" != ""; then
      sd_fft_cflags=""
      AC_MSG_NOTICE([resetting FFT C flags (implicit package)])
    fi
    if test "${sd_fft_fcflags}" != ""; then
      sd_fft_fcflags=""
      AC_MSG_NOTICE([resetting FFT Fortran flags (implicit package)])
    fi
    if test "${sd_fft_ldflags}" != ""; then
      sd_fft_ldflags=""
      AC_MSG_NOTICE([resetting FFT linker flags (implicit package)])
    fi
    if test "${sd_fft_libs}" != ""; then
      sd_fft_libs=""
      AC_MSG_NOTICE([resetting FFT library flags (implicit package)])
    fi
  fi

  # Clean-up
  unset tmp_fft_flavor_chk
  unset tmp_fft_invalid
]) # _SD_FFT_CHECK_CONFIG


AC_DEFUN([_SD_FFT_CHECK_DEFAULTS], [
  # Policy and status must be defined before initialisation
  # Note: this is a developer error, hence we abort unconditionally
  if test "${sd_fft_policy}" != "fail" -a \
          "${sd_fft_policy}" != "skip" -a \
          "${sd_fft_policy}" != "warn"; then
    AC_MSG_ERROR([invalid policy sd_fft_policy='${sd_fft_policy}'
                  Valid policies for broken configurations are:
                      'fail', 'skip', or 'warn'])
  fi
  if test "${sd_fft_status}" != "implicit" -a \
          "${sd_fft_status}" != "optional" -a \
          "${sd_fft_status}" != "required"; then
    AC_MSG_ERROR([invalid policy sd_fft_status='${sd_fft_status}'
                  Valid dependency statuses are:
                      'implicit', 'optional', or 'required'])
  fi

  # Default flavor must be valid
  tmp_fft_invalid="no"
  tmp_fft_flavor_chk=`echo "${sd_fft_flavors}" | grep "${sd_fft_flavor_def}"`
  if test "${tmp_fft_flavor_chk}" = ""; then
    case "${sd_fft_policy}" in
      fail)
        AC_MSG_ERROR([invalid default flavor: sd_fft_flavor_def = '${sd_fft_flavor_def}'])
        ;;
      skip)
        tmp_fft_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default flavor: sd_fft_flavor_def = '${sd_fft_flavor_def}'])
        tmp_fft_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong flavor default
  if test "${tmp_fft_invalid}" = "yes"; then
    sd_fft_flavor_def=`echo "${sd_fft_flavors}" | ${AWK} '{print [$]1}'`
    tmp_fft_invalid="no"
    AC_MSG_NOTICE([setting sd_fft_flavor_def to '${sd_fft_flavor_def}'])
  fi

  # Clean-up
  unset tmp_fft_flavor_chk
  unset tmp_fft_invalid
]) # _SD_FFT_CHECK_DEFAULTS


AC_DEFUN([_SD_FFT_DUMP_CONFIG], [
  AC_MSG_CHECKING([for FFT flavor])
  AC_MSG_RESULT([${sd_fft_flavor}])
  AC_MSG_CHECKING([for FFT C preprocessing flags])
  AC_MSG_RESULT([${sd_fft_cppflags}])
  AC_MSG_CHECKING([for FFT C flags])
  AC_MSG_RESULT([${sd_fft_cflags}])
  AC_MSG_CHECKING([for FFT Fortran flags])
  AC_MSG_RESULT([${sd_fft_fcflags}])
  AC_MSG_CHECKING([for FFT linker flags])
  AC_MSG_RESULT([${sd_fft_ldflags}])
  AC_MSG_CHECKING([for FFT library flags])
  AC_MSG_RESULT([${sd_fft_libs}])
]) # _SD_FFT_DUMP_CONFIG

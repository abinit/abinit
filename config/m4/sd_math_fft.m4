## Copyright (C) 2019 Yann Pouillon <devops@materialsevolution.es>

#
# Multi-flavor Fast Fourier Transform support
#


                    # ------------------------------------ #


AC_DEFUN([SD_FFT_INIT], [
  # Init
  sd_fft_cppflags=""
  sd_fft_cflags=""
  sd_fft_cxxflags=""
  sd_fft_fcflags=""
  sd_fft_ldflags=""
  sd_fft_libs=""
  sd_fft_enable=""
  sd_fft_flavor=""
  sd_fft_init="unknown"
  sd_fft_ok="unknown"
  sd_fft_selected_flavors=""

  # Set adjustable parameters
  _SD_FFT_INIT_FLAVORS

  # FIXME: policy and status are hard-coded for now
  sd_fft_policy="skip"
  sd_fft_status="optional"

  # Set reasonable defaults if not provided
  test -z "${sd_fft_policy}" && sd_fft_policy="fail"
  test -z "${sd_fft_status}" && sd_fft_status="optional"

  # Declare configure option
  AC_ARG_WITH([fft-flavor],
    [AS_HELP_STRING([--with-fft-flavor],
      [FFT flavor to select])],
    [ tmp_fft_flavor_ok=`echo "${sd_fft_selected_flavors}" | grep "${withval}"`
      if test "${tmp_fft_flavor_ok}" = ""; then
        AC_MSG_ERROR([invalid FFT flavor: '${withval}'])
      fi
      sd_fft_flavor="${withval}"
      sd_fft_enable="yes"
      sd_fft_init="kwd"
      sd_fft_selected_flavors="${sd_fft_flavor}"
      unset tmp_fft_flavor_ok],
    [ sd_fft_enable="unknown"; sd_fft_flavor="unknown"; sd_fft_init="def"])

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFT_CHECK_CONFIG
  fi

  # Display configuration
  _SD_FFT_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_fft_cppflags)
  AC_SUBST(sd_fft_cflags)
  AC_SUBST(sd_fft_cxxflags)
  AC_SUBST(sd_fft_fcflags)
  AC_SUBST(sd_fft_ldflags)
  AC_SUBST(sd_fft_libs)
  AC_SUBST(sd_fft_choices)
  AC_SUBST(sd_fft_enable)
  AC_SUBST(sd_fft_flavor)
  AC_SUBST(sd_fft_init)
  AC_SUBST(sd_fft_ok)
  AC_SUBST(with_fft_flavor)
]) # SD_FFT_INIT


                    # ------------------------------------ #


AC_DEFUN([SD_FFT_DETECT], [
  for sd_fft_flavor in ${sd_fft_selected_flavors}; do

    AC_MSG_CHECKING([for the FFT flavor to try])
    AC_MSG_RESULT([${sd_fft_flavor}])
    case "${sd_fft_flavor}" in
      dfti)
        _SD_DFTI_DETECT
        if test "${sd_dfti_ok}" = "yes"; then
          sd_fft_ok="yes"
        fi
        ;;
      fftw3|fftw3-mpi|fftw3-threads)
        SD_FFTW3_DETECT
        if test "${sd_fftw3_ok}" = "yes"; then
          sd_fft_cppflags="${sd_fftw3_cppflags}"
          sd_fft_cflags="${sd_fftw3_cflags}"
          sd_fft_fcflags="${sd_fftw3_fcflags}"
          sd_fft_ldflags="${sd_fftw3_ldflags}"
          sd_fft_libs="${sd_fftw3_libs}"
          sd_fft_ok="yes"
        fi
        ;;
      fftw3-mkl)
        # FIXME: linalg should be supported by Steredeg (abi_ -> sd_)
        sd_fftw3_enable="yes"
        sd_fftw3_init="mkl"
        sd_fftw3_cppflags="${abi_linalg_cppflags}"
        sd_fftw3_cflags="${abi_linalg_cflags}"
        sd_fftw3_cxxflags="${abi_linalg_cxxflags}"
        sd_fftw3_fcflags="${abi_linalg_fcflags}"
        sd_fftw3_ldflags="${abi_linalg_ldflags}"
        sd_fftw3_libs="${abi_linalg_libs}"
        SD_FFTW3_DETECT
        if test "${sd_fftw3_ok}" = "yes"; then
          sd_fft_ok="yes"
        fi
        sd_fftw3_cppflags=""
        sd_fftw3_cflags=""
        sd_fftw3_cxxflags=""
        sd_fftw3_fcflags=""
        sd_fftw3_ldflags=""
        sd_fftw3_libs=""
        ;;
      goedecker)
        AC_MSG_NOTICE([selecting the internal Goedecker FFT implementation])
        sd_fft_ok="yes"
        ;;
      pfft)
        SD_PFFT_DETECT
        if test "${sd_pfft_ok}" = "yes"; then
          sd_fft_cppflags="${sd_pfft_cppflags}"
          sd_fft_cflags="${sd_pfft_cflags}"
          sd_fft_fcflags="${sd_pfft_fcflags}"
          sd_fft_ldflags="${sd_pfft_ldflags}"
          sd_fft_libs="${sd_pfft_libs}"
          sd_fft_ok="yes"
        fi
        ;;
      *)
        AC_MSG_ERROR([unsupported FFT flavor: '${sd_fft_flavor}'])
        ;;
    esac

    test "${sd_fft_ok}" = "yes" && break

  done

  # Check that a working FFT implementation has been found
  if test "${sd_fft_ok}" = "yes"; then
    sd_fft_enable="yes"
    AC_MSG_CHECKING([for the actual FFT flavor to use])
    AC_MSG_RESULT([${sd_fft_flavor}])
  else
    AC_MSG_ERROR([invalid FFT configuration
                  Please adjust configure options to point to a working FFT
                  installation or change the requested FFT flavor through the
                  --with-fft-flavor option.])
  fi

  # FIXME: hard-coded FFTW3 options
  case "${sd_fft_flavor}" in
    fftw3-mpi)
      if test "${sd_mpi_enable}" = "yes"; then
        AC_DEFINE([HAVE_FFTW3_MPI], 1,
          [Define to 1 if you have a MPI-enabled FFTW3 library.])
      fi
      ;;
    fftw3-threads)
      AC_DEFINE([HAVE_FFTW3_THREADS], 1,
        [Define to 1 if you have a threads-enabled FFTW3 library.])
      ;;
  esac
]) # SD_FFT_DETECT


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_FFT_CHECK_CONFIG], [
  # Check consistency between trigger value and package status
  tmp_fft_invalid="no"
  tmp_fft_flavor_chk=`echo "${sd_fft_choices}" | grep "${sd_fft_flavor}"`
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
  tmp_fft_flavor_chk=`echo "${sd_fft_choices}" | grep "${sd_fft_flavor_def}"`
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
    sd_fft_flavor_def=`echo "${sd_fft_choices}" | ${AWK} '{print [$]1}'`
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


# FIXME: compiler vendors should be managed by Steredeg
# FIXME: linear algebra should be managed by Steredeg
AC_DEFUN([_SD_FFT_INIT_FLAVORS], [
  AC_MSG_CHECKING([which FFT flavors to enable])

  # Start from the internal implementation
  sd_fft_selected_flavors="goedecker"

  # Prepend PFFT if available
  if test "${sd_pfft_init}" != "" -a "${sd_pfft_enable}" != "no"; then
    if test "${sd_mpi_ok}" = "yes"; then
      sd_fft_selected_flavors="pfft ${sd_fft_selected_flavors}"
    fi
  fi

  # Prepend FFTW3 if available
  if test "${sd_fftw3_init}" != "" -a "${sd_fftw3_enable}" != "no"; then
    sd_fft_selected_flavors="fftw3 fftw3-threads ${sd_fft_selected_flavors}"
    if test "${sd_mpi_ok}" = "yes"; then
      sd_fft_selected_flavors="fftw3-mpi ${sd_fft_selected_flavors}"
    fi
  fi

  # Prepend FFTW3-in-MKL if MKL is present and FFTW3 is not set
  if test "${abi_linalg_flavor}" = "mkl" -a "${sd_fftw3_enable}" = "no"; then
    sd_fft_selected_flavors="fftw3-mkl ${sd_fft_selected_flavors}"
  fi

  # Prepend DFTI if linear algebra is MKL
  if test "${abi_linalg_flavor}" = "mkl"; then
    sd_fft_selected_flavors="dfti ${sd_fft_selected_flavors}"
  fi

  AC_MSG_RESULT([${sd_fft_selected_flavors}])
]) # _SD_FFT_INIT_FLAVORS


                    # ------------------------------------ #


#
# Private internal macros
#


# FIXME: linear algebra should be managed by Steredeg
AC_DEFUN([_SD_DFTI_DETECT], [
  # Init
  sd_dfti_cppflags="${abi_linalg_}"
  sd_dfti_cflags="${abi_linalg_}"
  sd_dfti_cxxflags="${abi_linalg_}"
  sd_dfti_fcflags="${abi_linalg_}"
  sd_dfti_ldflags="${abi_linalg_}"
  sd_dfti_libs="${abi_linalg_}"
  sd_dfti_enable="yes"
  sd_dfti_init="mkl"
  sd_dfti_ok="unknown"

  # Display configuration
  _SD_DFTI_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_dfti_enable}" = "auto" -o "${sd_dfti_enable}" = "yes"; then
    _SD_DFTI_CHECK_USE

    if test "${sd_dfti_ok}" = "yes"; then
      if test "${sd_dfti_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_dfti_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_dfti_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_dfti_ldflags}"

      AC_DEFINE([HAVE_DFTI], 1,
        [Define to 1 if you have the DFTI library.])
    else
      if test "${sd_dfti_status}" = "optional" -a \
              "${sd_dfti_init}" = "def"; then
        sd_dfti_enable="no"
        sd_dfti_cppflags=""
        sd_dfti_cflags=""
        sd_dfti_fcflags=""
        sd_dfti_ldflags=""
        sd_dfti_libs=""
      else
        AC_MSG_FAILURE([invalid DFTI configuration])
      fi
    fi
  else
    sd_dfti_enable="no"
    sd_dfti_cppflags=""
    sd_dfti_cflags=""
    sd_dfti_fcflags=""
    sd_dfti_ldflags=""
    sd_dfti_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_DFTI_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_dfti_init}" = "esl"; then
    AC_MSG_NOTICE([will look for DFTI in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_dfti_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_dfti_cppflags}"
    CFLAGS="${CFLAGS} ${sd_dfti_cflags}"
    FCFLAGS="${FCFLAGS} ${sd_dfti_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_dfti_ldflags}"
    LIBS="${sd_dfti_libs} ${LIBS}"
  fi

  # Check DFTI C API
  # FIXME: Very complex to have it work properly, would need a replacement
  #        of AC_LINK_IFELSE accepting prologues, because of the included
  #        mkl_dfti.f90 file.
  AC_MSG_CHECKING([whether the DFTI library works])
  sd_dfti_ok="yes"
  AC_MSG_RESULT([${sd_dfti_ok}])

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_DFTI_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_DFTI_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable DFTI])
  AC_MSG_RESULT([${sd_dfti_enable}])
  if test "${sd_dfti_enable}" != "no"; then
    AC_MSG_CHECKING([how DFTI parameters have been set])
    AC_MSG_RESULT([${sd_dfti_init}])
    AC_MSG_CHECKING([for DFTI C preprocessing flags])
    if test "${sd_dfti_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_cppflags}])
    fi
    AC_MSG_CHECKING([for DFTI C flags])
    if test "${sd_dfti_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_cflags}])
    fi
    AC_MSG_CHECKING([for DFTI Fortran flags])
    if test "${sd_dfti_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_fcflags}])
    fi
    AC_MSG_CHECKING([for DFTI linker flags])
    if test "${sd_dfti_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_ldflags}])
    fi
    AC_MSG_CHECKING([for DFTI library flags])
    if test "${sd_dfti_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_dfti_libs}])
    fi
  fi
]) # _SD_DFTI_DUMP_CONFIG

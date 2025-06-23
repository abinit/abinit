## Copyright (C) 2019-2025 ABINIT group (Yann Pouillon)

#
# Fastest Fourier Transform in the West library (FFTW3)
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_FFTW3_INIT], [
  # Init
  sd_fftw3_cppflags=""
  sd_fftw3_cflags=""
  sd_fftw3_cxxflags=""
  sd_fftw3_fcflags=""
  sd_fftw3_ldflags=""
  sd_fftw3_libs=""
  sd_fftw3_enable=""
  sd_fftw3_init="unknown"
  sd_fftw3_mpi_ok="unknown"
  sd_fftw3_ok="unknown"
  sd_fftw3_threads_ok="unknown"

  # Set adjustable parameters
  sd_fftw3_options="$1"
  sd_fftw3_libs_def="$2"
  sd_fftw3_cppflags_def="$3"
  sd_fftw3_cflags_def="$4"
  sd_fftw3_cxxflags_def="$5"
  sd_fftw3_fcflags_def="$6"
  sd_fftw3_ldflags_def="$7"
  sd_fftw3_enable_def=""
  sd_fftw3_policy=""
  sd_fftw3_status=""

  # Process options
  for kwd in ${sd_fftw3_options}; do
    case "${kwd}" in
      auto)
        sd_fftw3_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_fftw3_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_fftw3_policy="${kwd}"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg FFTW3 option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  if test -z "${sd_fftw3_libs_def}"; then
    sd_fftw3_libs_def="-lfftw3"
    if test "${abi_openmp_enable}" = "yes"; then 
      sd_fftw3_libs_def="-lfftw3_threads ${sd_fftw3_libs_def}"
    fi
    if test "${sd_mpi_enable}" = "yes"; then 
      sd_fftw3_libs_def="-lfftw3_mpi ${sd_fftw3_libs_def}"
    fi
  fi
  test -z "${sd_fftw3_policy}" && sd_fftw3_policy="fail"
  test -z "${sd_fftw3_status}" && sd_fftw3_status="optional"
  test -z "${sd_fftw3_enable_def}" && sd_fftw3_enable_def="no"
  case "${sd_fftw3_status}" in 
         implicit|required|optional)
      sd_fftw3_enable_def="yes"
      ;;
  esac

  AC_MSG_NOTICE([status for fftw3 TMP '${sd_fftw3_status}'])
  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([fftw3],
    [AS_HELP_STRING([--with-fftw3],
      [Install prefix of the FFTW3 library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_fftw3_enable="${withval}"
        sd_fftw3_init="yon"
      else
        sd_fftw3_enable="yes"
        sd_fftw3_init="dir"
      fi],
    [ sd_fftw3_enable="${sd_fftw3_enable_def}"; sd_fftw3_init="def"])

  # Declare environment variables
  AC_ARG_VAR([FFTW3_CPPFLAGS], [C preprocessing flags for FFTW3.])
  AC_ARG_VAR([FFTW3_CFLAGS], [C flags for FFTW3.])
  AC_ARG_VAR([FFTW3_FCFLAGS], [Fortran flags for FFTW3.])
  AC_ARG_VAR([FFTW3_FFLAGS], [Fortran flags for FFTW3 (better use FFTW3_FCFLAGS).])
  AC_ARG_VAR([FFTW3_LDFLAGS], [Linker flags for FFTW3.])
  AC_ARG_VAR([FFTW3_LIBS], [Library flags for FFTW3.])

  # Detect use of environment variables
  if test "${sd_fftw3_enable}" = "yes" -o "${sd_fftw3_enable}" = "auto"; then
    tmp_fftw3_vars="${FFTW3_CPPFLAGS}${FFTW3_CFLAGS}${FFTW3_FFLAGS}${FFTW3_FCFLAGS}${FFTW3_LDFLAGS}${FFTW3_LIBS}"
    if test "${sd_fftw3_init}" = "def" -a ! -z "${tmp_fftw3_vars}"; then
      sd_fftw3_enable="yes"
      sd_fftw3_init="env"
    fi
  fi


  # if mode is def and pkg_config exists and no prefix -> use pkg_config
  #
  if test "${sd_fftw3_enable}" = "yes" -a \( "${sd_fftw3_init}" = "def" -o "${sd_fftw3_init}" = "yon" \) -a "${sd_fftw3_prefix}" = ""; then
     #check if PKG_CONFIG exists (if not keep default mode)
     AC_MSG_NOTICE([setting for ${sd_fftw3_init} potential move to pkg])
     AC_CHECK_PROG([PKG_CONFIG], [pkg-config], [pkg-config], [no])
     
     AC_MSG_NOTICE([setting for ${sd_fftw3_init} potential move to pkg, PKG=${PKG_CONFIG}])
     if test "$PKG_CONFIG" != "no"; then
	 AC_MSG_CHECKING([for fftw3 via pkg-config])
         AC_PATH_TOOL(PKG_CONFIG,pkg-config)
         if "$PKG_CONFIG" --exists  fftw3; then
		AC_MSG_RESULT([yes])
                if test "${sd_fftw3_init}" = "yon";then
		   sd_fftw3_status="required"
		fi
		sd_fftw3_init="pkg"
	 else
		AC_MSG_RESULT([no])
		# sd_fftw3_init="def" or sd_fftw3_init="yon"
	 fi
     fi
  fi

  AC_MSG_NOTICE([setup to ${sd_fftw3_init} after checking for pkg])
  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_FFTW3_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_fftw3_enable}" = "yes" -o "${sd_fftw3_enable}" = "auto"; then

    # Set FFTW3-specific flags
    case "${sd_fftw3_init}" in

      def|yon)
        sd_fftw3_cppflags="${sd_fftw3_cppflags_def}"
        sd_fftw3_cflags="${sd_fftw3_cflags_def}"
        sd_fftw3_fcflags="${sd_fftw3_fcflags_def}"
        sd_fftw3_ldflags="${sd_fftw3_ldflags_def}"
        sd_fftw3_libs="${sd_fftw3_libs_def}"
        ;;

      dir)
        sd_fftw3_cppflags="-I${with_fftw3}/include"
        sd_fftw3_cflags="${sd_fftw3_cflags_def}"
        sd_fftw3_fcflags="${sd_fftw3_fcflags_def}"
        sd_fftw3_ldflags="${sd_fftw3_ldflags_def}"
        sd_fftw3_libs="-L${with_fftw3}/lib ${sd_fftw3_libs_def}"
        ;;

      env)
        sd_fftw3_cppflags="${sd_fftw3_cppflags_def}"
        sd_fftw3_cflags="${sd_fftw3_cflags_def}"
        sd_fftw3_fcflags="${sd_fftw3_fcflags_def}"
        sd_fftw3_ldflags="${sd_fftw3_ldflags_def}"
        sd_fftw3_libs="${sd_fftw3_libs_def}"
        test ! -z "${FFTW3_CPPFLAGS}" && sd_fftw3_cppflags="${FFTW3_CPPFLAGS}"
        test ! -z "${FFTW3_CFLAGS}" && sd_fftw3_cflags="${FFTW3_CFLAGS}"
        test ! -z "${FFTW3_FFLAGS}" && sd_fftw3_fcflags="${FFTW3_FFLAGS}"
        test ! -z "${FFTW3_FCFLAGS}" && sd_fftw3_fcflags="${FFTW3_FCFLAGS}"
        test ! -z "${FFTW3_LDFLAGS}" && sd_fftw3_ldflags="${FFTW3_LDFLAGS}"
        test ! -z "${FFTW3_LIBS}" && sd_fftw3_libs="${FFTW3_LIBS}"
        ;;

      pkg)
        TMP_FFTW_CPPFLAGS=`$PKG_CONFIG --cflags --keep-system-cflags fftw3`
        TMP_FFTW_FFFLAGS="${TMP_FFTW_CPPFLAGS}"
        TMP_FFTW_LIBS=`$PKG_CONFIG --libs  --keep-system-libs fftw3`
        TMP_FFTWF_CPPFLAGS=`$PKG_CONFIG --cflags --keep-system-cflags fftw3f`
        TMP_FFTWF_FFFLAGS="${TMP_FFTW_CPPFLAGS}"
        TMP_FFTWF_LIBS=`$PKG_CONFIG --libs  --keep-system-libs fftw3f`
        sd_fftw3_cppflags="${TMP_FFTW_CPPFLAGS} ${TMP_FFTWF_CPPFLAGS}"
        sd_fftw3_cflags="${TMP_FFTW_CPPFLAGS} ${TMP_FFTWF_CPPFLAGS}"
        #sd_fftw3_cxxflags="${TMP_FFTW_CPPFLAGS}"
        sd_fftw3_fcflags="${TMP_FFTW_FFLAGS} ${TMP_FFTWF_FFLAGS}"
        sd_fftw3_ldflags="${TMP_FFTW_LIBS} ${TMP_FFTWF_LIBS}"
        sd_fftw3_libs="${TMP_FFTW_LIBS} ${TMP_FFTWF_LIBS} "
        if test "${abi_openmp_enable}" = "yes"; then
           sd_fftw3_libs="-lfftw3_threads -lpthread -lfftw3f_threads ${sd_fftw3_libs}"
           sd_fftw3_ldflags="-lfftw3_threads -lpthread -lfftw3f_threads ${sd_fftw3_ldflags}"
        fi
        #if test "${sd_mpi_enable}" = "yes"; then
        #   sd_fftw3_libs="-lfftw3_mpi -lfftw3f_mpi -lm ${sd_fftw3_libs}"
        #   sd_fftw3_ldflags="-lfftw3_mpi -lm -lfftw3f_mpi ${sd_fftw3_ldflags}"
        #fi
        ;;

      *)
        AC_MSG_ERROR([invalid init type for FFTW3: '${sd_fftw3_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_fftw3_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_fftw3_init="esl"
    sd_fftw3_cppflags=""
    sd_fftw3_cflags=""
    sd_fftw3_fcflags=""
    sd_fftw3_ldflags=""
    sd_fftw3_libs=""
  fi

  # Display configuration
  _SD_FFTW3_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_fftw3_options)
  AC_SUBST(sd_fftw3_enable_def)
  AC_SUBST(sd_fftw3_policy)
  AC_SUBST(sd_fftw3_status)
  AC_SUBST(sd_fftw3_enable)
  AC_SUBST(sd_fftw3_init)
  AC_SUBST(sd_fftw3_ok)
  AC_SUBST(sd_fftw3_cppflags)
  AC_SUBST(sd_fftw3_cflags)
  AC_SUBST(sd_fftw3_fcflags)
  AC_SUBST(sd_fftw3_ldflags)
  AC_SUBST(sd_fftw3_libs)
  AC_SUBST(with_fftw3)

  # Clean-up
  unset tmp_fftw3_vars
]) # SD_FFTW3_INIT


AC_DEFUN([SD_FFTW3_DETECT], [
  # Display configuration
  _SD_FFTW3_DUMP_CONFIG

  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_fftw3_enable}" = "auto" -o "${sd_fftw3_enable}" = "yes"; then
    _SD_FFTW3_CHECK_USE

    if test "${sd_fftw3_ok}" = "yes"; then
      if test "${sd_fftw3_init}" = "esl"; then
        sd_esl_bundle_libs="${sd_fftw3_libs_def} ${sd_esl_bundle_libs}"
      else
        LIBS="${sd_fftw3_libs} ${LIBS}"
      fi
      LDFLAGS="${LDFLAGS} ${sd_fftw3_ldflags}"

      AC_DEFINE([HAVE_FFTW3], 1,
        [Define to 1 if you have the FFTW3 library.])
    else
      if test "${sd_fftw3_status}" = "optional" -a \
              "${sd_fftw3_init}" = "def"; then
        sd_fftw3_enable="no"
        sd_fftw3_cppflags=""
        sd_fftw3_cflags=""
        sd_fftw3_fcflags=""
        sd_fftw3_ldflags=""
        sd_fftw3_libs=""
      elif test "${sd_fftw3_status}" = "optional" -a \
              "${sd_fftw3_init}" = "pkg"; then
        sd_fftw3_enable="no"
        sd_fftw3_cppflags=""
        sd_fftw3_cflags=""
        sd_fftw3_fcflags=""
        sd_fftw3_ldflags=""
      else
        AC_MSG_FAILURE([invalid FFTW3 configuration])
      fi
    fi
  else
    sd_fftw3_enable="no"
    sd_fftw3_cppflags=""
    sd_fftw3_cflags=""
    sd_fftw3_fcflags=""
    sd_fftw3_ldflags=""
    sd_fftw3_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_FFTW3_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  if test "${sd_fftw3_init}" = "esl"; then
    AC_MSG_NOTICE([will look for FFTW3 in the installed ESL Bundle])
    SD_ESL_ADD_FLAGS
    SD_ESL_ADD_LIBS([${sd_fftw3_libs_def}])
  else
    CPPFLAGS="${CPPFLAGS} ${sd_fftw3_cppflags}"
    CFLAGS="${CFLAGS} ${sd_fftw3_cflags}"
    FCFLAGS="${FCFLAGS} ${sd_fftw3_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_fftw3_ldflags}"
    LIBS="${sd_fftw3_libs} ${LIBS}"
  fi

  # Check FFTW3 API
  AC_MSG_CHECKING([whether the FFTW3 library works])

  # ===== Fortran support
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        integer(8) :: plan
        complex(8) :: ff(10)
        call dfftw_execute_dft(plan,ff,ff)
    ]])], [sd_fftw3_ok="yes"], [sd_fftw3_ok="no"])
  AC_LANG_POP([Fortran])

  # ===== C support
  # AC_LANG_PUSH([C])
  # AC_LINK_IFELSE([AC_LANG_PROGRAM(
  #   [[
  #   #include <fftw3.h>
  #   ]],
  #   [[
  #     fftw_plan *plan;
  #     fftw_complex *a1, *a2;
  #     fftw_execute_dft(*plan, a1, a2);
  #   ]])], [sd_fftw3_ok="yes"], [sd_fftw3_ok="no"])
  # AC_LANG_POP([C])

  AC_MSG_RESULT([${sd_fftw3_ok}])

  # Check for FFTW3 threads support
  if test "${sd_fftw3_ok}" = "yes"; then
    AC_MSG_CHECKING([whether the FFTW3 library supports threads])

    # ===== Fortran support
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        integer iret
        call dfftw_init_threads(iret)
        call dfftw_plan_with_nthreads(8)
        call dfftw_cleanup_threads()
      ]])], [sd_fftw3_threads_ok="yes"], [sd_fftw3_threads_ok="no"])
    AC_LANG_POP([Fortran])

    # ===== C support
    # AC_LANG_PUSH([C])
    # AC_LINK_IFELSE([AC_LANG_PROGRAM(
    #    [[
    #    #include <fftw3.h>
    #    ]],
    #    [[
    #    fftw_init_threads;
    #    ]])], [sd_fftw3_ok="yes"], [sd_fftw3_ok="no"])
    # AC_LANG_POP([C])

    AC_MSG_RESULT([${sd_fftw3_threads_ok}])
  fi

  # Check for FFTW3 MPI support
  if test "${sd_fftw3_ok}" = "yes" -a "${sd_mpi_enable}" = "yes"; then
    AC_MSG_CHECKING([whether the FFTW3 MPI library works])

    # ===== Fortran support
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use,intrinsic :: iso_c_binding
        include 'fftw3-mpi.f03'
        call fftw_mpi_init()
      ]])], [sd_fftw3_mpi_ok="yes"], [sd_fftw3_mpi_ok="no"])
    AC_LANG_POP([Fortran])

    # ===== C support
    # AC_LANG_PUSH([C])
    # AC_LINK_IFELSE([AC_LANG_PROGRAM(
    #   [[
    #   #include <fftw3-mpi.h>
    #   ]],
    #   [[
    #     fftw_mpi_init();
    #   ]])], [sd_fftw3_mpi_ok="yes"], [sd_fftw3_mpi_ok="no"])
    # AC_LANG_POP([C])

    AC_MSG_RESULT([${sd_fftw3_mpi_ok}])
  else
    sd_fftw3_mpi_ok="no"
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_FFTW3_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_FFTW3_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_fftw3_invalid="no"
  if test "${sd_fftw3_enable_def}" != "auto" -a \
          "${sd_fftw3_enable_def}" != "no" -a \
          "${sd_fftw3_enable_def}" != "yes"; then
    case "${sd_fftw3_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_fftw3_enable_def = '${sd_fftw3_enable_def}'])
        ;;
      skip)
        tmp_fftw3_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_fftw3_enable_def = '${sd_fftw3_enable_def}'])
        tmp_fftw3_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_fftw3_invalid}" = "yes"; then
    if test "${sd_fftw3_status}" = "required"; then
      sd_fftw3_enable_def="yes"
    else
      sd_fftw3_enable_def="no"
    fi
    tmp_fftw3_invalid="no"
    AC_MSG_NOTICE([setting sd_fftw3_enable_def to '${sd_fftw3_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_fftw3_invalid="no"
  if test "${sd_fftw3_status}" = "implicit" -o \
          "${sd_fftw3_status}" = "required"; then
    if test "${sd_fftw3_enable}" = "no"; then
      case "${sd_fftw3_policy}" in
        fail)
          AC_MSG_ERROR([The FFTW3 package is required and cannot be disabled
                  See http://www.fftw.org/ for details on how to install it.])
          ;;
        skip)
          tmp_fftw3_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The FFTW3 package is required and cannot be disabled])
          tmp_fftw3_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_fftw3_enable}" = "auto"; then
      AC_MSG_NOTICE([setting FFTW3 trigger to yes])
      sd_fftw3_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_fftw3_invalid}" = "yes"; then
    case "${sd_fftw3_status}" in
      implicit|required)
        sd_fftw3_enable="yes"
        ;;
      optional)
        sd_fftw3_enable="no"
        ;;
    esac
    tmp_fftw3_invalid="no"
    AC_MSG_NOTICE([setting sd_fftw3_enable to '${sd_fftw3_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_fftw3_vars="${FFTW3_CPPFLAGS}${FFTW3_CFLAGS}${FFTW3_LDFLAGS}${FFTW3_LIBS}"
  tmp_fftw3_invalid="no"
  if test ! -z "${tmp_fftw3_vars}" -a ! -z "${with_fftw3}"; then
    case "${sd_fftw3_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for FFTW3
                  Please use FFTW3_CFLAGS + FFTW3_LIBS or --with-fftw3,
                  not both.])
        ;;
      skip)
        tmp_fftw3_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for FFTW3])
        tmp_fftw3_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_fftw3_vars}" -a ! "${sd_fftw3_init}" = "pkg"; then
    sd_fftw3_enable="yes"
    sd_fftw3_init="env"
    if test "${tmp_fftw3_invalid}" = "yes"; then
      tmp_fftw3_invalid="no"
      AC_MSG_NOTICE([overriding --with-fftw3 with FFTW3_{CPPFLAGS,CFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_fftw3_status}" = "implicit"; then
    if test "${sd_fftw3_ldflags}" != ""; then
      sd_fftw3_ldflags=""
      AC_MSG_NOTICE([resetting FFTW3 linker flags (implicit package)])
    fi
    if test "${sd_fftw3_libs}" != ""; then
      sd_fftw3_libs=""
      AC_MSG_NOTICE([resetting FFTW3 library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_fftw3_enable}" = "implicit"; then
    sd_fftw3_cppflags=""
    sd_fftw3_cflags=""
    sd_fftw3_fcflags=""
    sd_fftw3_ldflags=""
    sd_fftw3_libs=""
  fi

  # Clean-up
  unset tmp_fftw3_invalid
  unset tmp_fftw3_vars
]) # _SD_FFTW3_CHECK_CONFIG


AC_DEFUN([_SD_FFTW3_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable FFTW3])
  AC_MSG_RESULT([${sd_fftw3_enable}])
  if test "${sd_fftw3_enable}" != "no"; then
    AC_MSG_CHECKING([how FFTW3 parameters have been set])
    AC_MSG_RESULT([${sd_fftw3_init}])
    AC_MSG_CHECKING([for FFTW3 C preprocessing flags])
    if test "${sd_fftw3_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_fftw3_cppflags}])
    fi
    AC_MSG_CHECKING([for FFTW3 C flags])
    if test "${sd_fftw3_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_fftw3_cflags}])
    fi
    AC_MSG_CHECKING([for FFTW3 Fortran flags])
    if test "${sd_fftw3_fcflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_fftw3_fcflags}])
    fi
    AC_MSG_CHECKING([for FFTW3 linker flags])
    if test "${sd_fftw3_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_fftw3_ldflags}])
    fi
    AC_MSG_CHECKING([for FFTW3 library flags])
    if test "${sd_fftw3_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_fftw3_libs}])
    fi
  fi
]) # _SD_FFTW3_DUMP_CONFIG

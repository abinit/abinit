## Copyright (C) 2019 Yann Pouillon

#
# HDF5 I/O support
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_HDF5_INIT], [
  # Init
  sd_hdf5_h5cc=""
  sd_hdf5_h5fc=""
  sd_hdf5_cppflags=""
  sd_hdf5_cflags=""
  sd_hdf5_cxxflags=""
  sd_hdf5_fcflags=""
  sd_hdf5_ldflags=""
  sd_hdf5_libs=""
  sd_hdf5_libs_extra=""
  sd_hdf5_libs_hl=""
  sd_hdf5_prefix=""
  sd_hdf5_enable=""
  sd_hdf5_flavor=""
  sd_hdf5_init="unknown"
  sd_hdf5_has_fc="unknown"
  sd_hdf5_has_mpi="unknown"
  sd_hdf5_c_ok="unknown"
  sd_hdf5_fortran_ok="unknown"
  sd_hdf5_mpi_ok="unknown"
  sd_hdf5_ok="unknown"

  # Set adjustable parameters
  sd_hdf5_options="$1"
  sd_hdf5_libs_def="$2"
  sd_hdf5_cppflags_def="$3"
  sd_hdf5_cflags_def="$4"
  sd_hdf5_cxxflags_def="$5"
  sd_hdf5_fcflags_def="$6"
  sd_hdf5_ldflags_def="$7"
  sd_hdf5_enable_def=""
  sd_hdf5_enable_fc=""
  sd_hdf5_policy=""
  sd_hdf5_status=""
  sd_hdf5_cascade="no"

  # Process options
  for kwd in ${sd_hdf5_options}; do
    case "${kwd}" in
      auto)
        sd_hdf5_enable_def="${kwd}"
        ;;
      implicit|required|optional)
        sd_hdf5_status="${kwd}"
        ;;
      fail|skip|warn)
        sd_hdf5_policy="${kwd}"
        ;;
      no-fortran)
        sd_hdf5_enable_fc="no"
        ;;
      mandatory)
        sd_hdf5_enable="yes"
        sd_hdf5_enable_def="yes"
        ;;
      *)
        AC_MSG_ERROR([invalid Steredeg HDF5 option: '${kwd}'])
        ;;
    esac
  done

  # Set reasonable defaults if not provided
  test -z "${sd_hdf5_enable_fc}" && sd_hdf5_enable_fc="yes"
  if test "${sd_hdf5_enable_fc}" = "yes"; then
    test -z "${sd_hdf5_libs_def}" && sd_hdf5_libs_def="-lhdf5_fortran -lhdf5"
  else
    test -z "${sd_hdf5_libs_def}" && sd_hdf5_libs_def="-lhdf5"
  fi
  test -z "${sd_hdf5_policy}" && sd_hdf5_policy="fail"
  test -z "${sd_hdf5_status}" && sd_hdf5_status="optional"
  test -z "${sd_hdf5_enable_def}" && sd_hdf5_enable_def="no"
  case "${sd_hdf5_status}" in
    implicit|required)
      sd_hdf5_enable_def="yes"
      ;;
  esac

  # Declare configure option
  # TODO: make it switchable for the implicit case 
  AC_ARG_WITH([hdf5],
    [AS_HELP_STRING([--with-hdf5],
      [Install prefix of the HDF5 library (e.g. /usr/local).])],
    [ if test "${withval}" = "no" -o "${withval}" = "yes"; then
        sd_hdf5_enable="${withval}"
        sd_hdf5_init="yon"
      else
        sd_hdf5_prefix="${withval}"
        sd_hdf5_enable="yes"
        sd_hdf5_init="dir"
      fi],
    [ sd_hdf5_enable="${sd_hdf5_enable_def}"; sd_hdf5_init="def"])

  # Declare environment variables
  AC_ARG_VAR([H5CC], [HDF5-aware C compiler.])
  AC_ARG_VAR([HDF5_CPPFLAGS], [C preprocessing flags for HDF5.])
  AC_ARG_VAR([HDF5_CFLAGS], [C flags for HDF5.])
  AC_ARG_VAR([HDF5_CXXFLAGS], [C flags for HDF5.])
  AC_ARG_VAR([HDF5_FCFLAGS], [Fortran flags for HDF5.])
  AC_ARG_VAR([HDF5_LDFLAGS], [Linker flags for HDF5.])
  AC_ARG_VAR([HDF5_LIBS], [Library flags for HDF5.])

  # Detect use of HDF5-aware compiler
  if test "${sd_hdf5_init}" = "def" -a ! -z "${H5CC}"; then
    sd_hdf5_enable="yes"
    sd_hdf5_init="exe"
  fi

  # Detect use of environment variables
  if test "${sd_hdf5_enable}" != "no" -o "${sd_hdf5_init}" = "def"; then
    tmp_hdf5_vars="${HDF5_CPPFLAGS}${HDF5_CFLAGS}${HDF5_LDFLAGS}${HDF5_LIBS}"
    if test "${sd_hdf5_enable_fc}" = "yes"; then
      tmp_hdf5_vars="${tmp_hdf5_vars}${HDF5_FCFLAGS}"
    fi
    if test ! -z "${tmp_hdf5_vars}"; then
      case "${sd_hdf5_init}" in
        def)
          sd_hdf5_enable="yes"
          sd_hdf5_init="env"
          ;;
        exe)
          sd_hdf5_cascade="yes"
          ;;
      esac
    fi
  fi

  # Make sure configuration is correct
  if test "${STEREDEG_BYPASS_CONSISTENCY}" != "yes"; then
    _SD_HDF5_CHECK_CONFIG
  fi

  # Adjust configuration depending on init type
  if test "${sd_hdf5_enable}" = "yes" -o "${sd_hdf5_enable}" = "auto"; then

    # Set HDF5-specific flags
    case "${sd_hdf5_init}" in

      def|exe|yon)
        sd_hdf5_cppflags="${sd_hdf5_cppflags_def}"
        sd_hdf5_cflags="${sd_hdf5_cflags_def}"
        sd_hdf5_cxxflags="${sd_hdf5_cxxflags_def}"
        test "${sd_hdf5_enable_fc}" = "yes" && \
          sd_hdf5_fcflags="${sd_hdf5_fcflags_def}"
        sd_hdf5_ldflags="${sd_hdf5_ldflags_def}"
        sd_hdf5_libs="${sd_hdf5_libs_def}"
        ;;

      dir)
        sd_hdf5_cppflags="-I${sd_hdf5_prefix}/include"
        sd_hdf5_cflags="${sd_hdf5_cflags_def}"
        sd_hdf5_cxxflags="${sd_hdf5_cxxflags_def}"
        test "${sd_hdf5_enable_fc}" = "yes" && \
          sd_hdf5_fcflags="${sd_hdf5_fcflags_def} -I${sd_hdf5_prefix}/include"
        sd_hdf5_ldflags="${sd_hdf5_ldflags_def}"
        sd_hdf5_libs="-L${sd_hdf5_prefix}/lib ${sd_hdf5_libs_def}"
        ;;

      env)
        _SD_HDF5_SET_ENV
        ;;

      *)
        AC_MSG_ERROR([invalid init type for HDF5: '${sd_hdf5_init}'])
        ;;

    esac

  fi

  # Override the default configuration if the ESL Bundle is available
  # Note: The setting of the build flags is done once for all
  #       ESL-bundled packages
  if test "${sd_hdf5_init}" = "def" -a ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_hdf5_init="esl"
    sd_hdf5_cppflags=""
    sd_hdf5_cflags=""
    sd_hdf5_cxxflags=""
    sd_hdf5_fcflags=""
    sd_hdf5_ldflags=""
    sd_hdf5_libs=""
  fi

  # Display configuration
  _SD_HDF5_DUMP_CONFIG

  # Export configuration
  AC_SUBST(sd_hdf5_options)
  AC_SUBST(sd_hdf5_enable_def)
  AC_SUBST(sd_hdf5_enable_fc)
  AC_SUBST(sd_hdf5_policy)
  AC_SUBST(sd_hdf5_status)
  AC_SUBST(sd_hdf5_enable)
  AC_SUBST(sd_hdf5_init)
  AC_SUBST(sd_hdf5_ok)
  AC_SUBST(sd_hdf5_h5cc)
  AC_SUBST(sd_hdf5_h5fc)
  AC_SUBST(sd_hdf5_cppflags)
  AC_SUBST(sd_hdf5_cflags)
  AC_SUBST(sd_hdf5_cxxflags)
  AC_SUBST(sd_hdf5_fcflags)
  AC_SUBST(sd_hdf5_ldflags)
  AC_SUBST(sd_hdf5_libs)
  AC_SUBST(with_hdf5)

  # Clean-up
  unset tmp_hdf5_vars
]) # SD_HDF5_INIT


AC_DEFUN([SD_HDF5_DETECT], [
  # Check whether we can compile and link a simple program
  # and update build flags if successful
  if test "${sd_hdf5_enable}" = "auto" -o "${sd_hdf5_enable}" = "yes"; then
    if test "${sd_hdf5_init}" != "env"; then
      _SD_HDF5_CHECK_COMPILERS
    fi
    if test "${sd_hdf5_h5cc}" = ""; then
      if test "${sd_hdf5_cascade}" = "yes"; then
        _SD_HDF5_SET_ENV
      else
        AC_MSG_WARN([could not find the HDF5 C compiler
                    Please set the H5CC environment variable to a valid HDF5
                    C compiler (usually called 'h5cc' or 'h5pcc'), so that
                    essential information about your HDF5 installation will be
                    exposed to the build system of ABINIT.])
      fi
    fi

    _SD_HDF5_DUMP_CONFIG
    _SD_HDF5_CHECK_USE

    if test "${sd_hdf5_ok}" = "yes"; then
      CPPFLAGS="${CPPFLAGS} ${sd_hdf5_cppflags}"
      CFLAGS="${CFLAGS} ${sd_hdf5_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_hdf5_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_hdf5_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_hdf5_ldflags}"
      LIBS="${sd_hdf5_libs} ${LIBS}"

      AC_DEFINE([HAVE_HDF5], 1,
        [Define to 1 if you have the HDF5 library.])
      if test "${sd_hdf5_mpi_ok}" = "yes"; then
        AC_DEFINE([HAVE_HDF5_MPI], 1,
          [Define to 1 if you have a parallel HDF5 library.])
      fi
    else
      if test "${sd_hdf5_status}" = "optional" -a \
              "${sd_hdf5_init}" = "def"; then
        sd_hdf5_enable="no"
        sd_hdf5_cppflags=""
        sd_hdf5_cflags=""
        sd_hdf5_cxxflags=""
        sd_hdf5_fcflags=""
        sd_hdf5_ldflags=""
        sd_hdf5_libs=""
      else
        if test "${sd_hdf5_policy}" = "fail"; then
              AC_MSG_FAILURE([invalid HDF5 configuration])
        else
              AC_MSG_WARN([invalid HDF5 configuration])
        fi
      fi
    fi
  else
    sd_hdf5_enable="no"
    sd_hdf5_cppflags=""
    sd_hdf5_cflags=""
    sd_hdf5_cxxflags=""
    sd_hdf5_fcflags=""
    sd_hdf5_ldflags=""
    sd_hdf5_libs=""
  fi
])


                    # ------------------------------------ #


#
# Private macros
#


AC_DEFUN([_SD_HDF5_CHECK_USE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_hdf5_cppflags}"
  CFLAGS="${CFLAGS} ${sd_hdf5_cflags}"
  CXXFLAGS="${CFLAGS} ${sd_hdf5_cxxflags}"
  FCFLAGS="${FCFLAGS} ${sd_hdf5_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_hdf5_ldflags}"
  LIBS="${sd_hdf5_libs} ${LIBS}"

  # Look for hdf5.h
  AC_CHECK_HEADERS([hdf5.h], [sd_hdf5_hdr_ok="yes"], [sd_hdf5_hdr_ok="no"])
  if test "${sd_hdf5_hdr_ok}" != "yes"; then
    AC_MSG_WARN([could not find the HDF5 C header
                  Detecting HDF5 is a complex task due to the very high number
                  of possible configurations. If you do want to use HDF5,
                  please help the build system by pointing correctly the
                  HDF5_CPPFLAGS, HDF5_FCFLAGS and HDF5_LIBS parameters to
                  where the desired hdf5.h and libhdf5.* files are installed.])
  fi

  # Check HDF5 C API
  AC_MSG_CHECKING([whether the HDF5 C interface works])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <hdf5.h>
    ]],
    [[
      hid_t file_id;
      file_id = H5Fcreate("conftest.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    ]])], [sd_hdf5_c_ok="yes"], [sd_hdf5_c_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_hdf5_c_ok}])

  # Check HDF5 Fortran API
  if test "${sd_hdf5_c_ok}" = "yes" -a "${sd_hdf5_enable_fc}" = "yes"; then
    AC_MSG_CHECKING([whether the HDF5 Fortran interface works])
    for tmp_incs in "" "-I/usr/include"; do
      FCFLAGS="${FCFLAGS} ${tmp_incs}"
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          use hdf5
          integer(hid_t) :: file_id
          integer :: ierr
          call h5open_f(ierr)
          call h5fcreate_f("conftest.h5f", H5F_ACC_TRUNC_F, file_id, ierr)
        ]])], [sd_hdf5_fortran_ok="yes"], [sd_hdf5_fortran_ok="no"])
      AC_LANG_POP([Fortran])
      if test "${sd_hdf5_fortran_ok}" = "yes"; then
        test "${sd_fc_sys_incpath}" = "" && sd_fc_sys_incpath="${tmp_incs}"
        break
      fi
    done
    AC_MSG_RESULT([${sd_hdf5_fortran_ok}])
  fi
  unset tmp_incs

  # Check HDF5 parallel API
  if test "${sd_hdf5_c_ok}" = "yes" -a "${sd_mpi_ok}" = "yes" -a \
          "${sd_hdf5_has_mpi}" = "yes"; then
    AC_MSG_CHECKING([whether the parallel HDF5 API works])
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
      [[
#include <mpi.h>
#include <hdf5.h>
      ]],
      [[
        hid_t acc_tpl, file_id;
        herr_t ret;
        MPI_Comm comm = MPI_COMM_WORLD;
        MPI_Info info = MPI_INFO_NULL;
        acc_tpl = H5Pcreate(H5P_FILE_ACCESS);
        ret = H5Pset_fapl_mpio(acc_tpl, comm, info);
        file_id = H5Fcreate("conftest.h5", H5F_ACC_TRUNC, H5P_DEFAULT, acc_tpl);
      ]])], [sd_hdf5_mpi_ok="yes"], [sd_hdf5_mpi_ok="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${sd_hdf5_mpi_ok}])
  else
    sd_hdf5_has_mpi="no"
    sd_hdf5_mpi_ok="no"
  fi

  # Combine the available results
  sd_hdf5_ok="no"
  if test "${sd_hdf5_enable_fc}" = "yes"; then
    if test "${sd_hdf5_c_ok}" = "yes" -a \
            "${sd_hdf5_fortran_ok}" = "yes"; then
      sd_hdf5_ok="yes"
    fi
  else
    if test "${sd_hdf5_c_ok}" = "yes"; then
      sd_hdf5_ok="yes"
    fi
  fi

  # Restore environment
  SD_ESL_RESTORE_FLAGS
]) # _SD_HDF5_CHECK_USE


                    # ------------------------------------ #
                    # ------------------------------------ #


#
# Utility macros
#


AC_DEFUN([_SD_HDF5_CHECK_CONFIG], [
  # Default trigger must be yes, no, or auto
  tmp_hdf5_invalid="no"
  if test "${sd_hdf5_enable_def}" != "auto" -a \
          "${sd_hdf5_enable_def}" != "no" -a \
          "${sd_hdf5_enable_def}" != "yes"; then
    case "${sd_hdf5_policy}" in
      fail)
        AC_MSG_ERROR([invalid default value: sd_hdf5_enable_def = '${sd_hdf5_enable_def}'])
        ;;
      skip)
        tmp_hdf5_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([invalid default value: sd_hdf5_enable_def = '${sd_hdf5_enable_def}'])
        tmp_hdf5_invalid="yes"
        ;;
    esac
  fi

  # Fix wrong trigger default
  if test "${tmp_hdf5_invalid}" = "yes"; then
    if test "${sd_hdf5_status}" = "required"; then
      sd_hdf5_enable_def="yes"
    else
      sd_hdf5_enable_def="no"
    fi
    tmp_hdf5_invalid="no"
    AC_MSG_NOTICE([setting sd_hdf5_enable_def to '${sd_hdf5_enable_def}'])
  fi

  # Check consistency between trigger value and package status
  tmp_hdf5_invalid="no"
  if test "${sd_hdf5_status}" = "implicit" -o \
          "${sd_hdf5_status}" = "required"; then
    if test "${sd_hdf5_enable}" = "no"; then
      case "${sd_hdf5_policy}" in
        fail)
          AC_MSG_ERROR([The HDF5 package is required and cannot be disabled
                  See https://tddft.org/programs/hdf5/ for details on how to
                  install it.])
          ;;
        skip)
          tmp_hdf5_invalid="yes"
          ;;
        warn)
          AC_MSG_WARN([The HDF5 package is required and cannot be disabled])
          tmp_hdf5_invalid="yes"
          ;;
      esac
    fi
    if test "${sd_hdf5_enable}" = "auto"; then
      AC_MSG_NOTICE([setting HDF5 trigger to yes])
      sd_hdf5_enable="yes"
    fi
  fi

  # Fix wrong trigger value
  if test "${tmp_hdf5_invalid}" = "yes"; then
    case "${sd_hdf5_status}" in
      implicit|required)
        sd_hdf5_enable="yes"
        ;;
      optional)
        sd_hdf5_enable="no"
        ;;
    esac
    tmp_hdf5_invalid="no"
    AC_MSG_NOTICE([setting sd_hdf5_enable to '${sd_hdf5_enable}'])
  fi

  # Environment variables conflict with --with-* options
  tmp_hdf5_vars="${HDF5_CPPFLAGS}${HDF5_CFLAGS}${HDF5_FCFLAGS}${HDF5_LDFLAGS}${HDF5_LIBS}"
  tmp_hdf5_invalid="no"
  if test ! -z "${tmp_hdf5_vars}" -a ! -z "${sd_hdf5_prefix}"; then
    case "${sd_hdf5_policy}" in
      fail)
        # FIXME: use the new Steredeg specs
        AC_MSG_WARN([conflicting option settings for HDF5
                  Please use HDF5_FCFLAGS + HDF5_LIBS or --with-hdf5,
                  not both.])
        ;;
      skip)
        tmp_hdf5_invalid="yes"
        ;;
      warn)
        AC_MSG_WARN([conflicting option settings for HDF5])
        tmp_hdf5_invalid="yes"
        ;;
    esac
  fi

  # When using environment variables, triggers must be set to yes
  if test -n "${tmp_hdf5_vars}"; then
    sd_hdf5_enable="yes"
    sd_hdf5_init="env"
    if test "${tmp_hdf5_invalid}" = "yes"; then
      tmp_hdf5_invalid="no"
      AC_MSG_NOTICE([overriding --with-hdf5 with HDF5_{FCFLAGS,LDFLAGS,LIBS}])
    fi
  fi

  # Implicit status overrides everything
  if test "${sd_hdf5_status}" = "implicit"; then
    if test "${sd_hdf5_fcflags}" != ""; then
      sd_hdf5_fcflags=""
      AC_MSG_NOTICE([resetting HDF5 Fortran flags (implicit package)])
    fi
    if test "${sd_hdf5_ldflags}" != ""; then
      sd_hdf5_ldflags=""
      AC_MSG_NOTICE([resetting HDF5 linker flags (implicit package)])
    fi
    if test "${sd_hdf5_libs}" != ""; then
      sd_hdf5_libs=""
      AC_MSG_NOTICE([resetting HDF5 library flags (implicit package)])
    fi
  fi

  # Reset build parameters if disabled
  if test "${sd_hdf5_enable}" = "implicit"; then
    sd_hdf5_fcflags=""
    sd_hdf5_ldflags=""
    sd_hdf5_libs=""
  fi

  # Clean-up
  unset tmp_hdf5_invalid
  unset tmp_hdf5_vars
]) # _SD_HDF5_CHECK_CONFIG


AC_DEFUN([_SD_HDF5_CHECK_COMPILERS], [
  if test "${H5CC}" = ""; then
    tmp_h5cc_search="h5cc"
    if test "${sd_mpi_ok}" = "yes"; then
      tmp_h5cc_search="h5pcc ${tmp_h5cc_search}"
    fi
    case "${sd_hdf5_init}" in
      dir)
        AC_MSG_CHECKING([for a HDF5-aware compiler])
        for tmp_h5cc in ${tmp_h5cc_search}; do
          if test -x "${sd_hdf5_prefix}/bin/${tmp_h5cc}"; then
            sd_hdf5_h5cc="${sd_hdf5_prefix}/bin/${tmp_h5cc}"
            break
          fi
        done
        if test "${sd_hdf5_h5cc}" = ""; then
          AC_MSG_RESULT([none])
        else
          AC_MSG_RESULT([${sd_hdf5_h5cc}])
        fi
        ;;
      *)
        AC_CHECK_PROGS([sd_hdf5_h5cc], [${tmp_h5cc_search}])
        ;;
    esac
    unset tmp_h5cc_search
    unset tmp_h5cc
  else
    ${H5CC} -show >/dev/null 2>&1
    if test "${?}" != "0"; then
      AC_MSG_ERROR([invalid HDF5-aware compiler:
                  H5CC = '${H5CC}'])
    fi
    sd_hdf5_h5cc="${H5CC}"
  fi   # H5CC = ""

  if test "${sd_hdf5_h5cc}" != ""; then

    AC_MSG_CHECKING([for the HDF5 install prefix])
    if test "${sd_hdf5_prefix}" = ""; then
      sd_hdf5_prefix=`${sd_hdf5_h5cc} -showconfig | grep 'Installation point:' | awk '{print [$]NF}'`
    fi
    if test "${sd_hdf5_prefix}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_prefix}])
    fi

    AC_MSG_CHECKING([for the HDF5 build flavor])
    sd_hdf5_flavor=`${sd_hdf5_h5cc} -showconfig | grep 'Flavor name:' | awk '{print [$]NF}'`
    if test "${sd_hdf5_flavor}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_flavor}])
    fi

    if test "${sd_hdf5_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([whether HDF5 supports Fortran])
      sd_hdf5_has_fc=`${sd_hdf5_h5cc} -showconfig | grep 'Fortran:' | awk '{print [$]NF}'`
      test "${sd_hdf5_has_fc}" = "" && sd_hdf5_has_fc="no"
      AC_MSG_RESULT([${sd_hdf5_has_fc}])
    else
      sd_hdf5_has_fc="no"
    fi

    AC_MSG_CHECKING([for HDF5 high-level libraries])
    tmp_hdf5_hl=`${sd_hdf5_h5cc} -showconfig | grep 'High@<:@ \-@:>@@<:@Ll@:>@evel library:' | awk '{print [$]NF}'`
    if test "${tmp_hdf5_hl}" = "yes"; then
      sd_hdf5_libs_hl="-lhdf5_hl"
      test "${sd_hdf5_enable_fc}" = "yes" -a "${sd_hdf5_has_fc}" = "yes" && \
        sd_hdf5_libs_hl="-lhdf5hl_fortran ${sd_hdf5_libs_hl}"
    fi
    if test "${sd_hdf5_libs_hl}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_libs_hl}])
    fi
    unset tmp_hdf5_hl

    AC_MSG_CHECKING([for HDF5 extra dependencies])
    sd_hdf5_libs_extra=`${sd_hdf5_h5cc} -showconfig | grep 'Extra libraries: ' | sed -e 's/.*Extra libraries: //'`
    if test "${sd_hdf5_libs_extra}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_libs_extra}])
    fi

    if test "${sd_mpi_ok}" = "yes"; then
      AC_MSG_CHECKING([whether HDF5 supports MPI])
      sd_hdf5_has_mpi=`${sd_hdf5_h5cc} -showconfig | grep 'Parallel HDF5:' | awk '{print [$]NF}'`
      test "${sd_hdf5_has_mpi}" = "" && sd_hdf5_has_mpi="no"
      AC_MSG_RESULT([${sd_hdf5_has_mpi}])
    else
      sd_hdf5_has_mpi="no"
    fi

    tmp_hdf5_incs=""
    tmp_hdf5_libdirs=""
    for arg in `${sd_hdf5_h5cc} -show`; do
      case ${arg} in
        -I*)
          tmp_hdf5_incs="${tmp_hdf5_incs} ${arg}"
          ;;
        -L*)
          tmp_hdf5_libdirs="${tmp_hdf5_libdirs} ${arg}"
        ;;
      esac
    done
    sd_hdf5_cppflags="${sd_hdf5_cppflags} ${tmp_hdf5_incs}"
    test "${sd_hdf5_enable_fc}" = "yes" && \
      sd_hdf5_fcflags="${sd_hdf5_fcflags} ${tmp_hdf5_incs}"
    sd_hdf5_libs="${tmp_hdf5_libdirs} ${sd_hdf5_libs_hl} ${sd_hdf5_libs} ${sd_hdf5_libs_extra}"
    unset arg
    unset tmp_hdf5_incs
    unset tmp_hdf5_libdirs

  fi   # sd_hdf5_cc != ""
]) # _SD_HDF5_CHECK_COMPILERS


AC_DEFUN([_SD_HDF5_DUMP_CONFIG], [
  AC_MSG_CHECKING([whether to enable HDF5])
  AC_MSG_RESULT([${sd_hdf5_enable}])
  if test "${sd_hdf5_enable}" != "no"; then
    AC_MSG_CHECKING([how HDF5 parameters have been set])
    AC_MSG_RESULT([${sd_hdf5_init}])
    AC_MSG_CHECKING([for HDF5 C preprocessing flags])
    if test "${sd_hdf5_cppflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_cppflags}])
    fi
    AC_MSG_CHECKING([for HDF5 C flags])
    if test "${sd_hdf5_cflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_cflags}])
    fi
    AC_MSG_CHECKING([for HDF5 C++ flags])
    if test "${sd_hdf5_cxxflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_cxxflags}])
    fi
    if test "${sd_hdf5_enable_fc}" = "yes"; then
      AC_MSG_CHECKING([for HDF5 Fortran flags])
      if test "${sd_hdf5_fcflags}" = ""; then
        AC_MSG_RESULT([none])
      else
        AC_MSG_RESULT([${sd_hdf5_fcflags}])
      fi
    fi
    AC_MSG_CHECKING([for HDF5 linker flags])
    if test "${sd_hdf5_ldflags}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_ldflags}])
    fi
    AC_MSG_CHECKING([for HDF5 library flags])
    if test "${sd_hdf5_libs}" = ""; then
      AC_MSG_RESULT([none])
    else
      AC_MSG_RESULT([${sd_hdf5_libs}])
    fi
  fi
]) # _SD_HDF5_DUMP_CONFIG

AC_DEFUN([_SD_HDF5_SET_ENV], [
  sd_hdf5_cppflags="${sd_hdf5_cppflags_def}"
  sd_hdf5_cflags="${sd_hdf5_cflags_def}"
  sd_hdf5_cxxflags="${sd_hdf5_cxxflags_def}"
  test "${sd_hdf5_enable_fc}" = "yes" && \
    sd_hdf5_fcflags="${sd_hdf5_fcflags_def}"
  sd_hdf5_ldflags="${sd_hdf5_ldflags_def}"
  sd_hdf5_libs="${sd_hdf5_libs_def}"
  test ! -z "${HDF5_CPPFLAGS}" && sd_hdf5_cppflags="${HDF5_CPPFLAGS}"
  test ! -z "${HDF5_CFLAGS}" && sd_hdf5_cflags="${HDF5_CFLAGS}"
  test ! -z "${HDF5_CXXFLAGS}" && sd_hdf5_cxxflags="${HDF5_CXXFLAGS}"
  if test "${sd_hdf5_enable_fc}" = "yes"; then
    test ! -z "${HDF5_FCFLAGS}" && sd_hdf5_fcflags="${HDF5_FCFLAGS}"
  fi
  test ! -z "${HDF5_LDFLAGS}" && sd_hdf5_ldflags="${HDF5_LDFLAGS}"
  test ! -z "${HDF5_LIBS}" && sd_hdf5_libs="${HDF5_LIBS}"
]) # _SD_HDF5_SET_ENV

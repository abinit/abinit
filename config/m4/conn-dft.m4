# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for Density-Functional Theory libraries
#



# _ABI_DFT_CHECK_ATOMPAW()
# ------------------------
#
# Check whether the Atompaw library is working.
#
AC_DEFUN([_ABI_DFT_CHECK_ATOMPAW],[
  dnl Init
  abi_dft_atompaw_has_bins="no"
  abi_dft_atompaw_has_incs="no"
  abi_dft_atompaw_has_libs="no"
  abi_dft_atompaw_serial="no"
  abi_dft_atompaw_mpi="no"
  abi_dft_atompaw_fcflags=""
  abi_dft_atompaw_ldflags=""
  abi_dft_atompaw_bins="${with_atompaw_bins}"
  abi_dft_atompaw_incs="${with_atompaw_incs}"
  abi_dft_atompaw_libs="${with_atompaw_libs}"

  dnl Look for binaries
  if test "${abi_dft_atompaw_bins}" = ""; then
    AC_CHECK_PROGS([ATOMPAW_BIN],[atompaw])
    AC_CHECK_PROGS([GRAPHATOM_BIN],[graphatom])
    if test "${ATOMPAW_BIN}" != "" -a "${GRAPHATOM_BIN}" != ""; then
      abi_dft_atompaw_has_bins="yes"
    fi
  else
    ATOMPAW_BIN="${abi_dft_atompaw_bins}/atompaw"
    GRAPHATOM_BIN="${abi_dft_atompaw_bins}/graphatom"
    if test -x "${ATOMPAW_BIN}" -a -x "${GRAPHATOM_BIN}"; then
      abi_dft_atompaw_has_bins="yes"
    fi
  fi

  dnl Take final decision for the serial case
  if test "${abi_dft_atompaw_has_bins}" = "yes"; then
    abi_dft_atompaw_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_dft_atompaw_serial}" = "yes"; then
      abi_dft_atompaw_mpi="yes"
    fi
  fi
]) # _ABI_DFT_CHECK_ATOMPAW



# _ABI_DFT_CHECK_BIGDFT()
# -----------------------
#
# Check whether the BigDFT library is working.
#
AC_DEFUN([_ABI_DFT_CHECK_BIGDFT],[
  dnl Init
  abi_dft_bigdft_default_libs="-lbigdft-1 -labinit -lyaml"
  abi_dft_bigdft_has_incs="no"
  abi_dft_bigdft_has_libs="no"
  abi_dft_bigdft_serial="no"
  abi_dft_bigdft_mpi="no"
  abi_dft_bigdft_fcflags=""
  abi_dft_bigdft_ldflags=""
  abi_dft_bigdft_incs="${with_bigdft_incs}"
  abi_dft_bigdft_libs="${with_bigdft_libs}"

  dnl Prepare environment
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_linalg_incs} ${abi_trio_netcdf_incs} ${abi_trio_etsf_io_incs} ${abi_dft_libxc_incs} ${abi_dft_bigdft_incs}"
  LIBS="${lib_libxc_libs} ${lib_etsf_io_libs} ${lib_netcdf_libs} ${lib_linalg_libs} ${LIBS}"
  if test "${with_bigdft_libs}" = ""; then
    AC_MSG_CHECKING([for BigDFT libraries to try])
    LIBS="${abi_dft_bigdft_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_dft_bigdft_default_libs}])
  else
    LIBS="${abi_dft_bigdft_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([bigdft_api])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_dft_bigdft_has_incs="yes"
  fi

  dnl Look for libraries and routines
  AC_MSG_CHECKING([whether BigDFT libraries work])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use bigdft_api
      implicit none
      integer iproc
      type(input_variables) :: inputs
      type(atoms_data) :: atoms
      type(restart_objects) :: rst
      call init_restart_objects(iproc,inputs,atoms,rst)
    ]])], [abi_dft_bigdft_has_libs="yes"], [abi_dft_bigdft_has_libs="no"])
  AC_MSG_RESULT([${abi_dft_bigdft_has_libs}])

  dnl Take final decision for the serial case
  if test "${abi_dft_bigdft_has_incs}" = "yes" -a \
          "${abi_dft_bigdft_has_libs}" = "yes"; then
    abi_dft_bigdft_serial="yes"
    if test "${with_bigdft_libs}" = ""; then
      abi_dft_bigdft_libs="${abi_dft_bigdft_default_libs}"
    fi
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_dft_bigdft_serial}" = "yes"; then
      abi_dft_bigdft_mpi="yes"
    fi
  fi

  dnl Restore environment
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_DFT_CHECK_BIGDFT



# _ABI_DFT_CHECK_LIBXC(API_MAJOR_MIN, API_MINOR_MIN, API_MAJOR_MAX, API_MINOR_MAX)
# ---------------------------------------------------
#
# Check whether the LibXC library is working.
#
AC_DEFUN([_ABI_DFT_CHECK_LIBXC],[
  dnl Init
  abi_dft_libxc_has_incs="no"
  abi_dft_libxc_has_libs="no"
  abi_dft_libxc_version="no"
  abi_dft_libxc_serial="no"
  abi_dft_libxc_mpi="no"
  abi_dft_libxc_fcflags=""
  abi_dft_libxc_ldflags=""
  abi_dft_libxc_incs="${with_libxc_incs}"
  abi_dft_libxc_libs="${with_libxc_libs}"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_dft_libxc_incs}"
  FCFLAGS="${FCFLAGS} ${abi_dft_libxc_incs}"
  LIBS="${abi_dft_libxc_libs} ${LIBS}"
  AC_LANG_PUSH([C])

  dnl Look for C includes
  AC_CHECK_HEADERS([xc.h xc_funcs.h xc_version.h],[abi_dft_libxc_has_incs="yes"],[abi_dft_libxc_has_incs="no"])

  dnl Look for libraries and routines
  if test "${abi_dft_libxc_libs}" = ""; then
    AC_SEARCH_LIBS([xc_func_init],[xc dft_xc],[abi_dft_libxc_has_libs="yes"],[abi_dft_libxc_has_libs="yes"],[-lm])
    if test "${abi_dft_libxc_has_libs}" = "yes"; then
      if test "${ac_cv_search_xc_func_init}" != "none required"; then
        abi_dft_libxc_libs="${ac_cv_search_xc_func_init}"
      fi
    fi
    LIBS="${abi_dft_libxc_libs} ${LIBS}"
  else
    AC_CHECK_LIB([m],[main],abi_has_libm="yes",abi_has_libm="no")
    if test "${abi_has_libm}" = "yes"; then
      abi_dft_libxc_libs="${abi_dft_libxc_libs} -lm"
      with_libxc_libs="${with_libxc_libs} -lm"
      LIBS="${LIBS} -lm"
    fi
    AC_CHECK_LIB([xc],[xc_func_init],[abi_dft_libxc_has_libs="yes"],[abi_dft_libxc_has_libs="no"])
  fi

  dnl Check whether the C wrappers work
  if test "${abi_dft_libxc_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether LibXC is usable])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([[#include "xc.h"]],
      [[xc_func_type func;
        int func_id = 1;
        int i=xc_func_init(&func, func_id, XC_UNPOLARIZED);
      ]])], [abi_dft_libxc_has_libs="yes"], [abi_dft_libxc_has_libs="no"])
    AC_MSG_RESULT([${abi_dft_libxc_has_libs}])
  fi

  dnl Check that we have the correct LibXC version
  if test "${abi_dft_libxc_has_incs}" = "yes" -a \
          "${abi_dft_libxc_has_libs}" = "yes"; then
    AC_MSG_CHECKING([whether this is LibXC version $1.$2->$3.$4])
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
      [[#include "xc_version.h"]],
      [[int ver=100*XC_MAJOR_VERSION+XC_MINOR_VERSION;
        int ver_min=100*$1+$2,ver_max=100*$3+$4;
        if ( (ver<ver_min) || (ver>ver_max)) {return 1;}
      ]])], [abi_dft_libxc_version="yes"], [abi_dft_libxc_version="no"])
    AC_MSG_RESULT([${abi_dft_libxc_version}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_dft_libxc_has_incs}" = "yes" -a \
          "${abi_dft_libxc_has_libs}" = "yes" -a \
          "${abi_dft_libxc_version}" = "yes"; then
    abi_dft_libxc_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_dft_libxc_serial}" = "yes"; then
      abi_dft_libxc_mpi="yes"
    fi
  fi

  dnl Restore environment
  AC_LANG_POP([C])
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"

  dnl Make sure LIBS are properly set for the other packages
  if test "${abi_dft_libxc_serial}" = "yes"; then
    FCFLAGS="${FCFLAGS} ${abi_dft_libxc_incs}"
    LIBS="${abi_dft_libxc_libs} ${LIBS}"
  fi
]) # _ABI_DFT_CHECK_LIBXC



# _ABI_DFT_CHECK_WANNIER90()
# --------------------------
#
# Check whether the Wannier90 library is working.
#
AC_DEFUN([_ABI_DFT_CHECK_WANNIER90],[
  dnl Init
  abi_dft_wannier90_has_bins="no"
  abi_dft_wannier90_has_incs="no"
  abi_dft_wannier90_has_libs="no"
  abi_dft_wannier90_serial="no"
  abi_dft_wannier90_mpi="no"
  abi_dft_wannier90_v1="no"
  abi_dft_wannier90_fcflags=""
  abi_dft_wannier90_ldflags=""
  abi_dft_wannier90_bins="${with_wannier90_bins}"
  abi_dft_wannier90_incs="${with_wannier90_incs}"
  abi_dft_wannier90_libs="${with_wannier90_libs}"

  dnl Prepare environment
  tmp_saved_LIBS="${LIBS}"
  LIBS="${abi_dft_wannier90_libs} ${lib_linalg_libs} ${LIBS}"

  dnl Look for binaries
  if test "${abi_dft_wannier90_bins}" = ""; then
    AC_CHECK_PROGS([WANNIER90_X_BIN],[wannier.x wannier90.x])
    if test "${WANNIER90_X_BIN}" != ""; then
      abi_dft_wannier90_has_bins="yes"
    fi
    AC_CHECK_PROGS([W90CHK2CHK_X_BIN],[w90chk2chk.x])
    if test "${W90CHK2CHK_X_BIN}" = ""; then
      abi_dft_wannier90_v1="yes"
    fi
  else
    AC_MSG_CHECKING([for wannier.x])
    WANNIER90_X_BIN="${abi_dft_wannier90_bins}/wannier.x"
    if test -x "${WANNIER90_X_BIN}"; then
      abi_dft_wannier90_has_bins="yes"
    fi
    AC_MSG_RESULT([${abi_dft_wannier90_has_bins}])
    if test "${abi_dft_wannier90_has_bins}" != "yes"; then
      AC_MSG_CHECKING([for wannier90.x])
      WANNIER90_X_BIN="${abi_dft_wannier90_bins}/wannier90.x"
      if test -x "${WANNIER90_X_BIN}"; then
        abi_dft_wannier90_has_bins="yes"
      fi
      AC_MSG_RESULT([${abi_dft_wannier90_has_bins}])
    fi
    AC_MSG_CHECKING([for wannier90 v1.x (absence of w90chk2chk.x)])
    W90CHK2CHK_X_BIN="${abi_dft_wannier90_bins}/w90chk2chk.x"
    if test ! -x "${W90CHK2CHK_X_BIN}"; then
      abi_dft_wannier90_v1="yes"
    fi
    AC_MSG_RESULT([${abi_dft_wannier90_v1}])
  fi

  dnl Look for libraries and routines
  if test "${abi_dft_wannier90_libs}" = ""; then
    AC_SEARCH_LIBS([wannier_run],[wannier wannier90],
      [abi_dft_wannier90_has_incs="yes"; abi_dft_wannier90_has_libs="yes"],
      [abi_dft_wannier90_has_libs="no"])
    if test "${abi_dft_wannier90_has_libs}" = "yes"; then
      if test "${ac_cv_search_wannier_run}" != "none required"; then
        abi_dft_wannier90_libs="${ac_cv_search_wannier_run}"
      fi
    fi
  else
    AC_MSG_CHECKING([whether the specified Wannier90 library works])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        call wannier_run
      ]])],
      [abi_dft_wannier90_has_incs="yes"; abi_dft_wannier90_has_libs="yes"],
      [abi_dft_wannier90_has_libs="no"])
    AC_MSG_RESULT([${abi_dft_wannier90_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_dft_wannier90_has_bins}" = "yes" -a \
          "${abi_dft_wannier90_has_incs}" = "yes" -a \
          "${abi_dft_wannier90_has_libs}" = "yes"; then
    abi_dft_wannier90_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes"; then
    if test "${abi_dft_wannier_serial}" = "yes"; then
      abi_dft_wannier_mpi="yes"
    fi
  fi

  dnl Restore environment
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_DFT_CHECK_WANNIER90



# ABI_CONNECT_DFT()
# -----------------
#
# Sets all variables needed to handle the DFT libraries.
#
# FIXME: AtomPAW and BigDFT may or may not depend on LibXC. If they have
#        built externally, there is no means for now of controlling
#        whether they have LibXC support. They should output this info
#        along with their version number.
#
AC_DEFUN([ABI_CONNECT_DFT],[
  dnl Initial setup
  abi_test_atompaw="no"
  abi_test_bigdft="no"
  abi_test_libxc="no"
  abi_test_wannier90="no"
  lib_dft_flavor="${with_dft_flavor}"
  tmp_dft_flavor=""

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  LDFLAGS="${FC_LDFLAGS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested DFT support])
  AC_MSG_RESULT([${with_dft_flavor}])

  dnl Look for external DFT libraries
  if test "${with_dft_flavor}" != "none"; then

    dnl Make sure LibXC is looked for before the others
    abi_dft_iter=`echo "${lib_dft_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s",[$]1}'`
    abi_dft_tmp="${abi_dft_iter}"
    for abi_dft_flavor in ${abi_dft_iter}; do
      if test "${abi_dft_flavor}" = "libxc"; then
        abi_dft_tmp=`echo "${abi_dft_iter}" | sed -e 's/libxc//'`
        abi_dft_tmp="libxc ${abi_dft_tmp}"
      fi
      if test "${abi_dft_flavor}" = "libxc-fallback"; then
        abi_dft_tmp=`echo "${abi_dft_iter}" | sed -e 's/libxc-fallback//'`
        abi_dft_tmp="libxc-fallback ${abi_dft_tmp}"
      fi
    done
    abi_dft_iter="${abi_dft_tmp}"

    for abi_dft_flavor in ${abi_dft_iter}; do

      dnl Check if the user has requested a fallback
      tmp_dft_base_flavor=`echo "${abi_dft_flavor}" | cut -d- -f1`
      AC_MSG_CHECKING([whether to select a fallback for ${tmp_dft_base_flavor}])
      tmp_dft_fallback=`echo "${abi_dft_flavor}" | cut -s -d- -f2`
      if test "${tmp_dft_fallback}" = "fallback"; then
        tmp_dft_fallback="yes"
      else
        tmp_dft_fallback="no"
      fi
      AC_MSG_RESULT([${tmp_dft_fallback}])
      if test "${tmp_dft_fallback}" = "yes" -a \
              "${enable_fallbacks}" = "no"; then
        AC_MSG_ERROR([fallback requested while fallbacks have been globally disabled])
      fi

      dnl Look for DFT libraries
      case "${abi_dft_flavor}" in

        atompaw)
          if test "${abi_linalg_serial}" = "yes" -a \
                  "${abi_dft_libxc_serial}" = "yes"; then
            abi_dft_atompaw_prereqs="yes"
            _ABI_DFT_CHECK_ATOMPAW
          else
            if test "${abi_linalg_serial}" != "yes"; then
              AC_MSG_WARN([AtomPAW requires missing linear algebra support])
            fi
            if test "${abi_dft_libxc_serial}" != "yes"; then
              AC_MSG_WARN([AtomPAW requires missing LibXC support])
            fi
            if test \( "${abi_dft_libxc_serial}" = "no" -a \
                       "${abi_dft_libxc_fallback}" != "yes" \) -o \
                    \( "${abi_linalg_serial}" = "no" -a \
                       "${abi_dft_linalg_fallback}" != "yes" \); then
              abi_dft_atompaw_prereqs="no"
            fi
            abi_dft_atompaw_serial="no"
            abi_dft_atompaw_mpi="no"
          fi
          if test "${abi_dft_atompaw_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_ATOMPAW],1,
              [Define to 1 if you have the AtomPAW library.])
            abi_test_atompaw="yes"
          fi
          if test "${abi_dft_atompaw_serial}" = "yes"; then
            lib_atompaw_bins="${abi_dft_atompaw_bins}"
            lib_atompaw_incs="${abi_dft_atompaw_incs}"
            lib_atompaw_libs="${abi_dft_atompaw_libs}"
          fi
          ;;

        atompaw-fallback)
          if test "${abi_linalg_serial}" != "yes" -a \
                  "${abi_linalg_fallback}" != "yes"; then
            AC_MSG_WARN([AtomPAW requires missing linear algebra support])
            abi_dft_atompaw_prereqs="no"
            abi_dft_atompaw_serial="no"
            abi_dft_atompaw_mpi="no"
          fi
          if test "${abi_dft_libxc_serial}" != "yes" -a \
                  "${abi_dft_libxc_fallback}" != "yes"; then
            AC_MSG_WARN([AtomPAW requires missing LibXC support])
            abi_dft_atompaw_prereqs="no"
            abi_dft_atompaw_serial="no"
            abi_dft_atompaw_mpi="no"
          fi
          ;;

        bigdft)
          if test "${abi_dft_libxc_serial}" = "yes" -a \
                  "${abi_linalg_serial}" = "yes"; then
            abi_dft_bigdft_prereqs="yes"
            _ABI_DFT_CHECK_BIGDFT
          else
            if test "${abi_linalg_serial}" != "yes"; then
              AC_MSG_WARN([BigDFT requires missing linear algebra support])
            fi
            if test "${abi_dft_libxc_serial}" != "yes"; then
              AC_MSG_WARN([BigDFT requires missing LibXC support])
            fi
            if test \( "${abi_dft_libxc_serial}" = "no" -a \
                       "${abi_dft_libxc_fallback}" != "yes" \) -o \
                    \( "${abi_linalg_serial}" = "no" -a \
                       "${abi_dft_linalg_fallback}" != "yes" \); then
              abi_dft_bigdft_prereqs="no"
            fi
            abi_dft_bigdft_serial="no"
            abi_dft_bigdft_mpi="no"
          fi
          if test "${abi_dft_bigdft_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_BIGDFT],1,
              [Define to 1 if you have the BigDFT library.])
            abi_test_bigdft="yes"
          fi
          if test "${abi_dft_bigdft_serial}" = "yes"; then
            lib_bigdft_incs="${abi_dft_bigdft_incs}"
            lib_bigdft_libs="${abi_dft_bigdft_libs}"
          fi
          ;;

        bigdft-fallback)
          if test "${abi_linalg_serial}" != "yes" -a \
                  "${abi_linalg_fallback}" != "yes"; then
            AC_MSG_WARN([BigDFT requires missing linear algebra support])
            abi_dft_bigdft_prereqs="no"
            abi_dft_bigdft_serial="no"
            abi_dft_bigdft_mpi="no"
          fi
          if test "${abi_dft_libxc_serial}" != "yes" -a \
                  "${abi_dft_libxc_fallback}" != "yes"; then
            AC_MSG_WARN([BigDFT requires missing LibXC support])
            abi_dft_bigdft_prereqs="no"
            abi_dft_bigdft_serial="no"
            abi_dft_bigdft_mpi="no"
          fi
          ;;

        libxc)
          _ABI_DFT_CHECK_LIBXC(2,2,4,3)
          if test "${abi_dft_libxc_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_LIBXC],1,
              [Define to 1 if you have the LibXC library.])
            abi_test_libxc="yes"
          fi
          if test "${abi_dft_libxc_serial}" = "yes"; then
            lib_libxc_incs="${abi_dft_libxc_incs}"
            lib_libxc_libs="${abi_dft_libxc_libs}"
          elif test "${enable_fallbacks}" = "yes"; then
            abi_dft_libxc_fallback="yes"
          fi
          ;;

        wannier90)
          dnl Wannier90 requires linear algebra support
          if test "${abi_linalg_serial}" = "yes"; then
            abi_dft_wannier90_prereqs="yes"
            _ABI_DFT_CHECK_WANNIER90
          else
            AC_MSG_WARN([wannier90 requires missing linear algebra support])
            if test "${abi_dft_linalg_fallback}" != "yes"; then
              abi_dft_wannier90_prereqs="no"
            fi
            abi_dft_wannier90_serial="no"
            abi_dft_wannier90_mpi="no"
          fi
          if test "${abi_dft_wannier90_serial}" = "yes" -o \
                  "${enable_fallbacks}" = "yes"; then
            AC_DEFINE([HAVE_WANNIER90],1,
              [Define to 1 if you have the Wannier90 library.])
            abi_test_wannier90="yes"
          fi
          if test "${abi_dft_wannier90_serial}" = "yes"; then
            lib_wannier90_bins="${abi_dft_wannier90_bins}"
            lib_wannier90_incs="${abi_dft_wannier90_incs}"
            lib_wannier90_libs="${abi_dft_wannier90_libs}"
            if test "${abi_dft_wannier90_v1}" = "yes"; then
              AC_DEFINE([HAVE_WANNIER90_V1],1,
                [Define to 1 if you have the Wannier90 V1.x library.])
            fi
          fi
          ;;

        wannier90-fallback)
          if test "${abi_linalg_serial}" != "yes" -a \
                  "${abi_linalg_fallback}" != "yes"; then
            AC_MSG_WARN([Wannier90 requires missing linear algebra support])
            abi_dft_wannier90_prereqs="no"
            abi_dft_wannier90_serial="no"
            abi_dft_wannier90_mpi="no"
            AC_DEFINE([HAVE_WANNIER90_V1],0,
              [Define to 1 if you have the Wannier90 V1.x library.])
          fi
          ;;

        *)
          if test "${tmp_dft_fallback}" = "no"; then
            AC_MSG_ERROR([unknown DFT flavor '${abi_dft_flavor}'])
          fi
          ;;

      esac

      dnl Rebuild actual flavor
      if test "${tmp_dft_fallback}" = "yes"; then
        eval abi_dft_${tmp_dft_base_flavor}_fallback="yes"
        abi_fallbacks="${abi_fallbacks} ${tmp_dft_base_flavor}"
        tmp_dft_flavor="${tmp_dft_flavor}+${abi_dft_flavor}"
        tmp_dft_prereqs=`eval echo \$\{abi_dft_${tmp_dft_base_flavor}_prereqs\}`
        if test "${tmp_dft_prereqs}" = "no"; then
          ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
          AC_MSG_ERROR([prerequisites for ${abi_dft_flavor} not found])
        fi
      else
        tmp_dft_prereqs=`eval echo \$\{abi_dft_${abi_dft_flavor}_prereqs\}`
        tmp_dft_serial=`eval echo \$\{abi_dft_${abi_dft_flavor}_serial\}`
        tmp_dft_incs=`eval echo \$\{with_${abi_dft_flavor}_bins\}`
        tmp_dft_incs=`eval echo \$\{with_${abi_dft_flavor}_incs\}`
        tmp_dft_libs=`eval echo \$\{with_${abi_dft_flavor}_libs\}`
        AC_MSG_WARN([package: ${abi_dft_flavor} - preq=${tmp_dft_prereqs} - working=${tmp_dft_serial} - incs=${tmp_dft_incs} libs=${tmp_dft_libs}])
        if test "${tmp_dft_serial}" = "no"; then
          if test "${tmp_dft_bins}" = "" -a \
                  "${tmp_dft_incs}" = "" -a \
                  "${tmp_dft_libs}" = "" -a \
                  "${tmp_dft_prereqs}" != "no"; then
            AC_MSG_WARN([falling back to internal ${abi_dft_flavor} version])
            abi_fallbacks="${abi_fallbacks} ${abi_dft_flavor}"
            tmp_dft_flavor="${tmp_dft_flavor}+${abi_dft_flavor}-fallback"
          else
            ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
            AC_MSG_ERROR([external ${abi_dft_flavor} support does not work])
          fi
        else
          tmp_dft_flavor="${tmp_dft_flavor}+${abi_dft_flavor}"
        fi
      fi

    done

  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  if test "${tmp_dft_flavor}" != ""; then
    lib_dft_flavor=`echo "${tmp_dft_flavor}" | sed -e 's/^\+//;s/\+[$]//'`
  fi
  AC_MSG_CHECKING([for the actual DFT support])
  AC_MSG_RESULT([${lib_dft_flavor}])

  dnl Substitute variables needed for the use of the libraries
  AC_SUBST(lib_dft_flavor)
  AC_SUBST(lib_atompaw_bins)
  AC_SUBST(lib_atompaw_fcflags)
  AC_SUBST(lib_atompaw_incs)
  AC_SUBST(lib_atompaw_ldflags)
  AC_SUBST(lib_atompaw_libs)
  AC_SUBST(lib_bigdft_fcflags)
  AC_SUBST(lib_bigdft_incs)
  AC_SUBST(lib_bigdft_ldflags)
  AC_SUBST(lib_bigdft_libs)
  AC_SUBST(lib_libxc_fcflags)
  AC_SUBST(lib_libxc_incs)
  AC_SUBST(lib_libxc_ldflags)
  AC_SUBST(lib_libxc_libs)
  AC_SUBST(lib_wannier90_bins)
  AC_SUBST(lib_wannier90_fcflags)
  AC_SUBST(lib_wannier90_incs)
  AC_SUBST(lib_wannier90_ldflags)
  AC_SUBST(lib_wannier90_libs)
]) # ABI_CONNECT_DFT

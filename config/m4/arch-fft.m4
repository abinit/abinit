# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for FFT libraries
#



# _ABI_FFT_CHECK_LIB(LIBRARY,ROUTINE)
# -----------------------------------
#
# Check whether the specified library contains the specified routine.
#
AC_DEFUN([_ABI_FFT_CHECK_LIB],[
  dnl Init
  abi_fft_has_incs="no"
  abi_fft_has_libs="no"
  abi_fft_ok="no"
  abi_fft_mpi="no"

  dnl Prepare environment
  tmp_saved_LIBS="${LIBS}"
  LIBS="${abi_fft_libs} ${LIBS}"

  dnl FFT usually doesn't require any include (may change in the future)
  abi_fft_has_incs="yes"

  dnl Look for libraries and routines
  if test "${abi_fft_libs}" = ""; then
    AC_LANG_PUSH([Fortran])
    AC_SEARCH_LIBS($2,$1,[abi_fft_has_libs="yes"], [abi_fft_has_libs="no"])
    AC_LANG_POP([Fortran])
    if test "${abi_fft_has_libs}" = "yes"; then
      tmp_fft_libs=`eval echo \$\{ac_cv_search_$2\}`
      if test "${tmp_fft_libs}" != "none required"; then
        abi_fft_libs="${tmp_fft_libs}"
      fi
    fi
  else
    dnl Try the subroutine version
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [
        call $2
      ])], [abi_fft_has_libs="yes"], [abi_fft_has_libs="no"])
    AC_LANG_POP([Fortran])

    dnl Try the function version
    if test "${abi_fft_has_libs}" = "no"; then
      AC_LANG_PUSH([Fortran])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [
          integer :: xstatus
          xstatus = $2
        ])], [abi_fft_has_libs="yes"], [abi_fft_has_libs="no"])
      AC_LANG_POP([Fortran])
    fi
  fi

  dnl Take final decision
  if test "${abi_fft_has_incs}" = "yes" -a \
          "${abi_fft_has_libs}" = "yes"; then
    abi_fft_ok="yes"
  fi

  dnl Restore environment
  LIBS="${tmp_saved_LIBS}"
]) # _ABI_FFT_CHECK_LIB



# _ABI_FFTW3_CHECK_LIBS()
# -----------------------
#
# Check whether the FFTW3 library is working.
#
AC_DEFUN([_ABI_FFTW3_CHECK_LIBS],[
  dnl Init
  abi_fftw3_ok="unknown"
  abi_fftw3_has_hdrs="unknown"
  abi_fftw3_has_libs="unknown"
  abi_fftw3_has_mods="unknown"
  abi_fftw3_has_fort="unknown"
  abi_fftw3_has_mpi="unknown"
  abi_fftw3_has_omp="unknown"
  abi_fftw3_has_thr="unknown"
  abi_fftw3_fcflags=""
  abi_fftw3_ldflags=""
  abi_fftw3_incs="${with_fft_incs}"
  abi_fftw3_libs="${with_fft_libs}"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${abi_fftw3_incs}"
  FCFLAGS="${FCFLAGS} ${abi_fftw3_incs}"
  LIBS="${abi_fftw3_libs} ${LIBS}"

  dnl Look for C includes
  LDFLAGS="${CC_LDFLAGS}"
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([fftw3.h],
    [abi_fftw3_has_hdrs="yes"], [abi_fftw3_has_hdrs="no"])
  AC_LANG_POP([C])

  dnl Look for C libraries and routines
  if test "${abi_fftw3_has_hdrs}" = "yes"; then
    LDFLAGS="${CC_LDFLAGS}"
    AC_LANG_PUSH([C])
    if test "${abi_fftw3_libs}" = ""; then
      AC_SEARCH_LIBS([fftw_execute], [fftw3],
        [abi_fftw3_has_libs="c"], [abi_fftw3_has_libs="no"])
      if test "${ac_cv_search_fftw_execute}" != "no"; then
        if test "${ac_cv_search_fftw_execute}" != "none required"; then
          abi_fftw3_libs="${ac_cv_search_fftw_execute}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified FFTW3 C libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
#include <fftw3.h>
        ]],
        [[
          fftw_plan *plan;
          fftw_complex *a1, *a2;
          dfftw_execute_dft(plan, a1, a2);
          return 0;
        ]])], [abi_fftw3_has_libs="c"], [abi_fftw3_has_libs="no"])
      tmp_fftw3_has_libs=`echo "${abi_fftw3_has_libs}" | sed -e 's/^c$/yes/'`
      AC_MSG_RESULT([${tmp_fftw3_has_libs}])
      unset tmp_fftw3_has_libs
    fi
    AC_LANG_POP([C])
  fi

  dnl Look for Fortran libraries and routines
  if test "${abi_fftw3_has_libs}" = "c"; then
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    if test "${abi_fftw3_libs}" = ""; then
      AC_SEARCH_LIBS([sfftw_execute_dft], [fftw3f],
        [abi_fftw3_has_libs="yes"], [abi_fftw3_has_libs="no"])
      if test "${ac_cv_search_sfftw_execute_dft}" != "no"; then
        if test "${ac_cv_search_sfftw_execute_dft}" != "none required"; then
          abi_fftw3_libs="${ac_cv_search_sfftw_execute_dft} ${abi_fftw3_libs}"
        fi
      fi
    else
      AC_MSG_CHECKING([whether specified FFTW3 Fortran libraries work])
      AC_LINK_IFELSE([AC_LANG_PROGRAM([],
        [[
          call sfftw_execute_dft
        ]])], [abi_fftw3_has_libs="yes"], [abi_fftw3_has_libs="no"])
      AC_MSG_RESULT([${abi_fftw3_has_libs}])
    fi
    AC_LANG_POP([Fortran])
  fi
  # FIXME: _ABI_FFT_CHECK_LIB([fftw3_threads],[dfftw_init_threads])

  dnl Look for Fortran includes
  dnl Note: must be done after the libraries have been discovered
  if test "${abi_fftw3_has_libs}" = "yes"; then
    AC_MSG_CHECKING([whether we have FFTW3 Fortran include files])
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use iso_c_binding
        include "fftw3.f03"
      ]])], [abi_fftw3_has_mods="yes"], [abi_fftw3_has_mods="no"])
    if test "${abi_fftw3_has_mods}" = "no"; then
      FCFLAGS="${FCFLAGS} -I/usr/include"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
        [[
          use iso_c_binding
          include "fftw3.f03"
        ]])], [abi_fftw3_has_mods="yes"], [abi_fftw3_has_mods="no"])
      if test "${abi_fftw3_has_mods}" = "yes" -a \
              "${abi_fftw3_incs}" = ""; then
        abi_fftw3_incs="-I/usr/include"
      fi
    fi
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_fftw3_has_mods}])
  fi

  dnl Check Fortran support
  if test "${abi_fftw3_has_mods}" = "yes"; then
    AC_MSG_CHECKING([whether FFTW3 Fortran wrappers work])
    LDFLAGS="${FC_LDFLAGS}"
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use iso_c_binding
        include "fftw3.f03"
        integer, parameter :: N = 10
        complex :: a1(N), a2(N)
        integer :: plan
   
        call sfftw_plan_dft_1d(plan, N, a1, a2, FFTW_FORWARD, FFTW_ESTIMATE)
        call sfftw_execute_dft(plan, a1, a2)
        call sfftw_destroy_plan(plan)
      ]])], [abi_fftw3_has_fort="yes"], [abi_fftw3_has_fort="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_fftw3_has_fort}])
  fi

  dnl Set status of minimum FFTW3 support
  if test "${abi_fftw3_has_hdrs}" = "yes" -a \
          "${abi_fftw3_has_libs}" = "yes" -a \
          "${abi_fftw3_has_mods}" = "yes" -a \
          "${abi_fftw3_has_fort}" = "yes"; then
    abi_fftw3_ok="yes"
  else
    abi_fftw3_ok="no"
  fi

  dnl Handle the FFTW3 over MKL case
  AC_MSG_NOTICE([(DEBUG) vari = ${abi_fftw3_variant}])
  AC_MSG_NOTICE([(DEBUG) hdrs = ${abi_fftw3_has_hdrs}])
  AC_MSG_NOTICE([(DEBUG) libs = ${abi_fftw3_has_libs}])
  if test "${abi_fftw3_variant}" = "mkl" -a \
          "${abi_fftw3_has_hdrs}" = "yes" -a \
          "${abi_fftw3_has_libs}" = "yes"; then
    abi_fftw3_ok="yes"
    AC_MSG_WARN([FFT support detection may be incomplete with MKL
                    please check that the final configuration actually
                    corresponds to your situation])
  fi

  dnl Check for MPI support
  #if test "${abi_mpi_enable}" = "yes" -a \
  #        "${abi_fftw3_ok}" = "yes"; then
  #  AC_MSG_CHECKING([whether FFTW3 supports MPI])
  #  AC_LANG_PUSH([Fortran])
  #  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
  #    [[
  #      use fftw3
  #      character(len=*), parameter :: path = "dummy"
  #      integer :: cmode, comm, info, ncerr, ncid
  #      ncerr = sfftw_execute_dft_par(path, cmode, comm, info, ncid)
  #    ]])], [abi_fftw3_has_mpi="yes"], [abi_fftw3_has_mpi="no"])
  #  AC_LANG_POP([Fortran])
  #  AC_MSG_RESULT([${abi_fftw3_has_mpi}])
  #fi

  dnl Check for threads support
  if test "${abi_threads_enable}" = "yes" -a \
          "${abi_fftw3_ok}" = "yes"; then
    AC_MSG_CHECKING([whether FFTW3 supports threads])
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
#include <fftw3.h>
        ]],
      [[
        dfftw_init_threads;
      ]])], [abi_fftw3_has_thr="yes"], [abi_fftw3_has_thr="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${abi_fftw3_has_thr}])
  fi

  dnl Restore environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE
]) # _ABI_FFTW3_CHECK_LIBS



# ABI_FFT_DETECT(FLAVOR)
# ----------------------
#
# Detect FFT libraries according to the specified FFT flavor.
#
AC_DEFUN([ABI_FFT_DETECT],[
  dnl Initial setup
  abi_fft_fcflags=""
  abi_fft_ldflags=""
  abi_fft_incs=""
  abi_fft_libs=""
  abi_fft_ok="unknown"
  abi_fft_mpi="unknown"
  abi_fft_omp="unknown"
  abi_fft_thr="unknown"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${abi_fft_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${abi_fft_libs} ${LIBS}"

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested FFT support])
  AC_MSG_RESULT([${abi_fft_flavor}])

  dnl Look for external FFT libraries
  if test "${abi_fft_flavor}" != "none"; then

    case "${abi_fft_flavor}" in

      custom)
        if test "${abi_fft_libs}" == ""; then
          AC_MSG_ERROR([you must specify custom FFT libraries (--with-fft-libs)])
        fi
        abi_fft_ok="yes"
        abi_fft_serial="yes"
        abi_fft_mpi="yes"
        ;;

      dfti)
        _ABI_FFT_CHECK_LIB([dfti],[DftiCreateDescriptor])
        if test "${abi_fft_ok}" = "yes"; then
          AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
          abi_fft_serial="yes"
          abi_fft_mpi="no"
          abi_fft_fcflags="${abi_fft_fcflags}"
          abi_fft_ldflags="${abi_fft_ldflags}"
          abi_fft_incs="${abi_fft_incs}"
          abi_fft_libs="${abi_fft_libs}"
        fi
        ;;

      dfti-threads)
        _ABI_FFT_CHECK_LIB([dfti],[DftiCreateDescriptor])
        if test "${abi_fft_ok}" = "yes"; then
          AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
          AC_DEFINE([HAVE_FFT_DFTI_THREADS],1,[Define to 1 if you want to use the threaded DFTI library.])
          abi_fft_serial="yes"
          abi_fft_mpi="no"
          abi_fft_fcflags="${abi_fft_fcflags}"
          abi_fft_ldflags="${abi_fft_ldflags}"
          abi_fft_incs="${abi_fft_incs}"
          abi_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3|fftw3-*)
        abi_fftw3_variant=`echo "${abi_fft_flavor}" | cut -d- -f2-`
        _ABI_FFTW3_CHECK_LIBS
        abi_fft_ok="${abi_fftw3_ok}"
        abi_fft_mpi="${abi_fftw3_mpi}"
        abi_fft_omp="${abi_fftw3_omp}"
        abi_fft_thr="${abi_fftw3_thr}"
        if test "${abi_fftw3_ok}" = "yes"; then
          AC_DEFINE([HAVE_FFTW3],1,
            [Define to 1 if you want to use the FFTW3 library.])
          if test "${abi_fftw3_variant}" = "mpi" -a \
                  "${abi_fftw3_mpi}" = "yes"; then
            AC_DEFINE([HAVE_FFTW3_MPI], 1,
              [Define to 1 if you want to use the MPI FFTW3 features.])
          fi
          if test "${abi_fftw3_variant}" = "openmp" -a \
                  "${abi_fftw3_omp}" = "yes"; then
            AC_DEFINE([HAVE_FFTW3_OPENMP], 1,
              [Define to 1 if you want to use the OpenMP FFTW3 features.])
          fi
          if test "${abi_fftw3_variant}" = "threads" -a \
                  "${abi_fftw3_thr}" = "yes"; then
            AC_DEFINE([HAVE_FFTW3_THREADS], 1,
              [Define to 1 if you want to use the threaded FFTW3 features.])
          fi
          abi_fft_fcflags="${abi_fftw3_fcflags}"
          abi_fft_ldflags="${abi_fftw3_ldflags}"
          abi_fft_incs="${abi_fftw3_incs}"
          abi_fft_libs="${abi_fftw3_libs}"
        fi
        ;;

      fftw3-threads)
        _ABI_FFT_CHECK_LIB([fftw3],[dfftw_execute])	
        if test "${abi_fft_ok}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3],1,[Define to 1 if you want to use the FFTW3 library.])
          AC_DEFINE([HAVE_FFT_FFTW3_THREADS],1,[Define to 1 if you want to use the threaded FFTW3 library.])
          abi_fft_serial="yes"
          abi_fft_mpi="no"
          abi_fft_fcflags="${abi_fft_fcflags}"
          abi_fft_ldflags="${abi_fft_ldflags}"
          abi_fft_incs="${abi_fft_incs}"
          abi_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3-mpi)
	_ABI_FFT_CHECK_LIB([fftw3],[dfftw_execute])	
        _ABI_FFT_CHECK_LIB([fftw3_threads],[dfftw_init_threads])
        #_ABI_FFT_CHECK_LIB([fftw3_mpi],[fftw_mpi_init])
        if test "${abi_fft_ok}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3],1,[Define to 1 if you want to use the FFTW3 library.])
          AC_DEFINE([HAVE_FFT_FFTW3_THREADS],1,[Define to 1 if you want to use the threaded FFTW3 library.])	
          AC_DEFINE([HAVE_FFT_FFTW3_MPI],1,[Define to 1 if you want to use the distributed FFTW3 library.])	
          abi_fft_serial="yes"
          abi_fft_mpi="yes"
          abi_fft_fcflags="${abi_fft_fcflags}"
          abi_fft_ldflags="${abi_fft_ldflags}"
          abi_fft_incs="${abi_fft_incs}"
          abi_fft_libs="${abi_fft_libs}"
        fi
        ;;

      *)
        AC_MSG_ERROR([unknown FFT flavor '${abi_fft_flavor}'])
        ;;


    esac

  fi

  dnl Propagate build parameters
  if test "${abi_fft_ok}" = "yes"; then
    abi_fft_fcflags="${abi_fft_fcflags}"
    abi_fft_ldflags="${abi_fft_ldflags}"
    abi_fft_incs="${abi_fft_incs}"
    abi_fft_libs="${abi_fft_libs}"
  else
    abi_fft_flavor="broken"
  fi

  dnl Restore build environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual FFT support])
  AC_MSG_RESULT([${abi_fft_flavor}])
  if test "${abi_fft_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    if test "${with_fft_libs}" = ""; then
      AC_MSG_ERROR([the requested ${abi_fft_flavor} FFT flavor is not available])
    else
      AC_MSG_ERROR([the specified FFT libraries do not work])
    fi
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_fft_flavor)
  AC_SUBST(abi_fft_fcflags)
  AC_SUBST(abi_fft_ldflags)
  AC_SUBST(abi_fft_incs)
  AC_SUBST(abi_fft_libs)
]) # ABI_FFT_DETECT



# ABI_FFT_SET_FLAVOR(LINALG)
# --------------------------
#
# Select a FFT flavor according to the available linear algebra flavor.
#
AC_DEFUN([ABI_FFT_SET_FLAVOR],[
  case "$1" in
    mkl)
      abi_fft_flavor="fftw3-mkl"
      ;;
    *)
      if test "${abi_mpi_enable}" = "yes"; then
        abi_fft_flavor="fftw3-mpi"
      elif test "${abi_threads_enable}" = "yes"; then
        abi_fft_flavor="fftw3-threads"
      else
        abi_fft_flavor="fftw3"
      fi
      ;;
  esac
]) # ABI_FFT_SET_FLAVOR

# -*- Autoconf -*-
#
# Copyright (C) 2005-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for external FFT libraries
#



# _ABI_FFT_CHECK(LIBRARY,ROUTINE)
# -------------------------------
#
# Check whether the specified library contains the specified routine.
#
AC_DEFUN([_ABI_FFT_CHECK],[
  dnl Init
  abi_fft_has_incs="no"
  abi_fft_has_libs="no"
  abi_fft_serial="no"
  abi_fft_mpi="no"
  abi_fft_fcflags=""
  abi_fft_ldflags=""
  abi_fft_incs="${with_fft_incs}"
  abi_fft_libs="${with_fft_libs}"

  dnl FFT usually doesn't require any include (may change in the future)
  abi_fft_has_incs="yes"

  dnl Look for libraries and routines
  if test "${abi_fft_libs}" = ""; then
    AC_SEARCH_LIBS($2,$1,[abi_fft_has_libs="yes"])
    if test "${abi_fft_has_libs}" = "yes"; then
      tmp_fft_libs=`eval echo \$\{ac_cv_search_$2\}`
      if test "${tmp_fft_libs}" != "none required"; then
        abi_fft_libs="${tmp_fft_libs}"
      fi
    fi
  else
    tmp_saved_LIBS="${LIBS}"
    LIBS="${with_fft_libs} ${LIBS}"
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [
        call $2
      ])], [abi_fft_has_libs="yes"], [abi_fft_has_libs="no"])

    if test "${abi_fft_has_libs}" = "no"; then
    dnl Try the function version
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [
        integer :: xstatus
        xstatus = $2
      ])], [abi_fft_has_libs="yes"], [abi_fft_has_libs="no"])
    fi

    LIBS="${tmp_saved_LIBS}"
  fi

  dnl Take final decision for the serial case
  if test "${abi_fft_has_incs}" = "yes" -a \
          "${abi_fft_has_libs}" = "yes"; then
    abi_fft_serial="yes"
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes" -a \
          "${abi_fft_serial}" = "yes"; then
    abi_fft_mpi="yes"
  fi
]) # _ABI_FFT_CHECK



# ABI_CONNECT_FFT()
# -------------------
#
# Sets all variables needed to handle the FFT external libraries.
#
AC_DEFUN([ABI_CONNECT_FFT],[
  dnl Initial setup
  lib_fft_flavor="${with_fft_flavor}"
  lib_fft_fcflags=""
  lib_fft_ldflags=""
  lib_fft_incs=""
  lib_fft_libs=""
  abi_fft_serial="no"
  abi_fft_mpi="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  CPPFLAGS="${with_fft_incs} ${CPPFLAGS}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${with_fft_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Display requested flavor
  AC_MSG_CHECKING([for the requested FFT support])
  AC_MSG_RESULT([${with_fft_flavor}])

  dnl Look for external FFT libraries
  if test "${with_fft_flavor}" != "none"; then

    case "${with_fft_flavor}" in

      asl)
        _ABI_FFT_CHECK([asl],[zfc3fb])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_ASL],1,[Define to 1 if you want to use the ASL library for FFT.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      custom)
        if test "${with_fft_libs}" == ""; then
          AC_MSG_ERROR([you must specify custom FFT libraries (--with-fft-libs)])
        fi
        abi_fft_serial="yes"
        abi_fft_mpi="yes"
        lib_fft_incs="${with_fft_incs}"
        lib_fft_libs="${with_fft_libs}"
        ;;

      fftw2)
        _ABI_FFT_CHECK([fftw2],[fftw_execute])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW2],1,[Define to 1 if you want to use the FFTW2 library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw2-threads)
        _ABI_FFT_CHECK([fftw2],[fftw_init_threads])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW2_THREADS],1,[Define to 1 if you want to use the threaded FFTW2 library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3)
        _ABI_FFT_CHECK([fftw3],[dfftw_execute])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3],1,[Define to 1 if you want to use the FFTW3 library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3-mkl)
        _ABI_FFT_CHECK([fftw3],[dfftw_execute])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3_MKL],1,[Define to 1 if you want to use the threaded FFTW3 library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3-threads)
        _ABI_FFT_CHECK([fftw3],[dfftw_execute])	
        _ABI_FFT_CHECK([fftw3_threads],[dfftw_init_threads])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3],1,[Define to 1 if you want to use the FFTW3 library.])
          AC_DEFINE([HAVE_FFT_FFTW3_THREADS],1,[Define to 1 if you want to use the threaded FFTW3 library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      fftw3-mpi)
	_ABI_FFT_CHECK([fftw3],[dfftw_execute])	
        _ABI_FFT_CHECK([fftw3_threads],[dfftw_init_threads])
        #_ABI_FFT_CHECK([fftw3_mpi],[fftw_mpi_init])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_FFTW3],1,[Define to 1 if you want to use the FFTW3 library.])
          AC_DEFINE([HAVE_FFT_FFTW3_THREADS],1,[Define to 1 if you want to use the threaded FFTW3 library.])	
          AC_DEFINE([HAVE_FFT_FFTW3_MPI],1,[Define to 1 if you want to use the distributed FFTW3 library.])	
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      dfti)
        _ABI_FFT_CHECK([dfti],[DftiCreateDescriptor])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      dfti-threads)
        _ABI_FFT_CHECK([dfti],[DftiCreateDescriptor])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
          AC_DEFINE([HAVE_FFT_DFTI_THREADS],1,[Define to 1 if you want to use the threaded DFTI library.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      #dfti-mpi)
      #  _ABI_FFT_CHECK([dfti_mpi],[fftw_mpi_init])
      #  if test "${abi_fft_serial}" = "yes"; then
      #    AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
      #    AC_DEFINE([HAVE_FFT_DFTI_MPI],1,[Define to 1 if you want to use the distributed DFTI library.])	
      #    lib_fft_fcflags="${abi_fft_fcflags}"
      #    lib_fft_ldflags="${abi_fft_ldflags}"
      #    lib_fft_incs="${abi_fft_incs}"
      #    lib_fft_libs="${abi_fft_libs}"
      #  fi
      #  ;;

      #dfti-mpi-threads)
      #  _ABI_FFT_CHECK([dfti],[DftiCreateDescriptor])
      #  #_ABI_FFT_CHECK([dfti_mpi],[fftw_mpi_init])
      #  #_ABI_FFT_CHECK([dfti_threads],[dfftw_init_threads])
      #  if test "${abi_fft_serial}" = "yes"; then
      #    AC_DEFINE([HAVE_FFT_DFTI],1,[Define to 1 if you want to use the DFTI library.])
      #    AC_DEFINE([HAVE_FFT_DFTI_MPI],1,[Define to 1 if you want to use the distributed DFTI library.])	
      #    AC_DEFINE([HAVE_FFT_DFTI_THREADS],1,[Define to 1 if you want to use the threaded DFTI library.])	
      #    lib_fft_fcflags="${abi_fft_fcflags}"
      #    lib_fft_ldflags="${abi_fft_ldflags}"
      #    lib_fft_incs="${abi_fft_incs}"
      #    lib_fft_libs="${abi_fft_libs}"
      #  fi
      #  ;;

      mlib)
        _ABI_FFT_CHECK([veclib],[c1dfft])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_MLIB],1,[Define to 1 if you want to use the HP MLIB library for FFT.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      sgimath)
        _ABI_FFT_CHECK([complib.sgimath],[dfft1du])
        if test "${abi_fft_serial}" = "yes"; then
          AC_DEFINE([HAVE_FFT_SGIMATH],1,[Define to 1 if you want to use the SGIMATH library for FFT.])
          lib_fft_fcflags="${abi_fft_fcflags}"
          lib_fft_ldflags="${abi_fft_ldflags}"
          lib_fft_incs="${abi_fft_incs}"
          lib_fft_libs="${abi_fft_libs}"
        fi
        ;;

      *)
        AC_MSG_ERROR([unknown FFT flavor '${with_fft_flavor}'])
        ;;


    esac

  fi

  dnl Transmit serial status to the source code
  if test "${abi_fft_serial}" = "yes"; then
    AC_DEFINE([HAVE_FFT],1,[Define to 1 if you have an optimized FFT library.])
    AC_DEFINE([HAVE_FFT_SERIAL],1,[Define to 1 if you have an optimized serial FFT library.])
  elif test "${with_fft_flavor}" != "none"; then
    lib_fft_flavor="broken"
  fi

  dnl Transmit MPI status to the source code
  if test "${abi_fft_mpi}" = "yes"; then
    AC_DEFINE([HAVE_FFT_MPI],1,[Define to 1 if you have an optimized MPI-parallel FFT library.])
  fi

  dnl Restore build environment
  AC_LANG_POP([Fortran])
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Output final flavor
  AC_MSG_CHECKING([for the actual FFT support])
  AC_MSG_RESULT([${lib_fft_flavor}])
  if test "${lib_fft_flavor}" = "broken"; then
    ABI_MSG_NOTICE([connectors-failure],[Connector detection failure])
    if test "${with_fft_libs}" = ""; then
      AC_MSG_ERROR([the requested ${with_fft_flavor} FFT flavor is not available])
    else
      AC_MSG_ERROR([the specified FFT libraries do not work])
    fi
  fi

  dnl Substitute variables needed for the use of the library
  AC_SUBST(lib_fft_flavor)
  AC_SUBST(lib_fft_fcflags)
  AC_SUBST(lib_fft_ldflags)
  AC_SUBST(lib_fft_incs)
  AC_SUBST(lib_fft_libs)
]) # ABI_CONNECT_FFT

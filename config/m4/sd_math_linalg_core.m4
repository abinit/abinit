# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon, Marc Torrent)
#
# This file is part of the Steredeg software package. For license information,
# please see the COPYING file in the top-level directory of the source
# distribution.
#

#
# Support for external linear algebra libraries - Core macros
#


# _SD_LINALG_CHECK_LIBS()
# -----------------------
#
# Check whether the specified libraries are BLAS and LAPACK
# implementations.
#
AC_DEFUN([_SD_LINALG_CHECK_LIBS], [
  # Init
  sd_linalg_has_blas="no"
  sd_linalg_has_lapack="no"
  sd_linalg_has_lapacke="no"
  sd_linalg_has_blacs="no"
  sd_linalg_has_scalapack="no"
  sd_linalg_has_elpa="no"
  sd_linalg_has_elpa_2013="no"
  sd_linalg_has_elpa_2014="no"
  sd_linalg_has_elpa_2015="no"
  sd_linalg_has_elpa_2016="no"
  sd_linalg_has_plasma="no"
  sd_linalg_has_magma="no"

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags}"
  FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags}"
  LIBS="${sd_linalg_libs} ${sd_gpu_libs} ${sd_mpi_libs} ${LIBS}"
  AC_LANG_PUSH([Fortran])

  # BLAS?
  AC_MSG_CHECKING([for BLAS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zgemm
    ]])], [sd_linalg_has_blas="yes"], [sd_linalg_has_blas="no"])
  AC_MSG_RESULT([${sd_linalg_has_blas}])

  # BLAS extensions?
  _SD_LINALG_CHECK_BLAS_EXTS()

  # MKL BLAS extensions?
  _SD_LINALG_CHECK_BLAS_MKL_EXTS()

  # LAPACK?
  AC_MSG_CHECKING([for LAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call zhpev
    ]])], [sd_linalg_has_lapack="yes"], [sd_linalg_has_lapack="no"])
  AC_MSG_RESULT([${sd_linalg_has_lapack}])

  # LAPACKE?
  AC_MSG_CHECKING([for LAPACKE C API support in specified libraries])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [#include <lapacke.h>],
    [[
      zhpev_;
    ]])],[sd_linalg_has_lapacke="yes"], [sd_linalg_has_lapacke="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${sd_linalg_has_lapacke}])

  # BLACS?
  AC_MSG_CHECKING([for BLACS support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call blacs_gridinit
    ]])], [sd_linalg_has_blacs="yes"], [sd_linalg_has_blacs="no"])
  AC_MSG_RESULT([${sd_linalg_has_blacs}])

  # ScaLAPACK?
  AC_MSG_CHECKING([for ScaLAPACK support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call pzheevx
    ]])], [sd_linalg_has_scalapack="yes"], [sd_linalg_has_scalapack="no"])
  AC_MSG_RESULT([${sd_linalg_has_scalapack}])

  # ELPA
  AC_MSG_CHECKING([for ELPA support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      use elpa1
      integer,parameter :: n=1, comm=1
      integer :: comm1, comm2, success
      success = get_elpa_communicators(comm, n, n, comm1, comm2)
    ]])], [sd_linalg_has_elpa="yes"], [sd_linalg_has_elpa="no"])
  AC_MSG_RESULT([${sd_linalg_has_elpa}])
  if test "${sd_linalg_has_elpa}" = "yes"; then
    _SD_LINALG_CHECK_ELPA_2017()
    _SD_LINALG_CHECK_ELPA_2016()
    _SD_LINALG_CHECK_ELPA_2015()
    _SD_LINALG_CHECK_ELPA_2014()
    _SD_LINALG_CHECK_ELPA_2013()
  fi

  # PLASMA?
  AC_MSG_CHECKING([for PLASMA support in specified libraries])
  sd_linalg_chk_plasma="${sd_linalg_has_lapacke}"
  if test "${sd_linalg_chk_plasma}" = "yes"; then
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use plasma
        call plasma_zhegv
      ]])], [sd_linalg_has_plasma="yes"], [sd_linalg_has_plasma="no"])
  fi
  AC_MSG_RESULT([${sd_linalg_has_plasma}])

  # MAGMA?
  AC_MSG_CHECKING([for MAGMA (version>=1.1.0) support in specified libraries])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      call magmaf_zhegvd
    ]])], [sd_linalg_has_magma="yes"], [sd_linalg_has_magma="no"])
  AC_MSG_RESULT([${sd_linalg_has_magma}])
  if test "${sd_linalg_has_magma}" = "yes"; then
    _SD_LINALG_CHECK_MAGMA_15()
  fi

  # Restore environment
  AC_LANG_POP([Fortran])
  SD_ESL_RESTORE_FLAGS
]) # _SD_LINALG_CHECK_LIBS


                    # ------------------------------------ #


# _SD_LINALG_EXPLORE()
# --------------------
#
# Looks for linear algebra components by going through all the selected
# vendor sequences.
#
AC_DEFUN([_SD_LINALG_EXPLORE], [
  # Prepare environment
  SD_ESL_SAVE_FLAGS

  # Reset linear algebra build flags, i.e. discard defaults
  sd_linalg_cppflags=""
  sd_linalg_cflags=""
  sd_linalg_cxxflags=""
  sd_linalg_fcflags=""
  sd_linalg_ldflags=""
  sd_linalg_libs=""
  sd_linalg_flavor=""

  # Look for serial linear algebra support
  for tmp_linalg_vendor in ${sd_linalg_chk_serial}; do

    # Configure vendor libraries
    SD_ESL_RESTORE_FLAGS
    CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags}"
    CFLAGS="${CFLAGS} ${sd_linalg_cflags}"
    CXXFLAGS="${CXXFLAGS} ${sd_linalg_cxxflags}"
    FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags}"
    LDFLAGS="${LDFLAGS} ${sd_linalg_ldflags}"
    LIBS="${sd_linalg_libs} ${LIBS}"
    _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])

    # Look for BLAS
    tmp_linalg_blas_proceed=`echo "${sd_linalg_vendor_provided}" | grep "blas"`
    if test "${tmp_linalg_blas_proceed}" != "" -a \
            "${sd_linalg_has_blas}" != "yes"; then
      CPPFLAGS="${CPPFLAGS} ${sd_linalg_vendor_cppflags}"
      CFLAGS="${CFLAGS} ${sd_linalg_vendor_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_linalg_vendor_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_linalg_vendor_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_linalg_vendor_ldflags}"
      LIBS="${sd_linalg_vendor_blas_libs} ${LIBS}"
      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for BLAS])
      if test "${sd_linalg_vendor_blas_libs}" = ""; then
        AC_MSG_RESULT([none required])
      else
        AC_MSG_RESULT([${sd_linalg_vendor_blas_libs}])
      fi
      _SD_LINALG_CHECK_BLAS
      if test "${sd_linalg_has_blas}" = "yes"; then
         _SD_LINALG_CHECK_BLAS_EXTS
         if test "${tmp_linalg_vendor}" = "mkl"; then
           _SD_LINALG_CHECK_BLAS_MKL_EXTS
         fi
         sd_linalg_flavor="${tmp_linalg_vendor}"
         sd_linalg_blas_vendor="${tmp_linalg_vendor}"
         sd_linalg_provided="${sd_linalg_provided} blas"
         sd_linalg_cppflags="${sd_linalg_vendor_cppflags}"
         sd_linalg_cflags="${sd_linalg_vendor_cflags}"
         sd_linalg_cxxflags="${sd_linalg_vendor_cxxflags}"
         sd_linalg_fcflags="${sd_linalg_vendor_fcflags}"
         sd_linalg_ldflags="${sd_linalg_vendor_ldflags}"
         sd_linalg_libs="${sd_linalg_vendor_blas_libs}"
      fi
    fi

    # Look for LAPACK
    tmp_linalg_lapack_proceed=`echo "${sd_linalg_vendor_provided}" | grep "lapack"`
    if test "${tmp_linalg_lapack_proceed}" != "" -a \
            "${sd_linalg_has_blas}" = "yes" -a \
            "${sd_linalg_has_lapack}" != "yes"; then

      AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for LAPACK])
      if test "${sd_linalg_vendor_lapack_libs}" = ""; then
        AC_MSG_RESULT([none required])
      else
        AC_MSG_RESULT([${sd_linalg_vendor_lapack_libs}])
        LIBS="${sd_linalg_vendor_lapack_libs} ${LIBS}"
      fi
      _SD_LINALG_CHECK_LAPACK
      if test "${sd_linalg_has_lapack}" = "yes"; then
        sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
        sd_linalg_lapack_vendor="${tmp_linalg_vendor}"
        sd_linalg_provided="${sd_linalg_provided} lapack"
        if test "${sd_linalg_lapack_vendor}" != "${sd_linalg_blas_vendor}"; then
          test "${sd_linalg_vendor_cppflags}" != "" && \
            sd_linalg_cppflags="${sd_linalg_cppflags} ${sd_linalg_vendor_cppflags}"
          test "${sd_linalg_vendor_cflags}" != "" && \
            sd_linalg_cflags="${sd_linalg_cflags} ${sd_linalg_vendor_cflags}"
          test "${sd_linalg_vendor_cxxflags}" != "" && \
            sd_linalg_cxxflags="${sd_linalg_cxxflags} ${sd_linalg_vendor_cxxflags}"
          test "${sd_linalg_vendor_fcflags}" != "" && \
            sd_linalg_fcflags="${sd_linalg_fcflags} ${sd_linalg_vendor_fcflags}"
          test "${sd_linalg_vendor_ldflags}" != "" && \
            sd_linalg_ldflags="${sd_linalg_ldflags} ${sd_linalg_vendor_ldflags}"
        fi
        test "${sd_linalg_vendor_lapack_libs}" != "" && \
            sd_linalg_libs="${sd_linalg_vendor_lapack_libs} ${sd_linalg_libs}"
        break
      fi
    fi

  done

  # Checkpoint: validate the serial linear algebra support
  if test "${sd_linalg_has_blas}" = "yes" -a \
          "${sd_linalg_has_lapack}" = "yes"; then
    sd_linalg_serial_ok="yes"
  else
    sd_linalg_serial_ok="no"
  fi

  # Look for MPI linear algebra support
  if test "${sd_linalg_serial_ok}" = "yes" -a "${sd_mpi_enable}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_mpi}; do

      # Configure vendor libraries
      SD_ESL_RESTORE_FLAGS
      CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags}"
      CFLAGS="${CFLAGS} ${sd_linalg_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_linalg_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_linalg_ldflags}"
      LIBS="${sd_linalg_libs} ${LIBS}"
      _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])

      # FIXME: implement PLASMA first, with break

      # Look for BLACS
      tmp_linalg_blacs_proceed=`echo "${sd_linalg_vendor_provided}" | grep "blacs"`
      if test "${tmp_linalg_blacs_proceed}" != "" -a \
              "${sd_linalg_serial_ok}" = "yes" -a \
              "${sd_linalg_has_blacs}" != "yes"; then

        AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for BLACS])
        if test "${sd_linalg_vendor_blacs_libs}" = ""; then
          AC_MSG_RESULT([none required])
        else
          AC_MSG_RESULT([${sd_linalg_vendor_blacs_libs}])
          LIBS="${sd_linalg_vendor_blacs_libs} ${LIBS}"
        fi
        _SD_LINALG_CHECK_BLACS
        if test "${sd_linalg_has_blacs}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
          sd_linalg_blacs_vendor="${tmp_linalg_vendor}"
          sd_linalg_provided="${sd_linalg_provided} blacs"
          if test "${sd_linalg_blacs_vendor}" != "${sd_linalg_lapack_vendor}"; then
            test "${sd_linalg_vendor_cppflags}" != "" && \
              sd_linalg_cppflags="${sd_linalg_cppflags} ${sd_linalg_vendor_cppflags}"
            test "${sd_linalg_vendor_cflags}" != "" && \
              sd_linalg_cflags="${sd_linalg_cflags} ${sd_linalg_vendor_cflags}"
            test "${sd_linalg_vendor_cxxflags}" != "" && \
              sd_linalg_cxxflags="${sd_linalg_cxxflags} ${sd_linalg_vendor_cxxflags}"
            test "${sd_linalg_vendor_fcflags}" != "" && \
              sd_linalg_fcflags="${sd_linalg_fcflags} ${sd_linalg_vendor_fcflags}"
            test "${sd_linalg_vendor_ldflags}" != "" && \
              sd_linalg_ldflags="${sd_linalg_ldflags} ${sd_linalg_vendor_ldflags}"
          fi
          test "${sd_linalg_vendor_blacs_libs}" != "" && \
              sd_linalg_libs="${sd_linalg_vendor_blacs_libs} ${sd_linalg_libs}"
        fi
      fi

      # Look for ScaLAPACK
      tmp_linalg_scalapack_proceed=`echo "${sd_linalg_vendor_provided}" | grep "scalapack"`
      if test "${tmp_linalg_scalapack_proceed}" != "" -a \
              "${sd_linalg_serial_ok}" = "yes" -a \
              "${sd_linalg_has_scalapack}" != "yes"; then

        AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for ScaLAPACK])
        if test "${sd_linalg_vendor_scalapack_libs}" = ""; then
          AC_MSG_RESULT([none required])
        else
          AC_MSG_RESULT([${sd_linalg_vendor_scalapack_libs}])
          LIBS="${sd_linalg_vendor_scalapack_libs} ${LIBS}"
        fi
        _SD_LINALG_CHECK_SCALAPACK
        if test "${sd_linalg_has_scalapack}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
          sd_linalg_scalapack_vendor="${tmp_linalg_vendor}"
          sd_linalg_provided="${sd_linalg_provided} scalapack"
          if test "${sd_linalg_scalapack_vendor}" != "${sd_linalg_blacs_vendor}"; then
            test "${sd_linalg_vendor_cppflags}" != "" && \
              sd_linalg_cppflags="${sd_linalg_cppflags} ${sd_linalg_vendor_cppflags}"
            test "${sd_linalg_vendor_cflags}" != "" && \
              sd_linalg_cflags="${sd_linalg_cflags} ${sd_linalg_vendor_cflags}"
            test "${sd_linalg_vendor_cxxflags}" != "" && \
              sd_linalg_cxxflags="${sd_linalg_cxxflags} ${sd_linalg_vendor_cxxflags}"
            test "${sd_linalg_vendor_fcflags}" != "" && \
              sd_linalg_fcflags="${sd_linalg_fcflags} ${sd_linalg_vendor_fcflags}"
            test "${sd_linalg_vendor_ldflags}" != "" && \
              sd_linalg_ldflags="${sd_linalg_ldflags} ${sd_linalg_vendor_ldflags}"
          fi
          test "${sd_linalg_vendor_scalapack_libs}" != "" && \
              sd_linalg_libs="${sd_linalg_vendor_scalapack_libs} ${sd_linalg_libs}"
        fi
      fi

    done
  fi   # sd_mpi_enable = yes

  # FIXME: Look for GPU linear algebra support
  if test "${sd_gpu_enable}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_gpu}; do
      AC_MSG_WARN([GPU linear algebra exploration not implemented!])
    done
  fi

  # Reformat linear algebra flavor
  tmp_linalg_iter=`echo "${sd_linalg_flavor}" | tr '+' '\n' | sort -u | awk '{printf " %s", [$]1} END{printf "\n"}'`
  sd_linalg_flavor=""
  for tmp_linalg_vendor in ${tmp_linalg_iter}; do
    if test "${sd_linalg_flavor}" = ""; then
      sd_linalg_flavor="${tmp_linalg_vendor}"
    else
      sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
    fi
  done

  # Clean-up
  SD_ESL_RESTORE_FLAGS
  unset tmp_linalg_iter
  unset tmp_linalg_vendor
]) # SD_LINALG_EXPLORE


                    # ------------------------------------ #


# _SD_LINALG_SET_VENDOR_FLAGS(VENDOR)
# -----------------------------------
#
# Set libraries to look for depending on the specified flavor.
#
AC_DEFUN([_SD_LINALG_SET_VENDOR_FLAGS], [
  # Reset components
  sd_linalg_vendor_provided=""
  sd_linalg_vendor_cppflags=""
  sd_linalg_vendor_cflags=""
  sd_linalg_vendor_cxxflags=""
  sd_linalg_vendor_fcflags=""
  sd_linalg_vendor_ldflags=""
  sd_linalg_vendor_blacs_libs=""
  sd_linalg_vendor_blacs_prqs=""
  sd_linalg_vendor_blas_libs=""
  sd_linalg_vendor_blas_prqs=""
  sd_linalg_vendor_lapack_libs=""
  sd_linalg_vendor_lapack_prqs=""
  sd_linalg_vendor_lapacke_libs=""
  sd_linalg_vendor_lapacke_prqs=""
  sd_linalg_vendor_scalapack_libs=""
  sd_linalg_vendor_scalapack_prqs=""
  sd_linalg_vendor_elpa_libs=""
  sd_linalg_vendor_elpa_prqs=""
  sd_linalg_vendor_magma_libs=""
  sd_linalg_vendor_magma_prqs=""
  sd_linalg_vendor_plasma_libs=""
  sd_linalg_vendor_plasma_prqs=""

  # Update components according to specified vendor
  case "$1" in

    acml)
      sd_linalg_vendor_provided="blas lapack blacs scalapack"
      sd_linalg_vendor_blas_libs="-lacml -lacml_mv"
      ;;

    atlas)
      sd_linalg_vendor_provided="blas"
      sd_linalg_vendor_blas_libs="-lf77blas -lcblas -latlas"
      ;;

    debian-mpich)
      sd_linalg_vendor_provided="blas lapack blacs scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_blacs_libs="-lblacs-mpich -lblacsCinit-mpich -lblacsF77init-mpich"
      sd_linalg_vendor_scalapack_libs="-lscalapack-mpich"
      ;;

    debian-openmpi)
      sd_linalg_vendor_provided="blas lapack blacs scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_blacs_libs="-lblacs-openmpi -lblacsCinit-openmpi -lblacsF77init-openmpi"
      sd_linalg_vendor_scalapack_libs="-lscalapack-openmpi"
      ;;

    easybuild)
      sd_linalg_vendor_provided="blas lapack blacs scalapack"
      sd_linalg_vendor_blas_libs="-lopenblas"
      sd_linalg_vendor_blacs_libs="-lscalapack"
      ;;

    elpa)
      sd_linalg_vendor_provided="elpa"
      sd_linalg_vendor_elpa_libs="-lelpa"
      ;;

    essl)
      sd_linalg_vendor_provided="blas lapack lapacke blacs scalapack"
      sd_linalg_vendor_fcflags="-qessl"
      sd_linalg_vendor_ldflags="-qessl"
      sd_linalg_vendor_blas_libs="-lessl"
      ;;

    magma)
      sd_linalg_vendor_provided="magma"
      sd_linalg_vendor_magma_libs="-lmagma"
      ;;

    mkl)
      if test "${MKLROOT}" = ""; then
        AC_MSG_ERROR([MKLROOT is not set, which means that MKL is not
                  properly configured])
      fi
      sd_linalg_vendor_provided="blas lapack blacs scalapack"
      sd_linalg_vendor_cppflags="-I${MKLROOT}/include"
      sd_linalg_vendor_fcflags="-I${MKLROOT}/include"
      if test "${sd_mpi_enable}" = "yes"; then
        sd_linalg_vendor_ldflags="-mkl=cluster"
      else
        sd_linalg_vendor_ldflags="-mkl"
      fi
      ;;

    netlib)
      sd_linalg_vendor_provided="blas lapack lapacke blacs scalapack"
      sd_linalg_vendor_blas_libs="-lblas"
      sd_linalg_vendor_lapack_libs="-llapack"
      sd_linalg_vendor_lapacke_libs="-llapacke"
      sd_linalg_vendor_blacs_libs="-lblacs -lblacsCinit -lblacsF77init"
      sd_linalg_vendor_scalapack_libs="-lscalapack"
      ;;

    openblas)
      sd_linalg_vendor_provided="blas"
      sd_linalg_vendor_blas_libs="-lopenblas"
      ;;

    plasma)
      sd_linalg_vendor_provided="plasma"
      sd_linalg_vendor_plasma_libs="-lplasma -lcorelapack -lcoreblas"
      ;;

    *)
      AC_MSG_ERROR([no library settings for linear algebra flavor '$1'])
      ;;

  esac
]) # _SD_LINALG_SET_VENDOR_FLAGS

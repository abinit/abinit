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
  sd_linalg_has_blas="unknown"
  sd_linalg_has_lapack="unknown"
  sd_linalg_has_lapacke="unknown"
  sd_linalg_has_blacs="unknown"
  sd_linalg_has_scalapack="unknown"
  sd_linalg_has_elpa="unknown"
  sd_linalg_has_elpa_2013="unknown"
  sd_linalg_has_elpa_2014="unknown"
  sd_linalg_has_elpa_2015="unknown"
  sd_linalg_has_elpa_2016="unknown"
  sd_linalg_has_elpa_2017="unknown"
  sd_linalg_has_plasma="unknown"
  sd_linalg_has_magma="unknown"
  sd_linalg_provided=""

  # Prepare environment
  SD_ESL_SAVE_FLAGS
  CPPFLAGS="${CPPFLAGS} ${sd_linalg_cppflags}"
  CFLAGS="${CFLAGS} ${sd_linalg_cflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_linalg_cxxflags}"
  FCFLAGS="${FCFLAGS} ${sd_linalg_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_linalg_ldflags}"
  LIBS="${sd_linalg_libs} ${sd_gpu_libs} ${sd_mpi_libs} ${LIBS}"

  # BLAS?
  _SD_LINALG_CHECK_BLAS

  # BLAS extensions?
  _SD_LINALG_CHECK_BLAS_EXTS

  # MKL BLAS extensions?
  _SD_LINALG_CHECK_BLAS_MKL_EXTS

  # LAPACK?
  if test "${sd_linalg_has_blas}" = "yes"; then
    _SD_LINALG_CHECK_LAPACK
  fi

  # LAPACKE?
  if test "${sd_linalg_has_lapack}" = "yes"; then
    _SD_LINALG_CHECK_LAPACKE
  fi

  # Validate the serial linear algebra support
  if test "${sd_linalg_has_blas}" = "yes" -a \
          "${sd_linalg_has_lapack}" = "yes"; then
    sd_linalg_serial_ok="yes"
  else
    sd_linalg_serial_ok="no"
  fi

  if test "${sd_mpi_enable}" = "yes"; then

    # PLASMA?
    _SD_LINALG_CHECK_PLASMA

    if test "${sd_linalg_has_plasma}" != "yes"; then

      # BLACS?
      _SD_LINALG_CHECK_BLACS

      # ScaLAPACK?
      if test "${sd_linalg_has_blacs}" = "yes"; then
        _SD_LINALG_CHECK_SCALAPACK
      fi

    fi

    # ELPA
    _SD_LINALG_CHECK_ELPA

  fi

  # Validate the MPI linear algebra support
  if test "${sd_mpi_enable}" = "yes"; then
    if test \( "${sd_linalg_has_blacs}" = "yes" -a \
               "${sd_linalg_has_scalapack}" = "yes" \) -o \
               "${sd_linalg_has_plasma}" = "yes"; then
      sd_linalg_mpi_ok="yes"
    else
      sd_linalg_mpi_ok="no"
    fi
    if test "${sd_linalg_has_elpa}" = "yes"; then
      sd_linalg_mpiacc_ok="yes"
    else
      sd_linalg_mpiacc_ok="no"
    fi
  fi

  if test "${sd_gpu_enable}" = "yes"; then

    # MAGMA?
    _SD_LINALG_CHECK_MAGMA

  fi

  # Validate the GPU linear algebra support
  if test "${sd_linalg_has_magma}" = "yes"; then
    sd_linalg_gpu_ok="yes"
  else
    sd_linalg_gpu_ok="no"
  fi

  # Compose a flavor from the identified components
  if test "${sd_linalg_flavor_init}" = "def"; then
    if test "${sd_linalg_serial_ok}" = "yes"; then
      if test "${sd_linalg_mkl_has_imatcopy}" = "yes" -o \
              "${sd_linalg_mkl_has_omatcopy}" = "yes" -o \
              "${sd_linalg_mkl_has_omatadd}" = "yes"; then
        sd_linalg_flavor="mkl"
      else
        sd_linalg_flavor="netlib"
      fi
      if test "${sd_linalg_mpi_ok}" = "yes"; then
        if test "${sd_linalg_has_plasma}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+plasma"
        fi
      fi
      if test "${sd_linalg_mpiacc_ok}" = "yes"; then
        if test "${sd_linalg_has_elpa}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+elpa"
        fi
      fi
      if test "${sd_linalg_gpu_ok}" = "yes"; then
        if test "${sd_linalg_has_magma}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+magma"
        fi
      fi
    else
      sd_linalg_flavor="none"
    fi
  fi

  # Restore environment
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

  # FIXME: implement LAPACKE

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
      CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags} ${sd_linalg_cppflags}"
      CFLAGS="${CFLAGS} ${sd_mpi_cflags} ${sd_linalg_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_mpi_cxxflags} ${sd_linalg_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_mpi_fcflags} ${sd_linalg_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags} ${sd_linalg_ldflags}"
      LIBS="${sd_linalg_libs} ${sd_mpi_libs} ${LIBS}"
      _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])

      # Look for PLASMA
      tmp_linalg_plasma_proceed=`echo "${sd_linalg_vendor_provided}" | grep "plasma"`
      if test "${tmp_linalg_plasma_proceed}" != "" -a \
              "${sd_linalg_serial_ok}" = "yes" -a \
              "${sd_linalg_has_plasma}" != "yes"; then

        AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for PLASMA])
        if test "${sd_linalg_vendor_plasma_libs}" = ""; then
          AC_MSG_RESULT([none required])
        else
          AC_MSG_RESULT([${sd_linalg_vendor_plasma_libs}])
          LIBS="${sd_linalg_vendor_plasma_libs} ${LIBS}"
        fi
        _SD_LINALG_CHECK_PLASMA
        if test "${sd_linalg_has_plasma}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
          sd_linalg_plasma_vendor="${tmp_linalg_vendor}"
          sd_linalg_provided="${sd_linalg_provided} plasma"
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
          test "${sd_linalg_vendor_plasma_libs}" != "" && \
              sd_linalg_libs="${sd_linalg_vendor_plasma_libs} ${sd_linalg_libs}"
          break
        fi
      fi

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
          break
        fi
      fi

    done
  fi   # sd_mpi_enable = yes

  # Checkpoint: validate the MPI linear algebra support
  if test "${sd_mpi_enable}" = "yes"; then
    if test \( "${sd_linalg_has_blacs}" = "yes" -a \
               "${sd_linalg_has_scalapack}" = "yes" \) -o \
               "${sd_linalg_has_plasma}" = "yes"; then
      sd_linalg_mpi_ok="yes"
    else
      sd_linalg_mpi_ok="no"
    fi
  fi

  # Look for MPI acceleration linear algebra support
  if test "${sd_linalg_mpi_ok}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_mpiacc}; do

      # Configure vendor libraries
      SD_ESL_RESTORE_FLAGS
      CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags} ${sd_linalg_cppflags}"
      CFLAGS="${CFLAGS} ${sd_mpi_cflags} ${sd_linalg_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_mpi_cxxflags} ${sd_linalg_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_mpi_fcflags} ${sd_linalg_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags} ${sd_linalg_ldflags}"
      LIBS="${sd_linalg_libs} ${sd_mpi_libs} ${LIBS}"
      _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])

      # Look for ELPA
      tmp_linalg_elpa_proceed=`echo "${sd_linalg_vendor_provided}" | grep "elpa"`
      if test "${tmp_linalg_elpa_proceed}" != "" -a \
              "${sd_linalg_has_elpa}" != "yes"; then

        AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for ELPA])
        if test "${sd_linalg_vendor_elpa_libs}" = ""; then
          AC_MSG_RESULT([none required])
        else
          AC_MSG_RESULT([${sd_linalg_vendor_elpa_libs}])
          LIBS="${sd_linalg_vendor_elpa_libs} ${LIBS}"
        fi
        _SD_LINALG_CHECK_ELPA
        if test "${sd_linalg_has_elpa}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
          sd_linalg_elpa_vendor="${tmp_linalg_vendor}"
          sd_linalg_provided="${sd_linalg_provided} elpa"
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
          test "${sd_linalg_vendor_elpa_libs}" != "" && \
              sd_linalg_libs="${sd_linalg_vendor_elpa_libs} ${sd_linalg_libs}"
          break
        fi
      fi

    done
  fi   # sd_linalg_mpi_ok = yes

  # Checkpoint: validate the MPI acceleration linear algebra support
  if test "${sd_mpi_enable}" = "yes"; then
    if test "${sd_linalg_has_elpa}" = "yes"; then
      sd_linalg_mpiacc_ok="yes"
    else
      sd_linalg_mpiacc_ok="no"
    fi
  fi

  # FIXME: Look for GPU linear algebra support
  if test "${sd_gpu_enable}" = "yes"; then
    for tmp_linalg_vendor in ${sd_linalg_chk_gpu}; do

      # Configure vendor libraries
      SD_ESL_RESTORE_FLAGS
      CPPFLAGS="${CPPFLAGS} ${sd_mpi_cppflags} ${sd_gpu_cppflags} ${sd_linalg_cppflags}"
      CFLAGS="${CFLAGS} ${sd_mpi_cflags} ${sd_gpu_cflags} ${sd_linalg_cflags}"
      CXXFLAGS="${CXXFLAGS} ${sd_mpi_cxxflags} ${sd_gpu_cxxflags} ${sd_linalg_cxxflags}"
      FCFLAGS="${FCFLAGS} ${sd_mpi_fcflags} ${sd_gpu_fcflags} ${sd_linalg_fcflags}"
      LDFLAGS="${LDFLAGS} ${sd_mpi_ldflags} ${sd_gpu_ldflags} ${sd_linalg_ldflags}"
      LIBS="${sd_linalg_libs} ${sd_gpu_libs} ${sd_mpi_libs} ${LIBS}"
      _SD_LINALG_SET_VENDOR_FLAGS([${tmp_linalg_vendor}])

      # Look for MAGMA
      tmp_linalg_magma_proceed=`echo "${sd_linalg_vendor_provided}" | grep "magma"`
      if test "${tmp_linalg_magma_proceed}" != "" -a \
              "${sd_linalg_has_magma}" != "yes"; then

        AC_MSG_CHECKING([${tmp_linalg_vendor} libraries for MAGMA])
        if test "${sd_linalg_vendor_magma_libs}" = ""; then
          AC_MSG_RESULT([none required])
        else
          AC_MSG_RESULT([${sd_linalg_vendor_magma_libs}])
          LIBS="${sd_linalg_vendor_magma_libs} ${LIBS}"
        fi
        _SD_LINALG_CHECK_MAGMA
        if test "${sd_linalg_has_magma}" = "yes"; then
          sd_linalg_flavor="${sd_linalg_flavor}+${tmp_linalg_vendor}"
          sd_linalg_magma_vendor="${tmp_linalg_vendor}"
          sd_linalg_provided="${sd_linalg_provided} magma"
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
          test "${sd_linalg_vendor_magma_libs}" != "" && \
              sd_linalg_libs="${sd_linalg_vendor_magma_libs} ${sd_linalg_libs}"
          break
        fi
      fi

    done
  fi   # sd_linalg_mpi_ok = yes

  # Checkpoint: validate the GPU linear algebra support
  if test "${sd_linalg_has_magma}" = "yes"; then
    sd_linalg_gpu_ok="yes"
  else
    sd_linalg_gpu_ok="no"
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
      case "${sd_fc_vendor}" in
        gnu)
          sd_linalg_vendor_cppflags="-I${MKLROOT}/include"
          sd_linalg_vendor_cflags="-m64"
          sd_linalg_vendor_cxxflags="-m64"
          sd_linalg_vendor_fcflags="-m64 -I${MKLROOT}/include"
          sd_linalg_vendor_blas_libs="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl"
          if test "${sd_mpi_enable}" = "yes"; then
            sd_linalg_vendor_blas_libs="-L${MKLROOT}/lib/intel64 -Wl,--no-as-needed -lmkl_scalapack_lp64 -lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl"
          fi
          ;;
        intel)
          sd_linalg_vendor_cppflags="-I${MKLROOT}/include"
          sd_linalg_vendor_fcflags="-I${MKLROOT}/include"
          if test "${sd_mpi_enable}" = "yes"; then
            sd_linalg_vendor_ldflags="-mkl=cluster"
          else
            sd_linalg_vendor_ldflags="-mkl"
          fi
          ;;
        *)
          AC_MSG_ERROR([MKL is only supported with GNU and Intel compilers])
          ;;
      esac
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

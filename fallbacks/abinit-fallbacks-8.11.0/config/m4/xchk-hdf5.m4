# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_HDF5()
# ----------------
#
# Check whether the specified HDF5 library is working.
#
AC_DEFUN([AFB_CHECK_HDF5],[
  dnl Init
  afb_hdf5_default_libs="-lhdf5_hl -lhdf5"
  afb_hdf5_has_par="unknown"
  afb_hdf5_ok="unknown"

  dnl Prepare environment
  tmp_saved_with_hdf5="${with_hdf5}"
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_hdf5_incs}"
  AC_MSG_CHECKING([for HDF5 libraries to try])
  if test "${afb_hdf5_libs}" = ""; then
    LIBS="${afb_hdf5_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_hdf5_default_libs}])
  else
    LIBS="${afb_hdf5_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_hdf5_libs}])
  fi

  dnl Look for prerequisites
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([zlib.h curl/curl.h],
    [afb_hdf5_has_preh="yes"],[afb_hdf5_has_preh="no"])
  AC_SEARCH_LIBS([deflateInit], [z],
    [afb_hdf5_has_zlib="yes"], [afb_hdf5_has_prel="no"])
  AC_SEARCH_LIBS([curl_easy_init], [curl],
    [afb_hdf5_has_prel="yes"], [afb_hdf5_has_prel="no"])
  AC_LANG_POP([C])
  if test "${afb_hdf5_has_preh}" != "yes"; then
    AC_MSG_WARN([missing prerequisite C headers])
  fi
  if test "${afb_hdf5_has_prel}" != "yes"; then
    AC_MSG_WARN([missing prerequisite libraries])
  fi

  dnl Use existing HDF5 macro from Autoconf Archive,
  dnl first parallel, then serial if it fails
  dnl We have to do this because we prefer parallel, while the macro
  dnl prefers serial
  AX_LIB_HDF5([parallel])
  if test "${with_hdf5}" = "no"; then
    AC_MSG_NOTICE([no parallel HDF5 found, looking for a serial one])
    with_hdf5="${tmp_saved_with_hdf5}"
    AX_LIB_HDF5([serial])
  fi

  dnl Set internal variables according to the results of AX_LIB_HDF5
  dnl See ax_lib_hdf5.m4 for more information
  if test "${with_hdf5}" = "yes"; then
    AC_MSG_CHECKING([which HDF5 version we have])
    AC_MSG_RESULT([${HDF5_VERSION}])
    afb_hdf5_ok="yes"
  else
    afb_hdf5_ok="no"
  fi

  dnl Restore environment
  with_hdf5="${tmp_saved_with_hdf5}"
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_HDF5

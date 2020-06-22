# -*- Autoconf -*-
#
# Copyright (C) 2014-2020 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Python support
#



# _ABIBND_CHECK_NUMPY_HEADERS()
# --------------------------
#
# Checks for the existence of NumPy headers.
#
AC_DEFUN([_ABIBND_CHECK_NUMPY_HEADERS],[
  dnl Init
  abi_numpy_ok="no"

  dnl Look for a standard implementation
  if test "${abi_numpy_ok}" = "no"; then
    AC_MSG_CHECKING([for numpy/arrayobject.h])
    AC_LANG_PUSH([C])
    AC_COMPILE_IFELSE(
      [AC_LANG_PROGRAM(
        [[
#include <Python.h>
#include <numpy/arrayobject.h>
        ]],
        [[
        ]]
      )],
      [abi_numpy_ok="yes"],[abi_numpy_ok="no"])
    if test "${abi_numpy_ok}" = "yes"; then
      AC_DEFINE([HAVE_NUMPY],1,[Define to 1 if you have a standard implementation of NumPy.])
    fi
    AC_LANG_POP([C])
    AC_MSG_RESULT([${abi_numpy_ok}])
  fi
]) # _ABIBND_CHECK_NUMPY_HEADERS



# ABIBND_CHECK_PYTHON()
# ------------------
#
# Checks whether the Python environment satisfies the requirements of Abinit.
#
AC_DEFUN([ABIBND_CHECK_PYTHON],
[
  dnl Init
  abi_python_ok="no"
  abi_save_CPPFLAGS="${CPPFLAGS}"
  abi_save_CFLAGS="${CFLAGS}"
  abi_save_LDFLAGS="${LDFLAGS}"
  abi_save_LIBS="${LIBS}"

  dnl Look for Python interpreter
  AC_CHECK_PROGS(PYTHON,
    [python3.7 python3.6 python3.5 python3.4 python3 python2.7 python])
  AC_CHECK_PROGS(PYTHON_CONFIG,
    [python3.7-config python3.6-config python3.5-config python3.4-config python3-config python2.7-config python-config])

  dnl Get Python CPPFLAGS
  if test "${PYTHON}" != "" -a "${PYTHON_CONFIG}" != ""; then
    if test "${PYTHON_CPPFLAGS}" = ""; then
      PYTHON_CPPFLAGS=`${PYTHON_CONFIG} --includes`
    fi
    if test "${PYTHON_CFLAGS}" = ""; then
      PYTHON_CFLAGS=`${PYTHON_CONFIG} --cflags`
    fi
    if test "${PYTHON_LDFLAGS}" = ""; then
      PYTHON_LDFLAGS=`${PYTHON_CONFIG} --ldflags`
    fi
    if test "${PYTHON_LIBS}" = ""; then
      PYTHON_LIBS=`${PYTHON_CONFIG} --libs`
    fi
  fi
  AC_MSG_CHECKING([for Python CPPFLAGS])
  if test "${PYTHON_CPPFLAGS}" = ""; then
    AC_MSG_RESULT([none found])
  else
    AC_MSG_RESULT([${PYTHON_CPPFLAGS}])
  fi
  CPPFLAGS="${PYTHON_CPPFLAGS} ${CPPFLAGS}"
  AC_MSG_CHECKING([for Python CFLAGS])
  if test "${PYTHON_CFLAGS}" = ""; then
    AC_MSG_RESULT([none found])
  else
    AC_MSG_RESULT([${PYTHON_CFLAGS}])
  fi
  CFLAGS="${PYTHON_CFLAGS} ${CFLAGS}"
  AC_MSG_CHECKING([for Python LDFLAGS])
  if test "${PYTHON_LDFLAGS}" = ""; then
    AC_MSG_RESULT([none found])
  else
    AC_MSG_RESULT([${PYTHON_LDFLAGS}])
  fi
  LDFLAGS="${PYTHON_LDFLAGS} ${LDFLAGS}"
  AC_MSG_CHECKING([for Python LIBS])
  if test "${PYTHON_LIBS}" = ""; then
    AC_MSG_RESULT([none found])
  else
    AC_MSG_RESULT([${PYTHON_LIBS}])
  fi
  LIBS="${PYTHON_LIBS} ${LIBS}"

  dnl Preliminary Pyton tests
  AC_CHECK_HEADER([Python.h],[abi_python_ok="yes"])

  dnl Look for Python modules
  if test "${abi_python_ok}" = "yes"; then
    _ABIBND_CHECK_NUMPY_HEADERS
    if test "${abi_numpy_ok}" = "no"; then
      AC_MSG_NOTICE([adding "-I/usr/include/numpy" to CPPFLAGS])
      CPPFLAGS="${CPPFLAGS} -I/usr/include/numpy"
      _ABIBND_CHECK_NUMPY_HEADERS
      if test "${abi_numpy_ok}" = "yes"; then
        PYTHON_CPPFLAGS="${PYTHON_CPPFLAGS} -I/usr/include/numpy"
      fi
    fi
  else
    AC_MSG_WARN([your Python development environment is not working])
  fi

  dnl Restore environment
  CPPFLAGS="${abi_save_CPPFLAGS}"
  CFLAGS="${abi_save_CFLAGS}"
  LDFLAGS="${abi_save_LDFLAGS}"
  LIBS="${abi_save_LIBS}"
]) # ABIBND_CHECK_PYTHON

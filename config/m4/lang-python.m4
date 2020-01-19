# -*- Autoconf -*-
#
# Copyright (C) 2009-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Python support
#



# _ABI_CHECK_NUMPY_HEADERS()
# --------------------------
#
# Checks for the existence of NumPy headers.
#
AC_DEFUN([_ABI_CHECK_NUMPY_HEADERS],[
  # Init
  abi_numpy_ok="no"

  # Look for a standard implementation
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
]) # _ABI_CHECK_NUMPY_HEADERS



# ABI_PY_FEATURES()
# -----------------
#
# Checks whether scientific Python modules are available to other languages.
#
AC_DEFUN([ABI_PY_FEATURES],
[
  # Init
  abi_python_ok="no"

  # Preserve environment
  ABI_ENV_BACKUP

  # Get Python CPPFLAGS
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

  # Preliminary Pyton tests
  AC_CHECK_HEADER([Python.h],[abi_python_ok="yes"])

  # Look for Python modules
  if test "${abi_python_ok}" = "yes"; then
    _ABI_CHECK_NUMPY_HEADERS
    if test "${abi_numpy_ok}" = "no"; then
      AC_MSG_NOTICE([adding "-I/usr/include/numpy" to CPPFLAGS])
      CPPFLAGS="${CPPFLAGS} -I/usr/include/numpy"
      _ABI_CHECK_NUMPY_HEADERS
      if test "${abi_numpy_ok}" = "yes"; then
        PYTHON_CPPFLAGS="${PYTHON_CPPFLAGS} -I/usr/include/numpy"
      fi
    fi
  else
    AC_MSG_WARN([your Python development environment is not working])
  fi

  # Restore environment
  CPPFLAGS="${abi_env_CPPFLAGS}"
  CFLAGS="${abi_env_CFLAGS}"
  LDFLAGS="${abi_env_LDFLAGS}"
  LIBS="${abi_env_LIBS}"
]) # ABI_PY_FEATURES



# ABI_PROG_PYTHON()
# -----------------
#
# Looks for a suitable Python interpreter.
#
AC_DEFUN([ABI_PROG_PYTHON],
[
  # Look for a Python interpreter
  AC_CHECK_PROGS(PYTHON,
    [python3.7 python3.6 python3.5 python3.4 python3 python2.7 python])

  # Look for a Python configurator
  AC_CHECK_PROGS([PYTHON_CONFIG], [${PYTHON}-config])
  if test "${PYTHON_CONFIG}" = ""; then
    AC_CHECK_PROGS(PYTHON_CONFIG,
      [python3.7-config python3.6-config python3.5-config python3.4-config python3-config python2.7-config python-config])
  fi
]) # ABI_PROG_PYTHON

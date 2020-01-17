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
  CPPFLAGS="${PYTHON_CPPFLAGS} ${CPPFLAGS}"

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
]) # ABIBND_CHECK_PYTHON

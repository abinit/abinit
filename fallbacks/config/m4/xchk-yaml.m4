# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_YAML(API_MAJOR_MIN, API_MINOR_MIN)
# ------------------------------------------------------
#
# Check whether the specified YAML library is working.
#
AC_DEFUN([AFB_CHECK_YAML],[
  dnl Init
  afb_yaml_default_libs="-lyaml"
  afb_yaml_has_incs="unknown"
  afb_yaml_has_libs="unknown"
  afb_yaml_ext_ok="unknown"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_yaml_incs}"
  AC_MSG_CHECKING([for YAML libraries to try])
  if test "${afb_yaml_libs}" = ""; then
    LIBS="${afb_yaml_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_yaml_default_libs}])
  else
    LIBS="${afb_yaml_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_yaml_libs}])
  fi

  dnl Look for C includes
  AC_LANG_PUSH([C])
  AC_CHECK_HEADERS([yaml.h],[afb_yaml_has_incs="yes"],[afb_yaml_has_incs="no"])
  AC_LANG_POP([C])

  dnl Check that the library is working
  if test "${afb_yaml_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether YAML is working])
    AC_LANG_PUSH([C])
    AC_RUN_IFELSE([AC_LANG_PROGRAM(
      [[
#include <stdio.h>
#include "yaml.h"
      ]],
      [[
        yaml_parser_t parser;

        if(!yaml_parser_initialize(&parser))
          fputs("Failed to initialize parser!\n", stderr);
        yaml_parser_delete(&parser);
      ]])], [afb_yaml_has_libs="yes"], [afb_yaml_has_libs="no"])
    AC_LANG_POP([C])
    AC_MSG_RESULT([${afb_yaml_has_libs}])
  fi

  dnl Final adjustments
  if test "${afb_yaml_has_incs}" = "yes" -a \
          "${afb_yaml_has_libs}" = "yes"; then
    afb_yaml_ext_ok="yes"
  else
    afb_yaml_ext_ok="no"
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_YAML

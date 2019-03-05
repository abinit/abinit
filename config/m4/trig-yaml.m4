# -*- Autoconf -*-
#
# Copyright (C) 2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Support for the YAML library
#



# ABI_TRIGGER_YAML()
# ------------------
#
# Check whether the YAML library is working.
#
AC_DEFUN([ABI_TRIGGER_YAML],[
  dnl Init
  abi_yaml_default_libs="-lfyaml"
  abi_yaml_has_incs="no"
  abi_yaml_has_libs="no"
  abi_yaml_serial="no"
  abi_yaml_mpi="no"
  abi_yaml_fcflags=""
  abi_yaml_ldflags=""
  abi_yaml_incs="${with_yaml_incs}"
  abi_yaml_libs="${with_yaml_libs}"
  abi_test_yaml="no"

  dnl Prepare environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  FCFLAGS="${FCFLAGS} ${abi_yaml_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  if test "${with_yaml_libs}" = ""; then
    AC_MSG_CHECKING([for YAML libraries to try])
    LIBS="${abi_yaml_default_libs} ${LIBS}"
    AC_MSG_RESULT([${abi_yaml_default_libs}])
  else
    LIBS="${abi_yaml_libs} ${LIBS}"
  fi

  dnl Look for includes
  ABI_FC_MOD_INCS([yaml_output])
  FCFLAGS="${FCFLAGS} ${fc_mod_incs}"
  if test "${abi_fc_mod_incs_ok}" != "unknown"; then
    abi_yaml_has_incs="yes"
  fi

  dnl Look for libraries and routines
  if test "${abi_yaml_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether the YAML library works])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use yaml_output
        use dictionaries
        type(dictionary), pointer :: dict
        call yaml_new_document()
        call dict_init(dict)
        call set(dict//'tot_ncpus', 4)
      ]])], [abi_yaml_has_libs="yes"], [abi_yaml_has_libs="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_yaml_has_libs}])
  fi

  dnl Take final decision for the serial case
  if test "${abi_yaml_has_incs}" = "yes" -a \
          "${abi_yaml_has_libs}" = "yes"; then
    abi_yaml_serial="yes"
    AC_DEFINE([HAVE_YAML],1,
      [Define to 1 if you have the YAML library.])
    abi_test_yaml="yes"
    abi_yaml_incs="${abi_yaml_incs}"
    abi_yaml_libs="${abi_yaml_libs}"
  else
    ABI_MSG_NOTICE([connectors-failure],[YAML detection failure])
    AC_MSG_ERROR([the YAML library is absent or unusable])
  fi

  dnl Check for MPI support
  if test "${enable_mpi}" = "yes" -a \
          "${abi_yaml_serial}" = "yes"; then
    abi_yaml_mpi="yes"
  fi

  dnl Restore environment
  LIBS="${abi_saved_LIBS}"
  ABI_ENV_RESTORE

  dnl Substitute variables needed for the use of the library
  AC_SUBST(abi_yaml_fcflags)
  AC_SUBST(abi_yaml_incs)
  AC_SUBST(abi_yaml_ldflags)
  AC_SUBST(abi_yaml_libs)
]) # ABI_TRIGGER_YAML

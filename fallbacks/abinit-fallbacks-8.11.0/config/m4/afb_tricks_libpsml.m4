# -*- Autoconf -*-
#
# Copyright (C) 2006-2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Tricks for external packages
#



# AFB_TRICKS_LIBPSML(FC_VENDOR,FC_VERSION)
# -------------------------------------
#
# Applies tricks and workarounds to have the LibPSML library correctly
# linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_LIBPSML],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_libpsml_tricks="no"
  afb_libpsml_tricky_vars=""
  tmp_libpsml_num_tricks=0
  tmp_libpsml_cnt_tricks=0

  AC_MSG_NOTICE([no tricks to apply for LibPSML])

  dnl Count applied tricks
  case "${tmp_libpsml_cnt_tricks}" in
    0)
      afb_libpsml_tricks="no"
      ;;
    ${tmp_libpsml_num_tricks})
      afb_libpsml_tricks="yes"
      ;;
    *)
      afb_libpsml_tricks="partial"
      ;;
  esac
  unset tmp_libpsml_cnt_tricks
  unset tmp_libpsml_num_tricks
]) # AFB_TRICKS_LIBPSML

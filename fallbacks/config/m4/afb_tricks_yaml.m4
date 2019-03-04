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



# AFB_TRICKS_YAML(FC_VENDOR,FC_VERSION)
# --------------------------------------
#
# Applies tricks and workarounds to have the YAML library correctly
# linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_YAML],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], [], [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_yaml_tricks="no"
  afb_yaml_tricky_vars=""
  tmp_yaml_num_tricks=1
  tmp_yaml_cnt_tricks=0

  dnl Configure tricks
  if test "${afb_yaml_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying YAML tricks (vendor: $1, version: $2, flags: config)])
    dnl Internal YAML parameters
    CFGFLAGS_YAML="--enable-static --disable-shared"

    dnl Finish
    tmp_yaml_cnt_tricks=`expr ${tmp_yaml_cnt_tricks} \+ 1`
    afb_yaml_tricky_vars="${afb_yaml_tricky_vars} CFGFLAGS"
  else
    AC_MSG_NOTICE([CFGFLAGS_YAML set => skipping YAML config tricks])
  fi

  dnl Count applied tricks
  case "${tmp_yaml_cnt_tricks}" in
    0)
      afb_yaml_tricks="no"
      ;;
    ${tmp_yaml_num_tricks})
      afb_yaml_tricks="yes"
      ;;
    *)
      afb_yaml_tricks="partial"
      ;;
  esac
  unset tmp_yaml_cnt_tricks
  unset tmp_yaml_num_tricks
]) # AFB_TRICKS_YAML

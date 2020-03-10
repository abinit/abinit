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



# AFB_TRICKS_WANNIER90(FC_VENDOR,FC_VERSION)
# ------------------------------------------
#
# Applies tricks and workarounds to have the Wannier90 bindings correctly
# linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_WANNIER90],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_wannier90_tricks="no"
  afb_wannier90_tricky_vars=""
  tmp_wannier90_num_tricks=2
  tmp_wannier90_cnt_tricks=0

  dnl Configure tricks
  if test "${afb_wannier90_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying Wannier90 tricks (vendor: $1, version: $2, flags: config)])

    dnl Internal Wannier90 parameters
    dnl Note: Disable shared libraries, because linear algebra fallback is
    dnl only providing static libraries.
    CFGFLAGS_WANNIER90="${CFGFLAGS_WANNIER90} --disable-shared --enable-static"

    dnl Finish
    tmp_wannier90_cnt_tricks=`expr ${tmp_wannier90_cnt_tricks} \+ 1`
    afb_wannier90_tricky_vars="${afb_wannier90_tricky_vars} CFGFLAGS"
  else
    AC_MSG_NOTICE([CFGFLAGS_WANNIER90 set => skipping Wannier90 config tricks])
  fi

  dnl Libraries tricks
  if test "${afb_wannier90_libs_custom}" = "no"; then
    AC_MSG_NOTICE([applying Wannier90 tricks (vendor: $1, version: $2, flags: libraries)])

    dnl Linear algebra
    tmplibs_wannier90='$(afb_linalg_libs)'
    LIBS_WANNIER90="${tmplibs_wannier90} ${LIBS_WANNIER90}"
    unset tmplibs_wannier90

    case "$1" in
      intel)
        case "${target_cpu}" in
          ia64)
            # Do nothing
            ;;
          *)
            case "$2" in
              9.0|9.1|10.0|10.1|11.0|11.1)
              LIBS_WANNIER90="${LIBS_WANNIER90} -lsvml"
              ;;
            esac
            ;;
        esac
        ;;
    esac

    dnl Finish
    tmp_wannier90_cnt_tricks=`expr ${tmp_wannier90_cnt_tricks} \+ 1`
    afb_wannier90_tricky_vars="${afb_wannier90_tricky_vars} LIBS"
  else
    AC_MSG_NOTICE([LIBS_WANNIER90 set => skipping Wannier90 libraries tricks])
  fi

  dnl Count applied tricks
  case "${tmp_wannier90_cnt_tricks}" in
    0)
      afb_wannier90_tricks="no"
      ;;
    ${tmp_wannier90_num_tricks})
      afb_wannier90_tricks="yes"
      ;;
    *)
      afb_wannier90_tricks="partial"
      ;;
  esac
  unset tmp_wannier90_cnt_tricks
  unset tmp_wannier90_num_tricks
]) # AFB_TRICKS_WANNIER90

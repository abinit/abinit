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



# AFB_TRICKS_BIGDFT(FC_VENDOR,FC_VERSION)
# ---------------------------------------
#
# Applies tricks and workarounds to have the BigDFT library correctly
# linked to the binaries.
#
AC_DEFUN([AFB_TRICKS_BIGDFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  m4_if([$2], , [AC_FATAL([$0: missing argument 2])])dnl

  dnl Init
  afb_bigdft_tricks="no"
  afb_bigdft_tricky_vars=""
  tmp_bigdft_num_tricks=3
  tmp_bigdft_cnt_tricks=0

  dnl Configure tricks
  if test "${afb_bigdft_cfgflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying BigDFT tricks (vendor: $1, version: $2, flags: config)])

    dnl LibXC
    tmpflags_libxc='--disable-internal-libxc --with-libxc-incs="$(afb_libxc_incs)" --with-libxc-libs="$(afb_libxc_libs) ${LIBS}"'

    dnl YAML
    dnl FIXME: disabled internal YAML because PyYAML requires shared objects

    dnl Internal BigDFT parameters
    tmpflags_options='--without-archives --with-moduledir="$(prefix)/$(bigdft_pkg_inst)/include"'
    tmpflags_bigdft='--disable-binaries --disable-bindings --enable-libbigdft'
    CFGFLAGS_BIGDFT="${CFGFLAGS_BIGDFT} ${tmpflags_bigdft} ${tmpflags_options} ${tmpflags_libxc}"

    dnl Finish
    tmp_bigdft_cnt_tricks=`expr ${tmp_bigdft_cnt_tricks} \+ 1`
    afb_bigdft_tricky_vars="${afb_bigdft_tricky_vars} CFGFLAGS"
    unset tmpflags_libxc
    unset tmpflags_options
    unset tmpflags_bigdft
  else
    AC_MSG_NOTICE([CFGFLAGS_BIGDFT set => skipping BigDFT config tricks])
  fi

  dnl CPP tricks
  if test "${afb_bigdft_cppflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying BigDFT tricks (vendor: $1, version: $2, flags: C preprocessing)])

    CPPFLAGS_BIGDFT="${CPPFLAGS_BIGDFT} \$(afb_libxc_incs)"

    dnl Finish
    tmp_bigdft_cnt_tricks=`expr ${tmp_bigdft_cnt_tricks} \+ 1`
    afb_bigdft_tricky_vars="${afb_bigdft_tricky_vars} CPPFLAGS"
  else
    AC_MSG_NOTICE([CPPFLAGS_BIGDFT set => skipping BigDFT C preprocessing tricks])
  fi

  dnl Fortran tricks
  if test "${afb_bigdft_fcflags_custom}" = "no"; then
    AC_MSG_NOTICE([applying BigDFT tricks (vendor: $1, version: $2, flags: Fortran)])

    FCFLAGS_BIGDFT="${CPPFLAGS_BIGDFT} ${FCFLAGS_BIGDFT}"

    dnl Finish
    tmp_bigdft_cnt_tricks=`expr ${tmp_bigdft_cnt_tricks} \+ 1`
    afb_bigdft_tricky_vars="${afb_bigdft_tricky_vars} FCFLAGS"
  else
    AC_MSG_NOTICE([FCFLAGS_BIGDFT set => skipping BigDFT Fortran tricks])
  fi

  dnl Count applied tricks
  case "${tmp_bigdft_cnt_tricks}" in
    0)
      afb_bigdft_tricks="no"
      ;;
    ${tmp_bigdft_num_tricks})
      afb_bigdft_tricks="yes"
      ;;
    *)
      afb_bigdft_tricks="partial"
      ;;
  esac
  unset tmp_bigdft_cnt_tricks
  unset tmp_bigdft_num_tricks
]) # AFB_TRICKS_BIGDFT

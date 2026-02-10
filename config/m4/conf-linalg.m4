# -*- Autoconf -*-
#
# Copyright (C) 2005-2026 ABINIT Group (Marc Torrent)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Linear Algebra support for ABINIT
#


# ABI_LINALG_INIT()
# --------------
#
#
# Sets Linear Algebra parameters according to configure options.
#
#
AC_DEFUN([ABI_LINALG_INIT], [
  # Delegate most of the init stage to Steredeg
  SD_LINALG_INIT([required warn])
  SD_LINALG_INIT_FLAVOR

]) # ABI_LINALG_INIT


                    # ------------------------------------ #


# ABI_LINALG_DETECT()
# ----------------
#
# Tries first to determine whether the LinAlg implementation is usable,
# then takes appropriate actions.
#
AC_DEFUN([ABI_LINALG_DETECT], [
  # Delegate the actual detection to Steredeg
  SD_LINALG_DETECT
  abi_linalg_ok="${sd_linalg_enable}"

  if test "${abi_linalg_ok}" = "yes"; then

	# Test if BLAS library has buggy dot/norm interfaces
	if test "${abi_zdot_bugfix_enable}" = "no" -o "${abi_zdot_bugfix_enable}" = "auto"; then
	  abi_zdot_bugfix="${sd_linalg_has_buggy_zdot}"	
	  if test "${abi_zdot_bugfix}" = "yes" -a "${abi_zdot_bugfix_enable}" = "no"; then
		AC_MSG_ERROR([--enable-zdot-bugfix option is deactivated but the BLAS library has buggy interfaces!])
	  fi
	  if test "${abi_zdot_bugfix_enable}" = "auto"; then
		abi_zdot_bugfix_enable="${abi_zdot_bugfix}"
	  fi
	fi

	if test "${abi_zdot_bugfix_enable}" = "yes"; then
      AC_DEFINE([HAVE_LINALG_ZDOTC_BUG], 1,
        [Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.])
      AC_DEFINE([HAVE_LINALG_ZDOTU_BUG], 1,
        [Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.])
    fi

  fi # abi_linalg_ok

  AC_SUBST(abi_linalg_ok)

]) # ABI_LINALG_DETECT

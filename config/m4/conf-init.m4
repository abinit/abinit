# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the Abinit software package. For license information,
# please see the COPYING file in the top-level directory of the Abinit source
# distribution.
#

#
# Initialization
#



# ABI_INIT_ARCH()
# ---------------
#
# Initializes architecture-specific parameters.
#
AC_DEFUN([ABI_INIT_ARCH],[
  dnl Shared object extensions
  case "${target}" in
    *apple*)
      abi_so_ext="dylib"
      ;;
    *)
      abi_so_ext="so"
      ;;
  esac
]) # ABI_INIT_ARCH



# ABI_INIT_CPU_INFO()
# -------------------
#
# Sets architecture-related variables from the information given by the
# specified target. This is a helper for many other Abinit macros, that
# should be called quite early in the configure script.
#
# At present, the variables set are:
#
#  * abi_cpu_model  : CPU model, if guessed;
#  * abi_cpu_64bits : whether the CPU is 64 bits or not.
#
AC_DEFUN([ABI_INIT_CPU_INFO],[
  abi_cpu_platform=`echo "${target}" | cut -d- -f2`
  abi_cpu_vendor=""
  abi_cpu_model=""
  abi_cpu_spec=""
  abi_cpu_bits=""
  abi_cpu_64bits=""

  case "${target}" in

    alpha*)
      abi_cpu_vendor="dec"
      abi_cpu_model="${target_cpu}"
      abi_cpu_64bits=`echo "${abi_cpu_model}" | grep '64$'`
      if test "${abi_cpu_64bits}" = ""; then
        abi_cpu_64bits="no"
        abi_cpu_bits="32"
      else
        abi_cpu_64bits="yes"
        abi_cpu_bits="64"
      fi
      ;;

    powerpc*)
      abi_cpu_vendor="ibm"
      abi_cpu_model="${target_cpu}"
      abi_cpu_64bits=`echo "${abi_cpu_model}" | grep '64$'`
      if test "${abi_cpu_64bits}" = ""; then
        abi_cpu_64bits="no"
        abi_cpu_bits="32"
      else
        abi_cpu_64bits="yes"
        abi_cpu_bits="64"
      fi
      ;;

    i686-*linux*)
      dnl Athlon ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Athlon'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="amd"
          abi_cpu_model="athlon"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium 3 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Pentium III'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="pentium3"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium 4 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) 4'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="pentium4"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium 4M ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Pentium(R) M'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="pentium4"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Centrino ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2400'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="centrino"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium CoreDuo ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) CPU           T2050'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="coreduo"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium Core2 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Core(TM)2 CPU'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="core2"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Pentium Core2 Quad ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Intel(R) Core(TM)2 Quad CPU'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="core2"
          abi_cpu_64bits="no"
          abi_cpu_bits="32"
        fi
      fi
      dnl Unknown
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_vendor="unknown"
        abi_cpu_model="unknown"
        abi_cpu_64bits="unknown"
        abi_cpu_bits="32"
      fi
      ;;

    ia64-*linux*)
      dnl Itanium 1 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Itanium' | grep -v 'Itanium 2'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="itanium1"
        fi
      fi
      dnl Itanium 2 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Itanium 2'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="itanium2"
        fi
      fi
      dnl Madison ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat /proc/cpuinfo | grep 'Madison'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="itanium2"
        fi
      fi
      dnl Unknown
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_vendor="unknown"
        abi_cpu_model="unknown"
      fi
      dnl The processor is anyway 64-bit
      abi_cpu_64bits="yes"
      abi_cpu_bits="64"
      ;;

    mips*irix*)
      # Get processor type
      abi_cpu_vendor="mips"
      abi_cpu_model=`hinv 2> /dev/null | grep '^CPU: MIPS '`
      if test "${abi_cpu_model}" != ""; then
        abi_cpu_model=`echo "${abi_cpu_model}" | awk '{print tolower($3)}'`
      fi
      abi_cpu_64bits="yes"
      abi_cpu_bits="64"
      ;;

    x86_64-*linux*|x86_64*apple*)
      case "${target}" in
      x86_64-*linux*)
         cp /proc/cpuinfo cpuinfo
         ;;
      x86_64-*apple*)
         sysctl -A | grep 'machdep.cpu.brand_string' > cpuinfo
         ;;
      esac
      dnl Athlon64 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Athlon'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="amd"
          abi_cpu_model="athlon64"
        fi
      fi
      dnl Opteron ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Opteron'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="amd"
          abi_cpu_model="opteron"
        fi
      fi
      dnl Sempron ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Sempron'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="amd"
          abi_cpu_model="athlon64"
        fi
      fi
      dnl Pentium Core2 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Intel(R) Core(TM)2'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="core2"
        fi
      fi
      dnl Pentium Core2 Quad ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Intel(R) Core(TM)2 Quad'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="core2"
        fi
      fi
      dnl Pentium Core i3 ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Intel(R) Core(TM) i3'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="core_i3"
        fi
      fi
      dnl Xeon ?
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Intel(R) XEON(TM)'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="xeon"
        fi
      fi
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_model=`cat cpuinfo | grep 'Intel(R) Xeon(R)'`
        if test "${abi_cpu_model}" != ""; then
          abi_cpu_vendor="intel"
          abi_cpu_model="xeon"
        fi
      fi
      dnl Unknown
      if test "${abi_cpu_model}" = ""; then
        abi_cpu_vendor="unknown"
        abi_cpu_model="unknown"
      fi
      dnl The processor is anyway 64-bit
      abi_cpu_64bits="yes"
      abi_cpu_bits="64"
      rm -rf cpuinfo
      ;;

  esac

  dnl Generate CPU identifier
  abi_cpu_spec="${abi_cpu_vendor}_${abi_cpu_model}"

  dnl General system identifier
  abi_sys_spec=`echo "${target_os}" | cut -d- -f1 | sed -e 's/[[0-9]].*//'`
  abi_sys_spec="${abi_sys_spec}-${target_cpu}"

  AC_SUBST(abi_cpu_platform)
  AC_SUBST(abi_cpu_vendor)
  AC_SUBST(abi_cpu_model)
  AC_SUBST(abi_cpu_spec)
  AC_SUBST(abi_cpu_64bits)
  AC_SUBST(abi_cpu_bits)
  AC_SUBST(abi_sys_spec)
]) # ABI_INIT_CPU_INFO



# ABI_INIT_OS_INFO()
# ------------------
#
# Sets OS-related variables from the information given by the specified
# target.
#
AC_DEFUN([ABI_INIT_OS_INFO],[
  case "${target_os}" in

    *linux*)
      AC_DEFINE([HAVE_OS_LINUX],1,[Define to 1 if you are using Linux.])
      ;;

    *apple*)
      AC_DEFINE([HAVE_OS_MACOSX],1,[Define to 1 if you are using MacOS X.])
      ;;

    *cygwin*|*mingw*)
      AC_DEFINE([HAVE_OS_WINDOWS],1,[Define to 1 if you are using Windows.])
      ;;

  esac
]) # ABI_INIT_OS_INFO



# ABI_INIT_HEADER()
# -----------------
#
# Initializes the contents of the header file produced by Autoheader.
#
AC_DEFUN([ABI_INIT_HEADER],[
  dnl Set top of file ...
  AH_TOP([/*
 * Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
 *
 * This file is part of the Abinit software package. For license information,
 * please see the COPYING file in the top-level directory of the Abinit source
 * distribution.
 *
 */

/* Abinit configuration */

#ifndef _ABINIT_CONFIG_H
#define _ABINIT_CONFIG_H

#ifdef __INTEL_COMPILER
#define FC_INTEL 1
#endif

])

  dnl ... as well as bottom
  AH_BOTTOM([/* *** BEGIN sanity checks *** */

/* MPI options */
#if defined HAVE_MPI 

/* Check that one MPI level is actually defined */
#if ! defined HAVE_MPI1 && ! defined HAVE_MPI2 && ! defined HAVE_MPI3
#error "HAVE_MPI1, HAVE_MPI2, and HAVE_MPI3, are all undefined"
#endif

/* Check that only one MPI level has been defined */
#if defined HAVE_MPI1
#  if defined HAVE_MPI2
#    if defined HAVE_MPI3
#      error "HAVE_MPI1, Have_MPI2, and HAVE_MPI3, are all defined"
#    else
#error "HAVE_MPI1 and HAVE_MPI2 are both defined"
#    endif
#  else
#    if defined HAVE_MPI3
#      error "HAVE_MPI1 and HAVE_MPI3 are both defined"
#    endif
#  endif
#else
#  if defined HAVE_MPI2 && defined HAVE_MPI3
#    error "HAVE_MPI2 and HAVE_MPI3 are both defined"
#  endif
#endif

#else /* HAVE_MPI */

/* Check that no MPI level is defined */
#if defined HAVE_MPI1 || defined HAVE_MPI2 || HAVE_MPI3
#error "HAVE_MPI1, HAVE_MPI2, and HAVE_MPI3, must be undefined"
#endif

/* Check that MPI-IO is undefined */
#if defined HAVE_MPI_IO
#error "HAVE_MPI_IO must be undefined"
#endif

#endif /* HAVE_MPI */

/* *** END sanity checks *** */

#endif /* _ABINIT_CONFIG_H */])
]) # ABI_INIT_HEADER



# ABI_INIT_INSTALL_DIRS()
# -----------------------
#
# Sets installation directories.
#
AC_DEFUN([ABI_INIT_INSTALL_DIRS],[
  dnl Set-up prefix
  if test "${prefix}" = "NONE"; then
    abinit_prefix="${ac_default_prefix}"
  else
    abinit_prefix="${prefix}"
  fi

  dnl Set-up all directory names
  abinit_bindir="${abinit_prefix}/bin"
  abinit_chkdir="${abinit_prefix}/share/abinit/tests"
  abinit_datdir="${abinit_prefix}/share/abinit"
  abinit_docdir="${abinit_prefix}/doc/abinit"
  abinit_incdir="${abinit_prefix}/include"
  abinit_libdir="${abinit_prefix}/lib"
  abinit_mandir="${abinit_prefix}/share/man"

  dnl Substitute all variables
  AC_SUBST(abinit_prefix)
  AC_SUBST(abinit_bindir)
  AC_SUBST(abinit_chkdir)
  AC_SUBST(abinit_datdir)
  AC_SUBST(abinit_docdir)
  AC_SUBST(abinit_incdir)
  AC_SUBST(abinit_libdir)
  AC_SUBST(abinit_mandir)
]) # ABI_INIT_INSTALL_DIRS



# ABI_INIT_TARGET()
# -----------------
#
# Initializes the target name for the platform Abinit is about to be built on.
#
# Note: to be called after the detection of the Fortran compiler type.
#
AC_DEFUN([ABI_INIT_TARGET],[
  dnl Clean-up operating system name
  [abi_target_os=`echo ${target_os} | sed -e 's/-.*//'`]
  
  ABINIT_TARGET="${target_cpu}_${abi_target_os}_${abi_fc_vendor}${abi_fc_version}"
  AC_DEFINE_UNQUOTED(ABINIT_TARGET,"${ABINIT_TARGET}",
    [Abinit target description.])
  AC_SUBST(ABINIT_TARGET)
]) # ABI_INIT_TARGET



# ABI_INIT_VERSION()
# ------------------
#
# Sets all variables related to the current version of Abinit.
#
AC_DEFUN([ABI_INIT_VERSION],[
  dnl Get version from Autoconf
  ABINIT_VERSION="${PACKAGE_VERSION}"
  ABINIT_VERSION_MAJOR=`echo "${ABINIT_VERSION}" | cut -d. -s -f1`
  ABINIT_VERSION_MINOR=`echo "${ABINIT_VERSION}" | cut -d. -s -f2`
  ABINIT_VERSION_MICRO=`echo "${ABINIT_VERSION}" | cut -d. -s -f3`
  ABINIT_VERSION_MINOR=`echo "${ABINIT_VERSION_MINOR}" | sed -e 's/[a-z]//g'`
  if test "${ABINIT_VERSION_MICRO}" = ""; then
    ABINIT_VERSION_MICRO=`echo "${ABINIT_VERSION}" | cut -b4-`
  fi
  if test "${ABINIT_VERSION_MICRO}" = ""; then
    ABINIT_VERSION_MICRO="dev"
  fi
  ABINIT_VERSION_BUILD=`date '+%Y%m%d'`

  ABINIT_VERSION_BASE="${ABINIT_VERSION_MAJOR}.${ABINIT_VERSION_MINOR}"

  dnl Make numbers available to source files
  AC_DEFINE_UNQUOTED(ABINIT_VERSION,"${ABINIT_VERSION}",
    [Abinit whole version number.])
  AC_DEFINE_UNQUOTED(ABINIT_VERSION_MAJOR,"${ABINIT_VERSION_MAJOR}",
    [Abinit major version number.])
  AC_DEFINE_UNQUOTED(ABINIT_VERSION_MINOR,"${ABINIT_VERSION_MINOR}",
    [Abinit minor version number.])
  AC_DEFINE_UNQUOTED(ABINIT_VERSION_MICRO,"${ABINIT_VERSION_MICRO}",
    [Abinit micro version number (patch level).])
  AC_DEFINE_UNQUOTED(ABINIT_VERSION_BUILD,"${ABINIT_VERSION_BUILD}",
    [Abinit build date.])
  AC_DEFINE_UNQUOTED(ABINIT_VERSION_BASE,"${ABINIT_VERSION_BASE}",
    [Abinit base version number.])
  AC_SUBST(ABINIT_VERSION)
  AC_SUBST(ABINIT_VERSION_MAJOR)
  AC_SUBST(ABINIT_VERSION_MINOR)
  AC_SUBST(ABINIT_VERSION_MICRO)
  AC_SUBST(ABINIT_VERSION_BUILD)
  AC_SUBST(ABINIT_VERSION_BASE)
]) # ABI_INIT_VERSION

# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# C compilers support
#



# _ABI_CHECK_CC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the GNU C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
# Note: This macro should be called after AC_PROG_CC.
#
AC_DEFUN([_ABI_CHECK_CC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU C compiler])
  cc_info_string=`$1 --version 2>&1 | ${EGREP} '^g?cc' | head -n 1`
  if test "${ac_cv_c_compiler_gnu}" != "yes"; then
    abi_cc_vendor="unknown"
    abi_cc_version="unknown"
    abi_result="no"
  else
    AC_DEFINE([CC_GNU],1,[Define to 1 if you are using the GNU C compiler.])
    abi_cc_vendor="gnu"
    abi_cc_version=`echo ${cc_info_string} | sed -e 's/.*([[^)]]*) //; s/ .*//'`
    if test "${abi_cc_version}" = "${cc_info_string}"; then
      abi_result=`echo "${cc_info_string}" | grep ' '`
      if test "${abi_result}" != ""; then
        abi_cc_version="unknown"
      fi
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_GNU



# _ABI_CHECK_CC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the IBM XL C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL C compiler])
  cc_info_string=`$1 -qversion 2>&1 | head -n 1`
  cc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  abi_result=`echo "${cc_info_string}" | grep 'IBM XL C/C++'`
  if test "${abi_result}" = ""; then
    abi_result=`echo "${cc_info_string}" | grep 'IBM(R) XL C/C++'`
  fi
  if test "${abi_result}" = ""; then
    abi_result=`echo "${cc_info_string}" | grep 'C for AIX'`
  fi
  if test "${abi_result}" = ""; then
    abi_result="no"
    cc_info_string=""
    abi_cc_vendor="unknown"
    abi_cc_version="unknown"
    if test "${cc_garbage}" -gt 50; then
      AC_DEFINE([CC_IBM],1,[Define to 1 if you are using the IBM XL C compiler.])
      abi_cc_vendor="ibm"
      abi_cc_version="unknown"
      abi_result="yes"
    fi
  else
    AC_DEFINE([CC_IBM],1,[Define to 1 if you are using the IBM XL C compiler.])
    abi_cc_vendor="ibm"
    abi_cc_version=`echo "${abi_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${abi_cc_version}" = "${abi_result}"; then
      abi_cc_version=`echo "${abi_result}" | sed -e 's/C for AIX version //'`
    fi
    if test "${abi_cc_version}" = "${abi_result}"; then
      abi_cc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_IBM



# _ABI_CHECK_CC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the Intel C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel C compiler])
  cc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${cc_info_string}" | grep '^Intel(R) C'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cc_info_string=""
    abi_cc_vendor="unknown"
    abi_cc_version="unknown"
  else
    AC_DEFINE([CC_INTEL],1,[Define to 1 if you are using the Intel C compiler.])
    abi_cc_vendor="intel"
    abi_cc_version=`echo "${abi_result}" | sed -e 's/.*Version //; s/ .*//'`
    if test "${abi_cc_version}" = "${abi_result}"; then
      abi_cc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_INTEL



# _ABI_CHECK_CC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Portland Group C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PGI C compiler])
  cc_info_string=`$1 -V 2>&1 | sed -e '/^$/d' | head -n 1`
  abi_result=`echo "${cc_info_string}" | grep '^pgcc'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cc_info_string=""
    abi_cc_vendor="unknown"
    abi_cc_version="unknown"
  else
    AC_DEFINE([CC_PGI],1,[Define to 1 if you are using the Portland Group C compiler.])
    abi_cc_vendor="pgi"
    abi_cc_version=`echo "${abi_result}" | sed -e 's/.* //; s/-.*//'`
    if test "${abi_cc_version}" = "${abi_result}"; then
      abi_cc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CC_PGI



 #############################################################################



# _ABI_CHECK_CC_HEADERS()
# -----------------------
#
# Checks for header files.
#
AC_DEFUN([_ABI_CHECK_CC_HEADERS],[
  dnl Init AC_MSG_CHECKING([for C header files])
  
  dnl The following line causes trouble to aclocal
  dnl AC_HEADER_STDC
  AC_CHECK_HEADERS([stddef.h stdarg.h])
  dnl AC_CHECK_HEADERS([stdlib.h])
  AC_CHECK_HEADERS([stdio.h math.h termios.h])
  AC_CHECK_HEADERS([errno.h])
  AC_CHECK_HEADERS([malloc.h sys/malloc.h])
  AC_CHECK_HEADERS([mcheck.h])
  AC_CHECK_HEADERS([sys/time.h])
  AC_CHECK_HEADERS([sys/resource.h])
  dnl AC_CHECK_HEADERS([sys/ioctl.h sys/sysctl.h])
  dnl AC_CHECK_HEADERS([sys/stat.h])
  dnl AC_CHECK_HEADERS([string.h])
  dnl AC_CHECK_HEADERS([strings.h])
  dnl AC_CHECK_HEADERS([unistd.h])
  dnl AC_CHECK_HEADERS([limits.h])

]) # _ABI_CHECK_CC_HEADERS



# _ABI_CHECK_CC_FUNCTIONS()
# -------------------------
#
# Checks for library functions.
#
AC_DEFUN([_ABI_CHECK_CC_FUNCTIONS],[
  dnl Init AC_MSG_CHECKING([for library functions])

  dnl AC_CHECK_FUNCS([BSDgettimeofday gettimeofday gethrtime]) 
  AC_CHECK_FUNCS([abort])
  AC_CHECK_FUNCS([mallinfo])

]) # _ABI_CHECK_CC_FUNCTIONS



# _ABI_CHECK_CC_FEATURES()
# ------------------------
#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_DEFUN([_ABI_CHECK_CC_FEATURES],[
  dnl Init AC_MSG_CHECKING([for C compiler characteristics])

  AC_CHECK_SIZEOF(char)
  AC_CHECK_SIZEOF(short)
  AC_CHECK_SIZEOF(int)
  AC_CHECK_SIZEOF(long)
  AC_CHECK_SIZEOF(long long)
  AC_CHECK_SIZEOF(unsigned int)
  AC_CHECK_SIZEOF(unsigned long)
  AC_CHECK_SIZEOF(unsigned long long)
  AC_CHECK_SIZEOF(float)
  AC_CHECK_SIZEOF(double)
  AC_CHECK_SIZEOF(long double)
  AC_CHECK_SIZEOF(size_t)
  AC_CHECK_SIZEOF(ptrdiff_t)

  AC_C_CONST
  AC_TYPE_SIZE_T
  dnl AC_TYPE_PID_T

]) # _ABI_CHECK_CC_FEATURES



 #############################################################################



# ABI_CC_FEATURES()
# -----------------
#
# Explores the capabilities of the C compiler.
#
AC_DEFUN([ABI_CC_FEATURES],[
  dnl Explore compiler peculiarities
  _ABI_CHECK_CC_HEADERS
  _ABI_CHECK_CC_FUNCTIONS
  _ABI_CHECK_CC_FEATURES
]) # ABI_CC_FEATURES



# ABI_PROG_CC()
# -------------
#
# Tries to determine which type of C compiler is installed.
#
AC_DEFUN([ABI_PROG_CC],[
  dnl Init
  if test "${abi_cc_vendor}" = ""; then
    abi_cc_vendor="unknown"
  fi

  dnl Determine C compiler type (the order is important)
  AC_MSG_CHECKING([which type of compiler we have])

  dnl Always get rid of that one as early as possible
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CHECK_CC_IBM(${CC})
  fi

  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CHECK_CC_INTEL(${CC})
  fi
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CHECK_CC_PGI(${CC})
  fi

  dnl Check the GNU compiler last, because other compilers are cloning
  dnl its CLI
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CHECK_CC_GNU(${CC})
  fi

  dnl Fall back to generic when detection fails
  if test "${abi_cc_vendor}" = "unknown"; then
    abi_cc_vendor="generic"
    abi_cc_version="0.0"
  fi

  dnl Normalize C compiler version
  abi_cc_version=`echo ${abi_cc_version} | cut -d. -f1-2`

  dnl Display final result
  AC_MSG_RESULT([${abi_cc_vendor} ${abi_cc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(abi_cc_vendor)
  AC_SUBST(abi_cc_version)
  AC_SUBST(cc_info_string)
]) # ABI_PROG_CC

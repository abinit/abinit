# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# C++ compilers support
#



# _ABI_CHECK_CXX_ARM(COMPILER)
# ------------------------------
#
# Checks whether the specified C++ compiler is the ARMClang++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_ARM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ARM C++ compiler])
  cxx_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${cxx_info_string}" | grep '^Arm C/C++/Fortran Compiler'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_ARM],1,
      [Define to 1 if you are using the ARM C++ compiler.])
    abi_cxx_vendor="arm"
    abi_cxx_version=`echo ${abi_result} | sed -e 's/.*ersion //; s/ .*//'`
    if test "${abi_cxx_version}" = "${abi_result}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_ARM



# _ABI_CHECK_CXX_GNU(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the GNU C++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU C++ compiler])
  cxx_info_string=`$1 --version 2>&1 | head -n 1`
  if test "${ac_cv_cxx_compiler_gnu}" != "yes"; then
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
    abi_result="no"
  else
    AC_DEFINE([CXX_GNU],1,[Define to 1 if you are using the GNU C++ compiler.])
    abi_cxx_vendor="gnu"
    abi_cxx_version=`echo ${cxx_info_string} | sed -e 's/.*([[^)]]*) //; s/ .*//'`
    if test "${abi_cxx_version}" = "${cxx_info_string}"; then
      abi_result=`echo "${cxx_info_string}" | grep ' '`
      if test "${abi_result}" != ""; then
        abi_cxx_version="unknown"
      fi
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_GNU



# _ABI_CHECK_CXX_IBM(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the IBM XL C++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL C++ compiler])
  cxx_info_string=`$1 -qversion 2>&1 | head -n 1`
  cxx_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  abi_result=`echo "${cc_info_string}" | grep 'IBM XL C/C++'`
  if test "${abi_result}" = ""; then
    abi_result=`echo "${cxx_info_string}" | grep 'IBM(R) XL C/C++'`
  fi
  if test "${abi_result}" = ""; then
    abi_result=`echo "${cxx_info_string}" | grep 'C for AIX'`
  fi
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
    if test "${cxx_garbage}" -gt 50; then
      AC_DEFINE([CXX_IBM],1,[Define to 1 if you are using the IBM XL C++ compiler.])
      abi_cxx_vendor="ibm"
      abi_cxx_version="unknown"
      abi_result="yes"
    fi
  else
    AC_DEFINE([CXX_IBM],1,[Define to 1 if you are using the IBM XL C++ compiler.])
    abi_cxx_vendor="ibm"
    abi_cxx_version=`echo "${cxx_info_string}" | sed -e 's/.* V//; s/ .*//'`
    if test "${abi_cxx_version}" = "${cxx_info_string}"; then
      abi_cxx_version=`echo "${cxx_info_string}" | sed -e 's/C for AIX version //'`
    fi
    if test "${abi_cxx_version}" = "${cxx_info_string}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_IBM



# _ABI_CHECK_CXX_INTEL(COMPILER)
# ------------------------------
#
# Checks whether the specified C++ compiler is the Intel C++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel C++ compiler])
  cxx_info_string=`$1 -v -V 2>&1 | head -n 1`
  abi_result=`echo "${cxx_info_string}" | grep '^Intel(R) C++'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_INTEL],1,[Define to 1 if you are using the Intel C++ compiler.])
    abi_cxx_vendor="intel"
    abi_cxx_version=`echo "${abi_result}" | sed -e 's/.*Version //; s/ .*//'`
    if test "${abi_cxx_version}" = "${abi_result}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_INTEL


# _ABI_CHECK_CXX_LLVM(COMPILER)
# ------------------------------
#
# Checks whether the specified C compiler is the LLVM Clang++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_LLVM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the LLVM Clang++ C++ compiler])
  cxx_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${cxx_info_string}" | grep '[[Cc]]lang'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_LLVM],1,
      [Define to 1 if you are using the LLVM Clang++ C++ compiler.])
    abi_cxx_vendor="llvm"
    abi_cxx_version=`echo ${abi_result} | sed -e 's/.*ersion //; s/ .*//'`
    if test "${abi_cxx_version}" = "${abi_result}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_LLVM



# _ABI_CHECK_CXX_PGI(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the Portland Group C++
# compiler. If yes, tries to determine its version number and sets the
# abi_cxx_vendor and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_CXX_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Portland Group C++ compiler])
  cxx_info_string=`$1 -v -V 2>&1 | sed -e '/^$/d' | head -n 1`
  abi_result=`echo "${cxx_info_string}" | grep '^pgCC'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_PGI],1,[Define to 1 if you are using the Portland Group C++ compiler.])
    abi_cxx_vendor="pgi"
    abi_cxx_version=`echo "${abi_result}" | sed -e 's/.* //; s/-.*//'`
    if test "${abi_cxx_version}" = "${abi_result}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_CXX_PGI



# ABI_PROG_CXX()
# --------------
#
# Tries to determine which type of C++ compiler is installed.
#
AC_DEFUN([ABI_PROG_CXX],[
  dnl Init
  if test "${abi_cxx_vendor}" = ""; then
    abi_cxx_vendor="unknown"
  fi

  dnl Determine C++ compiler type (the order is important)
  AC_MSG_CHECKING([which type of C++ compiler we have])

  dnl Get rid of that one as early as possible
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_IBM(${CXX})
  fi

  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_ARM(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_INTEL(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_LLVM(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_PGI(${CXX})
  fi

  dnl Check the GNU compiler last, because other compilers are cloning
  dnl its CLI
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CHECK_CXX_GNU(${CXX})
  fi

  dnl Fall back to generic when detection fails
  if test "${abi_cxx_vendor}" = "unknown"; then
    abi_cxx_vendor="generic"
    abi_cxx_version="0.0"
  fi

  dnl Normalize C++ compiler version
  abi_cxx_version=`echo ${abi_cxx_version} | cut -d. -f1-2`

  dnl Display final result
  AC_MSG_RESULT([${abi_cxx_vendor} ${abi_cxx_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(abi_cxx_vendor)
  AC_SUBST(abi_cxx_version)
  AC_SUBST(cxx_info_string)
]) # ABI_PROG_CXX

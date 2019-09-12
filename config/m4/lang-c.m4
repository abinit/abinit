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



# _ABI_CC_CHECK_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the GNU C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
# Note: This macro should be called after AC_PROG_CC.
#
AC_DEFUN([_ABI_CC_CHECK_GNU],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  # AC_MSG_CHECKING([if we are using the GNU C compiler])
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
  # AC_MSG_RESULT(${abi_result})
]) # _ABI_CC_CHECK_GNU



# _ABI_CC_CHECK_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the IBM XL C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CC_CHECK_IBM],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  # AC_MSG_CHECKING([if we are using the IBM XL C compiler])
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
  # AC_MSG_RESULT(${abi_result})
]) # _ABI_CC_CHECK_IBM



# _ABI_CC_CHECK_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the Intel C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CC_CHECK_INTEL],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  # AC_MSG_CHECKING([if we are using the Intel C compiler])
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
  # AC_MSG_RESULT(${abi_result})
]) # _ABI_CC_CHECK_INTEL



# _ABI_CC_CHECK_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Portland Group C compiler.
# If yes, tries to determine its version number and sets the abi_cc_vendor
# and abi_cc_version variables accordingly.
#
AC_DEFUN([_ABI_CC_CHECK_PGI],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  # AC_MSG_CHECKING([if we are using the PGI C compiler])
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
  # AC_MSG_RESULT(${abi_result})
]) # _ABI_CC_CHECK_PGI



 #############################################################################



# _ABI_CC_CHECK_HEADERS()
# -----------------------
#
# Checks for header files.
#
AC_DEFUN([_ABI_CC_CHECK_HEADERS],[
  # Look for standard headers
  AC_HEADER_ASSERT
  AC_CHECK_HEADERS([stdio.h string.h termios.h unistd.h])
  AC_CHECK_HEADERS([errno.h])
  AC_CHECK_HEADERS([inttypes.h math.h stddef.h stdint.h])
  AC_CHECK_HEADERS([mcheck.h time.h])
  AC_CHECK_HEADERS([sys/ioctl.h sys/resource.h sys/stat.h sys/time.h sys/types.h])

  # Look for malloc.h
  CPPFLAGS_MALLOC=""
  AC_CHECK_HEADERS([sys/malloc.h])
  AC_CHECK_HEADERS([malloc.h], [abi_hdr_malloc="yes"], [abi_hdr_malloc="no"])
  if test "${abi_hdr_malloc}" = "no"; then
    AC_CHECK_HEADERS([malloc/malloc.h],
      [abi_hdr_malloc="yes"], [abi_hdr_malloc="no"])
  fi
]) # _ABI_CC_CHECK_HEADERS



# _ABI_CC_CHECK_FUNCTIONS()
# -------------------------
#
# Checks for library functions.
#
AC_DEFUN([_ABI_CC_CHECK_FUNCTIONS],[
  # Init AC_MSG_CHECKING([for library functions])

  # AC_CHECK_FUNCS([BSDgettimeofday gettimeofday gethrtime]) 
  AC_CHECK_FUNCS([abort])
  AC_CHECK_FUNCS([mallinfo])

]) # _ABI_CC_CHECK_FUNCTIONS



# _ABI_CC_CHECK_TYPES()
# ---------------------
#
# Checks for typedefs, structures, and compiler characteristics.
#
AC_DEFUN([_ABI_CC_CHECK_TYPES],[
  # Init AC_MSG_CHECKING([for C compiler characteristics])

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
  AC_C_BIGENDIAN
  # AC_TYPE_PID_T

]) # _ABI_CC_CHECK_TYPES



# _ABI_CC_CHECK_XML()
# -------------------
#
# Checks for LibXML2, a basic requirement on any decent Unix-like system.
#
AC_DEFUN([_ABI_CC_CHECK_XML],[
  abi_libxml2_cppflags=""
  abi_libxml2_fcflags=""
  abi_libxml2_ldflags=""
  abi_libxml2_libs=""
  if test "${abi_libxml2_enable}" = "yes" -o "${abi_libxml2_enable}" = "auto"; then
    AC_LANG_PUSH([C])
    AM_PATH_XML2([2.7.6], [abi_libxml2_ok="yes"], [abi_libxml2_ok="no"])
    AC_LANG_POP([C])
    if test "${abi_libxml2_ok}" = "yes"; then
      abi_libxml2_enable="yes"
      abi_libxml2_cppflags="${XML_CPPFLAGS}"
      abi_libxml2_fcflags="${XML_CPPFLAGS}"
      abi_libxml2_libs="${XML_LIBS}"
    else
      AC_MSG_WARN([could not find a usable XML library => disabling XML support])
      abi_libxml2_enable="no"
    fi
  fi
  AC_SUBST(abi_libxml2_cppflags)
  AC_SUBST(abi_libxml2_fcflags)
  AC_SUBST(abi_libxml2_ldflags)
  AC_SUBST(abi_libxml2_libs)
]) # _ABI_CC_CHECK_XML



 #############################################################################



# ABI_CC_FEATURES()
# -----------------
#
# Explores the capabilities of the C compiler.
#
AC_DEFUN([ABI_CC_FEATURES],[
  # Explore compiler peculiarities
  _ABI_CC_CHECK_HEADERS
  _ABI_CC_CHECK_FUNCTIONS
  _ABI_CC_CHECK_TYPES
  _ABI_CC_CHECK_XML
]) # ABI_CC_FEATURES



# ABI_PROG_CC()
# -------------
#
# Tries to determine which type of C compiler is installed.
#
AC_DEFUN([ABI_PROG_CC],[
  # Init
  if test "${abi_cc_vendor}" = ""; then
    abi_cc_vendor="unknown"
  fi

  # Preserve environment
  ABI_ENV_BACKUP

  # Look for the C compiler
  if test "${CC}" != "" -a ! -x "${CC}"; then
    abi_cc_probe=`echo "${CC}" | sed -e 's/ .*//'`
    if test ! -x "${abi_cc_probe}"; then
      AC_PATH_PROG([abi_cc_path], [${abi_cc_probe}])
      if test "${abi_cc_path}" = ""; then
        AC_MSG_ERROR([could not run C compiler "${CC}"])
      fi
    fi
  fi
  AC_PROG_CC

  # Fail if no C compiler is available
  if test "${CC}" = ""; then
    AC_MSG_ERROR([no C compiler available])
  fi

  # Look for the C preprocessor
  if test "${CPP}" != "" -a ! -x "${CPP}"; then
    AC_PATH_PROG([abi_cpp_path], [${CPP}])
    if test "${abi_cpp_path}" = ""; then
      AC_MSG_ERROR([could not run C preprocessor "${CPP}"])
    fi
  fi
  AC_PROG_CPP

  # Fail if no C preprocessor is available
  if test "${CPP}" = ""; then
    AC_MSG_ERROR([no C preprocessor available])
  fi

  # Determine C compiler vendor (the order is critical)
  AC_MSG_CHECKING([which type of compiler we have])
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CC_CHECK_IBM(${CC})
  fi
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CC_CHECK_INTEL(${CC})
  fi
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CC_CHECK_PGI(${CC})
  fi

  # Check the GNU compiler last, because other compilers are cloning
  # its CLI
  if test "${abi_cc_vendor}" = "unknown"; then
    _ABI_CC_CHECK_GNU(${CC})
  fi

  # Fall back to generic when detection fails
  if test "${abi_cc_vendor}" = "unknown"; then
    abi_cc_vendor="generic"
    abi_cc_version="0.0"
  fi

  # Normalize C compiler version
  abi_cc_version=`echo ${abi_cc_version} | cut -d. -f1-2`

  # Display final result
  AC_MSG_RESULT([${abi_cc_vendor} ${abi_cc_version}])

  # Restore back CPPFLAGS and CFLAGS
  CPPFLAGS="${abi_env_CPPFLAGS}"
  CFLAGS="${abi_env_CFLAGS}"

  # Schedule compiler info for substitution
  AC_SUBST(abi_cc_vendor)
  AC_SUBST(abi_cc_version)
  AC_SUBST(cc_info_string)
]) # ABI_PROG_CC

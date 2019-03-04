# -*- Autoconf -*-
#
# Copyright (C) 2005-2016 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# C compilers support
#



# _AFB_CHECK_CC_COMPAQ(COMPILER)
# ------------------------------
#
# Checks whether the specified C compiler is the COMPAQ C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_COMPAQ],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Compaq C compiler])
  cc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep '^Compaq C '`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_COMPAQ],1,[Define to 1 if you are using the COMPAQ C compiler.])
    afb_cc_vendor="compaq"
    cc_info_string=`$1 -V 2>&1 | grep '^Compiler Driver' | head -n 1`
    afb_cc_version=`echo "${cc_info_string}" | sed -e 's/Compiler Driver V//; s/ .*//'`
    if test "${afb_cc_version}" = "${cc_info_string}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_COMPAQ



# _AFB_CHECK_CC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the GNU C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
# Note: This macro should be called after AC_PROG_CC.
#
AC_DEFUN([_AFB_CHECK_CC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU C compiler])
  cc_info_string=`$1 --version 2>&1 | ${EGREP} '^g?cc' | head -n 1`
  if test "${ac_cv_c_compiler_gnu}" != "yes"; then
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
    afb_result="no"
  else
    AC_DEFINE([CC_GNU],1,[Define to 1 if you are using the GNU C compiler.])
    afb_cc_vendor="gnu"
    afb_cc_version=`echo ${cc_info_string} | sed -e 's/.*([[^)]]*) //; s/ .*//'`
    if test "${afb_cc_version}" = "${cc_info_string}"; then
      afb_result=`echo "${cc_info_string}" | grep ' '`
      if test "${afb_result}" != ""; then
        afb_cc_version="unknown"
      fi
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_GNU



# _AFB_CHECK_CC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the IBM XL C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL C compiler])
  cc_info_string=`$1 -qversion 2>&1 | head -n 1`
  cc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  afb_result=`echo "${cc_info_string}" | grep 'IBM XL C/C++'`
  if test "${afb_result}" = ""; then
    afb_result=`echo "${cc_info_string}" | grep 'IBM(R) XL C/C++'`
  fi
  if test "${afb_result}" = ""; then
    afb_result=`echo "${cc_info_string}" | grep 'C for AIX'`
  fi
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
    if test "${cc_garbage}" -gt 50; then
      AC_DEFINE([CC_IBM],1,[Define to 1 if you are using the IBM XL C compiler.])
      afb_cc_vendor="ibm"
      afb_cc_version="unknown"
      afb_result="yes"
    fi
  else
    AC_DEFINE([CC_IBM],1,[Define to 1 if you are using the IBM XL C compiler.])
    afb_cc_vendor="ibm"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version=`echo "${afb_result}" | sed -e 's/C for AIX version //'`
    fi
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_IBM



# _AFB_CHECK_CC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the Intel C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel C compiler])
  cc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep '^Intel(R) C'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_INTEL],1,[Define to 1 if you are using the Intel C compiler.])
    afb_cc_vendor="intel"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.*Version //; s/ .*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_INTEL



# _AFB_CHECK_CC_PATHSCALE(COMPILER)
# ---------------------------------
#
# Checks whether the specified C compiler is the PathScale C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_PATHSCALE],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PathScale C compiler])
  cc_info_string=`$1 --version 2>&1 | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep '^PathScale'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_PATHSCALE],1,[Define to 1 if you are using the PathScale C compiler.])
    afb_cc_vendor="pathscale"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_PATHSCALE


# _AFB_CHECK_CC_OPEN64(COMPILER)
# ---------------------------------
#
# Checks whether the specified C compiler is the Open64 C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_OPEN64],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Open64 C compiler])
  cc_info_string=`$1 --version 2>&1 | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep '^Open64'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_OPEN64],1,[Define to 1 if you are using the Open64 C compiler.])
    afb_cc_vendor="open64"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_OPEN64

# _AFB_CHECK_CC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Portland Group C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PGI C compiler])
  cc_info_string=`$1 -V 2>&1 | sed -e '/^$/d' | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep '^pgcc'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_PGI],1,[Define to 1 if you are using the Portland Group C compiler.])
    afb_cc_vendor="pgi"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.* //; s/-.*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_PGI



# _AFB_CHECK_CC_SUN(COMPILER)
# ---------------------------
#
# Checks whether the specified C compiler is the Sun C compiler.
# If yes, tries to determine its version number and sets the afb_cc_vendor
# and afb_cc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_CC_SUN],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Sun C compiler])
  cc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${cc_info_string}" | grep 'Sun' | grep ' C '`
  if test "${afb_result}" = ""; then
    afb_result="no"
    cc_info_string=""
    afb_cc_vendor="unknown"
    afb_cc_version="unknown"
  else
    AC_DEFINE([CC_SUN],1,[Define to 1 if you are using the Sun C compiler.])
    afb_cc_vendor="sun"
    afb_cc_version=`echo "${afb_result}" | sed -e 's/.* C //; s/ .*//'`
    if test "${afb_cc_version}" = "${afb_result}"; then
      afb_cc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_CC_SUN



 #############################################################################



# AFB_PROG_CC(CC)
# ---------------
#
# Tries to determine which type of C compiler is installed.
#
AC_DEFUN([AFB_PROG_CC],[
  dnl Init
  if test "${afb_cc_vendor}" = ""; then
    afb_cc_vendor="unknown"
  fi

  dnl Determine C compiler type (the order is important)
  AC_MSG_CHECKING([for compiler type of '$1'])

  dnl Always get rid of that one as early as possible
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_IBM($1)
  fi

  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_COMPAQ($1)
  fi
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_INTEL($1)
  fi
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_PATHSCALE($1)
  fi
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_OPEN64($1)
  fi
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_PGI($1)
  fi
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_SUN($1)
  fi

  dnl Check the GNU compiler last, because other compilers are cloning
  dnl its CLI
  if test "${afb_cc_vendor}" = "unknown"; then
    _AFB_CHECK_CC_GNU($1)
  fi

  dnl Fall back to generic when detection fails
  if test "${afb_cc_vendor}" = "unknown"; then
    afb_cc_vendor="generic"
    afb_cc_version="0.0"
  fi

  dnl Normalize C compiler version
  afb_cc_version=`echo ${afb_cc_version} | cut -d. -f1-2`

  dnl Display final result
  AC_MSG_RESULT([${afb_cc_vendor} ${afb_cc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(afb_cc_vendor)
  AC_SUBST(afb_cc_version)
  AC_SUBST(cc_info_string)
]) # AFB_PROG_CC

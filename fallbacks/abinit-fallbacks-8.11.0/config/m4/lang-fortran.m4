# -*- Autoconf -*-
#
# Copyright (C) 2005-2014 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _AFB_CHECK_FC_ABSOFT(COMPILER)
# ----------------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_ABSOFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^Pro Fortran'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_ABSOFT],1,
      [Define to 1 if you are using the ABSOFT Fortran compiler.])
    afb_fc_vendor="absoft"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/Pro Fortran //'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_ABSOFT



# _AFB_CHECK_FC_COMPAQ(COMPILER)
# ----------------------------------------
#
# Checks whether the specified Fortran compiler is the COMPAQ Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_COMPAQ],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Compaq Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^Compaq Fortran Compiler'`
  afb_result="${fc_info_string}"
  if test "${afb_result}" = ""; then
    fc_info_string=`$1 -version 2>&1 | sed -e 's/^	//' | grep '^HP Fortran Compiler'`
    afb_result="${fc_info_string}"
  fi
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_COMPAQ],1,
      [Define to 1 if you are using the COMPAQ Fortran compiler.])
    afb_fc_vendor="compaq"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.* V//;s/-.*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_COMPAQ



# _AFB_CHECK_FC_FUJITSU(COMPILER)
# -----------------------------------------
#
# Checks whether the specified Fortran compiler is the Fujitsu Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_FUJITSU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Fujitsu Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^Fujitsu Fortran'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_FUJITSU],1,
      [Define to 1 if you are using the Fujitsu Fortran compiler.])
    afb_fc_vendor="fujitsu"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_FUJITSU



# _AFB_CHECK_FC_G95(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the G95 Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_G95],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the G95 Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^G95'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_G95],1,
      [Define to 1 if you are using the G95 Fortran compiler.])
    afb_fc_vendor="g95"
    afb_fc_version=`echo ${afb_result} | sed -e 's/.*GCC //; s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_G95



# _AFB_CHECK_FC_GNU(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the GNU Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_GNU],1,
      [Define to 1 if you are using the GNU Fortran compiler.])
    AC_DEFINE([HAVE_FORTRAN2003],1,
      [Define to 1 if your Fortran compiler supports Fortran 2003.])
    afb_fc_vendor="gnu"
    afb_fc_version=`echo ${afb_result} | sed -e 's/^[[^(]]*([[^)]]*) //; s/ .*//'`
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_GNU



# _AFB_CHECK_FC_HITACHI(COMPILER)
# -----------------------------------------
#
# Checks whether the specified Fortran compiler is the Hitachi Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_HITACHI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Hitachi Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^Hitachi Fortran'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_HITACHI],1,
      [Define to 1 if you are using the Hitachi Fortran compiler.])
    afb_fc_vendor="hitachi"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.*Driver //;s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_HITACHI



# _AFB_CHECK_FC_IBM(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
  fc_info_string=`$1 -qversion 2>&1 | head -n 1`
  fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  afb_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
  if test "${afb_result}" = ""; then
    afb_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
  fi
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
    if test "${fc_garbage}" -gt 50; then
      AC_DEFINE([FC_IBM],1,
        [Define to 1 if you are using the IBM XL Fortran compiler.])
      afb_fc_vendor="ibm"
      afb_fc_version="unknown"
      afb_result="yes"
    fi
  else
    AC_DEFINE([FC_IBM],1,
      [Define to 1 if you are using the IBM XL Fortran compiler.])
    afb_fc_vendor="ibm"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_IBM



# _AFB_CHECK_FC_INTEL(COMPILER)
# ---------------------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_INTEL],1,
      [Define to 1 if you are using the Intel Fortran compiler.])
    afb_fc_vendor="intel"
    afb_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Version //;s/ .*//'`
    if test "${afb_fc_version}" = ""; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_INTEL



# _AFB_CHECK_FC_MIPSPRO(COMPILER)
# -----------------------------------------
#
# Checks whether the specified Fortran compiler is the MIPSpro Fortran
# compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_MIPSPRO],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the MIPSpro Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^MIPSpro'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_MIPSPRO],1,
      [Define to 1 if you are using the MIPSpro Fortran compiler.])
    afb_fc_vendor="mipspro"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.*Version //'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_MIPSPRO



# _AFB_CHECK_FC_NAG(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the NAGWare Fortran 95
# compiler. If yes, tries to determine its version number and sets the
# afb_fc_vendor and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_NAG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the NAGWare Fortran 95 compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^NAG'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_NAG],1,
      [Define to 1 if you are using the NAGWare Fortran 95 compiler.])
    afb_fc_vendor="nag"
    afb_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Release //;s/[[( ]].*//'`
    if test "${afb_fc_version}" = ""; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_NAG



# _AFB_CHECK_FC_OPEN64(COMPILER)
# ----------------------------------------
#
# Checks whether the specified Fortran compiler is the Open64
# Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_OPEN64],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Open64 Fortran compiler])
  fc_info_string=`$1 --version 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^Open64'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_OPEN64],1,
      [Define to 1 if you are using the Open64 Fortran compiler.])
    afb_fc_vendor="open64"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_OPEN64



# _AFB_CHECK_FC_PATHSCALE(COMPILER)
# -------------------------------------------
#
# Checks whether the specified Fortran compiler is the PathScale
# Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_PATHSCALE],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the PathScale Fortran compiler])
  fc_info_string=`$1 -version 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep '^PathScale'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_PATHSCALE],1,
      [Define to 1 if you are using the PathScale Fortran compiler.])
    afb_fc_vendor="pathscale"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.* Version //; s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_PATHSCALE



# _AFB_CHECK_FC_PGI(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 2 | grep -v "^$"`
  afb_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_PGI],1,
      [Define to 1 if you are using the Portland Group Fortran compiler.])
    afb_fc_vendor="pgi"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
    if test "${afb_fc_version}" = "${afb_result}"; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_PGI



# _AFB_CHECK_FC_SUN(COMPILER)
# -------------------------------------
#
# Checks whether the specified Fortran compiler is the Sun WorkShop Fortran compiler.
# If yes, tries to determine its version number and sets the afb_fc_vendor
# and afb_fc_version variables accordingly.
#
AC_DEFUN([_AFB_CHECK_FC_SUN],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Sun Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  afb_result=`echo "${fc_info_string}" | grep 'Sun' | grep 'Fortran 95'`
  if test "${afb_result}" = ""; then
    afb_result="no"
    fc_info_string=""
    afb_fc_vendor="unknown"
    afb_fc_version="unknown"
  else
    AC_DEFINE([FC_SUN],1,
      [Define to 1 if you are using the Sun Fortran compiler.])
    afb_fc_vendor="sun"
    afb_fc_version=`echo "${afb_result}" | sed -e 's/.* Fortran 95 //;s/ .*//'`
    if test "${afb_fc_version}" = "${afb_result}" -o "${afb_fc_version}" = ""; then
      afb_fc_version="unknown"
    fi
    afb_result="yes"
  fi
  dnl AC_MSG_RESULT(${afb_result})
]) # _AFB_CHECK_FC_SUN



 #############################################################################



# AFB_FC_EXTENSIONS()
# -----------------------------
#
# Sets the default extensions of Fortran source files and modules,
# whenever possible.
#
AC_DEFUN([AFB_FC_EXTENSIONS],[
  dnl Set Fortran module extension
  AX_F90_MODULE_EXTENSION
  if test "${ax_cv_f90_modext}" != ""; then
    MODEXT="${ax_cv_f90_modext}"
  else
    MODEXT="mod"
    AC_MSG_NOTICE([setting Fortran module extension to ".${MODEXT}"])
  fi

  dnl Change the default Fortran extension for tests
  AC_FC_SRCEXT(F90,[afb_fc_src_ok="yes"],[afb_fc_src_ok="no"])
  if test "${afb_fc_src_ok}" != "yes"; then
    AC_MSG_WARN([Fortran file extension could not be changed])
    AC_MSG_WARN([some advanced Fortran tests may fail])
  fi
]) # AFB_FC_EXTENSIONS



# AFB_FC_MOD_CASE()
# ---------------------------
#
# Checks whether the Fortran compiler creates upper-case or lower-case
# module files.
#
AC_DEFUN([AFB_FC_MOD_CASE],[
  AC_REQUIRE([AFB_FC_EXTENSIONS])

  dnl Init
  fc_mod_lowercase="yes"
  fc_mod_uppercase="no"
  AC_MSG_NOTICE([determining Fortran module case])

  dnl Compile a dummy module
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([[
    module conftest
    end module conftest
  ]],[],[AC_MSG_FAILURE([unable to compile a simple Fortran module])])
  AC_LANG_POP([Fortran])

  dnl Check module file existence
  if test -f "CONFTEST.${MODEXT}"; then
    fc_mod_lowercase="no"
    fc_mod_uppercase="yes"
  elif test ! -f "conftest.${MODEXT}"; then
    AC_MSG_WARN([conftest.${MODEXT} Fortran module could not be found])
  fi

  dnl Output final outcome
  AC_MSG_CHECKING([whether Fortran modules are upper-case])
  AC_MSG_RESULT([${fc_mod_uppercase}])
]) # AFB_FC_MOD_CASE



# AFB_FC_MOD_INCS(MODULE)
# ---------------------------------
#
# Checks whether the specified Fortran module is directly available, or
# if we need to add '-I/usr/include' to the compile flags. Returns the
# required includes.
#
AC_DEFUN([AFB_FC_MOD_INCS],[
  AC_MSG_CHECKING([for Fortran module includes])

  if test "${afb_fc_mod_incs_ok}" = "" -o \
          "${afb_fc_mod_incs_ok}" = "unknown"; then

    dnl Init
    fc_mod_incs=""

    dnl Prepare environment
    tmp_saved_FCFLAGS="${FCFLAGS}"
    AC_LANG_PUSH([Fortran])

    dnl Look for module without includes
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use $1
      ]])], [afb_fc_mod_incs_ok="none required"], [afb_fc_mod_incs_ok="unknown"])

    dnl Look for module with includes
    if test "${afb_fc_mod_incs_ok}" = "unknown"; then
      FCFLAGS="${FCFLAGS} -I/usr/include"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
        [[
          use $1
        ]])],
        [afb_fc_mod_incs_ok="-I/usr/include"; fc_mod_incs="-I/usr/include"],
        [afb_fc_mod_incs_ok="unknown"])
    fi
    AC_MSG_RESULT([${afb_fc_mod_incs_ok}])

    dnl Restore environment
    AC_LANG_POP([Fortran])
    FCFLAGS="${tmp_saved_FCFLAGS}"

  else

    AC_MSG_RESULT([${afb_fc_mod_incs_ok} (cached)])

  fi

  dnl Substitute variables
  AC_SUBST(fc_mod_incs)
]) # AFB_FC_MOD_INCS



# AFB_PROG_FC()
# -----------------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([AFB_PROG_FC],[
  dnl Init
  afb_fc_vendor="${with_fc_vendor}"
  afb_fc_version="${with_fc_version}"

  if test "${afb_fc_vendor}" = ""; then
    afb_fc_vendor="unknown"
  fi
  if test "${afb_fc_version}" = ""; then
    afb_fc_version="unknown"
  fi
  afb_fc_wrap="no"

  dnl Determine Fortran compiler type (the order is important)
  AC_MSG_CHECKING([which type of Fortran compiler we have])

  dnl Clear temporary info file

  dnl Get rid of that one as early as possible
  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_IBM(${FC})
  fi

  dnl Should be checked before gfortran because it mimics its behaviour
  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_INTEL(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_G95(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_GNU(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_PATHSCALE(${FC})
  fi

  #if test "${afb_fc_vendor}" = "unknown"; then
  #  _AFB_CHECK_FC_COMPAQ(${FC})
  #fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_ABSOFT(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_MIPSPRO(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_OPEN64(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_FUJITSU(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_SUN(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_HITACHI(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_NAG(${FC})
  fi

  if test "${afb_fc_vendor}" = "unknown"; then
    _AFB_CHECK_FC_PGI(${FC})
  fi

  dnl Fall back to generic when detection fails
  if test "${afb_fc_vendor}" = "unknown"; then
    afb_fc_vendor="generic"
  fi

  dnl Normalize Fortran compiler version
  if test "${afb_fc_version}" = "unknown"; then
    afb_fc_version="0.0"
  else
    afb_fc_version=`echo ${afb_fc_version} | cut -d. -f1-2`
  fi

  dnl Display final result
  AC_MSG_RESULT([${afb_fc_vendor} ${afb_fc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(afb_fc_vendor)
  AC_SUBST(afb_fc_version)
  AC_SUBST(afb_fc_wrap)
  AC_SUBST(fc_info_string)
]) # AFB_PROG_FC

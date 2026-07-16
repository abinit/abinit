# -*- Autoconf -*-
#
# Copyright (C) 2005-2026 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# C++ compilers support
#


# _ABI_CXX_CHECK_ARM(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the ARMClang++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_ARM],[
  # Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl

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
]) # _ABI_CXX_CHECK_ARM


# _ABI_CXX_CHECK_GNU(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the GNU C++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_GNU],[
  # Do some sanity checking of the arguments
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
]) # _ABI_CXX_CHECK_GNU


# _ABI_CXX_CHECK_IBM(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the IBM XL C++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_IBM],[
  # Do some sanity checking of the arguments
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
]) # _ABI_CXX_CHECK_IBM


# _ABI_CXX_CHECK_INTEL(COMPILER)
# ------------------------------
#
# Checks whether the specified C++ compiler is an Intel C++ compiler
# (either the classic "icpc" or the LLVM-based "icpx"/oneAPI compiler).
# If yes, tries to determine its version number and sets the abi_cxx_vendor,
# abi_cxx_version and abi_cxx_flavor variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_INTEL],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl
  dnl AC_MSG_CHECKING([if we are using an Intel C++ compiler])
  cxx_command="$1"

  # Capture the full output (a deprecation remark from icpc, e.g.
  # "icpc: remark #10441: ...", may appear before the real version banner)
  # and isolate the version banner line, wherever it is.
  cxx_output=`$cxx_command -V 2>&1`
  version_line=`echo "${cxx_output}" | grep -E '^Intel\(R\) C\+\+ Intel\(R\) 64 Compiler|^Intel\(R\) oneAPI DPC\+\+/C\+\+ Compiler' | head -n 1`
  intel_check="${version_line}"

  # If using mpiicpc/mpiicpx, it may crash with "usage: mpiicpc"/"usage: mpiicpx"
  if test "${intel_check}" = ""; then
    usage_line=`echo "${cxx_output}" | grep '^usage:' | head -n 1`
    if test "${usage_line}" != ""; then
      fallback_cxx=`echo "${usage_line}" | cut -d " " -f 2`
      if command -v "${fallback_cxx}" >/dev/null 2>&1; then
        cxx_output=`${fallback_cxx} -V 2>&1`
        version_line=`echo "${cxx_output}" | grep -E '^Intel\(R\) C\+\+ Intel\(R\) 64 Compiler|^Intel\(R\) oneAPI DPC\+\+/C\+\+ Compiler' | head -n 1`
        intel_check="${version_line}"
        cxx_command="${fallback_cxx}"
      fi
    fi
  fi

  if test "${intel_check}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
    abi_cxx_flavor="unknown"
  else
    cxx_info_string="${version_line}"
    AC_DEFINE([CXX_INTEL],1,[Define to 1 if you are using an Intel C++ compiler.])
    abi_cxx_vendor="intel"

    # Extract the version number from the isolated banner line only,
    # never from the full (possibly multi-line) cxx_output.
    abi_cxx_version=`echo "${version_line}" | sed -e 's/.*Version //; s/ .*//'`
    if test "${abi_cxx_version}" = ""; then
      abi_cxx_version="unknown"
    fi

    # Distinguish classic icpc from the LLVM-based icpx (oneAPI) driver.
    # icpc (any version) repeats "Intel(R)" twice right at the start:
    #   "Intel(R) C++ Intel(R) 64 Compiler [Classic ]for applications ..."
    # icpx has a completely different banner:
    #   "Intel(R) oneAPI DPC++/C++ Compiler for applications running on ..."
    classic_check=`echo "${version_line}" | grep '^Intel(R) C++ Intel(R) 64 Compiler'`
    if test "${classic_check}" != ""; then
      abi_cxx_flavor="classic"
    else
      abi_cxx_flavor="oneapi"
      AC_DEFINE([CXX_INTEL_ONEAPI], 1,
        [Define to 1 if you are using the LLVM-based Intel C++ compiler (icpx).])
    fi

    abi_cxx_vendor="${abi_cxx_vendor} ${abi_cxx_flavor}"
    abi_result="yes (${abi_cxx_flavor})"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CXX_CHECK_INTEL


# _ABI_CXX_CHECK_CRAY(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the CRAY Clang++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_CRAY],[
  # Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the CRAY Clang++ C++ compiler])
  cxx_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${cxx_info_string}" | grep 'Cray clang'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_CRAY],1,
      [Define to 1 if you are using the CRAY Clang++ C++ compiler.])
    abi_cxx_vendor="cray"
    abi_cxx_version=`echo ${abi_result} | sed -e 's/.*ersion //; s/ .*//'`
    if test "${abi_cxx_version}" = "${abi_result}"; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CXX_CHECK_CRAY


# _ABI_CXX_CHECK_LLVM(COMPILER)
# -----------------------------
#
# Checks whether the specified C compiler is the LLVM Clang++ compiler.
# If yes, tries to determine its version number and sets the abi_cxx_vendor
# and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_LLVM],[
  # Do some sanity checking of the arguments
  m4_if([$1], [], [AC_FATAL([$0: missing argument 1])])dnl

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
]) # _ABI_CXX_CHECK_LLVM


# _ABI_CXX_CHECK_NVHPC(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the NVIDIA HPC SDK C++
# compiler. If yes, tries to determine its version number and sets the
# abi_cxx_vendor and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_NVHPC],[
  # Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the NVIDIA HPC SDK C++ compiler])
  cxx_info_string=`$1 -V 2> /dev/null | grep "^nvc++"`
  abi_result=`echo "${cxx_info_string}"`
  if test "${abi_result}" = ""; then
    abi_result="no"
    cxx_info_string=""
    abi_cxx_vendor="unknown"
    abi_cxx_version="unknown"
  else
    AC_DEFINE([CXX_NVHPC],1,[Define to 1 if you are using the NVIDIA HPC SDK C++ compiler.])
    abi_cxx_vendor="nvhpc"
    abi_cxx_version=`echo "${abi_result}" | cut -f2 -d" "`
    if test "${abi_cxx_version}" = ""; then
      abi_cxx_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CXX_CHECK_NVHPC


# _ABI_CXX_CHECK_PGI(COMPILER)
# ----------------------------
#
# Checks whether the specified C++ compiler is the Portland Group C++
# compiler. If yes, tries to determine its version number and sets the
# abi_cxx_vendor and abi_cxx_version variables accordingly.
#
AC_DEFUN([_ABI_CXX_CHECK_PGI],[
  # Do some sanity checking of the arguments
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
]) # _ABI_CXX_CHECK_PGI


# ABI_PROG_CXX()
# --------------
#
# Tries to determine which type of C++ compiler is installed.
#
AC_DEFUN([ABI_PROG_CXX],[
  # Init
  if test "${abi_cxx_vendor}" = ""; then
    abi_cxx_vendor="unknown"
  fi

  # Preserve environment
  ABI_ENV_BACKUP

  # Look for the C++ compiler
  if test "${CXX}" != "" -a ! -x "${CXX}"; then
    abi_cxx_probe=`echo "${CXX}" | sed -e 's/ .*//'`
    if test ! -x "${abi_cxx_probe}"; then
      AC_PATH_PROG([abi_cxx_path],[${abi_cxx_probe}])
      if test "${abi_cxx_path}" = ""; then
        AC_MSG_ERROR([could not run C++ compiler "${CXX}"])
      fi
    fi
  fi
  AC_PROG_CXX([mpicxx mpiicpx cxx c++ icx icpx xlC CXX g++ nvc++ clang++])

  # Warn if no C++ compiler is available
  if test "${CXX}" = ""; then
    AC_MSG_WARN([no C++ compiler available])
  fi

  # Determine C++ compiler type (the order is important)
  AC_MSG_CHECKING([which type of C++ compiler we have])

  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_IBM(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_ARM(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_INTEL(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_CRAY(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_LLVM(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_NVHPC(${CXX})
  fi
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_PGI(${CXX})
  fi

  # Check the GNU compiler last, because other compilers are cloning
  # its CLI
  if test "${abi_cxx_vendor}" = "unknown"; then
    _ABI_CXX_CHECK_GNU(${CXX})
  fi

  # Fall back to generic when detection fails
  if test "${abi_cxx_vendor}" = "unknown"; then
    abi_cxx_vendor="generic"
    abi_cxx_version="0.0"
  fi

  # Normalize C++ compiler version
  abi_cxx_version=`echo ${abi_cxx_version} | cut -d. -f1-2`

  # Display final result
  AC_MSG_RESULT([${abi_cxx_vendor} ${abi_cxx_version}])

  # Restore back CXXFLAGS
  CXXFLAGS="${abi_env_CXXFLAGS}"

  # Schedule compiler info for substitution
  AC_SUBST(abi_cxx_vendor)
  AC_SUBST(abi_cxx_version)
  AC_SUBST(cxx_info_string)
]) # ABI_PROG_CXX

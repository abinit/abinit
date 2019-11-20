# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran compilers support
#



# _ABI_CHECK_FC_ABSOFT(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ABSoft Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ABSOFT],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ABSoft Fortran compiler])
  fc_info_string=`$1 -V 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep 'Pro Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_ABSOFT],1,
      [Define to 1 if you are using the ABSOFT Fortran compiler.])
    abi_fc_vendor="absoft"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.*Pro Fortran //'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ABSOFT



# _ABI_CHECK_FC_ARM(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the ARMFlang Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_ARM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the ARM Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Arm C/C++/Fortran Compiler'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_ARM],1,
      [Define to 1 if you are using the ARM Fortran compiler.])
    abi_fc_vendor="arm"
    abi_fc_version=`echo ${abi_result} | sed -e 's/.*ersion //; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_ARM



# _ABI_CHECK_FC_GNU(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the GNU Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_GNU],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the GNU Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^GNU Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_GNU],1,
      [Define to 1 if you are using the GNU Fortran compiler.])
    AC_DEFINE([HAVE_FORTRAN2003],1,
      [Define to 1 if your Fortran compiler supports Fortran 2003.])
    abi_fc_vendor="gnu"
    abi_fc_version=`echo ${abi_result} | sed -e 's/^[[^(]]*([[^)]]*) //; s/ .*//'`
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_GNU



# _ABI_CHECK_FC_IBM(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the IBM XL Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_IBM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the IBM XL Fortran compiler])
  fc_info_string=`$1 -qversion 2>&1 | head -n 1`
  fc_garbage=`$1 -qversion 2>&1 | wc -l | sed -e 's/ //g'`
  abi_result=`echo "${fc_info_string}" | grep 'IBM XL Fortran'`
  if test "${abi_result}" = ""; then
    abi_result=`echo "${fc_info_string}" | grep 'IBM(R) XL Fortran'`
  fi
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
    if test "${fc_garbage}" -gt 50; then
      AC_DEFINE([FC_IBM],1,
        [Define to 1 if you are using the IBM XL Fortran compiler.])
      abi_fc_vendor="ibm"
      abi_fc_version="unknown"
      abi_result="yes"
    fi
  else
    AC_DEFINE([FC_IBM],1,
      [Define to 1 if you are using the IBM XL Fortran compiler.])
    abi_fc_vendor="ibm"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/.* V//; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_IBM



# _ABI_CHECK_FC_INTEL(COMPILER)
# -----------------------------
#
# Checks whether the specified Fortran compiler is the Intel Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_INTEL],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Intel Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^Intel(R) Fortran'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_INTEL],1,
      [Define to 1 if you are using the Intel Fortran compiler.])
    abi_fc_vendor="intel"
    abi_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Version //;s/ .*//'`
    if test "${abi_fc_version}" = ""; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_INTEL



# _ABI_CHECK_FC_LLVM(COMPILER)
# ------------------------------
#
# Checks whether the specified Fortran compiler is the LLVM Flang compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_LLVM],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the LLVM Flang Fortran compiler])
  fc_info_string=`$1 --version 2>/dev/null | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep -e '^[[CcFf]]lang'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_LLVM],1,
      [Define to 1 if you are using the LLVM Flang Fortran compiler.])
    abi_fc_vendor="llvm"
    abi_fc_version=`echo ${abi_result} | sed -e 's/.*ersion //; s/ .*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_LLVM



# _ABI_CHECK_FC_NAG(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the NAGWare Fortran 95
# compiler. If yes, tries to determine its version number and sets the
# abi_fc_vendor and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_NAG],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the NAGWare Fortran 95 compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 1`
  abi_result=`echo "${fc_info_string}" | grep '^NAG'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_NAG],1,
      [Define to 1 if you are using the NAGWare Fortran 95 compiler.])
    abi_fc_vendor="nag"
    abi_fc_version=`echo "${fc_info_string}" | sed -e 's/.*Release //;s/[[( ]].*//'`
    if test "${abi_fc_version}" = ""; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_NAG



# _ABI_CHECK_FC_PGI(COMPILER)
# ---------------------------
#
# Checks whether the specified Fortran compiler is the Portland Group
# Fortran compiler.
# If yes, tries to determine its version number and sets the abi_fc_vendor
# and abi_fc_version variables accordingly.
#
AC_DEFUN([_ABI_CHECK_FC_PGI],[
  dnl Do some sanity checking of the arguments
  m4_if([$1], , [AC_FATAL([$0: missing argument 1])])dnl

  dnl AC_MSG_CHECKING([if we are using the Portland Group Fortran compiler])
  fc_info_string=`$1 -V 2>&1 | head -n 2 | grep -v "^$"`
  abi_result=`echo "${fc_info_string}" | grep '^pgf9[[05]]'`
  if test "${abi_result}" = ""; then
    abi_result="no"
    fc_info_string=""
    abi_fc_vendor="unknown"
    abi_fc_version="unknown"
  else
    AC_DEFINE([FC_PGI],1,
      [Define to 1 if you are using the Portland Group Fortran compiler.])
    abi_fc_vendor="pgi"
    abi_fc_version=`echo "${abi_result}" | sed -e 's/^pgf9[[05]] //' | sed -e 's/-.*//'`
    if test "${abi_fc_version}" = "${abi_result}"; then
      abi_fc_version="unknown"
    fi
    abi_result="yes"
  fi
  dnl AC_MSG_RESULT(${abi_result})
]) # _ABI_CHECK_FC_PGI



 ##############################################################################

# _ABI_CHECK_FC_ASYNC()
# --------------------
#
# Checks whether the Fortran compiler supports the ASYNCHRONOUS attribute (F2003).
#
AC_DEFUN([_ABI_CHECK_FC_ASYNC],[
  dnl Init
  fc_has_async="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts the ASYNCHRONOUS attribute])

  dnl Try to compile a program using asynchronous arrays (F2003)
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
     INTEGER, ASYNCHRONOUS :: int_array(10)
    ]])], [fc_has_async="yes"])
  AC_LANG_POP()

  if test "${fc_has_async}" = "yes"; then
    AC_DEFINE([HAVE_FC_ASYNC],1,
      [Define to 1 if your Fortran compiler supports the asynchronous attribute.])
  fi

  AC_MSG_RESULT(${fc_has_async})
]) # _ABI_CHECK_FC_ASYNC


 ##############################################################################

# _ABI_CHECK_FC_BACKTRACE()
# --------------------
#
# Checks whether the Fortran compiler supports BACKTRACE (gfortran extension added in 4.8)
#
AC_DEFUN([_ABI_CHECK_FC_BACKTRACE],[
  dnl Init
  fc_has_backtrace="no"

  AC_MSG_CHECKING([whether the Gfortran compiler supports BACKTRACE])

  dnl Try to compile a piece of code that calls BACKTRACE.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      call backtrace()
      
    ]])], [fc_has_backtrace="yes"])
  AC_LANG_POP()

  if test "${fc_has_backtrace}" = "yes"; then
    AC_DEFINE([HAVE_FC_BACKTRACE],1, 
      [Define to 1 if your Fortran compiler supports BACKTRACE.])
  fi

  AC_MSG_RESULT(${fc_has_backtrace})
]) # _ABI_CHECK_FC_BACKTRACE


 ##############################################################################

# _ABI_CHECK_FC_COMMAND_ARGUMENT()
# --------------------
#
# Checks whether the Fortran compiler supports GET_COMMAND_ARGUMENT.
#
AC_DEFUN([_ABI_CHECK_FC_COMMAND_ARGUMENT],[
  dnl Init
  fc_has_command_argument="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports GET_COMMAND_ARGUMENT])

  dnl Try to compile a piece of code that calls get_command_argument.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      integer :: ii
      character(len=500) :: arg

      call get_command(arg)

      do ii=1,command_argument_count()
        call get_command_argument(ii, arg)
      end do
      
    ]])], [fc_has_command_argument="yes"])
  AC_LANG_POP()

  if test "${fc_has_command_argument}" = "yes"; then
    AC_DEFINE([HAVE_FC_COMMAND_ARGUMENT],1, 
      [Define to 1 if your Fortran compiler supports GET_COMMAND_ARGUMENT.])
  fi

  AC_MSG_RESULT(${fc_has_command_argument})
]) # _ABI_CHECK_FC_COMMAND_ARGUMENT


 ##############################################################################

# _ABI_CHECK_FC_COMMAND_LINE()
# --------------------
#
# Checks whether the Fortran compiler supports EXECUTE_COMMAND_LINE.
#
AC_DEFUN([_ABI_CHECK_FC_COMMAND_LINE],[
  dnl Init
  fc_has_command_line="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports EXECUTE_COMMAND_LINE])

  dnl Try to compile a piece of code that calls execute_command_line.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      integer :: i
      
      call execute_command_line ("external_prog.exe", exitstat=i)
      print *, "Exit status of external_prog.exe was ", i
      
      call execute_command_line ("reindex_files.exe", wait=.false.)
      print *, "Now reindexing files in the background"
    ]])], [fc_has_command_line="yes"])
  AC_LANG_POP()

  if test "${fc_has_command_line}" = "yes"; then
    AC_DEFINE([HAVE_FC_COMMAND_LINE],1, 
      [Define to 1 if your Fortran compiler supports EXECUTE_COMMAND_LINE.])
  fi

  AC_MSG_RESULT(${fc_has_command_line})
]) # _ABI_CHECK_FC_COMMAND_LINE


 ##############################################################################

# _ABI_CHECK_FC_SYSTEM()
# --------------------
#
# Checks whether the Fortran compiler supports SYSTEM.
#
AC_DEFUN([_ABI_CHECK_FC_SYSTEM],[
  dnl Init
  fc_has_system="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports SYSTEM])

  dnl Try to compile a piece of code that calls system.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      call system ("ls -l")
    ]])], [fc_has_system="yes"])
  AC_LANG_POP()

  if test "${fc_has_system}" = "yes"; then
    AC_DEFINE([HAVE_FC_SYSTEM],1, 
      [Define to 1 if your Fortran compiler supports SYSTEM.])
  fi

  AC_MSG_RESULT(${fc_has_system})
]) # _ABI_CHECK_FC_SYSTEM


# _ABI_CHECK_FC_CONTIGUOUS()
# --------------------
#
# Checks whether the Fortran compiler supports the CONTIGUOUS attribute (F2008).
#
AC_DEFUN([_ABI_CHECK_FC_CONTIGUOUS],[
  dnl Init
  fc_has_contiguous="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts the CONTIGUOUS attribute])

  dnl Try to compile a program using contiguous (F2008)
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
     integer, parameter :: dp=kind(1.0d0)
     integer, parameter :: dpc=kind((1.0_dp,1.0_dp)) 

     integer,contiguous,pointer :: i_ptr(:)
     real(dp),contiguous,pointer :: r_ptr(:,:)
     complex(dpc),contiguous,pointer :: c_ptr(:,:,:)
    ]])], [fc_has_contiguous="yes"])
  AC_LANG_POP()

  if test "${fc_has_contiguous}" = "yes"; then
    AC_DEFINE([HAVE_FC_CONTIGUOUS],1,
      [Define to 1 if your Fortran compiler supports the contiguous attribute.])
  fi

  AC_MSG_RESULT(${fc_has_contiguous})
]) # _ABI_CHECK_FC_CONTIGUOUS


# _ABI_CHECK_FC_EXIT()
# --------------------
#
# Checks whether the Fortran compiler supports the exit() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_EXIT],[
  dnl Init
  fc_has_exit="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts exit()])

  dnl Try to compile a program calling exit()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call exit(1)
    ]])], [fc_has_exit="yes"])
  AC_LANG_POP()

  if test "${fc_has_exit}" = "yes"; then
    AC_DEFINE([HAVE_FC_EXIT],1,
      [Define to 1 if your Fortran compiler supports exit().])
  fi

  AC_MSG_RESULT(${fc_has_exit})
]) # _ABI_CHECK_FC_EXIT



# _ABI_CHECK_FC_FLUSH()
# ---------------------
#
# Checks whether the Fortran compiler supports the flush() subroutine.
#
AC_DEFUN([_ABI_CHECK_FC_FLUSH],[
  dnl Init
  fc_has_flush="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts flush()])

  dnl Try to compile a program calling flush()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call flush(1)
    ]])], [fc_has_flush="yes"])
  AC_LANG_POP()

  if test "${fc_has_flush}" = "yes"; then
    AC_DEFINE([HAVE_FC_FLUSH],1,
      [Define to 1 if your Fortran compiler supports flush().])
  fi

  AC_MSG_RESULT(${fc_has_flush})
]) # _ABI_CHECK_FC_FLUSH


# _ABI_CHECK_FC_FLUSH_()
# ----------------------
#   
# Checks whether the Fortran compiler supports the flush_() subroutine.
# 
AC_DEFUN([_ABI_CHECK_FC_FLUSH_],[
  dnl Init
  fc_has_flush_="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts flush_()])

  dnl Try to compile a program calling flush_()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            call flush_(1)
    ]])], [fc_has_flush_="yes"])
  AC_LANG_POP()

  if test "${fc_has_flush_}" = "yes"; then
    AC_DEFINE([HAVE_FC_FLUSH_],1,
      [Define to 1 if your Fortran compiler supports flush_().])
  fi

  AC_MSG_RESULT(${fc_has_flush_})
]) # _ABI_CHECK_FC_FLUSH_


# _ABI_CHECK_FC_GAMMA()
# ---------------------
#
# Checks whether the Fortran compiler supports the gamma() intrinsic
# (Fortran 2003 and later).
#
AC_DEFUN([_ABI_CHECK_FC_GAMMA],[
  dnl Init
  fc_has_gamma="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts gamma()])

  dnl Try to compile a program using gamma()
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
            real :: x
            x = gamma(1.5)
    ]])], [fc_has_gamma="yes"])
  AC_LANG_POP()

  if test "${fc_has_gamma}" = "yes"; then
    AC_DEFINE([HAVE_FC_GAMMA],1,
      [Define to 1 if your Fortran compiler supports gamma().])
  fi

  AC_MSG_RESULT(${fc_has_gamma})
]) # _ABI_CHECK_FC_GAMMA


# _ABI_CHECK_FC_SHIFTLR()
# ----------------------
#
# Checks whether the Fortran compiler supports SHIFTL/SHIFTR
# (Fortran 2008 and later).
#
AC_DEFUN([_ABI_CHECK_FC_SHIFTLR],[
  dnl Init
  fc_has_shiftlr="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts shiftl() and shiftr()])

  dnl Try to compile a call to cpu_time
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      integer :: ii,ishft,res
      res=shiftl(ii,ishft)
      res=shiftr(ii,ishft)

    ]])], [fc_has_shiftlr="yes"])
  AC_LANG_POP()

  if test "${fc_has_shiftlr}" = "yes"; then
    AC_DEFINE([HAVE_FC_SHIFTLR],1,
      [Define to 1 if your Fortran compiler supports shiftl() and shiftr().])
  fi

  AC_MSG_RESULT(${fc_has_shiftlr})
]) # _ABI_CHECK_FC_SHIFTLR


# _ABI_CHECK_FC_GETENV()
# ----------------------
#
# Checks whether the Fortran compiler supports GET_ENVIRONMENT_VARIABLE
# (Fortran 2003 and later).
#
AC_DEFUN([_ABI_CHECK_FC_GETENV],[
  dnl Init
  fc_has_getenv="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts getenv()])

  dnl Try to compile a call to getenv
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      character(len=255) :: homedir
      call getenv("HOME", homedir)

    ]])], [fc_has_getenv="yes"])
  AC_LANG_POP()

  if test "${fc_has_getenv}" = "yes"; then
    AC_DEFINE([HAVE_FC_GETENV],1, 
      [Define to 1 if your Fortran compiler supports getenv().])
  fi

  AC_MSG_RESULT(${fc_has_getenv})
]) # _ABI_CHECK_FC_GETENV



# _ABI_CHECK_FC_INT_QUAD()
# ------------------------
#
# Checks whether the Fortran compiler supports quadruple integers.
#
AC_DEFUN([_ABI_CHECK_FC_INT_QUAD],[
  dnl Init
  fc_has_int_quad="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts quadruple integers])

  dnl Try to compile a program defining a quadruple integer
  dnl Note: xlf "works around" the problem by changing the integer length
  dnl Note: need to test "integer*16" and "integer(kind=16)" (seems not equivalent)
  if test "${abi_fc_vendor}" != "ibm"; then
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
              integer*16 my_int
              integer(kind=16) :: my_int2
      ]])], [fc_has_int_quad="yes"])
    AC_LANG_POP()
  fi

  if test "${fc_has_int_quad}" = "yes"; then
    AC_DEFINE([HAVE_FC_INT_QUAD],1,
      [Define to 1 if your Fortran compiler accepts quadruple integers.])
  fi

  AC_MSG_RESULT(${fc_has_int_quad})
]) # _ABI_CHECK_FC_INT_QUAD


# _ABI_CHECK_FC_ISO_FORTRAN_2008()
# ------------------------
#
# Checks whether the Fortran compiler supports 2008 standard in ISO_FORTRAN_ENV.
#
AC_DEFUN([_ABI_CHECK_FC_ISO_FORTRAN_2008],[
  dnl Init
  fc_has_iso_fortran_2008="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports 2008 standard in ISO_FORTRAN_ENV])

  dnl Try to compile a program using ISO_FORTRAN_ENV with int16,int32,int64
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
          use ISO_FORTRAN_ENV, only : int16,int32,int64
      ]])], [fc_has_iso_fortran_2008="yes"])
  AC_LANG_POP()

  if test "${fc_has_iso_fortran_2008}" = "yes"; then
    AC_DEFINE([HAVE_FC_ISO_FORTRAN_2008],1,
      [Define to 1 if your Fortran compiler supports 2008 standard in ISO_FORTRAN_ENV.])
  fi

  AC_MSG_RESULT(${fc_has_iso_fortran_2008})
]) # _ABI_CHECK_FC_ISO_FORTRAN_2008


# _ABI_CHECK_FC_DTARRAYS()
# --------------------
#
# Checks whether the Fortran compiler supports allocatable arrays in Fortran datatypes.
#
AC_DEFUN([_ABI_CHECK_FC_DTARRAYS],[
  dnl Init
  fc_has_dtarrays="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports allocatable arrays in datatypes])

  dnl Try to compile a type with an allocatable array
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[

         integer, parameter :: dp=kind(1.0d0)
         integer, parameter :: dpc=kind((1.0_dp,1.0_dp))  

         type test_type
           integer,allocatable :: i(:) 
           real(dp),allocatable :: r(:,:) 
           complex(dpc),allocatable :: c(:,:,:) 
         end type test_type

    ]])], [fc_has_dtarrays="yes"])
  AC_LANG_POP()

  if test "${fc_has_dtarrays}" = "yes"; then
    AC_DEFINE([HAVE_FC_ALLOCATABLE_DTARRAYS],1, 
      [Define to 1 if your Fortran compiler supports allocatable arrays in datatypes.])
  fi

  AC_MSG_RESULT(${fc_has_dtarrays})
]) # _ABI_CHECK_FC_DTARRAYS


# _ABI_CHECK_FC_IEEE_ARITHMETIC()
# --------------------
#
# Checks whether the Fortran compiler supports the intrinsic module IEEE_ARITHMETIC
#
AC_DEFUN([_ABI_CHECK_FC_IEEE_ARITHMETIC],[
  dnl Init
  fc_has_ieee_arithmetic="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports IEEE_ARITHMETIC])

  dnl Try to compile a piece of code that uses the module.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      use, intrinsic :: ieee_arithmetic
      real :: val

      if (ieee_is_nan(val)) then  ! NaN
        write(*,*)"Hello NAN"
      end if

    ]])], [fc_has_ieee_arithmetic="yes"])
  AC_LANG_POP()

  if test "${fc_has_ieee_arithmetic}" = "yes"; then
    AC_DEFINE([HAVE_FC_IEEE_ARITHMETIC],1, 
      [Define to 1 if your Fortran compiler supports IEEE_ARITHMETIC module.])
  fi

  AC_MSG_RESULT(${fc_has_ieee_arithmetic})
]) # _ABI_CHECK_FC_IEEE_ARITHMETIC


 ##############################################################################

# _ABI_CHECK_FC_IEEE_EXCEPTIONS()
# --------------------
#
# Checks whether the Fortran compiler supports the intrinsic module IEEE_EXCEPTIONS
#
AC_DEFUN([_ABI_CHECK_FC_IEEE_EXCEPTIONS],[
  dnl Init
  fc_has_ieee_exceptions="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports IEEE_EXCEPTIONS])

  dnl Try to compile a piece of code that uses the module.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      use, intrinsic :: ieee_exceptions 
      type(ieee_status_type) :: status_value

      call ieee_get_status(status_value)   ! Get the flags
      call ieee_set_flag(ieee_all,.false.) ! Set the flags quiet
      call ieee_set_status(status_value)   ! Restore the flags
      
    ]])], [fc_has_ieee_exceptions="yes"])
  AC_LANG_POP()

  if test "${fc_has_ieee_exceptions}" = "yes"; then
    AC_DEFINE([HAVE_FC_IEEE_EXCEPTIONS],1, 
      [Define to 1 if your Fortran compiler supports IEEE_EXCEPTIONS.])
  fi

  AC_MSG_RESULT(${fc_has_ieee_exceptions})
]) # _ABI_CHECK_FC_IEEE_EXCEPTIONS


# _ABI_CHECK_FC_IOMSG()
# --------------------
#
# Checks whether the Fortran compiler supports IOMSG.
#
AC_DEFUN([_ABI_CHECK_FC_IOMSG],[
  dnl Init
  fc_has_iomsg="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports IOMSG])

  dnl Try to compile a piece of code that opens, reads, writes and closes a file using iomsg.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         IMPLICIT NONE
         CHARACTER(len=500) :: ERRMSG, DUMMY

         OPEN(UNIT=11, FILE="iomsg.demo", STATUS="NEW", IOMSG=ERRMSG)
         WRITE(11,IOMSG=ERRMSG) "first"
         READ(11, IOMSG=ERRMSG) DUMMY
         CLOSE(UNIT=11, IOMSG=ERRMSG)

    ]])], [fc_has_iomsg="yes"])
  AC_LANG_POP()

  if test "${fc_has_iomsg}" = "yes"; then
    AC_DEFINE([HAVE_FC_IOMSG],1, 
      [Define to 1 if your Fortran compiler supports IOMSG.])
  fi

  AC_MSG_RESULT(${fc_has_iomsg})
]) # _ABI_CHECK_FC_IOMSG


# _ABI_CHECK_FC_ISO_C_BINDING()
# -----------------------------
#
# Checks whether the Fortran compiler provides the intrinsic module ISO_C_BINDING.
#
AC_DEFUN([_ABI_CHECK_FC_ISO_C_BINDING],[
  dnl Init
  fc_has_iso_c_binding="no"

  AC_MSG_CHECKING([whether the Fortran compiler provides the iso_c_binding module])

  dnl Try to compile a simple piece of code using iso_c_binding
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         use iso_c_binding
         implicit none
         integer(c_int) :: ii
         logical :: lbool
         type(c_ptr) :: ptr
         ptr = c_null_ptr
         lbool = c_associated(ptr) 

    ]])], [fc_has_iso_c_binding="yes"])
  AC_LANG_POP()

  if test "${fc_has_iso_c_binding}" = "yes"; then
    AC_DEFINE([HAVE_FC_ISO_C_BINDING],1, 
      [Define to 1 if your Fortran compiler provides the iso_c_binding module.])
  else
    AC_MSG_ERROR([Fortran compiler does not provide iso_c_binding module. Use a more recent version or a different compiler])
  fi

  AC_MSG_RESULT(${fc_has_iso_c_binding})
]) # _ABI_CHECK_FC_ISO_C_BINDING



# _ABI_CHECK_FC_LONG_LINES()
# --------------------------
# 
# Checks whether the Fortran compiler supports long lines.
#
AC_DEFUN([_ABI_CHECK_FC_LONG_LINES],[
  dnl Init
  fc_has_long_lines="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts long lines])

  dnl Try to compile a single line exceeding 136 columns.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
         write(*,*)"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" !142
    ]])], [fc_has_long_lines="yes"])
  AC_LANG_POP()

  dnl This is not correctly implemented on LLVM
  if test "${abi_fc_vendor}" = "llvm" -o "${abi_fc_vendor}" = "arm" ; then
    fc_has_long_lines="no"
  fi

  if test "${fc_has_long_lines}" = "yes"; then
    AC_DEFINE([HAVE_FC_LONG_LINES],1, 
      [Define to 1 if your Fortran compiler supports long lines.])
  fi

  AC_MSG_RESULT(${fc_has_long_lines})

  dnl Official macro added in autoconf ??
  dnl AC_FC_LINE_LENGTH ([length], [action-if-success], [action-if-failure = AC_MSG_FAILURE])

]) # _ABI_CHECK_FC_LONG_LINES


# _ABI_CHECK_FC_MACRO_NEWLINE()
# --------------------------
#
# Checks whether the Fortran preprocessor supports \newline in macros.
#
AC_DEFUN([_ABI_CHECK_FC_MACRO_NEWLINE],[
  dnl Init
  fc_has_macro_newline="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports \newline in CPP macros])

  dnl Try to compile a piece of code that opens a file using \newline in a macro (has to use F90 ext).
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
#        define NEWLINE \newline
         print *,"foo1" NEWLINE print *,"foo2"

    ]])], [fc_has_macro_newline="yes"])
  AC_LANG_POP()

  if test "${fc_has_macro_newline}" = "yes"; then
    AC_DEFINE([HAVE_FC_MACRO_NEWLINE],1,
      [Define to 1 if your Fortran compiler supports \newline in a macros.])
  fi

  AC_MSG_RESULT(${fc_has_macro_newline})
]) # _ABI_CHECK_FC_MACRO_NEWLINE


# _ABI_CHECK_FC_MOVE_ALLOC()
# --------------------------
#
# Checks whether the Fortran compile supports MOVE_ALLOC.
#
AC_DEFUN([_ABI_CHECK_FC_MOVE_ALLOC],[
  dnl Init
  fc_has_move_alloc="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports MOVE_ALLOC (F2003)])

  dnl Try to compile a piece of code that uses move_alloc (F2003)
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
     integer, allocatable :: a(:), b(:)
     allocate(a(3))
     a = (/1, 2, 3/)
     call move_alloc(a, b)
    ]])], [fc_has_move_alloc="yes"])
  AC_LANG_POP()

  if test "${fc_has_move_alloc}" = "yes"; then
    AC_DEFINE([HAVE_FC_MOVE_ALLOC],1,
      [Define to 1 if your Fortran compiler supports MOVE_ALLOC (F2003).])
  fi

  AC_MSG_RESULT(${fc_has_move_alloc})
]) # _ABI_CHECK_FC_MOVE_ALLOC


# _ABI_CHECK_FC_PRIVATE()
# --------------------
#
# Checks whether the Fortran compiler supports the PRIVATE attribute (F2003).
#
AC_DEFUN([_ABI_CHECK_FC_PRIVATE],[
  dnl Init
  fc_has_private="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts the PRIVATE attribute])

  dnl Try to compile a program using private entities (F2003)
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([
     module foo
       type, public :: bar_t
         integer :: pub
         integer,private :: priv
       end type bar_t
     end module foo
    ], [fc_has_private="yes"])
  AC_LANG_POP()

  if test "${fc_has_private}" = "yes"; then
    AC_DEFINE([HAVE_FC_PRIVATE],1,
      [Define to 1 if your Fortran compiler supports the private attribute.])
  fi

  dnl remove the module.
  if test -f "FOO.${MODEXT}"; then
    rm -f "FOO.${MODEXT}"; 
  elif test -f "foo.${MODEXT}"; then
    rm -f "foo.${MODEXT}"; 
  fi

  AC_MSG_RESULT(${fc_has_private})
]) # _ABI_CHECK_FC_PRIVATE


 ##############################################################################

# _ABI_CHECK_FC_PROTECTED()
# --------------------
#
# Checks whether the Fortran compiler supports the PROTECTED attribute (F2003).
#
AC_DEFUN([_ABI_CHECK_FC_PROTECTED],[
  dnl Init
  fc_has_protected="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts the PROTECTED attribute])

  dnl Try to compile a program using protected module entities (F2003)
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([
    module foo
     real,save,protected :: aprot(10)
    end module foo
    ], [fc_has_protected="yes"])
  AC_LANG_POP()

  if test "${fc_has_protected}" = "yes"; then
    AC_DEFINE([HAVE_FC_PROTECTED],1,
      [Define to 1 if your Fortran compiler supports the protected attribute.])
  fi

  dnl remove the module.
  if test -f "FOO.${MODEXT}"; then
    rm -f "FOO.${MODEXT}"; 
  elif test -f "foo.${MODEXT}"; then
    rm -f "foo.${MODEXT}"; 
  fi

  AC_MSG_RESULT(${fc_has_protected})
]) # _ABI_CHECK_FC_PROTECTED


# _ABI_CHECK_FC_STREAM_IO()
# --------------------
#
# Checks whether the Fortran compiler supports stream IO (F2003)
#
AC_DEFUN([_ABI_CHECK_FC_STREAM_IO],[
  dnl Init
  fc_has_stream_io="no"

  AC_MSG_CHECKING([whether the Fortran compiler supports stream IO])

  dnl Try to compile a piece of code that opens a file using unformatted stream access.
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
         IMPLICIT NONE
         INTEGER :: myvalue = 12345, mypos
         OPEN(UNIT=11, FILE="ustream.demo", STATUS="NEW", ACCESS="STREAM")
         WRITE(11) "first"
         WRITE(11) "second"
         INQUIRE(UNIT=11, POS=mypos)
         PRINT *, "Myvalue will be written at position ", mypos
         WRITE(11) myvalue
         CLOSE(UNIT=11)

    ]])], [fc_has_stream_io="yes"])
  AC_LANG_POP()

  if test "${fc_has_stream_io}" = "yes"; then
    AC_DEFINE([HAVE_FC_STREAM_IO],1,
      [Define to 1 if your Fortran compiler supports stream IO.])
  fi

  AC_MSG_RESULT(${fc_has_stream_io})
]) # _ABI_CHECK_FC_STREAM_IO


# _ABI_CHECK_FC_TIMING()
# ----------------------
#
# Tries to determine which Fortran timing routines are available. (F95)
#
AC_DEFUN([_ABI_CHECK_FC_TIMING],[
  dnl Init
  fc_timing="standard"
  fc_has_etime="no"

  dnl Look for etime() support
  if test "${fc_timing}" = "standard"; then
    AC_MSG_CHECKING([whether the Fortran compiler accepts etime()])

    dnl Try to compile a program calling etime()
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
            call etime(1)
      ]])], [fc_has_etime="yes"])
    AC_LANG_POP()

    if test "${fc_has_etime}" = "yes"; then
      AC_DEFINE([HAVE_FC_ETIME],1,
        [Define to 1 if your Fortran compiler supports etime().])
    fi

    AC_MSG_RESULT(${fc_has_etime})

  fi

  dnl Determine whether to use C clock for timings
  AC_MSG_CHECKING([whether to use C clock for timings])
  AC_MSG_RESULT([${enable_cclock}])
  if test "${enable_cclock}" = "yes"; then
    AC_DEFINE([HAVE_CCLOCK],1,[Use C clock for timings.])
    fc_timing="cclock"
  fi
  AM_CONDITIONAL(DO_BUILD_CCLOCK,[test "${enable_cclock}" = "yes"])

  dnl Schedule info for substitution
  AC_SUBST(fc_timing)
]) # _ABI_CHECK_FC_TIMING


# _ABI_CHECK_FC_CPUTIME()
# ----------------------
#
# Checks whether the Fortran compiler supports CPU_TIME 
# (Fortran 95 and later).
#
AC_DEFUN([_ABI_CHECK_FC_CPUTIME],[
  dnl Init
  fc_has_cputime="no"

  AC_MSG_CHECKING([whether the Fortran compiler accepts cpu_time()])

  dnl Try to compile a call to cpu_time
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([], 
    [[
      real :: second
      call cpu_time(second)

    ]])], [fc_has_cputime="yes"])
  AC_LANG_POP()

  if test "${fc_has_cputime}" = "yes"; then
    AC_DEFINE([HAVE_FC_CPUTIME],1, 
      [Define to 1 if your Fortran compiler supports cpu_time().])
  fi

  AC_MSG_RESULT(${fc_has_cputime})
]) # _ABI_CHECK_FC_CPUTIME


# _ABI_CHECK_FC_GETPID()
# ----------------------
#
# Checks whether process IDs are available from Fortran. (F2003)
#
AC_DEFUN([_ABI_CHECK_FC_GETPID],[
  dnl Init
  fc_has_getpid="no"

  dnl Look for getpid() support
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
          call getpid()
    ]])], [fc_has_getpid="yes"])
  AC_LANG_POP([Fortran])

  dnl Determine whether to use getpid()
  AC_MSG_CHECKING([whether the Fortran compiler accepts getpid()])
  AC_MSG_RESULT([${fc_has_getpid}])
  if test "${fc_has_getpid}" = "yes"; then
    AC_DEFINE([HAVE_FC_GETPID],1,
      [Define to 1 if your Fortran compiler supports getpid().])
  fi
]) # _ABI_CHECK_FC_GETPID


# _ABI_CHECK_FC_ON_THE_FLY_SHAPE()
# --------------------------------
#
# Checks whether process IDs are available from Fortran. (F2003)
#
AC_DEFUN([_ABI_CHECK_FC_ON_THE_FLY_SHAPE],[
  dnl Init
  fc_has_getpid="no"

  dnl Look for getpid() support
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      integer :: x(3,4)
      call test_on_the_fly_shape(x)
      contains
      subroutine test_on_the_fly_shape(x)
        integer, intent(inout) :: x(:,:)
        integer :: y(product(shape(x)))
      end subroutine test_on_the_fly_shape
    ]])], [fc_has_on_the_fly_shape="yes"], [fc_has_on_the_fly_shape="no"])
  AC_LANG_POP([Fortran])

  dnl Determine whether to use getpid()
  AC_MSG_CHECKING([whether the Fortran compiler can shape arrays on-the-fly])
  AC_MSG_RESULT([${fc_has_on_the_fly_shape}])
  if test "${fc_has_on_the_fly_shape}" = "yes"; then
    AC_DEFINE([HAVE_FC_ON_THE_FLY_SHAPE],1,
      [Define to 1 if your Fortran compiler can shape arrays on-the-fly.])
  fi
]) # _ABI_CHECK_FC_ON_THE_FLY_SHAPE



 #############################################################################



# ABI_FC_EXTENSIONS()
# -------------------
#
# Sets the default extensions of Fortran source files and modules,
# whenever possible.
#
AC_DEFUN([ABI_FC_EXTENSIONS],[
  dnl Set Fortran module extension
  AX_F90_MODULE_EXTENSION
  if test "${ax_cv_f90_modext}" != ""; then
    MODEXT="${ax_cv_f90_modext}"
  else
    MODEXT="mod"
    AC_MSG_NOTICE([setting Fortran module extension to ".${MODEXT}"])
  fi

  dnl Change the default Fortran extension for tests
  AC_FC_SRCEXT(F90,[abi_fc_src_ok="yes"],[abi_fc_src_ok="no"])
  if test "${abi_fc_src_ok}" != "yes"; then
    AC_MSG_WARN([Fortran file extension could not be changed])
    AC_MSG_WARN([some advanced Fortran tests may fail])
  fi
]) # ABI_FC_EXTENSIONS



# ABI_FC_FEATURES()
# -----------------
#
# Explores the capabilities of the Fortran compiler.
#
AC_DEFUN([ABI_FC_FEATURES],[
  dnl Explore compiler peculiarities
  _ABI_CHECK_FC_ASYNC
  _ABI_CHECK_FC_BACKTRACE
  _ABI_CHECK_FC_COMMAND_ARGUMENT
  _ABI_CHECK_FC_COMMAND_LINE
  _ABI_CHECK_FC_SYSTEM
  _ABI_CHECK_FC_CONTIGUOUS
  _ABI_CHECK_FC_DTARRAYS
  _ABI_CHECK_FC_IEEE_ARITHMETIC
  _ABI_CHECK_FC_IEEE_EXCEPTIONS
  _ABI_CHECK_FC_IOMSG
  _ABI_CHECK_FC_ISO_C_BINDING
  _ABI_CHECK_FC_EXIT
  _ABI_CHECK_FC_FLUSH
  _ABI_CHECK_FC_FLUSH_
  _ABI_CHECK_FC_GAMMA
  _ABI_CHECK_FC_SHIFTLR
  _ABI_CHECK_FC_GETENV
  _ABI_CHECK_FC_GETPID
  _ABI_CHECK_FC_INT_QUAD
  _ABI_CHECK_FC_ISO_FORTRAN_2008
  _ABI_CHECK_FC_LONG_LINES
  _ABI_CHECK_FC_MACRO_NEWLINE
  _ABI_CHECK_FC_MOVE_ALLOC
  _ABI_CHECK_FC_PRIVATE
  _ABI_CHECK_FC_PROTECTED
  _ABI_CHECK_FC_STREAM_IO
  _ABI_CHECK_FC_CPUTIME
  _ABI_CHECK_FC_TIMING
  _ABI_CHECK_FC_ON_THE_FLY_SHAPE
]) # ABI_FC_FEATURES



# ABI_FC_MOD_CASE()
# -----------------
#
# Checks whether the Fortran compiler creates upper-case or lower-case
# module files.
#
AC_DEFUN([ABI_FC_MOD_CASE],[
  AC_REQUIRE([ABI_FC_EXTENSIONS])

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
]) # ABI_FC_MOD_CASE



# ABI_FC_MOD_INCS(MODULE)
# -----------------------
#
# Checks whether the specified Fortran module is directly available, or
# if we need to add '-I/usr/include' to the compile flags. Returns the
# required includes.
#
AC_DEFUN([ABI_FC_MOD_INCS],[
  AC_MSG_CHECKING([for Fortran module includes])

  if test "${abi_fc_mod_incs_ok}" = "" -o \
          "${abi_fc_mod_incs_ok}" = "unknown"; then

    dnl Init
    fc_mod_incs=""

    dnl Prepare environment
    tmp_saved_FCFLAGS="${FCFLAGS}"
    AC_LANG_PUSH([Fortran])

    dnl Look for module without includes
    AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
      [[
        use $1
      ]])], [abi_fc_mod_incs_ok="none required"], [abi_fc_mod_incs_ok="unknown"])

    dnl Look for module with includes
    if test "${abi_fc_mod_incs_ok}" = "unknown"; then
      FCFLAGS="${FCFLAGS} -I/usr/include"
      AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
        [[
          use $1
        ]])],
        [abi_fc_mod_incs_ok="-I/usr/include"; fc_mod_incs="-I/usr/include"],
        [abi_fc_mod_incs_ok="unknown"])
    fi
    AC_MSG_RESULT([${abi_fc_mod_incs_ok}])

    dnl Restore environment
    AC_LANG_POP([Fortran])
    FCFLAGS="${tmp_saved_FCFLAGS}"

  else

    AC_MSG_RESULT([${abi_fc_mod_incs_ok} (cached)])

  fi

  dnl Substitute variables
  AC_SUBST(fc_mod_incs)
]) # ABI_FC_MOD_INCS



# ABI_PROG_FC()
# -------------
#
# Tries to determine which type of Fortran compiler is installed.
#
AC_DEFUN([ABI_PROG_FC],[
  dnl Init
  abi_fc_vendor="${with_fc_vendor}"
  abi_fc_version="${with_fc_version}"
  tmp_fc_info_file="${abinit_builddir}/config.fc_info.tmp"

  if test "${abi_fc_vendor}" = ""; then
    abi_fc_vendor="unknown"
  fi
  if test "${abi_fc_version}" = ""; then
    abi_fc_version="unknown"
  fi
  abi_fc_wrap="no"

  dnl Determine Fortran compiler type (the order is important)
  AC_MSG_CHECKING([which type of Fortran compiler we have])

  dnl Clear temporary info file
  rm -f "${tmp_fc_info_file}"

  dnl Get rid of that one as early as possible
  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_IBM(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  dnl Should be checked before gfortran because it mimics its behaviour
  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_INTEL(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_GNU(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_ABSOFT(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_NAG(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_PGI(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_LLVM(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  if test "${abi_fc_vendor}" = "unknown"; then
    _ABI_CHECK_FC_ARM(${FC})
  fi
  echo "${fc_info_string}" >>"${tmp_fc_info_file}"

  dnl Fall back to generic when detection fails
  if test "${abi_fc_vendor}" = "unknown"; then
    abi_fc_vendor="generic"
  else
    rm -f "${tmp_fc_info_file}"
  fi

  dnl Normalize Fortran compiler version
  if test "${abi_fc_version}" = "unknown"; then
    abi_fc_version="0.0"
  else
    abi_fc_version=`echo ${abi_fc_version} | cut -d. -f1-2`
  fi

  dnl Display final result
  AC_MSG_RESULT([${abi_fc_vendor} ${abi_fc_version}])

  dnl Schedule compiler info for substitution
  AC_SUBST(abi_fc_vendor)
  AC_SUBST(abi_fc_version)
  AC_SUBST(abi_fc_wrap)
  AC_SUBST(fc_info_string)
]) # ABI_PROG_FC

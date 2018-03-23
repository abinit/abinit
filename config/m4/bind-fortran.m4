# -*- Autoconf -*-
#
# Copyright (C) 2010-2018 ABINIT Group (Damien Caliste, Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fortran bindings
#

# ABI_FC_MODULE_MANGLING()
# ------------------------
#
# Determines the Fortran module name-mangling scheme.
#
AC_DEFUN([ABI_FC_MODULE_MANGLING],[
  dnl Init
  abi_fc_module_ok="no"

  dnl Preserve environment
  ABI_ENV_BACKUP
  abi_saved_LIBS="${LIBS}"
  AC_LANG_PUSH([Fortran])

  dnl Compile a simple module
  AC_COMPILE_IFELSE(dnl
    [[
      module conftest
      contains
      subroutine foobar()
      return
      end subroutine foobar
      end module conftest
    ]],
    [abi_fc_module_ok="yes"; mv conftest.${ac_objext} cfortran_test.${ac_objext}])

  dnl Check that we got the information we need
  if test "${abi_fc_module_ok}" = "no"; then
    AC_MSG_FAILURE([cannot compile a simple Fortran program])
  fi

  dnl Extract information from the object file
  LIBS="cfortran_test.${ac_objext} ${LIBS} $[]_AC_LANG_PREFIX[]LIBS"

  dnl Look at the output of AC_FC_FUNC
  tmp_success="no"
  AC_LANG_PUSH([C])
  AC_FC_FUNC([foobar])
  for tmp_mod in "conftest" "CONFTEST" "foobar" "FOOBAR" "${foobar}" ; do
    for tmp_sub in "foobar" "FOOBAR" "${foobar}" "conftest" "CONFTEST" ; do
      for tmp_begin in "__" "" ; do
        for tmp_middle in "__" "_MOD_" "_MP_" "_mp_" ".in." "." "_" ; do
          tmp_func="${tmp_begin}${tmp_mod}${tmp_middle}${tmp_sub}"
          AC_LINK_IFELSE([AC_LANG_CALL([], [${tmp_func}])],
	    	         [tmp_success="yes"; break 4])
        done
      done
    done
  done
  AC_LANG_POP([C])

  dnl Define first concatenation symbol only if tmp_begin is not empty
  if test "${tmp_begin}" = ""; then
    tmp_bc=""
  else
    tmp_bc="##"
  fi

  dnl Define macro in config.h
  if test "${tmp_success}" = "yes"; then
    AH_TEMPLATE(ABI_FC_MOD,[Fortran module name mangling macro.])
    case "${tmp_mod}x${tmp_sub}" in
      conftestxfoobar)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} mod ## ${tmp_middle} ## sub])
        ;;
      conftestxFOOBAR)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} mod ## ${tmp_middle} ## SUB])
        ;;
      conftestx${foobar})
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} mod ## ${tmp_middle} ## FC_FUNC(sub,SUB)])
        ;;
      CONFTESTxfoobar)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} MOD ## ${tmp_middle} ## sub])
        ;;
      CONFTESTxFOOBAR)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} MOD ## ${tmp_middle} ## SUB])
        ;;
      CONFTESTx${foobar})
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} MOD ## ${tmp_middle} ## FC_FUNC(sub,SUB)])
        ;;
      foobarxconftest)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} sub ## ${tmp_middle} ## mod])
        ;;
      foobarxCONFTEST)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} sub ## ${tmp_middle} ## MOD])
        ;;
      ${foobar}xconftest)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} FC_FUNC(sub,SUB) ## ${tmp_middle} ## mod])
        ;;
      ${foobar}xCONFTEST)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} FC_FUNC(sub,SUB) ## ${tmp_middle} ## MOD])
        ;;
      FOOBARxconftest)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} SUB ## ${tmp_middle} ## mod])
        ;;
      FOOBARxCONFTEST)
        AC_DEFINE_UNQUOTED(ABI_FC_MOD(mod,MOD,sub,SUB), [${tmp_begin} ${tmp_bc} SUB ## ${tmp_middle} ## MOD])
        ;;
     esac
     abi_fc_mod_name="${tmp_begin}module${tmp_middle}subroutine"
  else
     abi_fc_mod_name="unknown"
  fi

  dnl Display final result
  AC_MSG_CHECKING([for the Fortran module name-mangling scheme])
  AC_MSG_RESULT([${abi_fc_mod_name}])

  dnl Restore environment
  AC_LANG_POP([Fortran])
  ABI_ENV_RESTORE
  LIBS="${abi_saved_LIBS}"
  rm -f cfortran_test* conftest*
]) # ABI_FC_MODULE_MANGLING

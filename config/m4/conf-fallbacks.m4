# -*- Autoconf -*-
#
# Copyright (C) 2016 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# Fallbacks support
#



# ABI_FALLBACKS_INIT()
# --------------------
#
# Sets fallbacks parameters according to configure options.
#
AC_DEFUN([ABI_FALLBACKS_INIT],[
  dnl Init
  abi_fbk_config=""
  abi_fbk_enable="no"
  abi_fbk_init="def"
  abi_fbk_list=""
  abi_fbk_prefix=""
  abi_fbk_required=""

  dnl Check fallbacks install prefix and look for the configurator
  if test "${with_fallbacks}" = "yes"; then
    AC_CHECK_PROGS([abi_fbk_config],[abinit-fallbacks-config])
    if test "${abi_fbk_config}" != ""; then
      abi_fbk_enable="yes"
      abi_fbk_init="yon"
    fi
  elif test -d "${with_fallbacks}"; then
    abi_fbk_prefix="${with_fallbacks}"
    abi_fbk_config="${abi_fbk_prefix}/bin/abinit-fallbacks-config"
    if test -x "${abi_fbk_config}"; then
      abi_fbk_enable="yes"
      abi_fbk_init="dir"
    fi
  elif test "${with_fallbacks}" != ""; then
    AC_MSG_ERROR([invalid fallbacks install prefix: '${with_fallbacks}'
                  please use --with-fallbacks without argument or make it
                  point to a readable directory actually containing
                  fallbacks (hint: it should contain an executable
                  program called 'bin/abinit-fallbacks-config')])
  fi

  dnl Decide whether to allow fallbacks
  if test "${abi_fbk_enable}" = "no"; then
    if test "${with_fallbacks}" != "" -a "${with_fallbacks}" != "no"; then
      AC_MSG_ERROR([fallbacks not found
                  please check that --with-fallbacks points to a readable
                  directory or that the abinit-fallbacks-config program
                  is available through your PATH environment variable])
    fi
  fi

  dnl Report final decision
  AC_MSG_CHECKING([whether we can use fallbacks])
  AC_MSG_RESULT([${abi_fbk_enable}])

  dnl Allow reporting of variables
  AC_SUBST(abi_fbk_config)
  AC_SUBST(abi_fbk_enable)
  AC_SUBST(abi_fbk_init)
  AC_SUBST(abi_fbk_list)
  AC_SUBST(abi_fbk_prefix)
  AC_SUBST(abi_fbk_required)
]) # ABI_FALLBACKS_INIT



# ABI_FALLBACKS_VALIDATE(FBK_LIST)
# --------------------------------
#
# Checks whether required fallbacks match those available.
#
AC_DEFUN([ABI_FALLBACKS_VALIDATE],[
  dnl Init
  abi_fbk_list="$1"

  dnl Display which fallbacks we need
  AC_MSG_CHECKING([which fallbacks to look for])
  if test "${abi_fbk_list}" = ""; then
    AC_MSG_RESULT([none])
  else
    AC_MSG_RESULT([${abi_fbk_list}])
  fi

  dnl Abort if fallbacks are disabled
  if test "${abi_fbk_enable}" = "no" -a "${abi_fbk_list}" != ""; then
    ABI_MSG_NOTICE([fallbacks-required],[Abinit Fallbacks required])
    AC_MSG_ERROR([please provide consistent options for optional features
                 or build fallbacks and reconfigure Abinit])
  fi

  dnl Look for each required fallback and set environment variables
  tmp_fbk_list_ok="yes"
  for pkg in ${abi_fbk_list}; do
    tmp_fbk_enabled=`${abi_fbk_config} --enabled ${pkg} 2>/dev/null`
    if test "${tmp_fbk_enabled}" = "yes"; then
      for item in `${abi_fbk_config} --avail ${pkg}`; do
        case "${item}" in
          incs)
            AC_MSG_NOTICE([setting ${item} flags for ${pkg}])
            eval abi_${pkg}_incs=\"`${abi_fbk_config} --incs ${pkg}`\"
            ;;
          libs)
            AC_MSG_NOTICE([setting ${item} flags for ${pkg}])
            eval abi_${pkg}_libs=\"`${abi_fbk_config} --libs ${pkg}`\"
            ;;
        esac
      done
    else
      AC_MSG_WARN([the available Abinit Fallbacks do not provide '${pkg}'])
      tmp_fbk_list_ok="no"
    fi
  done

  dnl Fail if there are missing fallbacks
  if test "${tmp_fbk_list_ok}" != "yes"; then
    ABI_MSG_NOTICE([fallbacks-required],[Abinit Fallbacks required])
    AC_MSG_ERROR([please provide a valid fallback installation prefix or
                 build the missing packages for which there is no fallback
                 (see above warnings)])
  fi
]) # ABI_FALLBACKS_VALIDATE

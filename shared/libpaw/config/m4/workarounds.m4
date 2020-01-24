#
# M4 macros for LibPAW (imported from Abinit)
#
# Copyright (C) 2006-2020 ABINIT Group (Yann Pouillon)



#
# Workarounds for portability issues
#



# AX_PROG_MKDIR_P()
# -----------------
#
# Wrapper for the bugged AC_PROG_MKDIR_P macro.
#
AC_DEFUN([AX_PROG_MKDIR_P],[
  _AC_SRCDIRS([.])
  AC_PROG_MKDIR_P
  ax_tmp_mkdir_p=`echo "${MKDIR_P}" | awk '{print [$]1}'`
  if test "${ax_tmp_mkdir_p}" = "config/gnu/install-sh"; then
    AC_MSG_NOTICE([fixing wrong path to mkdir replacement])
    MKDIR_P="${ac_abs_top_srcdir}/${MKDIR_P}"
  fi
  unset ax_tmp_mkdir_p
]) # AX_PROG_MKDIR_P

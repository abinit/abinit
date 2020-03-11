## Copyright (C) 2020 Yann Pouillon

#
#  Atompaw
#


                    # ------------------------------------ #


#
# Public macros
#

# ABI_CHECK_ATOMPAW_BINS()
# ------------------------
#
# Check whether the Atompaw bins are working.
#
AC_DEFUN([ABI_CHECK_ATOMPAW_BINS],[

  dnl Look for binaries
  AC_CHECK_PROGS([ATOMPAW_BIN],[atompaw])
  AC_CHECK_PROGS([GRAPHATOM_BIN],[graphatom])
  if test "${ATOMPAW_BIN}" != "" -a "${GRAPHATOM_BIN}" != ""; then
      AC_DEFINE([HAVE_ATOMPAW],1,[Define to 1 if you have the AtomPAW library.])
  fi
  

]) # ABI_CHECK_ATOMPAW_BINS

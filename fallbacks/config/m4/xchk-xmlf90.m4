# -*- Autoconf -*-
#
# Copyright (C) 2017 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#



# AFB_CHECK_XMLF90(API_MAJOR_MIN, API_MINOR_MIN)
# ----------------------------------------------
#
# Check whether the specified XMLF90 library is working.
#
AC_DEFUN([AFB_CHECK_XMLF90],[
  dnl Init
  afb_xmlf90_default_libs="-lxmlf90"
  afb_xmlf90_has_incs="unknown"
  afb_xmlf90_has_libs="unknown"
  afb_xmlf90_ext_ok="unknown"

  dnl Prepare environment
  tmp_saved_CPPFLAGS="${CPPFLAGS}"
  tmp_saved_FCFLAGS="${FCFLAGS}"
  tmp_saved_LIBS="${LIBS}"
  CPPFLAGS="${CPPFLAGS} ${afb_xmlf90_incs}"
  FCFLAGS="${FCFLAGS} ${afb_xmlf90_incs}"
  AC_MSG_CHECKING([for XMLF90 libraries to try])
  if test "${afb_xmlf90_libs}" = ""; then
    LIBS="${afb_xmlf90_default_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_xmlf90_default_libs}])
  else
    LIBS="${afb_xmlf90_libs} ${LIBS}"
    AC_MSG_RESULT([${afb_xmlf90_libs}])
  fi

  dnl Look for Fortran modules
  AC_LANG_PUSH([Fortran])
  AC_COMPILE_IFELSE([AC_LANG_PROGRAM([],
    [[
      use xmlf90_cml
      use xmlf90_dom
      use xmlf90_sax
      use xmlf90_wxml
      use xmlf90_xpath
    ]])], [afb_xmlf90_has_incs="yes"], [afb_xmlf90_has_incs="no"])
  AC_LANG_POP([Fortran])

  dnl Check that the library is working
  if test "${afb_xmlf90_has_incs}" = "yes"; then
    AC_MSG_CHECKING([whether XMLF90 is working])
    AC_LANG_PUSH([Fortran])
    AC_RUN_IFELSE([AC_LANG_PROGRAM([],
      [[
        use xmlf90_sax, only: xml_t, open_xmlfile, xml_parse, close_xmlfile

        integer :: ierr
        type(xml_t) :: fxml
        call open_xmlfile("config.xml", fxml, ierr)
        call close_xmlfile(fxml)
      ]])], [afb_xmlf90_has_libs="yes"], [afb_xmlf90_has_libs="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${afb_xmlf90_has_libs}])
  fi

  dnl Final adjustments
  if test "${afb_xmlf90_has_incs}" = "yes" -a \
          "${afb_xmlf90_has_libs}" = "yes"; then
    afb_xmlf90_ext_ok="yes"
  else
    afb_xmlf90_ext_ok="no"
  fi

  dnl Restore environment
  CPPFLAGS="${tmp_saved_CPPFLAGS}"
  FCFLAGS="${tmp_saved_FCFLAGS}"
  LIBS="${tmp_saved_LIBS}"
]) # AFB_CHECK_XMLF90

# ===========================================================================
#        ax_f90_module_case.m4
# ===========================================================================
#
# SYNOPSIS
#
#   AX_F90_MODULE_CASE
#
# DESCRIPTION
#
#   Find Fortran 90 modules file case. The module case is stored
#   in the cached variable ax_f90_mod_case, or "unknown" if the case
#   cannot be found. Two additional cache variables, containing "yes",
#   "no", or "unknown", are created as well: ax_f90_mod_lowercase and
#   ax_f90_mod_uppercase.
#
# LAST MODIFICATION
#
#   2011-08-31
#
# COPYLEFT
#
#   Copyright (C) 2011 Yann Pouillon <yann.pouillon@gmail.com>
#
#   Copying and distribution of this file, with or without modification, are
#   permitted in any medium without royalty provided the copyright notice
#   and this notice are preserved.

AC_DEFUN([AX_F90_MODULE_CASE],[
AC_REQUIRE([AX_F90_MODULE_EXTENSION])
AC_CACHE_CHECK([fortran 90 modules case],
ax_cv_f90_mod_case,
[AC_LANG_PUSH([Fortran])
i=0
while test \( -f tmpdir_$i \) -o \( -d tmpdir_$i \) ; do
    i=`expr $i + 1`
done
mkdir tmpdir_$i
cd tmpdir_$i
AC_COMPILE_IFELSE([
!234567
            module conftest_module
            contains
            subroutine conftest_routine
            write(*,'(a)') 'gotcha!'
            end subroutine conftest_routine
            end module conftest_module
    ],
    [ax_cv_f90_mod_case=`ls | sed -n "s,conftest_module\.$ax_cv_f90_modext,,p"`
      if test x$ax_cv_f90_mod_case = x ; then
          ax_cv_f90_mod_case=`ls | sed -n "s,CONFTEST_MODULE\.$ax_cv_f90_modext,,p"`
          if test x$ax_cv_f90_mod_case = x ; then
              ax_cv_f90_mod_case="unknown"
              ax_cv_f90_mod_lowercase="unknown"
              ax_cv_f90_mod_uppercase="unknown"
          else
            ax_cv_f90_mod_case="upper"
            ax_cv_f90_mod_lowercase="no"
            ax_cv_f90_mod_uppercase="yes"
          fi
      else
          ax_cv_f90_mod_case="lower"
          ax_cv_f90_mod_lowercase="yes"
          ax_cv_f90_mod_uppercase="no"
      fi
    ],
    [ax_cv_f90_mod_case=""; ax_cv_f90_mod_lowercase=""; ax_cv_mod_uppercase=""])
cd ..
rm -fr tmpdir_$i
AC_LANG_POP([Fortran])
])]) # AX_F90_MODULE_CASE

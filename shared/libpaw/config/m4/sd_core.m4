## Copyright (C) 2019 Yann Pouillon

#
# Steredeg core features
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_CORE_INIT], [
  # Init internal Fortran-related variables
  sd_sys_fcflags=""

  # Substitute Steredeg core variables
  AC_SUBST(sd_sys_fcflags)
]) # SD_CORE_INIT

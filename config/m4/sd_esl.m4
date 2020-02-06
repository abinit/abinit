## Copyright (C) 2019 Yann Pouillon

#
# Steredeg support for the ESL
#


                    # ------------------------------------ #


#
# Public macros
#


AC_DEFUN([SD_ESL_INIT], [
  # Declare relevant environment variables
  AC_ARG_VAR([ESL_BUNDLE_PREFIX],
    [Install prefix of the ESL Bundle])

  # Init internal ESL-related variables
  sd_esl_bundle_cppflags=""
  sd_esl_bundle_cflags=""
  sd_esl_bundle_cxxflags=""
  sd_esl_bundle_fcflags=""
  sd_esl_bundle_ldflags=""
  sd_esl_bundle_libdirs=""
  sd_esl_bundle_libs=""
  sd_esl_saved_CPPFLAGS=""
  sd_esl_saved_CFLAGS=""
  sd_esl_saved_CXXFLAGS=""
  sd_esl_saved_FCFLAGS=""
  sd_esl_saved_LDFLAGS=""
  sd_esl_saved_LIBS=""
  sd_esl_saved_bundle_status=""

  # Translate the available ESL configuration
  if test ! -z "${ESL_BUNDLE_PREFIX}"; then
    sd_esl_bundle_cppflags="-I${ESL_BUNDLE_PREFIX}/include"
    sd_esl_bundle_fcflags="-I${ESL_BUNDLE_PREFIX}/include"
    sd_esl_bundle_libdirs="-L${ESL_BUNDLE_PREFIX}/lib"
  fi
]) # SD_ESL_INIT


AC_DEFUN([SD_ESL_ADD_FLAGS], [
  CPPFLAGS="${CPPFLAGS} ${sd_esl_bundle_cppflags}"
  CFLAGS="${CFLAGS} ${sd_esl_bundle_cflags}"
  CXXFLAGS="${CXXFLAGS} ${sd_esl_bundle_cxxflags}"
  FCFLAGS="${FCFLAGS} ${sd_esl_bundle_fcflags}"
  LDFLAGS="${LDFLAGS} ${sd_esl_bundle_ldflags}"
  if test ! -z "$1"; then
    CPPFLAGS="${CPPFLAGS} -I${ESL_BUNDLE_PREFIX}/include/$1"
    FCFLAGS="${FCFLAGS} -I${ESL_BUNDLE_PREFIX}/include/$1"
  fi
  FCFLAGS="${FCFLAGS} -I/usr/include"
]) # SD_ESL_ADD_INCLUDES


AC_DEFUN([SD_ESL_ADD_LIBS], [
  LIBS="${sd_esl_bundle_libdirs} $1 ${sd_esl_bundle_libs} ${LIBS}"
]) # SD_ESL_ADD_LIBS


AC_DEFUN([SD_ESL_RESTORE_FLAGS], [
  CPPFLAGS="${sd_esl_saved_CPPFLAGS}"
  CFLAGS="${sd_esl_saved_CFLAGS}"
  CXXFLAGS="${sd_esl_saved_CXXFLAGS}"
  FCFLAGS="${sd_esl_saved_FCFLAGS}"
  LDFLAGS="${sd_esl_saved_LDFLAGS}"
  LIBS="${sd_esl_saved_LIBS}"
])


AC_DEFUN([SD_ESL_SAVE_FLAGS], [
  sd_esl_saved_CPPFLAGS="${CPPFLAGS}"
  sd_esl_saved_CFLAGS="${CFLAGS}"
  sd_esl_saved_CXXFLAGS="${CXXFLAGS}"
  sd_esl_saved_FCFLAGS="${FCFLAGS}"
  sd_esl_saved_LDFLAGS="${LDFLAGS}"
  sd_esl_saved_LIBS="${LIBS}"
])

# FindWANNIER.cmake
# -----------------
#
# Try to find wannier library (use env variable WANNIER_ROOT as a hint).
#
# Note a recent version of wannier (>= mid 2022) should a package-config file;
# when using a old version of wannier, use find_package(WANNIER) instead.
#
# Recommendation:
# if you use module-environment, just make sure variable WANNIER_ROOT is set to
# top-level directory where wannier was install.
#
# Result Variables
# ----------------
#
# This module defines the following variables::
#
#   WANNIER_FOUND          - True if WANNIER was found
#   WANNIER_INCLUDE_DIRS   - include directories for WANNIER
#   WANNIER_LIBRARIES      - link against this library to use WANNIER
#
# The module will also define two cache variables::
#
#   WANNIER_INCLUDE_DIR    - the WANNIER include directory
#   WANNIER_LIBRARY        - the path to the WANNIER library
#

if(NOT WANNIER_ROOT)
  if (DEFINED ENV{WANNIER_ROOT})
    set(WANNIER_ROOT $ENV{WANNIER_ROOT})
    message(STATUS "Env variable WANNIER_ROOT was defined : $ENV{WANNIER_ROOT}")
  else()
    message(STATUS "Env variable WANNIER_ROOT not defined...")
  endif()
endif()

find_path (WANNIER_INCLUDE_DIR
  NAMES w90_hamiltonian.mod
  HINTS ${WANNIER_ROOT})
mark_as_advanced (WANNIER_INCLUDE_DIR)
message(STATUS "WANNIER_INCLUDE_DIR : ${WANNIER_INCLUDE_DIR}")

find_library (WANNIER_LIBRARY
  NAMES wannier
  HINTS ${WANNIER_ROOT})
mark_as_advanced (WANNIER_LIBRARY)
message(STATUS "WANNIER_LIBRARY : ${WANNIER_LIBRARY}")


# handle the QUIETLY and REQUIRED arguments and set WANNIER_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (WANNIER
  REQUIRED_VARS WANNIER_LIBRARY WANNIER_INCLUDE_DIR)

if (WANNIER_FOUND)
  set (WANNIER_INCLUDE_DIRS ${WANNIER_INCLUDE_DIR})
else (WANNIER_FOUND)
  message("WANNIER nor found")
endif (WANNIER_FOUND)

if(WANNIER_FOUND AND NOT TARGET abinit::wannier)
  add_library(abinit::wannier SHARED IMPORTED)

  set_target_properties(abinit::wannier PROPERTIES
    INTERFACE_INCLUDE_DIRECTORIES "${WANNIER_INCLUDE_DIR}"
    IMPORTED_LOCATION "${WANNIER_LIBRARY}")
endif()

set(ABINIT_WANNIER_FOUND ${WANNIER_FOUND})
set(ABINIT_WANNIER_INCLUDE_DIR ${WANNIER_INCLUDE_DIR})
set(ABINIT_WANNIER_LIBRARY ${WANNIER_LIBRARY})

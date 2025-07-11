#
# 44_manage_kokkos
#

cmake_minimum_required(VERSION 3.18)

# Use variable to control path specific to autotools build integration
set(USE_AUTOTOOLS_BUILD_COMPAT OFF CACHE
  STRING "(Internal) Set on to use paths from autotools buildsystem")

#####################################################################
if(USE_AUTOTOOLS_BUILD_COMPAT)
  project(manage_kokkos_sublib
    LANGUAGES C CXX Fortran
    DESCRIPTION "manage_kokkos_sublib is a CMake project aiming to integrate Kokkos code within ABINIT autotools buildsystem.")

  if(NOT DEFINED ABINIT_SRC_DIR)
    message(FATAL_ERROR "Using CMake within ABINIT autotools buildsystem but ABINIT_SRC_DIR variable wasn't provided!")
  endif()
  # Including path of CMake modules to be used with autotools buildsystem
  list(INSERT CMAKE_MODULE_PATH 0 "${ABINIT_SRC_DIR}/config/cmake")
  include(detect_kokkos)
endif()
#####################################################################

set(44_manage_kokkos_SRC
  abi_common_kokkos.h
  add_array_kokkos.cpp
  assemble_energy_contrib_kokkos.cpp
  opernlc_ylm_allwf_kokkos.cpp
  xg_kokkos.cpp
  xg_kokkos.h
  m_manage_kokkos.F90
  m_xg_kokkos.F90
)

add_library(44_manage_kokkos STATIC
  ${44_manage_kokkos_SRC}
)

#####################################################################
# Set paths for integration with autotools buildsystem
if(USE_AUTOTOOLS_BUILD_COMPAT)

  set_target_properties(44_manage_kokkos
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods)

  target_include_directories(44_manage_kokkos
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/../..
    ${ABINIT_SRC_DIR}/shared/common/src/incs
    ${ABINIT_SRC_DIR}/src/incs
    ${ABINIT_SRC_DIR}/../..
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods
    )

  INSTALL(TARGETS 44_manage_kokkos DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

endif()
#####################################################################

# Set paths properly for CMake buildsystem
if(NOT USE_AUTOTOOLS_BUILD_COMPAT)

  set_target_properties(44_manage_kokkos
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

  target_include_directories(44_manage_kokkos
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/shared/common/src/incs
    ${CMAKE_SOURCE_DIR}/src/incs
    PUBLIC
    ${CMAKE_BINARY_DIR}/modules
    )

endif()

target_compile_definitions(44_manage_kokkos
  PRIVATE
  HAVE_CONFIG_H
)

target_link_libraries(44_manage_kokkos
  PUBLIC
  Kokkos::kokkos
)

add_library(abinit::44_manage_kokkos ALIAS 44_manage_kokkos)

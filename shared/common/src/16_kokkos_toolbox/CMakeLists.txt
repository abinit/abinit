#
# 16_kokkos_toolbox
#

cmake_minimum_required(VERSION 3.18)

# Use variable to control path specific to autotools build integration
set(USE_AUTOTOOLS_BUILD_COMPAT OFF CACHE
  STRING "(Internal) Set on to use paths from autotools buildsystem")

#####################################################################
if(USE_AUTOTOOLS_BUILD_COMPAT)
  project(kokkos_toolbox_sublib
    LANGUAGES C CXX Fortran
    DESCRIPTION "kokkos_toolbox_sublib is a CMake project aiming to integrate Kokkos code within ABINIT autotools buildsystem.")

  if(NOT DEFINED ABINIT_SRC_DIR)
    message(FATAL_ERROR "Using CMake within ABINIT autotools buildsystem but ABINIT_SRC_DIR variable wasn't provided!")
  endif()
  # Including path of CMake modules to be used with autotools buildsystem
  list(INSERT CMAKE_MODULE_PATH 0 "${ABINIT_SRC_DIR}/config/cmake")
  include(detect_kokkos)
endif()
#####################################################################

set(16_kokkos_toolbox_SRC
  kokkos_utils.cpp
  m_kokkos_utils.F90
  )

add_library(16_kokkos_toolbox STATIC
  ${16_kokkos_toolbox_SRC})

#####################################################################
# Set paths for integration with autotools buildsystem
if(USE_AUTOTOOLS_BUILD_COMPAT)

  set_target_properties(16_kokkos_toolbox
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods)

  target_include_directories(16_kokkos_toolbox
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_BINARY_DIR}/../../../..
    ${ABINIT_SRC_DIR}/shared/common/src/incs
    ${ABINIT_SRC_DIR}/src/incs
    ${ABINIT_SRC_DIR}/../../../..
    PUBLIC
    ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods
    )

  INSTALL(TARGETS 16_kokkos_toolbox DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

endif()
#####################################################################

# Set paths properly for CMake buildsystem
if(NOT USE_AUTOTOOLS_BUILD_COMPAT)

  set_target_properties(16_kokkos_toolbox
    PROPERTIES
    Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

  target_include_directories(16_kokkos_toolbox
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/shared/common/src/incs)

  target_include_directories(16_kokkos_toolbox
    PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules;${CMAKE_CURRENT_SOURCE_DIR}>"
    $<INSTALL_INTERFACE:include/libpaw>
    )

endif()

target_compile_definitions(16_kokkos_toolbox
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(16_kokkos_toolbox
  PUBLIC
  Kokkos::kokkos
  )

add_library(abinit::16_kokkos_toolbox ALIAS 16_kokkos_toolbox)

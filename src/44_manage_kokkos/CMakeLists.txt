#
# 44_manage_kokkos
#

# using 3.18 to have nice lapack imported target
cmake_minimum_required(VERSION 3.18)

include(autotools_integration.cmake)

set(44_manage_kokkos_SRC
  abi_common_kokkos.h
  add_array_kokkos.cpp
  assemble_energy_contrib_kokkos.cpp
  kokkos_utils.cpp
  opernlc_ylm_allwf_kokkos.cpp
  xg_kokkos.cpp
  xg_kokkos.h
  m_kokkos_utils.F90
  m_manage_kokkos.F90
  m_xg_kokkos.F90
  )

add_library(44_manage_kokkos STATIC
  ${44_manage_kokkos_SRC})

set_target_properties(44_manage_kokkos
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods)

target_include_directories(44_manage_kokkos
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/../../shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/../../src/incs
  ${CMAKE_SOURCE_DIR}/../..
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/../../src/mods
  )

target_compile_definitions(44_manage_kokkos
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(44_manage_kokkos
  PUBLIC
  Kokkos::kokkos
  )

INSTALL(TARGETS 44_manage_kokkos DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(abinit::44_manage_kokkos ALIAS 44_manage_kokkos)

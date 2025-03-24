#
# 62_cg_noabirule
#
add_library(62_cg_noabirule STATIC
  m_cgtk.F90
  m_frskerker1.F90
  m_frskerker2.F90
  )

set_target_properties(62_cg_noabirule
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(62_cg_noabirule
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(62_cg_noabirule
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(62_cg_noabirule
  PUBLIC
  abinit::52_fft_mpi_noabirule
  abinit::54_spacepar
  abinit::56_recipspace
  )

if(OpenMP_FOUND)
  target_link_libraries(62_cg_noabirule
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::62_cg_noabirule ALIAS 62_cg_noabirule)

#
# 56_recipspace
#
add_library(56_recipspace STATIC
  m_initylmg.F90
  m_ab7_kpoints.F90
  m_bz_mesh.F90
  m_double_grid.F90
  m_gsphere.F90
  m_getshell.F90
  m_kg.F90
  m_kpts.F90
  m_lgroup.F90
  m_nesting.F90
  m_skw.F90
  )

set_target_properties(56_recipspace
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(56_recipspace
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(56_recipspace
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(56_recipspace
  PUBLIC
  abinit::44_abitypes_defs
  abinit::54_spacepar
  )

if(OpenMP_FOUND)
  target_link_libraries(56_recipspace
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::56_recipspace ALIAS 56_recipspace)

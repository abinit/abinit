#
# 43_ptgroups
#
add_library(43_ptgroups STATIC
  m_defs_ptgroups.F90
  m_ptgroups.F90
  m_ptg_C1.F90
  m_ptg_Ci.F90
  m_ptg_C2.F90
  m_ptg_Cs.F90
  m_ptg_C2h.F90
  m_ptg_D2.F90
  m_ptg_C2v.F90
  m_ptg_D2h.F90
  m_ptg_C4.F90
  m_ptg_S4.F90
  m_ptg_C4h.F90
  m_ptg_D4.F90
  m_ptg_C4v.F90
  m_ptg_D2d.F90
  m_ptg_D4h.F90
  m_ptg_C3.F90
  m_ptg_C3i.F90
  m_ptg_D3.F90
  m_ptg_C3v.F90
  m_ptg_D3d.F90
  m_ptg_C6.F90
  m_ptg_C3h.F90
  m_ptg_C6h.F90
  m_ptg_D6.F90
  m_ptg_C6v.F90
  m_ptg_D3h.F90
  m_ptg_D6h.F90
  m_ptg_T.F90
  m_ptg_Th.F90
  m_ptg_O.F90
  m_ptg_Td.F90
  m_ptg_Oh.F90
  )

set_target_properties(43_ptgroups
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(43_ptgroups
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(43_ptgroups
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(43_ptgroups
  PUBLIC
  28_numeric_noabirule
  32_util)

add_library(abinit::43_ptgroups ALIAS 43_ptgroups)

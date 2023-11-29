#
# 80_tdep
#
add_library(80_tdep STATIC
  m_tdep_abitypes.F90
  m_tdep_constraints.F90
  m_tdep_latt.F90
  m_tdep_phdos.F90
  m_tdep_phi2.F90
  m_tdep_phi3.F90
  m_tdep_phi4.F90
  m_tdep_qpt.F90
  m_tdep_readwrite.F90
  m_tdep_shell.F90
  m_tdep_sym.F90
  m_tdep_utils.F90
  )

set_target_properties(80_tdep
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(80_tdep
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(80_tdep
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(80_tdep
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::44_abitypes_defs
  abinit::45_geomoptim
  abinit::51_manage_mpi
  abinit::72_response
  abinit::77_ddb
  )

add_library(abinit::80_tdep ALIAS 80_tdep)

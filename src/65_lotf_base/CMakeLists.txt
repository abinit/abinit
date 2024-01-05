#
# 65_lotf_base
#
add_library(65_lotf_base STATIC
  m_bond_lotf.F90
  defs_param_lotf.F90
  m_glue_lotf.F90
  m_pbc_lotf.F90
  m_tools_lotf.F90
  m_work_var_lotf.F90
  )

set_target_properties(65_lotf_base
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(65_lotf_base
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(65_lotf_base
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(65_lotf_base
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  )

add_library(abinit::65_lotf_base ALIAS 65_lotf_base)

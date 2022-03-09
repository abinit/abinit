#
# 68_lotf
#
add_library(68_lotf STATIC
  m_eval_lotf.F90
  m_lotf.F90
  m_pred_lotf.F90
  )

set_target_properties(68_lotf
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(68_lotf
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(68_lotf
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(68_lotf
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::45_geomoptim
  abinit::65_lotf_base
  )

add_library(abinit::68_lotf ALIAS 68_lotf)

#
# 68_recursion
#
add_library(68_recursion STATIC
  defs_rectypes.F90
  m_hidecudarec.F90
  m_pretty_rec.F90
  m_rec.F90
  m_rec_tools.F90
  m_vtorhorec.F90
  )

set_target_properties(68_recursion
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(68_recursion
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(68_recursion
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(68_recursion
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::44_abitypes_defs
  abinit::65_paw
  )

add_library(abinit::68_recursion ALIAS 68_recursion)

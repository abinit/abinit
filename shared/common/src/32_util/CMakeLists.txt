#
# 32_util
#
add_library(32_util STATIC
  m_atomdata.F90
  m_dummy.F90
  m_exp_mat.F90
  m_matrix.F90
  m_pptools.F90
  m_special_funcs.F90
  m_symtk.F90
  m_pstat.F90
  )

set_target_properties(32_util
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(32_util
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(32_util
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(32_util
  PUBLIC
  abinit::17_libtetra_ext
  abinit::28_numeric_noabirule
  )

add_library(abinit::32_util ALIAS 32_util)

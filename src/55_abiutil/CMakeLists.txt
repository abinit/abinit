#
# 55_abiutil
#
add_library(55_abiutil STATIC
  m_dens.F90
  m_results_out.F90
  m_scf_history.F90
  m_builtin_tests.F90
  m_timana.F90
  )

set_target_properties(55_abiutil
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(55_abiutil
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(55_abiutil
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(55_abiutil
  PUBLIC
  abinit::53_ffts
  )

add_library(abinit::55_abiutil ALIAS 55_abiutil)

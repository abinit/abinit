#
# 77_suscep
#
add_library(77_suscep STATIC
  m_suscep_stat.F90
  )

set_target_properties(77_suscep
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(77_suscep
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(77_suscep
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(77_suscep
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::44_abitypes_defs
  abinit::51_manage_mpi
  abinit::56_recipspace
  abinit::65_paw
  )

add_library(abinit::77_suscep ALIAS 77_suscep)

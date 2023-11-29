#
# 66_vdwxc
#
add_library(66_vdwxc STATIC
  m_evdw_wannier.F90
  m_vdw_dftd2.F90
  m_vdw_dftd3.F90
  m_vdw_dftd3_data.F90
  )

set_target_properties(66_vdwxc
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(66_vdwxc
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(66_vdwxc
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(66_vdwxc
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::32_util
  abinit::41_geometry
  )

add_library(abinit::66_vdwxc ALIAS 66_vdwxc)

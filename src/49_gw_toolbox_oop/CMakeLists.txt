#
# 49_gw_toolbox_oop
#
add_library(49_gw_toolbox_oop STATIC
  m_bs_defs.F90
  m_gwdefs.F90
  m_melemts.F90
  m_spectra.F90
  )

set_target_properties(49_gw_toolbox_oop
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(49_gw_toolbox_oop
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(49_gw_toolbox_oop
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(49_gw_toolbox_oop
  PUBLIC
  abinit::44_abitypes_defs)

add_library(abinit::49_gw_toolbox_oop ALIAS 49_gw_toolbox_oop)

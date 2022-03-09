#
# 27_toolbox_oop
#
add_library(27_toolbox_oop STATIC
  m_copy.F90
  m_iterators.F90
  m_mpiotk.F90
  m_nctk.F90
  )

set_target_properties(27_toolbox_oop
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(27_toolbox_oop
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(27_toolbox_oop
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(27_toolbox_oop
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(27_toolbox_oop
  PUBLIC
  abinit::10_defs
  abinit::14_hidewrite
  abinit::16_hideleave
  abinit::17_yaml_out
  )

add_library(abinit::27_toolbox_oop ALIAS 27_toolbox_oop)

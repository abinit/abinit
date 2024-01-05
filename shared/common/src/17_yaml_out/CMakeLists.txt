#
# 17_yaml_out
#
add_library(17_yaml_out STATIC
  c_pair_list.c
  m_type_pair_list.F90
  m_pair_list.F90
  m_stream_string.F90
  m_yaml.F90
)

set_target_properties(17_yaml_out
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(17_yaml_out
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(17_yaml_out
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(17_yaml_out
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(17_yaml_out
  PUBLIC
  abinit::16_hideleave
  )

add_library(abinit::17_yaml_out ALIAS 17_yaml_out)

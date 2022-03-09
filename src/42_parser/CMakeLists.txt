#
# 42_parser
#
add_library(42_parser STATIC
  m_parser.F90
  )

set_target_properties(42_parser
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(42_parser
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(42_parser
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(42_parser
  PUBLIC
  abinit::16_hideleave
  abinit::32_util
  abinit::41_geometry
  )

add_library(abinit::42_parser ALIAS 42_parser)

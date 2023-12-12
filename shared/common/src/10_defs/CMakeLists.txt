#
# 10_defs
#
add_library(10_defs STATIC
  defs_basis.F90
  defs_datatypes.F90)

set_target_properties(10_defs
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(10_defs
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(10_defs
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(10_defs
  PRIVATE
  HAVE_CONFIG_H)

add_library(abinit::10_defs ALIAS 10_defs)

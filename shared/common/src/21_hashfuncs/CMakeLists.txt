#
# 21_hashfuncs
#
add_library(21_hashfuncs STATIC
  m_hash_md5.F90
  )

set_target_properties(21_hashfuncs
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(21_hashfuncs
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(21_hashfuncs
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(21_hashfuncs
  PUBLIC
  abinit::16_hideleave
  )

add_library(abinit::21_hashfuncs ALIAS 21_hashfuncs)

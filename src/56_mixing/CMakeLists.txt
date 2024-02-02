#
# 56_mixing
#
add_library(56_mixing STATIC
  m_ab7_mixing.F90
  )

set_target_properties(56_mixing
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(56_mixing
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(56_mixing
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(56_mixing
  PUBLIC
  abinit::44_abitypes_defs
  )

if(OpenMP_FOUND)
  target_link_libraries(56_mixing
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::56_mixing ALIAS 56_mixing)

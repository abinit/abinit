#
# 54_spacepar
#
add_library(54_spacepar STATIC
  m_multipoles.F90
  m_spacepar.F90
  )

set_target_properties(54_spacepar
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(54_spacepar
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(54_spacepar
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(54_spacepar
  PUBLIC
  abinit::53_ffts
  )

if(OpenMP_FOUND)
  target_link_libraries(54_spacepar
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::54_spacepar ALIAS 54_spacepar)

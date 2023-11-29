#
# 44_abitools
#
add_library(44_abitools STATIC
  m_berrytk.F90
  m_cgtools.F90
  m_dynmat.F90
  m_integrals.F90
  m_macroave.F90
  m_wfutils.F90
  )

set_target_properties(44_abitools
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(44_abitools
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(44_abitools
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(44_abitools
  PUBLIC
  28_numeric_noabirule
  32_util
  39_libpaw
  41_geometry)

if(OpenMP_FOUND)
  target_link_libraries(44_abitools
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::44_abitools ALIAS 44_abitools)

#
# 42_nlstrain
#
add_library(42_nlstrain STATIC
  m_contistr01.F90
  m_contistr03.F90
  m_contistr12.F90
  m_contstr21.F90
  m_contstr23.F90
  m_contstr25.F90
  m_contstr25a.F90
  m_contstr26.F90
  )

set_target_properties(42_nlstrain
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(42_nlstrain
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(42_nlstrain
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(42_nlstrain
  PUBLIC
  abinit::14_hidewrite
  )

add_library(abinit::42_nlstrain ALIAS 42_nlstrain)

#
# 41_geometry
#
add_library(41_geometry STATIC
  m_ab7_symmetry.F90
  m_crystal.F90
  m_ewald.F90
  m_geometry.F90
  m_gtermcutoff.F90
  m_spgdata.F90
  m_spgbuilder.F90
  m_supercell.F90
  m_symlist.F90
  m_symsg.F90
  m_symfind.F90
  )

set_target_properties(41_geometry
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(41_geometry
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(41_geometry
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(41_geometry
  PUBLIC
  abinit::28_numeric_noabirule
  abinit::32_util
  abinit::39_libpaw
  )

add_library(abinit::41_geometry ALIAS 41_geometry)

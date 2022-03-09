#
# 29_kpoints
#
add_library(29_kpoints STATIC
  m_hashtable.F90
  m_htetra.F90
  m_octree.F90
  m_symkpt.F90
  )

set_target_properties(29_kpoints
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(29_kpoints
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(29_kpoints
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(29_kpoints
  PUBLIC
  abinit::17_libtetra_ext
  abinit::28_numeric_noabirule
  )

add_library(abinit::29_kpoints ALIAS 29_kpoints)

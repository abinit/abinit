#
# 44_abitypes_defs
#
add_library(44_abitypes_defs STATIC
  defs_abitypes.F90
  m_args_gs.F90
  m_distribfft.F90
  m_dtset.F90
  m_dtfil.F90
  m_efield.F90
  m_efmas_defs.F90
  m_energies.F90
  m_results_gs.F90
  m_results_img.F90
  m_results_respfn.F90
  )

set_target_properties(44_abitypes_defs
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(44_abitypes_defs
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(44_abitypes_defs
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(44_abitypes_defs
  PUBLIC
  29_kpoints
  32_util
  41_geometry
  42_parser)

add_library(abinit::44_abitypes_defs ALIAS 44_abitypes_defs)

#
# 77_ddb
#
add_library(77_ddb STATIC
  m_epweights.F90
  defs_elphon.F90
  m_ddb_diel.F90
  m_ddb_elast.F90
  m_ddb_interpolate.F90
  m_eliashberg_1d.F90
  m_elphon.F90
  m_ddb_flexo.F90
  m_harmonic_thermo.F90
  m_ddb_internalstr.F90
  m_a2ftr.F90
  m_anaddb_dataset.F90
  m_atprj.F90
  m_band2eps_dataset.F90
  m_gruneisen.F90
  m_ifc.F90
  m_lwf.F90
  m_phonons.F90
  m_primcell_ddb_info.F90
  m_sortph.F90
  m_ddb_piezo.F90
  m_raman.F90
  m_iogkk.F90
  m_relaxpol.F90
  m_thmeig.F90
  )

set_target_properties(77_ddb
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(77_ddb
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(77_ddb
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(77_ddb
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::17_libtetra_ext
  abinit::28_numeric_noabirule
  abinit::32_util
  abinit::41_geometry
  abinit::42_parser
  abinit::44_abitools
  abinit::71_wannier
  abinit::72_response
  )

add_library(abinit::77_ddb ALIAS 77_ddb)

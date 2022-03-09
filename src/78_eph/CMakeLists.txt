#
# 78_eph
#
add_library(78_eph STATIC
  m_cumulant.F90
  m_eph_double_grid.F90
  m_ephwg.F90
  m_frohlichmodel.F90
  m_gstore.F90
  m_gkk.F90
  m_migdal_eliashberg.F90
  m_phgamma.F90
  m_phpi.F90
  m_rta.F90
  m_sigmaph.F90
  m_berry_curvature.F90
  )

set_target_properties(78_eph
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(78_eph
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(78_eph
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(78_eph
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::17_libtetra_ext
  abinit::29_kpoints
  abinit::41_geometry
  abinit::44_abitypes_defs
  abinit::61_occeig
  abinit::62_iowfdenpot
  abinit::77_ddb
  )

add_library(abinit::78_eph ALIAS 78_eph)

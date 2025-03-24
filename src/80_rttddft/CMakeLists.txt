#
# 80_rttddft
#
add_library(80_rttddft STATIC
  m_rttddft.F90
  m_rttddft_exponential.F90
  m_rttddft_output.F90
  m_rttddft_propagate.F90
  m_rttddft_propagators.F90
  m_rttddft_properties.F90
  m_rttddft_tdks.F90
  )

set_target_properties(80_rttddft
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(80_rttddft
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(80_rttddft
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(80_rttddft
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::66_wfs
  abinit::67_common
  abinit::79_seqpar_mpi
  )

add_library(abinit::80_rttddft ALIAS 80_rttddft)

#
# 71_wannier
#
add_library(71_wannier STATIC
  defs_wannier90.F90
  m_mlwfovlp.F90
  m_mlwfovlp_qp.F90
  m_scdm_math.F90
  m_wannier_builder.F90
  m_wann_netcdf.F90
  )

set_target_properties(71_wannier
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(71_wannier
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(71_wannier
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(71_wannier
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::32_util
  abinit::52_fft_mpi_noabirule
  abinit::56_io_mpi
  abinit::56_recipspace
  abinit::61_occeig
  abinit::65_paw
  abinit::70_gw
  )

if(ABINIT_WANNIER_FOUND)
  target_link_libraries(71_wannier
    PUBLIC
    abinit::wannier)
endif()

add_library(abinit::71_wannier ALIAS 71_wannier)

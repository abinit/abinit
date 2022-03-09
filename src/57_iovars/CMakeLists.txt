#
# 57_iovars
#
add_library(57_iovars STATIC
  m_chkinp.F90
  m_ingeo.F90
  m_inkpts.F90
  m_invars1.F90
  m_invars2.F90
  m_memeval.F90
  m_mpi_setup.F90
  m_outvars.F90
  m_outvar_a_h.F90
  m_outvar_i_n.F90
  m_outvar_o_z.F90
  m_outxml.F90
  m_out_acknowl.F90
  )

set_target_properties(57_iovars
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(57_iovars
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(57_iovars
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(57_iovars
  PUBLIC
  abinit::44_abitypes_defs
  abinit::45_geomoptim
  abinit::49_gw_toolbox_oop
  abinit::52_fft_mpi_noabirule
  abinit::55_abiutil
  abinit::56_io_mpi
  abinit::56_recipspace
  abinit::56_xc
  )

add_library(abinit::57_iovars ALIAS 57_iovars)

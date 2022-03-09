#
# 71_bse
#
add_library(71_bse STATIC
  m_exc_build.F90
  m_exc_analyze.F90
  m_bse_io.F90
  m_exc_diago.F90
  m_bseinterp.F90
  m_exc_itdiago.F90
  m_exc_spectra.F90
  m_haydock.F90
  m_haydock_io.F90
  m_hexc.F90
  )

set_target_properties(71_bse
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(71_bse
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(71_bse
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(71_bse
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::49_gw_toolbox_oop
  abinit::56_io_mpi
  abinit::56_recipspace
  abinit::61_occeig
  abinit::67_common
  abinit::69_wfdesc
  abinit::70_gw
  )

if(OpenMP_FOUND)
  target_link_libraries(71_bse
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::71_bse ALIAS 71_bse)

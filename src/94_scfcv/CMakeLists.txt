#
# 94_scfcv
#
add_library(94_scfcv STATIC
  m_afterscfloop.F90
  m_elpolariz.F90
  m_gwr.F90
  m_scfcv_core.F90
  m_scfcv.F90
  m_outscfcv.F90
 )

set_target_properties(94_scfcv
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(94_scfcv
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(94_scfcv
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(94_scfcv
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::51_manage_mpi
  abinit::56_io_mpi
  abinit::68_dmft
  abinit::68_recursion
  abinit::68_rsprc
  abinit::70_gw
  abinit::71_wannier
  abinit::79_seqpar_mpi
  )

if(OpenMP_FOUND)
  target_link_libraries(94_scfcv
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::94_scfcv ALIAS 94_scfcv)

#
# 62_wvl_wfs
#
add_library(62_wvl_wfs STATIC
  m_wvl_psi.F90
  m_wvl_rwwf.F90
 )

set_target_properties(62_wvl_wfs
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(62_wvl_wfs
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(62_wvl_wfs
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(62_wvl_wfs
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::51_manage_mpi
  abinit::56_io_mpi
  )

add_library(abinit::62_wvl_wfs ALIAS 62_wvl_wfs)

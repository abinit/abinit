#
# 62_iowfdenpot
#
add_library(62_iowfdenpot STATIC
  m_io_screening.F90
  m_wfk.F90
  m_outqmc.F90
  m_outwant.F90
  )

set_target_properties(62_iowfdenpot
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(62_iowfdenpot
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(62_iowfdenpot
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(62_iowfdenpot
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::44_abitypes_defs
  abinit::49_gw_toolbox_oop
  abinit::56_io_mpi
  abinit::61_occeig
  abinit::62_cg_noabirule
  )

add_library(abinit::62_iowfdenpot ALIAS 62_iowfdenpot)

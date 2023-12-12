#
# 68_rsprc
#
add_library(68_rsprc STATIC
  mod_prc_memory.F90
  m_newrho.F90
  m_newvtr.F90
  m_prcref.F90
  )

set_target_properties(68_rsprc
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(68_rsprc
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(68_rsprc
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(68_rsprc
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_wvl_wrappers
  abinit::56_xc
  abinit::62_cg_noabirule
  abinit::67_common
  )

add_library(abinit::68_rsprc ALIAS 68_rsprc)

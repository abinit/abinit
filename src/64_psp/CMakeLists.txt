#
# 64_psp
#
add_library(64_psp STATIC
  m_atm2fft.F90
  m_psptk.F90
  m_psps.F90
  m_pspini.F90
  m_psp1.F90
  m_psp_hgh.F90
  m_psp5.F90
  m_psp6.F90
  m_psp8.F90
  m_psp9.F90
  m_upf2abinit.F90
  )

set_target_properties(64_psp
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(64_psp
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(64_psp
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(64_psp
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::53_ffts
  abinit::57_iovars
  abinit::57_iopsp_parser
  )

add_library(abinit::64_psp ALIAS 64_psp)

#
# 72_response
#
add_library(72_response STATIC
  m_dfpt_fef.F90
  m_dfptnl_pert.F90
  m_dfpt_cgwf.F90
  m_dfpt_elt.F90
  m_dfpt_lwwf.F90
  m_dfpt_mkrho.F90
  m_dfpt_mkvxcstr.F90
  m_dfpt_nstwf.F90
  m_dfpt_rhotov.F90
  m_dfpt_scfcv.F90
  m_dfpt_vtorho.F90
  m_dfpt_vtowfk.F90
  m_d2frnl.F90
  m_ddb.F90
  m_ddb_hdr.F90
  m_ddk.F90
  m_dvdb.F90
  m_efmas.F90
  m_eig2d.F90
  m_orbmag.F90
  m_paral_pert.F90
  m_rf2.F90
  m_rf2_init.F90
  m_strain.F90
  )

set_target_properties(72_response
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(72_response
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(72_response
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(72_response
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::44_abitools
  abinit::44_abitypes_defs
  abinit::56_io_mpi
  abinit::56_xc
  abinit::62_iowfdenpot
  abinit::64_psp
  abinit::66_nonlocal
  abinit::66_wfs
  abinit::68_rsprc
  abinit::69_wfdesc
  )

if(OpenMP_FOUND)
  target_link_libraries(72_response
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::72_response ALIAS 72_response)

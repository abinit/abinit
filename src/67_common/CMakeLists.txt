#
# 67_common
#
add_library(67_common STATIC
  m_berryphase.F90
  m_berryphase_new.F90
  m_conducti.F90
  m_dft_energy.F90
  m_extraprho.F90
  m_ephtk.F90
  m_forces.F90
  m_forstr.F90
  m_cgwf.F90
  m_cgwf_cprj.F90
  m_epjdos.F90
  m_ksdiago.F90
  mkcore_wvl.F90
  m_mklocl_realspace.F90
  m_mklocl.F90
  m_mkrho.F90
  m_nucprop.F90
  m_ioarr.F90
  m_iowf.F90
  m_eprenorms.F90
  m_sigtk.F90
  m_odamix.F90
  m_optic_tools.F90
  m_optics_vloc.F90
  m_rot_cg.F90
  m_rhotov.F90
  m_common.F90
  m_positron.F90
  m_setvtr.F90
  m_spin_current.F90
  m_stress.F90
  m_vtorhotf.F90
  m_wvl_rho.F90
  )

set_target_properties(67_common
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(67_common
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(67_common
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(67_common
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::44_abitools
  abinit::52_fft_mpi_noabirule
  abinit::56_io_mpi
  abinit::56_mixing
  abinit::62_iowfdenpot
  abinit::62_wvl_wfs
  abinit::66_nonlocal
  abinit::66_vdwxc
  abinit::66_wfs
  )

if(OpenMP_FOUND)
  target_link_libraries(67_common
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if (ABINIT_YAKL_WANTED)
  target_link_libraries(67_common
    PUBLIC
    abinit::yakl)
  target_include_directories(67_common
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::67_common ALIAS 67_common)

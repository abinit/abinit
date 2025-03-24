#
# 70_gw
#
add_library(70_gw STATIC
  m_calc_ucrpa.F90
  m_sigc.F90
  m_sigx.F90
  m_vhxc_me.F90
  m_chi0.F90
  m_cohsex.F90
  m_classify_bands.F90
  m_gwls_hamiltonian.F90
  m_gwls_wf.F90
  m_gwls_valenceWavefunctions.F90
  m_gwls_lineqsolver.F90
  m_gwls_polarisability.F90
  m_gwls_TimingLog.F90
  m_gwls_model_polarisability.F90
  m_gwls_QR_factorization.F90
  m_gwls_GWlanczos.F90
  m_gwls_GenerateEpsilon.F90
  m_gwls_LanczosBasis.F90
  m_gwls_DielectricArray.F90
  m_gwls_GWanalyticPart.F90
  m_gwls_LanczosResolvents.F90
  m_gwls_Projected_AT.F90
  m_gwls_Projected_BT.F90
  m_gwls_ComputePoles.F90
  m_gwls_ComputeCorrelationEnergy.F90
  m_gwls_utility.F90
  m_gwrdm.F90
  m_chi0tk.F90
  m_dyson_solver.F90
  m_fft_prof.F90
  m_model_screening.F90
  m_plowannier.F90
  m_ppmodel.F90
  m_prep_calc_ucrpa.F90
  m_qparticles.F90
  m_read_plowannier.F90
  m_screening.F90
  m_screen.F90
  m_sigma.F90
  )

set_target_properties(70_gw
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(70_gw
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(70_gw
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(70_gw
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::41_geometry
  abinit::44_abitypes_defs
  abinit::62_iowfdenpot
  abinit::66_wfs
  abinit::69_wfdesc
  )

if(OpenMP_FOUND)
  target_link_libraries(70_gw
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::70_gw ALIAS 70_gw)

#
# 95_drive
#
add_library(95_drive STATIC
  m_argparse.F90
  m_bader.F90
  m_cut3d.F90
  m_bethe_salpeter.F90
  m_dfptnl_loop.F90
  m_dfpt_looppert.F90
  m_driver.F90
  m_eph_driver.F90
  m_fold2block.F90
  m_gstate.F90
  m_gstateimg.F90
  m_gwls_sternheimer.F90
  m_gwr_driver.F90
  m_longwave.F90
  m_mover.F90
  m_mover_effpot.F90
  m_multibinit_driver.F90
  m_nonlinear.F90
  m_precpred_1geo.F90
  m_respfn_driver.F90
  m_rttddft_driver.F90
  m_pead_nl_loop.F90
  m_screening_driver.F90
  m_sigma_driver.F90
  m_unittests.F90
  m_wfk_analyze.F90
  )

set_target_properties(95_drive
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(95_drive
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(95_drive
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(95_drive
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::45_geomoptim
  abinit::49_gw_toolbox_oop
  abinit::56_io_mpi
  abinit::62_iowfdenpot
  abinit::67_common
  abinit::70_gw
  abinit::71_bse
  abinit::72_response
  abinit::77_ddb
  abinit::78_effpot
  abinit::78_eph
  abinit::79_seqpar_mpi
  abinit::80_rttddft
  abinit::80_tdep
  abinit::94_scfcv
  )

if (ABINIT_YAKL_WANTED)
  target_link_libraries(95_drive
    PUBLIC
    abinit::yakl)
  target_include_directories(95_drive
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::95_drive ALIAS 95_drive)

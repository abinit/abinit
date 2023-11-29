#
# 45_geomoptim
#
add_library(45_geomoptim STATIC
  m_predtk.F90
  m_abihist.F90
  m_abimover.F90
  m_abiimages.F90
  m_bfgs.F90
  m_ga.F90
  m_hmc.F90
  m_ipi.F90
  m_lbfgs.F90
  m_mep.F90
  m_m1geo.F90
  m_pimd.F90
  m_use_ga.F90
  m_pimd_langevin.F90
  m_pimd_nosehoover.F90
  m_pred_bfgs.F90
  m_pred_delocint.F90
  m_pred_diisrelax.F90
  m_pred_fire.F90
  m_pred_hmc.F90
  m_pred_isokinetic.F90
  m_pred_isothermal.F90
  m_pred_langevin.F90
  m_pred_moldyn.F90
  m_pred_nose.F90
  m_pred_simple.F90
  m_pred_srkna14.F90
  m_pred_steepdesc.F90
  m_pred_velverlet.F90
  m_pred_verlet.F90
  m_predict_neb.F90
  m_predict_pimd.F90
  m_predict_steepest.F90
  m_predict_string.F90
  m_xfpack.F90
  )

set_target_properties(45_geomoptim
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(45_geomoptim
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(45_geomoptim
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(45_geomoptim
  PUBLIC
  abinit::44_abitypes_defs)

add_library(abinit::45_geomoptim ALIAS 45_geomoptim)

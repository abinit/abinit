#
# 65_paw
#
add_library(65_paw STATIC
  m_fourier_interpol.F90
  m_paw_atomorb.F90
  m_paw_correlations.F90
  m_paw_denpot.F90
  m_paw_dfpt.F90
  m_paw_dfptnl.F90
  m_paw_dmft.F90
  m_paw_efield.F90
  m_paw_hr.F90
  m_paw_init.F90
  m_paw_lmn.F90
  m_paw_mkaewf.F90
  m_paw_mkrho.F90
  m_paw_nhat.F90
  m_paw_nmr.F90
  m_paw_occupancies.F90
  m_paw_optics.F90
  m_paw_overlap.F90
  m_paw_pwaves_lmn.F90
  m_paw_slater.F90
  m_paw_sym.F90
  m_paw_tools.F90
  m_paw_uj.F90
  m_pawfgr.F90
  m_pawpwij.F90
  )

set_target_properties(65_paw
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(65_paw
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(65_paw
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(65_paw
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::56_io_mpi
  abinit::56_xc
  abinit::56_recipspace
  abinit::61_occeig
  abinit::62_ctqmc
  abinit::64_psp
  )

if(OpenMP_FOUND)
  target_link_libraries(65_paw
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::65_paw ALIAS 65_paw)

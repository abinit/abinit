#
# 65_paw
#
add_library(65_paw STATIC
  m_fourier_interpol.F90
  m_paw_correlations.F90
  m_paw_denpot.F90
  m_paw_dfpt.F90
  m_paw_dfptnl.F90
  m_paw_dmft.F90
  m_paw_efield.F90
  m_paw_exactDC.F90
  m_paw_hr.F90
  m_paw_init.F90
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
  m_paw_yukawa.F90
  m_pawfgr.F90
  m_pawpwij.F90
  )

set_target_properties(65_paw
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(65_paw
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_BINARY_DIR}/modules
  )

target_compile_definitions(65_paw
  PRIVATE
  HAVE_CONFIG_H)

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(65_paw
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

target_link_libraries(65_paw
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::28_numeric_noabirule
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

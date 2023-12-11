#
# 62_ctqmc
#
add_library(62_ctqmc STATIC
  m_BathOperator.F90
  m_BathOperatoroffdiag.F90
  m_data4entropyDMFT.F90
  m_Ctqmc.F90
  m_Ctqmcoffdiag.F90
  m_CtqmcInterface.F90
  m_CtqmcoffdiagInterface.F90
  m_FFTHyb.F90
  m_GreenHyb.F90
  m_GreenHyboffdiag.F90
  m_ImpurityOperator.F90
  m_ListCdagC.F90
  m_MapHyb.F90
  m_MatrixHyb.F90
  m_OurRng.F90
  m_Stat.F90
  m_Vector.F90
  m_VectorInt.F90
  m_global.F90
  )

set_target_properties(62_ctqmc
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(62_ctqmc
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(62_ctqmc
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(62_ctqmc
  PUBLIC
  abinit::10_defs
  abinit::11_memory_mpi
  abinit::16_hideleave
  abinit::28_numeric_noabirule
  )

add_library(abinit::62_ctqmc ALIAS 62_ctqmc)

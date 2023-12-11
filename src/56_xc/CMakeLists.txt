#
# 56_xc
#
add_library(56_xc STATIC
  m_electronpositron.F90
  m_kxc.F90
  m_xcdata.F90
  m_xc_vdw.F90
  m_dfpt_mkvxc.F90
  m_gammapositron.F90
  m_mkcore.F90
  m_rhotoxc.F90
  m_xctk.F90
  m_xchybrid.F90
  )

set_target_properties(56_xc
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(56_xc
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(56_xc
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(56_xc
  PUBLIC
  abinit::44_abitypes_defs
  abinit::51_manage_mpi
  abinit::53_ffts
  abinit::54_spacepar
  )

if(OpenMP_FOUND)
  target_link_libraries(56_xc
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::56_xc ALIAS 56_xc)

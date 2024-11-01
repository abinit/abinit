#
# 62_poisson
#
add_library(62_poisson STATIC
  m_barevcoul.F90
  m_cutoff_sphere.F90
  m_cutoff_slab.F90
  m_cutoff_cylinder.F90
  m_jellium.F90
  m_psolver.F90
  m_qplusg.F90
  m_vcoul.F90
  )

set_target_properties(62_poisson
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(62_poisson
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(62_poisson
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(62_poisson
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::28_numeric_noabirule
  abinit::41_geometry
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::49_gw_toolbox_oop
  abinit::56_recipspace
  abinit::53_ffts
  )

add_library(abinit::62_poisson ALIAS 62_poisson)

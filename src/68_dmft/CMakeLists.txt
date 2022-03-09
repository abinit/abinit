#
# 68_dmft
#
add_library(68_dmft STATIC
  m_datafordmft.F90
  m_dmft.F90
  m_energy.F90
  m_entropyDMFT.F90
  m_forctqmc.F90
  m_green.F90
  m_hu.F90
  m_hubbard_one.F90
  m_dftu_self.F90
  m_matlu.F90
  m_oper.F90
  m_self.F90
  )

set_target_properties(68_dmft
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(68_dmft
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(68_dmft
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(68_dmft
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::28_numeric_noabirule
  abinit::39_libpaw
  abinit::41_geometry
  abinit::44_abitypes_defs
  abinit::65_paw
  )

if(DO_BUILD_67_PYTHON_INVOCATION_EXT)
  target_link_libraries(68_dmft
    PUBLIC
    abinit::67_python_invocation_ext)
endif()

add_library(abinit::68_dmft ALIAS 68_dmft)

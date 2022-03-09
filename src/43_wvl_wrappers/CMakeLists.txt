#
# 43_wvl_wrappers
#
add_library(43_wvl_wrappers STATIC
  defs_wvltypes.F90
  defs_PSolver.F90
  m_abi2big.F90
  m_paw2wvl.F90
  m_wvl_denspot.F90
  m_wvl_descr_psp.F90
  m_wvl_projectors.F90
  m_wvl_wfs.F90
  )

set_target_properties(43_wvl_wrappers
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(43_wvl_wrappers
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(43_wvl_wrappers
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(43_wvl_wrappers
  PUBLIC
  abinit::41_geometry)

if(ABINIT_BIGDFT_FOUND)
  target_link_libraries(43_wvl_wrappers
    PUBLIC
    abinit::bigdft)
endif()

add_library(abinit::43_wvl_wrappers ALIAS 43_wvl_wrappers)

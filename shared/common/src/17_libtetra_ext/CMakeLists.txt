#
# 17_libtetra_ext
#
add_library(17_libtetra_ext STATIC
  libtetra.h
  libtetrabz.F90
  libtetrabz_common.F90
  libtetrabz_dbldelta_mod.F90
  libtetrabz_dblstep_mod.F90
  libtetrabz_dos_mod.F90
  libtetrabz_fermigr_mod.F90
  libtetrabz_occ_mod.F90
  libtetrabz_polcmplx_mod.F90
  libtetrabz_polstat_mod.F90
  m_simtet.F90
  m_krank.F90
  m_tetrahedron.F90
)

set_target_properties(17_libtetra_ext
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(17_libtetra_ext
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(17_libtetra_ext
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(17_libtetra_ext
  PUBLIC
  abinit::14_hidewrite
  abinit::16_hideleave
  )

if(OpenMP_FOUND)
  target_link_libraries(17_libtetra_ext
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::17_libtetra_ext ALIAS 17_libtetra_ext)

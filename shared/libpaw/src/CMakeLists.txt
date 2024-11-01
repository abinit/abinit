#
# 39_libpaw
#
add_library(39_libpaw STATIC
  libpaw_libxc.c
  m_libpaw_mpi.F90
  m_libpaw_tools.F90
  m_libpaw_libxc.F90
  m_paral_atom.F90
  m_paw_numeric.F90
  m_paw_sphharm.F90
  m_paw_io.F90
  m_pawang.F90
  m_pawrad.F90
  m_pawtab.F90
  m_paw_an.F90
  m_paw_ij.F90
  m_pawfgrtab.F90
  m_paw_finegrid.F90
  m_paw_onsite.F90
  m_pawcprj.F90
  m_pawrhoij.F90
  m_pawdij.F90
  m_pawxc.F90
  m_paw_atom.F90
  m_paw_gaussfit.F90
  m_pawxmlps.F90
  m_pawpsp.F90
  )

set_target_properties(39_libpaw
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(39_libpaw
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs)

target_include_directories(39_libpaw
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(39_libpaw
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(39_libpaw
  PUBLIC
  abinit::33_xc_lowlevel
  )

if(OpenMP_FOUND)
  target_link_libraries(39_libpaw
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::39_libpaw ALIAS 39_libpaw)

#
# installable libpaw
#
set_target_properties(39_libpaw
  PROPERTIES OUTPUT_NAME paw)

# change export name, so that the user projet can import target abinit::paw
set_target_properties(39_libpaw
  PROPERTIES EXPORT_NAME paw)

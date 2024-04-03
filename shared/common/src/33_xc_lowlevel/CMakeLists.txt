#
# 33_xc_lowlevel
#
add_library(33_xc_lowlevel STATIC
  libxc_tools.c
  m_libxc_functionals.F90
  m_drivexc.F90
  m_xc_noncoll.F90
  m_xchcth.F90
  m_xciit.F90
  m_xcpbe.F90
  m_xcpositron.F90
  m_xclda.F90
  )

set_target_properties(33_xc_lowlevel
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(33_xc_lowlevel
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(33_xc_lowlevel
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(33_xc_lowlevel
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(33_xc_lowlevel
  PUBLIC
  abinit::28_numeric_noabirule
  )

if(OpenMP_FOUND)
  target_link_libraries(33_xc_lowlevel
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(ABINIT_LIBXC_FOUND)
  target_link_libraries(33_xc_lowlevel
    PUBLIC
    abinit::libxc)
endif()

add_library(abinit::33_xc_lowlevel ALIAS 33_xc_lowlevel)

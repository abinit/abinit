#
# 69_wfdesc
#
add_library(69_wfdesc STATIC
  m_wfd_optic.F90
  m_io_kss.F90
  m_oscillators.F90
  m_wfd.F90
  )

set_target_properties(69_wfdesc
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(69_wfdesc
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(69_wfdesc
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(69_wfdesc
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::27_toolbox_oop
  abinit::49_gw_toolbox_oop
  abinit::53_ffts
  abinit::56_io_mpi
  abinit::62_iowfdenpot
  abinit::65_paw
  abinit::66_nonlocal
  abinit::67_common
  )

if(OpenMP_FOUND)
  target_link_libraries(69_wfdesc
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::69_wfdesc ALIAS 69_wfdesc)

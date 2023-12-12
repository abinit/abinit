#
# 61_occeig
#
add_library(61_occeig STATIC
  m_occ.F90
  m_ebands.F90
  m_fstab.F90
  m_extfpmd.F90
  )

set_target_properties(61_occeig
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(61_occeig
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(61_occeig
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(61_occeig
  PUBLIC
  abinit::44_abitypes_defs
  abinit::56_recipspace
  abinit::56_io_mpi
  )

if(OpenMP_FOUND)
  target_link_libraries(61_occeig
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::61_occeig ALIAS 61_occeig)

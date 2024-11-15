#
# 56_io_mpi
#
add_library(56_io_mpi STATIC
  m_abi_etsf.F90
  m_exit.F90
  m_hdr.F90
  m_rwwf.F90
  )

set_target_properties(56_io_mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(56_io_mpi
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(56_io_mpi
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(56_io_mpi
  PUBLIC
  abinit::43_wvl_wrappers
  abinit::44_abitypes_defs
  abinit::51_manage_mpi
  )

add_library(abinit::56_io_mpi ALIAS 56_io_mpi)

#
# 11_memory_mpi
#
add_library(11_memory_mpi STATIC
  m_profiling_abi.F90)

set_target_properties(11_memory_mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(11_memory_mpi
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(11_memory_mpi
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(11_memory_mpi
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(11_memory_mpi
  PRIVATE
  abinit::02_clib
  abinit::10_defs
  MPI::MPI_Fortran)

add_library(abinit::11_memory_mpi ALIAS 11_memory_mpi)

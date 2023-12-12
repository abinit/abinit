#
# 16_hideleave
#
add_library(16_hideleave STATIC
  m_xieee.F90
  m_errors.F90
)

set_target_properties(16_hideleave
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(16_hideleave
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${NetCDF_MOD_INC}
  )

target_include_directories(16_hideleave
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(16_hideleave
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(16_hideleave
  PUBLIC
  abinit::10_defs
  abinit::11_memory_mpi
  abinit::12_hide_mpi
  abinit::14_hidewrite
  abinit::netcdf-fortran
  abinit::netcdf
  MPI::MPI_Fortran
  )

add_library(abinit::16_hideleave ALIAS 16_hideleave)

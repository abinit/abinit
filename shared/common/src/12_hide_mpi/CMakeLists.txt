#
# 12_hide_mpi
#
set(12_hide_mpi_SRC
  m_xmpi.F90
  m_xomp.F90
  xmpi_allgather.finc
  xmpi_iallgather.finc
  xmpi_allgatherv.finc
  xmpi_gather.finc
  xmpi_gatherv.finc
  xmpi_bcast.finc
  xmpi_min.finc
  xmpi_sum_master.finc
  xmpi_exch.finc
  xmpi_recv.finc
  xmpi_irecv.finc
  xmpi_sum.finc
  xmpi_ibcast.finc
  xmpi_isum.finc
  xmpi_alltoall.finc
  xmpi_ialltoall.finc
  xmpi_alltoallv.finc
  xmpi_ialltoallv.finc
  xmpi_max.finc
  xmpi_send.finc
  xmpi_isend.finc
  xmpi_scatterv.finc
  xmpi_land_lor.finc)

add_library(12_hide_mpi STATIC
  ${12_hide_mpi_SRC})

set_target_properties(12_hide_mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(12_hide_mpi
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(12_hide_mpi
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(12_hide_mpi
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(12_hide_mpi
  PUBLIC
  abinit::02_clib
  abinit::10_defs
  abinit::11_memory_mpi
  MPI::MPI_Fortran)

if(OpenMP_FOUND)
  target_link_libraries(12_hide_mpi
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::12_hide_mpi ALIAS 12_hide_mpi)

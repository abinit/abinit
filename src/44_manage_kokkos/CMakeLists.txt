#
# 44_manage_kokkos
#
set(44_manage_kokkos_SRC
  add_array_kokkos.cpp
  assemble_energy_contrib_kokkos.cpp
  opernlc_ylm_allwf_kokkos.cpp
  xg_kokkos.cpp
  xg_kokkos.h
  m_xg_kokkos.F90
  )

add_library(44_manage_kokkos STATIC
  ${44_manage_kokkos_SRC})

set_target_properties(44_manage_kokkos
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(44_manage_kokkos
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(44_manage_kokkos
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(44_manage_kokkos
  PUBLIC
  abinit::12_hide_mpi
  abinit::17_gpu_toolbox
  )

if (ABINIT_YAKL_WANTED)
  target_link_libraries(44_manage_kokkos
    PUBLIC
    abinit::yakl)
  target_include_directories(44_manage_kokkos
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

target_link_libraries(44_manage_kokkos
  PUBLIC
  Kokkos::kokkos
  )

add_library(abinit::44_manage_kokkos ALIAS 44_manage_kokkos)

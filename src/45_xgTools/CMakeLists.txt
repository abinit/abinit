#
# 45_xgTools
#
set(66_wfs_SRC
  m_xg.F90
  m_xgScalapack.F90
  m_xgTransposer.F90
  )

if (ABINIT_ENABLE_GPU_CUDA AND ABINIT_KOKKOS_WANTED)
  LIST(APPEND 66_wfs_SRC
    m_xg_kokkos.cpp
    m_xg_kokkos.h
    m_xg_kokkos.F90)
endif()

add_library(45_xgTools STATIC
  ${66_wfs_SRC})

set_target_properties(45_xgTools
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(45_xgTools
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(45_xgTools
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(45_xgTools
  PUBLIC
  abinit::44_abitypes_defs)

if(OpenMP_FOUND)
  target_link_libraries(45_xgTools
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(45_xgTools
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

if (ABINIT_YAKL_WANTED)
  target_link_libraries(45_xgTools
    PUBLIC
    abinit::yakl)
  target_include_directories(45_xgTools
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::45_xgTools ALIAS 45_xgTools)

#
# 46_ghc_omp
#
set(46_ghc_omp_SRC
  m_ompgpu_utils.F90
  m_ompgpu_fourwf.F90
  )

add_library(46_ghc_omp STATIC
  ${46_ghc_omp_SRC})

set_target_properties(46_ghc_omp
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(46_ghc_omp
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(46_ghc_omp
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(46_ghc_omp
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  )

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(46_ghc_omp
    PUBLIC
    abinit::46_manage_cuda
    )
endif()

if(OpenMP_FOUND)
  target_link_libraries(46_ghc_omp
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(46_ghc_omp
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

add_library(abinit::46_ghc_omp ALIAS 46_ghc_omp)

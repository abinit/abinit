#
# 53_ffts
#
add_library(53_ffts STATIC
  m_fft.F90
  )

set_target_properties(53_ffts
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(53_ffts
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(53_ffts
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(53_ffts
  PUBLIC
  abinit::46_ghc_omp
  abinit::52_fft_mpi_noabirule
  )

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(53_ffts
    PUBLIC
    abinit::46_manage_cuda
    )
endif()

if(OpenMP_FOUND)
  target_link_libraries(53_ffts
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::53_ffts ALIAS 53_ffts)

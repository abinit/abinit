#
# 17_gpu_toolbox
#
add_library(17_gpu_toolbox STATIC
  abi_gpu_header_common.h
  dev_spec.cpp
  gpu_fft.cpp
  gpu_linalg.cpp
  m_gpu_detect.F90
  m_gpu_toolbox.F90
  m_initcuda.F90
  m_nvtx.F90
  timing.cpp
)

set_target_properties(17_gpu_toolbox
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(17_gpu_toolbox
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs)

target_include_directories(17_gpu_toolbox
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/modules;${CMAKE_CURRENT_SOURCE_DIR}>"
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(17_gpu_toolbox
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(17_gpu_toolbox
  PUBLIC
  abinit::10_defs
  abinit::14_hidewrite
  )

if (ABINIT_KOKKOS_WANTED)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    abinit::16_kokkos_toolbox
  )
endif()

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    CUDA::cudart
    CUDA::cufft
    CUDA::cublas
    CUDA::cusolver
    )
endif()

if(ABINIT_ENABLE_GPU_HIP)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    hip::amdhip64
    hip::hipfft
    roc::hipblas
    roc::rocblas
    roc::rocfft
    roc::hipsolver
    )
endif()

if(HAVE_GPU_MARKERS)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    abinit::gpu_markers
    )
endif()

if (ABINIT_YAKL_WANTED)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    abinit::yakl)
  target_include_directories(17_gpu_toolbox
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::17_gpu_toolbox ALIAS 17_gpu_toolbox)

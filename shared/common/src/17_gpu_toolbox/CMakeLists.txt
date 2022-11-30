#
# 17_gpu_toolbox
#
add_library(17_gpu_toolbox STATIC
  abi_gpu_header_common.h
  dev_spec.cu
  gpu_fft.cu
  gpu_linalg.cu
  m_gpu_detect.F90
  m_gpu_toolbox.F90
  m_initcuda.F90
  m_nvtx.F90
  timing_cuda.cu
)

set_target_properties(17_gpu_toolbox
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(17_gpu_toolbox
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs)

target_include_directories(17_gpu_toolbox
  PUBLIC
  "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules;${CMAKE_CURRENT_SOURCE_DIR}>"
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

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(17_gpu_toolbox
    PUBLIC
    CUDA::cublas
    CUDA::cusolver
    CUDA::nvToolsExt
    )
endif()

add_library(abinit::17_gpu_toolbox ALIAS 17_gpu_toolbox)

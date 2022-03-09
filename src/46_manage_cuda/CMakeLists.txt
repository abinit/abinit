#
# 46_manage_cuda
#
add_library(46_manage_cuda STATIC
  m_manage_cuda.F90
  m_alloc_hamilt_gpu.F90
  m_nvtx_data.F90
  abi_gpu_header.h
  cuda_rec.cu
  cuda_rec_cut.cu
  density.cu
  gpu_apply_invovl_inner.cu
  gpu_apply_local_potential.cu
  gpu_compute_nl_hamiltonian.cu
  gpu_compute_nl_projections.cu
  gpu_density_accumulation.cu
  gpu_fourwf.cu
  gpu_gemm_nonlop.cu
  gpu_gemm_nonlop.h
  gpu_mkkpg.cu
  gpu_nonlop.cu
  gpu_sphere.cu
  gpu_fftbox.cu
  prt_utils_rec.cu
  rec_dens_calc.cu
  rec_kernels.cu
  recursion_bth.cu
  recursion_cut_bth.cu
  recursion_cut_no_bth.cu
  recursion_no_bth.cu
  )

set_target_properties(46_manage_cuda
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(46_manage_cuda
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_CURRENT_SOURCE_DIR}
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  ${CMAKE_SOURCE_DIR}/src/incs
  )

target_compile_definitions(46_manage_cuda
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(46_manage_cuda
  PUBLIC
  abinit::17_gpu_toolbox
  abinit::44_abitypes_defs
  )

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(46_manage_cuda
    PUBLIC
    CUDA::cufft
    CUDA::cublas
    )
endif()

add_library(abinit::46_manage_cuda ALIAS 46_manage_cuda)

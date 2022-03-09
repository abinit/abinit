#
# 28_numeric_noabirule
#
add_library(28_numeric_noabirule STATIC
  m_lib_four.F90
  m_abi_linalg.F90
  m_array.F90
  m_bessel.F90
  m_cplxtools.F90
  m_bspline.F90
  m_elpa.F90
  m_gaussian_quadrature.F90
  m_hide_blas.F90
  m_hide_lapack.F90
  m_lebedev.F90
  m_levenberg_marquardt.F90
  m_linalg_interfaces.F90
  m_numeric_tools.F90
  m_per_cond.F90
  m_random_xoroshiro128plus.F90
  m_random_zbq.F90
  m_slk.F90
  m_sort.F90
  m_splines.F90
  )

set_target_properties(28_numeric_noabirule
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(28_numeric_noabirule
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  )

target_include_directories(28_numeric_noabirule
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
  $<INSTALL_INTERFACE:include/libpaw>
  )

target_compile_definitions(28_numeric_noabirule
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(28_numeric_noabirule
  PUBLIC
  abinit::18_timing
  abinit::27_toolbox_oop
  BLAS::BLAS
  LAPACK::LAPACK
  )

if(ABINIT_SCALAPACK_ENABLED)
  target_link_libraries(28_numeric_noabirule
    PUBLIC
    abinit::scalapack
    )
endif()

if(OpenMP_FOUND)
  target_link_libraries(28_numeric_noabirule
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(28_numeric_noabirule
    PUBLIC
    abinit::17_gpu_toolbox
    )
endif()

add_library(abinit::28_numeric_noabirule ALIAS 28_numeric_noabirule)

#
# 44_abitools
#
add_library(44_abitools STATIC
  m_berrytk.F90
  m_cgtools.F90
  m_dynmat.F90
  m_integrals.F90
  m_macroave.F90
  m_nvtx_data.F90
  m_wfutils.F90
)

set_target_properties(44_abitools
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/modules)

target_include_directories(44_abitools
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  PUBLIC
  ${CMAKE_BINARY_DIR}/modules
  )

target_compile_definitions(44_abitools
  PRIVATE
  HAVE_CONFIG_H)

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(44_abitools
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

target_link_libraries(44_abitools
  PUBLIC
  28_numeric_noabirule
  32_util
  39_libpaw
  41_geometry
)

if(ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP OR ABINIT_ENABLE_GPU_MARKERS)
  target_link_libraries(44_abitools
    PUBLIC
    abinit::17_gpu_toolbox
    )
endif()

if(OpenMP_FOUND)
  target_link_libraries(44_abitools
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

add_library(abinit::44_abitools ALIAS 44_abitools)

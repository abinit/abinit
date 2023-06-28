#
# 66_nonlocal
#

set(66_nonlocal_SRC
  m_cgprj.F90
  m_contract.F90
  m_fock.F90
  m_gemm_nonlop_ompgpu.F90
  m_gemm_nonlop_gpu.F90
  m_gemm_nonlop.F90
  m_hamiltonian.F90
  m_metstr.F90
  m_mkffkg.F90
  m_mkffnl.F90
  m_nonlop.F90
  m_nonlop_pl.F90
  m_nonlop_test.F90
  m_nonlop_ylm.F90
  m_opernl.F90
  m_vkbr.F90
  m_opernla_ylm.F90
  m_opernla_ylm_mv.F90
  m_opernlb_ylm.F90
  m_opernlb_ylm_mv.F90
  m_opernlc_ylm.F90
  m_opernlc_ylm_allwf_cpu.F90
  m_opernlc_ylm_ompgpu.F90
  m_opernld_ylm.F90
  )

add_library(66_nonlocal STATIC
  ${66_nonlocal_SRC})


set_target_properties(66_nonlocal
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(66_nonlocal
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(66_nonlocal
  PRIVATE
  HAVE_CONFIG_H)

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(66_nonlocal
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

target_link_libraries(66_nonlocal
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::42_nlstrain
  abinit::62_poisson
  abinit::65_paw
  )

if(ABINIT_ENABLE_GPU_CUDA AND ABINIT_KOKKOS_WANTED)
  target_link_libraries(66_nonlocal
    PUBLIC
    abinit::44_manage_kokkos
  )
endif()

if(OpenMP_FOUND)
  target_link_libraries(66_nonlocal
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(ABINIT_ENABLE_GPU_CUDA)
  target_link_libraries(66_nonlocal
    PUBLIC
    abinit::46_manage_cuda
    )
endif()

if (ABINIT_YAKL_WANTED)
  target_link_libraries(66_nonlocal
    PUBLIC
    abinit::yakl)
  target_include_directories(66_nonlocal
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::66_nonlocal ALIAS 66_nonlocal)

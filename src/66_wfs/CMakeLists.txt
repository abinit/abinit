#
# 66_wfs
#
set(66_wfs_SRC
  m_esymm.F90
  m_bandfft_kpt.F90
  m_invovl.F90
  m_chebfi.F90
  m_cgcprj.F90
  m_fock_getghc.F90
  m_getchc.F90
  m_getghc.F90
  m_getghc_ompgpu.F90
  m_getgh1c.F90
  m_getgh2c.F90
  m_prep_kgb.F90
  m_rayleigh_ritz.F90
  m_rmm_diis.F90
  )

add_library(66_wfs STATIC
  ${66_wfs_SRC})

set_target_properties(66_wfs
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(66_wfs
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(66_wfs
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(66_wfs
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::43_ptgroups
  abinit::44_abitools
  abinit::44_abitypes_defs
  abinit::46_ghc_omp
  abinit::56_recipspace
  abinit::62_cg_noabirule
  abinit::66_nonlocal
  )

if(ABINIT_ENABLE_GPU_CUDA AND ABINIT_KOKKOS_WANTED)
  target_link_libraries(66_wfs
    PUBLIC
    abinit::44_manage_kokkos
  )
endif()

if(OpenMP_FOUND)
  target_link_libraries(66_wfs
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(HAVE_OPENMP_OFFLOAD)
  target_compile_options(66_wfs
    PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:${OPENMP_OFFLOAD_FLAGS}>)
endif()

add_library(abinit::66_wfs ALIAS 66_wfs)

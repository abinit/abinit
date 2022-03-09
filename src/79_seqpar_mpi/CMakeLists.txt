#
# 79_seqpar_mpi
#
add_library(79_seqpar_mpi STATIC
  m_chebfiwf.F90
  m_dfpt_lw.F90
  m_lobpcgwf_old.F90
  m_inwffil.F90
  m_lobpcgwf.F90
  m_lobpcg.F90
  m_tddft.F90
  m_vtowfk.F90
  m_vtorho.F90
  m_wvl_wfsinp.F90
  )

set_target_properties(79_seqpar_mpi
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(79_seqpar_mpi
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(79_seqpar_mpi
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(79_seqpar_mpi
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::44_abitypes_defs
  abinit::48_diago
  abinit::62_wvl_wfs
  abinit::66_wfs
  abinit::67_common
  abinit::68_dmft
  abinit::70_gw
  abinit::72_response
  abinit::77_suscep
  )

if(OpenMP_FOUND)
  target_link_libraries(79_seqpar_mpi
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if (ABINIT_YAKL_WANTED)
  target_link_libraries(79_seqpar_mpi
    PUBLIC
    abinit::yakl)
  target_include_directories(79_seqpar_mpi
    PUBLIC
    ${YAKL_BINARY_DIR})
endif()

add_library(abinit::79_seqpar_mpi ALIAS 79_seqpar_mpi)

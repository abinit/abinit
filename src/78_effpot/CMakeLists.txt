#
# 78_effpot
#
add_library(78_effpot STATIC
  m_compute_anharmonics.F90
  effpot_xml.c
  m_generate_training_set.F90
  m_init10.F90
  m_abstract_mover.F90
  m_abstract_potential.F90
  m_anharmonics_terms.F90
  m_dynarray.F90
  m_effective_potential.F90
  m_effective_potential_file.F90
  m_effpot_mpi.F90
  m_fit_data.F90
  m_fit_polynomial_coeff.F90
  m_harmonics_terms.F90
  m_hashtable_strval.F90
  m_lattice_berendsen_NPT_mover.F90
  m_lattice_berendsen_NVT_mover.F90
  m_lattice_dummy_mover.F90
  m_lattice_effpot.F90
  m_lattice_harmonic_potential.F90
  m_lattice_io_netcdf.F90
  m_lattice_langevin_mover.F90
  m_lattice_lwf_map.F90
  m_lattice_lwf_mover.F90
  m_lattice_mover.F90
  m_lattice_ncfile.F90
  m_lattice_verlet_mover.F90
  m_linked_list.F90
  m_lwf_berendsen_mover.F90
  m_lwf_dummy_mover.F90
  m_lwf_hist.F90
  m_lwf_mc_mover.F90
  m_lwf_mover.F90
  m_lwf_ncfile.F90
  m_lwf_observables.F90
  m_lwf_primitive.F90
  m_lwf_potential.F90
  m_mathfuncs.F90
  m_mergesort.F90
  m_mpi_scheduler.F90
  m_multibinit_cell.F90
  m_multibinit_dataset.F90
  m_multibinit_io_xml.F90
  m_multibinit_main.F90
  m_multibinit_manager.F90
  m_multibinit_unittest.F90
  m_observables.F90
  m_opt_effpot.F90
  m_polynomial.F90
  m_polynomial_coeff.F90
  m_polynomial_conf.F90
  m_polynomial_potential.F90
  m_polynomial_term.F90
  m_scup_dataset.F90
  m_potential_list.F90
  m_primitive_potential.F90
  m_primitive_potential_list.F90
  m_spin_hist.F90
  m_spin_mc_mover.F90
  m_spin_mover.F90
  m_spin_ncfile.F90
  m_spin_observables.F90
  m_spin_potential.F90
  m_spin_primitive_potential.F90
  m_spin_reciprocal.F90
  m_slc_dynamics.F90
  m_slc_potential.F90
  m_slc_primitive_potential.F90
  m_spmat_base.F90
  m_spmat_convert.F90
  m_spmat_coo.F90
  m_spmat_csr.F90
  m_spmat_dense.F90
  m_spmat_lco.F90
  m_spmat_lil.F90
  m_spmat_spvec.F90
  m_spmat_NDCOO.F90
  m_supercell_maker.F90
  m_tightbinding.F90
  m_twobody_interaction.F90
  )

set_target_properties(78_effpot
  PROPERTIES
  Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/modules)

target_include_directories(78_effpot
  PRIVATE
  ${CMAKE_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/shared/common/src/incs
  ${CMAKE_SOURCE_DIR}/src/incs
  PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}/modules
  )

target_compile_definitions(78_effpot
  PRIVATE
  HAVE_CONFIG_H)

target_link_libraries(78_effpot
  PUBLIC
  abinit::10_defs
  abinit::16_hideleave
  abinit::28_numeric_noabirule
  abinit::42_parser
  abinit::56_recipspace
  abinit::77_ddb
  )

if(OpenMP_FOUND)
  target_link_libraries(78_effpot
    PUBLIC
    OpenMP::OpenMP_Fortran)
endif()

if(ABINIT_XML2_FOUND)
  target_link_libraries(78_effpot
    PUBLIC
    abinit::xml2)
endif()

if(ABINIT_NETCDF_FOUND)
  target_link_libraries(78_effpot
    PUBLIC
    abinit::netcdf-fortran
    abinit::netcdf)
endif()


add_library(abinit::78_effpot ALIAS 78_effpot)

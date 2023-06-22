#
# 98_main
#

#
# Default cmake macro to build and link an executable using abinit libraries
#
macro(abi_build_exe exe_name)

  add_executable(${exe_name} ${exe_name}.F90)

  target_include_directories(${exe_name}
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/shared/common/src/incs
    ${CMAKE_SOURCE_DIR}/src/incs
    )

  target_compile_definitions(${exe_name}
    PRIVATE
    HAVE_CONFIG_H)

  target_link_libraries(${exe_name}
    PUBLIC
    abinit::95_drive
    )

  if (OpenMP_FOUND)
    target_link_libraries(${exe_name}
      PUBLIC
      OpenMP::OpenMP_Fortran)
  endif (OpenMP_FOUND)

  if(HAVE_OPENMP_OFFLOAD)
    target_link_options(${exe_name}
      PRIVATE
      ${OPENMP_OFFLOAD_FLAGS})
  endif()

  if (Libxc_FOUND)
    target_link_libraries(${exe_name}
      PUBLIC
      abinit::libxc)
  endif (Libxc_FOUND)

  if (ABINIT_NETCDF_FOUND)
    target_link_libraries(${exe_name}
      PUBLIC
      abinit::netcdf-fortran
      abinit::netcdf)
  endif (ABINIT_NETCDF_FOUND)

  set_target_properties(${exe_name}
    PROPERTIES
    LINKER_LANGUAGE "Fortran")

endmacro(abi_build_exe)

abi_build_exe(abinit)
abi_build_exe(abitk)
abi_build_exe(aim)
abi_build_exe(anaddb)
abi_build_exe(atdep)
abi_build_exe(band2eps)
abi_build_exe(conducti)
abi_build_exe(cut3d)
abi_build_exe(dummy_tests)
abi_build_exe(fftprof)
abi_build_exe(fold2Bloch)
abi_build_exe(ioprof)
abi_build_exe(lapackprof)
abi_build_exe(macroave)
abi_build_exe(mrgddb)
abi_build_exe(mrgdv)
abi_build_exe(mrggkk)
abi_build_exe(mrgscr)
abi_build_exe(multibinit)
abi_build_exe(optic)
abi_build_exe(testtransposer)
abi_build_exe(vdw_kernelgen)

# abinit library
add_library(abinit_lib STATIC
  $<TARGET_OBJECTS:02_clib>
  $<TARGET_OBJECTS:10_defs>
  $<TARGET_OBJECTS:11_memory_mpi>
  $<TARGET_OBJECTS:12_hide_mpi>
  $<TARGET_OBJECTS:14_hidewrite>
  $<TARGET_OBJECTS:16_hideleave>
  $<TARGET_OBJECTS:17_libtetra_ext>
  $<TARGET_OBJECTS:17_yaml_out>
  $<TARGET_OBJECTS:18_timing>
  $<TARGET_OBJECTS:21_hashfuncs>
  $<TARGET_OBJECTS:27_toolbox_oop>
  $<TARGET_OBJECTS:28_numeric_noabirule>
  $<TARGET_OBJECTS:29_kpoints>
  $<TARGET_OBJECTS:32_util>
  $<TARGET_OBJECTS:33_xc_lowlevel>
  $<TARGET_OBJECTS:39_libpaw>
  $<TARGET_OBJECTS:41_geometry>
  $<TARGET_OBJECTS:42_nlstrain>
  $<TARGET_OBJECTS:42_parser>
  $<TARGET_OBJECTS:43_ptgroups>
  $<TARGET_OBJECTS:43_wvl_wrappers>
  $<TARGET_OBJECTS:44_abitools>
  $<TARGET_OBJECTS:44_abitypes_defs>
  $<TARGET_OBJECTS:45_geomoptim>
  $<TARGET_OBJECTS:45_xgTools>
  $<TARGET_OBJECTS:48_diago>
  $<TARGET_OBJECTS:49_gw_toolbox_oop>
  $<TARGET_OBJECTS:51_manage_mpi>
  $<TARGET_OBJECTS:52_fft_mpi_noabirule>
  $<TARGET_OBJECTS:53_ffts>
  $<TARGET_OBJECTS:54_spacepar>
  $<TARGET_OBJECTS:55_abiutil>
  $<TARGET_OBJECTS:56_io_mpi>
  $<TARGET_OBJECTS:56_mixing>
  $<TARGET_OBJECTS:56_recipspace>
  $<TARGET_OBJECTS:56_xc>
  $<TARGET_OBJECTS:57_iopsp_parser>
  $<TARGET_OBJECTS:57_iovars>
  $<TARGET_OBJECTS:61_occeig>
  $<TARGET_OBJECTS:62_cg_noabirule>
  $<TARGET_OBJECTS:62_ctqmc>
  $<TARGET_OBJECTS:62_iowfdenpot>
  $<TARGET_OBJECTS:62_poisson>
  $<TARGET_OBJECTS:62_wvl_wfs>
  $<TARGET_OBJECTS:64_psp>
  $<TARGET_OBJECTS:65_lotf_base>
  $<TARGET_OBJECTS:65_paw>
  $<TARGET_OBJECTS:66_nonlocal>
  $<TARGET_OBJECTS:66_vdwxc>
  $<TARGET_OBJECTS:66_wfs>
  $<TARGET_OBJECTS:67_common>
  $<TARGET_OBJECTS:68_dmft>
  $<TARGET_OBJECTS:68_lotf>
  $<TARGET_OBJECTS:68_recursion>
  $<TARGET_OBJECTS:68_rsprc>
  $<TARGET_OBJECTS:69_wfdesc>
  $<TARGET_OBJECTS:70_gw>
  $<TARGET_OBJECTS:71_bse>
  $<TARGET_OBJECTS:71_wannier>
  $<TARGET_OBJECTS:72_response>
  $<TARGET_OBJECTS:77_ddb>
  $<TARGET_OBJECTS:77_suscep>
  $<TARGET_OBJECTS:78_effpot>
  $<TARGET_OBJECTS:78_eph>
  $<TARGET_OBJECTS:79_seqpar_mpi>
  $<TARGET_OBJECTS:80_rttddft>
  $<TARGET_OBJECTS:80_tdep>
  $<TARGET_OBJECTS:94_scfcv>
  $<TARGET_OBJECTS:95_drive>
  )

if(TARGET 17_gpu_toolbox)
  target_sources(abinit_lib PUBLIC $<TARGET_OBJECTS:17_gpu_toolbox>)
endif()

if(TARGET 44_manage_kokkos)
  target_sources(abinit_lib PUBLIC $<TARGET_OBJECTS:44_manage_kokkos>)
endif()

if(TARGET 46_manage_cuda)
  target_sources(abinit_lib PUBLIC $<TARGET_OBJECTS:46_manage_cuda>)
endif()

if (TARGET 67_python_invocation_ext)
  target_sources(abinit_lib PUBLIC $<TARGET_OBJECTS:67_python_invocation_ext>)
endif()

if (TARGET 67_triqs_ext)
  target_sources(abinit_lib PUBLIC $<TARGET_OBJECTS:67_triqs_ext>)
endif()

set_target_properties(abinit_lib
  PROPERTIES OUTPUT_NAME abinit)


#
# Add dummy example to illustrate how to build a yakl application.
#
if(ABINIT_YAKL_WANTED)

  add_executable(yakl_simple_demo yakl_simple_demo.F90)
  set_property(TARGET yakl_simple_demo PROPERTY CXX_STANDARD 17)
  set_target_properties(yakl_simple_demo PROPERTIES LINKER_LANGUAGE Fortran)

  target_include_directories(yakl_simple_demo
    PRIVATE
    ${CMAKE_BINARY_DIR}
    ${CMAKE_SOURCE_DIR}/shared/common/src/incs
    )

  target_compile_definitions(yakl_simple_demo
    PRIVATE
    HAVE_CONFIG_H)

  # we link with abinit::10_defs to avoid redefining std_out
  target_link_libraries(yakl_simple_demo
    PUBLIC
    abinit::yakl
    abinit::10_defs)

  target_include_directories(yakl_simple_demo
    PUBLIC
    ${YAKL_BINARY_DIR})

endif()

# gnu compatibility,
# see https://cmake.org/cmake/help/latest/module/GNUInstallDirs.html
include(GNUInstallDirs)

################################# EXPORT CONFIG #################################
include(CMakePackageConfigHelpers)

# setup some variables
set(version_config ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake)
set(project_config_src ${PROJECT_SOURCE_DIR}/cmake/${PROJECT_NAME}-config.cmake.in)
set(project_config_dst ${PROJECT_BINARY_DIR}/${PROJECT_NAME}-config.cmake)
set(targets_export_name ${PROJECT_NAME}-targets)

# important variables
set(INSTALL_BINDIR ${CMAKE_INSTALL_BINDIR} CACHE STRING
  "Installation directory for executables, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_LIBDIR ${CMAKE_INSTALL_LIBDIR} CACHE STRING
  "Installation directory for libraries, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_INCLUDEDIR ${CMAKE_INSTALL_INCLUDEDIR} CACHE STRING
  "Installation directory for include files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_PKGCONFIG_DIR ${CMAKE_INSTALL_LIBDIR}/pkgconfig CACHE PATH
  "Installation directory for pkgconfig (.pc) files, relative to ${CMAKE_INSTALL_PREFIX}.")

set(INSTALL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake CACHE STRING
  "Installation directory for cmake files, relative to ${CMAKE_INSTALL_PREFIX}.")

# Generate the version, config and target files into the build directory.
write_basic_package_version_file(
  ${version_config}
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Generate cmake my_package-config.cmake file
configure_package_config_file(
  ${project_config_src}
  ${project_config_dst}
  INSTALL_DESTINATION ${INSTALL_CMAKE_DIR})
  # Use a namespace because CMake provides better diagnostics
  # for namespaced imported targets.
export(
  TARGETS abinit_lib
  NAMESPACE abinit::
  FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

export(
  TARGETS 39_libpaw 33_xc_lowlevel 28_numeric_noabirule 27_toolbox_oop 18_timing 17_yaml_out 16_hideleave 14_hidewrite 12_hide_mpi 11_memory_mpi 10_defs 02_clib
  NAMESPACE abinit::
  APPEND FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)

if (ABINIT_ENABLE_GPU_CUDA)
  export(
    TARGETS 17_gpu_toolbox
    NAMESPACE abinit::
    APPEND FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)
endif()

if (ABINIT_ENABLE_GPU_CUDA AND ABINIT_YAKL_BUILD)
  export(
    TARGETS yakl
    NAMESPACE abinit::
    APPEND FILE ${PROJECT_BINARY_DIR}/${targets_export_name}.cmake)
endif()

# macro helper to generate pkg-config file abinit.pc
include(generate_pkgconfig)
generate_pkgconfig(abinit)

# install executables and libabinit
# archive => static libraries
# runtime => shared libraries
install(
  TARGETS abinit abitk aim anaddb atdep band2eps conducti cut3d dummy_tests fftprof
  fold2Bloch ioprof lapackprof macroave mrgddb mrgdv mrggkk mrgscr multibinit
  optic testtransposer ujdet vdw_kernelgen
  abinit_lib 39_libpaw 33_xc_lowlevel 28_numeric_noabirule 27_toolbox_oop 18_timing 17_yaml_out 16_hideleave 14_hidewrite 12_hide_mpi 11_memory_mpi 10_defs 02_clib
  EXPORT ${targets_export_name}
  ARCHIVE DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
  LIBRARY DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
  RUNTIME DESTINATION ${INSTALL_BINDIR} COMPONENT bin
  )

if (ABINIT_ENABLE_GPU_CUDA)
  install(
    TARGETS 17_gpu_toolbox
    EXPORT ${targets_export_name}
    ARCHIVE DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    LIBRARY DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    RUNTIME DESTINATION ${INSTALL_BINDIR} COMPONENT bin
    )
endif()

if (ABINIT_ENABLE_GPU_CUDA AND ABINIT_YAKL_BUILD)
  install(
    TARGETS yakl
    EXPORT ${targets_export_name}
    ARCHIVE DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    LIBRARY DESTINATION ${INSTALL_LIBDIR} COMPONENT lib
    RUNTIME DESTINATION ${INSTALL_BINDIR} COMPONENT bin
    )
endif()

# install cmake config and targets
install(
  FILES ${project_config_dst} ${version_config}
  DESTINATION ${INSTALL_CMAKE_DIR})

install(
  EXPORT ${targets_export_name}
  DESTINATION ${INSTALL_CMAKE_DIR}
  NAMESPACE abinit::)

# install pkgconfig
install(
  FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}.pc
  DESTINATION "${INSTALL_PKGCONFIG_DIR}")

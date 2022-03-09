# using 3.18 to have nice lapack imported target
cmake_minimum_required(VERSION 3.18)

#
# User resources to help writing CMake for fortran projects:
# - book: CMake CookBook + https://github.com/dev-cafe/cmake-cookbook
# - https://cliutils.gitlab.io/modern-cmake/
# - https://enccs.github.io/cmake-workshop/cxx-fortran/
# - https://github.com/scivision/fortran2018-examples
# - MKL + pkg-config : https://www.intel.com/content/www/us/en/developer/articles/technical/intel-math-kernel-library-intel-mkl-and-pkg-config-tool.html?wapkw=mkl%20pkgconfig

#
# Allow find_package() to use <PackageName>_ROOT variables,
# either cmake variable or environment variables
# see https://cmake.org/cmake/help/latest/policy/CMP0074.html
#
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

project(abinit
  LANGUAGES C CXX Fortran
  VERSION 9.9.0
  HOMEPAGE_URL https://www.abinit.org/
  DESCRIPTION "ABINIT is an atomic-scale simulation software suite.")

#
# Export compile command to json (for editors like emacs, clion, vim, etc...).
# It allows nice code editing features provided by LSP (Language Server Protocol)
#
set(CMAKE_EXPORT_COMPILE_COMMANDS on)

#
# default local cmake macro repository (located in abinit top level sources)
#
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
list(INSERT CMAKE_MODULE_PATH 0 "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules")

#
# Enforce user to use a build directory outside of source tree
#
include(prevent_build_in_source)
prevent_build_in_source()

#
# Set default compile optimization flag
#
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(DEFAULT_BUILD_TYPE "RelWithDebInfo")
  message(STATUS "Setting build type to '${DEFAULT_BUILD_TYPE}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${DEFAULT_BUILD_TYPE}" CACHE STRING
    "Choose the type of build, options are: Debug, Release, RelWithDebInfo and MinSizeRel." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif ()

#
# Use CPM project manager for external dependencies.
# See https://github.com/TheLartians/CPM.cmake for more info
#
#include(cmake/CPM.cmake)

#
# Mixing Fortran and C
#
include(FortranCInterface)
FortranCInterface_VERIFY()
FortranCInterface_HEADER(
  fc_mangle.h
  MACRO_NAMESPACE "FC_"
  )

string(TOLOWER ${CMAKE_BUILD_TYPE} build_type)
if (build_type STREQUAL debug)
  set(DEBUG_MODE 1)
endif()

# output dir for libraries
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    CACHE
    PATH "Single output directory for building all libraries.")

#------------------------------------------------------------------------------#
#
# Additional compiler flags, compiler specific
#
#------------------------------------------------------------------------------#
include(compiler_setup)

#------------------------------------------------------------------------------#
#
# Top level Options
#
#------------------------------------------------------------------------------#
include(options)

#------------------------------------------------------------------------------#
#
# Check external dependencies
#
#------------------------------------------------------------------------------#

#
# BLAS / LAPACK setup
#
# this is controlled by cmake var BLA_VENDOR
# available values slightly depends on cmake version
# currently in cmake 3.18, valid values are
#
#    Goto
#    OpenBLAS
#    FLAME
#    ATLAS PhiPACK
#    CXML
#    DXML
#    SunPerf
#    SCSL
#    SGIMATH
#    IBMESSL
#    Intel10_32 (intel mkl v10 32 bit)
#    Intel10_64lp (intel mkl v10+ 64 bit, threaded code, lp64 model)
#    Intel10_64lp_seq (intel mkl v10+ 64 bit, sequential code, lp64 model)
#    Intel10_64ilp (intel mkl v10+ 64 bit, threaded code, ilp64 model)
#    Intel10_64ilp_seq (intel mkl v10+ 64 bit, sequential code, ilp64 model)
#    Intel10_64_dyn (intel mkl v10+ 64 bit, single dynamic library)
#    Intel (obsolete versions of mkl 32 and 64 bit)
#    ACML
#    ACML_MP
#    ACML_GPU
#    Apple
#    NAS
#    Arm
#    Arm_mp
#    Arm_ilp64
#    Arm_ilp64_mp
#    Generic
#
# NOTE:
# - if variable BLA_VENDOR is not set, all vendors are considered.
# - if multiple implementation are available on your system, the first to be detected,
#   will be considered
#
# Recommended value for MKL : Intel10_32
find_package(BLAS)
find_package(LAPACK)

#
# Select FFT flavor
# available flavors are: FFTW3, MKL_DFTI, GOEDECKER, PFFT
#
set(ABINIT_FFT_FLAVOR "FFTW3" CACHE STRING
  "Possible FFT flavors are: FFTW3, MKL_DFTI, GOEDECKER or PFFT")
set_property(CACHE ABINIT_FFT_FLAVOR PROPERTY STRINGS FFTW3 MKL_DFTI GOEDECKER PFFT)

#
# Select SCALAPACK flavor
# available flavors are: NETLIB, MKL
#
# Important notes regarding MKL:
#
# If ABINIT_SCALAPACK_FLAVOR is MKL then  scalapack/mkl
#
# User must provide env variable MKL_ROOT if SCALAPACK/MKL
# is desired.
#
# Other variables of interest that are used to initialize
# MKL_BLAS_FLAVOR_DEFAULT and MKL_SCALAPACK_FLAVOR_DEFAULT
# - MKL_XX_BIT (32bit or 64bit)
# - MKL_THREADING_TYPE (omp or seq)
#
# User can override directly MKL_BLAS_FLAVOR and MKL_SCALAPACK_FLAVOR
#
set(ABINIT_SCALAPACK_FLAVOR "NETLIB" CACHE STRING
  "Possible SCALAPACK flavors are: NETLIB or MKL")
set_property(CACHE ABINIT_SCALAPACK_FLAVOR PROPERTY STRINGS NETLIB MKL)


#####################################################################
# checking for optional external dependencies :
# - MPI (OpenMPI, IntelMPI, MPICH, etc...)
# - OpenMP
# - MKL (Intel classic or newest OneAPI/MKL ?)
# - Scalapack (MKL or Netlib)
# - Libxc
#
# For each external dependency, you can :
# - instruct cmake to search for a particular directory
#   using either env variable or cmake variable <PackageName>_ROOT
# - (recommended) set env variable PKG_CONFIG_PATH to specify a list
#   of directories containing pkgconfig files (extension .pc)
# - (preferred) whenever possible, if dependency provides a cmake script
#   to detect it, just set env variable CMAKE_PREFIX_PATH
#
# On a supercomputer, most of the time, the modules environment will
# set either PKG_CONFIG_PATH and CMAKE_PREFIX_PATH
#####################################################################

# Enable MODULE mode (via pkg-config) to find packages
# see https://cmake.org/cmake/help/latest/command/find_package.html
find_package(PkgConfig REQUIRED)

#
# MPI
#
set(MPI_DETERMINE_LIBRARY_VERSION TRUE)
find_package(MPI COMPONENTS Fortran C CXX)
if(MPI_FOUND)
  add_compile_definitions(HAVE_MPI)
endif()

# determine mpi vendor (MPI_VENDOR)
# valid values are : OpenMPI, MPICH or IntelMPI
include(get_mpi_vendor)
get_mpi_vendor()

if (NOT DEFINED MPI_VERSION AND DEFINED MPI_Fortran_VERSION)
  set(MPI_VERSION ${MPI_Fortran_VERSION})
  set(MPI_VERSION_MAJOR ${MPI_Fortran_VERSION_MAJOR})
  set(MPI_VERSION_MINOR ${MPI_Fortran_VERSION_MINOR})
endif()

if (MPI_FOUND)

  # try to detect if our MPI implementation is cuda-aware
  # 1. If using OpenMPI, then try to parse output of ompi_info
  # 2. If using MPICH, TODO
  # 3. Other MPI impl, TODO
  if (ABINIT_ENABLE_GPU_CUDA)

    if (MPI_VENDOR STREQUAL "OpenMPI")
      find_program(OMPI_INFO
        NAMES ompi_info
        HINTS ${MPI_CXX_LIBRARIES}/../bin)

      # Full command line to probe if cuda support in MPI implementation is enabled
      # ompi_info --parsable --all | grep mpi_built_with_cuda_support:value
      if (OMPI_INFO)
        execute_process(COMMAND ${OMPI_INFO}
          OUTPUT_VARIABLE _output)
        if ( (_output MATCHES "smcuda") )
          message(STATUS "Found OpenMPI with CUDA-aware support available.")
          set(HAVE_GPU_MPI 1)
        else()
          message(WARNING "OpenMPI found, but it is not built with CUDA support.")
        endif()
      endif(OMPI_INFO)
    endif()

    if (MPI_VENDOR STREQUAL "MPICH")
      # TODO
    endif()

    if (MPI_VENDOR STREQUAL "IntelMPI")
      # TODO
    endif()

    # unconditionally set HAVE_GPU_MPI if the user say so
    # this may be need because our detection above is not bullet proof
    if (ABINIT_ENFORCE_CUDA_AWARE_MPI)
      set(HAVE_GPU_MPI 1)
    endif()

  endif(ABINIT_ENABLE_GPU_CUDA)

endif(MPI_FOUND)

#
# OpenMP
#
find_package(OpenMP)

# check for OpenMP collapse clause support
if(OpenMP_FOUND)
  set(HAVE_OPENMP 1)

  try_compile(HAVE_OMP_COLLAPSE_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_omp_collapse.F90
    LINK_LIBRARIES OpenMP::OpenMP_Fortran)

  if(HAVE_OMP_COLLAPSE_BOOL)
    set(HAVE_OMP_COLLAPSE 1)
  endif()

endif()

#
# Kokkos : https://github.com/kokkos/kokkos
#
include(build_or_find_kokkos)

#
# YAKL : https://github.com/mrnorman/YAKL
#
include(build_or_find_yakl)

#
# Scalapack setup
#
if (ABINIT_SCALAPACK_ENABLED)
  # setup abinit to use either netlib or mkl scalapack
  include(scalapack_setup)
endif()

#
# ELPA setup
#
if (ABINIT_ELPA_ENABLED)
  pkg_check_modules(ABINIT_ELPA QUIET IMPORTED_TARGET elpa)

  if(ABINIT_ELPA_FOUND)
    message("ELPA found via pkg-config")
    set(ELPA_FOUND TRUE)
    add_library(abinit::elpa ALIAS PkgConfig::ABINIT_ELPA)
  endif()
endif()

#
# Libxc
# 1. prefer using pkg-config,
# 2. use cmake target file if not found with pkg-config
#
pkg_check_modules(ABINIT_LIBXC QUIET IMPORTED_TARGET libxc)

if(ABINIT_LIBXC_FOUND)
  message("Libxc found via pkg-config")
  set(Libxc_FOUND TRUE)
  #set(Libxc_VERSION ${LibxcPkg_VERSION})
  #set(Libxc_INCLUDE_DIRS ${LibxcPkg_INCLUDE_DIRS})
  add_library(abinit::libxc ALIAS PkgConfig::ABINIT_LIBXC)
else()
  find_package(Libxc REQUIRED)
  message("Libxc found via cmake target")
endif()

if(Libxc_FOUND)
  add_compile_definitions(HAVE_LIBXC)
  set(HAVE_LIBXC 1)
endif()

#
# FFT FLAVOR :
#
if(ABINIT_FFT_FLAVOR STREQUAL "FFTW3")
  # use regular fftw3 library
  pkg_check_modules(ABINIT_FFTW3 QUIET IMPORTED_TARGET fftw3 fftw3f)
  if(ABINIT_FFTW3_FOUND)
    message("Selected fft flavor is FFTW3")
    add_library(abinit::fft ALIAS PkgConfig::ABINIT_FFTW3)
    set(HAVE_FFTW3 1)
  else()
    message(WARNING "Please adjust environment variable PKG_CONFIG_PATH to setup path to")
    message(WARNING "fftw3.pc or select another fft flavor.")
    message(FATAL_ERROR "fft flavor FFTW3 can't be found")
  endif()

  # TODO : search for fftw3_threads and create imported target
  # set(HAVE_FFTW3_THREADS 1)

  # TODO : search for fftw3_mpi and create imported target
  # set(HAVE_FFTW3_MPI 1)

elseif(ABINIT_FFT_FLAVOR STREQUAL "MKL_DFTI")

  # Note: OneMKL (version >= 2021.3) provides a MKLConfig.cmake, usually available
  # by setting env variable CMAKE_PREFIX_PATH to ${MKLROOT}/lib/cmake
  find_package(MKL)
  if(MKL_FOUND)
    message("Selected fft flavor is MKL_DFTI")
    # TODO : create imported target for mkl_intel_ilp64 / mkl_intel
    add_library(abinit::fft ALIAS MKL::MKL)
  else()
    message(WARNING "Please check your env variable CMAKE_PREFIX_PATH to contain path to")
    message(WARNING "MKLConfig.cmake (located into your intel OneAPI install path) or selected another flavor")
    message(FATAL_ERROR "MKL is not found, and thus MKL_DFTI fft flavor can't be set")
  endif()

elseif(ABINIT_FFT_FLAVOR STREQUAL "PFFT")

  pkg_check_modules(ABINIT_PFFT QUIET IMPORTED_TARGET pfft)
  if(ABINIT_PFFT_FOUND)
    message("pfft found via pkg-config")
    add_library(abinit::fft ALIAS PkgConfig::ABINIT_PFFT)
    set(HAVE_PFFT 1)
  else()
    message(FATAL_ERROR "pfft NOT FOUND, selected another flavor or adjust env variable PKG_CONFIG_PATH")
  endif()

  message("Selected fft flavor is PFFT")

elseif(ABINIT_FFT_FLAVOR STREQUAL "GOEDECKER")

  message("Selected fft flavor is GOEDECKER (internal implementation, no link flags required).")
  # WARNING: abinit::fft alias is not defined

endif()

#
# LEVMAR : use LevMar with cmake build-system
# see : https://github.com/pkestene/levmar-cmake
#
find_package(levmar)
if(levmar_FOUND)
  add_library(abinit::levmar ALIAS levmar::levmar)
  set(HAVE_LEVMAR 1)
else()
  message("levmar NOT FOUND")
endif()

#
# LIBPSML : TODO
#
pkg_check_modules(ABINIT_LIBPSML QUIET IMPORTED_TARGET libpsml)
if(ABINIT_LIBPSML_FOUND)
  message("libpsml found via pkg-config")
  add_library(abinit::libpsml ALIAS PkgConfig::ABINIT_LIBPSML)
  set(HAVE_LIBPSML 1)
else()
  message("libpsml NOT FOUND")
endif()

#
# HDF5
#
# HDF5 is not really used for now apparently
if (ABINIT_ENABLE_PARALLEL_HDF5)
  set(HDF5_PREFER_PARALLEL TRUE)
endif()

find_package(HDF5)
if(HDF5_FOUND)
  message("HDF5 found via find_package")
  add_library(abinit::hdf5 ALIAS HDF5::HDF5)
  set(HAVE_HDF5 1)

  if(HDF5_IS_PARALLEL)
    set(HAVE_HDF5_MPI 1)
  endif()
endif(HDF5_FOUND)


#
# NetCDF_C / NetCDF_Fortran
#
pkg_check_modules(ABINIT_NETCDF QUIET IMPORTED_TARGET netcdf)
if(ABINIT_NETCDF_FOUND)
  message("netcdf found via pkg-config")
  add_library(abinit::netcdf ALIAS PkgConfig::ABINIT_NETCDF)
else()
  message("netcdf NOT FOUND")
endif()

pkg_check_modules(ABINIT_NETCDF_FORTRAN QUIET IMPORTED_TARGET GLOBAL netcdf-fortran)
if(ABINIT_NETCDF_FORTRAN_FOUND)
  message("netcdf-fortran found via pkg-config : ${ABINIT_NETCDF_FORTRAN_PREFIX}")

  # kind of fix to make sure path PREFIX/include (where netcdf.mod is located)
  # is append in include dirs used in abinit::netcdf-fortran
  pkg_get_variable(ABINIT_NETCDF_FORTRAN_PREFIX netcdf-fortran prefix)
  get_target_property(ABINIT_NETCDF_FORTRAN_INC_TMP PkgConfig::ABINIT_NETCDF_FORTRAN INTERFACE_INCLUDE_DIRECTORIES)
  list(APPEND ABINIT_NETCDF_FORTRAN_INC_TMP "${ABINIT_NETCDF_FORTRAN_PREFIX}/include")
  set_target_properties(PkgConfig::ABINIT_NETCDF_FORTRAN PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ABINIT_NETCDF_FORTRAN_INC_TMP}")

  add_library(abinit::netcdf-fortran ALIAS PkgConfig::ABINIT_NETCDF_FORTRAN)

else()
  message("netcdf-fortran NOT FOUND")
endif()

if (MPI_FOUND)
  pkg_check_modules(ABINIT_NETCDF_MPI QUIET IMPORTED_TARGET netcdf-mpi)
  if(ABINIT_NETCDF_MPI_FOUND)
    message("netcdf-mpi found via pkg-config")
    add_library(abinit::netcdf-mpi ALIAS PkgConfig::ABINIT_NETCDF_MPI)
    #set(HAVE_NETCDF_MPI 1) => see generate_config.h
  else()
    message("netcdf-mpi NOT FOUND")
  endif()
endif()


#
# PAPI
#
pkg_check_modules(ABINIT_PAPI QUIET IMPORTED_TARGET papi)
if(ABINIT_PAPI_FOUND)

  # check if f90papi.h is available
  # note : can't use check_include_file macro here, because f90papi.h contains
  # actually fortran code
  find_path(PAPI_INCLUDE_F90PAPI NAMES f90papi.h
    HINTS ${ABINIT_PAPI_INCLUDE_DIRS})
  if (PAPI_INCLUDE_F90PAPI-NOTFOUND)
    message("f90papi.h not found")
  else()
    set(HAVE_F90PAPI_H 1)
  endif()

  # check if papi.h is available
  include(CheckIncludeFile)
  set(CMAKE_REQUIRED_INCLUDES ${ABINIT_PAPI_INCLUDEDIR})
  check_include_file(papi.h HAVE_PAPI_H_BOOL)
  set(CMAKE_REQUIRED_INCLUDES)
  if(HAVE_PAPI_H_BOOL)
    set(HAVE_PAPI_H 1)
  endif()

  add_library(abinit::papi ALIAS PkgConfig::ABINIT_PAPI)

  message("papi found via pkg-config")
else()
  message("papi NOT FOUND")
endif()

#
# TRIQS : TODO improve
# to help cmake detect TRIQS:
# - set env var or cmake var TRIQS_ROOT to the root path of triqs install directory
# - set env var CMAKE_PREFIX_PATH (usually by a modulefile on a supercomputer)
#
find_package(TRIQS)

if(TRIQS_FOUND)
  add_library(abinit::triqs ALIAS triqs)

  set(HAVE_TRIQS 1)

  # NOT sure about this, is TRIQS 3.x compatible with 2.0 ?
  # (TODO : check with someone who knows triqs)
  # Besides, it looks like currently abiniti is not compatible with triqs 3.x,
  # some porting is needed
  # if (TRIQS_VERSION_MAJOR VERSION_EQUAL "3")
  #   set(DO_BUILD_67_TRIQS_EXT ON)
  #   set(HAVE_TRIQS_v3 1)
  #   set(HAVE_TRIQS_v2_0 1)
  # endif()


  if (TRIQS_VERSION_MAJOR VERSION_EQUAL "2" AND TRIQS_VERSION_MINOR VERSION_EQUAL "0")
    set(DO_BUILD_67_TRIQS_EXT ON)
    set(HAVE_TRIQS_v2_0 1)
  endif()

  if (TRIQS_VERSION_MAJOR VERSION_EQUAL "1" AND TRIQS_VERSION_MINOR VERSION_EQUAL "4")
    set(DO_BUILD_67_TRIQS_EXT ON)
    set(HAVE_TRIQS_v1_4 1)
  endif()

  #get_property(TRIQS_COMPILE_DEFINITIONS TARGET triqs PROPERTY INTERFACE_COMPILE_DEFINITIONS)
  #get_property(TRIQS_COMPILE_FEATURES TARGET triqs PROPERTY INTERFACE_COMPILE_FEATURES)
  #get_property(TRIQS_COMPILE_OPTIONS TARGET triqs PROPERTY INTERFACE_COMPILE_OPTIONS)
  get_property(TRIQS_INCLUDE_DIRECTORIES TARGET triqs PROPERTY INTERFACE_INCLUDE_DIRECTORIES)
  get_property(TRIQS_LINK_LIBRARIES TARGET triqs PROPERTY INTERFACE_LINK_LIBRARIES)
else()
  message("triqs NOT FOUND")
endif()

#
# WANNIER
#
include(build_or_find_wannier90)

#
# libxml2
#
pkg_check_modules(ABINIT_XML2 QUIET IMPORTED_TARGET libxml-2.0)
if(ABINIT_XML2_FOUND)
  message("xml2 found via pkg-config")
  add_library(abinit::xml2 ALIAS PkgConfig::ABINIT_XML2)
  set(HAVE_XML 1)
else()
  message("xml2 NOT FOUND")
endif()

#
# XMLF90
#
pkg_check_modules(ABINIT_XMLF90 QUIET IMPORTED_TARGET xmlf90)
if(ABINIT_XMLF90_FOUND)
  message("xmlf90 found via pkg-config")
  add_library(abinit::xmlf90 ALIAS PkgConfig::ABINIT_XMLF90)
  set(HAVE_XMLF90 1)
else()
  message("xmlf90 NOT FOUND")
endif()


#
# AtomPAW
#
pkg_check_modules(ABINIT_ATOMPAW QUIET IMPORTED_TARGET atompaw)
if(ABINIT_ATOMPAW_FOUND)
  message("atompaw found via pkg-config")
  add_library(abinit::atompaw ALIAS PkgConfig::ABINIT_ATOMPAW)
  set(HAVE_ATOMPAW 1)
else()
  message("atompaw NOT FOUND")
endif()

#
# BigDFT
#
pkg_check_modules(ABINIT_BIGDFT QUIET IMPORTED_TARGET bigdft)
if(ABINIT_BIGDFT_FOUND)
  message("bigdft found via pkg-config")
  add_library(abinit::bigdft ALIAS PkgConfig::ABINIT_BIGDFT)
  set(HAVE_BIGDFT 1)
else()
  message("bigdft NOT FOUND")
endif()

#
# Generate config.h
#
include(generate_config_h)

##################### PRINT CONFIGURE STATUS ######################
message("//===================================================")
message("// ${PROJECT_NAME} build configuration:")
message("// ${PROJECT_NAME} version : ${PROJECT_VERSION}")
message("//===================================================")
message("  CMake version          : ${CMAKE_VERSION}")
if (NOT CMAKE_BUILD_TYPE)
  message("  CMake build type       : NOT SET !")
else()
  message("  CMake build type       : ${CMAKE_BUILD_TYPE}")
endif()
message("  CMake install prefix   : ${CMAKE_INSTALL_PREFIX}")
message("  CMake system processor : ${CMAKE_SYSTEM_PROCESSOR}")
message("  CMake system name (OS) : ${CMAKE_SYSTEM_NAME}")
message("")

message("  Fortran compiler Id      : ${CMAKE_Fortran_COMPILER_ID}")
message("  Fortran compiler version : ${CMAKE_Fortran_COMPILER_VERSION}")
message("  Fortran compiler exe     : ${CMAKE_Fortran_COMPILER}")
message("  Fortran flags            : ${CMAKE_Fortran_FLAGS}")
message("  Fortran compiler wrapper : ${CMAKE_Fortran_COMPILER_WRAPPER}")
message("")

message(STATUS "MPI config:")
message("    MPI found            : ${MPI_FOUND}")
message("    MPI standard version : ${MPI_Fortran_VERSION}")
# note : make sure you correctly your C, C++ and Fortran compiler
# through variable CC, CXX and FC
message("    MPI_VENDOR           : ${MPI_VENDOR}")
message("    MPI library version  : ${MPI_Fortran_LIBRARY_VERSION_STRING}")
message("    MPI fortran compiler : ${MPI_Fortran_COMPILER}")
message("    MPI headers          : ${MPI_Fortran_INCLUDE_DIRS}")
message("    MPI libraries        : ${MPI_Fortran_LIBRARIES}")
#message("    MPI link flags       : ${MPI_Fortran_LINK_FLAGS}")
message("")

message(STATUS "OpenMP config :")
message("    OpenMP found         : ${OpenMP_FOUND}")
message("    OpenMP version       : ${OpenMP_Fortran_VERSION}")
message("")


message(STATUS "Kokkos config :")
message("    Kokkos wanted         : ${ABINIT_KOKKOS_WANTED}")
if(ABINIT_KOKKOS_WANTED)
  message("    ABINIT_KOKKOS_BUILD   : ${ABINIT_KOKKOS_BUILD}")
  if (ABINIT_KOKKOS_BUILD)
    message("    ABINIT_KOKKOS_BACKEND : ${ABINIT_KOKKOS_BACKEND}")
    message("    Kokkos_ENABLE_HWLOC   : ${Kokkos_ENABLE_HWLOC}")
    message("    Kokkos_ENABLE_OPENMP  : ${Kokkos_ENABLE_OPENMP}")
    message("    Kokkos_ENABLE_CUDA    : ${Kokkos_ENABLE_CUDA}")
  else()
    message("    Using installed Kokkos : ${Kokkos_DIR}")
    message("    Kokkos version         : ${Kokkos_VERSION}")
    message("    Kokkos_CXX_COMPILER    : ${Kokkos_CXX_COMPILER}")
    message("    Kokkos_CXX_COMPILER_ID : ${Kokkos_CXX_COMPILER_ID}")
    message("    Kokkos_CXX_STANDARD    : ${Kokkos_CXX_STANDARD}")
    message("    Kokkos_OPTIONS         : ${Kokkos_OPTIONS}")
    message("    Kokkos_TPLS            : ${Kokkos_TPLS}")
    message("    Kokkos_DIRS            : ${Kokkos_DIR}")
  endif()
  message("    CMAKE_CXX_STANDARD     : ${CMAKE_CXX_STANDARD}")
endif()
message("")


message(STATUS "Yakl config :")
message("    YAKL wanted         : ${ABINIT_YAKL_WANTED}")
if(ABINIT_YAKL_WANTED)
  message("    YAKL_ARCH           : ${YAKL_ARCH}")
endif()
message("")


message(STATUS "BLAS config:")
message("  BLA_VENDOR : ${BLA_VENDOR}")
if (BLAS_FOUND)
  message("  BLAS_LINKER_FLAGS : ${BLAS_LINKER_FLAGS}")
  message("  BLAS_LIBRARIES    : ${BLAS_LIBRARIES}")
else()
  message("  BLAS not found ! Please recheck your environment variables.")
endif()
message("")

message(STATUS "LAPACK config:")
message("  BLA_VENDOR : ${BLA_VENDOR}")
if (LAPACK_FOUND)
  message("  LAPACK_LINKER_FLAGS : ${LAPACK_LINKER_FLAGS}")
  message("  LAPACK_LIBRARIES    : ${LAPACK_LIBRARIES}")
else()
  message("  LAPACK not found ! Please recheck your environment variables.")
endif()
message("")

message(STATUS "ScaLapack config:")
message("    ABINIT_SCALAPACK_ENABLED : ${ABINIT_SCALAPACK_ENABLED}")
if(ABINIT_SCALAPACK_ENABLED)
  message("    ABINIT_SCALAPACK_FLAVOR  : ${ABINIT_SCALAPACK_FLAVOR}")
  if(ABINIT_SCALAPACK_FLAVOR MATCHES "MKL")
    message("    MKL found                  : ${MKL_FOUND}")
    message("    MKL_BLAS_FLAVOR            : ${MKL_BLAS_FLAVOR}")
    message("    MKL_BLAS_FLAVOR_FOUND      : ${MKL_BLAS_FLAVOR_FOUND}")
    message("    MKL blas headers           : ${MKL_BLAS_HEADERS}")
    message("    MKL blas libraries         : ${MKL_BLAS_LIBRARIES}")
    message("    MKL_SCALAPACK_FLAVOR       : ${MKL_SCALAPACK_FLAVOR}")
    message("    MKL_SCALAPACK_FLAVOR_FOUND : ${MKL_SCALAPACK_FLAVOR_FOUND}")
    message("    MKL scalapack headers      : ${MKL_SCALAPACK_HEADERS}")
    message("    MKL scalapack libraries    : ${MKL_SCALAPACK_LIBRARIES}")
  elseif(ABINIT_SCALAPACK_FLAVOR MATCHES "NETLIB")
    message("    ScaLapack found            : ${scalapack_FOUND}")
    message("    Scalapack dir              : ${scalapack_DIR}")
    message("    Scalapack found via pkg-config : ${scalapack_found_using_pkg_config}")
    message("    Scalapack found via cmake cfg  : ${scalapack_found_using_cmake_target}")
  endif()
endif()
message("")

message(STATUS "ELPA config:")
message("    ABINIT_ELPA_ENABLED         : ${ABINIT_ELPA_ENABLED}")
if(ABINIT_ELPA_ENABLED)
  message("    ELPA found                  : ${ABINIT_ELPA_FOUND}")
  message("    ELPA version                : ${ABINIT_ELPA_VERSION}")
  message("    ELPA cflags                 : ${ABINIT_ELPA_CFLAGS}")
  message("    ELPA ldflags                : ${ABINIT_ELPA_LDFLAGS}")
endif()
message("")

#
# TODO PK : add the possibility to detect a version of libpaw already installed
#
message(STATUS "Libpaw config (only internal version is possible for now):")
message("    Libpaw (abinit)    : ${HAVE_LIBPAW_ABINIT}")
message("")

message(STATUS "Libxc config:")
message("    Libxc found        : ${Libxc_FOUND}")
message("    Libxc dir          : ${Libxc_DIR}")
message("    Libxc version      : ${Libxc_VERSION}")
message("    Libxc headers      : ${Libxc_INCLUDE_DIRS}")
message("    Libxc libraries    : ${Libxc_LIBRARIES}")
message("")

message(STATUS "fftw3 config")
message("    fftw3 found         : ${ABINIT_FFTW3_FOUND}")
message("    fftw3 version       : ${ABINIT_FFTW3_VERSION}")
message("    fftw3 cflags        : ${ABINIT_FFTW3_CFLAGS}")
message("    fftw3 ldflags       : ${ABINIT_FFTW3_LDFLAGS}")
message("")

message(STATUS "libpsml config")
message("    libpsml found      : ${ABINIT_LIBPSML_FOUND}")
message("    libpsml version    : ${ABINIT_LIBPSML_VERSION}")
message("    libpsml cflags     : ${ABINIT_LIBPSML_CFLAGS}")
message("    libpsml ldflags    : ${ABINIT_LIBPSML_LDFLAGS}")
message("")

message(STATUS "netcdf-c config")
message("    netcdf-c found         : ${ABINIT_NETCDF_FOUND}")
message("    netcdf-c prefix        : ${ABINIT_NETCDF_PREFIX}")
message("    netcdf-c version       : ${ABINIT_NETCDF_VERSION}")
message("    netcdf-c inc dirs      : ${ABINIT_NETCDF_INCLUDE_DIRS}")
message("    netcdf-c cflags        : ${ABINIT_NETCDF_CFLAGS}")
message("    netcdf-c cflags other  : ${ABINIT_NETCDF_CFLAGS_OTHER}")
message("    netcdf-c ldflags       : ${ABINIT_NETCDF_LDFLAGS}")
message("")

message(STATUS "netcdf-fortran config")
message("    netcdf-fortran found         : ${ABINIT_NETCDF_FORTRAN_FOUND}")
message("    netcdf-fortran prefix        : ${ABINIT_NETCDF_FORTRAN_PREFIX}")
message("    netcdf-fortran version       : ${ABINIT_NETCDF_FORTRAN_VERSION}")
message("    netcdf-fortran inc dirs      : ${ABINIT_NETCDF_FORTRAN_INCLUDE_DIRS}")
message("    netcdf-fortran cflags        : ${ABINIT_NETCDF_FORTRAN_CFLAGS}")
message("    netcdf-fortran cflags other  : ${ABINIT_NETCDF_FORTRAN_CFLAGS_OTHER}")
message("    netcdf-fortran ldflags       : ${ABINIT_NETCDF_FORTRAN_LDFLAGS}")
message("    netcdf C/MPI support         : ${HAVE_NETCDF_MPI}")
message("    netcdf fortran/MPI support   : ${HAVE_NETCDF_FORTRAN_MPI}")
message("")

message(STATUS "papi config")
message("    papi found         : ${ABINIT_PAPI_FOUND}")
message("    papi version       : ${ABINIT_PAPI_VERSION}")
message("    papi cflags        : ${ABINIT_PAPI_CFLAGS}")
message("    papi ldflags       : ${ABINIT_PAPI_LDFLAGS}")
message("")

message(STATUS "pfft config")
message("    pfft found         : ${ABINIT_PFFT_FOUND}")
message("    pfft version       : ${ABINIT_PFFT_VERSION}")
message("    pfft cflags        : ${ABINIT_PFFT_CFLAGS}")
message("    pfft ldflags       : ${ABINIT_PFFT_LDFLAGS}")
message("")

message(STATUS "wannier90 config")
if (ABINIT_WANNIER90_WANTED)
  message("    wannier90 builtin     : ${ABINIT_WANNIER90_BUILTIN}")
  message("    wannier90 found       : ${ABINIT_WANNIER90_FOUND}")
  message("    wannier90 root        : ${WANNIER_ROOT}")
  message("    wannier90 version     : ${ABINIT_WANNIER90_VERSION}")
  message("    wannier90 include dir : ${ABINIT_WANNIER90_INCLUDE_DIR}")
  message("    wannier90 library     : ${ABINIT_WANNIER90_LIBRARY}")
else()
  message("    wannier90 not wanted")
endif()
message("")

message(STATUS "triqs config")
message("    triqs found               : ${TRIQS_FOUND}")
message("    triqs version             : ${TRIQS_VERSION}")
message("    triqs root                : ${TRIQS_ROOT}")
message("    triqs python support      : ${TRIQS_WITH_PYTHON_SUPPORT}")
#message("    triqs compile definitions : ${TRIQS_COMPILE_DEFINITIONS}")
#message("    triqs compile features    : ${TRIQS_COMPILE_FEATURES}")
#message("    triqs compile options     : ${TRIQS_COMPILE_OPTIONS}")
message("    triqs include dir         : ${TRIQS_INCLUDE_DIRECTORIES}")
message("    triqs link libraries      : ${TRIQS_LINK_LIBRARIES}")
message("")

message(STATUS "xmlf90 config")
message("    xmlf90 found       : ${ABINIT_XMLF90_FOUND}")
message("    xmlf90 version     : ${ABINIT_XMLF90_VERSION}")
message("    xmlf90 cflags      : ${ABINIT_XMLF90_CFLAGS}")
message("    xmlf90 ldflags     : ${ABINIT_XMLF90_LDFLAGS}")
message("")

message(STATUS "xml-2.0 config")
message("    xml-2.0 found       : ${ABINIT_XML2_FOUND}")
message("    xml-2.0 version     : ${ABINIT_XML2_VERSION}")
message("    xml-2.0 cflags      : ${ABINIT_XML2_CFLAGS}")
message("    xml-2.0 ldflags     : ${ABINIT_XML2_LDFLAGS}")
message("")

message(STATUS "AtomPaw config")
message("    AtomPAW found       : ${ABINIT_ATOMPAW_FOUND}")
message("    AtomPAW version     : ${ABINIT_ATOMPAW_VERSION}")
message("    AtomPAW cflags      : ${ABINIT_ATOMPAW_CFLAGS}")
message("    AtomPAW ldflags     : ${ABINIT_ATOMPAW_LDFLAGS}")
message("")

message(STATUS "BigDFT config")
message("    BigDFT found       : ${ABINIT_BIGDFT_FOUND}")
message("    BigDFT version     : ${ABINIT_BIGDFT_VERSION}")
message("    BigDFT cflags      : ${ABINIT_BIGDFT_CFLAGS}")
message("    BigDFT ldflags     : ${ABINIT_BIGDFT_LDFLAGS}")
message("")

message(STATUS "GPU config")
message("    CUDA enabled       : ${ABINIT_ENABLE_GPU_CUDA}")
if (ABINIT_ENABLE_GPU_CUDA)
  message("  CUDA compiler ID      : ${CMAKE_CUDA_COMPILER_ID}")
  message("  CUDA compiler Version : ${CMAKE_CUDA_COMPILER_VERSION}")
  message("  C++ Compiler : ${CMAKE_CXX_COMPILER_ID} "
    "${CMAKE_CXX_COMPILER_VERSION} "
    "${CMAKE_CXX_COMPILER_WRAPPER}")
  message("    ${CMAKE_CXX_COMPILER}")
  message("  CUDA Compiler        : ${CMAKE_CUDA_COMPILER}")
  message("  CUDA Compiler exec   : ${CUDA_NVCC_EXECUTABLE}")
  message("  CUDA Compile flags   : ${CMAKE_CUDA_FLAGS}")
  message("  CUDA toolkit found   : ${CUDAToolkit_FOUND}")
  message("  CUDA toolkit version : ${CUDAToolkit_VERSION}")
  message("  CUDA toolkit nvcc    : ${CUDAToolkit_NVCC_EXECUTABLE}")
endif()
message("")

message(STATUS "ENV{PKG_CONFIG_PATH} was :")
message("    $ENV{PKG_CONFIG_PATH}")
message("")

message(STATUS "ENV{CMAKE_PREFIX_PATH} was :")
message("    $ENV{CMAKE_PREFIX_PATH}")
message("")

add_subdirectory(shared)
add_subdirectory(src)

# --- make git ignore build directory
if(NOT EXISTS ${PROJECT_BINARY_DIR}/.gitignore)
  file(WRITE ${PROJECT_BINARY_DIR}/.gitignore "*")
endif()

file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/tests)

#
# Robodoc
#
set(ROBODOC_DONE_FILE "robodoc_done")
add_custom_command(OUTPUT ${ROBODOC_DONE_FILE}
  COMMAND rm -rf tmp-robodoc robodoc-html
  COMMAND mkdir tmp-robodoc
  COMMAND cp -rfL ${CMAKE_SOURCE_DIR}/shared/common/src/[0-3]* tmp-robodoc
  COMMAND cp -rf ${CMAKE_SOURCE_DIR}/src/[4-9]* tmp-robodoc
  COMMAND cp ${CMAKE_SOURCE_DIR}/config/robodoc/robodoc-html.rc tmp-robodoc/robodoc.rc
  COMMAND cd tmp-robodoc && rm -f */*.in && rm -f */interfaces* && robodoc > ../robodoc.log 2> ../robodoc.err
  COMMAND mv -f tmp-robodoc/www/robodoc robodoc-html
  COMMAND ${CMAKE_COMMAND} -E tar "zcf" "robodoc-html-${CMAKE_PROJECT_VERSION}.tar.gz" robodoc-html
  COMMAND rm -rf robodoc-html tmp-robodoc
  COMMAND cat ${CMAKE_SOURCE_DIR}/doc/developers/robodoc.doc.txt >> robodoc.err
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  COMMENT "Generating documentation with ROBODOC"
  )

add_custom_target(robodoc DEPENDS ${ROBODOC_DONE_FILE})

#
# Install
#
include(abinit_install)

#
# Test
#
include(abinit_test)

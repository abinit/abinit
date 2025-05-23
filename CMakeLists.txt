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

# Retrieve ABINIT version
file (STRINGS ".current_version" ABI_VERSION)

project(abinit
  LANGUAGES C CXX Fortran
  VERSION ${ABI_VERSION}
  HOMEPAGE_URL https://www.abinit.org/
  DESCRIPTION "ABINIT is an atomic-scale simulation software suite.")

# will collect a list of all executables
Set(ALL_EXECUTABLES "" CACHE INTERNAL "")


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
# Top level Options
#
#------------------------------------------------------------------------------#
include(options)

#------------------------------------------------------------------------------#
#
# Additional compiler flags, compiler specific
#
#------------------------------------------------------------------------------#
include(compiler_setup)

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
# Check for buggy BLAS interfaces
# Set the ABINIT_ENABLE_ZDOT_BUGFIX flag
#
set(ABINIT_ENABLE_ZDOT_BUGFIX AUTO CACHE STRING "Activate workaround for bugged ZDOTC and ZDOTU (default: auto)")
set_property(CACHE ABINIT_ENABLE_ZDOT_BUGFIX PROPERTY STRINGS AUTO ON OFF)
if((ABINIT_ENABLE_ZDOT_BUGFIX STREQUAL "AUTO") OR (NOT ABINIT_ENABLE_ZDOT_BUGFIX))
  try_run(BLAS_BUGGY_RUN BLAS_BUGGY_COMPILE
    ${CMAKE_BINARY_DIR}/try_run ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_linalg_zdot_bug.F90
    LINK_LIBRARIES BLAS::BLAS
    RUN_OUTPUT_VARIABLE BLAS_BUGGY_OUTPUT)
  STRING(STRIP "${BLAS_BUGGY_OUTPUT}" BLAS_BUGGY_OUTPUT)
  if(BLAS_BUGGY_OUTPUT STREQUAL "OK")
    set(ABINIT_ENABLE_ZDOT_BUGFIX OFF CACHE STRING "Activate workaround for bugged ZDOTC and ZDOTU (default: auto)" FORCE)
  else()
    if(NOT ABINIT_ENABLE_ZDOT_BUGFIX)
      message(FATAL_ERROR "ABINIT_ENABLE_ZDOT_BUGFIX option is deactivated but the BLAS library has buggy interfaces!")
    endif()
    set(ABINIT_ENABLE_ZDOT_BUGFIX ON CACHE STRING "Activate workaround for bugged ZDOTC and ZDOTU (default: auto)" FORCE)
  endif()
endif()
if(ABINIT_ENABLE_ZDOT_BUGFIX)
  set(HAVE_LINALG_ZDOTC_BUG 1)
  set(HAVE_LINALG_ZDOTU_BUG 1)
endif()

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

  # try to detect if our MPI implementation supports MPI-IO
  set(ABINIT_ENABLE_MPI_IO AUTO CACHE STRING "Enable to have MPI I/O support (default: no)")
  set_property(CACHE ABINIT_ENABLE_MPI_IO PROPERTY STRINGS AUTO ON OFF)
  if(ABINIT_ENABLE_MPI_IO STREQUAL "AUTO" OR (ABINIT_ENABLE_MPI_IO))
    try_compile(MPI_HAS_IO
      ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_mpi_io.F90
      LINK_LIBRARIES MPI::MPI_Fortran)
    if(MPI_HAS_IO)
      set(ABINIT_ENABLE_MPI_IO ON CACHE STRING "Enable to have MPI I/O support (default: no)" FORCE)
    else()
      if(ABINIT_ENABLE_MPI_IO)
        message(FATAL_ERROR "ABINIT_ENABLE_MPI_IO option is activated but the MPI library does not support MPI-IO!")
      endif()
      set(ABINIT_ENABLE_MPI_IO OFF CACHE STRING "Enable to have MPI I/O support (default: no)" FORCE)
    endif()
  endif()
  if(ABINIT_ENABLE_MPI_IO)
    set(HAVE_MPI_IO 1)
  endif()

  # try to detect if our MPI implementation supports MPI_INPLACE
  set(ABINIT_ENABLE_MPI_INPLACE AUTO CACHE STRING "Enable the use of MPI_IN_PLACE (default: auto)")
  set_property(CACHE ABINIT_ENABLE_MPI_INPLACE PROPERTY STRINGS AUTO ON OFF)
  if((ABINIT_ENABLE_MPI_INPLACE STREQUAL "AUTO") OR (ABINIT_ENABLE_MPI_INPLACE))
    try_compile(MPI_HAS_INPLACE
      ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_mpi_inplace.F90
      LINK_LIBRARIES MPI::MPI_Fortran)
    if(MPI_HAS_INPLACE)
      set(ABINIT_ENABLE_MPI_INPLACE ON CACHE STRING "Enable the use of MPI_IN_PLACE (default: auto)" FORCE)
    else()
      if(ABINIT_ENABLE_MPI_INPLACE)
        message(FATAL_ERROR "ABINIT_ENABLE_MPI_INPLACE option is activated but the MPI library does not support MPI-INPLACE!")
      endif()
      set(ABINIT_ENABLE_MPI_INPLACE OFF CACHE STRING "Enable the use of MPI_IN_PLACE (default: auto)" FORCE)
    endif()
  endif()
  if(ABINIT_ENABLE_MPI_INPLACE)
    set(HAVE_MPI2_INPLACE 1)
  endif()

  # try to detect if our MPI implementation has buggy interfaces (mishandling scalars)
  set(ABINIT_ENABLE_MPI_INTERFACES_BUGFIX AUTO CACHE STRING "Enable a workaround for buggy MPI interfaces (mishandling scalars) (default: auto)")
  set_property(CACHE ABINIT_ENABLE_MPI_INTERFACES_BUGFIX PROPERTY STRINGS AUTO ON OFF)
  if((ABINIT_ENABLE_MPI_INTERFACES_BUGFIX STREQUAL "AUTO") OR (NOT ABINIT_ENABLE_MPI_INTERFACES_BUGFIX))
    try_compile(MPI_NOT_BUGGY
      ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_mpi_buggy_interfaces.F90
      LINK_LIBRARIES MPI::MPI_Fortran)
    if(MPI_NOT_BUGGY)
      set(ABINIT_ENABLE_MPI_INTERFACES_BUGFIX OFF CACHE STRING "Enable a workaround for buggy MPI interfaces (mishandling scalars) (default: auto)" FORCE)
    else()
      if(NOT ABINIT_ENABLE_MPI_INTERFACES_BUGFIX)
        message(FATAL_ERROR "ABINIT_ENABLE_MPI_INTERFACES_BUGFIX option is deactivated but the MPI library has buggy interfaces!")
      endif()
      set(ABINIT_ENABLE_MPI_INTERFACES_BUGFIX ON CACHE STRING "Enable a workaround for buggy MPI interfaces (mishandling scalars) (default: auto)" FORCE)
    endif()
  endif()
  if(ABINIT_ENABLE_MPI_INTERFACES_BUGFIX)
    set(HAVE_MPI_BUGGY_INTERFACES 1)
  endif()

  set(ABINIT_ENABLE_GPU_AWARE_MPI FALSE)
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
          set(ABINIT_ENABLE_GPU_AWARE_MPI TRUE)
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

  endif(ABINIT_ENABLE_GPU_CUDA)

  # unconditionally set HAVE_GPU_MPI if the user say so
  # this may be need because our detection above is not bullet proof
  if (ABINIT_ENFORCE_GPU_AWARE_MPI)
    set(HAVE_GPU_MPI 1)
    set(ABINIT_ENABLE_GPU_AWARE_MPI TRUE)
  endif()


endif(MPI_FOUND)

#
# OpenMP
#
find_package(OpenMP)

if(OpenMP_FOUND)
  set(HAVE_OPENMP 1)

  # check for OpenMP collapse clause support
  try_compile(HAVE_OMP_COLLAPSE_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_omp_collapse.F90
    LINK_LIBRARIES OpenMP::OpenMP_Fortran)
  if(HAVE_OMP_COLLAPSE_BOOL)
    set(HAVE_OMP_COLLAPSE 1)
  endif()

  if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "Cray")
    # Fix an inconsistency causing compilation problems despite '-fopenmp' and '-h omp' are supposed to be synonymous !!
    # Original flag for OpenMP in Cray FTN is : -h omp
    set(OpenMP_Fortran_FLAGS "-fopenmp")
    add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>>:-fopenmp>")
    add_link_options("$<$<AND:$<LINK_LANGUAGE:Fortran>>:-fopenmp>")
  endif()

  if (ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP)
    # check for OpenMP GPU offload support (OpenMP v5 onwards)
    try_compile(HAVE_OPENMP_OFFLOAD_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_openmp_offload.F90
      LINK_OPTIONS ${OpenMP_Fortran_FLAGS}
      LINK_LIBRARIES OpenMP::OpenMP_Fortran)

    if(HAVE_OPENMP_OFFLOAD_BOOL)
      set(HAVE_OPENMP_OFFLOAD 1)

      include(manage_openmp_target)
      #FIXME Experimental, only tested with NVHPC through Intel IFX, AMD Flang and GFortran might be handled.
      # Examples:
      #get_offload_flags(ARCH   gfx1030 RESULT OFFLOAD_FLAGS)  for AMD GPUs    (Flang > v16.0)
      #get_offload_flags(TARGET 80      RESULT OFFLOAD_FLAGS)  for NVIDIA GPUs (NVHPC)
      #get_offload_flags(ARCH   sm_80   RESULT OFFLOAD_FLAGS)  for NVIDIA GPUs (NVHPC > v23.x and Flang > v16.0)
      #get_offload_flags(TARGET spir64  RESULT OFFLOAD_FLAGS)  for Intel GPUs  (IntelLLVM a.k.a IFX)

      if (ABINIT_ENABLE_GPU_CUDA)
        get_offload_flags(TARGET ${CMAKE_CUDA_ARCHITECTURES}    RESULT OPENMP_OFFLOAD_FLAGS)
      elseif(ABINIT_ENABLE_GPU_HIP)
        get_offload_flags(ARCH   ${CMAKE_HIP_ARCHITECTURES}     RESULT OPENMP_OFFLOAD_FLAGS)
      endif()

      # check whether compiler implements omp_get_mapped_ptr
      try_compile(HAVE_OPENMP_GET_MAPPED_PTR_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_openmp_get_mapped_ptr.F90
      LINK_LIBRARIES OpenMP::OpenMP_Fortran)
      if(HAVE_OPENMP_GET_MAPPED_PTR_BOOL)
        set(HAVE_OPENMP_GET_MAPPED_PTR 1)
      endif()

      # some compilers do not correctly handle structured types in OMP OFFLOAD directives
      # at present only NVHPC can do it
      if ("${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC")
        set(HAVE_OPENMP_OFFLOAD_DATASTRUCTURE 1)
      endif()

    else()
      message(STATUS "OpenMP offload isn't supported and disabled.")
    endif()
  endif(ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP)

  if(ABINIT_ENABLE_NVIDIA_UNIFIED_MEM)
    if (NOT "${CMAKE_Fortran_COMPILER_ID}" MATCHES "NVHPC")
      message(FATAL_ERROR "ABINIT_ENABLE_NVIDIA_UNIFIED_MEM requires NVHPC as compiler !")
    endif()
    if(NOT HAVE_OPENMP_OFFLOAD_BOOL)
      message(FATAL_ERROR "Provided ABINIT_ENABLE_NVIDIA_UNIFIED_MEM but OpenMP offload isn't working !")
    endif()
    if(CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "24.3.0")
      message(FATAL_ERROR "Your NVHPC version is ${CMAKE_Fortran_COMPILER_VERSION}. Unified memory support starts from 24.3.0 !")
    else()
      set(UNIFIED_MEM_FLAG "-gpu=mem:unified")
      # Use older flag for NVHPC 24.3, deprecated in newer versions
      if(CMAKE_Fortran_COMPILER_VERSION VERSION_EQUAL "24.3.0")
        set(UNIFIED_MEM_FLAG "-gpu=unified")
      endif()
      list(APPEND OPENMP_OFFLOAD_FLAGS "${UNIFIED_MEM_FLAG}")
      add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:C>>:${UNIFIED_MEM_FLAG}>")
      add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:CXX>>:${UNIFIED_MEM_FLAG}>")
      add_compile_options("$<$<AND:$<COMPILE_LANGUAGE:Fortran>>:${UNIFIED_MEM_FLAG}>")
      add_link_options("$<$<AND:$<LINK_LANGUAGE:Fortran>>:${UNIFIED_MEM_FLAG}>")
    endif()
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
    message(STATUS "ELPA found via pkg-config")
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
  message(STATUS "Libxc found via pkg-config")
  set(Libxc_FOUND TRUE)
  #set(Libxc_VERSION ${LibxcPkg_VERSION})
  #set(Libxc_INCLUDE_DIRS ${LibxcPkg_INCLUDE_DIRS})
  add_library(abinit::libxc ALIAS PkgConfig::ABINIT_LIBXC)
else()
  find_package(Libxc REQUIRED)
  if(Libxc_FOUND)
    message(STATUS "Libxc found via cmake target")
  else()
    message(SEND_ERROR "Libxc NOT FOUND")
  endif()
endif()

if(Libxc_FOUND)
  set(HAVE_LIBXC 1)
endif()

#
# FFT FLAVOR :
#
if(ABINIT_FFT_FLAVOR STREQUAL "FFTW3")
  # use regular fftw3 library
  pkg_check_modules(ABINIT_FFTW3 QUIET IMPORTED_TARGET fftw3 fftw3f)
  if(ABINIT_FFTW3_FOUND)
    message(STATUS "Selected fft flavor is FFTW3")
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
    message(STATUS "Selected fft flavor is MKL_DFTI")
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
    message(STATUS "pfft found via pkg-config")
    add_library(abinit::fft ALIAS PkgConfig::ABINIT_PFFT)
    set(HAVE_PFFT 1)
  else()
    message(FATAL_ERROR "pfft NOT FOUND, selected another flavor or adjust env variable PKG_CONFIG_PATH")
  endif()

  message(STATUS "Selected fft flavor is PFFT")

elseif(ABINIT_FFT_FLAVOR STREQUAL "GOEDECKER")

  message(STATUS "Selected fft flavor is GOEDECKER (internal implementation, no link flags required).")
  # WARNING: abinit::fft alias is not defined

endif()

#
# LEVMAR : use LevMar with cmake build-system
# see : https://github.com/pkestene/levmar-cmake
#

find_library(LEVMAR_LIBRARY
  NAMES levmar
  PATHS ${LEVMAR_ROOT}/lib ${LEVMAR_ROOT}/lib64 /opt/levmar/lib
  DOC "Location of the levmar library"
)
check_include_file("levmar.h" LEVMAR_HEADER_FOUND)
if(LEVMAR_LIBRARY AND LEVMAR_HEADER_FOUND)
  target_link_libraries(abinit ${LEVMAR_LIBRARY})
  try_compile(LEVMAR_FOUND
    ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_levmar.c)
endif()
if (LEVMAR_FOUND)
  message(STATUS "levmar found: ${LEVMAR_LIBRARY}")
  set(HAVE_LEVMAR 1)
else()
  message(STATUS "levmar not found")
endif()

#
# LIBPSML : TODO
#
pkg_check_modules(ABINIT_LIBPSML QUIET IMPORTED_TARGET libpsml)
if(ABINIT_LIBPSML_FOUND)
  message(STATUS "libpsml found via pkg-config")
  add_library(abinit::libpsml ALIAS PkgConfig::ABINIT_LIBPSML)
  set(HAVE_LIBPSML 1)
else()
  message(STATUS "libpsml NOT FOUND")
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
  message(STATUS "HDF5 found via find_package")
  add_library(abinit::hdf5 ALIAS HDF5::HDF5)
  set(ABINIT_HDF5_FOUND TRUE)
else()
  pkg_check_modules(ABINIT_HDF5 QUIET IMPORTED_TARGET hdf5)
  if(NOT ABINIT_HDF5_FOUND)
    # Second try for Cray parallel version
    pkg_check_modules(ABINIT_HDF5 QUIET IMPORTED_TARGET hdf5_parallel)
  endif()
  if(ABINIT_HDF5_FOUND)
    add_library(abinit::hdf5 ALIAS PkgConfig::ABINIT_HDF5)
    message(STATUS "HDF5 found via pkg-config")
  endif()
endif()

if(ABINIT_HDF5_FOUND)
  set(HAVE_HDF5 1)

  if(HDF5_IS_PARALLEL)
    set(HAVE_HDF5_MPI 1)
  endif()
else()
  message(SEND_ERROR "HDF5 NOT FOUND")
endif()


#
# NetCDF_C / NetCDF_Fortran
#
pkg_check_modules(ABINIT_NETCDF QUIET IMPORTED_TARGET netcdf)
# Extra check for Cray env...
if(NOT ABINIT_NETCDF_FOUND)
  pkg_check_modules(ABINIT_NETCDF QUIET IMPORTED_TARGET netcdf_parallel)
endif()
if(ABINIT_NETCDF_FOUND)
  message(STATUS "netcdf found via pkg-config")
  add_library(abinit::netcdf ALIAS PkgConfig::ABINIT_NETCDF)
else()
  message(SEND_ERROR "netcdf NOT FOUND")
endif()

pkg_check_modules(ABINIT_NETCDF_FORTRAN QUIET IMPORTED_TARGET GLOBAL netcdf-fortran)
#Extra check for Cray env...
if(NOT ABINIT_NETCDF_FORTRAN_FOUND)
  pkg_check_modules(ABINIT_NETCDF_FORTRAN QUIET IMPORTED_TARGET GLOBAL netcdf-fortran_parallel)
endif()
if(ABINIT_NETCDF_FORTRAN_FOUND)
  message(STATUS "netcdf-fortran found via pkg-config : ${ABINIT_NETCDF_FORTRAN_PREFIX}")

  # kind of fix to make sure path PREFIX/include (where netcdf.mod is located)
  # is append in include dirs used in abinit::netcdf-fortran
  pkg_get_variable(ABINIT_NETCDF_FORTRAN_PREFIX netcdf-fortran prefix)
  if("${ABINIT_NETCDF_FORTRAN_PREFIX}" STREQUAL "")
    # Cray flavor of netcdf
    pkg_get_variable(ABINIT_NETCDF_FORTRAN_PREFIX netcdf-fortran_parallel netcdf_prefix)
  endif()
  get_target_property(ABINIT_NETCDF_FORTRAN_INC_TMP PkgConfig::ABINIT_NETCDF_FORTRAN INTERFACE_INCLUDE_DIRECTORIES)
  list(APPEND ABINIT_NETCDF_FORTRAN_INC_TMP "${ABINIT_NETCDF_FORTRAN_PREFIX}/include")
  set_target_properties(PkgConfig::ABINIT_NETCDF_FORTRAN PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${ABINIT_NETCDF_FORTRAN_INC_TMP}")

  add_library(abinit::netcdf-fortran ALIAS PkgConfig::ABINIT_NETCDF_FORTRAN)

else()
  message(SEND_ERROR "netcdf-fortran NOT FOUND")
endif()

if (MPI_FOUND)
  pkg_check_modules(ABINIT_NETCDF_MPI QUIET IMPORTED_TARGET netcdf-mpi)
  if(ABINIT_NETCDF_MPI_FOUND)
    message(STATUS "netcdf-mpi found via pkg-config")
    add_library(abinit::netcdf-mpi ALIAS PkgConfig::ABINIT_NETCDF_MPI)
    #set(HAVE_NETCDF_MPI 1) => see generate_config.h
  else()
    message(STATUS "netcdf-mpi NOT FOUND")
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
    message(STATUS "f90papi.h not found")
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

  message(STATUS "papi found via pkg-config")
else()
  message(STATUS "papi NOT FOUND")
endif()

#
# TRIQS : TODO improve
# to help cmake detect TRIQS:
# - set env var or cmake var TRIQS_ROOT to the root path of triqs install directory
# - set env var CMAKE_PREFIX_PATH (usually by a modulefile on a supercomputer)
#
if (ABINIT_ENABLE_TRIQS)
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
    message(SEND_ERROR "ABINIT_ENABLE_TRIQS activated but TRIQS has not been found!")
  endif()
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
  message(STATUS "xml2 found via pkg-config")
  add_library(abinit::xml2 ALIAS PkgConfig::ABINIT_XML2)
  set(HAVE_XML 1)
else()
  message(STATUS "xml2 NOT FOUND")
endif()

#
# XMLF90
#
pkg_check_modules(ABINIT_XMLF90 QUIET IMPORTED_TARGET xmlf90)
if(ABINIT_XMLF90_FOUND)
  message(STATUS "xmlf90 found via pkg-config")
  add_library(abinit::xmlf90 ALIAS PkgConfig::ABINIT_XMLF90)
  set(HAVE_XMLF90 1)
else()
  message(STATUS "xmlf90 NOT FOUND")
endif()


#
# AtomPAW
#
pkg_check_modules(ABINIT_ATOMPAW QUIET IMPORTED_TARGET atompaw)
if(ABINIT_ATOMPAW_FOUND)
  message(STATUS "atompaw found via pkg-config")
  add_library(abinit::atompaw ALIAS PkgConfig::ABINIT_ATOMPAW)
  set(HAVE_ATOMPAW 1)
else()
  message(STATUS "atompaw NOT FOUND")
endif()

#
# BigDFT
#
pkg_check_modules(ABINIT_BIGDFT QUIET IMPORTED_TARGET bigdft)
if(ABINIT_BIGDFT_FOUND)
  message(STATUS "bigdft found via pkg-config")
  add_library(abinit::bigdft ALIAS PkgConfig::ABINIT_BIGDFT)
  set(HAVE_BIGDFT 1)
else()
  message(STATUS "bigdft NOT FOUND")
endif()

#
# Generate config.h
#
include(generate_config_h)

##################### PRINT CONFIGURE STATUS ######################
message(STATUS "")
message(STATUS "=====================================================")
message(STATUS "== ${PROJECT_NAME} build configuration:")
message(STATUS "== ${PROJECT_NAME} version : ${PROJECT_VERSION}")
message(STATUS "=====================================================")
message(STATUS "")

message(STATUS "CMake config:")
message(STATUS "  > CMake version          : ${CMAKE_VERSION}")
if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "  > CMake build type       : NOT SET!")
else()
  message(STATUS "  > CMake build type       : ${CMAKE_BUILD_TYPE}")
endif()
message(STATUS "  > CMake install prefix   : ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  > CMake system processor : ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "  > CMake system name (OS) : ${CMAKE_SYSTEM_NAME}")

message(STATUS "Fortran compiler:")
message(STATUS "  > Fortran compiler Id      : ${CMAKE_Fortran_COMPILER_ID}")
message(STATUS "  > Fortran compiler version : ${CMAKE_Fortran_COMPILER_VERSION}")
message(STATUS "  > Fortran compiler exe     : ${CMAKE_Fortran_COMPILER}")
message(STATUS "  > Fortran flags            : ${CMAKE_Fortran_FLAGS}")
message(STATUS "  > Fortran compiler wrapper : ${CMAKE_Fortran_COMPILER_WRAPPER}")

if (MPI_FOUND)
  message(STATUS "MPI found:")
    # note : make sure you correctly set C, C++ and Fortran compilers
    # through variable CC, CXX and FC
  message(STATUS "  > MPI standard version : ${MPI_Fortran_VERSION}")
  message(STATUS "  > MPI_VENDOR           : ${MPI_VENDOR}")
  message(STATUS "  > MPI library version  : ${MPI_Fortran_LIBRARY_VERSION_STRING}")
  message(STATUS "  > MPI fortran compiler : ${MPI_Fortran_COMPILER}")
  message(STATUS "  > MPI headers          : ${MPI_Fortran_INCLUDE_DIRS}")
  message(STATUS "  > MPI libraries        : ${MPI_Fortran_LIBRARIES}")
 #message(STATUS "  > MPI link flags       : ${MPI_Fortran_LINK_FLAGS}")
  message(STATUS "  > MPI supports MPI-IO  : ${ABINIT_ENABLE_MPI_IO}")
  message(STATUS "  > MPI supports INPLACE : ${ABINIT_ENABLE_MPI_INPLACE}")
  message(STATUS "  > MPI has buggy interf.: ${ABINIT_ENABLE_MPI_INTERFACES_BUGFIX}")
else()
  message(STATUS "MPI not found!")
endif()

if (OpenMP_FOUND)
  message(STATUS "OpenMP found:")
  message(STATUS "  > OpenMP version : ${OpenMP_Fortran_VERSION}")
else()
  message(STATUS "OpenMP not found!")
endif()

if (BLAS_FOUND)
  message(STATUS "BLAS found:")
  message(STATUS "  > BLA_VENDOR        : ${BLA_VENDOR}")
  message(STATUS "  > BLAS_LINKER_FLAGS : ${BLAS_LINKER_FLAGS}")
  message(STATUS "  > BLAS_LIBRARIES    : ${BLAS_LIBRARIES}")
  message(STATUS "  > ZDOT workaround   : ${ABINIT_ENABLE_ZDOT_BUGFIX}")
else()
  message(STATUS "BLAS not found!")
  message(STATUS "  > Please recheck your environment variables!")
endif()

if (LAPACK_FOUND)
  message(STATUS "LAPACK found:")
  message(STATUS "  > BLA_VENDOR          : ${BLA_VENDOR}")
  message(STATUS "  > LAPACK_LINKER_FLAGS : ${LAPACK_LINKER_FLAGS}")
  message(STATUS "  > LAPACK_LIBRARIES    : ${LAPACK_LIBRARIES}")
else()
  message(STATUS "LAPACK not found!")
  message(STATUS "  > Please recheck your environment variables!")
endif()

message(STATUS "LAPACK config:")
message(STATUS "    BLA_VENDOR : ${BLA_VENDOR}")
if (LAPACK_FOUND)
  message(STATUS "    LAPACK_LINKER_FLAGS : ${LAPACK_LINKER_FLAGS}")
  message(STATUS "    LAPACK_LIBRARIES    : ${LAPACK_LIBRARIES}")
else()
  message(STATUS "    LAPACK not found ! Please recheck your environment variables.")
endif()

if (ABINIT_SCALAPACK_ENABLED)
  message(STATUS "Scalapack enabled:")
  message(STATUS "  > ABINIT_SCALAPACK_FLAVOR : ${ABINIT_SCALAPACK_FLAVOR}")
  if(ABINIT_SCALAPACK_FLAVOR MATCHES "MKL")
    message(STATUS "  > MKL found                  : ${MKL_FOUND}")
    message(STATUS "  > MKL_BLAS_FLAVOR            : ${MKL_BLAS_FLAVOR}")
    message(STATUS "  > MKL_BLAS_FLAVOR_FOUND      : ${MKL_BLAS_FLAVOR_FOUND}")
    message(STATUS "  > MKL blas headers           : ${MKL_BLAS_HEADERS}")
    message(STATUS "  > MKL blas libraries         : ${MKL_BLAS_LIBRARIES}")
    message(STATUS "  > MKL_SCALAPACK_FLAVOR       : ${MKL_SCALAPACK_FLAVOR}")
    message(STATUS "  > MKL_SCALAPACK_FLAVOR_FOUND : ${MKL_SCALAPACK_FLAVOR_FOUND}")
    message(STATUS "  > MKL scalapack headers      : ${MKL_SCALAPACK_HEADERS}")
    message(STATUS "  > MKL scalapack libraries    : ${MKL_SCALAPACK_LIBRARIES}")
  elseif(ABINIT_SCALAPACK_FLAVOR MATCHES "NETLIB")
    message(STATUS "  > ScaLapack found                : ${scalapack_FOUND}")
    message(STATUS "  > Scalapack dir                  : ${scalapack_DIR}")
    message(STATUS "  > Scalapack found via pkg-config : ${scalapack_found_using_pkg_config}")
    message(STATUS "  > Scalapack found via cmake cfg  : ${scalapack_found_using_cmake_target}")
  endif()
else()
  message(STATUS "Scalapack not enabled!")
endif()

if (ABINIT_ELPA_ENABLED)
  message(STATUS "ELPA enabled:")
  message(STATUS "  > ELPA found   : ${ABINIT_ELPA_FOUND}")
  message(STATUS "  > ELPA version : ${ABINIT_ELPA_VERSION}")
  message(STATUS "  > ELPA cflags  : ${ABINIT_ELPA_CFLAGS}")
  message(STATUS "  > ELPA ldflags : ${ABINIT_ELPA_LDFLAGS}")
else()
  message(STATUS "ELPA not enabled!")
endif()

if (Libxc_FOUND)
  message(STATUS "Libxc found:")
  message(STATUS "  > Libxc dir          : ${Libxc_DIR}")
  message(STATUS "  > Libxc version      : ${Libxc_VERSION}")
  message(STATUS "  > Libxc headers      : ${Libxc_INCLUDE_DIRS}")
  message(STATUS "  > Libxc libraries    : ${Libxc_LIBRARIES}")
else()
  message(STATUS "Libxc not found!")
endif()

if (ABINIT_NETCDF_FOUND)
  message(STATUS "NetCDF-c found:")
  message(STATUS "  > netcdf-c prefix       : ${ABINIT_NETCDF_PREFIX}")
  message(STATUS "  > netcdf-c version      : ${ABINIT_NETCDF_VERSION}")
  message(STATUS "  > netcdf-c inc dirs     : ${ABINIT_NETCDF_INCLUDE_DIRS}")
  message(STATUS "  > netcdf-c cflags       : ${ABINIT_NETCDF_CFLAGS}")
  message(STATUS "  > netcdf-c cflags other : ${ABINIT_NETCDF_CFLAGS_OTHER}")
  message(STATUS "  > netcdf-c ldflags      : ${ABINIT_NETCDF_LDFLAGS}")
else()
  message(STATUS "NetCDF-c not found!")
endif()

if (ABINIT_NETCDF_FORTRAN_FOUND)
  message(STATUS "NetCDF-Fortran found:")
  message(STATUS "  > netcdf-fortran prefix       : ${ABINIT_NETCDF_FORTRAN_PREFIX}")
  message(STATUS "  > netcdf-fortran version      : ${ABINIT_NETCDF_FORTRAN_VERSION}")
  message(STATUS "  > netcdf-fortran inc dirs     : ${ABINIT_NETCDF_FORTRAN_INCLUDE_DIRS}")
  message(STATUS "  > netcdf-fortran cflags       : ${ABINIT_NETCDF_FORTRAN_CFLAGS}")
  message(STATUS "  > netcdf-fortran cflags other : ${ABINIT_NETCDF_FORTRAN_CFLAGS_OTHER}")
  message(STATUS "  > netcdf-fortran ldflags      : ${ABINIT_NETCDF_FORTRAN_LDFLAGS}")
  message(STATUS "  > netcdf C/MPI support        : ${HAVE_NETCDF_MPI}")
  message(STATUS "  > netcdf fortran/MPI support  : ${HAVE_NETCDF_FORTRAN_MPI}")
else()
  message(STATUS "NetCDF-Fortran not found!")
endif()

#
# TODO PK : add the possibility to detect a version of libpaw already installed
#
if (HAVE_LIBPAW_ABINIT)
  message(STATUS "Internal libPAW selected!")
else()
  message(STATUS "External libPAW selected (not possible for now)!")
endif()

if (ABINIT_FFTW3_FOUND)
  message(STATUS "FFTW3 found:")
  message(STATUS "  > fftw3 version       : ${ABINIT_FFTW3_VERSION}")
  message(STATUS "  > fftw3 cflags        : ${ABINIT_FFTW3_CFLAGS}")
  message(STATUS "  > fftw3 ldflags       : ${ABINIT_FFTW3_LDFLAGS}")
else()
  message(STATUS "FFTW3 not found!")
endif()

if (ABINIT_PFFT_FOUND)
  message(STATUS "PFFT found:")
  message(STATUS "  > pfft version : ${ABINIT_PFFT_VERSION}")
  message(STATUS "  > pfft cflags  : ${ABINIT_PFFT_CFLAGS}")
  message(STATUS "  > pfft ldflags : ${ABINIT_PFFT_LDFLAGS}")
else()
  message(STATUS "PFFT not found!")
endif()

if (ABINIT_WANNIER90_WANTED)
  message(STATUS "Wannier90 enabled:")
  message(STATUS "  > Wannier90 builtin     : ${ABINIT_WANNIER90_BUILTIN}")
  message(STATUS "  > Wannier90 found       : ${ABINIT_WANNIER90_FOUND}")
  message(STATUS "  > Wannier90 root        : ${WANNIER_ROOT}")
  message(STATUS "  > Wannier90 version     : ${ABINIT_WANNIER90_VERSION}")
  message(STATUS "  > Wannier90 include dir : ${ABINIT_WANNIER90_INCLUDE_DIR}")
  message(STATUS "  > Wannier90 library     : ${ABINIT_WANNIER90_LIBRARY}")
else()
  message(STATUS "Wannier90 not enabled!")
endif()

if (ABINIT_XML2_FOUND)
  message(STATUS "xml-2.0 found:")
  message(STATUS "  > xml-2.0 version : ${ABINIT_XML2_VERSION}")
  message(STATUS "  > xml-2.0 cflags  : ${ABINIT_XML2_CFLAGS}")
  message(STATUS "  > xml-2.0 ldflags : ${ABINIT_XML2_LDFLAGS}")
else()
  message(STATUS "xml-2.0 not found!")
endif()

if (ABINIT_XMLF90_FOUND)
  message(STATUS "xmlf90 found:")
  message(STATUS "  > xmlf90 version : ${ABINIT_XMLF90_VERSION}")
  message(STATUS "  > xmlf90 cflags  : ${ABINIT_XMLF90_CFLAGS}")
  message(STATUS "  > xmlf90 ldflags : ${ABINIT_XMLF90_LDFLAGS}")
else()
  message(STATUS "xmlf90 not found!")
endif()

if (ABINIT_LIBPSML_FOUND)
  message(STATUS "libpsml found:")
  message(STATUS "  > libpsml version : ${ABINIT_LIBPSML_VERSION}")
  message(STATUS "  > libpsml cflags  : ${ABINIT_LIBPSML_CFLAGS}")
  message(STATUS "  > libpsml ldflags : ${ABINIT_LIBPSML_LDFLAGS}")
else()
  message(STATUS "libpsml not found!")
endif()

if (ABINIT_ATOMPAW_FOUND)
  message(STATUS "AtomPAW found:")
  message(STATUS "  > AtomPAW version : ${ABINIT_ATOMPAW_VERSION}")
  message(STATUS "  > AtomPAW cflags  : ${ABINIT_ATOMPAW_CFLAGS}")
  message(STATUS "  > AtomPAW ldflags : ${ABINIT_ATOMPAW_LDFLAGS}")
else()
  message(STATUS "AtomPAW not found!")
endif()

if (ABINIT_BIGDFT_FOUND)
  message(STATUS "BigDFT found:")
  message(STATUS "  > BigDFT version : ${ABINIT_BIGDFT_VERSION}")
  message(STATUS "  > BigDFT cflags  : ${ABINIT_BIGDFT_CFLAGS}")
  message(STATUS "  > BigDFT ldflags : ${ABINIT_BIGDFT_LDFLAGS}")
else()
  message(STATUS "BigDFT not found!")
endif()

if (LEVMAR_FOUND)
  message(STATUS "levmar found:")
  message(STATUS "  > levmar path : ${LEVMAR_LIBRARY}")
else()
  message(STATUS "levmar not found!")
endif()

if (TRIQS_FOUND)
  message(STATUS "TRIQS found:")
  message(STATUS "  > TRIQS version        : ${TRIQS_VERSION}")
  message(STATUS "  > TRIQS root           : ${TRIQS_ROOT}")
  message(STATUS "  > TRIQS python support : ${TRIQS_WITH_PYTHON_SUPPORT}")
  message(STATUS "  > TRIQS include dir    : ${TRIQS_INCLUDE_DIRECTORIES}")
  message(STATUS "  > TRIQS link libraries : ${TRIQS_LINK_LIBRARIES}")
# message(STATUS "  > TRIQS compile definitions : ${TRIQS_COMPILE_DEFINITIONS}")
# message(STATUS "  > TRIQS compile features    : ${TRIQS_COMPILE_FEATURES}")
# message(STATUS "  > TRIQS compile options     : ${TRIQS_COMPILE_OPTIONS}")
else()
  message(STATUS "TRIQS not found!")
endif()

if (ABINIT_PAPI_FOUND)
  message(STATUS "papi found:")
  message(STATUS "  > papi version : ${ABINIT_PAPI_VERSION}")
  message(STATUS "  > papi cflags  : ${ABINIT_PAPI_CFLAGS}")
  message(STATUS "  > papi ldflags : ${ABINIT_PAPI_LDFLAGS}")
else()
  message(STATUS "papi not found!")
endif()

if (ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP OR HAVE_OPENMP_OFFLOAD OR ABINIT_KOKKOS_WANTED OR ABINIT_YAKL_WANTED)
  message(STATUS "GPU enabled:")
  if (ABINIT_ENABLE_GPU_CUDA)
    message(STATUS "  > CUDA enabled:")
    message(STATUS "    - CUDA compiler ID      : ${CMAKE_CUDA_COMPILER_ID}")
    message(STATUS "    - CUDA compiler version : ${CMAKE_CUDA_COMPILER_VERSION}")
    message(STATUS "    - C++ compiler ID       : ${CMAKE_CXX_COMPILER_ID}")
    message(STATUS "    - C++ compiler version  : ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "    - CUDA Compiler         : ${CMAKE_CUDA_COMPILER}")
    message(STATUS "    - CUDA Compiler exec    : ${CUDA_NVCC_EXECUTABLE}")
    message(STATUS "    - CUDA Compile flags    : ${CMAKE_CUDA_FLAGS}")
    message(STATUS "    - CUDA toolkit found    : ${CUDAToolkit_FOUND}")
    message(STATUS "    - CUDA toolkit version  : ${CUDAToolkit_VERSION}")
    message(STATUS "    - CUDA toolkit nvcc     : ${CUDAToolkit_NVCC_EXECUTABLE}")
  else()
    message(STATUS "  > CUDA not enabled!")
  endif()
  if (ABINIT_ENABLE_GPU_HIP)
    message(STATUS "  > HIP enabled:")
    message(STATUS "    - HIP found   : ${HIP_FOUND}")
    message(STATUS "    - HIP version : ${HIP_VERSION}")
  else()
    message(STATUS "  > HIP not enabled!")
  endif()
  if (HAVE_OPENMP_OFFLOAD)
    message(STATUS "  > OpenMP offload enabled:")
    message(STATUS "    - OpenMP GPU offload flags : ${OPENMP_OFFLOAD_FLAGS}")
  else()
    message(STATUS "  > OpenMP offload not enabled!")
  endif()
  message(STATUS "  > GPU-aware MPI enabled : ${ABINIT_ENABLE_GPU_AWARE_MPI}")
  if (ABINIT_KOKKOS_WANTED)
    message(STATUS "  > Kokkos enabled:")
    if (ABINIT_KOKKOS_BUILD)
      message(STATUS "    - Building local Kokkos:")
      message(STATUS "      ABINIT_KOKKOS_BACKEND : ${ABINIT_KOKKOS_BACKEND}")
      message(STATUS "      Kokkos_ENABLE_HWLOC   : ${Kokkos_ENABLE_HWLOC}")
      message(STATUS "      Kokkos_ENABLE_OPENMP  : ${Kokkos_ENABLE_OPENMP}")
      message(STATUS "      Kokkos_ENABLE_CUDA    : ${Kokkos_ENABLE_CUDA}")
    else()
      message(STATUS "    - Using installed Kokkos:")
      message(STATUS "      Kokkos dir             : ${Kokkos_DIR}")
      message(STATUS "      Kokkos version         : ${Kokkos_VERSION}")
      message(STATUS "      Kokkos_CXX_COMPILER    : ${Kokkos_CXX_COMPILER}")
      message(STATUS "      Kokkos_CXX_COMPILER_ID : ${Kokkos_CXX_COMPILER_ID}")
      message(STATUS "      Kokkos_CXX_STANDARD    : ${Kokkos_CXX_STANDARD}")
      message(STATUS "      Kokkos_OPTIONS         : ${Kokkos_OPTIONS}")
      message(STATUS "      Kokkos_TPLS            : ${Kokkos_TPLS}")
      message(STATUS "      Kokkos_DIRS            : ${Kokkos_DIR}")
    endif()
  else()
    message(STATUS "  > Kokkos not enabled!")
  endif()
  if (ABINIT_YAKL_WANTED)
    message(STATUS "  > YAKL enabled:")
    message(STATUS"     - YAKL_ARCH : ${YAKL_ARCH}")
  else()
    message(STATUS "  > YAKL not enabled!")
  endif()
else()
  message(STATUS "GPU not enabled!")
endif()

message(STATUS "ENV{PKG_CONFIG_PATH} was:")
message(STATUS "  > $ENV{PKG_CONFIG_PATH}")

message(STATUS "ENV{CMAKE_PREFIX_PATH} was:")
message(STATUS "  > $ENV{CMAKE_PREFIX_PATH}")

message(STATUS "=====================================================")
message(STATUS "")

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

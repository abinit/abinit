#
# CMake recipes
# https://github.com/eth-cscs/cmake-recipes
#
# Copyright (c) 2018-2019, ETH Zurich
# BSD 3-Clause License. All rights reserved.
#
# Author: Teodor Nikolov (tnikolov@cscs.ch)
#
# Modified to define variables
# - MKL_SCALAPACK_FLAVOR_DEFAULT
# - MKL_BLAS_FLAVOR_DEFAULT
# according to what is available on execution platform
#
#[=======================================================================[.rst:
mkl_setup
---------

The following conventions are used:

intel / INTEL  - Bindings for everything except GNU Fortran
gf / GF        - GNU Fortran bindings
seq / SEQ      - sequential MKL
omp / OMP      - threaded MKL with OpenMP back end
tbb / TBB      - threaded MKL with TBB back end
32bit / 32BIT  - MKL 32 bit integer interface (used most often)
64bit / 64BIT  - MKL 64 bit integer interface
mpich / MPICH  - MPICH / IntelMPI BLACS back end
ompi / OMPI    - OpenMPI BLACS back end
st / ST        - static libraries
dyn / DYN      - dynamic libraries

The module attempts to define a target for each MKL configuration. The
configuration will not be available if there are missing library files or a
missing dependency.

MKL Link line advisor:
  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

Note: Mixing GCC and Intel OpenMP backends is a bad idea.

Search variables
^^^^^^^^^^^^^^^^

``MKLROOT``
  Environment variable set to MKL's root directory

``MKL_ROOT``
  CMake variable set to MKL's root directory

Example usage
^^^^^^^^^^^^^

To Find MKL:

  find_package(MKL REQUIRED)
  include(mkl_setup)

To check if target is available:

  if (TARGET mkl::scalapack_mpich_intel_32bit_omp_dyn)
    ...
  endif()

To link to an available target (see list below):

  target_link_libraries(... mkl::scalapack_mpich_intel_32bit_omp_dyn)

Note: dependencies are handled for you (MPI, OpenMP, ...)

Imported targets
^^^^^^^^^^^^^^^^

MKL (BLAS, LAPACK, FFT) targets:

  mkl::mkl_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::mkl_intel_32bit_omp_dyn

BLACS targets:

  mkl::blacs_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::blacs_intel_mpich_32bit_seq_st

ScaLAPACK targets:

  mkl::scalapack_[mpich|ompi]_[gf|intel]_[32bit|64bit]_[seq|omp|tbb]_[st|dyn] e.g.

  mkl::scalapack_mpich_intel_64bit_omp_dyn

Result variables
^^^^^^^^^^^^^^^^

MKL_SETUP_DONE

Not supported
^^^^^^^^^^^^^

- F95 interfaces

#]=======================================================================]

cmake_minimum_required(VERSION 3.12)

# check if compatible compiler is found
if(CMAKE_C_COMPILER_LOADED OR
   CMAKE_CXX_COMPILER_LOADED OR
   CMAKE_Fortran_COMPILER_LOADED)
    set(_mkl_compiler_found TRUE)
else()
    set(_mkl_compiler_found FALSE)
endif()

# Dependencies
#
find_package(Threads)
find_package(MPI COMPONENTS CXX)
find_package(OpenMP COMPONENTS CXX)

# If MKL_ROOT is not set, set it via the env variable MKLROOT.
#
if(NOT DEFINED MKL_ROOT)
    set(MKL_ROOT $ENV{MKLROOT})
endif()

# Determine MKL's library folder
#
set(_mkl_libpath_suffix "lib/intel64")
if(CMAKE_SIZEOF_VOID_P EQUAL 4) # 32 bit
    set(_mkl_libpath_suffix "lib/ia32")
endif()

if(WIN32)
    list(APPEND _mkl_libpath_suffix "${_mkl_libpath_suffix}_win")
    set(_mkl_libname_prefix "")
    set(_mkl_shared_lib "_dll.lib")
    set(_mkl_static_lib ".lib")
elseif(APPLE)
    list(APPEND _mkl_libpath_suffix "${_mkl_libpath_suffix}_mac")
    set(_mkl_libname_prefix "lib")
    set(_mkl_shared_lib ".dylib")
    set(_mkl_static_lib ".a")
else() # LINUX
    list(APPEND _mkl_libpath_suffix "${_mkl_libpath_suffix}_lin")
    set(_mkl_libname_prefix "lib")
    set(_mkl_shared_lib ".so")
    set(_mkl_static_lib ".a")
endif()
set(_mkl_search_paths "${MKL_ROOT}"
                      "${MKL_ROOT}/lib"
                      "${MKL_ROOT}/mkl"
                      "${MKL_ROOT}/compiler")

# Functions: finds both static and shared MKL libraries
#
function(__mkl_find_library _varname _libname)
    find_library(${_varname}_DYN
          NAMES ${_mkl_libname_prefix}${_libname}${_mkl_shared_lib}
          HINTS ${_mkl_search_paths}
          PATH_SUFFIXES ${_mkl_libpath_suffix})
    mark_as_advanced(${_varname}_DYN)
    find_library(${_varname}_ST
          NAMES ${_mkl_libname_prefix}${_libname}${_mkl_static_lib}
          HINTS ${_mkl_search_paths}
          PATH_SUFFIXES ${_mkl_libpath_suffix})
    mark_as_advanced(${_varname}_ST)
endfunction()

# Find MKL headers
#
find_path(MKL_INCLUDE_DIR mkl.h
    HINTS ${MKL_ROOT}/include
          ${MKL_ROOT}/mkl/include)
mark_as_advanced(MKL_INCLUDE_DIR)

# Group flags for static libraries on Linux (GNU, PGI, ICC -> same linker)
#
if(UNIX AND NOT APPLE)
    set(_mkl_linker_pre_flags_ST "-Wl,--start-group")
    set(_mkl_linker_post_flags_ST "-Wl,--end-group")
endif()

# Core MKL
#
__mkl_find_library(MKL_CORE_LIB mkl_core)

# Interface
#
__mkl_find_library(MKL_INTERFACE_INTEL_32BIT_LIB mkl_intel_lp64)
__mkl_find_library(MKL_INTERFACE_INTEL_64BIT_LIB mkl_intel_ilp64)
if(NOT APPLE AND CMAKE_Fortran_COMPILER_LOADED
             AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    __mkl_find_library(MKL_INTERFACE_GF_32BIT_LIB mkl_gf_lp64)
    __mkl_find_library(MKL_INTERFACE_GF_64BIT_LIB mkl_gf_ilp64)
endif()

# Threading
#
__mkl_find_library(MKL_SEQ_LIB mkl_sequential)
if(NOT APPLE AND (CMAKE_C_COMPILER_ID STREQUAL "GNU" OR
                  CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR
                  CMAKE_Fortran_COMPILER_ID STREQUAL "GNU"))
    __mkl_find_library(MKL_OMP_LIB mkl_gnu_thread)
else()
    __mkl_find_library(MKL_OMP_LIB mkl_intel_thread)
endif()
__mkl_find_library(MKL_TBB_LIB mkl_tbb_thread)

# BLACS
#
if(APPLE)
    __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_mpich_lp64)
    __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_mpich_ilp64)
else()
    __mkl_find_library(MKL_BLACS_MPICH_32BIT_LIB mkl_blacs_intelmpi_lp64)
    __mkl_find_library(MKL_BLACS_MPICH_64BIT_LIB mkl_blacs_intelmpi_ilp64)
endif()
__mkl_find_library(MKL_BLACS_OMPI_32BIT_LIB mkl_blacs_openmpi_lp64)
__mkl_find_library(MKL_BLACS_OMPI_64BIT_LIB mkl_blacs_openmpi_ilp64)

# ScaLAPACK
#
__mkl_find_library(MKL_SCALAPACK_32BIT_LIB mkl_scalapack_lp64)
__mkl_find_library(MKL_SCALAPACK_64BIT_LIB mkl_scalapack_ilp64)

# Check if core libs were found
#
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MKL REQUIRED_VARS MKL_INCLUDE_DIR
                                                    Threads_FOUND
                                                    _mkl_compiler_found)

# Sequential has no threading dependency. There is currently no TBB module
# shipped with CMake. The dependency is not accounted for.
#
set(_mkl_dep_found_SEQ TRUE)
set(_mkl_dep_found_TBB TRUE)
if (TARGET OpenMP::OpenMP_CXX)
  set(_mkl_dep_OMP OpenMP::OpenMP_CXX)
  set(_mkl_dep_found_OMP TRUE)
endif()

####################################################

#
# make default values for
# - MKL_BLAS_FLAVOR_DEFAULT
# - MKL_SCALAPACK_FLAVOR_DEFAULT
#
if(NOT DEFINED MKL_XX_BIT)
  set(MKL_XX_BIT 32bit)
endif()

if( (NOT MKL_XX_BIT STREQUAL "32bit") AND (NOT MKL_XX_BIT STREQUAL "64bit"))
  message(FATAL_ERROR "MKL_XX_BIT invalid value; valid values are 32bit and 64bit")
endif()

if (NOT DEFINED MKL_THREADING_TYPE)
  if(OPENMP_FOUND)
    set (MKL_THREADING_TYPE "omp")
    #elseif(TBB_FOUND)
    #  set (scalapack_threading "tbb")
  else()
    set (MKL_THREADING_TYPE "seq")
  endif()
endif()

include(get_mpi_vendor)
if(MPI_VENDOR STREQUAL "OpenMPI")
  set(mpi_type "ompi")
else()
  set(mpi_type "mpich")
endif()

if (CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  set(compiler_type gf)
else()
  set(compiler_type intel)
endif()

set(MKL_BLAS_FLAVOR_DEFAULT mkl_${compiler_type}_${MKL_XX_BIT}_${MKL_THREADING_TYPE}_dyn)
set(MKL_SCALAPACK_FLAVOR_DEFAULT scalapack_${mpi_type}_intel_${MKL_XX_BIT}_${MKL_THREADING_TYPE}_dyn)

####################################################

#
# Define all blas, blacs and scalapack
#
foreach(_libtype "ST" "DYN")
    set(_mkl_core_lib ${MKL_CORE_LIB_${_libtype}})
    foreach(_bits "32BIT" "64BIT")
        set(_mkl_scalapack_lib ${MKL_SCALAPACK_${_bits}_LIB_${_libtype}})
        foreach(_iface "INTEL" "GF")
            set(_mkl_interface_lib ${MKL_INTERFACE_${_iface}_${_bits}_LIB_${_libtype}})
            foreach(_threading "SEQ" "OMP" "TBB")
                set(_mkl_threading_lib ${MKL_${_threading}_LIB_${_libtype}})

                string(TOLOWER "${_iface}_${_bits}_${_threading}_${_libtype}" _tgt_config)
                set(_mkl_tgt mkl::mkl_${_tgt_config})

                if(MKL_FOUND
                   AND _mkl_interface_lib
                   AND _mkl_threading_lib
                   AND _mkl_core_lib
                   AND _mkl_dep_found_${_threading}
                   AND NOT TARGET ${_mkl_tgt})
                    set(_mkl_libs "${_mkl_linker_pre_flags_${_threading}}"
                                  "${_mkl_interface_lib}"
                                  "${_mkl_threading_lib}"
                                  "${_mkl_core_lib}"
                                  "${_mkl_linker_post_flags_${_threading}}"
                                  "${_mkl_dep_${_threading}}"
                                  "Threads::Threads")
                    add_library(${_mkl_tgt} INTERFACE IMPORTED)
                    set_target_properties(${_mkl_tgt} PROPERTIES
                      INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
                      INTERFACE_LINK_LIBRARIES "${_mkl_libs}")
                    message("Create target for ${_mkl_tgt}")
                endif()

                foreach(_mpi_impl "MPICH" "OMPI")
                    set(_mkl_blacs_lib ${MKL_BLACS_${_mpi_impl}_${_bits}_LIB_${_libtype}})

                    string(TOLOWER "${_mpi_impl}_${_iface}_${_bits}_${_threading}_${_libtype}" _tgt_config)
                    set(_blacs_tgt mkl::blacs_${_tgt_config})
                    set(_scalapack_tgt mkl::scalapack_${_tgt_config})

                    if(_mkl_blacs_lib
                        AND TARGET ${_mkl_tgt}
                        AND TARGET MPI::MPI_CXX
                        AND NOT TARGET ${_blacs_tgt})
                      set(_blacs_libs "${_mkl_linker_pre_flags_${_libtype}}"
                                      "${_mkl_interface_lib}"
                                      "${_mkl_threading_lib}"
                                      "${_mkl_core_lib}"
                                      "${_mkl_blacs_lib}"
                                      "${_mkl_linker_post_flags_${_libtype}}"
                                      "MPI::MPI_CXX"
                                      "${_mkl_dep_${_threading}}"
                                      "Threads::Threads")
                       add_library(${_blacs_tgt} INTERFACE IMPORTED)
                       set_target_properties(${_blacs_tgt} PROPERTIES
                         INTERFACE_INCLUDE_DIRECTORIES "${MKL_INCLUDE_DIR}"
                         INTERFACE_LINK_LIBRARIES "${_blacs_libs}")
                       message("Create target for ${_blacs_tgt}")
                    endif()

                    if(_mkl_scalapack_lib
                        AND TARGET ${_blacs_tgt}
                        AND NOT TARGET ${_scalapack_tgt})
                      set(_scalapack_libs "${_mkl_scalapack_lib}"
                        "${_blacs_tgt}")
                      add_library(${_scalapack_tgt} INTERFACE IMPORTED)
                      set_target_properties(${_scalapack_tgt} PROPERTIES
                        INTERFACE_LINK_LIBRARIES "${_scalapack_libs}")
                      message("Create target for ${_scalapack_tgt}")
                    endif()
                endforeach()
            endforeach()
        endforeach()
    endforeach()
endforeach()

#
# check  MKL_SCALAPACK_FLAVOR / MKL_BLAS_FLAVOR
#
if (NOT DEFINED MKL_SCALAPACK_FLAVOR)
  set(MKL_SCALAPACK_FLAVOR "${MKL_SCALAPACK_FLAVOR_DEFAULT}" CACHE STRING "MKL_SCALAPACK_FLAVOR, default : ${MKL_SCALAPACK_FLAVOR_DEFAULT}")
  message("[mkl_setup] setting MKL_SCALAPACK_FLAVOR to default value : ${MKL_SCALAPACK_FLAVOR_DEFAULT}")
else()
  message("[mkl_setup] MKL_SCALAPACK_FLAVOR already defined, value is ${MKL_SCALAPACK_FLAVOR}")
endif()

# list of valid value for MKL_SCALAPACK_FLAVOR
# used to cross-check that MKL_SCALAPACK_FLAVOR has a valid value
set_property(CACHE MKL_SCALAPACK_FLAVOR PROPERTY
  STRINGS
  scalapack_mpich_intel_32bit_seq_st
  scalapack_ompi_intel_32bit_seq_st
  scalapack_mpich_intel_32bit_omp_st
  scalapack_ompi_intel_32bit_omp_st
  scalapack_mpich_intel_32bit_tbb_st
  scalapack_ompi_intel_32bit_tbb_st
  scalapack_mpich_intel_64bit_seq_st
  scalapack_ompi_intel_64bit_seq_st
  scalapack_mpich_intel_64bit_omp_st
  scalapack_ompi_intel_64bit_omp_st
  scalapack_mpich_intel_64bit_tbb_st
  scalapack_ompi_intel_64bit_tbb_st
  scalapack_mpich_intel_32bit_seq_dyn
  scalapack_ompi_intel_32bit_seq_dyn
  scalapack_mpich_intel_32bit_omp_dyn
  scalapack_ompi_intel_32bit_omp_dyn
  scalapack_mpich_intel_32bit_tbb_dyn
  scalapack_ompi_intel_32bit_tbb_dyn
  scalapack_mpich_intel_64bit_seq_dyn
  scalapack_ompi_intel_64bit_seq_dyn
  scalapack_mpich_intel_64bit_omp_dyn
  scalapack_ompi_intel_64bit_omp_dyn
  scalapack_mpich_intel_64bit_tbb_dyn
  scalapack_ompi_intel_64bit_tbb_dyn)

# cross-check that mkl/scalapack target actually exists
if(TARGET mkl::${MKL_SCALAPACK_FLAVOR})
  set(MKL_SCALAPACK_FLAVOR_FOUND TRUE)
else()
  set(MKL_SCALAPACK_FLAVOR_FOUND FALSE)
  message(FATAL_ERROR "Warning: MKL_SCALAPACK_FLAVOR = ${MKL_SCALAPACK_FLAVOR} not found")
endif()

# set MKL_BLAS_FLAVOR according to MKL_SCALAPACK_FLAVOR
if (NOT DEFINED MKL_BLAS_FLAVOR)
  set(MKL_BLAS_FLAVOR "${MKL_BLAS_FLAVOR_DEFAULT}" CACHE STRING "MKL_BLAS_FLAVOR, default : ${MKL_BLAS_FLAVOR_DEFAULT}")
endif()

# list of valid value for MKL_SCALAPACK_FLAVOR
set_property(CACHE MKL_BLAS_FLAVOR PROPERTY
  STRINGS
  mkl_intel_32bit_seq_st
  mkl_intel_32bit_omp_st
  mkl_intel_32bit_tbb_st
  mkl_gf_32bit_seq_st
  mkl_gf_32bit_omp_st
  mkl_gf_32bit_tbb_st
  mkl_intel_64bit_seq_st
  mkl_intel_64bit_omp_st
  mkl_intel_64bit_tbb_st
  mkl_gf_64bit_seq_st
  mkl_gf_64bit_omp_st
  mkl_gf_64bit_tbb_st
  mkl_intel_32bit_seq_dyn
  mkl_intel_32bit_omp_dyn
  mkl_intel_32bit_tbb_dyn
  mkl_gf_32bit_seq_dyn
  mkl_gf_32bit_omp_dyn
  mkl_gf_32bit_tbb_dyn
  mkl_intel_64bit_seq_dyn
  mkl_intel_64bit_omp_dyn
  mkl_intel_64bit_tbb_dyn
  mkl_gf_64bit_seq_dyn
  mkl_gf_64bit_omp_dyn
  mkl_gf_64bit_tbb_dyn)

# cross-check that mkl/scalapack target actually exists
if(TARGET mkl::${MKL_BLAS_FLAVOR})
  set(MKL_BLAS_FLAVOR_FOUND TRUE)
else()
  set(MKL_BLAS_FLAVOR_FOUND FALSE)
  message("Warning: MKL_BLAS_FLAVOR = ${MKL_BLAS_FLAVOR} not found")
endif()


# get properties for print summary (see below)
get_target_property(MKL_BLAS_HEADERS mkl::${MKL_BLAS_FLAVOR} INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(MKL_BLAS_LIBRARIES mkl::${MKL_BLAS_FLAVOR} INTERFACE_LINK_LIBRARIES)

get_target_property(MKL_SCALAPACK_HEADERS mkl::${MKL_SCALAPACK_FLAVOR} INTERFACE_INCLUDE_DIRECTORIES)
get_target_property(MKL_SCALAPACK_LIBRARIES mkl::${MKL_SCALAPACK_FLAVOR} INTERFACE_LINK_LIBRARIES)

# usually MKL_SCALAPACK_HEADERS is empty
if (${MKL_SCALAPACK_HEADERS} STREQUAL "MKL_SCALAPACK_HEADERS-NOTFOUND")
  set(MKL_SCALAPACK_HEADERS "")
endif()

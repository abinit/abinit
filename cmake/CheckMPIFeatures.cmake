#
# This file is borrowed and slightly modified from
# https://github.com/eschnett/MPIwrapper/blob/main/cmake/CheckMPIFeatures.cmake
#
# function CheckMPIFeatures provides helper to check if MPI implementation
# has the runtime ability to probe GPU-awareness
# - cuda-aware (Nvidia GPU),
# - hip-aware (AMD GPU),
# - ze-aware (INTEL GPU)
#
#
# Apparently Intel MPI (as of version 2021.7.0) doesn't provide header mpi-ext.h, too bad.
#

include(CheckCSourceCompiles)
function(CheckMPIFeatures)
  if (NOT DEFINED HAVE_MPI_EXT OR NOT DEFINED MPI_HAS_QUERY_CUDA_SUPPORT)
    list(JOIN MPI_COMPILE_FLAGS " " CMAKE_REQUIRED_FLAGS)

    #set(CMAKE_REQUIRED_INCLUDES ${MPI_INCLUDE_PATH})
    #set(CMAKE_REQUIRED_LIBRARIES ${MPI_LIBRARIES})
    set(CMAKE_REQUIRED_LIBRARIES MPI::MPI_C)

    # We cannot use check_include_file here as <mpi.h> needs to be
    # included before <mpi-ext.h>, and check_include_file doesn't
    # support this.
    check_c_source_compiles(
      "
        #include <mpi.h>
        #include <mpi-ext.h>
        int main() {
          return 0;
        }
      "
      HAVE_MPI_EXT)

    if(NOT HAVE_MPI_EXT)
      set(HAVE_MPI_EXT 0)
    else()
      set(HAVE_MPI_EXT 1)
    endif()

    list(APPEND CMAKE_REQUIRED_DEFINITIONS -DHAVE_MPI_EXT=${HAVE_MPI_EXT})

    check_c_source_compiles(
      "
        #include <mpi.h>
        #if HAVE_MPI_EXT
        #include <mpi-ext.h>
        #endif
        int main() {
          int result = MPIX_Query_cuda_support();
          return 0;
        }
        "
      MPI_HAS_QUERY_CUDA_SUPPORT)

    if(NOT MPI_HAS_QUERY_CUDA_SUPPORT)
      set(MPI_HAS_QUERY_CUDA_SUPPORT 0)
    else()
      set(MPI_HAS_QUERY_CUDA_SUPPORT 1)
    endif()

    check_c_source_compiles(
      "
        #include <mpi.h>
        #if HAVE_MPI_EXT
        #include <mpi-ext.h>
        #endif
        int main() {
          int result = MPIX_Query_hip_support();
          return 0;
        }
        "
      MPI_HAS_QUERY_HIP_SUPPORT)

    if(NOT MPI_HAS_QUERY_HIP_SUPPORT)
      set(MPI_HAS_QUERY_HIP_SUPPORT 0)
    else()
      set(MPI_HAS_QUERY_HIP_SUPPORT 1)
    endif()

    check_c_source_compiles(
      "
        #include <mpi.h>
        #if HAVE_MPI_EXT
        #include <mpi-ext.h>
        #endif
        int main() {
          int result = MPIX_Query_ze_support();
          return 0;
        }
        "
      MPI_HAS_QUERY_ZE_SUPPORT)

    if(NOT MPI_HAS_QUERY_ZE_SUPPORT)
      set(MPI_HAS_QUERY_ZE_SUPPORT 0)
    else()
      set(MPI_HAS_QUERY_ZE_SUPPORT 1)
    endif()

    list(REMOVE_ITEM CMAKE_REQUIRED_DEFINITIONS -DHAVE_MPI_EXT)
  endif()
endfunction()

CheckMPIFeatures()

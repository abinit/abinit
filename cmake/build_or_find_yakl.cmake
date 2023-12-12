# currently YAKL must be used together with kokkos (im might change in the near future)
# currently abinit can only be used with YAKL_ARCH=Cuda; all other values are not yet supported

# ######################################################################################
# IMPORTANT NOTE: (January 2023):
# ######################################################################################
# We currently don't use official yakl sources, but the fork https://github.com/pkestene/yakl
# with the branch 'fix/install_target'
# A pull request has been made to official repository (https://github.com/mrnorman/YAKL), but
# until it has been merged we need to use the fork.

#
# Does abinit builds yakl ?
#
option(ABINIT_YAKL_BUILD "Turn ON if you want to build yakl (default: OFF)" OFF)

#
# Option to enable / disable yakl detection
#
option(ABINIT_YAKL_WANTED "Turn ON if you want to use yakl (default: OFF)" OFF)

# cross-check that kokkos is also enabled (yakl cannot be use with kokkos in abinit for now)
if(ABINIT_YAKL_WANTED AND (NOT ABINIT_KOKKOS_WANTED))
  message(FATAL_ERROR "You must also enable Kokkos to use YAKL")
endif()

#
# if yakl is wanted
# 1. we cross-check that yakl will use the same target architecture (Cuda, Hip, Sycl, etc...)
#    as kokkos
# 2. either we build yakl or we try to find yakl from environment
#
if(ABINIT_YAKL_WANTED)

  set(HAVE_YAKL 1)

  if (${ABINIT_KOKKOS_BACKEND} STREQUAL "Cuda")
    set(YAKL_ARCH_REQUIRED "CUDA")
    #elseif({ABINIT_KOKKOS_BACKEND} STREQUAL "HIP")
    #  set(YAKL_ARCH_REQUIRED "HIP")
    #elseif({ABINIT_KOKKOS_BACKEND} STREQUAL "SYCL")
    #  set(YAKL_ARCH_REQUIRED "SYCL")
  elseif({ABINIT_KOKKOS_BACKEND} STREQUAL "OpenMP")
    set(YAKL_ARCH_REQUIRED "OPENMP")
  else()
    message(FATAL_ERROR "Unknown or unsupported backend for YAKL: ${YAKL_ARCH_REQUIRED}")
  endif()

  # check if user requested a build of yakl
  # use carefully, it may strongly increase build time
  if(ABINIT_YAKL_BUILD)

    message("[abinit / yakl] Building yakl from source")

    add_compile_definitions(YAKL_ARCH=${YAKL_ARCH_REQUIRED})

    if ((YAKL_ARCH_REQUIRED STREQUAL "CUDA") OR
        (YAKL_ARCH_REQUIRED STREQUAL "HIP")  OR
        (YAKL_ARCH_REQUIRED STREQUAL "SYCL") )
      set(YAKL_MANAGED_MEMORY True)
    endif()

    find_package(MPI COMPONENTS Fortran C CXX)
    if (MPI_FOUND)
      set(YAKL_HAVE_MPI True)
    endif()

    set(YAKL_ARCH ${YAKL_ARCH_REQUIRED})
    message("[abinit / YAKL] Building yakl from source for YAKL_ARCH=${YAKL_ARCH}")

    set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/external)

    #
    # Yakl doesn't provide releases yet, so just use latest commit (end of Sept. 2022)
    #
    include (FetchContent)

    FetchContent_Declare(
      yakl_external
      #GIT_REPOSITORY https://github.com/pkestene/YAKL.git
      #GIT_TAG 67f5b71a3654c8db6600e243af0a8a9d5592444f
      GIT_REPOSITORY https://github.com/mrnorman/YAKL.git
      GIT_TAG e5ae8cbc1eadb5badab77e2a194a9d38f39f5df6
      )

    # Import yakl targets (download, and call add_subdirectory)
    FetchContent_MakeAvailable(yakl_external)

    # library alias
    add_library(abinit::yakl ALIAS yakl_fortran_interface)

  else(ABINIT_YAKL_BUILD)

    #
    # yakl is wanted, but we don't build it, so it must be found it in environment
    #

    # detect yakl from environment.
    # yakl-config.cmake must define YAKL_ARCH variable
    find_package(yakl REQUIRED)

    if(TARGET yakl::yakl)

      if (NOT DEFINED YAKL_ARCH)
        message(FATAL_ERROR "YAKL_ARCH is not defined. Please check you installed yakl correctly.")
      endif()

      if (${YAKL_ARCH_REQUIRED} STREQUAL "${YAKL_ARCH}")
        message(STATUS "Kokkos and YAKL target architecture match !")
      else()
        message(FATAL_ERROR "Kokkos and YAKL target architecture don't match !")
      endif()

      message("[abinit / yakl] Yakl found via find_package with YAKL_ARCH=${YAKL_ARCH}")
      set(ABINIT_YAKL_FOUND True)
      set(HAVE_YAKL 1)

      # library alias
      add_library(abinit::yakl ALIAS yakl::yakl_fortran_interface)

    else()

      message(FATAL_ERROR "[abinit / yakl] yakl is required but not found by find_package. Please adjust your env variable CMAKE_PREFIX_PATH (or yakl_ROOT) to where yakl is installed on your machine !")

    endif()

  endif(ABINIT_YAKL_BUILD)

endif(ABINIT_YAKL_WANTED)

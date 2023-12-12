# Two alternatives:
# 1. If ABINIT_KOKKOS_BUILD is ON, we download kokkos sources and build them using FetchContent (which actually uses add_subdirectory)
# 2. If ABINIT_KOKKOS_BUILD is OFF (default), we don't build kokkos, but use find_package for setup (you must have kokkos already installed)

# NOTE about required C++ standard
# we better chose to set the minimum C++ standard level if not already done:
# - when building kokkos <  4.0.00, it defaults to c++-14
# - when building kokkos >= 4.0.00, it defaults to c++-17
# - when using installed kokkos, we set C++ standard according to kokkos version

#
# Do we want to build kokkos (https://github.com/kokkos/kokkos) ?
#
option(ABINIT_KOKKOS_BUILD "Turn ON if you want to build kokkos (default: OFF)" OFF)

#
# Option to enable / disable kokkos detection
#
option(ABINIT_KOKKOS_WANTED "Turn ON if you want to use kokkos (default: OFF)" OFF)

#
# Option to use git (instead of tarball release) for downloading kokkos and kokkos-fortran-interop
#
option(ABINIT_KOKKOS_USE_GIT "Turn ON if you want to use git to download Kokkos sources (default: OFF)" OFF)

#
# Options to specify target device backend
#

# set default backend
set(ABINIT_KOKKOS_BACKEND "Undefined" CACHE STRING
  "Kokkos backend device")

# Set the possible values for kokkos backend device
set_property(CACHE ABINIT_KOKKOS_BACKEND PROPERTY STRINGS
  "OpenMP" "Cuda" "HIP" "Undefined")


if(ABINIT_KOKKOS_WANTED)

  # check if user requested a build of kokkos
  # use carefully, it may strongly increase build time
  if(ABINIT_KOKKOS_BUILD)

    message("[abinit / kokkos] Building kokkos from source")

    set_property(DIRECTORY PROPERTY EP_BASE ${CMAKE_BINARY_DIR}/external)

    # Kokkos default build options

    # set install path
    list (APPEND ABINIT_KOKKOS_CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${KOKKOS_INSTALL_DIR})

    # use predefined cmake args
    # can be override on the command line
    if (ABINIT_KOKKOS_BACKEND MATCHES "Cuda")

      if (NOT ABINIT_ENABLE_GPU_CUDA)
        message(FATAL_ERROR "[abinit / kokkos] You can't use Kokkos::Cuda backend if ABINIT_ENABLE_GPU_CUDA is OFF")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_HWLOC) OR (NOT Kokkos_ENABLE_HWLOC))
        set(Kokkos_ENABLE_HWLOC ON CACHE BOOL "")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_OPENMP) OR (NOT Kokkos_ENABLE_OPENMP))
        set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_CUDA) OR (NOT Kokkos_ENABLE_CUDA))
        set(Kokkos_ENABLE_CUDA ON CACHE BOOL "")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_CUDA_LAMBDA) OR (NOT Kokkos_ENABLE_CUDA_LAMBDA))
        set(Kokkos_ENABLE_CUDA_LAMBDA ON CACHE BOOL "")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_CUDA_CONSTEXPR) OR (NOT Kokkos_ENABLE_CUDA_CONSTEXPR))
        set(Kokkos_ENABLE_CUDA_CONSTEXPR ON CACHE BOOL "")
      endif()

      # Note : cuda architecture will probed by kokkos cmake configure

    elseif(ABINIT_KOKKOS_BACKEND MATCHES "HIP")

      message(FATAL_ERROR "[abinit / kokkos] HIP backend to supported yet")

      # uncomment the following lines when HIP will be supported

      # if (NOT ABINIT_ENABLE_GPU_HIP)
      #   message(FATAL_ERROR "[abinit / kokkos] You can't use Kokkos::HIP backend if ABINIT_ENABLE_GPU_HIP is OFF")
      # endif()

      # if ((NOT DEFINED Kokkos_ENABLE_HWLOC) OR (NOT Kokkos_ENABLE_HWLOC))
      #   set(Kokkos_ENABLE_HWLOC ON CACHE BOOL "")
      # endif()

      # if ((NOT DEFINED Kokkos_ENABLE_OPENMP) OR (NOT Kokkos_ENABLE_OPENMP))
      #   set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
      # endif()

      # if ((NOT DEFINED Kokkos_ENABLE_HIP) OR (NOT Kokkos_ENABLE_HIP))
      #   set(Kokkos_ENABLE_HIP ON CACHE BOOL "")
      # endif()

    elseif(ABINIT_KOKKOS_BACKEND MATCHES "OpenMP")

      if (ABINIT_ENABLE_GPU_CUDA)
        message(WARNING "[abinit / kokkos] ABINIT_ENABLE_GPU_CUDA is ON, you should consider setting ABINIT_KOKKOS_BACKEND=Cuda")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_HWLOC) OR (NOT Kokkos_ENABLE_HWLOC))
        set(Kokkos_ENABLE_HWLOC ON CACHE BOOL "")
      endif()

      if ((NOT DEFINED Kokkos_ENABLE_OPENMP) OR (NOT Kokkos_ENABLE_OPENMP))
        set(Kokkos_ENABLE_OPENMP ON CACHE BOOL "")
      endif()

    elseif(ABINIT_KOKKOS_BACKEND MATCHES "Undefined")

      message(FATAL_ERROR "[abinit / kokkos] You must chose a valid ABINIT_KOKKOS_BACKED !")

    endif()

    # we set c++ standard to c++-14 as kokkos < 4.0.00 requires at least c++-14
    # if later, we chose to upgrade to kokkos >= 4.0.00, we need to set it to c++-17
    if (NOT "${CMAKE_CXX_STANDARD}")
      set(CMAKE_CXX_STANDARD 14)
    endif()

    #find_package(Git REQUIRED)
    include (FetchContent)

    if (ABINIT_KOKKOS_USE_GIT)
      FetchContent_Declare( kokkos_external
        GIT_REPOSITORY https://github.com/kokkos/kokkos.git
        GIT_TAG 3.7.01
        )
    else()
      FetchContent_Declare( kokkos_external
        URL https://github.com/kokkos/kokkos/archive/refs/tags/3.7.01.tar.gz
        )
    endif()

    # Import kokkos targets (download, and call add_subdirectory)
    FetchContent_MakeAvailable(kokkos_external)

    if(TARGET Kokkos::kokkos)
      message("[abinit / kokkos] Kokkos found (using FetchContent)")
      set(ABINIT_KOKKOS_FOUND True)
      set(HAVE_KOKKOS 1)
    else()
      message("[abinit / kokkos] we shouldn't be here. We've just integrated kokkos build into abinit build !")
    endif()

    set(ABINIT_KOKKOS_BUILTIN TRUE)

  else()

    #
    # check if already installed kokkos exists
    #
    find_package(Kokkos 3.6.01 REQUIRED)

    if(TARGET Kokkos::kokkos)

      # set default c++ standard according to Kokkos version
      # Kokkos >= 4.0.00 requires c++-17
      if (NOT "${CMAKE_CXX_STANDARD}")
        if ( ${Kokkos_VERSION} VERSION_LESS 4.0.00)
          set(CMAKE_CXX_STANDARD 14)
        else()
          set(CMAKE_CXX_STANDARD 17)
        endif()
      endif()

      # kokkos_check is defined in KokkosConfigCommon.cmake
      kokkos_check( DEVICES "OpenMP" RETURN_VALUE KOKKOS_DEVICE_ENABLE_OPENMP)
      kokkos_check( DEVICES "Cuda" RETURN_VALUE KOKKOS_DEVICE_ENABLE_CUDA)
      kokkos_check( DEVICES "HIP" RETURN_VALUE KOKKOS_DEVICE_ENABLE_HIP)

      kokkos_check( TPLS "HWLOC" RETURN_VALUE Kokkos_TPLS_HWLOC_ENABLED)

      if(ABINIT_ENABLE_GPU_CUDA)
        set(ABINIT_KOKKOS_BACKEND "Cuda")
        kokkos_check( OPTIONS CUDA_LAMBDA RETURN_VALUE Kokkos_CUDA_LAMBDA_ENABLED)
        kokkos_check( OPTIONS CUDA_CONSTEXPR RETURN_VALUE Kokkos_CUDA_CONSTEXPR_ENABLED)
        kokkos_check( OPTIONS CUDA_UVM RETURN_VALUE Kokkos_CUDA_UVM_ENABLED)
      elseif(ABINIT_ENABLE_GPU_HIP)
        set(ABINIT_KOKKOS_BACKEND "HIP")
      elseif(KOKKOS_DEVICE_ENABLE_OPENMP)
        set(ABINIT_KOKKOS_BACKEND "OpenMP")
      endif()

      message("[abinit / kokkos] Kokkos found via find_package; default backend is ${ABINIT_KOKKOS_BACKEND}")
      set(ABINIT_KOKKOS_FOUND True)
      set(HAVE_KOKKOS 1)

    else()

      message(FATAL_ERROR "[abinit / kokkos] Kokkos is required but not found by find_package. Please adjust your env variable CMAKE_PREFIX_PATH (or Kokkos_ROOT) to where Kokkos is installed on your machine !")

    endif()

  endif()

else(ABINIT_KOKKOS_WANTED)

  message(NOTICE "[abinit / kokkos] kokkos is not wanted")

endif(ABINIT_KOKKOS_WANTED)

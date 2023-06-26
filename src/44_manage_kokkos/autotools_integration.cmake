#
# Allow find_package() to use <PackageName>_ROOT variables,
# either cmake variable or environment variables
# see https://cmake.org/cmake/help/latest/policy/CMP0074.html
#
if(NOT CMAKE_VERSION VERSION_LESS 3.12)
  cmake_policy(SET CMP0074 NEW)
endif()

project(manage_kokkos_sublib
  LANGUAGES C CXX Fortran
  DESCRIPTION "manage_kokkos_sublib is a CMake project aiming to integrate Kokkos code within abinit autotools buildsystem.")

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

  #FIXME Temporary, how should we detect HIP over CUDA ? Grepping config.h maybe ?
  set(ABINIT_ENABLE_GPU_CUDA YES)
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

endif()

message("[44_manage_kokkos] Kokkos found via find_package; default backend is ${ABINIT_KOKKOS_BACKEND}")
message("[44_manage_kokkos] C flags == ${CMAKE_C_FLAGS}")

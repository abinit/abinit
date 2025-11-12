#
# Does abinit use scalapack ?
#
option(ABINIT_SCALAPACK_ENABLED "Turn ON if you want abinit to use scalapack (default: OFF)" OFF)

#
# Does abinit use ELPA (https://elpa.mpcdf.mpg.de/) ?
#
option(ABINIT_ELPA_ENABLED "Turn ON if you want abinit to use ELPA (default: OFF)" OFF)

#
# Other options
#
option(ABINIT_DO_DEBUG_CONTRACT "Do you want to activate design-by-contract debugging tests (default OFF)." OFF)
if (ABINIT_DO_DEBUG_CONTRACT)
  set(DEBUG_CONTRACT 1)
endif()

option(ABINIT_DO_DEBUG_VERBOSE "Turn on verbose debug messages in the source code (default OFF)." OFF)
if (ABINIT_DO_DEBUG_VERBOSE)
  set(DEBUG_VERBOSE 1)
endif()

option(ABINIT_DO_DEBUG_VERBOSE_GPU "Turn on verbose debug messages in the GPU source code (default OFF)." OFF)
if (ABINIT_DO_DEBUG_VERBOSE_GPU)
  set(DEBUG_VERBOSE_GPU 1)
endif()

option(ABINIT_DO_MEM_PROFILING "Turn on memory profiling (default OFF)." OFF)
if (ABINIT_DO_MEM_PROFILING)
  set(HAVE_MEM_PROFILING 1)
endif()

option(ABINIT_AVX_SAFE_MODE "Disable vectorization in problematic procedures (default: OFF" OFF)
if (ABINIT_AVX_SAVE_MODE)
  set(HAVE_AVX_SAFE_MODE 1)
endif()

option(ABINIT_ENABLE_CCLOCK "Use C clock for timings (default: OFF)" OFF)
if(ABINIT_ENABLE_CCLOCK)
  set(HAVE_CCLOCK 1)
endif()

option(ABINIT_ENABLE_CRPA_OPTIM "Enable optimize cRPA calculations for ifort <= 17.0 (default: no)" OFF)
if(ABINIT_ENABLE_CRPA_OPTIM)
  set(HAVE_CRPA_OPTIM 1)
endif()

option(ABINIT_ENABLE_MPI_IO_DEFAULT "Enable to use MPI I/O as default I/O library (default: no)" OFF)
if(ABINIT_ENABLE_MPI_IO_DEFAULT)
  set(HAVE_MPI_IO_DEFAULT 1)
endif()

option(ABINIT_ENABLE_PARALLEL_HDF5 "Enable the use of parallel HDF5 (default: no)" OFF)

option(ABINIT_ENABLE_NETCDF_DEFAULT "Use this option if you want to use NetCDF as default I/O library (default: no)" OFF)
if(ABINIT_ENABLE_NETCDF_DEFAULT)
  set(HAVE_NETCDF_DEFAULT 1)
endif()

option(ABINIT_ENABLE_GW_DPC "Enable GW computing with double precision (default OFF)" OFF)
if(ABINIT_ENABLE_GW_DPC)
  set(HAVE_GW_DPC 1)
endif()

option(ABINIT_ENABLE_LIBTETRA "Enable internal support for libtetra (default ON)" ON)
if(ABINIT_ENABLE_LIBTETRA)
  set(HAVE_LIBTETRA_ABINIT 1)
endif()

option(ABINIT_ENABLE_TRIQS "Enable support for TRIQS (default OFF)" OFF)

option(ABINIT_ENABLE_PYTHON_INVOCATION "Enable python invocation (default OFF)" OFF)
if(ABINIT_ENABLE_PYTHON_INVOCATION)
  set(HAVE_PYTHON_INVOCATION 1)
  set(DO_BUILD_67_PYTHON_INVOCATION_EXT ON)
endif()

option (ABINIT_ENFORCE_GPU_AWARE_MPI "Some MPI GPU-aware implementations are not well detected; use this variable to enforce if you know that your MPI implementation is GPU-aware." OFF)

option(ABINIT_ENABLE_GPU_CUDA "Enable GPU build (using Nvidia CUDA backend, default OFF)" OFF)
if(ABINIT_ENABLE_GPU_CUDA)

  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    message(FATAL_ERROR
      "When using CUDA, you need to provide -DCMAKE_CUDA_ARCHITECTURES "
      "with the NVIDIA GPU compute capability matching your environment.\n"
      "For example, if targetting NVIDIA GPU A100:\n"
      "-DCMAKE_CUDA_ARCHITECTURES=80 \n"
    )
  endif()

  include(CheckLanguage)
  check_language(CUDA)
  if (CMAKE_CUDA_COMPILER)
    enable_language(CUDA)

    # make sure cuda libs are available (cufft, etc..)

    # FindCUDAToolkit.cmake is available since cmake 3.17 and is able
    # to correctly find cuda library either when using nvcc from cuda toolkit
    # or using nvc++ from nvhpc (but only in cmake 3.22.0)
    # if you want to use nvc++ you need at least cmake 3.22.0
    #
    # the following will make target like CUDA::cublas and CUDA::cufft available
    message(STATUS "Using CUDAToolkit macros")
    find_package(CUDAToolkit REQUIRED)

  endif()

  set(HAVE_GPU_CUDA 1)

  # enfore the use of double precision, as in many locations double are used
  # inconditionally
  set(HAVE_GPU_CUDA_DP 1)

  set(HAVE_GPU 1)
  set(HAVE_GPU_SERIAL 1)

endif()

option(ABINIT_ENABLE_GPU_HIP "Enable GPU build (using AMD HIP backend, default OFF)" OFF)
if(ABINIT_ENABLE_GPU_HIP)

  if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
    message(FATAL_ERROR
      "When using HIP, you need to provide -DCMAKE_HIP_ARCHITECTURES "
      "with the AMD GPU target matching your environment.\n"
      "For example, if targetting AMD GPU Instinct MI250:\n"
      "-DCMAKE_HIP_ARCHITECTURES=gfx90a \n"
    )
  endif()

  # For shutting annoying warning from ROCM Cmake module
  set(AMDGPU_TARGETS ${CMAKE_HIP_ARCHITECTURES})

  find_package(HIP)
  find_package(hipfft)
  find_package(rocfft)
  find_package(rocblas)
  find_package(hipblas)
  find_package(hipsolver)

  set(HAVE_GPU_HIP 1)

  set(HAVE_GPU 1)
  set(HAVE_GPU_SERIAL 1)

  add_compile_definitions("__HIP_PLATFORM_AMD__")

endif()

option(ABINIT_ENABLE_GPU_MARKERS "Enable GPU markers for profiling (requires NVTX3 or ROCM/HIP, default OFF)" OFF)
option(ABINIT_ENABLE_NVTX  "Enable NVTX markers for profiling (requires CUDA install, default OFF)" OFF)
option(ABINIT_ENABLE_ROCTX "Enable ROCTX markers for profiling (requires ROCM install, default OFF)" OFF)
if(ABINIT_ENABLE_GPU_MARKERS AND (ABINIT_ENABLE_NVTX OR ABINIT_ENABLE_ROCTX) OR (ABINIT_ENABLE_NVTX AND ABINIT_ENABLE_ROCTX))
  message(FATAL_ERROR "Please, activate only one of ABINIT_ENABLE_GPU_MARKERS, ABINIT_ENABLE_NVTX or ABINIT_ENABLE_ROCTX.")
endif()
if(ABINIT_ENABLE_GPU_MARKERS AND NOT (ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP))
  message(FATAL_ERROR "Please, activate ABINIT_ENABLE_GPU_MARKERS only along ABINIT_ENABLE_CUDA or ABINIT_ENABLE_HIP.")
endif()
if(ABINIT_ENABLE_NVTX OR (ABINIT_ENABLE_GPU_MARKERS AND ABINIT_ENABLE_GPU_CUDA))

  if(NOT ABINIT_ENABLE_GPU_CUDA)
    find_package(CUDAToolkit REQUIRED)
  endif()
  # check if NVTX library is available from imported CUDA install
  if (TARGET CUDA::nvToolsExt)
    set(HAVE_GPU_CUDA10 1)
    set(HAVE_GPU_MARKERS 1)
    set(HAVE_GPU_MARKERS_NVTX 1)
    add_library(abinit::gpu_markers ALIAS CUDA::nvToolsExt)
  # Else check manually for nvtx3interop (new name from CUDA 12.9)
  else()
    find_library(NVTX_LIBRARY
      NAMES libnvtx3interop.so
      PATHS ${CUDAToolkit_ROOT}/lib ${CUDAToolkit_ROOT}/lib64
      DOC "Location of the NVTX library"
    )
    if (NVTX_LIBRARY)
      set(HAVE_GPU_CUDA10 1)
      set(HAVE_GPU_MARKERS 1)
      set(HAVE_GPU_MARKERS_NVTX 1)
      add_library(nvtx3interop UNKNOWN IMPORTED ${NVTX_LIBRARY})
      set_target_properties(nvtx3interop PROPERTIES IMPORTED_LOCATION ${NVTX_LIBRARY})
      add_library(abinit::gpu_markers ALIAS nvtx3interop)
    else()
      message(SEND_ERROR "NVTX (GPU markers for NVIDIA tools) required but not found")
    endif()
  endif()

elseif(ABINIT_ENABLE_ROCTX OR (ABINIT_ENABLE_GPU_MARKERS AND ABINIT_ENABLE_GPU_HIP))
  # check ROCtx library is available
  find_library(ROCTX_LIBRARY
    NAMES libroctx64.so
    PATHS ${ROCM_ROOT}/roctracer/lib ${ROCM_PATH}/roctracer/lib ${ROCM_HOME}/roctracer/lib
    DOC "Location of the ROCTX library"
  )
  if (ROCTX_LIBRARY)
    message(STATUS "ROCTX found: ${ROCTX_LIBRARY}")
    add_library(roctx64 UNKNOWN IMPORTED ${ROCTX_LIBRARY})
    set_target_properties(roctx64 PROPERTIES IMPORTED_LOCATION ${ROCTX_LIBRARY})
    add_library(abinit::gpu_markers ALIAS roctx64)
    set(HAVE_GPU_MARKERS 1)
    set(HAVE_GPU_MARKERS_ROCTX 1)
  else()
    message(SEND_ERROR "ROCTX (GPU markers for AMD tools) required but not found")
  endif()

endif()

if (ABINIT_ENABLE_GPU_CUDA OR ABINIT_ENABLE_GPU_HIP OR ABINIT_ENABLE_GPU_MARKERS)
  set(DO_BUILD_17_GPU_TOOLBOX TRUE)
else()
  set(DO_BUILD_17_GPU_TOOLBOX FALSE)
endif()

option(ABINIT_ENABLE_NVIDIA_UNIFIED_MEM
  "Enable Unified Memory for NVIDIA GPU (requires OpenMP Offload & NVHPC compiler, default OFF)" OFF)
# Checks for this option occurs after compiler and OpenMP offload settings in root CMakeLists.txt

if (ABINIT_ENABLE_GPU_CUDA)
  set(DO_BUILD_46_MANAGE_CUDA TRUE)
else()
  set(DO_BUILD_46_MANAGE_CUDA FALSE)
endif()

if (ABINIT_ENABLE_GPU_CUDA AND ABINIT_KOKKOS_WANTED)
  set(DO_BUILD_16_KOKKOS_TOOLBOX TRUE)
  set(DO_BUILD_44_MANAGE_KOKKOS TRUE)
else()
  set(DO_BUILD_16_KOKKOS_TOOLBOX FALSE)
  set(DO_BUILD_44_MANAGE_KOKKOS FALSE)
endif()


option(ABINIT_ENABLE_LIBPAW_INTERNAL "Enable building libpaw as part of abinit build (default ON)" ON)
if (ABINIT_ENABLE_LIBPAW_INTERNAL)
  set(HAVE_LIBPAW_ABINIT 1)
endif()


set(ABI_DEBUG_FLAVOR "basic" CACHE STRING
  "Abinit C compiler debug flavor : basic, verbose, enhanced, paranoid, naughty")
set_property(CACHE ABI_DEBUG_FLAVOR PROPERTY STRINGS basic verbose enhanced paranoid naughty)

set(ABI_OPTIM_FLAVOR "standard" CACHE STRING
  "Abinit C/Fortran compiler optim flavor : safe, standard, aggressive")
set_property(CACHE ABI_OPTIM_FLAVOR PROPERTY STRINGS safe standard aggressive)

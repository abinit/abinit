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

option(ABINIT_ENABLE_LONG_LINES "Enable to have long lines in fortran source code (default: no)" OFF)
if(ABINIT_ENABLE_LONG_LINES)
  set(HAVE_FC_LONG_LINES 1)
endif()

option(ABINIT_ENABLE_MPI_IO "Enable to have MPI I/O support (default: no)" OFF)
if(ABINIT_ENABLE_MPI_IO)
  set(HAVE_MPI_IO 1)
endif()

option(ABINIT_ENABLE_MPI_IO_DEFAULT "Enable to use MPI I/O as default I/O library (default: no)" OFF)
if(ABINIT_ENABLE_MPI_IO_DEFAULT)
  set(HAVE_MPI_IO_DEFAULT 1)
endif()

option(ABINIT_ENABLE_MPI_INPLACE "Enable the use of MPI_IN_PLACE (default: no)" OFF)
if(ABINIT_ENABLE_MPI_INPLACE)
  set(HAVE_MPI2_INPLACE 1)
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

option(ABINIT_ENABLE_PYTHON_INVOCATION "Enable python invocation (default OFF)" OFF)
if(ABINIT_ENABLE_PYTHON_INVOCATION)
  set(HAVE_PYTHON_INVOCATION 1)
  set(DO_BUILD_67_PYTHON_INVOCATION_EXT ON)
endif()


option (ABINIT_ENFORCE_CUDA_AWARE_MPI "Some MPI cuda-aware implementation are not well detected; use this variable to enforce if you that your MPI implementation is Cuda-aware." OFF)

option(ABINIT_ENABLE_GPU_CUDA "Enable GPU build (using Nvidia CUDA backend, default OFF)" OFF)
if(ABINIT_ENABLE_GPU_CUDA)
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
    message("Using CUDAToolkit macros")
    find_package(CUDAToolkit REQUIRED)

  endif()

  set(HAVE_GPU_CUDA 1)

  # enfore the use of double precision, as in many locations double are used
  # inconditionally
  set(HAVE_GPU_CUDA_DP 1)

  set(HAVE_GPU 1)
  set(HAVE_GPU_SERIAL 1)

  # check nvtx library is available
  if (TARGET CUDA::nvToolsExt)
    set(HAVE_GPU_CUDA10 1)
    set(HAVE_GPU_NVTX_V3 1)
  endif()

endif()

if (ABINIT_ENABLE_GPU_CUDA)
  set(DO_BUILD_17_GPU_TOOLBOX TRUE)
  set(DO_BUILD_46_MANAGE_CUDA TRUE)
else()
  set(DO_BUILD_17_GPU_TOOLBOX FALSE)
  set(DO_BUILD_46_MANAGE_CUDA FALSE)
endif()

if (ABINIT_ENABLE_GPU_CUDA AND ABINIT_KOKKOS_WANTED)
  set(DO_BUILD_44_MANAGE_KOKKOS TRUE)
else()
  set(DO_BUILD_44_MANAGE_KOKKOS FALSE)
endif()


option(ABINIT_ENABLE_LIBPAW_INTERNAL "Enable building libpaw as part of abinit build (default ON)" ON)
if (ABINIT_ENABLE_LIBPAW_INTERNAL)
  set(HAVE_LIBPAW_ABINIT 1)
endif()


set(ABI_DEBUG_FLAVOR "basic" CACHE STRING
  "Abinit C compiler debug flavor : basic, verbose, enhanced, paranoid, naughty")
set_property(CACHE ABI_DEBUG_FLAVOR PROPERTY STRINGS basic verbose enhanced paranoid naughty)

set(ABI_OPTIM_FLAVOR "safe" CACHE STRING
  "Abinit C/Fortran compiler optim flavor : safe, standard, aggressive")
set_property(CACHE ABI_OPTIM_FLAVOR PROPERTY STRINGS safe standard aggressive)

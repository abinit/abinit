/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/*
 * Copyright (C) 2005-2022 ABINIT Group (Yann Pouillon)
 *
 * This file is part of the Abinit software package. For license information,
 * please see the COPYING file in the top-level directory of the Abinit source
 * distribution.
 *
 */

/* Abinit configuration */

#ifndef _ABINIT_CONFIG_H
#define _ABINIT_CONFIG_H

#ifdef __INTEL_COMPILER
#define FC_INTEL 1
#endif



/* Abinit target description. */
#define ABINIT_TARGET "x86_64_linux_gnu11.3"

/* Abinit whole version number. */
#define ABINIT_VERSION "9.9.0"

/* Abinit base version number. */
#define ABINIT_VERSION_BASE "9.9"

/* Abinit build date. */
#define ABINIT_VERSION_BUILD "20221117"

/* Abinit major version number. */
#define ABINIT_VERSION_MAJOR "9"

/* Abinit micro version number (patch level). */
#define ABINIT_VERSION_MICRO "0"

/* Abinit minor version number. */
#define ABINIT_VERSION_MINOR "9"

/* Define if building universal (internal helper macro) */
/* #undef AC_APPLE_UNIVERSAL_BUILD */

/* Define to 1 if you are using the ARM C compiler. */
/* #undef CC_ARM */

/* Define to 1 if you are using the GNU C compiler. */
#define CC_GNU 1

/* Define to 1 if you are using the IBM XL C compiler. */
/* #undef CC_IBM */

/* Define to 1 if you are using the Intel C compiler. */
/* #undef CC_INTEL */

/* Define to 1 if you are using the LLVM Clang C compiler. */
/* #undef CC_LLVM */

/* Define to 1 if you are using the Portland Group C compiler. */
/* #undef CC_PGI */

/* Define to 1 if you are using the ARM C++ compiler. */
/* #undef CXX_ARM */

/* Define to 1 if you are using the GNU C++ compiler. */
#define CXX_GNU 1

/* Define to 1 if you are using the IBM XL C++ compiler. */
/* #undef CXX_IBM */

/* Define to 1 if you are using the Intel C++ compiler. */
/* #undef CXX_INTEL */

/* Define to 1 if you are using the LLVM Clang++ C++ compiler. */
/* #undef CXX_LLVM */

/* Define to 1 if you are using the Portland Group C++ compiler. */
/* #undef CXX_PGI */

/* Define to 1 if you want to activate design-by-contract debugging tests. */
/* #undef DEBUG_CONTRACT */

/* Define to 1 to build debugging instructions in the source code. */
/* #undef DEBUG_MODE */

/* Define to 1 to turn on verbose debug messages in the source code. */
/* #undef DEBUG_VERBOSE */

/* Define to 1 if you are using the ABSOFT Fortran compiler. */
/* #undef FC_ABSOFT */

/* Define to 1 if you are using the LLVM Flang Fortran AOCC compiler. */
/* #undef FC_AOCC */

/* Define to 1 if you are using the ARM Fortran compiler. */
/* #undef FC_ARM */

/* Define to 1 if you are using the Cray Fortran compiler. */
/* #undef FC_CRAY */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if you are using the GNU Fortran compiler. */
#define FC_GNU 1

/* Define to 1 if you are using the IBM XL Fortran compiler. */
/* #undef FC_IBM */

/* Define to 1 if you are using the Intel Fortran compiler. */
/* #undef FC_INTEL */

/* Define to 1 if you are using the LLVM Flang Fortran compiler. */
/* #undef FC_LLVM */

/* Define to 1 if you are using the NAGWare Fortran 95 compiler. */
/* #undef FC_NAG */

/* Define to 1 if you are using the Portland Group Fortran compiler. */
/* #undef FC_PGI */

/* Define to 1 if you have the ABINIT Common library. */
/* #undef HAVE_ABINIT_COMMON */

/* Define to 1 if you have the `abort' function. */
#define HAVE_ABORT 1

/* Define to 1 if you have the AtomPAW library. */
/* #undef HAVE_ATOMPAW */

/* Define to 1 if you want to disable vectorization in problematic procedures.
   */
/* #undef HAVE_AVX_SAFE_MODE */

/* Define to 1 if you have the BigDFT library. */
/* #undef HAVE_BIGDFT */

/* Define to 1 if you want to activate Bethe-Salpeter unpacking. */
/* #undef HAVE_BSE_UNPACKED */

/* Use C clock for timings. */
/* #undef HAVE_CCLOCK */

/* Define to 1 if you have the `clock_gettime' function. */
/* #undef HAVE_CLOCK_GETTIME */

/* Define to 1 if you want to enable optimize cRPA calculations for ifort <=
   17.0. */
/* #undef HAVE_CRPA_OPTIM */

/* Define to 1 if you have the <cublas.h> header file. */
/* #undef HAVE_CUBLAS_H */

/* Define to 1 if you have the <cuda_runtime_api.h> header file. */
/* #undef HAVE_CUDA_RUNTIME_API_H */

/* Define to 1 if you have the <cufft.h> header file. */
/* #undef HAVE_CUFFT_H */

/* Define to 1 if you have the DFTI library. */
/* #undef HAVE_DFTI */

/* Define to 1 if you have the <errno.h> header file. */
#define HAVE_ERRNO_H 1

/* Define to 1 if you have the <f90papi.h> header file. */
/* #undef HAVE_F90PAPI_H */

/* Define to 1 if your Fortran compiler supports allocatable arrays in
   datatypes. */
#define HAVE_FC_ALLOCATABLE_DTARRAYS 1

/* Define to 1 if your Fortran compiler supports the asynchronous attribute.
   */
#define HAVE_FC_ASYNC 1

/* Define to 1 if your Fortran compiler supports BACKTRACE. */
#define HAVE_FC_BACKTRACE 1

/* Define to 1 if your Fortran compiler supports GET_COMMAND_ARGUMENT. */
#define HAVE_FC_COMMAND_ARGUMENT 1

/* Define to 1 if your Fortran compiler supports EXECUTE_COMMAND_LINE. */
#define HAVE_FC_COMMAND_LINE 1

/* Define to 1 if your Fortran compiler supports the contiguous attribute. */
#define HAVE_FC_CONTIGUOUS 1

/* Define to 1 if your Fortran compiler supports cpu_time(). */
#define HAVE_FC_CPUTIME 1

/* Define to 1 if your Fortran compiler supports etime(). */
/* #undef HAVE_FC_ETIME */

/* Define to 1 if your Fortran compiler supports exit(). */
#define HAVE_FC_EXIT 1

/* Define to 1 if your Fortran compiler supports flush(). */
#define HAVE_FC_FLUSH 1

/* Define to 1 if your Fortran compiler supports flush_(). */
/* #undef HAVE_FC_FLUSH_ */

/* Define to 1 if your Fortran compiler supports gamma(). */
#define HAVE_FC_GAMMA 1

/* Define to 1 if your Fortran compiler supports getenv(). */
#define HAVE_FC_GETENV 1

/* Define to 1 if your Fortran compiler supports getpid(). */
/* #undef HAVE_FC_GETPID */

/* Define to 1 if your Fortran compiler supports IEEE_ARITHMETIC module. */
#define HAVE_FC_IEEE_ARITHMETIC 1

/* Define to 1 if your Fortran compiler supports IEEE_EXCEPTIONS. */
#define HAVE_FC_IEEE_EXCEPTIONS 1

/* Define to 1 if your Fortran compiler accepts quadruple integers. */
#define HAVE_FC_INT_QUAD 1

/* Define to 1 if your Fortran compiler supports IOMSG. */
#define HAVE_FC_IOMSG 1

/* Define to 1 if your Fortran compiler provides the iso_c_binding module. */
#define HAVE_FC_ISO_C_BINDING 1

/* Define to 1 if your Fortran compiler supports 2008 standard in
   ISO_FORTRAN_ENV. */
#define HAVE_FC_ISO_FORTRAN_2008 1

/* Define to 1 if your Fortran compiler supports long lines. */
#define HAVE_FC_LONG_LINES 1

/* Define to 1 if your Fortran compiler supports \newline in a macros. */
/* #undef HAVE_FC_MACRO_NEWLINE */

/* Define to 1 if your Fortran compiler supports MOVE_ALLOC (F2003). */
#define HAVE_FC_MOVE_ALLOC 1

/* Define to 1 if your Fortran compiler can shape arrays on-the-fly. */
#define HAVE_FC_ON_THE_FLY_SHAPE 1

/* Define to 1 if your Fortran compiler supports the private attribute. */
#define HAVE_FC_PRIVATE 1

/* Define to 1 if your Fortran compiler supports the protected attribute. */
#define HAVE_FC_PROTECTED 1

/* Define to 1 if your Fortran compiler supports shiftl() and shiftr(). */
#define HAVE_FC_SHIFTLR 1

/* Define to 1 if your Fortran compiler supports stream IO. */
#define HAVE_FC_STREAM_IO 1

/* Define to 1 if your Fortran compiler supports SYSTEM. */
#define HAVE_FC_SYSTEM 1

/* Define to 1 if you have the FFTW3 library. */
/* #undef HAVE_FFTW3 */

/* Define to 1 if you have a MPI-enabled FFTW3 library. */
/* #undef HAVE_FFTW3_MPI */

/* Define to 1 if you have a threads-enabled FFTW3 library. */
/* #undef HAVE_FFTW3_THREADS */

/* Define to 1 if your Fortran compiler supports Fortran 2003. */
#define HAVE_FORTRAN2003 1

/* Define to 1 if you have a GPU library. */
/* #undef HAVE_GPU */

/* Define to 1 if you have the Cuda library. */
/* #undef HAVE_GPU_CUDA */

/* Define to 1 if you have a Cuda version >= 10 (for nvtx v3). */
/* #undef HAVE_GPU_CUDA10 */

/* Define to 1 if you have a Cuda version < 4. */
/* #undef HAVE_GPU_CUDA3 */

/* Define to 1 if you want to perform double-precision Cuda calculations. */
/* #undef HAVE_GPU_CUDA_DP */

/* Define to 1 if you want to perform single-precision Cuda calculations. */
/* #undef HAVE_GPU_CUDA_SP */

/* Define to 1 if you have a MPI-aware GPU library. */
/* #undef HAVE_GPU_MPI */

/* Define to 1 if you have library nvtx (v3). */
/* #undef HAVE_GPU_NVTX_V3 */

/* Define to 1 if you have a serial GPU library. */
/* #undef HAVE_GPU_SERIAL */

/* Define to 1 if you want to activate double-precision GW calculations. */
/* #undef HAVE_GW_DPC */

/* Define to 1 if you have the HDF5 library. */
#define HAVE_HDF5 1

/* Define to 1 if you have the <hdf5.h> header file. */
#define HAVE_HDF5_H 1

/* Define to 1 if you have a parallel HDF5 library. */
#define HAVE_HDF5_MPI 1

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the Levmar library. */
/* #undef HAVE_LEVMAR */

/* Define to 1 if you have the LibPAW library. */
/* #undef HAVE_LIBPAW */

/* Must be defined to 1 to have LibPAW properly working. */
#define HAVE_LIBPAW_ABINIT 1

/* Define to 1 if you have the LibPSML library. */
/* #undef HAVE_LIBPSML */

/* Define to 1 if you want to activate internal support for libtetra(hedron)
   in ABINIT. */
#define HAVE_LIBTETRA_ABINIT 1

/* Define to 1 if you have the LibXC library. */
#define HAVE_LIBXC 1

/* Define to 1 if you have the ASL linear algebra library. */
/* #undef HAVE_LINALG_ASL */

/* Define to 1 if you have an AXPBY BLAS1 extensions. */
#define HAVE_LINALG_AXPBY 1

/* Define to 1 if you have the ELPA linear algebra library. */
/* #undef HAVE_LINALG_ELPA */

/* Define to 1 if you have ELPA 2011-13 API support. */
/* #undef HAVE_LINALG_ELPA_2013 */

/* Define to 1 if you have ELPA 2014 API support. */
/* #undef HAVE_LINALG_ELPA_2014 */

/* Define to 1 if you have ELPA 2015.02 API support. */
/* #undef HAVE_LINALG_ELPA_2015_02 */

/* Define to 1 if you have ELPA 2015.11 API support. */
/* #undef HAVE_LINALG_ELPA_2015_11 */

/* Define to 1 if you have ELPA 2016 API support. */
/* #undef HAVE_LINALG_ELPA_2016 */

/* Define to 1 if you have ELPA 2017 API support. */
/* #undef HAVE_LINALG_ELPA_2017 */

/* Define to 1 if you have the ELPA Fortran 2008 API support. */
/* #undef HAVE_LINALG_ELPA_FORTRAN2008 */

/* Define to 1 if you have the ESSL linear algebra library. */
/* #undef HAVE_LINALG_ESSL */

/* Define to 1 if you have GEMM3M BLAS3 extensions. */
#define HAVE_LINALG_GEMM3M 1

/* Define to 1 if you have GEMMT BLAS3 extensions. */
/* #undef HAVE_LINALG_GEMMT */

/* Define to 1 if you have the MAGMA linear algebra library. */
/* #undef HAVE_LINALG_MAGMA */

/* Define to 1 if you have MAGMA >=1.5 API support */
/* #undef HAVE_LINALG_MAGMA_15 */

/* Define to 1 if you have MKL imatcopy extensions. */
/* #undef HAVE_LINALG_MKL_IMATCOPY */

/* Define to 1 if you have MKL omatadd extensions. */
/* #undef HAVE_LINALG_MKL_OMATADD */

/* Define to 1 if you have MKL omatcopy extensions. */
/* #undef HAVE_LINALG_MKL_OMATCOPY */

/* Define to 1 if you have mkl_*threads extensions. */
/* #undef HAVE_LINALG_MKL_THREADS */

/* Define to 1 if you have the PLASMA linear algebra library. */
/* #undef HAVE_LINALG_PLASMA */

/* Define to 1 if you have the ScaLAPACK linear algebra library. */
/* #undef HAVE_LINALG_SCALAPACK */

/* Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.
   */
/* #undef HAVE_LINALG_ZDOTC_BUG */

/* Define to 1 if you want to activate workaround for bugged ZDOTC and ZDOTU.
   */
/* #undef HAVE_LINALG_ZDOTU_BUG */

/* Define to 1 if you want to activate LOTF functionality (UNMAINTAINED). */
/* #undef HAVE_LOTF */

/* Define to 1 if you have the `mallinfo' function. */
#define HAVE_MALLINFO 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <malloc/malloc.h> header file. */
/* #undef HAVE_MALLOC_MALLOC_H */

/* Define to 1 if you have the <math.h> header file. */
#define HAVE_MATH_H 1

/* Define to 1 if you have the <mcheck.h> header file. */
#define HAVE_MCHECK_H 1

/* Define to 1 if you want to enable memory profiling. */
/* #undef HAVE_MEM_PROFILING */

/* Define to 1 if you have a working MPI installation. */
#define HAVE_MPI 1

/* Define to 1 if you have a MPI-1 implementation (obsolete, broken). */
/* #undef HAVE_MPI1 */

/* Define to 1 if you have a MPI-2 implementation. */
#define HAVE_MPI2 1

/* Define to 1 if you want to activate support for MPI_IN_PLACE. */
/* #undef HAVE_MPI2_INPLACE */

/* Define to 1 if you have a MPI-3 implementation. */
/* #undef HAVE_MPI3 */

/* Define to 1 if your MPI implementation provides MPI_Get_library_version. */
/* #undef HAVE_MPI_GET_LIBRARY_VERSION */

/* Define to 1 if your MPI library supports MPI_IALLGATHER. */
#define HAVE_MPI_IALLGATHER 1

/* Define to 1 if your MPI library supports MPI_IALLREDUCE. */
#define HAVE_MPI_IALLREDUCE 1

/* Define to 1 if your MPI library supports MPI_IALLTOALL. */
#define HAVE_MPI_IALLTOALL 1

/* Define to 1 if your MPI library supports MPI_IALLTOALLV. */
#define HAVE_MPI_IALLTOALLV 1

/* Define to 1 if your MPI library supports MPI_IBCAST. */
#define HAVE_MPI_IBCAST 1

/* Define to 1 if your MPI library supports MPI_IGATHERV. */
#define HAVE_MPI_IGATHERV 1

/* Define to 1 if you are using XLF. */
/* #undef HAVE_MPI_INCLUDED_ONCE */

/* Define to 1 if your MPI library supports MPI_INTEGER16. */
#define HAVE_MPI_INTEGER16 1

/* Define to 1 if you want MPI I/O support. */
#define HAVE_MPI_IO 1

/* Define to 1 if you want to use MPI I/O as default I/O library
   (maintainer-only option). */
/* #undef HAVE_MPI_IO_DEFAULT */

/* Define to 1 if your MPI library supports MPI_TYPE_CREATE_STRUCT. */
#define HAVE_MPI_TYPE_CREATE_STRUCT 1

/* Define to 1 if you have the NetCDF library. */
#define HAVE_NETCDF 1

/* Define to 1 if you want to use NetCDF as default I/O library
   (maintainer-only option). */
/* #undef HAVE_NETCDF_DEFAULT */

/* Define to 1 if you have the NetCDF Fortran interface library. */
#define HAVE_NETCDF_FORTRAN 1

/* Define to 1 if you have a parallel NetCDF Fortran interface library. */
#define HAVE_NETCDF_FORTRAN_MPI 1

/* Define to 1 if you have a parallel NetCDF library. */
#define HAVE_NETCDF_MPI 1

/* Define to 1 if you have a standard implementation of NumPy. */
/* #undef HAVE_NUMPY */

/* Set to 1 if OpenMP has a working implementation of COLLAPSE. */
/* #undef HAVE_OMP_COLLAPSE */

/* Define to 1 if you want to activate support for OpenMP. */
/* #undef HAVE_OPENMP */

/* Define to 1 if you are using Linux. */
#define HAVE_OS_LINUX 1

/* Define to 1 if you are using MacOS X. */
/* #undef HAVE_OS_MACOSX */

/* Define to 1 if you are using Windows. */
/* #undef HAVE_OS_WINDOWS */

/* Define to 1 if you have the PAPI library. */
/* #undef HAVE_PAPI */

/* Define to 1 if you have the <papi.h> header file. */
/* #undef HAVE_PAPI_H */

/* Define to 1 if you have the PFFT library. */
/* #undef HAVE_PFFT */

/* Define to 1 if you want to activate possibility to call python scripts
   externally by invoking a python interpreter. */
/* #undef HAVE_PYTHON_INVOCATION */

/* Define to 1 if you have the <stddef.h> header file. */
#define HAVE_STDDEF_H 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdio.h> header file. */
#define HAVE_STDIO_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the <sys/ioctl.h> header file. */
#define HAVE_SYS_IOCTL_H 1

/* Define to 1 if you have the <sys/malloc.h> header file. */
/* #undef HAVE_SYS_MALLOC_H */

/* Define to 1 if you have the <sys/resource.h> header file. */
#define HAVE_SYS_RESOURCE_H 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/time.h> header file. */
#define HAVE_SYS_TIME_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <termios.h> header file. */
#define HAVE_TERMIOS_H 1

/* Define to 1 if you want to use the Abinit timer. */
#define HAVE_TIMER_ABINIT 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the TRIQS library. */
/* #undef HAVE_TRIQS */

/* Define to 1 if you want to activate internal support for TRIQS 1.4. */
/* #undef HAVE_TRIQS_v1_4 */

/* Define to 1 if you want to activate internal support for TRIQS 2.0 (This
   option is dominant over the others versions). */
/* #undef HAVE_TRIQS_v2_0 */

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the Wannier90 library. */
/* #undef HAVE_WANNIER90 */

/* Define to 1 if you want to use Wanner90 1.x (awfully obsolete). */
/* #undef HAVE_WANNIER90_V1 */

/* Define to 1 if you want to use LibXML2-based XML I/O. */
/* #undef HAVE_XML */

/* Define to 1 if you have the XMLF90 library. */
/* #undef HAVE_XMLF90 */

/* Define to 1 if assertions should be disabled. */
/* #undef NDEBUG */

/* Name of package */
#define PACKAGE "abinit"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "https://bugs.launchpad.net/abinit/"

/* Define to the full name of this package. */
#define PACKAGE_NAME "ABINIT"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "ABINIT 9.9.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "abinit"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "9.9.0"

/* Define to 1 if you want to tell ABINIT to read file lists from standard
   input. */
/* #undef READ_FROM_FILE */

/* The size of `char', as computed by sizeof. */
#define SIZEOF_CHAR 1

/* The size of `double', as computed by sizeof. */
#define SIZEOF_DOUBLE 8

/* The size of `float', as computed by sizeof. */
#define SIZEOF_FLOAT 4

/* The size of `int', as computed by sizeof. */
#define SIZEOF_INT 4

/* The size of `long', as computed by sizeof. */
#define SIZEOF_LONG 8

/* The size of `long double', as computed by sizeof. */
#define SIZEOF_LONG_DOUBLE 16

/* The size of `long long', as computed by sizeof. */
#define SIZEOF_LONG_LONG 8

/* The size of `ptrdiff_t', as computed by sizeof. */
#define SIZEOF_PTRDIFF_T 8

/* The size of `short', as computed by sizeof. */
#define SIZEOF_SHORT 2

/* The size of `size_t', as computed by sizeof. */
#define SIZEOF_SIZE_T 8

/* The size of `unsigned int', as computed by sizeof. */
#define SIZEOF_UNSIGNED_INT 4

/* The size of `unsigned long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG 8

/* The size of `unsigned long long', as computed by sizeof. */
#define SIZEOF_UNSIGNED_LONG_LONG 8

/* Define to 1 if all of the C90 standard headers exist (not just the ones
   required in a freestanding environment). This macro is provided for
   backward compatibility; new code need not use it. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "9.9.0"

/* Define WORDS_BIGENDIAN to 1 if your processor stores words with the most
   significant byte first (like Motorola and SPARC, unlike Intel). */
#if defined AC_APPLE_UNIVERSAL_BUILD
# if defined __BIG_ENDIAN__
#  define WORDS_BIGENDIAN 1
# endif
#else
# ifndef WORDS_BIGENDIAN
/* #  undef WORDS_BIGENDIAN */
# endif
#endif

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */

/* *** BEGIN sanity checks *** */

/* MPI options */
#if defined HAVE_MPI 

/* Check that one MPI level is actually defined */
#if ! defined HAVE_MPI1 && ! defined HAVE_MPI2 && ! defined HAVE_MPI3
#error "HAVE_MPI1, HAVE_MPI2, and HAVE_MPI3, are all undefined"
#endif

/* Check that only one MPI level has been defined */
#if defined HAVE_MPI1
#  if defined HAVE_MPI2
#    if defined HAVE_MPI3
#      error "HAVE_MPI1, Have_MPI2, and HAVE_MPI3, are all defined"
#    else
#error "HAVE_MPI1 and HAVE_MPI2 are both defined"
#    endif
#  else
#    if defined HAVE_MPI3
#      error "HAVE_MPI1 and HAVE_MPI3 are both defined"
#    endif
#  endif
#else
#  if defined HAVE_MPI2 && defined HAVE_MPI3
#    error "HAVE_MPI2 and HAVE_MPI3 are both defined"
#  endif
#endif

#else /* HAVE_MPI */

/* Check that no MPI level is defined */
#if defined HAVE_MPI1 || defined HAVE_MPI2 || HAVE_MPI3
#error "HAVE_MPI1, HAVE_MPI2, and HAVE_MPI3, must be undefined"
#endif

/* Check that MPI-IO is undefined */
#if defined HAVE_MPI_IO
#error "HAVE_MPI_IO must be undefined"
#endif

#endif /* HAVE_MPI */

/* *** END sanity checks *** */

#endif /* _ABINIT_CONFIG_H */

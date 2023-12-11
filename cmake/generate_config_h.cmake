set(BUILD_WITH_CMAKE 1)

set(ABINIT_VERSION ${PROJECT_VERSION})
set(ABINIT_VERSION_BASE ${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR})
set(ABINIT_VERSION_MAJOR ${PROJECT_VERSION_MAJOR})
set(ABINIT_VERSION_MINOR ${PROJECT_VERSION_MINOR})
set(ABINIT_VERSION_MICRO ${PROJECT_VERSION_PATCH})
string(TIMESTAMP ABINIT_VERSION_BUILD "%Y%m%d")

#
# get compiler ID for C language
#
if (${CMAKE_C_COMPILER_ID} STREQUAL "ARMCC")
  set(CC_ARM 1)
endif()

if (${CMAKE_C_COMPILER_ID} STREQUAL "GNU")
  set(CC_GNU 1)
endif()

if (${CMAKE_C_COMPILER_ID} STREQUAL "XL")
  set(CC_IBM 1)
endif()

if ((${CMAKE_C_COMPILER_ID} STREQUAL "XLClang") OR
    (${CMAKE_C_COMPILER_ID} STREQUAL "IBMClang"))
  set(CC_IBMLLVM 1)
endif()

if ((${CMAKE_C_COMPILER_ID} STREQUAL "Intel") OR
    (${CMAKE_C_COMPILER_ID} STREQUAL "IntelLLVM"))
  set(CC_INTEL 1)
endif()

if (${CMAKE_C_COMPILER_ID} STREQUAL "Clang")
  set(CC_LLVM 1)
endif()

if (${CMAKE_C_COMPILER_ID} STREQUAL "PGI")
  set(CC_PGI 1)
endif()

if (${CMAKE_C_COMPILER_ID} STREQUAL "NVHPC")
  set(CC_NVHPC 1)
endif()

#
# get compiler ID for CXX language
#
if (${CMAKE_CXX_COMPILER_ID} STREQUAL "ARMCC")
  set(CXX_ARM 1)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU")
  set(CXX_GNU 1)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "XL")
  set(CXX_IBM 1)
endif()

if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "XLClang") OR
    (${CMAKE_CXX_COMPILER_ID} STREQUAL "IBMClang"))
  set(CXX_IBMLLVM 1)
endif()

if ((${CMAKE_CXX_COMPILER_ID} STREQUAL "Intel") OR
    (${CMAKE_CXX_COMPILER_ID} STREQUAL "IntelLLVM"))
  set(CXX_INTEL 1)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "Clang")
  set(CXX_LLVM 1)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "PGI")
  set(CXX_PGI 1)
endif()

if (${CMAKE_CXX_COMPILER_ID} STREQUAL "NVHPC")
  set(CXX_NVHPC 1)
endif()

#
# get compiler ID for Fortran language
#
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Absoft")
  set(FC_ABSOFT 1)
endif()

# TODO
if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Flang")
  # EXECUTE_PROCESS(COMMAND ${CMAKE_Fortran_COMPILER} --version 2>/dev/null | head -n 1
  #   OUTPUT_VARIABLE FC_VERSION)
  set(FC_AOCC 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "ARMCC")
  set(FC_ARM 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set(FC_GNU 1)
endif()

if ((${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel") OR
    (${CMAKE_Fortran_COMPILER_ID} STREQUAL "IntelLLVM"))
  set(FC_INTEL 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "NAG")
  set(FC_NAG 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "PGI")
  set(FC_PGI 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "NVHPC")
  set(FC_NVHPC 1)
endif()

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "XL")
  set(FC_IBM 1)
endif()

if ((${CMAKE_Fortran_COMPILER_ID} STREQUAL "XLClang") OR
    (${CMAKE_Fortran_COMPILER_ID} STREQUAL "IBMClang"))
  set(FC_IBMLLVM 1)
endif()

#################################################################

include(CheckTypeSize)

check_type_size(char                 SIZEOF_CHAR)
check_type_size(double               SIZEOF_DOUBLE)
check_type_size(float                SIZEOF_FLOAT)
check_type_size(int                  SIZEOF_INT)
check_type_size(long                 SIZEOF_LONG)
check_type_size("long double"        SIZEOF_LONG_DOUBLE)
check_type_size("long long"          SIZEOF_LONG_LONG)
check_type_size(ptrdiff_t            SIZEOF_PTRDIFF_T)
check_type_size(short                SIZEOF_SHORT)
check_type_size(size_t               SIZEOF_SIZE_T)
check_type_size("unsigned int"       SIZEOF_UNSIGNED_INT)
check_type_size("unsigned long"      SIZEOF_UNSIGNED_LONG)
check_type_size("unsigned long long" SIZEOF_UNSIGNED_LONG_LONG)

include(CheckFunctionExists)
include(CheckLibraryExists)

include(CheckIncludeFile)
include(CheckIncludeFileCXX)
include(CheckIncludeFiles)

include(CheckSymbolExists)

include(CheckFortranSourceRuns)

CHECK_INCLUDE_FILE(dlfcn.h HAVE_DLFCN_H)
CHECK_INCLUDE_FILE(inttypes.h HAVE_INTTYPES_H)
CHECK_INCLUDE_FILE(memory.h HAVE_MEMORY_H)
CHECK_INCLUDE_FILE(stddef.h HAVE_STDDEF_H)
CHECK_INCLUDE_FILE(stdint.h HAVE_STDINT_H)
CHECK_INCLUDE_FILE_CXX(cstdint HAVE_CSTDINT)
CHECK_INCLUDE_FILE(stdio.h HAVE_STDIO_H)
CHECK_INCLUDE_FILE(stdlib.h HAVE_STDLIB_H)
CHECK_INCLUDE_FILE(strings.h HAVE_STRINGS_H)
CHECK_INCLUDE_FILE(string.h HAVE_STRING_H)
CHECK_INCLUDE_FILE(sys/ioctl.h HAVE_SYS_IOCTL_H)
CHECK_INCLUDE_FILE(sys/resource.h HAVE_SYS_RESOURCE_H)
CHECK_INCLUDE_FILE(sys/stat.h HAVE_SYS_STAT_H)
CHECK_INCLUDE_FILE(sys/time.h HAVE_SYS_TIME_H)
CHECK_INCLUDE_FILE(sys/types.h HAVE_SYS_TYPES_H)
CHECK_INCLUDE_FILE(termios.h HAVE_TERMIOS_H)
CHECK_INCLUDE_FILE(time.h HAVE_TIME_H)
CHECK_INCLUDE_FILE(unistd.h HAVE_UNISTD_H)
CHECK_INCLUDE_FILE(errno.h HAVE_ERRNO_H)

#/* Define to 1 if you have the ANSI C header files. */
CHECK_INCLUDE_FILES("stdlib.h;stdarg.h;string.h;float.h" HAVE_STDC_HEADERS)

#
# check fortran compiler features
#

# TODO: I'm really not convinced that all this is really necessary, but for compatibility
# with the old autotools build system, let's check exhaustively all these features support

check_fortran_source_runs(
  "program test
     integer, parameter :: dp=kind(1.0d0)
     integer, parameter :: dpc=kind((1.0_dp,1.0_dp))

     type test_type
        integer,allocatable :: i(:)
        real(dp),allocatable :: r(:,:)
        complex(dpc),allocatable :: c(:,:,:)
     end type test_type
   end program test"
  HAVE_FC_ALLOCATABLE_DTARRAYS
  SRC_EXT f90)

# check async
try_compile(HAVE_FC_ASYNC_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_async.F90)
if (HAVE_FC_ASYNC_BOOL)
  set(HAVE_FC_ASYNC 1)
endif()

# check backtrace
try_compile(HAVE_FC_BACKTRACE_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_backtrace.F90)
if (HAVE_FC_BACKTRACE_BOOL)
  set(HAVE_FC_BACKTRACE 1)
endif()

# check command_argument
try_compile(HAVE_FC_COMMAND_ARGUMENT_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_command_argument.F90)
if (HAVE_FC_COMMAND_ARGUMENT_BOOL)
  set(HAVE_FC_COMMAND_ARGUMENT 1)
endif()

check_fortran_source_runs(
  "program test
     implicit none
     integer, parameter :: dp=kind(1.0d0)
     integer, parameter :: dpc=kind((1.0_dp,1.0_dp))

     integer,contiguous,pointer :: i_ptr(:)
     real(dp),contiguous,pointer :: r_ptr(:,:)
     complex(dpc),contiguous,pointer :: c_ptr(:,:,:)
   end program test"
  HAVE_FC_CONTIGUOUS
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
   end program test"
  HAVE_FC_CPUTIME
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     real, dimension(2) :: tarray
     real :: result
     call etime(tarray, result)
   end program test"
  HAVE_FC_ETIME
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     integer :: STATUS = 0
     call exit(STATUS)
   end program test"
  HAVE_FC_EXIT
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     call flush(1)
   end program test"
  HAVE_FC_FLUSH
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     call flush_(1)
   end program test"
  HAVE_FC_FLUSH_
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     real :: x
     x = gamma(1.5)
   end program test"
  HAVE_FC_GAMMA
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     character(len=255) :: homedir
     call getenv(\"HOME\", homedir)
   end program test"
  HAVE_FC_GETENV
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
   end program test"
  HAVE_FC_GETPID
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     use, intrinsic :: ieee_arithmetic
     implicit none
     real :: val

     if (ieee_is_nan(val)) then  ! NaN
       write(*,*) \"Hello NAN\"
     end if
   end program test"
  HAVE_FC_IEEE_ARITHMETIC
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     use, intrinsic :: ieee_exceptions
     implicit none
     type(ieee_status_type) :: status_value

     call ieee_get_status(status_value)   ! Get the flags
     call ieee_set_flag(ieee_all,.false.) ! Set the flags quiet
     call ieee_set_status(status_value)   ! Restore the flags
   end program test"
  HAVE_FC_IEEE_EXCEPTIONS
  SRC_EXT f90)

# TODO : not sure of the purpose of this one, if does not compile
# since -Werror=line-truncation
# check_fortran_source_runs(
#   "program test
#      implicit none
#      write(*,*) \"aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa\" !142
#    end program test"
#   HAVE_FC_LONG_LINES
#   SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
#        define NEWLINE \newline
     print *, \"foo1\" NEWLINE print *,\"foo2\"
   end program test"
  HAVE_FC_MACRO_NEWLINE
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     integer, allocatable :: a(:), b(:)
     allocate(a(3))
     a = (/1, 2, 3/)
     call move_alloc(a, b)
   end program test"
  HAVE_FC_MOVE_ALLOC
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
      implicit none
      integer :: x(3,4)
      call test_on_the_fly_shape(x)
      contains
      subroutine test_on_the_fly_shape(x)
        integer, intent(inout) :: x(:,:)
        integer :: y(product(shape(x)))
      end subroutine test_on_the_fly_shape
   end program test"
  HAVE_FC_ON_THE_FLY_SHAPE
  SRC_EXT f90)

check_fortran_source_runs(
  "module foo
    type, public :: bar_t
      integer :: pub
      integer,private :: priv
    end type bar_t
   end module foo

   program test
     use foo
     type(bar_t) :: data
   end program test"
  HAVE_FC_PRIVATE
  SRC_EXT f90)

check_fortran_source_runs(
  "module foo
     real,save,protected :: aprot(10)
   end module foo

   program test
    use foo
   end program test"
  HAVE_FC_PROTECTED
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     integer :: ii,ishft,res
     res=shiftl(ii,ishft)
     res=shiftr(ii,ishft)
   end program test"
  HAVE_FC_SHIFTLR
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
         IMPLICIT NONE
         INTEGER :: myvalue = 12345, mypos
         OPEN(UNIT=11, FILE=\"ustream.demo\", STATUS=\"NEW\", ACCESS=\"STREAM\")
         WRITE(11) \"first\"
         WRITE(11) \"second\"
         INQUIRE(UNIT=11, POS=mypos)
         PRINT *, \"Myvalue will be written at position \", mypos
         WRITE(11) myvalue
         CLOSE(UNIT=11)
   end program test"
  HAVE_FC_STREAM_IO
  SRC_EXT f90)

check_fortran_source_runs(
  "program test
     implicit none
     call system (\"ls -l\")
   end program test"
  HAVE_FC_SYSTEM
  SRC_EXT f90)

# check iso_c_binding
try_compile(HAVE_FC_ISO_C_BINDING_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_iso_c_binding.F90)
if (HAVE_FC_ISO_C_BINDING_BOOL)
  set(HAVE_FC_ISO_C_BINDING 1)
endif()

# check iso_c_fortran_env
try_compile(HAVE_FC_ISO_FORTRAN_2008_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_iso_fortran_2008.F90)
if (HAVE_FC_ISO_FORTRAN_2008_BOOL)
  set(HAVE_FC_ISO_FORTRAN_2008 1)
endif()

# check integer*16
try_compile(HAVE_FC_INT_QUAD_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_int_quad.F90)
if (HAVE_FC_INT_QUAD_BOOL)
  set(HAVE_FC_INT_QUAD 1)
endif()

# check IOMSG
try_compile(HAVE_FC_IOMSG_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_fc_iomsg.F90)
if (HAVE_FC_IOMSG_BOOL)
  set(HAVE_FC_IOMSG 1)
endif()

#
# MKL features check
#
if (MKL_FOUND)
  try_compile(HAVE_LINALG_MKL_IMATCOPY_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_linalg_mkl_imatcopy.F90
    LINK_LIBRARIES MKL::MKL)

  if (HAVE_LINALG_MKL_IMATCOPY_BOOL)
    set(HAVE_LINALG_MKL_IMATCOPY 1)
  endif()

  try_compile(HAVE_LINALG_MKL_OMATCOPY_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_linalg_mkl_omatcopy.F90
    LINK_LIBRARIES MKL::MKL)

  if (HAVE_LINALG_MKL_OMATCOPY_BOOL)
    set(HAVE_LINALG_MKL_OMATCOPY 1)
  endif()

  try_compile(HAVE_LINALG_MKL_OMATADD_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_linalg_mkl_omatadd.F90
    LINK_LIBRARIES MKL::MKL)

  if (HAVE_LINALG_MKL_OMATADD_BOOL)
    set(HAVE_LINALG_MKL_OMATADD 1)
  endif()

endif(MKL_FOUND)


#
# MPI features check
#
if (MPI_FOUND)

  set(HAVE_MPI 1)

  message("[MPI] MPI_Fortran_VERSION_MAJOR = ${MPI_Fortran_VERSION_MAJOR}")

  # the following reproduced the behavior existing in m4 macro
  # HAVE_MPI2 is used with the meaning MPI standard version (major) is greater if equal to 2
  if (MPI_Fortran_VERSION_MAJOR VERSION_GREATER_EQUAL "2")
    set(HAVE_MPI2 1)
  endif()

  if (MPI_Fortran_VERSION_MAJOR VERSION_GREATER_EQUAL "3")
    set(HAVE_MPI3 1)
  endif()

  # check mpi_integer16
  try_compile(HAVE_MPI_INTEGER16_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_mpi_integer16.F90
    LINK_LIBRARIES MPI::MPI_Fortran)
  if (HAVE_MPI_INTEGER16_BOOL)
    set(HAVE_MPI_INTEGER16 1)
  endif()

  try_compile(HAVE_MPI_GET_LIBRARY_VERSION_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_mpi_get_library_version.F90
    LINK_LIBRARIES MPI::MPI_Fortran)

  if(HAVE_MPI_GET_LIBRARY_VERSION_BOOL)
    set(HAVE_MPI_GET_LIBRARY_VERSION 1)
  endif()

  # check runtime GPU-awareness ability
  include(CheckMPIFeatures)

endif()

#
# NetCDF / MPI
#
if (MPI_FOUND)

  if (ABINIT_NETCDF_FOUND)

    # check NetCDF C / MPI
    # message(STATUS "Check netcdf/c/mpi support...")
    # try_compile(HAVE_NETCDF_C_MPI_BOOL ${CMAKE_BINARY_DIR}/try_compile ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_netcdf_c_mpi.c
    #   LINK_LIBRARIES MPI::MPI_C PkgConfig::ABINIT_NETCDF)
    # #PkgConfig::ABINIT_NETCDF_MPI)
    # if (HAVE_NETCDF_C_MPI_BOOL)
    #   message("NetCDF C/MPI support checked ok")
    #   set(HAVE_NETCDF_MPI 1)
    # else()
    #   message("NetCDF C/MPI not supported")
    # endif()

    message(STATUS "Check netcdf/fortran support...")
    if (NOT DEFINED HAVE_NETCDF_FORTRAN_BOOL)
      try_compile(HAVE_NETCDF_FORTRAN_BOOL
        ${CMAKE_BINARY_DIR}/try_compile_netcdf
        ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_netcdf_fortran.F90
        LINK_LIBRARIES PkgConfig::ABINIT_NETCDF_FORTRAN PkgConfig::ABINIT_NETCDF MPI::MPI_Fortran)
    endif()
    if (HAVE_NETCDF_FORTRAN_BOOL)
      message("NetCDF Fortran support checked ok")
      set(HAVE_NETCDF 1)
      set(HAVE_NETCDF_FORTRAN 1)
    else()
      message("NetCDF Fortran not supported")
    endif()

    # check NetCDF fortran / MPI
    message(STATUS "Check netcdf/fortran/mpi support...")
    if (NOT DEFINED HAVE_NETCDF_FORTRAN_MPI_BOOL)
      # try_run(
      #   HAVE_NETCDF_FORTRAN_MPI_BOOL_RUN
      #   HAVE_NETCDF_FORTRAN_MPI_BOOL_COMPILE
      #   ${CMAKE_BINARY_DIR}/try_compile
      #   ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_netcdf_fortran_mpi.F90
      #   LINK_LIBRARIES MPI::MPI_Fortran PkgConfig::ABINIT_NETCDF_FORTRAN
      #   RUN_OUTPUT_VARIABLE HAVE_NETCDF_FORTRAN_MPI_RUN_RES)
      try_compile(
        HAVE_NETCDF_FORTRAN_MPI_BOOL
        ${CMAKE_BINARY_DIR}/try_compile_netcdf_mpi
        ${CMAKE_SOURCE_DIR}/cmake/try_compile/have_netcdf_fortran_mpi.F90
        LINK_LIBRARIES MPI::MPI_Fortran PkgConfig::ABINIT_NETCDF_FORTRAN)
    endif()

    if (HAVE_NETCDF_FORTRAN_BOOL AND HAVE_NETCDF_FORTRAN_MPI_BOOL)
      message("NetCDF fortran/MPI support checked ok")
      set(HAVE_NETCDF_FORTRAN_MPI 1)

      # TO DO : evaluate if we should also set HAVE_NETCDF_MPI here, because actually
      # only HAVE_NETCDF_MPI is used e.g m_nctk.F90
      set(HAVE_NETCDF_MPI 1)
    else()
      message("NetCDF fortran/MPI not supported")
    endif()

  endif()

endif()

# TODO PK : evaluate if this really need to be checked, almost all MPI implementation
# have support for MPI2
set(HAVE_MPI_IALLGATHER 1)
set(HAVE_MPI_IALLREDUCE 1)
set(HAVE_MPI_IALLTOALL 1)
set(HAVE_MPI_IALLTOALLV 1)
set(HAVE_MPI_IBCAST 1)
set(HAVE_MPI_IGATHERV 1)
set(HAVE_MPI_TYPE_CREATE_STRUCT 1)

#
# System
#
if(CMAKE_SYSTEM_NAME STREQUAL Linux)
  set(HAVE_OS_LINUX 1)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL Darwin)
  set(HAVE_OS_MACOSX 1)
endif()

if(CMAKE_SYSTEM_NAME STREQUAL Windows)
  set(HAVE_OS_WINDOWS 1)
endif()

# TODO : other OS (FreeBSD, Android, MSYS ?)

#
# GPU features check
#

# check if cublas.h is available
if (TARGET CUDA::cublas)
  set(HAVE_CUBLAS_H 1)
endif()

if (TARGET CUDA::cudart)
  set(HAVE_CUDA_RUNTIME_API_H 1)
endif()

if (TARGET CUDA::cudafft)
  set(HAVE_CUFFT_H 1)
endif()

CHECK_SYMBOL_EXISTS(abort "stdlib.h" HAVE_ABORT)
CHECK_SYMBOL_EXISTS(clock_gettime "time.h" HAVE_CLOCK_GETTIME)

message(" Generating config.h ...")
configure_file(config.h.cmake config.h @ONLY)
message("")

# -*- Autoconf -*-
#
# Copyright (C) 2005-2019 ABINIT Group (Yann Pouillon)
#
# This file is part of the ABINIT software package. For license information,
# please see the COPYING file in the top-level directory of the ABINIT source
# distribution.
#

#
# MPI support for ABINIT
#


# _ABI_MPI_CHECK_CC()
# -------------------
#
# Checks whether the C compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_CC], [
  # Set default values
  abi_mpi_cc_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${CC_LDFLAGS}"
  LIBS="${CC_LIBS} ${abi_mpi_libs}"

  # Try to compile a C MPI program
  AC_MSG_CHECKING([whether the C compiler supports MPI])
  AC_LANG_PUSH([C])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
#include <stdlib.h>
#include <mpi.h>
    ]],
    [[
      int rc;

      MPI_Init(NULL, NULL);
      rc = MPI_Finalize();
    ]])], [abi_mpi_cc_ok="yes"], [abi_mpi_cc_ok="no"])
  AC_LANG_POP([C])
  AC_MSG_RESULT([${abi_mpi_cc_ok}])

  # Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_CC


                    # ------------------------------------ #


# _ABI_MPI_CHECK_CXX()
# --------------------
#
# Checks whether the C++ compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_CXX], [
  # Set default values
  abi_mpi_cxx_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${CXX_LDFLAGS}"
  LIBS="${CXX_LIBS} ${abi_mpi_libs}"

  # Try to compile a C++ MPI program
  AC_MSG_CHECKING([whether the C++ compiler supports MPI])
  AC_LANG_PUSH([C++])
  AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[@%:@include <mpi.h>]],
    [[
      MPI::Init();
      MPI::Finalize();
    ]])], [abi_mpi_cxx_ok="yes"], [abi_mpi_cxx_ok="no"])
  AC_LANG_POP([C++])
  AC_MSG_RESULT([${abi_mpi_cxx_ok}])

  # Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_CXX


                    # ------------------------------------ #


# _ABI_MPI_CHECK_FC()
# -------------------
#
# Checks whether the Fortran compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_FC], [
  # Set default values
  abi_mpi_fc_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  # Try to compile a Fortran MPI program
  AC_MSG_CHECKING([whether the Fortran Compiler supports MPI])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[
      include "mpif.h"
      integer :: ierr
      call mpi_init(ierr)
      call mpi_finalize(ierr)
    ]])], [abi_mpi_fc_ok="yes"], [abi_mpi_fc_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_fc_ok}])

  # Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_FC


                    # ------------------------------------ #


# _ABI_MPI_CHECK_FC_LEVEL()
# -------------------------
#
# Checks which MPI level is supported by the Fortran compiler.
#
AC_DEFUN([_ABI_MPI_CHECK_FC_LEVEL], [
  if test "${sd_mpi_fc_ok}" = "yes"; then

    # Back-up build environment
    ABI_ENV_BACKUP

    # Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    # Try to compile a MPI-2 Fortran program
    AC_MSG_CHECKING([which level of MPI is supported by the Fortran compiler])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
              use mpi
              integer :: ierr
              call mpi_init(ierr)
              call mpi_finalize(ierr)
      ]])], [abi_mpi_fc_level="2"], [abi_mpi_fc_level="1"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_mpi_fc_level}])

    # Restore build environment
    ABI_ENV_RESTORE

  fi # sd_mpi_fc_ok = yes
]) # _ABI_MPI_CHECK_FC_LEVEL


                    # ------------------------------------ #


# _ABI_MPI_CHECK_INTEGER16()
# --------------------------
#
# Checks whether the MPI library supports MPI_INTEGER16
#
AC_DEFUN([_ABI_MPI_CHECK_INTEGER16], [
  # Set default values
  abi_mpi_integer16_ok="no"


  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${lib_mpi_libs}"

  # Try to compile a Fortran program
  # Note: we assume a MPI implementation that provides the mpi module
  AC_MSG_CHECKING([whether the MPI library supports MPI_INTEGER16])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer, parameter :: ii = MPI_INTEGER16

    ]])], [abi_mpi_integer16_ok="yes"], [abi_mpi_integer16_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_integer16_ok}])

  # Restore build environment
  ABI_ENV_RESTORE

  # Forward information to the compiler
  if test "${abi_mpi_integer16_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_INTEGER16], 1,
      [Define to 1 if your MPI library supports MPI_INTEGER16.])
  fi
]) # _ABI_MPI_CHECK_INTEGER16


                    # ------------------------------------ #


# _ABI_MPI_CHECK_GET_LIBRARY_VERSION()
# ------------------------------------
#
# Checks whether the MPI library provides MPI_Get_library_version.
#
AC_DEFUN([_ABI_MPI_CHECK_GET_LIBRARY_VERSION], [
  # Set default values
  abi_mpi_get_library_version="unknown"

  if test "${abi_mpi_fc_ok}" = "yes"; then

    # Back-up build environment
    ABI_ENV_BACKUP

    # Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    AC_MSG_CHECKING([whether MPI provides MPI_Get_library_version])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
              use mpi
              character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: info
              integer :: ilen, ierr
              call mpi_init(ierr)
              call mpi_get_library_version(info, ilen, ierr)
              call mpi_finalize(ierr)
      ]])],
      [abi_mpi_get_library_version="yes"], [abi_mpi_get_library_version="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_mpi_get_library_version}])

    # Restore build environment
    ABI_ENV_RESTORE
  fi
]) # _ABI_MPI_CHECK_GET_LIBRARY_VERSION


                    # ------------------------------------ #


# _ABI_MPI_CHECK_CREATE_TYPE_STRUCT()
# -----------------------------------
#
# Checks whether the MPI library supports MPI_CREATE_TYPE_STRUCT (MPI2).
#
AC_DEFUN([_ABI_MPI_CHECK_CREATE_TYPE_STRUCT], [
  # Set default values
  abi_mpi_type_create_struct_ok="no"

  if test "${abi_mpi_fc_level}" -ge "2"; then

    # No problem should appear for MPI2 or MPI3 but we test it anyway.

    # Back-up build environment
    ABI_ENV_BACKUP

    # Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    # Try to compile a Fortran MPI program
    AC_MSG_CHECKING([whether the MPI library supports MPI_CREATE_TYPE_STRUCT])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        use mpi
        integer, parameter :: ncount=10
        integer :: ierr, new_type
        integer :: block_length(ncount), block_type(ncount)
        integer(MPI_ADDRESS_KIND) :: block_displ(ncount)
        call mpi_init(ierr)
        call MPI_TYPE_CREATE_STRUCT(ncount, block_length, block_displ, &
&         block_type, new_type, ierr)
        call mpi_finalize(ierr)
      ]])], [abi_mpi_type_create_struct_ok="yes"], [abi_mpi_type_create_struct_ok="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_mpi_type_create_struct_ok}])

    # Restore build environment
    ABI_ENV_RESTORE

  else

    # Back-up build environment
    ABI_ENV_BACKUP

    # Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    # Try to compile a Fortran MPI program
    AC_MSG_CHECKING([whether the MPI library supports MPI_CREATE_TYPE_STRUCT])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM([],
      [[
        include "mpif.h"
        integer, parameter :: ncount=10
        integer :: ierr, new_type
        integer :: block_length(ncount), block_type(ncount)
        integer(MPI_ADDRESS_KIND) :: block_displ(ncount)
        call mpi_init(ierr)
        call MPI_TYPE_CREATE_STRUCT(ncount, block_length, block_displ, &
&         block_type, new_type, ierr)
        call mpi_finalize(ierr)
      ]])], [abi_mpi_type_create_struct_ok="yes"], [abi_mpi_type_create_struct_ok="no"])
    AC_LANG_POP([Fortran])
    AC_MSG_RESULT([${abi_mpi_type_create_struct_ok}])

    # Restore build environment
    ABI_ENV_RESTORE

  fi

  # Forward information to the compiler
  if test "${abi_mpi_type_create_struct_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_TYPE_CREATE_STRUCT], 1,
      [Define to 1 if your MPI library supports MPI_TYPE_CREATE_STRUCT.])
  fi
]) # _ABI_MPI_CHECK_CREATE_TYPE_STRUCT


                    # ------------------------------------ #


# _ABI_MPI_CHECK_IBCAST()
# -----------------------
#
# Checks whether the MPI library supports MPI_IBCAST (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IBCAST],[
  dnl Set default values
  abi_mpi_ibcast_ok="no"

  dnl Try to compile a Fortran program
  AC_MSG_CHECKING([whether the MPI library supports MPI_IBCAST (MPI3)])

  dnl We assume a MPI implementation that provides the mpi module

  dnl Back-up build environment
  ABI_ENV_BACKUP
                                                                                            
  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${lib_mpi_libs}"

  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer,parameter :: siz=5
      integer :: buffer(siz)
      integer :: comm, root, request, ierr

      ! Prototype
      !  MPI_IBCAST(buffer, count, datatype, root, comm, request, ierr)
      !  INOUT buffer	starting address of buffer (choice)
      !  IN count	number of entries in buffer (non-negative integer)
      !  IN datatype	data type of buffer (handle)
      !  IN root	rank of broadcast root (integer)
      !  IN comm	communicator (handle)
      !  OUT request	communication request (handle)

      call MPI_IBCAST(buffer, siz, MPI_INTEGER, root, comm, request, ierr) 

    ]])], [abi_mpi_ibcast_ok="yes"], [abi_mpi_ibcast_ok="no"])
  AC_LANG_POP
                                                                                            
  dnl Restore build environment
  ABI_ENV_RESTORE

  AC_MSG_RESULT([${abi_mpi_ibcast_ok}])

  if test "${abi_mpi_ibcast_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IBCAST],1,
      [Define to 1 if your MPI library supports MPI_IBCAST.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi

]) # _ABI_MPI_CHECK_IBCAST     


                    ########################################


# _ABI_MPI_CHECK_IALLTOALL()
# --------------------------
#
# Checks whether the MPI library supports MPI_IALLTOALL (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IALLTOALL], [
  # Set default values
  abi_mpi_ialltoall_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  # Try to compile a Fortran program
  # Note: we assume a MPI implementation that provides the mpi module
  AC_MSG_CHECKING([whether the MPI library supports MPI_IALLTOALL (MPI3)])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer, parameter :: siz=5
      integer :: SENDBUF(siz), RECVBUF(siz)
      integer :: SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
      integer :: COMM, REQUEST, IERROR

      ! Prototype
      ! <type>    SENDBUF(*), RECVBUF(*)
      ! INTEGER    SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
      ! INTEGER    COMM, REQUEST, IERROR

      call MPI_IALLTOALL(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, &
      &  RECVTYPE, COMM, REQUEST, IERROR)

    ]])], [abi_mpi_ialltoall_ok="yes"], [abi_mpi_ialltoall_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_ialltoall_ok}])

  # Restore build environment
  ABI_ENV_RESTORE

  # Forward information to the compiler
  if test "${abi_mpi_ialltoall_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IALLTOALL], 1,
      [Define to 1 if your MPI library supports MPI_IALLTOALL.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi
]) # _ABI_MPI_CHECK_IALLTOALL     


                    # ------------------------------------ #


# _ABI_MPI_CHECK_IALLTOALLV()
# ---------------------------
#
# Checks whether the MPI library supports MPI_IALLTOALLV (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IALLTOALLV], [
  # Set default values
  abi_mpi_ialltoallv_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  # Try to compile a Fortran program
  # Note: we assume a MPI implementation that provides the mpi module
  AC_MSG_CHECKING([whether the MPI library supports MPI_IALLTOALLV (MPI3)])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer, parameter :: siz=5, group_size=3
      integer :: SENDBUF(siz), RECVBUF(siz)
      integer :: SENDCOUNTS(group_size), SDISPLS(group_size)
      integer :: RECVCOUNTS(group_size), RDISPLS(group_size)
      integer :: SENDTYPE, RECVTYPE
      integer :: COMM, REQUEST, IERROR

      call MPI_IALLTOALLV(SENDBUF, SENDCOUNTS, SDISPLS, SENDTYPE, &
&       RECVBUF, RECVCOUNTS, RDISPLS, RECVTYPE, COMM, REQUEST, IERROR)

      !call MPI_Ialltoallv(const void *sendbuf, const int sendcounts[], const int sdispls[],
      !MPI_Datatype sendtype, void *recvbuf, const int recvcounts[],
      !  const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm,
      !  MPI_Request *request)
    ]])], [abi_mpi_ialltoallv_ok="yes"], [abi_mpi_ialltoallv_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_ialltoallv_ok}])

  # Restore build environment
  ABI_ENV_RESTORE

  # Forward information to the compiler
  if test "${abi_mpi_ialltoallv_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IALLTOALLV], 1,
      [Define to 1 if your MPI library supports MPI_IALLTOALLV.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi
]) # _ABI_MPI_CHECK_IALLTOALLV


                    # ------------------------------------ #


# _ABI_MPI_CHECK_IGATHERV()
# -------------------------
#
# Checks whether the MPI library supports MPI_IGATHERV (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IGATHERV],[
  # Set default values
  abi_mpi_igatherv_ok="no"

  # Try to compile a Fortran program
  AC_MSG_CHECKING([whether the MPI library supports MPI_IGATHERV (MPI3)])

  # We assume a MPI implementation that provides the mpi module

  # Back-up build environment
  ABI_ENV_BACKUP
                                                                                            
  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${lib_mpi_libs}"

  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer,parameter :: siz=5,group_size=3
      integer :: SENDBUF(siz), RECVBUF(siz)
      integer :: SDISPLS(group_size), RECVCOUNTS(group_size),RDISPLS(group_size)
      integer :: SENDCOUNTS, SENDTYPE, RECVTYPE
      integer :: ROOT=0, COMM, REQUEST, IERROR

      call MPI_IGATHERV(SENDBUF,SENDCOUNTS,SENDTYPE,&
        RECVBUF,RECVCOUNTS,RDISPLS,RECVTYPE,ROOT,COMM,REQUEST,IERROR)

      !int MPI_Igatherv(const void *sendbuf, int sendcount, MPI_Datatype sendtype, void *recvbuf,
      !                 const int recvcounts[], const int displs[], MPI_Datatype recvtype, int root,
      !                 MPI_Comm comm, MPI_Request * request)

    ]])], [abi_mpi_igatherv_ok="yes"], [abi_mpi_igatherv_ok="no"])
  AC_LANG_POP
                                                                                            
  # Restore build environment
  ABI_ENV_RESTORE

  AC_MSG_RESULT([${abi_mpi_igatherv_ok}])

  if test "${abi_mpi_igatherv_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IGATHERV],1,
      [Define to 1 if your MPI library supports MPI_IGATHERV.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi

]) # _ABI_MPI_CHECK_IGATHERV


                    ########################################


# _ABI_MPI_CHECK_IALLREDUCE()
# ---------------------------
#
# Checks whether the MPI library supports MPI_IALLREDUCE (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IALLREDUCE], [
  # Set default values
  abi_mpi_iallreduce_ok="no"

  # Back-up build environment
  ABI_ENV_BACKUP

  # Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  # Try to compile a Fortran MPI program
  # Note: we assume a MPI implementation that provides the mpi module
  AC_MSG_CHECKING([whether the MPI library supports MPI_IALLREDUCE (MPI3)])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer, parameter :: count=5
      integer :: SENDBUF(count), RECVBUF(count)
      integer :: SENDCOUNT, SENDTYPE, RECVCOUNT, RECVTYPE
      integer :: COMM, REQUEST, IERROR

      ! Prototype
      !<type> SENDBUF(*), RECVBUF(*) INTEGER COUNT, DATATYPE, OP, COMM, REQUEST, IERROR

      call MPI_IALLREDUCE(sendbuf, recvbuf, count, MPI_INTEGER, MPI_SUM, comm, request, ierror)
    ]])], [abi_mpi_iallreduce_ok="yes"], [abi_mpi_iallreduce_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_iallreduce_ok}])

  # Restore build environment
  ABI_ENV_RESTORE

  # Forward information to the compiler
  if test "${abi_mpi_iallreduce_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IALLREDUCE], 1,
      [Define to 1 if your MPI library supports MPI_IALLREDUCE.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi
]) # _ABI_MPI_CHECK_IALLREDUCE     


                    # ------------------------------------ #


# _ABI_MPI_CREATE_WRAPPER(COMPILER_TYPE, SERIAL_COMPILER, MPI_COMPILER)
# ---------------------------------------------------------------------
#
# Creates a wrapper for MPI compilers when they can be configured to
# accept different sequential compilers.
#
# Note: this is impossible with the Autotools, because they require CC
#       CXX, and FC to be set to the actual compilers.
#
AC_DEFUN([_ABI_MPI_CREATE_WRAPPER], [
  # Init
  tmp_name=`echo "$1" | sed -e 's/.*/\L&/'`
  ${MKDIR_P} config/wrappers

  # Create file
  cat >${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name} <<EOF
#!/bin/sh

$1="$2"
export $1

$3 \[$]{*}
EOF

  # Fix permissions
  chmod u+x ${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}

  # Overwrite compiler setting
  eval $1="${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}"

  # Clean-up the mess
  unset tmp_name
]) # _ABI_MPI_CREATE_WRAPPER


                    # ------------------------------------ #


# ABI_MPI_INIT()
# --------------
#
# Looks for an implementation of MPI, using the provided prefix.
# Note 1: this is a convenience feature, purely for comfort.
# Note 2: should be run as early as possible.
#
AC_DEFUN([ABI_MPI_INIT], [
  # Delegate most of the init stage to Steredeg
  SD_MPI_INIT([optional fail], [-lmpi])

  # Allow MPI flavors
  AC_ARG_WITH([mpi-flavor],
    AC_HELP_STRING([--with-mpi-flavor],
      [Flavor of MPI to use (default: auto),
       see ~abinit/doc/build/config-template.ac9 for details]),
    [abi_mpi_flavor="${withval}"],
    [abi_mpi_flavor="auto"])

  # Init ABINIT MPI variables
  abi_mpi_enable="${sd_mpi_enable}"
  abi_mpi_fc_level="${abi_mpi_level}"
  abi_mpi_get_library_version="unknown"
  abi_mpi_init="${sd_mpi_init}"

  # Init ABINIT MPI build flags
  abi_mpi_cppflags=""
  abi_mpi_cflags=""
  abi_mpi_cxxflags=""
  abi_mpi_fcflags=""
  abi_mpi_ldflags=""
  abi_mpi_libs=""

  # Init ABINIT MPI build parameters
  if test "${abi_mpi_enable}" = "yes" -o "${abi_mpi_enable}" = "auto"; then

    if test "${with_mpi_level}" != ""; then
      AC_MSG_WARN([forcing MPI level to ${with_mpi_level} might make the build fail])
    fi

    if test "${sd_mpi_cc_set}" != "yes" -o "${sd_mpi_fc_set}" != "yes"; then
      abi_mpi_cppflags="${sd_mpi_cppflags}"
      abi_mpi_cflags="${sd_mpi_cflags}"
      abi_mpi_cxxflags="${sd_mpi_cxxflags}"
      abi_mpi_fcflags="${sd_mpi_fcflags}"
      abi_mpi_ldflags="${sd_mpi_ldflags}"
      abi_mpi_libs="${sd_mpi_libs}"
    fi

  else

    if test "${abi_mpi_init}" != "def"; then
      AC_MSG_NOTICE([MPI support disabled from command-line])
    fi
    abi_mpi_enable="no"
    abi_mpi_get_library_version="no"
    abi_mpi_init="cmd"
    abi_mpi_inplace_enable="no"
    abi_mpi_io_enable="no"
    abi_mpi_cppflags=""
    abi_mpi_cflags=""
    abi_mpi_cxxflags=""
    abi_mpi_fcflags=""
    abi_mpi_ldflags=""
    abi_mpi_level=""
    abi_mpi_incs=""
    abi_mpi_libs=""

  fi # abi_mpi_enable

  # Check that a permitted flavor value has been specified
  if test "${abi_mpi_enable}" = "yes"; then
    flavor_values=`echo "${abi_mpi_flavor}" | sed -e 's/+/ /g'`
    tmp_flavor_ok="no"
    AC_MSG_CHECKING([whether the '${abi_mpi_flavor}' MPI flavor is valid])
    for chk_flavor in auto double-wrap flags native prefix; do
      for set_flavor in ${flavor_values}; do
        if test "${set_flavor}" = "${chk_flavor}"; then
          tmp_flavor_ok="yes"
          break
        fi
      done
      test "${tmp_flavor_ok}" = "yes" && break
    done
    AC_MSG_RESULT([${tmp_flavor_ok}])
    if test "${tmp_flavor_ok}" = "no"; then
      AC_MSG_ERROR([invalid MPI flavor: '${abi_mpi_flavor}'])
    fi
  else
    abi_mpi_flavor="none"
  fi

  # Enable substitution
  AC_SUBST(abi_mpi_enable)
  AC_SUBST(abi_mpi_flavor)
  AC_SUBST(abi_mpi_cppflags)
  AC_SUBST(abi_mpi_cflags)
  AC_SUBST(abi_mpi_cxxflags)
  AC_SUBST(abi_mpi_fcflags)
  AC_SUBST(abi_mpi_ldflags)
  AC_SUBST(abi_mpi_level)
  AC_SUBST(abi_mpi_incs)
  AC_SUBST(abi_mpi_libs)
]) # ABI_MPI_INIT


                    # ------------------------------------ #


# ABI_MPI_DETECT()
# ----------------
#
# Tries first to determine whether the MPI implementation is usable,
# then takes appropriate actions.
#
AC_DEFUN([ABI_MPI_DETECT], [
  # Delegate the actual detection to Steredeg
  SD_MPI_DETECT
  abi_mpi_ok="${sd_mpi_ok}"

  if test "${sd_mpi_enable}" = "yes"; then

    if test "${sd_mpi_ok}" = "yes"; then

      # Force abi_mpi_enable to "yes", for clarity and to avoid having to
      # further test "auto"
      abi_mpi_enable="yes"

      # Propagate MPI I/O trigger
      AC_MSG_CHECKING([whether to build MPI I/O code])
      AC_MSG_RESULT([${abi_mpi_io_enable}])
      if test "${abi_mpi_io_enable}" = "yes" -o "${abi_mpi_io_enable}" = "auto"; then
        abi_mpi_io_enable="yes"
        AC_DEFINE([HAVE_MPI_IO], 1,
          [Define to 1 if you want MPI I/O support.])
      fi

      # Check MPI I/O trigger
      if test "${abi_mpi_enable}" = "yes" -a "${abi_mpi_io_enable}" = "no"; then
        AC_MSG_WARN([disabling MPI I/O is not recommended])
      fi

      # Check MPI level actually supported
      _ABI_MPI_CHECK_FC_LEVEL

      # Select MPI level
      if test "${abi_mpi_level}" = ""; then
        abi_mpi_level="${abi_mpi_fc_level}"
      else
        AC_MSG_NOTICE([forcing MPI-${abi_mpi_level} standard level support])
        if test "${abi_mpi_level}" != "${abi_mpi_fc_level}"; then
        AC_MSG_WARN([detected MPI-${abi_mpi_fc_level} support but using MPI-${abi_mpi_level} instructions])
        fi
      fi

      # Propagate MPI level
      case "${abi_mpi_level}" in
        1)
          AC_MSG_ERROR([prehistoric MPI implementations are not supported anymore])
          ;;
        2)
          AC_DEFINE([HAVE_MPI2], 1,
            [Define to 1 if you have a MPI-2 implementation.])
          ;;
        3)
          AC_DEFINE([HAVE_MPI3], 1,
            [Define to 1 if you have a MPI-3 implementation.])
          ;;
        *)
          AC_MSG_ERROR([invalid MPI level: MPI-${abi_mpi_level}])
          ;;
      esac

      # Test the availability of problematic MPI constants
      _ABI_MPI_CHECK_INTEGER16

      # Check the availability of useful MPI primitives
      _ABI_MPI_CHECK_GET_LIBRARY_VERSION
      if test "${abi_mpi_get_library_version}" = "yes"; then
        AC_DEFINE([HAVE_MPI_GET_LIBRARY_VERSION], 1,
          [Define to 1 if your MPI implementation provides MPI_Get_library_version.])
      fi

      # Check the availability of problematic MPI primitives
      _ABI_MPI_CHECK_CREATE_TYPE_STRUCT()

      # Check MPI3 extensions (very) important for HPC.
      _ABI_MPI_CHECK_IBCAST()
      _ABI_MPI_CHECK_IALLTOALL()
      _ABI_MPI_CHECK_IALLTOALLV()
      _ABI_MPI_CHECK_IGATHERV()
      _ABI_MPI_CHECK_IALLREDUCE()

    fi # sd_mpi_ok

  fi # sd_mpi_enable
]) # ABI_MPI_DETECT


                    # ------------------------------------ #


# ABI_MPI_DUMP()
# --------------
#
# Dumps all MPI parameters, which is useful to diagnose faulty MPI
# environments.
#
AC_DEFUN([ABI_MPI_DUMP], [
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([dumping all MPI parameters for diagnostics])
  AC_MSG_NOTICE([------------------------------------------])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([Configure options:])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([  * enable_mpi_inplace = '${enable_mpi_inplace}'])
  AC_MSG_NOTICE([  * enable_mpi_io      = '${enable_mpi_io}'])
  AC_MSG_NOTICE([  * with_mpi           = '${with_mpi}'])
  AC_MSG_NOTICE([  * with_mpi_level     = '${with_mpi_level}'])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([Internal parameters])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([  * MPI enabled (required)                       : ${abi_mpi_enable}])
  AC_MSG_NOTICE([  * MPI C compiler is set (required)             : ${sd_mpi_cc_set}])
  AC_MSG_NOTICE([  * MPI C compiler works (required)              : ${sd_mpi_cc_ok}])
  AC_MSG_NOTICE([  * MPI Fortran compiler is set (required)       : ${sd_mpi_fc_set}])
  AC_MSG_NOTICE([  * MPI Fortran compiler works (required)        : ${sd_mpi_fc_ok}])
  AC_MSG_NOTICE([  * MPI environment usable (required)            : ${abi_mpi_ok}])
  AC_MSG_NOTICE([  * MPI C++ compiler is set (optional)           : ${sd_mpi_cxx_set}])
  AC_MSG_NOTICE([  * MPI C++ compiler works (optional)            : ${sd_mpi_cxx_ok}])
  AC_MSG_NOTICE([  * MPI-in-place enabled (optional)              : ${abi_mpi_inplace_enable}])
  AC_MSG_NOTICE([  * MPI-IO enabled (optional)                    : ${abi_mpi_io_enable}])
  AC_MSG_NOTICE([  * MPI configuration type (computed)            : ${abi_mpi_init}])
  AC_MSG_NOTICE([  * MPI Fortran level supported (detected)       : ${abi_mpi_fc_level}])
  AC_MSG_NOTICE([  * MPI_Get_library_version available (detected) : ${abi_mpi_get_library_version}])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([All required parameters must be set to 'yes'.])
  AC_MSG_NOTICE([If not, the configuration and/or the build with])
  AC_MSG_NOTICE([MPI support will very likely fail.])
  AC_MSG_NOTICE([])
]) # ABI_MPI_DUMP

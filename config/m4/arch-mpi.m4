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
  dnl Set default values
  abi_mpi_cc_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${CC_LDFLAGS}"
  LIBS="${CC_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a C MPI program
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

  dnl Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_CC


                    # ------------------------------------ #


# _ABI_MPI_CHECK_CXX()
# --------------------
#
# Checks whether the C++ compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_CXX], [
  dnl Set default values
  abi_mpi_cxx_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${CXX_LDFLAGS}"
  LIBS="${CXX_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a C++ MPI program
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

  dnl Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_CXX


                    # ------------------------------------ #


# _ABI_MPI_CHECK_FC()
# -------------------
#
# Checks whether the Fortran compiler is able to produce MPI binaries.
#
AC_DEFUN([_ABI_MPI_CHECK_FC], [
  dnl Set default values
  abi_mpi_fc_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a Fortran MPI program
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

  dnl Restore build environment
  ABI_ENV_RESTORE
]) # _ABI_MPI_CHECK_FC


                    # ------------------------------------ #


# _ABI_MPI_CHECK_FC_LEVEL()
# -------------------------
#
# Checks which MPI level is supported by the Fortran compiler.
#
AC_DEFUN([_ABI_MPI_CHECK_FC_LEVEL], [
  if test "${abi_mpi_fc_ok}" = "yes"; then

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    dnl Try to compile a MPI-2 Fortran program
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

    dnl Restore build environment
    ABI_ENV_RESTORE

  fi dnl abi_mpi_fc_ok = yes
]) # _ABI_MPI_CHECK_FC_LEVEL


                    # ------------------------------------ #


# _ABI_MPI_CHECK_INTEGER16()
# --------------------------
#
# Checks whether the MPI library supports MPI_INTEGER16
#
AC_DEFUN([_ABI_MPI_CHECK_INTEGER16], [
  dnl Set default values
  abi_mpi_integer16_ok="no"


  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${lib_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${lib_mpi_libs}"

  dnl Try to compile a Fortran program
  dnl Note: we assume a MPI implementation that provides the mpi module
  AC_MSG_CHECKING([whether the MPI library supports MPI_INTEGER16])
  AC_LANG_PUSH([Fortran])
  AC_LINK_IFELSE([AC_LANG_PROGRAM([],
    [[  
      use mpi

      integer, parameter :: ii = MPI_INTEGER16

    ]])], [abi_mpi_integer16_ok="yes"], [abi_mpi_integer16_ok="no"])
  AC_LANG_POP([Fortran])
  AC_MSG_RESULT([${abi_mpi_integer16_ok}])

  dnl Restore build environment
  ABI_ENV_RESTORE

  dnl Forward information to the compiler
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
  dnl Set default values
  abi_mpi_get_library_version="unknown"

  if test "${abi_mpi_fc_ok}" = "yes"; then

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
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

    dnl Restore build environment
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
  dnl Set default values
  abi_mpi_type_create_struct_ok="no"

  if test "${abi_mpi_fc_level}" -ge "2"; then

    dnl No problem should appear for MPI2 or MPI3 but we test it anyway.

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    dnl Try to compile a Fortran MPI program
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

    dnl Restore build environment
    ABI_ENV_RESTORE

  else

    dnl Back-up build environment
    ABI_ENV_BACKUP

    dnl Prepare build environment
    CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
    LDFLAGS="${FC_LDFLAGS}"
    LIBS="${FC_LIBS} ${abi_mpi_libs}"

    dnl Try to compile a Fortran MPI program
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

    dnl Restore build environment
    ABI_ENV_RESTORE

  fi

  dnl Forward information to the compiler
  if test "${abi_mpi_type_create_struct_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_TYPE_CREATE_STRUCT], 1,
      [Define to 1 if your MPI library supports MPI_TYPE_CREATE_STRUCT.])
  fi
]) # _ABI_MPI_CHECK_CREATE_TYPE_STRUCT


                    # ------------------------------------ #

# _ABI_MPI_CHECK_IBCAST()
# -----------------------------------
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
  dnl Set default values
  abi_mpi_ialltoall_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a Fortran program
  dnl Note: we assume a MPI implementation that provides the mpi module
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

  dnl Restore build environment
  ABI_ENV_RESTORE

  dnl Forward information to the compiler
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
  dnl Set default values
  abi_mpi_ialltoallv_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a Fortran program
  dnl Note: we assume a MPI implementation that provides the mpi module
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

  dnl Restore build environment
  ABI_ENV_RESTORE

  dnl Forward information to the compiler
  if test "${abi_mpi_ialltoallv_ok}" = "yes"; then
    AC_DEFINE([HAVE_MPI_IALLTOALLV], 1,
      [Define to 1 if your MPI library supports MPI_IALLTOALLV.])
  else
    AC_MSG_WARN([Your MPI library does not support non-blocking communications. The wall time of certain algorithms will increase with the number of MPI processes. It is strongly suggested to use a more recent MPI2+ library!])
  fi
]) # _ABI_MPI_CHECK_IALLTOALLV


                    # ------------------------------------ #


# _ABI_MPI_CHECK_IALLREDUCE()
# -----------------------------------
#
# Checks whether the MPI library supports MPI_IALLREDUCE (MPI3)
#
AC_DEFUN([_ABI_MPI_CHECK_IALLREDUCE], [
  dnl Set default values
  abi_mpi_iallreduce_ok="no"

  dnl Back-up build environment
  ABI_ENV_BACKUP

  dnl Prepare build environment
  CPPFLAGS="${CPPFLAGS} ${abi_mpi_incs}"
  LDFLAGS="${FC_LDFLAGS}"
  LIBS="${FC_LIBS} ${abi_mpi_libs}"

  dnl Try to compile a Fortran MPI program
  dnl Note: we assume a MPI implementation that provides the mpi module
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

  dnl Restore build environment
  ABI_ENV_RESTORE

  dnl Forward information to the compiler
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
  dnl Init
  tmp_name=`echo "$1" | sed -e 's/.*/\L&/'`
  ${MKDIR_P} config/wrappers

  dnl Create file
  cat >${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name} <<EOF
#!/bin/sh

$1="$2"
export $1

$3 \[$]{*}
EOF

  dnl Fix permissions
  chmod u+x ${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}

  dnl Overwrite compiler setting
  eval $1="${abinit_builddir}/config/wrappers/wrap-mpi${tmp_name}"

  dnl Clean-up the mess
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
  dnl Init
  abi_mpi_cc_ok="unknown"
  abi_mpi_cxx_ok="unknown"
  abi_mpi_complete="unknown"
  abi_mpi_fc_level="none"
  abi_mpi_fc_ok="unknown"
  abi_mpi_get_library_version="unknown"
  abi_mpi_has_cc="no"
  abi_mpi_has_cxx="no"
  abi_mpi_has_fc="no"
  abi_mpi_inplace_enable="${enable_mpi_inplace}"
  abi_mpi_io_enable="${enable_mpi_io}"
  abi_mpi_level=""
  abi_mpi_usable="no"
  MPI_CC=""
  MPI_CXX=""
  MPI_FC=""

  dnl Note: abi_mpi_enable and abi_mpi_flavor have been set
  dnl       in ABI_TRIGGERS_INIT
  AC_MSG_CHECKING([for MPI support requested])
  AC_MSG_RESULT([${abi_mpi_enable} (flavor: ${abi_mpi_flavor})])

  if test "${abi_mpi_enable}" = "yes" -o "${abi_mpi_enable}" = "auto"; then

    dnl Banner
    AC_MSG_NOTICE([Initializing MPI support])

    dnl Check whether to look for generic files
    if test "${abi_mpi_prefix}" != ""; then
      AC_MSG_NOTICE([looking for MPI in ${abi_mpi_prefix}])

      dnl Look for incompatibilities
      if test "${CPP}" != ""; then
        AC_MSG_WARN([${CPP} might not be fully compatible with MPI])
      fi
      if test "${with_mpi_incs}" != ""; then
        AC_MSG_ERROR([use --with-mpi or --with-mpi-includes, not both])
      fi
      if test "${abi_mpi_level}" != ""; then
        AC_MSG_WARN([forcing MPI level to ${abi_mpi_level} might make the build fail])
      fi
      if test "${with_mpi_libs}" != ""; then
        AC_MSG_ERROR([use --with-mpi or --with-mpi-libs, not both])
      fi

      dnl Look for a C compiler
      AC_MSG_CHECKING([for a MPI C compiler])
      if test -x "${abi_mpi_prefix}/bin/mpicc"; then
        abi_mpi_has_cc="yes"
        MPI_CC="${abi_mpi_prefix}/bin/mpicc"
      fi
      if test "${MPI_CC}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_CC}])
      fi

      dnl Look for a C++ compiler
      AC_MSG_CHECKING([for a MPI C++ compiler])
      if test -x "${abi_mpi_prefix}/bin/mpicxx"; then
        abi_mpi_has_cxx="yes"
        MPI_CXX="${abi_mpi_prefix}/bin/mpicxx"
      elif test -x "${abi_mpi_prefix}/bin/mpic++"; then
        abi_mpi_has_cxx="yes"
        MPI_CXX="${abi_mpi_prefix}/bin/mpic++"
      fi
      if test "${MPI_CXX}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_CXX}])
      fi

      dnl Look for a Fortran 90 compiler
      AC_MSG_CHECKING([for a MPI Fortran compiler])
      if test -x "${abi_mpi_prefix}/bin/mpif90"; then
        abi_mpi_has_fc="yes"
        MPI_FC="${abi_mpi_prefix}/bin/mpif90"
      fi
      if test "${MPI_FC}" = ""; then
        AC_MSG_RESULT([none found])
      else
        AC_MSG_RESULT([${MPI_FC}])
      fi

      dnl Report whether generic MPI implementation is sufficiently complete
      if test "${abi_mpi_has_cc}" = "yes" -a \
              "${abi_mpi_has_fc}" = "yes"; then
        abi_mpi_complete="yes"
        abi_mpi_flavor="prefix"

        dnl Fool-proof warning for those who set compilers twice
        dnl FIXME: does not handle the case when full path is not
        dnl        provided, e.g. prefix=/usr, CC=mpicc.
        dnl Hint: TRUE_CC=`echo "${CC}" | cut -d' ' -f1`
        dnl       TRUE_CC=`type -t "${TRUE_CC}"`
        dnl       if "file" type -p ...
        tmp_chk_cc="no"
        tmp_chk_redundant="no"
        if test "${CC}" = "${MPI_CC}"; then
          tmp_chk_cc="yes"
          tmp_chk_redundant="yes"
          AC_MSG_WARN([redundant setting of MPI C compiler!
    Use --with-mpi preferably.])
        fi
        tmp_chk_cxx="no"
        if test "${CXX}" = "${MPI_CXX}"; then
          tmp_chk_cxx="yes"
          AC_MSG_WARN([redundant setting of MPI C++ compiler!
    Use --with-mpi preferably.])
        fi
        if test "${tmp_chk_cxx}" != "${tmp_chk_redundant}"; then
          AC_MSG_WARN([inconsistent compiler settings!
    Use --with-mpi or set (CC, CXX, FC), not both.])
        fi
        tmp_chk_fc="no"
        if test "${FC}" = "${MPI_FC}"; then
          tmp_chk_fc="yes"
          AC_MSG_WARN([redundant setting of MPI Fortran compiler
    Use --with-mpi preferably.])
        fi
        if test "${tmp_chk_fc}" != "${tmp_chk_redundant}"; then
          AC_MSG_ERROR([inconsistent compiler settings!
    Use --with-mpi or set (CC, CXX, FC), not both.])
        fi
        if test "${tmp_chk_redundant}" = "yes"; then
          CC=""
          CXX=""
          FC=""
          AC_MSG_NOTICE([ignoring CC, CXX, and FC settings])
        fi
        unset tmp_chk_cc tmp_chk_cxx tmp_chk_fc tmp_chk_redundant

        dnl Decide whether to wrap MPI compiler calls
        dnl FIXME: flavor does not reflect mixed situations
        if test "${CC}" = ""; then
          CC="${MPI_CC}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI C compiler])
          _ABI_MPI_CREATE_WRAPPER([CC], [${CC}], [${MPI_CC}])
          abi_mpi_flavor="double-wrap"
        fi
        if test "${CXX}" = ""; then
          CXX="${MPI_CXX}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI C++ compiler])
          _ABI_MPI_CREATE_WRAPPER([CXX], [${CXX}], [${MPI_CXX}])
          abi_mpi_flavor="double-wrap"
        fi
        if test "${FC}" = ""; then
          FC="${MPI_FC}"
        else
          AC_MSG_NOTICE([creating wrapper for MPI Fortran compiler])
          _ABI_MPI_CREATE_WRAPPER([FC], [${FC}], [${MPI_FC}])
          abi_mpi_flavor="double-wrap"
        fi
      else
        unset MPI_CC
        unset MPI_CXX
        unset MPI_FC
        abi_mpi_complete="no"
      fi

    else dnl abi_mpi_prefix

      case "${abi_mpi_flavor}" in

        auto|native)
          dnl Look for MPI wrappers in PATH
          if test -z "${CC}" -a -z "${CXX}" -a -z "${FC}" -a \
                  -z "${MPI_CC}" -a -z "${MPI_CXX}" -a -z "${MPI_FC}"; then
            AC_CHECK_PROGS([MPI_CC], [mpicc mpiicc])
            AC_CHECK_PROGS([MPI_CXX], [mpic++ mpicxx])
            AC_CHECK_PROGS([MPI_FC], [mpifort mpif90 mpiifort])
            if test -n "${MPI_CC}" -a -n "${MPI_CXX}" -a -n "${MPI_FC}"; then
              CC="${MPI_CC}"
              CXX="${MPI_CXX}"
              FC="${MPI_FC}"
              abi_mpi_complete="yes"
              abi_mpi_flavor="native"
            fi
          else
            if test "${abi_mpi_flavor}" != "auto"; then
              AC_MSG_WARN([could not find all required MPI compilers])
            fi
            abi_mpi_complete="no"
          fi
          ;;

        double-wrap)
          if test -n "${CC}" -a -n "${CXX}" -a -n "${FC}" -a \
                  -n "${MPI_CC}" -a -n "${MPI_CXX}" -a -n "${MPI_FC}"; then
            _ABI_MPI_CREATE_WRAPPER([CC], [${CC}], [${MPI_CC}])
            _ABI_MPI_CREATE_WRAPPER([CXX], [${CXX}], [${MPI_CXX}])
            _ABI_MPI_CREATE_WRAPPER([FC], [${FC}], [${MPI_FC}])
            abi_mpi_complete="yes"
          else
            AC_MSG_WARN([some of the serial and/or MPI compilers are not defined])
            abi_mpi_complete="no"
          fi
          ;;

      esac dnl abi_mpi_flavor

    fi dnl abi_mpi_prefix

  else dnl abi_mpi_enable

    if test "${abi_mpi_init}" != "def"; then
      AC_MSG_NOTICE([MPI support disabled from command-line])
    fi
    abi_mpi_enable="no"
    abi_mpi_inplace_enable="no"
    abi_mpi_io_enable="no"
    abi_mpi_fcflags=""
    abi_mpi_ldflags=""
    abi_mpi_level=""
    abi_mpi_incs=""
    abi_mpi_libs=""
    abi_mpi_prefix=""
    abi_mpi_complete="no"

  fi dnl abi_mpi_enable

  if test "${abi_mpi_enable}" != "no"; then
    AC_MSG_CHECKING([whether MPI support is complete])
    AC_MSG_RESULT([${abi_mpi_complete}])
    AC_MSG_CHECKING([for MPI flavor after init])
    AC_MSG_RESULT([${abi_mpi_flavor}])
  fi

  dnl Enable substitution
  AC_SUBST(abi_mpi_enable)
  AC_SUBST(abi_mpi_fcflags)
  AC_SUBST(abi_mpi_ldflags)
  AC_SUBST(abi_mpi_level)
  AC_SUBST(abi_mpi_incs)
  AC_SUBST(abi_mpi_libs)
  AC_SUBST(abi_mpi_prefix)
]) # ABI_MPI_INIT


                    # ------------------------------------ #


# ABI_MPI_DETECT()
# ----------------
#
# Tries first to determine whether the MPI implementation is usable,
# then takes appropriate actions.
#
AC_DEFUN([ABI_MPI_DETECT], [
  dnl Init
  AC_REQUIRE([ABI_MPI_INIT])

  dnl Report status
  AC_MSG_CHECKING([whether to enable MPI support])
  AC_MSG_RESULT([${abi_mpi_enable} (flavor: ${abi_mpi_flavor})])

  if test "${abi_mpi_enable}" = "yes" -o "${abi_mpi_enable}" = "auto"; then

    dnl Check whether MPI is usable
    if test "${abi_mpi_complete}" = "yes"; then
      _ABI_MPI_CHECK_CC
      _ABI_MPI_CHECK_CXX
      _ABI_MPI_CHECK_FC

      if test "${abi_mpi_cc_ok}" = "yes"; then
        abi_mpi_has_cc="yes"
      fi
      if test "${abi_mpi_cxx_ok}" = "yes"; then
        abi_mpi_has_cxx="yes"
      fi
      if test "${abi_mpi_fc_ok}" = "yes"; then
        abi_mpi_has_fc="yes"
      fi
      if test "${abi_mpi_cc_ok}" = "yes" -a \
              "${abi_mpi_fc_ok}" = "yes"; then
        abi_mpi_usable="yes"
      fi
    fi
    AC_MSG_CHECKING([whether MPI is usable])
    AC_MSG_RESULT([${abi_mpi_usable}])

    if test "${abi_mpi_usable}" = "yes"; then

      dnl Force abi_mpi_enable to "yes", for clarity and to avoid having to
      dnl further test "auto"
      abi_mpi_enable="yes"

      dnl Propagate main trigger
      AC_DEFINE([HAVE_MPI], 1,
        [Define to 1 if you want to enable MPI support.])

      dnl Propagate MPI I/O trigger
      AC_MSG_CHECKING([whether to build MPI I/O code])
      AC_MSG_RESULT([${abi_mpi_io_enable}])
      if test "${abi_mpi_io_enable}" = "yes" -o "${abi_mpi_io_enable}" = "auto"; then
        AC_DEFINE([HAVE_MPI_IO], 1,
          [Define to 1 if you want MPI I/O support.])
        abi_mpi_io_enable="yes"
      fi

      dnl Check MPI I/O trigger
      if test "${abi_mpi_enable}" = "yes" -a "${abi_mpi_io_enable}" = "no"; then
        AC_MSG_WARN([disabling MPI I/O is not recommended])
      fi

      dnl Check MPI level actually supported
      _ABI_MPI_CHECK_FC_LEVEL

      dnl Select MPI level
      if test "${abi_mpi_level}" = ""; then
        abi_mpi_level="${abi_mpi_fc_level}"
      else
        AC_MSG_NOTICE([forcing MPI-${abi_mpi_level} standard level support])
        if test "${abi_mpi_level}" != "${abi_mpi_fc_level}"; then
        AC_MSG_WARN([detected MPI-${abi_mpi_fc_level} support but using MPI-${abi_mpi_level} instructions])
        fi
      fi

      dnl Propagate MPI level
      case "${abi_mpi_level}" in
        1)
          AC_DEFINE([HAVE_MPI1], 1,
            [Define to 1 if you have a MPI-1 implementation.])
          ;;
        2)
          AC_DEFINE([HAVE_MPI2], 1,
            [Define to 1 if you have a MPI-2 implementation.])
          ;;
        3)
          AC_DEFINE([HAVE_MPI3], 1,
            [Define to 1 if you have a MPI-3 implementation.])
          ;;
      esac

      dnl Test the availability of problematic MPI constants
      _ABI_MPI_CHECK_INTEGER16()

      dnl Check the availability of useful MPI primitives
      _ABI_MPI_CHECK_GET_LIBRARY_VERSION()
      if test "${abi_mpi_get_library_version}" = "yes"; then
        AC_DEFINE([HAVE_MPI_GET_LIBRARY_VERSION], 1,
          [Define to 1 if your MPI implementation provides MPI_Get_library_version.])
      fi

      dnl Check the availability of problematic MPI primitives
      _ABI_MPI_CHECK_CREATE_TYPE_STRUCT()

      dnl Check MPI3 extensions (very) important for HPC.
      _ABI_MPI_CHECK_IBCAST()
      _ABI_MPI_CHECK_IALLTOALL()
      _ABI_MPI_CHECK_IALLTOALLV()
      _ABI_MPI_CHECK_IALLREDUCE()

    fi dnl mpi_usable

  fi dnl enable_mpi

  AM_CONDITIONAL(DO_TEST_MPI, [test "${abi_mpi_enable}" = "yes"])
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
  AC_MSG_NOTICE([  * enable_mpi_inplace = ${enable_mpi_inplace}])
  AC_MSG_NOTICE([  * enable_mpi_io      = ${enable_mpi_io}])
  AC_MSG_NOTICE([  * with_mpi           = ${with_mpi}])
  AC_MSG_NOTICE([  * with_mpi_incs      = ${with_mpi_incs}])
  AC_MSG_NOTICE([  * with_mpi_level     = ${with_mpi_level}])
  AC_MSG_NOTICE([  * with_mpi_libs      = ${with_mpi_libs}])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([Internal parameters])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([  * MPI enabled (required)                       : ${abi_mpi_enable}])
  AC_MSG_NOTICE([  * MPI in-place enabled (optional)              : ${abi_mpi_inplace_enable}])
  AC_MSG_NOTICE([  * MPI-IO enabled (optional)                    : ${abi_mpi_io_enable}])
  AC_MSG_NOTICE([  * MPI C compiler present (required)            : ${abi_mpi_has_cc}])
  AC_MSG_NOTICE([  * MPI C compiler works (required)              : ${abi_mpi_cc_ok}])
  AC_MSG_NOTICE([  * MPI C++ compiler present (optional)          : ${abi_mpi_has_cxx}])
  AC_MSG_NOTICE([  * MPI C++ compiler works (optional)            : ${abi_mpi_cxx_ok}])
  AC_MSG_NOTICE([  * MPI Fortran compiler present (required)      : ${abi_mpi_has_fc}])
  AC_MSG_NOTICE([  * MPI Fortran compiler works (required)        : ${abi_mpi_fc_ok}])
  AC_MSG_NOTICE([  * MPI environment complete (required)          : ${abi_mpi_complete}])
  AC_MSG_NOTICE([  * MPI environment usable (required)            : ${abi_mpi_usable}])
  AC_MSG_NOTICE([  * MPI actual flavor (computed)                 : ${abi_mpi_flavor}])
  AC_MSG_NOTICE([  * MPI configuration type (computed)            : ${abi_mpi_init}])
  AC_MSG_NOTICE([  * MPI Fortran level supported (detected)       : ${abi_mpi_fc_level}])
  AC_MSG_NOTICE([  * MPI_Get_library_version available (detected) : ${abi_mpi_get_library_version}])
  AC_MSG_NOTICE([])
  AC_MSG_NOTICE([All required parameters must be set to 'yes'.])
  AC_MSG_NOTICE([If not, the configuration and/or the build will])
  AC_MSG_NOTICE([very likely fail.])
  AC_MSG_NOTICE([])

  if test "${abi_mpi_usable}" = "no"; then

    case "${abi_mpi_enable}" in
      auto)
        AC_MSG_WARN([MPI support is broken - disabling MPI])
        ;;
      yes)
        AC_MSG_ERROR([MPI support is broken - please fix your config parameters and/or MPI installation])
        ;;
    esac

    abi_mpi_enable="no"
    abi_mpi_inplace_enable="no"
    abi_mpi_io_enable="no"
    abi_mpi_fcflags=""
    abi_mpi_ldflags=""
    abi_mpi_level=""
    abi_mpi_incs=""
    abi_mpi_libs=""
    abi_mpi_prefix=""
  fi
]) # ABI_MPI_DUMP

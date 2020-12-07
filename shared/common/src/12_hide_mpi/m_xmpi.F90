!!****m* ABINIT/m_xmpi
!! NAME
!!  m_xmpi
!!
!! FUNCTION
!!  This module provides MPI named constants, tools for inquiring the MPI environment
!!  and a set of generic interfaces wrapping the most commonly used MPI primitives.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2020 ABINIT group (MG, MB, XG, YP, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! TODO
!!  Get rid of xmpi_paral. Sequential code is the **exception**. Developers should code parallel
!!  code or code that is compatible both with MPI and seq (thanks to the wrappers provided by this module)
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

MODULE m_xmpi

 use defs_basis
 use m_profiling_abi
 use iso_c_binding
#ifdef HAVE_FC_ISO_FORTRAN_2008
 use ISO_FORTRAN_ENV, only : int16,int32,int64
#endif
#ifdef HAVE_MPI2
 use mpi
#endif
#ifdef FC_NAG
 use f90_unix_proc
#endif

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif
#ifndef HAVE_FC_ISO_FORTRAN_2008
 integer,parameter :: int16=2,int32=4,int64=8
#endif

#ifdef HAVE_MPI
 ! MPI constants used in abinit. Make sure that a corresponding fake value is provided for the sequential version.
 integer,public,parameter :: xmpi_paral          = 1
 integer,public,parameter :: xmpi_world          = MPI_COMM_WORLD
 integer,public,parameter :: xmpi_comm_self      = MPI_COMM_SELF
 integer,public,parameter :: xmpi_undefined      = MPI_UNDEFINED
 integer,public,parameter :: xmpi_undefined_rank = MPI_UNDEFINED  ! MPI_UNDEFINED_RANK is not portable.
 integer,public,parameter :: xmpi_comm_null      = MPI_COMM_NULL
 integer,public,parameter :: xmpi_group_null     = MPI_GROUP_NULL
 integer,public,parameter :: xmpi_any_source     = MPI_ANY_SOURCE
 integer,public,parameter :: xmpi_request_null   = MPI_REQUEST_NULL
 integer,public,parameter :: xmpi_msg_len        = MPI_MAX_ERROR_STRING ! Length of fortran string used to store MPI error strings.
 integer,public,parameter :: xmpi_info_null      = MPI_INFO_NULL
 integer,public,parameter :: xmpi_success        = MPI_SUCCESS
#else
 ! Fake replacements for the sequential version. Values are taken from
 ! http://www.mit.edu/course/13/13.715/sun-hpc-ct-8.2.1/Linux/sun/include/mpif-common.h
 ! Please use these conventions when adding new replacements in order to avoid collisions between values.
 integer,public,parameter :: xmpi_paral          = 0
 integer,public,parameter :: xmpi_world          = 0
 integer,public,parameter :: xmpi_comm_self      = 1
 integer,public,parameter :: xmpi_undefined      =-32766
 integer,public,parameter :: xmpi_undefined_rank =-32766
 integer,public,parameter :: xmpi_comm_null      = 2
 integer,public,parameter :: xmpi_group_null     = 0
 integer,public,parameter :: xmpi_any_source     = -1
 integer,public,parameter :: xmpi_request_null   = 0
 integer,public,parameter :: xmpi_msg_len        = 1000
 integer,public,parameter :: xmpi_info_null      = 0
 integer,public,parameter :: xmpi_success        = 0
#endif

#ifdef HAVE_MPI
 integer,save,private  :: xmpi_tag_ub=32767
 ! The tag upper bound value must be at least 32767. An MPI implementation is free to make
 ! the value of MPI_TAG_UB larger than this hence xmpi_tag_ub is redefined when MPI is init in xmpi_init.
#endif

 ! Size in bytes of the entries used in MPI datatypes.
 integer,save, public ABI_PROTECTED:: xmpi_bsize_ch  = 0
 integer,save, public ABI_PROTECTED:: xmpi_bsize_int = 0
 integer,save, public ABI_PROTECTED:: xmpi_bsize_sp  = 0
 integer,save, public ABI_PROTECTED:: xmpi_bsize_dp  = 0
 integer,save, public ABI_PROTECTED:: xmpi_bsize_spc = 0
 integer,save, public ABI_PROTECTED:: xmpi_bsize_dpc = 0

 ! kind of the offset used for MPI-IO.
#ifdef HAVE_MPI_IO
 integer,public,parameter :: xmpi_offset_kind  = MPI_OFFSET_KIND
 integer,public,parameter :: xmpi_address_kind = MPI_ADDRESS_KIND
 integer,public,parameter :: xmpi_mpiio = 1
#else
 integer,public,parameter :: xmpi_offset_kind = i8b
 integer,public,parameter :: xmpi_address_kind = i8b
 integer,public,parameter :: xmpi_mpiio = 0
#endif

 ! The byte size and the MPI type of the Fortran record marker.
 ! These quantities are compiler-dependent and are initalized here
 ! for selected compilers or in xmpio_get_info_frm that is called by xmpi_init (only if MPI-IO is on).
#if defined HAVE_MPI && (defined FC_INTEL || defined FC_GNU || defined FC_IBM)
 integer,save,public ABI_PROTECTED :: xmpio_bsize_frm   = 4
 integer,save,public ABI_PROTECTED :: xmpio_mpi_type_frm= MPI_INTEGER4
#else
 integer,save,public ABI_PROTECTED :: xmpio_bsize_frm    = 0
 integer,save,public ABI_PROTECTED :: xmpio_mpi_type_frm = 0
#endif

 integer,save, public ABI_PROTECTED :: xmpio_info = xmpi_info_null
 ! Global variable used to pass hints to the MPI-IO routines.

 integer(XMPI_OFFSET_KIND),public,parameter :: xmpio_chunk_bsize = 2000 * (1024.0_dp**2)
 ! Defines the chunk size (in bytes) used to (read|write) data in a single MPI-IO call.
 ! MPI-IO, indeed, crashes if we try to do the IO of a large array with a single call.
 ! We use a value <= 2  Gb to avoid wraparound errors with standard integers.

 ! Options used for the MPI-IO wrappers used in abinit.
 integer,public,parameter :: xmpio_single     = 1  ! Individual IO.
 integer,public,parameter :: xmpio_collective = 2  ! Collective IO.

 integer,save, public ABI_PROTECTED :: xmpi_count_requests = 0
 ! Count number of requests (+1 for each call to non-blocking API, -1 for each call to xmpi_wait)
 ! This counter should be zero at the end of the run if all requests have been released)

 logical,save, private :: xmpi_use_inplace_operations = .False.
 ! Enable/disable usage of MPI_IN_PLACE in e.g. xmpi_sum

 ! For MPI <v4, collective communication routines accept only a 32bit integer as data count.
 ! To exchange more than 2^32 data we need to create specific user-defined datatypes
 ! For this, we need some parameters:
 integer(KIND=int32),public,parameter :: xmpi_maxint32 = huge(0_int32)
 integer(KIND=int64),public,parameter :: xmpi_maxint32_64 = int(xmpi_maxint32,kind=int64)
 ! Max. integer that can be represented with 32 bits
 integer(KIND=int64),public,save :: xmpi_largetype_size = 0
 ! Number of data to be used in user-defined operations related to user-defined "largetype" type
!!***

!----------------------------------------------------------------------

!!****t* m_xmpi/xcomm_t
!! NAME
!! xcomm_t
!!
!! FUNCTION
!!  A small object storing the MPI communicator, the rank of the processe and the size of the communicator.
!!  Provides helper functions to perform typical operations and parallelize loops.
!!  The datatype is initialized with xmpi_comm_self
!!
!! SOURCE

 type, public :: xcomm_t
   integer :: value = xmpi_comm_self
   integer :: nproc = 1
   integer :: me = 0
 contains
   ! procedure :: iam_master => xcomm_iam_master
   procedure :: skip => xcomm_skip
   procedure :: set_to_null => xcomm_set_to_null
   procedure :: set_to_self => xcomm_set_to_self
   procedure :: free => xcomm_free
 end type xcomm_t
!!***

! Public procedures.
 public :: xmpi_init                  ! Initialize the MPI environment.
 public :: xmpi_set_inplace_operations! Set internal flag to use MPI_IN_PLACE whenever possible.
 public :: xmpi_end                   ! Terminate the MPI environment.
 public :: xmpi_abort                 ! Hides MPI_ABORT from MPI library.
 public :: xmpi_show_info             ! Printout of the basic variables stored in this module (useful for debugging).
 public :: xmpi_group_free            ! Hides MPI_GROUP_FREE from MPI library.
 public :: xmpi_group_incl            ! Hides MPI_GROUP_INCL from MPI library.
 public :: xmpi_group_translate_ranks ! Hides MPI_GROUP_TRANSLATE_RANKS from MPI library.
 public :: xmpi_comm_create           ! Hides MPI_COMM_CREATE from MPI library.
 public :: xmpi_comm_rank             ! Hides MPI_COMM_RANK from MPI library.
 public :: xmpi_comm_size             ! Hides MPI_COMM_SIZE from MPI library.
 public :: xmpi_comm_free             ! Hides MPI_COMM_FREE from MPI library.
 public :: xmpi_comm_group            ! Hides MPI_COMM_GROUP from MPI library.
 public :: xmpi_comm_translate_ranks  ! Hides MPI_GROUP_TRANSLATE_RANKS from MPI library.
 public :: xmpi_comm_split            ! Hides MPI_COMM_SPLIT from MPI library.
 public :: xmpi_subcomm               ! Creates a sub-communicator from an input communicator.
 public :: xmpi_barrier               ! Hides MPI_BARRIER from MPI library.
 public :: xmpi_name                  ! Hides MPI_NAME from MPI library.
 public :: xmpi_iprobe                ! Hides MPI_IPROBE from MPI library.
 public :: xmpi_wait                  ! Hides MPI_WAIT from MPI library.
 public :: xmpi_waitall               ! Hides MPI_WAITALL from MPI library.
 public :: xmpi_request_free          ! Hides MPI_REQUEST_FREE from MPI library.
 public :: xmpi_comm_set_errhandler   ! Hides MPI_COMM_SET_ERRHANDLER from MPI library.
 public :: xmpi_error_string          ! Return a string describing the error from ierr.
 public :: xmpi_split_work            ! Splits tasks inside communicator using blocks
 public :: xmpi_split_block           ! Splits tasks inside communicator using block distribution.
 public :: xmpi_split_cyclic          ! Splits tasks inside communicator using cyclic distribution.
 public :: xmpi_split_list            ! Splits list of indices inside communicator using block distribution.
 public :: xmpi_distab                ! Fill table defining the distribution of the tasks according to the # of processors
 public :: xmpi_distrib_with_replicas ! Distribute tasks among MPI ranks (replicas are allowed)

! Private procedures.
 private :: xmpi_largetype_create      ! Build a large-count contiguous datatype (to handle a very large # of data)
 private :: xmpi_largetype_free        ! Release a large-count contiguous datatype

 interface xmpi_comm_free
   module procedure xmpi_comm_free_0D
   module procedure xmpi_comm_free_1D
   module procedure xmpi_comm_free_2D
   module procedure xmpi_comm_free_3D
 end interface xmpi_comm_free

 interface xmpi_waitall
   module procedure xmpi_waitall_1d
   module procedure xmpi_waitall_2d
 end interface xmpi_waitall

 interface xmpi_split_work
   module procedure xmpi_split_work_i4b
 end interface xmpi_split_work

 public :: xmpi_split_work2_i4b
 public :: xmpi_split_work2_i8b
 !public :: xmpi_split_work2
 !
 ! g95@green v0.93 is not able to resolve the interface.
 ! For the time being, this generic interface has been disabled.
 !interface xmpi_split_work2
 !  module procedure xmpi_split_work2_i4b
 !  module procedure xmpi_split_work2_i8b
 !end interface xmpi_split_work2

 interface xmpi_distab
   module procedure xmpi_distab_4D
 end interface xmpi_distab

 ! MPI generic interfaces.
 public :: xmpi_allgather
 public :: xmpi_iallgather
 public :: xmpi_allgatherv
 public :: xmpi_alltoall
 public :: xmpi_ialltoall
 public :: xmpi_alltoallv
 public :: xmpi_ialltoallv
 public :: xmpi_bcast
 public :: xmpi_ibcast
 public :: xmpi_exch
 public :: xmpi_gather
 public :: xmpi_gatherv
 public :: xmpi_max
 public :: xmpi_min
 public :: xmpi_recv
 public :: xmpi_irecv
 public :: xmpi_scatterv
 public :: xmpi_send
 public :: xmpi_isend
 public :: xmpi_sum_master
 public :: xmpi_sum
 public :: xmpi_isum
 public :: xmpi_isum_ip
 public :: xmpi_land              ! allreduce with MPI_LAND
 public :: xmpi_lor               ! allreduce with MPI_LOR

#ifdef HAVE_MPI_IO
 public :: xmpio_max_address      !  Returns .TRUE. if offset cannot be stored in integer(kind=XMPI_ADDRESS_KIND).
 public :: xmpio_type_struct
 public :: xmpio_get_info_frm
 public :: xmpio_check_frmarkers
 public :: xmpio_read_frm
 public :: xmpio_read_int
 public :: xmpio_read_dp
 public :: xmpio_write_frm
 public :: xmpio_write_frmarkers

 public :: xmpio_create_fstripes
 public :: xmpio_create_fsubarray_2D
 public :: xmpio_create_fsubarray_3D
 public :: xmpio_create_fsubarray_4D
 public :: xmpio_create_fherm_packed
 public :: xmpio_create_coldistr_from_fpacked
 public :: xmpio_create_coldistr_from_fp3blocks

!interface xmpio_read
!  module procedure xmpio_read_int
!  module procedure xmpio_read_dp
!end interface xmpio_read
!
!interface xmpio_write
!  module procedure xmpio_write_int
!  module procedure xmpio_write_dp
!end interface xmpio_write
#endif

!----------------------------------------------------------------------

interface xmpi_allgather
  module procedure xmpi_allgather_int
  module procedure xmpi_allgather_char
  module procedure xmpi_allgather_int1d_1b
  module procedure xmpi_allgather_int1d
  module procedure xmpi_allgather_int2d
  module procedure xmpi_allgather_dp1d
  module procedure xmpi_allgather_dp2d
  module procedure xmpi_allgather_dp3d
  module procedure xmpi_allgather_dp4d
end interface xmpi_allgather

!----------------------------------------------------------------------

! non-blocking version (requires MPI3)
! Prototype:
!
!   call xmpi_iallgather(SENDBUF, SENDCOUNT, SENDTYPE, RECVBUF, RECVCOUNT, RECVTYPE, COMM, REQUEST, IERROR)
!
! If the MPI library does not provide ialltoall, we call the blocking version and
! we return xmpi_request_null (see xmpi_iallgather.finc)
! Client code should always test/wait the request so that code semantics is preserved.

interface xmpi_iallgather
  module procedure xmpi_iallgather_dp4d
end interface xmpi_iallgather

interface xmpi_allgatherv
  module procedure xmpi_allgatherv_int2d
  module procedure xmpi_allgatherv_int
  module procedure xmpi_allgatherv_int1_dp1
  module procedure xmpi_allgatherv_dp
  module procedure xmpi_allgatherv_dp2d
  module procedure xmpi_allgatherv_dp3d
  module procedure xmpi_allgatherv_dp4d
  module procedure xmpi_allgatherv_dp5d
  module procedure xmpi_allgatherv_dp6d
  module procedure xmpi_allgatherv_coeff2d
  module procedure xmpi_allgatherv_coeff2d_indx
end interface xmpi_allgatherv

!----------------------------------------------------------------------

! blocking
interface xmpi_alltoall
  module procedure xmpi_alltoall_int
  module procedure xmpi_alltoall_dp2d
  module procedure xmpi_alltoall_dp4d
end interface xmpi_alltoall

! non-blocking version (requires MPI3)
! Prototype:
!
!   call xmpi_ialltoall(xval, sendsize, recvbuf, recvsize, comm, request)
!
! If the MPI library does not provide ialltoall, we call the blocking version and
! we return xmpi_request_null (see xmpi_ialltoall.finc)
! Client code should always test/wait the request so that code semantics is preserved.

interface xmpi_ialltoall
  module procedure xmpi_ialltoall_dp4d
end interface xmpi_ialltoall

!----------------------------------------------------------------------

interface xmpi_alltoallv
  module procedure xmpi_alltoallv_dp2d
  module procedure xmpi_alltoallv_int2d
  module procedure xmpi_alltoallv_dp1d
  module procedure xmpi_alltoallv_dp1d2
end interface xmpi_alltoallv

!----------------------------------------------------------------------

! non-blocking version (requires MPI3)
! Prototype:
!
!   call xmpi_ialltoallv(xval, sendcnts, sdispls, recvbuf, recvcnts, rdispls, comm, request)
!
! If the MPI library does not provide ialltoallv, we call the blocking version and
! we return xmpi_request_null (see xmpi_ialltoallv.finc)
! Client code should always test/wait the request so that code semantics is preserved.

interface xmpi_ialltoallv
  module procedure xmpi_ialltoallv_dp2d
  module procedure xmpi_ialltoallv_int2d
  module procedure xmpi_ialltoallv_dp1d2
end interface xmpi_ialltoallv

!----------------------------------------------------------------------

interface xmpi_bcast
  module procedure xmpi_bcast_intv
  module procedure xmpi_bcast_int1d
  module procedure xmpi_bcast_int2d
  module procedure xmpi_bcast_int3d
  module procedure xmpi_bcast_int4d
  module procedure xmpi_bcast_dpv
  module procedure xmpi_bcast_dp1d
  module procedure xmpi_bcast_dp2d
  module procedure xmpi_bcast_dp3d
  module procedure xmpi_bcast_dp4d
  module procedure xmpi_bcast_dp5d
  module procedure xmpi_bcast_dp6d
  module procedure xmpi_bcast_spv
  module procedure xmpi_bcast_sp1d
  module procedure xmpi_bcast_sp2d
  module procedure xmpi_bcast_sp3d
  module procedure xmpi_bcast_sp4d
  module procedure xmpi_bcast_cplxv
  module procedure xmpi_bcast_cplx1d
  module procedure xmpi_bcast_cplx2d
  module procedure xmpi_bcast_cplx3d
  module procedure xmpi_bcast_cplx4d
  module procedure xmpi_bcast_dcv
  module procedure xmpi_bcast_dc1d
  module procedure xmpi_bcast_dc2d
  module procedure xmpi_bcast_dc3d
  module procedure xmpi_bcast_dc4d
  module procedure xmpi_bcast_ch0d
  module procedure xmpi_bcast_ch1d
  module procedure xmpi_bcast_log0d
  module procedure xmpi_bcast_coeffi2_1d
  module procedure xmpi_bcast_coeff2_1d
end interface xmpi_bcast

!----------------------------------------------------------------------

interface xmpi_ibcast
  module procedure xmpi_ibcast_int1d
  module procedure xmpi_ibcast_int4d
  module procedure xmpi_ibcast_dp1d
  module procedure xmpi_ibcast_dp2d
  module procedure xmpi_ibcast_dp3d
  module procedure xmpi_ibcast_dp4d
end interface xmpi_ibcast

!----------------------------------------------------------------------

interface xmpi_exch
  module procedure xmpi_exch_intn
  module procedure xmpi_exch_int2d
  module procedure xmpi_exch_dpn
  module procedure xmpi_exch_dp2d
  module procedure xmpi_exch_dp3d
  module procedure xmpi_exch_dp4d_tag
  module procedure xmpi_exch_dp5d_tag
  module procedure xmpi_exch_spc_1d
  module procedure xmpi_exch_dpc_1d
  module procedure xmpi_exch_dpc_2d
end interface xmpi_exch

!----------------------------------------------------------------------

interface xmpi_gather
  module procedure xmpi_gather_int
  module procedure xmpi_gather_int2d
  module procedure xmpi_gather_dp
  module procedure xmpi_gather_dp2d
  module procedure xmpi_gather_dp3d
  module procedure xmpi_gather_dp4d
end interface xmpi_gather

!----------------------------------------------------------------------

interface xmpi_gatherv
  module procedure xmpi_gatherv_int
  module procedure xmpi_gatherv_int1_dp1
  module procedure xmpi_gatherv_int2d
  module procedure xmpi_gatherv_dp
  module procedure xmpi_gatherv_dp2d
  module procedure xmpi_gatherv_dp3d
  module procedure xmpi_gatherv_dp4d
  module procedure xmpi_gatherv_dp5d
  module procedure xmpi_gatherv_dp6d
end interface xmpi_gatherv

!----------------------------------------------------------------------

interface xmpi_max
  module procedure xmpi_max_int0d_i4b
  module procedure xmpi_max_int0d_i8b
  module procedure xmpi_max_int
  module procedure xmpi_max_dpv
  module procedure xmpi_max_dp0d_ip
end interface xmpi_max

!----------------------------------------------------------------------

interface xmpi_min
  module procedure xmpi_min_intv
  module procedure xmpi_min_dpv
end interface xmpi_min

!----------------------------------------------------------------------

!interface xmpi_min_max
!  module procedure xmpi_min_max_int0d_i4b
!end interface xmpi_min_max

!----------------------------------------------------------------------

interface xmpi_recv
  module procedure xmpi_recv_char
  module procedure xmpi_recv_intv
  module procedure xmpi_recv_int1d
  module procedure xmpi_recv_int2d
  module procedure xmpi_recv_int3d
  module procedure xmpi_recv_dp
  module procedure xmpi_recv_dp1d
  module procedure xmpi_recv_dp2d
  module procedure xmpi_recv_dp3d
end interface xmpi_recv

!----------------------------------------------------------------------

interface xmpi_irecv
  module procedure xmpi_irecv_intv
  module procedure xmpi_irecv_int1d
  module procedure xmpi_irecv_dp1d
  module procedure xmpi_irecv_dp2d
end interface xmpi_irecv

!----------------------------------------------------------------------

interface xmpi_scatterv
  module procedure xmpi_scatterv_int
  module procedure xmpi_scatterv_int2d
  module procedure xmpi_scatterv_dp
  module procedure xmpi_scatterv_dp2d
  module procedure xmpi_scatterv_dp3d
  module procedure xmpi_scatterv_dp4d
end interface xmpi_scatterv

!----------------------------------------------------------------------

interface xmpi_isend
  module procedure xmpi_isend_int1d
  module procedure xmpi_isend_dp1d
  module procedure xmpi_isend_dp2d
end interface xmpi_isend

!----------------------------------------------------------------------

interface xmpi_send
  module procedure xmpi_send_char
  module procedure xmpi_send_intv
  module procedure xmpi_send_int1d
  module procedure xmpi_send_int2d
  module procedure xmpi_send_int3d
  module procedure xmpi_send_dp
  module procedure xmpi_send_dp1d
  module procedure xmpi_send_dp2d
  module procedure xmpi_send_dp3d
end interface xmpi_send

!----------------------------------------------------------------------

interface xmpi_sum_master
  module procedure xmpi_sum_master_int
  module procedure xmpi_sum_master_int2d
  module procedure xmpi_sum_master_int4d
  module procedure xmpi_sum_master_dp
  module procedure xmpi_sum_master_dp1d
  module procedure xmpi_sum_master_dp2d
  module procedure xmpi_sum_master_dp3d
  module procedure xmpi_sum_master_dp4d
  module procedure xmpi_sum_master_dp5d
  module procedure xmpi_sum_master_dp6d
  module procedure xmpi_sum_master_dp7d
  module procedure xmpi_sum_master_c1cplx
  module procedure xmpi_sum_master_c2cplx
  module procedure xmpi_sum_master_c3cplx
  module procedure xmpi_sum_master_c4cplx
  module procedure xmpi_sum_master_c5cplx
  module procedure xmpi_sum_master_c1dpc
  module procedure xmpi_sum_master_c2dpc
  module procedure xmpi_sum_master_c3dpc
  module procedure xmpi_sum_master_c4dpc
  module procedure xmpi_sum_master_c5dpc
end interface xmpi_sum_master

!----------------------------------------------------------------------

!MG:TODO procedure marked with !? are considered obsolete.
!   and will be removed in future versions.
!   Please use interfaces where array dimensions are not passed explicitly.
!   Rationale: The array descriptor is already passed to the routine
!   so it does not make sense to pass the dimension explicitly.

interface xmpi_sum
  module procedure xmpi_sum_int
  module procedure xmpi_sum_intv
  module procedure xmpi_sum_intv2
  module procedure xmpi_sum_intn   !?
  module procedure xmpi_sum_int2t  !?
  module procedure xmpi_sum_int2d
  module procedure xmpi_sum_int3d
  module procedure xmpi_sum_int4d
  module procedure xmpi_sum_dp
  module procedure xmpi_sum_dpvt
  module procedure xmpi_sum_dpv
  module procedure xmpi_sum_dpn    !?
  module procedure xmpi_sum_sp2d
  module procedure xmpi_sum_sp3d
  module procedure xmpi_sum_sp4d
  module procedure xmpi_sum_sp5d
  module procedure xmpi_sum_sp6d
  module procedure xmpi_sum_sp7d
  module procedure xmpi_sum_dp2d
  module procedure xmpi_sum_dp3d
  module procedure xmpi_sum_dp4d
  module procedure xmpi_sum_dp5d
  module procedure xmpi_sum_dp6d
  module procedure xmpi_sum_dp7d
  module procedure xmpi_sum_dp2t   !?
  module procedure xmpi_sum_dp2d2t
  module procedure xmpi_sum_dp3d2t !?
  module procedure xmpi_sum_dp4d2t !?
  module procedure xmpi_sum_c0dc
  module procedure xmpi_sum_c1dc
  module procedure xmpi_sum_c2dc
  module procedure xmpi_sum_c3dc
  module procedure xmpi_sum_c4dc
  module procedure xmpi_sum_c5dc
  module procedure xmpi_sum_c6dc
  module procedure xmpi_sum_c7dc
  module procedure xmpi_sum_c1cplx
  module procedure xmpi_sum_c2cplx
  module procedure xmpi_sum_c3cplx
  module procedure xmpi_sum_c4cplx
  module procedure xmpi_sum_c5cplx
  module procedure xmpi_sum_c6cplx
end interface xmpi_sum
!!***

! Non-blocking version
interface xmpi_isum
  module procedure xmpi_isum_int0d
end interface xmpi_isum
!!***

! Non-blocking in-place version
interface xmpi_isum_ip
  module procedure xmpi_isum_ip_dp2d
  module procedure xmpi_isum_ip_dp4d
end interface xmpi_isum_ip
!!***

interface xmpi_land
  module procedure xmpi_land_log0d
end interface xmpi_land
!!***

interface xmpi_lor
  module procedure xmpi_lor_log1d
  module procedure xmpi_lor_log2d
  module procedure xmpi_lor_log3d
end interface xmpi_lor
!!!***


!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_xmpi/xmpi_init
!! NAME
!!  xmpi_init
!!
!! FUNCTION
!!  Hides MPI_INIT from MPI library. Perform the initialization of some basic variables
!!  used by the MPI routines employed in abinit.
!!
!! INPUTS
!!  None
!!
!! PARENTS
!!      abinit,abitk,aim,anaddb,atdep,band2eps,conducti,cut3d,dummy_tests
!!      fftprof,fold2Bloch,ioprof,lapackprof,macroave,mrgddb,mrgdv,mrggkk
!!      mrgscr,multibinit,optic,testtransposer,ujdet,vdw_kernelgen
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_init()

!Local variables-------------------
 integer :: mpierr,ierr,unt
 logical :: exists
#ifdef HAVE_MPI
 integer :: attribute_val
 logical :: lflag
#ifdef HAVE_OPENMP
 integer :: required,provided
#endif
#endif

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI

#ifndef HAVE_OPENMP
 call MPI_INIT(mpierr)
#else
 required = MPI_THREAD_SINGLE
 !required = MPI_THREAD_FUNNELED
 !required = MPI_THREAD_SERIALIZED
 !required = MPI_THREAD_MULTIPLE
 call MPI_INIT_THREAD(required,provided,mpierr)
 if (provided /= required) then
   call xmpi_abort(msg="MPI_INIT_THREADS: provided /= required")
 end if
#endif

 !%comm_world = xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
 !%call xmpi_comm_set_errhandler(comm_world, MPI_ERRORS_RETURN, err_handler_sav, mpierr)

 ! Deprecated in MPI2 but not all MPI2 implementations provide MPI_Comm_get_attr !
 call MPI_ATTR_GET(xmpi_world, MPI_TAG_UB, attribute_val, lflag, mpierr)
 !call MPI_Comm_get_attr(xmpi_world, MPI_TAG_UB, attribute_val, lflag, mpierr)

 if (lflag) xmpi_tag_ub = attribute_val

!  Define type values.
 call MPI_TYPE_SIZE(MPI_CHARACTER,xmpi_bsize_ch,mpierr)
 call MPI_TYPE_SIZE(MPI_INTEGER,xmpi_bsize_int,mpierr)
 call MPI_TYPE_SIZE(MPI_REAL,xmpi_bsize_sp,mpierr)
 call MPI_TYPE_SIZE(MPI_DOUBLE_PRECISION,xmpi_bsize_dp,mpierr)
 call MPI_TYPE_SIZE(MPI_COMPLEX,xmpi_bsize_spc,mpierr)
 call MPI_TYPE_SIZE(MPI_DOUBLE_COMPLEX,xmpi_bsize_dpc,mpierr)

 ! Find the byte size of Fortran record marker used in MPI-IO routines.
 if (xmpio_bsize_frm == 0) then
   call xmpio_get_info_frm(xmpio_bsize_frm, xmpio_mpi_type_frm, xmpi_world)
 end if
#endif

 ! Master Removes the ABI_MPIABORTFILE if present so that we start with a clean environment
 if (xmpi_comm_rank(xmpi_world) == 0) then
    inquire(file=ABI_MPIABORTFILE, exist=exists)
    if (exists) then
       ! Get free unit (emulate F2008 newunit for portability reasons)
       unt = xmpi_get_unit()
       if (unt == -1) call xmpi_abort(msg="Cannot find free unit!!")
       open(unit=unt, file=trim(ABI_MPIABORTFILE), status="old", iostat=ierr)
       if (ierr == 0) close(unit=unt, status="delete", iostat=ierr)
       if (ierr /= 0) call xmpi_abort(msg="Cannot remove ABI_MPIABORTFILE")
    end if
 end if

end subroutine xmpi_init
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_set_inplace_operations
!! NAME
!!  xmpi_set_inplace_operations
!!
!! FUNCTION
!!  Set internal flag to use MPI_IN_PLACE whenever possible.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine xmpi_set_inplace_operations(bool)

!Local variables-------------------
 logical :: bool

! *************************************************************************

 xmpi_use_inplace_operations = bool

end subroutine xmpi_set_inplace_operations
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_get_unit
!! NAME
!!  xmpi_get_unit
!!
!! FUNCTION
!! Get free unit (emulate F2008 newunit for portability reasons)
!! Return -1 if no unit is found.
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

integer function xmpi_get_unit() result(unt)

!Local variables-------------------
 logical :: isopen

! *************************************************************************

 do unt=1024,-1,-1
   inquire(unit=unt, opened=isopen)
   if (.not.isopen) exit
 end do

end function xmpi_get_unit
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_end
!! NAME
!!  xmpi_end
!!
!! FUNCTION
!!  Hides MPI_FINALIZE from MPI library.
!!
!! INPUTS
!!  None
!!
!! PARENTS
!!      aim,atdep,band2eps,conducti,cut3d,fold2Bloch,lapackprof
!!      m_multibinit_driver,macroave,mrggkk,optic,testtransposer,ujdet
!!      vdw_kernelgen
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_end()

!Local variables-------------------
 integer :: mpierr

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI
 call MPI_BARRIER(MPI_COMM_WORLD,mpierr)  !  Needed by some HPC architectures (MT, 20110315)
 call MPI_FINALIZE(mpierr)
#endif

#ifndef FC_IBM
 ! IBM8 returns 260. 320 ...
 call sys_exit(0)
#endif

end subroutine xmpi_end
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_abort
!! NAME
!!  xmpi_abort
!!
!! FUNCTION
!!  Hides MPI_ABORT from MPI library.
!!
!! INPUTS
!!  [comm]=communicator of tasks to abort.
!!  [mpierr]=Error code to return to invoking environment.
!!  [msg]=User message
!!  [exit_status]=optional, shell return code, default 1
!!
!! PARENTS
!!      m_errors,m_initcuda,m_libpaw_tools,m_mpinfo,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_abort(comm,mpierr,msg,exit_status)

!Arguments-------------------------
 integer,optional,intent(in) :: comm,mpierr,exit_status
 character(len=*),optional,intent(in) :: msg

!Local variables-------------------
 integer :: ierr,my_comm,my_errorcode,ilen,ierr2
 logical :: testopen
 character(len=xmpi_msg_len) :: mpi_msg_error

! *************************************************************************

 ierr=0
 my_comm = xmpi_world; if (PRESENT(comm)) my_comm = comm

 if (PRESENT(msg)) then
   write(std_out,'(2a)')"User message: ",TRIM(msg)
 end if

 ! Close std_out and ab_out and flush units.
 ! Note that flush does not guarantee that the data is committed to disk.
 ! This is rather annoying because we may end up with incomplete log files
 ! that cannot be parsed by Abinit.
 ! For a possible approach based on fsync, see
 ! https://gcc.gnu.org/onlinedocs/gcc-4.7.4/gfortran/FLUSH.html

 inquire(std_out, opened=testopen)
 if (testopen) then
#if defined HAVE_FC_FLUSH
   call flush(std_out)
#endif
   close(std_out)
 end if

 inquire(ab_out,opened=testopen)
 if (testopen) then
#if defined HAVE_FC_FLUSH
   call flush(ab_out)
#endif
   close(ab_out)
 end if

#ifdef HAVE_MPI
 my_errorcode=MPI_ERR_UNKNOWN; if (PRESENT(mpierr)) my_errorcode=mpierr

 call MPI_ERROR_STRING(my_errorcode, mpi_msg_error, ilen, ierr2)

 !if (ilen>xmpi_msg_len) write(std_out,*)" WARNING: MPI message has been truncated!"
 !if (ierr2/=MPI_SUCCESS) then
 !  write(std_out,'(a,i0)')" WARNING: MPI_ERROR_STRING returned ierr2= ",ierr2
 !else
 !  write(std_out,'(2a)')" MPI_ERROR_STRING: ",TRIM(mpi_msg_error)
 !end if

 call MPI_ABORT(my_comm, my_errorcode, ierr)
#endif

 if (present(exit_status)) then
   call sys_exit(exit_status)
 else
   call sys_exit(1)
 end if

end subroutine xmpi_abort
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/sys_exit
!! NAME
!! sys_exit
!!
!! FUNCTION
!! Routine for clean exit of f90 code by one processor
!!
!! INPUTS
!!   exit_status:
!!     return code.
!!
!! NOTES
!!  By default, it uses "call exit(1)", that is not completely portable.
!!
!! PARENTS
!!      m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine sys_exit(exit_status)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: exit_status

! **********************************************************************

#if defined FC_NAG
 call exit(exit_status)
#elif defined HAVE_FC_EXIT
 call exit(exit_status)
#else
 ! stop with exit_status
 ! MT 06-2013:stop function only accept parameters !
 if (exit_status== 0) stop  "0"
 if (exit_status== 1) stop  "1"
 if (exit_status==-1) stop "-1"
#endif
 stop 1

end subroutine sys_exit
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_show_info
!! NAME
!!  xmpi_show_info
!!
!! FUNCTION
!!  Printout of the most important variables stored in this module (useful for debugging).
!!
!! INPUTS
!!  unt=Unit number for formatted output.
!!
!! PARENTS
!!      abinit,m_argparse,m_errors
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_show_info(unit)

!Arguments-------------------------
 integer,optional,intent(in) :: unit

!Local variables-------------------
 integer :: my_unt

! *************************************************************************

 !@m_xmpi
 my_unt = std_out; if (PRESENT(unit)) my_unt=unit

#ifdef HAVE_MPI1
  write(my_unt,*)" ==== Using MPI-1 specifications ==== "
#endif
#ifdef HAVE_MPI2
  write(my_unt,*)" ==== Using MPI-2 specifications ==== "
#endif

#ifdef HAVE_MPI_IO
  write(my_unt,*)" MPI-IO support is ON"
#else
  write(my_unt,*)" MPI-IO support is OFF"
#endif

#ifdef HAVE_MPI
 write(my_unt,*)" xmpi_tag_ub ................ ",xmpi_tag_ub
 write(my_unt,*)" xmpi_bsize_ch .............. ",xmpi_bsize_ch
 write(my_unt,*)" xmpi_bsize_int ............. ",xmpi_bsize_int
 write(my_unt,*)" xmpi_bsize_sp .............. ",xmpi_bsize_sp
 write(my_unt,*)" xmpi_bsize_dp .............. ",xmpi_bsize_dp
 write(my_unt,*)" xmpi_bsize_spc ............. ",xmpi_bsize_spc
 write(my_unt,*)" xmpi_bsize_dpc ............. ",xmpi_bsize_dpc
 write(my_unt,*)" xmpio_bsize_frm ............ ",xmpio_bsize_frm
 write(my_unt,*)" xmpi_address_kind .......... ",xmpi_address_kind
 write(my_unt,*)" xmpi_offset_kind ........... ",xmpi_offset_kind
 write(my_unt,*)" MPI_WTICK .................. ",MPI_WTICK()
#endif

end subroutine xmpi_show_info
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_rank
!! NAME
!!  xmpi_comm_rank
!!
!! FUNCTION
!!  Hides MPI_COMM_RANK from MPI library.
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  xmpi_comm_rank=The rank of the node inside comm
!!
!! PARENTS
!!
!! SOURCE

function xmpi_comm_rank(comm)

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: xmpi_comm_rank

!Local variables-------------------
 integer :: mpierr

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI
 xmpi_comm_rank=-1  ! Return non-sense value if the proc does not belong to the comm
 if (comm/=xmpi_comm_null) then
   call MPI_COMM_RANK(comm,xmpi_comm_rank,mpierr)
 end if
#else
 xmpi_comm_rank=0
#endif

end function xmpi_comm_rank
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_size
!! NAME
!!  xmpi_comm_size
!!
!! FUNCTION
!!  Hides MPI_COMM_SIZE from MPI library.
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  xmpi_comm_size=The number of processors inside comm.
!!
!! PARENTS
!!
!! SOURCE

function xmpi_comm_size(comm)

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: xmpi_comm_size

!Local variables-------------------------------
!scalars
 integer :: mpierr

! *************************************************************************

 mpierr=0; xmpi_comm_size=1
#ifdef HAVE_MPI
 if (comm/=xmpi_comm_null) then
   call MPI_COMM_SIZE(comm,xmpi_comm_size,mpierr)
 end if
#endif

end function xmpi_comm_size
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_free_0D
!! NAME
!!  xmpi_comm_free_0D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library.
!!  Does not abort MPI in case of an invalid communicator
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_free_0D(comm)

!Arguments-------------------------
 integer,intent(inout) :: comm

!Local variables-------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: comm_world,err_handler_dum,err_handler_sav,ierr,mpierr,mpierr_class

! *************************************************************************

 if (comm/=xmpi_comm_null.and.comm/=xmpi_world.and.comm/=xmpi_comm_self) then

   comm_world=xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
   call xmpi_comm_set_errhandler(comm_world,MPI_ERRORS_RETURN,err_handler_sav,ierr)
   call MPI_COMM_FREE(comm,mpierr)
   call xmpi_comm_set_errhandler(comm_world,err_handler_sav,err_handler_dum,ierr)

   if (mpierr/=MPI_SUCCESS) then
     call MPI_ERROR_CLASS(mpierr,mpierr_class,ierr)
     if (mpierr_class/=MPI_ERR_COMM) then
       write(std_out,*)" WARNING: MPI_COMM_FREE returned ierr= ",mpierr
     end if
   end if

 end if

#else
 if (.false.) write(std_out,*) comm
#endif

end subroutine xmpi_comm_free_0D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_free_1D
!! NAME
!!  xmpi_comm_free_1D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 1D arrays
!!  Does not abort MPI in case of an invalid communicator
!!
!! INPUTS
!!  comms(:)=MPI communicators
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_free_1D(comms)

!Arguments-------------------------
 integer,intent(inout) :: comms(:)

!Local variables-------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: comm_world,err_handler_dum,err_handler_sav,ii,mpierr

! *************************************************************************

 comm_world=xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
 call xmpi_comm_set_errhandler(comm_world,MPI_ERRORS_RETURN,err_handler_sav,mpierr)

 do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
   if (comms(ii)/=xmpi_comm_null.and.comms(ii)/=xmpi_world.and.comms(ii)/=xmpi_comm_self) then
     call MPI_COMM_FREE(comms(ii),mpierr)
   end if
 end do

 call xmpi_comm_set_errhandler(comm_world,err_handler_sav,err_handler_dum,mpierr)

#else
 if (.false.) write(std_out,*) comms(1)
#endif

end subroutine xmpi_comm_free_1D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_free_2D
!! NAME
!!  xmpi_comm_free_2D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 2D arrays
!!  Does not abort MPI in case of an invalid communicator
!!
!! INPUTS
!!  comms=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_free_2D(comms)

!Arguments-------------------------
 integer,intent(inout) :: comms(:,:)

!Local variables-------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: comm_world,err_handler_dum,err_handler_sav,ii,jj,mpierr

! *************************************************************************

 comm_world=xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
 call xmpi_comm_set_errhandler(comm_world,MPI_ERRORS_RETURN,err_handler_sav,mpierr)

 do jj=LBOUND(comms,DIM=2),UBOUND(comms,DIM=2)
   do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
     if (comms(ii,jj)/=xmpi_comm_null.and.comms(ii,jj)/=xmpi_world.and. &
&        comms(ii,jj)/=xmpi_comm_self) then
       call MPI_COMM_FREE(comms(ii,jj),mpierr)
     end if
   end do
 end do

 call xmpi_comm_set_errhandler(comm_world,err_handler_sav,err_handler_dum,mpierr)

#else
 if (.false.) write(std_out,*) comms(1,1)
#endif

end subroutine xmpi_comm_free_2D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_free_3D
!! NAME
!!  xmpi_comm_free_3D
!!
!! FUNCTION
!!  Hides MPI_COMM_FREE from MPI library. Target 3D arrays
!!  Does not abort MPI in case of an invalid communicator
!!
!! INPUTS
!!  comms=MPI communicator.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_free_3D(comms)

!Arguments-------------------------
 integer,intent(inout) :: comms(:,:,:)

!Local variables-------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: comm_world,err_handler_dum,err_handler_sav,ii,jj,kk,mpierr

! *************************************************************************

 comm_world=xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
 call xmpi_comm_set_errhandler(comm_world,MPI_ERRORS_RETURN,err_handler_sav,mpierr)

 do kk=LBOUND(comms,DIM=3),UBOUND(comms,DIM=3)
   do jj=LBOUND(comms,DIM=2),UBOUND(comms,DIM=2)
     do ii=LBOUND(comms,DIM=1),UBOUND(comms,DIM=1)
       if (comms(ii,jj,kk)/=xmpi_comm_null.and.comms(ii,jj,kk)/=xmpi_world.and. &
&          comms(ii,jj,kk)/=xmpi_comm_self) then
         call MPI_COMM_FREE(comms(ii,jj,kk),mpierr)
       end if
     end do
   end do
 end do

 call xmpi_comm_set_errhandler(comm_world,err_handler_sav,err_handler_dum,mpierr)

#else
 if (.false.) write(std_out,*) comms(1,1,1)
#endif

end subroutine xmpi_comm_free_3D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_group_free
!! NAME
!!  xmpi_group_free
!!
!! FUNCTION
!!  Hides MPI_GROUP_FREE from MPI library.
!!  Does not abort MPI in case of an invalid group
!!
!! INPUTS
!!  spaceGroup=MPI group
!!
!! PARENTS
!!      m_paw_tools,m_wfd,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_group_free(spaceGroup)

!Arguments-------------------------
 integer,intent(inout) :: spaceGroup

!Local variables-------------------------------
!scalars
#ifdef HAVE_MPI
 integer :: comm_world,err_handler_dum,err_handler_sav,ierr,mpierr,mpierr_class

! *************************************************************************

 if (spaceGroup/=xmpi_group_null) then

   comm_world=xmpi_world ! Needed to bypass a bug in some OMPI implementations (intent(inout))
   call xmpi_comm_set_errhandler(comm_world,MPI_ERRORS_RETURN,err_handler_sav,ierr)
   call MPI_GROUP_FREE(spaceGroup,mpierr)
   call xmpi_comm_set_errhandler(comm_world,err_handler_sav,err_handler_dum,ierr)

   if (mpierr/=MPI_SUCCESS) then
     call MPI_ERROR_CLASS(mpierr,mpierr_class,ierr)
     if (mpierr_class/=MPI_ERR_GROUP) write(std_out,*)" WARNING: MPI_GROUP_FREE returned ierr= ",mpierr
   end if

 end if

#else
 if (.false.) write(std_out,*) spaceGroup
#endif

end subroutine xmpi_group_free
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_group_incl
!! NAME
!!  xmpi_group_incl
!!
!! FUNCTION
!!  Hides MPI_GROUP_INCL from MPI library.
!!
!! INPUTS
!!  group=input group
!!  nrank=number of elements in array ranks (size of newgroup)
!!  ranks=ranks of processes in group to appear in newgroup
!!
!! OUTPUT
!!  newgroup= new group derived from above, in the order defined by ranks
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_group_incl(group,nranks,ranks,newgroup,mpierr)

!Arguments-------------------------
!scalars
 integer,intent(in) :: group,nranks
 integer,intent(out) :: mpierr
 integer,intent(inout) :: newgroup
!arrays
 integer,intent(in) :: ranks(nranks)

! *************************************************************************

 mpierr=0 ; newgroup=xmpi_group_null
#ifdef HAVE_MPI
 if (group/=xmpi_group_null) then
   call MPI_GROUP_INCL(group,nranks,ranks,newgroup,mpierr)
 end if
#endif

end subroutine xmpi_group_incl
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_create
!! NAME
!!  xmpi_comm_create
!!
!! FUNCTION
!!  Hides MPI_COMM_CREATE from MPI library.
!!
!! INPUTS
!!  comm=communicator
!!  group=group, which is a subset of the group of comm
!!
!! OUTPUT
!!  newcomm=new communicator
!!
!! PARENTS
!!      m_wfd
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_create(comm,group,newcomm,mpierr)

!Arguments-------------------------
!scalars
 integer,intent(in) :: comm,group
 integer,intent(out) :: mpierr
 integer,intent(inout) :: newcomm

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI
 if (group/=xmpi_group_null) then
   call MPI_comm_create(comm,group,newcomm,mpierr)
 else
   newcomm=xmpi_comm_null
 end if
#else
  newcomm=xmpi_comm_self
#endif

end subroutine xmpi_comm_create
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_subcomm
!! NAME
!!  xmpi_subcomm
!!
!! FUNCTION
!!  Return a sub-communicator from an input communicator and a given proc. ranks set.
!!  (hides subgroup creation/destruction)
!!
!! INPUTS
!!  comm=input communicator
!!  nrank=number of elements in array ranks (size of subcomm)
!!  ranks=ranks of processes in group to appear in subcomm
!!
!! OUTPUT
!!  [my_rank_in_group]=optional: my rank in the group of new sub-communicator
!!  xmpi_subcomm=new (sub-)communicator
!!
!! PARENTS
!!
!! SOURCE

function xmpi_subcomm(comm,nranks,ranks,my_rank_in_group)

!Arguments-------------------------
!scalars
 integer,intent(in) :: comm,nranks
 integer,intent(out),optional :: my_rank_in_group
 integer :: xmpi_subcomm
!arrays
 integer,intent(in) :: ranks(nranks)

!Local variables-------------------------------
#ifdef HAVE_MPI
 integer :: group,ierr,subgroup
#endif

! *************************************************************************

 xmpi_subcomm=xmpi_comm_null
 if (present(my_rank_in_group)) my_rank_in_group=xmpi_undefined

#ifdef HAVE_MPI
 if (comm/=xmpi_comm_null.and.nranks>=0) then
   call MPI_COMM_GROUP(comm,group,ierr)
   call MPI_GROUP_INCL(group,nranks,ranks,subgroup,ierr)
   call MPI_COMM_CREATE(comm,subgroup,xmpi_subcomm,ierr)
   if ( nranks == 0 )xmpi_subcomm=xmpi_comm_self
   if (present(my_rank_in_group)) then
     call MPI_Group_rank(subgroup,my_rank_in_group,ierr)
   end if
   call MPI_GROUP_FREE(subgroup,ierr)
   call MPI_GROUP_FREE(group,ierr)
 end if
#else
 if (nranks>0) then
   if (ranks(1)==0) then
     xmpi_subcomm=xmpi_comm_self
     if (present(my_rank_in_group)) my_rank_in_group=0
   end if
 end if
#endif

end function xmpi_subcomm
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_group
!! NAME
!!  xmpi_comm_group
!!
!! FUNCTION
!!  Hides MPI_COMM_GROUP from MPI library.
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  spaceGroup=The group associated to comm.
!!  mpierr=error code returned
!!
!! PARENTS
!!      m_paw_tools,m_wfd,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_group(comm,spaceGroup,mpierr)

!Arguments-------------------------
 integer,intent(in) :: comm
 integer,intent(out) :: mpierr,spaceGroup

! *************************************************************************

 mpierr=0; spaceGroup=xmpi_group_null
#ifdef HAVE_MPI
 if (comm/=xmpi_comm_null) then
   call MPI_COMM_GROUP(comm,spaceGroup,mpierr)
 end if
#endif

end subroutine xmpi_comm_group
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_split
!! NAME
!!  xmpi_comm_split
!!
!! FUNCTION
!!  Hides MPI_COMM_SPLIT from MPI library.
!!
!! INPUTS
!!  input_comm=Input MPI communicator (to be splitted)
!!  color=Control of subset assignment (nonnegative integer).
!!        Processes with the same color are in the same new communicator
!!  key=Control of rank assigment (integer)
!!
!! OUTPUT
!!  mpierr=error code returned
!!  output_comm=new splitted communicator
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_split(input_comm,color,key,output_comm,mpierr)

!Arguments-------------------------
!scalars
 integer,intent(in) :: color,input_comm,key
 integer,intent(out) :: mpierr,output_comm

! *************************************************************************

 mpierr=0; output_comm=input_comm
#ifdef HAVE_MPI
 if (input_comm/=xmpi_comm_null.and.input_comm/=xmpi_comm_self) then
   call MPI_COMM_SPLIT(input_comm,color,key,output_comm,mpierr)
 end if
#endif

end subroutine xmpi_comm_split
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_group_translate_ranks
!! NAME
!!  xmpi_group_translate_ranks
!!
!! FUNCTION
!!  Hides MPI_GROUP_TRANSLATE_RANKS from MPI library.
!!
!! INPUTS
!!  nrank=number of ranks in ranks1 and ranks2 arrays
!!  ranks1(nrank)=array of zero or more valid ranks in group1
!!  spaceGroup1=group1
!!  spaceGroup2=group2
!!
!! OUTPUT
!!  mpierr=error code returned
!!  ranks2(nrank)=array of corresponding ranks in group2,
!!                xmpi_undefined when no correspondence exists
!!
!! PARENTS
!!      m_paw_tools,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_group_translate_ranks(spaceGroup1,nrank,ranks1,&
&                                     spaceGroup2,ranks2,mpierr)

!Arguments-------------------------
!scalars
 integer,intent(in) :: nrank,spaceGroup1,spaceGroup2
 integer,intent(out) :: mpierr
!arrays
 integer,intent(in) :: ranks1(nrank)
 integer,intent(out) :: ranks2(nrank)

! *************************************************************************

 mpierr=0; ranks2(:)=xmpi_undefined
#ifdef HAVE_MPI
 if (spaceGroup1/=xmpi_group_null.and.spaceGroup2/=xmpi_group_null) then
   call MPI_GROUP_TRANSLATE_RANKS(spaceGroup1,nrank,ranks1, spaceGroup2,ranks2,mpierr)
 end if
#else
 ranks2(1)=0
#endif

end subroutine xmpi_group_translate_ranks
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_translate_ranks
!! NAME
!!  xmpi_comm_translate_ranks
!!
!! FUNCTION
!!  Helper function that translate the ranks from a communicator to another one.
!!  Wraps xmpi_group_translate_ranks but provides a more user-friendly interface
!!
!! INPUTS
!!  from_comm=MPI communicator where from_ranks are defined.
!!  nrank=number of ranks in from_ranks and to_ranks arrays
!!  from_ranks(nrank)=array of zero or more valid ranks in from_comm
!!
!! OUTPUT
!!  to_ranks(nrank)=array of corresponding ranks in to_comm
!!                xmpi_undefined when no correspondence exists
!!
!! PARENTS
!!      m_paral_pert
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_comm_translate_ranks(from_comm,nrank,from_ranks,to_comm,to_ranks)

!Arguments-------------------------
!scalars
 integer,intent(in) :: nrank,from_comm,to_comm
!arrays
 integer,intent(in) :: from_ranks(nrank)
 integer,intent(out) :: to_ranks(nrank)

!Local variables-------------------------------
!scalars
 integer :: ierr,from_group,to_group

! *************************************************************************

 ! Get the groups
 call xmpi_comm_group(from_comm,from_group,ierr)
 call xmpi_comm_group(to_comm,to_group,ierr)

 call xmpi_group_translate_ranks(from_group,nrank,from_ranks,to_group,to_ranks,ierr)

 ! Release the groups
 call xmpi_group_free(from_group)
 call xmpi_group_free(to_group)

end subroutine xmpi_comm_translate_ranks
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_barrier
!! NAME
!!  xmpi_barrier
!!
!! FUNCTION
!!  Hides MPI_BARRIER from MPI library.
!!
!! INPUTS
!!  comm=MPI communicator
!!
!! PARENTS
!!      m_Ctqmcoffdiag,m_abihist,m_alloc_hamilt_gpu,m_bse_io,m_calc_ucrpa
!!      m_chebfi,m_datafordmft,m_ddk,m_dfpt_looppert,m_dfpt_nstwf,m_dfpt_scfcv
!!      m_dtfil,m_dvdb,m_errors,m_exc_build,m_exc_diago,m_exc_itdiago
!!      m_exc_spectra,m_fit_polynomial_coeff,m_forctqmc,m_green,m_gstateimg
!!      m_haydock,m_hdr,m_io_kss,m_io_redirect,m_ioarr,m_iowf,m_ksdiago,m_mkrho
!!      m_mlwfovlp,m_mover_effpot,m_paw_mkaewf,m_paw_mkrho,m_plowannier
!!      m_polynomial_coeff,m_precpred_1geo,m_primitive_potential_list
!!      m_rf2_init,m_sigma_driver,m_sigmaph,m_slk,m_spmat_csr,m_tddft,m_vtorho
!!      m_vtorhorec,m_wfd,m_wfd_optic,m_wffile,m_wfk,m_wfk_analyze
!!      testtransposer
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_barrier(comm)

!Arguments-------------------------
 integer,intent(in) :: comm

!Local variables-------------------
 integer   :: ier
#ifdef HAVE_MPI
 integer :: nprocs
#endif

! *************************************************************************

 ier = 0
#ifdef HAVE_MPI
 if (comm/=xmpi_comm_null) then
   call MPI_COMM_SIZE(comm,nprocs,ier)
   if(nprocs>1) call MPI_BARRIER(comm,ier)
 end if
#endif

end subroutine xmpi_barrier
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_name
!! NAME
!!  xmpi_name
!!
!! FUNCTION
!!  Hides MPI_GET_PROCESSOR_NAME from MPI library.
!!
!! OUTPUT
!!  name= the host name transformed to integer variable.
!!  mpierr=Status error.
!!
!! PARENTS
!!      m_gpu_detect
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_name(name_ch, mpierr)

!Arguments-------------------------
 integer,intent(out) ::  mpierr
 character(20),intent(out) :: name_ch

!Local variables-------------------
 integer :: name,len
! character(len=MPI_MAX_PROCESSOR_NAME) :: name_ch

! *************************************************************************
!Get the name of this processor (usually the hostname)

 name   = 0
 mpierr = 0

#ifdef HAVE_MPI
 call MPI_GET_PROCESSOR_NAME(name_ch, len, mpierr)
 name_ch = trim(name_ch)

#else
 name_ch ='0'
#endif

end subroutine xmpi_name
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_iprobe
!! NAME
!!  xmpi_iprobe
!!
!! FUNCTION
!!  Hides MPI_IPROBE from MPI library.
!!  Nonblocking test for a message.
!!
!! INPUTS
!!  source= source processes
!!  tag= tag value
!!  mpicomm= communicator
!!
!! OUTPUT
!!  flag= True if a message with the specified source, tag, and communicator is available
!!  mpierr= status error
!!
!! PARENTS
!!      m_paw_an,m_paw_ij,m_pawfgrtab,m_pawrhoij
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_iprobe(source,tag,mpicomm,flag,mpierr)

!Arguments-------------------------
 integer,intent(in) :: mpicomm,source,tag
 integer,intent(out) :: mpierr
 logical,intent(out) :: flag

!Local variables-------------------
#ifdef HAVE_MPI
 integer :: ier,status(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 mpierr = 0
#ifdef HAVE_MPI
  call MPI_IPROBE(source,tag,mpicomm,flag,status,ier)
  mpierr=ier
#endif

end subroutine xmpi_iprobe
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_wait
!! NAME
!!  xmpi_wait
!!
!! FUNCTION
!!  Hides MPI_WAIT from MPI library.
!!  Waits for an MPI request to complete.
!!
!! INPUTS
!!  request= MPI request handle to wait for
!!
!! OUTPUT
!!  mpierr= status error
!!
!! PARENTS
!!      m_dfpt_scfcv,m_dvdb,m_fftw3,m_mover,m_paw_an,m_paw_ij,m_paw_occupancies
!!      m_pawfgrtab,m_pawrhoij,m_scfcv_core,m_sg2002
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_wait(request, mpierr)

!Arguments-------------------------
 integer,intent(inout) :: request
 integer,intent(out) :: mpierr

!Local variables-------------------
#ifdef HAVE_MPI
 integer :: ier,status(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 mpierr = 0
#ifdef HAVE_MPI
 if (request /= xmpi_request_null) xmpi_count_requests = xmpi_count_requests - 1
 call MPI_WAIT(request,status,ier)
 mpierr=ier
#endif

end subroutine xmpi_wait
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_waitall_1d
!! NAME
!!  xmpi_waitall_1d
!!
!! FUNCTION
!!  Hides MPI_WAITALL from MPI library.
!!  Waits for all given MPI Requests to complete.
!!
!! INPUTS
!!  array_of_requests= array of request handles
!!
!! OUTPUT
!!  mpierr= status error
!!
!! PARENTS
!!      m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_waitall_1d(array_of_requests, mpierr)

!Arguments-------------------------
 integer,intent(inout) :: array_of_requests(:)
 integer,intent(out) :: mpierr

!Local variables-------------------
#ifdef HAVE_MPI
 integer :: ier,status(MPI_STATUS_SIZE,size(array_of_requests))
#endif

! *************************************************************************

 mpierr = 0
#ifdef HAVE_MPI
 xmpi_count_requests = xmpi_count_requests - count(array_of_requests /= xmpi_request_null)
 call MPI_WAITALL(size(array_of_requests), array_of_requests, status, ier)
 mpierr=ier
#endif

end subroutine xmpi_waitall_1d
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_waitall_2d
!! NAME
!!  xmpi_waitall_2d
!!
!! FUNCTION
!!  Hides MPI_WAITALL from MPI library.
!!  Waits for all given MPI Requests to complete.
!!
!! INPUTS
!!  array_of_requests= array of request handles
!!
!! OUTPUT
!!  mpierr= status error
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_waitall_2d(array_of_requests, mpierr)

!Arguments-------------------------
 integer,intent(inout) :: array_of_requests(:,:)
 integer,intent(out) :: mpierr

!Local variables-------------------
 integer :: flat_requests(product(shape(array_of_requests)))

! *************************************************************************

 ! MPI_WAITALL is a Fortran interface so cannot pass count and base address a la C
 ! so flat 2d array and copy in-out. See https://github.com/open-mpi/ompi/issues/587
 flat_requests = pack(array_of_requests, mask=.True.)
 call xmpi_waitall_1d(flat_requests, mpierr)
 array_of_requests = reshape(flat_requests, shape(array_of_requests))

end subroutine xmpi_waitall_2d
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_request_free
!! NAME
!!  xmpi_request_free
!!
!! FUNCTION
!!  Hides MPI_REQUEST_FREE from MPI library.
!!  Frees an array of communication request objects.
!!
!! INPUTS
!!  requests(:)= communication request array (array of handles)
!!
!! OUTPUT
!!  mpierr= status error
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_request_free(requests,mpierr)

!Arguments-------------------------
 integer,intent(inout) :: requests(:)
 integer,intent(out)  :: mpierr

!Local variables-------------------
#ifdef HAVE_MPI
 integer :: ier,ii
#endif

! *************************************************************************

 mpierr = 0
#ifdef HAVE_MPI
 do ii=1,size(requests)
   if (requests(ii) /= xmpi_request_null) xmpi_count_requests = xmpi_count_requests - 1
   call MPI_REQUEST_FREE(requests(ii),ier)
 end do
 mpierr=ier
#endif

end subroutine xmpi_request_free
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_error_string
!! NAME
!!  xmpi_error_string
!!
!! FUNCTION
!!  Hides MPI_ERROR_STRING from MPI library.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_error_string(mpierr,err_string,ilen,ierror)

!Arguments-------------------------
 integer,intent(in) :: mpierr
 integer,intent(out) :: ilen,ierror
 character(len=*),intent(out) :: err_string

! *************************************************************************

 ilen=0
#ifdef HAVE_MPI
 call MPI_Error_string(mpierr,err_string,ilen,ierror)
#else
 ierror=1
 err_string="Sorry, no MPI_Error_string routine is available to interpret the error message"
#endif

end subroutine xmpi_error_string
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_comm_set_errhandler
!! NAME
!!  xmpi_set_errhandler
!!
!! FUNCTION
!!  Hides MPI_COMM_SET_ERRHANDLER from MPI library.
!!
!! INPUTS
!!  new_err_handler= new error handler
!!
!! OUTPUT
!!  ierror=error code
!!  old_err_handler= old error handler
!!
!! SIZE EFFECTS
!!  comm= communicator (should be intent(in) but is intent(inout) in some
!!             OMPI implementation ; known as a bug)
!!
!! PARENTS
!!
!! SOURCE

subroutine xmpi_comm_set_errhandler(comm,new_err_handler,old_err_handler,ierror)

!Arguments-------------------------
 integer,intent(in) :: new_err_handler
 integer,intent(in) :: comm
 integer,intent(out) :: ierror,old_err_handler

!Local variables-------------------------
 integer :: mpierr1,mpierr2,my_comm

! *************************************************************************

 ierror=0
 my_comm = comm  !should be intent(in) but is intent(inout) in some OMPI implementation ; known as a bug)

#if defined HAVE_MPI

 mpierr1=MPI_SUCCESS; mpierr2=MPI_SUCCESS

#if defined HAVE_MPI1
   call MPI_Errhandler_get(my_comm,old_err_handler,mpierr1)
   call MPI_Errhandler_set(my_comm,new_err_handler,mpierr2)
#endif
#if defined HAVE_MPI2
   call MPI_comm_get_Errhandler(my_comm,old_err_handler,mpierr1)
   call MPI_comm_set_Errhandler(my_comm,new_err_handler,mpierr2)
#endif

 ierror=MPI_SUCCESS
 if (mpierr1/=MPI_SUCCESS) then
   ierror=mpierr1
 else if (mpierr2/=MPI_SUCCESS) then
   ierror=mpierr2
 end if
#endif

end subroutine xmpi_comm_set_errhandler
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work_i4b
!! NAME
!!  xmpi_split_work_i4b
!!
!! FUNCTION
!!  Splits the number of tasks, ntasks, among nprocs processors.
!!  Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks=number of tasks
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  my_start,my_stop= indices defining the initial and final task for this processor
!!
!! NOTES
!!  If nprocs > ntasks then:
!!
!!    my_start = ntasks + 1
!!    my_stop = ntask
!!
!!  In this particular case, loops of the form
!!
!!  do ii=my_start,my_stop
!!   ...
!!  end do
!!
!!  are not executed. Moreover allocation such as foo(my_start:my_stop) will generate a zero-sized array.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_work_i4b(ntasks, comm, my_start, my_stop)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,comm
 integer,intent(out) :: my_start, my_stop

!Local variables-------------------------------
 integer :: res,nprocs,my_rank,block_p1,block

! *************************************************************************

 nprocs  = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 block   = ntasks/nprocs
 res     = MOD(ntasks,nprocs)
 block_p1= block+1

 if (my_rank < res) then
   my_start =  my_rank   *block_p1+1
   my_stop  = (my_rank+1)*block_p1
 else
   my_start = res*block_p1 + (my_rank-res  )*block + 1
   my_stop  = res*block_p1 + (my_rank-res+1)*block
 end if

end subroutine xmpi_split_work_i4b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_block
!! NAME
!!  xmpi_split_block
!!
!! FUNCTION
!!  Splits tasks inside communicator using cyclic distribution.
!!  Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks: number of tasks
!!  comm: MPI communicator.
!!
!! OUTPUT
!!  my_ntasks: Number of tasks received by this rank. May be zero if ntasks > nprocs.
!!  my_inds(my_ntasks): List of tasks treated by this rank. Allocated by the routine. May be zero-sized.
!!
!! PARENTS
!!      m_sigmaph
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_block(ntasks, comm, my_ntasks, my_inds)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks, comm
 integer,intent(out) :: my_ntasks
 integer,allocatable,intent(out) :: my_inds(:)

!Local variables-------------------------------
 integer :: ii, istart, istop

! *************************************************************************

 call xmpi_split_work(ntasks, comm, istart, istop)
 my_ntasks = istop - istart + 1
 ABI_MALLOC(my_inds, (my_ntasks))
 if (my_ntasks > 0) my_inds = [(istart + (ii - 1), ii=1, my_ntasks)]

end subroutine xmpi_split_block
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_cyclic
!! NAME
!!  xmpi_split_cyclic
!!
!! FUNCTION
!!  Splits tasks inside communicator using cyclic distribution.
!!  Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks: number of tasks
!!  comm: MPI communicator.
!!
!! OUTPUT
!!  my_ntasks: Number of tasks received by this rank. May be zero if ntasks > nprocs.
!!  my_inds(my_ntasks): List of tasks treated by this rank. Allocated by the routine. May be zero-sized.
!!
!! PARENTS
!!      m_phgamma,m_sigmaph
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_cyclic(ntasks, comm, my_ntasks, my_inds)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks, comm
 integer,intent(out) :: my_ntasks
 integer,allocatable,intent(out) :: my_inds(:)

!Local variables-------------------------------
 integer :: ii, cnt, itask, my_rank, nprocs

! *************************************************************************

 nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 do ii=1,2
   if (ii == 2) then
     ABI_MALLOC(my_inds, (my_ntasks))
   end if
   cnt = 0
   do itask=1,ntasks
     if (mod(itask, nprocs) == my_rank) then
       cnt = cnt + 1
       if (ii == 2) my_inds(cnt) = itask
     end if
   end do
   if (ii == 1) my_ntasks = cnt
 end do

end subroutine xmpi_split_cyclic
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_list
!! NAME
!!  xmpi_split_list
!!
!! FUNCTION
!!  Splits list of itmes inside communicator using block distribution.
!!  Used for the MPI parallelization of simple loops.
!!
!! INPUTS
!!  ntasks:Number of items in list (global)
!!  list(ntasks): List of indices
!!  comm: MPI communicator.
!!
!! OUTPUT
!!  my_ntasks: Number of tasks received by this rank. May be zero if ntasks > nprocs.
!!  my_inds(my_ntasks): List of tasks treated by this rank. Allocated by the routine. May be zero-sized.
!!
!! PARENTS
!!      m_phgamma
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_list(ntasks, list, comm, my_ntasks, my_inds)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks, comm
 integer,intent(out) :: my_ntasks
 integer,intent(in) :: list(ntasks)
 integer,allocatable,intent(out) :: my_inds(:)

!Local variables-------------------------------
 integer :: my_start, my_stop

! *************************************************************************

 call xmpi_split_work(ntasks, comm, my_start, my_stop)

 my_ntasks = my_stop - my_start + 1

 if (my_stop >= my_start) then
   ABI_MALLOC(my_inds, (my_ntasks))
   my_inds = list(my_start:my_stop)
 else
   my_ntasks = 0
   ABI_MALLOC(my_inds, (0))
 end if

end subroutine xmpi_split_list
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work2_i4b
!! NAME
!!  xmpi_split_work2_i4b
!!
!! FUNCTION
!!  Splits a number of tasks, ntasks, among nprocs processors.
!!  The output arrays istart(1:nprocs) and istop(1:nprocs)
!!  report the starting and final task index for each CPU.
!!  Namely CPU with rank ii has to perform all the tasks between
!!  istart(ii+1) and istop(ii+1). Note the Fortran convention of using
!!  1 as first index of the array.
!!  Note, moreover, that if a proc has rank > ntasks then:
!!   istart(rank+1)=ntasks+1
!!   istop(rank+1)=ntask
!!
!!  In this particular case, loops of the form
!!
!!  do ii=istart(rank),istop(rank)
!!   ...
!!  end do
!!
!!  are not executed. Moreover allocation such as foo(istart(rank):istop(rank))
!!  will generate a zero-sized array
!!
!! INPUTS
!!  ntasks= number of tasks
!!  nprocs=Number of processors.
!!
!! OUTPUT
!!  istart(nprocs),istop(nprocs)= indices defining the initial and final task for each processor
!!
!! PARENTS
!!      m_exc_build,m_screening,m_screening_driver,m_skw
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_work2_i4b(ntasks, nprocs, istart, istop)

!Arguments ------------------------------------
 integer,intent(in)  :: ntasks,nprocs
 integer,intent(inout) :: istart(nprocs), istop(nprocs)

!Local variables-------------------------------
 integer :: res,irank,block,block_tmp

! *************************************************************************

 block_tmp = ntasks/nprocs
 res       = MOD(ntasks,nprocs)
 block     = block_tmp+1

 do irank=0,nprocs-1
   if (irank<res) then
     istart(irank+1) = irank    *block+1
     istop (irank+1) = (irank+1)*block
   else
     istart(irank+1) = res*block + (irank-res  )*block_tmp+1
     istop (irank+1) = res*block + (irank-res+1)*block_tmp
   end if
 end do

end subroutine xmpi_split_work2_i4b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_split_work2_i8b
!! NAME
!!  xmpi_split_work2_i8b
!!
!! FUNCTION
!!  Same as xmpi_split_work2_i8b but accepts 8 bytes integer.
!!
!! INPUTS
!!  ntasks= number of tasks
!!  nprocs=Number of processors.
!!
!! OUTPUT
!!  istart(nprocs),istop(nprocs)= indices defining the initial and final task for each processor
!!
!! PARENTS
!!      m_exc_build
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_split_work2_i8b(ntasks,nprocs,istart,istop)

!Arguments ------------------------------------
 integer,intent(in)  :: nprocs
 integer(i8b),intent(in)  :: ntasks
 integer(i8b),intent(inout) :: istart(nprocs),istop(nprocs)

!Local variables-------------------------------
 integer(i8b) :: res,irank,block,block_tmp

! *************************************************************************

 block_tmp = ntasks/nprocs
 res       = MOD(ntasks,INT(nprocs,KIND=i8b))
 block     = block_tmp+1

 do irank=0,nprocs-1
   if (irank<res) then
     istart(irank+1)= irank   *block+1
     istop (irank+1)=(irank+1)*block
   else
     istart(irank+1)=res*block+(irank-res  )*block_tmp+1
     istop (irank+1)=res*block+(irank-res+1)*block_tmp
   end if
 end do

end subroutine xmpi_split_work2_i8b
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_distab_4D
!! NAME
!!  xmpi_distab_4D
!!
!! FUNCTION
!!  Fill table defining the distribution of the tasks according to the number of processors involved in the
!!  calculation. For each set of indeces, the table contains the rank of the node in the MPI communicator.
!!
!! INPUTS
!!  nprocs=The number of processors performing the calculation in parallel.
!!
!! OUTPUT
!!  task_distrib(:,:,:,:) = Contains the rank of the node that is taking care of this particular set of loop indeces.
!!  Tasks are distributed across the nodes in column-major order.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_distab_4D(nprocs, task_distrib)

!Arguments ------------------------------------
 integer,intent(in) :: nprocs
!arrays
 integer,intent(inout) :: task_distrib(:,:,:,:)

!Local variables ------------------------------
!scalars
 integer :: ii,jj,n1,n2,n3,n4,ntasks,irank,remainder,ntpblock
 integer,allocatable :: list(:)

!************************************************************************

 n1= SIZE(task_distrib,DIM=1)
 n2= SIZE(task_distrib,DIM=2)
 n3= SIZE(task_distrib,DIM=3)
 n4= SIZE(task_distrib,DIM=4)
 ntasks = n1*n2*n3*n4

 ABI_MALLOC(list, (ntasks))
 list=-999

 ntpblock  = ntasks/nprocs
 remainder = MOD(ntasks,nprocs)

 if (ntpblock==0) then ! nprocs > ntasks
   do ii=1,ntasks
     list(ii) = ii-1
   end do
 else
   ii=1
   do irank=nprocs-1,0,-1 ! If remainder/=0, master will get less tasks.
     jj = ii+ntpblock-1
     if (remainder>0) then
       jj=jj+1
       remainder = remainder-1
     end if
     list(ii:jj)=irank
     ii=jj+1
   end do
 end if

 task_distrib = RESHAPE(list, [n1,n2,n3,n4])

 if (ANY(task_distrib==-999)) call xmpi_abort(msg="task_distrib == -999")

 ABI_FREE(list)

end subroutine xmpi_distab_4D
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_distrib_with_replicas
!! NAME
!!  xmpi_distrib_with_replicas
!!
!! FUNCTION
!!  This function distributes the i-th task among `nprocs` inside a MPI communicator.
!!  If nprocs > ntasks, multiple MPI ranks will be assigned to a given task.
!!
!! INPUTS
!!  itask=Index of the task (must be <= ntasks)
!!  ntasks= number of tasks
!!  rank=MPI Rank of this processor
!!  nprocs=Number of processors in the MPI communicator.
!!
!! OUTPUT
!!  True if this node will treat itask (replicas are possible if nprocs > ntasks)
!!
!! PARENTS
!!
!! SOURCE

pure function xmpi_distrib_with_replicas(itask,ntasks,rank,nprocs) result(bool)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: itask,rank,nprocs,ntasks
 logical :: bool

!Local variables-------------------------------
!scalars
 integer :: ii,mnp_pool,rk_base

! *************************************************************************

 ! If the number of processors is less than ntasks, we have max one task per processor,
 ! else we replicate the tasks inside a pool of max size mnp_pool
 if (nprocs <= ntasks) then
   bool = (MODULO(itask-1, nprocs)==rank)
 else
   mnp_pool = (nprocs / ntasks)
   !write(std_out,*)"Will duplicate itasks"
   !write(std_out,*)"mnp_pool",mnp_pool,"nprocs, ntasks",nprocs,ntasks

   rk_base = MODULO(itask-1, nprocs)
   bool = .False.
   do ii=1,mnp_pool+1
     if (rank == rk_base + (ii-1) * ntasks) then
        bool = .True.; exit
     end if
   end do
 end if

end function xmpi_distrib_with_replicas
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_largetype_create
!! NAME
!!  xmpi_largetype_create
!!
!! FUNCTION
!!  This function builds a large-count contiguous datatype made of "small" adjacent
!!  chunks (of same original type). The new type can then be used in MPI
!!  routines when the number of elements to communicate exceeds a 32bit integer.
!!
!! INPUTS
!!  largecount= total number of elements expressed as a 64bit integer
!!  inputtype= (INTEGER) input type (typically INTEGER, REAL(dp), ...)
!!  op_type= type of operation that will be applied during collective comms
!!           At present, MPI_SUM, MPI_LOR, MPI_LAND are implemented
!!
!! OUTPUT
!!  largetype= (INTEGER) new MPI type made of a serie of adjacent chunks
!!  largetype_op= (INTEGER) MPI user-defined operation associated to largetype type
!!
!! NOTE
!!  This routine is partially inspired by https://github.com/jeffhammond/BigMPI
!!  See: J.R. Hammond. A. Schafer, R. Latham,
!!       "ToINT_MAX. . . and beyond. Exploring large-count support in MPI",
!!       2014 Workshop on Exascale MPI at Supercomputing Conference
!!       MIT License (MIT)
!!       Permission is hereby granted, free of charge, to any person obtaining a copy
!!       of this software and associated documentation files (the "Software"), to deal
!!       in the Software without restriction, including without limitation the rights
!!       to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
!!       copies of the Software, and to permit persons to whom the Software is
!!       furnished to do so.
!!
!!  From MPI4 specification, this routine is useless as large-count MPI communications
!!    can be called with the use of the MPI_count datatype (instead of INTEGER).
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_largetype_create(largecount,inputtype,largetype,largetype_op,op_type)

!Arguments ------------------------------------
!scalars
 integer(KIND=int64),intent(in) :: largecount
 integer,intent(in) :: inputtype,op_type
 integer,intent(out) :: largetype,largetype_op

!Local variables-------------------------------
#ifdef HAVE_MPI
!scalars
 integer,parameter :: INT_MAX=max(1,xmpi_maxint32/2)
 integer(KIND=int32) :: cc,rr,ierr
 integer(KIND=XMPI_ADDRESS_KIND) :: extent,lb,remdisp
 integer :: chunks,remainder
!arrays
 integer(KIND=int32) :: blklens(2)
 integer(KIND=XMPI_ADDRESS_KIND) :: disps(2)
 integer :: types(2)
#endif

! *************************************************************************

#ifdef HAVE_MPI
 if (XMPI_ADDRESS_KIND<int64) &
&  call xmpi_abort(msg="Too much data to communicate for this architecture!")

!Divide data in chunks
 cc=int(largecount,kind=int32)/INT_MAX
 rr=int(largecount,kind=int32)-cc*INT_MAX

!Create user-defined datatype
 if (rr==0) then
   call MPI_TYPE_VECTOR(cc,INT_MAX,INT_MAX,inputtype,largetype,ierr)
   if (ierr==0) call MPI_TYPE_COMMIT(largetype,ierr)
 else
   call MPI_TYPE_VECTOR(cc,INT_MAX,INT_MAX,inputtype,chunks,ierr)
   call MPI_TYPE_CONTIGUOUS(rr,inputtype,remainder,ierr)
   if (ierr==0) then
     call MPI_TYPE_GET_EXTENT(inputtype,lb,extent,ierr)
     remdisp=cc*INT_MAX*extent
     blklens(1:2)=1
     disps(1)=0;disps(2)=remdisp
     types(1)=chunks;types(2)=remainder
#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
     call MPI_TYPE_CREATE_STRUCT(2,blklens,disps,types,largetype,ierr)
#else
     call MPI_TYPE_STRUCT(2,blklens,disps,types,largetype,ierr)
#endif
     if (ierr==0) then
       call MPI_TYPE_COMMIT(largetype,ierr)
       call MPI_TYPE_FREE(chunks,ierr)
       call MPI_TYPE_FREE(remainder,ierr)
     end if
   end if
 end if
 if (ierr/=0) call xmpi_abort(msg="Cannot remove ABI_MPIABORTFILE")

!Associate user-defined MPI operation
 xmpi_largetype_size=largecount ; largetype_op=-1111
 if (op_type==MPI_SUM) then
   select case(inputtype)
     case(MPI_INTEGER)
      call MPI_OP_CREATE(largetype_sum_int  ,.true.,largetype_op,ierr)
     case(MPI_REAL)
      call MPI_OP_CREATE(largetype_sum_real ,.true.,largetype_op,ierr)
     case(MPI_DOUBLE_PRECISION)
      call MPI_OP_CREATE(largetype_sum_dble ,.true.,largetype_op,ierr)
     case(MPI_COMPLEX)
      call MPI_OP_CREATE(largetype_sum_cplx ,.true.,largetype_op,ierr)
     case(MPI_DOUBLE_COMPLEX)
      call MPI_OP_CREATE(largetype_sum_dcplx,.true.,largetype_op,ierr)
   end select
 else if (op_type==MPI_LOR) then
   select case(inputtype)
     case(MPI_LOGICAL)
      call MPI_OP_CREATE(largetype_lor_log,.true.,largetype_op,ierr)
   end select
 else if (op_type==MPI_LAND) then
   select case(inputtype)
     case(MPI_LOGICAL)
      call MPI_OP_CREATE(largetype_land_log,.true.,largetype_op,ierr)
   end select
 else if (op_type==MPI_OP_NULL) then
   largetype_op=-1111
 end if
#else
 ABI_UNUSED(largecount)
 ABI_UNUSED(inputtype)
 ABI_UNUSED(largetype)
 ABI_UNUSED(largetype_op)
 ABI_UNUSED(op_type)
#endif

end subroutine xmpi_largetype_create
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_sum_int
!! NAME
!!  largetype_sum_int
!!
!! FUNCTION
!!  Routine used to overload MPI_SUM for integers
 subroutine largetype_sum_int(invec,inoutvec,len,datatype)
  integer :: len,datatype
  integer :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk)+invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_sum_int
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_sum_real
!! NAME
!!  largetype_sum_real
!!
!! FUNCTION
!!  Routine used to overload MPI_SUM for reals
 subroutine largetype_sum_real(invec,inoutvec,len,datatype)
  integer :: len,datatype
  real(sp) :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk)+invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_sum_real
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_sum_dble
!! NAME
!!  largetype_sum_dble
!!
!! FUNCTION
!!  Routine used to overload MPI_SUM for double precision reals
 subroutine largetype_sum_dble(invec,inoutvec,len,datatype)
  integer :: len,datatype
  real(dp) :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk)+invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_sum_dble
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_sum_cplx
!! NAME
!!  largetype_sum_cplx
!!
!! FUNCTION
!!  Routine used to overload MPI_SUM for complex
 subroutine largetype_sum_cplx(invec,inoutvec,len,datatype)
  integer :: len,datatype
  complex(spc) :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk)+invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_sum_cplx
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_sum_dcplx
!! NAME
!!  largetype_sum_dcplx
!!
!! FUNCTION
!!  Routine used to overload MPI_SUM for double commplex
 subroutine largetype_sum_dcplx(invec,inoutvec,len,datatype)
  integer :: len,datatype
  complex(dpc) :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk)+invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_sum_dcplx
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_lor_log
!! NAME
!!  largetype_lor_log
!!
!! FUNCTION
!!  Routine used to overload MPI_LOR for logicals
 subroutine largetype_lor_log(invec,inoutvec,len,datatype)
  integer :: len,datatype
  logical :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk).or.invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_lor_log
!!***
!--------------------------------------
!!****f* m_xmpi/largetype_lang_log
!! NAME
!!  largetype_lang_log
!!
!! FUNCTION
!!  Routine used to overload MPI_LANG for logicals
 subroutine largetype_land_log(invec,inoutvec,len,datatype)
  integer :: len,datatype
  logical :: invec(len*xmpi_largetype_size),inoutvec(len*xmpi_largetype_size)
  integer(KIND=int64) :: ii,jj,kk
  kk=0
  do ii=1,len
    do jj=1,xmpi_largetype_size
      kk=kk+1
      inoutvec(kk)=inoutvec(kk).and.invec(kk)
    end do
  end do
  ABI_UNUSED(datatype)
 end subroutine largetype_land_log
!!***

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpi_largetype_free
!! NAME
!!  xmpi_largetype_free
!!
!! FUNCTION
!!  This function release a large-count contiguous datatype.
!!
!! SIDE EFFECTS
!!  largetype= (INTEGER) MPI type to release
!!  largetype_op= (INTEGER) MPI user-defined operation associated to largetype type
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

subroutine xmpi_largetype_free(largetype,largetype_op)

!Arguments ------------------------------------
!scalars
 integer,intent(inout) :: largetype,largetype_op
!Local variables-------------------------------
#ifdef HAVE_MPI
 integer :: ierr
#endif

! *************************************************************************

#ifdef HAVE_MPI
   xmpi_largetype_size=0
   if (largetype_op/=-1111) call MPI_OP_FREE(largetype_op,ierr)
   call MPI_TYPE_FREE(largetype,ierr)
#else
 ABI_UNUSED(largetype)
 ABI_UNUSED(largetype_op)
#endif

end subroutine xmpi_largetype_free
!!***

!----------------------------------------------------------------------

! Include files providing wrappers for some of the most commonly used MPI primitives.

#include "xmpi_iallgather.finc"
#include "xmpi_allgather.finc"
#include "xmpi_allgatherv.finc"
#include "xmpi_alltoall.finc"
#include "xmpi_ialltoall.finc"
#include "xmpi_alltoallv.finc"
#include "xmpi_ialltoallv.finc"
#include "xmpi_bcast.finc"
#include "xmpi_ibcast.finc"
#include "xmpi_exch.finc"
#include "xmpi_gather.finc"
#include "xmpi_gatherv.finc"
#include "xmpi_max.finc"
#include "xmpi_min.finc"
#include "xmpi_recv.finc"
#include "xmpi_irecv.finc"
#include "xmpi_scatterv.finc"
#include "xmpi_send.finc"
#include "xmpi_isend.finc"
#include "xmpi_sum_master.finc"
#include "xmpi_sum.finc"
#include "xmpi_isum.finc"
#include "xmpi_land_lor.finc"

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_type_struct
!! NAME
!!  xmpio_type_struct
!!
!! FUNCTION
!!  Some highly non-standard MPI implementations support MPI-IO without
!!  implementing the full set of MPI-2 extensions.
!!  This wrapper will call the obsolete MPI_TYPE_STRUCT if MPI_TYPE_CREATE_STRUCT
!!  is not supported. Note that MPI_TYPE_STRUCT requires the displacement arrays
!!  to be an array of default integers whereas the argument block_displ is an array of kind XMPI_ADDRESS_KIND.
!!  The routine will abort if the displacement cannot be represented with a default integer.
!!
!! INPUTS
!! ncount= number of blocks (integer) --- also number of entries in arrays array_of_types, array_of_displacements and array_of_blocklengths
!! array_of_blocklength(ncount)=number of elements in each block (array of integer)
!! array_of_displacements(ncount)=byte displacement of each block (array of integer)
!! array_of_types(ncount)=type of elements in each block (array of handles to datatype objects)
!!
!! OUTPUT
!! new_type=new datatype (handle)
!! mpierr=MPI status error
!!
!! PARENTS
!!      m_slk,m_wffile,m_wfk,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_type_struct(ncount,block_length,block_displ,block_type,new_type,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncount
 integer,intent(out) :: new_type,mpierr
!arrays
 integer,intent(in) :: block_length(ncount),block_type(ncount)
 integer(XMPI_ADDRESS_KIND),intent(in) :: block_displ(ncount)

!Local variables-------------------
#ifndef HAVE_MPI_TYPE_CREATE_STRUCT
 integer,allocatable :: tmp_displ(:)
#endif

!************************************************************************

#ifdef HAVE_MPI_TYPE_CREATE_STRUCT
 call MPI_TYPE_CREATE_STRUCT(ncount,block_length,block_displ,block_type,new_type,mpierr)
#else

 ABI_MALLOC(tmp_displ,(ncount))
 tmp_displ = block_displ
 if (ANY(block_displ > HUGE(tmp_displ(1)) ))then
   call xmpi_abort(msg=" byte displacement cannot be represented with a default integer")
 end if

 call MPI_TYPE_STRUCT(ncount,block_length,block_displ,block_type,new_type,mpierr)
 ABI_FREE(tmp_displ)
#endif

end subroutine xmpio_type_struct
!!***

#endif

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpio_get_info_frm
!! NAME
!!  xmpio_marker_info
!!
!! FUNCTION
!!  Return the byte size of the Fortran record and its corresponding MPI_type (compiler-dependent).
!!  These two values are needed to access sequential binary Fortran files with MPI/IO routines where
!!  C-streams are used.
!!
!! INPUTS
!! comm=MPI communicator. Only master will find the values for the record marker. The results
!! are then broadcast to all the other nodes in comm.
!!
!! OUTPUT
!!  bsize_frm=Byte size of the Fortran record marker.
!!  mpi_type_frm=MPI type of the marker.
!!
!! PARENTS
!!
!! SOURCE

subroutine xmpio_get_info_frm(bsize_frm,mpi_type_frm,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 integer,intent(out) :: mpi_type_frm,bsize_frm

!Local variables-------------------------------
 integer :: my_rank
#ifdef HAVE_MPI_IO
!scalars
 integer,parameter :: master=0
 integer :: spt,ept,ii
 integer :: f90_unt,iimax,mpio_fh,bsize_int,mpierr
 integer(XMPI_OFFSET_KIND) :: offset,rml
 character(len=fnlen) :: fname
 character(len=500) :: errmsg
 logical :: file_exists
!arrays
 integer :: xvals(2),ivals(100),read_5ivals(5),ref_5ivals(5)
 integer :: rm_lengths(4)=(/4,8,2,16/)
 integer :: statux(MPI_STATUS_SIZE)
 real(dp) :: xrand(fnlen)
#endif

!************************************************************************

 bsize_frm=0; mpi_type_frm=0

 my_rank = xmpi_comm_rank(comm) !; RETURN

#ifdef HAVE_MPI_IO
 if ( my_rank == master ) then
   ! Fortran scratch files cannot have a name so have to generate a random one.
   ! cannot use pick_aname since it is higher level.
   fname = "__MPI_IO_FRM__"
   spt=LEN(trim(fname))+1; ept=spt

   inquire(file=trim(fname),exist=file_exists)

   do while (file_exists)
     call RANDOM_NUMBER(xrand(spt:ept))
     xrand(spt:ept) = 64+xrand(spt:ept)*26
     do ii=spt,ept
       fname(ii:ii) = ACHAR(NINT(xrand(ii)))
     end do
     ept = MIN(ept+1,fnlen)
     inquire(file=trim(fname),exist=file_exists)
   end do
   !
   ! Write five integers on the binary file open in Fortran mode, then try
   ! to reread the values with MPI-IO using different offsets for the record marker.
   !
   f90_unt = xmpi_get_unit()
   if (f90_unt == -1) call xmpi_abort(msg="Cannot find free unit!!")
   ! MT dec 2013: suppress the new attribute: often cause unwanted errors
   !              and theoretically useless because of the previous inquire
   open(unit=f90_unt,file=trim(fname),form="unformatted",err=10, iomsg=errmsg)

   ref_5ivals = (/(ii, ii=5,9)/)
   ivals = HUGE(1); ivals(5:9)=ref_5ivals
   write(f90_unt, err=10, iomsg=errmsg) ivals
   close(f90_unt, err=10, iomsg=errmsg)

   call MPI_FILE_OPEN(xmpi_comm_self, trim(fname), MPI_MODE_RDONLY, MPI_INFO_NULL, mpio_fh,mpierr)

   iimax=3 ! Define number of INTEGER types to be tested
#ifdef HAVE_FC_INT_QUAD
   iimax=4
#endif
   !
   ! Try to read ivals(5:9) from file.
   ii=0; bsize_frm=-1
   call MPI_TYPE_SIZE(MPI_INTEGER,bsize_int,mpierr)

   do while (bsize_frm<=0 .and. ii<iimax)
     ii=ii+1
     rml = rm_lengths(ii)
     offset = rml + 4 * bsize_int
     call MPI_FILE_READ_AT(mpio_fh,offset,read_5ivals,5,MPI_INTEGER,statux,mpierr)
     !write(std_out,*)read_5ivals
     if (mpierr==MPI_SUCCESS .and. ALL(read_5ivals==ref_5ivals) ) bsize_frm=rml
   end do

   if (ii==iimax.and.bsize_frm<=0) then
     write(std_out,'(7a)') &
       'Error during FORTRAN file record marker detection:',ch10,&
       'It was not possible to read/write a small file!',ch10,&
       'ACTION: check your access permissions to the file system.',ch10,&
       'Common sources of this problem: quota limit exceeded, R/W incorrect permissions, ...'
     call xmpi_abort()
   else
     !write(std_out,'(a,i0)')' Detected FORTRAN record mark length: ',bsize_frm
   end if

   call MPI_FILE_CLOSE(mpio_fh, mpierr)
   !
   ! Select MPI datatype corresponding to the Fortran marker.
   SELECT CASE (bsize_frm)
   CASE (4)
     mpi_type_frm=MPI_INTEGER4
   CASE (8)
     mpi_type_frm=MPI_INTEGER8
#if defined HAVE_FC_INT_QUAD && defined HAVE_MPI_INTEGER16
   CASE (16)
     mpi_type_frm=MPI_INTEGER16
#endif
   CASE (2)
     mpi_type_frm=MPI_INTEGER2
   CASE DEFAULT
     write(std_out,'(a,i0)')" Wrong bsize_frm: ",bsize_frm
     call xmpi_abort()
   END SELECT

   open(unit=f90_unt,file=trim(fname), err=10, iomsg=errmsg)
   close(f90_unt,status="delete", err=10, iomsg=errmsg)
 end if
 !
 ! Broadcast data.
 xvals = (/bsize_frm,mpi_type_frm/)
 call xmpi_bcast(xvals,master,comm,mpierr)

 bsize_frm    = xvals(1)
 mpi_type_frm = xvals(2)

 return

!HANDLE IO ERROR
10 continue
 call xmpi_abort(msg=errmsg)
#endif

end subroutine xmpio_get_info_frm
!!***

!----------------------------------------------------------------------

!!****f* m_wffile/xmpio_read_frm
!! NAME
!!  xmpio_read_frm
!!
!! FUNCTION
!!  Read the content of a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  the file pointer is modified according to the value of advance.
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  sc_mode=
!!    xmpio_single     ==> for reading by current proc.
!!    xmpio_collective ==> for collective reading.
!!  offset=MPI/IO file pointer
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next read will continue picking information
!!    off of the currect record.
!!
!! OUTPUT
!!  fmarker=Content of the Fortran record marker.
!!  mpierr= MPI error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: file pointer used to access the Fortran marker.
!!     output: new offset updated after the reading, depending on advance.
!!
!! PARENTS
!!      m_bse_io,m_exc_diago,m_exc_itdiago,m_hdr,m_io_screening,m_xmpi
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_read_frm(fh,offset,sc_mode,fmarker,mpierr,advance)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,sc_mode
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer(XMPI_OFFSET_KIND),intent(out) :: fmarker
 integer,intent(out) :: mpierr
 logical,optional,intent(in) :: advance

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,mpi_type_frm,myfh
 integer(kind=int16) :: delim_record2
 integer(kind=int32) :: delim_record4
 integer(kind=int64) :: delim_record8
#if defined HAVE_FC_INT_QUAD
 integer*16 :: delim_record16
#endif
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 !Workaround for XLF.
 myfh = fh

 bsize_frm    = xmpio_bsize_frm    ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm ! MPI type of the record marker.

 SELECT CASE (sc_mode)

 CASE (xmpio_single)

   if (bsize_frm==4) then
     call MPI_FILE_READ_AT(myfh,offset,delim_record4,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record4
   else if (bsize_frm==8) then
     call MPI_FILE_READ_AT(myfh,offset,delim_record8,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     call MPI_FILE_READ_AT(myfh,offset,delim_record16,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record16
#endif
   else if (bsize_frm==2) then
     call MPI_FILE_READ_AT(myfh,offset,delim_record2 ,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record2
   else
     call xmpi_abort(msg='Wrong record marker length!')
   end if

 CASE (xmpio_collective)

   if (bsize_frm==4) then
     call MPI_FILE_READ_AT_ALL(myfh,offset,delim_record4 ,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record4
   else if (bsize_frm==8) then
     call MPI_FILE_READ_AT_ALL(myfh,offset,delim_record8 ,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record8
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     call MPI_FILE_READ_AT_ALL(myfh,offset,delim_record16,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record16
#endif
   else if (bsize_frm==2) then
     call MPI_FILE_READ_AT_ALL(myfh,offset,delim_record2 ,1,mpi_type_frm,statux,mpierr)
     fmarker = delim_record2
   else
     call xmpi_abort(msg='Wrong record marker length!')
   end if

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for sc_mode: ",sc_mode
   call xmpi_abort(msg=msg)
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm  ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_read_frm
!!***

#endif

!------------------------------------------------------------------------------------

!!****f* m_wffile/xmpio_write_frm
!! NAME
!!  xmpio_write_frm
!!
!! FUNCTION
!!  Write a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  The file pointer is modified according to the value of advance.
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  sc_mode=
!!         xmpio_single     ==> for reading by current proc.
!!         xmpio_collective ==> for collective reading.
!!  fmarker=The content of the Fortran marker i.e. the size of the record in bytes.
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next write will continue writing data
!!    on the currect record.
!!
!! OUTPUT
!!  mpierr= error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: offset of  the Fortran marker.
!!     output: new offset updated after the writing, depending on advance.
!!
!! PARENTS
!!      m_ioarr
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_write_frm(fh,offset,sc_mode,fmarker,mpierr,advance)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: fmarker
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer,intent(out) :: mpierr
 logical,optional,intent(in) :: advance

!Local variables-------------------------------
!scalars
 integer :: myfh,bsize_frm,mpi_type_frm
 integer(XMPI_OFFSET_KIND) :: last
 integer(kind=int16)  :: delim_record2
 integer(kind=int32)  :: delim_record4
 integer(kind=int64)  :: delim_record8
#if defined HAVE_FC_INT_QUAD
 integer*16 :: delim_record16
#endif
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 bsize_frm    = xmpio_bsize_frm      ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm   ! MPI type of the record marker.
 last = offset + bsize_frm + fmarker ! position of the end marker

 SELECT CASE (sc_mode)

 CASE (xmpio_single)
   if (bsize_frm==4) then
     delim_record4 = fmarker
     call MPI_FILE_WRITE_AT(myfh,offset,delim_record4 ,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT(myfh,last,delim_record4 ,1,mpi_type_frm,statux,mpierr)

   else if (bsize_frm==8) then
     delim_record8 = fmarker
     call MPI_FILE_WRITE_AT(myfh,offset,delim_record8 ,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT(myfh,last,delim_record8 ,1,mpi_type_frm,statux,mpierr)
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     delim_record16 = fmarker
     call MPI_FILE_WRITE_AT(myfh,offset,delim_record16,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT(myfh,last,delim_record16 ,1,mpi_type_frm,statux,mpierr)
#endif
   else if (bsize_frm==2) then
     delim_record2 = fmarker
     call MPI_FILE_WRITE_AT(myfh,offset,delim_record2, 1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT(myfh,last,delim_record2 ,1,mpi_type_frm,statux,mpierr)
   else
     call xmpi_abort(msg='Wrong record marker length!')
   end if

 CASE (xmpio_collective)
   if (bsize_frm==4) then
     delim_record4 = fmarker
     call MPI_FILE_WRITE_AT_ALL(myfh,offset,delim_record4 ,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT_ALL(myfh,last,delim_record4 ,1,mpi_type_frm,statux,mpierr)
   else if (bsize_frm==8) then
     delim_record8 = fmarker
     call MPI_FILE_WRITE_AT_ALL(myfh,offset,delim_record8 ,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT_ALL(myfh,last,delim_record8 ,1,mpi_type_frm,statux,mpierr)
#if defined HAVE_FC_INT_QUAD
   else if (bsize_frm==16) then
     delim_record16 = fmarker
     call MPI_FILE_WRITE_AT_ALL(myfh,offset,delim_record16,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT_ALL(myfh,last,delim_record16 ,1,mpi_type_frm,statux,mpierr)
#endif
   else if (bsize_frm==2) then
     delim_record2 = fmarker
     call MPI_FILE_WRITE_AT_ALL(myfh,offset,delim_record2 ,1,mpi_type_frm,statux,mpierr)
     call MPI_FILE_WRITE_AT_ALL(myfh,last,delim_record2 ,1,mpi_type_frm,statux,mpierr)
   else
     call xmpi_abort(msg='Wrong record marker length!')
   end if

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for sc_mode: ",sc_mode
   call xmpi_abort(msg=msg)
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm  ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm              ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_write_frm
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fstripes
!! NAME
!!  xmpio_create_fstripes
!!
!! FUNCTION
!!  Return a MPI type that can be used to (read|write) a set of interleaved Fortran records.
!!
!!  <FRM> type(1), type(1), ... <FRM>  ! size(1) elements
!!  <FRM> type(2), type(2), ... <FRM>  ! size(2) elements
!!  <FRM> type(1), type(1), ... <FRM>  ! size(1) elements
!!  ....
!!
!! INPUTS
!!  ncount = Number of records with elements of type types(1) to (read|write)
!!  sizes(1:2) = Number of elements of each type in the two sets of record
!!  type(1:2) = MPI Type of the elements in the first and in the second record.
!!
!! OUTPUT
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker individuating the beginning of the matrix. (lets call it "base").
!!    Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  new_type=New MPI type.
!!  mpierr= MPI error code
!!
!! PARENTS
!!      m_wfk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_fstripes(ncount,sizes,types,new_type,my_offpad,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ncount
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
 integer,intent(out) :: new_type,mpierr
!arrays
 integer,intent(in) :: types(2),sizes(2)

!Local variables-------------------------------
!scalars
 integer :: type_x,type_y,bsize_frm,bsize_x,bsize_y,nx,ny,column_type
 integer(MPI_ADDRESS_KIND) :: stride

!************************************************************************

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Number of elements in the two stripes.
 nx = sizes(1)
 ny = sizes(2)

 type_x = types(1)
 type_y = types(2)

 ! Byte size of type_x and type_y
 call MPI_TYPE_SIZE(type_x,bsize_x,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_SIZE(type_y,bsize_y,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! The view starts at the first element of the first stripe.
 my_offpad = xmpio_bsize_frm

 call MPI_Type_contiguous(nx,type_x,column_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! Byte size of the Fortran record + the two markers.
 stride = nx*bsize_x + 2*bsize_frm  + ny*bsize_y + 2*bsize_frm

 ! ncount colum_type separated by stride bytes
 call MPI_Type_create_hvector(ncount,1,stride,column_type,new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(column_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine xmpio_create_fstripes
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_2D
!! NAME
!!  xmpio_create_fsubarray_2D
!!
!! FUNCTION
!!  Return a MPI type that can be used to (read|write) a 2D matrix of elements of type old_type stored in a Fortran file.
!!
!! INPUTS
!!  sizes(2)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  subsizes(2)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  array_of_starts(2)=starting coordinates of the subarray in each dimension (array of nonnegative integers >=1, <=sizes)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker individuating the beginning of the matrix. (lets call it "base").
!!    Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  new_type=New MPI type.
!!  mpierr= MPI error code
!!
!! PARENTS
!!      m_exc_build,m_exc_itdiago,m_mpiotk,m_wfk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_fsubarray_2D(sizes,subsizes,array_of_starts,old_type,new_type,my_offpad,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
 integer,intent(out) :: mpierr,new_type
!arrays
 integer,intent(in) :: sizes(2),subsizes(2),array_of_starts(2)
!Local variables-------------------------------
!scalars
 integer :: bsize_frm,bsize_old,nx,ny
 integer :: column_type,ldx
 integer(XMPI_OFFSET_KIND) :: st_x,st_y
 integer(MPI_ADDRESS_KIND) :: stride_x
 !character(len=500) :: msg

!************************************************************************

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)
 ABI_HANDLE_MPIERR(mpierr)
 !
 ! Number of columns and rows of the submatrix.
 nx = subsizes(1)
 ny = subsizes(2)

 ldx = sizes(1)
 st_x = array_of_starts(1)
 st_y = array_of_starts(2)

 ! The view starts at the first element of the submatrix.
 my_offpad = (st_x-1)*bsize_old + (st_y-1)*(ldx*bsize_old+2*xmpio_bsize_frm) + xmpio_bsize_frm

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(ny,1,stride_x,column_type,new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_COMMIT(new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(column_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine xmpio_create_fsubarray_2D
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_3D
!! NAME
!!  xmpio_create_fsubarray_3D
!!
!! FUNCTION
!!  Return a MPI type that can be used to (read|write) a 3D matrix of elements of type old_type stored in a Fortran file.
!!
!! INPUTS
!!  sizes(3)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  subsizes(3)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  array_of_starts(3)=starting coordinates of the subarray in each dimension (array of nonnegative integers >=1, <=sizes)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker individuating the beginning of the matrix. (lets call it "base").
!!    Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  new_type=New MPI type.
!!  mpierr= MPI error code
!!
!! PARENTS
!!      m_mpiotk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_fsubarray_3D(sizes,subsizes,array_of_starts,old_type,new_type,my_offpad,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: mpierr,new_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: sizes(3),subsizes(3),array_of_starts(3)
!Local variables-------------------------------
!scalars
 integer :: bsize_frm,bsize_old,nx,ny,nz
 integer :: column_type,plane_type,ldx,ldy,ldz
 integer(XMPI_OFFSET_KIND) :: st_x,st_y,st_z
 integer(MPI_ADDRESS_KIND) :: stride_x
 !character(len=500) :: msg

!************************************************************************

 bsize_frm = xmpio_bsize_frm    ! Byte size of the Fortran record marker.

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)
 ABI_HANDLE_MPIERR(mpierr)
 !
 ! Number of columns and rows of the submatrix.
 nx = subsizes(1)
 ny = subsizes(2)
 nz = subsizes(3)

 ldx = sizes(1)
 ldy = sizes(2)
 ldz = sizes(3)

 st_x = array_of_starts(1)
 st_y = array_of_starts(2)
 st_z = array_of_starts(3)

 ! The view starts at the first element of the submatrix.
 my_offpad = (st_x-1)*bsize_old + &
&            (st_y-1)*    (ldx*bsize_old+2*xmpio_bsize_frm) + &
&            (st_z-1)*ldy*(ldx*bsize_old+2*xmpio_bsize_frm) + &
&             xmpio_bsize_frm

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(ny,1,stride_x,column_type,plane_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(nz,1,ldy*stride_x,plane_type,new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! Commit the datatype
 call MPI_TYPE_COMMIT(new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! Free memory
 call MPI_TYPE_FREE(plane_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine xmpio_create_fsubarray_3D
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fsubarray_4D
!! NAME
!!  xmpio_create_fsubarray_4D
!!
!! FUNCTION
!!  Return a MPI type that can be used to (read|write) a 2D matrix of elements of type old_type stored in a Fortran file.
!!
!! INPUTS
!!  sizes(4)=number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  subsizes(4)=number of elements of type old_type in each dimension of the subarray (array of positive integers)
!!  array_of_starts(4)=starting coordinates of the subarray in each dimension (array of nonnegative integers >=1, <=sizes)
!!  old_type=Old MPI type.
!!
!! OUTPUT
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker individuating the beginning of the matrix. (lets call it "base").
!!    Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  new_type=New MPI type.
!!  mpierr= MPI error code
!!
!! PARENTS
!!      m_mpiotk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_fsubarray_4D(sizes,subsizes,array_of_starts,old_type,new_type,my_offpad,mpierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: mpierr,new_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: sizes(4),subsizes(4),array_of_starts(4)

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,bsize_old,nx,ny,nz,na
 integer :: column_type,plane_type,ldx,ldy,ldz,lda,vol_type
 integer(XMPI_OFFSET_KIND) :: st_x,st_y,st_z,st_a
 integer(MPI_ADDRESS_KIND) :: stride_x

!************************************************************************

 bsize_frm = xmpio_bsize_frm    ! Byte size of the Fortran record marker.

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)
 ABI_HANDLE_MPIERR(mpierr)
 !
 ! Number of columns and rows of the submatrix.
 nx = subsizes(1)
 ny = subsizes(2)
 nz = subsizes(3)
 na = subsizes(4)

 ldx = sizes(1)
 ldy = sizes(2)
 ldz = sizes(3)
 lda = sizes(4)

 st_x = array_of_starts(1)
 st_y = array_of_starts(2)
 st_z = array_of_starts(3)
 st_a = array_of_starts(4)

 ! The view starts at the first element of the submatrix.
 my_offpad = (st_x-1)*bsize_old + &
&            (st_y-1)*        (ldx*bsize_old+2*xmpio_bsize_frm) + &
&            (st_z-1)*ldy*    (ldx*bsize_old+2*xmpio_bsize_frm) + &
&            (st_a-1)*lda*ldy*(ldx*bsize_old+2*xmpio_bsize_frm) + &
&             xmpio_bsize_frm

 ! Byte size of the Fortran record + the two markers.
 stride_x = ldx*bsize_old + 2*bsize_frm

 call MPI_Type_contiguous(nx,old_type,column_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(ny,1,stride_x,column_type,plane_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(nz,1,ldy*stride_x,plane_type,vol_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_Type_create_hvector(na,1,ldz*ldy*stride_x,vol_type,new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! Commit the datatype
 call MPI_TYPE_COMMIT(new_type,mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 ! Free memory
 call MPI_TYPE_FREE(column_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(plane_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

 call MPI_TYPE_FREE(vol_type, mpierr)
 ABI_HANDLE_MPIERR(mpierr)

end subroutine xmpio_create_fsubarray_4D
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_check_frmarkers
!! NAME
!!  xmpio_check_frmarkers
!!
!! FUNCTION
!!  Check a set of Fortran record markers starting at a given offset using MPI-IO.
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  sc_mode=Option for individual or collective reading.
!!  nfrec=Number of Fortran records to be checked.
!!  bsize_frecord(nfrec)=Byte size of the Fortran records (markers are NOT included)
!!    These values will be compared with the markers reported in the file.
!!
!! OUTPUT
!!  ierr=A non-zero error code signals failure.
!!
!! PARENTS
!!      m_bse_io,m_exc_itdiago,m_slk,m_wfk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_check_frmarkers(fh,offset,sc_mode,nfrec,bsize_frecord,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nfrec,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: ierr
!arrays
 integer(XMPI_OFFSET_KIND),intent(in) :: bsize_frecord(nfrec)

!Local variables-------------------------------
!scalars
 integer :: nb,irec,frmarkers_type,jj,bsize_frm,mpi_type_frm,mpierr,myfh
 integer(XMPI_OFFSET_KIND) :: displ
!arrays
 integer(kind=int16),allocatable :: bufdelim2(:)
 integer(kind=int32),allocatable :: bufdelim4(:)
 integer(kind=int64),allocatable :: bufdelim8(:)
#ifdef HAVE_FC_INT_QUAD
 integer*16,allocatable :: bufdelim16(:)
#endif
!integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND),allocatable :: delim_record(:)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 ierr=0

 bsize_frm    = xmpio_bsize_frm     ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm  ! MPI type of the record marker.
 !
 ! Define the view for the file.
 nb=2*nfrec
 ABI_MALLOC(block_length,(nb+2))
 ABI_MALLOC(block_displ,(nb+2))
 ABI_MALLOC(block_type,(nb+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 jj=2; displ=0
 do irec=1,nfrec
   block_type (jj:jj+1) =mpi_type_frm
   block_length(jj:jj+1)=1
   block_displ(jj  )     = displ
   block_displ(jj+1)     = bsize_frm + displ + bsize_frecord(irec)
   jj=jj+2
   displ = displ + bsize_frecord(irec) + 2*bsize_frm ! Move to the beginning of the next column.
   if (xmpio_max_address(displ)) ierr=-1  ! Check for wraparound.
 end do

 block_length(nb+2)=1
 block_displ (nb+2)=displ
 block_type  (nb+2)=MPI_UB

 call xmpio_type_struct(nb+2,block_length,block_displ,block_type,frmarkers_type,mpierr)
 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 call MPI_TYPE_COMMIT(frmarkers_type,mpierr)
 call MPI_FILE_SET_VIEW(myfh,offset,MPI_BYTE,frmarkers_type,"native",MPI_INFO_NULL,mpierr)

 jj=1
 ABI_MALLOC(delim_record,(nb))
 do irec=1,nfrec
   delim_record(jj:jj+1)=bsize_frecord(irec)
   jj=jj+2
 end do

 ! Read markers according to the MPI type of the Fortran marker.
 SELECT CASE (bsize_frm)

 CASE (4)
   ABI_MALLOC(bufdelim4,(nb))
   if (sc_mode==xmpio_single) then
     call MPI_FILE_READ    (myfh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   if (ANY(bufdelim4/=delim_record)) ierr=1
   !if (ierr==1) then
   !  do irec=1,2*nfrec
   !    write(std_out,*)"irec, bufdelim4, delim_record: ",irec,bufdelim4(irec),delim_record(irec)
   !  end do
   !end if
   ABI_FREE(bufdelim4)

 CASE (8)
   ABI_MALLOC(bufdelim8,(nb))
   if (sc_mode==xmpio_single) then
     call MPI_FILE_READ    (myfh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   if (ANY(bufdelim8/=delim_record)) ierr=1
   ABI_FREE(bufdelim8)

#ifdef HAVE_FC_INT_QUAD
 CASE (16)
   ABI_MALLOC(bufdelim16,(nb))
   if (sc_mode==xmpio_single) then
     call MPI_FILE_READ    (myfh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   if (ANY(bufdelim16/=delim_record)) ierr=1
   ABI_FREE(bufdelim16)
#endif

 CASE (2)
   ABI_MALLOC(bufdelim2,(nb))
   if (sc_mode==xmpio_single) then
     call MPI_FILE_READ    (myfh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_READ_ALL(myfh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   if (ANY(bufdelim2/=delim_record)) ierr=1
   ABI_FREE(bufdelim2)

 CASE DEFAULT
   ierr=-2
 END SELECT

 ! Free memory
 call MPI_TYPE_FREE(frmarkers_type,mpierr)
 ABI_FREE(delim_record)

end subroutine xmpio_check_frmarkers
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpio_read_int
!! NAME
!!  xmpio_read_int
!!
!! FUNCTION
!!  Read the content of a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  the file pointer is modified according to the value of advance.
!!  target: integer array
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  sc_mode=
!!         xmpio_single     ==> for reading by current proc.
!!         xmpio_collective ==> for collective reading.
!!  ncount=Number of elements in the buffer
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next read will continue picking information
!!    off of the currect record.
!!
!! OUTPUT
!!  buf(ncount)=array with the values read from file
!!  fmarker=Content of the Fortran record marker.
!!  mpierr= MPI error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: file pointer used to access the Fortran marker.
!!     output: new offset updated after the reading, depending on advance.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_read_int(fh,offset,sc_mode,ncount,buf,fmarker,mpierr,advance)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,sc_mode,ncount
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer(XMPI_OFFSET_KIND),intent(out) :: fmarker
 integer,intent(out) :: mpierr
 logical,optional,intent(in) :: advance
!arrays
 integer,intent(out) :: buf(ncount)

!Local variables-------------------------------
!scalars
 integer :: myfh,bsize_frm
 integer(XMPI_OFFSET_KIND) :: my_offset
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 my_offset = offset
 bsize_frm = xmpio_bsize_frm  ! Byte size of the Fortran record marker.

 call xmpio_read_frm(myfh,my_offset,sc_mode,fmarker,mpierr,advance=.FALSE.)

 SELECT CASE (sc_mode)
 CASE (xmpio_single)
   call MPI_FILE_READ_AT(myfh, my_offset, buf, ncount, MPI_INTEGER, statux, mpierr)

 CASE (xmpio_collective)
   call MPI_FILE_READ_AT_ALL(myfh, my_offset, buf, ncount, MPI_INTEGER, statux, mpierr)

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for sc_mode: ",sc_mode
   call xmpi_abort(msg=msg)
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm  ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_read_int
!!***
#endif

!----------------------------------------------------------------------

!!****f* m_xmpi/xmpio_read_dp
!! NAME
!!  xmpio_read_dp
!!
!! FUNCTION
!!  Read the content of a single record marker in a FORTRAN file at a given offset using MPI-IO.
!!  the file pointer is modified according to the value of advance.
!!  targer: double precision real array
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  sc_mode=
!!         xmpio_single     ==> for reading by current proc.
!!         xmpio_collective ==> for collective reading.
!!  ncount=Number of elements in the buffer
!!  [advance]=By default the routine will move the file pointer to the next record.
!!    advance=.FALSE. can be used so that the next read will continue picking information
!!    off of the currect record.
!!
!! OUTPUT
!!  buf(ncount)=array with the values read from file
!!  fmarker=Content of the Fortran record marker.
!!  mpierr= MPI error code
!!
!! SIDE EFFECTS
!!  offset=
!!     input: file pointer used to access the Fortran marker.
!!     output: new offset updated after the reading, depending on advance.
!!
!! PARENTS
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_read_dp(fh,offset,sc_mode,ncount,buf,fmarker,mpierr,advance)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,sc_mode,ncount
 integer(XMPI_OFFSET_KIND),intent(inout) :: offset
 integer(XMPI_OFFSET_KIND),intent(out) :: fmarker
 integer,intent(out) :: mpierr
 logical,optional,intent(in) :: advance
!arrays
 real(dp),intent(out) :: buf(ncount)

!Local variables-------------------------------
!scalars
 integer :: bsize_frm,myfh
 integer(XMPI_OFFSET_KIND) :: my_offset
 character(len=500) :: msg
!arrays
 integer :: statux(MPI_STATUS_SIZE)

!************************************************************************

 ! Workaround for XLF
 myfh = fh

 my_offset = offset
 bsize_frm = xmpio_bsize_frm  ! Byte size of the Fortran record marker.

 call xmpio_read_frm(myfh,my_offset,sc_mode,fmarker,mpierr,advance=.FALSE.)

 SELECT CASE (sc_mode)
 CASE (xmpio_single)
   call MPI_FILE_READ_AT(myfh, my_offset, buf, ncount, MPI_DOUBLE_PRECISION, statux, mpierr)

 CASE (xmpio_collective)
   call MPI_FILE_READ_AT_ALL(myfh, my_offset, buf, ncount, MPI_DOUBLE_PRECISION, statux, mpierr)

 CASE DEFAULT
   write(msg,"(a,i0)")" Wrong value for sc_mode: ",sc_mode
   call xmpi_abort(msg=msg)
 END SELECT

 if (PRESENT(advance)) then
   if (advance) then
     offset = offset + fmarker + 2*bsize_frm ! Move the file pointer to the next record.
   else
     offset = offset + bsize_frm  ! Move the pointer after the marker.
   end if
 else
   offset = offset + fmarker + 2*bsize_frm
 end if

end subroutine xmpio_read_dp
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_max_address
!! NAME
!!  xmpio_max_address
!!
!! FUNCTION
!!  Returns .TRUE. if offset cannot be stored in a Fortran integer of kind XMPI_ADDRESS_KIND.
!!
!! PARENTS
!!
!! SOURCE

#ifdef HAVE_MPI_IO

function xmpio_max_address(offset)

!Arguments ------------------------------------
!scalars
 logical :: xmpio_max_address
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
!arrays

!Local variables-------------------------------
!scalars
 integer(XMPI_ADDRESS_KIND) :: address
 integer(XMPI_OFFSET_KIND),parameter :: max_address=HUGE(address)-100

!************************************************************************

 xmpio_max_address = (offset >= max_address)

end function xmpio_max_address
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_write_frmarkers
!! NAME
!!  xmpio_write_frmarkers
!!
!! FUNCTION
!!  Write a set of Fortran record markers starting at a given offset using MPI-IO.
!!
!! INPUTS
!!  fh=MPI-IO file handler.
!!  offset=MPI-IO file pointer
!!  sc_mode=Option for individual or collective reading.
!!  nfrec=Number of Fortran records to be written.
!!  bsize_frecord(nfrec)=Byte size of the Fortran records to be written (markers are NOT included in the size)
!!
!! OUTPUT
!!  ierr=A non-zero error code signals failure.
!!
!! PARENTS
!!      m_exc_build,m_exc_itdiago,m_ioarr,m_slk,m_wfk
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_write_frmarkers(fh,offset,sc_mode,nfrec,bsize_frecord,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: fh,nfrec,sc_mode
 integer(XMPI_OFFSET_KIND),intent(in) :: offset
 integer,intent(out) :: ierr
!arrays
 integer(XMPI_OFFSET_KIND),intent(in) :: bsize_frecord(nfrec)

!Local variables-------------------------------
!scalars
 integer :: nb,irec,frmarkers_type,jj,bsize_frm,mpi_type_frm,mpierr,myfh
 integer(XMPI_OFFSET_KIND) :: displ
!integer(XMPI_OFFSET_KIND) :: my_offset
!character(len=500) :: msg
!arrays
 integer(kind=int16),allocatable :: bufdelim2(:)
 integer(kind=int32),allocatable :: bufdelim4(:)
 integer(kind=int64),allocatable :: bufdelim8(:)
#ifdef HAVE_FC_INT_QUAD
 integer*16,allocatable :: bufdelim16(:)
#endif
!integer :: statux(MPI_STATUS_SIZE)
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND),allocatable :: delim_record(:)

!************************************************************************

 ! Workaround for XLF
 myfh = fh; ierr=0

 !my_offset = offset
 !do irec=1,nfrec
 !  call xmpio_write_frm(myfh,my_offset,sc_mode,bsize_frecord(irec),mpierr)
 !end do
 !return

 ! FIXME: This is buggy
 bsize_frm    = xmpio_bsize_frm     ! Byte size of the Fortran record marker.
 mpi_type_frm = xmpio_mpi_type_frm  ! MPI type of the record marker.

 ! Define the view for the file
 nb=2*nfrec
 ABI_MALLOC(block_length,(nb+2))
 ABI_MALLOC(block_displ,(nb+2))
 ABI_MALLOC(block_type,(nb+2))
 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 jj=2; displ=0
 do irec=1,nfrec
   block_type (jj:jj+1)  = mpi_type_frm
   block_length(jj:jj+1) = 1
   block_displ(jj  )     = displ
   block_displ(jj+1)     = displ + bsize_frm + bsize_frecord(irec)
   jj=jj+2
   displ = displ + bsize_frecord(irec) + 2*bsize_frm ! Move to the beginning of the next column.
   if (xmpio_max_address(displ)) then ! Check for wraparound.
      ierr = -1; return
   end if
 end do

 block_length(nb+2) = 1
 block_displ (nb+2) = displ
 block_type  (nb+2) = MPI_UB

 call xmpio_type_struct(nb+2,block_length,block_displ,block_type,frmarkers_type,mpierr)

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 call MPI_TYPE_COMMIT(frmarkers_type,mpierr)
 call MPI_FILE_SET_VIEW(myfh,offset,MPI_BYTE,frmarkers_type,"native",MPI_INFO_NULL,mpierr)

 jj=1
 ABI_MALLOC(delim_record,(nb))
 do irec=1,nfrec
   delim_record(jj:jj+1)=bsize_frecord(irec)
   jj=jj+2
 end do

 ! Write all markers according to the MPI type of the Fortran marker.
 SELECT CASE (bsize_frm)

 CASE (4)
   ABI_MALLOC(bufdelim4,(nb))
   bufdelim4=delim_record
   if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE    (myfh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,bufdelim4,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   ABI_FREE(bufdelim4)

 CASE (8)
   ABI_MALLOC(bufdelim8,(nb))
   bufdelim8=delim_record
   if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE    (myfh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,bufdelim8,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   ABI_FREE(bufdelim8)

#ifdef HAVE_FC_INT_QUAD
 CASE (16)
   ABI_MALLOC(bufdelim16,(nb))
   bufdelim16=delim_record
   if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE    (myfh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,bufdelim16,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   ABI_FREE(bufdelim16)
#endif

 CASE (2)
   ABI_MALLOC(bufdelim2,(nb))
   bufdelim2=delim_record
   if (sc_mode==xmpio_single) then
     call MPI_FILE_WRITE    (myfh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else if (sc_mode==xmpio_collective) then
     call MPI_FILE_WRITE_ALL(myfh,bufdelim2,2*nfrec,mpi_type_frm,MPI_STATUS_IGNORE,mpierr)
   else
     ierr=2
   end if
   ABI_FREE(bufdelim2)

 CASE DEFAULT
   ierr=-2
 END SELECT

 ! Free memory
 call MPI_TYPE_FREE(frmarkers_type,mpierr)
 ABI_FREE(delim_record)

end subroutine xmpio_write_frmarkers
#endif
!!***

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_fherm_packed
!! NAME
!!  xmpio_create_fherm_packed
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to (read|write) with MPI-IO the columns of an
!!  Hermitian matrix whose upper triangle is written on a Fortran binary file.
!!  Note that the view assumes that the file pointer used to create the MPI-IO view
!!  points to the first element of the first column. In other words,the first Fortran record marker
!!  (if any) is not taken into account in the calculation of the displacements.
!!
!! INPUTS
!!  array_of_starts(2)=starting coordinates in the global Hermitian matrix
!!     (array of positive integers with jj>=ii, Fortran convention)
!!  array_of_ends(2)=final coordinates in the global Hermitian matrix
!!     (array of positive integers, jj>=ii, Fortran convention)
!!  is_fortran_file=.FALSE. is C stream is used. .TRUE. for writing Fortran binary files.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  my_offset=Offset relative to the beginning of the matrix in the file.
!!  hmat_type=New MPI type.
!!  offset_err= error code
!!
!! NOTES
!!  The matrix on file is written in the following FORTRAN format (let us assume a 3x3 matrix for simplicity)
!!
!!    m (1,1)             m
!!    m (1,2) (2,2)       m
!!    m (1,3) (2,3) (3,3) m
!!
!!  each Fortran record stores a column of the packed Hermitian matrix, "m" denotes the Fortran
!!  record marker that introduces holes in the MPI-IO file view.
!!  To read the columns from (1,2) up to (2,2) one should use array_of_starts=(1,2) and array_of_ends=(2,2).
!!  The MPI-IO file view should be created by moving the file pointer so that it points to the elements (1,2).
!!
!!  File views for C-streams is not optimal since one can use a single slice of contigous data.
!!
!! PARENTS
!!      m_exc_build
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_fherm_packed(array_of_starts,array_of_ends,is_fortran_file,my_offset,old_type,hmat_type,offset_err)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: offset_err,hmat_type
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offset
 logical,intent(in) :: is_fortran_file
!arrays
 integer,intent(in) :: array_of_starts(2),array_of_ends(2)

!Local variables-------------------------------
!scalars
 integer :: nrow,my_ncol,ii,bsize_old,col,jj_glob,bsize_frm,prev_col,mpierr
 integer(XMPI_OFFSET_KIND) :: col_displ
!arrays
 integer,allocatable :: col_type(:),block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

 offset_err=0

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)

 bsize_frm=0; if (is_fortran_file) bsize_frm = xmpio_bsize_frm

 my_ncol = array_of_ends(2) - array_of_starts(2) + 1
 !
 ! Calculate my offset relative to the beginning of the matrix in the file.
 prev_col = array_of_starts(2)-1
 my_offset = (prev_col*(prev_col+1)/2)*bsize_old + (array_of_starts(1)-1)*bsize_old + 2*prev_col*bsize_frm + bsize_frm
 !
 ! col_type(col) describes the col-th column of the packed matrix.
 ! block_displ(col+1) stores its displacement taking into account the Fortran marker.
 ABI_MALLOC(col_type,(my_ncol))
 ABI_MALLOC(block_displ,(my_ncol+2))

 if (my_ncol>1) then
   col_displ=0
   do col=1,my_ncol
    jj_glob = (col-1) + array_of_starts(2)
    nrow = jj_glob
    if (jj_glob==array_of_starts(2)) nrow = jj_glob - array_of_starts(1) + 1 ! First column treated by me.
    if (jj_glob==array_of_ends(2))   nrow = array_of_ends(1)                 ! Last column treated by me.
    call MPI_Type_contiguous(nrow,old_type,col_type(col),mpierr)
    !
    if (xmpio_max_address(col_displ)) offset_err=1  ! Test for wraparounds
    block_displ(col+1) = col_displ
    col_displ = col_displ + nrow * bsize_old + 2 * bsize_frm  ! Move to the next column.
   end do

 else if (my_ncol==1) then  ! The case of a single column is treated separately.
    block_displ(2) = 0
    nrow = array_of_ends(1) - array_of_starts(1) + 1
    call MPI_Type_contiguous(nrow,old_type,col_type(2),mpierr)
    col_displ= nrow*bsize_old
    if (xmpio_max_address(col_displ)) offset_err=1  ! Test for wraparounds
 else
   call xmpi_abort(msg="my_ncol cannot be negative!")
 end if

 ABI_MALLOC(block_length,(my_ncol+2))
 ABI_MALLOC(block_type,(my_ncol+2))

 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB

 do ii=2,my_ncol+1
   block_length(ii)=1
   block_type(ii)  =col_type(ii-1)
   !write(std_out,*)" ii-1, depl, length, type: ",ii-1,block_displ(ii),block_length(ii),block_type(ii)
 end do

 block_length(my_ncol+2)= 1
 block_displ (my_ncol+2)= col_displ
 block_type  (my_ncol+2)= MPI_UB

 call xmpio_type_struct(my_ncol+2,block_length,block_displ,block_type,hmat_type,mpierr)

 call MPI_TYPE_COMMIT(hmat_type,mpierr)

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

 do col=1,my_ncol
   call MPI_TYPE_FREE(col_type(col),mpierr)
 end do

 ABI_FREE(col_type)

end subroutine xmpio_create_fherm_packed
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_coldistr_from_fpacked
!! NAME
!!  xmpio_create_coldistr_from_fpacked
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to MPI-IO (read|write) the columns of an
!!  (Hermitian|Symmetric) matrix whose upper triangle is written on a Fortran binary file.
!!  Note that the view assumes that the file pointer used to instanciate the MPI-IO view
!!  points to the first element of the first column. In other words,the first Fortran record marker
!!  (if any) is not taken into account in the calculation of the displacements.
!!
!! INPUTS
!!  sizes(2)=Number of elements of type old_type in each dimension of the full array (array of positive integers)
!!  my_cols(2)=initial and final column to (read|write). Array of positive integers, Fortran convention.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  new_type=New MPI type that can be used to instanciate the MPI-IO view for the Fortran file.
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker (lets call it "base"). Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  offset_err=Error code. A non-zero returned value signals that the global matrix is tool large
!!    for a single MPI-IO access (see notes below).
!!
!! NOTES
!!  1) The matrix on file is written in the following FORTRAN format (let us assume a 3x3 matrix for simplicity)
!!
!!      m (1,1)             m
!!      m (1,2) (2,2)       m
!!      m (1,3) (2,3) (3,3) m
!!
!!     each Fortran record stores a column of the packed matrix, "m" denotes the Fortran
!!     record marker that introduces holes in the file view.
!!
!!  2) With (signed) Fortran integers, the maximum size of the file that
!!     that can be read in one-shot is around 2Gb when etype is set to byte.
!!     Using a larger etype might create portability problems (real data on machines using
!!     integer*16 for the marker) since etype must be a multiple of the Fortran record marker
!!     Due to the above reason, block_displ is given in bytes but it has to be defined as Fortran
!!     integer. If the displacement cannot be stored in a Fortran integer, the routine returns
!!     offset_err=1 so that the caller will know that several MPI-IO reads are nedded to
!!     read the file.
!!
!! PARENTS
!!      m_bse_io
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_coldistr_from_fpacked(sizes,my_cols,old_type,new_type,my_offpad,offset_err)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: new_type,offset_err
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: sizes(2),my_cols(2)

!Local variables-------------------------------
!scalars
 integer :: my_ncol,bsize_old,my_col
 integer :: my_nels,my_el,row_glob,ii_hpk,jj_hpk,col_glob,bsize_frm,mpierr
 integer(XMPI_OFFSET_KIND) :: my_offset,ijp_glob
 !character(len=500) :: msg
!arrays
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)

!************************************************************************

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)

 ! my number of columns and total numer of elements to be read.
 my_ncol = my_cols(2) - my_cols(1) + 1
 my_nels = my_ncol*sizes(1)
 !
 ! block_displ(el+1) stores the displacement of the local element el taking into account the Fortran marker.
 ABI_MALLOC(block_displ,(my_nels+2))
 ABI_MALLOC(block_length,(my_nels+2))
 ABI_MALLOC(block_type,(my_nels+2))

 block_length(1)=1
 block_displ (1)=0
 block_type  (1)=MPI_LB
 !
 ! * the view assumes that the file pointer used to instanciate the MPI-IO view
 !   points to the first element of the first column. In other words,the first Fortran record marker
 !   is not taken into account in the calculation of the displacements.
 my_offpad=xmpio_bsize_frm

 ! * Some matrix elements are read twice. This part has to be tested.
 offset_err=0; my_el=0
 do my_col=1,my_ncol
   col_glob = (my_col-1) + my_cols(1)
   do row_glob=1,sizes(1)
     if (col_glob>=row_glob) then
       ii_hpk = row_glob
       jj_hpk = col_glob
       ijp_glob = row_glob + col_glob*(col_glob-1)/2  ! Index for packed form
     else ! Exchange the indeces as (jj,ii) will be read.
       ii_hpk = col_glob
       jj_hpk = row_glob
       ijp_glob = col_glob + row_glob*(row_glob-1)/2  ! Index for packed form
     end if
     my_el = my_el+1
     my_offset = (ijp_glob-1)* bsize_old + (jj_hpk-1)*2*bsize_frm
     if (xmpio_max_address(my_offset)) offset_err=1   ! Check for wraparounds.
     block_displ (my_el+1)=my_offset
     block_length(my_el+1)=1
     block_type  (my_el+1)=old_type
     !write(std_out,*)" my_el, displ: ",my_el,block_displ(my_el+1)
   end do
 end do

 block_length(my_nels+2)=1
 block_displ (my_nels+2)=my_offset
 block_type  (my_nels+2)=MPI_UB

 call xmpio_type_struct(my_nels+2,block_length,block_displ,block_type,new_type,mpierr)

 call MPI_TYPE_COMMIT(new_type,mpierr)

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

end subroutine xmpio_create_coldistr_from_fpacked
!!***
#endif

!------------------------------------------------------------------------------------

!!****f* m_xmpi/xmpio_create_coldistr_from_fp3blocks
!! NAME
!!  xmpio_create_coldistr_from_fp3blocks
!!
!! FUNCTION
!!  Returns an MPI datatype that can be used to MPI-IO (read|write) the columns of a
!!  matrix of the form  M = (S1    F3)
!!                          (F3^H  S2)
!!  where S1 and S2 are square (symmetric|Hermitian) matrices whose upper triangle is stored on file
!!  while F3 is a generic matrix (not necessarily square) stored in full mode.
!!  The Fortran file contains the blocks in the following order.
!!      upper(S1)
!!      upper(S2)
!!      F3
!! INPUTS
!!  sizes(2)=Number of elements of type old_type in each dimension of the full array M (array of positive integers)
!!  my_cols(2)=initial and final column to (read|write). Array of positive integers, Fortran convention.
!!  block_sizes(2,3)=The sizes of S1, S2, F.
!!  old_type=MPI datatype of the elements of the matrix.
!!
!! OUTPUT
!!  new_type=New MPI type that can be used to instanciate the MPI-IO view for the Fortran file.
!!  my_offpad=Offset to be added to the file pointer giving the position of the first Fortran record
!!    marker (lets call it "base"). Each node should (read|write) using my_offset = base + my_offpad.
!!    my_offpad is used so that one can safely change the way the fileview is generated (for example
!!    to make it more efficient) without having to change the client code.
!!  offset_err=Error code. A non-zero returned value signals that the global matrix is tool large
!!    for a single MPI-IO access (see notes below).
!!
!! NOTES
!!  1) block_displ is given in bytes due to the presence of the marker.
!!     If the displacement of an element is too large, the routine returns
!!     offset_err=1 so that the caller knows that several MPI-IO reads are required to (read| write) the file.
!!
!! PARENTS
!!      m_bse_io
!!
!! CHILDREN
!!      xmpi_comm_free
!!
!! SOURCE

#ifdef HAVE_MPI_IO

subroutine xmpio_create_coldistr_from_fp3blocks(sizes,block_sizes,my_cols,old_type,new_type,my_offpad,offset_err)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: old_type
 integer,intent(out) :: new_type,offset_err
 integer(XMPI_OFFSET_KIND),intent(out) :: my_offpad
!arrays
 integer,intent(in) :: sizes(2),my_cols(2),block_sizes(2,3)

!Local variables-------------------------------
!scalars
 integer :: my_ncol,bsize_old,my_col,which_block,uplo,swap
 integer :: my_nels,my_el,row_glob,ii_hpk,jj_hpk,ii,jj
 integer :: col_glob,bsize_frm,mpierr,row_shift,col_shift,n1,n2
 integer(XMPI_OFFSET_KIND) :: my_offset,ijp,bsize_tot,max_displ,min_displ
 integer(XMPI_ADDRESS_KIND) :: address
!arrays
 integer,allocatable :: block_length(:),block_type(:)
 integer(XMPI_ADDRESS_KIND),allocatable :: block_displ(:)
 integer(XMPI_OFFSET_KIND) :: bsize_mat(2)

!************************************************************************

 if (sizes(1) /= SUM(block_sizes(1,1:2)) .or. &
     sizes(2) /= SUM(block_sizes(2,1:2)) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: Inconsistency between block_sizes ans sizes "
   call xmpi_abort()
 end if

 if (block_sizes(1,1) /= block_sizes(2,1) .or.&
     block_sizes(1,2) /= block_sizes(2,2) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: first two blocks must be square"
   call xmpi_abort()
 end if

 if (block_sizes(2,3) /= block_sizes(2,2) .or.&
     block_sizes(1,3) /= block_sizes(1,1) ) then
   write(std_out,*)" xmpio_create_coldistr_from_fp3blocks: Full matrix must be square"
   call xmpi_abort()
 end if

 write(std_out,*)" xmpio_create_coldistr_from_fp3blocks is still under testing"
 !call xmpi_abort()

 ! Byte size of the Fortran record marker.
 bsize_frm = xmpio_bsize_frm

 ! Byte size of old_type.
 call MPI_TYPE_SIZE(old_type,bsize_old,mpierr)

 ! my number of columns and total numer of elements to be read.
 my_ncol = my_cols(2) - my_cols(1) + 1
 my_nels = sizes(1)*my_ncol
 !
 ! block_displ(el+1) stores the displacement of the local element el taking into account the Fortran marker.
 ABI_MALLOC(block_displ,(my_nels+2))
 ABI_MALLOC(block_length,(my_nels+2))
 ABI_MALLOC(block_type,(my_nels+2))
 !
 ! * the view assumes that the file pointer used to instanciate the MPI-IO view
 !   points to the first element of the first column. In other words,the first Fortran record marker
 !   is not taken into account in the calculation of the displacements.
 my_offpad=xmpio_bsize_frm
 !
 ! Byte size of the first two blocks including the markers.
 n1=block_sizes(1,1)
 bsize_mat(1) = (n1*(n1+1)/2)*bsize_old + 2*n1*bsize_frm

 n2=block_sizes(1,2)
 bsize_mat(2) = (n2*(n2+1)/2)*bsize_old + 2*n2*bsize_frm

 bsize_tot=SUM(bsize_mat) +  PRODUCT(block_sizes(:,3))*bsize_old + block_sizes(2,3)*2*bsize_frm - bsize_frm
 write(std_out,*)"bsize_mat",bsize_mat,"bsize_tot",bsize_tot
 !
 ! * Some matrix elements are read twice. This part has to be tested.
 offset_err=0; my_el=0; max_displ=0; min_displ=HUGE(address)
 do my_col=1,my_ncol
   col_glob = (my_col-1) + my_cols(1)
   do row_glob=1,sizes(1)
     !
     which_block=3
     if (row_glob<=block_sizes(1,1).and.col_glob<=block_sizes(2,1)) which_block=1
     if (row_glob >block_sizes(1,1).and.col_glob >block_sizes(2,1)) which_block=2

     if ( ANY(which_block == (/1,2/)) ) then ! S1 or S2
       !
       row_shift=(which_block-1)*block_sizes(1,1)
       col_shift=(which_block-1)*block_sizes(2,1)

       ii_hpk = row_glob - row_shift
       jj_hpk = col_glob - col_shift
       if (jj_hpk<ii_hpk) then ! Exchange the indeces so that the symmetric is read.
         swap   = jj_hpk
         jj_hpk = ii_hpk
         ii_hpk = swap
       end if
       ijp = ii_hpk + jj_hpk*(jj_hpk-1)/2  ! Index for packed form
       my_offset = (ijp-1)*bsize_old + (jj_hpk-1)*2*bsize_frm
       if (which_block==2) my_offset=my_offset+bsize_mat(1)    ! Shift the offset to account for S1.
       !my_offset=4
       !
     else
       ! The element belongs either to F3 of F3^H.
       ! Now find whether it is the upper or the lower block since only F3 is stored on file.
       uplo=1; if (row_glob>block_sizes(1,1)) uplo=2

       if (uplo==1) then
         row_shift=0
         col_shift=block_sizes(2,1)
       else
         row_shift=block_sizes(1,1)
         col_shift=0
       end if
       ii = row_glob - row_shift
       jj = col_glob - col_shift

       if (uplo==2) then ! Exchange the indeces since the symmetric element will be read.
         swap=jj
         jj  =ii
         ii  =swap
       end if

       my_offset = (ii-1)*bsize_old + (jj-1)*block_sizes(1,3)*bsize_old + (jj-1)*2*bsize_frm
       my_offset = my_offset + SUM(bsize_mat)
       !if (uplo==1) my_offset=my_offset + bsize_mat(1)
       !my_offset=0
       !if (ii==1.and.jj==1) write(std_out,*)" (1,1) offset = ",my_offset
       !if (ii==block_sizes(1,3).and.jj==block_sizes(2,3)) write(std_out,*)" (n,n) offset =", my_offset
       if (my_offset>=bsize_tot-1*bsize_old) then
         write(std_out,*)"WARNING (my_offset>bsize_tot-bsize_old),",ii,jj,my_offset,bsize_tot
       end if
     end if

     if (xmpio_max_address(my_offset)) offset_err=1   ! Check for wraparounds.
     my_el = my_el+1
     block_displ (my_el+1)=my_offset
     block_length(my_el+1)=1
     block_type  (my_el+1)=old_type
     max_displ = MAX(max_displ,my_offset)
     min_displ = MIN(min_displ,my_offset)
     !if (which_block==3) write(std_out,*)" my_el, which, displ: ",my_el,which_block,block_displ(my_el+1)
   end do
 end do

 write(std_out,*)" MAX displ = ",max_displ," my_nels = ",my_nels
 write(std_out,*)" MIN displ = ",MINVAL(block_displ(2:my_nels+1))

 !block_displ (1)=max_displ ! Do not change this value.
 !if (min_displ>0) block_displ (1)=min_displ ! Do not change this value.

 block_displ (1)=min_displ
 block_displ (1)=0
 block_length(1)=0
 block_type  (1)=MPI_LB

 block_length(my_nels+2)=0
 !block_displ (my_nels+2)=bsize_tot
 block_displ (my_nels+2)=max_displ
 block_type  (my_nels+2)=MPI_UB

 call xmpio_type_struct(my_nels+2,block_length,block_displ,block_type,new_type,mpierr)
 !call xmpio_type_struct(my_nels,block_length(2:),block_displ(2:),block_type(2:),new_type,mpierr)

 !call MPI_TYPE_CREATE_INDEXED_BLOCK(my_nels, block_length(2:), block_displ(2:), old_type, new_type, mpierr)

 call MPI_TYPE_COMMIT(new_type,mpierr)

 ABI_FREE(block_length)
 ABI_FREE(block_displ)
 ABI_FREE(block_type)

end subroutine xmpio_create_coldistr_from_fp3blocks
!!***
#endif

 !type(xcomm_t) function from_mpi_int(comm_value) result(new)
 !  new%value = comm_value
 !  new%nproc  xmpi_comm_size(comm_value)
 !  new%me  xmpi_comm_rank(comm_value)
 !end function from_mpi_int
 !pure logical function xcomm_iam_master(self)
 !  class(xcomm_t),intent(in) :: self
 !  xcomm_iam_master = self%me == 0
 !end function xcomm_iam_master
 pure logical function xcomm_skip(self, iter)
   class(xcomm_t),intent(in) :: self
   integer,intent(in) :: iter
   xcomm_skip = mod(iter, self%nproc) /= self%me
 end function xcomm_skip
 subroutine xcomm_set_to_self(self)
   class(xcomm_t),intent(inout) :: self
   call self%free()
   self%value = xmpi_comm_self; self%me = 0; self%nproc = 1
 end subroutine xcomm_set_to_self
 subroutine xcomm_set_to_null(self)
   class(xcomm_t),intent(inout) :: self
   call self%free()
   self%value = xmpi_comm_null
 end subroutine xcomm_set_to_null
 subroutine xcomm_free(self)
   class(xcomm_t),intent(inout) :: self
   call xmpi_comm_free(self%value)
   self%me = -1; self%nproc = 0
 end subroutine xcomm_free

END MODULE m_xmpi
!!***

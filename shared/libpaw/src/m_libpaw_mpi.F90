!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_libpaw_mpi
!! NAME
!!  m_libpaw_mpi
!!
!! FUNCTION
!!  libPAW wrappers for MPI library.
!!  Provides MPI named constants or tools as well as
!!  a set of generic interfaces wrapping MPI primitives.
!!
!! COPYRIGHT
!!  Copyright (C) 2014-2020 ABINIT group (MT, MG, ...)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  This file comes directly from m_xmpi.F90 module delivered with ABINIT.
!!
!!  FOR DEVELOPPERS: in order to preserve the portability of libPAW library,
!!  please consult ~abinit/src/??_libpaw/libpaw-coding-rules.txt
!!
!! SOURCE

#include "libpaw.h"

module m_libpaw_mpi
    
 USE_DEFS

#ifdef HAVE_MPI2
 use mpi
#endif

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

!----------------------------------------------------------------------
!Some MPI constants

#ifdef HAVE_MPI
 integer,public,parameter :: xpaw_mpi_world     = MPI_COMM_WORLD
 integer,public,parameter :: xpaw_mpi_comm_self = MPI_COMM_SELF
 integer,public,parameter :: xpaw_mpi_comm_null = MPI_COMM_NULL
#else
 integer,public,parameter :: xpaw_mpi_world     = 0
 integer,public,parameter :: xpaw_mpi_comm_self = 0
 integer,public,parameter :: xpaw_mpi_comm_null = 0
#endif

!----------------------------------------------------------------------
!MPI public procedures

 public :: xpaw_mpi_abort                 ! Wrapper for MPI_ABORT
 public :: xpaw_mpi_comm_rank             ! Wrapper for MPI_COMM_RANK
 public :: xpaw_mpi_comm_size             ! Wrapper for MPI_COMM_SIZE
 public :: xpaw_mpi_barrier               ! Wrapper for MPI_BARRIER
 public :: xpaw_mpi_wait                  ! Wrapper for MPI_WAIT
 public :: xpaw_mpi_waitall               ! Wrapper for MPI_WAITALL
 public :: xpaw_mpi_iprobe                ! Wrapper for MPI_IPROBE

!----------------------------------------------------------------------
!MPI private procedures
 private :: xpaw_mpi_get_tag_ub           ! Get MPI_TAG_UB attribute

!----------------------------------------------------------------------
!MPI generic interfaces

 public :: xpaw_mpi_allgather             ! Wrapper for MPI_ALLGATHER
 interface xpaw_mpi_allgather
  module procedure xpaw_mpi_allgather_int1d
  module procedure xpaw_mpi_allgather_dp1d
  module procedure xpaw_mpi_allgather_dp2d
  module procedure xpaw_mpi_allgather_dp3d
 end interface xpaw_mpi_allgather

 public :: xpaw_mpi_allgatherv            ! Wrapper for MPI_ALLGATHERV
 interface xpaw_mpi_allgatherv
  module procedure xpaw_mpi_allgatherv_int1d
  module procedure xpaw_mpi_allgatherv_dp1d
  module procedure xpaw_mpi_allgatherv_dp2d
 end interface xpaw_mpi_allgatherv

 public :: xpaw_mpi_scatterv              ! Wrapper for MPI_SCATTERV
 interface xpaw_mpi_scatterv
  module procedure xpaw_mpi_scatterv_int1d
  module procedure xpaw_mpi_scatterv_dp1d
  module procedure xpaw_mpi_scatterv_dp2d
 end interface xpaw_mpi_scatterv

 public :: xpaw_mpi_alltoall              ! Wrapper for MPI_ALLTOALL
 interface xpaw_mpi_alltoall
  module procedure xpaw_mpi_alltoall_int1d
  module procedure xpaw_mpi_alltoall_dp1d
  module procedure xpaw_mpi_alltoall_dp2d
 end interface xpaw_mpi_alltoall

 public :: xpaw_mpi_alltoallv             ! Wrapper for MPI_ALLTOALLV
 interface xpaw_mpi_alltoallv
  module procedure xpaw_mpi_alltoallv_int1d
  module procedure xpaw_mpi_alltoallv_dp1d
  module procedure xpaw_mpi_alltoallv_dp2d
 end interface xpaw_mpi_alltoallv

 public :: xpaw_mpi_bcast                 ! Wrapper for MPI_BCAST
 interface xpaw_mpi_bcast
  module procedure xpaw_mpi_bcast_int
  module procedure xpaw_mpi_bcast_int1d
  module procedure xpaw_mpi_bcast_dp1d
  module procedure xpaw_mpi_bcast_dp2d
  module procedure xpaw_mpi_bcast_dp3d
 end interface xpaw_mpi_bcast

 public :: xpaw_mpi_gather                ! Wrapper for MPI_GATHER
 interface xpaw_mpi_gather
  module procedure xpaw_mpi_gather_int1d
  module procedure xpaw_mpi_gather_dp1d
  module procedure xpaw_mpi_gather_dp2d
 end interface xpaw_mpi_gather

 public :: xpaw_mpi_gatherv               ! Wrapper for MPI_GATHERV
 interface xpaw_mpi_gatherv
  module procedure xpaw_mpi_gatherv_int1d
  module procedure xpaw_mpi_gatherv_dp1d
  module procedure xpaw_mpi_gatherv_dp2d
 end interface xpaw_mpi_gatherv

 public :: xpaw_mpi_recv                  ! Wrapper for MPI_RECV
 interface xpaw_mpi_recv
  module procedure xpaw_mpi_recv_int1d
  module procedure xpaw_mpi_recv_dp1d
  module procedure xpaw_mpi_recv_dp2d
  module procedure xpaw_mpi_recv_dp3d
 end interface xpaw_mpi_recv

 public :: xpaw_mpi_irecv                 ! Wrapper for MPI_IRECV
 interface xpaw_mpi_irecv
  module procedure xpaw_mpi_irecv_int1d
  module procedure xpaw_mpi_irecv_dp1d
  module procedure xpaw_mpi_irecv_dp2d
 end interface xpaw_mpi_irecv

 public :: xpaw_mpi_send                  ! Wrapper for MPI_SEND
 interface xpaw_mpi_send
  module procedure xpaw_mpi_send_int1d
  module procedure xpaw_mpi_send_dp1d
  module procedure xpaw_mpi_send_dp2d
  module procedure xpaw_mpi_send_dp3d
 end interface xpaw_mpi_send

 public :: xpaw_mpi_isend                 ! Wrapper for MPI_ISEND
 interface xpaw_mpi_isend
  module procedure xpaw_mpi_isend_int1d
  module procedure xpaw_mpi_isend_dp1d
  module procedure xpaw_mpi_isend_dp2d
 end interface xpaw_mpi_isend

 public :: xpaw_mpi_exch                  ! Wrapper for MPI_SEND/MPI_RECV
 interface xpaw_mpi_exch
  module procedure xpaw_mpi_exch_int1d
  module procedure xpaw_mpi_exch_dp1d
  module procedure xpaw_mpi_exch_dp2d
  module procedure xpaw_mpi_exch_dp3d
 end interface xpaw_mpi_exch

 public :: xpaw_mpi_sum                   ! Wrapper for MPI_ALLREDUCE(SUM)
 interface xpaw_mpi_sum
  module procedure xpaw_mpi_sum_int
  module procedure xpaw_mpi_sum_int1d
  module procedure xpaw_mpi_sum_dp1d
  module procedure xpaw_mpi_sum_dp2d
  module procedure xpaw_mpi_sum_dp3d
 end interface xpaw_mpi_sum
!!***

!----------------------------------------------------------------------

CONTAINS  !===========================================================
!!***

!!****f* m_libpaw_mpi/xpaw_mpi_abort
!! NAME
!!  xpaw_mpi_abort
!!
!! FUNCTION
!!  Wrapper for MPI_ABORT
!!
!! INPUTS
!!  [comm]=communicator of tasks to abort.
!!  [mpierr]=Error code to return to invoking environment.
!!  [msg]=User message
!!  [exit_status]=optional, shell return code, default 1
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_abort(comm,mpierr,msg,exit_status)

!Arguments-------------------------
 integer,optional,intent(in) :: comm,mpierr,exit_status
 character(len=*),optional,intent(in) :: msg

!Local variables-------------------
 integer :: ilen,ierr,ierr2,my_comm,my_errorcode,my_exit_status
 logical :: testopen
#ifdef HAVE_MPI
 character(len=MPI_MAX_ERROR_STRING) :: mpi_msg_error
#endif

! *************************************************************************

 ierr=0
 my_comm = xpaw_mpi_world; if (PRESENT(comm)) my_comm = comm
 my_exit_status = 1; if (PRESENT(exit_status)) my_exit_status=exit_status

 if (PRESENT(msg)) then
   write(std_out,'(2a)')"User message: ",TRIM(msg)
 end if

 ! Close std_out and ab_out
 inquire(std_out,opened=testopen)
 if (testopen) close(std_out)
 inquire(ab_out,opened=testopen)
 if (testopen) close(ab_out)

#ifdef HAVE_MPI
 my_errorcode=MPI_ERR_UNKNOWN; if (PRESENT(mpierr)) my_errorcode=mpierr
 call MPI_ERROR_STRING(my_errorcode, mpi_msg_error, ilen, ierr2)
 call MPI_ABORT(my_comm,my_errorcode,ierr)
#endif

#if defined FC_NAG
 call exit(exit_status)
#elif defined HAVE_FC_EXIT
 call exit(exit_status)
#else
 if (exit_status== 0) stop  "0"
 if (exit_status== 1) stop  "1"
 if (exit_status==-1) stop "-1"
#endif
 stop "1"

end subroutine xpaw_mpi_abort
!!***


!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_comm_rank
!! NAME
!!  xpaw_mpi_comm_rank
!!
!! FUNCTION
!!  Wrapper for MPI_COMM_RANK
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  xpaw_mpi_comm_rank=The rank of the node inside comm
!!
!! SOURCE

function xpaw_mpi_comm_rank(comm)

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: xpaw_mpi_comm_rank

!Local variables-------------------
 integer :: mpierr

! *************************************************************************

 mpierr=0
#ifdef HAVE_MPI
 xpaw_mpi_comm_rank=-1  ! Return non-sense value if the proc does not belong to the comm
 if (comm/=xpaw_mpi_comm_null) then
   call MPI_COMM_RANK(comm,xpaw_mpi_comm_rank,mpierr)
 end if
#else
 xpaw_mpi_comm_rank=0
#endif

end function xpaw_mpi_comm_rank
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_comm_size
!! NAME
!!  xpaw_mpi_comm_size
!!
!! FUNCTION
!!  Wrapper for MPI_COMM_SIZE
!!
!! INPUTS
!!  comm=MPI communicator.
!!
!! OUTPUT
!!  xpaw_mpi_comm_size=The number of processors inside comm.
!!
!! SOURCE

function xpaw_mpi_comm_size(comm)

!Arguments-------------------------
 integer,intent(in) :: comm
 integer :: xpaw_mpi_comm_size

!Local variables-------------------------------
!scalars
 integer :: mpierr

! *************************************************************************

 mpierr=0; xpaw_mpi_comm_size=1
#ifdef HAVE_MPI
 if (comm/=xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,xpaw_mpi_comm_size,mpierr)
 end if
#endif

end function xpaw_mpi_comm_size
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_barrier
!! NAME
!!  xpaw_mpi_barrier
!!
!! FUNCTION
!!  Wrapper for MPI_BARRIER
!!
!! INPUTS
!!  comm=MPI communicator
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_barrier(comm)

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
 if (comm/=xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nprocs,ier)
   if(nprocs>1)then
     call MPI_BARRIER(comm,ier)
   end if
 end if
#endif

end subroutine xpaw_mpi_barrier
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_wait
!! NAME
!!  xpaw_mpi_wait
!!
!! FUNCTION
!!  Wrapper for MPI_WAIT
!!
!! INPUTS
!!  request= MPI request handle to wait for
!!
!! OUTPUT
!!  mpierr= status error
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_wait(request,mpierr)

!Arguments-------------------------
 integer,intent(out) :: mpierr
 integer,intent(inout) :: request

!Local variables-------------------
#ifdef HAVE_MPI
 integer :: ier,status(MPI_STATUS_SIZE)
#endif

! *************************************************************************

 mpierr = 0
#ifdef HAVE_MPI
  call MPI_WAIT(request,status,ier)
  mpierr=ier
#endif

end subroutine xpaw_mpi_wait
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_waitall
!! NAME
!!  xpaw_mpi_waitall
!!
!! FUNCTION
!!  Wrapper for MPI_WAITALL
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
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_waitall(array_of_requests,mpierr)

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
  call MPI_WAITALL(size(array_of_requests),array_of_requests,status,ier)
  mpierr=ier
#endif

end subroutine xpaw_mpi_waitall
!!***

!----------------------------------------------------------------------

!!****f* m_libpaw_mpi/xpaw_mpi_iprobe
!! NAME
!!  xpaw_mpi_iprobe
!!
!! FUNCTION
!!  Wrapper for MPI_IPROBE
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
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_iprobe(source,tag,mpicomm,flag,mpierr)

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

end subroutine xpaw_mpi_iprobe
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_allgather_int1d
!! NAME
!!  xpaw_mpi_allgather_int1d
!!
!! FUNCTION
!!  MPI_ALLGATHER for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received elements
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgather_int1d(xval,nelem,recvbuf,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHER(xval,nelem,MPI_INTEGER,recvbuf,nelem,MPI_INTEGER,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf(1:nelem)=xval(1:nelem)
 end if
#else
 recvbuf(1:nelem)=xval(1:nelem)
#endif
end subroutine xpaw_mpi_allgather_int1d
!!***

!!****f* ABINIT/xpaw_mpi_allgather_dp1d
!! NAME
!!  xpaw_mpi_allgather_dp1d
!!
!! FUNCTION
!!  MPI_ALLGATHER for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received elements
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgather_dp1d(xval,nelem,recvbuf,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout) :: recvbuf(:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf(1:nelem)=xval(1:nelem)
 end if
#else
 recvbuf(1:nelem)=xval(1:nelem)
#endif
end subroutine xpaw_mpi_allgather_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_allgather_dp2d
!! NAME
!!  xpaw_mpi_allgather_dp2d
!!
!! FUNCTION
!!  MPI_ALLGATHER for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received elements
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgather_dp2d(xval,nelem,recvbuf,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf(:,:)=xval(:,:)
 end if
#else
 recvbuf(:,:)=xval(:,:)
#endif
end subroutine xpaw_mpi_allgather_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_allgather_dp3d
!! NAME
!!  xpaw_mpi_allgather_dp3d
!!
!! FUNCTION
!!  MPI_ALLGATHER for 3D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received elements
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgather_dp3d(xval,nelem,recvbuf,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:,:)
 real(dp), intent(inout) :: recvbuf(:,:,:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHER(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,nelem,MPI_DOUBLE_PRECISION,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf(:,:,:)=xval(:,:,:)
 end if
#else
 recvbuf(:,:,:)=xval(:,:,:)
#endif
end subroutine xpaw_mpi_allgather_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_allgatherv_int1d
!! NAME
!!  xpaw_mpi_allgatherv_int1d
!!
!! FUNCTION
!!  MPI_ALLGATHERV for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgatherv_int1d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: cc,dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHERV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_allgatherv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_allgatherv_dp1d
!! NAME
!!  xpaw_mpi_allgatherv_dp1d
!!
!! FUNCTION
!!  MPI_ALLGATHERV for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgatherv_dp1d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout) :: recvbuf(:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
 integer :: cc,dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_allgatherv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_allgatherv_dp2d
!! NAME
!!  xpaw_mpi_allgatherv_dp2d
!!
!! FUNCTION
!!  MPI_ALLGATHERV for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_allgatherv_dp2d(xval,nelem,recvbuf,recvcounts,displs,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
 integer :: cc,dd,sz1
 
! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLGATHERV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   sz1=size(xval,1)
   dd=0;if (size(displs)>0) dd=displs(1)/sz1
   cc=size(xval,2);if (size(recvcounts)>0) cc=recvcounts(1)/sz1
   recvbuf(:,dd+1:dd+cc)=xval(:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_allgatherv_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_scatterv_int1d
!! NAME
!!  xpaw_mpi_scatterv_int1d
!!
!! FUNCTION
!!  MPI_SCATTERV for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcount= number of received elements
!!  displs= relative offsets for incoming data (array)
!!  sendcounts= number of sent elements (array)
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_scatterv_int1d(xval,sendcounts,displs,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: sendcounts(:),displs(:)
 integer, intent(in) :: recvcount,root,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_SCATTERV(xval,sendcounts,displs,MPI_INTEGER,recvbuf,recvcount,&
&   MPI_INTEGER,root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   recvbuf(1:recvcount)=xval(dd+1:dd+recvcount)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_scatterv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_scatterv_dp1d
!! NAME
!!  xpaw_mpi_scatterv_dp1d
!!
!! FUNCTION
!!  MPI_SCATTERV for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcount= number of received elements
!!  displs= relative offsets for incoming data (array)
!!  sendcounts= number of sent elements (array)
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_scatterv_dp1d(xval,sendcounts,displs,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout)   :: recvbuf(:)
 integer, intent(in) :: sendcounts(:),displs(:)
 integer, intent(in) :: recvcount,root,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_SCATTERV(xval,sendcounts,displs,MPI_DOUBLE_PRECISION,recvbuf,recvcount,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   recvbuf(1:recvcount)=xval(dd+1:dd+recvcount)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_scatterv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_scatterv_dp2d
!! NAME
!!  xpaw_mpi_scatterv_dp2d
!!
!! FUNCTION
!!  MPI_SCATTERV for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcount= number of received elements
!!  displs= relative offsets for incoming data (array)
!!  sendcounts= number of sent elements (array)
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_scatterv_dp2d(xval,sendcounts,displs,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout)   :: recvbuf(:,:)
 integer, intent(in) :: sendcounts(:),displs(:)
 integer, intent(in) :: recvcount,root,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: cc,dd,sz1 

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_SCATTERV(xval,sendcounts,displs,MPI_DOUBLE_PRECISION,recvbuf,recvcount,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
#endif
   sz1=size(recvbuf,1);cc=recvcount/sz1
   dd=0;if (size(displs)>0) dd=displs(1)/sz1
   recvbuf(:,1:cc)=xval(:,dd+1:dd+cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_scatterv_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_alltoall_int1d
!! NAME
!!  xpaw_mpi_alltoall_int1d
!!
!! FUNCTION
!!  MPI_ALLTOALL for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendsize= size of sent buffer
!!  recvsize= size of received buffer
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoall_int1d(xval,sendsize,recvbuf,recvsize,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: sendsize, recvsize
 integer, intent(in) :: spaceComm
 integer, intent(out) :: ier

!Local variables-------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLTOALL(xval, sendsize, MPI_INTEGER, recvbuf, &
&   recvsize, MPI_INTEGER, spaceComm, ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_alltoall_int1d
!!***

!!****f* ABINIT/xpaw_mpi_alltoall_dp1d
!! NAME
!!  xpaw_mpi_alltoall_dp1d
!!
!! FUNCTION
!!  MPI_ALLTOALL for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendsize= size of sent buffer
!!  recvsize= size of received buffer
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoall_dp1d(xval,sendsize,recvbuf,recvsize,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in)    :: xval(:)
 real(dp), intent(inout) :: recvbuf(:)
 integer, intent(in) :: sendsize, recvsize
 integer, intent(in) :: spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLTOALL(xval, sendsize, MPI_DOUBLE_PRECISION, recvbuf, &
&   recvsize, MPI_DOUBLE_PRECISION, spaceComm, ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_alltoall_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_alltoall_dp2d
!! NAME
!!  xpaw_mpi_alltoall_dp2d
!!
!! FUNCTION
!!  MPI_ALLTOALL for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendsize= size of sent buffer
!!  recvsize= size of received buffer
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoall_dp2d(xval,sendsize,recvbuf,recvsize,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in)    :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: sendsize, recvsize
 integer, intent(in) :: spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_ALLTOALL(xval, sendsize, MPI_DOUBLE_PRECISION, recvbuf, &
&   recvsize, MPI_DOUBLE_PRECISION, spaceComm, ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_alltoall_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_alltoallv_int1d
!! NAME
!!  xpaw_mpi_alltoallv_int1d
!!
!! FUNCTION
!!  MPI_ALLTOALLV for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of elements to send to each processor
!!  sdispls= displacements from which to take the outgoing data
!!  recvcnts= number of elements that can be received from each processor 
!!  rdispls= displacement at which to place the incoming data from each processor
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoallv_int1d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,comm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: sendcnts(:),sdispls(:),recvcnts(:)
 integer, intent(in) :: comm, rdispls
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: sc,sds,sdr
#if defined HAVE_MPI
 integer :: rdispls_on(size(sendcnts))
#endif

! *********************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   rdispls_on = 0
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_INTEGER,recvbuf,&
&   recvcnts,rdispls_on,MPI_INTEGER,comm,ier)
 else if (comm == xpaw_mpi_comm_self) then
#endif
   sdr=rdispls;sds=0;if (size(sdispls)>0) sds=sdispls(1)
   sc=size(xval);if (size(sendcnts)>0) sc=sendcnts(1)
   recvbuf(1:sc)=xval(sds+1:sds+sc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_alltoallv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_alltoallv_dp1d
!! NAME
!!  xpaw_mpi_alltoallv_dp1d
!!
!! FUNCTION
!!  MPI_ALLTOALLV for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of elements to send to each processor
!!  sdispls= displacements from which to take the outgoing data
!!  recvcnts= number of elements that can be received from each processor 
!!  rdispls= displacement at which to place the incoming data from each processor
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoallv_dp1d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,comm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout) :: recvbuf(:)
 integer, intent(in) :: sendcnts(:),sdispls(:),recvcnts(:),rdispls(:)
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: sc,sds,sdr

! *********************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&   recvcnts,rdispls,MPI_DOUBLE_PRECISION,comm,ier)
 else if (comm == MPI_COMM_SELF) then
#endif
   sds=0;if (size(sdispls)>0) sds=sdispls(1)
   sdr=0;if (size(rdispls)>0) sdr=rdispls(1)
   sc=size(xval);if (size(sendcnts)>0) sc=sendcnts(1)
   recvbuf(sdr+1:sdr+sc)=xval(sds+1:sds+sc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_alltoallv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_alltoallv_dp2d
!! NAME
!!  xpaw_mpi_alltoallv_dp2d
!!
!! FUNCTION
!!  MPI_ALLTOALLV for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcnts= number of elements to send to each processor
!!  sdispls= displacements from which to take the outgoing data
!!  recvcnts= number of elements that can be received from each processor 
!!  rdispls= displacement at which to place the incoming data from each processor
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_alltoallv_dp2d(xval,sendcnts,sdispls,recvbuf,recvcnts,rdispls,comm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: sendcnts(:),sdispls(:),rdispls(:),recvcnts(:)
 integer,intent(in) :: comm
 integer,intent(out) :: ier

!Local variables-------------------
 integer :: sc,sds,sdr,sz1
 
! *********************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_ALLTOALLV(xval,sendcnts,sdispls,MPI_DOUBLE_PRECISION,recvbuf,&
&   recvcnts,rdispls,MPI_DOUBLE_PRECISION,comm,ier)
 else if (comm == xpaw_mpi_comm_self) then
#endif
   sz1=size(xval,1)
   sds=0;if (size(sdispls)>0) sds=sdispls(1)/sz1
   sdr=0;if (size(rdispls)>0) sdr=rdispls(1)/sz1
   sc=size(xval,2);if (size(sendcnts)>0) sc=sendcnts(1)/sz1
   recvbuf(:,sdr+1:sdr+sc)=xval(:,sds+1:sds+sc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_alltoallv_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_bcast_int
!! NAME
!!  xpaw_mpi_bcast_int
!!
!! FUNCTION
!!  MPI_BCAST for integers
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_bcast_int(xval,master,spaceComm,ier)

!Arguments ------------------------------------
 integer, intent(inout) :: xval
 integer, intent(in) :: spaceComm,master
 integer, intent(out) :: ier

!Local variables-------------------------------

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_BCAST(xval,1,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_bcast_int
!!***

!!****f* ABINIT/xpaw_mpi_bcast_int1d
!! NAME
!!  xpaw_mpi_bcast_int1d
!!
!! FUNCTION
!!  MPI_BCAST for 1D integer arrays
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_bcast_int1d(xval,master,spaceComm,ier)

!Arguments ------------------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: spaceComm,master
 integer, intent(out) :: ier

!Local variables-------------------------------
 integer :: n

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n=size(xval)
   call MPI_BCAST(xval,n,MPI_INTEGER,master,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_bcast_int1d
!!***

!!****f* ABINIT/xpaw_mpi_bcast_dp1d
!! NAME
!!  xpaw_mpi_bcast_dp1d
!!
!! FUNCTION
!!  MPI_BCAST for 1D double precision arrays
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE
subroutine xpaw_mpi_bcast_dp1d(xval,master,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: spaceComm,master
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n=size(xval,dim=1)
   call MPI_BCAST(xval,n,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_bcast_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_bcast_dp2d
!! NAME
!!  xpaw_mpi_bcast_dp2d
!!
!! FUNCTION
!!  MPI_BCAST for 2D double precision arrays
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_bcast_dp2d(xval,master,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: spaceComm,master
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2)
   call MPI_BCAST(xval,n1*n2,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_bcast_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_bcast_dp3d
!! NAME
!!  xpaw_mpi_bcast_dp3d
!!
!! FUNCTION
!!  MPI_BCAST for 3D double precision arrays
!!
!! INPUTS
!!  spaceComm= MPI communicator
!!  master= master MPI node
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_bcast_dp3d(xval,master,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:,:)
 integer, intent(in) :: spaceComm,master
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2) ; n3=size(xval,dim=3)
   call MPI_BCAST(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,master,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_bcast_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_gather_int1d
!! NAME
!!  xpaw_mpi_gather_int1d
!!
!! FUNCTION
!!  MPI_GATHER for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcont= number of sent elements
!!  recvcount= number of received elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gather_int1d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: sendcount,recvcount
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: root,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_gather(xval,sendcount,MPI_INTEGER,recvbuf,recvcount,MPI_INTEGER,root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_gather_int1d
!!***

!!****f* ABINIT/xpaw_mpi_gather_dp1d
!! NAME
!!  xpaw_mpi_gather_dp1d
!!
!! FUNCTION
!!  MPI_GATHER for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcont= number of sent elements
!!  recvcount= number of received elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gather_dp1d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: sendcount,recvcount
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout)   :: recvbuf(:)
 integer, intent(in) :: root,spaceComm
 integer, intent(out) :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_gather(xval,sendcount,MPI_DOUBLE_PRECISION,recvbuf,recvcount,MPI_DOUBLE_PRECISION,&
&   root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_gather_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_gather_dp2d
!! NAME
!!  xpaw_mpi_gather_dp2d
!!
!! FUNCTION
!!  MPI_GATHER for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  sendcont= number of sent elements
!!  recvcount= number of received elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gather_dp2d(xval,sendcount,recvbuf,recvcount,root,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: sendcount,recvcount
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: root,spaceComm
 integer, intent(out)   :: ier

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   call MPI_gather(xval,sendcount,MPI_DOUBLE_PRECISION,recvbuf,recvcount,MPI_DOUBLE_PRECISION,&
&   root,spaceComm,ier)
 else if (spaceComm == xpaw_mpi_comm_self) then
   recvbuf=xval
 end if
#else
 recvbuf=xval
#endif
end subroutine xpaw_mpi_gather_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_gatherv_int1d
!! NAME
!!  xpaw_mpi_gatherv_int1d
!!
!! FUNCTION
!!  MPI_GATHERV for 1D integer arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gatherv_int1d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)

!Arguments-------------------------
 integer, intent(in) :: xval(:)
 integer, intent(inout) :: recvbuf(:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,root,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
 integer :: cc,dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /=xpaw_mpi_comm_self .and. spaceComm /=xpaw_mpi_comm_null) then
   call MPI_gatherV(xval,nelem,MPI_INTEGER,recvbuf,recvcounts,displs,&
&   MPI_INTEGER,root,spaceComm,ier)
 else if (spaceComm ==xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_gatherv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_gatherv_dp1d
!! NAME
!!  xpaw_mpi_gatherv_dp1d
!!
!! FUNCTION
!!  MPI_GATHERV for 1D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gatherv_dp1d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:)
 real(dp), intent(inout) :: recvbuf(:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,root,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
 integer :: cc,dd

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /=xpaw_mpi_comm_self .and. spaceComm /=xpaw_mpi_comm_null) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm ==xpaw_mpi_comm_self) then
#endif
   dd=0;if (size(displs)>0) dd=displs(1)
   cc=size(xval);if (size(recvcounts)>0) cc=recvcounts(1)
   recvbuf(dd+1:dd+cc)=xval(1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_gatherv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_gatherv_dp2d
!! NAME
!!  xpaw_mpi_gatherv_dp2d
!!
!! FUNCTION
!!  MPI_GATHERV for 2D double precision arrays
!!
!! INPUTS
!!  xval= buffer array
!!  recvcounts= number of received elements
!!  displs= relative offsets for incoming data
!!  nelem= number of elements
!!  root= rank of receiving process
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  recvbuf= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_gatherv_dp2d(xval,nelem,recvbuf,recvcounts,displs,root,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(in) :: xval(:,:)
 real(dp), intent(inout) :: recvbuf(:,:)
 integer, intent(in) :: recvcounts(:),displs(:)
 integer, intent(in) :: nelem,root,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
 integer :: cc,dd,sz1

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /=xpaw_mpi_comm_self .and. spaceComm /=xpaw_mpi_comm_null) then
   call MPI_gatherV(xval,nelem,MPI_DOUBLE_PRECISION,recvbuf,recvcounts,displs,&
&   MPI_DOUBLE_PRECISION,root,spaceComm,ier)
 else if (spaceComm ==xpaw_mpi_comm_self) then
#endif
   sz1=size(xval,1)
   dd=0;if (size(displs)>0) dd=displs(1)/sz1
   cc=size(xval,2);if (size(recvcounts)>0) cc=recvcounts(1)/sz1
   recvbuf(:,dd+1:dd+cc)=xval(:,1:cc)
#if defined HAVE_MPI
 end if
#endif
end subroutine xpaw_mpi_gatherv_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_recv_int1d
!! NAME
!!  xpaw_mpi_recv_int1d
!!
!! FUNCTION
!!  MPI_RECV for 1D integer arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_recv_int1d(xval,source,tag,spaceComm,ier)

!Arguments-------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: my_tag, n1
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_RECV(xval,n1,MPI_INTEGER,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif
 end subroutine xpaw_mpi_recv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_recv_dp1d
!! NAME
!!  xpaw_mpi_recv_dp1d
!!
!! FUNCTION
!!  MPI_RECV for 1D double precision arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_recv_dp1d(xval,source,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_RECV(xval,n1,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif
end subroutine xpaw_mpi_recv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_recv_dp2d
!! NAME
!!  xpaw_mpi_recv_dp2d
!!
!! FUNCTION
!!  MPI_RECV for 2D double precision arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_recv_dp2d(xval,source,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_RECV(xval,n1*n2,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif
end subroutine xpaw_mpi_recv_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_recv_dp3d
!! NAME
!!  xpaw_mpi_recv_dp3d
!!
!! FUNCTION
!!  MPI_RECV for 3D double precision arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_RECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_recv_dp3d(xval,source,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:,:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2) ; n3=size(xval,dim=3)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_RECV(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,MPI_STATUS_IGNORE,ier)
 end if
#endif
end subroutine xpaw_mpi_recv_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_irecv_int1d
!! NAME
!!  xpaw_mpi_irecv_int1d
!!
!! FUNCTION
!!  MPI_IRECV for 1D integer arrays
!!
!! INPUTS
!!  dest :: rank of destination process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFETS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_IRECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_irecv_int1d(xval,source,tag,spaceComm,request,ierr)

!Arguments-------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: request,ierr
!Local variables-------------------
#if defined HAVE_MPI
  integer :: ier,n1,my_tag
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_IRECV(xval,n1,MPI_INTEGER,source,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
 end subroutine xpaw_mpi_irecv_int1d
!!***

!!****f* ABINIT/xpaw_mpi_irecv_dp1d
!! NAME
!!  xpaw_mpi_irecv_dp1d
!!
!! FUNCTION
!!  MPI_IRECV for 1D double precision arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_IRECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_irecv_dp1d(xval,source,tag,spaceComm,request,ierr)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: request, ierr

!Local variables-------------------
#if defined HAVE_MPI
 integer :: ier,my_tag,n1
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_IRECV(xval,n1,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
end subroutine xpaw_mpi_irecv_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_irecv_dp2d
!! NAME
!!  xpaw_mpi_irecv_dp2d
!!
!! FUNCTION
!!  MPI_IRECV for 2D double precision arrays
!!
!! INPUTS
!!  source :: rank of source process
!!  tag :: integer message tag
!!  spaceComm :: MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! NOTES
!!  status of MPI_IRECV is explicitly ignored
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_irecv_dp2d(xval,source,tag,spaceComm,request,ierr)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: source,tag,spaceComm
 integer, intent(out) :: request, ierr

!Local variables-------------------
#if defined HAVE_MPI
 integer :: ier,my_tag,n1,n2
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1);n2=size(xval,dim=2)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_IRECV(xval,n1*n2,MPI_DOUBLE_PRECISION,source,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
end subroutine xpaw_mpi_irecv_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_send_int1d
!! NAME
!!  xpaw_mpi_send_int1d
!!
!! FUNCTION
!!  MPI_ISEND for 1D integer arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_send_int1d(xval,dest,tag,spaceComm,ier)

!Arguments-------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: my_tag, n1
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_SEND(xval,n1,MPI_INTEGER,dest,my_tag,spaceComm,ier)
 end if
#endif
 end subroutine xpaw_mpi_send_int1d
!!***

!!****f* ABINIT/xpaw_mpi_send_dp1d
!! NAME
!!  xpaw_mpi_send_dp1d
!!
!! FUNCTION
!!  MPI_SEND for 1D double precision arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_send_dp1d(xval,dest,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_SEND(xval,n1,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_send_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_send_dp2d
!! NAME
!!  xpaw_mpi_send_dp2d
!!
!! FUNCTION
!!  MPI_SEND for 2D double precision arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_send_dp2d(xval,dest,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_SEND(xval,n1*n2,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_send_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_send_dp3d
!! NAME
!!  xpaw_mpi_send_dp3d
!!
!! FUNCTION
!!  MPI_SEND for 3D double precision arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_send_dp3d(xval,dest,tag,spaceComm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:,:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,my_tag
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n2=size(xval,dim=2) ; n3=size(xval,dim=3)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_SEND(xval,n1*n2*n3,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_send_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_isend_int1d
!! NAME
!!  xpaw_mpi_isend_int1d
!!
!! FUNCTION
!!  MPI_ISEND for 1D integer arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_isend_int1d(xval,dest,tag,spaceComm,request,ierr)

!Arguments-------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: request,ierr

!Local variables-------------------
#if defined HAVE_MPI
 integer :: ier,my_tag,n1
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_ISEND(xval,n1,MPI_INTEGER,dest,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
 end subroutine xpaw_mpi_isend_int1d
!!***

!!****f* ABINIT/xpaw_mpi_isend_dp1d
!! NAME
!!  xpaw_mpi_isend_dp1d
!!
!! FUNCTION
!!  MPI_ISEND for 1D double precision arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_isend_dp1d(xval,dest,tag,spaceComm,request,ierr)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: request,ierr

!Local variables-------------------
#if defined HAVE_MPI
 integer :: ier,my_tag,n1
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_ISEND(xval,n1,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
end subroutine xpaw_mpi_isend_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_isend_dp2d
!! NAME
!!  xpaw_mpi_isend_dp2d
!!
!! FUNCTION
!!  MPI_ISEND for 2D double precision arrays
!!
!! INPUTS
!!  dest= rank of destination process
!!  tag= integer message tag
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ierr= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_isend_dp2d(xval,dest,tag,spaceComm,request,ierr)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: dest,tag,spaceComm
 integer, intent(out) :: request,ierr

!Local variables-------------------
#if defined HAVE_MPI
 integer :: ier,my_tag,n1,n2
#endif

! *************************************************************************
 ierr=0
#if defined HAVE_MPI
 if (spaceComm /= xpaw_mpi_comm_self .and. spaceComm /= xpaw_mpi_comm_null) then
   n1=size(xval,dim=1) ; n1=size(xval,dim=2)
   my_tag = MOD(tag,xpaw_mpi_get_tag_ub(spaceComm))
   call MPI_ISEND(xval,n1*n2,MPI_DOUBLE_PRECISION,dest,my_tag,spaceComm,request,ier)
   ierr=ier
 end if
#endif
end subroutine xpaw_mpi_isend_dp2d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_exch_int1d
!! NAME
!!  xpaw_mpi_exch_int1d
!!
!! FUNCTION
!!  MPI_SEND/MPI_RECV for 1D integer arrays
!!
!! INPUTS
!!  n1= vector length
!!  vsend= sent buffer
!!  sender= node sending the data
!!  recever= node receiving the data
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  vrecv= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_exch_int1d(vsend,n1,sender,vrecv,recever,spaceComm,ier)

!Arguments----------------
 integer, intent(in) :: n1
 integer, intent(in) :: vsend(:)
 integer, intent(inout) :: vrecv(:)
 integer, intent(in) :: sender,recever,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer :: status(MPI_STATUS_SIZE)
 integer :: tag,me
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (sender==recever.or.spaceComm==xpaw_mpi_comm_null.or.(n1==0)) return
 call MPI_COMM_RANK(spaceComm,me,ier)
 tag = MOD(n1,xpaw_mpi_get_tag_ub(spaceComm))
 if (recever==me) then
   call MPI_RECV(vrecv,n1,MPI_INTEGER,sender,tag,spaceComm,status,ier)
 end if
 if (sender==me) then
   call MPI_SEND(vsend,n1,MPI_INTEGER,recever,tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_exch_int1d
!!***

!!****f* ABINIT/xpaw_mpi_exch_dp1d
!! NAME
!!  xpaw_mpi_exch_dp1d
!!
!! FUNCTION
!!  MPI_SEND/MPI_RECV for 1D double precision arrays
!!
!! INPUTS
!!  n1= first dimension of the array
!!  vsend= send buffer
!!  sender= node sending the data
!!  recever= node receiving the data
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  vrecv= receive buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_exch_dp1d(vsend,n1,sender,vrecv,recever,spaceComm,ier)

!Arguments----------------
 integer,intent(in) :: n1
 real(dp), intent(in) :: vsend(:)
 real(dp), intent(inout) :: vrecv(:)
 integer, intent(in) :: sender,recever,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer :: status(MPI_STATUS_SIZE)
 integer :: tag,me
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (sender==recever.or.spaceComm==xpaw_mpi_comm_null.or.(n1==0)) return
 call MPI_COMM_RANK(spaceComm,me,ier)
 tag = MOD(n1,xpaw_mpi_get_tag_ub(spaceComm))
 if (recever==me) then
   call MPI_RECV(vrecv,n1,MPI_DOUBLE_PRECISION,sender,tag,spaceComm,status,ier)
 end if
 if (sender==me) then
   call MPI_SEND(vsend,n1,MPI_DOUBLE_PRECISION,recever,tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_exch_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_exch_dp2d
!! NAME
!!  xpaw_mpi_exch_dp2d
!!
!! FUNCTION
!!  MPI_SEND/MPI_RECV for 2D double precision arrays
!!
!! INPUTS
!!  nt= vector length
!!  vsend= sent buffer
!!  sender= node sending the data
!!  recever= node receiving the data
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  vrecv= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_exch_dp2d(vsend,nt,sender,vrecv,recever,spaceComm,ier)

!Arguments----------------
 integer,intent(in) :: nt
 real(dp), intent(in) :: vsend(:,:)
 real(dp), intent(inout) :: vrecv(:,:)
 integer, intent(in) :: sender,recever,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer :: status(MPI_STATUS_SIZE)
 integer :: tag,me
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (sender==recever.or.spaceComm==xpaw_mpi_comm_null.or.(nt==0)) return
 call MPI_COMM_RANK(spaceComm,me,ier)
 tag = MOD(nt,xpaw_mpi_get_tag_ub(spaceComm))
 if (recever==me) then
   call MPI_RECV(vrecv,nt,MPI_DOUBLE_PRECISION,sender,tag,spaceComm,status,ier)
 end if
 if (sender==me) then
   call MPI_SEND(vsend,nt,MPI_DOUBLE_PRECISION,recever,tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_exch_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_exch_dp3d
!! NAME
!!  xpaw_mpi_exch_dp3d
!!
!! FUNCTION
!!  MPI_SEND/MPI_RECV for 3D double precision arrays
!!
!! INPUTS
!!  nt= vector length
!!  vsend= sent buffer
!!  sender= node sending the data
!!  recever= node receiving the data
!!  spaceComm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  vrecv= received buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_exch_dp3d(vsend,nt,sender,vrecv,recever,spaceComm,ier)

!Arguments----------------
 integer,intent(in) :: nt
 real(dp), intent(in) :: vsend(:,:,:)
 real(dp), intent(inout) :: vrecv(:,:,:)
 integer, intent(in) :: sender,recever,spaceComm
 integer, intent(out) :: ier

!Local variables--------------
#if defined HAVE_MPI
 integer :: status(MPI_STATUS_SIZE)
 integer :: tag,me
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (sender==recever.or.spaceComm==xpaw_mpi_comm_null.or.(nt==0)) return
 call MPI_COMM_RANK(spaceComm,me,ier)
 tag = MOD(nt,xpaw_mpi_get_tag_ub(spaceComm))
 if (recever==me) then
   call MPI_RECV(vrecv,nt,MPI_DOUBLE_PRECISION,sender,tag,spaceComm,status,ier)
 end if
 if (sender==me) then
   call MPI_SEND(vsend,nt,MPI_DOUBLE_PRECISION,recever,tag,spaceComm,ier)
 end if
#endif
end subroutine xpaw_mpi_exch_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_sum_int
!! NAME
!!  xpaw_mpi_sum_int
!!
!! FUNCTION
!!  MPI_ALLREDUCE(SUM) for integers
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_sum_int(xval,comm,ier)

!Arguments-------------------------
 integer, intent(inout) :: xval
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: nproc
#if !defined HAVE_MPI2 || !defined HAVE_MPI2_INPLACE
 integer :: xsum
#endif
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nproc,ier)
   if (nproc /= 1) then
#if defined HAVE_MPI2 && defined HAVE_MPI2_INPLACE
     call MPI_ALLREDUCE(MPI_IN_PLACE,xval,1,MPI_INTEGER,MPI_SUM,comm,ier)
#else
     call MPI_ALLREDUCE(xval,xsum,1,MPI_INTEGER,MPI_SUM,comm,ier)
     xval=xsum
#endif
   end if
 end if
#endif
end subroutine xpaw_mpi_sum_int
!!***

!!****f* ABINIT/xpaw_mpi_sum_int1d
!! NAME
!!  xpaw_mpi_sum_int1d
!!
!! FUNCTION
!!  MPI_ALLREDUCE(SUM) for 1D integer arrays
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_sum_int1d(xval,comm,ier)

!Arguments-------------------------
 integer, intent(inout) :: xval(:)
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,nproc
#if !defined HAVE_MPI2 || !defined HAVE_MPI2_INPLACE
 integer :: xsum(size(xval,dim=1))
#endif
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nproc,ier)
   if (nproc /= 1) then
     n1=size(xval,dim=1)
#if defined HAVE_MPI2 && defined HAVE_MPI2_INPLACE
     call MPI_ALLREDUCE(MPI_IN_PLACE,xval,n1,MPI_INTEGER,MPI_SUM,comm,ier)
#else
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_INTEGER,MPI_SUM,comm,ier)
     xval(:)=xsum(:)
#endif
   end if
 end if
#endif
end subroutine xpaw_mpi_sum_int1d
!!***

!!****f* ABINIT/xpaw_mpi_sum_dp1d
!! NAME
!!  xpaw_mpi_sum_dp1d
!!
!! FUNCTION
!!  MPI_ALLREDUCE(SUM) for 1D double precision arrays
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_sum_dp1d(xval,comm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:)
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,nproc
#if !defined HAVE_MPI2 || !defined HAVE_MPI2_INPLACE
 real(dp) :: xsum(size(xval,dim=1))
#endif
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nproc,ier)
   if (nproc /= 1) then
     n1=size(xval,dim=1)
#if defined HAVE_MPI2 && defined HAVE_MPI2_INPLACE
     call MPI_ALLREDUCE(MPI_IN_PLACE,xval,n1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
#else
     call MPI_ALLREDUCE(xval,xsum,n1,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
     xval(:)=xsum(:)
#endif
   end if
 end if
#endif
end subroutine xpaw_mpi_sum_dp1d
!!***

!!****f* ABINIT/xpaw_mpi_sum_dp2d
!! NAME
!!  xpaw_mpi_sum_dp2d
!!
!! FUNCTION
!!  MPI_ALLREDUCE(SUM) for 2D double precision arrays
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_sum_dp2d(xval,comm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:)
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,nproc
#if !defined HAVE_MPI2 || !defined HAVE_MPI2_INPLACE
 real(dp) :: xsum(size(xval,dim=1),size(xval,dim=2))
#endif
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nproc,ier)
   if (nproc /= 1) then
     n1=size(xval,dim=1) ; n2=size(xval,dim=2)
#if defined HAVE_MPI2 && defined HAVE_MPI2_INPLACE
     call MPI_ALLREDUCE(MPI_IN_PLACE,xval,n1*n2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
#else
     call MPI_ALLREDUCE(xval,xsum,n1*n2,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
     xval(:,:)=xsum(:,:)
#endif
   end if
 end if
#endif
end subroutine xpaw_mpi_sum_dp2d
!!***

!!****f* ABINIT/xpaw_mpi_sum_dp3d
!! NAME
!!  xpaw_mpi_sum_dp3d
!!
!! FUNCTION
!!  MPI_ALLREDUCE(SUM) for 3D double precision arrays
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  ier= exit status, a non-zero value meaning there is an error
!!
!! SIDE EFFECTS
!!  xval= buffer array
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_allreduce,mpi_attr_get,mpi_comm_size
!!
!! SOURCE

subroutine xpaw_mpi_sum_dp3d(xval,comm,ier)

!Arguments-------------------------
 real(dp), intent(inout) :: xval(:,:,:)
 integer, intent(in) :: comm
 integer, intent(out) :: ier

!Local variables-------------------
#if defined HAVE_MPI
 integer :: n1,n2,n3,nproc
#if !defined HAVE_MPI2 || !defined HAVE_MPI2_INPLACE
 real(dp) :: xsum(size(xval,dim=1),size(xval,dim=2),size(xval,dim=3))
#endif
#endif

! *************************************************************************
 ier=0
#if defined HAVE_MPI
 if (comm /= xpaw_mpi_comm_self .and. comm /= xpaw_mpi_comm_null) then
   call MPI_COMM_SIZE(comm,nproc,ier)
   if (nproc /= 1) then
     n1=size(xval,dim=1) ; n2=size(xval,dim=2) ; n3=size(xval,dim=3)
#if defined HAVE_MPI2 && defined HAVE_MPI2_INPLACE
     call MPI_ALLREDUCE(MPI_IN_PLACE,xval,n1*n2*n3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
#else
     call MPI_ALLREDUCE(xval,xsum,n1*n2*n3,MPI_DOUBLE_PRECISION,MPI_SUM,comm,ier)
     xval(:,:,:)=xsum(:,:,:)
#endif
   end if
 end if
#endif
end subroutine xpaw_mpi_sum_dp3d
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/xpaw_mpi_get_tag_ub
!! NAME
!!  xpaw_mpi_get_tag_ub
!!
!! FUNCTION
!!  Get MPI_TAG_UB attribute
!!
!! INPUTS
!!  comm= MPI communicator
!!
!! OUTPUT
!!  xpaw_mpi_get_tag_ub=value for the MPI_TAG_UB attribute attached to comm
!!
!! PARENTS
!!
!! CHILDREN
!!      mpi_attr_get
!!
!! SOURCE

function xpaw_mpi_get_tag_ub(comm)

!Arguments-------------------------
 integer, intent(in) :: comm
 integer :: xpaw_mpi_get_tag_ub

!Local variables-------------------
#if defined HAVE_MPI
 integer :: attribute_val,ier
 logical :: lflag
#endif

! *************************************************************************

#if defined HAVE_MPI
 !Deprecated in MPI2 but not all MPI2 implementations provide MPI_Comm_get_attr !
 call MPI_ATTR_GET(comm,MPI_TAG_UB,attribute_val,lflag,ier)
!call MPI_Comm_get_attr(comm MPI_TAG_UB,attribute_val,lflag,ier)

 if (lflag) xpaw_mpi_get_tag_ub = attribute_val

#else
 xpaw_mpi_get_tag_ub=32767
#endif

end function xpaw_mpi_get_tag_ub
!!***

!----------------------------------------------------------------------

end module m_libpaw_mpi
!!***

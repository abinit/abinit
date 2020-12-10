!!****m* ABINIT/m_xredistribute
!! NAME
!!  m_xredistribute
!!
!! FUNCTION
!!  This module contains a function to re-dristribute
!!  results on different procs.
!!
!! COPYRIGHT
!! Copyright (C) 2008-2020 ABINIT group (MMANCINI)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!
!! OUTPUT
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_xredistribute

 use defs_basis
 use m_errors
 use m_abicore
#if defined HAVE_MPI2
 use mpi
#endif

 implicit none

#if defined HAVE_MPI1
 include 'mpif.h'
#endif

 private

 public :: xredistribute         ! Redistribute the work


 interface xredistribute
  module procedure xredistribute_mpi_dp
  module procedure xredistribute_mpi_2d_dp
 end interface xredistribute


CONTAINS  !===========================================================
!!***

!!****f* ABINIT/xredistribute_mpi_dp
!! NAME
!!  xredistribute_mpi_dp
!!
!! FUNCTION
!!
!! INPUTS
!!  xval= buffer array
!!  send_counts= number of sent elements (initial distribution)
!!  send_displs= postions of values sent by the processor (initial positions)
!!  rec_counts= number of received elements (final distribution)
!!  rec_displs= positions of values received by the processors (final position)
!!  nproc=number of processor
!!  me=proc me
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
!!      mpi_allgatherv,mpi_scatterv
!!
!! SOURCE
subroutine xredistribute_mpi_dp(xval,send_counts,send_displs,recvbuf,&
                                rec_counts,rec_displs,me,nproc,spaceComm,ier)

!Arguments-------------------------
 integer ,intent(in) :: me,nproc
 real(dp),intent(in) :: xval(:)
 real(dp),intent(inout) :: recvbuf(:)
 integer ,intent(in) :: send_displs(0:nproc-1),send_counts(0:nproc-1)
 integer ,intent(in) :: rec_displs(0:nproc-1),rec_counts(0:nproc-1)
 integer ,intent(in) :: spaceComm
 integer ,intent(out):: ier

!Local variables-------------------
 integer :: size
 real(dp),allocatable :: totbuff(:)
 character(500) :: msg
! *********************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
 size = sum(send_counts)
 if(size /=sum(rec_counts))then
   msg = 'the total sizes of sent and receved msg are not equal'
   ABI_ERROR(msg)
 endif

 ABI_ALLOCATE(totbuff,(size))
   !--join all the vector in to a single one
   call MPI_ALLGATHERV(xval,send_counts(me),MPI_DOUBLE_PRECISION,totbuff,&
&                      send_counts,send_displs,MPI_DOUBLE_PRECISION,spaceComm,ier)


   !--now distribute the total vector on the procs
   call MPI_SCATTERV(totbuff,rec_counts,rec_displs,MPI_DOUBLE_PRECISION,&
&                    recvbuf,rec_counts(me),MPI_DOUBLE_PRECISION,&
&                    0,spaceComm,ier)

   ABI_DEALLOCATE(totbuff)
 end if
#endif
end subroutine xredistribute_mpi_dp
!!***


!!****f* ABINIT/xredistribute_mpi_2d_dp
!! NAME
!!  xredistribute_mpi_2d_dp
!!
!! FUNCTION
!!
!! INPUTS
!!  xval= buffer array
!!  send_counts= number of sent elements (initial distribution)
!!  send_displs= postions of values sent by the processor (initial positions)
!!  rec_counts= number of received elements (final distribution)
!!  rec_displs= positions of values received by the processors (final position)
!!  nproc=number of processor
!!  me=proc me
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
!!      mpi_allgatherv,mpi_scatterv
!!
!! SOURCE
subroutine xredistribute_mpi_2d_dp(xval,send_counts,send_displs,recvbuf,&
                                   rec_counts,rec_displs,me,nproc,spaceComm,ier)

!Arguments-------------------------
 integer ,intent(in) :: me,nproc
 real(dp),intent(in) :: xval(:,:)
 real(dp),intent(inout) :: recvbuf(:,:)
 integer ,intent(in) :: send_displs(0:nproc-1),send_counts(0:nproc-1)
 integer ,intent(in) :: rec_displs(0:nproc-1),rec_counts(0:nproc-1)
 integer ,intent(in) :: spaceComm
 integer ,intent(out):: ier

!Local variables-------------------
 integer :: size
 real(dp),allocatable :: totbuff(:)
 character(500) :: msg
! *********************************************************************

 ier=0
#if defined HAVE_MPI
 if (spaceComm /= MPI_COMM_SELF .and. spaceComm /= MPI_COMM_NULL) then
 size = sum(send_counts)
 if(size /=sum(rec_counts))then
   msg = 'the total sizes of sent and receved msg are not equal'
   ABI_ERROR(msg)
 endif

 ABI_ALLOCATE(totbuff,(size))
   !--join all the vector in to a single one
   call MPI_ALLGATHERV(xval,send_counts(me),MPI_DOUBLE_PRECISION,totbuff,&
&                      send_counts,send_displs,MPI_DOUBLE_PRECISION,spaceComm,ier)


   !--now distribute the total vector on the procs
   call MPI_SCATTERV(totbuff,rec_counts,rec_displs,MPI_DOUBLE_PRECISION,&
&                    recvbuf,rec_counts(me),MPI_DOUBLE_PRECISION,&
&                    0,spaceComm,ier)
   ABI_DEALLOCATE(totbuff)
 end if
#endif
end subroutine xredistribute_mpi_2d_dp
!!***

end module m_xredistribute
!!***

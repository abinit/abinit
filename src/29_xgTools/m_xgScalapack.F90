!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xgScalapack
!! NAME
!!  m_xgScalapack
!! 
!! FUNCTION 
!! 
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xgScalapack
  use defs_basis, only : std_err, std_out
  use m_profiling_abi
  use m_xmpi 
  use m_errors
  use m_slk
  use m_xg

  implicit none

  private

  integer, parameter :: _SLK = 1
  integer, parameter :: _ROW = 2
  integer, parameter :: _COL = 3
  integer, parameter :: _UNUSED = 4
  integer, parameter :: _WORLD = 5
  integer, parameter :: _NDATA = 5

  type, public :: xgScalapack_t
    integer :: comms(_NDATA)
    integer :: rank(_NDATA)
    integer :: size(_NDATA)


  end type xgScalapack_t

  contains 
!!***

!!****f* m_xgScalapack/xgScalapack_init
!! NAME
!!  xgScalapack_init
!!
!! FUNCTION
!!  Init the scalapack communicator for next operations.
!!  If the comm has to many cpus, then take only a subgroup of this comm
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
  subroutine  xgScalapack_init(xgScalapack,comm,maxData)
    type(xgScalapack_t), intent(inout) :: xgScalapack
    integer            , intent(in   ) :: comm
    integer            , intent(in   ) :: maxData
    integer :: maxProc
    integer :: nproc
    integer :: subgroup
    integer :: mycomm(2)
    integer :: ierr

    xgScalapack%comms = xmpi_comm_null
    xgScalapack%rank = MPI_UNDEFINED_RANK

    nproc = xmpi_comm_size(comm)
    xgScalapack%comms(_WORLD) = comm
    xgScalapack%rank(_WORLD) = xmpi_comm_rank(comm)
    xgScalapack%size(_WORLD) = nproc

    maxProc = MAX(1,maxData / 1000000) ! ( 1000 x 1000 matrice per MPI )
    if ( maxProx < nproc ) then
      if ( xgScalapack%rank(_WORLD) < maxProc ) then
        subgroup = 0
        mycomm(1) = _SLK
        mycomm(2) = _UNUSED
      else
        subgroup = 1
        mycomm(1) = _UNUSED
        mycomm(2) = _SLK
      end if
       call MPI_Comm_split(comm, subgroup, xgScalapack%rank(_WORLD), xgScalapack%comms(mycomm(1)),ierr);
       xgScalapack%comms(mycomm(2)) = xmpi_comm_null
       xgScalapack%rank(mycomm(1)) = xmpi_comm_rank(xgScalapack%comms(mycomm(1))
       xgScalapack%rank(mycomm(2)) = MPI_UNDEFINED_RANK
       xgScalapack%size(mycomm(1)) = xmpi_comm_size(xgScalapack%comms(mycomm(1))
       xgScalapack%size(mycomm(2)) = MPI_UNDEFINED
    else 
       call MPI_Comm_dup(comm,xgScalapack%comms(_SLK))
       xgScalapack%rank(_SLK)) = xmpi_comm_rank(xgScalapack%comms(_SLK))
       xgScalapack%size(_SLK) = nproc
    end if

    call init_scalapack(xgScalapack%processor,xgScalapack%comms(_SLK))

  end subroutine xgScalapack_init

  subroutine  xgScalapack_free(xgScalapack)
    if ( xgScalapack%comms(_SLK) /= xmpi_comm_null ) then
      call MPI_Comm_free(xgScalapack%comms(_SLK))
    end if
    if ( xgScalapack%comms(_UNUSED) /= xmpi_comm_null ) then
      call MPI_Comm_free(xgScalapack%comms(_UNUSED))
    end if 
  end subroutine xgScalapack_free

end module m_xgScalapack

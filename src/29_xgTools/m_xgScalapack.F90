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
  use defs_basis, only : std_err, std_out, dp
  use m_profiling_abi
  use m_xmpi 
  use m_errors
  use m_slk
  use m_xg

  implicit none

  private

  integer, parameter :: M_SLK = 1
  integer, parameter :: M_ROW = 1
  integer, parameter :: M_COL = 2
  integer, parameter :: M_UNUSED = 4
  integer, parameter :: M_WORLD = 5
  integer, parameter :: M_NDATA = 5

  type, public :: xgScalapack_t
    integer :: comms(M_NDATA)
    integer :: rank(M_NDATA)
    integer :: size(M_NDATA)
    integer :: coords(2)
    type(grid_scalapack) :: grid
  end type xgScalapack_t

  public :: xgScalapack_init
  public :: xgScalapack_free
  public :: xgScalapack_heev
  public :: xgScalapack_hegv
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
  subroutine  xgScalapack_init(xgScalapack,comm,maxData,usable)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgScalapack_init'
!End of the abilint section

    type(xgScalapack_t), intent(inout) :: xgScalapack
    integer            , intent(in   ) :: comm
    integer            , intent(in   ) :: maxData
    logical            , intent(  out) :: usable
    integer :: maxProc
    integer :: nproc
    integer :: subgroup
    integer :: mycomm(2)
    integer :: ierr
    integer :: test_row
    integer :: test_col

#ifdef HAVE_SCALAPACK
    xgScalapack%comms = xmpi_comm_null
    xgScalapack%rank = xmpi_undefined_rank

    nproc = xmpi_comm_size(comm)
    xgScalapack%comms(M_WORLD) = comm
    xgScalapack%rank(M_WORLD) = xmpi_comm_rank(comm)
    xgScalapack%size(M_WORLD) = nproc

    maxProc = MAX(1,maxData / 1000000) ! ( 1000 x 1000 matrice per MPI )
    if ( maxProc == 1 ) then
      usable = .false.
      return
    else
      usable = .true.
    end if

    if ( maxProc < nproc ) then
      if ( xgScalapack%rank(M_WORLD) < maxProc ) then
        subgroup = 0
        mycomm(1) = M_SLK
        mycomm(2) = M_UNUSED
      else
        subgroup = 1
        mycomm(1) = M_UNUSED
        mycomm(2) = M_SLK
      end if
       call MPI_Comm_split(comm, subgroup, xgScalapack%rank(M_WORLD), xgScalapack%comms(mycomm(1)),ierr)
       xgScalapack%comms(mycomm(2)) = xmpi_comm_null
       xgScalapack%rank(mycomm(1)) = xmpi_comm_rank(xgScalapack%comms(mycomm(1)))
       xgScalapack%rank(mycomm(2)) = xmpi_undefined_rank
       xgScalapack%size(mycomm(1)) = xmpi_comm_size(xgScalapack%comms(mycomm(1)))
       xgScalapack%size(mycomm(2)) = xmpi_undefined
    else 
       call MPI_Comm_dup(comm,xgScalapack%comms(M_SLK))
       xgScalapack%rank(M_SLK) = xmpi_comm_rank(xgScalapack%comms(M_SLK))
       xgScalapack%size(M_SLK) = nproc
    end if

    if ( xgScalapack%comms(M_SLK) /= xmpi_comm_null ) then
      call build_grid_scalapack(xgScalapack%grid, xgScalapack%size(M_SLK), xgScalapack%comms(M_SLK))
      call BLACS_GridInfo(xgScalapack%grid%ictxt, &
        xgScalapack%grid%dims(M_ROW), xgScalapack%grid%dims(M_COL),&
        xgScalapack%coords(M_ROW), xgScalapack%coords(M_COL))
     
     !These values are the same as those computed by BLACS_GRIDINFO
     !except in the case where the myproc argument is not the local proc
      test_row = INT((xgScalapack%rank(M_SLK)) / xgScalapack%grid%dims(2))
      test_col = MOD((xgScalapack%rank(M_SLK)), xgScalapack%grid%dims(2))
      if ( test_row /= xgScalapack%coords(M_ROW) ) then
        MSG_WARNING("Row id mismatc")
      end if
      if ( test_col /= xgScalapack%coords(m_COL) ) then
        MSG_WARNING("Col id mismatc")
      end if
    end if
#else
    usable = .false.
#endif
  end subroutine xgScalapack_init

  function toProcessorScalapack(xgScalapack) result(processor)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'toProcessorScalapack'
!End of the abilint section

    type(xgScalapack_t), intent(in) :: xgScalapack
    type(processor_scalapack) :: processor

    processor%myproc = xgScalapack%rank(M_SLK)
    processor%comm = xgScalapack%comms(M_SLK)
    processor%coords = xgScalapack%coords
  end function toProcessorScalapack

  !This is for testing purpose.
  !May not be optimal since I do not control old implementation but at least gives a reference.
  subroutine xgScalapack_heev(xgScalapack,matrixA,eigenvalues)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgScalapack_heev'
!End of the abilint section

    type(xgScalapack_t), intent(inout) :: xgScalapack
    type(xgBlock_t)    , intent(inout) :: matrixA
    type(xgBlock_t)    , intent(inout) :: eigenvalues
    real(dp), pointer :: matrix(:,:) !(cplex*nbli_global,nbco_global)
    real(dp), pointer :: eigenvalues_tmp(:,:)
    real(dp), pointer :: vector(:)
    integer :: cplex
    integer :: istwf_k
    integer :: nbli_global, nbco_global
    type(c_ptr) :: cptr

    ! Keep only working processors 
    if ( xgScalapack%comms(M_SLK) == xmpi_comm_null ) return

    call xgBlock_getSize(eigenvalues,nbli_global,nbco_global)
    if ( cols(matrixA) /= nbli_global ) then
      MSG_ERROR("Number of eigen values differ from number of vectors")
    end if

    if ( space(matrixA) == SPACE_C ) then
      cplex = 2
      istwf_k = 1
    else
      cplex = 1
      istwf_k = 2
    endif

    call xgBlock_getSize(matrixA,nbli_global,nbco_global)

    call xgBlock_reverseMap(matrixA,matrix,nbli_global,nbco_global)
    call xgBlock_reverseMap(eigenvalues,eigenvalues_tmp,nbco_global,1)
    cptr = c_loc(eigenvalues_tmp)
    call c_f_pointer(cptr,vector,(/ nbco_global /))

    call compute_eigen1(xgScalapack%comms(M_SLK), &
      toProcessorScalapack(xgScalapack), &
      cplex,nbli_global,nbco_global,matrix,vector,istwf_k)

  end subroutine xgScalapack_heev

  !This is for testing purpose.
  !May not be optimal since I do not control old implementation but at least gives a reference.
  subroutine xgScalapack_hegv(xgScalapack,matrixA,matrixB,eigenvalues)
    use iso_c_binding

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgScalapack_hegv'
!End of the abilint section

    type(xgScalapack_t), intent(inout) :: xgScalapack
    type(xgBlock_t)    , intent(inout) :: matrixA
    type(xgBlock_t)    , intent(inout) :: matrixB
    type(xgBlock_t)    , intent(inout) :: eigenvalues
    real(dp), pointer :: matrix1(:,:) !(cplex*nbli_global,nbco_global)
    real(dp), pointer :: matrix2(:,:) !(cplex*nbli_global,nbco_global)
    real(dp), pointer :: eigenvalues_tmp(:,:)
    real(dp), pointer :: vector(:)
    integer :: cplex
    integer :: istwf_k
    integer :: nbli_global, nbco_global
    type(c_ptr) :: cptr

    ! Keep only working processors 
    if ( xgScalapack%comms(M_SLK) == xmpi_comm_null ) return

    call xgBlock_getSize(eigenvalues,nbli_global,nbco_global)
    if ( cols(matrixA) /= cols(matrixB) ) then
      MSG_ERROR("Matrix A and B don't have the same number of vectors")
    end if
    if ( cols(matrixA) /= nbli_global ) then
      MSG_ERROR("Number of eigen values differ from number of vectors")
    end if

    if ( space(matrixA) == SPACE_C ) then
      cplex = 2
      istwf_k = 1
    else
      cplex = 1
      istwf_k = 2
    endif

    call xgBlock_getSize(matrixA,nbli_global,nbco_global)

    call xgBlock_reverseMap(matrixA,matrix1,nbli_global,nbco_global)
    call xgBlock_reverseMap(matrixB,matrix2,nbli_global,nbco_global)
    call xgBlock_reverseMap(eigenvalues,eigenvalues_tmp,nbco_global,1)
    cptr = c_loc(eigenvalues_tmp)
    call c_f_pointer(cptr,vector,(/ nbco_global /))

    call compute_eigen2(xgScalapack%comms(M_SLK), &
      toProcessorScalapack(xgScalapack), &
      cplex,nbli_global,nbco_global,matrix1,matrix2,vector,istwf_k)

  end subroutine xgScalapack_hegv

  subroutine  xgScalapack_free(xgScalapack)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgScalapack_free'
!End of the abilint section

    type(xgScalapack_t), intent(inout) :: xgScalapack

    call BLACS_GridExit(xgScalapack%grid%ictxt)
    if ( xgScalapack%comms(M_SLK) /= xmpi_comm_null ) then
      call MPI_Comm_free(xgScalapack%comms(M_SLK))
    end if
    if ( xgScalapack%comms(M_UNUSED) /= xmpi_comm_null ) then
      call MPI_Comm_free(xgScalapack%comms(M_UNUSED))
    end if 
  end subroutine xgScalapack_free

end module m_xgScalapack

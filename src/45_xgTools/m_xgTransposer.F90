!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_xgTransposer
!! NAME
!!  m_xgTransposer
!! 
!! FUNCTION 
!! This module is to be user to go to "KGB" representation and to "linear
!! algebra representation" It will replace most of prep_* subroutine
!! This should really help to do the transposition operataion
!! 
!! COPYRIGHT
!!  Copyright (C) 2017 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
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

module m_xgTransposer
  use defs_basis, only : std_err, std_out
  use m_profiling_abi
  use m_xmpi
  use m_errors
  use m_xg

#ifdef HAVE_MPI2
 use mpi
#endif

  implicit none

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

  private

  integer, parameter, public :: STATE_LINALG = 1
  integer, parameter, public :: STATE_COLSROWS  = 2
  integer, parameter, public :: STATE_UNKNOW    = 3
  integer, parameter :: MPI_LINALG = 1
  integer, parameter :: MPI_ROWS   = 2
  integer, parameter :: MPI_COLS   = 3
  integer, parameter :: MPI_2DCART = 4
  integer, parameter, public :: TRANS_ALL2ALL = 1
  integer, parameter, public :: TRANS_GATHER = 2

  integer, parameter :: tim_toColsRows = 1670

  type, private :: mpiData_t
    integer :: comm
    integer :: rank
    integer :: size
  end type mpiData_t

  type, private :: ptr_t
    double precision, pointer :: ptr(:,:) => null()
  end type ptr_t

  type, public :: xgTransposer_t
    type(xgBlock_t), pointer, private :: xgBlock_linalg => null()
    type(xgBlock_t), pointer, private :: xgBlock_colsrows => null()
    integer :: state
    type(mpiData_t), private :: mpiData(4)
    integer, allocatable, private :: lookup(:)
    integer, allocatable, private :: nrows_ind(:)
    integer :: mpiAlgo
  end type xgTransposer_t

  public :: xgTransposer_init
  public :: xgTransposer_transpose
  public :: xgTransposer_free


  contains
!!***

!!****f* m_xgTransposer/xgTransposer_init
!!
!! NAME
!! xgTransposer_init

  subroutine xgTransposer_init(xgTransposer,xgBlock_linalg,xgBlock_colsrows,ncpuRows,ncpuCols,state,algo)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_init'
!End of the abilint section

    type(xgTransposer_t)   , intent(inout) :: xgTransposer
    type(xgBlock_t), target, intent(in   ) :: xgBlock_linalg
    type(xgBlock_t), target, intent(in   ) :: xgBlock_colsrows
    integer                , intent(inout) :: ncpuRows
    integer                , intent(inout) :: ncpuCols
    integer                , intent(in   ) :: state
    integer                , intent(in   ) :: algo
    integer :: commLinalg
    integer :: ncols
    integer :: nrows
    integer :: ierr
    integer :: icol
    integer :: icplx
    character(len=500) :: message

    xgTransposer%xgBlock_linalg => xgBlock_linalg
    xgTransposer%xgBlock_colsrows => xgBlock_colsrows
    xgTransposer%state = state
    commLinalg = comm(xgBlock_linalg)
    xgTransposer%mpiData(MPI_LINALG)%comm = commLinalg
    xgTransposer%mpiData(MPI_LINALG)%rank = xmpi_comm_rank(commLinalg)
    xgTransposer%mpiData(MPI_LINALG)%size = xmpi_comm_size(commLinalg)

!    if ( space(xgBlock_linalg) /= space(xgBlock_colsrows) ) then
!      MSG_ERROR("Linalg xgBlock and ColsRows xgBlocks are not in the same space")
!    end if

    if ( xgTransposer%mpiData(MPI_LINALG)%size < ncpuCols*ncpuRows ) then
      write(message,'(a,i6,a,i6,a)') "There is not enough MPI processes in the communcation (", &
        xgTransposer%mpiData(MPI_LINALG)%size, "). Need at least ", ncpuCols*ncpuRows, " processes"
      MSG_ERROR(message)
    end if

    if ( ( algo == TRANS_ALL2ALL .or. algo == TRANS_GATHER ) ) then
      xgTransposer%mpiAlgo = algo
    else
      xgTransposer%mpiAlgo = TRANS_ALL2ALL
      MSG_COMMENT("Bad value for transposition MPI_algo")
    end if

    if ( xgTransposer%mpiAlgo == TRANS_ALL2ALL ) then
      MSG_COMMENT("Using mpi_ialltoall for transposition")
    else
      MSG_COMMENT("Using mpi_igatherv for transposition")
    end if

    select case (state) 
    case (STATE_LINALG)
      ! We are in the linalg representation. 
      ! We need to construct the colsrows parallelization
      call xgBlock_getSize(xgBlock_linalg,nrows,ncols)
      call xmpi_sum(nrows,commLinalg,ierr)

      if ( MOD(ncols,ncpuCols) /=0 ) then
        if ( ncols > ncpuCols ) then
          write(message,'(a,i6,a,i6,a)') "Unbalanced parallelization : ", ncols, " columns for ", ncpuCols, " MPI"
          MSG_ERROR(message)
        else
          write(message,'(i6,a)') (ncpuCols-ncols)*ncpuRows, " MPI will not be used  because of the number of columns!!"
          MSG_ERROR(message)
          ncpuCols = ncols
        end if
      end if

      select case ( space(xgBlock_linalg) )
        case (SPACE_CR, SPACE_R)
          write(*,*) "SPACE R"
          icplx = 1
        case (SPACE_C)
          write(*,*) "SPACE C"
          icplx = 2
        case default
          MSG_ERROR("Space value unknown !")
      end select

      write(*,*) "There is ", nrows, "rows of real or complex for ", ncpuRows, "fft cpu"
      if ( MOD(nrows,ncpuRows) /=0 ) then
        if ( nrows > ncpuRows ) then
          write(message,'(a,i6,a,i6,a)') "Unbalanced parallelization : ", nrows, " rows for ", ncpuRows, " MPI"
          MSG_ERROR(message)
        else
          write(message,'(i6,a)') (ncpuRows-nrows)*ncpuCols, " MPI will not be used because of the number of rows!!"
          MSG_ERROR(message)
          ncpuRows = nrows
        end if
      end if

      ! Build the lookup table
      ABI_MALLOC(xgTransposer%lookup,(1:ncols))
      do icol = 0, ncols-1
        xgTransposer%lookup(icol) = MOD(icol,ncpuCols)
      end do

      call xgTransposer_makeComm(xgTransposer,ncpuRows,ncpuCols)

    case (STATE_COLSROWS)
      MSG_BUG("Not yet implemented")
    case default
      MSG_ERROR("State is undefined")
    end select

  end subroutine xgTransposer_init
!!***

!!****f* m_xgTransposer/xgTransposer_makeComm
!!
!! NAME
!! xgTransposer_makeComm

  subroutine xgTransposer_makeComm(xgTransposer,ncpuRows,ncpuCols)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_makeComm'
!End of the abilint section

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer             , intent(in   ) :: ncpuRows
    integer             , intent(in   ) :: ncpuCols
#if defined HAVE_MPI
    integer :: coordInGrid(2),sizeGrid(2)
    logical :: periodic(2), selectDim(2), reorder
    integer :: ierr
    integer :: commColsRows

    sizeGrid(1) = ncpuRows
    sizeGrid(2) = ncpuCols
    periodic = (/ .false., .false. /)
    reorder  = .false.
    call mpi_cart_create(xgTransposer%mpiData(MPI_LINALG)%comm,2,sizeGrid,periodic,reorder,commColsRows,ierr)
    if ( ierr /= MPI_SUCCESS ) then
      MSG_ERROR("xgTransposer failed to creat cartesian grid")
    end if

    selectDim = (/ .true., .false. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_ROWS)%comm,ierr)
    if ( ierr /= MPI_SUCCESS ) then
      MSG_ERROR("xgTransposer failed to creat rows communicator")
    end if
    selectDim = (/ .false., .true. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_COLS)%comm,ierr)
    if ( ierr /= MPI_SUCCESS ) then
      MSG_ERROR("xgTransposer failed to creat columns communicator")
    end if
#else 
    commColsRows = xmpi_comm_null
    xgTransposer%mpiData(MPI_ROWS)%comm = xmpi_comm_null
    xgTransposer%mpiData(MPI_COLS)%comm = xmpi_comm_null
#endif

    xgTransposer%mpiData(MPI_2DCART)%comm = commColsRows

    xgTransposer%mpiData(MPI_ROWS)%rank = xmpi_comm_rank(xgTransposer%mpiData(MPI_ROWS)%comm)
    xgTransposer%mpiData(MPI_ROWS)%size = xmpi_comm_size(xgTransposer%mpiData(MPI_ROWS)%comm)

    xgTransposer%mpiData(MPI_COLS)%rank = xmpi_comm_rank(xgTransposer%mpiData(MPI_COLS)%comm)
    xgTransposer%mpiData(MPI_COLS)%size = xmpi_comm_size(xgTransposer%mpiData(MPI_COLS)%comm)

    call xgBlock_setComm(xgTransposer%xgBlock_colsrows,xgTransposer%mpiData(MPI_ROWS)%comm)

  end subroutine xgTransposer_makeComm
!!***

!!****f* m_xgTransposer/xgTransposer_transpose
!!
!! NAME
!! xgTransposer_transpose

  subroutine xgTransposer_transpose(xgTransposer,toState)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_transpose'
!End of the abilint section

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer             , intent(in   ) :: toState

    if ( toState /= STATE_LINALG .and. toState /= STATE_COLSROWS ) then
      MSG_ERROR("Bad value for toState")
    end if

    select case (toState)
    case (STATE_LINALG)
      if ( xgTransposer%state == STATE_LINALG ) then
        MSG_WARNING("Array linalg has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_LINALG)%size == 1 ) then
        call xgBlock_copy(xgTransposer%xgBlock_colsrows,xgTransposer%xgBlock_linalg)
      else
        call xgTransposer_toLinalg(xgTransposer)
      end if
    case (STATE_COLSROWS)
      if ( xgTransposer%state == STATE_COLSROWS ) then
        MSG_WARNING("Array colsrows has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_LINALG)%size == 1 ) then
        call xgBlock_copy(xgTransposer%xgBlock_linalg,xgTransposer%xgBlock_colsrows)
      else
        call xgTransposer_toColsRows(xgTransposer)
      end if
    end select

  end subroutine xgTransposer_transpose
!!***

!!****f* m_xgTransposer/xgTransposer_toLinalg
!!
!! NAME
!! xgTransposer_toLinalg

  subroutine xgTransposer_toLinalg(xgTransposer)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_toLinalg'
!End of the abilint section

   type(xgTransposer_t), intent(inout) :: xgTransposer
   double precision, allocatable, target :: sendbuf(:,:)
   double precision, pointer :: recvbuf(:,:)
   double precision, pointer :: buffer(:,:)
   integer, allocatable :: sendcounts(:), recvcounts(:)
   integer, allocatable :: sdispls(:), rdispls(:)
   integer :: ncpu, comm, me
   integer :: icplx
   integer :: i
   integer :: nrowsTotal, ncolsInd
   integer :: col, icpu, shiftCpu, bufiCol
   integer :: myrequest
   integer, allocatable :: nrowsInd(:)
   integer, allocatable :: request(:), status(:)
   type(xgBlock_t) :: xgBlock_toTransposed
   type(ptr_t), allocatable :: sendptrbuf(:)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me   = xgTransposer%mpiData(MPI_COLS)%rank
   
   select case ( space(xgTransposer%xgBlock_colsrows) )
     case (SPACE_CR, SPACE_R)
       icplx = 2
     case (SPACE_C)
       icplx = 1
     case default
       MSG_ERROR("Space value unknown !")
   end select
   
   ABI_MALLOC(nrowsInd,(ncpu))
   nrowsTotal = rows(xgTransposer%xgBlock_linalg)/icplx

   call xmpi_allgather(nrowsTotal,nrowsInd,comm,i)
   if ( i /= MPI_SUCCESS ) then
     MSG_ERROR("Error while gathering number of rows")
   end if 
   nrowsTotal = sum(nrowsInd)
   ncolsInd = cols(xgTransposer%xgBlock_colsrows)

   if ( nrowsTotal /= rows(xgTransposer%xgBlock_colsrows)/icplx ) then
     MSG_ERROR("Missmatch number of rows")
   end if
   if ( ncolsInd*ncpu /= cols(xgTransposer%xgBlock_linalg) ) then
     MSG_ERROR("Missmatch number of rows")
   end if

   ABI_MALLOC(sendbuf,(2,nrowsTotal*ncolsInd))
   call xgBlock_reverseMap(xgTransposer%xgBlock_colsrows,buffer,icplx,nrowsTotal*ncolsInd)
!$omp parallel do private(shiftCpu), collapse(2)
   do col = 1, ncolsInd
     do icpu = 1, ncpu
       shiftCpu = ncolsInd*sum(nrowsInd(1:icpu-1))
       sendbuf(:,(shiftCpu+(col-1)*nrowsInd(icpu)+1):(shiftCpu+col*nrowsInd(icpu))) = &

       buffer(:,((col-1)*nrowsTotal+sum(nrowsInd(1:icpu-1))+1):((col-1)*nrowsTotal+sum(nrowsInd(1:icpu))))
     end do
   end do

   ABI_MALLOC(recvcounts,(ncpu))
   ABI_MALLOC(rdispls,(ncpu))
   recvcounts(:) = 2*nrowsInd(me+1)*ncolsInd !! Thank you fortran for not starting at 0 !
   rdispls(1) = 0
   do i = 2, ncpu
     rdispls(i) = rdispls(i-1)+recvcounts(i-1)
   end do

   call xgBlock_reverseMap(xgTransposer%xgBlock_linalg,recvbuf,icplx,cols(xgTransposer%xgBlock_linalg)*nrowsInd(me+1))
   select case(xgTransposer%mpiAlgo)
   case (TRANS_ALL2ALL)
     ABI_MALLOC(sendcounts,(ncpu))
     ABI_MALLOC(sdispls,(ncpu))
     ABI_MALLOC(request,(1))
     myrequest = 1

     sendcounts(:) = 2*nrowsInd(:)*ncolsInd 
     sdispls(1) = 0
     do i = 2, ncpu
       sdispls(i) = sdispls(i-1)+sendcounts(i-1)
     end do

     !call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
     !                    recvbuf, recvcounts, rdispls, & 
     !                    comm, i)
     call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
                         recvbuf, recvcounts, rdispls, & 
                         comm, request(myrequest))
   case (TRANS_GATHER)
     ABI_MALLOC(request,(ncpu))
     myrequest = me+1

     ABI_ALLOCATE(sendptrbuf,(1:ncpu))
     do icpu = 1, ncpu
       sendptrbuf(me+1)%ptr => sendbuf(:,(icpu-1)*ncolsInd*nrowsInd(icpu)+1:icpu*ncolsInd*nrowsInd(icpu+1))
       !call xmpi_gatherv(sendptrbuf(me+1)%ptr,2*ncolsInd*nrowsInd(icpu),recvbuf,recvcounts,rdispls,icpu-1,comm,i)
       call mpi_igatherv(sendptrbuf(me+1)%ptr,2*ncolsInd*nrowsInd(icpu),MPI_DOUBLE_PRECISION,&
         recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,icpu-1,comm,request(icpu),i)
     end do

   case default
     MSG_BUG("This algo does not exist")
   end select

   ABI_FREE(nrowsInd)

   if ( cols(xgTransposer%xgBlock_colsrows)*rows(xgTransposer%xgBlock_colsrows)/icplx /= nrowsTotal*ncolsInd ) then
     MSG_ERROR("Wrong size for xgBlock linalg")
   end if

   xgTransposer%state = STATE_LINALG


   ABI_MALLOC(status,(MPI_STATUS_SIZE))
   call mpi_wait(request(myrequest),status,i)
   if ( i /= MPI_SUCCESS ) then
     MSG_ERROR("Error while waiting for mpi")
   end if 

   if ( allocated(sendcounts) ) then
     ABI_FREE(sendcounts)
   end if
   if ( allocated(sdispls) ) then
     ABI_FREE(sdispls)
   end if

   ABI_FREE(recvcounts)
   ABI_FREE(rdispls)


   if ( allocated(sendptrbuf) ) then
     ABI_FREE(sendptrbuf)
   end if

   do icpu = 1, size(request)
     if ( icpu /= myrequest ) then
       call mpi_wait(request(icpu),status,i)
       if ( i /= MPI_SUCCESS ) then
         MSG_ERROR("Error while waiting for other mpi")
       end if 
     end if
   end do
   ABI_FREE(sendbuf)
   ABI_FREE(status)
   ABI_FREE(request)

  end subroutine xgTransposer_toLinalg
!!***

!!****f* m_xgTransposer/xgTransposer_toColsRows
!!
!! NAME
!! xgTransposer_toColsRows

  subroutine xgTransposer_toColsRows(xgTransposer)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_toColsRows'
!End of the abilint section

   type(xgTransposer_t), intent(inout) :: xgTransposer
   double precision, pointer :: sendbuf(:,:)
   double precision, allocatable :: recvbuf(:,:)
   double precision, allocatable :: buffer(:,:)
   integer, allocatable :: sendcounts(:), recvcounts(:)
   integer, allocatable :: sdispls(:), rdispls(:)
   integer :: ncpu, comm, me
   integer :: icplx
   integer :: i
   integer :: nrowsTotal, ncolsInd
   integer :: col, icpu, shiftCpu, bufiCol
   integer :: licpu
   integer :: myrequest
   integer, allocatable :: nrowsInd(:)
   integer, allocatable :: request(:), status(:)
   type(xgBlock_t) :: xgBlock_toTransposed
   type(ptr_t), allocatable :: sendptrbuf(:)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me   = xgTransposer%mpiData(MPI_COLS)%rank
   
   select case ( space(xgTransposer%xgBlock_linalg) )
     case (SPACE_CR, SPACE_R)
       icplx = 2
     case (SPACE_C)
       icplx = 1
     case default
       MSG_ERROR("Space value unknown !")
   end select
   
   ABI_MALLOC(nrowsInd,(xgTransposer%mpiData(MPI_LINALG)%size))
   nrowsTotal = rows(xgTransposer%xgBlock_linalg)/icplx !number of pair of reals

   call xmpi_allgather(nrowsTotal,nrowsInd,xgTransposer%mpiData(MPI_LINALG)%comm,i)
   if ( i /= MPI_SUCCESS ) then
     MSG_ERROR("Error while gathering number of rows")
   end if 
   icpu = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu
   nrowsTotal = sum(nrowsInd(icpu+1:icpu+ncpu))
   ncolsInd = cols(xgTransposer%xgBlock_linalg)/xgTransposer%mpiData(MPI_COLS)%size
   write(*,*) xgTransposer%mpiData(MPI_ROWS)%rank, "in transpose to col_row there is a total of ", nrowsTotal, "real in rows to gather"
   write(*,*) nrowsInd

   ABI_MALLOC(recvbuf,(2,nrowsTotal*ncolsInd))
   ABI_MALLOC(recvcounts,(ncpu))
   ABI_MALLOC(rdispls,(ncpu))

   recvcounts(:) = 2*nrowsInd(icpu+1:icpu+ncpu)*ncolsInd
   rdispls(1) = 0
   do icpu = 2, ncpu
     rdispls(icpu) = rdispls(icpu-1)+recvcounts(icpu-1)
   end do

   select case(xgTransposer%mpiAlgo)
   case (TRANS_ALL2ALL)
     ABI_MALLOC(sendcounts,(ncpu))
     ABI_MALLOC(sdispls,(ncpu))
     ABI_MALLOC(request,(1))
     myrequest = 1

     sendcounts(:) = 2*nrowsInd(xgTransposer%mpiData(MPI_LINALG)%rank+1)*ncolsInd !! Thank you fortran for not starting at 0 !
     sdispls(1) = 0
     do i = 2, ncpu
       sdispls(i) = sdispls(i-1)+sendcounts(i-1)
     end do

     call xgBlock_reverseMap(xgTransposer%xgBlock_linalg,sendbuf,icplx,cols(xgTransposer%xgBlock_linalg)*nrowsInd(me+1))
     !call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
     !                    recvbuf, recvcounts, rdispls, & 
     !                    comm, i)
     call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
                         recvbuf, recvcounts, rdispls, & 
                         comm, request(myrequest))

   case (TRANS_GATHER)
     ABI_MALLOC(request,(ncpu))
     myrequest = me+1

     ABI_ALLOCATE(sendptrbuf,(1:ncpu))
     do icpu = 0, ncpu-1
       call xgBlock_setBlock(xgTransposer%xgBlock_linalg,xgBlock_toTransposed,icpu*ncolsInd+1,rows(xgTransposer%xgBlock_linalg),ncolsInd)
       call xgBlock_reverseMap(xgBlock_toTransposed,sendptrbuf(me+1)%ptr,icplx,ncolsInd*nrowsInd(me+1))
       !call xmpi_gatherv(sendptrbuf(me+1)%ptr,2*ncolsInd*nrowsInd(me+1),recvbuf,recvcounts,rdispls,icpu,comm,i)
       call mpi_igatherv(sendptrbuf(me+1)%ptr,2*ncolsInd*nrowsInd(me+1),MPI_DOUBLE_PRECISION,&
         recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,icpu,comm,request(icpu+1),i)
     end do

   case default
     MSG_BUG("This algo does not exist")
   end select

!call xgBlock_map(xgeigen,l_gvnlc,space,l_icplx*l_npw*l_nspinor,nband/dtset%npband,l_mpi_enreg%comm_fft)
!   call xgBlock_reverseMap(xgTransposer%xgBlock_colsrows,buffer,icplx,nrowsTotal*ncolsInd)
   xgTransposer%state = STATE_COLSROWS


   ABI_MALLOC(status,(MPI_STATUS_SIZE))
   call mpi_wait(request(myrequest),status,i)
   if ( i /= MPI_SUCCESS ) then
     MSG_ERROR("Error while waiting for mpi")
   end if 

   licpu = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu
   ABI_MALLOC(buffer,(2,nrowsTotal*ncolsInd))
   !!!!!!$omp parallel do private(shiftCpu), collapse(2)
   do col = 1, ncolsInd
     do icpu = 1, ncpu
       shiftCpu = ncolsInd*sum(nrowsInd(licpu+1:licpu+icpu-1))

       buffer(:,((col-1)*nrowsTotal+sum(nrowsInd(licpu+1:licpu+icpu-1))+1):((col-1)*nrowsTotal+sum(nrowsInd(licpu+1:licpu+icpu)))) = &
       recvbuf(:,(shiftCpu+(col-1)*nrowsInd(licpu+icpu)+1):(shiftCpu+col*nrowsInd(licpu+icpu)))
     end do
   end do
   call xgBlock_map(xgTransposer%xgBlock_colsrows,buffer,space(xgTransposer%xgBlock_linalg),icplx*nrowsTotal,ncolsInd,xgTransposer%mpiData(MPI_ROWS)%comm)

   ABI_FREE(nrowsInd)

   if ( allocated(sendcounts) ) then
     ABI_FREE(sendcounts)
   end if
   if ( allocated(sdispls) ) then
     ABI_FREE(sdispls)
   end if

   ABI_FREE(recvcounts)
   ABI_FREE(rdispls)

   ABI_FREE(recvbuf)

   if ( allocated(sendptrbuf) ) then
     ABI_FREE(sendptrbuf)
   end if

   do icpu = 1, size(request)
     if ( icpu /= myrequest ) then
       call mpi_wait(request(icpu),status,i)
       if ( i /= MPI_SUCCESS ) then
         MSG_ERROR("Error while waiting for other mpi")
       end if 
     end if
   end do
   ABI_FREE(status)
   ABI_FREE(request)

  end subroutine xgTransposer_toColsRows
!!***


!!****f* m_xgTransposer/xgTransposer_free
!!
!! NAME
!! xgTransposer_free

  subroutine xgTransposer_free(xgTransposer)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'xgTransposer_free'
!End of the abilint section

    type(xgTransposer_t), intent(inout) :: xgTransposer
#ifdef HAVE_MPI
    integer :: i
    call mpi_comm_free(xgTransposer%mpiData(MPI_ROWS)%comm,i)
    call mpi_comm_free(xgTransposer%mpiData(MPI_COLS)%comm,i)
    call mpi_comm_free(xgTransposer%mpiData(MPI_2DCART)%comm,i)
#endif

    if ( allocated(xgTransposer%lookup) ) then
      ABI_FREE(xgTransposer%lookup)
    end if


  end subroutine xgTransposer_free
!!***


end module m_xgTransposer
!!***

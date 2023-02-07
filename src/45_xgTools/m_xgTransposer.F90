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
!!  Copyright (C) 2017-2022 ABINIT group (J. Bieder)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_xgTransposer

  use, intrinsic :: iso_c_binding, only: c_double, c_size_t, c_loc

  use defs_basis, only : std_err, std_out, dp
  use m_xomp
  use m_profiling_abi
  use m_xmpi
  use m_errors
  use m_xg
  use m_time

#if defined HAVE_YAKL
  use gator_mod
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize, gpu_data_prefetch_async
#endif

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
  integer, parameter, public :: MPI_LINALG = 1
  integer, parameter, public :: MPI_ROWS   = 2
  integer, parameter, public :: MPI_COLS   = 3
  integer, parameter, public :: MPI_2DCART = 4
  integer, parameter, public :: TRANS_ALL2ALL = 1
  integer, parameter, public :: TRANS_GATHER = 2
  integer, parameter         :: TRANS_TYPE_CONSTRUCTED = 1
  integer, parameter         :: TRANS_TYPE_COPIED = 2

  integer, parameter :: tim_toColsRows  = 1662
  integer, parameter :: tim_toLinalg    = 1663
  integer, parameter :: tim_all2allv    = 1664
  integer, parameter :: tim_gatherv     = 1665
  integer, parameter :: tim_reorganize  = 1666
  integer, parameter :: tim_init        = 1667
  integer, parameter :: tim_free        = 1668
  integer, parameter :: tim_transpose   = 1669

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
    integer, pointer, private :: nrowsLinalg(:) => null()
    integer :: nrowsColsRows
    integer :: ncolsColsRows
    integer :: mpiAlgo
    integer :: type
    integer :: perPair
    integer :: use_gpu = 0
    integer :: gpu_num_openmp_threads = 1
#if defined HAVE_GPU && defined HAVE_YAKL
    real(kind=c_double), ABI_CONTIGUOUS pointer:: buffer(:,:) => null()
#else
    double precision, allocatable :: buffer(:,:)
#endif
  end type xgTransposer_t

  public :: xgTransposer_constructor
  public :: xgTransposer_copyConstructor
  public :: xgTransposer_transpose
  public :: xgTransposer_getRank
  public :: xgTransposer_getComm
  public :: xgTransposer_free


  contains
!!***

!!****f* m_xgTransposer/xgTransposer_constructor
!!
!! NAME
!! xgTransposer_constructor

  subroutine xgTransposer_constructor(xgTransposer,xgBlock_linalg,xgBlock_colsrows,ncpuRows,ncpuCols,state,algo)

    type(xgTransposer_t)   , intent(inout) :: xgTransposer
    type(xgBlock_t), target, intent(in   ) :: xgBlock_linalg
    type(xgBlock_t), target, intent(in   ) :: xgBlock_colsrows
    integer                , intent(in   ) :: ncpuRows
    integer                , intent(in   ) :: ncpuCols
    integer                , intent(in   ) :: state
    integer                , intent(in   ) :: algo
    integer :: commLinalg
    integer :: ncols
    integer :: nrows
    integer :: ierr
    integer :: icol
    character(len=500) :: message
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)

    xgTransposer%type = TRANS_TYPE_CONSTRUCTED

    xgTransposer%xgBlock_linalg => xgBlock_linalg
    xgTransposer%xgBlock_colsrows => xgBlock_colsrows
    xgTransposer%state = state
    commLinalg = comm(xgBlock_linalg)
    xgTransposer%mpiData(MPI_LINALG)%comm = commLinalg
    xgTransposer%mpiData(MPI_LINALG)%rank = xmpi_comm_rank(commLinalg)
    xgTransposer%mpiData(MPI_LINALG)%size = xmpi_comm_size(commLinalg)

    if ( xgTransposer%mpiData(MPI_LINALG)%size < ncpuCols*ncpuRows ) then
      write(message,'(a,i6,a,i6,a)') "There is not enough MPI processes in the communcation (", &
        xgTransposer%mpiData(MPI_LINALG)%size, "). Need at least ", ncpuCols*ncpuRows, " processes"
      ABI_ERROR(message)
    end if

    if ( ( algo == TRANS_ALL2ALL .or. algo == TRANS_GATHER ) ) then
      xgTransposer%mpiAlgo = algo
    else
      xgTransposer%mpiAlgo = TRANS_ALL2ALL
      ABI_COMMENT("Bad value for transposition MPI_algo. Will use ALLTOALL")
    end if

    !if ( xgTransposer%mpiAlgo == TRANS_ALL2ALL ) then
    !  ABI_COMMENT("Using mpi_alltoall for transposition")
    !else
    !  ABI_COMMENT("Using mpi_gatherv for transposition")
    !end if

    select case (state)
    case (STATE_LINALG)
      ! We are in the linalg representation.
      ! We need to construct the colsrows parallelization
      call xgBlock_getSize(xgBlock_linalg,nrows,ncols)
      call xmpi_sum(nrows,commLinalg,ierr)

      if ( MOD(ncols,ncpuCols) /=0 ) then
        if ( ncols > ncpuCols ) then
          write(message,'(a,i6,a,i6,a)') "Unbalanced parallelization : ", ncols, " columns for ", ncpuCols, " MPI"
          ABI_ERROR(message)
        else
          write(message,'(i6,a)') (ncpuCols-ncols)*ncpuRows, " MPI will not be used  because of the number of columns!!"
          ABI_ERROR(message)
          !ncpuCols = ncols
        end if
      end if

      select case ( space(xgBlock_linalg) )
        case (SPACE_CR, SPACE_R)
          xgTransposer%perPair = 2
        case (SPACE_C)
          xgTransposer%perPair = 1
        case default
          ABI_ERROR("Space value unknown !")
      end select

      !write(*,*) "There is a total of ", nrows/xgTransposer%perPair, "rows of real pairs for ", ncpuRows, "fft cpu"

      ! Build the lookup table
      ABI_MALLOC(xgTransposer%lookup,(1:ncols))
      !ABI_MALLOC(xgTransposer%me_g0_lookup,())
      do icol = 0, ncols-1
        xgTransposer%lookup(icol+1) = MOD(icol,ncpuCols)
      end do

      call xgTransposer_makeComm(xgTransposer,ncpuRows,ncpuCols)
      call xgTransposer_computeDistribution(xgTransposer)
      call xgTransposer_makeXgBlock(xgTransposer)

    case (STATE_COLSROWS)
      ABI_BUG("Not yet implemented")
    case default
      ABI_ERROR("State is undefined")
    end select

    call timab(tim_init,2,tsec)

  end subroutine xgTransposer_constructor
!!***

!!****f* m_xgTransposer/xgTransposer_copyConstructor
!!
!! NAME
!! xgTransposer_copyConstructor

  subroutine xgTransposer_copyConstructor(xgTransposer,xgTransposerInitialized,xgBlock_linalg,xgBlock_colsrows,state)

    type(xgTransposer_t)   , intent(inout) :: xgTransposer
    type(xgTransposer_t)   , intent(in   ) :: xgTransposerInitialized
    type(xgBlock_t), target, intent(in   ) :: xgBlock_linalg
    type(xgBlock_t), target, intent(in   ) :: xgBlock_colsrows
    integer                , intent(in   ) :: state
    integer :: commLinalg
    integer :: ncols, ncpuCols
    integer :: nrows, ncpuRows
    integer :: ierr
    integer :: icol
    character(len=500) :: message
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)

    xgTransposer%type = TRANS_TYPE_COPIED

    xgTransposer%xgBlock_linalg => xgBlock_linalg
    xgTransposer%xgBlock_colsrows => xgBlock_colsrows
    xgTransposer%state = state
    commLinalg = comm(xgBlock_linalg)

    if ( commLinalg /= xgTransposerInitialized%mpiData(MPI_LINALG)%comm ) then
      ABI_ERROR("Linalg communicators are different for the two transposers, this is not allowed.")
    end if

    xgTransposer%mpiData(:)%comm = xgTransposerInitialized%mpiData(:)%comm
    xgTransposer%mpiData(MPI_LINALG)%rank = xmpi_comm_rank(commLinalg)
    xgTransposer%mpiData(MPI_LINALG)%size = xmpi_comm_size(commLinalg)
    xgTransposer%mpiData(MPI_COLS)%rank = xmpi_comm_rank(xgTransposer%mpiData(MPI_COLS)%comm)
    xgTransposer%mpiData(MPI_COLS)%size = xmpi_comm_size(xgTransposer%mpiData(MPI_COLS)%comm)
    xgTransposer%mpiData(MPI_ROWS)%rank = xmpi_comm_rank(xgTransposer%mpiData(MPI_ROWS)%comm)
    xgTransposer%mpiData(MPI_ROWS)%size = xmpi_comm_size(xgTransposer%mpiData(MPI_ROWS)%comm)

    xgTransposer%mpiAlgo = xgTransposerInitialized%mpiAlgo

    ncpuCols = xgTransposer%mpiData(MPI_COLS)%size
    ncpuRows = xgTransposer%mpiData(MPI_ROWS)%size

    select case (state)
    case (STATE_LINALG)
      call xgBlock_getSize(xgBlock_linalg,nrows,ncols)
      call xmpi_sum(nrows,commLinalg,ierr)

      if ( MOD(ncols,ncpuCols) /=0 ) then
        if ( ncols > ncpuCols ) then
          write(message,'(a,i6,a,i6,a)') "Unbalanced parallelization : ", ncols, " columns for ", ncpuCols, " MPI"
          ABI_ERROR(message)
        else
          write(message,'(i6,a)') (ncpuCols-ncols)*ncpuRows, " MPI will not be used  because of the number of columns!!"
          ABI_ERROR(message)
        end if
      end if

      select case ( space(xgBlock_linalg) )
        case (SPACE_CR, SPACE_R)
          xgTransposer%perPair = 2
        case (SPACE_C)
          xgTransposer%perPair = 1
        case default
          ABI_ERROR("Space value unknown !")
      end select

      !write(*,*) "There is a total of ", nrows/xgTransposer%perPair, "rows of real pairs for ", ncpuRows, "fft cpu"

      ! Build the lookup table
      ABI_MALLOC(xgTransposer%lookup,(1:ncols))
      if ( cols(xgTransposerInitialized%xgBlock_linalg) /= ncols ) then
        do icol = 0, ncols-1
        xgTransposer%lookup(icol+1) = MOD(icol,ncpuCols)
        end do
      else
        xgTransposer%lookup(:) = xgTransposerInitialized%lookup(:)
      endif

      call xgTransposer_computeDistribution(xgTransposer)
      call xgTransposer_makeXgBlock(xgTransposer)

    case (STATE_COLSROWS)
      ABI_BUG("Not yet implemented")
    case default
      ABI_ERROR("State is undefined")
    end select

    call timab(tim_init,2,tsec)

  end subroutine xgTransposer_copyConstructor
!!***

!!****f* m_xgTransposer/xgTransposer_makeComm
!!
!! NAME
!! xgTransposer_makeComm

  subroutine xgTransposer_makeComm(xgTransposer,ncpuRows,ncpuCols)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer             , intent(in   ) :: ncpuRows
    integer             , intent(in   ) :: ncpuCols
    integer :: commColsRows
#if defined HAVE_MPI
    integer :: sizeGrid(2) !coordInGrid(2),
    logical :: periodic(2), selectDim(2), reorder
    integer :: ierr

    sizeGrid(1) = ncpuRows
    sizeGrid(2) = ncpuCols
    periodic = (/ .false., .false. /)
    reorder  = .false.
    call mpi_cart_create(xgTransposer%mpiData(MPI_LINALG)%comm,2,sizeGrid,periodic,reorder,commColsRows,ierr)
    if ( ierr /= xmpi_success ) then
      ABI_ERROR("xgTransposer failed to creat cartesian grid")
    end if

    selectDim = (/ .true., .false. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_ROWS)%comm,ierr)
    if ( ierr /= xmpi_success ) then
      ABI_ERROR("xgTransposer failed to creat rows communicator")
    end if
    selectDim = (/ .false., .true. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_COLS)%comm,ierr)
    if ( ierr /= xmpi_success ) then
      ABI_ERROR("xgTransposer failed to creat columns communicator")
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

!!****f* m_xgTransposer/xgTransposer_computeDistribution
!!
!! NAME
!! xgTransposer_computeDistribution

  subroutine xgTransposer_computeDistribution(xgTransposer)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer :: nRealPairs
    integer :: ierr
    integer :: icpu
    integer :: ncpuCols

    ABI_MALLOC(xgTransposer%nrowsLinalg,(xgTransposer%mpiData(MPI_LINALG)%size))
    nRealPairs = rows(xgTransposer%xgBlock_linalg)/xgTransposer%perPair !number of pair of reals

    call xmpi_allgather(nRealPairs,xgTransposer%nrowsLinalg,xgTransposer%mpiData(MPI_LINALG)%comm,ierr)
    if ( ierr /= xmpi_success ) then
      ABI_ERROR("Error while gathering number of rows in linalg")
    end if

    ncpuCols = xgTransposer%mpiData(MPI_COLS)%size
    icpu = xgTransposer%mpiData(MPI_ROWS)%rank*ncpuCols
    xgTransposer%nrowsColsRows = sum(xgTransposer%nrowsLinalg(icpu+1:icpu+ncpuCols))
    xgTransposer%ncolsColsRows = cols(xgTransposer%xgBlock_linalg)/ncpuCols

    !write(*,*) "In linalg, # of real pairs:", xgTransposer%nrowsLinalg
    !write(*,*) "In rows, # of real pairs for proc ", xgTransposer%mpiData(MPI_ROWS)%rank, ":", xgTransposer%nrowsColsRows
    !write(*,*) "In cols, # of cols for proc ", xgTransposer%mpiData(MPI_COLS)%rank, ":", xgTransposer%ncolsColsRows

  end subroutine xgTransposer_computeDistribution
!!***

!!****f* m_xgTransposer/xgTransposer_makeXgBlock
!!
!! NAME
!! xgTransposer_makeXgBlock

  subroutine xgTransposer_makeXgBlock(xgTransposer)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    !integer :: cols, rows

    select case (xgTransposer%state)
    case (STATE_LINALG)
      ! Assume xgBlock_colsrows is empty and not constructed because user cannot
      ! predict the size
#if defined HAVE_GPU && defined HAVE_YAKL
      if ( associated(xgTransposer%buffer) ) then
        ABI_FREE_MANAGED(xgTransposer%buffer)
      end if
#else
      if ( allocated(xgTransposer%buffer) ) then
        ABI_FREE(xgTransposer%buffer)
      end if
#endif

      if ( xgTransposer%mpiData(MPI_COLS)%size == 1 ) then
        xgTransposer%xgBlock_colsrows = xgTransposer%xgBlock_linalg
      else
#if defined HAVE_GPU && defined HAVE_YAKL
        ABI_MALLOC_MANAGED(xgTransposer%buffer,(/2,xgTransposer%ncolsColsRows*xgTransposer%nrowsColsRows/))
#else
        ABI_MALLOC(xgTransposer%buffer,(2,xgTransposer%ncolsColsRows*xgTransposer%nrowsColsRows))
#endif
        call xgBlock_map(xgTransposer%xgBlock_colsrows,xgTransposer%buffer,space(xgTransposer%xgBlock_linalg),&
          xgTransposer%perPair*xgTransposer%nrowsColsRows,&
          xgTransposer%ncolsColsRows,xgTransposer%mpiData(MPI_ROWS)%comm)
      end if
    case (STATE_COLSROWS)
      ABI_ERROR("Not yet implemented")
    case default
      ABI_ERROR("State unknown")
    end select
  end subroutine xgTransposer_makeXgBlock
!!***

!!****f* m_xgTransposer/xgTransposer_transpose
!!
!! NAME
!! xgTransposer_transpose

  subroutine xgTransposer_transpose(xgTransposer,toState)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer             , intent(in   ) :: toState
    double precision :: tsec(2)

    call timab(tim_transpose,1,tsec)

    if ( toState /= STATE_LINALG .and. toState /= STATE_COLSROWS ) then
      ABI_ERROR("Bad value for toState")
    end if

    !write(std_out,*) "linalg", rows(xgTransposer%xgBlock_linalg)*cols(xgTransposer%xgBlock_linalg)
    !write(std_out,*) "colsrows", rows(xgTransposer%xgBlock_colsrows)*cols(xgTransposer%xgBlock_colsrows)
    select case (toState)
    case (STATE_LINALG)
      if ( xgTransposer%state == STATE_LINALG ) then
        ABI_WARNING("Array linalg has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_COLS)%size > 1 ) then
        call xgTransposer_toLinalg(xgTransposer)
      else
        xgTransposer%state = STATE_LINALG
      end if
    case (STATE_COLSROWS)
      if ( xgTransposer%state == STATE_COLSROWS ) then
        ABI_WARNING("Array colsrows has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_COLS)%size > 1 ) then
        call xgTransposer_toColsRows(xgTransposer)
      else
        xgTransposer%state = STATE_COLSROWS
      end if
    end select

    call timab(tim_transpose,2,tsec)

  end subroutine xgTransposer_transpose
!!***

!!****f* m_xgTransposer/xgTransposer_toLinalg
!!
!! NAME
!! xgTransposer_toLinalg

  subroutine xgTransposer_toLinalg(xgTransposer)

   type(xgTransposer_t), intent(inout) :: xgTransposer
   double precision, allocatable, target :: sendbuf(:,:)
   double precision, pointer :: recvbuf(:,:)
   !double precision, pointer :: buffer(:,:)
   integer, allocatable :: sendcounts(:), recvcounts(:)
   integer, allocatable :: sdispls(:), rdispls(:)
   integer :: ncpu, comm, me
   integer :: nrowsColsRows
   integer :: ncolsColsRows
   integer :: nrowsLinalgMe
   integer :: icpu, ierr
   !integer :: myrequest
   !integer, allocatable :: request(:), status(:)
   type(ptr_t), allocatable :: sendptrbuf(:)
   integer, pointer :: nrowsLinalg(:)
   double precision :: tsec(2)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
   double precision, allocatable :: recvbuf_mpi(:,:)
   integer(c_size_t) :: buffer_size
#endif

   call timab(tim_toLinalg,1,tsec)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu

   nrowsColsRows = xgTransposer%nrowsColsRows
   ncolsColsRows = xgTransposer%ncolsColsRows

   nrowsLinalg => xgTransposer%nrowsLinalg
   nrowsLinalgMe = nrowsLinalg(xgTransposer%mpiData(MPI_LINALG)%rank+1)

   ABI_MALLOC(sendbuf,(2,nrowsColsRows*ncolsColsRows))
   call xgTransposer_reorganizeData(xgTransposer,sendbuf)

   ABI_MALLOC(recvcounts,(ncpu))
   ABI_MALLOC(rdispls,(ncpu))
   recvcounts(:) = 2*nrowsLinalgMe*ncolsColsRows !! Thank you fortran for not starting at 0 !
   rdispls(1) = 0
   do icpu = 2, ncpu
     rdispls(icpu) = rdispls(icpu-1)+recvcounts(icpu-1)
   end do

   call xgBlock_reverseMap(xgTransposer%xgBlock_linalg,recvbuf,xgTransposer%perPair,cols(xgTransposer%xgBlock_linalg)*nrowsLinalgMe)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
   ! just for debug
   if( xgTransposer%use_gpu == 1) then
      !call gpu_managed_ptr_status(C_LOC(recvbuf))
      buffer_size = size(recvbuf) * dp
      call gpu_data_prefetch_async(C_LOC(recvbuf), buffer_size, CPU_DEVICE_ID)
      call gpu_device_synchronize()

      ABI_MALLOC(recvbuf_mpi, (size(recvbuf,1), size(recvbuf,2)) )

   end if
#endif

   ABI_MALLOC(sendcounts,(ncpu))
   ABI_MALLOC(sdispls,(ncpu))
   sendcounts(:) = 2*nrowsLinalg(me+1:me+ncpu)*ncolsColsRows
   sdispls(1) = 0
   do icpu = 2, ncpu
   sdispls(icpu) = sdispls(icpu-1)+sendcounts(icpu-1)
   end do

   select case(xgTransposer%mpiAlgo)
   case (TRANS_ALL2ALL)
     !ABI_MALLOC(request,(1))
     !myrequest = 1

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
     if( xgTransposer%use_gpu == 1) then

       call timab(tim_all2allv,1,tsec)
       call xmpi_alltoallv(sendbuf,     sendcounts, sdispls, &
                           recvbuf_mpi, recvcounts, rdispls, &
                           comm, ierr)
       call timab(tim_all2allv,2,tsec)
       !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
       !                    recvbuf, recvcounts, rdispls, &
       !                    comm, request(myrequest))

       ! copy back recvbuf_mpi into recvbuf
       recvbuf(:,:) = recvbuf_mpi(:,:)

       ABI_FREE(recvbuf_mpi)

     end if
#endif


    if( xgTransposer%use_gpu == 0) then

      call timab(tim_all2allv,1,tsec)
      call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
                          recvbuf, recvcounts, rdispls, &
                          comm, ierr)
      call timab(tim_all2allv,2,tsec)
      !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
      !                    recvbuf, recvcounts, rdispls, &
      !                    comm, request(myrequest))
   end if

   case (TRANS_GATHER)
     !ABI_MALLOC(request,(ncpu))
     me = xgTransposer%mpiData(MPI_COLS)%rank
     !myrequest = me+1

     ABI_MALLOC(sendptrbuf,(1:ncpu))
     do icpu = 1, ncpu
       sendptrbuf(me+1)%ptr => sendbuf(:,sdispls(icpu)/2+1:sdispls(icpu)/2+sendcounts(icpu))
       call timab(tim_gatherv,1,tsec)
       call xmpi_gatherv(sendptrbuf(me+1)%ptr,sendcounts(icpu),recvbuf,recvcounts,rdispls,icpu-1,comm,ierr)
       call timab(tim_gatherv,2,tsec)
       !call mpi_igatherv(sendptrbuf(me+1)%ptr,sendcounts(icpu),MPI_DOUBLE_PRECISION,&
       !  recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,icpu-1,comm,request(icpu),ierr)
     end do

   case default
     ABI_BUG("This algo does not exist")
   end select

   xgTransposer%state = STATE_LINALG

   !ABI_MALLOC(status,(MPI_STATUS_SIZE))
   !call mpi_wait(request(myrequest),status,ierr)
   if ( ierr /= xmpi_success ) then
     ABI_ERROR("Error while waiting for mpi")
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

   !do icpu = 1, size(request)
   !  if ( icpu /= myrequest ) then
   !    call mpi_wait(request(icpu),status,ierr)
   !    if ( ierr /= MPI_SUCCESS ) then
   !      ABI_ERROR("Error while waiting for other mpi")
   !    end if
   !  end if
   !end do
   ABI_FREE(sendbuf)
   !ABI_FREE(status)
   !ABI_FREE(request)

    call timab(tim_toLinalg,2,tsec)

  end subroutine xgTransposer_toLinalg
!!***

!!****f* m_xgTransposer/xgTransposer_toColsRows
!!
!! NAME
!! xgTransposer_toColsRows

  subroutine xgTransposer_toColsRows(xgTransposer)

   type(xgTransposer_t), intent(inout) :: xgTransposer
   double precision, pointer :: sendbuf(:,:)
   double precision, allocatable :: recvbuf(:,:)
   !double precision, allocatable :: buffer(:,:)
   integer, allocatable :: sendcounts(:), recvcounts(:)
   integer, allocatable :: sdispls(:), rdispls(:)
   integer :: ncpu, comm, me
   integer :: nrowsColsRows
   integer :: ncolsColsRows
   integer :: nrowsLinalgMe
   integer :: icpu,ierr
   !integer :: myrequest
   !integer, allocatable :: request(:), status(:)
   type(xgBlock_t) :: xgBlock_toTransposed
   type(ptr_t), allocatable :: sendptrbuf(:)
   integer, pointer :: nrowsLinalg(:)
   double precision :: tsec(2)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
   double precision, allocatable :: sendbuf_mpi(:,:)
   integer(c_size_t) :: buffer_size
#endif

    call timab(tim_toColsRows,1,tsec)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu

   nrowsColsRows = xgTransposer%nrowsColsRows
   ncolsColsRows = xgTransposer%ncolsColsRows

   nrowsLinalg => xgTransposer%nrowsLinalg
   nrowsLinalgMe = nrowsLinalg(xgTransposer%mpiData(MPI_LINALG)%rank+1)

   ABI_MALLOC(recvbuf,(2,nrowsColsRows*ncolsColsRows))
   ABI_MALLOC(recvcounts,(ncpu))
   ABI_MALLOC(rdispls,(ncpu))

   recvcounts(:) = 2*nrowsLinalg(me+1:me+ncpu)*ncolsColsRows
   rdispls(1) = 0
   do icpu = 2, ncpu
     rdispls(icpu) = rdispls(icpu-1)+recvcounts(icpu-1)
   end do

   select case(xgTransposer%mpiAlgo)
   case (TRANS_ALL2ALL)
     ABI_MALLOC(sendcounts,(ncpu))
     ABI_MALLOC(sdispls,(ncpu))
     !ABI_MALLOC(request,(1))
     !myrequest = 1

     sendcounts(:) = 2*nrowsLinalgMe*ncolsColsRows !! Thank you fortran for not starting at 0 !
     sdispls(1) = 0
     do icpu = 2, ncpu
       sdispls(icpu) = sdispls(icpu-1)+sendcounts(icpu-1)
     end do

     call xgBlock_reverseMap(xgTransposer%xgBlock_linalg,sendbuf, &
&      xgTransposer%perPair,cols(xgTransposer%xgBlock_linalg)*nrowsLinalgMe)
     !write(*,*) "Before ialltoall"

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
    ! if gpu is enabled, data are located in GPU memory, so we copy them on a host buffer
    if( xgTransposer%use_gpu == 1) then
      ABI_MALLOC(sendbuf_mpi, (size(sendbuf,1), size(sendbuf,2)) )

      ! sync sendbuf on host and then copy to sendbuf_mpi
      buffer_size = size(sendbuf) * dp
      call gpu_data_prefetch_async(C_LOC(sendbuf), buffer_size, CPU_DEVICE_ID)
      call gpu_device_synchronize()

      sendbuf_mpi(:,:) = sendbuf(:,:)

      call timab(tim_all2allv,1,tsec)
      ! this is a cpu mpi comm
      call xmpi_alltoallv(sendbuf_mpi, sendcounts, sdispls, &
                          recvbuf,     recvcounts, rdispls, &
                          comm, ierr)
      call timab(tim_all2allv,2,tsec)
      !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
      !                    recvbuf, recvcounts, rdispls, &
      !                    comm, request(myrequest))
      !write(*,*) "After ialltoall"

      ABI_FREE(sendbuf_mpi)
    end if
#endif

    if( xgTransposer%use_gpu == 0) then
      call timab(tim_all2allv,1,tsec)
      call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
                          recvbuf, recvcounts, rdispls, &
                          comm, ierr)
      call timab(tim_all2allv,2,tsec)
      !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
      !                    recvbuf, recvcounts, rdispls, &
      !                    comm, request(myrequest))
      !write(*,*) "After ialltoall"
    end if

   case (TRANS_GATHER)
     !ABI_MALLOC(request,(ncpu))
     me = xgTransposer%mpiData(MPI_COLS)%rank
     !myrequest = me+1

     ABI_MALLOC(sendptrbuf,(1:ncpu))
     !call flush(6)
     !call xmpi_barrier(xgTransposer%mpiData(MPI_LINALG)%comm)
     do icpu = 0, ncpu-1
       !write(*,*) me, "->", icpu, "from col ",icpu*ncolsColsRows+1, " number of rows:", nrowsLinalgMe
       call xgBlock_setBlock(xgTransposer%xgBlock_linalg,xgBlock_toTransposed,icpu*ncolsColsRows+1,nrowsLinalgMe,ncolsColsRows)
       call xgBlock_reverseMap(xgBlock_toTransposed,sendptrbuf(me+1)%ptr,xgTransposer%perPair,ncolsColsRows*nrowsLinalgMe)
       call timab(tim_gatherv,1,tsec)
       call xmpi_gatherv(sendptrbuf(me+1)%ptr,2*ncolsColsRows*nrowsLinalgMe,recvbuf,recvcounts,rdispls,icpu,comm,ierr)
       call timab(tim_gatherv,2,tsec)
       !call mpi_igatherv(sendptrbuf(me+1)%ptr,2*ncolsColsRows*nrowsLinalgMe,MPI_DOUBLE_PRECISION,&
       !  recvbuf,recvcounts,rdispls,MPI_DOUBLE_PRECISION,icpu,comm,request(icpu+1),ierr)
     end do
     !call xmpi_barrier(xgTransposer%mpiData(MPI_LINALG)%comm)
     !call flush(6)
     !write(*,*) me, request

   case default
     ABI_BUG("This algo does not exist")
   end select

   !ABI_MALLOC(status,(MPI_STATUS_SIZE))
   !call mpi_wait(request(myrequest),status,ierr)
   !write(*,*) "Request ended"
   if ( ierr /= xmpi_success ) then
     ABI_ERROR("Error while waiting for mpi")
   end if
   !write(*,*) "with success"

   call xgTransposer_reorganizeData(xgTransposer,recvbuf)

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

   xgTransposer%state = STATE_COLSROWS

   !do icpu = 1, size(request)
   !  if ( icpu /= myrequest ) then
   !    call mpi_wait(request(icpu),status,ierr)
   !    if ( ierr /= MPI_SUCCESS ) then
   !      ABI_ERROR("Error while waiting for other mpi")
   !    end if
   !  end if
   !end do
   !ABI_FREE(status)
   !ABI_FREE(request)

    call timab(tim_toColsRows,2,tsec)

  end subroutine xgTransposer_toColsRows
!!***

!!****f* m_xgTransposer/xgTransposer_reorganizeData
!!
!! NAME
!! xgTransposer_reorganizeData

  subroutine xgTransposer_reorganizeData(xgTransposer,bufferMess)

    type(xgTransposer_t), intent(inout)          :: xgTransposer
    double precision    , intent(inout), target  :: bufferMess(:,:)
    double precision,                    pointer :: bufferOrdered(:,:) => null()
    integer :: shiftCpu
    integer :: nrowsColsRows
    integer :: ncolsColsRows
    integer :: tos,toe,froms,frome
    integer :: col, icpu
    integer :: me
    integer :: nPair
    integer, pointer :: nrowsLinalg(:)
    double precision :: tsec(2)
#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
    integer(c_size_t) :: buffer_size
#endif

    call timab(tim_reorganize,1,tsec)

    me = xgTransposer%mpiData(MPI_ROWS)%rank*xgTransposer%mpiData(MPI_COLS)%size
    nrowsColsRows = xgTransposer%nrowsColsRows
    ncolsColsRows = xgTransposer%ncolsColsRows
    nPair = nrowsColsRows*ncolsColsRows

    call xgBlock_reverseMap(xgTransposer%xgBlock_colsrows,bufferOrdered,xgTransposer%perPair,nPair)
    nrowsLinalg => xgTransposer%nrowsLinalg

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
    ! if gpu is enabled, data are located in GPU memory, so we prefetch them on host
    ! to do the following reorganization.
    ! Alternatively, we should provide a GPU implementation of this data layout reorganization
    if( xgTransposer%use_gpu == 1) then
       buffer_size = size(bufferOrdered) * dp
       call gpu_data_prefetch_async(C_LOC(bufferOrdered), buffer_size, CPU_DEVICE_ID)
       call gpu_device_synchronize()
    end if

    ! if gpu enabled increase locally OpenMP num threads
    if (xgTransposer%use_gpu == 1) then
       call xomp_set_num_threads(xgTransposer%gpu_num_openmp_threads)
    end if
#endif

    select case (xgTransposer%state)
    case (STATE_LINALG)
      ! We are going to STATE_COLSROWS so we are after all2all
      !$omp parallel do private(shiftCpu,toe,tos,frome,froms), collapse(2)
      do col = 1, ncolsColsRows
        do icpu = 1, xgTransposer%mpiData(MPI_COLS)%size
          shiftCpu = ncolsColsRows*sum(nrowsLinalg(me+1:me+icpu-1))
          tos=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu-1))+1)
          toe=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu)))
          froms=(shiftCpu+(col-1)*nrowsLinalg(me+icpu)+1)
          frome=(shiftCpu+col*nrowsLinalg(me+icpu))
          bufferOrdered(:,tos:toe) = bufferMess(:,froms:frome)
        end do
      end do
    case (STATE_COLSROWS)
      ! We are going to STATE_LINALG so we are before all2all
      !$omp parallel do private(shiftCpu,toe,tos,frome,froms), collapse(2)
      do col = 1, ncolsColsRows
        do icpu = 1, xgTransposer%mpiData(MPI_COLS)%size
          shiftCpu = ncolsColsRows*sum(nrowsLinalg(me+1:me+icpu-1))
          tos=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu-1))+1)
          toe=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu)))
          froms=(shiftCpu+(col-1)*nrowsLinalg(me+icpu)+1)
          frome=(shiftCpu+col*nrowsLinalg(me+icpu))
          bufferMess(:,froms:frome) = bufferOrdered(:,tos:toe)
        end do
      end do
    end select

#if defined(HAVE_GPU_CUDA) && defined(HAVE_KOKKOS) && defined(HAVE_YAKL)
    ! if gpu enable restore OpenMP num threads to 1
    if (xgTransposer%use_gpu == 1) then
       ! restore OMP_NUM_THREADS=1
       call xomp_set_num_threads(1)
    end if

    ! if gpu is enabled, transfer back data on GPU
    if (xgTransposer%use_gpu == 1) then
       buffer_size = size(bufferOrdered) * dp
       call gpu_data_prefetch_async(C_LOC(bufferOrdered), buffer_size)
       call gpu_device_synchronize()
    end if
#endif

    call timab(tim_reorganize,2,tsec)

  end subroutine xgTransposer_reorganizeData
!!***

!!****f* m_xgTransposer/xgTransposer_getRank
!!
!! NAME
!! xgTransposer_getRank

  function xgTransposer_getRank(xgTransposer, comm) result(rank)
    type(xgTransposer_t), intent(in   ) :: xgTransposer
    integer             , intent(in   ) :: comm
    integer :: rank
    if ( (comm > ubound(xgTransposer%mpiData,1)) .or.  (comm < lbound(xgTransposer%mpiData,1)) ) then
      ABI_ERROR("Value for communicator is wrong")
    end if
    rank = xgTransposer%mpiData(comm)%rank
  end function xgTransposer_getRank
!!***
 
!!****f* m_xgTransposer/xgTransposer_getComm
!!
!! NAME
!! xgTransposer_getComm

  function xgTransposer_getComm(xgTransposer, comm1) result(communicator)
    type(xgTransposer_t), intent(in   ) :: xgTransposer
    integer             , intent(in   ) :: comm1
    integer :: communicator
    if ( (comm1 > ubound(xgTransposer%mpiData,1)) .or.  (comm1 < lbound(xgTransposer%mpiData,1)) ) then
      ABI_ERROR("Value for communicator is wrong")
    end if
    communicator = xgTransposer%mpiData(comm1)%comm
  end function xgTransposer_getComm
!!***

!!****f* m_xgTransposer/xgTransposer_free
!!
!! NAME
!! xgTransposer_free

  subroutine xgTransposer_free(xgTransposer)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    double precision :: tsec(2)
    integer :: i

    call timab(tim_free,1,tsec)
#ifdef HAVE_MPI
    if ( xgTransposer%type == TRANS_TYPE_CONSTRUCTED ) then
      call mpi_comm_free(xgTransposer%mpiData(MPI_ROWS)%comm,i)
      call mpi_comm_free(xgTransposer%mpiData(MPI_COLS)%comm,i)
      call mpi_comm_free(xgTransposer%mpiData(MPI_2DCART)%comm,i)
    end if
#else
    ABI_UNUSED(i)
#endif

    if ( allocated(xgTransposer%lookup) ) then
      ABI_FREE(xgTransposer%lookup)
    end if

    if ( associated(xgTransposer%nrowsLinalg) ) then
      ABI_FREE(xgTransposer%nrowsLinalg)
    end if

#if defined HAVE_GPU && defined HAVE_YAKL
    if ( associated(xgTransposer%buffer) ) then
      ABI_FREE_MANAGED(xgTransposer%buffer)
    end if
#else
    if ( allocated(xgTransposer%buffer) ) then
      ABI_FREE(xgTransposer%buffer)
    end if
#endif

    call timab(tim_free,2,tsec)

  end subroutine xgTransposer_free
!!***

end module m_xgTransposer
!!***

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
!!  Copyright (C) 2017-2019 ABINIT group (J. Bieder)
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
  use m_time

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
    !integer, allocatable, private :: me_g0_lookup(:,:)
    integer, pointer, private :: nrowsLinalg(:) => null()
    integer :: nrowsColsRows
    integer :: ncolsColsRows
    integer :: mpiAlgo
    integer :: perPair
    double precision, allocatable :: buffer(:,:)
    integer :: debug_rank
  end type xgTransposer_t

  public :: xgTransposer_init
  public :: xgTransposer_transpose
  public :: xgTransposer_free
  public :: xgTransposer_getCPURow


  contains
!!***

!!****f* m_xgTransposer/xgTransposer_init
!!
!! NAME
!! xgTransposer_init

  subroutine xgTransposer_init(xgTransposer,xgBlock_linalg,xgBlock_colsrows,ncpuRows,ncpuCols,state,algo, d_rank)

    type(xgTransposer_t)   , intent(inout) :: xgTransposer
    type(xgBlock_t), target, intent(in   ) :: xgBlock_linalg
    type(xgBlock_t), target, intent(in   ) :: xgBlock_colsrows
    integer                , intent(in   ) :: ncpuRows
    integer                , intent(in   ) :: ncpuCols
    integer                , intent(in   ) :: state
    integer                , intent(in   ) :: algo
    integer                , intent(in   ) :: d_rank
    integer :: commLinalg
    integer :: ncols
    integer :: nrows
    integer :: ierr
    integer :: icol
    character(len=500) :: message
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)

    xgTransposer%xgBlock_linalg => xgBlock_linalg
    xgTransposer%xgBlock_colsrows => xgBlock_colsrows
    xgTransposer%state = state
    xgTransposer%debug_rank = d_rank
    commLinalg = comm(xgBlock_linalg)
    !print *, "commLinalg", commLinalg
    !stop
    xgTransposer%mpiData(MPI_LINALG)%comm = commLinalg
    xgTransposer%mpiData(MPI_LINALG)%rank = xmpi_comm_rank(commLinalg)
    xgTransposer%mpiData(MPI_LINALG)%size = xmpi_comm_size(commLinalg)
!    if (xmpi_comm_rank(xmpi_world) == 0) then
!      print *, "xcgColsRows UNUTAR TRANSPOSERA:"
!      call xgBlock_print(xgBlock_colsrows,6)  
!    end if
!    if ( space(xgBlock_linalg) /= space(xgBlock_colsrows) ) then
!      MSG_ERROR("Linalg xgBlock and ColsRows xgBlocks are not in the same space")
!    end if
!    print *, "ncpuRows", ncpuRows
!    print *, "ncpuCols", ncpuCols 
!    print *, "xgTransposer%mpiData(MPI_LINALG)%rank", xgTransposer%mpiData(MPI_LINALG)%rank
!    print *, "xgTransposer%mpiData(MPI_LINALG)%size", xgTransposer%mpiData(MPI_LINALG)%size   
!    stop
    
    if ( xgTransposer%mpiData(MPI_LINALG)%size < ncpuCols*ncpuRows ) then
      write(message,'(a,i6,a,i6,a)') "There is not enough MPI processes in the communcation (", &
        xgTransposer%mpiData(MPI_LINALG)%size, "). Need at least ", ncpuCols*ncpuRows, " processes"
      MSG_ERROR(message)
    end if

    if ( ( algo == TRANS_ALL2ALL .or. algo == TRANS_GATHER ) ) then
      xgTransposer%mpiAlgo = algo
    else
      xgTransposer%mpiAlgo = TRANS_ALL2ALL
      MSG_COMMENT("Bad value for transposition MPI_algo. Will use ALLTOALL")
    end if
    
    !print *, "xgTransposer%mpiAlgo", xgTransposer%mpiAlgo
    !stop

    !if ( xgTransposer%mpiAlgo == TRANS_ALL2ALL ) then
    !  MSG_COMMENT("Using mpi_alltoall for transposition")
    !else
    !  MSG_COMMENT("Using mpi_gatherv for transposition")
    !end if

    select case (state)
    case (STATE_LINALG)
      ! We are in the linalg representation.
      ! We need to construct the colsrows parallelization
      call xgBlock_getSize(xgBlock_linalg,nrows,ncols)
      
      !print *, "NROWS 2", nrows
      !print *, "NCLOS 2", ncols
      !stop

      !call xmpi_sum(nrows,commLinalg,ierr) !suma  = Npw * broj procesora po redovima

      !print *, "NROWS 3", nrows
      !print *, "NCLOS 3", ncols

      if ( MOD(ncols,ncpuCols) /=0 ) then
        if ( ncols > ncpuCols ) then
          write(message,'(a,i6,a,i6,a)') "Unbalanced parallelization : ", ncols, " columns for ", ncpuCols, " MPI"
          MSG_ERROR(message)
        else
          write(message,'(i6,a)') (ncpuCols-ncols)*ncpuRows, " MPI will not be used  because of the number of columns!!"
          MSG_ERROR(message)
          !ncpuCols = ncols
        end if
      end if

      select case ( space(xgBlock_linalg) )
        case (SPACE_CR, SPACE_R)
          xgTransposer%perPair = 2
        case (SPACE_C)
          xgTransposer%perPair = 1
        case default
          MSG_ERROR("Space value unknown !")
      end select
      
      !print *, "xgTransposer%perPair", xgTransposer%perPair
      !stop

      !write(*,*) "There is a total of ", nrows/xgTransposer%perPair, "rows of real pairs for ", ncpuRows, "fft cpu"

      ! Build the lookup table
      !print *, "NCOLS", ncols
      !print *, "NCPUCOLS", ncpuCols
      !stop
      ABI_MALLOC(xgTransposer%lookup,(1:ncols))
      !ABI_MALLOC(xgTransposer%me_g0_lookup,())
      do icol = 0, ncols-1
        xgTransposer%lookup(icol+1) = MOD(icol,ncpuCols)
      end do
      
      !print *, "xgTransposer%lookup", xgTransposer%lookup
      !stop

      call xgTransposer_makeComm(xgTransposer,ncpuRows,ncpuCols)
      !print *, "PROSAO"
      call xgTransposer_computeDistribution(xgTransposer)
      !print *, "xgTransposer_makeXgBlock"
      call xgTransposer_makeXgBlock(xgTransposer) !ovo valjda pravi send-receive buffer

    case (STATE_COLSROWS)
      MSG_BUG("Not yet implemented")
    case default
      MSG_ERROR("State is undefined")
    end select

    call timab(tim_init,2,tsec)

  end subroutine xgTransposer_init
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
      MSG_ERROR("xgTransposer failed to creat cartesian grid")
    end if

    selectDim = (/ .true., .false. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_ROWS)%comm,ierr)
    if ( ierr /= xmpi_success ) then
      MSG_ERROR("xgTransposer failed to creat rows communicator")
    end if
    selectDim = (/ .false., .true. /)
    call mpi_cart_sub(commColsRows, selectDim, xgTransposer%mpiData(MPI_COLS)%comm,ierr)
    if ( ierr /= xmpi_success ) then
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
    
    print *, "xgTransposer%mpiData(MPI_ROWS)%rank", xgTransposer%mpiData(MPI_ROWS)%rank
    print *, "xgTransposer%mpiData(MPI_ROWS)%size", xgTransposer%mpiData(MPI_ROWS)%size 
    !stop

    xgTransposer%mpiData(MPI_COLS)%rank = xmpi_comm_rank(xgTransposer%mpiData(MPI_COLS)%comm)
    xgTransposer%mpiData(MPI_COLS)%size = xmpi_comm_size(xgTransposer%mpiData(MPI_COLS)%comm)

    print *, "xgTransposer%mpiData(MPI_COLS)%rank", xgTransposer%mpiData(MPI_COLS)%rank
    print *, "xgTransposer%mpiData(MPI_COLS)%size", xgTransposer%mpiData(MPI_COLS)%size 
    !stop

    call xgBlock_setComm(xgTransposer%xgBlock_colsrows,xgTransposer%mpiData(MPI_ROWS)%comm)

  end subroutine xgTransposer_makeComm
!!***

  subroutine xgTransposer_computeDistribution(xgTransposer)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    integer :: nRealPairs
    integer :: ierr
    integer :: icpu
    integer :: ncpuCols

    ABI_MALLOC(xgTransposer%nrowsLinalg,(xgTransposer%mpiData(MPI_LINALG)%size))
    !print *, "xgTransposer%mpiData(MPI_LINALG)%size", xgTransposer%mpiData(MPI_LINALG)%size
    !stop
    nRealPairs = rows(xgTransposer%xgBlock_linalg)/xgTransposer%perPair !number of pair of reals
     
    !print *, "rows(xgTransposer%xgBlock_linalg)", rows(xgTransposer%xgBlock_linalg)
    !print *, "nRealPairs", nRealPairs
    !stop
    
    call xmpi_allgather(nRealPairs,xgTransposer%nrowsLinalg,xgTransposer%mpiData(MPI_LINALG)%comm,ierr)
    if ( ierr /= xmpi_success ) then
      MSG_ERROR("Error while gathering number of rows in linalg")
    end if
    
    !print *, "xgTransposer%nrowsLinalg", xgTransposer%nrowsLinalg
    !stop

    ncpuCols = xgTransposer%mpiData(MPI_COLS)%size
    
    !print *, "ncpuCols", ncpuCols
    !stop
    
    icpu = xgTransposer%mpiData(MPI_ROWS)%rank*ncpuCols  
    
!    print *, "xgTransposer%mpiData(MPI_ROWS)%rank", xgTransposer%mpiData(MPI_ROWS)%rank
!    print *, "xgTransposer%mpiData(MPI_COLS)%rank", xgTransposer%mpiData(MPI_COLS)%rank
!    print *, "ncpuCols", ncpuCols
!    print *, "icpu", icpu
   
    
    !print *, "xgTransposer%nrowsColsRows", xgTransposer%nrowsColsRows
    xgTransposer%nrowsColsRows = sum(xgTransposer%nrowsLinalg(icpu+1:icpu+ncpuCols))
    !print *, "xgTransposer%nrowsColsRows", xgTransposer%nrowsColsRows
    !stop
    
!    if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then 
!      print *, "icpu, ncpuCols", icpu, ncpuCols
!      print *, "xgTransposer%nrowsLinalg(icpu+1)", xgTransposer%nrowsLinalg(icpu+1)
!      print *, "xgTransposer%nrowsLinalg(icpu+ncpuCols)", xgTransposer%nrowsLinalg(icpu+ncpuCols)
!    end if
    
    !print *, "nrowsColsRows", " icpu+1 ", "icpu+ncpuCols",  xgTransposer%nrowsColsRows, icpu+1, icpu+ncpuCols
    !stop
    
    xgTransposer%ncolsColsRows = cols(xgTransposer%xgBlock_linalg)/ncpuCols
    
    !print *, "xgTransposer%ncolsColsRows", xgTransposer%ncolsColsRows
    !stop

    !write(*,*) "In linalg, # of real pairs:", xgTransposer%nrowsLinalg
    !write(*,*) "In rows, # of real pairs for proc ", xgTransposer%mpiData(MPI_ROWS)%rank, ":", xgTransposer%nrowsColsRows
    !write(*,*) "In cols, # of cols for proc ", xgTransposer%mpiData(MPI_COLS)%rank, ":", xgTransposer%ncolsColsRows

  end subroutine xgTransposer_computeDistribution

  subroutine xgTransposer_makeXgBlock(xgTransposer)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    !integer :: cols, rows

    select case (xgTransposer%state)
    case (STATE_LINALG)
      ! Assume xgBlock_colsrows is empty and not constructed because user cannot
      ! predict the size
      if ( allocated(xgTransposer%buffer) ) then
        ABI_FREE(xgTransposer%buffer)
      end if
      if ( xgTransposer%mpiData(MPI_COLS)%size == 1 ) then
        xgTransposer%xgBlock_colsrows = xgTransposer%xgBlock_linalg
        
        !print *, "xgTransposer%ncolsColsRows", xgTransposer%ncolsColsRows
        !print *, "xgTransposer%nrowsColsRows", xgTransposer%nrowsColsRows
        !stop
        
      else
        !print *, "HEREEEEEEE"
        !print *, "xgTransposer%ncolsColsRows", xgTransposer%ncolsColsRows
        !print *, "xgTransposer%nrowsColsRows", xgTransposer%nrowsColsRows
        !stop
 
        !stop
        ABI_MALLOC(xgTransposer%buffer,(2,xgTransposer%ncolsColsRows*xgTransposer%nrowsColsRows))
        !print *, "xgTransposer%ncolsColsRows*xgTransposer%nrowsColsRows", xgTransposer%ncolsColsRows*xgTransposer%nrowsColsRows
        
        !print *, "xgTransposer%perPair*xgTransposer%nrowsColsRows", xgTransposer%perPair*xgTransposer%nrowsColsRows
        !print *, "xgTransposer%ncolsColsRows", xgTransposer%ncolsColsRows
        !stop
        !OVO JE VALJDA SEND RECEIVE BUFFER
!        if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then
!          print *, "xcgColsRows PRE MAPIRANJA:"
!          call xgBlock_print(xgTransposer%xgBlock_colsrows,6)  
!        end if
        call xgBlock_map(xgTransposer%xgBlock_colsrows,xgTransposer%buffer,space(xgTransposer%xgBlock_linalg),&
          xgTransposer%perPair*xgTransposer%nrowsColsRows,&
          xgTransposer%ncolsColsRows,xgTransposer%mpiData(MPI_ROWS)%comm)
!        if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then
!          print *, "xcgColsRows NAKON MAPIRANJA:"
!          call xgBlock_print(xgTransposer%xgBlock_colsrows,6)  
!        end if
      end if
    case (STATE_COLSROWS)
      MSG_ERROR("Not yet implemented")
    case default
      MSG_ERROR("State unknown")
    end select
  end subroutine xgTransposer_makeXgBlock

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
      MSG_ERROR("Bad value for toState")
    end if

    !write(std_out,*) "linalg", rows(xgTransposer%xgBlock_linalg)*cols(xgTransposer%xgBlock_linalg)
    !write(std_out,*) "colsrows", rows(xgTransposer%xgBlock_colsrows)*cols(xgTransposer%xgBlock_colsrows)
    
    !print *, "TRANSPOSEEEEEEEEEEE"
    select case (toState)
    case (STATE_LINALG)
      if ( xgTransposer%state == STATE_LINALG ) then
        MSG_WARNING("Array linalg has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_COLS)%size > 1 ) then
        call xgTransposer_toLinalg(xgTransposer)
      else
        !print *, "SAMO PROMENIO STANJE"
        xgTransposer%state = STATE_LINALG
      end if
    case (STATE_COLSROWS)
      !print *, "OVDEEEEEEEEEEEEEEEE"
      if ( xgTransposer%state == STATE_COLSROWS ) then
        MSG_WARNING("Array colsrows has already been transposed")
      end if
      if ( xgTransposer%mpiData(MPI_COLS)%size > 1 ) then
        !print *, "LUDIRANJE"
        !stop
        !print *, "1111111111111"
        call xgTransposer_toColsRows(xgTransposer)
        !print *, "2222222222222"
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

    call timab(tim_toLinalg,1,tsec)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu  !!TODO TEST TEST TEST

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
     !print *, "TRANS ALLTOALLV"  !!TODO OVDE SE BAFERI NESTO IZDEBILISU SKONTATI ZASTO :(((
     call timab(tim_all2allv,1,tsec)
     call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
                         recvbuf, recvcounts, rdispls, &
                         comm, ierr)
     call timab(tim_all2allv,2,tsec)
     !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
     !                    recvbuf, recvcounts, rdispls, &
     !                    comm, request(myrequest))
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
     MSG_BUG("This algo does not exist")
   end select

   xgTransposer%state = STATE_LINALG


   !ABI_MALLOC(status,(MPI_STATUS_SIZE))
   !call mpi_wait(request(myrequest),status,ierr)
   if ( ierr /= xmpi_success ) then
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

   !do icpu = 1, size(request)
   !  if ( icpu /= myrequest ) then
   !    call mpi_wait(request(icpu),status,ierr)
   !    if ( ierr /= MPI_SUCCESS ) then
   !      MSG_ERROR("Error while waiting for other mpi")
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

  subroutine xgTransposer_toColsRows(xgTransposer) !!TODO DEBUGGGGGGGGGGGGGGGGGGGGGGGGGGGGg

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

    call timab(tim_toColsRows,1,tsec)

   ncpu = xgTransposer%mpiData(MPI_COLS)%size
   comm = xgTransposer%mpiData(MPI_COLS)%comm
   me = xgTransposer%mpiData(MPI_ROWS)%rank*ncpu
   
   !print *, "NCPU", ncpu
   !print *, "COMM", comm
   !print *, "ME", me
   !stop

   nrowsColsRows = xgTransposer%nrowsColsRows
   ncolsColsRows = xgTransposer%ncolsColsRows
   
   !print *, "nrowsColsRows", nrowsColsRows
   !print *, "ncolsColsRows", ncolsColsRows
   !stop

   !print *, "SIZE", xgTransposer%mpiData(MPI_LINALG)%size
   !stop

   nrowsLinalg => xgTransposer%nrowsLinalg  !size = 4?
   nrowsLinalgMe = nrowsLinalg(xgTransposer%mpiData(MPI_LINALG)%rank+1)
   
   !print *, "nrowsLinalgMe", nrowsLinalgMe
   !stop
   

   ABI_MALLOC(recvbuf,(2,nrowsColsRows*ncolsColsRows)) !2x8000x1000
   ABI_MALLOC(recvcounts,(ncpu))
   ABI_MALLOC(rdispls,(ncpu))

   recvcounts(:) = 2*nrowsLinalg(me+1:me+ncpu)*ncolsColsRows
   !print *, "me+1, me+ncpu, nrowsLinalg(me+1:me+ncpu), ncolsColsRows, recvcounts(:)", me+1, me+ncpu, nrowsLinalg(me+1:me+ncpu), ncolsColsRows, recvcounts(:) 
   !stop
   rdispls(1) = 0
   do icpu = 2, ncpu
     rdispls(icpu) = rdispls(icpu-1)+recvcounts(icpu-1)
   end do
   !print *, "       "
!   if (xgTransposer%mpiData(MPI_LINALG)%rank == xgTransposer%debug_rank) then
!     print *, "RANK", xgTransposer%mpiData(MPI_LINALG)%rank
!     print *, "nrowsColsRows", nrowsColsRows 
!     print *, "ncolsColsRows", ncolsColsRows
!     write(*,*) "RECVCOUNTS", recvcounts
!     write(*,*) "RDISPL", rdispls
!     !stop
!   end if

   select case(xgTransposer%mpiAlgo)
   case (TRANS_ALL2ALL)
     !print *, "OVDE"
     !stop
     ABI_MALLOC(sendcounts,(ncpu))
     ABI_MALLOC(sdispls,(ncpu))
     !ABI_MALLOC(request,(1))
     !myrequest = 1

     sendcounts(:) = 2*nrowsLinalgMe*ncolsColsRows !! Thank you fortran for not starting at 0 !
     !write(*,*) "SENDCOUNTS", sendcounts
     !print *, "KRAJ"
     !stop
     sdispls(1) = 0
     do icpu = 2, ncpu
       sdispls(icpu) = sdispls(icpu-1)+sendcounts(icpu-1)
     end do
     
!     if (xgTransposer%mpiData(MPI_LINALG)%rank == xgTransposer%debug_rank) then
!       print *, "RANK", xgTransposer%mpiData(MPI_LINALG)%rank
!       write(*,*) "SENDCOUNTS", sendcounts
!       write(*,*) "SDISPL", sdispls
!     end if
     
     !print *, "xgTransposer%perPair", xgTransposer%perPair
     !print *, "cols(xgTransposer%xgBlock_linalg)*nrowsLinalgMe", cols(xgTransposer%xgBlock_linalg)*nrowsLinalgMe
     !stop

!     if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then
!       print *, "SEND BUFF BEFORE ALLTOALL"
!       call xgBlock_print(xgTransposer%xgBlock_linalg,6)
!       !stop
!     end if
     !print *, "BEFORE ALLTOALL"
     call xgBlock_reverseMap(xgTransposer%xgBlock_linalg,sendbuf, &
&      xgTransposer%perPair,cols(xgTransposer%xgBlock_linalg)*nrowsLinalgMe)
     !write(*,*) "Before ialltoall"
     call timab(tim_all2allv,1,tsec)
     call xmpi_alltoallv(sendbuf, sendcounts, sdispls, &
                         recvbuf, recvcounts, rdispls, &
                         comm, ierr)
     call timab(tim_all2allv,2,tsec)
     !print *, "AFTER ALLTOALL"
     !call xmpi_ialltoallv(sendbuf, sendcounts, sdispls, &
     !                    recvbuf, recvcounts, rdispls, &
     !                    comm, request(myrequest))
     !write(*,*) "After ialltoall"

   case (TRANS_GATHER)
     !print *, "USAO TRANDZA"
     !stop
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
     MSG_BUG("This algo does not exist")
   end select

   !ABI_MALLOC(status,(MPI_STATUS_SIZE))
   !call mpi_wait(request(myrequest),status,ierr)
   !write(*,*) "Request ended"
   if ( ierr /= xmpi_success ) then
     MSG_ERROR("Error while waiting for mpi")
   end if
   !write(*,*) "with success"

   !print *, "REORGANIZER"
   !stop
   call xgTransposer_reorganizeData(xgTransposer,recvbuf)

   !print *, "AFTER REORGANIZER"
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
   !      MSG_ERROR("Error while waiting for other mpi")
   !    end if
   !  end if
   !end do
   !ABI_FREE(status)
   !ABI_FREE(request)

    call timab(tim_toColsRows,2,tsec)

  end subroutine xgTransposer_toColsRows
!!***

  subroutine xgTransposer_reorganizeData(xgTransposer,bufferMess)

    type(xgTransposer_t), intent(inout) :: xgTransposer
    double precision    , intent(inout) :: bufferMess(:,:)
    double precision, pointer :: bufferOrdered(:,:) => null()
    integer :: shiftCpu
    integer :: nrowsColsRows
    integer :: ncolsColsRows
    integer :: tos,toe,froms,frome
    integer :: col, icpu
    integer :: me
    integer :: nPair
    integer, pointer :: nrowsLinalg(:)
    double precision :: tsec(2)

    call timab(tim_reorganize,1,tsec)

    me = xgTransposer%mpiData(MPI_ROWS)%rank*xgTransposer%mpiData(MPI_COLS)%size
    nrowsColsRows = xgTransposer%nrowsColsRows
    ncolsColsRows = xgTransposer%ncolsColsRows
    nPair = nrowsColsRows*ncolsColsRows
    
    !print *, "ME", me
    !print *, "nrowsColsRows", nrowsColsRows
    !print *, "ncolsColsRows", ncolsColsRows
    !print *, "nPair", nPair
    !print *, "xgTransposer%mpiData(MPI_LINALG)%rank", xgTransposer%mpiData(MPI_LINALG)%rank
    !stop
!    if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then
!      print *, "xcgColsRows PRE REVERSE MAPIRANJA:"
!      call xgBlock_print(xgTransposer%xgBlock_colsrows,6)  
!    end if  
    call xgBlock_reverseMap(xgTransposer%xgBlock_colsrows,bufferOrdered,xgTransposer%perPair,nPair)
!    if (xmpi_comm_rank(xmpi_world) == xgTransposer%debug_rank) then
!      print *, "xcgColsRows NAKON REVERSE MAPIRANJA:"
!      call xgBlock_print(xgTransposer%xgBlock_colsrows,6)  
!    end if  
    !bufferOrdered je samo prijemni bafer nula

    nrowsLinalg => xgTransposer%nrowsLinalg

    select case (xgTransposer%state)
    case (STATE_LINALG)
      !print *, "STATE_LINALG"
      !print *, "ncolsColsRows", ncolsColsRows
      !print *, "xgTransposer%mpiData(MPI_COLS)%size", xgTransposer%mpiData(MPI_COLS)%size
      !stop
      ! We are going to STATE_COLSROWS so we are after all2all
      !$omp parallel do private(shiftCpu,toe,tos,frome,froms), collapse(2)
!      if (xgTransposer%mpiData(MPI_LINALG)%rank == 1) then 
!        print *, "ncolsColsRows", ncolsColsRows
!      end if
      !TODO TODO TODO TODO vidi ovde na malom primeru
      do col = 1, ncolsColsRows
        do icpu = 1, xgTransposer%mpiData(MPI_COLS)%size
          shiftCpu = ncolsColsRows*sum(nrowsLinalg(me+1:me+icpu-1))
          tos=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu-1))+1)
          toe=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu)))
          froms=(shiftCpu+(col-1)*nrowsLinalg(me+icpu)+1)
          frome=(shiftCpu+col*nrowsLinalg(me+icpu))
!          if (xgTransposer%mpiData(MPI_LINALG)%rank == xgTransposer%debug_rank) then 
!            if (col < 5 .and. icpu == 1) then
!              print *, "col", col
!              print *, "shiftCpu", shiftCpu
!              print *, "tos", tos
!              print *, "toe", toe
!              print *, "froms", froms
!              print *, "frome", frome
!            end if
!          end if
          !print *, "loc BORDERED", loc(bufferOrdered)
          bufferOrdered(:,tos:toe) = bufferMess(:,froms:frome)
        end do
      end do
!      if (xgTransposer%mpiData(MPI_LINALG)%rank == xgTransposer%debug_rank) then 
!        call xgBlock_print(xgTransposer%xgBlock_colsrows,6)     
!        print *, "bufferOrdered", bufferOrdered(1,:)
!        print *, "bufferMess", bufferMess(1,:)
!      end if
    case (STATE_COLSROWS)
      !print *, "xgTransposer_reorganizeData"
      ! We are going to STATE_LINALG so we are before all2all
      !$omp parallel do private(shiftCpu,toe,tos,frome,froms), collapse(2)
      do col = 1, ncolsColsRows
        do icpu = 1, xgTransposer%mpiData(MPI_COLS)%size
          shiftCpu = ncolsColsRows*sum(nrowsLinalg(me+1:me+icpu-1))
          tos=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu-1))+1)
          toe=((col-1)*nrowsColsRows+sum(nrowsLinalg(me+1:me+icpu)))
          froms=(shiftCpu+(col-1)*nrowsLinalg(me+icpu)+1)
          frome=(shiftCpu+col*nrowsLinalg(me+icpu))
          !print *, "loc BMESS", loc(bufferOrdered)
          bufferMess(:,froms:frome) = bufferOrdered(:,tos:toe)
        end do
      end do
    end select

    call timab(tim_reorganize,2,tsec)

  end subroutine xgTransposer_reorganizeData


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
    call mpi_comm_free(xgTransposer%mpiData(MPI_ROWS)%comm,i)
    call mpi_comm_free(xgTransposer%mpiData(MPI_COLS)%comm,i)
    call mpi_comm_free(xgTransposer%mpiData(MPI_2DCART)%comm,i)
#else
    ABI_UNUSED(i)
#endif

    if ( allocated(xgTransposer%lookup) ) then
      ABI_FREE(xgTransposer%lookup)
    end if

    if ( associated(xgTransposer%nrowsLinalg) ) then
      ABI_FREE(xgTransposer%nrowsLinalg)
    end if

    if ( allocated(xgTransposer%buffer) ) then
      ABI_FREE(xgTransposer%buffer)
    end if
    call timab(tim_free,2,tsec)


  end subroutine xgTransposer_free
!!***

  subroutine xgTransposer_getCPURow(xgTransposer, row)
  
    type(xgTransposer_t), intent(inout) :: xgTransposer
    type(integer), intent(inout) :: row
    
    row = xgTransposer%mpiData(MPI_ROWS)%rank
  
  end subroutine xgTransposer_getCPURow

end module m_xgTransposer
!!***

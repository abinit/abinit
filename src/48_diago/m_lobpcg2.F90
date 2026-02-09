!!****f* ABINIT/m_lobpcg2
!! NAME
!! m_lobpcg2
!!
!! FUNCTION
!! This module contains the types and routines used to apply the
!! LOBPCG method (second version introduced by J. Bieder), using the xg_tools.
!!
!! COPYRIGHT
!! Copyright (C) 2015-2026 ABINIT group (J. Bieder, L. Baguet)
!! This file is distributed under the terms of the
!! gnu general public license, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! for the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

! nvtx related macro definition
#include "nvtx_macros.h"

module m_lobpcg2

  use m_xg
  use m_xgTransposer
  use m_xg_ortho_RR
  use defs_basis
  use m_abicore
  use m_errors
  use m_xomp
#ifdef HAVE_OPENMP
  use omp_lib
#endif
  use m_xmpi
 use, intrinsic :: iso_c_binding, only: c_size_t

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

  use m_time, only : timab

  implicit none

  private

  integer, parameter :: VAR_X   = 1000
  integer, parameter :: VAR_W   = 1001
  integer, parameter :: VAR_P   = 1002
  integer, parameter :: VAR_XW  = 1010
  integer, parameter :: VAR_WP  = 1011
  integer, parameter :: VAR_XWP = 1100

  integer, parameter :: tim_init     = 1651
  integer, parameter :: tim_free     = 1652
  integer, parameter :: tim_copy     = 1653
  integer, parameter :: tim_getAX_BX = 1654
  integer, parameter :: tim_ortho    = 1655
  integer, parameter :: tim_nbdbuf   = 1656
!  integer, parameter :: tim_RR       = 1657
  integer, parameter :: tim_maxres   = 1658
  integer, parameter :: tim_ax_bx    = 1659
  integer, parameter :: tim_pcond    = 1660
!  integer, parameter :: tim_hegv     = 1661

  integer, parameter :: tim_Bortho_X    = 1641
  integer, parameter :: tim_Bortho_XW   = 1642
  integer, parameter :: tim_Bortho_XWP  = 1643
  integer, parameter :: tim_Bortho_Xall = 1644
  integer, parameter :: tim_RR_X        = 1645
  integer, parameter :: tim_RR_XW       = 1646
  integer, parameter :: tim_RR_XWP      = 1647
  integer, parameter :: tim_RR_Xall     = 1648

  integer, parameter :: tim_transpose   = 1649

  type, public :: lobpcg_t
    logical :: is_nested                     ! For OpenMP nested region
    integer :: spacedim                      ! Space dimension for one vector
    integer :: neigenpairs                   ! Number of eigen values/vectors we want
    integer :: nblock                        ! Number of block in the space of dim spacedim
    integer :: blockdim                      ! Number of vectors in one block
    integer :: nline                         ! Number of line to perform
    integer :: spacecom                      ! Communicator for MPI
    integer :: paral_kgb                     ! paral_kgb formalism or not
    integer :: comm_rows                     ! communicator for rows
    integer :: comm_cols                     ! communicator for cols
    integer :: me_g0                         ! =1 if the processor have G=0 (linalg representation)
    integer :: me_g0_fft                     ! =1 if the processor have G=0 (fft_representation)
    integer :: gpu_option                    ! Which GPU version is used (0=none)
    integer :: gpu_thread_limit              ! When GPU is enabled, how many CPU threads to use in sensitive areas
    double precision :: tolerance            ! Tolerance on the residu to stop the minimization
    integer :: prtvol
    type(xgBlock_t) :: AllX0 ! Block of initial and final solution.
    type(xgBlock_t) :: X0 ! Block of initial and final solution.
                                             ! Dim is (cplx*spacedim,neigenpair)
    !double precision, allocatable :: eig(:)  ! Actual eigen values

    ! Some denomination here
    ! X, W, P are the shifts where their block vectors starts.
    ! X is the current block we are solving
    ! W is the residu of the pencil A-lambda*B ie W=A*X-lambdaB*X
    ! P correspond to X_{iline} - X_{iline-1}
    ! if prefixed by A (AX, AW, AP) it means the result of A*X (A*W, A*P ! respectively)
    ! Idem for B (BX, BW, BP) which is the result of B*X (B*W, B*P)
    ! An extend to this is AXWP which is the vector [AX,AW,AP]
    ! and BXWP which is [BX,BW,BP] columnwise
    ! Now we only allocate the triple XWP matrix and each individual X, W and P
    ! Example : XWP=[x1,x2,x3,w1,w2,w3,p1,p2,p3] X points to the indice of x1-1, W points to w1-1 and
    ! P to p1-1
    type(xg_t) :: XWP
    type(xg_t) :: AXWP
    type(xg_t) :: BXWP

    type(xgBlock_t) :: XColsRows
    type(xgBlock_t) :: AXColsRows
    type(xgBlock_t) :: BXColsRows

    type(xgTransposer_t) :: xgTransposerX
    type(xgTransposer_t) :: xgTransposerAX
    type(xgTransposer_t) :: xgTransposerBX

    type(xgBlock_t) :: WColsRows
    type(xgBlock_t) :: AWColsRows
    type(xgBlock_t) :: BWColsRows

    type(xgTransposer_t) :: xgTransposerW
    type(xgTransposer_t) :: xgTransposerAW
    type(xgTransposer_t) :: xgTransposerBW

    type(xgBlock_t) :: X ! Shift to apply to start reading the X values
    type(xgBlock_t) :: W ! Shift to apply to start reading the W values
    type(xgBlock_t) :: P ! Shift to apply to start reading the P values
    type(xgBlock_t) :: XW ! Shift to apply to start reading the P values
    type(xgBlock_t) :: WP ! Shift to apply to start reading the P values

    type(xgBlock_t) :: AX ! Shift to apply to start reading the X values
    type(xgBlock_t) :: AW ! Shift to apply to start reading the W values
    type(xgBlock_t) :: AP ! Shift to apply to start reading the P values
    type(xgBlock_t) :: AXW ! Shift to apply to start reading the P values
    type(xgBlock_t) :: AWP ! Shift to apply to start reading the P values

    type(xgBlock_t) :: BX ! Shift to apply to start reading the X values
    type(xgBlock_t) :: BW ! Shift to apply to start reading the W values
    type(xgBlock_t) :: BP ! Shift to apply to start reading the P values
    type(xgBlock_t) :: BXW ! Shift to apply to start reading the P values
    type(xgBlock_t) :: BWP ! Shift to apply to start reading the P values

    type(xg_t) :: AllBX0
    type(xg_t) :: AllAX0
    type(xgBlock_t) :: BX0

    type(xgBlock_t) :: AllX0ColsRows
    type(xgBlock_t) :: AllAX0ColsRows
    type(xgBlock_t) :: AllBX0ColsRows

    type(xgTransposer_t) :: xgTransposerAllX0
    type(xgTransposer_t) :: xgTransposerAllAX0
    type(xgTransposer_t) :: xgTransposerAllBX0

    ! Variable for work for lapack
  end type lobpcg_t

  public :: lobpcg_init
  public :: lobpcg_memInfo
  public :: lobpcg_run
  public :: lobpcg_free

  contains


  subroutine lobpcg_init(lobpcg, neigenpairs, spacedim, blockdim, tolerance, nline, &
&      space, spacecom, paral_kgb, comm_rows, comm_cols, me_g0, me_g0_fft, gpu_option, &
&      gpu_thread_limit)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: neigenpairs
    integer         , intent(in   ) :: spacedim
    integer         , intent(in   ) :: blockdim
    double precision, intent(in   ) :: tolerance
    integer         , intent(in   ) :: nline
    integer         , intent(in   ) :: comm_rows,comm_cols
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: spacecom
    integer         , intent(in   ) :: paral_kgb
    integer         , intent(in   ) :: me_g0
    integer         , intent(in   ) :: me_g0_fft
    integer         , intent(in   ) :: gpu_option
    integer,optional, intent(in   ) :: gpu_thread_limit
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)
    lobpcg%neigenpairs = neigenpairs
    lobpcg%spacedim    = spacedim
    lobpcg%blockdim    = blockdim
    if (tolerance > 0.0) then
      lobpcg%tolerance = tolerance
    else
      lobpcg%tolerance = 1.0e-20
    end if
    lobpcg%nline       = nline
    lobpcg%spacecom    = spacecom
    lobpcg%nblock      = neigenpairs / blockdim
    lobpcg%paral_kgb   = paral_kgb
    lobpcg%comm_rows   = comm_rows
    lobpcg%comm_cols   = comm_cols
    lobpcg%me_g0       = me_g0
    lobpcg%me_g0_fft   = me_g0_fft
    lobpcg%gpu_option  = gpu_option
    lobpcg%gpu_thread_limit  = 0
    if(present(gpu_thread_limit)) lobpcg%gpu_thread_limit  = gpu_thread_limit

    call lobpcg_allocateAll(lobpcg,space,me_g0)
    call timab(tim_init,2,tsec)

  end subroutine lobpcg_init


  subroutine lobpcg_allocateAll(lobpcg,space,me_g0)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: me_g0
    integer :: spacedim
    integer :: blockdim

    spacedim = lobpcg%spacedim
    blockdim = lobpcg%blockdim

    call lobpcg_free(lobpcg) ! Make sure everything is not allocated and
    ! pointer point to null()

    call xg_init(lobpcg%XWP,space,spacedim,3*blockdim,lobpcg%spacecom,me_g0=me_g0,gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%XWP,lobpcg%X,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%W,spacedim,blockdim,fcol=blockdim+1)
    call xg_setBlock(lobpcg%XWP,lobpcg%P,spacedim,blockdim,fcol=2*blockdim+1)
    call xg_setBlock(lobpcg%XWP,lobpcg%XW,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%WP,spacedim,2*blockdim,fcol=blockdim+1)

    call xg_init(lobpcg%AXWP,space,spacedim,3*blockdim,lobpcg%spacecom,me_g0=me_g0,gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AX,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AW,spacedim,blockdim,fcol=blockdim+1)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AP,spacedim,blockdim,fcol=2*blockdim+1)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AXW,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AWP,spacedim,2*blockdim,fcol=blockdim+1)

    call xg_init(lobpcg%BXWP,space,spacedim,3*blockdim,lobpcg%spacecom,me_g0=me_g0,gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BX,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BW,spacedim,blockdim,fcol=blockdim+1)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BP,spacedim,blockdim,fcol=2*blockdim+1)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BXW,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BWP,spacedim,2*blockdim,fcol=blockdim+1)

    if ( lobpcg%nblock /= 1 ) then
      call xg_init(lobpcg%AllBX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom,me_g0=me_g0,gpu_option=lobpcg%gpu_option)
      call xg_init(lobpcg%AllAX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom,me_g0=me_g0,gpu_option=lobpcg%gpu_option)
    else
      lobpcg%AllBX0%self = lobpcg%BX
      lobpcg%AllAX0%self = lobpcg%AX
    end if


  end subroutine lobpcg_allocateAll


  function lobpcg_memInfo(neigenpairs, spacedim, space, paral_kgb, blockdim) result(arraymem)

    integer         , intent(in   ) :: neigenpairs
    integer         , intent(in   ) :: spacedim
    integer         , intent(in   ) :: blockdim
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: paral_kgb
    integer(kind=c_size_t) :: memXWP
    integer(kind=c_size_t) :: memAXWP
    integer(kind=c_size_t) :: memBXWP
    integer(kind=c_size_t) :: mem_xgTransposer
    integer(kind=c_size_t) :: memAllBX0
    integer(kind=c_size_t) :: memAllAX0
    integer(kind=c_size_t) :: memeigenvalues3N
    integer(kind=c_size_t) :: membufferOrtho
    integer(kind=c_size_t) :: membufferBOrtho
    integer(kind=c_size_t) :: memsubA
    integer(kind=c_size_t) :: memsubB
    integer(kind=c_size_t) :: memsubBtmp
    integer(kind=c_size_t) :: memvec
    integer(kind=c_size_t) :: memRR
    integer(kind=c_size_t) :: maxmemTmp
    integer(kind=c_size_t) :: cplx
    integer :: nblock
    integer(kind=c_size_t) :: arraymem(2)

    cplx = 1 ; if ( space == SPACE_C ) cplx = 2
    nblock = neigenpairs/blockdim

    ! Permanent in lobpcg
    memXWP  = int(cplx,c_size_t)* kind(1.d0) * spacedim * 3*blockdim
    memAXWP = int(cplx,c_size_t)* kind(1.d0) * spacedim * 3*blockdim
    memBXWP = int(cplx,c_size_t)* kind(1.d0) * spacedim * 3*blockdim
    if(paral_kgb > 1) then
      ! xgtransposer *_ColRows buffers for X, AX, BX, W, AW, BW + internal send/recvbuf
      mem_xgTransposer = int(cplx,c_size_t)* kind(1.d0) * spacedim * 7*blockdim
    else
      mem_xgTransposer = 0
    end if
    if ( nblock > 1 ) then
      memAllAX0 = int(cplx,c_size_t) * kind(1.d0) * spacedim * 3*blockdim
      memAllBX0 = int(cplx,c_size_t) * kind(1.d0) * spacedim * 3*blockdim
      membufferOrtho = int(cplx,c_size_t) * kind(1.d0) * blockdim * (nblock-1) * blockdim
      mem_xgTransposer = mem_xgTransposer + int(cplx,c_size_t) * kind(1.d0) * spacedim * 3*blockdim
    else
      memAllAX0 = 0
      memAllBX0 = 0
      membufferOrtho = 0
    endif
    memeigenvalues3N = kind(1.d0) * 3*blockdim

    ! Temporary arrays

    ! For the moment being, only Bortho with X or WP at the same time
    membufferBOrtho = int(cplx,c_size_t) * kind(1.d0) * 2*blockdim * 2*blockdim
    memsubA = int(cplx,c_size_t) * kind(1.d0) * 3*blockdim * 3*blockdim
    memsubB = int(cplx,c_size_t) * kind(1.d0) * 3*blockdim * 3*blockdim
    memsubBtmp = int(cplx,c_size_t) * kind(1.d0) * spacedim * blockdim
    memvec = int(cplx,c_size_t) * kind(1.d0) * 3*blockdim * blockdim
    memRR = max(memsubA+memsubB+memvec,memsubBtmp+memvec)

    maxmemTmp = max( membufferBOrtho,memRR,membufferOrtho )

    arraymem(1) = memXWP+memAXWP+memBXWP+memAllAX0+memAllBX0+memeigenvalues3N+mem_xgTransposer
    arraymem(2) = maxmemTmp

  end function lobpcg_memInfo


  subroutine lobpcg_run(lobpcg, X0, getAX_BX, pcond, eigen, occ, residu, prtvol, nspinor, isppol, ikpt, inonsc, istep, nbdbuf)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: X0   ! Full initial vectors
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: occ
    type(xgBlock_t), intent(inout) :: residu
    type(xgBlock_t), intent(in)    :: pcond
    integer        , intent(in   ) :: prtvol
    integer        , intent(in   ) :: nspinor
    integer        , intent(in   ) :: isppol,ikpt,inonsc,istep,nbdbuf

    type(xg_t) :: eigenvalues3N   ! eigen values for Rayleight-Ritz
    type(xg_t) :: residu_eff
    type(xgBlock_t) :: eigenvaluesN   ! eigen values for Rayleight-Ritz
    type(xgBlock_t) :: eigenvalues2N   ! eigen values for Rayleight-Ritz
    logical :: skip,compute_residu
    integer :: blockdim, blockdim3, blockdim2
    integer :: spacedim
    integer :: iblock, nblock
    integer :: iline, nline
    integer :: rows_tmp, cols_tmp, nband_eff, iband_min, iband_max
    type(xgBlock_t) :: eigenBlock   !
    type(xgBlock_t) :: residuBlock,occBlock
    double precision :: maxResidu, minResidu, dummy
    double precision :: dlamch,tolerance
    integer :: ierr = 0
    integer :: nrestart
    double precision :: tsec(2)
    character(len=500) :: msg

    interface
      subroutine getAX_BX(X,AX,BX)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
      end subroutine getAX_BX
    end interface

!    call timab(tim_run,1,tsec)

    lobpcg%prtvol = prtvol

    tolerance=2*dlamch('E')

    blockdim = lobpcg%blockdim
    blockdim2 = 2*blockdim
    blockdim3 = 3*blockdim
    spacedim = lobpcg%spacedim

    nblock = lobpcg%nblock
    nline = lobpcg%nline

    if (nbdbuf>0) then
       nband_eff = lobpcg%neigenpairs - nbdbuf
    else
       nband_eff = lobpcg%neigenpairs
    end if

    call xgBlock_getSize(eigen,rows_tmp, cols_tmp)
    if ( rows_tmp /= lobpcg%neigenpairs .and. cols_tmp /= 1 ) then
      ABI_ERROR("Error eigen size")
    endif
    call xgBlock_getSize(X0,rows_tmp, cols_tmp)
    if ( rows_tmp /= lobpcg%spacedim ) then
      ABI_ERROR("Error X0 spacedim")
    endif
    if ( cols_tmp /= lobpcg%neigenpairs ) then
      ABI_ERROR("Error X0 npairs")
    endif

    if (isppol==1.and.ikpt==1.and.inonsc==1.and.istep==1) then
      write(msg,'(a,es16.6)') ' lobpcg%tolerance(tolwfr_diago)=',lobpcg%tolerance
      call wrtout(std_out,msg,'COLL')
    end if

    call xg_init(eigenvalues3N,SPACE_R,blockdim3,1, gpu_option=lobpcg%gpu_option)
    call xg_setBlock(eigenvalues3N,eigenvaluesN,blockdim,1)
    call xg_setBlock(eigenvalues3N,eigenvalues2N,blockdim2,1)

    call xgBlock_reshape(eigen,blockdim,nblock)
    call xgBlock_reshape(residu,blockdim,nblock)
    call xgBlock_reshape(occ,blockdim,nblock)

    lobpcg%AllX0 = X0

    call xg_init(residu_eff,SPACE_R,blockdim,1,gpu_option=ABI_GPU_DISABLED)

    if ( lobpcg%paral_kgb == 1 ) then
      call timab(tim_transpose,1,tsec)
      call xgTransposer_constructor(lobpcg%xgTransposerX,lobpcg%X,lobpcg%XColsRows,nspinor,&
        STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0,lobpcg%me_g0_fft,&
        gpu_option=lobpcg%gpu_option,gpu_thread_limit=lobpcg%gpu_thread_limit)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerAX,lobpcg%xgTransposerX,&
        lobpcg%AX,lobpcg%AXColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerBX,lobpcg%xgTransposerX,&
        lobpcg%BX,lobpcg%BXColsRows,STATE_LINALG)

      call xgTransposer_copyConstructor(lobpcg%xgTransposerW,lobpcg%xgTransposerX,&
        lobpcg%W,lobpcg%WColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerAW,lobpcg%xgTransposerX,&
        lobpcg%AW,lobpcg%AWColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerBW,lobpcg%xgTransposerX,&
        lobpcg%BW,lobpcg%BWColsRows,STATE_LINALG)
      call timab(tim_transpose,2,tsec)
    else
      call xgBlock_setBlock(lobpcg%X, lobpcg%XColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AX, lobpcg%AXColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%BX, lobpcg%BXColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%W, lobpcg%WColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AW, lobpcg%AWColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%BW, lobpcg%BWColsRows, spacedim, blockdim)
    end if

    !! Start big loop over blocks
    do iblock = 1, nblock
      ABI_NVTX_START_RANGE(NVTX_LOBPCG2_BLOCK)
      nrestart = 0

      call lobpcg_getX0(lobpcg,iblock)
      call xgBlock_setBlock(residu,residuBlock,blockdim,1,fcol=iblock)
      call xgBlock_setBlock(occ,   occBlock,   blockdim,1,fcol=iblock)

      if ( iblock > 1 ) then
        call lobpcg_setPreviousX0_BX0(lobpcg,iblock)

        ! Orthogonalize current iblock X block With Respect To previous Blocks in B-basis
        call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%X,iblock)
      end if

      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerX,STATE_COLSROWS)
        lobpcg%xgTransposerAX%state=STATE_COLSROWS
        lobpcg%xgTransposerBX%state=STATE_COLSROWS
        call timab(tim_transpose,2,tsec)
      end if
      ! Initialize some quantitites (AX and BX)
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(lobpcg%XColsRows,lobpcg%AXColsRows,lobpcg%BXColsRows)
      call xgBlock_zero_im_g0(lobpcg%AXColsRows)
      call xgBlock_zero_im_g0(lobpcg%BXColsRows)
      call timab(tim_ax_bx,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerX,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAX,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerBX,STATE_LINALG)
        call timab(tim_transpose,2,tsec)
      end if

      ! B-orthonormalize X, BX and AX
      call xg_Borthonormalize(lobpcg%X,lobpcg%BX,ierr,tim_Bortho_X,lobpcg%gpu_option,AX=lobpcg%AX)

      ! Do first RR on X to get the first eigen values
      call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvaluesN,ierr,lobpcg%prtvol,tim_RR_X,lobpcg%gpu_option)

      compute_residu = .true.

      do iline = 1, nline
        ABI_NVTX_START_RANGE(NVTX_LOBPCG2_LINE)

        if ( ierr /= 0 ) then
          !ABI_COMMENT("Consider using more bands and nbdbuf if necessary.")
          ierr = 0
        end if

        !write(*,*) "    -> Iteration ", iline

        ! Compute AX-Lambda*BX
        call lobpcg_getResidu(lobpcg,eigenvaluesN)

        ! Compute residu norm here !
        call timab(tim_maxres,1,tsec)
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock)
        call timab(tim_maxres,2,tsec)

        ! Apply preconditioner
        call timab(tim_pcond,1,tsec)
        call xgBlock_apply_diag(lobpcg%W,pcond,nspinor)
        call timab(tim_pcond,2,tsec)

        call timab(tim_nbdbuf,1,tsec)
        if (nbdbuf>=0) then
          ! There is a transfer from GPU to CPU in this copy
          call xgBlock_copy(residuBlock,residu_eff%self)
          iband_min = 1 + blockdim*(iblock-1)
          iband_max = blockdim*iblock
          if (iband_max<=nband_eff) then ! all bands of this block are below nband_eff
            call xgBlock_minmax(residu_eff%self,minResidu,maxResidu)
          else if (iband_min<=nband_eff) then ! some bands of this block are below nband_eff
            call xgBlock_minmax(residu_eff%self,minResidu,maxResidu,row_bound=(nband_eff-iband_min+1))
          else ! all bands of this block are above nband_eff
            minResidu = 0.0
            maxResidu = 0.0
          end if
        else if (nbdbuf==-101) then
          call xgBlock_minmax(residuBlock,minResidu,dummy) ! Get minimum of true residuals
          ! Compute effective residuals : res_eff = res * occ
          call xgBlock_apply_diag(residuBlock,occBlock,1,Y=residu_eff%self)
          call xgBlock_minmax(residu_eff%self,dummy,maxResidu) ! Get maximum of effective residuals
        else
          ABI_ERROR('Bad value of nbdbuf')
        end if
        call timab(tim_nbdbuf,2,tsec)
        if ( maxResidu < lobpcg%tolerance ) then
          compute_residu = .false.
          ABI_NVTX_END_RANGE()
          exit
        end if

        ! Orthonormalize with respect to previous blocks
        if ( iblock > 1 ) then
          call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%W,iblock)
        end if

        if (lobpcg%paral_kgb == 1) then
          call timab(tim_transpose,1,tsec)
          call xgTransposer_transpose(lobpcg%xgTransposerW,STATE_COLSROWS)
          lobpcg%xgTransposerAW%state=STATE_COLSROWS
          lobpcg%xgTransposerBW%state=STATE_COLSROWS
          call timab(tim_transpose,2,tsec)
        end if
        ! Apply A and B on W
        call timab(tim_ax_bx,1,tsec)
        call getAX_BX(lobpcg%WColsRows,lobpcg%AWColsRows,lobpcg%BWColsRows)
        call xgBlock_zero_im_g0(lobpcg%AWColsRows)
        call xgBlock_zero_im_g0(lobpcg%BWColsRows)
        call timab(tim_ax_bx,2,tsec)
        if (lobpcg%paral_kgb == 1) then
          call timab(tim_transpose,1,tsec)
          call xgTransposer_transpose(lobpcg%xgTransposerW,STATE_LINALG)
          call xgTransposer_transpose(lobpcg%xgTransposerAW,STATE_LINALG)
          call xgTransposer_transpose(lobpcg%xgTransposerBW,STATE_LINALG)
          call timab(tim_transpose,2,tsec)
        end if

        ! DO RR in the correct subspace
        ! if residu starts to be too small, there is an accumulation error in
        ! P with values such as 1e-29 that make the eigenvectors diverge
        if ( iline == 1 .or. minResidu < 1e-27) then
          ! Do RR on XW to get the eigen vectors
          call xg_Borthonormalize(lobpcg%XW,lobpcg%BXW,ierr,tim_Bortho_XW,lobpcg%gpu_option,AX=lobpcg%AXW) ! Do rotate AW
          call xgBlock_zero(lobpcg%P)
          call xgBlock_zero(lobpcg%AP)
          call xgBlock_zero(lobpcg%BP)
          if ( ierr /= 0 ) then
            ABI_COMMENT("B-orthonormalization (XW) did not work.")
          end if
          call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvalues2N,ierr,lobpcg%prtvol,tim_RR_XW,lobpcg%gpu_option,&
           & tolerance=tolerance,&
           & XW=lobpcg%XW,AW=lobpcg%AW,BW=lobpcg%BW,P=lobpcg%P,AP=lobpcg%AP,BP=lobpcg%BP,WP=lobpcg%WP,&
           & AWP=lobpcg%AWP,BWP=lobpcg%BWP)
          if ( ierr /= 0 ) then
            ABI_WARNING("RayleighRitz (XW) did not work, but continue anyway.")
            ABI_NVTX_END_RANGE()
            exit
          end if
        else
          ! B-orthonormalize P, BP
          call xg_Borthonormalize(lobpcg%XWP%self,lobpcg%BXWP%self,ierr,tim_Bortho_XWP,lobpcg%gpu_option,AX=lobpcg%AXWP%self) ! Do rotate AW
          ! Do RR on XWP to get the eigen vectors
          if ( ierr == 0 ) then
            call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvalues3N%self,ierr,lobpcg%prtvol,tim_RR_XWP,lobpcg%gpu_option,&
           & tolerance=tolerance,XW=lobpcg%XW,AW=lobpcg%AW,BW=lobpcg%BW,P=lobpcg%P,AP=lobpcg%AP,BP=lobpcg%BP,WP=lobpcg%WP,&
           & AWP=lobpcg%AWP,BWP=lobpcg%BWP,XWP=lobpcg%XWP%self)
            if ( ierr /= 0 ) then
              ABI_WARNING("RayleighRitz (XWP) did not work, but continue anyway.")
              ABI_NVTX_END_RANGE()
              exit
            end if
          else
            ABI_COMMENT("B-orthonormalization (XWP) did not work, try on XW.")
            call xg_Borthonormalize(lobpcg%XW,lobpcg%BXW,ierr,tim_Bortho_XW,lobpcg%gpu_option,AX=lobpcg%AXW) ! Do rotate AW
            if ( ierr /= 0 ) then
              ABI_COMMENT("B-orthonormalization (XW) did not work.")
            end if
            call xgBlock_zero(lobpcg%P)
            call xgBlock_zero(lobpcg%AP)
            call xgBlock_zero(lobpcg%BP)
            nrestart = nrestart + 1
            call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvalues2N,ierr,lobpcg%prtvol,tim_RR_XW,lobpcg%gpu_option,&
           & tolerance=tolerance,&
           & XW=lobpcg%XW,AW=lobpcg%AW,BW=lobpcg%BW,P=lobpcg%P,AP=lobpcg%AP,BP=lobpcg%BP,WP=lobpcg%WP,&
           & AWP=lobpcg%AWP,BWP=lobpcg%BWP)
            if ( ierr /= 0 ) then
              ABI_WARNING("RayleighRitz (XWP) did not work, but continue anyway.")
              ABI_NVTX_END_RANGE()
              exit
            end if
          end if
        end if

        ABI_NVTX_END_RANGE()
      end do

      if ( compute_residu ) then
        ! Recompute AX-Lambda*BX for the last time
        call lobpcg_getResidu(lobpcg,eigenvaluesN)
        ! Recompute residu norm here !
        call timab(tim_maxres,1,tsec)
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock)
        call timab(tim_maxres,2,tsec)
        ! Apply preconditioner
        call timab(tim_pcond,1,tsec)
        call xgBlock_apply_diag(lobpcg%W,pcond,nspinor)
        call timab(tim_pcond,2,tsec)

        call timab(tim_nbdbuf,1,tsec)
        if(lobpcg%gpu_option==ABI_GPU_OPENMP) call xgBlock_copy_from_gpu(residuBlock)
        if (nbdbuf>=0) then
          call xgBlock_copy(residuBlock,residu_eff%self)
          iband_min = 1 + blockdim*(iblock-1)
          iband_max = blockdim*iblock
          if (iband_max<=nband_eff) then ! all bands of this block are below nband_eff
            call xgBlock_minmax(residu_eff%self,minResidu,maxResidu)
          else if (iband_min<=nband_eff) then ! some bands of this block are below nband_eff
            call xgBlock_minmax(residu_eff%self,minResidu,maxResidu,row_bound=(nband_eff-iband_min+1))
          else ! all bands of this block are above nband_eff
            minResidu = 0.0
            maxResidu = 0.0
          end if
        else if (nbdbuf==-101) then
          call xgBlock_minmax(residuBlock,minResidu,dummy) ! Get minimum of true residuals
          ! Compute effective residuals : res_eff = res * occ
          call xgBlock_apply_diag(residuBlock,occBlock,1,Y=residu_eff%self)
          call xgBlock_minmax(residu_eff%self,dummy,maxResidu) ! Get maximum of effective residuals
        else
          ABI_ERROR('Bad value of nbdbuf')
        end if
        call timab(tim_nbdbuf,2,tsec)
      end if

      if (prtvol==5.and.xmpi_comm_rank(lobpcg%spacecom)==0) then
        write(msg,'(6(a,i4),2(a,es16.6))') 'lobpcg | istep=',istep,'| isppol=',isppol,'| ikpt=',ikpt,&
          & '| inonsc=',inonsc,'| iblock=',iblock,'| nline_done=',iline-1,'| minRes=',minResidu,'| maxRes=',maxResidu
        call wrtout(std_out,msg,'PERS')
      end if

      ! Save eigenvalues
      call timab(tim_copy,1,tsec)
      call xgBlock_setBlock(eigen,eigenBlock,blockdim,1,fcol=iblock)
      call xgBlock_copy(eigenvaluesN,eigenBlock)
      call timab(tim_copy,2,tsec)

      ! Save new X in X0
      call lobpcg_setX0(lobpcg,iblock)

      ! Copy previous BX into BX0 for previous block
      if ( nblock > 1 ) then
        call lobpcg_transferAX_BX(lobpcg,iblock)
      end if

      ABI_NVTX_END_RANGE()
    end do !! End iblock loop

    call xgBlock_reshape(eigen,blockdim*nblock,1)
    call xgBlock_reshape(residu,blockdim*nblock,1)
    call xgBlock_reshape(occ,blockdim*nblock,1)

    call xg_free(eigenvalues3N)
    call xg_free(residu_eff)

    skip = .false.
    if ( ierr /= 0 ) then
      ABI_COMMENT("Some errors happened, so H|Psi> and S|Psi> are computed before leaving")
      if ( lobpcg%paral_kgb == 1 ) then
        call xgTransposer_constructor(lobpcg%xgTransposerAllX0,lobpcg%AllX0,lobpcg%AllX0ColsRows,nspinor,&
          STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0,lobpcg%me_g0_fft,&
          gpu_option=lobpcg%gpu_option)
        call xgTransposer_copyConstructor(lobpcg%xgTransposerAllAX0,lobpcg%xgTransposerAllX0,&
          lobpcg%AllAX0%self,lobpcg%AllAX0ColsRows,STATE_LINALG)
        call xgTransposer_copyConstructor(lobpcg%xgTransposerAllBX0,lobpcg%xgTransposerAllX0,&
          lobpcg%AllBX0%self,lobpcg%AllBX0ColsRows,STATE_LINALG)
      else
        call xgBlock_setBlock(lobpcg%AllX0      , lobpcg%AllX0ColsRows , spacedim, lobpcg%neigenpairs)
        call xgBlock_setBlock(lobpcg%AllAX0%self, lobpcg%AllAX0ColsRows, spacedim, lobpcg%neigenpairs)
        call xgBlock_setBlock(lobpcg%AllBX0%self, lobpcg%AllBX0ColsRows, spacedim, lobpcg%neigenpairs)
      end if
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerAllX0,STATE_COLSROWS)
        lobpcg%xgTransposerAllAX0%state=STATE_COLSROWS
        lobpcg%xgTransposerAllBX0%state=STATE_COLSROWS
        call timab(tim_transpose,2,tsec)
      end if
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(lobpcg%AllX0ColsRows,lobpcg%AllAX0ColsRows,lobpcg%AllBX0ColsRows)
      call xgBlock_zero_im_g0(lobpcg%AllAX0ColsRows)
      call xgBlock_zero_im_g0(lobpcg%AllBX0ColsRows)
      call timab(tim_ax_bx,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerAllX0,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAllAX0,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAllBX0,STATE_LINALG)
        call timab(tim_transpose,2,tsec)
      end if
      call xgTransposer_free(lobpcg%xgTransposerAllX0)
      call xgTransposer_free(lobpcg%xgTransposerAllAX0)
      call xgTransposer_free(lobpcg%xgTransposerAllBX0)
      skip = .true.
    end if

    if (.not.skip) then
      if ( nblock > 1 ) then
        call xg_Borthonormalize(X0,lobpcg%AllBX0%self,ierr,tim_Bortho_Xall,&
          & lobpcg%gpu_option,AX=lobpcg%AllAX0%self) ! Do rotate AX
        call xg_RayleighRitz(X0,lobpcg%AllAX0%self,lobpcg%AllBX0%self,eigen,ierr,lobpcg%prtvol,tim_RR_Xall,&
          & lobpcg%gpu_option,tolerance=tolerance)
      end if
    end if

    if ( lobpcg%paral_kgb == 1 ) then
      call xgTransposer_free(lobpcg%xgTransposerX)
      call xgTransposer_free(lobpcg%xgTransposerAX)
      call xgTransposer_free(lobpcg%xgTransposerBX)
      call xgTransposer_free(lobpcg%xgTransposerW)
      call xgTransposer_free(lobpcg%xgTransposerAW)
      call xgTransposer_free(lobpcg%xgTransposerBW)
    end if

!    call timab(tim_run,2,tsec)

  end subroutine lobpcg_run


  subroutine lobpcg_getX0(lobpcg,iblock)

    type(lobpcg_t), intent(inout) :: lobpcg
    integer       , intent(in   ) :: iblock
    integer :: blockdim
    integer :: spacedim
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    !lobpcg%XWP(:,X+1:X+blockdim) = lobpcg%X0(:,(iblock-1)*blockdim+1:iblock*blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,lobpcg%X0,spacedim,blockdim,fcol=(iblock-1)*blockdim+1)
    call xgBlock_copy(lobpcg%X0,lobpcg%X)

    call timab(tim_copy,2,tsec)

  end subroutine lobpcg_getX0


  subroutine lobpcg_setPreviousX0_BX0(lobpcg,iblock)

    type(lobpcg_t) , intent(inout) :: lobpcg
    integer        , intent(in   ) :: iblock

    if (iblock<2) then
      ABI_ERROR("iblock<2")
    end if
    call xg_setBlock(lobpcg%AllBX0,lobpcg%BX0,lobpcg%spacedim,(iblock-1)*lobpcg%blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,lobpcg%X0,lobpcg%spacedim,(iblock-1)*lobpcg%blockdim)
  end subroutine lobpcg_setPreviousX0_BX0


  subroutine lobpcg_orthoXwrtBlocks(lobpcg,var,iblock)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: var
    integer        , intent(in   ) :: iblock
    integer :: previousBlock
    integer :: blockdim
    integer :: spacedim
    integer :: space_buf
    type(xg_t) :: buffer
    double precision :: tsec(2)

    call timab(tim_ortho,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_LOBPCG2_ORTHO_X_WRT)

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim
    previousBlock = (iblock-1)*lobpcg%blockdim

    space_buf = space(var)
    if (space(var)==SPACE_CR) then
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,previousBlock,blockdim,comm=lobpcg%spacecom,gpu_option=lobpcg%gpu_option)

    ! buffer = BX0^T*X
    call xgBlock_gemm('t','n',1.0d0,lobpcg%BX0,var,0.d0,buffer%self,comm=lobpcg%spacecom)

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm('n','n',-1.0d0,lobpcg%X0,buffer%self,1.0d0,var)

    call xg_free(buffer)

    ABI_NVTX_END_RANGE()
    call timab(tim_ortho,2,tsec)

  end subroutine lobpcg_orthoXwrtBlocks

  subroutine lobpcg_getResidu(lobpcg,eigenvalues)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(in   ) :: eigenvalues
    double precision :: tsec(2)

    call timab(tim_maxres,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_LOBPCG2_RESIDUE)
      !lobpcg%XWP(1:spacedim,shiftW+iblock) = lobpcg%AXWP(:,shiftX+iblock) - lobpcg%BXWP(:,shiftX+iblock)*eigenvalues(iblock)
    call xgBlock_colwiseCymax(lobpcg%W,eigenvalues,lobpcg%BX,lobpcg%AX)
    ABI_NVTX_END_RANGE()
    call timab(tim_maxres,2,tsec)
  end subroutine lobpcg_getResidu

  subroutine lobpcg_setX0(lobpcg,iblock)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: iblock
    type(xgBlock_t) :: Xtmp
    integer :: blockdim
    integer :: spacedim
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)
    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    !X0(:,(iblock-1)*blockdim+1:iblock*blockdim) = lobpcg%XWP(:,lobpcg%X+1:lobpcg%X+blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,Xtmp,spacedim,blockdim,fcol=(iblock-1)*blockdim+1)
    call xgBlock_copy(lobpcg%X,Xtmp)
    call timab(tim_copy,2,tsec)

  end subroutine lobpcg_setX0


  subroutine lobpcg_transferAX_BX(lobpcg,iblock)

    type(lobpcg_t), intent(inout) :: lobpcg
    integer       , intent(in   ) :: iblock
    type(xgBlock_t) :: CXtmp
    integer :: firstcol
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)

    if (iblock<1) then
      ABI_ERROR("iblock<1")
    end if

    ! iblock goes from 1 to nblock-1 included
    firstcol = (iblock-1)*lobpcg%blockdim+1  ! Start of each block

    ! BX
    call xg_setBlock(lobpcg%AllBX0,CXtmp,lobpcg%spacedim,lobpcg%blockdim,fcol=firstcol)
    call xgBlock_copy(lobpcg%BX,CXtmp)

    ! AX
    call xg_setBlock(lobpcg%AllAX0,CXtmp,lobpcg%spacedim,lobpcg%blockdim,fcol=firstcol)
    call xgBlock_copy(lobpcg%AX,CXtmp)

    call timab(tim_copy,2,tsec)

  end subroutine lobpcg_transferAX_BX

  subroutine lobpcg_free(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

    call xg_free(lobpcg%XWP)
    call xg_free(lobpcg%AXWP)
    call xg_free(lobpcg%BXWP)
    call xg_free(lobpcg%AllAX0)
    call xg_free(lobpcg%AllBX0)
  end subroutine lobpcg_free

end module m_lobpcg2
!!***

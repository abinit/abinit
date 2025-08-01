
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lobpcg2_cprj

  use defs_basis
  use m_xg
  use m_xg_nonlop
  use m_xg_ortho_RR
  use m_xgTransposer
  use m_xgScalapack
  use defs_basis
  use m_abicore
  use m_errors
  use m_xomp
#ifdef HAVE_OPENMP
  use omp_lib
#endif

  use m_xmpi
  use m_time, only : timab

  implicit none

  private

  integer, parameter :: tim_init     = 2040
  integer, parameter :: tim_free     = 2041
  integer, parameter :: tim_cprj     = 2043
  integer, parameter :: tim_ortho    = 2044
  integer, parameter :: tim_transpose = 2039
  integer, parameter :: tim_maxres   = 2046
  integer, parameter :: tim_ax_k     = 2048
  integer, parameter :: tim_ax_v     = 2049
  integer, parameter :: tim_ax_nl    = 2050
  integer, parameter :: tim_pcond    = 2047

  integer, parameter :: tim_Bortho_X    = 2031
  integer, parameter :: tim_Bortho_XW   = 2032
  integer, parameter :: tim_Bortho_XWP  = 2033
  integer, parameter :: tim_Bortho_Xall = 2034
  integer, parameter :: tim_RR_X        = 2035
  integer, parameter :: tim_RR_XW       = 2036
  integer, parameter :: tim_RR_XWP      = 2037
  integer, parameter :: tim_RR_Xall     = 2038

  integer, parameter :: tim_copy     = 2042
  integer, parameter :: tim_nbdbuf   = 2045

  integer, parameter :: tim_enl = 2051

  type, public :: lobpcg_t
    logical :: is_nested                     ! For OpenMP nested region
    integer :: spacedim                      ! Space dimension for one vector (WFs)
    integer :: space                         ! Space type for one vector
    integer :: space_cprj                    ! Space type for one projector
    integer :: cprjdim                       ! Space dimension for one vector (cprj)
    integer :: neigenpairs                   ! Number of eigen values/vectors we want
    integer :: nblock                        ! Number of block in the space of dim spacedim
    integer :: blockdim                      ! Number of vectors in one block
    integer :: blockdim_cprj                 ! Number of vectors in one block (for cprj array)
    integer :: nline                         ! Number of line to perform
    integer :: spacecom                      ! Communicator for MPI
    integer :: paral_kgb                     ! paral_kgb formalism or not
    integer :: comm_rows                     ! communicator for rows
    integer :: comm_cols                     ! communicator for cols
    integer :: me_g0
    integer :: me_g0_fft
    double precision :: tolerance            ! Tolerance on the residu to stop the minimization
    integer :: prtvol
    type(xg_nonlop_t) :: xg_nonlop ! contains data for nonlocal operations
    type(xgBlock_t) :: AllX0 ! Block of initial and final solution.
    type(xgBlock_t) :: AllcprjX0 ! Block of initial and final solution.
    type(xgBlock_t) :: X0
    type(xgBlock_t) :: cprjX0
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
    type(xg_t) :: cprjXWP

    type(xg_t) :: cprj_work
    type(xg_t) :: cprj_work2

    type(xgBlock_t) :: XColsRows
    type(xgBlock_t) :: AXColsRows

    type(xgTransposer_t) :: xgTransposerX
    type(xgTransposer_t) :: xgTransposerAX

    type(xgBlock_t) :: WColsRows
    type(xgBlock_t) :: AWColsRows

    type(xgTransposer_t) :: xgTransposerW
    type(xgTransposer_t) :: xgTransposerAW

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

    type(xgBlock_t) :: cprjX ! Shift to apply to start reading the X values
    type(xgBlock_t) :: cprjW ! Shift to apply to start reading the W values
    type(xgBlock_t) :: cprjP ! Shift to apply to start reading the P values
    type(xgBlock_t) :: cprjXW ! Shift to apply to start reading the P values
    type(xgBlock_t) :: cprjWP ! Shift to apply to start reading the P values

    type(xg_t) :: AllAX0

    type(xgBlock_t) :: AllX0ColsRows
    type(xgBlock_t) :: AllAX0ColsRows

    type(xgTransposer_t) :: xgTransposerAllX0
    type(xgTransposer_t) :: xgTransposerAllAX0

  end type lobpcg_t

  public :: lobpcg_init
  public :: lobpcg_memInfo
  public :: lobpcg_run_cprj
  public :: lobpcg_free

  contains


  subroutine lobpcg_init(lobpcg, bandpp, neigenpairs, spacedim, cprjdim, blockdim, tolerance, nline, &
&      space, space_cprj, spacecom, paral_kgb, xg_nonlop, comm_rows, comm_cols, me_g0, me_g0_fft)

    type(lobpcg_t)   , intent(inout) :: lobpcg
    type(xg_nonlop_t), intent(in   ) :: xg_nonlop
    integer          , intent(in   ) :: neigenpairs
    integer          , intent(in   ) :: spacedim
    integer          , intent(in   ) :: cprjdim
    integer          , intent(in   ) :: blockdim
    integer          , intent(in   ) :: bandpp
    double precision , intent(in   ) :: tolerance
    integer          , intent(in   ) :: nline
    integer          , intent(in   ) :: comm_rows,comm_cols
    integer          , intent(in   ) :: space
    integer          , intent(in   ) :: space_cprj
    integer          , intent(in   ) :: spacecom
    integer          , intent(in   ) :: paral_kgb
    integer          , intent(in   ) :: me_g0
    integer          , intent(in   ) :: me_g0_fft
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)
    lobpcg%neigenpairs   = neigenpairs
    lobpcg%spacedim      = spacedim
    lobpcg%space         = space
    lobpcg%space_cprj    = space_cprj
    lobpcg%cprjdim       = cprjdim
    lobpcg%nblock        = neigenpairs / blockdim
    lobpcg%blockdim      = blockdim
    lobpcg%blockdim_cprj = bandpp*xg_nonlop%nspinor
    if (tolerance > 0.0) then
      lobpcg%tolerance   = tolerance
    else
      lobpcg%tolerance   = 1e-20
    end if
    lobpcg%nline         = nline
    lobpcg%spacecom      = spacecom
    lobpcg%paral_kgb     = paral_kgb
    lobpcg%comm_rows     = comm_rows
    lobpcg%comm_cols     = comm_cols
    lobpcg%xg_nonlop     = xg_nonlop
    lobpcg%me_g0         = me_g0
    lobpcg%me_g0_fft     = me_g0_fft

    call lobpcg_allocateAll(lobpcg)
    call timab(tim_init,2,tsec)

  end subroutine lobpcg_init


  subroutine lobpcg_allocateAll(lobpcg)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer :: spacedim
    integer :: cprjdim
    integer :: blockdim
    integer :: blockdim_cprj

    spacedim = lobpcg%spacedim
    blockdim = lobpcg%blockdim

    cprjdim       = lobpcg%cprjdim
    blockdim_cprj = lobpcg%blockdim_cprj

    call lobpcg_free(lobpcg) ! Make sure everything is not allocated and
    ! pointer point to null()

    call xg_init(lobpcg%XWP,lobpcg%space,spacedim,3*blockdim,lobpcg%spacecom,me_g0=lobpcg%me_g0)
    call xg_setBlock(lobpcg%XWP,lobpcg%X,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%W,spacedim,blockdim,fcol=blockdim+1)
    call xg_setBlock(lobpcg%XWP,lobpcg%P,spacedim,blockdim,fcol=2*blockdim+1)
    call xg_setBlock(lobpcg%XWP,lobpcg%XW,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%WP,spacedim,2*blockdim,fcol=blockdim+1)

    call xg_init(lobpcg%AXWP,lobpcg%space,spacedim,3*blockdim,lobpcg%spacecom,me_g0=lobpcg%me_g0)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AX,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AW,spacedim,blockdim,fcol=blockdim+1)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AP,spacedim,blockdim,fcol=2*blockdim+1)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AXW,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AWP,spacedim,2*blockdim,fcol=blockdim+1)

    call xg_init(lobpcg%cprjXWP,lobpcg%space_cprj,cprjdim,3*blockdim_cprj,lobpcg%spacecom)
    call xg_setBlock(lobpcg%cprjXWP,lobpcg%cprjX, cprjdim,  blockdim_cprj)
    call xg_setBlock(lobpcg%cprjXWP,lobpcg%cprjW, cprjdim,  blockdim_cprj,fcol=   blockdim_cprj+1)
    call xg_setBlock(lobpcg%cprjXWP,lobpcg%cprjP, cprjdim,  blockdim_cprj,fcol= 2*blockdim_cprj+1)
    call xg_setBlock(lobpcg%cprjXWP,lobpcg%cprjXW,cprjdim,2*blockdim_cprj)
    call xg_setBlock(lobpcg%cprjXWP,lobpcg%cprjWP,cprjdim,2*blockdim_cprj,fcol=   blockdim_cprj+1)

    call xg_init(lobpcg%cprj_work ,lobpcg%space_cprj,cprjdim,blockdim_cprj,lobpcg%spacecom)
    call xg_init(lobpcg%cprj_work2,lobpcg%space_cprj,cprjdim,blockdim_cprj,lobpcg%spacecom)

    if ( lobpcg%nblock /= 1 ) then
      call xg_init(lobpcg%AllAX0,lobpcg%space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom,me_g0=lobpcg%me_g0)
    else
      lobpcg%AllAX0%self = lobpcg%AX
    end if

  end subroutine lobpcg_allocateAll


  function lobpcg_memInfo(neigenpairs, spacedim, blockdim, space) result(arraymem)

    integer         , intent(in   ) :: neigenpairs
    integer         , intent(in   ) :: spacedim
    integer         , intent(in   ) :: blockdim
    integer         , intent(in   ) :: space
    double precision :: memXWP
    double precision :: memAXWP
    double precision :: memBXWP
    double precision :: memAllBX0
    double precision :: memAllAX0
    double precision :: memeigenvalues3N
    double precision :: membufferOrtho
    double precision :: membufferBOrtho
    double precision :: memsubA
    double precision :: memsubB
    double precision :: memsubBtmp
    double precision :: memvec
    double precision :: memRR
    double precision :: maxmemTmp
    double precision :: cplx
    integer :: nblock
    double precision :: arraymem(2)

    cplx = 1 ; if ( space == SPACE_C ) cplx = 2
    nblock = neigenpairs/blockdim

    ! Permanent in lobpcg
    memXWP  = cplx* kind(1.d0) * spacedim * 3*blockdim
    memAXWP = cplx* kind(1.d0) * spacedim * 3*blockdim
    memBXWP = cplx* kind(1.d0) * spacedim * 3*blockdim
    if ( nblock > 1 ) then
      memAllAX0 = cplx * kind(1.d0) * spacedim * 3*blockdim
      memAllBX0 = cplx * kind(1.d0) * spacedim * 3*blockdim
      membufferOrtho = cplx * kind(1.d0) * blockdim * (nblock-1) * blockdim
    else
      memAllAX0 = 0
      memAllBX0 = 0
      membufferOrtho = 0
    endif
    memeigenvalues3N = kind(1.d0) * 3*blockdim

    ! Temporary arrays
    membufferBOrtho = cplx * kind(1.d0) * 2*blockdim * 2*blockdim ! For the moment being, only Bortho with X or WP at the same time
    memsubA = cplx * kind(1.d0) * 3*blockdim * 3*blockdim
    memsubB = cplx * kind(1.d0) * 3*blockdim * 3*blockdim
    memsubBtmp = cplx * kind(1.d0) * spacedim * blockdim
    memvec = cplx * kind(1.d0) * 3*blockdim * blockdim
    memRR = max(memsubA+memsubB+memvec,memsubBtmp+memvec)

    maxmemTmp = max( membufferBOrtho,memRR,membufferOrtho )

    arraymem(1) = memXWP+memAXWP+memBXWP+memAllAX0+memAllBX0+memeigenvalues3N
    arraymem(2) = maxmemTmp

  end function lobpcg_memInfo


  subroutine lobpcg_run_cprj(lobpcg, X0, cprjX0, getAX, kin, pcond, eigen, occ, residu, enl, prtvol, nspinor, &
      isppol, ikpt, inonsc, istep, nbdbuf)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: X0      ! Full initial vectors
    type(xgBlock_t), intent(inout) :: cprjX0  ! Full initial cprj
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: occ
    type(xgBlock_t), intent(inout) :: residu
    type(xgBlock_t), intent(inout) :: enl
    type(xgBlock_t), intent(in   ) :: kin
    type(xgBlock_t), intent(in   ) :: pcond
    integer        , intent(in   ) :: prtvol
    integer        , intent(in   ) :: nspinor
    integer        , intent(in   ) :: isppol,ikpt,inonsc,istep,nbdbuf

    type(xg_t) :: eigenvalues3N   ! eigen values for Rayleight-Ritz
    type(xg_t) :: residu_eff
    type(xgBlock_t) :: eigenvaluesN   ! eigen values for Rayleight-Ritz
    type(xgBlock_t) :: eigenvalues2N   ! eigen values for Rayleight-Ritz
    logical :: skip,compute_residu
    integer :: blockdim, blockdim3, blockdim2, blockdim_cprj
    integer :: spacedim
    integer :: iblock, nblock
    integer :: iline, nline
    integer :: rows_tmp, cols_tmp, nband_eff, iband_min, iband_max
    integer,parameter :: gpu_option=ABI_GPU_DISABLED
    type(xgBlock_t) :: eigenBlock   !
    type(xgBlock_t) :: residuBlock,occBlock
    type(xg_t):: cprj_work_all
    double precision :: maxResidu, minResidu, dummy
    double precision :: dlamch,tolerance
    integer :: ierr = 0
    integer :: nrestart
    double precision :: tsec(2)
    character(len=500) :: msg

    interface
      subroutine getAX(X,AX)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
      end subroutine getAX
    end interface

    lobpcg%prtvol = prtvol

    tolerance=2*dlamch('E')

    blockdim = lobpcg%blockdim
    blockdim2 = 2*blockdim
    blockdim3 = 3*blockdim
    blockdim_cprj = lobpcg%blockdim_cprj
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

    call xg_init(eigenvalues3N,SPACE_R,blockdim3,1)
    call xg_setBlock(eigenvalues3N,eigenvaluesN,blockdim,1)
    call xg_setBlock(eigenvalues3N,eigenvalues2N,blockdim2,1)

    call xgBlock_reshape(eigen,blockdim,nblock)
    call xgBlock_reshape(residu,blockdim,nblock)
    call xgBlock_reshape(occ,blockdim,nblock)

    lobpcg%AllX0     = X0
    lobpcg%AllcprjX0 = cprjX0

    call xg_init(residu_eff,SPACE_R,blockdim,1,gpu_option=ABI_GPU_DISABLED)

    if ( lobpcg%paral_kgb == 1 ) then
      call timab(tim_transpose,1,tsec)
      call xgTransposer_constructor(lobpcg%xgTransposerX,lobpcg%X,lobpcg%XColsRows,nspinor,&
        STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0,lobpcg%me_g0_fft)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerAX,lobpcg%xgTransposerX,&
        lobpcg%AX,lobpcg%AXColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerW,lobpcg%xgTransposerX,&
        lobpcg%W,lobpcg%WColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(lobpcg%xgTransposerAW,lobpcg%xgTransposerX,&
        lobpcg%AW,lobpcg%AWColsRows,STATE_LINALG)
      call timab(tim_transpose,2,tsec)
    else
      call xgBlock_setBlock(lobpcg%X, lobpcg%XColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AX, lobpcg%AXColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%W, lobpcg%WColsRows, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AW, lobpcg%AWColsRows, spacedim, blockdim)
    end if

    !! Start big loop over blocks
    do iblock = 1, nblock
      nrestart = 0

      call lobpcg_getX0(lobpcg,iblock)
      call xgBlock_setBlock(residu,residuBlock,blockdim,1,fcol=iblock)
      call xgBlock_setBlock(occ,   occBlock,   blockdim,1,fcol=iblock)

      call timab(tim_cprj,1,tsec)
      call xg_nonlop_getcprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%cprj_work%self)
      call timab(tim_cprj,2,tsec)

      if ( iblock > 1 ) then
        call lobpcg_setPreviousX0(lobpcg,iblock)
        ! Orthogonalize current iblock X block With Respect To previous Blocks in B-basis
        call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%X,lobpcg%cprjX,iblock,lobpcg%cprj_work%self)
      end if

      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerX,STATE_COLSROWS)
        lobpcg%xgTransposerAX%state=STATE_COLSROWS
        call timab(tim_transpose,2,tsec)
      end if

      ! Initialize some quantitites (AX)
      call timab(tim_ax_v,1,tsec)
      call getAX(lobpcg%XColsRows,lobpcg%AXColsRows)
      call xgBlock_zero_im_g0(lobpcg%AXColsRows)
      call timab(tim_ax_v,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        lobpcg%xgTransposerX%state=STATE_LINALG
        call xgTransposer_transpose(lobpcg%xgTransposerAX,STATE_LINALG)
        call timab(tim_transpose,2,tsec)
      end if

      call timab(tim_ax_k,1,tsec)
      call xgBlock_add_diag(lobpcg%X,kin,nspinor,lobpcg%AX)
      call timab(tim_ax_k,2,tsec)

      ! B-orthonormalize X, BX and AX
      call xg_Borthonormalize_cprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,ierr,tim_Bortho_X,gpu_option,AX=lobpcg%AX)

      ! Do first RR on X to get the first eigen values
      call xg_RayleighRitz_cprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%AX,eigenvaluesN,ierr,&
        & lobpcg%prtvol,tim_RR_X,gpu_option,add_Anl=.True.)

      compute_residu = .true.

      do iline = 1, nline

        if ( ierr /= 0 ) then
          !ABI_COMMENT("Consider using more bands and nbdbuf if necessary.")
          ierr = 0
        end if

        ! Compute AX-Lambda*BX
        call timab(tim_copy,1,tsec)
        call xgBlock_copy(lobpcg%AX,lobpcg%W)
        call timab(tim_copy,2,tsec)
        ! Add the non-local part: (A_nl - Lambda*B)X
        call timab(tim_ax_nl,1,tsec)
        if (lobpcg%xg_nonlop%paw) then
          call xg_nonlop_getHmeSX(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%W,eigenvaluesN,&
            lobpcg%cprj_work%self,lobpcg%cprj_work2%self)
        else
          call xg_nonlop_getHX(lobpcg%xg_nonlop,lobpcg%W,lobpcg%cprjX,&
            lobpcg%cprj_work%self,lobpcg%cprj_work2%self)
          call xgBlock_yxmax(lobpcg%W,eigenvaluesN,lobpcg%X)
        end if
        call timab(tim_ax_nl,2,tsec)

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
          exit
        end if

        call timab(tim_cprj,1,tsec)
        call xg_nonlop_getcprj(lobpcg%xg_nonlop,lobpcg%W,lobpcg%cprjW,lobpcg%cprj_work%self)
        call timab(tim_cprj,2,tsec)

        ! Orthonormalize with respect to previous blocks
        if ( iblock > 1 ) then
          call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%W,lobpcg%cprjW,iblock,lobpcg%cprj_work%self)
        end if

        if (lobpcg%paral_kgb == 1) then
          call timab(tim_transpose,1,tsec)
          call xgTransposer_transpose(lobpcg%xgTransposerW,STATE_COLSROWS)
          lobpcg%xgTransposerAW%state=STATE_COLSROWS
          call timab(tim_transpose,2,tsec)
        end if
        ! Apply A and B on W
        call timab(tim_ax_v,1,tsec)
        call getAX(lobpcg%WColsRows,lobpcg%AWColsRows)
        call xgBlock_zero_im_g0(lobpcg%AWColsRows)
        call timab(tim_ax_v,2,tsec)
        if (lobpcg%paral_kgb == 1) then
          call timab(tim_transpose,1,tsec)
          lobpcg%xgTransposerW%state=STATE_LINALG
          call xgTransposer_transpose(lobpcg%xgTransposerAW,STATE_LINALG)
          call timab(tim_transpose,2,tsec)
        end if

        call timab(tim_ax_k,1,tsec)
        call xgBlock_add_diag(lobpcg%W,kin,nspinor,lobpcg%AW)
        call timab(tim_ax_k,2,tsec)

        ! DO RR in the correct subspace
        ! if residu starts to be too small, there is an accumulation error in
        ! P with values such as 1e-29 that make the eigenvectors diverge
        if ( iline == 1 .or. minResidu < 1e-27) then
          ! Do RR on XW to get the eigen vectors
          call xg_Borthonormalize_cprj(lobpcg%xg_nonlop,lobpcg%XW,lobpcg%cprjXW,&
            & ierr,tim_Bortho_XW,gpu_option,AX=lobpcg%AXW,blockdim_cprj=blockdim_cprj)
          call xgBlock_zero(lobpcg%P)
          call xgBlock_zero(lobpcg%AP)
          call xgBlock_zero(lobpcg%cprjP)
          if ( ierr /= 0 ) then
            ABI_COMMENT("B-orthonormalization (XW) did not work.")
          end if
          call xg_RayleighRitz_cprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%AX,eigenvaluesN,ierr,&
           & lobpcg%prtvol,tim_RR_XW,gpu_option,tolerance=tolerance,&
           & XW=lobpcg%XW,W=lobpcg%W,cprjXW=lobpcg%cprjXW,cprjW=lobpcg%cprjW,AW=lobpcg%AW,&
           & P=lobpcg%P,cprjP=lobpcg%cprjP,AP=lobpcg%AP,WP=lobpcg%WP,cprjWP=lobpcg%cprjWP,&
           & AWP=lobpcg%AWP,blockdim_cprj=blockdim_cprj,add_Anl=.True.)
          if ( ierr /= 0 ) then
            ABI_WARNING("RayleighRitz (XW) did not work, but continue anyway.")
            exit
          end if
        else
          ! B-orthonormalize P, BP
          call xg_Borthonormalize_cprj(lobpcg%xg_nonlop,lobpcg%XWP%self,lobpcg%cprjXWP%self,&
            & ierr,tim_Bortho_XWP,gpu_option,AX=lobpcg%AXWP%self,blockdim_cprj=blockdim_cprj)
          ! Do RR on XWP to get the eigen vectors
          if ( ierr == 0 ) then
            call xg_RayleighRitz_cprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%AX,eigenvaluesN,ierr,&
             & lobpcg%prtvol,tim_RR_XWP,gpu_option,tolerance=tolerance,&
             & XW=lobpcg%XW,W=lobpcg%W,cprjXW=lobpcg%cprjXW,cprjW=lobpcg%cprjW,AW=lobpcg%AW,&
             & P=lobpcg%P,cprjP=lobpcg%cprjP,AP=lobpcg%AP,WP=lobpcg%WP,cprjWP=lobpcg%cprjWP,&
             & AWP=lobpcg%AWP,XWP=lobpcg%XWP%self,cprjXWP=lobpcg%cprjXWP%self,blockdim_cprj=blockdim_cprj,add_Anl=.True.)
             if ( ierr /= 0 ) then
               ABI_WARNING("RayleighRitz (XWP) did not work, but continue anyway.")
               exit
             end if
          else
            ABI_COMMENT("B-orthonormalization (XWP) did not work, try on XW.")
            call xg_Borthonormalize_cprj(lobpcg%xg_nonlop,lobpcg%XW,lobpcg%cprjXW,&
              & ierr,tim_Bortho_XW,gpu_option,AX=lobpcg%AXW,blockdim_cprj=blockdim_cprj)
            if ( ierr /= 0 ) then
              ABI_COMMENT("B-orthonormalization (XW) did not work.")
            end if
            call xgBlock_zero(lobpcg%P)
            call xgBlock_zero(lobpcg%AP)
            call xgBlock_zero(lobpcg%cprjP)
            nrestart = nrestart + 1
            call xg_RayleighRitz_cprj(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%AX,eigenvaluesN,ierr,&
             & lobpcg%prtvol,tim_RR_XW,gpu_option,tolerance=tolerance,&
             & XW=lobpcg%XW,W=lobpcg%W,cprjXW=lobpcg%cprjXW,cprjW=lobpcg%cprjW,AW=lobpcg%AW,&
             & P=lobpcg%P,cprjP=lobpcg%cprjP,AP=lobpcg%AP,WP=lobpcg%WP,cprjWP=lobpcg%cprjWP,&
             & AWP=lobpcg%AWP,blockdim_cprj=blockdim_cprj,add_Anl=.True.)
             if ( ierr /= 0 ) then
               ABI_WARNING("RayleighRitz (XW) did not work, but continue anyway.")
               exit
             end if
            end if
        end if

      end do

      if ( compute_residu ) then
        ! Recompute AX-Lambda*BX for the last time
        call timab(tim_copy,1,tsec)
        call xgBlock_copy(lobpcg%AX,lobpcg%W)
        call timab(tim_copy,2,tsec)
        call timab(tim_ax_nl,1,tsec)
        if (lobpcg%xg_nonlop%paw) then
          call xg_nonlop_getHmeSX(lobpcg%xg_nonlop,lobpcg%X,lobpcg%cprjX,lobpcg%W,eigenvaluesN,&
            lobpcg%cprj_work%self,lobpcg%cprj_work2%self)
        else
          call xg_nonlop_getHX(lobpcg%xg_nonlop,lobpcg%W,lobpcg%cprjX,&
            lobpcg%cprj_work%self,lobpcg%cprj_work2%self)
          call xgBlock_yxmax(lobpcg%W,eigenvaluesN,lobpcg%X)
        end if
        call timab(tim_ax_nl,2,tsec)
        ! Recompute residu norm here !
        call timab(tim_maxres,1,tsec)
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock)
        call timab(tim_maxres,2,tsec)
        ! Apply preconditioner
        call timab(tim_pcond,1,tsec)
        call xgBlock_apply_diag(lobpcg%W,pcond,nspinor)
        call timab(tim_pcond,2,tsec)

        call timab(tim_nbdbuf,1,tsec)
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

    end do !! End iblock loop

    call xgBlock_reshape(eigen,blockdim*nblock,1)
    call xgBlock_reshape(residu,blockdim*nblock,1)
    call xgBlock_reshape(occ,blockdim*nblock,1)

    call xg_free(eigenvalues3N)
    call xg_free(residu_eff)

    skip = .false.
    if ( ierr /= 0 ) then
      ABI_COMMENT("But before that, I want to recalculate H|Psi> and S|Psi>")
      if ( lobpcg%paral_kgb == 1 ) then
        call xgTransposer_constructor(lobpcg%xgTransposerAllX0,lobpcg%AllX0,lobpcg%AllX0ColsRows,nspinor,&
          STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0,lobpcg%me_g0_fft)
        call xgTransposer_copyConstructor(lobpcg%xgTransposerAllAX0,lobpcg%xgTransposerAllX0,&
          lobpcg%AX,lobpcg%AXColsRows,STATE_LINALG)
      else
        call xgBlock_setBlock(lobpcg%AllX0      , lobpcg%AllX0ColsRows , spacedim, lobpcg%neigenpairs)
        call xgBlock_setBlock(lobpcg%AllAX0%self, lobpcg%AllAX0ColsRows, spacedim, lobpcg%neigenpairs)
      end if
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        call xgTransposer_transpose(lobpcg%xgTransposerAllX0,STATE_COLSROWS)
        lobpcg%xgTransposerAllAX0%state=STATE_COLSROWS
        call timab(tim_transpose,2,tsec)
      end if
      call timab(tim_ax_v,1,tsec)
      call getAX(lobpcg%AllX0,lobpcg%AllAX0%self)
      call xgBlock_zero_im_g0(lobpcg%AllAX0%self)
      call timab(tim_ax_v,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call timab(tim_transpose,1,tsec)
        lobpcg%xgTransposerAllX0%state=STATE_LINALG
        call xgTransposer_transpose(lobpcg%xgTransposerAllAX0,STATE_LINALG)
        call timab(tim_transpose,2,tsec)
      end if

      call timab(tim_ax_k,1,tsec)
      call xgBlock_add_diag(lobpcg%AllX0,kin,nspinor,lobpcg%AllAX0%self)
      call timab(tim_ax_k,2,tsec)

      call xgTransposer_free(lobpcg%xgTransposerAllX0)
      call xgTransposer_free(lobpcg%xgTransposerAllAX0)
      skip = .true.
    end if

    if (.not.skip) then

      if ( nblock > 1 ) then
        call timab(tim_cprj,1,tsec)
        call xg_init(cprj_work_all,space(cprjX0),rows(cprjX0),cols(cprjX0),comm(cprjX0))
        ! cprj computation could be avoided, with MPI comm instead
        call xg_nonlop_getcprj(lobpcg%xg_nonlop,X0,cprjX0,cprj_work_all%self)
        call xg_free(cprj_work_all)
        call timab(tim_cprj,2,tsec)
        call xg_Borthonormalize_cprj(lobpcg%xg_nonlop,X0,cprjX0,&
          & ierr,tim_Bortho_Xall,gpu_option,AX=lobpcg%AllAX0%self)
        call xg_RayleighRitz_cprj(lobpcg%xg_nonlop,X0,cprjX0,lobpcg%AllAX0%self,eigen,ierr,&
          & lobpcg%prtvol,tim_RR_Xall,gpu_option,add_Anl=.True.)
      end if

    end if

    if ( lobpcg%paral_kgb == 1 ) then
      call xgTransposer_free(lobpcg%xgTransposerX)
      call xgTransposer_free(lobpcg%xgTransposerAX)
      call xgTransposer_free(lobpcg%xgTransposerW)
      call xgTransposer_free(lobpcg%xgTransposerAW)
    end if

    if (.not.lobpcg%xg_nonlop%paw) then
      call timab(tim_enl,1,tsec)
      call xg_init(cprj_work_all,space(cprjX0),rows(cprjX0),cols(cprjX0),comm(cprjX0))
      call xg_nonlop_colwiseXHX(lobpcg%xg_nonlop,cprjX0,cprj_work_all%self,enl)
      call xg_free(cprj_work_all)
      call timab(tim_enl,2,tsec)
    end if

  end subroutine lobpcg_run_cprj

  subroutine lobpcg_getX0(lobpcg,iblock)

    type(lobpcg_t) , intent(inout)  :: lobpcg
    integer        , intent(in   )  :: iblock
    type(xgBlock_t) :: Xtmp,cprjXtmp
    integer :: blockdim
    integer :: spacedim
    integer :: cprjdim
    integer :: blockdim_cprj
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    cprjdim       = lobpcg%cprjdim
    blockdim_cprj = lobpcg%blockdim_cprj

    call xgBlock_setBlock(lobpcg%AllX0,Xtmp,spacedim,blockdim,fcol=(iblock-1)*blockdim+1)
    call xgBlock_copy(Xtmp,lobpcg%X)

    call xgBlock_setBlock(lobpcg%AllcprjX0,cprjXtmp,cprjdim,blockdim_cprj,fcol=(iblock-1)*blockdim_cprj+1)
    call xgBlock_copy(cprjXtmp,lobpcg%cprjX)

    call timab(tim_copy,2,tsec)

  end subroutine lobpcg_getX0


  subroutine lobpcg_setPreviousX0(lobpcg,iblock)

    type(lobpcg_t) , intent(inout) :: lobpcg
    integer        , intent(in   ) :: iblock

    if (iblock<2) then
      ABI_ERROR("iblock<2")
    end if
    call xgBlock_setBlock(lobpcg%AllX0    ,lobpcg%X0    ,lobpcg%spacedim,(iblock-1)*lobpcg%blockdim)
    call xgBlock_setBlock(lobpcg%AllcprjX0,lobpcg%cprjX0,lobpcg%cprjdim ,(iblock-1)*lobpcg%blockdim_cprj)
  end subroutine lobpcg_setPreviousX0

  subroutine lobpcg_orthoXwrtBlocks(lobpcg,var,cprjvar,iblock,cprj_work)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: var
    type(xgBlock_t), intent(inout) :: cprjvar
    type(xgBlock_t), intent(inout) :: cprj_work
    integer        , intent(in   ) :: iblock
    integer :: previousBlock
    integer :: blockdim
    integer :: spacedim
    integer :: space_buf
    integer :: nspinor,cprjdim
    type(xg_t) :: buffer
    double precision :: tsec(2)
    type(xgBlock_t) :: cprjX0_spinor,cprj_work_spinor

    call timab(tim_ortho,1,tsec)

    if (iblock<2) then
      ABI_ERROR("iblock<2")
    end if

    if (cols(cprjvar)/=cols(cprj_work)) then
      ABI_ERROR("cprjvar and cprj_work should have same number of columns")
    end if
    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim
    previousBlock = (iblock-1)*lobpcg%blockdim

    cprjdim = lobpcg%xg_nonlop%cprjdim
    nspinor = lobpcg%xg_nonlop%nspinor

    space_buf = space(var)
    if (space(var)==SPACE_CR) then
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,previousBlock,blockdim,lobpcg%spacecom)

    ! buffer = X0^T*X
    call xgBlock_gemm('t','n',1.0d0,lobpcg%X0,var,0.d0,buffer%self,comm=lobpcg%spacecom)

    ! Add the nonlocal part if paw
    if (lobpcg%xg_nonlop%paw) then
      call xg_nonlop_getXSX(lobpcg%xg_nonlop,lobpcg%cprjX0,cprjvar,cprj_work,buffer%self)
    end if

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm('n','n',-1.0d0,lobpcg%X0,buffer%self,1.0d0,var)

    call xgBlock_zero(cprj_work)
    call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,nspinor,COLS2ROWS)
    call xgBlock_reshape_spinor(lobpcg%cprjX0,cprjX0_spinor,nspinor,COLS2ROWS)
    call xgBlock_gemm_mpi_cyclic_permutation(cprjX0_spinor,buffer%self,cprj_work_spinor,&
      & lobpcg%xg_nonlop%me_band,blocksize=lobpcg%blockdim_cprj/nspinor,comm=lobpcg%xg_nonlop%comm_band)
    call xgBlock_saxpy(cprjvar,-1.0d0,cprj_work)

    call xg_free(buffer)

    call timab(tim_ortho,2,tsec)

  end subroutine lobpcg_orthoXwrtBlocks

  subroutine lobpcg_setX0(lobpcg,iblock)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: iblock
    type(xgBlock_t) :: Xtmp,cprjXtmp
    integer :: blockdim
    integer :: spacedim
    integer :: cprjdim
    integer :: blockdim_cprj
    double precision :: tsec(2)

    call timab(tim_copy,1,tsec)
    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    cprjdim       = lobpcg%cprjdim
    blockdim_cprj = lobpcg%blockdim_cprj

    !X0(:,(iblock-1)*blockdim+1:iblock*blockdim) = lobpcg%XWP(:,lobpcg%X+1:lobpcg%X+blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,Xtmp,spacedim,blockdim,fcol=(iblock-1)*blockdim+1)
    call xgBlock_copy(lobpcg%X,Xtmp)

    call xgBlock_setBlock(lobpcg%AllcprjX0,cprjXtmp,cprjdim,blockdim_cprj,fcol=(iblock-1)*blockdim_cprj+1)
    call xgBlock_copy(lobpcg%cprjX,cprjXtmp)
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

    ! AX
    call xg_setBlock(lobpcg%AllAX0,CXtmp,lobpcg%spacedim,lobpcg%blockdim,fcol=firstcol)
    call xgBlock_copy(lobpcg%AX,CXtmp)

    call timab(tim_copy,2,tsec)

  end subroutine lobpcg_transferAX_BX

  subroutine lobpcg_free(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

    call xg_free(lobpcg%XWP)
    call xg_free(lobpcg%AXWP)
    call xg_free(lobpcg%cprjXWP)
    call xg_free(lobpcg%cprj_work)
    call xg_free(lobpcg%cprj_work2)
    call xg_free(lobpcg%AllAX0)

  end subroutine lobpcg_free

end module m_lobpcg2_cprj
!!***

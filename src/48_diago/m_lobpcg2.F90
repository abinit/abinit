
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

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
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
!  integer, parameter :: tim_run      = 1653
  integer, parameter :: tim_getAX_BX = 1654
  integer, parameter :: tim_ortho    = 1655
!  integer, parameter :: tim_Bortho   = 1656
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
    integer :: gpu_option                    ! Which GPU version is used (0=none)
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
&      space, spacecom, paral_kgb, comm_rows, comm_cols, gpu_option)

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
    integer         , intent(in   ) :: gpu_option
    double precision :: tsec(2)
    double precision :: advice
    double precision :: advice_target
    !character(len=255) :: linalg_threads
    !integer :: ierr
    integer :: iadvice, nthread
#ifdef HAVE_LINALG_MKL_THREADS
    integer :: mkl_get_max_threads
#endif
#ifdef HAVE_LINALG_OPENBLAS_THREADS
    integer :: openblas_get_num_threads
#endif

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
    lobpcg%comm_rows  = comm_rows
    lobpcg%comm_cols  = comm_cols
    lobpcg%gpu_option = gpu_option

    nthread = 1
#ifdef HAVE_LINALG_MKL_THREADS
    nthread =  mkl_get_max_threads()
#elif HAVE_LINALG_OPENBLAS_THREADS
    nthread =  openblas_get_num_threads()
#else
!#elif defined HAVE_FC_GETENV
    !call getenv("OMP_NUM_THREADS",linalg_threads)
    nthread = xomp_get_num_threads(open_parallel=.true.)
    !read(linalg_threads,'(i5)',iostat=ierr) nthread
    !if ( ierr /= 0 ) nthread = 1
    if ( nthread == 0 ) nthread = 1
#endif

    advice_target = 2.5d6*dble(nthread)
    advice = advice_target/dble(spacedim) ! assume npband*npfft = cst and we adjust bandpp obtain the correct blocksize
    iadvice = 1+int(dble(neigenpairs)/advice) ! get the int so that advice is a divisor or neigenpairs
    do while (iadvice >= 1)
      if ( mod(neigenpairs,iadvice) == 0 ) then
        exit
      end if
      iadvice = iadvice - 1
    end do

!    if ( abs(dble(spacedim * blockdim)/advice_target-1.d0) > 0.5 ) then
!      if ( neigenpairs /= blockdim*iadvice ) then
!        write(std_out,'(1x,A,i5)') "You should try to get npband*bandpp=", neigenpairs/iadvice
!        write(std_out,'(1x,A,i8)',advance="no") "For information matrix size is ", spacedim*blockdim
!        if ( nthread > 1 ) then
!          write(std_out,'(1x,A,i3,1x,A)') "and linalg will use", nthread, "threads"
!        else
!          write(std_out,*)
!        end if
!      end if
!    end if

    call lobpcg_allocateAll(lobpcg,space)
    call timab(tim_init,2,tsec)
  end subroutine lobpcg_init


  subroutine lobpcg_allocateAll(lobpcg,space)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: space
    integer :: spacedim
    integer :: blockdim

    spacedim = lobpcg%spacedim
    blockdim = lobpcg%blockdim

    call lobpcg_free(lobpcg) ! Make sure everything is not allocated and
    ! pointer point to null()

    call xg_init(lobpcg%XWP,space,spacedim,3*blockdim,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%XWP,lobpcg%X,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%W,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%P,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%XW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%WP,blockdim+1,spacedim,2*blockdim)

    call xg_init(lobpcg%AXWP,space,spacedim,3*blockdim,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AX,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AW,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AP,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AXW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AWP,blockdim+1,spacedim,2*blockdim)

    call xg_init(lobpcg%BXWP,space,spacedim,3*blockdim,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BX,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BW,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BP,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BXW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BWP,blockdim+1,spacedim,2*blockdim)

    if ( lobpcg%nblock /= 1 ) then
      call xg_init(lobpcg%AllBX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)
      call xg_init(lobpcg%AllAX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)
    else
      lobpcg%AllBX0%self = lobpcg%BX
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


  subroutine lobpcg_run(lobpcg, X0, getAX_BX, pcond, eigen, occ, residu, prtvol, nspinor, isppol, ikpt, inonsc, istep, nbdbuf)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: X0   ! Full initial vectors
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: occ
    type(xgBlock_t), intent(inout) :: residu
    integer        , intent(in   ) :: prtvol
    integer        , intent(in   ) :: nspinor
    integer        , intent(in   ) :: isppol,ikpt,inonsc,istep,nbdbuf

    type(xg_t) :: eigenvalues3N   ! eigen values for Rayleight-Ritz
    type(xg_t) :: residu_eff
    type(xgBlock_t) :: eigenvaluesN   ! eigen values for Rayleight-Ritz
    type(xgBlock_t) :: eigenvalues2N   ! eigen values for Rayleight-Ritz
    integer :: blockdim, blockdim3, blockdim2
!    integer :: comm_fft_save,comm_band_save
    integer :: spacedim
    integer :: iblock, nblock
    integer :: iline, nline
    integer :: rows_tmp, cols_tmp, nband_eff, iband_min, iband_max
    type(xgBlock_t) :: eigenBlock   !
    type(xgBlock_t) :: residuBlock,occBlock
    double precision :: maxResidu, minResidu
    double precision :: dlamch,tolerance
    integer :: ierr = 0
    integer :: nrestart
    double precision :: tsec(2)
    logical :: compute_residu
    character(len=500) :: msg

    interface
      subroutine getAX_BX(X,AX,BX,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
        type(xgTransposer_t), intent(inout) :: transposer
      end subroutine getAX_BX
    end interface
    interface
      subroutine pcond(W,gpu_option)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
        integer, intent(in) :: gpu_option
      end subroutine pcond
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
      write(msg,'(a,es16.6)') 'lobpcg%tolerance(tolwfr_diago)=',lobpcg%tolerance
      call wrtout(std_out,msg,'COLL')
    end if

    call xg_init(eigenvalues3N,SPACE_R,blockdim3,1, gpu_option=lobpcg%gpu_option)
    call xg_setBlock(eigenvalues3N,eigenvaluesN,1,blockdim,1)
    call xg_setBlock(eigenvalues3N,eigenvalues2N,1,blockdim2,1)

    call xgBlock_reshape(eigen,(/ blockdim, nblock /))
    call xgBlock_reshape(residu,(/ blockdim, nblock /))
    call xgBlock_reshape(occ,(/ blockdim, nblock /))

    lobpcg%AllX0 = X0

    call xg_init(residu_eff,SPACE_R,blockdim,1)

    if ( lobpcg%paral_kgb == 1 ) then
      call xgTransposer_constructor(lobpcg%xgTransposerX,lobpcg%X,lobpcg%XColsRows,nspinor,&
        STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0)
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
    else
      call xgBlock_setBlock(lobpcg%X, lobpcg%XColsRows, 1, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AX, lobpcg%AXColsRows, 1, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%BX, lobpcg%BXColsRows, 1, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%W, lobpcg%WColsRows, 1, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%AW, lobpcg%AWColsRows, 1, spacedim, blockdim)
      call xgBlock_setBlock(lobpcg%BW, lobpcg%BWColsRows, 1, spacedim, blockdim)
    end if

    !! Start big loop over blocks
    do iblock = 1, nblock
      ABI_NVTX_START_RANGE(NVTX_LOBPCG2_BLOCK)
      nrestart = 0

      call lobpcg_getX0(lobpcg,iblock)
      call xgBlock_setBlock(residu,residuBlock,iblock,blockdim,1)
      call xgBlock_setBlock(occ,occBlock,iblock,blockdim,1)

      if ( iblock > 1 ) then
        call lobpcg_setPreviousX0_BX0(lobpcg,iblock)

        ! Orthogonalize current iblock X block With Respect To previous Blocks in B-basis
        call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%X,iblock)
      end if

      if (lobpcg%paral_kgb == 1) then
        call xgTransposer_transpose(lobpcg%xgTransposerX,STATE_COLSROWS)
        lobpcg%xgTransposerAX%state=STATE_COLSROWS
        lobpcg%xgTransposerBX%state=STATE_COLSROWS
      end if
      ! Initialize some quantitites (AX and BX)
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(lobpcg%XColsRows,lobpcg%AXColsRows,lobpcg%BXColsRows,lobpcg%xgTransposerX)
      call timab(tim_ax_bx,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call xgTransposer_transpose(lobpcg%xgTransposerX,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAX,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerBX,STATE_LINALG)
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

        ! Apply preconditioner
        call timab(tim_pcond,1,tsec)
        call pcond(lobpcg%W,lobpcg%gpu_option)
        call timab(tim_pcond,2,tsec)

        ! Compute residu norm here !
        call timab(tim_maxres,1,tsec)
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock,gpu_option=lobpcg%gpu_option)
        call timab(tim_maxres,2,tsec)

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
          call xgBlock_apply_diag_nospin(residuBlock,occBlock,1,Y=residu_eff%self)
          call xgBlock_minmax(residu_eff%self,minResidu,maxResidu)
        else
          ABI_ERROR('Bad value of nbdbuf')
        end if
        if ( maxResidu < lobpcg%tolerance ) then
          compute_residu = .false.
          exit
        end if

        ! Orthonormalize with respect to previous blocks
        if ( iblock > 1 ) then
          call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%W,iblock)
        end if

        if (lobpcg%paral_kgb == 1) then
          call xgTransposer_transpose(lobpcg%xgTransposerW,STATE_COLSROWS)
          lobpcg%xgTransposerAW%state=STATE_COLSROWS
          lobpcg%xgTransposerBW%state=STATE_COLSROWS
        end if
        ! Apply A and B on W
        call timab(tim_ax_bx,1,tsec)
        call getAX_BX(lobpcg%WColsRows,lobpcg%AWColsRows,lobpcg%BWColsRows,lobpcg%xgTransposerX)
        call timab(tim_ax_bx,2,tsec)
        if (lobpcg%paral_kgb == 1) then
          call xgTransposer_transpose(lobpcg%xgTransposerW,STATE_LINALG)
          call xgTransposer_transpose(lobpcg%xgTransposerAW,STATE_LINALG)
          call xgTransposer_transpose(lobpcg%xgTransposerBW,STATE_LINALG)
        end if

        ! DO RR in the correct subspace
        ! if residu starts to be too small, there is an accumulation error in
        ! P with values such as 1e-29 that make the eigenvectors diverge
        if ( iline == 1 .or. minResidu < 1e-27) then
          ! Do RR on XW to get the eigen vectors
          call xg_Borthonormalize(lobpcg%XW,lobpcg%BXW,ierr,tim_Bortho_XW,lobpcg%gpu_option,AX=lobpcg%AXW) ! Do rotate AW
          call xgBlock_zero(lobpcg%P, gpu_option=lobpcg%gpu_option)
          call xgBlock_zero(lobpcg%AP, gpu_option=lobpcg%gpu_option)
          call xgBlock_zero(lobpcg%BP, gpu_option=lobpcg%gpu_option)
          if ( ierr /= 0 ) then
            ABI_COMMENT("B-orthonormalization (XW) did not work.")
          end if
          call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvalues2N,ierr,lobpcg%prtvol,tim_RR_XW,lobpcg%gpu_option,tolerance=tolerance,&
           & XW=lobpcg%XW,AW=lobpcg%AW,BW=lobpcg%BW,P=lobpcg%P,AP=lobpcg%AP,BP=lobpcg%BP,WP=lobpcg%WP,&
           & AWP=lobpcg%AWP,BWP=lobpcg%BWP)
          if ( ierr /= 0 ) then
            ABI_WARNING("RayleighRitz (XW) did not work, but continue anyway.")
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
              exit
            end if
          else
            ABI_COMMENT("B-orthonormalization (XWP) did not work, try on XW.")
            call xg_Borthonormalize(lobpcg%XW,lobpcg%BXW,ierr,tim_Bortho_XW,lobpcg%gpu_option,AX=lobpcg%AXW) ! Do rotate AW
            if ( ierr /= 0 ) then
              ABI_COMMENT("B-orthonormalization (XW) did not work.")
            end if
            call xgBlock_zero(lobpcg%P, gpu_option=lobpcg%gpu_option)
            call xgBlock_zero(lobpcg%AP, gpu_option=lobpcg%gpu_option)
            call xgBlock_zero(lobpcg%BP, gpu_option=lobpcg%gpu_option)
            nrestart = nrestart + 1
            call xg_RayleighRitz(lobpcg%X,lobpcg%AX,lobpcg%BX,eigenvalues2N,ierr,lobpcg%prtvol,tim_RR_XW,lobpcg%gpu_option,tolerance=tolerance,&
           & XW=lobpcg%XW,AW=lobpcg%AW,BW=lobpcg%BW,P=lobpcg%P,AP=lobpcg%AP,BP=lobpcg%BP,WP=lobpcg%WP,&
           & AWP=lobpcg%AWP,BWP=lobpcg%BWP)
            if ( ierr /= 0 ) then
              ABI_WARNING("RayleighRitz (XWP) did not work, but continue anyway.")
              exit
            end if
          end if
        end if

        ABI_NVTX_END_RANGE()
      end do

      if ( compute_residu ) then
        ! Recompute AX-Lambda*BX for the last time
        call lobpcg_getResidu(lobpcg,eigenvaluesN)
        ! Apply preconditioner
        call pcond(lobpcg%W,lobpcg%gpu_option)
        ! Recompute residu norm here !
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock,gpu_option=lobpcg%gpu_option)
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
          call xgBlock_apply_diag_nospin(residuBlock,occBlock,1,Y=residu_eff%self)
          call xgBlock_minmax(residu_eff%self,minResidu,maxResidu)
        else
          ABI_ERROR('Bad value of nbdbuf')
        end if
      end if

      if (prtvol==5.and.xmpi_comm_rank(lobpcg%spacecom)==0) then
        write(msg,'(6(a,i4),2(a,es16.6))') 'lobpcg | istep=',istep,'| isppol=',isppol,'| ikpt=',ikpt,&
          & '| inonsc=',inonsc,'| iblock=',iblock,'| nline_done=',iline-1,'| minRes=',minResidu,'| maxRes=',maxResidu
        call wrtout(std_out,msg,'PERS')
      end if

      ! Save eigenvalues
      call xgBlock_setBlock(eigen,eigenBlock,iblock,blockdim,1)
      call xgBlock_copy(eigenvaluesN,eigenBlock,gpu_option=lobpcg%gpu_option)

      ! Save new X in X0
      call lobpcg_setX0(lobpcg,iblock)

      ! Copy previous BX into BX0 for previous block
      if ( nblock > 1 ) then
        call lobpcg_transferAX_BX(lobpcg,iblock)
      end if

      ABI_NVTX_END_RANGE()
    end do !! End iblock loop

    call xgBlock_reshape(eigen,(/ blockdim*nblock, 1 /))
    call xgBlock_reshape(residu,(/ blockdim*nblock, 1 /))
    call xgBlock_reshape(occ,(/ blockdim*nblock, 1 /))

    call xg_free(eigenvalues3N)
    call xg_free(residu_eff)

    if ( ierr /= 0 ) then
      ABI_COMMENT("Some errors happened, so H|Psi> and S|Psi> are computed before leaving")
      if ( lobpcg%paral_kgb == 1 ) then
        call xgTransposer_constructor(lobpcg%xgTransposerAllX0,lobpcg%AllX0,lobpcg%AllX0ColsRows,nspinor,&
          STATE_LINALG,TRANS_ALL2ALL,lobpcg%comm_rows,lobpcg%comm_cols,0,0)
        call xgTransposer_copyConstructor(lobpcg%xgTransposerAllAX0,lobpcg%xgTransposerAllX0,&
          lobpcg%AllAX0%self,lobpcg%AllAX0ColsRows,STATE_LINALG)
        call xgTransposer_copyConstructor(lobpcg%xgTransposerAllBX0,lobpcg%xgTransposerAllX0,&
          lobpcg%AllBX0%self,lobpcg%AllBX0ColsRows,STATE_LINALG)
      else
        call xgBlock_setBlock(lobpcg%AllX0      , lobpcg%AllX0ColsRows , 1, spacedim, lobpcg%neigenpairs)
        call xgBlock_setBlock(lobpcg%AllAX0%self, lobpcg%AllAX0ColsRows, 1, spacedim, lobpcg%neigenpairs)
        call xgBlock_setBlock(lobpcg%AllBX0%self, lobpcg%AllBX0ColsRows, 1, spacedim, lobpcg%neigenpairs)
      end if
      if (lobpcg%paral_kgb == 1) then
        call xgTransposer_transpose(lobpcg%xgTransposerAllX0,STATE_COLSROWS)
        lobpcg%xgTransposerAllAX0%state=STATE_COLSROWS
        lobpcg%xgTransposerAllBX0%state=STATE_COLSROWS
      end if
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(lobpcg%AllX0ColsRows,lobpcg%AllAX0ColsRows,lobpcg%AllBX0ColsRows,lobpcg%xgTransposerAllX0)
      call timab(tim_ax_bx,2,tsec)
      if (lobpcg%paral_kgb == 1) then
        call xgTransposer_transpose(lobpcg%xgTransposerAllX0,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAllAX0,STATE_LINALG)
        call xgTransposer_transpose(lobpcg%xgTransposerAllBX0,STATE_LINALG)
      end if
      call xgTransposer_free(lobpcg%xgTransposerAllX0)
      call xgTransposer_free(lobpcg%xgTransposerAllAX0)
      call xgTransposer_free(lobpcg%xgTransposerAllBX0)
      nblock = 1 ! Avoid the next RR
    end if

    if ( nblock > 1 ) then
      call xg_Borthonormalize(X0,lobpcg%AllBX0%self,ierr,tim_Bortho_Xall,&
        & lobpcg%gpu_option,AX=lobpcg%AllAX0%self) ! Do rotate AX
      call xg_RayleighRitz(X0,lobpcg%AllAX0%self,lobpcg%AllBX0%self,eigen,ierr,lobpcg%prtvol,tim_RR_Xall,&
        & lobpcg%gpu_option,tolerance=tolerance)
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

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    !lobpcg%XWP(:,X+1:X+blockdim) = lobpcg%X0(:,(iblock-1)*blockdim+1:iblock*blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,lobpcg%X0,(iblock-1)*blockdim+1,spacedim,blockdim)
    call xgBlock_copy(lobpcg%X0,lobpcg%X, gpu_option=lobpcg%gpu_option)

  end subroutine lobpcg_getX0


  subroutine lobpcg_setPreviousX0_BX0(lobpcg,iblock)

    type(lobpcg_t) , intent(inout) :: lobpcg
    integer        , intent(in   ) :: iblock
    call xg_setBlock(lobpcg%AllBX0,lobpcg%BX0,1,lobpcg%spacedim,(iblock-1)*lobpcg%blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,lobpcg%X0,1,lobpcg%spacedim,(iblock-1)*lobpcg%blockdim)
  end subroutine lobpcg_setPreviousX0_BX0


  subroutine lobpcg_orthoXwrtBlocks(lobpcg,var,iblock)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: var
    integer        , intent(in   ) :: iblock
    integer :: previousBlock
    integer :: blockdim
    integer :: spacedim
    !integer :: shift
    type(xg_t) :: buffer
    double precision :: tsec(2)

    call timab(tim_ortho,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_LOBPCG2_ORTHO_X_WRT)

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim
    previousBlock = (iblock-1)*lobpcg%blockdim

    call xg_init(buffer,space(var),previousBlock,blockdim,lobpcg%spacecom, gpu_option=lobpcg%gpu_option)

    ! buffer = BX0^T*X
    call xgBlock_gemm(lobpcg%BX0%trans,lobpcg%X%normal,1.0d0,lobpcg%BX0,var,0.d0,buffer%self,&
        gpu_option=lobpcg%gpu_option)

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm(lobpcg%X0%normal,lobpcg%X0%normal,-1.0d0,lobpcg%X0,buffer%self,1.0d0,&
        var,gpu_option=lobpcg%gpu_option)

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
    call xgBlock_colwiseCymax(lobpcg%W,eigenvalues,lobpcg%BX,lobpcg%AX, gpu_option=lobpcg%gpu_option)
    ABI_NVTX_END_RANGE()
    call timab(tim_maxres,2,tsec)
  end subroutine lobpcg_getResidu


  subroutine lobpcg_setX0(lobpcg,iblock)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: iblock
    type(xgBlock_t) :: Xtmp
    integer :: blockdim
    integer :: spacedim

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    !X0(:,(iblock-1)*blockdim+1:iblock*blockdim) = lobpcg%XWP(:,lobpcg%X+1:lobpcg%X+blockdim)
    call xgBlock_setBlock(lobpcg%AllX0,Xtmp,(iblock-1)*blockdim+1,spacedim,blockdim)
    call xgBlock_copy(lobpcg%X,Xtmp, gpu_option=lobpcg%gpu_option)
  end subroutine lobpcg_setX0


  subroutine lobpcg_transferAX_BX(lobpcg,jblock)

    type(lobpcg_t), intent(inout) :: lobpcg
    integer       , intent(in   ) :: jblock
    type(xgBlock_t) :: CXtmp
    integer :: firstcol
    ! jblock goes from 1 to nblock-1 included
    firstcol = (jblock-1)*lobpcg%blockdim+1  ! Start of each block

    ! BX
    call xg_setBlock(lobpcg%AllBX0,CXtmp,firstcol,lobpcg%spacedim,lobpcg%blockdim)
    call xgBlock_copy(lobpcg%BX,CXtmp, gpu_option=lobpcg%gpu_option)

    ! AX
    call xg_setBlock(lobpcg%AllAX0,CXtmp,firstcol,lobpcg%spacedim,lobpcg%blockdim)
    call xgBlock_copy(lobpcg%AX,CXtmp, gpu_option=lobpcg%gpu_option)
  end subroutine lobpcg_transferAX_BX


  subroutine lobpcg_allowNested(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

!#ifdef HAVE_OPENMP
!    lobpcg%is_nested = omp_get_nested()
!    call omp_set_nested(.true.)
!#else
    lobpcg%is_nested = .false.
!#endif
  end subroutine lobpcg_allowNested

  subroutine lobpcg_restoreNested(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

!#ifdef HAVE_OPENMP
!    call omp_set_nested(lobpcg%is_nested)
!#else
    lobpcg%is_nested = .false.
!#endif
  end subroutine lobpcg_restoreNested


  subroutine lobpcg_free(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

    call xg_free(lobpcg%XWP)
    call xg_free(lobpcg%AXWP)
    call xg_free(lobpcg%BXWP)
    call xg_free(lobpcg%AllAX0)
    call xg_free(lobpcg%AllBX0)
  end subroutine lobpcg_free

end module m_lobpcg2

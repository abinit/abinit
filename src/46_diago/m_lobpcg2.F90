
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_lobpcg2

  use m_xg
  use m_xgScalapack
  use defs_basis, only : std_err, std_out
  use m_abicore
  use m_errors
  use m_xomp
#ifdef HAVE_OPENMP
  use omp_lib
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

  integer, parameter :: EIGENVX = 1
  integer, parameter :: EIGENVD = 2
  integer, parameter :: EIGENV = 3
  integer, parameter :: EIGENPVX = 4
  integer, parameter :: EIGENPVD = 5
  integer, parameter :: EIGENPV = 6
  integer, parameter :: EIGENEVD = 7
  integer, parameter :: EIGENEV = 8
  integer, parameter :: EIGENPEVD = 9
  integer, parameter :: EIGENPEV = 10
  integer, parameter :: EIGENSLK = 11
  logical, parameter :: EIGPACK(11) = &
    (/ .false.,.false.,.false., &
       .true. ,.true. ,.true. ,&
       .false.,.false.,&
       .true. ,.true., .false.  /)

  integer, parameter :: tim_init     = 1651
  integer, parameter :: tim_free     = 1652
  integer, parameter :: tim_run      = 1653
  integer, parameter :: tim_getAX_BX = 1654
  integer, parameter :: tim_ortho    = 1655
  integer, parameter :: tim_Bortho   = 1656
  integer, parameter :: tim_RR       = 1657
  integer, parameter :: tim_maxres   = 1658
  integer, parameter :: tim_ax_bx    = 1659
  integer, parameter :: tim_pcond    = 1660
  integer, parameter :: tim_hegv     = 1661

#ifdef HAVE_OPENMP
  integer, save :: eigenSolver = EIGENVD     ! Type of eigen solver to use
#else
  integer, save :: eigenSolver = EIGENVX     ! Type of eigen solver to use
#endif
  !double precision, save :: eigenSolverTime(10) = tiny(0.d0)      ! Store time for each solver
  !double precision, save :: eigenSolverCount(10) = 0      ! Store time for each solver

  type, public :: lobpcg_t
    logical :: is_nested                     ! For OpenMP nested region
    integer :: spacedim                      ! Space dimension for one vector
    integer :: neigenpairs                   ! Number of eigen values/vectors we want
    integer :: nblock                        ! Number of block in the space of dim spacedim
    integer :: blockdim                      ! Number of vectors in one block
    integer :: nline                         ! Number of line to perform
    integer :: spacecom                      ! Communicator for MPI
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

    ! Variable for work for lapack
  end type lobpcg_t

  public :: lobpcg_init
  public :: lobpcg_memInfo
  public :: lobpcg_run
  public :: lobpcg_free
  public :: lobpcg_getAX_BX

  contains


  subroutine lobpcg_init(lobpcg, neigenpairs, spacedim, blockdim, tolerance, nline, space, spacecom)

    type(lobpcg_t)  , intent(inout) :: lobpcg
    integer         , intent(in   ) :: neigenpairs
    integer         , intent(in   ) :: spacedim
    integer         , intent(in   ) :: blockdim
    double precision, intent(in   ) :: tolerance
    integer         , intent(in   ) :: nline
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: spacecom
    double precision :: tsec(2)
    double precision :: advice
    double precision :: advice_target
    !character(len=255) :: linalg_threads
    !integer :: ierr
    integer :: iadvice, nthread
#ifdef HAVE_LINALG_MKL_THREADS
    integer :: mkl_get_max_threads
#endif

    call timab(tim_init,1,tsec)
    lobpcg%neigenpairs = neigenpairs
    lobpcg%spacedim    = spacedim
    lobpcg%blockdim    = blockdim
    lobpcg%tolerance   = tolerance
    lobpcg%nline       = nline
    lobpcg%spacecom    = spacecom
    lobpcg%nblock      = neigenpairs / blockdim

    nthread = 1
#ifdef HAVE_LINALG_MKL_THREADS
    nthread =  mkl_get_max_threads()
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

    if ( abs(dble(spacedim * blockdim)/advice_target-1.d0) > 0.5 ) then
      if ( neigenpairs /= blockdim*iadvice ) then
        write(std_out,'(1x,A,i5)') "You should try to get npband*bandpp=", neigenpairs/iadvice
        write(std_out,'(1x,A,i8)',advance="no") "For information matrix size is ", spacedim*blockdim
        if ( nthread > 1 ) then
          write(std_out,'(1x,A,i3,1x,A)') "and linalg will use", nthread, "threads"
        else
          write(std_out,*)
        end if
      end if
    end if

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

    call xg_init(lobpcg%XWP,space,spacedim,3*blockdim,lobpcg%spacecom)
    call xg_setBlock(lobpcg%XWP,lobpcg%X,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%W,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%P,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%XW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%XWP,lobpcg%WP,blockdim+1,spacedim,2*blockdim)

    call xg_init(lobpcg%AXWP,space,spacedim,3*blockdim,lobpcg%spacecom)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AX,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AW,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AP,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AXW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%AXWP,lobpcg%AWP,blockdim+1,spacedim,2*blockdim)

    call xg_init(lobpcg%BXWP,space,spacedim,3*blockdim,lobpcg%spacecom)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BX,1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BW,blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BP,2*blockdim+1,spacedim,blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BXW,1,spacedim,2*blockdim)
    call xg_setBlock(lobpcg%BXWP,lobpcg%BWP,blockdim+1,spacedim,2*blockdim)

    if ( lobpcg%nblock /= 1 ) then
      call xg_init(lobpcg%AllBX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom)
      call xg_init(lobpcg%AllAX0,space,spacedim,lobpcg%neigenpairs,lobpcg%spacecom)
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


  subroutine lobpcg_run(lobpcg, X0, getAX_BX, pcond, eigen, residu, prtvol)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(inout) :: X0   ! Full initial vectors
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: residu
    integer        , intent(in   ) :: prtvol

    type(xg_t) :: eigenvalues3N   ! eigen values for Rayleight-Ritz
    type(xgBlock_t) :: eigenvaluesN   ! eigen values for Rayleight-Ritz
    type(xgBlock_t) :: eigenvalues2N   ! eigen values for Rayleight-Ritz
    integer :: blockdim, blockdim3, blockdim2
    integer :: spacedim
    integer :: iblock, nblock
    integer :: iline, nline
    integer :: rows_tmp, cols_tmp
    integer :: RR_var
    type(xgBlock_t) :: eigenBlock   !
    type(xgBlock_t) :: residuBlock
    type(xgBlock_t):: RR_eig ! Will be eigenvaluesXN
    double precision :: maxResidu, minResidu, average, deviation
    double precision :: prevMaxResidu
    double precision :: dlamch
    integer :: eigResiduMax, eigResiduMin
    integer :: ierr = 0
    integer :: nrestart
    double precision :: tsec(2)

    interface
      subroutine getAX_BX(X,AX,BX)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
      end subroutine getAX_BX
    end interface
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface

    call timab(tim_run,1,tsec)

    lobpcg%prtvol = prtvol

    blockdim = lobpcg%blockdim
    blockdim2 = 2*blockdim
    blockdim3 = 3*blockdim
    spacedim = lobpcg%spacedim

    nblock = lobpcg%nblock
    nline = lobpcg%nline
    prevMaxResidu = huge(1d0)/100.d0 ! Divide by 100 to avoid 10*huge at the first iteration  which is a FPE

    call xgBlock_getSize(eigen,rows_tmp, cols_tmp)
    if ( rows_tmp /= lobpcg%neigenpairs .and. cols_tmp /= 1 ) then
      MSG_ERROR("Error eigen size")
    endif
    call xgBlock_getSize(X0,rows_tmp, cols_tmp)
    if ( rows_tmp /= lobpcg%spacedim ) then
      MSG_ERROR("Error X0 spacedim")
    endif
    if ( cols_tmp /= lobpcg%neigenpairs ) then
      MSG_ERROR("Error X0 npairs")
    endif

    call xg_init(eigenvalues3N,SPACE_R,blockdim3,1)
    call xg_setBlock(eigenvalues3N,eigenvaluesN,1,blockdim,1)
    call xg_setBlock(eigenvalues3N,eigenvalues2N,1,blockdim2,1)

    call xgBlock_reshape(eigen,(/ blockdim, nblock /))
    call xgBlock_reshape(residu,(/ blockdim, nblock /))

    lobpcg%AllX0 = X0

    !! Start big loop over blocks
    do iblock = 1, nblock
      if ( prtvol == 4 ) write(std_out,*) "  -- Block ", iblock
      nrestart = 0

      call lobpcg_getX0(lobpcg,iblock)
      call xgBlock_setBlock(residu,residuBlock,iblock,blockdim,1)

      if ( iblock > 1 ) then
        call lobpcg_setPreviousX0_BX0(lobpcg,iblock)

        ! Orthogonalize current iblock X block With Respect To previous Blocks in B-basis
        call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%X,iblock)
      end if

      ! Initialize some quantitites (AX and BX)
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(lobpcg%X,lobpcg%AX,lobpcg%BX)
      call timab(tim_ax_bx,2,tsec)

      ! B-orthonormalize X, BX and AX
      call lobpcg_Borthonormalize(lobpcg,VAR_X,.true.,ierr) ! true to rotate AX as well

      ! Do first RR on X to get the first eigen values
      call lobpcg_rayleighRitz(lobpcg,VAR_X,eigenvaluesN,ierr)

      do iline = 1, nline

        if ( ierr /= 0 ) then
          !MSG_COMMENT("Consider using more bands and nbdbuf if necessary.")
          ierr = 0
        end if

        !write(*,*) "    -> Iteration ", iline

        ! Compute AX-Lambda*BX
        call lobpcg_getResidu(lobpcg,eigenvaluesN)

        ! Apply preconditioner
        call timab(tim_pcond,1,tsec)
        call pcond(lobpcg%W)
        call timab(tim_pcond,2,tsec)

        ! Compute residu norm here !
        call timab(tim_maxres,1,tsec)
        call xgBlock_colwiseNorm2(lobpcg%W,residuBlock,max_val=maxResidu,max_elt=eigResiduMax,&
                                                       min_val=minResidu,min_elt=eigResiduMin)
        call timab(tim_maxres,2,tsec)
        if ( prtvol == 4 ) then
          write(std_out,'(2x,a1,es10.3,a1,es10.3,a,i4,a,i4,a)') &
            "(",minResidu,",",maxResidu, ") for eigen vectors (", &
            eigResiduMin,",",eigResiduMax,")"
          call xgBlock_average(residuBlock,average)
          call xgBlock_deviation(residuBlock,deviation)
          write(std_out,'(a,es21.14,a,es21.14)') "Average : ", average, " +/-", deviation
          if ( maxResidu < lobpcg%tolerance ) then
            write(std_out,*) "Block ", iblock, "converged at iline =", iline
            exit
          else if ( 10.d0*prevMaxResidu < maxResidu .and. iline > 1) then
            write(std_out,*) "Block ", iblock, "stopped at iline =", iline
            exit
          endif
        end if
        prevMaxResidu = maxResidu

        ! Orthonormalize with respect to previous blocks
        if ( iblock > 1 ) then
          call lobpcg_orthoXwrtBlocks(lobpcg,lobpcg%W,iblock)
        end if

        ! Apply A and B on W
        call timab(tim_ax_bx,1,tsec)
        call getAX_BX(lobpcg%W,lobpcg%AW,lobpcg%BW)
        call timab(tim_ax_bx,2,tsec)

        ! B-orthonormalize W, BW
        !call lobpcg_Borthonormalize(lobpcg,VAR_XW,.true.,ierr) ! Do rotate AW
        !call lobpcg_Borthonormalize(lobpcg,VAR_W,.true.,ierr) ! Do rotate AW

        ! DO RR in the correct subspace
        ! if residu starts to be too small, there is an accumulation error in
        ! P with values such as 1e-29 that make the eigenvectors diverge
        if ( iline == 1 .or. minResidu < 1e-27) then
          ! Do RR on XW to get the eigen vectors
          call lobpcg_Borthonormalize(lobpcg,VAR_XW,.true.,ierr) ! Do rotate AW
          RR_var = VAR_XW
          call xgBlock_zero(lobpcg%P)
          call xgBlock_zero(lobpcg%AP)
          call xgBlock_zero(lobpcg%BP)
          RR_eig = eigenvalues2N
          if ( ierr /= 0 ) then
            MSG_COMMENT("This is embarrasing. Let's pray")
          end if
        else
          ! B-orthonormalize P, BP
          call lobpcg_Borthonormalize(lobpcg,VAR_XWP,.true.,ierr) ! Do rotate AWP
          ! Do RR on XWP to get the eigen vectors
          if ( ierr == 0 ) then
            RR_var = VAR_XWP
            RR_eig = eigenvalues3N%self
          else
            call lobpcg_Borthonormalize(lobpcg,VAR_XW,.true.,ierr) ! Do rotate AW
            RR_var = VAR_XW
            RR_eig = eigenvalues2N
            call xgBlock_zero(lobpcg%P)
            call xgBlock_zero(lobpcg%AP)
            call xgBlock_zero(lobpcg%BP)
            nrestart = nrestart + 1
          end if
          !RR_eig = eigenvalues3N%self
        end if
        !RR_eig = eigenvalues3N%self
        call lobpcg_rayleighRitz(lobpcg,RR_var,RR_eig,ierr,2*dlamch('E'))
        if ( ierr /= 0 ) then
          MSG_ERROR_NOSTOP("I'm so so sorry I could not make it, I did my best but I failed. Sorry. I'm gonna suicide",ierr)
          exit
        end if

      end do

      ! Recompute AX-Lambda*BX for the last time
      call lobpcg_getResidu(lobpcg,eigenvaluesN)
      ! Apply preconditioner
      call pcond(lobpcg%W)
      ! Recompute residu norm here !
      call xgBlock_colwiseNorm2(lobpcg%W,residuBlock,max_val=maxResidu,max_elt=eigResiduMax, &
                                                     min_val=minResidu,min_elt=eigResiduMin)

      ! Save eigenvalues
      call xgBlock_setBlock(eigen,eigenBlock,iblock,blockdim,1)
      call xgBlock_copy(eigenvaluesN,eigenBlock)

      ! Save new X in X0
      call lobpcg_setX0(lobpcg,iblock)

      ! Copy previous BX into BX0 for previous block
      if ( nblock > 1 ) then
        call lobpcg_transferAX_BX(lobpcg,iblock)
      end if

    end do !! End iblock loop

    call xgBlock_reshape(eigen,(/ blockdim*nblock, 1 /))
    call xgBlock_reshape(residu,(/ blockdim*nblock, 1 /))

    call xg_free(eigenvalues3N)

    if ( ierr /= 0 ) then
      MSG_COMMENT("But before that, I want to recalculate H|Psi> and S|Psi>")
      call timab(tim_ax_bx,1,tsec)
      call getAX_BX(X0,lobpcg%AllAX0%self,lobpcg%AllBX0%self)
      call timab(tim_ax_bx,2,tsec)
      nblock = 1 ! Avoid the next RR
    end if

    if ( nblock > 1 ) then
      lobpcg%X = X0
      lobpcg%AX = lobpcg%AllAX0%self
      lobpcg%BX = lobpcg%AllBX0%self
      lobpcg%blockdim = blockdim*nblock
      call lobpcg_Borthonormalize(lobpcg,VAR_X,.true.,ierr) ! Do rotate AX
      call lobpcg_rayleighRitz(lobpcg,VAR_X,eigen,ierr,2*dlamch('E'))
    end if

    call timab(tim_run,2,tsec)

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
    call xgBlock_copy(lobpcg%X0,lobpcg%X)

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

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim
    previousBlock = (iblock-1)*lobpcg%blockdim

    call xg_init(buffer,space(var),previousBlock,blockdim,lobpcg%spacecom)

    ! buffer = BX0^T*X
    call xgBlock_gemm(lobpcg%BX0%trans,lobpcg%X%normal,1.0d0,lobpcg%BX0,var,0.d0,buffer%self)

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm(lobpcg%X0%normal,lobpcg%X0%normal,-1.0d0,lobpcg%X0,buffer%self,1.0d0,var)

    call xg_free(buffer)

    call timab(tim_ortho,2,tsec)

  end subroutine lobpcg_orthoXwrtBlocks


  subroutine lobpcg_Borthonormalize(lobpcg,var,BorthoA,info)

    type(lobpcg_t), intent(inout) :: lobpcg
    integer       , intent(in   ) :: var
    logical       , intent(in   ) :: BorthoA
    integer       , intent(  out) :: info
    type(xg_t) :: buffer
    type(xgBlock_t) :: X
    type(xgBlock_t) :: BX
    type(xgBlock_t) :: AX
    double precision :: tsec(2)

    call timab(tim_Bortho,1,tsec)

    select case (var)
    case (VAR_X) ! Select X vectors
      X = lobpcg%X
      AX = lobpcg%AX
      BX = lobpcg%BX
    case (VAR_W) ! Select W vectors
      X = lobpcg%W
      AX = lobpcg%AW
      BX = lobpcg%BW
    case (VAR_P) ! Select P vectors
      X = lobpcg%P
      AX = lobpcg%AP
      BX = lobpcg%BP
    case (VAR_WP) ! Select W vectors
      X = lobpcg%WP
      AX = lobpcg%AWP
      BX = lobpcg%BWP
    case (VAR_XW) ! Select W vectors
      X = lobpcg%XW
      AX = lobpcg%AXW
      BX = lobpcg%BXW
    case (VAR_XWP) ! Select W vectors
      X = lobpcg%XWP%self
      AX = lobpcg%AXWP%self
      BX = lobpcg%BXWP%self
    case default
      MSG_ERROR("Bortho")
    end select

    call xg_init(buffer,space(X),cols(X),cols(X),lobpcg%spacecom)

    ! Compute X^TBX
    call xgBlock_gemm(X%trans,BX%normal,1.d0,X,BX,0.d0,buffer%self)

    ! Compute Cholesky decomposition (Upper part)
    call xgBlock_potrf(buffer%self,'u',info)

    if ( info /= 0 ) then
      MSG_COMMENT("An old style abi_xorthonormalize happened but now I'll try to continue ;-)")
      call xg_free(buffer)
      return
    end if

!!$omp parallel default(shared)
!!$omp single
!!$omp task
    ! Solve YU=X
    call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,X)
!!$omp end task

!!$omp task
    ! Solve BYU=BX
    call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,BX)
!!$omp end task

!!$omp task
    if ( BorthoA .eqv. .true. ) then
      ! Solve AYU=AX
      call xgBlock_trsm('r','u',buffer%normal,'n',1.d0,buffer%self,AX)
    end if
!!$omp end task
!!$omp end single nowait
!!$omp end parallel

    call xg_free(buffer)

    call timab(tim_Bortho,2,tsec)

  end subroutine lobpcg_Borthonormalize


  subroutine lobpcg_rayleighRitz(lobpcg,var,eigenvalues,info,tolerance)

    use m_time
    type(lobpcg_t) , intent(inout) :: lobpcg
    integer        , intent(in   ) :: var
    type(xgBlock_t), intent(inout) :: eigenvalues
    integer        , intent(  out) :: info
    double precision, optional, intent(in) :: tolerance
    integer :: blockdim
    integer :: spacedim
    integer :: subdim
    !integer :: neigen
    double precision :: abstol
#ifdef HAVE_LINALG_SCALAPACK
    logical :: use_slk
#endif
    type(xg_t) :: vec
    type(xg_t) :: subA
    type(xg_t) :: subB
    type(xgBlock_t) :: subsub
    type(xgBlock_t) :: X
    type(xgBlock_t) :: AX
    type(xgBlock_t) :: BX
    type(xgBlock_t) :: Cwp
    type(xgBlock_t) :: WP
    type(xgBlock_t) :: AWP
    type(xgBlock_t) :: BWP
    type(xgScalapack_t) :: scalapack
    double precision :: tsec(2)
#ifdef HAVE_LINALG_MKL_THREADS
    integer :: mkl_get_max_threads
#endif

    call timab(tim_RR, 1, tsec)

    blockdim = lobpcg%blockdim
    spacedim = lobpcg%spacedim

    select case(var)
    case(VAR_X)
      subdim = blockdim
      X = lobpcg%X
      AX = lobpcg%AX
      BX = lobpcg%BX
      !eigenSolver = minloc(eigenSolverTime(7:10), dim=1) + 6
#ifdef HAVE_LINALG_MKL_THREADS
      if ( mkl_get_max_threads() > 1 ) then
        eigenSolver = EIGENEVD
      else
#endif
        eigenSolver = EIGENEV
#ifdef HAVE_LINALG_MKL_THREADS
      end if
#endif
    case(VAR_XW)
      subdim = blockdim*2
      X = lobpcg%XW
      AX = lobpcg%AXW
      BX = lobpcg%BXW
      WP = lobpcg%W
      AWP = lobpcg%AW
      BWP = lobpcg%BW
      !eigenSolver = minloc(eigenSolverTime(1:6), dim=1)
#ifdef HAVE_LINALG_MKL_THREADS
      if ( mkl_get_max_threads() > 1 ) then
        eigenSolver = EIGENVD
      else
#endif
        eigenSolver = EIGENVX
#ifdef HAVE_LINALG_MKL_THREADS
      end if
#endif
    case (VAR_XWP)
      subdim = blockdim*3
      X = lobpcg%XWP%self
      AX = lobpcg%AXWP%self
      BX = lobpcg%BXWP%self
      WP = lobpcg%WP
      AWP = lobpcg%AWP
      BWP = lobpcg%BWP
      !eigenSolver = minloc(eigenSolverTime(1:6), dim=1)

#ifdef HAVE_LINALG_MKL_THREADS
      if ( mkl_get_max_threads() > 1 ) then
        eigenSolver = EIGENVD
      else
#endif
        eigenSolver = EIGENVX
#ifdef HAVE_LINALG_MKL_THREADS
      end if
#endif
    case default
      MSG_ERROR("RR")
    end select

#ifdef HAVE_LINALG_SCALAPACK
    call xgScalapack_init(scalapack,lobpcg%spacecom,subdim,lobpcg%prtvol-2,use_slk)
    if ( use_slk) then
      eigenSolver = EIGENSLK
    end if
#endif

    ! Select diago algorithm

    abstol = 0d0 ; if ( present(tolerance) ) abstol = tolerance

    call xg_init(subA,space(X),subdim,subdim,lobpcg%spacecom)
    !call xgBlock_zero(subA%self)
    if ( var /= VAR_X ) then
      call xg_init(subB,space(X),subdim,subdim,lobpcg%spacecom)
      !call xgBlock_zero(subB%self)
      !call xgBlock_one(subB%self)
    end if

    if ( eigenSolver == EIGENVX .or. eigenSolver == EIGENPVX ) then
      call xg_init(vec,space(x),subdim,blockdim)
    else if ( EIGPACK(eigenSolver) ) then
      call xg_init(vec,space(x),subdim,subdim)
    else
      call xg_setBlock(subA,vec%self,1,subdim,blockdim)
    endif

     ! Compute subA and subB by part
    !--- begin
    ! |  E  |  XAW  | XAP |  |  I  |  XBW  | XBP |
    ! |  *  |  WAW  | WAP |  |  *  |   I   | WBP |
    ! |  *  |   *   | PAP |  |  *  |   *   |  I  |

    call xg_setBlock(subA,subsub,1,blockdim,blockdim)
    call xgBlock_gemm(X%trans,AX%normal,1.0d0,X,AX,0.d0,subsub)

    if ( var /= VAR_X ) then
      call xg_setBlock(subB,subsub,1,blockdim,blockdim)
      call xgBlock_gemm(X%trans,BX%normal,1.0d0,X,BX,0.d0,subsub)
    endif

    if ( var == VAR_XW .or. var == VAR_XWP ) then
      ! subA
      call xg_setBlock(subA,subsub,blockdim+1,2*blockdim,blockdim)
      call xgBlock_gemm(lobpcg%XW%trans,lobpcg%AW%normal,1.0d0,lobpcg%XW,lobpcg%AW,0.d0,subsub)

      ! subB
      call xg_setBlock(subB,subsub,blockdim+1,2*blockdim,blockdim)
      call xgBlock_gemm(lobpcg%XW%trans,lobpcg%BW%normal,1.0d0,lobpcg%XW,lobpcg%BW,0.d0,subsub)
    end if

    if ( var == VAR_XWP ) then
      ! subA
      call xg_setBlock(subA,subsub,2*blockdim+1,3*blockdim,blockdim)
      call xgBlock_gemm(lobpcg%XWP%trans,lobpcg%AP%normal,1.0d0,lobpcg%XWP%self,lobpcg%AP,0.d0,subsub)

      ! subB
      call xg_setBlock(subB,subsub,2*blockdim+1,3*blockdim,blockdim)
      call xgBlock_gemm(lobpcg%XWP%trans,lobpcg%BP%normal,1.0d0,lobpcg%XWP%self,lobpcg%BP,0.d0,subsub)
    end if

    if ( EIGPACK(eigenSolver) ) then
      call xgBlock_pack(subA%self,subA%self,'u')
      if ( var /= VAR_X ) then
        call xgBlock_pack(subB%self,subB%self,'u')
      end if
    end if

    ! Compute X*AX subspace matrix
    !call xgBlock_gemm(X%trans,AX%normal,1.0d0,X,AX,0.d0,subA%self)

    ! Compute X*BX subspace matrix
    !call xgBlock_gemm(X%trans,BX%normal,1.0d0,X,BX,0.d0,subB%self)
    !---end

    call timab(tim_hegv,1,tsec)
    tsec(2) = abi_wtime()
    if ( var == VAR_X ) then
    ! Solve Hermitian eigen problem
      select case (eigenSolver)
      case (EIGENEVD)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heevd"
        call xgBlock_heevd('v','u',subA%self,eigenvalues,info)
      case (EIGENEV)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using heev"
        call xgBlock_heev('v','u',subA%self,eigenvalues,info)
      case (EIGENPEVD)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpevd"
        call xgBlock_hpevd('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENPEV)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpev"
        call xgBlock_hpev('v','u',subA%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using pheev"
        call xgScalapack_heev(scalapack,subA%self,eigenvalues)
        info = 0 ! No error code returned for the moment
      case default
        MSG_ERROR("Error for Eigen Solver HEEV")
      end select
    else
      ! Solve Hermitian general eigen problem only for first blockdim eigenvalues
      select case (eigenSolver)
      case (EIGENVX)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvx"
        call xgBlock_hegvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENVD)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegvd"
        call xgBlock_hegvd(1,'v','u',subA%self,subB%self,eigenvalues,info)
      case (EIGENV)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hegv"
        call xgBlock_hegv(1,'v','u',subA%self,subB%self,eigenvalues,info)
      case (EIGENPVX)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvx"
        call xgBlock_hpgvx(1,'v','i','u',subA%self,subB%self,0.d0,0.d0,1,blockdim,abstol,&
          eigenvalues,vec%self,info)
      case (EIGENPVD)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgvd"
        call xgBlock_hpgvd(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENPV)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using hpgv"
        call xgBlock_hpgv(1,'v','u',subA%self,subB%self,eigenvalues,vec%self,info)
      case (EIGENSLK)
        if ( lobpcg%prtvol == 4 ) write(std_out,'(A,1x)',advance="no") "Using phegv"
        call xgScalapack_hegv(scalapack,subA%self,subB%self,eigenvalues)
        info = 0 ! No error code returned for the moment
      case default
        MSG_ERROR("Error for Eigen Solver HEGV")
      end select
    end if
    if ( eigenSolver == EIGENSLK ) then
      call xgScalapack_free(scalapack)
    end if
    tsec(2) = abi_wtime() - tsec(2)
!    if ( var /= VAR_XW ) then
!      eigenSolverTime(eigenSolver) = (eigenSolverTime(eigenSolver)*eigenSolverCount(eigenSolver) + tsec(2))/(eigenSolverCount(eigenSolver)+1)
!      eigenSolverCount(eigenSolver) = eigenSolverCount(eigenSolver)+1
!    end if
    if ( lobpcg%prtvol == 4 ) write(std_out,*) tsec(2)
    call timab(tim_hegv,2,tsec)

    if ( eigenSolver == EIGENVX .or. EIGPACK(eigenSolver)) then
      call xg_free(subA)
    end if
    call xg_free(subB)

    if ( info == 0 ) then
      call xg_init(subB,space(X),spacedim,blockdim)

      !/* Easy basic solution */
      !/* Compute first part of X here */
      ! Use subB as buffer
      !lobpcg%XWP (:,X+1:X+blockdim) = matmul(lobpcg%XWP (:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_setBlock(vec%self,Cwp,1,blockdim,blockdim)
      call xgBlock_gemm(lobpcg%X%normal,Cwp%normal,1.0d0,lobpcg%X,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,lobpcg%X)

      !lobpcg%AXWP(:,X+1:X+blockdim) = matmul(lobpcg%AXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm(lobpcg%AX%normal,Cwp%normal,1.0d0,lobpcg%AX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,lobpcg%AX)

      !lobpcg%BXWP(:,X+1:X+blockdim) = matmul(lobpcg%BXWP(:,X+1:X+blockdim),vec(1:blockdim,1:blockdim))
      call xgBlock_gemm(lobpcg%BX%normal,Cwp%normal,1.0d0,lobpcg%BX,Cwp,0.d0,subB%self)
      call xgBlock_copy(subB%self,lobpcg%BX)

      if ( var /= VAR_X ) then
        ! Cost to pay to avoid temporary array in xgemm
        call xgBlock_cshift(vec%self,blockdim,1) ! Bottom 2*blockdim lines are now at the top
        call xgBlock_setBlock(vec%self,Cwp,1,subdim-blockdim,blockdim)

        !lobpcg%XWP (:,P+1:P+blockdim) = matmul(lobpcg%XWP (:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(WP%normal,Cwp%normal,1.0d0,WP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,lobpcg%P)

        !lobpcg%AXWP(:,P+1:P+blockdim) = matmul(lobpcg%AXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(AWP%normal,Cwp%normal,1.0d0,AWP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,lobpcg%AP)

        !lobpcg%BXWP(:,P+1:P+blockdim) = matmul(lobpcg%BXWP(:,W+1:W+subdim-blockdim),vec(1:subdim-blockdim,1:blockdim))
        call xgBlock_gemm(BWP%normal,Cwp%normal,1.0d0,BWP,Cwp,0.d0,subB%self)
        call xgBlock_copy(subB%self,lobpcg%BP)

        !/* Maybe faster solution
        ! * Sum previous contribution plus P direction
        ! */
        call xgBlock_add(lobpcg%X,lobpcg%P)
        call xgBlock_add(lobpcg%AX,lobpcg%AP)
        call xgBlock_add(lobpcg%BX,lobpcg%BP)
      end if
    end if

    ! Doing free on an already free object does not doe anything
    call xg_free(vec)
    call xg_free(subA)
    call xg_free(subB)

    call timab(tim_RR, 2, tsec)

  end subroutine lobpcg_rayleighRitz


  subroutine lobpcg_getResidu(lobpcg,eigenvalues)

    type(lobpcg_t) , intent(inout) :: lobpcg
    type(xgBlock_t), intent(in   ) :: eigenvalues
    double precision :: tsec(2)

    call timab(tim_maxres,1,tsec)
      !lobpcg%XWP(1:spacedim,shiftW+iblock) = lobpcg%AXWP(:,shiftX+iblock) - lobpcg%BXWP(:,shiftX+iblock)*eigenvalues(iblock)
    call xgBlock_colwiseCymax(lobpcg%W,eigenvalues,lobpcg%BX,lobpcg%AX)
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
    call xgBlock_copy(lobpcg%X,Xtmp)
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
    call xgBlock_copy(lobpcg%BX,CXtmp)

    ! AX
    call xg_setBlock(lobpcg%AllAX0,CXtmp,firstcol,lobpcg%spacedim,lobpcg%blockdim)
    call xgBlock_copy(lobpcg%AX,CXtmp)
  end subroutine lobpcg_transferAX_BX


  subroutine lobpcg_getAX_BX(lobpcg,AX,BX)

    type(lobpcg_t) , intent(in   ) :: lobpcg
    type(xgBlock_t), intent(  out) :: AX
    type(xgBlock_t), intent(  out) :: BX

    AX = lobpcg%AllAX0%self
    BX = lobpcg%AllBX0%self
  end subroutine lobpcg_getAX_BX


  subroutine lobpcg_allowNested(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

#ifdef HAVE_OPENMP
    lobpcg%is_nested = omp_get_nested()
    call omp_set_nested(.true.)
#else
    lobpcg%is_nested = .false.
#endif
  end subroutine lobpcg_allowNested

  subroutine lobpcg_restoreNested(lobpcg)

    type(lobpcg_t), intent(inout) :: lobpcg

#ifdef HAVE_OPENMP
    call omp_set_nested(lobpcg%is_nested)
#else
    lobpcg%is_nested = .false.
#endif
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

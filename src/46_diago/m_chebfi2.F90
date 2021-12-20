
#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_chebfi2

  use m_xg
  use m_xgTransposer
  use m_xgScalapack
  use defs_basis, only : std_err, std_out
  use defs_abitypes
  use m_abicore
  use m_errors
  use m_xomp
  use m_xmpi
  use m_cgtools
#ifdef HAVE_OPENMP
  use omp_lib
  !use mkl_service
#endif

  use m_time, only : timab

  implicit none
  
  private
  
  integer, parameter :: EIGENV = 1
  integer, parameter :: EIGENVD = 2
  integer, parameter :: DEBUG_ROWS = 5
  integer, parameter :: DEBUG_COLUMNS = 5
  
#ifdef HAVE_OPENMP 
  integer, save :: eigenSolver = EIGENVD     ! Type of eigen solver to use
#else
  integer, save :: eigenSolver = EIGENV     ! Type of eigen solver to use
#endif
  
  integer, parameter :: tim_init         = 1751
  integer, parameter :: tim_free         = 1752
  integer, parameter :: tim_run          = 1753
  integer, parameter :: tim_getAX_BX     = 1754
  integer, parameter :: tim_invovl       = 1755
  integer, parameter :: tim_residu       = 1756
  integer, parameter :: tim_RR           = 1757
  integer, parameter :: tim_pcond        = 1758
  integer, parameter :: tim_RR_q         = 1759
  integer, parameter :: tim_next_p       = 1760
  integer, parameter :: tim_swap         = 1761
  integer, parameter :: tim_amp_f        = 1762
  integer, parameter :: tim_alltoall     = 1763
  integer, parameter :: tim_RR_hegv      = 1764
  integer, parameter :: tim_RR_scale     = 1765
  integer, parameter :: tim_RR_XNP_reset = 1766
  integer, parameter :: tim_RR_gemm_1    = 1767
  integer, parameter :: tim_RR_gemm_2    = 1768
  integer, parameter :: tim_X_NP_init    = 1769
  integer, parameter :: tim_AX_BX_init   = 1770

  type, public :: chebfi_t
    integer :: space
    integer :: spacedim                      ! Space dimension for one vector
    integer :: total_spacedim                ! Maybe not needed
    integer :: neigenpairs                   ! Number of eigen values/vectors we want
    integer :: nline                         ! Number of line to perform
    integer :: spacecom                      ! Communicator for MPI
    double precision :: tolerance            ! Tolerance on the residu to stop the minimization
    double precision :: ecut                 ! Ecut for Chebfi oracle
    
    integer :: paral_kgb                     ! MPI parallelization variables
    integer :: nproc_band
    integer :: bandpp
    integer :: nproc_fft
    
    logical :: paw
    integer :: eigenProblem   !1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
    integer :: istwf_k
    integer :: me_g0

    !ARRAYS
    type(xgBlock_t) :: X 
   
    type(xg_t) :: X_NP
    type(xgBlock_t) :: X_next
    type(xgBlock_t) :: X_prev
        
    type(xg_t) :: AX  
    type(xg_t) :: BX  
        
    type(xgBlock_t) :: xXColsRows
    type(xgBlock_t) :: xAXColsRows
    type(xgBlock_t) :: xBXColsRows
    
    type(xgTransposer_t) :: xgTransposerX
    type(xgTransposer_t) :: xgTransposerAX
    type(xgTransposer_t) :: xgTransposerBX

    type(xgBlock_t) :: eigenvalues
    
    !SWAP POINTERS
    type(xgBlock_t) :: X_swap
    type(xgBlock_t) :: AX_swap   
    type(xgBlock_t) :: BX_swap     

  end type chebfi_t

  public :: chebfi_init
  public :: chebfi_memInfo
  public :: chebfi_run
  public :: chebfi_free

  contains

  subroutine chebfi_init(chebfi, neigenpairs, spacedim, tolerance, ecut, paral_kgb, nproc_band, bandpp, nproc_fft, &
                         nline, space, eigenProblem, istwf_k, spacecom, me_g0, paw)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_init'
!End of the abilint section

    type(chebfi_t)  , intent(inout) :: chebfi
    integer         , intent(in   ) :: neigenpairs
    integer         , intent(in   ) :: spacedim
    double precision, intent(in   ) :: tolerance
    double precision, intent(in   ) :: ecut
    integer         , intent(in   ) :: paral_kgb
    integer         , intent(in   ) :: nproc_band
    integer         , intent(in   ) :: bandpp
    integer         , intent(in   ) :: nproc_fft
    integer         , intent(in   ) :: nline
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: eigenProblem
    integer         , intent(in   ) :: istwf_k
    integer         , intent(in   ) :: spacecom
    integer         , intent(in   ) :: me_g0
    logical         , intent(in   ) :: paw
    double precision :: tsec(2)

    call timab(tim_init,1,tsec)
    
    chebfi%space = space
    chebfi%neigenpairs = neigenpairs  
    chebfi%spacedim    = spacedim
    chebfi%tolerance   = tolerance
    chebfi%ecut        = ecut
    chebfi%paral_kgb   = paral_kgb
    chebfi%nproc_band  = nproc_band
    chebfi%bandpp      = bandpp
    chebfi%nproc_fft   = nproc_fft
    chebfi%nline       = nline
    chebfi%spacecom    = spacecom
    chebfi%eigenProblem = eigenProblem
    chebfi%istwf_k = istwf_k
    chebfi%me_g0 = me_g0
    chebfi%paw = paw
    
    call chebfi_allocateAll(chebfi)
    
    call timab(tim_init,2,tsec)
    
  end subroutine chebfi_init

  subroutine chebfi_allocateAll(chebfi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_allocateAll'
!End of the abilint section

    type(chebfi_t)  , intent(inout) :: chebfi
    integer :: spacedim
    integer :: neigenpairs
    integer :: space

    integer :: row, col
    integer :: total_spacedim, ierr 
    double precision :: tsec(2)
    
    space = chebfi%space
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs

    call chebfi_free(chebfi) 
       
    call timab(tim_X_NP_init,1,tsec)
    if (chebfi%paral_kgb == 0) then
      chebfi%total_spacedim = spacedim
      call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs, chebfi%spacecom) !regular arrays
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs)  
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)  
    else     
      total_spacedim = spacedim
      call xmpi_sum(total_spacedim,chebfi%spacecom,ierr)
      chebfi%total_spacedim = total_spacedim
      call xg_init(chebfi%X_NP,space,total_spacedim,2*chebfi%bandpp,chebfi%spacecom) !transposed arrays
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,total_spacedim,chebfi%bandpp)  
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,chebfi%bandpp+1,total_spacedim,chebfi%bandpp)  
    end if
    call timab(tim_X_NP_init,2,tsec) 

    call timab(tim_AX_BX_init,1,tsec)
    !transposer will handle these arrays automatically
    call xg_init(chebfi%AX,space,spacedim,neigenpairs,chebfi%spacecom)
    call xg_init(chebfi%BX,space,spacedim,neigenpairs,chebfi%spacecom)
    call timab(tim_AX_BX_init,2,tsec)

  end subroutine chebfi_allocateAll

  function chebfi_memInfo(neigenpairs, spacedim, space, paral_kgb, total_spacedim, bandpp) result(arraymem)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_memInfo'
!End of the abilint section

    integer, intent(in   ) :: neigenpairs
    integer, intent(in   ) :: spacedim
    integer, intent(in   ) :: space
    integer, intent(in   ) :: paral_kgb
    integer, intent(in   ) :: total_spacedim
    integer, intent(in   ) :: bandpp

    double precision :: memX_Re
    double precision :: memX_Im
    double precision :: memX

    double precision :: memX_next
    double precision :: memX_prev

    double precision :: memAX
    double precision :: memBX

    !Transposer variables 
    double precision :: memX_CR
    double precision :: memAX_CR
    double precision :: memBX_CR
    
    !chebfi_rayleighRitz function variables
    double precision :: memA_und_X
    double precision :: memB_und_X
    double precision :: memEigenvalues
    double precision :: cplx

    double precision :: arraymem(2)

    cplx = 1 
    if ( space == SPACE_C ) cplx = 2 !for now only complex

    ! Permanent in chebfi
    memX = cplx * kind(1.d0) * spacedim * neigenpairs 

    if (paral_kgb == 0) then
      memX_next = cplx * kind(1.d0) * spacedim * neigenpairs 
      memX_prev = cplx * kind(1.d0) * spacedim * neigenpairs 
    else 
      memX_next = cplx * kind(1.d0) * total_spacedim * bandpp 
      memX_prev = cplx * kind(1.d0) * total_spacedim * bandpp 
    end if

    memAX = cplx * kind(1.d0) * spacedim * neigenpairs
    memBX = cplx * kind(1.d0) * spacedim * neigenpairs

    !Transposer colrow array 
    if (paral_kgb == 1) then
      memX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
      memAX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
      memBX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
    else 
      memX_CR = 0
      memAX_CR = 0
      memBX_CR = 0
    end if    
    
    !chebfi_rayleighRitz function variables
    memA_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
    memB_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
    memEigenvalues = kind(1.d0) * neigenpairs

    arraymem(1) = memX + memX_next + memX_prev + &
                  memAX + memBX + memX_CR + memAX_CR + memBX_CR
    arraymem(2) = memA_und_X + memB_und_X + memEigenvalues

  end function chebfi_memInfo

  subroutine chebfi_run(chebfi, X0, getAX_BX, getBm1X, pcond, eigen, residu,  mpi_enreg)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_run'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0   ! Full initial vectors
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: residu
    type(mpi_type),  intent(inout) :: mpi_enreg
          
    integer :: spacedim
    integer :: neigenpairs
    integer :: nline, nline_max, nline_decrease, nline_tolwfr 
    integer :: shift
    double precision :: tolerance
    double precision :: maxeig, maxeig_global
    double precision :: mineig, mineig_global
    
    type(xg_t)::DivResults

    integer :: iline, iband    
    integer :: nCpuCols, nCpuRows
        
    !Pointers similar to old Chebfi    
    integer, allocatable :: nline_bands(:)      !Oracle variable
    double precision, pointer :: eig(:,:)

    double precision :: lambda_minus
    double precision :: lambda_plus
    double precision :: maximum
    double precision :: one_over_r
    double precision :: two_over_r
    double precision :: center
    double precision :: radius
    
    !FFT and BAND MPI communicators from rest of ABinit, to be saved
    integer :: comm_fft_save
    integer :: comm_band_save
    integer :: ierr
    
    character(len=15) :: str
   
    double precision :: tsec(2)
                    
    interface
      subroutine getAX_BX(X,AX,BX,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
        type(xgTransposer_t), optional, intent(inout) :: transposer
      end subroutine getAX_BX
    end interface
    interface
      subroutine getBm1X(X,Bm1X,iline_t,xXColsRows,X1,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
        type(integer), intent(inout) :: iline_t
        type(xgBlock_t), intent(inout) :: xXColsRows
        type(xgBlock_t), intent(inout) :: X1
        type(xgTransposer_t), optional, intent(inout) :: transposer
      end subroutine getBm1X
    end interface 
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface
    
    call timab(tim_run,1,tsec)

    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs
    nline = chebfi%nline
    chebfi%eigenvalues = eigen 
    
    if (chebfi%paral_kgb == 0) then
      allocate(nline_bands(neigenpairs))
      call xg_init(DivResults, chebfi%space, neigenpairs, 1)
    else
      allocate(nline_bands(chebfi%bandpp))
      call xg_init(DivResults, chebfi%space, chebfi%bandpp, 1)
    end if
        
    tolerance = chebfi%tolerance
    lambda_plus = chebfi%ecut
    chebfi%X = X0
       
    ! Transpose
    if (chebfi%paral_kgb == 1) then      
      nCpuRows = chebfi%nproc_fft
      nCpuCols = chebfi%nproc_band
      
      call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nCpuRows,nCpuCols,STATE_LINALG,TRANS_ALL2ALL)
     
      !save existing ABinit communicators
      comm_fft_save = mpi_enreg%comm_fft
      comm_band_save = mpi_enreg%comm_band
      
      !set new communicators from Transposer so it can interact with getghc
      !transpose correctly
      mpi_enreg%comm_fft = xgTransposer_getComm(chebfi%xgTransposerX, 2)
      mpi_enreg%comm_band = xgTransposer_getComm(chebfi%xgTransposerX, 3)
      
      call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,chebfi%AX%self,chebfi%xAXColsRows,STATE_LINALG)
      call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,chebfi%BX%self,chebfi%xBXColsRows,STATE_LINALG)
           
      call xgTransposer_transpose(chebfi%xgTransposerX,STATE_COLSROWS)
      call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_COLSROWS) 
      call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_COLSROWS)    
    else
      call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, 1, spacedim, neigenpairs)   !use xXColsRows instead of X notion
      call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, 1, spacedim, neigenpairs)   !use xAXColsRows instead of AX notion
      call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, 1, spacedim, neigenpairs)
    end if
       
    call timab(tim_getAX_BX,1,tsec)         
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows,chebfi%xgTransposerX)
    call timab(tim_getAX_BX,2,tsec)
       
    if (chebfi%paral_kgb == 1) then   
      call xmpi_barrier(chebfi%spacecom)
    end if
               
    !********************* Compute Rayleigh quotients for every band, and set λ − equal to the largest one *****!
    call timab(tim_RR_q, 1, tsec)
    call chebfi_rayleighRitzQuotiens(chebfi, maxeig, mineig, DivResults%self) !OK
    call timab(tim_RR_q, 2, tsec)
     
    if (chebfi%paral_kgb == 1) then
      call xmpi_max(maxeig,maxeig_global,chebfi%spacecom,ierr)
      call xmpi_min(mineig,mineig_global,chebfi%spacecom,ierr)
    else
      maxeig_global = maxeig
      mineig_global = mineig
    end if
                
    lambda_minus = maxeig_global
                                
    nline_max = cheb_oracle1(mineig_global, lambda_minus, lambda_plus, 1D-16, 40)
     
    if (chebfi%paral_kgb == 0) then
      call xgBlock_reverseMap(DivResults%self,eig,1,neigenpairs)
      do iband=1, neigenpairs !TODO TODO
  !      ! nline necessary to converge to tolerance
  !      !nline_tolwfr = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, tolerance / resids_filter(iband), nline)
  !      ! nline necessary to decrease residual by a constant factor
  !      !nline_decrease = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, 0.1D, dtset%nline)
  !      !nline_bands(iband) = MAX(MIN(nline_tolwfr, nline_decrease, nline_max, chebfi%nline), 1)
        nline_bands(iband) = nline ! fiddle with this to use locking
      end do
    else 
      call xgBlock_reverseMap(DivResults%self,eig,1,chebfi%bandpp)
      do iband=1, chebfi%bandpp !TODO TODO
        nline_bands(iband) = nline ! fiddle with this to use locking
      end do
    end if

    center = (lambda_plus + lambda_minus)*0.5  
    radius = (lambda_plus - lambda_minus)*0.5
                
    one_over_r = 1/radius
    two_over_r = 2/radius
       
    do iline = 0, nline - 1 
       
      call timab(tim_next_p,1,tsec)         
      call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X) 
      call timab(tim_next_p,2,tsec)  
          
      call timab(tim_swap,1,tsec)  
      if (chebfi%paral_kgb == 0) then             
        call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
      else 
        call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
      end if
      call timab(tim_swap,2,tsec) 
                
      !A * ψ  
      call timab(tim_getAX_BX,1,tsec)              
      call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows,chebfi%xgTransposerX)  
      call timab(tim_getAX_BX,2,tsec)  
    end do
    
    if (chebfi%paral_kgb == 1) then  
      call xmpi_barrier(chebfi%spacecom)
    end if 
   
    call timab(tim_amp_f,1,tsec)
    call chebfi_ampfactor(chebfi, eig, lambda_minus, lambda_plus, nline_bands)
    call timab(tim_amp_f,2,tsec)
    
    if (chebfi%paral_kgb == 1) then
      call xmpi_barrier(chebfi%spacecom)
     
      call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) 
      call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) 
      call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG) 
      
      if (xmpi_comm_size(chebfi%spacecom) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
        call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, 1, spacedim, neigenpairs)  
        call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, 1, spacedim, neigenpairs)  
     end if  
    else 
      call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
      call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, 1, spacedim, neigenpairs)  
      call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, 1, spacedim, neigenpairs)  
    end if
        
    call timab(tim_RR, 1, tsec)
    call chebfi_rayleighRitz(chebfi, nline)
    call timab(tim_RR, 2, tsec) 

    call timab(tim_residu, 1, tsec)
    maximum =  chebfi_computeResidue(chebfi, residu, pcond)
    call timab(tim_residu, 2, tsec)
                    
    deallocate(nline_bands)
    
    if (xmpi_comm_size(chebfi%spacecom) > 1) then
      call xgBlock_copy(chebfi%X_swap,chebfi%X, 1, 1)    !copy cannot be avoided :(
    end if
    
    call xg_free(DivResults)
       
    if (chebfi%paral_kgb == 1) then
      call xgTransposer_free(chebfi%xgTransposerX)
      call xgTransposer_free(chebfi%xgTransposerAX)
      call xgTransposer_free(chebfi%xgTransposerBX)
      
      !Reset communicators to original ABinit values for rest of ABinit
      mpi_enreg%comm_fft = comm_fft_save
      mpi_enreg%comm_band = comm_band_save
    end if
    
    call timab(tim_run,2,tsec)
        
  end subroutine chebfi_run

  subroutine chebfi_free(chebfi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_free'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    
    call xg_free(chebfi%X_NP)
    
    call xg_free(chebfi%AX)
    call xg_free(chebfi%BX)
    
    !call xg_finalize()

  end subroutine chebfi_free

  subroutine chebfi_rayleighRitzQuotiens(chebfi, maxeig, mineig, DivResults)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_rayleighRitzQuotiens'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    double precision, intent(inout) :: maxeig
    double precision, intent(inout) :: mineig
    type(xgBlock_t) , intent(inout) :: DivResults

    integer :: maxeig_pos(2)
    integer :: mineig_pos(2)
    integer :: info

    type(xg_t)::Results1
    type(xg_t)::Results2

    if (chebfi%paral_kgb == 0) then
      call xg_init(Results1, chebfi%space, chebfi%neigenpairs, 1)
      call xg_init(Results2, chebfi%space, chebfi%neigenpairs, 1)
    else
      call xg_init(Results1, chebfi%space, chebfi%bandpp, 1)
      call xg_init(Results2, chebfi%space, chebfi%bandpp, 1)
    end if
    
    call xgBlock_colwiseDotProduct(chebfi%xXColsRows,chebfi%xAXColsRows,Results1%self)
        
    !PAW
    call xgBlock_colwiseDotProduct(chebfi%xXColsRows,chebfi%xBXColsRows,Results2%self)

    !PAW
    call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, maxeig, maxeig_pos, mineig, mineig_pos)
         
    call xg_free(Results1)
    call xg_free(Results2)
    
  end subroutine chebfi_rayleighRitzQuotiens

  subroutine chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_computeNextOrderChebfiPolynom'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(integer) , intent(in) :: iline
    type(double precision ) , intent(in) :: center
    type(double precision ) , intent(in) :: one_over_r
    type(double precision ) , intent(in) :: two_over_r
    
    double precision :: tsec(2)   
    character(15) :: str
    
    integer :: iline_t
             
    interface
      subroutine getBm1X(X,Bm1X,iline_t,xXColsRows,X1,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
        type(integer), intent(inout) :: iline_t
        type(xgBlock_t), intent(inout) :: xXColsRows
        type(xgBlock_t), intent(inout) :: X1
        type(xgTransposer_t), optional, intent(inout) :: transposer    
      end subroutine getBm1X
    end interface

    if (chebfi%paw) then
      call timab(tim_invovl, 1, tsec)
          
      iline_t = iline
      call getBm1X(chebfi%xAXColsRows, chebfi%X_next, iline_t, chebfi%xXColsRows, chebfi%X, chebfi%xgTransposerX)
      call timab(tim_invovl, 2, tsec)
    else
      call xgBlock_copy(chebfi%xAXColsRows,chebfi%X_next, 1, 1) 
    end if  
    
    call xgBlock_scale(chebfi%xXColsRows, center, 1) !scale by c
         
    !(B-1 * A * ψi-1 - c * ψi-1)
    call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%xXColsRows)  
         
    !ψi-1  = 1/c * ψi-1 
    call xgBlock_scale(chebfi%xXColsRows, 1/center, 1) !counter scale by c
   
    if (iline == 0) then  
      call xgBlock_scale(chebfi%X_next, one_over_r, 1) 
   else
      call xgBlock_scale(chebfi%X_next, two_over_r, 1)      
      
      call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X_prev)
   end if
    
    !print *, "ZAVRSIO POLY"
    !stop

  end subroutine chebfi_computeNextOrderChebfiPolynom

  subroutine chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_swapInnerBuffers'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(integer) , intent(in) :: spacedim
    type(integer) , intent(in) :: neigenpairs
       
    call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !X_swap = X_prev     
    call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = xXColsRows
    call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, spacedim, neigenpairs) !xXColsRows = X_next
    call xgBlock_setBlock(chebfi%X_swap, chebfi%X_next, 1, spacedim, neigenpairs) !X_next = X_swap
    
  end subroutine chebfi_swapInnerBuffers
    
  subroutine chebfi_rayleighRitz(chebfi, nline)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_rayleighRitz'
!End of the abilint section
    
    use m_time,  only : timab
    
    type(chebfi_t) , intent(inout) :: chebfi
    type(integer) , intent(in) :: nline
    
    type(xg_t) :: A_und_X !H_UND_PSI
    type(xg_t) :: B_und_X !S_UND_PSI
        
    type(xgBlock_t) :: X_first_row
    type(xgBlock_t) :: AX_first_row
    type(xgBlock_t) :: BX_first_row
                
    double precision, pointer :: evec(:,:)
    double precision, pointer :: sub_pointer(:,:)
    double precision, allocatable :: edummy(:,:)
    
    integer :: neigenpairs
    integer :: spacedim
    integer :: i
    integer :: info
    integer :: space
    integer :: eigenProblem
    integer :: me_g0
    integer :: remainder
    integer :: ierr
    double precision :: tsec(2)
    
    integer :: multiply, size_remainder
    integer :: X_NP_TEST = 1
    character(len=15) :: str

    
    space = chebfi%space
    eigenProblem = chebfi%eigenProblem
    me_g0 = chebfi%me_g0
    
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs   !remains whole nband domain since it is after transpose
    
    call xg_init(A_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)
    call xg_init(B_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)
    
    call timab(tim_RR_gemm_1,1,tsec)
    !print *, "ZGEMM(T,N,M,N,K)", "M=", spacedim, "N=", neigenpairs, "K=", neigenpairs
    call xgBlock_gemm(chebfi%AX%self%trans, chebfi%X%normal, 1.0d0, chebfi%AX%self, chebfi%X, 0.d0, A_und_X%self) 
    
    !print *, "ZGEMM(T,N,M,N,K)", "M=", spacedim, "N=", neigenpairs, "K=", neigenpairs      
    if (chebfi%paw) then
      call xgBlock_gemm(chebfi%BX%self%trans, chebfi%X%normal, 1.0d0, chebfi%BX%self, chebfi%X, 0.d0, B_und_X%self) 
    else
      call xgBlock_gemm(chebfi%X%trans, chebfi%X%normal, 1.0d0, chebfi%X, chebfi%X, 0.d0, B_und_X%self) 
    end if
    call timab(tim_RR_gemm_1,2,tsec)
          
    call timab(tim_RR_scale,1,tsec)
    if(chebfi%istwf_k == 2) then    
       call xgBlock_scale(chebfi%X, 1/sqrt2, 1)   
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%X, X_first_row, 1, 2, neigenpairs) !has to be 2 rows in SPACE_CR
         call xgBlock_scale(X_first_row, sqrt2, 1)
       end if              
       call xgBlock_scale(chebfi%AX%self, 1/sqrt2, 1)
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%AX%self, AX_first_row, 1, 2, neigenpairs)
         call xgBlock_scale(AX_first_row, sqrt2, 1)
       end if  
       if (chebfi%paw) then
         call xgBlock_scale(chebfi%BX%self, 1/sqrt2, 1)
         if (me_g0 == 1)  then 
           call xgBlock_setBlock(chebfi%BX%self, BX_first_row, 1, 2, neigenpairs)
           call xgBlock_scale(BX_first_row, sqrt2, 1)
         end if
       end if 
    end if
    call timab(tim_RR_scale,2,tsec)
    
    call timab(tim_RR_hegv,1,tsec)
    !tsec(2) = abi_wtime()
    select case (eigenSolver)
    case (EIGENVD)
      call xgBlock_hegvd(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case (EIGENV)
      call xgBlock_hegv(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case default
       MSG_ERROR("Error for Eigen Solver HEGV")
    end select    
    !tsec(2) = abi_wtime() - tsec(2)
    call timab(tim_RR_hegv,2,tsec)
               
    remainder = mod(nline, 3) !3 buffer swap, keep the info which one contains X_data at the end of loop
        
    call timab(tim_RR_XNP_reset,1,tsec)
    !resize X_NP from colrwos to linalg since it will be used in RR
    if (chebfi%paral_kgb == 1 .and. xmpi_comm_size(chebfi%spacecom) > 1) then
      call xg_free(chebfi%X_NP) 
      
      call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs,chebfi%spacecom)
      
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs) 
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)   
      call xgBlock_zero(chebfi%X_NP%self)
    end if 
    call timab(tim_RR_XNP_reset,2,tsec)
    
    !print *, "ZGEMM(N,N,M,N,K)", "M=", spacedim, "N=", neigenpairs, "K=", neigenpairs
    call timab(tim_RR_gemm_2,1,tsec)      
    if (remainder == 1) then        
      call xgBlock_setBlock(chebfi%X_next, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) 
      
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_swap)                           
                            
    else if (remainder == 2) then 
      call xgBlock_setBlock(chebfi%X_prev, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_next, chebfi%X_swap, 1, spacedim, neigenpairs)
      
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_swap)    
     
    else if (remainder == 0) then  
      call xgBlock_setBlock(chebfi%X_prev, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_next, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%X_swap, 1, spacedim, neigenpairs)
           
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_next) 
      
      call xgBlock_copy(chebfi%X_next,chebfi%X_swap, 1, 1)    !copy cannot be avoided :(      
    end if
    
    !print *, "ZGEMM(N,N,M,N,K)", "M=", spacedim, "N=", neigenpairs, "K=", neigenpairs
    call xgBlock_gemm(chebfi%AX%self%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%AX%self, A_und_X%self, 0.d0, chebfi%AX_swap)
                                                              
    !print *, "ZGEMM(N,N,M,N,K)", "M=", spacedim, "N=", neigenpairs, "K=", neigenpairs
    if (chebfi%paw) then         
      call xgBlock_gemm(chebfi%BX%self%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%BX%self, A_und_X%self, 0.d0, chebfi%BX_swap)                   
    end if  
    call timab(tim_RR_gemm_2,2,tsec)
    
    call xg_free(A_und_X)
    call xg_free(B_und_X)
    
    
  end subroutine chebfi_rayleighRitz
  
  
  double precision function chebfi_computeResidue(chebfi, residu, pcond)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_computeResidue'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t) , intent(inout) :: residu
    
    double precision :: maxResidu, minResidu
    integer :: eigResiduMax, eigResiduMin
    
    double precision :: tsec(2)
    
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface
     
    if (chebfi%paw) then
      call xgBlock_colwiseCymax(chebfi%AX_swap,chebfi%eigenvalues,chebfi%BX_swap,chebfi%AX_swap)
    else
      call xgBlock_colwiseCymax(chebfi%AX_swap,chebfi%eigenvalues,chebfi%X_swap,chebfi%AX_swap)
    end if
   
    !pcond call
    call timab(tim_pcond,1,tsec)
    call pcond(chebfi%AX_swap)
    call timab(tim_pcond,2,tsec)
    
    call xgBlock_colwiseNorm2(chebfi%AX_swap, residu, max_val=maxResidu, max_elt=eigResiduMax,&
                                                      min_val=minResidu, min_elt=eigResiduMin) 
    
    chebfi_computeResidue = maxResidu
    
  end function chebfi_computeResidue  
  
  subroutine chebfi_ampfactor(chebfi, eig, lambda_minus, lambda_plus, nline_bands)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_ampfactor'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    double precision, pointer , intent(in) :: eig(:,:)
    double precision, intent(in) :: lambda_minus
    double precision, intent(in) :: lambda_plus
    integer, intent(in) :: nline_bands(:) 
    
    type(xgBlock_t) :: X_part
    type(xgBlock_t) :: AX_part
    type(xgBlock_t) :: BX_part
    
    integer :: iband, shift, nbands
    double precision :: ampfactor
    double precision eig_per_band
    
    if (chebfi%paral_kgb == 0) then
      nbands = chebfi%neigenpairs
    else
      nbands = chebfi%bandpp
    end if
      
    do iband = 1, nbands
          
      if(chebfi%istwf_k == 2) then
        eig_per_band = dble(eig(iband,1))
      else
        eig_per_band = dble(eig(iband*2-1,1))
      end if
      !cheb_poly1(x, n, a, b)
      ampfactor = cheb_poly1(eig_per_band, nline_bands(iband), lambda_minus, lambda_plus)
           
      if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much
              
      call xgBlock_setBlock(chebfi%xXColsRows, X_part, iband, chebfi%total_spacedim, 1) 
      call xgBlock_setBlock(chebfi%xAXColsRows, AX_part, iband, chebfi%total_spacedim, 1) 
      
      call xgBlock_scale(X_part, 1/ampfactor, 1)    
      call xgBlock_scale(AX_part, 1/ampfactor, 1)
      
      if(chebfi%paw) then
        call xgBlock_setBlock(chebfi%xBXColsRows, BX_part, iband, chebfi%total_spacedim, 1) 
        call xgBlock_scale(BX_part, 1/ampfactor, 1)
      end if
    end do
    
  end subroutine chebfi_ampfactor

  function cheb_oracle1(x, a, b, tol, nmax) result(n)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cheb_oracle1'
!End of the abilint section

   implicit none

   double precision :: tol

   integer :: nmax
   integer :: n, i
   double precision, intent(in) :: x, a, b
   double precision :: y, xred, temp
   double precision :: yim1

! *************************************************************************

   xred = (x-(a+b)/2)/(b-a)*2
   y = xred
   yim1 = 1 !ONE

   n = nmax
   if(1/(y**2) < tol) then
     n = 1
   else
     do i=2, nmax-1
       temp = y
       y = 2*xred*y - yim1
       yim1 = temp
       if(1/(y**2) < tol) then
         n = i
         exit
       end if
     end do
   end if

  end function cheb_oracle1
  !!***
  
  function cheb_poly1(x, n, a, b) result(y)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'cheb_poly1'
!End of the abilint section

   implicit none

   integer, intent(in) :: n
   integer :: i
   double precision, intent(in) :: x, a, b
   double precision :: y, xred, temp
   double precision :: yim1

! *************************************************************************

   xred = (x-(a+b)/2)/(b-a)*2
   y = xred
   yim1 = 1
   do i=2, n
     temp = y
     y = 2*xred*y - yim1
     yim1 = temp
   end do

  end function cheb_poly1
 
 
end module m_chebfi2


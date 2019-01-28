#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"


module m_chebfi2

  use m_xg
  use m_xgScalapack
  use defs_basis, only : std_err, std_out
  use defs_abitypes
  use m_abicore
  use m_errors
  use m_xomp
  use m_cgtools
#ifdef HAVE_OPENMP
  use omp_lib
#endif

  use m_time, only : timab

  implicit none
  
  private
  
  integer, parameter :: EIGENV = 1
  integer, parameter :: EIGENVD = 2
  
#ifdef HAVE_OPENMP 
  integer, save :: eigenSolver = EIGENVD     ! Type of eigen solver to use
#else
  integer, save :: eigenSolver = EIGENV     ! Type of eigen solver to use
#endif

  type, public :: chebfi_t
    integer :: space
    integer :: spacedim                      ! Space dimension for one vector
    integer :: neigenpairs                   ! Number of eigen values/vectors we want
    integer :: nline                         ! Number of line to perform
    integer :: spacecom                      ! Communicator for MPI
    double precision :: tolerance            ! Tolerance on the residu to stop the minimization
    double precision :: ecut                 ! Ecut for Chebfi oracle
    
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

    type(xgBlock_t) :: eigenvalues

    !SWAP POINTERS
    type(xgBlock_t) :: SWP_X  

  end type chebfi_t

  public :: chebfi_init
  public :: chebfi_memInfo
  public :: chebfi_run
  public :: chebfi_free

  contains

  subroutine chebfi_init(chebfi, neigenpairs, spacedim, tolerance, ecut, nline, space, eigenProblem, istwf_k, spacecom, me_g0)


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
    integer         , intent(in   ) :: nline
    integer         , intent(in   ) :: space
    integer         , intent(in   ) :: eigenProblem
    integer         , intent(in   ) :: istwf_k
    integer         , intent(in   ) :: spacecom
    integer         , intent(in   ) :: me_g0
    integer :: nthread

    chebfi%space = space
    chebfi%neigenpairs = neigenpairs
    chebfi%spacedim    = spacedim
    chebfi%tolerance   = tolerance
    chebfi%ecut        = ecut
    chebfi%nline       = nline
    chebfi%spacecom    = spacecom
    chebfi%paw = .false.
    chebfi%eigenProblem = eigenProblem
    chebfi%istwf_k = istwf_k
    chebfi%me_g0 = me_g0

    nthread = 16

    !call OMP_SET_NUM_THREADS(nthread) !REMOVE THIS IN THE END
    
    !SHOULD HERE I DO ADVICE TARGET AND BLOCK CALCULATIONS???

    call chebfi_allocateAll(chebfi)

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
    
    space = chebfi%space
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs

    call chebfi_free(chebfi) 

    call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs)
    
    call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs)  
    call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)  

    call xg_init(chebfi%AX,space,spacedim,neigenpairs)
    call xg_init(chebfi%BX,space,spacedim,neigenpairs)

  end subroutine chebfi_allocateAll


  function chebfi_memInfo(neigenpairs, spacedim, space) result(arraymem)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_memInfo'
!End of the abilint section

    integer, intent(in   ) :: neigenpairs
    integer, intent(in   ) :: spacedim
    integer, intent(in   ) :: space

    double precision :: memX_Re
    double precision :: memX_Im
    double precision :: memX

    double precision :: memX_next
    double precision :: memX_prev

    double precision :: memAX
    double precision :: memBX
    
    !chebfi_rayleightRitz function variables
    double precision :: memA_und_X
    double precision :: memB_und_X
    double precision :: memEigenvalues
    double precision :: cplx

    double precision :: arraymem(2)

    cplx = 1 ; if ( space == SPACE_C ) cplx = 2 !for now only complex

    ! Permanent in chebfi
    memX = cplx * kind(1.d0) * spacedim * neigenpairs

    memX_next = cplx * kind(1.d0) * spacedim * neigenpairs
    memX_prev = cplx * kind(1.d0) * spacedim * neigenpairs

    memAX = cplx * kind(1.d0) * spacedim * neigenpairs
    memBX = cplx * kind(1.d0) * spacedim * neigenpairs
    
    !chebfi_rayleightRitz function variables
    memA_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
    memB_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
    memEigenvalues = kind(1.d0) * neigenpairs

    arraymem(1) = memX + memX_next + memX_prev + &
                  memAX + memBX
    arraymem(2) = memA_und_X + memB_und_X + memEigenvalues

  end function chebfi_memInfo

  subroutine chebfi_run(chebfi, X0, getAX_BX, getBm1X, pcond, eigen, residu)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_run'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0   ! Full initial vectors
    type(xgBlock_t), intent(inout) :: eigen   ! Full initial eigen values
    type(xgBlock_t), intent(inout) :: residu
          
    integer :: spacedim
    integer :: neigenpairs
    integer :: nline, nline_max, nline_decrease, nline_tolwfr 
    integer :: shift
    double precision :: tolerance
    !double precision :: ecut
    double precision :: maxeig
    double precision :: mineig
    double precision :: ampfactor
    
    type(xg_t)::DivResults

    integer :: test, iline, iband
        
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
    
    integer :: info

    interface
      subroutine getAX_BX(X,AX,BX)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
      end subroutine getAX_BX
    end interface
    interface
      subroutine getBm1X(X,Bm1X)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
      end subroutine getBm1X
    end interface 
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface

    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs
    nline = chebfi%nline
    chebfi%eigenvalues = eigen
    
    allocate(nline_bands(neigenpairs))
    
    call xg_init(DivResults, chebfi%space, neigenpairs, 1)
        
    tolerance = chebfi%tolerance
    lambda_plus = chebfi%ecut
    chebfi%X = X0
    
    call getAX_BX(chebfi%X,chebfi%AX%self,chebfi%BX%self) 
        
    !********************* Compute Rayleigh quotients for every band, and set λ − equal to the largest one *****!

    call chebfi_rayleighRitzQuotiens(chebfi, neigenpairs, maxeig, mineig, DivResults%self) !OK
        
    lambda_minus = maxeig
    
    nline_max = cheb_oracle1(mineig, lambda_minus, lambda_plus, 1D-16, 40)

    call xgBlock_reverseMap(DivResults%self,eig,1,neigenpairs)
    
    do iband=1, neigenpairs !TODO
!      ! nline necessary to converge to tolerance
!      !nline_tolwfr = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, tolerance / resids_filter(iband), nline)
!      ! nline necessary to decrease residual by a constant factor
!      !nline_decrease = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, 0.1D, dtset%nline)
!      !nline_bands(iband) = MAX(MIN(nline_tolwfr, nline_decrease, nline_max, chebfi%nline), 1)
      nline_bands(iband) = nline ! fiddle with this to use locking
    end do

    center = (lambda_plus + lambda_minus)*0.5  
    radius = (lambda_plus - lambda_minus)*0.5
        
    one_over_r = 1/radius
    two_over_r = 2/radius
    
    do iline = 0, nline - 1 

      call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)
      
      call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
   
      !A * ψ      
      call getAX_BX(chebfi%X,chebfi%AX%self,chebfi%BX%self)
           
    end do
        
    call chebfi_ampfactor(chebfi, eig, lambda_minus, lambda_plus, nline_bands)    !ampfactor
        
    call chebfi_rayleightRitz(chebfi)

    maximum =  chebfi_computeResidue(chebfi, residu, pcond)
       
    call xgBlock_copy(chebfi%X,chebfi%X_prev, 1, 1) !copy back to X_prev (has to be done otherwise CG in big loop has bad data)
        
    deallocate(nline_bands)
    
    call xg_free(DivResults)
 
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
    
    call xg_finalize()

  end subroutine chebfi_free

  subroutine chebfi_rayleighRitzQuotiens(chebfi, neigenpairs, maxeig, mineig, DivResults)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_rayleighRitzQuotiens'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    integer, intent(in) :: neigenpairs
    double precision, intent(inout) :: maxeig
    double precision, intent(inout) :: mineig
    type(xgBlock_t) , intent(inout) :: DivResults

    integer :: maxeig_pos(2)
    integer :: mineig_pos(2)
    integer :: info

    type(xg_t)::Results1
    type(xg_t)::Results2

    call xg_init(Results1, chebfi%space, neigenpairs, 1)
    call xg_init(Results2, chebfi%space, neigenpairs, 1)
    
    call xgBlock_colwiseDotProduct(chebfi%X,chebfi%AX%self,Results1%self)
        
    !PAW
    call xgBlock_colwiseDotProduct(chebfi%X,chebfi%BX%self,Results2%self)

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
        
    interface
      subroutine getBm1X(X,Bm1X)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
      end subroutine getBm1X
    end interface

    !B-1 * A * ψ    
    if (chebfi%paw) then
      call getBm1X(chebfi%AX%self, chebfi%X_next)
    else
      !call xg_setBlock(chebfi%AX,chebfi%X_next, 1, chebfi%spacedim, chebfi%neigenpairs) !references problem
      call xgBlock_copy(chebfi%AX%self,chebfi%X_next, 1, 1)
    end if   
    
    !ψi-1 = c * ψi-1   
    call xgBlock_scale(chebfi%X, center, 1) !scale by c
            
    !(B-1 * A * ψi-1 - c * ψi-1)
    call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X)
        
    !ψi-1  = 1/c * ψi-1 
    call xgBlock_scale(chebfi%X, 1/center, 1) !counter scale by c
    
    if (iline == 0) then  
      call xgBlock_scale(chebfi%X_next, one_over_r, 1) !scale by one_over_r
    else
      call xgBlock_scale(chebfi%X_next, two_over_r, 1) !scale by two_over_r
      
      call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X_prev)
    end if
        
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

    !swap buffers
    call xgBlock_setBlock(chebfi%X_prev, chebfi%SWP_X, 1, spacedim, neigenpairs) !*swap_psi = X_prev;	 
    call xgBlock_setBlock(chebfi%X, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = X;   !OK
    call xgBlock_setBlock(chebfi%X_next, chebfi%X, 1, spacedim, neigenpairs) !X = X_next; !OK
    call xgBlock_setBlock(chebfi%SWP_X, chebfi%X_next, 1, spacedim, neigenpairs) !X_next = *swap_psi;

  end subroutine chebfi_swapInnerBuffers
  
  subroutine chebfi_rayleightRitz(chebfi)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_rayleightRitz'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    
    type(xg_t) :: A_und_X !H_UND_PSI
    type(xg_t) :: B_und_X !S_UND_PSI
    
    type(xgBlock_t) :: X_first_row
    type(xgBlock_t) :: AX_first_row
    
    double precision, pointer :: evec(:,:)
    double precision, allocatable :: edummy(:,:)
    
    integer :: neigenpairs
    integer :: spacedim
    integer :: i
    integer :: info
    integer :: space
    integer :: eigenProblem
    integer :: me_g0

    space = chebfi%space
    eigenProblem = chebfi%eigenProblem
    me_g0 = chebfi%me_g0
    
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs
    
    call xg_init(A_und_X,space,neigenpairs,neigenpairs)
    call xg_init(B_und_X,space,neigenpairs,neigenpairs)
    
    if(chebfi%istwf_k == 2) then
       call xgBlock_scale(chebfi%X, sqrt2, 1)
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%X, X_first_row, 1, 1, neigenpairs)
         call xgBlock_scale(X_first_row, 1/sqrt2, 1)
       end if
       call xgBlock_scale(chebfi%AX%self, sqrt2, 1)
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%AX%self, AX_first_row, 1, 1, neigenpairs)
         call xgBlock_scale(AX_first_row, 1/sqrt2, 1)
       end if  
    end if
        
    !call xgBlock_gemm(chebfi%AX%self%trans, chebfi%X%normal, 1.0d0, chebfi%AX%self, chebfi%X, 0.d0, A_und_X%self)
    
    call xgBlock_gemm(chebfi%X%trans, chebfi%AX%self%normal, 1.0d0, chebfi%X, chebfi%AX%self, 0.d0, A_und_X%self)
        
    call xgBlock_gemm(chebfi%X%trans, chebfi%X%normal, 1.0d0, chebfi%X, chebfi%X, 0.d0, B_und_X%self) !changed from original
      
     if(chebfi%istwf_k == 2) then
       call xgBlock_scale(chebfi%X, 1/sqrt2, 1)
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%X, X_first_row, 1, 1, neigenpairs)
         call xgBlock_scale(X_first_row, sqrt2, 1)
       end if
       call xgBlock_scale(chebfi%AX%self, 1/sqrt2, 1)
       if (me_g0 == 1)  then 
         call xgBlock_setBlock(chebfi%AX%self, AX_first_row, 1, 1, neigenpairs)
         call xgBlock_scale(AX_first_row, sqrt2, 1)
       end if  
    end if
    !*********************************** EIGENVALUE A Ψ X = B Ψ XΛ **********************************************!
    
    select case (eigenSolver)
    case (EIGENVD) 
      call xgBlock_hegvd(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case (EIGENV)
      call xgBlock_hegv(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case default
       MSG_ERROR("Error for Eigen Solver HEGV")
    end select
        
    call xgBlock_reverseMap(A_und_X%self,evec,1,neigenpairs*neigenpairs)
    
    !ABI_ALLOCATE(edummy, (2*neigenpairs, neigenpairs)) !2 = cplx
    !call fxphas_seq(evec,edummy,0,0,1,neigenpairs*neigenpairs,neigenpairs*neigenpairs,neigenpairs,neigenpairs,0)  
    !ABI_DEALLOCATE(edummy)
                              
    call xgBlock_gemm(chebfi%AX%self%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%AX%self, A_und_X%self, 0.d0, chebfi%X_prev)
                      
    call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_next)                      
     
    call xgBlock_copy(chebfi%X_prev,chebfi%AX%self, 1, 1) !AX = X_prev;
    call xgBlock_copy(chebfi%X_next,chebfi%X, 1, 1)  !X = X_next;
     
    call xg_free(A_und_X)
    call xg_free(B_und_X)

  end subroutine chebfi_rayleightRitz
  
  
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
    
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface
  
    !PAW
    call xgBlock_colwiseCaxmy(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
    
    !call xgBlock_colwiseCaxmy(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
        
    !pcond call
    call pcond(chebfi%AX%self)
    
    call xgBlock_colwiseNorm2(chebfi%AX%self, residu, max_val=maxResidu, max_elt=eigResiduMax,&
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
    
    integer :: iband, shift
    double precision :: ampfactor
    
    !ampfactor
    do iband = 1, chebfi%neigenpairs
      ampfactor = cheb_poly1(dble(eig(iband*2-1,1)), nline_bands(iband), lambda_minus, lambda_plus) !OK

      if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much
      !shift = chebfi%spacedim*(iband-1)
      
      call xgBlock_setBlock(chebfi%X, X_part, iband, chebfi%spacedim, 1)
      call xgBlock_setBlock(chebfi%AX%self, AX_part, iband, chebfi%spacedim, 1) 
      
      call xgBlock_scale(X_part, 1/ampfactor, 1)
      call xgBlock_scale(AX_part, 1/ampfactor, 1)	 	 
      
!   if(paw) then
!     gsc_filter(:, shift+1:shift+npw_filter*nspinor) = gsc_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
!   else
!     gvnlc_filter(:, shift+1:shift+npw_filter*nspinor) = gvnlc_filter(:, shift+1:shift+npw_filter*nspinor) / ampfactor
!   end if
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


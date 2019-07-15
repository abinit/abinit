
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

  integer, parameter :: tim_init     = 1751
  integer, parameter :: tim_free     = 1752
  integer, parameter :: tim_run      = 1753
  integer, parameter :: tim_getAX_BX = 1754
  integer, parameter :: tim_invovl   = 1755
  integer, parameter :: tim_residu   = 1756
  integer, parameter :: tim_RR       = 1757
  integer, parameter :: tim_pcond    = 1758
  integer, parameter :: tim_RR_q     = 1759
  integer, parameter :: tim_next_p   = 1760
  integer, parameter :: tim_swap     = 1761
  integer, parameter :: tim_amp_f    = 1762
  integer, parameter :: tim_alltoall = 1763

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
    
    !type(xg_t) :: X_all_to_all1, AX_all_to_all1, BX_all_to_all1
    !type(xg_t) :: X_all_to_all2, AX_all_to_all2, BX_all_to_all2
    
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
    !type(xgBlock_t) :: XXX  

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
    !integer :: nthread


    call timab(tim_init,1,tsec)
    
    chebfi%space = space
    !if (paral_kgb == 0) then
      chebfi%neigenpairs = neigenpairs  !transposer sam ovo postavlja
    !else
      !chebfi%neigenpairs = bandpp
    !end if
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
    
    !print *, "chebfi%nproc_band", chebfi%nproc_band
!    print *, "chebfi%spacedim", chebfi%spacedim
!    print *, "chebfi%neigenpairs", chebfi%neigenpairs
!    print *, "chebfi%bandpp", chebfi%bandpp
!    stop
    
    !SHOULD HERE I DO ADVICE TARGET AND BLOCK CALCULATIONS???

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
    integer :: remainder = 0

    integer :: row, col
    integer :: total_spacedim, ierr !don't know why but have to have different var for MPI spacedim
    
    space = chebfi%space
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs

    call chebfi_free(chebfi) 
    !stop
!    print *, "spacedim", spacedim
!    print *, "neigenpairs", neigenpairs
!    print *, "chebfi%bandpp", chebfi%bandpp
!    stop
    
      !print *, "neigenpairs", neirgenpairs
      !print *, "chebfi%nproc_band", chebfi%nproc_band
      !call xmpi_barrier(chebfi%spacecom)
      !stop
    
    if (chebfi%paral_kgb == 0) then
      chebfi%total_spacedim = spacedim
      call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs, chebfi%spacecom) !regular arrays
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs)  
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)  
    else
      !print *, "chebfi%nproc_band", chebfi%nproc_band
      !print *, "spacedim*chebfi%nproc_band", spacedim*chebfi%nproc_band
      !print *, "2*chebfi%bandpp", 2*chebfi%bandpp
      !stop
!      print *, "spacedim", spacedim
      !print *, "spacedim*chebfi%nproc_band", spacedim*chebfi%nproc_band
      !print *, "chebfi%bandpp", chebfi%bandpp
      !print *, "spacedim*blockdim 1", spacedim*chebfi%nproc_band*chebfi%bandpp
!      dummy = 3887
!      call xg_init(chebfi%X_NP,space,dummy,2*chebfi%bandpp) !transposed arrays
!      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,dummy,chebfi%bandpp)  
!      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,chebfi%bandpp+1,dummy,chebfi%bandpp) 
      !stop
      total_spacedim = spacedim
      call xmpi_sum(total_spacedim,chebfi%spacecom,ierr)
      !print *, "spacedim", spacedim
      !stop
      !ovde ako nije deljivo sa 2 zauzeti za X_NP jedan red vise pa onda posle prepakovati
      chebfi%total_spacedim = total_spacedim
      !chebfi%bandpp = neigenpairs/chebfi%nproc_band
      !ovako nesto mora da bi se dobio tacan broj redova
!      if (MOD(total_spacedim,xmpi_comm_size(xmpi_world)) /=0) then
!        remainder = MOD(total_spacedim,xmpi_comm_size(xmpi_world))
!      end if
      !print *, "chebfi%bandpp", chebfi%bandpp
      !stop
      !ne moze se dodati jos jedan dummy red jer posle ldim propadne za x_next i x_prev
      call xg_init(chebfi%X_NP,space,total_spacedim,2*chebfi%bandpp,chebfi%spacecom) !transposed arrays
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,total_spacedim,chebfi%bandpp)  
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,chebfi%bandpp+1,total_spacedim,chebfi%bandpp)  
      !print *, "prosao"
      !stop
    end if
    !print *, "spacedim", spacedim
    !stop
    
    !transposer will handle these arrays automatically
    call xg_init(chebfi%AX,space,spacedim,neigenpairs,chebfi%spacecom)
    call xg_init(chebfi%BX,space,spacedim,neigenpairs,chebfi%spacecom)
    !stop

    !print *, "AX spacedim", spacedim
    !stop
    !neigenpairs * spacedim * spacedim * neigenpairs
  end subroutine chebfi_allocateAll

  !!TODO TODO TODO TODO TODO fix meminfo it is for OpenMP only
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

    cplx = 1 
    if ( space == SPACE_C ) cplx = 2 !for now only complex

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
    double precision :: maxeig, maxeig_global
    double precision :: mineig, mineig_global
    double precision :: ampfactor
    
    type(xg_t)::DivResults

    integer :: test, iline, iband
    
    integer :: ikpt_this_proc, node_spacedim
    
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
    
    double precision :: errmax
    double precision, allocatable, target :: cg(:,:)
    double precision, pointer :: cg0(:,:)
    
    logical :: status_p
    
    integer :: remainder
    
    integer :: info
    integer :: ierr
    
    integer :: cols, rows
    
    type(xgBlock_t) :: HELPER
    type(xg_t) :: xCOPY
   
    double precision :: tsec(2)
                    
    interface
      subroutine getAX_BX(X,AX,BX,npband,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: AX
        type(xgBlock_t), intent(inout) :: BX
        type(integer), intent(inout) :: npband
        type(xgTransposer_t), optional, intent(inout) :: transposer
      end subroutine getAX_BX
    end interface
    interface
      !call getBm1X(chebfi%xAXColsRows, chebfi%X_next, iline_t, xXColsRows, X, chebfi%xgTransposerX) !ovde pobrljavi X skroz za
      subroutine getBm1X(X,Bm1X,iline_t,xXColsRows,X1,npband,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
        type(integer), intent(inout) :: iline_t
        type(xgBlock_t), intent(inout) :: xXColsRows
        type(xgBlock_t), intent(inout) :: X1
        type(integer), intent(inout) :: npband
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
    
    !stop
    
    !if (xmpi_comm_rank(chebfi%spacecom) == 1) then
      !call xgBlock_print(chebfi%eigenvalues, 6)
    !end if
    
    !call xmpi_barrier(chebfi%spacecom)
    
    !stop
!  call omp_set_num_threads(4)

!  print *, "Native", omp_get_num_threads()
!  write(*,*) "Num threads", OMP_GET_MAX_THREADS()
!  WRITE(*,*) "Num procs", omp_get_num_procs()
!  
!  write ( *, * ) 'A sequential hello to you!'
!  !$omp parallel
!  print *, "Native", omp_get_num_threads()
!  write ( *, * ) ' Parallel hello''s to you!'
!  !$omp end parallel
!  
!  stop
!    print *, "Native", omp_get_num_threads()
!    print *, "Current", xomp_get_num_threads(open_parallel=.true.)
!    write(*,*) "Num threads", OMP_GET_MAX_THREADS()
!    WRITE(*,*) "Num procs", omp_get_num_procs()
!    stop

    !print *, "KRENUO"
    !stop
    
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
    
!    if (associated(cg0)) then
!      print *, "CGA"
!    else
!      print *, "CGNA"
!    end if
!    stop
    
    !print *, "xmpi_comm_rank(comm)", xmpi_comm_rank(chebfi%spacecom)
    !print *, "xmpi_comm_size(chebfi%spacecom)", xmpi_comm_size(chebfi%spacecom)
    !print *, "xmpi_comm_size(xmpi_world)", xmpi_comm_size(xmpi_world)
    !stop
   
    !call xg_getPointer(chebfi%X)
    ! Initialize the _filter pointers. Depending on paral_kgb, they might point to the actual arrays or to _alltoall variables

    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
    !stop
    
    !call debug_me_g0_LINALG(chebfi%X, chebfi)  !OK 1-2 ISTWFK2 MPI
    !stop
    
    !call debug_helper_colrows(chebfi%X_prev, chebfi) 
    !call debug_helper_colrows(chebfi%X_next, chebfi) 
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !stop  
    
    !print *, "AJDE" 
   
    !call xgBlock_getSize(chebfi%X,rows,cols) 
    !print *, "rows", rows
    !print *, "cols", cols
    !stop  
    
    
    !call debug_helper_colrows(chebfi%X, chebfi) 
    !print *, "STA SE DESAVA"
    !stop
    !stop
    
    !if (iline == 0) then
      !call debug_helper_linalg(chebfi%X, chebfi) 
      !call xg_getPointer(chebfi%X)
      !call xg_getPointer(chebfi%X_next)
      !stop 
    !end if
    !stop
     
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
    !stop
    !print *, "RANK", xmpi_comm_rank(chebfi%spacecom)
    !stop
    ! Transpose
    if (chebfi%paral_kgb == 1) then
      
      nCpuRows = chebfi%nproc_fft
      nCpuCols = chebfi%nproc_band
  
      !print *, "spacedim", spacedim
      !print *, "neigenpairs", neigenpairs
      
      !call xg_init(xCOPY, chebfi%space, spacedim, neigenpairs)
      !call xgBlock_reverseMap(chebfi%X,cg,spacedim,neigenpairs)
      !call xgBlock_copy(chebfi%X,xCOPY%self, 1, 1)
      !call xgBlock_reverseMap(xCOPY%self,cg0,spacedim,neigenpairs)
      !kopirati X u poseban array pa to posle porediti nakon trans
  
      !print *, "nCpuRows", nCpuRows
      !print *, "nCpuCols", nCpuCols
      !stop
      
      !!TODO promesati ovde X transponovati ga i nazad i opet koristiti isti debag. Deluje da je 
      ! problem negde u transponovanju ili sta vec, ali AX kasnije ne fukcionise kao X
      !call debug_helper_linalg(chebfi%X, chebfi, 1)
      !stop
    
      
      !all stride arrays are hidden inside transposer
      call xgTransposer_init(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nCpuRows,nCpuCols,STATE_LINALG,1,0)
      !print *, "PROSAO INIT"
      !call xgTransposer_init(chebfi%xgTransposerX,X0,chebfi%X,nCpuRows,nCpuCols,STATE_LINALG,1,0)
      
 !     call xgBlock_getSize(chebfi%X,rows,cols)
      
!      print *, "CCCCCCCCCCCCCCC"
!      write (100+xmpi_comm_rank(xmpi_world),*), "ROWS 3", rows
!      write (100+xmpi_comm_rank(xmpi_world),*), "COLS 3", cols
!      call xgBlock_getSize(xXColsRows,rows,cols)
!      write (100+xmpi_comm_rank(xmpi_world),*), "ROWS 4", rows
!      write (100+xmpi_comm_rank(xmpi_world),*), "COLS 4", cols
!      flush (100+xmpi_comm_rank(xmpi_world))
!      stop
      !print *, "EEEEEEEEEEEEEEEEE"
      call xgTransposer_init(chebfi%xgTransposerAX,chebfi%AX%self,chebfi%xAXColsRows,nCpuRows,nCpuCols,STATE_LINALG,1,0)
      
      !print *, "DDDDDDDDDDDDDDDD"
      call xgTransposer_init(chebfi%xgTransposerBX,chebfi%BX%self,chebfi%xBXColsRows,nCpuRows,nCpuCols,STATE_LINALG,1,0)
      
      !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
      !stop
      !print *, "SSSSSSSSSSSSSS"
      
      call xgTransposer_transpose(chebfi%xgTransposerX,STATE_COLSROWS) !all_to_all
      call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_COLSROWS) !all_to_all
      call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_COLSROWS) 
      
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi)   
      !stop 
      !DEBUG DEBUG DEBUG
      !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all 
      !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
      !stop
!      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc  
!        call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
!        !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
!      else
!        call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
!        !call debug_helper_colrows(chebfi%AX%self, chebfi) 
!      end if
!      
!      stop
      
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
      !stop
      
      
      !errmax = (sum(cg0-cg))/neigenpairs
      !print *, "loc cg0", loc(cg0)
      !print *, "loc cg", loc(cg)   
      !print *, " Difference: ",errmax
      
      !stop
      
      !print *, "BBBBBBBBBBBBBB"
                        
!      write (100+xmpi_comm_rank(xmpi_world),*), "ROWS 1", rows
!      write (100+xmpi_comm_rank(xmpi_world),*), "COLS 1", cols
!      call xgBlock_getSize(xXColsRows,rows,cols)
!      write (100+xmpi_comm_rank(xmpi_world),*), "ROWS 2", rows
!      write (100+xmpi_comm_rank(xmpi_world),*), "COLS 2", cols
!      flush (100+xmpi_comm_rank(xmpi_world))

      !call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, chebfi%bandpp) !ovo treba jos nekako promeniti
      !call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, 1, spacedim, chebfi%bandpp)
      !call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, 1, spacedim, chebfi%bandpp)
      !stop
    else
      call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, 1, spacedim, neigenpairs)   !use xXColsRows instead of X notion
      call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, 1, spacedim, neigenpairs)   !use xAXColsRows instead of AX notion
      call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, 1, spacedim, neigenpairs)
    end if
    
    !print *, "xXColsRows after transpose"
    !call xg_getPointer(chebfi%xXColsRows)
  
    
    !stop
    
    !call debug_me_g0_COLROWS(chebfi%xXColsRows, chebfi) !OK 1-2 ISTWFK2 MPI
    !stop
    
    !print *, "AAAAAAAAAAAAAAAAAAAA"
    !stop
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !stop
   
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
    
    call timab(tim_getAX_BX,1,tsec)         
    !call getAX_BX(chebfi%X,chebfi%AX%self,chebfi%BX%self) 
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows,chebfi%nproc_band,chebfi%xgTransposerX)  !OVO SAD MORA SVUDA DA SE MENJA
    call timab(tim_getAX_BX,2,tsec)
    
    
    print *, "PROSAO getAX_BX"
    !stop
    !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
    !call debug_helper_linalg(chebfi%X, chebfi, 1) 
    !stop
    
    !!on nproc 2 for some reason
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
    !stop
    
    !print *, "getAXBX"
    !call debug_me_g0_COLROWS(chebfi%xXColsRows, chebfi)
    !stop
    
      
    if (chebfi%paral_kgb == 1) then   
      call xmpi_barrier(chebfi%spacecom)
    end if
    !stop
               
    !********************* Compute Rayleigh quotients for every band, and set λ − equal to the largest one *****!
    call timab(tim_RR_q, 1, tsec)
    call chebfi_rayleighRitzQuotiens(chebfi, maxeig, mineig, DivResults%self) !OK
    call timab(tim_RR_q, 2, tsec)
    
    !print *, "PROSAO RRQ"
    !stop
    !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
    !call debug_helper_linalg(chebfi%X, chebfi, 1) 
    !stop
    
    !stop
    if (chebfi%paral_kgb == 1) then
      call xmpi_max(maxeig,maxeig_global,chebfi%spacecom,ierr)
      call xmpi_min(mineig,mineig_global,chebfi%spacecom,ierr)
    else
      maxeig_global = maxeig
      mineig_global = mineig
    end if
    
    !print *, "maxeig_global", maxeig_global  
    !print *, "mineig_global", mineig_global    
    !stop
                
    lambda_minus = maxeig_global
    
    !print *, "lambda_minus", lambda_minus
    !stop
        
    nline_max = cheb_oracle1(mineig, lambda_minus, lambda_plus, 1D-16, 40)
    
    !print *, "nline_max", nline_max
    !stop
   
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
    
    !print *, "center", center
    !print *, "radius", radius
    !stop
               
    one_over_r = 1/radius
    two_over_r = 2/radius
    

    !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
    !stop
     !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)    
   ! print *, "AAAAAAAAAAAAAAAAAAAA"
   
    !if (xmpi_comm_rank(xmpi_world) == 0) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
    print *, "************************"
    print *, "BITNI POINTERI PRE LOOPA"
    print *, "xXcolsRows, rank", xmpi_comm_rank(xmpi_world)
    call xg_getPointer(chebfi%xXColsRows)
    print *, "X, rank", xmpi_comm_rank(xmpi_world)
    call xg_getPointer(chebfi%X)
    print *, "X_next, rank", xmpi_comm_rank(xmpi_world)
    call xg_getPointer(chebfi%X_next)
    print *, "X_prev, rank", xmpi_comm_rank(xmpi_world)
    call xg_getPointer(chebfi%X_prev)
    print *, "X_NP, rank", xmpi_comm_rank(xmpi_world)
    call xg_getPointer(chebfi%X_NP%self)
    !stop
    print *, "***********************"
    !end if
    
    !stop
    !stop
    do iline = 0, nline - 1 
    
            
      if (iline == 2) then
!        print *, "2 before chebfi_computeNextOrderChebfiPolynom"
!        call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!        call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!        call xmpi_barrier(chebfi%spacecom)
!        call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!        stop    
      end if
      
      !print *, "X_next BEFORE chebfi_computeNextOrderChebfiPolynom, rank", xmpi_comm_rank(xmpi_world)
      !call xg_getPointer(chebfi%X_next)
      !call debug_helper_colrows(chebfi%X_next, chebfi)
      !stop 
      
      call timab(tim_next_p,1,tsec)         
      call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X) !inline 2 problem istwfk 2 MPI
      call timab(tim_next_p,2,tsec)   
      
      !print *, "X_next AFTER chebfi_computeNextOrderChebfiPolynom, rank", xmpi_comm_rank(xmpi_world)
      !call xg_getPointer(chebfi%X_next)
      !call debug_helper_colrows(chebfi%X_next, chebfi)
      !stop 
      
      !if (iline == 2) then
!        print *, "2 after chebfi_computeNextOrderChebfiPolynom"
!        call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!        call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!        call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!        call xmpi_barrier(chebfi%spacecom)
!        stop    
      !end if
      
      
      !stop
      !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
      !call debug_helper_linalg(chebfi%X, chebfi, 1) 
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
      !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
      !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
      !stop
      
      !stop      
      !call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)  !OVO SAD MORA SVUDA DA SE MENJA  
      !print *, "PROSAO SVE"
      !stop
      !print *, "spacedim swap", spacedim
      !stop
      
      !call xg_getPointer(chebfi%X)
      !call xg_getPointer(chebfi%xXColsRows)
      
      !print *, "chebfi%bandpp SIB", chebfi%bandpp
      !print *, "chebfi%total_spacedim", chebfi%total_spacedim
      !stop
      
      !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all 
      !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
      !stop
      
      !if (iline == 2) then
        !print *, "2 before swap"
        !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
        !stop
        !call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
        !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
        !call xmpi_barrier(chebfi%spacecom)
        !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
        !stop    
      !end if

      !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)  
      !stop
      
      !call debug_helper_colrows(chebfi%X_prev, chebfi) 
      !call debug_helper_colrows(chebfi%X_next, chebfi) 
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
      !stop  
      
      call timab(tim_swap,1,tsec)  
      if (chebfi%paral_kgb == 0) then             
        call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
      else 
        !print *, "chebfi%bandpp SIB", chebfi%bandpp
        !print *, "chebfi%total_spacedim", chebfi%total_spacedim
        !stop
        !call debug_helper_colrows(chebfi%X_next, chebfi)  !!TODO MYSTIRIOUS BUG HERE
        !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
        !stop 
        call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
        !call debug_helper_colrows(chebfi%X_next, chebfi)
        !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
        !stop 
      end if
      call timab(tim_swap,2,tsec) 
      
      !print *, "iline", iline
      !stop     
      
!      if (iline == 2) then
!        print *, "2 after swap"
!        
!        if (xmpi_comm_rank(xmpi_world) == 0) then !only one MPI
!        
!        print *, "************************"
!        print *, "BITNI POINTERI POSLE SWAPA"
!        print *, "xXcolsRows"
!        call xg_getPointer(chebfi%xXColsRows)
!        print *, "X"
!        call xg_getPointer(chebfi%X)
!        print *, "X_next"
!        call xg_getPointer(chebfi%X_next)
!        print *, "X_prev"
!        call xg_getPointer(chebfi%X_prev)
!        print *, "X_NP"
!        call xg_getPointer(chebfi%X_NP%self)
!        print *, "***********************"
!        
!        end if
!        !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
!        !stop
!        call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!        call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!        call xmpi_barrier(chebfi%spacecom)
!        call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!        stop    
!      end if
            
      !if (iline == 2) then
!          if (xmpi_comm_rank(xmpi_world) == 0) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
!            print *, "************************"
!            print *, "BITNI POINTERI POSLE SWAPA"
!            print *, "xXcolsRows"
!            call xg_getPointer(chebfi%xXColsRows)
!            print *, "X"
!            call xg_getPointer(chebfi%X)
!            print *, "X_next"
!            call xg_getPointer(chebfi%X_next)
!            print *, "X_prev"
!            call xg_getPointer(chebfi%X_prev)
!            print *, "X_NP"
!            call xg_getPointer(chebfi%X_NP%self)
!            print *, "***********************"
!        end if
        
        !if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps) 
          !call xgBlock_setBlock(chebfi%X_prev, chebfi%X, 1, spacedim, neigenpairs) 
          !call xgBlock_copy(chebfi%X_prev, chebfi%X, 1, 1)   
        !end if
        !print *, "USAOASASASA"
        !call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
        !call xgBlock_copy(chebfi%X_prev, chebfi%xXColsRows, 1, 1)  
        !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
        !call debug_helper_colrows(chebfi%X, chebfi, 1, 1)
        !call xmpi_barrier(chebfi%spacecom)
        !stop
      
      !end if
      
!      call xgBlock_setBlock(chebfi%X_prev, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp) 
!      call xgBlock_copy(chebfi%X_prev, chebfi%xXColsRows, 1, 1)  
!      call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!      call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!      stop
      
      
      !call xg_getPointer(chebfi%X)
      !print *, "xXColsRows before LINE 905"
      !call xg_getPointer(chebfi%xXColsRows)
      !stop
      
 
      !print *, "BUFFERS SWAPPED"
      !call xmpi_barrier(chebfi%spacecom)
      !stop  
            
      !A * ψ  
      call timab(tim_getAX_BX,1,tsec)             
      !call getAX_BX(chebfi%X,chebfi%AX%self,chebfi%BX%self)
      !stop
      !print *, "getAX_BX before"
      call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows,chebfi%nproc_band,chebfi%xgTransposerX)  !OVO SAD MORA SVUDA DA SE MENJA
      !print *, "getAX_BX after"
      
      print *, "xXColsRows before LINE 905, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%xXColsRows)
      !call xgBlock_reverseMap(chebfi%xXColsRows,cg,spacedim,chebfi%bandpp)
      if (xg_associated(chebfi%xXColsRows)) then
        print *, "xXColsRows associated!!!, rank", xmpi_comm_rank(xmpi_world) 
      end if
      call xmpi_barrier(chebfi%spacecom)
      !stop
      !stop
      !call debug_helper_colrows(chebfi%X_prev, chebfi) 
      !call debug_helper_colrows(chebfi%X_next, chebfi) 
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
      !stop
      !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
      !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
      call timab(tim_getAX_BX,2,tsec) 
      
      !stop  
      
!      if (iline == 2) then
!        print *, "2 after getAX_BX"
!        
!        if (xmpi_comm_rank(xmpi_world) == 0) then !only one MPI
!        
!        print *, "************************"
!        print *, "BITNI POINTERI POSLE SWAPA"
!        print *, "xXcolsRows"
!        call xg_getPointer(chebfi%xXColsRows)
!        print *, "X"
!        call xg_getPointer(chebfi%X)
!        print *, "X_next"
!        call xg_getPointer(chebfi%X_next)
!        print *, "X_prev"
!        call xg_getPointer(chebfi%X_prev)
!        print *, "X_NP"
!        call xg_getPointer(chebfi%X_NP%self)
!        print *, "***********************"
!        
!        end if
!        !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
!        !stop
!        !call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!        call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!        call xmpi_barrier(chebfi%spacecom)
!        call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!        stop    
!      end if
      
    end do
    
    !print *, "TRANSPOSE AFTER LOOP"
    
!    remainder = mod(nline-1, 3) 3 buffer swap, keep the info which one contains X_data at the end of loop
!    
!    if (chebfi%paral_kgb == 1) then             
!      if (remainder == 0) then
!        call xgBlock_setBlock(chebfi%X_prev, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!      else if (remainder == 1) then
!        call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
!      else 
!        do nothing 
!      end if
!    end if
    
!    call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!    call xmpi_barrier(chebfi%spacecom)
!    call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!    stop    
!    call xgBlock_copy(chebfi%X_next, chebfi%xXColsRows, 1, 1)  
!    call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp) 
!    call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!    call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
!    stop
!    
 
!    print *, "************************"
!    print *, "BITNI POINTERI POSLE LOOPA"
!    print *, "xXcolsRows"
!    call xg_getPointer(chebfi%xXColsRows)
!    print *, "X"
!    call xg_getPointer(chebfi%X)
!    print *, "X_next"
!    call xg_getPointer(chebfi%X_next)
!    print *, "X_prev"
!    call xg_getPointer(chebfi%X_prev)
!    print *, "X_NP"
!    call xg_getPointer(chebfi%X_NP%self)
!    !stop
!    print *, "***********************"
    !stop

    
!    call xgBlock_setBlock(chebfi%xXColsRows, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
!    call xgBlock_print(HELPER, 200 + chebfi%paral_kgb) 
!    stop
     
!    call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!    call debug_helper_linalg(chebfi%X, chebfi, 1) 
!    stop
    
    !call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) !all_to_all
    !call debug_helper_linalg(chebfi%AX%self, chebfi, 1) 
    !stop
    

    !print *, "LOOP FINISHED"
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
    !stop
    
    !print *, "LOOP FINISHED"
    !stop

    if (chebfi%paral_kgb == 1) then  
      call xmpi_barrier(chebfi%spacecom)
    end if 
    
   !call xgBlock_setBlock(chebfi%X_prev, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
    
!   if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
!     call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) !all_to_all
!     !call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, chebfi%total_spacedim, chebfi%bandpp)  
!     !call xgBlock_copy(chebfi%xXColsRows, chebfi%X, 1, 1)  
!     !call xgBlock_copy(chebfi%xAXColsRows, chebfi%AX%self, 1, 1)  
!     !call xgBlock_copy(chebfi%xBXColsRows, chebfi%BX%self, 1, 1)  
!    else
!      !call xgBlock_copy(chebfi%X_prev, chebfi%xXColsRows, 1, 1)  
!      !call xgBlock_setBlock(chebfi%X_prev, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp) 
!      !call xgBlock_copy(chebfi%X_next, chebfi%xXColsRows, 1, 1)  
!      call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) !all_to_all
!    end if
!    
!    call xmpi_barrier(chebfi%spacecom)    
!    print *, "************************"
!    print *, "BITNI POINTERI POSLE TRANSPOSE"
!    print *, "xXcolsRows"
!    call xg_getPointer(chebfi%xXColsRows)
!    print *, "X"
!    call xg_getPointer(chebfi%X)
!    print *, "X_next"
!    call xg_getPointer(chebfi%X_next)
!    print *, "X_prev"
!    call xg_getPointer(chebfi%X_prev)
!    print *, "X_NP"
!    call xg_getPointer(chebfi%X_NP%self)
!    !stop
!    print *, "***********************"
!    
!    call debug_helper_linalg(chebfi%AX%self, chebfi, 1, 1)
!    
!    stop
    !print *, "EIG", eig
!    if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
!      call xgBlock_print(DivResults%self, 200)
!    else
!      call xgBlock_print(DivResults%self, 100+xmpi_comm_rank(chebfi%spacecom))
!    end if
!    call xmpi_barrier(chebfi%spacecom)
!    stop

!    print *, "DEBUG REACH"
    !if (xmpi_comm_size(xmpi_world) == 1 .or. chebfi%paral_kgb == 0) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
      !print *, "VRATIO ADRESU"
      !stop
      !call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
    !end if
!      
!    call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!    call debug_helper_linalg(chebfi%X, chebfi, 1) 
    !stop
    
    !call xgBlock_setBlock(chebfi%xXColsRows, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
    !call xgBlock_print(HELPER, 200 + chebfi%paral_kgb) 
    !stop
    !BAG RESEN, SPACEDIM PROBLEM BIO U AMPFACTORU
    !eig should be OK here DivResults are split in the right way over MPI processes  
    call timab(tim_amp_f,1,tsec)
    call chebfi_ampfactor(chebfi, eig, lambda_minus, lambda_plus, nline_bands)    !ampfactor
    call timab(tim_amp_f,2,tsec)
    
    !call xgBlock_setBlock(chebfi%xXColsRows, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
    !call xgBlock_print(HELPER, 200 + chebfi%paral_kgb) 
    !stop
   
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi)
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi)
    !stop 
   
!    print *, "DEBUG REACH"
!    if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
!      call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
!    end if
!      
!    call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
!    call debug_helper_linalg(chebfi%X, chebfi, 1) 
!    stop
!    
!    call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) !all_to_all
!    call debug_helper_linalg(chebfi%AX%self, chebfi, 1) 
!    stop
!   
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
    !call debug_helper_colrows(chebfi%xBXColsRows, chebfi) 
    !print *, "AMPFACTOR FINISHED" 
    !stop
   
    !call xgBlock_setBlock(chebfi%xXColsRows, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
    !call xgBlock_print(HELPER, 200) 
    
    !print *, "xXColsRows before backtranspose"
    !call xg_getPointer(chebfi%xXColsRows)
    !stop
    !Filtering done transpose back 
    if (chebfi%paral_kgb == 1) then
    
      call xmpi_barrier(chebfi%spacecom)
     
      !WHICHEVER POINTER IS INSIDE xXColsRows it will be used for transpose. It doesn't have to be same as in beginning
      !call debug_helper_linalg(chebfi%X, chebfi) 
      call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
      call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG) !all_to_all !no need to transpose (only X)
      call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG) !all_to_all
      
      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
        !call xgBlock_copy(chebfi%xXColsRows, chebfi%X, 1, 1)  
        !call xgBlock_copy(chebfi%xAXColsRows, chebfi%AX%self, 1, 1)  
        !call xgBlock_copy(chebfi%xBXColsRows, chebfi%BX%self, 1, 1)  
      end if
     
    else
      !if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)  !OVO JE OK ONO IZNAD NIJE
        !call xg_getPointer(chebfi%xXColsRows)
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, 1, spacedim, neigenpairs)  
        !call xgBlock_copy(chebfi%xXColsRows, chebfi%X, 1, 1)  
      !end if
    end if
    
!    print *, "************************"
!    print *, "BITNI POINTERI"
!    print *, "xXcolsRows"
!    call xg_getPointer(chebfi%xXColsRows)
!    print *, "X"
!    call xg_getPointer(chebfi%X)
!    print *, "X_next"
!    call xg_getPointer(chebfi%X_next)
!    print *, "X_prev"
!    call xg_getPointer(chebfi%X_prev)
!    print *, "X_NP"
!    call xg_getPointer(chebfi%X_NP%self)
!    !stop
!    print *, "***********************"
!    
!    stop
    
    
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
    !stop
  
    
    !print *, "PROSAO TRANSPOSE"
    !stop
    
    !!TODO TRANSPOSE IS NOT WORKING CORRECTLY FOR ISTWFK2 (or some other thing)
    !print *, "AFTER BACKTRANSPOSE"
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)  !!TODO FROM HERE
    !call debug_helper_linalg(chebfi%AX%self, chebfi, 1, 1) 
    !call debug_helper_linalg(chebfi%BX%self, chebfi, 1, 1) 
    !stop 
    
    !call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_COLSROWS) !proba ali ne radi bas kako treba :(
    !stop
   
    call timab(tim_RR, 1, tsec)
    call chebfi_rayleightRitz(chebfi, nline)
    call timab(tim_RR, 2, tsec) 

    !print *, "PRE RESIDUE"
    !stop
    
    call timab(tim_residu, 1, tsec)
    maximum =  chebfi_computeResidue(chebfi, residu, pcond)
    call timab(tim_residu, 2, tsec)
                       
    deallocate(nline_bands)
    
    call xg_free(DivResults)
    
    if (chebfi%paral_kgb == 1) then
      call xgTransposer_free(chebfi%xgTransposerX)
      call xgTransposer_free(chebfi%xgTransposerAX)
      call xgTransposer_free(chebfi%xgTransposerBX)
    end if
    
    call timab(tim_run,2,tsec)
    
    
    !X_swap izgleda ne mora da se koristi kod MPI > 1 jer je X vec poseban bafer
    !if (xmpi_comm_rank(chebfi%spacecom) == 0) then
      print *, "X_swap AT THE END"
      call xg_getPointer(chebfi%X_swap)
    
      print *, "************************"
      print *, "BITNI POINTERI"
      print *, "xXcolsRows, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%xXColsRows)
      print *, "X, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%X)
      print *, "X_next, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%X_next)
      print *, "X_prev, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%X_prev)
      print *, "X_NP, rank", xmpi_comm_rank(xmpi_world)
      call xg_getPointer(chebfi%X_NP%self)
      !stop
      print *, "***********************"
    !end if
    
    if (xmpi_comm_size(chebfi%spacecom) > 1) then
      call xgBlock_copy(chebfi%X_swap,chebfi%X, 1, 1)    !copy cannot be avoided :(
    end if
    
    print *, "xgx0 INSIDE CB2 AT THE END"
    call xg_getPointer(chebfi%X) 
    call debug_helper_linalg(chebfi%X, chebfi, 1, 1) !MPI 1,2 OK
    !stop   
     
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

    !print *, "neigenpairs", neigenpairs
    !stop

    if (chebfi%paral_kgb == 0) then
      !print *, "chebfi%neigenpairs", chebfi%neigenpairs
      !stop
      call xg_init(Results1, chebfi%space, chebfi%neigenpairs, 1)
      call xg_init(Results2, chebfi%space, chebfi%neigenpairs, 1)
    else
      call xg_init(Results1, chebfi%space, chebfi%bandpp, 1)
      call xg_init(Results2, chebfi%space, chebfi%bandpp, 1)
    end if
    
    !call xgBlock_colwiseDotProduct(chebfi%X,chebfi%AX%self,Results1%self)
    !call xgBlock_colwiseDotProduct(chebfi%X,chebfi%xAXColsRows,Results1%self)
    call xgBlock_colwiseDotProduct(chebfi%xXColsRows,chebfi%xAXColsRows,Results1%self)
        
    !PAW
    !call xgBlock_colwiseDotProduct(chebfi%X,chebfi%BX%self,Results2%self)
    !call xgBlock_colwiseDotProduct(chebfi%X,chebfi%xBXColsRows,Results2%self)
    call xgBlock_colwiseDotProduct(chebfi%xXColsRows,chebfi%xBXColsRows,Results2%self)

    !PAW
    call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, maxeig, maxeig_pos, mineig, mineig_pos)
    
    !print *, "maxeig", maxeig
    !stop
        
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
    
    integer :: iline_t
             
    interface
      subroutine getBm1X(X,Bm1X,iline_t,xXColsRows,X1,npband,transposer)
        use m_xg, only : xgBlock_t
        use m_xgTransposer !, only: xgTransposer_t
        type(xgBlock_t), intent(inout) :: X
        type(xgBlock_t), intent(inout) :: Bm1X
        type(integer), intent(inout) :: iline_t
        type(xgBlock_t), intent(inout) :: xXColsRows
        type(xgBlock_t), intent(inout) :: X1
        type(integer), intent(inout) :: npband
        type(xgTransposer_t), optional, intent(inout) :: transposer    
      end subroutine getBm1X
    end interface


    !B-1 * A * ψ    
    !stop
    if (chebfi%paw) then
      call timab(tim_invovl, 1, tsec)
      !call getBm1X(chebfi%AX%self, chebfi%X_next)
      !print *, "SSSSSSSSS"
      !stop
!      if (iline == 2) then
!        call debug_helper_linalg(chebfi%X, chebfi) 
!        call xg_getPointer(chebfi%X)
!        call xg_getPointer(chebfi%X_next)
!        stop 
!      end if
      !stop
      !if (xmpi_comm_rank(chebfi%spacecom) == 0) then
        print *, "ILINE", iline
        print *, "X_next before getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%X_next)
        print *, "X_prev before getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%X_prev)
        print *, "X_xXColsRows before getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%xXColsRows)
      !end if
      if (iline == 2) then
     
!        call debug_helper_colrows(chebfi%X_next, chebfi)
        print *, "chebfi%total_spacedim",  chebfi%total_spacedim
        print *, "chebfi%bandpp", chebfi%bandpp
        !call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
        !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
        !call xmpi_barrier(chebfi%spacecom)
        !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
        !stop    
      end if
      
      iline_t = iline
      !!TODO UCI U WRAPPER SA ILINE=2 i transposerom i probati transponovati unutra
      call getBm1X(chebfi%xAXColsRows, chebfi%X_next, iline_t, chebfi%xXColsRows, chebfi%X, chebfi%nproc_band, chebfi%xgTransposerX) !ovde pobrljavi X skroz za iline 2 MPI istwfk 2
      
      !if (xmpi_comm_rank(chebfi%spacecom) == 0) then
        print *, "ILINE", iline
        print *, "X_next after getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%X_next)
        print *, "X_prev after getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%X_prev)
        print *, "X_xXColsRows after getBm1X, rank", xmpi_comm_rank(xmpi_world)
        call xg_getPointer(chebfi%xXColsRows)
      !end if
      if (iline == 2) then
        !call debug_helper_colrows(chebfi%X_next, chebfi)
        !call xmpi_barrier(chebfi%spacecom)
        !call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, chebfi%total_spacedim, chebfi%bandpp)  
        !call xgTransposer_transpose(chebfi%xgTransposerX,STATE_LINALG) !all_to_all
        !call xmpi_barrier(chebfi%spacecom)
        !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)
        !stop    
      end if
      
      !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
      !call debug_helper_colrows(chebfi%X_next, chebfi) 
      !stop
!      if (iline == 2) then
!        call debug_helper_linalg(chebfi%X, chebfi) 
!        call xg_getPointer(chebfi%X)
!        call xg_getPointer(chebfi%X_next)
!        stop 
!      end if

      !print *, "PPPPP"
      !stop
      !call debug_helper_colrows(chebfi%xAXColsRows, chebfi) 
      !call debug_helper_colrows(chebfi%X_next, chebfi) 
      !stop
      !call xgBlock_setBlock(chebfi%X_next, chebfi%AX_swap, 1, chebfi%spacedim, chebfi%neigenpairs) !AX_swap = X_next;
      call timab(tim_invovl, 2, tsec)
    else
      !call xgBlock_copy(chebfi%AX%self,chebfi%X_next, 1, 1)
      call xgBlock_copy(chebfi%xAXColsRows,chebfi%X_next, 1, 1)
      !!TODO try to swap buffers in a way that last copy is avoided
      !call xgBlock_setBlock(chebfi%AX%self, chebfi%AX_swap, 1, chebfi%spacedim, chebfi%neigenpairs) !AX_swap = AX;
      !call xgBlock_setBlock(chebfi%X_next, chebfi%XXX, 1, chebfi%spacedim, chebfi%neigenpairs) !AX_swap = AX;
    end if  
    

    
    !print *, "PROSAO APPLY"
    !stop 
            
    !ψi-1 = c * ψi-1   
    !call xgBlock_scale(chebfi%X, center, 1) !scale by c
    call xgBlock_scale(chebfi%xXColsRows, center, 1) !scale by c
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !print *, "PROSAO SCALE"
    !stop
            
    !(B-1 * A * ψi-1 - c * ψi-1)
    !call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X)
    call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%xXColsRows)  
    !call xgBlock_saxpy(chebfi%AX_swap, dble(-1.0), chebfi%X)
    !call xgBlock_saxpy(chebfi%XXX, dble(-1.0), chebfi%X)
    !call debug_helper_colrows(chebfi%X_next, chebfi) 
    !print *, "PROSAO SAXPY"
    !stop
        
    !ψi-1  = 1/c * ψi-1 
    !call xgBlock_scale(chebfi%X, 1/center, 1) !counter scale by c
    call xgBlock_scale(chebfi%xXColsRows, 1/center, 1) !counter scale by c
    !call debug_helper_colrows(chebfi%xXColsRows, chebfi) 
    !stop
    
    if (iline == 0) then  
      call xgBlock_scale(chebfi%X_next, one_over_r, 1) 
      !call xgBlock_scale(chebfi%AX_swap, one_over_r, 1) !scale by one_over_r
      !call xgBlock_scale(chebfi%XXX, one_over_r, 1) !scale by one_over_r
      !call debug_helper_colrows(chebfi%X_next, chebfi) 
      !print *, "SCALE NEXT"
      !stop
    else
      call xgBlock_scale(chebfi%X_next, two_over_r, 1) 
      !call xgBlock_scale(chebfi%AX_swap, two_over_r, 1) !scale by one_over_r
      !call xgBlock_scale(chebfi%XXX, two_over_r, 1) !scale by one_over_r
      
      call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X_prev)
      !call xgBlock_saxpy(chebfi%AX_swap, dble(-1.0), chebfi%X_prev)
      !call xgBlock_saxpy(chebfi%XXX, dble(-1.0), chebfi%X_prev)
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
    
    !print *, "SWAP START"
    !stop
    !WORKING OPENMP
!    call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !X_swap = X_prev	      
!    call xgBlock_setBlock(chebfi%X, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = X
!    call xgBlock_setBlock(chebfi%X_next, chebfi%X, 1, spacedim, neigenpairs) !X = X_next
!    call xgBlock_setBlock(chebfi%X_swap, chebfi%X_next, 1, spacedim, neigenpairs) !X_next = X_swap
    
    !MPI
    !print *, "spacedim", spacedim
    !print *, "neigenpairs", neigenpairs
    !stop
    call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !X_swap = X_prev     
    call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = xXColsRows
    call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, spacedim, neigenpairs) !xXColsRows = X_next
    call xgBlock_setBlock(chebfi%X_swap, chebfi%X_next, 1, spacedim, neigenpairs) !X_next = X_swap
    
    !!TODO try to swap buffers in a way that last copy is avoided
    !swap buffers
!    if (chebfi%paw) then
!      call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !X_swap = X_prev	
!      call xgBlock_setBlock(chebfi%X, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = X
!      call xgBlock_setBlock(chebfi%AX_swap, chebfi%X, 1, spacedim, neigenpairs) !X = AX_swap = X_next
!      call xgBlock_setBlock(chebfi%X_swap, chebfi%AX_swap, 1, spacedim, neigenpairs) !AX_swap = X_next = X_swap
!    else 
!      if (mod(chebfi%nline, 3) == 1) then
!        call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !*X_swap = X_prev
!      else !if ((mod(chebfi%nline, 3) == 2) .or. (mod(chebfi%nline, 3) == 0)) 
!        call xgBlock_setBlock(chebfi%X_next, chebfi%X_swap, 1, spacedim, neigenpairs) !*X_swap = X_prev
!      end if
!      
!      call xgBlock_setBlock(chebfi%X, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = X 
!      call xgBlock_setBlock(chebfi%AX_swap, chebfi%X, 1, spacedim, neigenpairs) !X = AX_swap = AX
!      
!      if (mod(chebfi%nline, 3) == 1) then
!        call xgBlock_setBlock(chebfi%X_swap, chebfi%AX%self, 1, spacedim, neigenpairs) !X = AX_swap = AX
!      else !if (mod(chebfi%nline, 3) == 2 .or. mod(chebfi%nline, 3) == 0) 
!        call xgBlock_setBlock(chebfi%X_swap, chebfi%AX%self, 1, spacedim, neigenpairs) !*X_swap = X_prev
!      end if     
!    end if

  end subroutine chebfi_swapInnerBuffers
    
  subroutine chebfi_rayleightRitz(chebfi, nline)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'chebfi_rayleightRitz'
!End of the abilint section

    type(chebfi_t) , intent(inout) :: chebfi
    type(integer) , intent(in) :: nline
    
    type(xg_t) :: A_und_X !H_UND_PSI
    type(xg_t) :: B_und_X !S_UND_PSI
    
!    type(xg_t) :: dummy_AX
!    type(xg_t) :: dummy_X
!    type(xg_t) :: dummy_A_und_X
!    
    type(xgBlock_t) :: X_first_row
    type(xgBlock_t) :: AX_first_row
    type(xgBlock_t) :: BX_first_row
    
    type(xgBlock_t) :: HELPER
                
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
    
    integer :: multiply, size_remainder
    integer :: X_NP_TEST = 1
    
    

    space = chebfi%space
    eigenProblem = chebfi%eigenProblem
    me_g0 = chebfi%me_g0
    
    spacedim = chebfi%spacedim
    !if (chebfi%paral_kgb == 0) then
      neigenpairs = chebfi%neigenpairs   !remains whole nband domain since it is after transpose
    !else 
      !neigenpairs = chebfi%bandpp
    !end if
    
    !print *, "spacedim", spacedim
    !print *, "neigenpairs", neigenpairs
    !stop
   
    call xg_init(A_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)
    call xg_init(B_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)
    call xgBlock_zero(A_und_X%self)
    
    
    !call xg_init(dummy_AX, space, 5, 3, chebfi%spacecom)
    !call xg_init(dummy_X, space, 5, 3, chebfi%spacecom)
    !call xg_init(dummy_A_und_X, space, 3, 3, chebfi%spacecom)
    !print *, "neigenpairs", neigenpairs
    !stop
    !call xgBlock_zero(dummy_AX%self)
    !call xgBlock_one(dummy_AX%self)
    !call xgBlock_zero(dummy_X%self)
    !call xgBlock_one(dummy_X%self)
    !if (xmpi_comm_rank(chebfi%spacecom) == 0) then
      !call xgBlock_scale(dummy_AX%self, dble(2), 1)
      !call xgBlock_scale(chebfi%X, dble(2), 1)
    !else
      !call xgBlock_scale(dummy_X%self, dble(4), 1)
    !end if
       
    !ispada da je prvi deo AX jednak za proc 0 i 1?!?!?!
   ! call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
    !call debug_helper_linalg(chebfi%AX%self, chebfi, 1, 1) 
    !call debug_helper_linalg(chebfi%BX%self, chebfi, 1, 1) 
    !stop
    
    !!TODO TODO TODO TODO A_und_X se ne izracuna dobro iz nekog razloga. Ulazni nizovi deluju dobro DEBUG
    
    !call xg_getPointer(chebfi%AX%self)
    !call xg_getPointer(chebfi%X)
    !stop
    
    !if (xmpi_comm_size(xmpi_world) == 2) then
      !if (xmpi_comm_rank(xmpi_world) == 0) then
        !print *, "RANK 0"
        !call xgBlock_one(chebfi%AX%self)
      !else
        !print *, "RANK 1"
        !call xgBlock_zero(chebfi%AX%self)
      !end if
    !end if
    !pointers are not important they seem ok
   
    !print *, "AAAAAAAAAAAAAA"
    !stop
    !call xgBlock_print(A_und_X%self, 100+xmpi_comm_rank(chebfi%spacecom))
    !stop
    
    !print *, "BEFORE GEMM"
    !stop
    
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(A_und_X%self, 1000+xmpi_comm_rank(chebfi%spacecom))
!    else
!      call xgBlock_print(A_und_X%self, 2000+xmpi_comm_rank(chebfi%spacecom))
!    end if
    
    call xgBlock_gemm(chebfi%AX%self%trans, chebfi%X%normal, 1.0d0, chebfi%AX%self, chebfi%X, 0.d0, A_und_X%self) 
    !call xgBlock_gemm(dummy_AX%self%trans, dummy_X%self%normal, 1.0d0, dummy_AX%self, dummy_X%self, 0.d0, dummy_A_und_X%self) 
    
    !call xmpi_barrier(chebfi%spacecom)
    !stop
     
    !BELOW NOT NECESSARY IT'S AUTOMATIC INSIDE GEMM!!!
!    if (chebfi%paral_kgb == 1) then
!      print *, "SPACE", space
!      !stop
!    
!      call xgBlock_reverseMap(A_und_X%self,sub_pointer,neigenpairs,neigenpairs)
!      
!      print *, "PROSAO 1"
!      !stop
!      
!      !call xmpi_sum_master(sub_pointer, 0, chebfi%spacecom, ierr)
!      call xmpi_sum(sub_pointer, chebfi%spacecom, ierr)
!    
!      print *, "PROSAO 2"    
!    end if
    
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(A_und_X%self, 1000+xmpi_comm_rank(chebfi%spacecom)+1)
!    else
!      call xgBlock_print(A_und_X%self, 2000+xmpi_comm_rank(chebfi%spacecom)+1)
!    end if
!    stop
    !call debug_helper_linalg(A_und_X%self, chebfi, 1, 1)
    
    !print *, "PROSAO GEMM"
    !stop
            
    if (chebfi%paw) then
      call xgBlock_gemm(chebfi%BX%self%trans, chebfi%X%normal, 1.0d0, chebfi%BX%self, chebfi%X, 0.d0, B_und_X%self) 
      !call xgBlock_gemm(chebfi%xBXColsRows%trans, chebfi%X%normal, 1.0d0, chebfi%xBXColsRows, chebfi%X, 0.d0, B_und_X%self) 
    else
      call xgBlock_gemm(chebfi%X%trans, chebfi%X%normal, 1.0d0, chebfi%X, chebfi%X, 0.d0, B_und_X%self) 
    end if
   
         
!    print *, "PROSAO GEMM 2"
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(B_und_X%self, 1000+xmpi_comm_rank(chebfi%spacecom)+1)
!    else
!      call xgBlock_print(B_und_X%self, 2000+xmpi_comm_rank(chebfi%spacecom)+1)
!    end if
!    stop
          
    if(chebfi%istwf_k == 2) then    
       call xgBlock_scale(chebfi%X, 1/sqrt2, 1)   
       if (me_g0 == 1)  then   !in linalg no need to check cpurow, first one is always me_g0 (hopefully)
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
    
    print *, "PROSLO SKALIRANJE"
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
    !call debug_helper_linalg(chebfi%AX%self, chebfi, 1, 1) 
    !call debug_helper_linalg(chebfi%BX%self, chebfi, 1, 1) 
    !stop
        
    !*********************************** EIGENVALUE A Ψ X = B Ψ XΛ **********************************************!
    
    !if (xmpi_comm_rank(chebfi%spacecom) == 1) then
      !print *, "USAO"
      !call xgBlock_print(chebfi%eigenvalues, 6)
      !print *, "IZASAO"
    !end if
    
    call xmpi_barrier(chebfi%spacecom)
    !print *, "info", inf
    !print *, "PROSAO HEGV"
    !stop
    
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(A_und_X%self, 1000+xmpi_comm_rank(chebfi%spacecom))
!    else
!      call xgBlock_print(A_und_X%self, 2000+xmpi_comm_rank(chebfi%spacecom))
!    end if
!    stop
    
    
    !MPI OK 1-2 processes
    select case (eigenSolver)
    case (EIGENVD) 
      !print *, "EIGENVD"
      call xgBlock_hegvd(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case (EIGENV)
      call xgBlock_hegv(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
    case default
       MSG_ERROR("Error for Eigen Solver HEGV")
    end select
    
    !print *, "comm", chebfi%spacecom
    !print *, "RANK", xmpi_comm_rank(chebfi%spacecom)
    !print *, "SIZE", xmpi_comm_size(chebfi%spacecom)
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(chebfi%eigenvalues, 1000+xmpi_comm_rank(chebfi%spacecom)+1)  !EV GOOD FOR 1 and 2 MPI
!    else 
!      call xgBlock_print(chebfi%eigenvalues, 2000+xmpi_comm_rank(chebfi%spacecom)+1)  !EV GOOD FOR 1 and 2 MPI
!    end if
    
    !call xmpi_barrier(chebfi%spacecom)
    !print *, "info", info
    !print *, "PROSAO HEGV"
    !stop
    
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)  !X OK HERE FOR 1 and 2 MPI
    !stop
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(A_und_X%self, 1000+xmpi_comm_rank(chebfi%spacecom)+1)
!    else
!      call xgBlock_print(A_und_X%self, 2000+xmpi_comm_rank(chebfi%spacecom)+1)
!    end if
    
!    if (xmpi_comm_size(xmpi_world) == 1) then
!      call xgBlock_print(chebfi%eigenvalues, 1000+xmpi_comm_rank(chebfi%spacecom)+1)
!    else
!      call xgBlock_print(chebfi%eigenvalues, 2000+xmpi_comm_rank(chebfi%spacecom)+1)
!    end if
!    
!    stop
    
        
!    call xgBlock_reverseMap(A_und_X%self,evec,1,neigenpairs*neigenpairs)
!    
!    ABI_ALLOCATE(edummy, (2*neigenpairs, neigenpairs)) !2 = cplx
!    call fxphas_seq(evec,edummy,0,0,1,neigenpairs*neigenpairs,neigenpairs*neigenpairs,neigenpairs,neigenpairs,0)  
!    ABI_DEALLOCATE(edummy)
!    stop
                
    !print *, "PRE SWAPA"
    !stop            
                   
    remainder = mod(nline, 3) !3 buffer swap, keep the info which one contains X_data at the end of loop
        
    !print *, "remainder", remainder
    !print *, "spacedim", spacedim
    !print *, "neigenpairs", neigenpairs
    !stop   
    
    !print *, "chebfi%paral_kgb", chebfi%paral_kgb
    !stop
     !if (X_NP_TEST == 1) then  
    if (chebfi%paral_kgb == 1 .and. xmpi_comm_size(xmpi_world) > 1) then  !!TODO TODO TODO ovde je problem za 1, izgube se adrese ili sta vec...
!      multiply = chebfi%total_spacedim * (chebfi%bandpp * 2 + 1)
!      size_remainder = multiply - (spacedim * 2 * neigenpairs)
!      print *, "X_NP old size", multiply
!      print *, "X_NP new size", spacedim * 2 * neigenpairs
!      print *, "size_remainder", size_remainder
!      print *, "size_remainder/spacedim", size_remainder/spacedim
!      print *, "MOD", MOD(multiply, spacedim)
!      stop
      !call xgBlock_reshape(chebfi%X_next,(/ spacedim, neigenpairs /))
      !call xgBlock_reshape(chebfi%X_prev,(/ spacedim, neigenpairs /))
      !call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,total_spacedim,chebfi%bandpp)  
      !call xgBlock_reshape(chebfi%X_NP%self,(/ spacedim, 2*neigenpairs /))
      call xg_free(chebfi%X_NP)   !I don't know how to avoid this
      
!      print *, "chebfi%X"
!      call xg_getPointer(chebfi%X)
!      print *, "chebfi%X_next"
!      call xg_getPointer(chebfi%X_next) 
!      print *, "chebfi%X_prev"
!      call xg_getPointer(chebfi%X_prev) 
!      print *, "chebfi%xXColsRows"
!      call xg_getPointer(chebfi%xXColsRows)
!      print *, "chebfi%X_NP%self"
!      call xg_getPointer(chebfi%X_NP%self)
!      
!      call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
!      stop    

      !print *, "spacedim", spacedim
      !stop
      !!TODO DEBUG FOR MPI > 1 now. 1 works (I don't know exactly how)
      call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs,chebfi%spacecom) !transposed arrays
      call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs) 
      call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)   
      !print *, "spacedim", spacedim
      !print *, "neigenpairs", neigenpairs
      !stop
      !print *, "PROSAO AASASASASASA"
      !stop
      !call xg_setBlock(chebfi%X_NP,chebfi%X_prev,chebfi%bandpp+1,total_spacedim,chebfi%bandpp)  
      
      !call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs) !regular arrays
      !call xg_setBlock(chebfi%X_NP,chebfi%X_next,1,spacedim,neigenpairs)  
      !call xg_setBlock(chebfi%X_NP,chebfi%X_prev,neigenpairs+1,spacedim,neigenpairs)  
      !!TODO UOPSTE NE ZNAM DA LI OVO MOZE DA SE ODRADI SA POSTOJECIM FUNKCIJAMA
      !print *, "USAAAAAAAAAAAAAAAAAAAAAAAAA"
    end if 
    !end if
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1)  !X OK HERE FOR 1 and 2 MPI
    !stop
 
    !call xgBlock_setBlock(chebfi%X, HELPER, 1, 2, 5) 
    !call xgBlock_print(HELPER, 201)
    !stop
    
    !call debug_helper_linalg(chebfi%X, chebfi, 1, 1) 
    !stop
    !print *, "REMAINDER", remainder      
    if (remainder == 1) then  !save original cg address into temp
      !print *, "USAO REMEINDER 1"      
      call xgBlock_setBlock(chebfi%X_next, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) 
      
    
      !print *, "SAD RADI"
      !stop
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_swap)      !put X into expected buffer directly
                      
      !call debug_helper_linalg(chebfi%X_swap, chebfi, 1, 1)  !X_swap OK HERE FOR 1 and 2 MPI
      !stop
                            
    else if (remainder == 2) then 
      call xgBlock_setBlock(chebfi%X_prev, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_next, chebfi%X_swap, 1, spacedim, neigenpairs)
      
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_swap)     !put X into expected buffer directly
     
    else if (remainder == 0) then  
      call xgBlock_setBlock(chebfi%X_prev, chebfi%AX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X_next, chebfi%BX_swap, 1, spacedim, neigenpairs) 
      call xgBlock_setBlock(chebfi%X, chebfi%X_swap, 1, spacedim, neigenpairs)
      
      call xgBlock_gemm(chebfi%X%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%X, A_und_X%self, 0.d0, chebfi%X_next)   !put X into expected buffer directly
            
      call xgBlock_copy(chebfi%X_next,chebfi%X_swap, 1, 1)    !copy cannot be avoided :(
            
    end if
    
    !call debug_helper_linalg(chebfi%X_swap, chebfi, 1, 1) 
    !print *, "X_SWAP print"
    !stop
                          
    call xgBlock_gemm(chebfi%AX%self%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%AX%self, A_und_X%self, 0.d0, chebfi%AX_swap)
        
    !print *, "AX_swap"              
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1)  !AX_swap OK HERE FOR 1 and 2 MPI
    !stop
      
    !call xgBlock_setBlock(chebfi%AX_swap, HELPER, 1, 2, 5) 
    !call xgBlock_print(HELPER, 200 + chebfi%paral_kgb)
    !stop
                                                        
    if (chebfi%paw) then     
      !print *, "PAW"    
      call xgBlock_gemm(chebfi%BX%self%normal, A_und_X%self%normal, 1.0d0, & 
                      chebfi%BX%self, A_und_X%self, 0.d0, chebfi%BX_swap)                   
    end if  
    
    !print *, "BX_swap"
    !call debug_helper_linalg(chebfi%BX_swap, chebfi, 1, 1)  !BX_swap OK HERE FOR 1 and 2 MPI
    !stop
    
!    print *, "AX_swap"
!    call xg_getPointer(chebfi%AX_swap)
!    print *, "BX_swap"
!    call xg_getPointer(chebfi%BX_swap)
!    print *, "X pointer"
!    call xg_getPointer(chebfi%X)
!    print *, "X colwise pointer"
!    call xg_getPointer(chebfi%xXColsRows)
    
    !call xgBlock_setBlock(chebfi%AX_swap, HELPER, 1, 2, 5) 
    !call xgBlock_print(HELPER, 200 + chebfi%paral_kgb)
    !stop

    
    !call debug_helper_linalg(chebfi%BX_swap, chebfi, 1, 1) 
    !print *, "BX_swap print"
    !stop              
     
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
    
    type(xgBlock_t) :: HELPER
    
    double precision :: maxResidu, minResidu
    integer :: eigResiduMax, eigResiduMin
    
    double precision :: tsec(2)
    
    interface
      subroutine pcond(W)
        use m_xg, only : xgBlock_t
        type(xgBlock_t), intent(inout) :: W
      end subroutine pcond
    end interface
    
    !!TODO ODAVDE NA DALJE
    !print *, "chebfi_computeResidue AX_SWAP print"
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) !MPI 1,2 OK
    !stop       
    
    !PAW
    !print *, "PRE COLWISE"
    !stop
    
    
    !call xgBlock_setBlock(chebfi%AX_swap, HELPER, 1, 2, 5) 
    !call xgBlock_print(HELPER, 200 + chebfi%paral_kgb)
    !stop
    
    !print *, "AX_swap residue"
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) 
    !stop  
    
    if (chebfi%paw) then
      call xgBlock_colwiseCymax(chebfi%AX_swap,chebfi%eigenvalues,chebfi%BX_swap,chebfi%AX_swap)
    else
      call xgBlock_colwiseCymax(chebfi%AX_swap,chebfi%eigenvalues,chebfi%X_swap,chebfi%AX_swap)
    end if
    !print *, "POSLE COLWISE"
    !stop
    
    !print *, "AFTER COLWISE"
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) !MPI 1,2 OK
    !stop      
    
    !pcond call
    call timab(tim_pcond,1,tsec)
    call pcond(chebfi%AX%self)
    call timab(tim_pcond,2,tsec)
    
    !print *, "AX_swap after pcond"
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) !MPI 1,2 OK
    !stop
    
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) 
    !print *, "pcond AX_SWAP print"
    !stop 
    !print *, "BEFORE NORM"
    !stop
    call xgBlock_colwiseNorm2(chebfi%AX_swap, residu, max_val=maxResidu, max_elt=eigResiduMax,&
                                                      min_val=minResidu, min_elt=eigResiduMin) 
    
    !print *, "AX_swap after xgBlock_colwiseNorm2"
    !call debug_helper_linalg(chebfi%AX_swap, chebfi, 1, 1) !MPI 1,2 OK
    !stop
    
    
    !print *, "AFTER NORM"
    !stop
    chebfi_computeResidue = maxResidu
    
    !print *, "MAXRES", maxResidu   
    !stop

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
    
    !print *, "nbands", nbands
    !print *, "chebfi%spacedim", chebfi%spacedim
    !stop
    !print *, "chebfi%total_spacedim", chebfi%total_spacedim
    !stop
    
    do iband = 1, nbands
      !print *, "iband", iband
      
    
      if(chebfi%istwf_k == 2) then
        eig_per_band = dble(eig(iband,1))
      else
        eig_per_band = dble(eig(iband*2-1,1))
      end if
      !cheb_poly1(x, n, a, b)
      ampfactor = cheb_poly1(eig_per_band, nline_bands(iband), lambda_minus, lambda_plus) !ovo je valjda OK
      
      !if (xmpi_comm_rank(xmpi_world) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        !print *, "AMPFACTOR, iband", ampfactor, iband
      !end if
      !ovi ampfactori bi trebali da budu OK s obzirom da su izvedeni iz DivResults (eig) koji je ok podeljen
      
      if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much
      !shift = chebfi%spacedim*(iband-1)
        
      !!TODO OVDE JE BAG
      call xgBlock_setBlock(chebfi%xXColsRows, X_part, iband, chebfi%total_spacedim, 1)
      !call xgBlock_setBlock(chebfi%AX%self, AX_part, iband, chebfi%spacedim, 1) 
      call xgBlock_setBlock(chebfi%xAXColsRows, AX_part, iband, chebfi%total_spacedim, 1) 
      
      call xgBlock_scale(X_part, 1/ampfactor, 1)    
      !call debug_helper_colrows(chebfi%xXColsRows, chebfi) X_next
      !stop
     
      call xgBlock_scale(AX_part, 1/ampfactor, 1)	 	 
      
      if(chebfi%paw) then
        !call xgBlock_setBlock(chebfi%BX%self, BX_part, iband, chebfi%spacedim, 1) 
        call xgBlock_setBlock(chebfi%xBXColsRows, BX_part, iband, chebfi%total_spacedim, 1) 
        call xgBlock_scale(BX_part, 1/ampfactor, 1)
      end if !missing gvnlc???
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
  
  subroutine debug_me_g0_LINALG(debugBlock, chebfi)
  
    type(xgBlock_t) , intent(inout) :: debugBlock
    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t) :: HELPER
    
    call xmpi_barrier(chebfi%spacecom)
    
    !if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
    if (xmpi_comm_rank(chebfi%spacecom) == 0) then
      call xgBlock_setBlock(debugBlock, HELPER, 1, 2, 5) 
      call xgBlock_print(HELPER, 200+xmpi_comm_size(xmpi_world)) 
    end if
    !else 
    
    !end if
    call xmpi_barrier(chebfi%spacecom)   
    
  end subroutine
  
  subroutine debug_me_g0_COLROWS(debugBlock, chebfi)
  
    type(xgBlock_t) , intent(inout) :: debugBlock
    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t) :: HELPER
    
    call xmpi_barrier(chebfi%spacecom)
    
    if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
      if (xmpi_comm_rank(chebfi%spacecom) == 0) then
        call xgBlock_setBlock(debugBlock, HELPER, 1, 2, 5) 
        call xgBlock_print(HELPER, 100) 
      
        call xgBlock_setBlock(debugBlock, HELPER, chebfi%bandpp/2+1, 2, 5) 
        call xgBlock_print(HELPER, 101) 
      end if
    else 
      call xgBlock_setBlock(debugBlock, HELPER, 1, 2, 5) 
      call xgBlock_print(HELPER, 200 + xmpi_comm_rank(chebfi%spacecom)) 
    end if
    
    call xmpi_barrier(chebfi%spacecom)   
    
  end subroutine debug_me_g0_COLROWS
  
  subroutine debug_helper_colrows(debugBlock, chebfi)
      
    type(xgBlock_t) , intent(inout) :: debugBlock
    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t) :: HELPER
    
    call xmpi_barrier(chebfi%spacecom)
    
    if (xmpi_comm_size(chebfi%spacecom) == 1) then !only one MPI proc
      print *, "size 1 -> chebfi%bandpp/2+1", chebfi%bandpp/2+1
      call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
      call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom)) 
      
      call xgBlock_setBlock(debugBlock, HELPER, chebfi%bandpp/2+1, DEBUG_ROWS, DEBUG_COLUMNS) 
      call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom)+1) 
    else
      !print *, "2222222222222222222"
      call xgBlock_setBlock(debugBlock, HELPER, 1, DEBUG_ROWS, DEBUG_COLUMNS) 
      call xgBlock_print(HELPER, 200+xmpi_comm_rank(chebfi%spacecom)) 
    end if
    
    !print *, "debugBlock%rows", rows(debugBlock)
    !print *, "debugBlock%cols", cols(debugBlock)
    call xmpi_barrier(chebfi%spacecom)

  end subroutine debug_helper_colrows
  
  subroutine debug_helper_linalg(debugBlock, chebfi, set_option, first_row)
      
    type(xgBlock_t) , intent(inout) :: debugBlock
    type(chebfi_t) , intent(inout) :: chebfi
    type(integer) , intent(in) :: set_option
    type(integer), intent(in) :: first_row
    type(xgBlock_t) :: HELPER
    
    integer, parameter :: FROW = 1, FCOL = 1, DROWS = 1, DCOLS = 10
    
    call xmpi_barrier(chebfi%spacecom) 
    
    !print *, "chebfi%bandpp", chebfi%bandpp
    print *, "spacedim", chebfi%spacedim
    !stop
    !stop
    
    if (set_option == 0) then
      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
        ! *, "PISE 1"
        call xgBlock_setBlock(debugBlock, HELPER, 1, DROWS, DCOLS) 
        call xgBlock_print(HELPER, 200+xmpi_comm_rank(chebfi%spacecom)) 
      
        call xgBlock_setBlock(debugBlock, HELPER, chebfi%neigenpairs/2+1, DROWS, DCOLS) 
        call xgBlock_print(HELPER, 200+xmpi_comm_rank(chebfi%spacecom)+1) 
      else
        !print *, "2222222222222222222"
        !print *, "xmpi_comm_rank(chebfi%spacecom)", xmpi_comm_rank(chebfi%spacecom)
        if (xmpi_comm_rank(chebfi%spacecom) == 0) then
          call xgBlock_setBlock(debugBlock, HELPER, 1, DROWS, DCOLS) 
          call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom))
        
          call xgBlock_setBlock(debugBlock, HELPER, chebfi%neigenpairs/2+1, DROWS, DCOLS) 
          call xgBlock_print(HELPER, 100+xmpi_comm_rank(chebfi%spacecom)+1)  
        end if
      end if
    else
      if (xmpi_comm_size(xmpi_world) == 1) then !only one MPI proc
        print *, "STAMPAM MPI1"
        print *, "chebfi%spacedim/2", chebfi%spacedim/2
        call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
        call xgBlock_print(HELPER, 100) 
        
        call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim/2+1, FCOL, DROWS, DCOLS) 
        call xgBlock_print(HELPER, 101)
      
        !call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim/2+2, FCOL, DROWS, DCOLS) 
        !call xgBlock_print(HELPER, 101)
        !!TODO OVO IZNAD JE OK OVO ISPOD NE VALJA NE ZNAM ZASTO POLUDECU!!!!
        
        !call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim/2+2, FCOL, DROWS, DCOLS) 
        !call xgBlock_print(HELPER, 110) 
      
        !call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim/2+3, FCOL, DROWS, DCOLS) 
        !call xgBlock_print(HELPER, 111)
         
      else
        !!TODO STOP DEBUGGING PARTS OF AX AND TRY TO PRINT THE SUM
        print *, "STAMPAM MPI2"
        print *, "rank, spacedim", xmpi_comm_rank(chebfi%spacecom), chebfi%spacedim
        if (xmpi_comm_rank(chebfi%spacecom) == 0) then
          !print *, "chebfi%spacedim-1", chebfi%spacedim-1
          call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
          call xgBlock_print(HELPER, 200)
          
          !call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim, FCOL, DROWS, DCOLS) 
          !call xgBlock_print(HELPER, 201)
          
          !print *, "chebfi%spacedim", chebfi%spacedim
          !print *, "FCOL", FCOL
          
          !call xgBlock_setBlock1(debugBlock, HELPER, chebfi%spacedim, FCOL, DROWS, DCOLS) 
          !call xgBlock_print(HELPER, 201)  
        else !!TODO OVO IZNAD JE OK OVO ISPOD NE VALJA NE ZNAM ZASTO POLUDECU!!!!
          !print *, "FROW, FCOL, DROWS, DCOLS", FROW, FCOL, DROWS, DCOLS
          call xgBlock_setBlock1(debugBlock, HELPER, 1, FCOL, DROWS, DCOLS) 
          call xgBlock_print(HELPER, 201)
        
          !call xgBlock_setBlock1(debugBlock, HELPER, FROW+1, FCOL, DROWS, DCOLS) 
          !call xgBlock_print(HELPER, 211)
        end if
      end if
    end if
    
    !print *, "debugBlock%rows", rows(debugBlock)
    !print *, "debugBlock%cols", cols(debugBlock)
    call xmpi_barrier(chebfi%spacecom)

  end subroutine debug_helper_linalg
  
end module m_chebfi2


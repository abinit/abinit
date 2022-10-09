!!****f* ABINIT/m_chebfi2
!! NAME
!! m_chebfi2
!!
!! FUNCTION
!! This module contains the types and routines used to apply the
!! Chebyshev filtering method (2021 implementation using xG abstraction layer)
!! It mainly defines a 'chebfi' datatypes and associated methods.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2022 ABINIT group (BS)
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

module m_chebfi2

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_time, only : timab

 use m_cgtools
 use m_xg
 use m_xgTransposer
 use m_xgScalapack
 use m_cgtools

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_GPU_NVTX_V3)
 use m_nvtx_data
#endif

 implicit none

 private

!Several (private) parameters
!-------------------------------------------------

 integer, parameter :: EIGENV = 1
 integer, parameter :: EIGENVD = 2
 integer, parameter :: DEBUG_ROWS = 5
 integer, parameter :: DEBUG_COLUMNS = 5

!Type of eigen solver to use
#ifdef HAVE_OPENMP
 integer, save :: eigenSolver = EIGENVD
#else
 integer, save :: eigenSolver = EIGENV
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

!Public 'chebfi' datatype
!-------------------------------------------------

 type, public :: chebfi_t
   integer :: space
   integer :: spacedim                      ! Space dimension for one vector
   integer :: total_spacedim                ! Maybe not needed
   integer :: neigenpairs                   ! Number of eigen values/vectors we want
   integer :: nline                         ! Number of line to perform
   integer :: spacecom                      ! Communicator for MPI
   real(dp) :: tolerance            ! Tolerance on the residu to stop the minimization
   real(dp) :: ecut                 ! Ecut for Chebfi oracle

   integer :: paral_kgb                     ! MPI parallelization variables
   integer :: nproc_band
   integer :: bandpp
   integer :: nproc_fft

   logical :: paw
   integer :: eigenProblem   !1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
   integer :: istwf_k
   integer :: use_gpu_cuda
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

!Public methods associated to 'chebfi' datatype
!-------------------------------------------------
 public :: chebfi_init
 public :: chebfi_free
 public :: chebfi_memInfo
 public :: chebfi_run

 CONTAINS  !========================================================================================
!!***

!!****f* m_chebfi2/chebfi_init
!! NAME
!! chebfi_init
!!
!! FUNCTION
!! Initialize a 'chebfi' datastructure.
!!
!! INPUTS
!!  bandpp= number of 'bands' handled by a processor
!!  eigenProblem= type of eigenpb: 1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
!!  istwf_k= parameter that describes the storage of wfs
!!  me_g0= 1 if this processors treats G=0, 0 otherwise
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  nline= Chebyshev polynomial level (.i.e. number of H applications)
!!  nproc_band= size of "band" communicator
!!  nproc_fft= size of "FFT" communicator
!!  paral_kgb= flag controlling (k,g,bands) parallelization
!!  space= defines in which space we are (columns, rows, etc.)
!!  spacecom= MPI communicator
!!  spacedim= space dimension for one vector
!!  paw= flag. TRUE if current calculation ses the PAW approach
!!  ecut= plane-wave cut-off energy
!!  tolerance= tolerance criterion on the residu to stop the minimization
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_init(chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,nproc_band,bandpp,nproc_fft, &
                       nline,space,eigenProblem,istwf_k,spacecom,me_g0,paw,use_gpu_cuda)

 implicit none

!Arguments ------------------------------------
 integer       , intent(in   ) :: bandpp
 integer       , intent(in   ) :: eigenProblem
 integer       , intent(in   ) :: istwf_k
 integer       , intent(in   ) :: me_g0
 integer       , intent(in   ) :: neigenpairs
 integer       , intent(in   ) :: nline
 integer       , intent(in   ) :: nproc_band
 integer       , intent(in   ) :: nproc_fft
 integer       , intent(in   ) :: paral_kgb
 integer       , intent(in   ) :: space
 integer       , intent(in   ) :: spacecom
 integer       , intent(in   ) :: spacedim
 integer       , intent(in   ) :: use_gpu_cuda
 logical       , intent(in   ) :: paw
 real(dp)      , intent(in   ) :: ecut
 real(dp)      , intent(in   ) :: tolerance
 type(chebfi_t), intent(inout) :: chebfi

!Local variables-------------------------------
 real(dp)                    :: tsec(2)

! *********************************************************************

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
 chebfi%use_gpu_cuda = use_gpu_cuda

 call chebfi_allocateAll(chebfi)

 call timab(tim_init,2,tsec)

end subroutine chebfi_init
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_allocateAll
!! NAME
!! chebfi_allocateAll
!!
!! FUNCTION
!! Allocate all memory spaces in a 'chebfi' datastructure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_allocateAll(chebfi)

 implicit none

!Arguments ------------------------------------
 type(chebfi_t)  , intent(inout) :: chebfi

!Local variables-------------------------------
!scalars
 integer :: neigenpairs
 integer :: space
 integer :: spacedim
 integer :: total_spacedim, ierr
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

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
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_free
!! NAME
!! chebfi_free
!!
!! FUNCTION
!! Destroy a 'chebfi' datastructure.
!!
!! INPUTS
!!
!! OUTPUT
!!  arraymem(2)= memory information
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_free(chebfi)

 implicit none

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi

! *********************************************************************

 call xg_free(chebfi%X_NP)

 call xg_free(chebfi%AX)
 call xg_free(chebfi%BX)

!call xg_finalize()

end subroutine chebfi_free
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_memInfo
!! NAME
!! chebfi_memInfo
!!
!! FUNCTION
!! Provides memory information about a 'chebfi' datastructure.
!!
!! INPUTS
!!  bandpp= number of 'bands' handled by a processor
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  paral_kgb= flag controlling (k,g,bands) parallelization
!!  space= defines in which space we are (columns, rows, etc.)
!!  spacedim= dimension of MPI communicator
!!  total_spacedim= size of global KGB communicator (typically 'banspinorfft' comm.)
!!
!! OUTPUT
!!  arraymem(2)= memory information
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

function chebfi_memInfo(neigenpairs,spacedim,space,paral_kgb,total_spacedim,bandpp) result(arraymem)

 implicit none

!Arguments ------------------------------------
 integer, intent(in   ) :: bandpp
 integer, intent(in   ) :: neigenpairs
 integer, intent(in   ) :: paral_kgb
 integer, intent(in   ) :: space
 integer, intent(in   ) :: spacedim
 integer, intent(in   ) :: total_spacedim

!Local variables-------------------------------
!scalars
 real(dp) :: memX
 real(dp) :: memX_next
 real(dp) :: memX_prev
 real(dp) :: memAX
 real(dp) :: memBX
!Transposer variables
 real(dp) :: memX_CR
 real(dp) :: memAX_CR
 real(dp) :: memBX_CR
!chebfi_rayleighRitz function variables
 real(dp) :: memA_und_X
 real(dp) :: memB_und_X
 real(dp) :: memEigenvalues
 real(dp) :: cplx
!arrays
 real(dp) :: arraymem(2)

! *********************************************************************
 cplx = 1
 if ( space == SPACE_C ) cplx = 2 !for now only complex

 !Permanent in chebfi
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
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_run
!! NAME
!! chebfi_run
!!
!! FUNCTION
!! Apply the Chebyshev Filtering algorithm on a set of vectors.
!!
!! INPUTS
!!  mpi_enreg = information about MPI parallelization
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X= pointer to the function giving B^-1|X>
!!           B is typically the overlap operator S
!!  pcond= pointer to the function used to apply the preconditioning
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine chebfi_run(chebfi,X0,getAX_BX,getBm1X,pcond,eigen,residu,mpi_enreg)

 implicit none

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi
 type(mpi_type),  intent(inout) :: mpi_enreg
 type(xgBlock_t), intent(inout) :: X0
 type(xgBlock_t), intent(inout) :: eigen
 type(xgBlock_t), intent(inout) :: residu
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
   subroutine getBm1X(X,Bm1X,transposer)
     use m_xg, only : xgBlock_t
     use m_xgTransposer !, only: xgTransposer_t
     type(xgBlock_t), intent(inout) :: X
     type(xgBlock_t), intent(inout) :: Bm1X
     type(xgTransposer_t), optional, intent(inout) :: transposer
   end subroutine getBm1X
 end interface
 interface
   subroutine pcond(W)
     use m_xg, only : xgBlock_t
     type(xgBlock_t), intent(inout) :: W
   end subroutine pcond
 end interface

!Local variables-------------------------------
!scalars
 integer :: spacedim
 integer :: neigenpairs
 integer :: nline,nline_max
 integer :: iline, iband, ierr
 integer :: nCpuCols,nCpuRows
 integer :: comm_fft_save,comm_band_save !FFT and BAND MPI communicators from rest of ABinit, to be saved
 real(dp) :: tolerance
 real(dp) :: maxeig, maxeig_global
 real(dp) :: mineig, mineig_global
 real(dp) :: lambda_minus
 real(dp) :: lambda_plus
 real(dp) :: maximum
 real(dp) :: one_over_r
 real(dp) :: two_over_r
 real(dp) :: center
 real(dp) :: radius
 type(xg_t) :: DivResults
!arrays
 real(dp) :: tsec(2)
 !Pointers similar to old Chebfi
 integer,allocatable :: nline_bands(:) !Oracle variable
 real(dp),pointer :: eig(:,:)

! *********************************************************************

 call timab(tim_run,1,tsec)

 spacedim = chebfi%spacedim
 neigenpairs = chebfi%neigenpairs
 nline = chebfi%nline
 chebfi%eigenvalues = eigen

 if (chebfi%paral_kgb == 0) then
   ABI_MALLOC(nline_bands,(neigenpairs))
   call xg_init(DivResults, chebfi%space, neigenpairs, 1)
 else
   ABI_MALLOC(nline_bands,(chebfi%bandpp))
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

!********************* Compute Rayleigh quotients for every band, and set lambda equal to the largest one *****
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RRQ)

 ! NOTICE : the following lines are kept for reference
 ! they are no longer necessary as chebfi_rayleighRitzQuotients can fully run on GPU
 ! it is no longer necessary to issue a data prefetch, data are already present on device

!  if (chebfi%use_gpu_cuda == 1) then
! #if defined(HAVE_GPU_CUDA)
!    call xgBlock_prefetch_async(chebfi%xXColsRows,  CPU_DEVICE_ID)
!    call xgBlock_prefetch_async(chebfi%xAXColsRows, CPU_DEVICE_ID)
!    call xgBlock_prefetch_async(chebfi%xBXColsRows, CPU_DEVICE_ID)
! #endif
!  end if

 call timab(tim_RR_q, 1, tsec)
 call chebfi_rayleighRitzQuotients(chebfi, maxeig, mineig, DivResults%self) !OK
 call timab(tim_RR_q, 2, tsec)

 if (chebfi%paral_kgb == 1) then
   call xmpi_max(maxeig,maxeig_global,chebfi%spacecom,ierr)
   call xmpi_min(mineig,mineig_global,chebfi%spacecom,ierr)
 else
   maxeig_global = maxeig
   mineig_global = mineig
 end if
 ABI_NVTX_END_RANGE()

 lambda_minus = maxeig_global

 nline_max = cheb_oracle1(mineig_global, lambda_minus, lambda_plus, 1D-16, 40)

 if (chebfi%paral_kgb == 0) then
   call xgBlock_reverseMap(DivResults%self,eig,1,neigenpairs)
   do iband=1, neigenpairs !TODO TODO
  ! !nline necessary to converge to tolerance
  ! !nline_tolwfr = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, tolerance / resids_filter(iband), nline)
  ! !nline necessary to decrease residual by a constant factor
  ! !nline_decrease = cheb_oracle1(dble(eig(iband*2-1,1)), lambda_minus, lambda_plus, 0.1D, dtset%nline)
  ! !nline_bands(iband) = MAX(MIN(nline_tolwfr, nline_decrease, nline_max, chebfi%nline), 1)
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

 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
 do iline = 0, nline - 1

   call timab(tim_next_p,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
   call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)
   ABI_NVTX_END_RANGE()
   call timab(tim_next_p,2,tsec)

   call timab(tim_swap,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
   if (chebfi%paral_kgb == 0) then
     call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
   else
     call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
   end if
   ABI_NVTX_END_RANGE()
   call timab(tim_swap,2,tsec)

   !A * Psi
   call timab(tim_getAX_BX,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
   call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows,chebfi%xgTransposerX)
   ABI_NVTX_END_RANGE()
   call timab(tim_getAX_BX,2,tsec)
 end do
 ABI_NVTX_END_RANGE()

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

 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
 call timab(tim_RR, 1, tsec)
 call chebfi_rayleighRitz(chebfi, nline)
 call timab(tim_RR, 2, tsec)
 ABI_NVTX_END_RANGE()

 call timab(tim_residu, 1, tsec)
 maximum =  chebfi_computeResidue(chebfi, residu, pcond)
 call timab(tim_residu, 2, tsec)

 ABI_FREE(nline_bands)

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
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_rayleighRitzQuotients
!! NAME
!! chebfi_rayleighRitzQuotients
!!
!! FUNCTION
!! Compute the Rayleigh-Ritz quotients.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!  maxeig= highest eigenvalue
!!  mineig= lowest eigenvalue
!!  DivResults= Rayleigh-Ritz quotients
!!
!! SOURCE

subroutine chebfi_rayleighRitzQuotients(chebfi,maxeig,mineig,DivResults)

 implicit none

!Arguments ------------------------------------
 real(dp), intent(inout) :: maxeig
 real(dp), intent(inout) :: mineig
 type(chebfi_t), intent(inout) :: chebfi
 type(xgBlock_t), intent(inout) :: DivResults

!Local variables-------------------------------
!scalars
 type(xg_t)::Results1
 type(xg_t)::Results2
!arrays
 integer :: maxeig_pos(2)
 integer :: mineig_pos(2)

! *********************************************************************

!Doesnt work with npfft (ncols=1 in the formula below) ???

 if (chebfi%paral_kgb == 0) then
   call xg_init(Results1, chebfi%space, chebfi%neigenpairs, 1)
   call xg_init(Results2, chebfi%space, chebfi%neigenpairs, 1)
 else
   call xg_init(Results1, chebfi%space, chebfi%bandpp, 1)
   call xg_init(Results2, chebfi%space, chebfi%bandpp, 1)
 end if

 call xgBlock_colwiseDotProduct(chebfi%xXColsRows, chebfi%xAXColsRows, Results1%self, &
   & use_gpu_cuda=chebfi%use_gpu_cuda)

!PAW
 call xgBlock_colwiseDotProduct(chebfi%xXColsRows, chebfi%xBXColsRows, Results2%self, &
   & use_gpu_cuda=chebfi%use_gpu_cuda)

 call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, &
   & maxeig, maxeig_pos, mineig, mineig_pos, &
   & use_gpu_cuda=chebfi%use_gpu_cuda)

 call xg_free(Results1)
 call xg_free(Results2)

end subroutine chebfi_rayleighRitzQuotients
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_computeNextOrderChebfiPolynom
!! NAME
!! chebfi_computeNextOrderChebfiPolynom
!!
!! FUNCTION
!! From P_n(B-^1.A)|X> (where P_n is the Chebyshev polynom of order n),
!!   computes P_n+1(B-^1.A)|X>
!!
!! INPUTS
!!  iline=order of polynom
!!  center=
!!  one_over_r,two_over_r=
!!  getBm1X= pointer to the function giving B^-1|X>
!!           B is typically the overlap operator S
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_computeNextOrderChebfiPolynom(chebfi,iline,center,one_over_r,two_over_r,getBm1X)

 implicit none

!Arguments ------------------------------------
 real(dp)       , intent(in) :: center
 integer        , intent(in) :: iline
 real(dp)       , intent(in) :: one_over_r
 real(dp)       , intent(in) :: two_over_r
 type(chebfi_t) , intent(inout) :: chebfi
 interface
   subroutine getBm1X(X,Bm1X,transposer)
     use m_xg, only : xgBlock_t
     use m_xgTransposer !, only: xgTransposer_t
     type(xgBlock_t), intent(inout) :: X
     type(xgBlock_t), intent(inout) :: Bm1X
     type(xgTransposer_t), optional, intent(inout) :: transposer
   end subroutine getBm1X
 end interface

!Local variables-------------------------------
 real(dp) :: tsec(2)

! *********************************************************************

 if (chebfi%paw) then
   call timab(tim_invovl, 1, tsec)
   ABI_NVTX_START_RANGE(NVTX_INVOVL)
   call getBm1X(chebfi%xAXColsRows, chebfi%X_next, chebfi%xgTransposerX)
   ABI_NVTX_END_RANGE()
   call timab(tim_invovl, 2, tsec)
 else
   call xgBlock_copy(chebfi%xAXColsRows,chebfi%X_next, 1, 1)
 end if

 call xgBlock_scale(chebfi%xXColsRows, center, 1) !scale by c

 !(B-1 * A * Psi^i-1 - c * Psi^i-1)
 call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%xXColsRows)

 !Psi^i-1  = 1/c * Psi^i-1
 call xgBlock_scale(chebfi%xXColsRows, 1/center, 1) !counter scale by c

 if (iline == 0) then
   call xgBlock_scale(chebfi%X_next, one_over_r, 1)
else
   call xgBlock_scale(chebfi%X_next, two_over_r, 1)

   call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X_prev)
end if

end subroutine chebfi_computeNextOrderChebfiPolynom
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_swapInnerBuffers
!! NAME
!! chebfi_swapInnerBuffers
!!
!! FUNCTION
!! Swap buffers inside a 'chebfi' datastructure.
!!
!! INPUTS
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  spacedim= space dimension for one vector
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_swapInnerBuffers(chebfi,spacedim,neigenpairs)

 implicit none

!Arguments ------------------------------------
 integer        , intent(in) :: spacedim
 integer        , intent(in) :: neigenpairs
 type(chebfi_t) , intent(inout) :: chebfi

! *********************************************************************

 call xgBlock_setBlock(chebfi%X_prev, chebfi%X_swap, 1, spacedim, neigenpairs) !X_swap = X_prev
 call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X_prev, 1, spacedim, neigenpairs) !X_prev = xXColsRows
 call xgBlock_setBlock(chebfi%X_next, chebfi%xXColsRows, 1, spacedim, neigenpairs) !xXColsRows = X_next
 call xgBlock_setBlock(chebfi%X_swap, chebfi%X_next, 1, spacedim, neigenpairs) !X_next = X_swap

end subroutine chebfi_swapInnerBuffers
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_rayleighRitz
!! NAME
!! chebfi_rayleighRitz
!!
!! FUNCTION
!! Apply Rayleigh-Ritz procedure
!!
!! INPUTS
!!  nline= order of Chebyshev polynom
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_rayleighRitz(chebfi,nline)

 implicit none

!Arguments ------------------------------------
 integer        , intent(in) :: nline
 type(chebfi_t) , intent(inout) :: chebfi
!Local variables-------------------------------
!scalars
 integer :: eigenProblem
 integer :: info
 integer :: me_g0
 integer :: neigenpairs
 integer :: space
 integer :: spacedim
 integer :: remainder
 type(xg_t) :: A_und_X !H_UND_PSI
 type(xg_t) :: B_und_X !S_UND_PSI
 type(xgBlock_t) :: X_first_row
 type(xgBlock_t) :: AX_first_row
 type(xgBlock_t) :: BX_first_row
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

 space = chebfi%space
 eigenProblem = chebfi%eigenProblem
 me_g0 = chebfi%me_g0

 spacedim = chebfi%spacedim
 neigenpairs = chebfi%neigenpairs   !remains whole nband domain since it is after transpose

 call xg_init(A_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)
 call xg_init(B_und_X,space,neigenpairs,neigenpairs,chebfi%spacecom)

 call timab(tim_RR_gemm_1,1,tsec)
 call xgBlock_gemm(chebfi%AX%self%trans, chebfi%X%normal, 1.0d0, chebfi%AX%self, chebfi%X, 0.d0, A_und_X%self)

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
 select case (eigenSolver)
 case (EIGENVD)
   call xgBlock_hegvd(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
 case (EIGENV)
   call xgBlock_hegv(eigenProblem, 'v','u', A_und_X%self, B_und_X%self, chebfi%eigenvalues, info)
 case default
    ABI_ERROR("Error for Eigen Solver HEGV")
 end select
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
&                     chebfi%X, A_und_X%self, 0.d0, chebfi%X_next)
   call xgBlock_copy(chebfi%X_next,chebfi%X_swap, 1, 1)    !copy cannot be avoided :(
 end if

 call xgBlock_gemm(chebfi%AX%self%normal, A_und_X%self%normal, 1.0d0, &
                   chebfi%AX%self, A_und_X%self, 0.d0, chebfi%AX_swap)

 if (chebfi%paw) then
   call xgBlock_gemm(chebfi%BX%self%normal, A_und_X%self%normal, 1.0d0, &
                     chebfi%BX%self, A_und_X%self, 0.d0, chebfi%BX_swap)
 end if

 call timab(tim_RR_gemm_2,2,tsec)

 call xg_free(A_und_X)
 call xg_free(B_und_X)

end subroutine chebfi_rayleighRitz
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_computeResidue
!! NAME
!! chebfi_computeResidue
!!
!! FUNCTION
!! Compute vector of residuals (AX-e.BX), including predconditionning
!!
!! INPUTS
!!  pcond= pointer to the function used to apply the preconditioning
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  residu<type(xgBlock_t)>= vector of residuals
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

 real(dp) function chebfi_computeResidue(chebfi, residu, pcond)

 implicit none

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi
 type(xgBlock_t) , intent(inout) :: residu
 interface
   subroutine pcond(W)
     use m_xg, only : xgBlock_t
     type(xgBlock_t), intent(inout) :: W
   end subroutine pcond
 end interface

!Local variables-------------------------------
!scalars
 integer :: eigResiduMax,eigResiduMin
 real(dp) :: maxResidu,minResidu
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

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
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_ampfactor
!! NAME
!! chebfi_ampfactor
!!
!! FUNCTION
!! Compute amplification factor
!!
!! INPUTS
!! eig (:,:)= eigenvalues
!! lambda_minus,lambda_plus=
!! nline_bands(:)= order of Chebyshev polynom for each band
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  residu<type(xgBlock_t)>= vector of residuals
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_ampfactor(chebfi,eig,lambda_minus,lambda_plus,nline_bands)

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nline_bands(:)
 real(dp), pointer, intent(in) :: eig(:,:)
 real(dp), intent(in) :: lambda_minus
 real(dp), intent(in) :: lambda_plus
 type(chebfi_t) , intent(inout) :: chebfi

!Local variables-------------------------------
!scalars
 integer :: iband,nbands
 real(dp) :: ampfactor
 real(dp) :: eig_per_band
 type(xgBlock_t) :: X_part
 type(xgBlock_t) :: AX_part
 type(xgBlock_t) :: BX_part

! *********************************************************************

 if (chebfi%paral_kgb == 0) then
   nbands = chebfi%neigenpairs
 else
   nbands = chebfi%bandpp
 end if

 do iband = 1, nbands

   eig_per_band = eig(1,iband)

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
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_oracle1
!! NAME
!! chebfi_oracle1
!!
!! FUNCTION
!! Compute order of Chebyshev polynom necessary to converge to a given tol
!!
!! INPUTS
!!  xx= input variable
!!  aa= left bound of the interval
!!  bb= right bound of the interval
!!  tol= needed precision
!!  nmax= max number of iterations
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

function cheb_oracle1(xx,aa,bb,tol,nmax) result(nn)

 implicit none

!Arguments ------------------------------------
 integer :: nn
 integer,intent(in) :: nmax
 real(dp), intent(in) :: xx,aa,bb
 real(dp) :: tol

!Local variables-------------------------------
 integer :: ii
 real(dp) :: yy,yim1,xred,temp

! *************************************************************************

 xred = (xx-(aa+bb)/2)/(bb-aa)*2
 yy = xred
 yim1 = 1 !ONE

 nn = nmax
 if(1/(yy**2) < tol) then
   nn = 1
 else
   do ii=2, nmax-1
     temp = yy
     yy = 2*xred*yy - yim1
     yim1 = temp
     if(1/(yy**2) < tol) then
       nn = ii
       exit
     end if
   end do
 end if

end function cheb_oracle1
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_poly1
!! NAME
!! chebfi_poly1
!!
!! FUNCTION
!! Compute Chebyshev polynomial???
!!
!! INPUTS
!!  xx= input variable
!!  aa= left bound of the interval
!!  bb= right bound of the interval
!!  nn=
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

function cheb_poly1(xx,nn,aa,bb) result(yy)

 implicit none

!Arguments ------------------------------------
 integer, intent(in) :: nn
 real(dp), intent(in) :: xx, aa, bb
 real(dp) :: yy

!Local variables-------------------------------
 integer :: ii
 real(dp) :: xred,yim1,temp

! *************************************************************************

 xred = (xx-(aa+bb)/2)/(bb-aa)*2
 yy = xred
 yim1 = 1
 do ii= 2, nn
   temp = yy
   yy = 2*xred*yy - yim1
   yim1 = temp
 end do

end function cheb_poly1
!!***

end module m_chebfi2
!!***

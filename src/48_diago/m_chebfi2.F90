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
!! Copyright (C) 2018-2026 ABINIT group (BS, L. Baguet)
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
 use m_xg_ortho_RR

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif
 use, intrinsic :: iso_c_binding, only: c_size_t

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU_MARKERS)
 use m_nvtx_data
#endif

 implicit none

 private

!Several (private) parameters
!-------------------------------------------------

 integer, parameter :: tim_init         = 1751
 integer, parameter :: tim_free         = 1752
 !                                        1753 is used by chebfi2_nonlop
 integer, parameter :: tim_getAX_BX     = 1754
 integer, parameter :: tim_invovl       = 1755
 integer, parameter :: tim_residu       = 1756
 integer, parameter :: tim_RR           = 1757
 integer, parameter :: tim_transpose    = 1758
 integer, parameter :: tim_RR_q         = 1759
 integer, parameter :: tim_postinvovl   = 1760
 integer, parameter :: tim_swap         = 1761
 integer, parameter :: tim_amp_f        = 1762
 integer, parameter :: tim_oracle       = 1763
 integer, parameter :: tim_barrier      = 1764
 integer, parameter :: tim_copy         = 1765

!Public 'chebfi' datatype
!-------------------------------------------------

 type, public :: chebfi_t
   integer :: space
   integer :: spacedim                      ! Space dimension for one vector
   integer :: total_spacedim                ! Maybe not needed
   integer :: neigenpairs                   ! Number of eigen values/vectors we want
   integer :: ndeg_filter                   ! Degree of the polynomial filter
   integer :: nbdbuf                        ! Number of bands in the buffer
   integer :: spacecom                      ! Communicator for MPI
   integer :: oracle                        ! Option to compute ndeg_filter from residuals
   real(dp) :: tolerance            ! Tolerance on the residu to stop the minimization
   real(dp) :: ecut                 ! Ecut for Chebfi oracle
   real(dp) :: oracle_factor                ! factor used to decrease residuals
   real(dp) :: oracle_min_occ               ! threshold on occupancies used for nbdbuf=-101

   integer :: paral_kgb                     ! MPI parallelization variables
   integer :: bandpp
   integer :: comm_cols
   integer :: comm_rows
   integer :: me_g0
   integer :: me_g0_fft

   logical :: paw
   integer :: eigenProblem   !1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)

   ! when GPU with Kokkos is enabled, currently OpenMP is not fully supported, abinit is launched
   ! with OMP_NUM_THREADS=1, but we may locally increase the number of OpenMP threads
   ! wherever it is safe to do; in that case we use gpu_kokkos_nthrd to specify
   ! the number of OpenMP threads. This value is controlled by dtset variable
   ! dtset%gpu_kokkos_nthrd
   integer :: gpu_option
   integer :: gpu_kokkos_nthrd = 1 ! only used if GPU Kokkos is enabled, number of OpenMP threads used
   integer :: gpu_thread_limit = 1 ! only used if GPU is enabled, max number of OpenMP threads used in sensitive areas

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
!!  me_g0= 1 if this processors treats G=0, 0 otherwise
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  ndeg_filter= polynomial degree of the Chebyshev filter (.i.e. number of H applications)
!!  comm_rows= "rows" communicator
!!  comm_cols= "cols" communicator
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

subroutine chebfi_init(chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp, &
                       ndeg_filter,nbdbuf,space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols, &
                       oracle,oracle_factor,oracle_min_occ,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit)

 ! Arguments ------------------------------------
 integer       , intent(in   ) :: bandpp
 integer       , intent(in   ) :: eigenProblem
 integer       , intent(in   ) :: me_g0
 integer       , intent(in   ) :: me_g0_fft
 integer       , intent(in   ) :: neigenpairs
 integer       , intent(in   ) :: ndeg_filter
 integer       , intent(in   ) :: nbdbuf
 integer       , intent(in   ) :: comm_cols
 integer       , intent(in   ) :: comm_rows
 integer       , intent(in   ) :: paral_kgb
 integer       , intent(in   ) :: space
 integer       , intent(in   ) :: spacecom
 integer       , intent(in   ) :: spacedim
 integer       , intent(in   ) :: gpu_option
 integer       , intent(in   ) :: oracle
 logical       , intent(in   ) :: paw
 real(dp)      , intent(in   ) :: ecut
 real(dp)      , intent(in   ) :: tolerance
 real(dp)      , intent(in   ) :: oracle_factor
 real(dp)      , intent(in   ) :: oracle_min_occ
 type(chebfi_t), intent(inout) :: chebfi
 integer       , intent(in   ), optional :: gpu_kokkos_nthrd
 integer       , intent(in   ), optional :: gpu_thread_limit

 ! Local variables-------------------------------
 real(dp)                      :: tsec(2)

 ! *********************************************************************

 call timab(tim_init,1,tsec)

 chebfi%space = space
 chebfi%neigenpairs = neigenpairs
 chebfi%spacedim    = spacedim
 if (tolerance > 0.0) then
   chebfi%tolerance = tolerance
 else
   chebfi%tolerance = 1.0e-20
 end if
 chebfi%ecut        = ecut
 chebfi%paral_kgb   = paral_kgb
 chebfi%comm_cols   = comm_cols
 chebfi%bandpp      = bandpp
 chebfi%comm_rows   = comm_rows
 chebfi%ndeg_filter = ndeg_filter
 chebfi%nbdbuf      = nbdbuf
 chebfi%spacecom    = spacecom
 chebfi%eigenProblem = eigenProblem
 chebfi%me_g0        = me_g0
 chebfi%me_g0_fft    = me_g0_fft
 chebfi%paw          = paw
 chebfi%gpu_option  = gpu_option
 chebfi%oracle      = oracle
 chebfi%oracle_factor = oracle_factor
 chebfi%oracle_min_occ = oracle_min_occ

 chebfi%gpu_kokkos_nthrd = 1
 if (present(gpu_kokkos_nthrd)) chebfi%gpu_kokkos_nthrd = gpu_kokkos_nthrd
 chebfi%gpu_thread_limit = 0
 if (present(gpu_thread_limit)) chebfi%gpu_thread_limit = gpu_thread_limit

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

 ! Arguments ------------------------------------
 type(chebfi_t)  , intent(inout) :: chebfi

 ! Local variables-------------------------------
 ! scalars
 integer  :: neigenpairs
 integer  :: space
 integer  :: spacedim
 integer  :: total_spacedim, ierr

! *********************************************************************

 space       = chebfi%space
 spacedim    = chebfi%spacedim
 neigenpairs = chebfi%neigenpairs

 call chebfi_free(chebfi)

 if (chebfi%paral_kgb == 0) then
   chebfi%total_spacedim = spacedim
   call xg_init(chebfi%X_NP,space,spacedim,2*neigenpairs,chebfi%spacecom,me_g0=chebfi%me_g0,gpu_option=chebfi%gpu_option) !regular arrays
   call xg_setBlock(chebfi%X_NP, chebfi%X_next,spacedim, neigenpairs)
   call xg_setBlock(chebfi%X_NP, chebfi%X_prev,spacedim, neigenpairs, fcol=neigenpairs+1)
 else
   total_spacedim = spacedim
   call xmpi_sum(total_spacedim,chebfi%spacecom,ierr)
   chebfi%total_spacedim = total_spacedim
   call xg_init(chebfi%X_NP,space,total_spacedim,2*chebfi%bandpp,chebfi%spacecom,me_g0=chebfi%me_g0_fft,&
     & gpu_option=chebfi%gpu_option) !transposed arrays
   call xg_setBlock(chebfi%X_NP, chebfi%X_next, total_spacedim, chebfi%bandpp)
   call xg_setBlock(chebfi%X_NP, chebfi%X_prev, total_spacedim, chebfi%bandpp, fcol=chebfi%bandpp+1)
 end if

 !transposer will handle these arrays automatically
 call xg_init(chebfi%AX,space,spacedim,neigenpairs,chebfi%spacecom,me_g0=chebfi%me_g0,gpu_option=chebfi%gpu_option)
 call xg_init(chebfi%BX,space,spacedim,neigenpairs,chebfi%spacecom,me_g0=chebfi%me_g0,gpu_option=chebfi%gpu_option)

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

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi

! *********************************************************************

 call xg_free(chebfi%X_NP)

 call xg_free(chebfi%AX)
 call xg_free(chebfi%BX)

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

!Arguments ------------------------------------
 integer, intent(in   ) :: bandpp
 integer, intent(in   ) :: neigenpairs
 integer, intent(in   ) :: paral_kgb
 integer, intent(in   ) :: space
 integer, intent(in   ) :: spacedim
 integer, intent(in   ) :: total_spacedim

!Local variables-------------------------------
!scalars
 integer(kind=c_size_t) :: memX
 integer(kind=c_size_t) :: memX_next
 integer(kind=c_size_t) :: memX_prev
 integer(kind=c_size_t) :: memAX
 integer(kind=c_size_t) :: memBX
!Transposer variables
 integer(kind=c_size_t) :: memX_CR
 integer(kind=c_size_t) :: memAX_CR
 integer(kind=c_size_t) :: memBX_CR
 integer(kind=c_size_t) :: mem_sendrecv_CR
!chebfi_rayleighRitz function variables
 integer(kind=c_size_t) :: memA_und_X
 integer(kind=c_size_t) :: memB_und_X
 integer(kind=c_size_t) :: memEigenvalues
 integer(kind=c_size_t) :: cplx
!arrays
 integer(kind=c_size_t) :: arraymem(2)

! *********************************************************************
 cplx = 1
 if ( space == SPACE_C ) cplx = 2 !for now only complex

 !Permanent in chebfi
 memX = int(cplx,c_size_t) * kind(1.d0) * spacedim * neigenpairs

 if (paral_kgb == 0) then
   memX_next = int(cplx,c_size_t) * kind(1.d0) * spacedim * neigenpairs
   memX_prev = int(cplx,c_size_t) * kind(1.d0) * spacedim * neigenpairs
 else
   memX_next = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
   memX_prev = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
 end if

 memAX = int(cplx,c_size_t) * kind(1.d0) * spacedim * neigenpairs
 memBX = int(cplx,c_size_t) * kind(1.d0) * spacedim * neigenpairs

 !Transposer colrow array
 if (paral_kgb == 1) then
   memX_CR = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
   memAX_CR = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
   memBX_CR = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
   mem_sendrecv_CR = int(cplx,c_size_t) * kind(1.d0) * total_spacedim * bandpp
 else
   memX_CR = 0
   memAX_CR = 0
   memBX_CR = 0
 end if

 !chebfi_rayleighRitz function variables
 memA_und_X = int(cplx,c_size_t) * kind(1.d0) * neigenpairs * neigenpairs
 memB_und_X = int(cplx,c_size_t) * kind(1.d0) * neigenpairs * neigenpairs
 memEigenvalues = int(kind(1.d0),c_size_t) * neigenpairs

 arraymem(1) = memX + memX_next + memX_prev + &
               memAX + memBX + memX_CR + memAX_CR + memBX_CR + mem_sendrecv_CR
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

subroutine chebfi_run(chebfi,X0,getAX_BX,getBm1X,eigen,occ,residu,nspinor)

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi
 integer,         intent(in)    :: nspinor
 type(xgBlock_t), intent(inout) :: X0
 type(xgBlock_t), intent(inout) :: eigen
 type(xgBlock_t), intent(in)    :: occ
 type(xgBlock_t), intent(inout) :: residu
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

!Local variables-------------------------------
!scalars
 integer :: spacedim
 integer :: space_res
 integer :: neigenpairs,bandpp
 integer :: ndeg_filter,ndeg_filter_max
 integer :: ideg, ierr
 real(dp) :: tolerance
 real(dp) :: maxeig, maxeig_global
 real(dp) :: mineig, mineig_global
 real(dp) :: lambda_minus
 real(dp) :: lambda_plus
 real(dp) :: one_over_r
 real(dp) :: two_over_r
 real(dp) :: center
 real(dp) :: radius
 type(xg_t) :: DivResults
!arrays
 real(dp) :: tsec(2)
 !Pointers similar to old Chebfi
 integer,allocatable :: ndeg_filter_bands(:) !Oracle variable

! *********************************************************************

! call timab(tim_run,1,tsec)

 spacedim = chebfi%spacedim
 neigenpairs = chebfi%neigenpairs
 bandpp = chebfi%bandpp
 ndeg_filter = chebfi%ndeg_filter
 chebfi%eigenvalues = eigen

 if (chebfi%space==SPACE_C) then
   space_res = SPACE_C
 else if (chebfi%space==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if

 if (chebfi%paral_kgb == 0) then
   ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
   call xg_init(DivResults, space_res, neigenpairs, 1, gpu_option=chebfi%gpu_option)
 else
   ABI_MALLOC(ndeg_filter_bands,(bandpp))
   call xg_init(DivResults, space_res, bandpp, 1, gpu_option=chebfi%gpu_option)
 end if

 tolerance = chebfi%tolerance
 lambda_plus = chebfi%ecut
 chebfi%X = X0

 ! Transpose
 if (chebfi%paral_kgb == 1) then

   call timab(tim_transpose,1,tsec)
   call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
     STATE_LINALG,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0_fft,&
     gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit)

   call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,chebfi%AX%self,chebfi%xAXColsRows,STATE_LINALG)
   call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,chebfi%BX%self,chebfi%xBXColsRows,STATE_LINALG)

   chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
   chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
   chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
   call xgTransposer_transpose(chebfi%xgTransposerX,STATE_COLSROWS)
   chebfi%xgTransposerAX%state = STATE_COLSROWS
   chebfi%xgTransposerBX%state = STATE_COLSROWS
   ABI_NVTX_END_RANGE()
   call timab(tim_transpose,2,tsec)
 else
   call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)   !use xXColsRows instead of X notion
   call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)   !use xAXColsRows instead of AX notion
   call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
 end if

 call timab(tim_getAX_BX,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
 call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
 call xgBlock_zero_im_g0(chebfi%xAXColsRows)
 call xgBlock_zero_im_g0(chebfi%xBXColsRows)
 ABI_NVTX_END_RANGE()
 call timab(tim_getAX_BX,2,tsec)

 if (chebfi%paral_kgb == 1) then
   call timab(tim_barrier,1,tsec)
   call xmpi_barrier(chebfi%spacecom)
   call timab(tim_barrier,2,tsec)
 end if

!********************* Compute Rayleigh quotients for every band, and set lambda equal to the largest one *****
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RRQ)

 ! NOTICE : the following lines are kept for reference
 ! they are no longer necessary as chebfi_rayleighRitzQuotients can fully run on GPU
 ! it is no longer necessary to issue a data prefetch, data are already present on device

!  if (chebfi%gpu_option == ABI_GPU_KOKKOS) then
! #if defined(HAVE_GPU_CUDA)
!    call xgBlock_prefetch_async(chebfi%xXColsRows,  CPU_DEVICE_ID)
!    call xgBlock_prefetch_async(chebfi%xAXColsRows, CPU_DEVICE_ID)
!    call xgBlock_prefetch_async(chebfi%xBXColsRows, CPU_DEVICE_ID)
! #endif
!  end if

 call timab(tim_RR_q, 1, tsec)
 call chebfi_rayleighRitzQuotients(chebfi, maxeig, mineig, DivResults%self)

 if (chebfi%paral_kgb == 1) then
   call xmpi_max(maxeig,maxeig_global,chebfi%spacecom,ierr)
   call xmpi_min(mineig,mineig_global,chebfi%spacecom,ierr)
 else
   maxeig_global = maxeig
   mineig_global = mineig
 end if
 call timab(tim_RR_q, 2, tsec)
 ABI_NVTX_END_RANGE()

 lambda_minus = maxeig_global

 call timab(tim_oracle,1,tsec)

 ! ndeg_filter_max limits the reduction of the residual of the smallest eigenvalue (i.e. the most amplified one by the filter) by a factor 1e8.
 ! Also, the maximal value of ndeg_filter_max is 40.
 ndeg_filter_max = cheb_oracle1(mineig_global, lambda_minus, lambda_plus, 1D-16, 40)
 ndeg_filter = MIN(ndeg_filter_max,chebfi%ndeg_filter)
 if (chebfi%oracle>0) then
   call chebfi_set_ndeg_from_residu(chebfi,lambda_minus,lambda_plus,occ,DivResults%self,ndeg_filter_max,ndeg_filter)
 end if
 ndeg_filter_bands(:) = ndeg_filter

 call timab(tim_oracle,2,tsec)

 center = (lambda_plus + lambda_minus)*0.5
 radius = (lambda_plus - lambda_minus)*0.5

 one_over_r = 1/radius
 two_over_r = 2/radius

 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
 do ideg = 0, ndeg_filter - 1

   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
   call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
   ABI_NVTX_END_RANGE()

   call timab(tim_swap,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
   if (chebfi%paral_kgb == 0) then
     call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
   else
     call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, bandpp)
   end if
   ABI_NVTX_END_RANGE()
   call timab(tim_swap,2,tsec)

   !A * Psi
   call timab(tim_getAX_BX,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
   call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
   call xgBlock_zero_im_g0(chebfi%xAXColsRows)
   call xgBlock_zero_im_g0(chebfi%xBXColsRows)
   ABI_NVTX_END_RANGE()
   call timab(tim_getAX_BX,2,tsec)

 end do ! ideg
 ABI_NVTX_END_RANGE()

 if (chebfi%paral_kgb == 1) then
   call timab(tim_barrier,1,tsec)
   call xmpi_barrier(chebfi%spacecom)
   call timab(tim_barrier,2,tsec)
 end if

 call timab(tim_amp_f,1,tsec)
 call chebfi_ampfactor(chebfi, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)
 call timab(tim_amp_f,2,tsec)

 call xg_free(DivResults)
 ABI_FREE(ndeg_filter_bands)

 call timab(tim_transpose,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
 if (chebfi%paral_kgb == 1) then
   call xmpi_barrier(chebfi%spacecom)

   call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG)

   if (xmpi_comm_size(chebfi%spacecom) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
     call xgBlock_setBlock(chebfi%xXColsRows,  chebfi%X,       spacedim, neigenpairs)
     call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
     call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
   end if
 else
   call xgBlock_setBlock(chebfi%xXColsRows,  chebfi%X,       spacedim, neigenpairs)
   call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
   call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
 end if
 ABI_NVTX_END_RANGE()
 call timab(tim_transpose,2,tsec)

 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
 call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,chebfi%eigenvalues,ierr,0,tim_RR,chebfi%gpu_option,&
 &    solve_ax_bx=.true.)
 ABI_NVTX_END_RANGE()

 call timab(tim_residu, 1, tsec)
 if (chebfi%paw) then
   call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
 else
   call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
 end if

 call xgBlock_colwiseNorm2(chebfi%AX%self, residu)
 call timab(tim_residu, 2, tsec)

 call timab(tim_copy, 1, tsec)
 call xgBlock_copy(chebfi%X,X0)
 call timab(tim_copy, 2, tsec)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
   if (chebfi%gpu_option==ABI_GPU_KOKKOS) then
     call gpu_device_synchronize()
   end if
#endif

 if (chebfi%paral_kgb == 1) then
   call xgTransposer_free(chebfi%xgTransposerX)
   call xgTransposer_free(chebfi%xgTransposerAX)
   call xgTransposer_free(chebfi%xgTransposerBX)
 end if

! call timab(tim_run,2,tsec)

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
 integer :: space_res

! *********************************************************************

 if (chebfi%space==SPACE_C) then
   space_res = SPACE_C
 else if (chebfi%space==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if

!Doesnt work with npfft (ncols=1 in the formula below) ???
 if (chebfi%paral_kgb == 0) then
   call xg_init(Results1, space_res, chebfi%neigenpairs, 1, gpu_option=chebfi%gpu_option)
   call xg_init(Results2, space_res, chebfi%neigenpairs, 1, gpu_option=chebfi%gpu_option)
 else
   call xg_init(Results1, space_res, chebfi%bandpp, 1, gpu_option=chebfi%gpu_option)
   call xg_init(Results2, space_res, chebfi%bandpp, 1, gpu_option=chebfi%gpu_option)
 end if

 ! <Psi|H|Psi>
 call xgBlock_colwiseDotProduct(chebfi%xXColsRows, chebfi%xAXColsRows, Results1%self, comm_loc=xmpi_comm_null)

 ! <Psi|S|Psi>
 call xgBlock_colwiseDotProduct(chebfi%xXColsRows, chebfi%xBXColsRows, Results2%self, comm_loc=xmpi_comm_null)

 ! eig = <Psi|H|Psi> / <Psi|S|Psi>
 call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, &
   & maxeig, maxeig_pos, mineig, mineig_pos)

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
!!  ideg=current degree of polynom
!!  center=filter center
!!  one_over_r,two_over_r=1/R, 2/R, R being the radius of the filter
!!  getBm1X= pointer to the function giving B^-1|X>
!!           B is typically the overlap operator S
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_computeNextOrderChebfiPolynom(chebfi,ideg,center,one_over_r,two_over_r,getBm1X)

!Arguments ------------------------------------
 real(dp)       , intent(in) :: center
 integer        , intent(in) :: ideg
 real(dp)       , intent(in) :: one_over_r
 real(dp)       , intent(in) :: two_over_r
 type(chebfi_t) , intent(inout) :: chebfi
 interface
   subroutine getBm1X(X,Bm1X)
     use m_xg, only : xgBlock_t
     type(xgBlock_t), intent(inout) :: X
     type(xgBlock_t), intent(inout) :: Bm1X
   end subroutine getBm1X
 end interface

 !Local variables-------------------------------
 real(dp) :: tsec(2)

 ! *********************************************************************

 if (chebfi%paw) then
   call timab(tim_invovl, 1, tsec)
   ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_BM1X)
   call getBm1X(chebfi%xAXColsRows, chebfi%X_next)
   ABI_NVTX_END_RANGE()
   call timab(tim_invovl, 2, tsec)
 else
   call timab(tim_copy, 1, tsec)
   call xgBlock_copy(chebfi%xAXColsRows,chebfi%X_next)
   call timab(tim_copy, 2, tsec)
 end if

 call timab(tim_postinvovl, 1, tsec)
 ABI_NVTX_START_RANGE(NVTX_INVOVL_POST3)
 call xgBlock_scale(chebfi%xXColsRows, center, 1) !scale by center

 !(B-1 * A * Psi^i-1 - c * Psi^i-1)
 call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%xXColsRows)

 !Psi^i-1  = 1/c * Psi^i-1
 call xgBlock_scale(chebfi%xXColsRows, 1/center, 1) !counter scale by 1/center

 if (ideg == 0) then
   call xgBlock_scale(chebfi%X_next, one_over_r, 1)
 else
   call xgBlock_scale(chebfi%X_next, two_over_r, 1)

   call xgBlock_saxpy(chebfi%X_next, dble(-1.0), chebfi%X_prev)
 end if

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
 if (chebfi%gpu_option==ABI_GPU_KOKKOS) then
   call gpu_device_synchronize()
 end if
#endif
 ABI_NVTX_END_RANGE()
 call timab(tim_postinvovl, 2, tsec)

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

  ! Arguments ------------------------------------
  integer        , intent(in   ) :: spacedim
  integer        , intent(in   ) :: neigenpairs
  type(chebfi_t) , intent(inout) :: chebfi

  ! *********************************************************************

  call xgBlock_setBlock(chebfi%X_prev,     chebfi%X_swap,     spacedim, neigenpairs) !X_swap = X_prev
  call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X_prev,     spacedim, neigenpairs) !X_prev = xXColsRows
  call xgBlock_setBlock(chebfi%X_next,     chebfi%xXColsRows, spacedim, neigenpairs) !xXColsRows = X_next
  call xgBlock_setBlock(chebfi%X_swap,     chebfi%X_next,     spacedim, neigenpairs) !X_next = X_swap

end subroutine chebfi_swapInnerBuffers
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
!! ndeg_filter_bands(:)= degree of Chebyshev polynomial filter for each band
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  residu<type(xgBlock_t)>= vector of residuals
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_ampfactor(chebfi,DivResults,lambda_minus,lambda_plus,ndeg_filter_bands)

  ! Arguments ------------------------------------
  integer,           intent(in   ) :: ndeg_filter_bands(:)
  type(xgBlock_t),   intent(in   ) :: DivResults
  real(dp),          intent(in   ) :: lambda_minus
  real(dp),          intent(in   ) :: lambda_plus
  type(chebfi_t),    intent(inout) :: chebfi

  ! Local variables-------------------------------
  ! scalars
  integer         :: iband,nbands
  real(dp)        :: ampfactor
  real(dp)        :: eig_per_band
  type(xgBlock_t) :: X_part
  type(xgBlock_t) :: AX_part
  type(xgBlock_t) :: BX_part
  real(dp),pointer :: eig(:,:)

  ! *********************************************************************

  if (chebfi%paral_kgb == 0) then
    nbands = chebfi%neigenpairs
  else
    nbands = chebfi%bandpp
  end if

  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=chebfi%bandpp)

  do iband = 1, nbands

    eig_per_band = eig(1,iband)

    !cheb_poly1(x, n, a, b)
    ampfactor = cheb_poly1(eig_per_band, ndeg_filter_bands(iband), lambda_minus, lambda_plus)

    if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much

    call xgBlock_setBlock(chebfi%xXColsRows, X_part, chebfi%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(chebfi%xAXColsRows, AX_part, chebfi%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(chebfi%xBXColsRows, BX_part, chebfi%total_spacedim, 1, fcol=iband)

    call xgBlock_scale(X_part, 1/ampfactor, 1)
    call xgBlock_scale(AX_part, 1/ampfactor, 1)
    call xgBlock_scale(BX_part, 1/ampfactor, 1)

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

  ! Arguments ------------------------------------
  integer              :: nn
  integer,  intent(in) :: nmax
  real(dp), intent(in) :: xx,aa,bb
  real(dp), intent(in) :: tol

  ! Local variables-------------------------------
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

  ! Arguments ------------------------------------
  integer,  intent(in) :: nn
  real(dp), intent(in) :: xx, aa, bb
  real(dp)             :: yy

  ! Local variables-------------------------------
  integer  :: ii
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

!!****f* m_chebfi2/chebfi_set_ndeg_from_residu
!! NAME
!! chebfi_set_ndeg_from_residu
!!
!! FUNCTION
!! Compute ndeg_filter using the oracle and residuals.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! SOURCE

subroutine chebfi_set_ndeg_from_residu(chebfi,lambda_minus,lambda_plus,occ,DivResults,ndeg_filter_max,ndeg_filter)

 integer,intent(in) :: ndeg_filter_max
 integer,intent(out) :: ndeg_filter
 type(chebfi_t), intent(inout) :: chebfi
 type(xgBlock_t), intent(in)    :: occ
 type(xgBlock_t), intent(in)    :: DivResults
 real(dp), intent(in) :: lambda_minus, lambda_plus

 logical :: test1,test2,test3
 integer :: iband_tot,iband
 integer :: bandpp,ierr,ndeg_filter_tolwfr,ndeg_filter_decrease,nbdbuf,ndeg_filter_all,shift
 integer,allocatable :: ndeg_filter_bands(:)
 type(xgBlock_t) :: occBlock,occ_reshaped
 type(xg_t) :: residu
 real(dp),pointer :: residu_(:,:),occ_(:,:)
 real(dp) :: eig_iband,res_iband,occ_iband
 real(dp),pointer :: eig(:,:)
! character(len=500) :: msg

 bandpp = chebfi%bandpp

 !Compute residu here for oracle, use X_next as a work space
 ! X_next = S|Psi>
 call xgBlock_copy(chebfi%xBXColsRows,chebfi%X_next)
 ! X_next = - eig * S|Psi>
 call xgBlock_ymax(chebfi%X_next,DivResults,0,1)
 ! X_next = H|Psi> - eig * S|Psi>
 call xgBlock_add(chebfi%X_next,chebfi%xAXColsRows)
 ! resid = |X_next|^2
 call xg_init(residu,SPACE_R,bandpp,1)
 call xgBlock_colwiseNorm2(chebfi%X_next, residu%self,comm_loc=xmpi_comm_null)

 occ_reshaped = occ
 shift=xmpi_comm_rank(chebfi%comm_cols)*bandpp
 call xgBlock_reshape(occ_reshaped,1,chebfi%neigenpairs)
 call xgBlock_setBlock(occ_reshaped,occBlock,1,bandpp,fcol=1+shift)
 call xgBlock_reshape(occBlock,bandpp,1)
 if (chebfi%nbdbuf==-101) then
   call xgBlock_apply_diag(residu%self,occBlock,1)
 end if

 ABI_MALLOC(ndeg_filter_bands,(bandpp))

 ! DivResults could be complex (with null imaginary part), so bandpp has to be in cols, not rows
 call xgBlock_reverseMap(DivResults,eig,rows=1,cols=bandpp)
 call xgBlock_reverseMap(residu%self,residu_,rows=1,cols=bandpp)
 call xgBlock_reverseMap(occBlock,occ_,rows=1,cols=bandpp)

 if (chebfi%nbdbuf>0) then
   nbdbuf = chebfi%nbdbuf
 else if (chebfi%nbdbuf==-101) then
   nbdbuf = 0
 end if

 do iband=1, bandpp
   eig_iband = eig(1,iband)
   res_iband = residu_(1,iband)
   occ_iband = occ_(1,iband)
   iband_tot = iband + shift
   test1 = res_iband<chebfi%tolerance ! band already converged
   test2 = iband_tot>chebfi%neigenpairs-nbdbuf ! band in the buffer
   test3 = chebfi%nbdbuf==-101.and.occ_iband<chebfi%oracle_min_occ ! occupancy is too low
   if (test1.or.test2.or.test3) then
     ndeg_filter_bands(iband) = 0
   else
     !ndeg_filter necessary to converge to tolerance
     ndeg_filter_tolwfr = cheb_oracle1(eig_iband, lambda_minus, lambda_plus, chebfi%tolerance / res_iband, 1000)
     if (chebfi%oracle==1) then
       ndeg_filter_bands(iband) = MIN(ndeg_filter_max, ndeg_filter_tolwfr, chebfi%ndeg_filter)
     else if (chebfi%oracle==2) then
       !ndeg_filter necessary to decrease residual by a constant factor
       ndeg_filter_decrease = cheb_oracle1(eig_iband, lambda_minus, lambda_plus, chebfi%oracle_factor, 15)
       ndeg_filter_bands(iband) = MIN(ndeg_filter_max, ndeg_filter_tolwfr, ndeg_filter_decrease)
     else
       ABI_ERROR('Wrong value for chebfi%oracle')
     end if
   end if
 end do
 ndeg_filter = MAXVAL(ndeg_filter_bands)
 call xmpi_max(ndeg_filter,ndeg_filter_all,chebfi%comm_cols,ierr)
 ndeg_filter=ndeg_filter_all

 call xg_free(residu)
 ABI_FREE(ndeg_filter_bands)

end subroutine chebfi_set_ndeg_from_residu
!!***

end module m_chebfi2
!!***

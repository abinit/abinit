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
!! Copyright (C) 2018-2026 ABINIT group (BS, L. Baguet, IML)
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
   real(dp) :: tolerance                    ! Tolerance on the residu to stop the minimization
   real(dp) :: ecut                         ! Ecut for Chebfi oracle
   real(dp) :: oracle_factor                ! factor used to decrease residuals
   real(dp) :: oracle_min_occ               ! threshold on occupancies used for nbdbuf=-101

   integer :: paral_kgb                     ! MPI parallelization variables
   integer :: bandpp
   integer :: comm_cols
   integer :: comm_rows
   integer :: me_g0
   integer :: me_g0_fft

   logical :: from_linalg    ! Transposer allocates representation ColsRows if true, Linalg if false
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
   type(xg_t) :: xAXColsRows_W ! only used if from_linalg is false
   type(xg_t) :: xBXColsRows_W ! only used if from linalg is false

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

  ! Partition column vectors X=(X_lock X_active)
  ! Function: handles MPI distribution and transitions between them
  !-------------------------------------------------
  type, private :: bandPartition_t

    ! MPI sizes in linalg
    integer :: n_locked
    integer :: n_active

    ! MPI sizes in colsrows
    integer :: my_rank
    logical :: rank_active ! true if process treats active bands
    integer :: n_active_bandpp ! rank specific in colsrows

    ! communicator for processes treating active bands
    integer :: comm_active 

    type(xgBlock_t) :: linalg_active
    type(xgBlock_t) :: colsrows_locked
    type(xgBlock_t) :: colsrows_active
    
    type(xgTransposer_t) :: transposer_active

  end type bandPartition_t

!Public methods associated to 'chebfi' datatype
!-------------------------------------------------
 public :: chebfi_init
 public :: chebfi_free
 public :: chebfi_memInfo
 public :: chebfi_run
 public :: chebfi_runSlice
 public :: chebfi_runSubspaceIteration
 public :: chebfi_rayleighRitzQuotients
 public :: chebfi_computeNextOrderChebfiPolynom
 public :: chebfi_swapInnerBuffers
 public :: bandpassIndicator_sca    ! polynomial bandpass filter at point x

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
!!  me_g0_fft= 1 if this processors treats G=0 in FFT, 0 otherwise
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
!!  nbdbuf= number of bands in the buffer
!!  oracle= option compute ndeg_filter from residuals
!!  oracle_factor= factor used to decrease residuals
!!  oracle_min_occ= threshold on occupancies used for nbdbuf=-101
!!  gpu_option= flag. Enable GPU if true
!!  gpu_kokkos_nthrd= number of OpenMP offloaded threads used
!!  gpu_thread_limit= maximum number of OpenMP offloaded threads
!!  from_linalg= flag. Transposer allocates representation ColsRows if true, Linalg if false
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_init(chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp, &
                       ndeg_filter,nbdbuf,space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols, &
                       oracle,oracle_factor,oracle_min_occ,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit,from_linalg)

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
 logical       , intent(in   ), optional :: from_linalg

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
 chebfi%from_linalg = .true.
 if (present(from_linalg)) chebfi%from_linalg = from_linalg

 call chebfi_allocateAll(chebfi)

 call timab(tim_init,2,tsec)

 write(std_out,*) 'At chebfi_init:'
 write(std_out,*) 'chebfi%spacecom id=', chebfi%spacecom 
 write(std_out,*) 'chebfi%spacecom size=', xmpi_comm_size(chebfi%spacecom)
 write(std_out,*) 'chebfi%spacedim=', chebfi%spacedim
 write(std_out,*) 'chebfi%bandpp=', chebfi%bandpp

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
   if (chebfi%from_linalg) then
    total_spacedim = spacedim
    call xmpi_sum(total_spacedim,chebfi%spacecom,ierr)
    chebfi%total_spacedim = total_spacedim
   else
    chebfi%total_spacedim = chebfi%spacedim
   end if
   call xg_init(chebfi%X_NP,space,chebfi%total_spacedim,2*chebfi%bandpp,chebfi%spacecom,me_g0=chebfi%me_g0_fft,&
     & gpu_option=chebfi%gpu_option) !transposed arrays
   call xg_setBlock(chebfi%X_NP, chebfi%X_next, chebfi%total_spacedim, chebfi%bandpp)
   call xg_setBlock(chebfi%X_NP, chebfi%X_prev, chebfi%total_spacedim, chebfi%bandpp, fcol=chebfi%bandpp+1)
 end if

 !transposer will handle these arrays automatically
 if (chebfi%from_linalg) then
    call xg_init(chebfi%AX,space,spacedim,neigenpairs,chebfi%spacecom,me_g0=chebfi%me_g0,gpu_option=chebfi%gpu_option)
    call xg_init(chebfi%BX,space,spacedim,neigenpairs,chebfi%spacecom,me_g0=chebfi%me_g0,gpu_option=chebfi%gpu_option)
 else
    call xg_init(chebfi%xAXColsRows_W,space,chebfi%total_spacedim,chebfi%bandpp,chebfi%comm_rows,&
        me_g0=chebfi%me_g0_fft,gpu_option=chebfi%gpu_option)
    call xg_init(chebfi%xBXColsRows_W,space,chebfi%total_spacedim,chebfi%bandpp,chebfi%comm_rows,&
        me_g0=chebfi%me_g0_fft,gpu_option=chebfi%gpu_option)
    chebfi%xAXColsRows = chebfi%xAXColsRows_W%self
    chebfi%xBXColsRows = chebfi%xBXColsRows_W%self
 end if

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

 if (chebfi%from_linalg) then
    call xg_free(chebfi%AX)
    call xg_free(chebfi%BX)
 else 
    call xg_free(chebfi%xAXColsRows_W)
    call xg_free(chebfi%xBXColsRows_W)
 end if

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
   mem_sendrecv_CR = 0
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
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X= pointer to the function giving B^-1|X>
!!           B is typically the overlap operator S
!!  nspinor= number of spinorial components of the wavefunctions
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry) distributed in Linalg representation
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
 ! IML
 real(dp), pointer :: thetas(:,:) => null()

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

 ! IML debug
 write(std_out,*) 'maxeig_global=', maxeig_global
 !write(std_out,*) 'divresults=', xgBlock_getid(DivResults%self)
 !write(std_out,*) 'X0=', xgBlock_getid(X0)
 !write(std_out,*) 'xX=', xgBlock_getid(chebfi%xXColsRows)
 !call xgBlock_print(DivResults%self, std_out)
 flush(std_out)
 ! IML debug

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
 ABI_SFREE(ndeg_filter_bands)
 
 call timab(tim_transpose,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
 if (chebfi%paral_kgb == 1) then
   call xmpi_barrier(chebfi%spacecom)

   call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG)

   !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
   if (xmpi_comm_size(chebfi%spacecom) == 1) then 
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

 ! Apply Rayleigh-Ritz
 ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
 call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,chebfi%eigenvalues,ierr,0,tim_RR,&
&                     chebfi%gpu_option,solve_ax_bx=.true.)
 ABI_NVTX_END_RANGE()

 ! Start IML write converged eigenvalues here
 if (xmpi_comm_rank(chebfi%spacecom)==0) then
    if (chebfi%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(chebfi%eigenvalues)
    end if
    open(unit=1201, file='converged_eigenvalues.csv', status='replace')
    call xgBlock_reverseMap(chebfi%eigenvalues, thetas, rows=1, cols=chebfi%neigenpairs)
    write(1201,'(A)') "eig"
    do ideg=1, chebfi%neigenpairs
        write(1201,'(F15.5)') thetas(1,ideg)
    end do
    close(1201)
 end if 
 ! End IML

 ! Compute residual
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
 call xgBlock_scale(chebfi%xXColsRows, dble(1.0)/center, 1) !counter scale by 1/center

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

!!****f* m_chebfi2/chebfi_prepareAmpfactor
!! NAME
!! chebfi_prepAmpfactor
!! 
!! FUNCTION
!! Utility function to distribute correctly eigenvalues to MPI procs.
!! Performs MPI communication
!!
!! INPUTS
!! chebfi=
!! eigen= eigenvalues not distributed
!! 
!! OUTPUT
!! DivResults= prepared eigenvalues same array as in chebfi_run
!! 
!! SOURCE

subroutine chebfi_prepAmpfactor(chebfi, eigen, DivResults)

    implicit none

    ! Arguments ------------------------------------
    type(xg_t), intent(inout) :: DivResults
    type(xgBlock_t), intent(inout) :: eigen
    type(chebfi_t),  intent(inout) :: chebfi

    ! Local variables-------------------------------
    ! scalars
    integer :: space_res
    integer :: my_rank, num_proc, shift, ierr
    type(xgBlock_t) :: eigen_block
    ! Arrays
    integer, allocatable, target :: allbandpp(:)
    real(dp), allocatable, target :: theta_reshaped(:,:)
    integer, pointer :: allbandpp_ptr(:) => null()
    real(dp), pointer :: theta_reshaped_ptr(:,:) => null()
    real(dp), pointer :: theta(:,:) => null()

    ! *********************************************************************

    if (chebfi%space==SPACE_C) then
        space_res = SPACE_C
    else if (chebfi%space==SPACE_CR) then
        space_res = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if

    if (chebfi%paral_kgb == 0) then
        call xg_init(DivResults, space_res, rows=chebfi%neigenpairs, cols=1, gpu_option=chebfi%gpu_option)
        ! Fill DivResults with full eigenvalues
        ! TODO fix for workaround copy between space_res and SPACE_R
        if (space_res==SPACE_R) then
            call xgBlock_copy(eigen, DivResults%self)
        else
            ! workaround to copy from SPACE_R to SPACE_C
            ABI_MALLOC_IFNOT(theta_reshaped,(2,chebfi%neigenpairs))
            theta_reshaped_ptr => theta_reshaped
            call xgBlock_reverseMap(eigen, theta, rows=1, cols=chebfi%neigenpairs)
            theta_reshaped = 0.d0
            theta_reshaped(1,1:chebfi%neigenpairs) = theta(1,1:chebfi%neigenpairs)
#ifdef HAVE_OPENMP_OFFLOAD
            !$OMP TARGET ENTER DATA MAP(to:theta_reshaped) IF(chebfi%gpu_option==ABI_GPU_OPENMP)
#endif
            call xgBlock_map(eigen, theta_reshaped_ptr, space_res, rows=1, &
                cols=chebfi%neigenpairs, gpu_option=chebfi%gpu_option)
            call xgBlock_copy(eigen, DivResults%self)
#ifdef HAVE_OPENMP_OFFLOAD
            !$OMP TARGET EXIT DATA MAP(delete:theta_reshaped) IF(chebfi%gpu_option==ABI_GPU_OPENMP)
#endif
            ABI_SFREE(theta_reshaped)
        end if
    else
        call xg_init(DivResults, space_res, chebfi%bandpp, 1, gpu_option=chebfi%gpu_option) 
        if (xmpi_comm_size(chebfi%spacecom) > 1) then
            my_rank = xmpi_comm_rank(chebfi%spacecom)
            !shift = my_rank * chebfi%bandpp ! FIXME not working for different bandpp per rank
            num_proc = xmpi_comm_size(chebfi%spacecom)
            ABI_MALLOC_IFNOT(allbandpp,(num_proc))
            allbandpp_ptr => allbandpp
            call xmpi_allgather(chebfi%bandpp, allbandpp_ptr, chebfi%spacecom, ierr)
            if ( ierr /= xmpi_success ) then
                ABI_ERROR("Error while gathering number of bandpp for spacecom")
            end if
            if (my_rank==0) then
                shift = 0
            else
                shift = sum(allbandpp(1:my_rank)) ! fixed
            end if
            ABI_SFREE(allbandpp)
        else
            shift = 0
        end if
        ! Fill DivResults(bandpp,1) with block of eigen(neigenpairs,1) of size bandpp
        ! reshape to access column range
        call xgBlock_reshape(DivResults%self, 1, chebfi%bandpp)
        call xgBlock_reshape(eigen, 1, chebfi%neigenpairs)
        if (space_res==SPACE_R) then
            call xgBlock_setBlock(eigen, eigen_block, rows=1, cols=chebfi%bandpp, fcol=1+shift)
            call xgBlock_copy(eigen_block, DivResults%self)
        else
            ! workaround to copy from SPACE_R to SPACE_C
            ABI_MALLOC_IFNOT(theta_reshaped,(2,chebfi%bandpp))
            theta_reshaped_ptr => theta_reshaped
            call xgBlock_setBlock(eigen, eigen_block, rows=1, cols=chebfi%bandpp, fcol=1+shift)
            call xgBlock_reverseMap(eigen_block, theta, rows=1, cols=chebfi%bandpp)
            theta_reshaped = 0.d0
            theta_reshaped(1,1:chebfi%bandpp) = theta(1,1:chebfi%bandpp)
#ifdef HAVE_OPENMP_OFFLOAD
            !$OMP TARGET ENTER DATA MAP(to:theta_reshaped) IF(chebfi%gpu_option==ABI_GPU_OPENMP)
#endif
            call xgBlock_map(eigen_block, theta_reshaped_ptr, space_res, rows=1, &
                cols=chebfi%bandpp, gpu_option=chebfi%gpu_option)
            call xgBlock_copy(eigen_block, DivResults%self)
#ifdef HAVE_OPENMP_OFFLOAD
            !$OMP TARGET EXIT DATA MAP(delete:theta_reshaped) IF(chebfi%gpu_option==ABI_GPU_OPENMP)
#endif
            ABI_SFREE(theta_reshaped)
        end if
        ! restore dimensions
        call xgBlock_reshape(eigen, chebfi%neigenpairs, 1) 
        call xgBlock_reshape(DivResults%self, chebfi%bandpp, 1) 
    end if

    ! DivResults must be on CPU for ampfactor routine
    if (chebfi%gpu_option==ABI_GPU_OPENMP) then
        call xgBlock_copy_from_gpu(DivResults%self)
    end if

end subroutine chebfi_prepAmpfactor
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

!!****f* m_chebfi2/chebfi_ampfactorBandpass
!! NAME
!! chebfi_ampfactorBandpass
!!
!! FUNCTION
!! Compute amplification factor for bandpass polynomial
!! Numerical zero is 1e-3. Assumes prepAmpfactor prior to this.
!!
!! INPUTS
!! eig (:,:)= eigenvalues
!! lambda_minus,lambda_plus=
!! center,radius= used to rescale bandpass filter to interval
!! ndeg_filter= degree of bandpass polynomial filter
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  residu<type(xgBlock_t)>= vector of residuals
!!  chebfi <type(chebfi_t)>=all data used to apply Chebyshev Filtering algorithm
!!
!! SOURCE

subroutine chebfi_ampfactorBandpass(chebfi,DivResults,lambda_minus,lambda_plus,center,radius,ndeg_filter)

  implicit none

  ! Arguments ------------------------------------
  integer,           intent(in   ) :: ndeg_filter
  type(xgBlock_t),   intent(in   ) :: DivResults
  real(dp),          intent(in   ) :: lambda_minus
  real(dp),          intent(in   ) :: lambda_plus
  real(dp),          intent(in   ) :: center
  real(dp),          intent(in   ) :: radius
  type(chebfi_t),    intent(inout) :: chebfi

  ! Local variables-------------------------------
  ! scalars
  integer         :: iband,nbands
  real(dp)        :: ampfactor
  real(dp)        :: xred, ls, us
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
  ls = (lambda_minus-center)/radius
  us = (lambda_plus-center)/radius

  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=chebfi%bandpp)

  do iband = 1, nbands

    eig_per_band = eig(1,iband)

    !poly(x, a, b, n), where x,a,b are scaled!!!
    xred = (eig_per_band-center)/radius
    ampfactor = bandpassIndicator_sca(xred, ls, us, ndeg_filter)

    if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much
    
    call xgBlock_setBlock(chebfi%xXColsRows, X_part, chebfi%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(chebfi%xAXColsRows, AX_part, chebfi%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(chebfi%xBXColsRows, BX_part, chebfi%total_spacedim, 1, fcol=iband)

    !write(std_out,*) 'ampfactor, eig, iband=', eig_per_band, ampfactor, iband

    call xgBlock_scale(X_part, 1/ampfactor, 1)
    call xgBlock_scale(AX_part, 1/ampfactor, 1)
    call xgBlock_scale(BX_part, 1/ampfactor, 1)

  end do

end subroutine chebfi_ampfactorBandpass
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi/chebfi_runSlice
!! NAME
!! chebfi_runSlice
!! 
!! FUNCTION
!! Apply Polynomial filtering to set of vectors.
!! 
!! NOTES
!! List of differences with chebfi_run:
!! - X0 input and output is in ColsRows MPI representation (not Linalg!).
!! - the filter polynomial can be lowpass or bandpass.
!! - Rayleigh Quotients are read from eigen.
!! 
!! INPUT
!! chebfi <type(chebfi_t)>= all data used to apply Polynomial Filtering algorithm
!! X0= eigenvector guess distributed in ColsRows representation
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!! eigen= Rayleigh quotients associated to X0
!! residu= empty array
!! nspinor= number of spinorial components of the wavefunctions
!! lambda_minus= lower interval to amplify/diminish for bandpass/lowpass
!! lambda_plus= upper interval to amplify/diminish for bandpass/lowpass
!! mineig_global= guaranteed lower bound for entire spectrum
!! maxeig_global= guaranteed upper bound for entire spectrum
!! is_lowpass= flag. True if Chebyshev otherwise use bandpass Chebyshev-Jackson
!! nrows_blockrows= number of rows per MPI block in Linalg representation
!! 
!! SIDE EFFECTS
!! chebfi= workspaces used
!! X0= full set of vectors distributed in ColsRows representation 
!! eigen= full eigenvalues
!! residu= residuals, i.e. norm of (A-lambdaB)|X>
!! 
!! SOURCE

subroutine chebfi_runSlice(chebfi,X0,getAX_BX,getBm1X,eigen,residu,nspinor,&
        mineig_global,maxeig_global,lambda_minus,lambda_plus,is_lowpass,k_rank,nrows_blockrows)

    implicit none

    !Arguments ------------------------------------    
    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
    integer        , intent(in   ) :: nspinor
    integer        , intent(in   ) :: k_rank
    integer, pointer, intent(in  ) :: nrows_blockrows(:)
    real(dp)       , intent(in   ) :: mineig_global 
    real(dp)       , intent(in   ) :: maxeig_global
    real(dp)       , intent(in   ) :: lambda_minus
    real(dp)       , intent(in   ) :: lambda_plus
    logical        , intent(in   ) :: is_lowpass
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
    integer :: spacedim, neigenpairs, num_proc, ierr
    type(xg_t) :: X_k
    ! Arrays
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    real(dp) :: tsec(2)

    ! *********************************************************************
 
    if (chebfi%from_linalg) then
        ABI_ERROR("chebfi should be from colsrows")
    end if

    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs
    num_proc = xmpi_comm_size(chebfi%spacecom)
    chebfi%eigenvalues = eigen

    ABI_MALLOC_IFNOT(nrowsLinalg,(num_proc))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg = nrows_blockrows
    
    !write(std_out,*) 'getid inside runSlice xXColsRows', xgBlock_getId(chebfi%xXColsRows)
    !write(std_out,*) 'getid inside runSlice xAXColsRows', xgBlock_getId(chebfi%xAXColsRows)
    !flush(std_out)

    !A * Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    !write(std_out,*) 'getid inside runSlice xXColsRows (filtered 1)', xgBlock_getId(chebfi%xXColsRows)
    !write(std_out,*) 'getid inside runSlice xAXColsRows (filtered 1)', xgBlock_getId(chebfi%xAXColsRows)
    
    write(std_out,*) 'starting filter in proc', xmpi_comm_rank(chebfi%spacecom)
    flush(std_out)

    ! Apply polynomial filtering to active MPI ColsRows block-column
    if (is_lowpass) then
        ! [lambda_minus,lambda_plus) is diminished using Chebyshev
        write(std_out,*) 'lambda_minus=', lambda_minus
        write(std_out,*) 'lambda_plus=', lambda_plus
        flush(std_out)
        call chebfi_lowpassFilter(chebfi,eigen,lambda_minus,lambda_plus,getAX_BX,getBm1X)
    else
        write(std_out,*) 'lambda_minus=', lambda_minus
        write(std_out,*) 'lambda_plus=', lambda_plus
        write(std_out,*) 'mineig_global=', mineig_global
        write(std_out,*) 'maxeig_global=', maxeig_global
        flush(std_out)
        call chebfi_bandpassFilter(chebfi,eigen,lambda_minus,lambda_plus,mineig_global,&
            maxeig_global,getAX_BX,getBm1X)
    end if

    !write(std_out,*) 'getid inside runSlice xXColsRows (filtered N)', xgBlock_getId(chebfi%xXColsRows)
    !write(std_out,*) 'getid inside runSlice xAXColsRows (filtered N)', xgBlock_getId(chebfi%xAXColsRows)
    !flush(std_out)

    ! MPI transpose to linalg state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb==1) then

        ! Allocate chebfi%X
        call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
            STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0,&
            gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
            custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! Note: bandpp is custom because it is created from resource allocator

        ! Allocate chebfi%AX, chebfi%BX
        call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,&
            chebfi%AX%self,chebfi%xAXColsRows,STATE_COLSROWS)
        call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,&
            chebfi%BX%self,chebfi%xBXColsRows,STATE_COLSROWS)
        ! Note: at this point chebfi%AX and chebfi%BX are empty. Must transpose
        !       to fill with correct values.

        ! todo use copy constructor to create an object by copying an existing object
        ! actually copy constructor *allocates* memory for chebfi%AX. Write a version
        ! that does not allocate memory and only reassigns pointers. Can pointers be reassigned
        ! directly then used in the global constructor? Perform tests.

        !write(std_out,*) 'getid before transpose AX', xgBlock_getId(chebfi%AX%self)
        !write(std_out,*) 'getid before transpose xAX', xgBlock_getId(chebfi%xAXColsRows)
        !write(std_out,*) 'getid before transpose xX', xgBlock_getId(chebfi%xXColsRows)
        !flush(std_out)

        chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
        chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
        chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
        call xgTransposer_transpose(chebfi%xgTransposerAX, STATE_LINALG)
        call xgTransposer_transpose(chebfi%xgTransposerBX, STATE_LINALG)
        call xmpi_barrier(chebfi%spacecom)

        !write(std_out,*) 'getid after transpose AX', xgBlock_getId(chebfi%AX%self)
        !flush(std_out)

    else
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
    end if
    call timab(tim_transpose,2,tsec)
    ABI_NVTX_END_RANGE()

    if (chebfi%paral_kgb==1 .and. rows(chebfi%X) /= nrowsLinalg(xmpi_comm_rank(chebfi%spacecom)+1)) then
        ABI_ERROR("wrong linalg representation")
    end if
    write(std_out,'(a,i6,i6)') 'local # proc has # rows ', xmpi_comm_rank(chebfi%spacecom), rows(chebfi%X)

    !write(std_out,*) 'chebfi%eigenvalues before RR'
    !call xgBlock_print(chebfi%eigenvalues,std_out)
    
    !write(std_out,*) 'id of X, (before RR) ncols=', xgBlock_getId(chebfi%X), cols(chebfi%X)
    !flush(std_out)

    !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,chebfi%gpu_option,AX=chebfi%xAXColsRows)
    
    ! Apply Rayleigh-Ritz to active MPI Linalg row-block
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
    call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,eigen,ierr,0,tim_RR,&
        chebfi%gpu_option,solve_ax_bx=.true.)
    ABI_NVTX_END_RANGE()
    
    !write(std_out,*) 'id of X, (after RR) ncols=', xgBlock_getId(chebfi%X), cols(chebfi%X)

    if ( ierr /= 0 ) then
        ABI_WARNING("RayleighRitz did not work")
    else
        !write(std_out,*) 'is lowpass=', is_lowpass
        !write(std_out,*) 'chebfi%eigenvalues after RR'
        !call xgBlock_print(chebfi%eigenvalues,std_out)
        !flush(std_out)
    end if

    ! Compute residual norm *squared*
    if (chebfi%paw) then
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
    else
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
    end if
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu) ! performs MPI comm

    !write(std_out,*) 'max colwise residual norm squared='; call xgBlock_print(residu, std_out)
    !flush(std_out)

    ! Copy in Linalg representation (see chebfi_run, kept for reference)
    ! call xgBlock_copy(chebfi%X,X0)
    
    ! MPI Transpose to recover colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)
        call xgTransposer_transpose(chebfi%xgTransposerAX, STATE_COLSROWS)
        call xgTransposer_transpose(chebfi%xgTransposerBX, STATE_COLSROWS)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)
 
    ! Copy in ColsRows representation
    call xgBlock_copy(chebfi%xXColsRows, X0)

    if (cols(X0) /= chebfi%bandpp) then
        ABI_ERROR('wrong colsrows representation')
    end if
    write(std_out,'(a,i6,i6,i6)') 'local # proc has # rows cols ', xmpi_comm_rank(chebfi%spacecom), rows(X0), cols(X0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (chebfi%gpu_option==ABI_GPU_KOKKOS) then
        call gpu_device_synchronize()
    end if
#endif

    ! Free transposer objects
    if (chebfi%paral_kgb == 1) then
        call xgTransposer_free(chebfi%xgTransposerX)
        call xgTransposer_free(chebfi%xgTransposerAX)
        call xgTransposer_free(chebfi%xgTransposerBX)
    end if
    
    ! Free temporary memory
    ABI_SFREE(nrowsLinalg)

end subroutine chebfi_runSlice
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi/chebfi_runSubspaceIteration
!! NAME
!! chebfi_runSubspaceIteration 
!! 
!! FUNCTION
!! Applies subspace iteration by Chebyshev polynomial filtering on slice
!! then extracts eigenvectors using Rayleigh-Ritz. This execution is completely
!! local on slice processors. Each slice does not see others. 
!! 
!! NOTES
!! Restart logic in subspace iteration:
!!   Allocate constructor (state=colsrows)
!!   while (convergence not reached):
!!    |  if (state=linalg) Transpose
!!    |  Filter
!!    |  Transpose (state=linalg)
!!    |  Orthogonalize
!!   Rayleigh-Ritz
!! 
!! INPUTS
!! X0 input vectors distibuted by colsrows along slice processors
!! 
!! OUTPUTS
!! X0 (in-place) converged vectors distributed by linalg along slice processors
!! 
!! dev=============
!! - version 1: count number of converged vectors within iteration
!!
!! SOURCE

subroutine chebfi_runSubspaceIteration(chebfi,X0,getAX_BX,getBm1X,eigen,residu,nspinor,&
        mineig_global,maxeig_global,lambda_minus,lambda_plus,is_lowpass,k_rank,nrows_blockrows)

    implicit none

    !Arguments ------------------------------------    
    type(chebfi_t) , intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
    integer        , intent(in   ) :: nspinor
    integer        , intent(in   ) :: k_rank
    integer, pointer, intent(in  ) :: nrows_blockrows(:)
    real(dp)       , intent(in   ) :: mineig_global 
    real(dp)       , intent(in   ) :: maxeig_global
    real(dp)       , intent(in   ) :: lambda_minus
    real(dp)       , intent(in   ) :: lambda_plus
    logical        , intent(in   ) :: is_lowpass
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
    integer :: iter_subspace, niter_subspace_max, n_locked
    integer :: spacedim, neigenpairs, num_proc, ierr
    type(xg_t) :: X_k
    type(xg_t) :: resid_active
    type(bandPartition_t) :: part_X
    type(bandPartition_t) :: part_AX
    type(bandPartition_t) :: part_BX
    ! Arrays
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    real(dp) :: tsec(2)

    ! *********************************************************************
 
    if (chebfi%from_linalg) then
        ABI_ERROR("chebfi should be from colsrows")
    end if

    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs
    num_proc = xmpi_comm_size(chebfi%spacecom)
    chebfi%eigenvalues = eigen

    ABI_MALLOC_IFNOT(nrowsLinalg,(num_proc))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg = nrows_blockrows
    
    call xg_init(resid_active, SPACE_R, neigenpairs, 1, gpu_option=chebfi%gpu_option)
    
    !write(std_out,*) 'getid inside runSlice xXColsRows', xgBlock_getId(chebfi%xXColsRows)
    !write(std_out,*) 'getid inside runSlice xAXColsRows', xgBlock_getId(chebfi%xAXColsRows)
    !flush(std_out)

    ! Initialize values of AX and BX
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ! Allocate memory for linalg distribution from colsrows
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb==1) then

        ! Allocate chebfi%X
        call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
            STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0,&
            gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
            custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! Note: bandpp is custom because it is created from resource allocator

        ! Allocate chebfi%AX, chebfi%BX
        call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,&
            chebfi%AX%self,chebfi%xAXColsRows,STATE_COLSROWS)
        call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,&
            chebfi%BX%self,chebfi%xBXColsRows,STATE_COLSROWS)
        ! Note: at this point chebfi%AX and chebfi%BX are empty. Must transpose
        !       to fill with correct values.

        chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
        chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
        chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

    else
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xAXColsRows, chebfi%AX%self, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%xBXColsRows, chebfi%BX%self, spacedim, neigenpairs)
    end if
    call timab(tim_transpose,2,tsec)
    ABI_NVTX_END_RANGE()
    
    write(std_out,*) 'starting filter in proc', xmpi_comm_rank(chebfi%spacecom)
    flush(std_out)

    niter_subspace_max = 10
    n_locked = 0

    ! chebfi%X    contains active vectors
    ! X_lock      contains locked vectors

    ! todo initialize partitions

    xX_active = chebfi%xXColsRows
    xAX_active = chebfi%xAXColsRows
    xBX_active = chebfi%xBXColsRows
    
    call bandPartition_setColsRowsActive(X_part , chebfi%xXColsRows, n_locked)
    call bandPartition_setColsRowsActive(AX_part, chebfi%AXColsRows, n_locked)
    call bandPartition_setColsRowsActive(BX_part, chebfi%BXColsRows, n_locked)

    do iter_subspace=1, niter_subspace_max

        write(std_out,*) 'subspace iteration no=', iter_subspace
        flush(std_out)

        ! si n_locked_bandpp est trop petit, on ne va pas pouvoir distribuer
        ! cas minimal: si n_active < bandpp alors utilise 1 seul proc bandpp=n_active pas de distr
        ! en gros il faut mettre à jour comm_cols
        ! comm_cols peut être le communicateur de p procs (p=size(comm_rows))
        ! ou comm_cols peut être un sous-communicateur. Il nous faut une fonction qui décide
        ! combien de processus les vecteurs active ont besoin.
        ! pour cela on fait p_active =
        ! min p = 1 et max p = nproc slice
        ! avec capacité en bandes : min 1 et max bandpp

        call bandPartition_getActiveDistribution(X_part)
        call bandPartition_copyActiveDistribution(X_part, AX_part)
        call bandPartition_copyActiveDistribution(X_part, BX_part)

        ! Construct transposer only for active vectors from linalg distribution
        call chebfi_constructActiveTransposers(chebfi, n_locked, &
            xgTransposerXactive, xgTransposerAXactive, xgTransposerBXactive)
        
        ! Transpose active vectors to colsrows distribution
        call timab(tim_transpose,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
        call xgTransposer_transpose(xgTransposerXactive, STATE_COLSROWS)
        call xgTransposer_transpose(xgTransposerAXactive, STATE_COLSROWS)
        call xgTransposer_transpose(xgTransposerBXactive, STATE_COLSROWS)
        call timab(tim_transpose,2,tsec)
        ABI_NVTX_END_RANGE()

        if (iter_subspace>1 .and. chebfi%paral_kgb==1) then
            call timab(tim_transpose,1,tsec)
            ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
            call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)
            call xgTransposer_transpose(chebfi%xgTransposerAX, STATE_COLSROWS)
            call xgTransposer_transpose(chebfi%xgTransposerBX, STATE_COLSROWS)
            call timab(tim_transpose,2,tsec)
            ABI_NVTX_END_RANGE()

            if (n_locked>1) then
                ! n_locked_mpi each rank has different number of locked vectors todo
                call xgBlock_setBlock(chebfi%xXColsRows , xX_active , spacedim, n_locked_mpi)
                call xgBlock_setBlock(chebfi%xAXColsRows, xAX_active, spacedim, n_locked_mpi)
                call xgBlock_setBlock(chebfi%xBXColsRows, xBX_active, spacedim, n_locked_mpi)
            end if
        end if
        ! todo this will transform the active+locked in colsrows. Locked is not necessary
        ! therefore try to communicate less data is possible.
    
        ! ############################ Filter active  ##############################
        ! ############################ column vectors ##############################
        if (is_lowpass) then
            ! todo set active pointers
            call chebfi_lowpassFilterActive(chebfi,xX_active,xAX_active,xBX_active,eigen,&
                lambda_minus,lambda_plus,getAX_BX,getBm1X)
        else
            call chebfi_bandpassFilterActive(chebfi,xX_active,xAX_active,xBX_active,eigen,&
                lambda_minus,lambda_plus,mineig_global,maxeig_global,getAX_BX,getBm1X)
        end if

        call bandPartition_transposeActive(X_part)
        call bandPartition_transposeActive(AX_part)
        call bandPartition_transposeActive(BX_part)
        
        write(std_out,'(a,i6,i6)') 'chebfi%xXColsRows # rows # cols ', &
            rows(chebfi%xXColsRows), cols(chebfi%xXColsRows)
        flush(std_out)

        if (chebfi%paral_kgb==1) then
            call timab(tim_transpose,1,tsec)
            ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
            call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
            call xgTransposer_transpose(chebfi%xgTransposerAX, STATE_LINALG)
            call xgTransposer_transpose(chebfi%xgTransposerBX, STATE_LINALG)
            call timab(tim_transpose,2,tsec)
            ABI_NVTX_END_RANGE()
        end if
    
        write(std_out,'(a,i6,i6)') 'chebfi%X # rows # cols ', rows(chebfi%X), cols(chebfi%X)
        flush(std_out)

        if (n_locked>0) then
            call chebfi_deflateWrtLocked(chebfi, n_locked)
        end if

        call xg_Borthonormalize(chebfi%X,chebfi%BX%self,ierr,1,chebfi%gpu_option,AX=chebfi%AX%self)

        call chebfi_getSubspaceResidual(chebfi, resid_active%self)

        ! todo change to swapConvergedVectors
        call chebfi_lockConvergedVectors(chebfi, resid_active%self, 1e-3_dp, n_locked) 

        ! Definir les espace en utilisant n_locked, n_active
        call bandPartition_setLinalg(chebfi, X_part, AX_part, BX_part)


        ! todo debug by recomputing the getSubspaceResidual for locked vectors and verify it is smaller than tol
        
        ! lock vectors in small residual
        ! actually redefine pointers Xactive and Xlock pointing to X that's all
    
    end do

    ! Apply Rayleigh-Ritz to active MPI Linalg row-block
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
    call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,eigen,ierr,0,tim_RR,&
        chebfi%gpu_option,solve_ax_bx=.true.)
    ABI_NVTX_END_RANGE()
    
    !write(std_out,*) 'id of X, (after RR) ncols=', xgBlock_getId(chebfi%X), cols(chebfi%X)

    if ( ierr /= 0 ) then
        ABI_WARNING("RayleighRitz did not work")
    end if

    ! Compute residual norm *squared*
    if (chebfi%paw) then
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
    else
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
    end if
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu) ! performs MPI comm

    ! MPI Transpose to recover colsrows state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)
        call xgTransposer_transpose(chebfi%xgTransposerAX, STATE_COLSROWS)
        call xgTransposer_transpose(chebfi%xgTransposerBX, STATE_COLSROWS)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)
            call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)
        call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)
 
    ! Copy in ColsRows representation
    call xgBlock_copy(chebfi%xXColsRows, X0)

    if (cols(X0) /= chebfi%bandpp) then
        ABI_ERROR('wrong colsrows representation')
    end if
    write(std_out,'(a,i6,i6,i6)') 'local # proc has # rows cols ', xmpi_comm_rank(chebfi%spacecom), rows(X0), cols(X0)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (chebfi%gpu_option==ABI_GPU_KOKKOS) then
        call gpu_device_synchronize()
    end if
#endif

    ! Free transposer objects
    if (chebfi%paral_kgb == 1) then
        call xgTransposer_free(chebfi%xgTransposerX)
        call xgTransposer_free(chebfi%xgTransposerAX)
        call xgTransposer_free(chebfi%xgTransposerBX)
    end if
    
    ! Free temporary memory
    ABI_SFREE(nrowsLinalg)
    call xg_free(resid_active)

end subroutine chebfi_runSubspaceIteration
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi/chebfi_lowpassFilter
!! NAME
!! chebfi_lowpassFilter
!!
!! FUNCTION
!! Apply Lowpass filter using Chebyshev polynomial on a set of vectors.
!! Amplifies interval [-oo, lambda_minus) and diminishes [lambda_minus,lambda_plus).
!!
!! INPUTS
!! chebfi <type(chebfi_t)>=memory workspace used to apply filter
!! eigen= Rayleigh quotients to use in amplification of chebfi%xXColsRows
!! lambda_minus= lower bound of interval to diminish
!! lambda_plus= upper bound of interval to diminish
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  chebfi%xXColsRows= Filtered vectors to use in Subspace iteration
!!
!! SOURCE

subroutine chebfi_lowpassFilter(chebfi,eigen,lambda_minus,lambda_plus,getAX_BX,getBm1X)

    implicit none

    ! Arguments ------------------------------------
    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: eigen
    real(dp), intent(in) :: lambda_minus
    real(dp), intent(in) :: lambda_plus
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
    integer :: ideg
    real(dp) :: center, radius, one_over_r, two_over_r
    type(xg_t) :: DivResults ! Rayleigh quotients
    ! Arrays
    real(dp) :: tsec(2)
    integer, allocatable :: ndeg_filter_bands(:)

    ! *********************************************************************
   
    if (chebfi%paral_kgb == 0) then
        ABI_MALLOC_IFNOT(ndeg_filter_bands,(chebfi%neigenpairs))
    else    
        ABI_MALLOC_IFNOT(ndeg_filter_bands,(chebfi%bandpp))
    end if
    
    ! [lambda_minus,lambda_plus) is diminished using Chebyshev
    write(std_out,*) '@lowpass lambda_minus=', lambda_minus
    write(std_out,*) '@lowpass lambda_plus=', lambda_plus
    flush(std_out)

    ! Filter parameters
    ndeg_filter_bands(:) = chebfi%ndeg_filter
    center = (lambda_plus + lambda_minus)*0.5
    radius = (lambda_plus - lambda_minus)*0.5
    one_over_r = 1/radius
    two_over_r = 2/radius

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
    do ideg = 0, chebfi%ndeg_filter - 1

        ! X_next=2/r*(AX_next-c*X_next)-X_prev
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, chebfi%spacedim, chebfi%neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
        end if
        ABI_NVTX_END_RANGE()

        !A * Psi (=AX_next=A*X_next)
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)

    end do ! ideg
    ABI_NVTX_END_RANGE()

    ! Avoid overflow rescale
    call chebfi_prepAmpfactor(chebfi, eigen, DivResults)
    call chebfi_ampfactor(chebfi, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)

    ! Free temporary memory
    call xg_free(DivResults)
    ABI_SFREE(ndeg_filter_bands)

end subroutine chebfi_lowpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi/chebfi_lowpassFilterActive
!! NAME
!! chebfi_lowpassFilterActive
!!
!! FUNCTION
!! Apply Lowpass filter on Active vectors only
!! TODO this function will replace lowpassFilter
!!
!! SOURCE

subroutine chebfi_lowpassFilterActive(chebfi,X_active,AX_active,BX_active,eigen,&
        lambda_minus,lambda_plus,getAX_BX,getBm1X)

    implicit none

    ! Arguments ------------------------------------
    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X_active
    type(xgBlock_t), intent(inout) :: AX_active
    type(xgBlock_t), intent(inout) :: BX_active
    type(xgBlock_t), intent(inout) :: eigen
    real(dp), intent(in) :: lambda_minus
    real(dp), intent(in) :: lambda_plus
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
    integer :: ideg
    real(dp) :: center, radius, one_over_r, two_over_r
    type(xg_t) :: DivResults ! Rayleigh quotients
    ! Arrays
    real(dp) :: tsec(2)
    integer, allocatable :: ndeg_filter_bands(:)

    ! *********************************************************************
   
    if (chebfi%paral_kgb == 0) then
        ABI_MALLOC_IFNOT(ndeg_filter_bands,(chebfi%neigenpairs))
    else    
        ABI_MALLOC_IFNOT(ndeg_filter_bands,(chebfi%bandpp))
    end if
    
    ! [lambda_minus,lambda_plus) is diminished using Chebyshev
    write(std_out,*) '@lowpass lambda_minus=', lambda_minus
    write(std_out,*) '@lowpass lambda_plus=', lambda_plus
    flush(std_out)

    ! Filter parameters
    ndeg_filter_bands(:) = chebfi%ndeg_filter
    center = (lambda_plus + lambda_minus)*0.5
    radius = (lambda_plus - lambda_minus)*0.5
    one_over_r = 1/radius
    two_over_r = 2/radius

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
    do ideg = 0, chebfi%ndeg_filter - 1

        ! X_next=2/r*(AX_next-c*X_next)-X_prev
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, chebfi%spacedim, chebfi%neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
        end if
        ABI_NVTX_END_RANGE()

        !A * Psi (=AX_next=A*X_next)
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
        call getAX_BX(X_active, AX_active, BX_active)
        call xgBlock_zero_im_g0(AX_active)
        call xgBlock_zero_im_g0(BX_active)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)

    end do ! ideg
    ABI_NVTX_END_RANGE()

    ! Avoid overflow rescale
    call chebfi_prepAmpfactor(chebfi, eigen, DivResults)
    call chebfi_ampfactor(chebfi, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)

    ! Free temporary memory
    call xg_free(DivResults)
    ABI_SFREE(ndeg_filter_bands)

end subroutine chebfi_lowpassFilterActive
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi/chebfi_bandpassFilter
!! NAME
!! chebfi_bandpassFilter
!!
!! FUNCTION
!! Apply Bandpass filter using Chebyshev-Jackson polynomial, that is an 
!! approximation of Heaviside step function by a Chebyshev expansion plus
!! Jackson damping to reduce oscillations, applied on a set of vectors. 
!! Amplifies interval [lambda_minus, lambda_plus).
!! 
!! INPUTS
!! chebfi <type(chebfi_t)>=memory workspace used to apply filter
!! eigen= Rayleigh quotients to use in amplification of chebfi%xXColsRows
!! lambda_minus= lower bound of interval to amplify
!! lambda_plus= upper bound of interval to amplify
!! mineig_global= used to rescale to [-1,1), will le -1
!! maxeig_global= used to rescale to [-1,1), will be 1
!! getAX_BX= pointer to the function giving A|X> and B|X>
!!           A is typically the Hamiltonian H, and B the overlap operator S
!! getBm1X= pointer to the function giving B^-1|X>
!!          B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  chebfi%xXColsRows= Filtered vectors to use in Subspace iteration
!!
!! SOURCE

subroutine chebfi_bandpassFilter(chebfi,eigen,lambda_minus,lambda_plus,mineig_global,&
        maxeig_global,getAX_BX,getBm1X)

    implicit none

    ! Arguments ------------------------------------
    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: eigen
    real(dp), intent(in) :: lambda_minus
    real(dp), intent(in) :: lambda_plus
    real(dp), intent(in) :: mineig_global
    real(dp), intent(in) :: maxeig_global
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

    ! Local variables-------------------------------
    integer :: ndeg, n
    real(dp) :: center, radius, one_over_r, two_over_r
    real(dp) :: ls, us, cdeg, mu, damp
    type(xg_t) :: Heaviside
    type(xg_t) :: DivResults ! Rayleigh quotients
    real(dp) :: tsec(2)

    ! *********************************************************************

    ndeg = chebfi%ndeg_filter 
  
    ! Allocate memory for Chebyshev expansion of Heaviside step function
    call xg_init(Heaviside, chebfi%space, chebfi%total_spacedim, chebfi%bandpp, chebfi%spacecom, &
        gpu_option=chebfi%gpu_option)
    !call xgBlock_zero(Heaviside%self)

    write(std_out,*) '@bandpass lambda_minus=', lambda_minus
    write(std_out,*) '@bandpass lambda_plus=', lambda_plus
    write(std_out,*) '@bandpass mineig_global=', mineig_global
    write(std_out,*) '@bandpass maxeig_global=', maxeig_global
    flush(std_out)
 
    ! Filter parameters
    center = (maxeig_global + mineig_global)*0.5
    radius = (maxeig_global - mineig_global)*0.5
    one_over_r = 1/radius
    two_over_r = 2/radius

    ! Scaled slice bounds to be amplified
    ls = (lambda_minus - center) / radius
    us = (lambda_plus - center) / radius
    
    ! TODO IL 10/3/2025 Deflate vectors to reduce linear dependence: Y=X-(B-projection)
    !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,chebfi%gpu_option,AX=chebfi%xAXColsRows)

    ! Initialize expansion: Heaviside = mu(0)*damp(0)*X + Heaviside
    cdeg = Pi/(ndeg+2)
    mu = 1/Pi*(ACOS(ls)-ACOS(us))
    damp = 1.d0 ! Jackson damping
    call xgBlock_saxpy(Heaviside%self, mu*damp, chebfi%xXColsRows)

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
    do n = 0, ndeg - 1  

        ! X_next=2/r*(AX_next-c*X_next)-X_prev
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, n, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, chebfi%spacedim, chebfi%neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, chebfi%bandpp)
        end if
        ABI_NVTX_END_RANGE()

        ! Update expansion: Heaviside = damp(i+1)*mu(i+1)*X_next + Heaviside
        mu = 2/Pi * (SIN((n+1)*ACOS(ls)) - SIN((n+1)*ACOS(us)))/(n+1)
        damp = ((1 - (n+1)/(ndeg+2))*SIN(cdeg)*COS((n+1)*cdeg) + 1/(ndeg+2)*COS(cdeg)*SIN((n+1)*cdeg))/SIN(cdeg)
        call xgBlock_saxpy(Heaviside%self, mu*damp, chebfi%xXColsRows)
   
        ! Store final expansion before exit, X_next=Heaviside
        if (n==ndeg-1) then
            call xgBlock_copy(Heaviside%self, chebfi%xXColsRows)
        end if

        !A * Psi (=AX_next=A*X_next)
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows, chebfi%xAXColsRows, chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
        call timab(tim_getAX_BX,2,tsec)
    
    end do ! end n
    ABI_NVTX_END_RANGE()

    ! Avoid overflow
    !call chebfi_prepAmpfactor(chebfi, eigen, DivResults)
    !call chebfi_ampfactorBandpass(chebfi, DivResults%self, lambda_minus, lambda_plus, center, radius, ndeg)

    ! Free temporary memory
    call xg_free(DivResults)
    call xg_free(Heaviside)

end subroutine chebfi_bandpassFilter
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

!----------------------------------------------------------------------

!!****f* m_chebfi2/bandpassIndicator_sca
!! NAME
!! bandpassIndicator_sca
!!
!! FUNCTION
!! Scalar Chebyshev-Jackson polynomial filter f(x) approximating an 
!! indicator function, using degree deg evaluated at point x=t
!!
!! INPUTS
!!  a,b=    interval to amplify included in -1,1
!!  t=      scalar to evaluate filter on
!!  deg=    order of Chebyshev expansion
!!
!! OUTPUT
!!  res
!!
!! SOURCE

function bandpassIndicator_sca(t,a,b,deg) result(f_t)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: t,a,b
    integer , intent(in ) :: deg

    real(dp) :: f_t
    
    !Local variables-------------------------------
    real(dp) :: yt0,yt,yt_swap,ck,mu,damp
    integer  :: i
    
    ! *********************************************************************

    ! init cheby of deg=0,1 eval at t
    yt0 = 1.d0
    yt = t

    ! init filter for deg=0
    ck = Pi/(deg+2)
    mu = 1/Pi*(ACOS(a)-ACOS(b))
    damp = 1.d0
    f_t = mu * damp * yt0

    do i=1,deg 
        
        ! Update damping and expansion coefficient
        mu = 2/Pi * (SIN(i*ACOS(a)) - SIN(i*ACOS(b)))/i
        damp = ((1 - i/(deg+2))*SIN(ck)*COS(i*ck) + 1/(deg+2)*COS(ck)*SIN(i*ck))/SIN(ck)

        ! Sum terms
        f_t = f_t + mu * damp * yt

        ! Update Chebyshev polynomial
        yt_swap = yt
        yt = 2 * t * yt - yt0
        yt0 = yt_swap
        
    end do

end function bandpassIndicator_sca
!!***

!----------------------------------------------------------------------

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
 ABI_SFREE(ndeg_filter_bands)

end subroutine chebfi_set_ndeg_from_residu
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_Bdeflate_active
!! NAME
!! chebfi_Bdeflate_active
!!
!! FUNCTION
!! Orthogonalize active block with respect to locked block in B-basis
!!
!! INPUT
!! Xlock assumed to be B-orthonormal
!! 
!! OUTPUT
!! in-place
!! Xactive = Xactive - Xlock*(BXlock'*Xactive)

!  subroutine chebfi_Bdeflate_active(lobpcg,var,iblock)
!
!    type(lobpcg_t) , intent(inout) :: lobpcg
!    type(xgBlock_t), intent(inout) :: var
!    integer        , intent(in   ) :: iblock
!    integer :: previousBlock
!    integer :: blockdim
!    integer :: spacedim
!    integer :: space_buf
!    type(xg_t) :: buffer
!    double precision :: tsec(2)
!
!    call timab(tim_ortho,1,tsec)
!    ABI_NVTX_START_RANGE(NVTX_LOBPCG2_ORTHO_X_WRT)
!
!    blockdim = lobpcg%blockdim
!    spacedim = lobpcg%spacedim
!    previousBlock = (iblock-1)*lobpcg%blockdim
!
!    ! replace var by chebfi%Xactive
!
!    space_buf = space(var)
!    if (space(var)==SPACE_CR) then
!      space_buf = SPACE_R
!    end if
!   call xg_init(buffer,space_buf,previousBlock,blockdim,comm=lobpcg%spacecom,gpu_option=lobpcg%gpu_option)
!
!    ! buffer = BX0^T*X
!    call xgBlock_gemm('t','n',1.0d0,lobpcg%BX0,var,0.d0,buffer%self,comm=lobpcg%spacecom)
!
!    ! sum all process contribution of X
!    ! X = - X0*(BX0^T*X) + X 
!   call xgBlock_gemm('n','n',-1.0d0,lobpcg%X0,buffer%self,1.0d0,var)
!
!    call xg_free(buffer)
!
!   ABI_NVTX_END_RANGE()
!   call timab(tim_ortho,2,tsec)
!
!  end subroutine lobpcg_orthoXwrtBlocks
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_getSubspaceResidual
!! NAME
!! chebfi_getSubspaceResidual
!! 
!! FUNCTION
!! || (I-QQ'B)Aq_j ||^2 for j=1,..,bandpp
!! stored to resid

subroutine chebfi_getSubspaceResidual(chebfi, resid)

    implicit none

    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: resid
    
    type(xg_t) :: M, R 
    integer :: nrows, ncols
    integer :: space_buf

! *********************************************************************

    nrows = rows(chebfi%X) ! todo this should be chebfi%spacedim
    ncols = chebfi%neigenpairs ! todo active part only

    if (chebfi%space==SPACE_C) then
        space_buf = SPACE_C
    else if (chebfi%space==SPACE_CR) then
        space_buf = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if

    call xg_init(M, space_buf, ncols, ncols, comm=chebfi%spacecom, gpu_option=chebfi%gpu_option)
    call xg_init(R, space_buf, nrows, ncols, comm=chebfi%spacecom, gpu_option=chebfi%gpu_option)

    ! Compute M = X^T*AX 
    ! sum all process contribution
    call xgBlock_gemm('t','n',1.0d0,chebfi%X,chebfi%AX%self,0.d0,M%self,comm=chebfi%spacecom)

    ! Compute R = AX - BX*M
    call xgBlock_copy(chebfi%AX%self, R%self)
    call xgBlock_gemm('n','n',-1.0d0,chebfi%BX%self,M%self,1.0d0,R%self)

    call xgBlock_colwiseNorm2(R%self, resid, comm_loc=xmpi_comm_null)

    call xg_free(M)
    call xg_free(R)

end subroutine chebfi_getSubspaceResidual
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_lockConvergedVectors
!! NAME
!! chebfi_lockConvergedVectors
!! 
!! FUNCTION
!! In-place contiguous-in-memory locking of column vectors in row distribution.
!! Algorithm is based on swapping: at the end first k are locked last are active.
!! | locked columns | active columns |
!! |   1 ... k      |  k+1 ... m     |
!! TODO logic could also be applied to divide Xext to slices.
!! 

subroutine chebfi_lockConvergedVectors(chebfi, resid, tol, n_locked)

    implicit none

    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(in) :: resid
    real(dp), intent(in) :: tol
    integer, intent(out) :: n_locked

    integer :: i, nrows, n_active
    integer :: j, left
    real(dp), pointer :: resid_vals(:,:)
    logical, allocatable :: mask(:)
    logical, allocatable :: is_locked(:)
    integer, allocatable :: idx(:)
    integer :: space_buf
    type(xg_t) :: T_swap ! temporary column buffer

! *********************************************************************

    nrows = rows(resid)

    if (chebfi%space==SPACE_C) then
        space_buf = SPACE_C
    else if (chebfi%space==SPACE_CR) then
        space_buf = SPACE_R
    else
        ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
    end if
    
    ABI_MALLOC(mask, (nrows))
    ABI_MALLOC(is_locked, (nrows))
    
    call xg_init(T_swap, space_buf, nrows, 1, comm=chebfi%spacecom, gpu_option=chebfi%gpu_option)
    
    call xgBlock_reverseMap(resid, resid_vals, rows=nrows, cols=1)

    mask = resid_vals(:,1) < tol
    idx = pack([(i, i=1,nrows)], mask)
    
    n_locked = size(idx)
    is_locked = .false.
    is_locked(idx) = .true.

    ! partition via swapping (in-place)
    left = 1
    do j = 1, nrows
        if (is_locked(j)) then
            if (j /= left) then
                call xgBlock_colwiseSwap(chebfi%X, j, left, T_swap%self)
                call xgBlock_colwiseSwap(chebfi%AX%self, j, left, T_swap%self)
                call xgBlock_colwiseSwap(chebfi%BX%self, j, left, T_swap%self)

                ! keep mask consistent after swap
                is_locked(j)    = is_locked(left)
                is_locked(left) = .true.
            end if
            left = left + 1
        end if
    end do

    ! update pointers
    n_active = chebfi%neigenpairs - n_locked
    call xgBlock_setBlock(chebfi%X      , X_active , rows(chebfi%X), n_active, fcol=n_locked+1)
    call xgBlock_setBlock(chebfi%AX%self, AX_active, rows(chebfi%X), n_active, fcol=n_locked+1)
    call xgBlock_setBlock(chebfi%BX%self, BX_active, rows(chebfi%X), n_active, fcol=n_locked+1)

    ! what I do is that I use X_lock vectors and update chebfi%X

    ! todo this should be updated automatically after transposing...
    ! ongoing modif in xgTransposer to be able to inverse on a subset of data
    ! also test dimensions what happens if chebfi%X is smaller
    !call xgBlock_setBlock(chebfi%xXColsRows , xX_active , chebfi%spacedim, n_active, f
    !call xgBlock_setBlock(chebfi%xAXColsRows, xAX_active, 
    !call xgBlock_setBlock(chebfi%xBXColsRows, xBX_active, 

    ABI_FREE(mask)
    ABI_FREE(is_locked)
    call xg_free(T_swap)

end subroutine chebfi_lockConvergedVectors
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_deflateWrtLocked
!! NAME
!! chebfi_deflateWrtLocked

  subroutine chebfi_deflateWrtLocked(chebfi, X_lock, AX_lock, BX_lock, n_locked)

    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(in) :: X_lock, AX_lock, BX_lock
    integer, intent(in) :: n_locked

    integer :: previousBlock
    integer :: blockdim
    integer :: spacedim
    integer :: space_buf
    type(xg_t) :: buffer
    type(xgBlock_t) :: 

    blockdim = chebfi%blockdim
    spacedim = chebfi%spacedim
    previousBlock = (iblock-1)*lobpcg%blockdim

    n_active = chebfi%neigenpairs - n_locked
    call xgBlock_setBlock(chebfi%X      , X_active , rows(chebfi%X), n_active, fcol=n_locked+1)
    call xgBlock_setBlock(chebfi%AX%self, AX_active, rows(chebfi%X), n_active, fcol=n_locked+1)
    call xgBlock_setBlock(chebfi%BX%self, BX_active, rows(chebfi%X), n_active, fcol=n_locked+1)


    space_buf = space(var)
    if (space(var)==SPACE_CR) then
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,previousBlock,blockdim,comm=chebfi%spacecom,gpu_option=chebfi%gpu_option)

    ! buffer = BX0^T*X
    call xgBlock_gemm('t','n',1.0d0,chebfi%BX0,var,0.d0,buffer%self,comm=chebfi%spacecom)

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm('n','n',-1.0d0,chebfi%X0,buffer%self,1.0d0,var)

    call xg_free(buffer)

    ! todo AX and BX
    ! X = X - X0*buffer
    ! call xgBlock_gemm('n','n',-1.0d0,X_locked,buffer%self,1.0d0,chebfi%AX%self)
    ! AX = AX - X0*buffer
    ! call xgBlock_gemm('n','n',-1.0d0,X_locked,buffer%self,1.0d0,chebfi%AX%self)
    ! BX = BX - X0*buffer
    ! call xgBlock_gemm('n','n',-1.0d0,X_locked,buffer%self,1.0d0,chebfi%BX%self)


  end subroutine chebfi_deflateWrtLocked
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/bandPartition_setColsRowsActive
!! NAME
!! bandPartition_setColsRowsActive

  subroutine bandPartition_setColsRowsActive(bpart, X, n_locked)

      implicit none

      type(bandPartition_t), intent(inout) :: bpart
      type(xgBlock_t), intent(in) :: X
      integer, intent(in) :: n_locked
      integer :: nrows
  
  ! *********************************************************************

      nrows = rows(X)
      bpart%n_locked = 

      call xgBlock_setBlock(X, bpart%colsrows_locked, nrows, bpart%n_locked_mpi)
    
  end subroutine bandPartition_setColsRowsActive
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/bandPartition_setColsRowsActive
!! NAME
!! bandPartition_setColsRowsActive

  subroutine bandPartition_setColsRowsActive(bpart, X, n_locked)

      implicit none

      type(bandPartition_t), intent(inout) :: bpart
      type(xgBlock_t), intent(in) :: X
      integer, intent(in) :: n_locked
      integer :: nrows
  
  ! *********************************************************************

      nrows = rows(X)
      bpart%n_locked = 

      call xgBlock_setBlock(X, bpart%colsrows_locked, nrows, bpart%n_locked_bandpp)
    
  end subroutine bandPartition_setColsRowsActive
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/bandPartition_getActiveDistribution
!! NAME
!! bandPartition_getActiveDistribution
!! 
!! FUNCTION
!! Every process can treat between 1 and bandpp bands
!! with p processes we cover p x bandpp bands at maximum.
!! Solve problem find minimal p to cover n_active bands.
!! Constraints: use at least 1 process and at max nproc.
!! 
!! OUTPUT
!! bpart%n_active_bandpp    stores the result per mpi rank
!! bpart%rank_active        true if process treats active bands

  subroutine bandPartition_getActiveDistribution(bpart, n_locked)

      implicit none

      type(bandPartition_t), intent(inout) :: bpart
      type(xgBlock_t), intent(in) :: X
      integer, intent(in) :: n_locked
      integer :: nrows
  
  ! *********************************************************************

      nrows = rows(X)
      bpart%n_locked =

      min_p = ceiling(n_active / bandpp)
      p = minval(nproc, maxval(1, min_p))

      bpart%rank_active = ..
      bpart%n_active_bandpp = ..
      
    
  end subroutine bandPartition_getActiveDistribution
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/bandPartition_copyActiveDistribution
!! NAME
!! bandPartition_copyActiveDistribution
!! 

  subroutine bandPartition_copyActiveDistribution(bpart_in, bpart_out)

      implicit none

      type(bandPartition_t), intent(in   ) :: bpart_in
      type(bandPartition_t), intent(inout) :: bpart_out
  
  ! *********************************************************************

      bpart_out%rank_active = bpart_in%rank_active
      bpart_out%n_active_bandpp = bpart_in%n_active_bandpp
      bpart_out%comm_active = bpart_in%comm_active
    
  end subroutine bandPartition_copyActiveDistribution
!!***

!----------------------------------------------------------------------

!!****f* m_chebfi2/chebfi_constructActiveTransposer
!! NAME
!! chebfi_constructActiveTransposer

  subroutine chebfi_constructActiveTransposer(chebfi, bpart)

      implicit none

      type(chebfi_t), intent(inout) :: chebfi
      type(bandPartition_t), intent(inout) :: bpart

      integer :: n_active

  ! *********************************************************************

      n_active = chebfi%neigenpairs - bpart%n_locked

      ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
      call timab(tim_transpose,1,tsec)

      ! todo use n_active
      call xgTransposer_constructor(bpart%transposer_active,bpart%linalg_active,bpart%colsrows_active,&
          nspinor,STATE_LINALG,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0_fft,&
          gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit)

      bpart%transposer_active%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd

      ABI_NVTX_END_RANGE()
      call timab(tim_transpose,2,tsec)

  end subroutine chebfi_constructActiveTransposer
!!***

end module m_chebfi2
!!***

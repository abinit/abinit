!!****f* ABINIT/m_slice_cprj
!! NAME
!! m_slice_cprj
!!
!! FUNCTION
!! This module contains the types and routines used to apply the
!! Spectrum Slicing method (2021 implementation using xG abstraction layer)
!! It mainly defines a 'slice' datatypes and associated methods.
!!
!! COPYRIGHT
!! Copyright (C) 2023-2025 ABINIT group (LB,IML)
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

module m_slice_cprj

 use defs_basis
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_time, only : timab
 use m_sort, only: sort_dp

 use m_cgtools
 use m_xg
 use m_xgTransposer
 use m_xg_ortho_RR
 use m_xg_nonlop

 use m_xmpi
 use m_xomp
#ifdef HAVE_OPENMP
 use omp_lib
#endif

 implicit none

 private

!Several (private) parameters
!-------------------------------------------------

 integer, parameter :: tim_init         = 2171
 integer, parameter :: tim_free         = 2172
 integer, parameter :: tim_cprj         = 2173
 integer, parameter :: tim_invovl       = 2174
 integer, parameter :: tim_residu       = 2175
 integer, parameter :: tim_RR           = 2176
 integer, parameter :: tim_transpose    = 2177
 integer, parameter :: tim_RR_q         = 2178
 integer, parameter :: tim_postinvovl   = 2179
 integer, parameter :: tim_swap         = 2180
 integer, parameter :: tim_amp_f        = 2181
 integer, parameter :: tim_barrier      = 2182
 integer, parameter :: tim_copy         = 2183
 integer, parameter :: tim_ax_k         = 2184
 integer, parameter :: tim_ax_v         = 2185
 integer, parameter :: tim_ax_nl        = 2186
 integer, parameter :: tim_enl          = 2187
 integer, parameter :: tim_ortho        = 2188

!Public 'slice' datatype
!-------------------------------------------------

 type, public :: slice_t

   integer :: space                         
   integer :: space_cprj
   integer :: spacedim                      ! Space dimension for one vector
   integer :: cprjdim                       ! cprj dimension
   integer :: blockdim_cprj                 ! ncols of cprj for slice only
   integer :: all_blockdim_cprj             ! ncols of cprj for all vectors
   integer :: neigenpairs                   ! Number of eigen values/vectors we want
   integer :: ndeg_filter                   ! Degree of the polynomial filter as input
   integer :: nslice                        ! Number of spectral slices
   integer :: slicedim                      ! Number of eigen values/vectors we want in one slice
   integer :: paral_slice                   ! slice parallelization strategy (off)
   integer :: spectral_cut                  ! how to cut the eigenvalue spectral interval (off)
   real(dp) :: tolerance                    ! Tolerance on the residu to stop the minimization
   real(dp) :: ecut                         ! Ecut for Chebfi oracle
   real(dp) :: tolfilter                    ! Polynomial filter wanted amplification

   ! MPI-related
   integer :: bandpp 
   integer :: total_spacedim                ! Maybe not needed
   integer :: spacecom                      ! Communicator for MPI

   logical :: paw
   integer :: eigenProblem   !1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
   integer :: me_g0
   integer :: me_g0_fft

   type(xg_nonlop_t) :: xg_nonlop

   !ARRAYS for entire spectum
   type(xgBlock_t) :: AllX ! Block of initial and final solution
   type(xg_t) :: AllAX     ! space to save AX Hamiltonian application

   ! cprj for entire spectrum
   type(xgBlock_t) :: AllcprjX
   type(xg_t) :: Allcprj_work
   
   !ARRAYS on slice
   type(xg_t) :: X_SLICE
   type(xgBlock_t) :: X
   type(xgBlock_t) :: AX

   ! space to hold X_next, X_prev for Chebyshev recursion, slice only
   type(xg_t) :: X_NP
   type(xgBlock_t) :: X_next
   type(xgBlock_t) :: X_prev
   !SWAP POINTERS on slice, also for Chebyshev recursion
   type(xgBlock_t) :: X_swap
   type(xgBlock_t) :: AX_swap

   ! cprj for slice only
   type(xgBlock_t) :: cprjX ! pointer to AllcprjX
   type(xgBlock_t) :: cprj_work ! pointer to Allcprj_work
   type(xg_t) :: cprj_work2

   ! the following is independent of number of vectors, common to slice and spectrum
   type(xg_t) :: proj_work

   type(xg_t) :: X_PROBE

   type(xgBlock_t) :: eigenvalues

  end type slice_t

!Public methods associated to 'slice' datatype
!-------------------------------------------------
 public :: slice_init
 public :: slice_free
 public :: slice_memInfo
 public :: slice_run_cprj

 CONTAINS  !========================================================================================
!!***

!!****f* m_slice_cprj/slice_init
!! NAME
!! slice_init
!!
!! FUNCTION
!! Initialize a 'slice' datastructure.
!!
!! INPUTS
!!  bandpp= number of 'bands' handled by a processor
!!  eigenProblem= type of eigenpb: 1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
!!  me_g0= 1 if this processors treats G=0, 0 otherwise
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  ndeg_filter= polynomial degree of the Chebyshev filter (.i.e. number of H applications)
!!  space= defines in which space we are (columns, rows, etc.)
!!  spacecom= MPI communicator
!!  spacedim= space dimension for one vector
!!  paw= flag. TRUE if current calculation ses the PAW approach
!!  ecut= plane-wave cut-off energy
!!  tolerance= tolerance criterion on the residu to stop the minimization
!!  nslice= number of slices
!!  tolfilter= tolerance for filter amplification (a priori)
!!  paral_slice= slice parallelization
!!  spectral_cut= how to cut the eigenvalue spectrum (a priori)
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_init(slice,neigenpairs,spacedim,cprjdim,tolerance,ecut,bandpp, &
                      ndeg_filter,space,space_cprj,eigenProblem,spacecom,me_g0,paw,&
                      nslice,tolfilter,paral_slice,spectral_cut,&
                      xg_nonlop,me_g0_fft)

!Arguments ------------------------------------
 integer          , intent(in   ) :: bandpp
 integer          , intent(in   ) :: eigenProblem
 integer          , intent(in   ) :: me_g0
 integer          , intent(in   ) :: me_g0_fft
 integer          , intent(in   ) :: neigenpairs
 integer          , intent(in   ) :: ndeg_filter
 integer          , intent(in   ) :: space
 integer          , intent(in   ) :: space_cprj
 integer          , intent(in   ) :: spacecom
 integer          , intent(in   ) :: spacedim
 integer          , intent(in   ) :: cprjdim
 integer          , intent(in   ) :: nslice
 integer          , intent(in   ) :: paral_slice
 integer          , intent(in   ) :: spectral_cut
 logical          , intent(in   ) :: paw
 real(dp)         , intent(in   ) :: ecut
 real(dp)         , intent(in   ) :: tolerance
 real(dp)         , intent(in   ) :: tolfilter
 type(xg_nonlop_t), intent(in   ) :: xg_nonlop
 type(slice_t)   , intent(inout) :: slice

!Local variables-------------------------------
 real(dp)                    :: tsec(2)

! *********************************************************************

 call timab(tim_init,1,tsec)

 slice%space         = space
 slice%space_cprj    = space_cprj
 slice%neigenpairs   = neigenpairs
 slice%spacedim      = spacedim
 slice%bandpp        = bandpp
 slice%spacecom      = spacecom
 slice%cprjdim       = cprjdim
 if (tolerance > 0.0) then
   slice%tolerance = tolerance
 else
   slice%tolerance = 1.0e-20
 end if
 slice%ecut          = ecut
 slice%ndeg_filter   = ndeg_filter
 slice%eigenProblem  = eigenProblem
 slice%me_g0         = me_g0
 slice%me_g0_fft     = me_g0_fft
 slice%paw           = paw
 slice%xg_nonlop     = xg_nonlop

 ! slice specific
 slice%nslice        = nslice
 slice%tolfilter     = tolfilter
 slice%paral_slice   = paral_slice
 slice%spectral_cut  = spectral_cut

 !!!!! HARDCODED !!!!!
 slice%slicedim = 20 ! hard-coded
 !!!!!!!!!!!!!!!!!!!!!
 
 ABI_CHECK(slice%slicedim == neigenpairs/2, "must change hard-coded quantity to continue")

 ! cprj depending on number of bands specific
 !slice%blockdim_cprj = slice%slicedim*xg_nonlop%nspinor
 slice%blockdim_cprj = slice%neigenpairs*xg_nonlop%nspinor
 slice%all_blockdim_cprj = neigenpairs*xg_nonlop%nspinor

 call slice_allocateAll(slice)

 call timab(tim_init,2,tsec)

end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_allocateAll
!! NAME
!! slice_allocateAll
!!
!! FUNCTION
!! Allocate all memory spaces in a 'slice' datastructure.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_allocateAll(slice)

 ! Arguments ------------------------------------
 type(slice_t)  , intent(inout) :: slice

 ! Local variables-------------------------------
 ! scalars
 integer  :: neigenpairs
 integer  :: space,space_cprj
 integer  :: spacedim
 integer  :: slicedim
 integer  :: total_spacedim, ierr
 integer  :: nspinor

! *********************************************************************

 space       = slice%space
 space_cprj  = slice%space_cprj
 spacedim    = slice%spacedim
 !slicedim    = slice%slicedim
 slicedim    = slice%neigenpairs
 neigenpairs = slice%neigenpairs
 nspinor     = slice%xg_nonlop%nspinor

 call slice_free(slice)

 total_spacedim = spacedim
 call xmpi_sum(total_spacedim,slice%spacecom,ierr)
 slice%total_spacedim = total_spacedim

 ! transposed arrays (for slice only)
 call xg_init(slice%X_NP,space,total_spacedim,2*slicedim,xmpi_comm_self,me_g0=slice%me_g0_fft)
 call xg_setBlock(slice%X_NP,slice%X_next,total_spacedim,slicedim)
 call xg_setBlock(slice%X_NP,slice%X_prev,total_spacedim,slicedim,fcol=slicedim+1)

 call xg_init(slice%X_SLICE,space,total_spacedim,2*slicedim,xmpi_comm_self,me_g0=slice%me_g0_fft)
 call xg_setBlock(slice%X_SLICE,slice%X,total_spacedim,slicedim)
 call xg_setBlock(slice%X_SLICE,slice%AX,total_spacedim,slicedim,fcol=slicedim+1)
 ! note: the last two are necessary to set space and me_g0 for slice%X and slice%AX
 
 call xg_init(slice%X_PROBE,space,total_spacedim,slicedim,xmpi_comm_self,me_g0=slice%me_g0_fft)

 ! cprj workspaces (for entire spectrum)
 call xg_init(slice%AllAX,space,spacedim,neigenpairs,slice%spacecom,me_g0=slice%me_g0)
 call xg_init(slice%Allcprj_work,space_cprj,slice%cprjdim,slice%all_blockdim_cprj,slice%spacecom)

 ! cprj workspaces (for slice only) 
 call xg_init(slice%cprj_work2,space_cprj,slice%cprjdim,slice%blockdim_cprj,slice%spacecom)
 
! cprj workspaces common to slices and entire spectrum
 call xg_init(slice%proj_work,space,slice%xg_nonlop%max_npw_k,slice%xg_nonlop%cprjdim,slice%spacecom,me_g0=slice%me_g0)

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_free
!! NAME
!! slice_free
!!
!! FUNCTION
!! Destroy a 'slice' datastructure.
!!
!! INPUTS
!!
!! OUTPUT
!!  arraymem(2)= memory information
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_free(slice)

!Arguments ------------------------------------
 type(slice_t) , intent(inout) :: slice

! *********************************************************************

 call xg_free(slice%X_NP)
 call xg_free(slice%X_SLICE)
 call xg_free(slice%X_PROBE)
 call xg_free(slice%Allcprj_work)
 call xg_free(slice%cprj_work2)
 call xg_free(slice%proj_work)
 call xg_free(slice%AllAX)
 call xg_free(slice%Allcprj_work)

end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_memInfo
!! NAME
!! slice_memInfo
!!
!! FUNCTION
!! Provides memory information about a 'slice' datastructure.
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
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

function slice_memInfo(neigenpairs,spacedim,space,total_spacedim,bandpp) result(arraymem)

!Arguments ------------------------------------
 integer, intent(in   ) :: bandpp
 integer, intent(in   ) :: neigenpairs
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
!slice_rayleighRitz function variables
 real(dp) :: memA_und_X
 real(dp) :: memB_und_X
 real(dp) :: memEigenvalues
 real(dp) :: cplx
!arrays
 real(dp) :: arraymem(2)

! *********************************************************************
 cplx = 1
 if ( space == SPACE_C ) cplx = 2 !for now only complex

 !Permanent in slice
 memX = cplx * kind(1.d0) * spacedim * neigenpairs

! if (paral_kgb == 0) then
!   memX_next = cplx * kind(1.d0) * spacedim * neigenpairs
!   memX_prev = cplx * kind(1.d0) * spacedim * neigenpairs
! else
   memX_next = cplx * kind(1.d0) * total_spacedim * bandpp
   memX_prev = cplx * kind(1.d0) * total_spacedim * bandpp
! end if

 memAX = cplx * kind(1.d0) * spacedim * neigenpairs
 memBX = cplx * kind(1.d0) * spacedim * neigenpairs

 !Transposer colrow array
! if (paral_kgb == 1) then
   memX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
   memAX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
   memBX_CR = cplx * kind(1.d0) * total_spacedim * bandpp
! else
!   memX_CR = 0
!   memAX_CR = 0
!   memBX_CR = 0
! end if

 !slice_rayleighRitz function variables
 memA_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
 memB_und_X = cplx * kind(1.d0) * neigenpairs * neigenpairs
 memEigenvalues = kind(1.d0) * neigenpairs

 arraymem(1) = memX + memX_next + memX_prev + &
               memAX + memBX + memX_CR + memAX_CR + memBX_CR
 arraymem(2) = memA_und_X + memB_und_X + memEigenvalues

end function slice_memInfo
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_run_cprj
!! NAME
!! slice_run_cprj
!!
!! FUNCTION
!! Apply the Spectrum Slicing algorithm on a set of vectors.
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
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_run_cprj(slice,X0,cprjX0,getAX,kin,eigen,occ,residu,enl,nspinor)

!Arguments ------------------------------------
 type(slice_t) , intent(inout) :: slice
 integer,         intent(in)    :: nspinor
 type(xgBlock_t), intent(inout) :: X0
 type(xgBlock_t), intent(inout) :: cprjX0
 type(xgBlock_t), intent(inout) :: eigen
 type(xgBlock_t), intent(in)    :: occ
 type(xgBlock_t), intent(inout) :: residu
 type(xgBlock_t), intent(inout) :: enl
 type(xgBlock_t), intent(in   ) :: kin
 interface
   subroutine getAX(X,AX)
     use m_xg, only : xgBlock_t
     type(xgBlock_t), intent(inout) :: X
     type(xgBlock_t), intent(inout) :: AX
   end subroutine getAX
 end interface

!Local variables-------------------------------
!scalars
 integer :: spacedim
 integer :: space_res
 integer :: neigenpairs
 integer :: nslice
 integer :: spectral_cut
 integer :: paral_slice
 integer :: ndeg_filter,ndeg_filter_max
 integer :: ideg, ierr
 integer :: iband, islice
 integer :: slicedim
 integer :: shift_x,shift_cprj
 integer :: niter, max_niter_restart
 integer :: nband_slice, nband_slice_buf
 integer :: count_mask
 integer :: icount
 real(dp) :: tolerance
 real(dp) :: tolfilter
 real(dp) :: maxeig, maxeig_global
 real(dp) :: mineig, mineig_global
 real(dp) :: lambda_minus
 real(dp) :: lambda_plus
 real(dp) :: one_over_r
 real(dp) :: two_over_r
 real(dp) :: center
 real(dp) :: radius
 real(dp) :: ls, us, cdeg, mu, damp
 real(dp) :: max_dist2
 real(dp) :: prev_mineig
 real(dp) :: tol12 = 1.0e-12
 real(dp) :: tol_probe = 1.0d0 ! this depends on the residual
 type(xg_t) :: Xsum
 type(xg_t) :: DivResults
 type(xg_t) :: dist1, dist2, dist12
 type(xgBlock_t) :: DivResults_part
 type(xgBlock_t) :: X_prev
 type(xgBlock_t) :: cprjX_prev
 type(xgBlock_t) :: cprj_work_prev
 type(xgBlock_t) :: eig_part
 type(xgBlock_t) :: X_col, AX_col
 type(xg_t) :: X_kept, AX_kept
 type(xgBlock_t) :: X_kept_col, AX_kept_col
 type(xg_t) :: X_part
 type(xg_t) :: AX_part
 type(xg_t) :: cprjX_part
 integer,parameter :: gpu_option=ABI_GPU_DISABLED
!arrays
 real(dp) :: tsec(2)
 integer, allocatable :: permute_cols(:)
 integer, allocatable :: sorted_idx(:) ! same as permute_cols but used elsewhere
 real(dp), allocatable :: rayleigh_quotients(:)
 real(dp), pointer :: probe(:)
 real(dp), pointer :: lambda_apost(:)
 real(dp), pointer :: dist2_array(:) => null()
 real(dp), pointer :: theta_(:,:) => null()
 !Pointers similar to old Chebfi
 integer,allocatable :: ndeg_filter_slice(:) !Slice variable
 integer,allocatable :: ndeg_filter_bands(:) !Oracle variable
 type(xg_nonlop_t) :: xg_nonlop

! *********************************************************************

 ! ITEST 
 write(901,*) 'inside slice_run'
 write(901,*) 'nslice=', slice%nslice
 flush(901)
 ! ITEST

 ! Warning; the entire code assumes this for simplicity and debugging purposes
 ABI_CHECK(slice%bandpp == slice%neigenpairs, "slice_cprj not implemented in MPI")

 ! Read scalar variables
 nslice = slice%nslice
 slicedim = slice%slicedim
 tolfilter = slice%tolfilter
 spectral_cut = slice%spectral_cut
 paral_slice = slice%paral_slice
 spacedim = slice%spacedim
 neigenpairs = slice%neigenpairs
 ndeg_filter = slice%ndeg_filter
 tolerance = slice%tolerance
 xg_nonlop = slice%xg_nonlop

 ! Set space of results (with symmetry or not)..
 if (slice%space==SPACE_C) then
   space_res = SPACE_C
 else if (slice%space==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if

 ! Allocations
 ! DivResults stores Rayleigh quotients
 call xg_init(DivResults, space_res, neigenpairs, 1)
 ABI_MALLOC(ndeg_filter_slice,(nslice))
 ABI_MALLOC(permute_cols, (neigenpairs))
 ABI_MALLOC(rayleigh_quotients, (neigenpairs))

 ! Set initial vectors from input guess
 slice%eigenvalues = eigen

 ! workspace named "All" will store input and output solution
 slice%AllX = X0
 slice%AllcprjX = cprjX0

 ! Compute cprjX for all X (in colsrows)
 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%AllX,slice%AllcprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

 ! ITEST
 write(901,*) 'slice%X0 before filter=', xgBlock_getid(X0)
 write(901,*) 'slice%AllX before filter=', xgBlock_getid(slice%AllX) 
 flush(901)
 ! ITEST

 !Compute A * Psi
 call timab(tim_AX_v,1,tsec)
 call getAX(slice%AllX,slice%AllAX%self)
 call timab(tim_AX_v,2,tsec)
 call timab(tim_AX_k,1,tsec)
 call xgBlock_add_diag(slice%AllX,kin,nspinor,slice%AllAX%self)
 call timab(tim_AX_k,2,tsec)
 call timab(tim_AX_nl,1,tsec)
 call xg_nonlop_getHX(xg_nonlop,slice%AllAX%self,slice%AllcprjX,slice%Allcprj_work%self,slice%proj_work%self)
 call timab(tim_AX_nl,2,tsec)

 ! ITEST
 write(901,*) 'slice%AllX init guess=', xgBlock_getid(slice%AllX) 
 write(901,*) 'slice%AllAX init guess=', xgBlock_getid(slice%AllAX%self) 
 flush(901)
 ! ITEST

 ! B-orthonormalize X and AX for all bands(assuming linalg)
 !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
 !   gpu_option,AX=slice%AllAX%self)

 ! ITEST
 !write(901,*) 'slice%AllX after Bortho no1=', xgBlock_getid(slice%AllX) 
 !write(901,*) 'slice%AllAX after Bortho no1=', xgBlock_getid(slice%AllAX%self) 
 !flush(901)
 ! ITEST

 ! B-orthonormalize X and AX for all bands(assuming linalg)
 !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
 !    gpu_option,AX=slice%AllAX%self)

 ! ITEST
 !write(901,*) 'slice%AllX after Bortho no2=', xgBlock_getid(slice%AllX) 
 !write(901,*) 'slice%AllAX after Bortho no2=', xgBlock_getid(slice%AllAX%self) 
 !flush(901)
 ! ITEST

 ! Compute Rayleigh quotients for every band
 call timab(tim_RR_q, 1, tsec)
 call slice_rayleighRitzQuotients(slice, maxeig, mineig, DivResults%self)
 call xmpi_max(maxeig,maxeig_global,slice%spacecom,ierr)
 call xmpi_min(mineig,mineig_global,slice%spacecom,ierr)
 call timab(tim_RR_q, 2, tsec)

 ! Recover column indices of sorted Rayleigh quotients in increasing order
 call xgBlock_reverseMap(DivResults%self, theta_, rows=1, cols=neigenpairs)
 rayleigh_quotients(1:neigenpairs) = theta_(1,1:neigenpairs)
 permute_cols(1:neigenpairs) = (/ (iband, iband=1,neigenpairs) /)
 call sort_dp(neigenpairs, rayleigh_quotients, permute_cols, tol12)
 ! store order into theta_
 theta_(1,1:neigenpairs) = rayleigh_quotients(1:neigenpairs)

 ! ITEST
 write(901,*) 'rayleigh quotients='
 call xgBlock_print(DivResults%self, 901)
 flush(901)
 ! ITEST

 ! Permute columns 1,..,neigenpairs according to order
 call xgBlock_permuteCols(slice%AllX, slice%spacedim, neigenpairs, permute_cols)
 call xgBlock_permuteCols(slice%AllAX%self, slice%spacedim, neigenpairs, permute_cols)
 ! ITEST
 write(901,*) 'slice%AllX after perm=', xgBlock_getid(slice%AllX) 
 write(901,*) 'slice%AllAX after perm=', xgBlock_getid(slice%AllAX%self) 
 flush(901)
 ! ITEST 

 ! TODO insert probe here
 ! 1) Sequential loop on slices. Start with first slice. The first slice
 !    only needs an upper bound. This is the rayleigh quotient (overestimation)
 !    at the i-th eigenvalue, where i is the maximum slice width.
 
 ! What we don't want to do; 
 ! Apply RR to 1,..,nband and choose nband_slice
 ! What we expect;
 ! Find vectors of number nband_slice_pad (<nband_slice) to apply RR.
 ! What we hope;
 ! That the gain from applying RR to fewer vectors is less than the loss from 
 ! filtering more vectors, when compared to Chebyshev where everything is on nband.
 islice = 1
 nband_slice_buf = 15 ! maximum number of vectors we afford to RR
 
 ! ITEST
 write(901,*)
 write(901,*) '====================Slice=================', islice
 write(901,*) 'Filter maximum buffer of size=', nband_slice_buf
 flush(901)
 ! ITEST

 ! Compute |X|^2 colwise L2-norm (before any filter)
 call xg_init(dist1,SPACE_R,neigenpairs,1)
 call xgBlock_colwiseNorm2(slice%AllX,dist1%self,comm_loc=xmpi_comm_null)
 ! preallocate, will be updated during ideg iterations
 call xg_init(dist2,SPACE_R,neigenpairs,1)

 ! ITEST
 write(901,*) 'probe dist1='
 call xgBlock_print(dist1%self,901)
 flush(901)
 ! ITEST

 ! Start filter up to degree
 !lambda_minus = rayleigh_quotients(nband_slice_buf)
 lambda_minus = 0.51404d0
 lambda_plus = slice%ecut
 center = (lambda_plus + lambda_minus)*0.5
 radius = (lambda_plus - lambda_minus)*0.5
 one_over_r = 1.0/radius
 two_over_r = 2.0/radius
 
 ! ITEST
 write(901,*) 'unwanted part of the spectrum=', lambda_minus, lambda_plus
 flush(901)
 ! ITEST

 ! Temporarily all slice pointers coincide with entire pointers
 !slice%X = slice%AllX
 !slice%AX = slice%AllAX%self
 call xgBlock_copy(slice%AllX, slice%X)
 call xgBlock_copy(slice%AllAX%self, slice%AX)
 slice%cprjX = slice%AllcprjX
 slice%cprj_work = slice%Allcprj_work%self

 ! ITEST
 write(901,*) 
 write(901,*) 'X id at slice2 init', xgBlock_getid(slice%X)
 write(901,*) 'AX id at slice2 init', xgBlock_getid(slice%AX)
 write(901,*) 
 flush(901)
 ! ITEST
 
 ! Start loop on degree
 ndeg_filter = slice%ndeg_filter ! take the degree from abi
 do ideg = 0, ndeg_filter - 1

    ! ITEST
    write(901,*) 'polynomial degree=', ideg
    flush(901)
    ! ITEST

    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop,slice%AX,slice%cprjX,slice%proj_work%self)
    call timab(tim_cprj,2,tsec)

    call slice_computeNextOrderChebfiPolynom(slice, ideg, center, one_over_r, two_over_r)

    call timab(tim_swap,1,tsec)
    call slice_swapInnerBuffers(slice, slice%total_spacedim, neigenpairs)
    call timab(tim_swap,2,tsec)

    ! Amplify, use for probe only for degree ideg
    call xgBlock_copy(slice%X, slice%X_PROBE%self)
    !call slice_ampfactorProbe(slice, DivResults%self, lambda_minus, lambda_plus, ideg+1)
    
    ! ######### Compute probe
    ! dist2 = |X_PROBE|^2 colwise L2-norm
    call xgBlock_colwiseNorm2(slice%X_PROBE%self, dist2%self, max_dist2, comm_loc=xmpi_comm_null)
    !call xgBlock_scale(dist2%self, 1/max_dist2, 1)

    ! ITEST
    write(901,*) 'maximum dist2=', max_dist2
    write(901,*) 'probe dist2(scaled by max)='
    call xgBlock_print(dist2%self,901)
    flush(901)
    ! ITEST
            
    !A * Psi for next iteration
    call timab(tim_AX_v,1,tsec)
    call getAX(slice%X,slice%AX)
    call timab(tim_AX_v,2,tsec)
    call timab(tim_AX_k,1,tsec)
    call xgBlock_add_diag(slice%X,kin,nspinor,slice%AX)
    call timab(tim_AX_k,2,tsec)
    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
    call timab(tim_cprj,2,tsec)
    call timab(tim_AX_nl,1,tsec)
    call xg_nonlop_getHX(xg_nonlop,slice%AX,slice%cprjX,slice%cprj_work,slice%proj_work%self)
    call timab(tim_AX_nl,2,tsec)

 end do

 ! Amplify to finalize filter application
 call timab(tim_amp_f,1,tsec)
 ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
 ndeg_filter_bands(:) = ndeg_filter
 call slice_ampfactor(slice, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)
 ABI_FREE(ndeg_filter_bands)
 call timab(tim_amp_f,2,tsec)
 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

 ! Use dist2%self as probe to select indices
 ! first count
 call xgBlock_reverseMap_1d(dist2%self,probe)
 write(901,*) 
 write(901,*) 'probe values used to discard indices=', probe(:)
 flush(901)

 count_mask = count(probe > tol_probe, 1)

 write(901,*) 
 write(901,*) 'Keep count_mask(should be 15)= out of neigenpairs=', count_mask, neigenpairs
 flush(901)

 if (nband_slice_buf < count_mask) then
     write(901,*) 'Number of kept is more than initial guess! Initial guess missed eigenpairs'
     flush(901)
 end if

 ! Procedure for keeping selected to new workspaces
 call xg_init(X_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
 call xg_init(AX_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
 icount = 1
 do iband=1,neigenpairs
    if (probe(iband) > tol_probe) then
        call xgBlock_setBlock(X_kept%self, X_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(AX_kept%self, AX_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(slice%X, X_col, slice%total_spacedim, 1, fcol=iband)
        call xgBlock_setBlock(slice%AX, AX_col, slice%total_spacedim, 1, fcol=iband)
        call xgBlock_copy(X_col, X_kept_col)
        call xgBlock_copy(AX_col, AX_kept_col)
        icount = icount + 1
    end if
 end do

 ! reset pointers to temporary (kept)
 slice%X = X_kept%self
 slice%AX = AX_kept%self
 ! content is not important, but dimensions 
 call xgBlock_setBlock(slice%AllcprjX, slice%cprjX, slice%cprjdim, count_mask*nspinor)

 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

 ! Apply Rayleigh Ritz on slice (refinement)
 call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,slice%eigenvalues,&
     slice%blockdim_cprj,ierr,0,tim_RR,ABI_GPU_DISABLED,solve_ax_bx=.true.)

 if ( ierr /= 0 ) then
    ABI_BUG("TEST slice1 RayleighRitz did not work")
 end if

 write(901,*) 'converged eigenval1='
 call xgBlock_print(slice%eigenvalues, 901)
 flush(901)

 ! reset and free
 call xg_free(dist1)
 call xg_free(dist2)
 call xg_free(X_kept)
 call xg_free(AX_kept)

 ! continue testing .. slice 2
 islice = 2
 !nband_slice_buf = 15

 lambda_minus = 0.51404
 lambda_plus = 1.25348
 prev_mineig = -1.14360
 center = (slice%ecut + prev_mineig)*0.5
 radius = (slice%ecut - prev_mineig)*0.5 
 ls = (lambda_minus - center) / radius
 us = (lambda_plus - center) / radius
 one_over_r = 1.0/radius
 two_over_r = 2.0/radius

 ! ITEST
 write(901,*)
 write(901,*) '====================Slice=================', islice
 write(901,*) 'Probe for vector space of wanted interval=', lambda_minus, lambda_plus
 flush(901)
 ! ITEST

 ! Reset pointers for dimensions of nband
 call xg_setBlock(slice%X_SLICE,slice%X,slice%total_spacedim,neigenpairs)
 call xg_setBlock(slice%X_SLICE,slice%AX,slice%total_spacedim,neigenpairs,fcol=neigenpairs+1)
 ! Init workspaces
 call xgBlock_copy(slice%AllX, slice%X)
 call xgBlock_copy(slice%AllAX%self, slice%AX)
 ! the two following ones will be recomputed so whatever
 slice%cprjX = slice%AllcprjX
 slice%cprj_work = slice%Allcprj_work%self

 ! ITEST
 write(901,*) 
 write(901,*) 'X id at slice2 init', xgBlock_getid(slice%X)
 write(901,*) 'AX id at slice2 init', xgBlock_getid(slice%AX)
 write(901,*) 
 flush(901)
 ! ITEST

 ! preallocate, will be updated during ideg iterations
 call xg_init(dist2,SPACE_R,neigenpairs,1)
 
 ! Initialize Chebyshev expansion of indicator function of order ndeg_filter
 ndeg_filter = 60
 call xg_init(Xsum,slice%space,slice%total_spacedim,neigenpairs,slice%spacecom,gpu_option=gpu_option)
 cdeg = Pi/(ndeg_filter+2)
 mu = 1.d0/Pi*(ACOS(ls)-ACOS(us))
 damp = 1.d0 ! Jackson damping
 call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)

 ! Start loop on degree
 do ideg = 0, ndeg_filter - 1

   ! ITEST
   write(901,*) 'polynomial degree=', ideg
   flush(901)
   ! ITEST

   call timab(tim_cprj,1,tsec)
   call xg_nonlop_getcprj(xg_nonlop,slice%AX,slice%cprjX,slice%proj_work%self)
   call timab(tim_cprj,2,tsec)

   call slice_computeNextOrderChebfiPolynom(slice, ideg, center, one_over_r, two_over_r)

   call timab(tim_swap,1,tsec)
   call slice_swapInnerBuffers(slice, slice%total_spacedim, neigenpairs)
   call timab(tim_swap,2,tsec)

   ! Accumulate X with weight in Xsum for bandpass filters
   mu = 2/Pi * (SIN((ideg+1)*ACOS(ls)) - SIN((ideg+1)*ACOS(us)))/(ideg+1)
   damp = ((1 - (ideg+1)/(ndeg_filter+2))*SIN(cdeg)*COS((ideg+1)*cdeg) + &
            1/(ndeg_filter+2)*COS(cdeg)*SIN((ideg+1)*cdeg))/SIN(cdeg)
   call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)
   
   ! ######### Compute probe
   ! dist2 = |Xsum|^2 colwise L2-norm
   call xgBlock_colwiseNorm2(Xsum%self, dist2%self, max_dist2, comm_loc=xmpi_comm_null)
   call xgBlock_scale(dist2%self, 1/max_dist2, 1)
   
   ! ITEST
   write(901,*) 'maximum dist2=', max_dist2
   write(901,*) 'probe dist2(scaled by max)='
   !write(901,*) 'probe dist2='
   call xgBlock_print(dist2%self,901)
   flush(901)
   ! ITEST
   
   ! Count kept: n-maximum values until ordered indices are stabilized
   call xgBlock_reverseMap_1d(dist2%self, probe)
   ABI_MALLOC(sorted_idx, (neigenpairs))
   sorted_idx = (/ (iband, iband=1,neigenpairs) /)
   call sort_dp(neigenpairs, probe, sorted_idx, tol12)

   ! ITEST
   write(901,*) 'sorted indices=', sorted_idx(neigenpairs:1:-1)
   flush(901)
   ! ITEST

   ABI_FREE(sorted_idx)

   ! store final expansion Xsum to X
   if (ideg==ndeg_filter - 1) then 
     call xgBlock_copy(Xsum%self, slice%X)
   end if

   !A * Psi
   call timab(tim_AX_v,1,tsec)
   call getAX(slice%X,slice%AX)
   call timab(tim_AX_v,2,tsec)
   call timab(tim_AX_k,1,tsec)
   call xgBlock_add_diag(slice%X,kin,nspinor,slice%AX)
   call timab(tim_AX_k,2,tsec)
   call timab(tim_cprj,1,tsec)
   call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
   call timab(tim_cprj,2,tsec)
   call timab(tim_AX_nl,1,tsec)
   call xg_nonlop_getHX(xg_nonlop,slice%AX,slice%cprjX,slice%cprj_work,slice%proj_work%self)
   call timab(tim_AX_nl,2,tsec)

 end do

 ! Use dist2%self as probe to select indices
 ! first count
 call xgBlock_reverseMap_1d(dist2%self,probe)
 write(901,*) 
 write(901,*) 'probe values used to discard indices=', probe(:)
 flush(901)

 tol_probe = 0.25d0
 count_mask = count(probe > tol_probe, 1)

 write(901,*) 
 write(901,*) 'Keep count_mask(should be 15)= out of neigenpairs=', count_mask, neigenpairs
 flush(901)

 if (nband_slice_buf < count_mask) then
     write(901,*) 'Number of kept is more than initial guess! Initial guess missed eigenpairs'
     flush(901)
 end if

 ! Procedure for keeping selected to new workspaces
 call xg_init(X_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
 call xg_init(AX_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
 icount = 1
 do iband=1,neigenpairs
    if (probe(iband) > tol_probe) then
        call xgBlock_setBlock(X_kept%self, X_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(AX_kept%self, AX_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(slice%X, X_col, slice%total_spacedim, 1, fcol=iband)
        call xgBlock_setBlock(slice%AX, AX_col, slice%total_spacedim, 1, fcol=iband)
        call xgBlock_copy(X_col, X_kept_col)
        call xgBlock_copy(AX_col, AX_kept_col)
        icount = icount + 1
    end if
 end do

 ! reset pointers to temporary (kept)
 slice%X = X_kept%self
 slice%AX = AX_kept%self
 ! content is not important, but dimensions 
 call xgBlock_setBlock(slice%AllcprjX, slice%cprjX, slice%cprjdim, count_mask*nspinor)

 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

 ! Apply Rayleigh Ritz on slice (refinement)
 call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,slice%eigenvalues,&
     slice%blockdim_cprj,ierr,0,tim_RR,ABI_GPU_DISABLED,solve_ax_bx=.true.)

 if ( ierr /= 0 ) then
    ABI_BUG("TEST slice2 RayleighRitz did not work")
 end if

 ! Count how many eigenvalues converged in wanted interval
 call xgBlock_reverseMap_1d(slice%eigenvalues, lambda_apost)
 count_mask = count(lambda_minus < lambda_apost .and. lambda_apost < lambda_plus, 1)

 ! ITEST
 write(901,*) 'What we want: 15 vectors in ', lambda_minus, lambda_plus
 write(901,*) 'what we get: converged eigenval2 of number=', count_mask
 call xgBlock_print(slice%eigenvalues, 901)
 flush(901)
 ! ITEST

 ! Compute H-eSX
 if (slice%paw) then
    call timab(tim_AX_nl,1,tsec)
    call xg_nonlop_getHmeSX(xg_nonlop,slice%X,slice%cprjX,slice%AX,slice%eigenvalues,&
        slice%Allcprj_work%self,slice%cprj_work2%self,no_H=.True.)
    call timab(tim_AX_nl,2,tsec)
 end if

 ! Compute residual norm squared
 call timab(tim_residu, 1, tsec)
 if (.not.slice%paw) then
   call xgBlock_yxmax(slice%AX,slice%eigenvalues,slice%X)
 end if
 call xgBlock_colwiseNorm2(slice%AX, residu)
 call timab(tim_residu, 2, tsec)

 ! ITEST
 write(901,*) 'residu='
 call xgBlock_print(residu, 901)
 flush(901)
 ! ITEST

 call xg_free(Xsum)
 call xg_free(dist2)
 call xg_free(X_kept)
 call xg_free(AX_kept)

 !!!! ===============


 ! Compare this probe to rayleigh quotients. Actually count the number of
 ! rayleigh quotients in the slice interval. We expect to see that it is 
 ! not the same as the final converged number


 ! TODO oracle must be replaced by the bandpass filter degree computation...

 !!!!! HARDCODED !!!!!
 ndeg_filter_slice(1) = ndeg_filter
 ndeg_filter_slice(2) = 60 ! hard-coded
 !!!!!!!!!!!!!!!!!!!!!

 ! Sequential loop on slices
 do islice=2, nslice

    ! ITEST
    write(901,*)
    write(901,*) '====================Slice=================', islice
    write(901,*) 'Apply slicing for slicedim(fixed)=', slicedim 
    write(901,*) 'computing slicedim from probe...'
    flush(901)
    ! ITEST

    ! Init: at first all bands are into slice
    slicedim = neigenpairs

    ! In latest version, *_part is a separate object that does not impact AllX
    call xg_init(X_part, slice%space, spacedim, slicedim, slice%spacecom, me_g0=slice%me_g0_fft)
    call xg_init(AX_part, slice%space, spacedim, slicedim, slice%spacecom, me_g0=slice%me_g0_fft)
    call xg_init(cprjX_part, slice%space_cprj, slice%cprjdim, slicedim*nspinor, comm=slice%spacecom)

    ! old shifts assuming all slices have the same slicedim kept for reference
    shift_x = (islice-1)*slice%slicedim
    !shift_cprj = (islice-1)*slicedim*xg_nonlop%nspinor

    ! All bands are copied to current slice workspace
    call xgBlock_copy(slice%AllX, X_part%self)
    call xgBlock_copy(slice%AllAX%self, AX_part%self)
    call xgBlock_copy(slice%AllcprjX, cprjX_part%self)
    
    ! Working pointers to slice workspace
    slice%X = X_part%self
    slice%AX = AX_part%self
    slice%cprjX = cprjX_part%self
    call xgBlock_setBlock(slice%Allcprj_work%self, slice%cprj_work, slice%cprjdim, slice%blockdim_cprj) 
    call xgBlock_setBlock(DivResults%self, DivResults_part, slicedim, 1)
    call xgBlock_reshape(DivResults_part, 1, neigenpairs)

    ! Old version where AllX was modified
    ! Pointers to global objects common to all slices
    !call xgBlock_setBlock(slice%AllX, X_part%self, spacedim, slicedim, fcol=shift_x+1)
    !call xgBlock_setBlock(slice%AllAX%self, AX_part%self, spacedim, slicedim, fcol=shift_x+1)
    !call xgBlock_setBlock(slice%AllcprjX, cprjX_part%self, slice%cprjdim, slice%blockdim_cprj, fcol=shift_cprj+1)

    ! Get part of eigenvalues
    !call xgBlock_reshape(DivResults%self, 1, neigenpairs)
    !call xgBlock_setBlock(DivResults%self, DivResults_part, 1, slicedim, fcol=shift_x+1)
    !call xgBlock_reshape(DivResults%self, neigenpairs, 1)

    !! ######### Deactivate orthogonalization wrt Blocks
    !! Orthogonalize current X_part With Respect To previous blocks in B-basis
    !if ( islice > 1 ) then
    !    call xgBlock_setBlock(slice%AllX             , X_prev         , spacedim     , shift_x   , fcol=1)
    !    call xgBlock_setBlock(slice%AllcprjX         , cprjX_prev     , slice%cprjdim, shift_cprj, fcol=1)
    !    call xgBlock_setBlock(slice%Allcprj_work%self, cprj_work_prev , slice%cprjdim, shift_cprj, fcol=1)

    !    ! ITEST
    !    write(901,*) 'X prev=', xgBlock_getid(X_prev) 
    !    write(901,*) 'cprjX prev=', xgBlock_getid(cprjX_prev) 
    !    write(901,*) 'slice%X before orthoXwrt prev=', xgBlock_getid(slice%X) 
    !    write(901,*) 'slice%AX before orthoXwrt prev=', xgBlock_getid(slice%AX) 
    !    flush(901)
    !    ! ITEST

    !    call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)

    !    ! ITEST
    !    write(901,*) 'slice%X after orthoXwrt prev no1=', xgBlock_getid(slice%X) 
    !    write(901,*) 'slice%AX after orthoXwrt prev no1=', xgBlock_getid(slice%AX) 
    !    flush(901)
    !    ! ITEST
        
    !    !call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)

    !    ! ITEST
    !    !write(901,*) 'slice%X after orthoXwrt prev no2=', xgBlock_getid(slice%X) 
    !    !write(901,*) 'slice%AX after orthoXwrt prev no2=', xgBlock_getid(slice%AX) 
    !    !flush(901)
    !    ! ITEST

    !end if

    if (islice==1) then
        !*** hardcoded
        lambda_minus = 0.81501785305036
        ! unwanted spectrum
        !lambda_minus = rayleigh_quotients(slicedim) ! largest eigenvalue from current slice
        lambda_plus = slice%ecut
        center = (lambda_plus + lambda_minus)*0.5
        radius = (lambda_plus - lambda_minus)*0.5

        ! ITEST
        write(901,*) 'unwanted part of the spectrum=', lambda_minus, lambda_plus
        flush(901)
        ! ITEST

    else if (islice>1) then
 
        !write(901,*) 'rayleigh_quotients=', rayleigh_quotients
        !flush(901)

        !call timab(tim_RR_q, 1, tsec)
        !call slice_rayleighRitzQuotientsOnSlice(slice, maxeig, mineig, DivResults%self)
        !call timab(tim_RR_q, 2, tsec)
        
        !call xgBlock_reverseMap(DivResults%self, theta_, rows=1, cols=slicedim)
        !rayleigh_quotients(shift_x+1:shift_x+slicedim) = theta_(1,1:slicedim)
        
        !write(901,*) 'rayleigh_quotients=', rayleigh_quotients
        !flush(901)

        !*** hardcoded
        lambda_minus = 0.50
        lambda_plus = 0.9
        ! wanted second slice
        !lambda_minus = rayleigh_quotients(slicedim) ! largest eigenvalue from previous slice
        !lambda_plus = maxeig-0.3
        !center = (slice%ecut + rayleigh_quotients(1))*0.5
        !radius = (slice%ecut - rayleigh_quotients(1))*0.5
        center = (slice%ecut - 0.15)*0.5 !** hard-coded
        radius = (slice%ecut + 0.15)*0.5 !** hard-coded
        ls = (lambda_minus - center) / radius
        us = (lambda_plus - center) / radius

        ! ITEST
        write(901,*) 'wanted part of the spectrum [a,b)=', lambda_minus, lambda_plus
        write(901,*) 'scaled in [-1,1)=', ls, us
        flush(901)
        ! ITEST
    end if

    one_over_r = 1.0/radius
    two_over_r = 2.0/radius
    ndeg_filter = ndeg_filter_slice(islice)
 
    !!!!! HARDCODED !!!!!
    if (islice==1) then
        max_niter_restart=1
    else
        max_niter_restart = 1! number of inner iterations
    end if
    !!!!!!!!!!!!!!!!!!!!!

    do niter=1, max_niter_restart

        write(901,*) '%%%%%%%%%%%%%%%%%%%%%%%% inner restart=', niter
        flush(901)

        ! Initialize Chebyshev expansion
        if (islice>1) then
   
            ! # deactivate update for this version
            !if (niter>1) then
            !    ! adapt limits of wanted interval
            !    ! TODO lambda_plus should take into account occupancies. Put to last occupied state
            !    lambda_minus = rayleigh_quotients(shift_x+1) ! smalles eigenvalue in slice
            !    lambda_plus = rayleigh_quotients(shift_x+slicedim) ! largest eigenvalue in slice
            !    ls = (lambda_minus - center) / radius
            !    us = (lambda_plus - center) / radius
            !end if

            if (niter==1) then
                call xg_init(Xsum,slice%space,slice%total_spacedim,slicedim,slice%spacecom,&
                    gpu_option=gpu_option)
            else
                call xgBlock_zero(Xsum%self)
            end if
            cdeg = Pi/(ndeg_filter+2)
            mu = 1/Pi*(ACOS(ls)-ACOS(us))
            damp = 1.d0 ! Jackson damping
            call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)
        end if

        if (niter==1) then
            ! dist1 = |slice%X|^2 colwise L2-norm
            call xg_init(dist1,SPACE_R,slicedim,1)
            call xgBlock_colwiseNorm2(slice%X,dist1%self,comm_loc=xmpi_comm_null)
            ! preallocate
            call xg_init(dist2,SPACE_R,slicedim,1)
            call xg_init(dist12,SPACE_R,slicedim,1)
        end if

        ! ITEST
        write(901,*) 'restarted interval=', lambda_minus, lambda_plus
        flush(901)
        ! ITEST

        do ideg = 0, ndeg_filter - 1

            ! ITEST
            write(901,*) 'inner polynomial degree=', ideg
            flush(901)
            ! ITEST

            call timab(tim_cprj,1,tsec)
            call xg_nonlop_getcprj(xg_nonlop,slice%AX,slice%cprjX,slice%proj_work%self)
            call timab(tim_cprj,2,tsec)

            call slice_computeNextOrderChebfiPolynom(slice, ideg, center, one_over_r, two_over_r)

            call timab(tim_swap,1,tsec)
            call slice_swapInnerBuffers(slice, slice%total_spacedim, slicedim)
            call timab(tim_swap,2,tsec)

            ! Accumulate X with weight in Xsum for bandpass filters
            if (islice==2) then
                mu = 2/Pi * (SIN((ideg+1)*ACOS(ls)) - SIN((ideg+1)*ACOS(us)))/(ideg+1)
                damp = ((1 - (ideg+1)/(ndeg_filter+2))*SIN(cdeg)*COS((ideg+1)*cdeg) + &
                    1/(ndeg_filter+2)*COS(cdeg)*SIN((ideg+1)*cdeg))/SIN(cdeg)
                call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)
            end if

            ! ######### Compute probe
            ! dist2 = |Xsum%self|^2 colwise L2-norm
            if (islice==1) then
                call xgBlock_colwiseNorm2(slice%X, dist2%self, comm_loc=xmpi_comm_null)
            else
                call xgBlock_colwiseNorm2(Xsum%self,dist2%self,comm_loc=xmpi_comm_null)
            end if
            !call xgBlock_reverseMap_1d(dist2%self, dist2_array)
            !call xgBlock_colwiseMul(dist2%self, dist2_array) ! dist2 = dist2^2
            !call xgBlock_colwiseDivision(dist2%self,dist1%self,dist12%self) ! dist2/dist
            ! ITEST
            write(901,*) 'probe dist1='
            call xgBlock_print(dist1%self,901)
            write(901,*) 'probe dist2='
            call xgBlock_print(dist2%self,901)
            flush(901)
            ! ITEST
            
            if (islice==2) then
                ! store final expansion Xsum to X
                if (ideg==ndeg_filter - 1) then 
                    call xgBlock_copy(Xsum%self, slice%X)
                end if
            end if

            !A * Psi
            call timab(tim_AX_v,1,tsec)
            call getAX(slice%X,slice%AX)
            call timab(tim_AX_v,2,tsec)
            call timab(tim_AX_k,1,tsec)
            call xgBlock_add_diag(slice%X,kin,nspinor,slice%AX)
            call timab(tim_AX_k,2,tsec)
            call timab(tim_cprj,1,tsec)
            call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
            call timab(tim_cprj,2,tsec)
            call timab(tim_AX_nl,1,tsec)
            call xg_nonlop_getHX(xg_nonlop,slice%AX,slice%cprjX,slice%cprj_work,slice%proj_work%self)
            call timab(tim_AX_nl,2,tsec)

        end do

        ! Amplify for first slice only
        if (islice==1) then
            call timab(tim_amp_f,1,tsec)
            ABI_MALLOC(ndeg_filter_bands,(slicedim))
            ndeg_filter_bands(:) = ndeg_filter
            call slice_ampfactor(slice, DivResults_part, lambda_minus, lambda_plus, ndeg_filter_bands)
            ABI_FREE(ndeg_filter_bands)
            call timab(tim_amp_f,2,tsec)

            call timab(tim_cprj,1,tsec)
            call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
            call timab(tim_cprj,2,tsec)
        end if

        ! ITEST
        write(901,*) 'slice%AllX before ortho=', xgBlock_getid(slice%AllX) 
        write(901,*) 'slice%X before ortho=', xgBlock_getid(slice%X) 
        write(901,*) 'slice%AX before ortho=', xgBlock_getid(slice%AX) 
        write(901,*) 'slice%AllcprjX before ortho=', xgBlock_getid(slice%AllcprjX) 
        flush(901)
        ! ITEST

        ! B-orthonormalize X and AX for all bands(assuming linalg)
        call xg_Borthonormalize_cprj(xg_nonlop,slice%blockdim_cprj,slice%X,slice%cprjX,ierr,tim_ortho,&
            gpu_option,AX=slice%AX)

        ! ITEST
        write(901,*) 'slice%AllX after Bortho no1=', xgBlock_getid(slice%AllX) 
        write(901,*) 'slice%X after Bortho no1=', xgBlock_getid(slice%X) 
        write(901,*) 'slice%AX after Bortho no1=', xgBlock_getid(slice%AX) 
        write(901,*) 'slice%AllcprjX after Bortho no1=', xgBlock_getid(slice%AllcprjX) 
        flush(901)
        ! ITEST

        !! B-orthonormalize X and AX for all bands(assuming linalg)
        call xg_Borthonormalize_cprj(xg_nonlop,slice%blockdim_cprj,slice%X,slice%cprjX,ierr,tim_ortho,&
            gpu_option,AX=slice%AX)

        ! ITEST
        write(901,*) 'slice%AllX after Bortho no2=', xgBlock_getid(slice%AllX) 
        write(901,*) 'slice%X after Bortho no2=', xgBlock_getid(slice%X) 
        write(901,*) 'slice%AX after Bortho no2=', xgBlock_getid(slice%AX) 
        write(901,*) 'slice%AllcprjX after Bortho no2=', xgBlock_getid(slice%AllcprjX) 
        flush(901)
        ! ITEST

        !! ######### Deactivate orthogonalization wrt Blocks
        !! Orthogonalize current X_part With Respect To previous blocks in B-basis
        !if ( islice > 1 ) then
        !    call xgBlock_setBlock(slice%AllX             , X_prev         , spacedim     , shift_x   , fcol=1)
        !    call xgBlock_setBlock(slice%AllcprjX         , cprjX_prev     , slice%cprjdim, shift_cprj, fcol=1)
        !    call xgBlock_setBlock(slice%Allcprj_work%self, cprj_work_prev , slice%cprjdim, shift_cprj, fcol=1)

        !    ! ITEST
        !    write(901,*) 'X prev=', xgBlock_getid(X_prev) 
        !    write(901,*) 'cprjX prev=', xgBlock_getid(cprjX_prev) 
        !    write(901,*) 'slice%X before orthoXwrt prev=', xgBlock_getid(slice%X) 
        !    write(901,*) 'slice%AX before orthoXwrt prev=', xgBlock_getid(slice%AX) 
        !    flush(901)
        !    ! ITEST

        !    call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)

        !    ! ITEST
        !    write(901,*) 'slice%X after orthoXwrt prev no1=', xgBlock_getid(slice%X) 
        !    write(901,*) 'slice%AX after orthoXwrt prev no1=', xgBlock_getid(slice%AX) 
        !    flush(901)
        !    ! ITEST
        
        !    !call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)

        !    ! ITEST
        !    !write(901,*) 'slice%X after orthoXwrt prev no2=', xgBlock_getid(slice%X) 
        !    !write(901,*) 'slice%AX after orthoXwrt prev no2=', xgBlock_getid(slice%AX) 
        !    !flush(901)
        !    ! ITEST

        !end if
    
        ! Get part of eigenvalues
        !call xgBlock_reshape(slice%eigenvalues, 1, neigenpairs)
        !call xgBlock_setBlock(slice%eigenvalues, eig_part, 1, slicedim, fcol=shift_x+1)
        !call xgBlock_reshape(slice%eigenvalues, neigenpairs, 1)
        !call xgBlock_reshape(eig_part, slicedim, 1)

        ! Apply Rayleigh Ritz on slice (refinement)
        call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,slice%eigenvalues,&
            slice%blockdim_cprj,ierr,0,tim_RR,ABI_GPU_DISABLED,solve_ax_bx=.true.)

        ! ### current version does no update
        ! Extract eigenvalues to array
        !call xgBlock_reverseMap(eig_part, theta_, rows=1, cols=slicedim)
        !rayleigh_quotients(shift_x+1:shift_x+slicedim) = theta_(1,1:slicedim)

        !write(901,*) 'converged slice eigenvalues='
        !call xgBlock_print(eig_part, 901)
        !flush(901)

        if ( ierr /= 0 ) then
            ABI_BUG("RayleighRitz did not work")
        end if

        write(901,*) '%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
        flush(901)

    end do ! end restart 

    if (islice>1) then
        call xg_free(Xsum)
    end if
    
    if (niter == max_niter_restart) then
        call xg_free(dist1)
        call xg_free(dist2)
        call xg_free(dist12)
    end if

    ! ** do not update collective space
    ! Store to all-workspace
    !call xgBlock_copy(slice%X, X_part%self)
    !call xgBlock_copy(slice%AX, AX_part%self)
    !call xgBlock_copy(slice%cprjX, cprjX_part%self)

    ! ITEST
    write(901,*) 'On partial output - slice: Bortho test'
    write(901,*) 'slice%X before ortho=', xgBlock_getid(slice%X) 
    flush(901)
    call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%X,slice%cprjX,ierr,tim_ortho,&
        gpu_option,AX=slice%AX)
    write(901,*) 'slice%X after Bortho no1=', xgBlock_getid(slice%X) 
    flush(901)
    call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%X,slice%cprjX,ierr,tim_ortho,&
        gpu_option,AX=slice%AX)
    write(901,*) 'slice%X after Bortho no2=', xgBlock_getid(slice%X) 
    flush(901)
    ! ITEST

    ! ### current version does no update
    ! Update eigenvalues
    !call xgBlock_reverseMap(eig_part, theta_, rows=1, cols=slicedim)
    !rayleigh_quotients(shift_x+1:shift_x+slicedim) = theta_(1,1:slicedim)

    ! ITEST
    !write(901,*) 'eigenvalues after diago on slice=', rayleigh_quotients(:)
    write(901,*) 'slice%X after diago=', xgBlock_getid(slice%X) 
    write(901,*) 'slice%AllX after diago=', xgBlock_getid(slice%AllX) 
    write(901,*) 'slice%cprjX after diago=', xgBlock_getid(slice%cprjX)
    write(901,*) 'eigenvalues (slice merged)='
    call xgBlock_print(slice%eigenvalues,901)
    flush(901)
    ! ITEST

    ! free xg on slice
    call xg_free(X_part)
    call xg_free(AX_part)
    call xg_free(cprjX_part)

 end do

 ! Deallocate
 call xg_free(DivResults)
 ABI_FREE(ndeg_filter_slice)

 ! ITEST
 !write(901,*) 'On output: Bortho test'
 !write(901,*) 'slice%AllX before ortho=', xgBlock_getid(slice%AllX) 
 !flush(901)
 !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
 !   gpu_option,AX=slice%AllAX%self)
 !write(901,*) 'slice%AllX after Bortho no1=', xgBlock_getid(slice%AllX) 
 !flush(901)
 !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
 !    gpu_option,AX=slice%AllAX%self)
 !write(901,*) 'slice%AllX after Bortho no2=', xgBlock_getid(slice%AllX) 
 !flush(901)
 ! ITEST
        
 ! Apply Rayleigh Ritz on slice (refinement)
 !call xg_RayleighRitz_cprj(xg_nonlop,slice%AllX,slice%AllcprjX,slice%AllAX%self,slice%eigenvalues,&
 !       slice%all_blockdim_cprj,ierr,0,tim_RR,ABI_GPU_DISABLED,solve_ax_bx=.true.)

 !if ( ierr /= 0 ) then
 !   ABI_BUG("RayleighRitz did not work")
 !end if

 ! ITEST
 !write(901,*) 'eigenvalues(global RR)='
 !call xgBlock_print(slice%eigenvalues,901)
 !flush(901)
 ! ITEST

 ! Compute H-eSX
 if (slice%paw) then
    call timab(tim_AX_nl,1,tsec)
    call xg_nonlop_getHmeSX(xg_nonlop,slice%AllX,slice%AllcprjX,slice%AllAX%self,slice%eigenvalues,&
        slice%Allcprj_work%self,slice%cprj_work2%self,no_H=.True.)
    call timab(tim_AX_nl,2,tsec)
 end if

 ! Compute residual norm squared
 call timab(tim_residu, 1, tsec)
 if (.not.slice%paw) then
   call xgBlock_yxmax(slice%AllAX%self,slice%eigenvalues,slice%AllX)
 end if
 call xgBlock_colwiseNorm2(slice%AllAX%self, residu)
 call timab(tim_residu, 2, tsec)

 ! ITEST
 write(901,*) 'residu='
 call xgBlock_print(residu, 901)
 flush(901)
 ! ITEST

 ! Store result to output
 call timab(tim_copy, 1, tsec)
 call xgBlock_copy(slice%AllX,X0)
 call xgBlock_copy(slice%AllcprjX,cprjX0)
 call timab(tim_copy, 2, tsec)

 if (.not.slice%paw) then
   call timab(tim_enl,1,tsec)
   call xg_nonlop_colwiseXHX(xg_nonlop,slice%AllcprjX,slice%Allcprj_work%self,enl)
   call timab(tim_enl,2,tsec)
 end if

end subroutine slice_run_cprj
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_orthoXwrtBlocks
!! NAME
!! slice_orthoXwrtBlocks
!! 
!! FUNCTION
!! same as lobpcg_orthoXwrtBlocks but X0 and cprjX0 is given in input

subroutine slice_orthoXwrtBlocks(slice,X0,cprjX0,var,cprjvar,islice,cprj_work)

    type(slice_t) , intent(inout) :: slice
    type(xgBlock_t), intent(in) :: X0
    type(xgBlock_t), intent(inout) :: cprjX0
    type(xgBlock_t), intent(inout) :: var
    type(xgBlock_t), intent(inout) :: cprjvar
    type(xgBlock_t), intent(inout) :: cprj_work
    integer        , intent(in   ) :: islice
    integer :: previousBlock
    integer :: slicedim
    integer :: spacedim
    integer :: space_buf
    integer :: nspinor,cprjdim
    type(xg_t) :: buffer
    double precision :: tsec(2)
    type(xgBlock_t) :: cprjX0_spinor,cprj_work_spinor

    call timab(tim_ortho,1,tsec)

    if (islice<2) then
      ABI_ERROR("islice<2")
    end if

    if (cols(cprjvar)/=cols(cprj_work)) then
      ABI_ERROR("cprjvar and cprj_work should have same number of columns")
    end if
    slicedim = slice%slicedim
    spacedim = slice%spacedim
    previousBlock = (islice-1)*slice%slicedim

    cprjdim = slice%xg_nonlop%cprjdim
    nspinor = slice%xg_nonlop%nspinor

    space_buf = space(var)
    if (space(var)==SPACE_CR) then
      space_buf = SPACE_R
    end if
    call xg_init(buffer,space_buf,previousBlock,slicedim,slice%spacecom)

    ! buffer = X0^T*X
    call xgBlock_gemm('t','n',1.0d0,X0,var,0.d0,buffer%self,comm=slice%spacecom)

    ! Add the nonlocal part if paw
    if (slice%xg_nonlop%paw) then
      call xg_nonlop_getXSX(slice%xg_nonlop,cprjX0,cprjvar,cprj_work,buffer%self,slice%blockdim_cprj)
    end if

    ! sum all process contribution
    ! X = - X0*(BX0^T*X) + X
    call xgBlock_gemm('n','n',-1.0d0,X0,buffer%self,1.0d0,var)

    call xgBlock_zero(cprj_work)
    call xgBlock_reshape_spinor(cprj_work,cprj_work_spinor,nspinor,COLS2ROWS)
    call xgBlock_reshape_spinor(cprjX0,cprjX0_spinor,nspinor,COLS2ROWS)
    call xgBlock_gemm_mpi_cyclic_permutation(cprjX0_spinor,buffer%self,cprj_work_spinor,&
      & slice%xg_nonlop%me_band,slice%blockdim_cprj/nspinor,comm=slice%xg_nonlop%comm_band)
    call xgBlock_saxpy(cprjvar,-1.0d0,cprj_work)

    call xg_free(buffer)

    call timab(tim_ortho,2,tsec)

end subroutine slice_orthoXwrtBlocks
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_rayleighRitzQuotients
!! NAME
!! slice_rayleighRitzQuotients
!!
!! FUNCTION
!! Compute the Rayleigh-Ritz quotients.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!  maxeig= highest eigenvalue
!!  mineig= lowest eigenvalue
!!  DivResults= Rayleigh-Ritz quotients
!!
!! SOURCE

subroutine slice_rayleighRitzQuotients(slice,maxeig,mineig,DivResults)

!Arguments ------------------------------------
 real(dp), intent(inout) :: maxeig
 real(dp), intent(inout) :: mineig
 type(slice_t), intent(inout) :: slice
 type(xgBlock_t), intent(inout) :: DivResults

!Local variables-------------------------------
!scalars
 type(xg_t)::Results1
 type(xg_t)::Results2
 type(xg_t)::Results_work
!arrays
 integer :: maxeig_pos(2)
 integer :: mineig_pos(2)
 integer :: space_res

! *********************************************************************

 if (space(slice%AllX)==SPACE_C) then
   space_res = SPACE_C
 else if (space(slice%AllX)==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if
 call xg_init(Results1, space_res, slice%bandpp, 1)
 call xg_init(Results2, space_res, slice%bandpp, 1)

 call xgBlock_colwiseDotProduct(slice%AllX,slice%AllAX%self,Results1%self,comm_loc=xmpi_comm_null)

 call xgBlock_colwiseDotProduct(slice%AllX,slice%AllX,Results2%self,comm_loc=xmpi_comm_null)
 if (slice%xg_nonlop%paw) then
   call xg_init(Results_work, space_res, slice%bandpp, 1)
   call xg_nonlop_colwiseXAX(slice%xg_nonlop,slice%xg_nonlop%Sij%self,slice%AllcprjX,&
       slice%Allcprj_work%self,Results_work%self)
   call xgBlock_add(Results2%self,Results_work%self)
   call xg_free(Results_work)
 end if

 call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, maxeig, maxeig_pos, mineig, mineig_pos)

 call xg_free(Results1)
 call xg_free(Results2)

end subroutine slice_rayleighRitzQuotients
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_rayleighRitzQuotientsOnSlice
!! NAME
!! slice_rayleighRitzQuotientsOnSlice
!!
!! FUNCTION
!! Compute the Rayleigh-Ritz quotients.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!  maxeig= highest eigenvalue
!!  mineig= lowest eigenvalue
!!  DivResults= Rayleigh-Ritz quotients
!!
!! SOURCE

subroutine slice_rayleighRitzQuotientsOnSlice(slice,maxeig,mineig,DivResults)

!Arguments ------------------------------------
 real(dp), intent(inout) :: maxeig
 real(dp), intent(inout) :: mineig
 type(slice_t), intent(inout) :: slice
 type(xgBlock_t), intent(inout) :: DivResults

!Local variables-------------------------------
!scalars
 type(xg_t)::Results1
 type(xg_t)::Results2
 type(xg_t)::Results_work
!arrays
 integer :: maxeig_pos(2)
 integer :: mineig_pos(2)
 integer :: space_res

! *********************************************************************

 if (space(slice%X)==SPACE_C) then
   space_res = SPACE_C
 else if (space(slice%X)==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if
 call xg_init(Results1, space_res, slice%neigenpairs, 1)
 call xg_init(Results2, space_res, slice%neigenpairs, 1)

 call xgBlock_colwiseDotProduct(slice%X,slice%AX,Results1%self,comm_loc=xmpi_comm_null)

 call xgBlock_colwiseDotProduct(slice%X,slice%X,Results2%self,comm_loc=xmpi_comm_null)
 if (slice%xg_nonlop%paw) then
   call xg_init(Results_work, space_res, slice%neigenpairs, 1)
   call xg_nonlop_colwiseXAX(slice%xg_nonlop,slice%xg_nonlop%Sij%self,slice%cprjX,&
       slice%cprj_work,Results_work%self)
   call xgBlock_add(Results2%self,Results_work%self)
   call xg_free(Results_work)
 end if

 call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, maxeig, maxeig_pos, mineig, mineig_pos)

 call xg_free(Results1)
 call xg_free(Results2)

end subroutine slice_rayleighRitzQuotientsOnSlice
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_computeNextOrderChebfiPolynom
!! NAME
!! slice_computeNextOrderChebfiPolynom
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
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_computeNextOrderChebfiPolynom(slice,ideg,center,one_over_r,two_over_r)

!Arguments ------------------------------------
 real(dp)       , intent(in) :: center
 integer        , intent(in) :: ideg
 real(dp)       , intent(in) :: one_over_r
 real(dp)       , intent(in) :: two_over_r
 type(slice_t) , intent(inout) :: slice

 !Local variables-------------------------------
 real(dp) :: tsec(2)

 ! *********************************************************************

 call timab(tim_copy, 1, tsec)
 call xgBlock_copy(slice%AX,slice%X_next)
 call timab(tim_copy, 2, tsec)

 if (slice%paw) then
   call timab(tim_invovl, 1, tsec)
   call xg_nonlop_getSm1X(slice%xg_nonlop,slice%X_next,slice%cprjX,&
     & slice%cprj_work,slice%cprj_work2%self,slice%proj_work%self)
   call timab(tim_invovl, 2, tsec)
 else
   call timab(tim_copy, 1, tsec)
   call xgBlock_copy(slice%AX,slice%X_next)
   call timab(tim_copy, 2, tsec)
 end if

 call timab(tim_postinvovl, 1, tsec)
 call xgBlock_scale(slice%X, center, 1) !scale by center

 !(B-1 * A * Psi^i-1 - c * Psi^i-1)
 call xgBlock_saxpy(slice%X_next, dble(-1.0), slice%X)

 !Psi^i-1  = 1/c * Psi^i-1
 call xgBlock_scale(slice%X, 1/center, 1) !counter scale by 1/center

 if (ideg == 0) then
   call xgBlock_scale(slice%X_next, one_over_r, 1)
 else
   call xgBlock_scale(slice%X_next, two_over_r, 1)

   call xgBlock_saxpy(slice%X_next, dble(-1.0), slice%X_prev)
 end if

 call timab(tim_postinvovl, 2, tsec)

end subroutine slice_computeNextOrderChebfiPolynom
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_swapInnerBuffers
!! NAME
!! slice_swapInnerBuffers
!!
!! FUNCTION
!! Swap buffers inside a 'slice' datastructure.
!!
!! INPUTS
!!  ncols= number of requested eigenvectors/eigenvalues
!!  spacedim= space dimension for one vector
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_swapInnerBuffers(slice,spacedim,ncols)

  ! Arguments ------------------------------------
  integer        , intent(in   ) :: spacedim
  integer        , intent(in   ) :: ncols
  type(slice_t) , intent(inout) :: slice

  ! *********************************************************************

  call xgBlock_setBlock(slice%X_prev, slice%X_swap, spacedim, ncols) !X_swap = X_prev
  call xgBlock_setBlock(slice%X,      slice%X_prev, spacedim, ncols) !X_prev = X
  call xgBlock_setBlock(slice%X_next, slice%X,      spacedim, ncols) !X = X_next
  call xgBlock_setBlock(slice%X_swap, slice%X_next, spacedim, ncols) !X_next = X_swap

end subroutine slice_swapInnerBuffers
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_ampfactor
!! NAME
!! slice_ampfactor
!!
!! FUNCTION
!! Compute amplification factor
!!
!! INPUTS
!! eig (:,:)= eigenvalues
!! lambda_minus,lambda_plus=
!! ndeg_filter_bands(:)= degree of Spectrum Slicing filter for each band
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  residu<type(xgBlock_t)>= vector of residuals
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_ampfactor(slice,DivResults,lambda_minus,lambda_plus,ndeg_filter_bands)

  ! Arguments ------------------------------------
  integer,           intent(in   ) :: ndeg_filter_bands(:)
  type(xgBlock_t),   intent(in   ) :: DivResults
  real(dp),          intent(in   ) :: lambda_minus
  real(dp),          intent(in   ) :: lambda_plus
  type(slice_t),    intent(inout) :: slice

  ! Local variables-------------------------------
  ! scalars
  integer         :: iband
  real(dp)        :: ampfactor
  real(dp)        :: eig_per_band
  type(xgBlock_t) :: X_part
  type(xgBlock_t) :: AX_part
  real(dp),pointer :: eig(:,:)

  ! *********************************************************************

  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=cols(DivResults))

  do iband = 1, cols(DivResults)

    eig_per_band = eig(1,iband)

    !cheb_poly1(x, n, a, b)
    ampfactor = cheb_poly1(eig_per_band, ndeg_filter_bands(iband), lambda_minus, lambda_plus)

    if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much

    call xgBlock_setBlock(slice%X, X_part, slice%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(slice%AX, AX_part, slice%total_spacedim, 1, fcol=iband)

    call xgBlock_scale(X_part, 1/ampfactor, 1)
    call xgBlock_scale(AX_part, 1/ampfactor, 1)

  end do

end subroutine slice_ampfactor
!!***

subroutine slice_ampfactorProbe(slice,DivResults,lambda_minus,lambda_plus,ndeg_filter)

  ! Arguments ------------------------------------
  integer,           intent(in   ) :: ndeg_filter
  type(xgBlock_t),   intent(in   ) :: DivResults
  real(dp),          intent(in   ) :: lambda_minus
  real(dp),          intent(in   ) :: lambda_plus
  type(slice_t),    intent(inout) :: slice

  ! Local variables-------------------------------
  ! scalars
  integer         :: iband
  real(dp)        :: ampfactor
  real(dp)        :: eig_per_band
  type(xgBlock_t) :: X_col
  real(dp),pointer :: eig(:,:)

  ! *********************************************************************

  ! Apply amplification factor to f(Bm1A)X
  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=cols(DivResults))

  do iband = 1, cols(DivResults)

    eig_per_band = eig(1,iband)

    !cheb_poly1(x, n, a, b)
    ampfactor = cheb_poly1(eig_per_band, ndeg_filter, lambda_minus, lambda_plus)

    if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much

    call xgBlock_setBlock(slice%X_PROBE%self, X_col, slice%total_spacedim, 1, fcol=iband)
    call xgBlock_scale(X_col, 1/ampfactor, 1)

  end do

end subroutine slice_ampfactorProbe

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_probeProximity
!! NAME
!! slice_probeProximity
!! 
!! FUNCTION
!! Compute the principal angles to measure the subspace convergence
!! regarding the slice subspace associated to spectral interval [a,b).
!! Essentially computes the distance between two subspaces.
!! Note that A=Span(a1,..ak) and B=Span(b1,...bk) where families
!! are assumed to be linearly independent set of vectors.
!! They do not need to be orthogonal!!
!! 
!! IML debug version 08/07:
!! probe an interior slice, works for test case only

subroutine slice_probeProximity(slice,islice)

    implicit none
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: islice


end subroutine slice_probeProximity
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/slice_oracle1
!! NAME
!! slice_oracle1
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

!!****f* m_slice_cprj/cheb_poly1
!! NAME
!! cheb_poly1
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

end module m_slice_cprj
!!***

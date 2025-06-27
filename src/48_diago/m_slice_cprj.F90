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

!Public 'slice' datatype
!-------------------------------------------------

 type, public :: slice_t
   integer :: space
   integer :: space_cprj
   integer :: spacedim                      ! Space dimension for one vector
   integer :: cprjdim                       ! cprj dimension
   integer :: blockdim_cprj                 !
   integer :: total_spacedim                ! Maybe not needed
   integer :: neigenpairs                   ! Number of eigen values/vectors we want
   integer :: ndeg_filter                   ! Degree of the polynomial filter
   integer :: spacecom                      ! Communicator for MPI
   integer :: nslice                        ! Number of spectral slices
   integer :: paral_slice                   ! slice parallelization strategy
   integer :: spectral_cut                  ! how to cut the eigenvalue spectral interval
   real(dp) :: tolerance                    ! Tolerance on the residu to stop the minimization
   real(dp) :: ecut                         ! Ecut for Chebfi oracle
   real(dp) :: tolfilter                    ! Polynomial filter wanted amplification

   integer :: bandpp

   logical :: paw
   integer :: eigenProblem   !1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)
   integer :: me_g0
   integer :: me_g0_fft

   type(xg_nonlop_t) :: xg_nonlop

   !ARRAYS
   type(xgBlock_t) :: X

   type(xg_t) :: X_NP
   type(xgBlock_t) :: X_next
   type(xgBlock_t) :: X_prev

   type(xg_t) :: AX
   type(xgBlock_t) :: cprjX
   type(xg_t) :: cprj_work
   type(xg_t) :: cprj_work2
   type(xg_t) :: proj_work

   type(xgBlock_t) :: xXColsRows
   type(xgBlock_t) :: xAXColsRows

   type(xgTransposer_t) :: xgTransposerX
   type(xgTransposer_t) :: xgTransposerAX

   type(xgBlock_t) :: eigenvalues

   !SWAP POINTERS
   type(xgBlock_t) :: X_swap
   type(xgBlock_t) :: AX_swap

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
 slice%cprjdim       = cprjdim
 slice%blockdim_cprj = bandpp*xg_nonlop%nspinor
 if (tolerance > 0.0) then
   slice%tolerance = tolerance
 else
   slice%tolerance = 1.0e-20
 end if
 slice%ecut          = ecut
 slice%bandpp        = bandpp
 slice%ndeg_filter   = ndeg_filter
 slice%spacecom      = spacecom
 slice%eigenProblem  = eigenProblem
 slice%me_g0         = me_g0
 slice%me_g0_fft     = me_g0_fft
 slice%paw           = paw
 slice%xg_nonlop     = xg_nonlop
 slice%nslice        = nslice
 slice%tolfilter     = tolfilter
 slice%paral_slice   = paral_slice
 slice%spectral_cut  = spectral_cut

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
 integer  :: total_spacedim, ierr
 integer  :: nspinor

! *********************************************************************

 space       = slice%space
 space_cprj  = slice%space_cprj
 spacedim    = slice%spacedim
 neigenpairs = slice%neigenpairs
 nspinor = slice%xg_nonlop%nspinor

 call slice_free(slice)

 total_spacedim = spacedim
 call xmpi_sum(total_spacedim,slice%spacecom,ierr)
 slice%total_spacedim = total_spacedim
 call xg_init(slice%X_NP,space,total_spacedim,2*slice%bandpp,xmpi_comm_self,me_g0=slice%me_g0_fft) !transposed arrays
 call xg_setBlock(slice%X_NP,slice%X_next,total_spacedim,slice%bandpp)
 call xg_setBlock(slice%X_NP,slice%X_prev,total_spacedim,slice%bandpp,fcol=slice%bandpp+1)

 call xg_init(slice%AX,space,spacedim,neigenpairs,slice%spacecom,me_g0=slice%me_g0)
 call xg_init(slice%cprj_work ,space_cprj,slice%cprjdim,slice%blockdim_cprj,slice%spacecom)
 call xg_init(slice%cprj_work2,space_cprj,slice%cprjdim,slice%blockdim_cprj,slice%spacecom)

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

 call xg_free(slice%AX)
 call xg_free(slice%cprj_work)
 call xg_free(slice%cprj_work2)
 call xg_free(slice%proj_work)

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

!!****f* m_slice_cprj/slice_run
!! NAME
!! slice_run
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
 type(xg_t) :: DivResults
!arrays
 real(dp) :: tsec(2)
 !Pointers similar to old Chebfi
 integer,allocatable :: ndeg_filter_bands(:) !Oracle variable
 type(xg_nonlop_t) :: xg_nonlop

! *********************************************************************

 nslice = slice%nslice
 tolfilter = slice%tolfilter
 spectral_cut = slice%spectral_cut
 paral_slice = slice%paral_slice
 spacedim = slice%spacedim
 neigenpairs = slice%neigenpairs
 ndeg_filter = slice%ndeg_filter
 xg_nonlop = slice%xg_nonlop
 slice%eigenvalues = eigen

 write(std_out,*) 'inside slice_run'
 write(std_out,*) 'nslice=', nslice

 ABI_MALLOC(ndeg_filter_bands,(slice%bandpp))
 if (slice%space==SPACE_C) then
   space_res = SPACE_C
 else if (slice%space==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if
 call xg_init(DivResults, space_res, slice%bandpp, 1)

 tolerance = slice%tolerance
 lambda_plus = slice%ecut
 slice%X = X0
 slice%cprjX = cprjX0

! Transpose
 call timab(tim_transpose,1,tsec)
 call xgTransposer_constructor(slice%xgTransposerX,slice%X,slice%xXColsRows,nspinor,&
   STATE_LINALG,TRANS_ALL2ALL,xmpi_comm_self,slice%spacecom,0,0,slice%me_g0_fft)

 call xgTransposer_copyConstructor(slice%xgTransposerAX,slice%xgTransposerX,slice%AX%self,slice%xAXColsRows,STATE_LINALG)

 call xgTransposer_transpose(slice%xgTransposerX,STATE_COLSROWS)
 slice%xgTransposerAX%state = STATE_COLSROWS
 call timab(tim_transpose,2,tsec)

 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%xXColsRows,slice%cprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)
 call timab(tim_AX_v,1,tsec)
 call getAX(slice%xXColsRows,slice%xAXColsRows)
 call timab(tim_AX_v,2,tsec)
 call timab(tim_AX_k,1,tsec)
 call xgBlock_add_diag(slice%xXColsRows,kin,nspinor,slice%xAXColsRows)
 call timab(tim_AX_k,2,tsec)
 call timab(tim_AX_nl,1,tsec)
 call xg_nonlop_getHX(xg_nonlop,slice%xAXcolsRows,slice%cprjX,slice%cprj_work%self,slice%proj_work%self)
 call timab(tim_AX_nl,2,tsec)

 call timab(tim_barrier,1,tsec)
 call xmpi_barrier(slice%spacecom)
 call timab(tim_barrier,2,tsec)

!********************* Compute Rayleigh quotients for every band, and set lambda equal to the largest one *****
 call timab(tim_RR_q, 1, tsec)
 call slice_rayleighRitzQuotients(slice, maxeig, mineig, DivResults%self)

 call xmpi_max(maxeig,maxeig_global,slice%spacecom,ierr)
 call xmpi_min(mineig,mineig_global,slice%spacecom,ierr)
 call timab(tim_RR_q, 2, tsec)

 lambda_minus = maxeig_global
 ndeg_filter = slice%ndeg_filter

 ! TODO oracle must be replaced by the bandpass filter degree computation...
 call slice_set_ndeg_from_residu(slice,lambda_minus,lambda_plus,occ,DivResults%self,ndeg_filter_max,ndeg_filter)
 ndeg_filter_bands(:) = ndeg_filter

 center = (lambda_plus + lambda_minus)*0.5
 radius = (lambda_plus - lambda_minus)*0.5

 one_over_r = 1/radius
 two_over_r = 2/radius

 do ideg = 0, ndeg_filter - 1

   call timab(tim_cprj,1,tsec)
   call xg_nonlop_getcprj(xg_nonlop,slice%xAXcolsrows,slice%cprjX,slice%proj_work%self)
   call timab(tim_cprj,2,tsec)
   call slice_computeNextOrderChebfiPolynom(slice, ideg, center, one_over_r, two_over_r)

   call timab(tim_swap,1,tsec)
   call slice_swapInnerBuffers(slice, slice%total_spacedim, slice%bandpp)
   call timab(tim_swap,2,tsec)

   !A * Psi
   call timab(tim_AX_v,1,tsec)
   call getAX(slice%xXColsRows,slice%xAXColsRows)
   call timab(tim_AX_v,2,tsec)
   call timab(tim_AX_k,1,tsec)
   call xgBlock_add_diag(slice%xXColsRows,kin,nspinor,slice%xAXColsRows)
   call timab(tim_AX_k,2,tsec)
   call timab(tim_cprj,1,tsec)
   call xg_nonlop_getcprj(xg_nonlop,slice%xXColsRows,slice%cprjX,slice%proj_work%self)
   call timab(tim_cprj,2,tsec)
   call timab(tim_AX_nl,1,tsec)
   call xg_nonlop_getHX(xg_nonlop,slice%xAXcolsRows,slice%cprjX,slice%cprj_work%self,slice%proj_work%self)
   call timab(tim_AX_nl,2,tsec)

 end do

 call timab(tim_barrier,1,tsec)
 call xmpi_barrier(slice%spacecom)
 call timab(tim_barrier,2,tsec)

 call timab(tim_amp_f,1,tsec)
 call slice_ampfactor(slice, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)
 call timab(tim_amp_f,2,tsec)

 call xg_free(DivResults)
 ABI_FREE(ndeg_filter_bands)

 call timab(tim_transpose,1,tsec)
 call xmpi_barrier(slice%spacecom)

 call xgTransposer_transpose(slice%xgTransposerX,STATE_LINALG)
 call xgTransposer_transpose(slice%xgTransposerAX,STATE_LINALG)

 if (xmpi_comm_size(slice%spacecom) == 1) then !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
   call xgBlock_setBlock(slice%xXColsRows , slice%X      , spacedim, neigenpairs)
   call xgBlock_setBlock(slice%xAXColsRows, slice%AX%self, spacedim, neigenpairs)
 end if
 call timab(tim_transpose,2,tsec)

 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%cprj_work%self)
 call timab(tim_cprj,2,tsec)
 call xg_RayleighRitz_cprj(slice%xg_nonlop,slice%X,slice%cprjX,slice%AX%self,slice%eigenvalues,slice%blockdim_cprj,ierr,0,&
   tim_RR,ABI_GPU_DISABLED,solve_ax_bx=.true.)

 if (slice%paw) then
   call timab(tim_AX_nl,1,tsec)
   call xg_nonlop_getHmeSX(xg_nonlop,slice%X,slice%cprjX,slice%AX%self,slice%eigenvalues,slice%cprj_work%self,&
   & slice%cprj_work2%self,no_H=.True.)
   call timab(tim_AX_nl,2,tsec)
 end if

 call timab(tim_residu, 1, tsec)

 if (.not.slice%paw) then
   call xgBlock_yxmax(slice%AX%self,slice%eigenvalues,slice%X)
 end if

 call xgBlock_colwiseNorm2(slice%AX%self, residu)
 call timab(tim_residu, 2, tsec)

 call timab(tim_copy, 1, tsec)
 call xgBlock_copy(slice%X,X0)
 call timab(tim_copy, 2, tsec)

 call xgTransposer_free(slice%xgTransposerX)
 call xgTransposer_free(slice%xgTransposerAX)

 if (.not.slice%paw) then
   call timab(tim_enl,1,tsec)
   call xg_nonlop_colwiseXHX(slice%xg_nonlop,slice%cprjX,slice%cprj_work%self,enl)
   call timab(tim_enl,2,tsec)
 end if

end subroutine slice_run_cprj
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

 if (space(slice%xXcolsRows)==SPACE_C) then
   space_res = SPACE_C
 else if (space(slice%xXcolsRows)==SPACE_CR) then
   space_res = SPACE_R
 else
   ABI_ERROR('space(X) should be SPACE_C or SPACE_CR')
 end if
 call xg_init(Results1, space_res, slice%bandpp, 1)
 call xg_init(Results2, space_res, slice%bandpp, 1)

 call xgBlock_colwiseDotProduct(slice%xXColsRows,slice%xAXColsRows,Results1%self,comm_loc=xmpi_comm_null)

 call xgBlock_colwiseDotProduct(slice%xXColsRows,slice%xXColsRows,Results2%self,comm_loc=xmpi_comm_null)
 if (slice%xg_nonlop%paw) then
   call xg_init(Results_work, space_res, slice%bandpp, 1)
   call xg_nonlop_colwiseXAX(slice%xg_nonlop,slice%xg_nonlop%Sij%self,slice%cprjX,slice%cprj_work%self,Results_work%self)
   call xgBlock_add(Results2%self,Results_work%self)
   call xg_free(Results_work)
 end if

 call xgBlock_colwiseDivision(Results1%self, Results2%self, DivResults, maxeig, maxeig_pos, mineig, mineig_pos)

 call xg_free(Results1)
 call xg_free(Results2)

end subroutine slice_rayleighRitzQuotients
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
 call xgBlock_copy(slice%xAXColsRows,slice%X_next)
 call timab(tim_copy, 2, tsec)

 if (slice%paw) then
   call timab(tim_invovl, 1, tsec)
   call xg_nonlop_getSm1X(slice%xg_nonlop,slice%X_next,slice%cprjX,&
     & slice%cprj_work%self,slice%cprj_work2%self,slice%proj_work%self)
   call timab(tim_invovl, 2, tsec)
 else
   call timab(tim_copy, 1, tsec)
   call xgBlock_copy(slice%xAXColsRows,slice%X_next)
   call timab(tim_copy, 2, tsec)
 end if

 call timab(tim_postinvovl, 1, tsec)
 call xgBlock_scale(slice%xXColsRows, center, 1) !scale by center

 !(B-1 * A * Psi^i-1 - c * Psi^i-1)
 call xgBlock_saxpy(slice%X_next, dble(-1.0), slice%xXColsRows)

 !Psi^i-1  = 1/c * Psi^i-1
 call xgBlock_scale(slice%xXColsRows, 1/center, 1) !counter scale by 1/center

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
!!  neigenpairs= number of requested eigenvectors/eigenvalues
!!  spacedim= space dimension for one vector
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm
!!
!! SOURCE

subroutine slice_swapInnerBuffers(slice,spacedim,neigenpairs)

  ! Arguments ------------------------------------
  integer        , intent(in   ) :: spacedim
  integer        , intent(in   ) :: neigenpairs
  type(slice_t) , intent(inout) :: slice

  ! *********************************************************************

  call xgBlock_setBlock(slice%X_prev,     slice%X_swap,     spacedim, neigenpairs) !X_swap = X_prev
  call xgBlock_setBlock(slice%xXColsRows, slice%X_prev,     spacedim, neigenpairs) !X_prev = xXColsRows
  call xgBlock_setBlock(slice%X_next,     slice%xXColsRows, spacedim, neigenpairs) !xXColsRows = X_next
  call xgBlock_setBlock(slice%X_swap,     slice%X_next,     spacedim, neigenpairs) !X_next = X_swap

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

  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=slice%bandpp)

  do iband = 1, slice%bandpp

    eig_per_band = eig(1,iband)

    !cheb_poly1(x, n, a, b)
    ampfactor = cheb_poly1(eig_per_band, ndeg_filter_bands(iband), lambda_minus, lambda_plus)

    if(abs(ampfactor) < 1e-3) ampfactor = 1e-3 !just in case, avoid amplifying too much

    call xgBlock_setBlock(slice%xXColsRows, X_part, slice%total_spacedim, 1, fcol=iband)
    call xgBlock_setBlock(slice%xAXColsRows, AX_part, slice%total_spacedim, 1, fcol=iband)

    call xgBlock_scale(X_part, 1/ampfactor, 1)
    call xgBlock_scale(AX_part, 1/ampfactor, 1)

  end do

end subroutine slice_ampfactor
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

!!****f* m_slice/slice_set_ndeg_from_residu
!! NAME
!! slice_set_ndeg_from_residu
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

subroutine slice_set_ndeg_from_residu(slice,lambda_minus,lambda_plus,occ,DivResults,ndeg_filter_max,ndeg_filter)

 integer,intent(in) :: ndeg_filter_max
 integer,intent(out) :: ndeg_filter
 type(slice_t), intent(inout) :: slice
 type(xgBlock_t), intent(in)    :: occ
 type(xgBlock_t), intent(in)    :: DivResults
 real(dp), intent(in) :: lambda_minus, lambda_plus

 logical :: test1,test2,test3
 integer :: iband_tot,iband
 integer :: bandpp,ierr,ndeg_filter_tolwfr,ndeg_filter_decrease,ndeg_filter_all,shift
 integer,allocatable :: ndeg_filter_bands(:)
 type(xgBlock_t) :: occBlock,occ_reshaped
 type(xg_t) :: residu
 real(dp),pointer :: residu_(:,:),occ_(:,:)
 real(dp) :: eig_iband,res_iband,occ_iband
 real(dp),pointer :: eig(:,:)

 bandpp = slice%bandpp

 !Compute residu here for oracle, use X_next as a work space
 ! X_next = S|Psi>
 call xgBlock_copy(slice%xXColsRows,slice%X_next)
 if (slice%paw) then
   call xg_nonlop_getSX(slice%xg_nonlop,slice%X_next,slice%cprjX,slice%cprj_work%self,slice%proj_work%self)
 end if
 ! X_next = - eig * S|Psi>
 call xgBlock_ymax(slice%X_next,DivResults,0,1)
 ! X_next = H|Psi> - eig * S|Psi>
 call xgBlock_add(slice%X_next,slice%xAXColsRows)
 ! resid = |X_next|^2
 call xg_init(residu,SPACE_R,bandpp,1)
 call xgBlock_colwiseNorm2(slice%X_next, residu%self,comm_loc=xmpi_comm_null)

 occ_reshaped = occ
 shift=xmpi_comm_rank(slice%spacecom)*bandpp
 call xgBlock_reshape(occ_reshaped,1,slice%neigenpairs)
 call xgBlock_setBlock(occ_reshaped,occBlock,1,bandpp,fcol=1+shift)
 call xgBlock_reshape(occBlock,bandpp,1)
 !if (slice%nbdbuf==-101) then
 !  call xgBlock_apply_diag(residu%self,occBlock,1)
 !end if

 ABI_MALLOC(ndeg_filter_bands,(bandpp))

 ! DivResults could be complex (with null imaginary part), so bandpp has to be in cols, not rows
 call xgBlock_reverseMap(DivResults,eig,rows=1,cols=bandpp)
 call xgBlock_reverseMap(residu%self,residu_,rows=1,cols=bandpp)
 call xgBlock_reverseMap(occBlock,occ_,rows=1,cols=bandpp)

 do iband=1, bandpp
   eig_iband = eig(1,iband)
   res_iband = residu_(1,iband)
   occ_iband = occ_(1,iband)
   iband_tot = iband + shift
   !ndeg_filter necessary to converge to tolerance
   ndeg_filter_tolwfr = cheb_oracle1(eig_iband, lambda_minus, lambda_plus, slice%tolerance / res_iband, 1000)
   if (slice%tolfilter<0) then
     ndeg_filter_bands(iband) = MIN(ndeg_filter_max, ndeg_filter_tolwfr, slice%ndeg_filter)
   else if (slice%tolfilter>0) then
     !ndeg_filter necessary to decrease residual by a constant factor
     ndeg_filter_decrease = cheb_oracle1(eig_iband, lambda_minus, lambda_plus, slice%tolfilter, 15)
     ndeg_filter_bands(iband) = MIN(ndeg_filter_max, ndeg_filter_tolwfr, ndeg_filter_decrease)
   else
     ABI_ERROR('Wrong value for slice%tolfilter')
   end if
 end do
 ndeg_filter = MAXVAL(ndeg_filter_bands)
 call xmpi_max(ndeg_filter,ndeg_filter_all,slice%spacecom,ierr)
 ndeg_filter=ndeg_filter_all

 call xg_free(residu)
 ABI_FREE(ndeg_filter_bands)

end subroutine slice_set_ndeg_from_residu
!!***

end module m_slice_cprj
!!***

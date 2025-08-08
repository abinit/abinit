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
 
 integer, parameter :: tim_slice1_fi    = 2190
 integer, parameter :: tim_slice1_rr    = 2191
 integer, parameter :: tim_slice1_pr    = 2192
 integer, parameter :: tim_slice2_fi    = 2193
 integer, parameter :: tim_slice2_rr    = 2194
 integer, parameter :: tim_slice2_pr    = 2195

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
   integer :: nbdbuf                        ! Number of bands in the buffer
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
                      ndeg_filter,nbdbuf,space,space_cprj,eigenProblem,spacecom,me_g0,paw,&
                      nslice,tolfilter,paral_slice,spectral_cut,&
                      xg_nonlop,me_g0_fft)

!Arguments ------------------------------------
 integer          , intent(in   ) :: bandpp
 integer          , intent(in   ) :: eigenProblem
 integer          , intent(in   ) :: me_g0
 integer          , intent(in   ) :: me_g0_fft
 integer          , intent(in   ) :: neigenpairs
 integer          , intent(in   ) :: ndeg_filter
 integer          , intent(in   ) :: nbdbuf
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
 slice%nbdbuf        = nbdbuf
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

 ! This is temporary slice workspace that is the same for all slice initially
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
 integer :: count_rr
 integer :: count_merge
 integer :: count_slice
 integer :: icount
 integer :: blockdim_cprj
 integer :: fcol_in, lcol_in
 integer :: fcol_global
 integer :: ndeg, ndeg_max
 integer :: tim_slice_fi
 integer :: tim_slice_rr
 integer :: tim_slice_pr
 logical :: is_close_to_V
 real(dp) :: tol_step
 real(dp) :: tol_probe
 real(dp) :: tolerance
 real(dp) :: tolfilter
 real(dp) :: ramp
 real(dp) :: maxeig, maxeig_global
 real(dp) :: mineig, mineig_global
 real(dp) :: lambda_minus, alpha_minus
 real(dp) :: lambda_plus, alpha_plus
 real(dp) :: ls_in, us_in
 real(dp) :: f_u, f_l, f_uw, f_lw
 real(dp) :: overlap_width
 real(dp) :: one_over_r
 real(dp) :: two_over_r
 real(dp) :: center
 real(dp) :: radius
 real(dp) :: ls, us, cdeg, mu, damp
 real(dp) :: max_dist2
 real(dp) :: slice1_mineig, slice1_maxeig
 real(dp) :: amp_ideg
 real(dp) :: ein_ideg, eout_ideg
 real(dp) :: tol12 = 1.0e-12
 type(xg_t) :: Xsum
 type(xg_t) :: DivResults
 type(xg_t) :: dist1, dist2, dist3, dist12
 type(xg_t) :: X0_out, eigen_out
 type(xgBlock_t) :: X0_out_part, eigen_out_part
 type(xgBlock_t) :: X_in, eigen_in
 type(xgBlock_t) :: DivResults_part
 type(xgBlock_t) :: X_prev
 type(xgBlock_t) :: cprjX_prev
 type(xgBlock_t) :: cprj_work_prev
 type(xgBlock_t) :: eig_part
 type(xgBlock_t) :: X_col, AX_col
 type(xg_t) :: X_kept, AX_kept
 type(xg_t) :: cprj_work_slice, cprj_work2_slice
 type(xgBlock_t) :: X_kept_col, AX_kept_col
 type(xgBlock_t) :: eigenvalues_slice, residu_slice
 type(xg_t) :: X_part
 type(xg_t) :: AX_part
 type(xg_t) :: cprjX_part
 integer,parameter :: gpu_option=ABI_GPU_DISABLED
!arrays
 real(dp) :: tsec(2)
 integer, allocatable :: permute_cols(:)
 integer, allocatable :: sorted_idx(:) ! same as permute_cols but used elsewhere
 real(dp), allocatable :: rayleigh_quotients(:)
 real(dp), pointer :: probe(:) => null()
 real(dp), pointer :: X0_norm2(:) => null()
 real(dp), pointer :: lambda_apost(:) => null()
 real(dp), pointer :: lambda_apost_slice(:) => null()
 real(dp), pointer :: dist2_array(:) => null()
 real(dp), pointer :: dist3_array(:) => null()
 real(dp), pointer :: theta_(:,:) => null()
 real(dp), pointer :: resid(:) => null()
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

 ! Memory used to store results of slice merging. 
 call xg_init(X0_out, slice%space, spacedim, neigenpairs, slice%spacecom, me_g0=slice%me_g0)
 call xg_init(eigen_out, SPACE_R, 1, neigenpairs)

 ! preallocate, will be updated during ideg iterations
 call xg_init(dist1,SPACE_R,neigenpairs,1)
 call xg_init(dist2,SPACE_R,neigenpairs,1)
 call xg_init(dist3,SPACE_R,neigenpairs,1)
 
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
 !theta_(1,1:neigenpairs) = rayleigh_quotients(1:neigenpairs)

 ! ITEST
 write(901,*) 'rayleigh quotients='
 call xgBlock_print(DivResults%self, 901)
 flush(901)
 ! ITEST

 ! Compute |X|^2 colwise L2-norm (before any filter)
 call xgBlock_colwiseNorm2(slice%AllX,dist1%self,comm_loc=xmpi_comm_null)
 call xgBlock_reverseMap_1d(dist1%self,X0_norm2)

 ! ITEST
 write(901,*) 'norm2(squared) ||X||='
 call xgBlock_print(dist1%self,901)
 flush(901)
 ! ITEST
 
 ! Permute columns 1,..,neigenpairs according to order
 !call xgBlock_permuteCols(slice%AllX, slice%spacedim, neigenpairs, permute_cols)
 !call xgBlock_permuteCols(slice%AllAX%self, slice%spacedim, neigenpairs, permute_cols)
 ! ITEST
 !write(901,*) 'slice%AllX after perm=', xgBlock_getid(slice%AllX) 
 !write(901,*) 'slice%AllAX after perm=', xgBlock_getid(slice%AllAX%self) 
 !flush(901)
 ! ITEST 

 ! Compare this probe to rayleigh quotients. Actually count the number of
 ! rayleigh quotients in the slice interval. We expect to see that it is 
 ! not the same as the final converged number


 !! ------------------------------------------------------------
 !! 
 !! -                      Main slice loop                     -
 !! 
 !! ------------------------------------------------------------

 count_merge = 0
 fcol_global = 1

 do islice=1, nslice

    ! ITEST
    write(901,*)
    write(901,*) '====================Slice=================', islice
    flush(901)
    ! ITEST

    !! ------------------------------------------------------------
    !! 
    !! -                Global to Slice operation                 -
    !! 
    !! ------------------------------------------------------------

    if (islice>1) then
        
        ! Reset pointers to dimensions of nband
        call xg_setBlock(slice%X_SLICE,slice%X,slice%total_spacedim,neigenpairs)
        call xg_setBlock(slice%X_SLICE,slice%AX,slice%total_spacedim,neigenpairs,fcol=neigenpairs+1)

    end if

    ! Initial slice workspaces are independent entire arrays
    call xgBlock_copy(slice%AllX, slice%X)
    call xgBlock_copy(slice%AllAX%self, slice%AX)
 
    ! the two following ones will be recomputed so whatever
    slice%cprjX = slice%AllcprjX
    slice%cprj_work = slice%Allcprj_work%self

    ! ITEST
    write(901,*) 
    write(901,*) 'X id at slice init', xgBlock_getid(slice%X)
    write(901,*) 'AX id at slice init', xgBlock_getid(slice%AX)
    write(901,*) 
    flush(901)
    ! ITEST

    !! ------------------------------------------------------------
    !! 
    !! -                Scalar polynomial tuning                  -
    !! 
    !! ------------------------------------------------------------

    if (islice==1) then

        nband_slice_buf = neigenpairs/nslice
        lambda_minus = rayleigh_quotients(nband_slice_buf)
        !lambda_minus = (maxval(rayleigh_quotients) - minval(rayleigh_quotients)) / 2.0

        !lambda_minus = 0.51404d0
        lambda_plus = slice%ecut
        center = (lambda_plus + lambda_minus)*0.5
        radius = (lambda_plus - lambda_minus)*0.5

        ndeg_filter = slice%ndeg_filter

        ! ITEST
        write(901,*) 'unwanted part of the spectrum=', lambda_minus, lambda_plus
        flush(901)
        ! ITEST

    else

        ! Assumes sequential, =smallest converged value from first slice
        call xgBlock_minmax(eigen_out%self, slice1_mineig, slice1_maxeig)

        ! ITEST
        write(901,*) 'detected previous slice, mineig=', slice1_mineig
        write(901,*) '                         maxeig=', slice1_maxeig
        flush(901)
        ! ITEST

        !alpha_minus = 0.51404 ! assumes sequential, =largest previously converged value
        alpha_minus = slice1_maxeig
        alpha_plus = rayleigh_quotients(neigenpairs) ! can be in parallel, only depends on Rayleigh value
        ! FIXME only works for slice=1
        ABI_WARNING('present code only works for nslice=2')

        overlap_width = (alpha_plus - alpha_minus)/8.0
        lambda_minus = alpha_minus-overlap_width 
        lambda_plus = alpha_plus+overlap_width

        ! ITEST
        write(901,*) 'wanted part of the spectrum=', alpha_minus, alpha_plus
        write(901,*) '               with overlap=', lambda_minus, lambda_plus
        flush(901)
        ! ITEST
 
        ! TODO add overlap width to filter in [a-w,b+w) for wanted [a,b) for convergence reasons
        center = (slice%ecut + slice1_mineig)*0.5
        radius = (slice%ecut - slice1_mineig)*0.5 
        ls = (lambda_minus - center) / radius
        us = (lambda_plus - center) / radius

        ! Optimize polynomial degree with given tolerance
        ramp = slice%tolfilter
        ndeg_max = 200

        ls_in = (alpha_minus - center) / radius ! scaled point inside slice
        us_in = (alpha_plus - center) / radius ! scaled point inside slice
        ndeg = 8
        f_lw = 1.d0; f_uw = 1.d0; f_l = 1.d0; f_u = 1.d0

        ! Amplification is f(l-w)/f(l) and f(u+w)/f(u) (between 0 and 1)
        do while ( (f_lw/f_l > ramp .or. f_uw/f_u > ramp) .and. ndeg < ndeg_max )
            ndeg = ndeg + 1
            f_l  = bandpassIndicator_sca(ls_in,ls,us,ndeg)
            f_lw = bandpassIndicator_sca(ls   ,ls,us,ndeg)
            f_u  = bandpassIndicator_sca(us_in,ls,us,ndeg)
            f_uw = bandpassIndicator_sca(us   ,ls,us,ndeg)
        end do

        ndeg_filter = ndeg
       
        ! ITEST
        write(901,*) 'left/right amplif factor f(out)/f(in)=', f_lw/f_l, f_uw/f_u
        write(901,*) 'minimal polynomial degree=', ndeg_filter
        flush(901)
        ! ITEST

    end if

    one_over_r = 1.0/radius
    two_over_r = 2.0/radius

    !! ------------------------------------------------------------
    !! 
    !! -                Polynomial degree loop                    -
    !! 
    !! ------------------------------------------------------------
 
    if (islice==1) then
        tim_slice_fi = tim_slice1_fi
    else
        tim_slice_fi = tim_slice2_fi
    end if

    call timab(tim_slice_fi,1,tsec)
    
    if (islice>1) then
    
        ! Initialize Chebyshev expansion of indicator function of order ndeg_filter
        call xg_init(Xsum,slice%space,slice%total_spacedim,neigenpairs,slice%spacecom,gpu_option=gpu_option)
        cdeg = Pi/(ndeg_filter+2)
        mu = 1.d0/Pi*(ACOS(ls)-ACOS(us))
        damp = 1.d0 ! Jackson damping
        call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)

    end if
        
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

        if (islice>1) then

            ! Accumulate X with weight in Xsum for bandpass filters
            mu = 2/Pi * (SIN((ideg+1)*ACOS(ls)) - SIN((ideg+1)*ACOS(us)))/(ideg+1)
            damp = ((1 - (ideg+1)/(ndeg_filter+2))*SIN(cdeg)*COS((ideg+1)*cdeg) + &
                     1/(ndeg_filter+2)*COS(cdeg)*SIN((ideg+1)*cdeg))/SIN(cdeg)
            call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)

        end if

        if (islice>1 .and. ideg==ndeg_filter - 1) then 
        
            ! store final expansion Xsum to X
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

    end do ! End polynomial degree loop

    if (islice==1) then

        ! Amplify to finalize filter application
        call timab(tim_amp_f,1,tsec)
        ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
        ndeg_filter_bands(:) = ndeg_filter
        !call slice_ampfactorMax(slice, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)
        call slice_ampfactor(slice, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)
        ABI_FREE(ndeg_filter_bands)
        call timab(tim_amp_f,2,tsec)
        
        call timab(tim_cprj,1,tsec)
        call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
        call timab(tim_cprj,2,tsec) 

    end if
    
    call timab(tim_slice_fi,2,tsec)

    !! ------------------------------------------------------------
    !! 
    !! -           Construct approximated subspace                -
    !! 
    !! ------------------------------------------------------------

    if (islice==1) then
        tim_slice_pr = tim_slice1_pr
    else
        tim_slice_pr = tim_slice2_pr
    end if

    call timab(tim_slice_pr,1,tsec)

    ! Compute probe
    ! dist2 = |X_PROBE|^2 colwise L2-norm
    call xgBlock_copy(slice%X, slice%X_PROBE%self)
    call xgBlock_colwiseNorm2(slice%X_PROBE%self, dist2%self, max_dist2, comm_loc=xmpi_comm_null)
    
    if (slice%spectral_cut == 1) then

        ! TODO rename variables
        call xgBlock_reverseMap_1d(dist2%self,probe)
        !call xgBlock_reverseMap_1d(dist3%self, dist3_array)

        ! Normalize by norm of initial X
        probe(:) = probe(:) / sqrt(X0_norm2(:))

        ! Scale probe to get a pivot between 0 and 1
        ! useful for absolute probe
        probe(:) = probe(:) / maxval(probe)

        ! ITEST
        write(901,*) 
        write(901,*) 'indicator for components in wanted eigenspace=', probe(:)
        !write(901,*) 'indicator for components in orthogonal complement=', dist3_array(:)
        !write(901,*) 'eout_ideg=', eout_ideg
        !write(901,*) 'ein_ideg=', ein_ideg
        !write(901,*) 'estimator=', eout_ideg**2 + ((1+ein_ideg)*0.1d0)**2
        flush(901)
        ! ITEST

    end if

    ! Criterion is deactivated by default
    count_mask = neigenpairs
   
    ! Old criterion kept for reference
    !count_mask = count(probe > dist3_array .or. ( probe > 1 ))

    !tol_probe = 0.3
    ! Initialize with quantity that is very small
    tol_probe = sum(probe)/neigenpairs ! initialize with average value (a little less)
    tol_step = tol_probe * 0.05 ! step is 10%
    write(901,*) 'initial tolerance (average)=', tol_probe
    write(901,*) 'starting refinement with step=', tol_step
    flush(901)

    ! Adaptive refinement
    if (slice%spectral_cut == 1) then
        
        !if (islice==1) then
        !    
        !    !count_mask = count( probe > eout_ideg**2 + (ein_ideg*0.1d0)**2 )
        !
        !else
        !
        !    ! Criterion to take into account error of f
        !    ! FIXME norm of X? Tolerance? slice=1?
        !    !count_mask = count( probe > eout_ideg**2 + ((1+ein_ideg)*0.1d0)**2 )
        !    count_mask = count( probe > tol_probe) 
        !
        !end if
        count_mask = count( probe > tol_probe ) 

    end if
    
    write(901,*) 
    write(901,*) 'initial count_mask=', count_mask
    flush(901)

    if (islice==1) then

        ! Discard vectors: Increase tolerance for probe if too many vectors in slice
!        icount = 1
!        do while(count_mask > 1.4*neigenpairs/nslice) ! maximum columns in block
!            tol_probe = tol_probe + tol_step
!            count_mask = count( probe > tol_probe)
!            write(901,*) '#icount, tol_probe=, count_mask=', icount, tol_probe, count_mask
!            icount = icount + 1
!        end do
        !tol_probe = sum(probe)/neigenpairs*0.3
        tol_probe = -1
        !count_mask = count( probe > tol_probe )
        count_mask = neigenpairs - 3*slice%nbdbuf
        ! TODO perform Alternating method where we adjust sizes upper and lower by alternating
        ! between the two
        write(901,*) 'fixed tolerance using nbdbuf=', 3*slice%nbdbuf
        flush(901)
 
    else
        
        ! Add vectors: Decrease tolerance for probe if not enough vectors after merge
        icount = 1
        do while (count_mask < 1.3*neigenpairs/nslice .or. count_mask + count_merge < neigenpairs)
            tol_probe = tol_probe - tol_step
            count_mask = count( probe > tol_probe)
            write(901,*) '#icount, tol_probe=, count_mask=', icount, tol_probe, count_mask
            icount = icount + 1
        end do
        write(901,*) 'refined tolerance, #iterations=', tol_probe, icount
        flush(901)

    end if

    ! ITEST
    write(901,*) 'Keep count_mask= out of neigenpairs=', count_mask, neigenpairs
    write(901,*) 
    flush(901)
    ! ITEST

    if (slice%spectral_cut == 1) then

        ! Allocate slice subspace memory 
        call xg_init(X_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
        call xg_init(AX_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
   
        ! Copy data
        icount = 1
        do iband=1,neigenpairs
            is_close_to_V = .true.
            if (islice==1) then
            !    !is_close_to_V = probe(iband) > eout_ideg**2 + (ein_ideg*0.1d0)**2
                 is_close_to_V = iband < neigenpairs - 3*slice%nbdbuf + 1
            else
            !    !is_close_to_V = probe(iband) > eout_ideg**2 + ((1+ein_ideg)*0.1d0)**2
                 is_close_to_V = probe(iband) > tol_probe
            end if
            if (is_close_to_V) then
                call xgBlock_setBlock(X_kept%self, X_kept_col, slice%total_spacedim, 1, fcol=icount)
                call xgBlock_setBlock(AX_kept%self, AX_kept_col, slice%total_spacedim, 1, fcol=icount)
                call xgBlock_setBlock(slice%X, X_col, slice%total_spacedim, 1, fcol=iband)
                call xgBlock_setBlock(slice%AX, AX_col, slice%total_spacedim, 1, fcol=iband)
    
                call timab(tim_copy, 1, tsec)
                call xgBlock_copy(X_col, X_kept_col)
                call xgBlock_copy(AX_col, AX_kept_col)
                call timab(tim_copy, 2, tsec)

                icount = icount + 1
            end if
        end do

        ! reset pointers to temporary (kept)
        slice%X = X_kept%self
        slice%AX = AX_kept%self

    end if
 
    ! content is not important, but dimensions 
    call xgBlock_setBlock(slice%AllcprjX, slice%cprjX, slice%cprjdim, count_mask*nspinor)

    ! Recompute cprj to be sure
    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
    call timab(tim_cprj,2,tsec)
    
    call timab(tim_slice_pr,2,tsec)

    !! ------------------------------------------------------------
    !! 
    !! -               Apply Rayleigh-Ritz step                   -
    !! 
    !! ------------------------------------------------------------

    if (islice==1) then
        tim_slice_rr = tim_slice1_rr
    else
        tim_slice_rr = tim_slice2_rr
    end if

    ! Number of vectors on which we apply rr
    count_rr = count_mask

    ! Restrict eigenvalue array for size consistancy (essentially keep nonzero entries)
    call xgBlock_reshape(slice%eigenvalues, 1, neigenpairs) 
    call xgBlock_setBlock(slice%eigenvalues, eigenvalues_slice, rows=1, cols=count_rr)
    call xgBlock_reshape(eigenvalues_slice, count_rr, 1)
    call xgBlock_reshape(slice%eigenvalues, neigenpairs, 1)

    ! Orthonormalize
    !call xg_Borthonormalize_cprj(xg_nonlop,slice%blockdim_cprj,slice%X,slice%cprjX,&
    !    ierr,tim_ortho,gpu_option,AX=slice%AX)

    ! Apply Rayleigh Ritz on slice (refinement)
    ! prtvol = 15015015 to print condition number of overlap matrix
    call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,eigenvalues_slice,&
        slice%blockdim_cprj,ierr,15015015,tim_slice_rr,ABI_GPU_DISABLED,solve_ax_bx=.true.)
    
    if ( ierr /= 0 ) then
        ABI_BUG("RayleighRitz did not work")
    end if

    ! restart!
!    call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,eigenvalues_slice,&
!        slice%blockdim_cprj,ierr,15015015,tim_slice_rr,ABI_GPU_DISABLED,solve_ax_bx=.true.)

!    if ( ierr /= 0 ) then
!        ABI_BUG("RayleighRitz did not work")
!    end if

    ! ITEST
    write(901,*) 'converged eigenval='
    call xgBlock_print(eigenvalues_slice, 901)
    flush(901)
    ! ITEST

    ! Restrict dimension of residual array
    call xgBlock_reshape(residu, 1, neigenpairs) 
    call xgBlock_setBlock(residu, residu_slice, rows=1, cols=count_rr)
    call xgBlock_reshape(residu_slice, count_rr, 1)
    call xgBlock_reshape(residu, neigenpairs, 1)
    
    !! ------------------------------------------------------------
    !! 
    !! -            Compute residual HX-eSX on slice              -
    !! 
    !! ------------------------------------------------------------

    blockdim_cprj = count_rr*xg_nonlop%nspinor
    call xg_init(cprj_work_slice,slice%space_cprj,slice%cprjdim,blockdim_cprj,slice%spacecom)
    call xg_init(cprj_work2_slice,slice%space_cprj,slice%cprjdim,blockdim_cprj,slice%spacecom)
 
    ! Compute H-eSX
    if (slice%paw) then
        call timab(tim_AX_nl,1,tsec)
        call xg_nonlop_getHmeSX(xg_nonlop,slice%X,slice%cprjX,slice%AX,eigenvalues_slice,&
            cprj_work_slice%self,cprj_work2_slice%self,no_H=.True.)
        call timab(tim_AX_nl,2,tsec)
    end if

    call xg_free(cprj_work_slice)
    call xg_free(cprj_work2_slice)

    ! Compute residual norm squared
    call timab(tim_residu, 1, tsec)
    if (.not.slice%paw) then
        call xgBlock_yxmax(slice%AX,eigenvalues_slice,slice%X)
    end if
    call xgBlock_colwiseNorm2(slice%AX, residu_slice)
    call timab(tim_residu, 2, tsec)

    call xgBlock_reverseMap_1d(residu_slice, resid)
 
    ! ITEST
    write(901,*) 'Slice 2: colwiseNorm2 residu='
    call xgBlock_print(residu_slice, 901)
    write(901,*) 'Frobenius norm=', sqrt(sum(resid))
    flush(901)
    ! ITEST
    
    !! ------------------------------------------------------------
    !! 
    !! -                Slice to Global operation                 -
    !! 
    !! ------------------------------------------------------------
    
    ! Merge to X0_out, eigen_out using Rayleigh value as criterion
    call xgBlock_reverseMap_1d(eigenvalues_slice, lambda_apost_slice)
    fcol_in = 1
    lcol_in = count_rr
    if (islice==1) then

        lcol_in = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < lambda_minus))

        !! If too many columns and not enough for second slice, then reduce
        !! can happen if guess is too bad
        !if (lcol_in - fcol_in + 1 > neigenpairs / nslice + 20) then
        !    lcol_in = neigenpairs / nslice
        !end if

        if (lcol_in < count_rr) then
        
            ! if multiple eigenvalue is at endpoint, backwards !remove! its multiplicities
            write(901,*) 'last      =', lambda_apost_slice(lcol_in), lcol_in
            write(901,*) 'after last=', lambda_apost_slice(lcol_in+1)
            flush(901)
            do while (lambda_apost_slice(lcol_in+1) - lambda_apost_slice(lcol_in) < 1.0e-3)
                lcol_in = lcol_in - 1
                write(901,*) 'is multiple eigenvalue, force to include it', lcol_in
                flush(901)
            end do
            lcol_in = lcol_in - 2 
            ! plus some extra space
            write(901,*) 'safety:', lcol_in
            flush(901)

        end if

    else

        fcol_in = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < alpha_minus)) + 1

        if (islice<nslice) then

            lcol_in = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < alpha_plus))

        else

            write(901,*) 'pass', lcol_in
            flush(901)
            if (count_merge + lcol_in-fcol_in+1 > neigenpairs) then
                lcol_in = fcol_in + neigenpairs - count_merge - 1
                write(901,*) 'pass', lcol_in, fcol_in, neigenpairs, count_merge
                flush(901)
            end if

        end if

    end if

    count_slice = lcol_in - fcol_in + 1

    ! ITEST
    if (islice==1) then
        write(901,*) 'inside indices [a,b)=', fcol_in, lcol_in, lambda_minus
    else
        write(901,*) 'inside indices [a,b)=', fcol_in, lcol_in, alpha_minus, alpha_plus
    end if
    write(901,*) 'Frobenius norm (inside only)=', sqrt(sum(resid(fcol_in:lcol_in)))
    flush(901)
    ! ITEST

    ! ITEST
    write(901,*) 
    write(901,*) '======================================'
    write(901,*) 'count_merged=', count_merge
    write(901,*) 'count_slice =', count_slice
    write(901,*) 'count_rr    =', count_rr
    write(901,*) 'fcol_in     =', fcol_in
    write(901,*) 'lcol_in     =', lcol_in
    write(901,*) 'lambda_fcol =', lambda_apost_slice(fcol_in)
    write(901,*) 'lambda_lcol =', lambda_apost_slice(lcol_in)
    write(901,*) '======================================'
    write(901,*)
    flush(901)
    ! ITEST

    count_merge = count_merge + count_slice
    if (count_merge > neigenpairs) then

        ABI_WARNING("Attempting to merge more bands than possible")

    end if

    if (count_merge < neigenpairs .and. islice==nslice) then
    
        write(901,*) 'missing eigenvalues!'
        flush(901)

    end if
    
    call xgBlock_setBlock(slice%X, X_in, spacedim, count_slice, fcol=fcol_in)
    call xgBlock_reshape(eigenvalues_slice, 1, count_rr)
    call xgBlock_setBlock(eigenvalues_slice, eigen_in, 1, count_slice, fcol=fcol_in)
    call xgBlock_reshape(eigenvalues_slice, count_rr, 1)

    call xgBlock_setBlock(X0_out%self, X0_out_part, spacedim, count_slice, fcol=fcol_global)
    call xgBlock_setBlock(eigen_out%self, eigen_out_part, 1, count_slice, fcol=fcol_global)
 
    call timab(tim_copy, 1, tsec)
    call xgBlock_copy(X_in, X0_out_part)
    call xgBlock_copy(eigen_in, eigen_out_part)
    call timab(tim_copy, 2, tsec)

    fcol_global = count_merge + 1

    !! Free slice memory whose size depends on count_mask, different for every slice
    if (islice>1) then
        
        call xg_free(Xsum)

    end if

    if (slice%spectral_cut==1) then

        call xg_free(X_kept)
        call xg_free(AX_kept)

    end if

    !! ------------------------------------------------------------
    !! 
    !! -     B-ortho X_part With Respect To previous blocks       -
    !! -                 kept for reference                       -
    !! 
    !! ------------------------------------------------------------
    
    ! call xgBlock_setBlock(slice%AllX             , X_prev         , spacedim     , shift_x   , fcol=1)
    ! call xgBlock_setBlock(slice%AllcprjX         , cprjX_prev     , slice%cprjdim, shift_cprj, fcol=1)
    ! call xgBlock_setBlock(slice%Allcprj_work%self, cprj_work_prev , slice%cprjdim, shift_cprj, fcol=1)

    ! ITEST
    ! write(901,*) 'X prev=', xgBlock_getid(X_prev) 
    ! write(901,*) 'cprjX prev=', xgBlock_getid(cprjX_prev) 
    ! write(901,*) 'slice%X before orthoXwrt prev=', xgBlock_getid(slice%X) 
    ! write(901,*) 'slice%AX before orthoXwrt prev=', xgBlock_getid(slice%AX) 
    ! flush(901)
    ! ITEST

    ! call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)

    ! ITEST
    ! write(901,*) 'slice%X after orthoXwrt prev no1=', xgBlock_getid(slice%X) 
    ! write(901,*) 'slice%AX after orthoXwrt prev no1=', xgBlock_getid(slice%AX) 
    ! flush(901)
    ! ITEST
        
    ! call slice_orthoXwrtBlocks(slice, X_prev, cprjX_prev, slice%X, slice%cprjX, islice, cprj_work_prev)
    
    ! ITEST
    ! write(901,*) 'slice%X after orthoXwrt prev no2=', xgBlock_getid(slice%X) 
    ! write(901,*) 'slice%AX after orthoXwrt prev no2=', xgBlock_getid(slice%AX) 
    ! flush(901)
    ! ITEST

    !! ------------------------------------------------------------
    !! 
    !! -               B-ortho kept for reference                 -
    !! 
    !! ------------------------------------------------------------
   
    ! ITEST
    !write(901,*) 'On output: Bortho test'
    !write(901,*) 'slice%AllX before ortho=', xgBlock_getid(slice%AllX) 
    !flush(901)
    ! ITEST
    
    !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
    !   gpu_option,AX=slice%AllAX%self)
    
    ! ITEST
    !write(901,*) 'slice%AllX after Bortho no1=', xgBlock_getid(slice%AllX) 
    !flush(901)
    ! ITEST
    
    !call xg_Borthonormalize_cprj(xg_nonlop,slice%all_blockdim_cprj,slice%AllX,slice%AllcprjX,ierr,tim_ortho,&
    !    gpu_option,AX=slice%AllAX%self)
    
    ! ITEST
    !write(901,*) 'slice%AllX after Bortho no2=', xgBlock_getid(slice%AllX) 
    !flush(901)
    ! ITEST

 end do ! End loop on slices

 !! ------------------------------------------------------------
 !! 
 !! -           Final computation of cprjX, residu, enl        -
 !! 
 !! ------------------------------------------------------------

 call xgBlock_reshape(eigen_out%self, neigenpairs, 1)

 call timab(tim_copy, 1, tsec)
 call xgBlock_copy(X0_out%self, X0)
 call xgBlock_copy(eigen_out%self, eigen)
 call timab(tim_copy, 2, tsec)

 slice%eigenvalues = eigen
 slice%AllX = X0
 slice%AllcprjX = cprjX0

 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%AllX,slice%AllcprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

! A * Psi
 call timab(tim_AX_v,1,tsec)
 call getAX(slice%AllX,slice%AllAX%self)
 call timab(tim_AX_v,2,tsec)
 call timab(tim_AX_k,1,tsec)
 call xgBlock_add_diag(slice%AllX,kin,nspinor,slice%AllAX%self)
 call timab(tim_AX_k,2,tsec)
 call timab(tim_AX_nl,1,tsec)
 call xg_nonlop_getHX(xg_nonlop,slice%AllAX%self,slice%AllcprjX,slice%Allcprj_work%self,slice%proj_work%self)
 call timab(tim_AX_nl,2,tsec)

 ! restart!
! call xg_RayleighRitz_cprj(xg_nonlop,slice%AllX,slice%AllcprjX,slice%AllAX%self,eigen,&
!     slice%all_blockdim_cprj,ierr,15015015,tim_slice_rr,ABI_GPU_DISABLED,solve_ax_bx=.true.)

! if ( ierr /= 0 ) then
!     ABI_BUG("RayleighRitz did not work")
! end if

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
 call xgBlock_reverseMap_1d(residu, resid)
 write(901,*) 'Frobenius norm (merged slices)=', sqrt(sum(resid(1:slice%nbdbuf)))
 write(901,*) 'resid (merged slices)='
 call xgBlock_print(residu, 901)
 write(901,*) 'eigen (merged slices)='
 call xgBlock_print(eigen, 901)
 flush(901)
 ! ITEST

 if (.not.slice%paw) then
   call timab(tim_enl,1,tsec)
   call xg_nonlop_colwiseXHX(xg_nonlop,slice%AllcprjX,slice%Allcprj_work%self,enl)
   call timab(tim_enl,2,tsec)
 end if

 ! Free memory used for slice merging
 call xg_free(X0_out)
 call xg_free(eigen_out)

 call xg_free(dist1)
 call xg_free(dist2)
 call xg_free(dist3)
 ABI_FREE(permute_cols)
 ABI_FREE(rayleigh_quotients)

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

subroutine slice_ampfactorMax(slice,DivResults,lambda_minus,lambda_plus,ndeg_filter_bands)

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
  type(xgBlock_t) :: X_part
  type(xgBlock_t) :: AX_part
  real(dp),pointer :: eig(:,:)

  ! *********************************************************************

  call xgBlock_reverseMap(DivResults,eig,rows=1,cols=cols(DivResults))

  !cheb_poly1(x, n, a, b)
  ampfactor = maxval( (/ (cheb_poly1(eig(1,iband), ndeg_filter_bands(iband), lambda_minus, lambda_plus),& 
      iband=1,cols(DivResults)) /) )

  call xgBlock_scale(slice%X, 1/ampfactor, 1)
  call xgBlock_scale(slice%AX, 1/ampfactor, 1)

end subroutine slice_ampfactorMax
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

end module m_slice_cprj
!!***

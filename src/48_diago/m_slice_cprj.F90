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

 use m_polynomial_filter, only: buildChebyshevJacksonCoeffs
 use m_trace_estimation, only: smallestTridiagEigenpair 

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

   !ARRAYS for entire spectum
   type(xgBlock_t) :: AllX ! Block of initial and final solution
   type(xg_t) :: AllAX     ! space to save AX Hamiltonian application

   ! cprj for entire spectrum
   type(xgBlock_t) :: AllcprjX
   type(xg_t) :: Allcprj_work
   
   !ARRAYS on slice using slice-safe memory
   type(xg_t) :: X_SLICE ! memory independent of the entire spectrum
   type(xgBlock_t) :: X  ! pointers to slice memory
   type(xgBlock_t) :: AX ! 

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

   ! Independent of number of vectors, common to slice and spectrum
   type(xg_t) :: proj_work
   type(xg_nonlop_t) :: xg_nonlop

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
                      nslice,tolfilter,paral_slice,spectral_cut,xg_nonlop,me_g0_fft)

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
 slice%tolerance     = 1.0e-20
 if (tolerance > 0.0) then
   slice%tolerance = tolerance
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
 integer :: deg_i
 integer :: ideg, ierr, j
 integer :: iband, islice
 integer :: slicedim
 integer :: shift_x,shift_cprj
 integer :: niter, max_niter_restart
 integer :: nband_slice
 integer :: count_mask
 integer :: count_rr
 integer :: count_merge
 integer :: count_slice
 integer :: icount
 integer :: blockdim_cprj
 integer :: fcol_in, lcol_in
 integer :: fcol_dummy, lcol_dummy, ncount_out
 integer :: fcol_global
 integer :: ndeg, ndeg_max
 integer :: tim_slice_fi
 integer :: tim_slice_rr
 integer :: tim_slice_pr
 integer :: my_rank
 integer :: trace_degree, trace_rank
 integer :: nstep_spectrum
 integer :: kmax
 integer :: nfilters, ifilter, nstep_bisect, ishift
 integer :: nvec_approx
 integer :: ideg_shift
 integer :: nm
 real(dp) :: balance_prev, balance_this
 real(dp) :: conf_tol
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
 real(dp) :: min_low_est
 real(dp) :: amp_ideg
 real(dp) :: ein_ideg, eout_ideg
 real(dp) :: trace_est
 real(dp) :: trace_est_slice1, trace_est_slice2
 real(dp) :: low_bound, upp_bound, min_low_bound
 real(dp) :: min_upp_bound, max_upp_bound
 real(dp) :: lambda_min, res_norm ! lanczos
 real(dp) :: lower_i, upper_i, mid_i, width
 complex(dp) :: meanz
 real(dp)    :: norm2_, variance
 real(dp) :: tol12 = 1.0e-12
 type(xg_t) :: Xsum
 type(xg_t) :: DivResults
 type(xg_t) :: norm2_X
 type(xg_t) :: dot_XfX, norm2_fX
 type(xg_t) :: res_temp
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
 integer, allocatable :: probe_idx(:)
 integer, allocatable :: nb_vec_slices(:)
 real(dp), allocatable :: upper_bound_slices(:)
 real(dp), allocatable :: rayleigh_quotients(:)
 real(dp), allocatable :: energy_filters_left(:,:)
 real(dp), allocatable :: energy_filters_right(:,:)
 real(dp), allocatable :: lower_bounds(:)
 real(dp), allocatable :: upper_bounds(:)
 real(dp), allocatable, target :: cja(:)
 real(dp), allocatable, target :: energy_interval(:)
 real(dp), allocatable :: confi_interval_left(:), confi_interval_right(:)
 real(dp), pointer :: probe(:) => null()
 real(dp), pointer :: probe_XfX(:,:) => null()
 real(dp), pointer :: X0_norm2(:) => null()
 real(dp), pointer :: lambda_apost(:) => null()
 real(dp), pointer :: lambda_apost_slice(:) => null()
 real(dp), pointer :: theta_(:,:) => null()
 real(dp), pointer :: cheby_moments(:,:) => null()
 real(dp), pointer :: resid(:) => null()
 !Pointers similar to old Chebfi
 type(xg_nonlop_t) :: xg_nonlop
 type(xg_t) :: xgX
 complex(dp), pointer :: Xprobe(:,:) => null() 

! *********************************************************************

 ! ITEST 
 write(901,*) 'inside slice_run'
 write(901,*) 'nslice=', slice%nslice
 flush(901)
 ! ITEST

 ! Warning; the entire code assumes this for simplicity and debugging purposes
 ABI_CHECK(slice%bandpp == slice%neigenpairs, "slice_cprj not implemented in MPI")
 
 ABI_WARNING("experimental code is not debugged")

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
 ABI_MALLOC(permute_cols, (neigenpairs))
 ABI_MALLOC(probe_idx, (neigenpairs))
 ABI_MALLOC(rayleigh_quotients, (neigenpairs))
 ABI_MALLOC(confi_interval_left, (neigenpairs))
 ABI_MALLOC(confi_interval_right, (neigenpairs))
 ABI_MALLOC(nb_vec_slices, (nslice))
 ABI_MALLOC(upper_bound_slices, (nslice))
 ABI_MALLOC(lower_bounds, (nslice))
 ABI_MALLOC(upper_bounds, (nslice))

 ! Memory used to store results of slice merging. 
 call xg_init(X0_out, slice%space, spacedim, neigenpairs, slice%spacecom, me_g0=slice%me_g0)
 call xg_init(eigen_out, SPACE_R, 1, neigenpairs)

 call xg_init(norm2_fX, SPACE_R, neigenpairs, 1)
 call xg_init(norm2_X,SPACE_R,neigenpairs,1)
 call xg_init(dot_XfX, space_res, neigenpairs, 1)

 ! Set initial vectors from input guess
 slice%eigenvalues = eigen

 ! workspace named "All" will store input and output solution
 slice%AllX = X0
 slice%AllcprjX = cprjX0

 ! Compute cprjX for all X (in colsrows)
 call timab(tim_cprj,1,tsec)
 call xg_nonlop_getcprj(xg_nonlop,slice%AllX,slice%AllcprjX,slice%proj_work%self)
 call timab(tim_cprj,2,tsec)

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

 !
 ! ------------------------------------------------------------
 !          Compute Rayleigh quotient for every band
 ! ------------------------------------------------------------
 !

 call timab(tim_RR_q, 1, tsec)
 call slice_rayleighRitzQuotients(slice, maxeig, mineig, DivResults%self)
 call xmpi_max(maxeig,maxeig_global,slice%spacecom,ierr)
 call xmpi_min(mineig,mineig_global,slice%spacecom,ierr)
 call timab(tim_RR_q, 2, tsec)

 ! Recover column indices of sorted Rayleigh quotients in increasing order
 call xgBlock_reverseMap(DivResults%self, theta_, rows=1, cols=neigenpairs)
 rayleigh_quotients(1:neigenpairs) = theta_(1,1:neigenpairs)

 ! ITEST
 write(901,*) 'rayleigh quotients='
 call xgBlock_print(DivResults%self, 901)
 flush(901)
 ! ITEST

 ! Compute |X|^2 colwise L2-norm (before any filter)
 call xgBlock_colwiseNorm2(slice%AllX,norm2_X%self,comm_loc=xmpi_comm_null)
 call xgBlock_reverseMap_1d(norm2_X%self,X0_norm2)

 ! ITEST
 !write(901,*) 'norm2(squared) ||X||='
 !call xgBlock_print(norm2_X%self,901)
 !flush(901)
 ! ITEST
 ! Results: |X|=1 so we don't have to normalize everything after..

 !
 ! ------------------------------------------------------------
 !              Approximate lowest eigenvalue
 ! ------------------------------------------------------------
 !

 kmax = 100
 call computeBLanczos(slice, getAX, kin, spacedim, kmax, lambda_min, res_norm, gpu_option)
 min_low_bound = lambda_min - res_norm

 write(901,*) 'Lanczos lambda_min=', lambda_min
 write(901,*) 'Lanczos res_norm  =', res_norm
 write(901,*) 'Lanczos guarantee =', min_low_bound
 flush(901)

 if (res_norm > 0.1d0) then
     ABI_WARNING("Lanczos has residual > 0.1 may need greater kmax to guarantee lower bound")
 end if

 !
 ! ------------------------------------------------------------
 !        Split spectrum to slices based on spectral gaps
 ! ------------------------------------------------------------
 !
 
 nstep_spectrum = 30
 min_upp_bound = maxval(rayleigh_quotients)
 my_rank = xmpi_comm_rank(slice%spacecom)
 max_upp_bound = slice%ecut
 min_low_est = min_low_bound - min(5*res_norm,0.1d0) ! experimental
 trace_rank = neigenpairs ! FIXME for the moment changing this produces a bug
 trace_degree = 50

 write(901,*) 'splitting spectrum to slices within global', min_low_est, max_upp_bound
 write(901,*) '                                    wanted', min_low_bound, min_upp_bound
 flush(901)
! call splitSpectrumToSlices(slice, nslice, trace_rank, trace_degree, nstep_spectrum, &
!     min_low_bound, min_upp_bound, min_low_est, max_upp_bound, my_rank, getAX, kin, &
!     nb_vec_slices, upper_bound_slices, & ! output
!     gpu_option=gpu_option)
 
 ! Au31 system
 !upper_bound_slices(1) = 0.10d0 ! TODO hardcoded!
 ! Alu system
 upper_bound_slices(1) = 1.37d0 ! hardcoded

 ! Slice 1: 
 my_rank = xmpi_comm_rank(slice%spacecom)
 low_bound = min_low_bound  ! [a,b)
 upp_bound = upper_bound_slices(1)
 trace_degree = 80

 !call computeTraceEstimation(slice, trace_rank, trace_degree, low_bound, upp_bound,&
 !    min_low_est, max_upp_bound, trace_est_slice1, getAX, kin, my_rank, gpu_option=gpu_option)
 write(901,*) 'trace estimation for slice1, deg=', trace_est_slice1, trace_degree
 flush(901)

 ! Slice 2: 
 my_rank = xmpi_comm_rank(slice%spacecom)
 low_bound = upper_bound_slices(1) ! [a,b)
 upp_bound = maxval(rayleigh_quotients)
 overlap_width = (upp_bound - low_bound)/10.0
 low_bound = low_bound - overlap_width
 min_low_est = min_low_est
 max_upp_bound = slice%ecut
 trace_degree = 50
 trace_rank = neigenpairs ! FIXME for the moment changing this produces a bug

 !call computeTraceEstimation(slice, trace_rank, trace_degree, low_bound, upp_bound,&
 !    min_low_est, max_upp_bound, trace_est_slice2, getAX, kin, my_rank, gpu_option=gpu_option)
 !write(901,*) 'trace estimation for slice2, deg=', trace_est_slice2, trace_degree
 !flush(901)
 
 ! attention slice%X est modifié n'est plus X0..
 !! normalement il faut remettre slice%X à valeurs de AllX --->


  !! ------------------------------------------------------------
  !! 
  !! -                    Spectral trees                        -
  !!
  !! Probe nband with low degree just to prune and reduce to n. 
  !! Then actually filter in parallel using high optimal degree. 
  !! Spirit of two-step multiresolution.
  !! 
  !! Spectral multigrid method. Tree traversal method.
  !! Try to detect number of vectors with spectral trees (hierarchy).
  !! Pay high filter degrees only where spectrum mass lies.
  !! 
  !! ------------------------------------------------------------

  ! Reset pointers to dimensions of nband (assumes sequential slices)
  call xg_setBlock(slice%X_SLICE,slice%X,slice%total_spacedim,neigenpairs)
  call xg_setBlock(slice%X_SLICE,slice%AX,slice%total_spacedim,neigenpairs,fcol=neigenpairs+1)

  ! Initial slice workspaces are independent entire arrays
  call xgBlock_copy(slice%AllX, slice%X)
  call xgBlock_copy(slice%AllAX%self, slice%AX)
    
  ! the two following ones will be recomputed so whatever
  slice%cprjX = slice%AllcprjX
  slice%cprj_work = slice%Allcprj_work%self

  ! reinitialize pointers to workspaces ...
  call xg_setBlock(slice%X_NP,slice%X_next,slice%total_spacedim,slice%neigenpairs)
  call xg_setBlock(slice%X_NP,slice%X_prev,slice%total_spacedim,slice%neigenpairs,fcol=slice%neigenpairs+1)

  ndeg_filter_max = 80 ! low degree
  ! if energy oscillates instead of being monotonous this means we have to increase
  ! degree 

  write(901,*) 'Spectral trees %%%%%%%%%'
  write(901,*) 'ndeg_filter_max=', ndeg_filter_max
  write(901,*) 'min_low_est=', min_low_est
  write(901,*) 'max_upp_bound=', max_upp_bound
  write(901,*) ' ********************** '
  flush(901) 

  call xg_init(xgX, slice%space, spacedim, neigenpairs, slice%spacecom, me_g0=slice%me_g0)
  call xgBlock_reverseMap(xgX%self, Xprobe, spacedim, neigenpairs)
  call generateRademacherMatrix(Xprobe, spacedim, neigenpairs, my_rank)
  
  nm = spacedim * neigenpairs
  meanz = sum(Xprobe) / dcmplx(nm,0.0d0)
  norm2_ = sum(abs(Xprobe)**2) / nm
  variance = norm2_ - abs(meanz)**2

  write(901,*) "mean = ", meanz
  write(901,*) "E|z|^2 = ", norm2_
  write(901,*) "variance = ", variance
  call xgBlock_copy(xgX%self, slice%X)
  write(901,*) xgBlock_getid(slice%X)
  flush(901)

  call computeChebyshevMoments(slice, getAX, kin, min_low_est, max_upp_bound, &
      ndeg_filter_max, cheby_moments, gpu_option)

  call xg_free(xgX)

  ! Ugly loop to set slice intervals
  do islice=1, nslice
    lambda_minus = upper_bound_slices(1)
    alpha_minus = upper_bound_slices(1)
        if (islice==1) then
            alpha_minus = min_low_bound
            alpha_plus = lambda_minus
        else
            alpha_plus = maxval(rayleigh_quotients)
        end if
        overlap_width = (alpha_plus - alpha_minus)/10.0
        lambda_plus = alpha_plus+overlap_width
        lambda_minus = alpha_minus-overlap_width
        if (islice==1) then
            lambda_minus = alpha_minus ! otherwise it is outside the center
            lambda_plus = alpha_plus
        end if
        lower_bounds(islice) = lambda_minus
        upper_bounds(islice) = lambda_plus
  end do

  nfilters = nslice
  nstep_bisect = 8 ! number of bisection steps

  ABI_MALLOC(energy_filters_left, (neigenpairs, nstep_bisect))
  ABI_MALLOC(energy_filters_right, (neigenpairs, nstep_bisect))
  ABI_MALLOC(energy_interval, (neigenpairs))
  ABI_MALLOC(cja, (ndeg_filter_max+1))

  upper_bounds(2) = 3.0d0 ! hardcoded FIXME auto
  center = (max_upp_bound + min_low_est)*0.5
  radius = (max_upp_bound - min_low_est)*0.5
  lower_i = lower_bounds(1)
  upper_i = upper_bounds(2)
  width = (upper_i - lower_i) / (nstep_bisect + 1)

  write(901,*)
  write(901,*) 'lower_i=', lower_i
  write(901,*) 'upper_i=', upper_i

  ! total mass
  call buildChebyshevJacksonCoeffs((lower_i-center)/radius, (upper_i-center)/radius, &
      ndeg_filter_max, cja)
  call computeFilterEnergy(neigenpairs, ndeg_filter_max, cja, cheby_moments, &
      energy_interval, nvec_approx)
  write(901,*) 'nvec estimate in total', lower_i, upper_i
  write(901,*) nvec_approx
  flush(901)

  ! shifted bisection
  do ishift = 1, nstep_bisect

    mid_i = lower_i + ishift * width
    deg_i = ndeg_filter_max
     
    write(901,*) '========================================'
    write(901,*) 'ishift=', ishift
    write(901,*) 'mid_i  =', mid_i

    ! Slice Left [a,b)
    call buildChebyshevJacksonCoeffs((lower_i-center)/radius, (mid_i-center)/radius, deg_i, cja)
    call computeFilterEnergy(neigenpairs, deg_i, cja, cheby_moments, &
          energy_interval, nvec_approx)

    energy_filters_left(:,ishift) = energy_interval(:)

    write(901,*) 'nvec estimate left from X0 probe=', nvec_approx
    flush(901)

    ! Slice Right [b,c)
    call buildChebyshevJacksonCoeffs((mid_i-center)/radius, (upper_i-center)/radius, deg_i, cja)
    call computeFilterEnergy(neigenpairs, deg_i, cja, cheby_moments, &
          energy_interval, nvec_approx)

    energy_filters_right(:,ishift) = energy_interval(:)

    write(901,*) 'nvec estimate right from X0 probe=', nvec_approx
    flush(901)

    !! Uniform mass bisection
    !! si le degré est assez elevé alors la masse totale est
    !! - constante
    !! - égale à nband
    !!
    !! Actuellement la masse totale est 88<192. Ça veut dire qu'on pert de la masse. 
    !! Pour ça il faut augmenter l'overlap à droite et à gauche à mon avis.
    !! todo nvec_approx sign flip means that we found important spectral
    !! mass. Do not split there. also follow how the sign flip moves with degree
    !! if it tends to go left or right. Take convergence into account by taking the limit.

  end do 

  ! for every band print its nfilter scores
  write(901,*)
  write(901,*) 'j=    ', nstep_bisect, 'scores'
  do j = 1, neigenpairs
    write(901,*) 'eigenvalue', j
    do ishift = 1, nstep_bisect
        balance_this = energy_filters_left(j,ishift) - energy_filters_right(j,ishift)
        if (ishift > 1) then
            if (balance_prev * balance_this < 0) then
                write(901,*) 'sign flip! significant mass between=', ishift-1, ishift
            end if
        end if
        write(901,*) energy_filters_left(j,ishift), energy_filters_right(j,ishift)
        balance_prev = balance_this
    end do
  end do
  write(901,*)
  flush(901)

  ! sign flip indicates presence of importance spectral mass
  ! we should NOT cut in the (bi,bi+1) that contains the largest spectral mass
  !! adding more nstep_bisect refines the interval (bi,bi+1) containing the large 
  !! spectral mass

  ABI_FREE(energy_filters_left)
  ABI_FREE(energy_filters_right)
  ABI_FREE(cja)

  !! TODO do the tree traversal
  !! like subdivide [A,B) and compute energy is each step

 !! ------------------------------------------------------------
 !! 
 !! -                      Main slice loop                     -
 !! 
 !! How to optimize degree using energies:
 !! Increase degree until the band becomes “energetically isolated” 
 !! and further degree increase changes nothing important.
 !! 
 !! TODO compute Chebyshev recursion for all filters at the same
 !! time with ONE recursion and MULTIPLE Xsum per slice. Then
 !! proceed to Rayleigh-Ritz if slice degree reached. Local operation.
 !! 
 !! ------------------------------------------------------------

 count_merge = 0
 fcol_global = 1

 do islice=1, nslice

    write(901,*)
    write(901,*) '====================Slice=================', islice
    flush(901)

    !! ------------------------------------------------------------
    !! 
    !! -     Initialize slice workspaces from global workspaces
    !!                  Global to Slice operation                 -
    !! 
    !! ------------------------------------------------------------

    ! Reset pointers to dimensions of nband (assumes sequential slices)
    call xg_setBlock(slice%X_SLICE,slice%X,slice%total_spacedim,neigenpairs)
    call xg_setBlock(slice%X_SLICE,slice%AX,slice%total_spacedim,neigenpairs,fcol=neigenpairs+1)

    ! Initial slice workspaces are independent entire arrays
    call xgBlock_copy(slice%AllX, slice%X)
    call xgBlock_copy(slice%AllAX%self, slice%AX)
    
    ! the two following ones will be recomputed so whatever
    slice%cprjX = slice%AllcprjX
    slice%cprj_work = slice%Allcprj_work%self

    ! reinitialize pointers to workspaces ...
    call xg_setBlock(slice%X_NP,slice%X_next,slice%total_spacedim,slice%neigenpairs)
    call xg_setBlock(slice%X_NP,slice%X_prev,slice%total_spacedim,slice%neigenpairs,fcol=slice%neigenpairs+1)

    !! ------------------------------------------------------------
    !! 
    !! -                Scalar polynomial tuning                  -
    !! 
    !! Defining the degree, interval bounds, overlap..
    !! TODO move outside and before islice loop and store parameters
    !! into arrays
    !! 
    !! ------------------------------------------------------------

    ! ongoing, hardcoded depends on previous code
    lambda_minus = upper_bound_slices(1)
    alpha_minus = upper_bound_slices(1)

    if (islice==1) then
        alpha_minus = min_low_bound
        alpha_plus = lambda_minus
    else
        alpha_plus = maxval(rayleigh_quotients)
    end if

    ! overlapping between slices
    overlap_width = (alpha_plus - alpha_minus)/10.0
    lambda_plus = alpha_plus+overlap_width
    lambda_minus = alpha_minus-overlap_width
    if (islice==1) then
        lambda_minus = alpha_minus ! otherwise it is outside the center
        lambda_plus = alpha_plus
    else if (islice==nslice) then
        lambda_plus = alpha_plus ! there is nothing righwise of max anyway..
        ! lambda_plus too wide slows down convergence
    end if

    ! ITEST
    write(901,*) 'global spectrum=', min_low_est, slice%ecut
    write(901,*) 'wanted slice=', alpha_minus, alpha_plus
    write(901,*) 'with overlap=', lambda_minus, lambda_plus
    flush(901)
    ! ITEST

    lower_bounds(islice) = lambda_minus
    upper_bounds(islice) = lambda_plus
 
    center = (slice%ecut + min_low_est)*0.5
    radius = (slice%ecut - min_low_est)*0.5 
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
    if (islice==1) then
        do while ( f_uw/f_u > ramp .and. ndeg < ndeg_max )
            ndeg = ndeg + 1
            f_u  = bandpassIndicator_sca(us_in,ls,us,ndeg)
            f_uw = bandpassIndicator_sca(us   ,ls,us,ndeg)
        end do
    else
        do while ( (f_lw/f_l > ramp .or. f_uw/f_u > ramp) .and. ndeg < ndeg_max )
            ndeg = ndeg + 1
            f_l  = bandpassIndicator_sca(ls_in,ls,us,ndeg)
            f_lw = bandpassIndicator_sca(ls   ,ls,us,ndeg)
            f_u  = bandpassIndicator_sca(us_in,ls,us,ndeg)
            f_uw = bandpassIndicator_sca(us   ,ls,us,ndeg)
        end do
    end if

    ndeg = 90

    ndeg_filter = ndeg
       
    ! ITEST
    write(901,*) 'left/right amplif factor f(out)/f(in)=', f_lw/f_l, f_uw/f_u
    write(901,*) 'minimal polynomial degree=', ndeg_filter
    flush(901)
    ! ITEST

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
    
    ! Initialize Chebyshev expansion of indicator function of order ndeg_filter
    call xg_init(Xsum,slice%space,slice%total_spacedim,neigenpairs,slice%spacecom,gpu_option=gpu_option)
    cdeg = Pi/(ndeg_filter+2)
    mu = 1.d0/Pi*(ACOS(ls)-ACOS(us))
    damp = 1.d0 ! Jackson damping
    one_over_r = 1.d0/radius
    two_over_r = 2.d0/radius
    call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)
    
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

        if (ideg==ndeg_filter - 1) then 
        
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
    
    call timab(tim_slice_fi,2,tsec)

    !! ------------------------------------------------------------
    !! 
    !! -           Vector pruning for subspace basis              -
    !!       
    !! The purpose of this step is to select the nvec vectors that 
    !! have the greatest energy norm for the slice. No tolerance.
    !! We start from all vectors then we keep only nvec of them.
    !! For every column vector x with npw rows,
    !! * f(x) = ||f(M)x|| (option 1)
    !! * f(x) = x^T f(M) x (option 2)
    !! 
    !! ------------------------------------------------------------

    if (islice==1) then
        tim_slice_pr = tim_slice1_pr
    else
        tim_slice_pr = tim_slice2_pr
    end if

    call timab(tim_slice_pr,1,tsec)

    !    Step 1
    ! =============
    ! Compute probe
    ! Note: option 1 gives betten results than option 2 so far
    ! =============

    if (slice%spectral_cut == 1) then
        ! norm2_fX = <fX,fX> colwise L2-norm
        call xgBlock_colwiseNorm2(slice%X, norm2_fX%self, comm_loc=xmpi_comm_null)
        call xgBlock_reverseMap_1d(norm2_fX%self, probe)
    else if (slice%spectral_cut == 2) then
        ! dot_XfX = <X,fX> colwise L2-dot product
        call xgBlock_colwiseDotProduct(X0, slice%X, dot_XfX%self, comm_loc=xmpi_comm_null)
        call xgBlock_reverseMap(dot_XfX%self, probe_XfX, rows=1, cols=neigenpairs)
        probe => probe_XfX(1,1:neigenpairs)
    end if
    
    !write(901,*) 
    !write(901,*) 'probe=', probe(:)
    !flush(901)

    !    Step 2
    ! =============
    ! Probe pruning
    ! =============

    if (islice==1) then
        !count_mask = min(ceiling(trace_est_slice1), neigenpairs)
        !count_mask = ceiling(trace_est_slice1) + 20
        count_mask = 150
    else if (islice==2) then
        !count_mask = ceiling(neigenpairs - trace_est_slice1*0.8)
        count_mask = 150
    end if

    probe = -probe
    probe_idx(1:neigenpairs) = (/ (iband, iband=1,neigenpairs) /)
    call sort_dp(neigenpairs, probe, probe_idx, tol12)
    probe = -probe

    write(901,*) 'kept probes', probe(1:count_mask)

    ! TODO 
    ! deal with extra vectors: if great probes are found outside the kept ones maybe include them
    ! appending vectors should be within a loop here
        
    !    Step 3
    ! ==============================================
    ! Store kept vectors in contiguous memory layout
    ! ==============================================

    ! Allocate slice subspace memory, this is contiguous !!  
    call xg_init(X_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
    call xg_init(AX_kept,slice%space,slice%total_spacedim,count_mask,xmpi_comm_self,me_g0=slice%me_g0_fft)
   
    do icount=1,count_mask

        iband = probe_idx(icount)
        call xgBlock_setBlock(X_kept%self, X_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(AX_kept%self, AX_kept_col, slice%total_spacedim, 1, fcol=icount)
        call xgBlock_setBlock(slice%X, X_col, slice%total_spacedim, 1, fcol=iband)
        call xgBlock_setBlock(slice%AX, AX_col, slice%total_spacedim, 1, fcol=iband)
    
        call timab(tim_copy, 1, tsec)
        call xgBlock_copy(X_col, X_kept_col)
        call xgBlock_copy(AX_col, AX_kept_col)
        call timab(tim_copy, 2, tsec)

    end do
 
    ! reset pointers to temporary (kept)
    slice%X = X_kept%self
    slice%AX = AX_kept%self
    
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

    ! restrict eigenvalue array for size consistency (essentially keep nonzero entries)
    call xgblock_reshape(slice%eigenvalues, 1, neigenpairs) 
    call xgblock_setblock(slice%eigenvalues, eigenvalues_slice, rows=1, cols=count_rr)
    call xgblock_reshape(eigenvalues_slice, count_rr, 1)
    call xgblock_reshape(slice%eigenvalues, neigenpairs, 1)

    ! Orthonormalize (ça fait aucune différence)
    !call xg_Borthonormalize_cprj(xg_nonlop,slice%X,slice%cprjX,ierr,tim_ortho,&
    !    gpu_option,count_rr*xg_nonlop%nspinor,AX=slice%AX)
                  !=blocksize_cprj

    ! Apply Rayleigh Ritz on slice (refinement)
    ! prtvol = 15015015 to print condition number of overlap matrix
    call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,eigenvalues_slice,&
        ierr,15015015,tim_slice_rr,ABI_GPU_DISABLED,solve_ax_bx=.true.)
    
    if ( ierr /= 0 ) then
        ABI_BUG("RayleighRitz did not work")
    end if

    ! restart!
!    call xg_RayleighRitz_cprj(xg_nonlop,slice%X,slice%cprjX,slice%AX,eigenvalues_slice,&
!        ierr,15015015,tim_slice_rr,ABI_GPU_DISABLED,solve_ax_bx=.true.)

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
    write(901,*) 'Slice ', islice, ': colwiseNorm2 residu='
    call xgBlock_print(residu_slice, 901)
    flush(901)
    ! ITEST

    !! ------------------------------------------------------------
    !! 
    !! -              Diagnostic for full slices                  -
    !! 
    !! ------------------------------------------------------------

    call xgBlock_reverseMap_1d(eigenvalues_slice, lambda_apost_slice)

    ! count the confidence intervals within the slice and outside the slice
   
    confi_interval_left(1:count_rr) = lambda_apost_slice + sqrt(resid)
    confi_interval_right(1:count_rr) = lambda_apost_slice - sqrt(resid)
   
    if (islice==1) then
        write(901,*) 'wanted is <', lambda_minus
    else
        write(901,*) 'wanted is ', alpha_minus, alpha_plus
    end if
    flush(901)

    !! TODO étape suivante: une fois qu'on a diagnostiquer une mauvaise convergence
    !! dans une slice on pourrait faire une procédure de restart pour corriger
    !! l'erreur soit en ajoutant plus de vecteurs soit jsp à réflechir

    
    !! ------------------------------------------------------------
    !! 
    !! -                Slice to Global operation                 -
    !!
    !! Merge to X0_out, eigen_out using Rayleigh value as criterion
    !! 
    !! ------------------------------------------------------------
    
    lcol_in = count_rr
    fcol_in = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < alpha_minus)) + 1
   
    conf_tol = maxval(sqrt(resid))
    fcol_dummy = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice+conf_tol < alpha_minus)) + 1
    write(901,*) 'mergeD: fcol slice    ', islice, '        ', fcol_in
    write(901,*) 'mergeD: fcol slice    ', islice, 'interval', fcol_dummy
    ncount_out = count(lambda_apost_slice+conf_tol < alpha_minus)
    write(901,*) 'mergeD: count out left', islice, '        ', ncount_out, 'out of', count_rr
    ncount_out = count(lambda_apost_slice-conf_tol > alpha_plus)
    write(901,*) 'mergeD: count outright', islice, '        ', ncount_out, 'out of', count_rr
    write(901,*) 'mergeD-----------------------'
    flush(901)
    
    if (islice<nslice) then

        lcol_in = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < alpha_plus))
        write(901,*) 'merge: interval on slice   ', islice, ':', lcol_in; flush(901)

    else
        
        write(901,*) 'merge: maximized on slice  ', islice, ':', lcol_in
        lcol_dummy = maxloc(lambda_apost_slice, dim=1, mask=(lambda_apost_slice < alpha_plus))
        write(901,*) 'merge: interval on slice***', islice, ':', lcol_in
        write(901,*) 'merge: resid***************', islice, ':', sqrt(sum(resid(fcol_in:lcol_dummy)))
        flush(901)

        if (count_merge + lcol_in-fcol_in+1 > neigenpairs) then
            write(901,*) 'merge: pass on slice       ', islice, ':', lcol_in; flush(901)
            lcol_in = fcol_in + neigenpairs - count_merge - 1
            write(901,*) 'merge: pass after on slice ', islice, ':', lcol_in, fcol_in, neigenpairs, count_merge
            write(901,*) 'merge: residBBBBBBBBBBBBBBB', islice, ':', sqrt(sum(resid(fcol_in:lcol_in)))
            flush(901)
        end if

    end if

    count_slice = lcol_in - fcol_in + 1

    ! ITEST
    write(901,*) 'Frobenius norm (inside slice', islice, 'only)=', sqrt(sum(resid(fcol_in:lcol_in)))
    flush(901)
    ! ITEST

    ! ITEST
    write(901,*) 
    write(901,*) '====================================== Slice', islice
    write(901,*) 'interval bounds', lambda_minus, alpha_minus, alpha_plus
    !write(901,*) 'count_merged=', count_merge 
    write(901,*) 'count_slice =', count_slice
    write(901,*) 'count_rr    =', count_rr
    write(901,*) 'fcol_in,val =', fcol_in, lambda_apost_slice(fcol_in)
    write(901,*) 'lcol_in,val =', lcol_in, lambda_apost_slice(lcol_in)
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
    call xg_free(Xsum)
    call xg_free(X_kept)
    call xg_free(AX_kept)

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
 write(901,*) 'Frobenius norm (merged slices)=', sqrt(sum(resid))
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

 call xg_free(norm2_X)
 call xg_free(norm2_fX)
 call xg_free(dot_XfX)
 call xg_free(norm2_fX)
 ABI_FREE(permute_cols)
 ABI_FREE(rayleigh_quotients)
 ABI_FREE(probe_idx)
 ABI_FREE(confi_interval_left)
 ABI_FREE(confi_interval_right)
 ABI_FREE(nb_vec_slices)
 ABI_FREE(upper_bound_slices)
 ABI_FREE(lower_bounds)
 ABI_FREE(upper_bounds)

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
      call xg_nonlop_getSX(slice%xg_nonlop,cprjX0,cprjvar,cprj_work,buffer%self)
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
  integer         :: npw, nband
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

!----------------------------------------------------------------------

!!****f* m_slice_cprj/generateRademacherMatrix
!! NAME
!! generateRademacherMatrix
!!
!! SOURCE
subroutine generateRademacherMatrix(V, n, m, rank)
  
    implicit none
    
    ! input/output
    integer, intent(in) :: n, m, rank
    complex(dp), intent(out) :: V(n,m)
    ! local arguments
    integer :: nseed, i, j, k
    integer :: base_seed
    real(dp) :: R(n,m)
    integer, allocatable :: seed(:)
    
    ! *********************************************************************

    ! MPI-safe seed: deterministic way to generate a unique seed per MPI rank
    call random_seed(size=nseed) ! runtime value of nseed
    
    ABI_MALLOC(seed, (nseed))
    base_seed = 123456789
    seed = mod( base_seed + rank*73856093 + [(i*19349663, i=1,nseed)], 2147483647 )
    call random_seed(put=seed)
  
    call random_number(R)
    do j = 1, m
        do i = 1, n
            k = int(4.0d0 * R(i,j))      ! 0,1,2,3
            V(i,j) = dcmplx( cos(0.5d0*pi*k), sin(0.5d0*pi*k))
        end do
    end do

    ABI_FREE(seed)

end subroutine generateRademacherMatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/generateGaussianMatrix
!! NAME
!! generateGaussianMatrix
!!
!! SOURCE
subroutine generateGaussianMatrix(V, n, m, rank)

    implicit none

    ! input/output
    integer, intent(in)  :: n, m, rank
    real(dp), intent(out) :: V(2, n*m)

    ! local arguments
    integer :: nseed, i, j, nm
    integer :: base_seed
    integer, allocatable :: seed(:)
    real(dp) :: u1, u2

    ! *********************************************************************

    nm = n * m

    ! MPI-safe seed: deterministic way to generate a unique seed per MPI rank
    call random_seed(size = nseed)

    ABI_MALLOC(seed, (nseed))
    base_seed = 123456789
    seed = mod( base_seed + rank*73856093 + [(i*19349663, i=1,nseed)], 2147483647 )
    call random_seed(put = seed)

    ! Generate i.i.d. N(0,1) entries (real-valued)
    i = 1
    do while (i <= nm)
        call random_number(u1)
        call random_number(u2)

        ! Box–Muller transform
        V(1, i) = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * Pi * u2)

        if (i + 1 <= nm) then
            V(1, i+1) = sqrt(-2.0_dp * log(u1)) * sin(2.0_dp * Pi * u2)
        end if

        i = i + 2
    end do

    ! Imaginary part = 0 (consistent with your Rademacher routine)
    V(2, 1:nm) = 0.0_dp

    ABI_FREE(seed)

end subroutine generateGaussianMatrix
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/applyLowpassFilter
!! NAME
!! applyLowpassFilter
!!
!! SOURCE
subroutine applyLowpassFilter(slice, getAX, ampl_factors, kin, low_bound, upp_bound, &
        ndeg_filter, gpu_option)

    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    type(xgBlock_t), intent(in) :: ampl_factors
    integer, intent(in) :: ndeg_filter
    real(dp), intent(in) :: low_bound, upp_bound
    integer, optional, intent(in) :: gpu_option
    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    integer :: neigenpairs
    integer :: ideg
    integer :: nspinor
    integer :: l_gpu_option
    real(dp) :: center, radius
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: tsec(2)
    integer, allocatable :: ndeg_filter_bands(:) !Oracle variable
    type(xg_nonlop_t) :: xg_nonlop

    ! *********************************************************************

    ! Initialize values
    neigenpairs = cols(slice%X)
    xg_nonlop = slice%xg_nonlop
    nspinor = slice%xg_nonlop%nspinor
    l_gpu_option = ABI_GPU_DISABLED
    
    ! Process input arguments
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    ! Spectral interval to be amplified scaled to [-1,1)
    center = (upp_bound + low_bound)*0.5
    radius = (upp_bound - low_bound)*0.5 
    one_over_r = 1.0/radius
    two_over_r = 2.0/radius
    
    ! Loop on degree 
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

    ! Amplify to finalize filter application
    call timab(tim_amp_f,1,tsec)
    ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
    ndeg_filter_bands(:) = ndeg_filter
    call slice_ampfactor(slice, ampl_factors, low_bound, upp_bound, ndeg_filter_bands)
    ABI_FREE(ndeg_filter_bands)
    call timab(tim_amp_f,2,tsec)

end subroutine applyLowpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/applyBandpassFilter
!! NAME
!! applyBandpassFilter
!!
!! SOURCE
subroutine applyBandpassFilter(slice, getAX, kin, low_bound, upp_bound, &
        min_low_bound, max_upp_bound, ndeg_filter, gpu_option)

    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: ndeg_filter
    real(dp), intent(in) :: low_bound, upp_bound
    real(dp), intent(in) :: min_low_bound, max_upp_bound
    integer, optional, intent(in) :: gpu_option
    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    integer :: neigenpairs
    integer :: ideg
    integer :: nspinor
    integer :: l_gpu_option
    real(dp) :: cdeg
    real(dp) :: center, radius
    real(dp) :: ls, us
    real(dp) :: mu, damp
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: tsec(2)
    type(xg_t) :: Xsum
    type(xg_nonlop_t) :: xg_nonlop

    ! *********************************************************************

    ! Initialize values
    neigenpairs = cols(slice%X)
    xg_nonlop = slice%xg_nonlop
    nspinor = slice%xg_nonlop%nspinor
    l_gpu_option = ABI_GPU_DISABLED
    
    ! Process input arguments
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    ! Allocate space
    call xg_init(Xsum,slice%space,slice%total_spacedim,neigenpairs,&
        slice%spacecom,gpu_option=l_gpu_option) 

    ! Spectral interval to be amplified scaled to [-1,1)
    center = (max_upp_bound + min_low_bound)*0.5
    radius = (max_upp_bound - min_low_bound)*0.5 
    ls = (low_bound - center) / radius
    us = (upp_bound - center) / radius
    one_over_r = 1.0/radius
    two_over_r = 2.0/radius
    
    ! Initialize Chebyshev expansion of indicator function of order ndeg_filter
    cdeg = Pi/(ndeg_filter+2)
    mu = 1.d0/Pi*(ACOS(ls)-ACOS(us))
    damp = 1.d0 ! Jackson damping
    call xgBlock_saxpy(Xsum%self, mu*damp, slice%X)

    ! Loop on degree 
    do ideg = 0, ndeg_filter - 1
        
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

        if (ideg==ndeg_filter - 1) then 
        
            ! store final expansion Xsum to X (inplace rewrites X)
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

    ! Free memory
    call xg_free(Xsum)
    
end subroutine applyBandpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/computeChebyshevMoments
!! NAME
!! computeChebyshevMoments
!! 
!! FUNCTION
!! Compute Chebyshev moments up to maximal degree all centered in [A,B)
!!
!! OUTPUT
!! M_n = <X, f_n(B^{-1}AX) X> for n=1,..,ndeg_filter_max
!! 
!! SOURCE
subroutine computeChebyshevMoments(slice, getAX, kin, min_low_bound, max_upp_bound, &
        ndeg_filter_max, cheby_moments, gpu_option)

    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: ndeg_filter_max
    real(dp), intent(in) :: min_low_bound, max_upp_bound
    real(dp), pointer, intent(inout) :: cheby_moments(:,:)
    integer, optional, intent(in) :: gpu_option
    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    integer :: neigenpairs
    integer :: ideg
    integer :: nspinor
    integer :: nrows
    integer :: l_gpu_option
    real(dp) :: center, radius
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: tsec(2)
    type(xg_t) :: Moments, X0
    type(xgBlock_t) :: moment_ideg
    type(xg_nonlop_t) :: xg_nonlop

    ! *********************************************************************

    ! Initialize values
    neigenpairs = cols(slice%X)
    nrows = slice%total_spacedim
    xg_nonlop = slice%xg_nonlop
    nspinor = slice%xg_nonlop%nspinor
    l_gpu_option = ABI_GPU_DISABLED
    
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    ! Allocate space
    call xg_init(Moments, slice%space, neigenpairs, ndeg_filter_max+1)        ! M_n=<X0,f_n(A)X0>
    call xg_init(X0, slice%space, nrows, neigenpairs, slice%spacecom, me_g0=slice%me_g0)  ! X0

    call xgBlock_copy(slice%X, X0%self)

    ! Spectral interval to be amplified scaled to [-1,1)
    center = (max_upp_bound + min_low_bound)*0.5
    radius = (max_upp_bound - min_low_bound)*0.5 
    one_over_r = 1.0/radius
    two_over_r = 2.0/radius

    ! Initialize
    call xgBlock_setBlock(Moments%self, moment_ideg, neigenpairs, 1) 
    call xgBlock_colwiseDotProduct(slice%X, slice%X, moment_ideg, comm_loc=xmpi_comm_null)
   
    ! Loop on degree 
    do ideg = 0, ndeg_filter_max - 1
        
        call timab(tim_cprj,1,tsec)
        call xg_nonlop_getcprj(xg_nonlop,slice%AX,slice%cprjX,slice%proj_work%self)
        call timab(tim_cprj,2,tsec)

        call slice_computeNextOrderChebfiPolynom(slice, ideg, center, one_over_r, two_over_r)

        ! slice%X = f_ideg X0
        call timab(tim_swap,1,tsec)
        call slice_swapInnerBuffers(slice, slice%total_spacedim, neigenpairs)
        call timab(tim_swap,2,tsec)

        ! M_ideg = < X0, f_ideg X0 > in C^nband for every ideg
        call xgBlock_setBlock(Moments%self, moment_ideg, neigenpairs, 1, fcol=ideg+2) 
        call xgBlock_colwiseDotProduct(X0%self, slice%X, moment_ideg, comm_loc=xmpi_comm_null)

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

    call xgBlock_reverseMap(Moments%self, cheby_moments, neigenpairs, ndeg_filter_max+1)
  
    ! Free memory
    call xg_free(X0)
    call xg_free(Moments)
    
end subroutine computeChebyshevMoments
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/computeTraceEstimation
!! NAME
!! computeTraceEstimation
!! 
!! FUNCTION
!! Compute Girard-Hutchinson trace estimator
!!
!! SOURCE
subroutine computeTraceEstimation(slice, m_vecs, trace_degree, low_bound, upp_bound, &
        min_low_bound, max_upp_bound, trace_est, getAX, kin, my_rank, gpu_option)
  
    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: my_rank
    integer, intent(in) :: trace_degree
    integer, intent(in) :: m_vecs
    real(dp), intent(out) :: trace_est
    real(dp), intent(in) :: low_bound, upp_bound
    real(dp), intent(in) :: min_low_bound, max_upp_bound
    integer, optional, intent(in) :: gpu_option

    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    type(xg_t) :: dot_XfX
    type(xg_t) :: xgX
    type(xg_nonlop_t) :: xg_nonlop
    integer :: cprjdim
    integer :: spacedim, i, j, k
    integer :: blockdim_cprj
    integer :: nspinor
    integer :: l_gpu_option
    integer :: idx_i, idx_j
    real(dp) :: tolerance
    real(dp) :: trace_tmp
    real(dp) :: normX
    real(dp), pointer :: accum(:,:) => null()
    complex(dp), pointer :: X(:,:) => null() 
    real(dp) :: tsec(2)

    ! *********************************************************************

    spacedim = slice%total_spacedim
    tolerance = slice%tolerance
    cprjdim = slice%cprjdim
    nspinor = slice%xg_nonlop%nspinor
    blockdim_cprj = m_vecs*nspinor
    xg_nonlop = slice%xg_nonlop

    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    ! Allocate contiguous memory blocks
    call xg_init(dot_XfX, slice%space, m_vecs, 1)                                        ! <X,f(A)X>
    call xg_init(xgX, slice%space, spacedim, m_vecs, slice%spacecom, me_g0=slice%me_g0)  ! X

    call xgBlock_reverseMap(xgX%self, X, spacedim, m_vecs)

    ! Fill entries of random matrices
    call generateRademacherMatrix(X, spacedim, m_vecs, my_rank)
    !call generateGaussianMatrix(X, spacedim, m_vecs, my_rank)
    
    ! Initialize X (dimensions then values)
    call xg_setBlock(slice%X_SLICE, slice%X, spacedim, m_vecs)
    call xgBlock_copy(xgX%self, slice%X)
    
    ! Initialize Xprev, Xnext (dimensions)
    call xg_setBlock(slice%X_NP,slice%X_next,spacedim,m_vecs)
    call xg_setBlock(slice%X_NP,slice%X_prev,spacedim,m_vecs,fcol=m_vecs+1)

    ! Initialize cprj (dimensions and values)
    call xgBlock_setBlock(slice%AllcprjX, slice%cprjX, slice%cprjdim, blockdim_cprj)
    call xgBlock_setBlock(slice%Allcprj_work%self, slice%cprj_work, slice%cprjdim, blockdim_cprj)

    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop,slice%X,slice%cprjX,slice%proj_work%self)
    call timab(tim_cprj,2,tsec)

    ! Initialize A * Psi (dimensions and values)
    call xg_setBlock(slice%X_SLICE, slice%AX, spacedim, m_vecs, fcol=m_vecs+1)
    
    call timab(tim_AX_v,1,tsec)
    call getAX(slice%X,slice%AX)
    call timab(tim_AX_v,2,tsec)
    call timab(tim_AX_k,1,tsec)
    call xgBlock_add_diag(slice%X,kin,nspinor,slice%AX)
    call timab(tim_AX_k,2,tsec)
    call timab(tim_AX_nl,1,tsec)
    call xg_nonlop_getHX(xg_nonlop,slice%AX,slice%cprjX,slice%cprj_work,slice%proj_work%self)
    call timab(tim_AX_nl,2,tsec)

    ! Compute f(A) * X
    call applyBandpassFilter(slice, getAX, kin, low_bound, upp_bound, &
        min_low_bound, max_upp_bound, trace_degree, l_gpu_option)
    
    call xgBlock_colwiseDotProduct(xgX%self, slice%X, dot_XfX%self, comm_loc=xmpi_comm_null)
    call xgBlock_reverseMap(dot_XfX%self, accum, m_vecs, 1)
   
    trace_est = sum(accum) / m_vecs
    
    ! Free memory
    call xg_free(dot_XfX)
    call xg_free(xgX)

end subroutine computeTraceEstimation
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/splitSpectrumToSlices
!! NAME
!! splitSpectrumToSlices
!!
!! SOURCE
subroutine splitSpectrumToSlices( &
        slice, nslice, trace_rank, trace_degree, nstep_spectrum, &
        low_bound_wanted, upp_bound_wanted, min_bound, max_bound, my_rank, &
        getAX, kin, nb_vec_slices, upp_bound_slices, gpu_option)

    implicit none
 
    ! Arguments
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: nslice
    integer, intent(in) :: trace_rank
    integer, intent(in) :: trace_degree
    integer, intent(in) :: nstep_spectrum
    integer, intent(in) :: my_rank
    real(dp), intent(in) :: low_bound_wanted
    real(dp), intent(in) :: upp_bound_wanted
    real(dp), intent(in) :: min_bound
    real(dp), intent(in) :: max_bound
    integer, intent(inout) :: nb_vec_slices(nslice)
    real(dp), intent(inout) :: upp_bound_slices(nslice)
    integer, optional, intent(in) :: gpu_option
    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    ! Local variables
    integer, parameter :: trace_crit = 2
    integer :: l_gpu_option
    integer :: trace_sum
    integer :: neigenpairs
    integer :: ipart
    integer :: mid, idx, k, n
    logical :: found_gap(nstep_spectrum)
    real(dp) :: width
    real(dp) :: trace_est
    real(dp) :: lower_i
    real(dp) :: upper_i
    real(dp) :: trace_estim_spectrum(nstep_spectrum)
    real(dp) :: upper_bound_spectrum(nstep_spectrum)
    
    ! *********************************************************************

    neigenpairs = slice%neigenpairs
    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    if (low_bound_wanted < min_bound) then
        ABI_ERROR('wanted spectrum falls outside given lower bound')
    end if
    if (upp_bound_wanted > max_bound) then
        ABI_ERROR('wanted spectrum falls outside given upper bound')
    end if
 
    ! Split working spectrum into N intervals (N=nstep_spectrum)
    width = (upp_bound_wanted - low_bound_wanted) / nstep_spectrum
    trace_sum = 0
 
    do ipart=1, nstep_spectrum

         lower_i = min(max_bound, low_bound_wanted + (ipart - 1) * width)
         upper_i = min(max_bound, low_bound_wanted + ipart * width)

         !call computeTraceEstimation(slice, trace_rank, trace_degree, lower_i, upper_i,&
         !   min_bound, max_bound, trace_est, getAX, kin, my_rank, gpu_option=l_gpu_option)
        
        ! TODO keep this version and add second version that reuses Chebyshev 
        !      recursion to avoid repeated calculations... 

        found_gap(ipart) = ( ipart>1 .and. trace_est < trace_crit )

        write(901,*) 'spectrum part : nbvecs', ipart, lower_i, upper_i, ceiling(trace_est)
        flush(901)
        
        upper_bound_spectrum(ipart) = upper_i
        trace_estim_spectrum(ipart) = ceiling(trace_est)
        
        trace_sum = trace_sum + ceiling(trace_est)

    end do

    if (trace_sum < slice%neigenpairs) then
        ABI_WARNING("trace estimation missed eigenvalues")
    end if

    ! Chose the spectral gap with the index closer to the middle
    ! same for interval bound closer to the midpoint
    write(901,*) "found_gap=", found_gap
    write(901,*) "upper_bound=", upper_bound_spectrum
    flush(901)

    found_gap = ( trace_estim_spectrum == minval(trace_estim_spectrum) )
    
    write(901,*) "found_gap(mod)=", found_gap

    n = nstep_spectrum
    mid = (n + 1) / 2
    idx = 0   ! 0 = not found

    if (found_gap(mid)) then
        idx = mid
    else
        do k = 1, max(mid-1, n-mid)
            if (mid-k >= 1) then
                if (found_gap(mid-k)) then
                    idx = mid-k
                    exit
                end if
            end if
            if (mid+k <= n) then
                if (found_gap(mid+k)) then
                    idx = mid+k
                    exit
                end if
            end if
        end do
    end if

    nb_vec_slices(1) = min(ceiling(sum(trace_estim_spectrum(1:idx+1))), neigenpairs)
    nb_vec_slices(2) = min(ceiling(sum(trace_estim_spectrum(idx+2:nstep_spectrum))), neigenpairs)
    write(901,*) 'gap idx    =', idx
    write(901,*) 'gap value  =', upper_bound_spectrum(idx)
    write(901,*) 'trace left =', nb_vec_slices(1)
    write(901,*) 'trace right=', nb_vec_slices(2)
    flush(901)

    ! output the gap
    upp_bound_slices(1) = upper_bound_spectrum(idx)


end subroutine splitSpectrumToSlices
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/computeBLanczos
!! NAME
!! computeBLanczos
!! 
!! FUNCTION
!! B-Lanczos three-term recurrence (using B-inner product)
!!
!! SOURCE
  
  subroutine computeBLanczos(slice, getAX, kin, n, k, lambda_min, res_norm, gpu_option)

    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: n, k
    real(dp), intent(out) :: lambda_min, res_norm
    integer, optional, intent(in) :: gpu_option
    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface
    
    type(xg_nonlop_t) :: xg_nonlop
    type(xg_t) :: xg_Bv
    real(dp) :: Bv(n), Bv2(2,n)
    real(dp) :: q(n), q2(2,n), v(n), v2(2,n)
    real(dp) :: Bm1v2(2,n)
    real(dp) :: alpha(k), beta(k-1)
    real(dp) :: q_prev(n)
    real(dp) :: v_min(n)
    real(dp) :: beta_prev
    real(dp) :: normB

    integer :: i, j
    integer :: space
    integer :: spacedim
    integer :: l_gpu_option

    ! *********************************************************************

    xg_nonlop = slice%xg_nonlop
    space = slice%space
    spacedim = slice%spacedim

    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if

    ! ONGOING
    ! unit tests per elementary operation
    
    call random_number(q)
    do i = 1, n
        q2(1, i) = q(i)
        q2(2, i) = 0.0d0
    end do
    v2(:,:) = 0.0d0; v(:) = 0.0d0
    Bv2(:,:) = 0.0d0; Bv(:) = 0.0d0
    Bm1v2(:,:) = 0.0d0
    
    ! Bv = B * q
    call matmul_op_B(slice, n, q2, Bv2, l_gpu_option)
    do i = 1, n
        Bv(i) = Bv2(1, i) ! avoids temporary
    end do

    normB = sqrt(dot_product(q, Bv))
    q = q / normB

    q_prev = 0.0_dp
    beta_prev = 0.0_dp

    do j = 1, k
        ! v = A * q    
        do i = 1, n
            q2(1, i) = q(i)
        end do
        call matmul_op_A(slice, getAX, kin, n, q2, v2, l_gpu_option)
        do i = 1, n
            v(i) = v2(1, i)
        end do

        ! alpha_j = q^T * Aq
        alpha(j) = dot_product(q, v)

        ! v = B^{-1} * A * q
        call matmul_op_Binv(slice, n, v2, Bm1v2, l_gpu_option) 
        do i = 1, n
            v(i) = Bm1v2(1, i)
        end do

        ! v = B^{-1} A q - alpha q - beta_prev q_prev
        v = v - alpha(j)*q
        if (j > 1) v = v - beta_prev*q_prev

        ! Compute beta_j if j<k
        if (j < k) then
            ! Bv = B * v        
            do i = 1, n
                v2(1, i) = v(i)
            end do
            call matmul_op_B(slice, n, v2, Bv2, l_gpu_option)
            do i = 1, n
                Bv(i) = Bv2(1, i)
            end do

            beta(j) = sqrt(dot_product(v, Bv))

            ! Update q_prev, q, beta_prev
            q_prev = q
            q = v / beta(j)
            beta_prev = beta(j)
        end if
    end do

    ! Diagonalize T
    call smallestTridiagEigenpair(k, alpha, beta, lambda_min, v_min)

    ! residual norm using Lanczos shortcut
    res_norm = abs(beta(k-1)*v_min(k))

  end subroutine computeBLanczos
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/matmul_op_A
!! NAME
!! matmul_op_A
!! 
!! FUNCTION
!! Apply A to Fortran 1D array q2 and store result to Fortran 1D array Aq2
!! All in complex.
!! 
!! SOURCE

subroutine matmul_op_A(slice, getAX, kin, n, q2, Aq2, gpu_option)
    
    implicit none

    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(in) :: kin
    integer, intent(in) :: n
    real(dp), intent(in) :: q2(2,n) ! 1d complex array (in)
    real(dp), intent(out) :: Aq2(2,n) ! 1d complex array (out)
    integer, optional, intent(in) :: gpu_option

    interface
        subroutine getAX(X,AX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
        end subroutine getAX
    end interface

    type(xg_t) :: cprjW, cprj_work
    type(xgBlock_t) :: W, AW, proj_work
    type(xg_nonlop_t) :: xg_nonlop
    integer :: l_gpu_option, space, spacecom, space_cprj, cprjdim, nspinor
    integer :: blockdim_cprj
    real(dp) :: tsec(2)
    
    ! *********************************************************************

    space = slice%space
    spacecom = slice%spacecom
    space_cprj = slice%space_cprj
    cprjdim = slice%cprjdim
    nspinor = slice%xg_nonlop%nspinor
    xg_nonlop = slice%xg_nonlop
    proj_work = slice%proj_work%self

    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if
    
    ! TODO make this optional
    blockdim_cprj = nspinor ! = number_columns * nspinor
    call xg_init(cprjW, space_cprj, cprjdim, blockdim_cprj, spacecom)
    call xg_init(cprj_work, space_cprj, cprjdim, blockdim_cprj, spacecom)
   
    ! TODO make this optional otherwise W=slice%X and AW=slice%AX
    call xgBlock_map(W, q2, space, n, 1, spacecom, gpu_option=gpu_option)
    call xgBlock_map(AW, Aq2, space, n, 1, spacecom, gpu_option=gpu_option)

    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop, W, cprjW%self, proj_work)
    call timab(tim_cprj,2,tsec)
        
    call timab(tim_ax_v,1,tsec)
    call getAX(W,AW)
    call timab(tim_ax_v,2,tsec)

    call timab(tim_ax_k,1,tsec)
    call xgBlock_add_diag(W,kin,nspinor,AW)
    call timab(tim_ax_k,2,tsec)

    call timab(tim_AX_nl,1,tsec)
    call xg_nonlop_getHX(xg_nonlop, AW, cprjW%self, cprj_work%self, proj_work)
    call timab(tim_AX_nl,2,tsec)

    call xg_free(cprjW)
    call xg_free(cprj_work)

end subroutine matmul_op_A
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/matmul_op_B
!! NAME
!! matmul_op_B
!! 
!! FUNCTION
!! Apply B to Fortran 1D array q2 and store result to Fortran 1D array Bq2
!! All in complex.
!! 
!! SOURCE

subroutine matmul_op_B(slice, n, q2, Bq2, gpu_option)
    
    implicit none

    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: n
    real(dp), intent(in) :: q2(2,n) ! 1d complex array (in)
    real(dp), intent(out) :: Bq2(2,n) ! 1d complex array (out)
    integer, optional, intent(in) :: gpu_option

    type(xg_t) :: cprjW, cprj_work
    type(xgBlock_t) :: W, BW, proj_work
    type(xg_nonlop_t) :: xg_nonlop
    integer :: l_gpu_option, space, spacecom, space_cprj, cprjdim, nspinor
    integer :: blockdim_cprj
    real(dp) :: tsec(2)
    
    ! *********************************************************************

    space = slice%space
    spacecom = slice%spacecom
    space_cprj = slice%space_cprj
    cprjdim = slice%cprjdim
    nspinor = slice%xg_nonlop%nspinor
    xg_nonlop = slice%xg_nonlop
    proj_work = slice%proj_work%self

    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if 
   
    ! TODO make this optional
    blockdim_cprj = nspinor ! = number_columns * nspinor
    call xg_init(cprjW, space_cprj, cprjdim, blockdim_cprj, spacecom)
    call xg_init(cprj_work, space_cprj, cprjdim, blockdim_cprj, spacecom)
    
    ! TODO make this optional otherwise W=slice%X and AW=slice%AX
    call xgBlock_map(W, q2, space, n, 1, spacecom, gpu_option=gpu_option)
    call xgBlock_map(BW, Bq2, space, n, 1, spacecom, gpu_option=gpu_option)

    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop, W, cprjW%self, proj_work)
    call timab(tim_cprj,2,tsec)

    ! BW = S|Psi>
    call timab(tim_copy, 1, tsec)
    call xgBlock_copy(W,BW)
    call timab(tim_copy, 2, tsec)

    if (slice%paw) then
        call xg_nonlop_getSX(xg_nonlop, BW, cprjW%self, cprj_work%self, proj_work)
    end if

    call xg_free(cprjW)
    call xg_free(cprj_work)

end subroutine matmul_op_B
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/matmul_op_Binv
!! NAME
!! matmul_op_Binv
!! 
!! FUNCTION
!! Apply Binv to Fortran 1D array q2 and store result to Fortran 1D array Binvq2
!! All in complex.
!! 
!! SOURCE

subroutine matmul_op_Binv(slice, n, q2, Binvq2, gpu_option)
    
    implicit none

    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: n
    real(dp), intent(in) :: q2(2,n) ! 1d complex array (in)
    real(dp), intent(out) :: Binvq2(2,n) ! 1d complex array (out)
    integer, optional, intent(in) :: gpu_option

    type(xg_t) :: cprjW, cprj_work, cprj_work2
    type(xgBlock_t) :: W, BinvW, proj_work
    type(xg_nonlop_t) :: xg_nonlop
    integer :: l_gpu_option, space, spacecom, space_cprj, cprjdim, nspinor
    integer :: blockdim_cprj
    real(dp) :: tsec(2)
    
    ! *********************************************************************

    space = slice%space
    spacecom = slice%spacecom
    space_cprj = slice%space_cprj
    cprjdim = slice%cprjdim
    nspinor = slice%xg_nonlop%nspinor
    xg_nonlop = slice%xg_nonlop
    proj_work = slice%proj_work%self

    l_gpu_option = ABI_GPU_DISABLED
    if (present(gpu_option)) then
      l_gpu_option = gpu_option
    end if
    
    blockdim_cprj = nspinor ! = number_columns * nspinor
    call xg_init(cprjW, space_cprj, cprjdim, blockdim_cprj, spacecom)
    call xg_init(cprj_work, space_cprj, cprjdim, blockdim_cprj, spacecom)
    call xg_init(cprj_work2, space_cprj, cprjdim, blockdim_cprj, spacecom)
    
    call xgBlock_map(W, q2, space, n, 1, spacecom, gpu_option=gpu_option)
    call xgBlock_map(BinvW, Binvq2, space, n, 1, spacecom, gpu_option=gpu_option)

    call timab(tim_cprj,1,tsec)
    call xg_nonlop_getcprj(xg_nonlop, W, cprjW%self, proj_work)
    call timab(tim_cprj,2,tsec)

    ! BinvW = S^{-1}|Psi>
    call timab(tim_copy, 1, tsec)
    call xgBlock_copy(W,BinvW)
    call timab(tim_copy, 2, tsec)

    if (slice%paw) then
        call timab(tim_invovl, 1, tsec)
        call xg_nonlop_getSm1X(xg_nonlop, BinvW, cprjW%self, cprj_work%self, &
            cprj_work2%self,proj_work)
        call timab(tim_invovl, 2, tsec)
    end if

    call xg_free(cprjW)
    call xg_free(cprj_work)
    call xg_free(cprj_work2)

end subroutine matmul_op_Binv
!!***

!----------------------------------------------------------------------

!!****f* m_slice_cprj/computeFilterEnergy
!! NAME
!! computeFilterEnergy
!! 
!! SOURCE

subroutine computeFilterEnergy(neigenpairs, ndeg_filter, cja, cheby_moments, &
        energy_interval, nvec_approx)

      implicit none

      integer, intent(in) :: ndeg_filter
      integer, intent(in) :: neigenpairs
      integer, intent(out) :: nvec_approx
      real(dp), pointer, intent(in) :: cja(:)
      real(dp), pointer, intent(in) :: cheby_moments(:,:)
      real(dp), pointer, intent(in) :: energy_interval(:)
 
      real(dp) :: E_re(neigenpairs)
      real(dp) :: E_im(neigenpairs)
      real(dp) :: E_tot
      real(dp) :: energy_magn_j
      integer :: ideg
      integer :: j

      ! *********************************************************************

      E_re = 0.0d0
      E_im = 0.0d0
      E_tot = 0.d0
      do ideg = 1, ndeg_filter+1
          do j = 1, neigenpairs
              E_re(j) = E_re(j) + cja(ideg) * cheby_moments(2*j-1, ideg)
              E_im(j) = E_im(j) + cja(ideg) * cheby_moments(2*j,   ideg)
          end do
      end do
      do j = 1, neigenpairs
          !energy_magn_j = hypot(E_re(j), E_im(j))
          energy_magn_j = real(E_re(j))
          energy_interval(j) = energy_magn_j
          E_tot = E_tot + energy_magn_j 
      end do
      nvec_approx = ceiling(E_tot)

end subroutine computeFilterEnergy
!!***

end module m_slice_cprj
!!***

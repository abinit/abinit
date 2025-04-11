!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply 
!! the Spectrum Slicing method. It mainly defines 'slice' 
!! datatypes and associated methods. Features:
!! - uses xgTools implementation as matrix data structure.
!! - based on 'chebfi' data structure for most vector routines.
!! - polynomial filters are Chebyshev for first slice and
!!   Chebyshev-Jackson expansion of indicator otherwise.
!! - implements new parallel level between slices.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (IML, LB)
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

module m_slice

    use defs_basis
    use defs_abitypes
    use m_abicore
    use m_errors
    use m_time, only : timab
    use m_sort, only: sort_dp
    use m_io_tools, only : flush_unit

    use m_cgtools
    use m_xg
    use m_xgTransposer
    use m_xg_ortho_RR
    use m_chebfi2

    use m_xmpi
    use m_xomp
#ifdef HAVE_OPENMP
    use omp_lib
#endif

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    use m_gpu_toolbox, only : CPU_DEVICE_ID, gpu_device_synchronize
#endif

#if defined(HAVE_GPU) && defined(HAVE_GPU_MARKERS)
    use m_nvtx_data
#endif

    implicit none

    private

    ! Load balance criterion for slice distribution
    !-------------------------------------------------
    integer, parameter :: BALANCE_BANDPP      = 0 ! balanced bandpp        across processes
    integer, parameter :: BALANCE_BANDPP_NDEG = 1 ! balanced bandpp x ndeg across processes

    ! Private 'parallelEnv' datatype
    ! Parameters of slice parallelisation environment
    !------------------------------------------------
    type, private :: parallelEnv_t

        ! MPI-related
        integer :: rank
        integer :: size
        integer :: ierr
        integer :: comm
        integer :: comm_sub
        integer :: group
        integer :: status(MPI_STATUS_SIZE)

        ! Flags
        logical :: on_host = .false.
        logical :: on_device = .false.
        logical :: is_row = .false.
        logical :: is_col = .false.

        ! Arrays
        integer, allocatable :: ncolsColsRows(:)
        integer, allocatable :: nrowsLinalg(:)
        integer, allocatable :: group_lookup(:)

    end type parallelEnv_t

    ! Public 'slice' datatype
    ! [IML 9/4 work in progress] I don't think this is necessary
    ! Parameters specific to individual slices
    !-------------------------------------------------
    type, public :: sliceWorker_t

        integer :: islice                       ! slice index in 1,..,nslice
        integer :: nband                        ! number of bands in slice
        integer :: i1                           ! first band index in global
        integer :: i2                           ! last band index in global
        integer :: degree                       ! filter degree
        integer :: bandpp                       ! number of slice bands per process
        real(dp) :: ramp                        ! filter amplification factor
        real(dp) :: low                         ! filter support lower bound
        real(dp) :: upp                         ! filter support upper bound
        real(dp) :: glb                         ! spectrum lower bound
        real(dp) :: gub                         ! spectrum upper bound

        ! Memory allocated for a slice
        type(xg_t) :: XW                        ! input/ output eigenvectors
        type(xg_t) :: EW                        ! output eigenvalues
        type(xg_t) :: RW                        ! output residuals
        type(xg_t) :: OCCW                      ! unused memory space in Slicing
        type(chebfi_t) :: chebfi                ! workspace for slice diago
   
        ! Pointers to slice memory
        type(xgBlock_t) :: xgx0
        type(xgBlock_t) :: xgeigen
        type(xgBlock_t) :: xgresidu
        type(xgBlock_t) :: xgocc ! not used at all

    end type sliceWorker_t

    ! Public 'slice' datatype
    ! Parameters common to all slices
    !-------------------------------------------------
    type, public :: sliceAll_t

        integer :: nslice                        ! Total number of slices
        integer :: space
        integer :: space_res
        integer :: spacedim                      ! Space dimension for one vector
        integer :: total_spacedim                ! Maybe not needed
        integer :: neigenpairs                   ! Number of eigen values/vectors we want
        integer :: mdeg_filter                   ! Degree of the polynomial filter
        integer :: spacecom                      ! Communicator for MPI
        real(dp) :: tolerance                    ! Tolerance on the residu to stop the minimization
        real(dp) :: ecut                         ! Ecut used for Polynomial filtering
       
        ! Switches that allow to follow the state
        logical :: buffer_mem = .false.
        logical :: buffer_rows = .false.
        logical :: buffer_cols = .false.

        ! Variables deactivated for chebfi
        integer :: oracle = 0                        
        integer :: nbdbuf = 0                       
        real(dp) :: oracle_factor = 1.d0                
        real(dp) :: oracle_min_occ = 0.d0              

        integer :: paral_kgb                     ! allow MPI distribution of k-points, bands or PW
        integer :: paral_slice                   ! distribute slices over MPI processes
        integer :: bandpp                        ! nbands per MPI process
        integer :: comm_cols                     ! MPI column communicator
        integer :: comm_rows                     ! MPI row communicator
        integer :: me_g0           
        integer :: me_g0_fft

        ! General slicing params
        integer :: nband_ovlp                        ! number of columns in overlap-free space
        integer :: npband                            ! number of MPI processes
        integer :: spectral_cut                         ! spectral partition technique
        real(dp) :: glb                              ! global lower spectral bound
        real(dp) :: gub                              ! global upper spectral bound
        real(dp) :: ramp                             ! amplification factor for slices
        
        logical :: paw
        integer :: eigenProblem   ! 1 (A*x = (lambda)*B*x), 2 (A*B*x = (lambda)*x), 3 (B*A*x = (lambda)*x)

        ! when GPU is enabled, currently OpenMP is not fully supported, abinit is launched
        ! with OMP_NUM_THREADS=1, but we may locally increase the number of OpenMP threads
        ! wherever it is safe to do; in that case we use gpu_kokkos_nthrd to specify
        ! the number of OpenMP threads. This value is controlled by dtset variable
        ! dtset%gpu_kokkos_nthrd
        integer :: gpu_option
        integer :: gpu_kokkos_nthrd = 1 ! only used if gpu is enabled, number of OpenMP threads used
        integer :: gpu_thread_limit = 1 ! only used if GPU is enabled, max number of OpenMP threads used in sensitive areas

        ! DOS (all eigenvalues+residuals) before Slicing. Not distributed across MPI proc. 
        ! Useful for convergence test
        type(xg_t) :: Eig0                           ! Rayleigh quotients before Slicing
        type(xg_t) :: Res0                           ! residual norms before Slicing

        ! Initial memory space (in/out)
        type(xgBlock_t) :: X
        type(xgBlock_t) :: eigen
        type(xgBlock_t) :: residu

        ! Extended memory space
        type(xg_t) :: Xext_                        ! input/ output eigenvectors
        type(xg_t) :: Eext_                        ! output eigenvalues
        type(xg_t) :: Rext_                        ! output residuals

        ! Pointers to extended memory (linalg representation)
        type(xgBlock_t) :: Xext_linalg
        type(xgBlock_t) :: Eext_linalg
        type(xgBlock_t) :: Rext_linalg
        ! colsrows representation
        type(xgBlock_t) :: Xext

        ! Transposer
        type(xgBlock_t) :: Transposer_Xext

        ! Pointers to all slice parameters
        integer, pointer :: pband(:) => NULL()        ! Eigenvector column permutation
        integer, pointer :: idx(:,:) => NULL()        ! start and end indices per slice
        integer, pointer :: idx_merge(:,:) => NULL()  ! start and end converged indices per slice
        integer, pointer :: ndeg(:) => NULL()         ! filter degree per slice 
        integer, pointer :: npbandSlice(:) => NULL()
        real(dp), pointer :: sbound(:,:) => NULL()    ! (slice low bound,slice upp bound,
                                                      !  filter low support,filter upp support)

    end type slice_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: slice_init                        ! initiate slice data type object
    public :: slice_run                         ! diagonalize individual slices
    public :: slice_free                        ! free     slice data type object
    public :: slice_getSpectrum                 ! compute spectrum as Rayleigh quotients
    public :: slice_distributeSpectrum          ! distribute spectral slices to processes
    public :: slice_mergeConverged              ! merge converged slices by removing duplicates
    public :: slice_unitTest                    ! used for debugging

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!!
!! FUNCTION
!! Initialize a 'slice' datastructure.
!! See chebfi_init.
!!
!! SOURCE

subroutine slice_init(slice,neigenpairs,spacedim,tolerance,ecut,paral_kgb,paral_slice,&
        bandpp,mdeg_filter,space,eigenProblem,spacecom,me_g0,me_g0_fft,&
        paw,comm_rows,comm_cols,nslice,ramp,balance,gpu_option,&
        gpu_kokkos_nthrd,gpu_thread_limit)

    implicit none

    ! Arguments ------------------------------------
    integer      , intent(in   ) :: bandpp
    integer      , intent(in   ) :: npband
    integer      , intent(in   ) :: nslice
    integer      , intent(in   ) :: eigenProblem
    integer      , intent(in   ) :: me_g0
    integer      , intent(in   ) :: me_g0_fft
    integer      , intent(in   ) :: neigenpairs
    integer      , intent(in   ) :: mdeg_filter
    integer      , intent(in   ) :: comm_cols
    integer      , intent(in   ) :: comm_rows
    integer      , intent(in   ) :: paral_kgb
    integer      , intent(in   ) :: paral_slice
    integer      , intent(in   ) :: space
    integer      , intent(in   ) :: spacecom
    integer      , intent(in   ) :: spacedim
    integer      , intent(in   ) :: balance
    integer      , intent(in   ) :: gpu_option
    logical      , intent(in   ) :: paw
    real(dp)     , intent(in   ) :: ramp
    real(dp)     , intent(in   ) :: ecut
    real(dp)     , intent(in   ) :: tolerance
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ), optional :: gpu_kokkos_nthrd
    integer      , intent(in   ), optional :: gpu_thread_limit

    ! *********************************************************************

    slice%space        = space
    slice%neigenpairs  = neigenpairs
    slice%spacedim     = spacedim
    slice%tolerance    = tolerance
    slice%ecut         = ecut
    slice%paral_kgb    = paral_kgb
    slice%paral_slice  = paral_slice
    slice%comm_cols    = comm_cols
    slice%bandpp       = bandpp
    slice%comm_rows    = comm_rows
    slice%mdeg_filter  = mdeg_filter
    slice%spacecom     = spacecom
    slice%eigenProblem = eigenProblem
    slice%me_g0        = me_g0
    slice%me_g0_fft    = me_g0_fft
    slice%paw          = paw
    slice%gpu_option   = gpu_option
    slice%nslice       = nslice
    slice%npband       = xmpi_comm_size(comm_cols)
    slice%ramp         = ramp
    slice%spectral_cut = balance
    slice%nband_ovlp   = neigenpairs ! see initExtended
    slice%nbdbuf       = nbdbuf

    slice%gpu_kokkos_nthrd = 1
    if (present(gpu_kokkos_nthrd)) slice%gpu_kokkos_nthrd = gpu_kokkos_nthrd
    slice%gpu_thread_limit = 0
    if (present(gpu_thread_limit)) slice%gpu_thread_limit = gpu_thread_limit

    ! Space of eigenvalues
    if (space==SPACE_C) then
        slice%space_res = SPACE_C
    else if (space==SPACE_CR) then
        slice%space_res = SPACE_R
    end if

    call slice_allocateAll(slice)
    
end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateAll
!! NAME
!! slice_allocateAll
!! 
 
subroutine slice_allocateAll(slice)

    implicit none
    
    type(slice_t), intent(inout) :: slice
    integer :: space_res
    integer :: neigenpairs
    integer :: comm_cols
    integer :: gpu_option

    call slice_free(slice)

    space_res = slice%space_res
    neigenpairs = slice%neigenpairs
    comm_cols = slice%comm_cols
    !gpu_option = slice%gpu_option
    ! FIXME forced CPU
    gpu_option = ABI_GPU_DISABLED

    ! Eigenvalues and residuals before slicing
    ! Every MPI process contains this array
    ! IML 9/4 FIXME this is useful if we compute convergence rates aposteriori
    ! otherwise, I am not sure how to use this information
    call xg_init(slice%Eig0,space_res,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)
    call xg_init(slice%Res0,SPACE_R,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)

    ! FIXME incorporate? Add to mpiData? Add to slice?
    ! Memory allocations of size depending on fixed nslice
    ABI_MALLOC(pband, (nband)); pband_ptr => pband                     ! band permutation
    ABI_MALLOC(idx, (nslice,2)); idx_ptr => idx                        ! idx in overlapping mem
    ABI_MALLOC(idx_ovlp, (nslice,2)); idx_ovlp_ptr => idx_ovlp         ! idx in overlap-free mem
    ABI_MALLOC(idx_merge, (nslice,2)); idx_merge_ptr => idx_merge      ! converged idx in overlap-free mem
    ABI_MALLOC(ndeg, (nslice)); ndeg_ptr => ndeg                       ! filter degree per slice
    ABI_MALLOC(sbound, (nslice,4)); sbound_ptr => sbound               ! eigenvalue bounds per slice
    ABI_MALLOC(npbandSlice, (nslice)); npbandSlice_ptr => npbandSlice  ! number of mpi processes per slice

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! NAME
!! slice_free
!! 
 
subroutine slice_free(slice)

    implicit none
    
    type(slice_t), intent(inout) :: slice

    call slice_freeOverlapFree(slice)
    call xg_free(slice%Eig0)
    call xg_free(slice%Res0)
 
    ! Free slice parameters
    if (allocated(pband)) ABI_FREE(pband)
    if (allocated(idx)) ABI_FREE(idx)
    if (allocated(idx_ovlp)) ABI_FREE(idx_ovlp)
    if (allocated(idx_merge)) ABI_FREE(idx_merge)
    if (allocated(ndeg)) ABI_FREE(ndeg)
    if (allocated(sbound)) ABI_FREE(sbound)
    if (allocated(npbandSlice)) ABI_FREE(npbandSlice)
    
end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_computeDos
!! NAME
!! slice_computeDos
!! 
!! FUNCTION
!! Compute density of states (DOS) for all bands.
!! 
!! SOURCE

subroutine slice_dos(slice,X0,getAX_BX,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: X0
    integer, intent(in) :: nspinor
    interface
        subroutine getAX_BX(X,AX,BX)
            use m_xg, only : xgBlock_t
            type(xgBlock_t), intent(inout) :: X
            type(xgBlock_t), intent(inout) :: AX
            type(xgBlock_t), intent(inout) :: BX
        end subroutine getAX_BX
    end interface
    ! Local variables-------------------------------    
    integer :: neigenpairs,spacedim,paral_kgb,space_res
    integer :: bandpp,mdeg_filter,space,eigenProblem,spacecom
    integer :: me_g0,me_g0_fft,comm_rows,comm_cols
    integer :: gpu_option,gpu_kokkos_nthrd,gpu_thread_limit
    integer :: nrows
    integer :: nbdbuf,oracle
    integer :: my_rank,shift_row ! for MPI all to all
    logical :: paw
    real(dp) :: tolerance,ecut
    real(dp) :: oracle_factor,oracle_min_occ
    type(chebfi_t) :: chebfi
    type(xg_t) :: Eig0_XW, Res0_XW
    type(xgBlock_t) :: Eig0_paral, Res0_paral
 
    ! *********************************************************************

    ! TODO examine if this short alternative is sufficient
    !if (chebfi%paral_kgb == 0) then
    !    call xg_init(DivResults, space_res, neigenpairs, 1, gpu_option=chebfi%gpu_option)
    !else
    !   call xg_init(DivResults, space_res, bandpp, 1, gpu_option=chebfi%gpu_option)
    !end if
    !call chebfi_rayleighRitzQuotients(chebfi, maxeig, mineig, DivResults%self)

    space          = slice%space
    spacedim       = slice%spacedim
    neigenpairs    = slice%neigenpairs
    space_res      = slice%space_res
    bandpp         = slice%bandpp
    mdeg_filter    = slice%mdeg_filter
    tolerance      = slice%tolerance
    ecut           = slice%ecut
    paral_kgb      = slice%paral_kgb
    spacecom       = slice%spacecom
    me_g0          = slice%me_g0
    me_g0_fft      = slice%me_g0_fft
    eigenproblem   = slice%eigenproblem
    paw            = slice%paw
    comm_rows      = slice%comm_rows
    comm_cols      = slice%comm_cols
    nbdbuf         = slice%nbdbuf
    oracle         = slice%oracle
    oracle_factor  = slice%oracle_factor
    oracle_min_occ = slice%oracle_min_occ

    !gpu_option       = slice%gpu_option
    ! FIXME forced CPU, to be done in GPU
    gpu_option = ABI_GPU_DISABLED
    gpu_kokkos_nthrd = slice%gpu_kokkos_nthrd
    gpu_thread_limit = slice%gpu_thread_limit

    ! ==================== INITIALIZATION ===========================
    
    ! Initialize DOS workspace
    write(std_out,'(a,i0)') '-----> allocate chebfi for DOS, bandpp=', bandpp
    call chebfi_init(chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
&                    mdeg_filter,nbdbuf,space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,&
&                    comm_cols,oracle,oracle_factor,oracle_min_occ,gpu_option,&
&                    gpu_kokkos_nthrd=gpu_kokkos_nthrd,gpu_thread_limit=gpu_thread_limit)

    ! Number of vectors per MPI process
    nrows = neigenpairs
    if (paral_kgb == 1) nrows = bandpp

    ! Initialize result space eigenvalues and residuals before slicing, MPI distributed
    call xg_init(Eig0_XW,space_res,rows=nrows,cols=1,comm=comm_cols,gpu_option=gpu_option)
    call xg_init(Res0_XW,SPACE_R,rows=nrows,cols=1,comm=comm_cols,gpu_option=gpu_option)
    Eig0_paral = Eig0_XW%self
    Res0_paral = Res0_XW%self

    ! ======================== RUN ===============================

    ! Compute Rayleigh values and residuals of all bands (using MPI)
    call chebfi_RayleighValues(chebfi,X0,Eig0_paral,Res0_paral,getAX_BX,nspinor)
    
    ! Debug; problem is -7 is too low
    !write(std_out,*) 'Rayleigh values (before MPI row coms)'
    !call xgBlock_print(Eig0_paral, std_out)

    ! ===================== MPI ROW COMS =========================

    ! Store to slice object without distribution across MPI processes
    ! All processors have all eigenvalues and residuals
    ! FIXME remove if condition after partialcopy is debugged on GPU
    !if (gpu_option==ABI_GPU_OPENMP) then
    !    write(std_out,'(a)') 'Copy (Eig0,Res0) from gpu'
    !    call xgBlock_copy_from_gpu(Eig0_paral)
    !    call xgBlock_copy_from_gpu(Res0_paral)
    !end if
    if (xmpi_comm_size(comm_cols)>1) then

        ! Initialize entire array with zeros, every MPI contains this array
        call xgBlock_zero(slice%Eig0%self)
        call xgBlock_zero(slice%Res0%self)

        ! Recover rank of current MPI process
        my_rank = xmpi_comm_rank(comm_cols)
        shift_row = my_rank * bandpp

        ! Reshape in order to use blockCopy on columns (requires same number of rows)
        call xgBlock_reshape(Eig0_paral,(/1,bandpp/))     
        call xgBlock_reshape(Res0_paral,(/1,bandpp/))     
 
        ! Fill entire array with entry=(current MPI part | zero otherwise)
        ! FIXME replace blockCopy by a simple xgBlock_copy
        call slice_blockCopy(Eig0_paral,slice%Eig0%self,1,shift_row+1,bandpp,shift_row+bandpp)
        call slice_blockCopy(Res0_paral,slice%Res0%self,1,shift_row+1,bandpp,shift_row+bandpp)
        
        ! Sum entire object across MPI, every MPI contains the same entire array
        call xgBlock_mpi_sum(slice%Eig0%self,comm=comm_cols)
        call xgBlock_mpi_sum(slice%Res0%self,comm=comm_cols)

        ! Undo reshape
        call xgBlock_reshape(slice%Eig0%self,(/neigenpairs,1/))     
        call xgBlock_reshape(slice%Res0%self,(/neigenpairs,1/))

    else
        call xgBlock_copy(Eig0_paral,slice%Eig0%self)
        call xgBlock_copy(Res0_paral,slice%Res0%self)
    end if
    
    ! Debug; problem is -7 is too low
    !write(std_out,*) 'Rayleigh values (after MPI row coms)'
    !call xgBlock_print(slice%Eig0%self, std_out)

    ! Free workspace
    write(std_out,'(a)') 'free chebfi for DOS <-----'
    call chebfi_free(chebfi)
    call xg_free(Eig0_XW)
    call xg_free(Res0_XW)
 
end subroutine slice_dos
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_distributeSpectrum
!! NAME
!! slice_distributeSpectrum
!! 
!! FUNCTION
!! Split then distribute spectrum across processes.
!! Then create extended memory buffers and distribute them.
!! 
!! SOURCE

subroutine slice_distributeSpectrum(slice,spectral_cut,paral_slice)

    ! Compute nband per slice
    call slice_cutSpectrum(slice,pband_ptr,spectral_cut)
 
    ! Compute target mpi distribution
    ! in here assume that mpiData contains the rule_nband
    ! FIXME workinprogress
    call parallelEnv_init(env,slice%mpi_slice,paral_slice)

    ! Create the extended workspaces on CPU
    call slice_initExtended(slice,pband_ptr)

    ! Distribute extended space across MPI processes
    call slice_distributeExtended(slice)

end subroutine slice_distributeSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_cutSpectrum
!! NAME
!! slice_cutSpectrum
!! 
!! FUNCTION
!! Split spectrum into overlapping slices. 
!!
!! INPUT
!! sliceAll=         parameters common to all slices, such as
!!                   balance_option=balance number of vectors(1)
!!                                  balance filter degrees(2)
!!                   ramp=filter convergence threshold, greater than 1
!! idxAll=           start and end index per slice
!! ndegAll=          filter degree per slice
!! sboundAll=        various spectral bounds
!! pband=            eigenvector permutation
!!
!! SOURCE

subroutine slice_cutSpectrum(slice,nband,idxAll,ndegAll,sboundAll,pband,npbandSlice,plot_filter)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nband
    integer, pointer, intent(inout) :: idxAll(:,:)
    integer, pointer, intent(inout) :: ndegAll(:)
    integer, pointer, intent(inout) :: pband(:)
    integer, pointer, intent(inout) :: npbandSlice(:)
    real(dp), pointer, intent(inout) :: sboundAll(:,:)
    logical, optional, intent(in) :: plot_filter
    
    !Local variables-------------------------------
    integer :: j,k,jmax,spos,nvec_ovlp,nvec,k1,k2
    integer :: nv_pad,ndeg,npband,comm_cols
    integer :: paral_slice
    integer :: nslice,neigenpairs,nline
    integer :: ndeg_max = 200
    integer :: balance_option
    integer :: ipt,npt,iptL,iptR
    logical :: plot_filter_
    real(dp) :: ramp
    real(dp) :: tol12 = 1.0e-12
    real(dp) :: ecut,low,upp,glb,gub,c,r
    real(dp) :: lj,uj,wj,finL,finR,foutL,foutR
    real(dp) :: wlj,wuj
    real(dp) :: fun_pt,pt
    real(dp) :: a_,b_ ! target interval scaled in -1,1
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    integer :: jperm(nband-1)
    real(dp) :: consdiff(nband-1)
    real(dp), allocatable :: slice_cut(:)
    real(dp), allocatable :: resid_cut(:)
    real(dp), pointer :: resid_(:,:)
    real(dp), pointer :: theta_(:,:)
    real(dp), allocatable, target :: resid(:)
    real(dp), allocatable, target :: theta(:)

    ! *********************************************************************

    plot_filter_ = .false.
    if (present(plot_filter)) plot_filter_ = plot_filter

    ! Interval slicing parameters
    npband         = slice%npband ! number of MPI processes
    bandpp         = slice%bandpp ! fixed number of bands per MPI process
    comm_cols      = slice%comm_cols
    neigenpairs    = slice%neigenpairs
    nslice         = slice%nslice
    ecut           = slice%ecut
    nline          = slice%mdeg_filter
    paral_slice    = slice%paral_slice
    balance_option = slice%spectral_cut
    ramp           = slice%ramp

    ! Set pointers to eigenvalue and residual memory
    Eig0_all = slice%Eig0%self
    Res0_all = slice%Res0%self

    ! Results could be complex, so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(Eig0_all,theta_,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(Res0_all,resid_,rows=1,cols=neigenpairs)

    ! Save theta_,resid_ first row in theta,resid
    ABI_MALLOC(theta,(neigenpairs))
    ABI_MALLOC(resid,(neigenpairs))
    theta(1:neigenpairs) = theta_(1,1:neigenpairs)
    resid(1:neigenpairs) = resid_(1,1:neigenpairs)

    ! Sort thetas in increasing order and store permutation
    call sort_dp(neigenpairs,theta,pband,tol12)
    
    !write(std_out,*) 'Rayleigh Values after sort:'
    !write(std_out,*) theta(:)

    ! Spectrum bounds
    low = theta(1)
    upp = theta(neigenpairs)
    ! Bounds that contain entire spectrum (guaranteed)
    glb = low - sqrt(resid(pband(1)))
    gub = ecut 
    ! ecut too large resulting in fine slices, narrow down:
    !gub = upp + sqrt(resid(pband(neigenpairs))) + 2.d0
    ! this results in eigenvalues between gub and ecut

    ! Center and radius of the entire spectrum mapped to -1,1
    c = (glb + gub) / 2.d0
    r = (gub - glb) / 2.d0
    
    write(std_out,*) ' '
    write(std_out,*) '------------ A priori spectrum '
    write(std_out,*) 'Spectral bounds=', glb,gub
    write(std_out,*) ' Minmax eigvals=', low,upp
    write(std_out,*) '           ecut=', ecut
    write(std_out,*) '         center=', c
    write(std_out,*) '         radius=', r
    write(std_out,*) ' '

    ! Store to slicing object
    slice%gub = gub
    slice%glb = glb
   
    ! Initialize spectral cuts
    ABI_MALLOC(slice_cut,(nslice+1))
    slice_cut(:) = 0.d0
    slice_cut(1) = low
    slice_cut(nslice+1) = upp
    
    ! Compute spectral cuts on interior slices
    select case(balance_option)
    case(1)
        ! Balance interval widths
        slice_cut(2:nslice) = (/ (low+(upp-low)/nslice*j, j=1,nslice-1) /)
    case(2)
        ! Balance number of vectors
        slice_cut(2:nslice) = (/ (theta(neigenpairs/nslice*j), j=1,nslice-1) /)
    case(3)
        ! Cut on spectral gaps=where eigenvalues are less concentrated
        ! Sort consecutive differences by increasing order
        jperm = (/ (k, k=1,nband-1) /)
        consdiff = (/ (theta(k+1) - theta(k), k=1,nband-1) /)
        call sort_dp(nband-1,consdiff,jperm,tol12)
        ! Take median of largest gaps
        do i=1,nslice
            jmax = jperm(nband-i)
            slice_cut(i+1) = (theta(jmax) + theta(jmax+1)) / 2.d0
        end do
    end select

    ! FIXME function starts becoming too big.
    ! Split into two parts 
    ! * slice_cutSpectrum. * slice_setFilters

    ! Define spectral subintervals and optimize degrees for individual slices
    j1 = 1 ! band index ! declare and change
    do j=1,nslice
            
        ! Slice position, first (1), interior (2), last (3)
        spos = 2
        if (j==1) spos = 1
        if (j==nslice) spos = 3

        ! Spectral slice of interest is [l,u)
        lj = slice_cut(j)
        uj = slice_cut(j+1)
        
        ! Overlap width *****************************************************
        ! TODO investigate if we also need to add residual
        !                                                    (IML 28/01/2025)
        ! for (u-l)/10: faster because less vectors, but ratio > ramp sometimes
        !               SCF cycle errors are quite large (deltaE,res2)
        ! for (u-l)/8 : slower because more vectors, but ratio < ramp always
        !               also improves all three errors in SCF cycles
        ! Conclusion: there is a critical w above which convergence ramp is achieved
        ! Idea: use residual intervals to find critical w

        ! This is critical w fixed empirical
        ! Heuristic rule: if nvec and nvec_ovlp is very close, update wj
        !wj = (uj - lj) / 8.d0
        wj = (uj - lj) / 10.d0
        ! ndeg is ok (30)

        write(std_out,*) '------------ /Divide/ Slice ',j
        
        ! Uncomment following two lines to use wj fixed overlap width
        wlj = wj
        wuj = wj
        write(std_out,*) '       overlap widths=', wlj, wuj
        ! *******************************************************************
        
        ndeg = nline
       
        ! Count eigenvalues in slice cut plus overlap 
        call count_values(lj-wlj,uj+wuj,theta,spos,neigenpairs,k1,k2,nvec)
        nvec_ovlp = nvec
        write(std_out,*) '    after overlap', k1,k2

        ! How to compute optimal degree in slice with overlap:
        ! ********* control amplification ratio ********
        ! The convergence ratio r0/rN depends on the amplification
        ! ratios f(l)/f(l-w) and f(u)/f(u+w). Increasing w should
        ! improve the convergence ratio.

        ! Filter support [l-w,u+w) scaled to [-1,1)
        a_ = (lj-wlj-c)/r
        b_ = (uj+wuj-c)/r

        if (j==1) then
            ! uj,gub is the interval mapped to -1,1
            ! in this interval Chebyshev poly is bounded by 1
            ! remember uj,gub is the interval to ignore
            ndeg = 4
            do while
            !write(std_out,*) 'first slice apriori=', 1.d0/cheb_poly1(lj,12,uj+wuj,ecut)
        else
            ndeg = 4
            finL = 0.d0; finR = 0.d0; foutL = 1.d0; foutR = 1.d0
            do while ( (finL/foutL < ramp) .and. (finR/foutR < ramp) .and. (ndeg<ndeg_max) )
                finL  = bandpassIndicator_sca((lj-c)/r,a_,b_,ndeg)
                foutL = bandpassIndicator_sca(a_      ,a_,b_,ndeg)
                finR  = bandpassIndicator_sca((uj-c)/r,a_,b_,ndeg)
                foutR = bandpassIndicator_sca(b_      ,a_,b_,ndeg)
                ndeg = ndeg + 1
            end do
        end if

        ! Plot filter
        if (plot_filter_) then
            write(std_out,*) ' '
            write(std_out,*) 'Plot filter ==== x | f(x)'
            npt = 100
            if (j==1) then
                do ipt=1,npt
                    pt = lj + (ipt-1)*(uj-lj)/npt
                    fun_pt = cheb_poly1(pt,ndeg,uj+wuj,gub)
                    write(std_out,*) pt, fun_pt
                end do
            else
                do ipt=1,npt
                    pt = (lj + (ipt-1)*(uj-lj)/npt - c)/r
                    fun_pt = bandpassIndicator_sca(pt,a_,b_,ndeg)
                    write(std_out,*) pt, fun_pt
                end do
            end if
            write(std_out,*) ' '
        end if

        ! Print slice interval info
        write(std_out,*) 'Without overlap=', lj, uj
        write(std_out,*) '          width=', uj-lj
        write(std_out,*) '      scaled to=', (lj-c)/r,(uj-c)/r
        write(std_out,*) 'With    overlap=', lj-wlj, uj+wuj
        write(std_out,*) '          width=', uj+wuj-lj+wlj
        write(std_out,*) '      scaled to=', (lj-wlj-c)/r,(uj+wuj-c)/r
        write(std_out,*) 'With    overlap=', lj-wlj, uj+wuj
        write(std_out,*) '  nvec(balance)=', neigenpairs/nslice
        write(std_out,*) '      nvec_ovlp=', nvec_ovlp
        write(std_out,*) '           nvec=', nvec
        write(std_out,*) '           ndeg=', ndeg
        write(std_out,*) ' '
        ! end print

        ! Compute last index in extended memory (without ovlp)
        j2 = j1 + nvec_ovlp - 1

        ! Store slice parameters
        idxAll(j,1:2) = (/k1,k2/)
        idx_extAll(j,1:2) = (/j1,j2/)
        ! FIXME add nband (nband per slice with overlap)
        nbandAll(j) = k2 - k1 + 1
        ndegAll(j) = ndeg
        sboundAll(j,1:4) = (/lj,uj,lj-wlj,uj+wuj/) 

        slice%fcol(islice) = k1 
        slice%lcol(islice) = k2
        
        slice%fcol_ext(islice) = j1 ! FIXME declare this 
        slice%lcol_ext(islice) = j2 ! FIXME declare this 
       
        slice%nband(islice) = nvec_ovlp
        
        ! Update first index in next extended
        j1 = j2 + 1
           
    end do

    ! Update pointers for all slices
    ! TODO rename idx to icol
    slice%pband => pband
    slice%ndeg => ndegAll
    slice%sbound => sboundAll
    !slice%npbandSlice => npbandSlice

    mpi_slice%rule_nband => nbandAll
    mpi_slice%rule_ndeg => ndegAll

    ! Free workspace not needed
    if (allocated(slice_cut)) ABI_FREE(slice_cut)
    if (allocated(resid_cut)) ABI_FREE(resid_cut)
    if (allocated(theta)) ABI_FREE(theta)
    if (allocated(resid)) ABI_FREE(resid)
   
end subroutine slice_cutSpectrum
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_mergeConverged
!! NAME
!! slice_mergeConverged
!! 
!! FUNCTION
!! Find converged eigenvalues in each slice using a criterion.
!! [^Not true. For the moment we brutally merge using interval limits]
!! Return range of first and last index to merge per slice, both 
!! computed from a residual criterion on eigenvalues.
!! slice_mergeConverged is done locally on slice
!! while on linalg representation
!! so it is after Rayleigh-Ritz and before the last transposition
!! 
!! SOURCE

subroutine slice_mergeConverged(sliceAll,idx_ovlp,idx_merge,merge_option)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t) , intent(inout) :: sliceAll
    integer          , intent(in   ) :: merge_option
    integer, pointer , intent(inout) :: idx_merge(:,:)
    integer, pointer , intent(in   ) :: idx_ovlp(:,:)

    ! Local variables-------------------------------    
    integer :: neigenpairs,nslice,islice,i1,i2,j1,j2
    integer :: iband,nband_slice,nband_ovlp
    integer :: k,k1,k2,nvec_count,spos
    integer :: nband_conv,ndeg,k_out
    real(dp) :: ramp,slow,supp,flow,fupp
    real(dp) :: glb,gub
    real(dp) :: c,r ! center, radius
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    integer, pointer :: pband(:)
    real(dp), pointer :: theta0(:,:), resid0(:,:)
    real(dp), pointer :: thetaN(:,:), residN(:,:)
    !real(dp), pointer :: theta0_(:), resid0_(:)
    !real(dp), pointer :: thetaN_(:), residN_(:)
    real(dp), allocatable :: resid0_slice(:)
    real(dp), allocatable :: residN_slice(:)
    real(dp), allocatable :: theta_slice(:)
 
    ! *********************************************************************

    ! Various variables
    nslice = sliceAll%nslice
    neigenpairs = sliceAll%neigenpairs
    nband_ovlp = sliceAll%nband_ovlp
    glb = sliceAll%glb
    gub = sliceAll%gub
    ramp = sliceAll%ramp
    r = (gub - glb)/2.d0
    c = (gub + glb)/2.d0
    
    ! Memory pointers
    Eig0_all = sliceAll%Eig0%self
    Res0_all = sliceAll%Res0%self
    pband => sliceAll%pband ! sorting theta0 in increasing order 
    
    ! Results could be complex (with null imaginary part), so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(Eig0_all,theta0,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(Res0_all,resid0,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(sliceAll%xgeigen_ovlp,thetaN,rows=1,cols=nband_ovlp)
    call xgBlock_reverseMap(sliceAll%xgresidu_ovlp,residN,rows=1,cols=nband_ovlp)

    ! Use index maps to recover bands associated to a slice
    nband_conv = 0
    do islice=1,nslice
        slow = sliceAll%sbound(islice,1) ! slice low
        supp = sliceAll%sbound(islice,2) ! slice upp
        flow = sliceAll%sbound(islice,3) ! filter low
        fupp = sliceAll%sbound(islice,4) ! filter upp
        ndeg = sliceAll%ndeg(islice)

        ! Indices in overlapping objects theta0,resid0
        i1 = sliceAll%idx(islice,1)
        i2 = sliceAll%idx(islice,2)
        nband_slice = i2 - i1 + 1

        ! Indices in overlap-free objects thetaN,residN
        j1 = idx_ovlp(islice,1)
        j2 = idx_ovlp(islice,2)

        ABI_MALLOC(theta_slice,(nband_slice))
        
        ABI_MALLOC(residN_slice,(nband_slice))
        ABI_MALLOC(resid0_slice,(nband_slice))

        ! Assumes row distribution (process has all nbands_slice in memory)
        theta_slice(1:nband_slice) = thetaN(1, j1:j2) ! after slicing
        residN_slice(1:nband_slice) = sqrt(residN(1,j1:j2))
        resid0_slice(1:nband_slice) = sqrt(resid0(1,pband(i1:i2)))

        ! Slice position, first (1), interior (2), last (3)
        spos = 2
        if (islice==1) spos = 1
        if (islice==nslice) spos = 3

        write(std_out,*) ''
        write(std_out,*) '-------------- /Mark/ Slice', islice
        write(std_out,*) '      nband_slice=', nband_slice
        write(std_out,*) ''

        ! Number of eigenvalues in slice interval
        ! This is the ones we keep
        call count_values(slow,supp,theta_slice,spos,nband_slice,k1,k2,nvec_count)
        write(std_out,*) 'Partition:', slow, supp
        write(std_out,*) 'number of eigenvalues in Partition =', nvec_count
        write(std_out,*) 'max residual rN       in Partition =', maxval(sqrt(residN_slice(k1:k2)))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Partition (kept) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Partition (kept) =', maxval(theta_slice(k1:k2))
 
        ! Store limited indices in overlap-free memory, attention shift k1,k2 by j1
        if (nband_conv+k2-k1+1>neigenpairs) then
            ! excess of eigenvalues, ignore last ones
            write(std_out,'(a,i0)') 'last one is ', k2
            k2 = k2 - (nband_conv+k2-k1+1-neigenpairs)
            write(std_out,'(a,i0)') 'due to excess, shift last one to ', k2
            idx_merge(islice,1) = j1 + k1 - 1
            idx_merge(islice,2) = j1 + k2 - 1
            nband_conv = nband_conv + k2 - k1 + 1
        else
            idx_merge(islice,1) = j1 + k1 - 1
            idx_merge(islice,2) = j1 + k2 - 1
            nband_conv = nband_conv + k2 - k1 + 1
        end if 
        write(std_out,*) ''

        ! Number of eigenvalues in filter support
        call count_values(flow,fupp,theta_slice,spos,nband_slice,k1,k2,nvec_count)
        write(std_out,*) 'Support:', flow, fupp
        write(std_out,*) 'number of eigenvalues in Support      =', nvec_count
        write(std_out,*) 'max residual rN       in Support      =', maxval(sqrt(residN_slice(k1:k2)))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Support (converged) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Support (converged) =', maxval(theta_slice(k1:k2))

        if (allocated(theta_slice)) ABI_FREE(theta_slice)
        if (allocated(residN_slice)) ABI_FREE(residN_slice)
        if (allocated(resid0_slice)) ABI_FREE(resid0_slice)

    end do
    
    ! Report missing or extra eigenvalues
    if (nband_conv<neigenpairs) then
        ABI_ERROR("Not enough converged eigenvalues")
    else if (nband_conv>neigenpairs) then
        ABI_ERROR("Too many converged eigenvalues")
    end if

    ! TODO 
    ! * count how many thetaN_ are in low,upp for every slice
    ! * count how many thetaN_ are outside current slice, and if they converged
    ! * count how many are in overlap region
    ! This will help diagnostic convergence "slice full"

    ! FIXME actually do the copy from extended to io
    ! replace blockCopy by xgBlock_copy
    call xgBlock_reshape(xgeigen, (/1,nband/))
    call xgBlock_reshape(xgresidu, (/1,nband/))
    j1 = 1                      ! start copy to overlapping mem (cg,eig,resid)
    do islice=1,nslice
        i1 = idx_merge(islice,1) ! start read from overlap-free mem
        i2 = idx_merge(islice,2)
        nband_merge = i2 - i1 + 1
        j2 = j1 + nband_merge - 1
        call slice_blockCopy(sliceAll%xgx0_ovlp,xgx0,i1,j1,i2,j2)
        call slice_blockCopy(sliceAll%xgeigen_ovlp,xgeigen,i1,j1,i2,j2)
        call slice_blockCopy(sliceAll%xgresidu_ovlp,xgresidu,i1,j1,i2,j2)
        j1 = j2 + 1
    end do
    ! Write clean as this:
    !do islice=1,nslice
    !    ncols = spsl%nband_slice(islice)
    !    fcol = spsl%fcol_slice_merge(islice)
    !    fcol_buf = spsl%fcol_buf_merge(islice)
    !    call xgBlock_setBlock(X0,xgcols_out,spacedim,ncols,fcol=fcol)
    !    call xgBlock_setBlock(spsl%Bufr,xgcols_in,spacedim,ncols,fcol=fcol_buf)
    !    call xgBlock_copy(xgcols_in,xgcols_out)
    !end do
    call xgBlock_reshape(xgeigen, (/nband,1/))
    call xgBlock_reshape(xgresidu, (/nband,1/))
 
end subroutine slice_mergeConverged
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_init
!! NAME
!! slice_init
!! 
!! SOURCE

subroutine slice_init(sliceAll,slice,islice)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t), intent(in   ) :: sliceAll
    type(slice_t   ), intent(inout) :: slice
    integer         , intent(in   ) :: islice

    ! Local variables-------------------------------    
    integer :: bandpp_slice,i1,i2,npband_slice,nband_slice
    integer :: tim_slice_init
    real(dp) :: tsec(2)

    ! *********************************************************************

    tim_slice_init = tim_slice2_init + islice - 1
    call timab(tim_slice_init,1,tsec)
    
    ! Input variables
    slice%islice = islice
    
    ! Variables inherited from sliceAll
    slice%glb = sliceAll%glb
    slice%gub = sliceAll%gub
    slice%low = sliceAll%sbound(islice,3)
    slice%upp = sliceAll%sbound(islice,4)
    slice%degree = sliceAll%ndeg(islice)
    npband_slice = sliceAll%npbandSlice(islice)
    i1 = sliceAll%idx(islice,1)
    i2 = sliceAll%idx(islice,2)

    ! Number of bands in slice
    nband_slice = i2 - i1 + 1
    slice%nband = nband_slice

    ! Number of bands per MPI process
    bandpp_slice = nband_slice
    if (sliceAll%paral_kgb==1) bandpp_slice=nband_slice/npband_slice
    slice%bandpp = bandpp_slice

    call slice_allocateAll(sliceAll,slice)

    ! Make dimension compatible with xgeigen,xgresidu of slicewf
    call xgBlock_reshape(slice%EW%self, (/nband_slice,1/))
    call xgBlock_reshape(slice%RW%self, (/nband_slice,1/))

    ! Attention T_n<0 for n odd
    ! TODO add a warning for fist slice

    ! Set pointers to slice memory
    slice%xgx0 = slice%XW%self
    slice%xgeigen = slice%EW%self
    slice%xgresidu = slice%RW%self
    slice%xgocc = slice%OCCW%self
    
    call timab(tim_slice_init,2,tsec)

end subroutine slice_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_allocateAll
!! NAME
!! slice_allocateAll
!! 
!! FUNCTION
!! Allocate slice memory spaces on CPU/GPU
!! 
!! SOURCE

subroutine slice_allocateAll(sliceAll,slice)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t), intent(in   ) :: sliceAll
    type(slice_t)   , intent(inout) :: slice

    ! Local variables-------------------------------    
    integer :: spacedim,paral_kgb,nband,bandpp
    integer :: ndeg,space,eigenProblem,spacecom
    integer :: me_g0,me_g0_fft,comm_rows,comm_cols
    integer :: gpu_option,gpu_kokkos_nthrd,gpu_thread_limit
    integer :: nbdbuf, oracle
    integer  :: total_spacedim, ierr
    logical :: paw
    real(dp) :: tolerance,ecut
    real(dp) :: oracle_factor, oracle_min_occ

    ! *********************************************************************
    
    ! Various parameters common to all slices
    space            = sliceAll%space
    spacedim         = sliceAll%spacedim
    spacecom         = sliceAll%spacecom
    me_g0            = sliceAll%me_g0
    me_g0_fft        = sliceAll%me_g0_fft
    tolerance        = sliceAll%tolerance
    ecut             = sliceAll%ecut
    paral_kgb        = sliceAll%paral_kgb
    comm_rows        = sliceAll%comm_rows
    comm_cols        = sliceAll%comm_cols
    paw              = sliceAll%paw
    gpu_option       = sliceAll%gpu_option
    gpu_kokkos_nthrd = sliceAll%gpu_kokkos_nthrd
    gpu_thread_limit = sliceAll%gpu_thread_limit

    ! not used but passed to chebfi with deactivated values
    oracle           = 0
    nbdbuf           = 0
    oracle_factor    = 1.d0
    oracle_min_occ   = 0.d0

    ! Various parameters of current slice
    nband  = slice%nband
    bandpp = slice%bandpp
    ndeg   = slice%degree
  
    call slice_free(slice)

    ! With current def each xg is distributed along spacecom=plane-wave MPI distr
    ! TODO Every slice MPI proc has entire space
    ! FIXME EW,RW devrait être alloués exactement comme xgx0 et xgresidu in chebfiwf 
    !! FIXME this is not necessary because extended spaces does the pointer
    write(std_out,'(a,i0)') '-----> Allocating slice result memory nband_slice=', nband
    call xg_init(slice%XW,space,spacedim,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
        
    ! TODO each column MPI has all rows and bandpp columns
    ! The subcomm_band should be defined outside this function. So in 
    ! slice_init 
    call xg_init(slice%XW,space,total_spacedim,bandpp,subcomm_band,&
&                me_g0=me_g0,gpu_option=gpu_option)
    
    call xg_init(slice%EW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%RW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%OCCW,SPACE_R,nband,1,gpu_option=gpu_option)! fill with zero

    ! Get total number of rows that are distributed along MPI processes
    total_spacedim = spacedim
    if (paral_kgb == 1) then
        call xmpi_sum(total_spacedim,spacecom,ierr)
    end if
    ! TODO paral_slice values should be encoded
    if (paral_slice == 0) then
        ! the communicator of chebfi is the large one
    else if (paral_slice == 1) then
        ! the communicator of chebfi is the subcomm one
    end if

    ! Every slice MPI proc has a part (bandpp) of this space
    write(std_out,'(a,i0)') '-----> Allocating slice working memory bandpp=', bandpp
    call chebfi_init(slice%chebfi,nband,spacedim,tolerance,ecut,paral_kgb,bandpp,&
&                    ndeg,nbdbuf,space,1,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
&                    oracle,oracle_factor,oracle_min_occ,gpu_option,&
&                    gpu_kokkos_nthrd=gpu_kokkos_nthrd,gpu_thread_limit=gpu_thread_limit)

end subroutine slice_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_free
!! NAME
!! slice_free
!! 
!! SOURCE

subroutine slice_free(slice)

    implicit none

    type(slice_t), intent(inout) :: slice
    integer :: tim_slice_free
    real(dp) :: tsec(2)

    tim_slice_free = tim_slice2_free + slice%islice - 1
    call timab(tim_slice_free,1,tsec)
    
    call xg_free(slice%XW)
    call xg_free(slice%EW)
    call xg_free(slice%RW)
    call xg_free(slice%OCCW)
    call chebfi_free(slice%chebfi)
    
    call timab(tim_slice_free,2,tsec)

end subroutine slice_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/bandpass_sca
!! NAME
!! bandpass_sca
!!
!! FUNCTION
!! Computes delta-Dirac polynomial filter f(x) approximated by a Chebyshev
!! expansion of order deg centered at gamma, evaluated at point x=t
!!
!! INPUTS
!!  t=      scalar to evaluate filter on
!!  deg=    order of Chebyshev expansion
!!  gam=    center of Chebyshev expansion
!!
!! OUTPUT
!!  res
!!
!! SOURCE

function bandpass_sca(t, deg, gam) result(f_t)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: t, gam
    integer , intent(in ) :: deg

    real(dp) :: f_t
    
    !Local variables-------------------------------
    real(dp) :: yt0, yt, yg0, yg, yt_swap, yg_swap
    real(dp) :: mu, damp, rho, rhog, theta
    integer  :: i
    
! *********************************************************************

    ! init cheby of deg=0,1 eval at t,gamma
    yt0 = 1.d0
    yt = t
    
    yg0 = 1.d0
    yg = gam

    ! init delta-Dirac filters of deg=0
    theta = Pi/(deg + 1)
    damp = SIN(theta) / theta
    rho = 0.5d0 + gam * damp * yt
    rhog = 0.5d0 + gam * damp * yg

    do i=2,deg 

        ! Update Chebyshev polynomials
        yt_swap = yt
        yt = 2 * t * yt - yt0
        yt0 = yt_swap
        
        yg_swap = yg
        yg = 2 * gam * yg - yg0
        yg0 = yg_swap

        ! Update delta-Dirac filters
        mu = COS(i * ACOS(gam))
        damp = SIN(i * theta) / (i * theta)
        rho = rho + mu * damp * yt
        rhog = rhog + mu * damp * yg
        
    end do

    f_t = rho / rhog

end function bandpass_sca
!!***

!----------------------------------------------------------------------

!!****f* m_slice/bandpassIndicator_sca
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

!!****f* m_slice/count_values
!! NAME
!! count_values
!! 
!! FUNCTION
!! Return number of consecutive theta values in interval [a,b)
!! and their index range (first and last indices).
!! Slice position (spos) takes into account extremal indices.
!! 
!! SOURCE

subroutine count_values(low,upp,theta,spos,nband,i1,i2,nvec)

    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: spos,nband
    integer, intent(inout) :: i1,i2,nvec
    real(dp), intent(in) :: low,upp
    real(dp), intent(in) :: theta(1:nband)
    !Local variables-------------------------------

! *********************************************************************

    i1 = maxloc(theta, dim=1, mask=(theta < low)) + 1
    i2 = maxloc(theta, dim=1, mask=(theta < upp))            
    if (spos == 1) i1 = 1
    if (spos == 3) i2 = nband
    nvec = i2 - i1 + 1
    nvec = i2 - i1 + 1

end subroutine count_values
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_blockCopy
!! NAME
!! slice_blockCopy
!! 
!! FUNCTION
!! Read from A deep copy to B, given column ranges.
!! Handles GPU / CPU data location by calling xgBlock_copy as well as
!! enter data map OpenMP ranges for execution on CPU or GPU.
!! 
!! SOURCE

subroutine slice_blockCopy(A_in,B_out,a1,b1,a2,b2)

    implicit none

    !Arguments ------------------------------------
    type(xgBlock_t)  , intent(in   ) :: A_in
    type(xgBlock_t)  , intent(inout) :: B_out
    integer          , intent(in   ) :: a1,b1 ! first columns
    integer, optional, intent(in   ) :: a2,b2 ! last columns
    
    !Local variables-------------------------------
    integer :: nrowsA,ncolsA
    integer :: nrowsB,ncolsB
    integer :: a2_,b2_
    logical :: on_host 
    integer :: gpua,gpub ! gpu_options
    integer :: ncolsA_block, ncolsB_block
    type(xgBlock_t) :: A_block, B_block
    real(dp) :: tsec(2)

! *********************************************************************

    call timab(tim_slice_Acopy,1,tsec)
    
    call xgBlock_getSize(A_in ,nrowsA,ncolsA)
    call xgBlock_getSize(B_out,nrowsB,ncolsB)

    if (nrowsA/=nrowsB) ABI_ERROR("A and B should have the same number of rows")

    ! DEBUG print gpu_option of in/out objects
    ! OMP query: are we on CPU or not?
    call xgBlock_get_gpu_option(A_in ,gpua)
    call xgBlock_get_gpu_option(B_out,gpub)
    write(std_out,'(a,i0,a,i0,a,i0)') 'Memcopy: gpu_optionA ',gpua,' gpu_optionB ',gpub
    if (gpua==ABI_GPU_OPENMP .or. gpub==ABI_GPU_OPENMP) then 
        on_host = xomp_is_initial_device()
        write(std_out,*) 'Memcopy: on_host', on_host
    end if
    !! IML 14/03 TODO decide whether we perform actions of copy from to gpu in here ..
    ! first detect cases of incompatiblity
    ! print on_host has no point because only the CPU prints anyway .. always true
    ! ..

    ! Get last index in block
    a2_ = ncolsA
    b2_ = ncolsB
    if (present(a2)) a2_ = a2
    if (present(b2)) b2_ = b2
    
    ! Number of block columns in range
    ncolsA_block = a2 - a1 + 1
    ncolsB_block = b2 - b1 + 1

    ! Deep copy from X to Y
    ! ============================
    ! xgBlock_copy will do FIXME also depends on on_host!! 
    ! if A CPU and B GPU then: copy from B GPU to B CPU, copy from A CPU to B CPU
    ! so the copy is performed on CPU if one of A and B is on CPU.
    call xgBlock_setBlock(A_in ,A_block,rows=nrowsA,cols=ncolsA_block,fcol=a1)
    call xgBlock_setBlock(B_out,B_block,rows=nrowsB,cols=ncolsB_block,fcol=b1)
    call xgBlock_copy(A_block,B_block)
    
    call timab(tim_slice_Acopy,2,tsec)

end subroutine slice_blockCopy
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_initExtended
!! NAME
!! slice_initExtended
!! 
!! FUNCTION
!! Allocate extended buffer by replicating overlapping slice data.
!! This creates a buffer without data overlap on which we can safely
!! read and write data avoiding concurrent memory access. The 
!! overlapping data has two independent copies for adjacent slices.
!! The extended memory space is bigger than original one.
!! 
!! SOURCE

subroutine slice_initExtended(slice,pband_ptr)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, pointer, intent(in) :: pband_ptr
    
    ! Local variables-------------------------------    
    integer :: me_g0,nband,comm
    integer :: nrow,ncol
    integer :: islice,nslice
    integer :: fcol,fcol_ext
    type(xgBlock_t) :: xgcols_in, xgcols_out 
    
    ! *********************************************************************

    nslice = slice%nslice
    nrow = slice%spacedim
    comm = slice%spacecom
    nband = slice%nband
    me_g0 = slice%me_g0
    ncol = sum(slice%nband_slice)

    ! Permute columns of X in linalg representation
    call xgBlock_permuteCols(slice%X,nrow,nband,pband_ptr)

    if (nslice==1) then
        slice%Xext_linalg = slice%X
        slice%Eext_linalg = slice%eigen
        slice%Rext_linalg = slice%residu
    else
        ! Allocate extended spaces in linalg representation (on CPU)
        call xg_init(slice%Xext_,space,nrow,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
        call xg_init(slice%Eext_,SPACE_R,1,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
        call xg_init(slice%Rext_,SPACE_R,1,ncol,comm,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
    
        ! Define pointers
        slice%Xext_linalg = slice%Xext_%self
        slice%Eext_linalg = slice%Eext_%self
        slice%Rext_linalg = slice%Rext_%self

        ! Copy X to extended by blocks
        ! Reminder: xgBlock_copy is always on CPU expect if both blocks are on GPU
        do islice=1,nslice
            ncol = slice%nband_slice(islice)
            fcol = slice%fcol(islice)
            fcol_ext = slice%fcol_ext(islice)
            call xgBlock_setBlock(slice%X,xgcols_in,nrow,ncol,fcol=fcol)
            call xgBlock_setBlock(slice%Xext_linalg,xgcols_out,nrow,ncol,fcol=fcol_ext)
            call xgBlock_copy(xgcols_in,xgcols_out)
        end do
    end if

    ! Unitary test
    ABI_CHECK(cols(slice%Xext_linalg)==ncol,'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # cols of Xext_linalg ', xmpi_comm_rank(comm), cols(slice%Xext_linalg)

end subroutine slice_initExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_freeExtended
!! NAME
!! slice_freeExtended
!! 
!! SOURCE

subroutine slice_freeExtended(slice)

    implicit none

    type(slice_t), intent(inout) :: slice
   
    if (slice%nslice>1) then
        call xg_free(slice%Xext_)
        call xg_free(slice%Eext_)
        call xg_free(slice%Rext_)
    end if

end subroutine slice_freeExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_distributeExtended
!! NAME
!! slice_distributeExtended
!!
!! FUNCTION
!! Distribute Xext_linalg across *all* MPI processes and allocate 
!! its distributed version Xext on individual processes.
!! This routine applied transposition across *all* MPI processes.
!! After the transposition each process contains the correct
!! bandpp corresponding to the slice so that no additional communication
!! has to be performed in order to bring band slices to processes.
!!
!! SOURCE

subroutine slice_distributeExtended(slice,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nspinor

    ! Local variables -------------------------------
    integer :: specedim,neigenpairs
    integer, pointer :: ncolsColsRows_ptr(:) => null()

    ! *********************************************************************
 
    spacedim = slice%spacedim
    neigenpairs = slice%neigenpairs

    if (chebfi%paral_kgb == 1) then

        nprocs = xmpi_comm_size(comm(X0))

        ! Rule for number of bands per process
        ncolsColsRows_ptr => slice%mpiData%ncolsColsRows

        ! Allocate slice%Xext according to the target MPI distribution for slices
        call xgTransposer_constructor(slice%xgTransposerX,slice%Xext_linalg,slice%Xext,nspinor,&
            STATE_LINALG,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0_fft,&
            gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
            custom_ncolsColsRows=.true.,ncolsColsRows_sub=ncolsColsRows_ptr)
   
        slice%xgTransposerX%gpu_kokkos_nthrd  = slice%gpu_kokkos_nthrd
   
        ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE_XEXPAND)
        call xgTransposer_transpose(slice%xgTransposerX,STATE_COLSROWS)
        ABI_NVTX_END_RANGE()

    else
        call xgBlock_setBlock(slice%Xext, slice%Xext_linalg, spacedim, neigenpairs)
    end if

    ! Unitary test
    ABI_CHECK(cols(slice%Xext)==ncolsColsRows(xmpi_comm_rank(spacecom)),'wrong colsrows representation')
    write(*,'(a,i6,i6)') '# proc has # cols of Xext ', xmpi_comm_rank(spacecom), cols(slice%Xext)
    ! cols(slice%X_expand) should be equal to ncolsColsRows(i) where i is the rank of MPI process
    
end subroutine slice_distributeExtended
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run
!! NAME
!! slice_run
!! 
!! FUNCTION
!! Diagonalize slices in parallel. Input/output is the 
!! extended buffer in colsrows representation.
!! Notice that initial objects xgx0,eigen,residu are not input
!!
!! IML TODO add count
!! In early SCF iterations, count the number of Ritz values that fall within the
!! perturbed spectral interval of each slice, where the size of the perturbation
!! is related to the residual norm of each Ritz pair. We could therefore terminate
!! the subspace iterations when the counts no longer change.
!!
!! SOURCE

subroutine slice_run(slice,getAX_BX,getBm1X,nspinor)

    implicit none

    !Arguments ------------------------------------    
    type(slice_t), intent(inout) :: slice
    integer      , intent(in   ) :: nspinor
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
    integer :: tot_nrows,ncols_slice,spacecom
    integer :: color,my_rank,slice_comm,ierr
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    type(xgBlock_t) :: X0

    ! *********************************************************************

    call xgBlock_getSize(slice%X_expand,tot_nrows,ncols_slice)
    spacecom = comm(xgBlock_colsrows) ! global communicator
    ! spacedim is smaller than tot_nrows

    ! slice colsrows workspace points to column range of distributed X_expand
    call xgBlock_setBlock(slice%X_expand,X0,rows=tot_nrows,cols=ncols_slice)
    slice%X = X0

    ! if gpu: This is important! Because slice%Xext is on CPU
    ! call xgBlock_copy_to_gpu(slice%X)
    ! call xgBlock_set_gpu_option(slice%X)
    
    ! Synchronize before creating a subcomm
    call xmpi_barrier(spacecom)
    
    ! Split global communicator so that only procs with the same color communicate
    my_rank = xmpi_comm_rank(spacecom)
    color = mpiSlice_getSliceMe(mpi_slice,my_rank)
    call xmpi_comm_split(spacecom,color,my_rank,slice_comm,ierr)

    ! Restrict all communications to slice subcommunicator
    call xgBlock_setComm(slice%X,slice_comm) ! colsrows representation
    call xgBlock_setComm(slice%X_linalg,slice_comm) ! linalg representation

    ! Compute distribution of rows (plane waves) in the subcommunicator
    call mpiSlice_distributeRows(mpi_slice,slice_comm)
    nrowsLinalg_ptr => mpi_slice%lookup
    nproc = xmpi_comm_size(slice_comm)
    ABI_MALLOC(nrowsLinalg, (nproc)) ! FIXME size known before slice construction. Move out
    nrowsLinalg_ptr => nrowsLinalg
    call compute_uniform_distribution(nrowsLinalg_ptr,nproc,tot_nrows)

    ! Create slice transposer using subcommunicator and allocate slice%X_linalg
    call xgTransposer_constructor(slice%xgTransposerX,slice%X_linalg,slice%X,nspinor,&
        STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,slice_comm,0,0,chebfi%me_g0_fft,&
        gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
        custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! true to allow different bandpp (avoid pad)

    ncols_slice = cols(slice%X)

    ! Do the same for AX and BX ..
    call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,&
        chebfi%AX%self,chebfi%xAXColsRows,STATE_COLSROWS)
    call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,&
        chebfi%BX%self,chebfi%xBXColsRows,STATE_COLSROWS)

    chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
    chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

    ! Body of computation
    ! ===============
    
    ! Apply polynomial filtering (requires colsrows state)
    if (islice==0) then
        !call slice_applyLowpassFilter(slice,getAX_BX,getBm1X,nspinor)
        call slice_applyLowpassFilter(chebfi,getAX_BX,getBm1X,nspinor)
    else
        !call slice_applyBandpassFilter(slice,getAX_BX,getBm1X,nspinor)
        call slice_applyBandpassFilter(chebfi,getAX_BX,getBm1X,nspinor)
    end if

    ! Transpose to linalg state
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

        ! All MPI columns wait to finish
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

    ! Unitary test
    ABI_CHECK(rows(slice%X_linalg)==nrowsLinalg(xmpi_comm_rank(spacecom)),'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # rows of slice X ', xmpi_comm_rank(spacecom), rows(slice%X_linalg)

    ! Perform Rayleigh-Ritz and compute residuals (requires linalg state)
    call slice_RayleighRitz(slice,eigen,residu)
    
    ! Transpose to colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_TRANSPOSE)
    if (chebfi%paral_kgb == 1) then

        ! All MPI rows wait to finish
        call xmpi_barrier(chebfi%spacecom)

        call xgTransposer_transpose(chebfi%xgTransposerX, STATE_COLSROWS)

        !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
        if (xmpi_comm_size(chebfi%spacecom) == 1) then 
            call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
        end if
    else
        call xgBlock_setBlock(chebfi%xXColsRows, chebfi%X, spacedim, neigenpairs)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_transpose,2,tsec)

    ! TODO Deal with xgeigen and xgresidu
    !call xgBlock_reshape(slice%xgeigen, (/1,nband_slice/))
    !call xgBlock_reshape(slice%xgresidu, (/1,nband_slice/))
    !call xgBlock_copy_from_gpu(slice%xgeigen)
    !call xgBlock_copy_from_gpu(slice%xgresidu)
    !call slice_blockCopy(slice%xgeigen,sliceAll%xgeigen_ovlp,1,j1,nband_slice,j2)
    !call slice_blockCopy(slice%xgresidu,sliceAll%xgresidu_ovlp,1,j1,nband_slice,j2)

    ! Unitary test
    ABI_CHECK(cols(slice%X)==ncols_slice,'wrong colsrows representation')
    write(*,'(a,i6,i6)') '# proc has # cols of slice X ', xmpi_comm_rank(spacecom), cols(slice%X)

    ! Copy slice solution to the extended buffer (requires colsrows state)
    call xgBlock_copy(slice%X,X0)
    ! FIXME same for eigen, residu?

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
    if (gpu_option==ABI_GPU_KOKKOS) then
        call gpu_device_synchronize()
    end if
#endif

    ! Free transposer objects
    if (chebfi%paral_kgb == 1) then
        call xgTransposer_free(chebfi%xgTransposerX)
        call xgTransposer_free(chebfi%xgTransposerAX)
        call xgTransposer_free(chebfi%xgTransposerBX)
    end if
    
    ! Free memory
    if (allocated(nrowsLinalg)) ABI_FREE(nrowsLinalg)
    call xgTransposer_free(slice%xgTransposerX)

end subroutine slice_run
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_applyLowpassFilter
!! NAME
!! slice_applyLowpassFilter
!!
!! FUNCTION
!! Apply Lowpass filter using Chebyshev polynomial on a set of vectors.
!!
!! INPUTS
!!  slice   = spectral slice parameters
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm 
!!  on a single spectral slice
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_applyLowpassFilter(slice,getAX_BX,getBm1X,nspinors)

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nspinor
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
    type(chebfi_t) :: chebfi

    ! *********************************************************************

    chebfi = slice%chebfi

    ! When entering DivResult%self should contain all eigenvalues in slice
    ! in the form (bandpp,1) as a ROW vector
    comm_slice = slice%mpiData%comm_sub

    ! This is the maximum and the minimum eigenvalue in the slice
    if (slice%paral_kgb == 1) then
        call xmpi_max(slice%maxeig,maxeig_global,comm_slice,ierr)
        call xmpi_min(slice%mineig,mineig_global,comm_slice,ierr)
    else
        call xmpi_max(slice%maxeig,maxeig_global,comm,ierr)
        call xmpi_min(slice%mineig,mineig_global,comm,ierr)
    end if

    eigenvalues = DivResults%self !! ....
    ! DivResults%self is already constructed from previous routines
    ! this routine should also set the maxeig_global, mineig_global then

    if (chebfi%paral_kgb == 0) then
        ABI_MALLOC(ndeg_filter_bands,(neigenpairs))
    else
        ABI_MALLOC(ndeg_filter_bands,(bandpp))
    end if
    ndeg_filter_bands(:) = ndeg_filter
    
    ! Spectral interval to amplify is [-oo, lambda_minus)
    lambda_minus = maxeig_global
    lambda_plus = slice%%ecut

    center = (lambda_plus + lambda_minus)*0.5
    radius = (lambda_plus - lambda_minus)*0.5

    one_over_r = 1/radius
    two_over_r = 2/radius

    !A * Psi
    call timab(tim_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_CORE)
    do ideg = 0, ndeg_filter - 1

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, bandpp)
        end if
        ABI_NVTX_END_RANGE()

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

    ! Scale X,AX,BX by amplification factor to reduce large values
    call chebfi_ampfactor(chebfi, DivResults%self, lambda_minus, lambda_plus, ndeg_filter_bands)

    call xg_free(DivResults) ! en fait ne pas faire ça ici
    ! car on garde DivResults à l'extérieur des slices aussi pour
    ! comparer les convergences
    ABI_FREE(ndeg_filter_bands)

end subroutine slice_applyLowpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_applyBandpassFilter
!! NAME
!! slice_applyBandpassFilter
!!
!! FUNCTION
!! Apply Bandpass filter using Chebyshev-Jackson polynomial on a set of vectors.
!!
!! INPUTS
!!  slice   = spectral slice parameters
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  slice <type(slice_t)>=all data used to apply Spectrum Slicing algorithm 
!!  on a single spectral slice
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_applyBandpassFilter(slice,getAX_BX,getBm1X,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(slice_t), intent(inout) :: slice
    integer, intent(in) :: nspinor
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
    type(chebfi_t) :: chebfi
    integer :: space
    integer :: spacedim
    integer :: neigenpairs
    integer :: nline
    integer :: gpu_option
    integer :: nrows, ncols
    integer :: iline, ilinep1, iband, ierr
    real(dp) :: tolerance
    real(dp) :: one_over_r
    real(dp) :: two_over_r
    real(dp) :: center, radius
    real(dp) :: ck, mu, damp, tau  ! bandpass filter parameters
    real(dp) :: alow,bupp,low,upp  ! slice interval 
    type(xg_t) :: ChebyExpansion   ! Chebyshev expansion for vectors

    ! *********************************************************************

    chebfi = slice%chebfi

    ! Initialize solutions for slice using chebfi 
    space = chebfi%space
    spacedim = chebfi%spacedim
    neigenpairs = chebfi%neigenpairs 
    tolerance = chebfi%tolerance
    gpu_option = chebfi%gpu_option
    chebfi%eigenvalues = eigen
    chebfi%X = X0
    nrows = spacedim
    ncols = neigenpairs
    if (chebfi%paral_kgb==1) then
        nrows = chebfi%total_spacedim
        ncols = chebfi%bandpp
    end if
    
    ! Global spectral interval
    radius = (slice%gub - slice%glb)/2.d0   ! entire spectrum radius
    center = (slice%gub + slice%glb)/2.d0   ! entire spectrum center
    one_over_r = 1/radius
    two_over_r = 2/radius
    ! Target local to amplify scaled to [-1,1)
    nline = slice%degree                    ! polynomial filter degree
    low = slice%low                         ! filter support low bound
    upp = slice%upp                         ! filter support upper bound
    alow = (low - center) / radius          ! scaled filter support low
    bupp = (upp - center) / radius          ! scaled filter support upp

    ! AX_next=A*X -> 1 Hamiltonian application
    call chebfi_getAX_BX(chebfi, getAX_BX)

    ! B-orthonormalize X, BX and AX
    !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,gpu_option,AX=chebfi%xAXColsRows)
    ! IL TODO Deflate vectors 10/03/2025

    ! Why do this here?
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    write(std_out,*) 'TRACE initialize Chebyshev expansion (hopefuly on GPU)'
    ! Compute Chebyshev polynomial expansion on X iteratively on iline=0,nline
    ! Initialize Xsum = 0 (bands are distributed)
    ABI_NVTX_START_RANGE(NVTX_SLICE_EXPANSION)
    call xg_init(ChebyExpansion, chebfi%space, nrows, ncols, chebfi%spacecom, gpu_option=gpu_option)
    call xgBlock_zero(ChebyExpansion%self)
    ! X_next=X -> iline=0 Hamiltonian applications
    ck = Pi/(nline+2)
    mu = 1/Pi*(ACOS(alow)-ACOS(bupp))
    damp = 1.d0
    !Xsum = mu(0)*damp(0)*X_next + Xsum
    call xgBlock_saxpy(ChebyExpansion%self, mu*damp, chebfi%xXColsRows)
    ABI_NVTX_END_RANGE()

    write(std_out,*) 'TRACE start Slice core'
    ABI_NVTX_START_RANGE(NVTX_BANDPASS_CORE)
    do iline = 0, nline - 1  

        ! X_next=2/r*(AX_next-c*X_next)-X_prev, -> iline+1 Hamiltonian applications
        ABI_NVTX_START_RANGE(NVTX_CHEBFI2_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ! xXColsRows=X_next
        call chebfi_swapInnerBuffers(chebfi, nrows, ncols)

        ! Add new term to the Chebyshev expansion
        !Xsum = damp(i+1)*mu(i+1)*X_next + Xsum
        ABI_NVTX_START_RANGE(NVTX_CHEBFI_EXPANSION)
        ilinep1 = iline + 1
        mu = 2/Pi * (SIN(ilinep1*ACOS(alow)) - SIN(ilinep1*ACOS(bupp)))/ilinep1
        damp = ((1 - ilinep1/(nline+2))*SIN(ck)*COS(ilinep1*ck) + 1/(nline+2)*COS(ck)*SIN(ilinep1*ck))/SIN(ck)
        call xgBlock_saxpy(ChebyExpansion%self, mu*damp, chebfi%xXColsRows)
   
        ! Store term before exit
        ! AX_next=A*X_next -> iline+2 Hamiltonian applications
        if (iline==nline-1) then
            ! X_next=Xsum (copy Xsum to X_next)
            call xgBlock_copy(ChebyExpansion%self, chebfi%xXColsRows)
        end if
        ABI_NVTX_END_RANGE()

        ! Apply A and B to X
        call chebfi_getAX_BX(chebfi, getAX_BX)
    
    end do ! end iline
    ABI_NVTX_END_RANGE()

    ! All slice processes wait for filter done
    ! FIXME Is this necessary? No communication happens actually
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Free Chebyshev expansion workspace
    call xg_free(ChebyExpansion)

end subroutine slice_applyBandpassFilter
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_RayleighRitz
!! NAME
!! slice_RayleighRitz
!!
!! FUNCTION
!! Wrapper for xg_RayleighRitz + compute residual.
!! Assumes (X,AX,BX) have MPI row distribution.
!!
!! SOURCE

subroutine slice_RayleighRitz(slice,eigen,residu)

    ! Arguments ***
    type(slice_t), intent(inout) :: slice
    type(xgBlock_t), intent(inout) :: eigen
    type(xgBlock_t), intent(inout) :: residu
    ! Local variables ***
    type(chebfi_t) :: chebfi
    integer :: ierr

    chebfi = slice%chebfi
    chebfi%eigenvalues = eigen

    ! Apply Rayleigh-Ritz for each MPI row
    ABI_NVTX_START_RANGE(NVTX_CHEBFI2_RR)
    call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,chebfi%eigenvalues,ierr,0,tim_RR,&
        chebfi%gpu_option,solve_ax_bx=.true.)
    ABI_NVTX_END_RANGE()
    
    if ( ierr /= 0 ) then
        ABI_WARNING("RayleighRitz did not work, but continue anyway.")
    end if

    ! Compute residual for each MPI row and store it to AX
    if (chebfi%paw) then
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
    else
        call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
    end if

    ! Wait until all MPI rows have computed their residual
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Communicate MPI rows to compute residual norm squared
    call xgBlock_colwiseNorm2(chebfi%AX%self,residu)

end subroutine slice_RayleighRitz
!!***

!----------------------------------------------------------------------

!!****f* m_slice/mpiSlice_getSliceMe
!! NAME
!! mpiSlice_getSliceMe
!! 
!! SOURCE

integer function mpiSlice_getSliceMe(mpi_slice,my_rank) result(slice_me)

    type(mpiSlice_t), intent(in) :: mpi_slice
    integer, intent(in) :: my_rank

    slice_me = mpi_slice%lookup_slice_me(my_rank+1)

end function mpiSlice_getSliceMe
!!***

!----------------------------------------------------------------------

!!****f* m_slice/parallelEnv_init
!! NAME
!! parallelEnv_init
!! 
!! FUNCTION
!! Set parameters of slice parallelisation environment.
!! Assumes that spectral slices have already been splitted
!! (uses degree).
!! 
!! SOURCE

subroutine parallelEnv_init(env,distribution)

    implicit none

    type(parallelEnv_t), intent(inout) :: env
    

    if (.not. allocated(env%data_buffer)) then
        ABI_ALLOC(env%data_buffer(100))
    end if


    ! color MPI processes with the slice TODO automatize for arbitrary number of slices
    slice_color = -1
    if (my_rank==0) slice_color=0
    if (my_rank==1) slice_color=1
    if (my_rank==2) slice_color=1
    if (my_rank==3) slice_color=1
    slice_color = lookup(my_rank)

    !rule_nband(islice)= number of bands contained in slice
    !rule_ndeg(islice)= computational load for slice (eg degree)
    !rule_nproc(islice)= number of processes working for slice
    !lookup_bandpp_me(iproc)= number of bands contained in proc
    !lookup_slice_me(iproc)= slice index contained in proc

    ABI_MALLOC(mpiSlice%rule_nband, (nslice)) 
    ABI_MALLOC(mpiSlice%rule_ndeg, (nslice))
    do islice=1,nslice
        ABI_MALLOC(slice%mpiData(islice)%rule_nproc, (nslice))
    end do
    ABI_MALLOC(mpiSlice%lookup_bandpp_me, (nproc))
    ABI_MALLOC(mpiSlice%lookup_slice_me, (nproc))

    ! rule_nband: defined by the spectral splitter based on DOS
    ! rule_nproc: defined by the load balance strategy
    ! lookups: deduced

    !comm_slice(islice)=comm id working for slice
    ! Attention this must be constructed at the same time for all processes
    ! before doing any communication

    ! Used in mpi_enreg...
    bandpp = slice%bandpp

    ! Solve allocation problem to find the amount of resource allocated to each slice
    !In the first case, we are minimizing the maximum of mini/ximi​ni​/xi​ 
    !with xixi​ summing to pp, which aims to ensure none of the ratios is excessively large. 
    !In the second, the goal is to distribute resources so that the ratios mi/ximi​/xi​ 
    !are as equal as possible across all tasks. So it's about choice in balancing extremes versus uniformity.

    !Goal: We need to allocate pp resources to nn groups, where each group ii 
    !has a weight wiwi​ and size mimi​. The allocation xixi​ should be 
    !proportional to the weighted size wi×miwi​×mi​ for each group.

    nband_ptr => slice%nband_sub
    select case(resource_allocation)
    case(FAIR_ALLOCATION)
        ! fair share, equal division. Uses the same bandpp per process regardless ndeg
        ABI_MALLOC(ones, (nslice))
        ones(:) = 1 ! weight is 1
        ones_ptr => ones
        nproc_ptr = optimize_allocation(nband_ptr, ones_ptr, nproc)
        ABI_FREE(ones)
    case(WEIGHTED_FAIR_ALLOCATION)
        ! applies load balancing per process (=ndeg*bandpp)
        ! This ensures that each slice receives resources in proportion to its degree 
        ! and the number of vectors it has.
        ndeg_ptr => slice%ndeg_sub ! weight is degree
        nproc_ptr = optimize_allocation(nband_ptr, ndeg_ptr, nproc) 
    end select
    slice%mpiData%rule_nproc(:) = nproc_ptr(:)
    
    ! 'bandpp_slice'= uniformly distribute 'nband_slice' across 'nproc_slice' processes
    do islice=1,nslice
        bandpp_ptr => slice%mpiData(islice)%rule_nband
        nband_slice = slice%nband_sub(islice)
        nproc_slice = slice%mpiData(islice)%rule_nproc
        call uniform_distribution(bandpp_ptr,nproc_slice,nband_slice)
    end do

    ! Assumes every process has fixed capacity of bandpp
    !my_rank = 0
    !do i=1,nslice
    !    my_rank = my_rank + nband_slice(i) / bandpp
    !    proc_id = my_rank + 1 ! fortran indices start from 1!
    !    my_slice(proc_id) = i
    !end do

    ! For nrows linalg
    ABI_MALLOC(mpiSlice%nrowsLinalg, (nproc))

    map_slice_to_proc

    map_proc_to_slice

    ! lookup: my_rank+1 -> color
    ! example of number of vectors per process
    vectors = [3, 2, 5, 7, 1, 4, 8, 3, 6, 5]

    if (.not.allocated(env%group_lookup)) then
        ABI_MALLOC(env%group_lookup, (nproc))
    end if

    ! Initialization with invalid values
    env%group_lookup = -1
  
    total_vectors = 0
    current_slice = 1
    vectors_in_current_slice = 0

    ! Distribution of processes at groups/slices
    do i = 1, n
    
        ! Include vectors of current process
        vectors_in_current_slice = vectors_in_current_slice + vectors(i)
    
        ! If vectors of current process exceed limits of current slice
        ! TODO rename slice%nband number of bands per slice to slice_size
        if (vectors_in_current_slice > slice_size) then
            ! Update moving to next slice
            current_slice = current_slice + 1
        vectors_in_current_slice = vectors(i) ! initialize new slice with vectors of current process
        end if
    
        ! Assign process to the group/slice
        env%group_lookup(i) = current_slice
    end do


    device_id = xomp_get_device_num()
    env%use_host = (device_id < 1) ! outside target (-1), inside target host (0)
    env%use_device = (device_id > 0) ! inside target not host (device number)

    do i = 1, n
        write(*,'(a,i5,a,i5)') "Process ", i, " is in group ", group_assignment(i)
    end do

    lookup(1) = 0
    lookup(2) = 1
    lookup(3) = 1
    lookup(4) = 1
    ABI_MALLOC(lookup,(nproc))
    nband_prev = 0
    do i=1,nproc
        lookup(i) = nband_prev + mod((i-1) ! NO
        nband_prev = nband_prev + nband(i)
    end do

end subroutine parallelEnv_init
!***

subroutine parallelEnv_free(env)

    implicit none 

    type(parallelEnv_t), intent(inout) :: env

    if (allocated(env%group_lookup)) then
        ABI_FREE(env%group_lookup)
    end if

end subroutine parallelEnv_free

!----------------------------------------------------------------------

!!****f* m_slice/uniform_distribution
!! NAME
!! uniform_distribution
!!
!! FUNCTION
!  Distribute 'm_tot' into 'nproc' processes as uniformly as possible.
!! If 'mod(m_tot,nproc) =/= 0' give less charge to the last process.
!! Return 'm_distr' array with 'nproc' values. Note what we *do not* do:
!! pad m so that each process has a multiple of m_uni.
!! 
!! SOURCE

subroutine uniform_distribution(m_distr,nproc,m_tot)

    ! Arguments
    pointer, integer, intent(inout) :: m_distr
    integer, intent(in) :: nproc
    integer, intent(in) :: m_tot
    ! Local variables
    integer :: m_uni
    integer :: m_rem

    if (nproc > 1) then
        if (modulo(m_tot,nproc) == 0) then
            m_distr(1:nproc) = m_tot / nproc
        else
            m_uni = ceiling(real(m_tot) / real(nproc))
            m_rem = m_tot - (nproc-1)*m_uni
            m_distr(1:nproc-1) = m_uni
            m_distr(nproc) = m_rem
        end if
    else
        m_distr(1) = m_tot
    end if

end subroutine uniform_distribution
!!***

!----------------------------------------------------------------------

!!****f* m_slice/fair_allocation
!! NAME
!! fair_allocation
!! 
!! FUNCTION
!! Allocate resources so that m(i)/x(i) is approximately equal across i
!! 
!! SOURCE

function fair_charge_allocation(m, p) result(x)
    
    ! Arguments
    real(dp), intent(in) :: m(:)
    integer, intent(in) :: p
    ! Local variables
    integer :: x(size(m))
    integer :: s, i, remaining, max_i
    real(dp), allocatable :: ideal_x(:), frac_part(:)
    real(dp) :: total_m

    s = size(m)
    ABI_MALLOC(ideal_x,(s))
    ABI_MALLOC(frac_part,(s))

    total_m = sum(m)
    ideal_x = (m * real(p, dp)) / total_m

    do i = 1, s
        x(i) = floor(ideal_x(i))
        frac_part(i) = ideal_x(i) - real(x(i), dp)
    end do

    remaining = p - sum(x)

    ! Distribute remaining units to highest fractional parts
    do while (remaining > 0)
        max_i = maxloc(frac_part, 1)
        x(max_i) = x(max_i) + 1
        frac_part(max_i) = 0.0_dp  ! mark as used
        remaining = remaining - 1
    end do

    if (allocated(ideal_x)) ABI_FREE(ideal_x)
    if (allocated(frac_part)) ABI_FREE(frac_part)

end function fair_allocation
!!***

!----------------------------------------------------------------------

!!****f* m_slice/optimize_allocation
!! NAME
!! optimize_allocation
!! 
!! FUNCTION
!! Solve integer optimization problem under constraint: 
!! 
!!     min_{x_1,..,x_s} max_{1,..,s} f_i(x_i)
!!     subject to:   x_1 + .. + x_s = p
!!                   x_i integers
!! 
!! with objective cost function f_i(x)=m_i*n_i/x.
!! The solution x_i is the amount of resource allocated to the i-th task.
!! The algorithm uses binary search to deal with integer rounding.
!! 
!! INPUTS
!! arrays m and n (length s), and integer p
!! m can be the group size, n can be another measure or need ..
!! p in the number of total resources
!! 
!! OUTPUT
!! integer array x of size s such that sum(x) = p and max(m_i*n_i/x_i) is minimized
!!
!! SOURCE

function optimize_allocation(m, n, p) result(x)

    ! Arguments
    real(dp), intent(in) :: m(:), n(:)
    integer, intent(in) :: p
    ! Local variables
    integer :: x(size(m))
    integer :: s, i, total
    real(dp) :: t_low, t_high, t_mid, eps
    integer :: max_iter, iter
    real(dp), allocatable :: temp(:)
        
    s = size(m)
    ABI_MALLOC(temp,(s))    

    ! Binary search parameters
    t_low = 1.0e-6_dp
    t_high = maxval(m * n)
    eps = 1.0e-6_dp
    max_iter = 100

    do iter = 1, max_iter
        t_mid = (t_low + t_high) / 2.0_dp
        total = 0
        do i = 1, s
            x(i) = ceiling(m(i) * n(i) / t_mid)
            total = total + x(i)
        end do

        if (total > p) then
            t_low = t_mid
        else
            t_high = t_mid
        end if

        if (abs(t_high - t_low) < eps) exit
    end do

    ! Final rounding and optional adjustment
    total = 0
    do i = 1, s
        x(i) = ceiling(m(i) * n(i) / t_high)
        total = total + x(i)
    end do

    ! Distribute leftover units
    do while (total < p)
        real(dp) :: min_increase, current, next
        integer :: best_i
        min_increase = 1.0e9_dp
        best_i = -1
        do i = 1, s
            current = m(i) * n(i) / real(x(i), dp)
            next = m(i) * n(i) / real(x(i) + 1, dp)
            if (next - current < min_increase) then
                min_increase = next - current
                best_i = i
            end if
        end do
        x(best_i) = x(best_i) + 1
        total = total + 1
    end do

    if (allocated(temp)) ABI_FREE(temp)

end function optimize_allocation
!!***

! This is the same as before but uses greedy to deal with integer rounding
! I think this is not optimal

    ! Subroutine to perform weighted fair allocation with integer results
    subroutine greedy_allocation(w, m, p, n, x)
        real, dimension(n), intent(in) :: w  ! Weights of the groups
        integer, dimension(n), intent(in) :: m  ! Sizes of the groups (number of vectors)
        real, intent(in) :: p  ! Total resources available
        integer, intent(in) :: n  ! Number of groups
        integer, dimension(n), intent(out) :: x  ! Allocated resources for each group (integers)
        
        real :: total_weighted_size  ! Total weighted size (sum of w_i * m_i)
        real :: allocated_real(n)  ! Real-valued proportional allocation before rounding
        integer :: i, total_allocated, leftover

        ! Calculate the total weighted size (sum of w_i * m_i)
        total_weighted_size = 0.0
        do i = 1, n
            total_weighted_size = total_weighted_size + w(i) * m(i)
        end do

        ! Perform the proportional allocation: x_i = (w_i * m_i) / (sum(w_j * m_j)) * p
        total_allocated = 0
        do i = 1, n
            allocated_real(i) = (w(i) * m(i) / total_weighted_size) * p
            x(i) = nint(allocated_real(i))  ! Round to the nearest integer
            total_allocated = total_allocated + x(i)
        end do

        ! If the total allocation doesn't sum to p, adjust the allocations
        leftover = p - total_allocated
        do i = 1, leftover
            x(i) = x(i) + 1  ! Distribute the leftover resources
        end do

    end subroutine greedy_allocation

!----------------------------------------------------------------------

!!****f* m_slice/slice_unitTest
!! NAME
!! slice_unitTest
!! 
!! FUNCTION
!! Given a set of vectors X computes the columnwise norm squared, on:
!! - CPU sequential,
!! - CPU MPI,
!! - GPU sequential (if GPU available),
!! - GPU MPI (if GPU available).
!! Checks that all four are the same.
!! 
!! SOURCE

function slice_unitTest(X) result(ierr)

    implicit none

    !Arguments ------------------------------------
    type(xgBlock_t), intent(in) :: X
    integer                     :: ierr         
    
    !Local variables-------------------------------
    integer :: nrows,ncols
    integer :: comm
    integer :: gpu_option
    real(dp) :: id_cpu, id_from_gpu, id_to_gpu

! *********************************************************************

    ! Get nrows, nrcols, gpu_option of X
    call xgBlock_getSize(X,nrows,ncols)
    call xgBlock_get_gpu_option(X,gpu_option)
    call xgBlock_get_communicator(X,comm)
    ! above this is just comm(X)
    write(std_out,*) '<--*--> Start unitary test for matrix'
    write(std_out,*) 'nrows       =', nrows
    write(std_out,*) 'ncols       =', ncols
    write(std_out,*) 'gpu_option  =', gpu_option
    write(std_out,*) 'communicator=', comm

    id_cpu = -1.0d0
    id_from_gpu = -1.0d0
    id_to_gpu = -1.0d0

    if (gpu_option==ABI_GPU_DISABLED) then
        id_cpu = xgBlock_getid(X, comm)
        write(std_out,*) 'id_cpu     =', id_cpu
    end if        

#if defined HAVE_GPU && defined HAVE_OPENMP_OFFLOAD
    if (gpu_option/=ABI_GPU_DISABLED) then
        call xgBlock_copy_from_gpu(X)
        id_from_gpu = xgBlock_getid(X, comm)
        call xgBlock_copy_to_gpu(X)
        write(std_out,*) 'id_from_gpu=', id_from_gpu
    else
        call xgBlock_copy_to_gpu(X)
        call xgBlock_copy_from_gpu(X)
        id_to_gpu = xgBlock_getid(X, comm)
        call xgBlock_copy_to_gpu(X)
        call xgBlock_copy_from_gpu(X)
        write(std_out,*) 'id_to_gpu  =', id_to_gpu
    end if
#endif
    
    ! TODO getid on xmpi_comm_null
    ! then sum across comm
    ! should give the same as
    ! getid on comm
    ! this is useful for MPI

    ierr = 0
    if (id_cpu < 0) ierr = -1
    if (id_from_gpu < 0) ierr = -1
    if (id_to_gpu < 0) ierr = -1

end function slice_unitTest
!!***

!----------------------------------------------------------------------

!!****f* m_slice/chebfi_poly1
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

end module m_slice
!!***



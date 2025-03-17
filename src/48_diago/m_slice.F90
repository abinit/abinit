!!****f* ABINIT/m_slice
!! NAME
!! m_slice
!!
!! FUNCTION
!! This module contains the types and routines used to apply the
!! Spectrum Slicing method (uses xG tools). It mainly defines a 
!! 'slice' datatype and associated methods.
!!
!! COPYRIGHT
!! Copyright (C) 2018-2025 ABINIT group (Ioanna-Maria Lygatsika, Lucas Baguet)
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

    ! Several (private) parameters
    !-------------------------------------------------
    integer, parameter :: tim_sliceAll_dos          = 2161
    integer, parameter :: tim_sliceAll_init         = 2162
    integer, parameter :: tim_sliceAll_split        = 2163
    integer, parameter :: tim_sliceAll_merge        = 2165
    integer, parameter :: tim_sliceAll_free         = 2166
    integer, parameter :: tim_slice2_RR             = 2169
    integer, parameter :: tim_slice2_invovl         = 2170
    integer, parameter :: tim_slice2_barrier        = 2171
    integer, parameter :: tim_slice2_getAX_BX       = 2172
    integer, parameter :: tim_slice2_copy           = 2173
    integer, parameter :: tim_slice2_swap           = 2174
    integer, parameter :: tim_slice2_RR_q           = 2175
    integer, parameter :: tim_slice2_transpose      = 2176
    integer, parameter :: tim_slice2_residu         = 2177
    integer, parameter :: tim_slice2_postinvovl     = 2178
    integer, parameter :: tim_slice2_init           = 2179
    integer, parameter :: tim_slice2_free           = 2180
    integer, parameter :: tim_slice2_expansion      = 2181
    integer, parameter :: tim_slice_blockCopy       = 2182
    integer, parameter :: tim_slice_Acopy           = 2185

    ! Public 'slice' datatype
    ! Parameters specific to current slice, manages isolated slice memory
    !-------------------------------------------------
    type, public :: slice_t

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
        type(xg_t) :: OCCW                      ! useless memory space in Slicing
        type(chebfi_t) :: chebfi                ! workspace for slice diago
    
        ! Pointers to slice memory
        type(xgBlock_t) :: xgx0
        type(xgBlock_t) :: xgeigen
        type(xgBlock_t) :: xgresidu
        type(xgBlock_t) :: xgocc ! useless

    end type slice_t

    ! Public 'sliceAll' datatype
    ! Parameters common to all slices, manages overlap-free memory
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
        
        ! Variables not used in slicing but passed to chebfi (with default values cause not used)
        integer :: oracle                        ! Option to compute mdeg_filter from residuals
        integer :: nbdbuf                        ! Number of bands in the buffer
        real(dp) :: oracle_factor                ! factor used to decrease residuals
        real(dp) :: oracle_min_occ               ! threshold on occupancies used for nbdbuf=-101
        ! end chebfi variables 

        integer :: paral_kgb                     ! MPI parallelization variables
        integer :: bandpp
        integer :: comm_cols
        integer :: comm_rows
        integer :: me_g0
        integer :: me_g0_fft

        ! General slicing params
        integer :: nband_ovlp                        ! number of columns in overlap-free space
        integer :: npband                            ! number of MPI processes
        integer :: balfilter                         ! spectral partition technique
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

        ! DOS (all eigenvalues+residuals) before Slicing. Not distributed across MPI proc. 
        ! Useful for convergence test
        type(xg_t) :: Eig0                           ! Rayleigh quotients before Slicing
        type(xg_t) :: Res0                           ! residual norms before Slicing

        ! Overlap-free memory allocated for all slices
        type(xg_t) :: XW_ovlp                        ! input/ output eigenvectors
        type(xg_t) :: EW_ovlp                        ! output eigenvalues
        type(xg_t) :: RW_ovlp                        ! output residuals

        ! Pointers to overlap-free memory
        type(xgBlock_t) :: xgx0_ovlp
        type(xgBlock_t) :: xgeigen_ovlp
        type(xgBlock_t) :: xgresidu_ovlp

        ! Pointers to all slice parameters
        integer, pointer :: pband(:) => NULL()        ! Eigenvector column permutation
        integer, pointer :: idx(:,:) => NULL()        ! start and end indices per slice
        integer, pointer :: idx_merge(:,:) => NULL()  ! start and end converged indices per slice
        integer, pointer :: ndeg(:) => NULL()         ! filter degree per slice 
        integer, pointer :: npbandSlice(:) => NULL()
        real(dp), pointer :: sbound(:,:) => NULL()    ! (slice low bound,slice upp bound,
                                                      !  filter low support,filter upp support)

    end type sliceAll_t

    ! Public methods associated to 'slice' datatype
    !-------------------------------------------------
    public :: sliceAll_init
    public :: sliceAll_initOverlapFree
    public :: sliceAll_free
    public :: sliceAll_freeOverlapFree
    public :: sliceAll_dos
    public :: sliceAll_split
    public :: sliceAll_merge
    public :: slice_blockCopy
    public :: slice_findOptimalNumMpiProcs
    public :: slice_init
    public :: slice_run
    public :: slice_free
    public :: slice_unitTest

    CONTAINS  
!=====================================================================
!!***

!!****f* m_slice/sliceAll_init
!! NAME
!! sliceAll_init
!!
!! FUNCTION
!! Initialize a 'sliceAll' datastructure.
!!
!! INPUTS
!! Same as chebfi_init, plus nslice total number of slices.
!!
!! SOURCE

subroutine sliceAll_init(sliceAll,neigenpairs,spacedim,tolerance,ecut,paral_kgb,&
&                        bandpp,mdeg_filter,space,eigenProblem,spacecom,me_g0,me_g0_fft,&
&                        paw,comm_rows,comm_cols,nslice,npband,ramp,balance,&
&                        nbdbuf,oracle,oracle_factor,oracle_min_occ,gpu_option,gpu_kokkos_nthrd)

 implicit none

 ! Arguments ------------------------------------
 integer         , intent(in   ) :: bandpp
 integer         , intent(in   ) :: npband
 integer         , intent(in   ) :: nslice
 integer         , intent(in   ) :: eigenProblem
 integer         , intent(in   ) :: me_g0
 integer         , intent(in   ) :: me_g0_fft
 integer         , intent(in   ) :: neigenpairs
 integer         , intent(in   ) :: mdeg_filter
 integer         , intent(in   ) :: comm_cols
 integer         , intent(in   ) :: comm_rows
 integer         , intent(in   ) :: paral_kgb
 integer         , intent(in   ) :: space
 integer         , intent(in   ) :: spacecom
 integer         , intent(in   ) :: spacedim
 integer         , intent(in   ) :: balance
 integer         , intent(in   ) :: nbdbuf
 integer         , intent(in   ) :: oracle
 integer         , intent(in   ) :: gpu_option
 logical         , intent(in   ) :: paw
 real(dp)        , intent(in   ) :: ramp
 real(dp)        , intent(in   ) :: ecut
 real(dp)        , intent(in   ) :: tolerance
 real(dp)        , intent(in   ) :: oracle_factor
 real(dp)        , intent(in   ) :: oracle_min_occ
 type(sliceAll_t), intent(inout) :: sliceAll
 integer         , intent(in   ), optional :: gpu_kokkos_nthrd
 real(dp) :: tsec(2)

 ! *********************************************************************

 call timab(tim_sliceAll_init,1,tsec)
 
 sliceAll%space        = space
 sliceAll%neigenpairs  = neigenpairs
 sliceAll%spacedim     = spacedim
 sliceAll%tolerance    = tolerance
 sliceAll%ecut         = ecut
 sliceAll%paral_kgb    = paral_kgb
 sliceAll%comm_cols    = comm_cols
 sliceAll%bandpp       = bandpp
 sliceAll%comm_rows    = comm_rows
 sliceAll%mdeg_filter  = mdeg_filter
 sliceAll%spacecom     = spacecom
 sliceAll%eigenProblem = eigenProblem
 sliceAll%me_g0        = me_g0
 sliceAll%me_g0_fft    = me_g0_fft
 sliceAll%paw          = paw
 sliceAll%gpu_option   = gpu_option
 sliceAll%nslice       = nslice
 sliceAll%npband       = npband
 sliceAll%ramp         = ramp
 sliceAll%balfilter    = balance
 sliceAll%nband_ovlp   = neigenpairs ! see initOverlapFree
 sliceAll%nbdbuf       = nbdbuf
 sliceAll%oracle       = oracle
 sliceAll%oracle_factor = oracle_factor
 sliceAll%oracle_min_occ = oracle_min_occ

 if (present(gpu_kokkos_nthrd)) then
    sliceAll%gpu_kokkos_nthrd = gpu_kokkos_nthrd
 else
    sliceAll%gpu_kokkos_nthrd = 1
 end if

 ! Space of eigenvalues
 if (space==SPACE_C) then
    sliceAll%space_res = SPACE_C
 else if (space==SPACE_CR) then
    sliceAll%space_res = SPACE_R
 end if

 call sliceAll_allocateAll(sliceAll)
    
 call timab(tim_sliceAll_init,2,tsec)

end subroutine sliceAll_init
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_allocateAll
!! NAME
!! sliceAll_allocateAll
!! 
 
subroutine sliceAll_allocateAll(sliceAll)

    implicit none
    
    type(sliceAll_t), intent(inout) :: sliceAll
    integer :: space_res
    integer :: neigenpairs
    integer :: comm_cols
    integer :: gpu_option

    call sliceAll_free(sliceAll)

    space_res = sliceAll%space_res
    neigenpairs = sliceAll%neigenpairs
    comm_cols = sliceAll%comm_cols
    !gpu_option = sliceAll%gpu_option
    ! FIXME forced CPU 
    gpu_option = ABI_GPU_DISABLED

    ! Eigenvalues and residuals before slicing
    ! Every MPI process contains this array
    call xg_init(sliceAll%Eig0,space_res,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)
    call xg_init(sliceAll%Res0,SPACE_R,rows=1,cols=neigenpairs,comm=comm_cols,gpu_option=gpu_option)

end subroutine sliceAll_allocateAll
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_free
!! NAME
!! sliceAll_free
!! 
 
subroutine sliceAll_free(sliceAll)

    implicit none
    
    type(sliceAll_t), intent(inout) :: sliceAll
    real(dp) :: tsec(2)

    call timab(tim_sliceAll_free,1,tsec)

    call sliceAll_freeOverlapFree(sliceAll)
    call xg_free(sliceAll%Eig0)
    call xg_free(sliceAll%Res0)
    
    call timab(tim_sliceAll_free,2,tsec)
    
end subroutine sliceAll_free
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_initOverlapFree
!! NAME
!! sliceAll_initOverlapFree
!! 
 
subroutine sliceAll_initOverlapFree(sliceAll)

    implicit none
    
    ! Arguments ------------------------------------
    type(sliceAll_t), intent(inout) :: sliceAll
    ! Local variables-------------------------------    
    integer :: nslice,space,spacedim,spacecom,me_g0
    integer :: nband_ovlp,islice
    real(dp) :: tsec(2)
    
    ! *********************************************************************

    call timab(tim_sliceAll_init,1,tsec)
    
    nslice   = sliceAll%nslice
    space    = sliceAll%space
    spacedim = sliceAll%spacedim
    spacecom = sliceAll%spacecom
    me_g0    = sliceAll%me_g0

    ! This space is bigger than original one.
    ! Purpose: * Prevent race condition using overlap-safe memory space for concurrent 
    !            read-write on slices simultaneously.
    !          * Perform MPI coms and allocations only between CPU hosts. 
 
    nband_ovlp = sum((/ (sliceAll%idx(islice,2)-sliceAll%idx(islice,1)+1, islice=1,nslice) /))
    sliceAll%nband_ovlp = nband_ovlp
    
    write(std_out,'(a,i0)') '-----> Allocating overlap-free memory nband_ovlp=', nband_ovlp
    call xg_init(sliceAll%XW_ovlp,space,spacedim,nband_ovlp,spacecom,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
    call xg_init(sliceAll%EW_ovlp,SPACE_R,1,nband_ovlp,spacecom,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)
    call xg_init(sliceAll%RW_ovlp,SPACE_R,1,nband_ovlp,spacecom,me_g0=me_g0,gpu_option=ABI_GPU_DISABLED)

    ! Pointers to overlap-free memory
    sliceAll%xgx0_ovlp = sliceAll%XW_ovlp%self
    sliceAll%xgeigen_ovlp = sliceAll%EW_ovlp%self
    sliceAll%xgresidu_ovlp = sliceAll%RW_ovlp%self
    
    call timab(tim_sliceAll_init,2,tsec)
    
end subroutine sliceAll_initOverlapFree
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_freeOverlapFree
!! NAME
!! sliceAll_freeOverlapFree
!! 
 
subroutine sliceAll_freeOverlapFree(sliceAll)

    implicit none
    
    type(sliceAll_t), intent(inout) :: sliceAll
    real(dp) :: tsec(2)

    call timab(tim_sliceAll_free,1,tsec)
    
    ! Free overlap-free memory
    call xg_free(sliceAll%XW_ovlp)
    call xg_free(sliceAll%EW_ovlp) 
    call xg_free(sliceAll%RW_ovlp) 
    
    call timab(tim_sliceAll_free,2,tsec)
    
end subroutine sliceAll_freeOverlapFree
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_dos
!! NAME
!! sliceAll_dos
!! 
!! FUNCTION
!! Compute Density Of States (DOS) for all bands.
!! 
!! SOURCE

subroutine sliceAll_dos(sliceAll,X0,getAX_BX,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(sliceAll_t), intent(inout) :: sliceAll
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
    integer :: gpu_option,gpu_kokkos_nthrd,nrows
    integer :: nbdbuf,oracle
    integer :: my_rank,shift_row ! for MPI all to all
    logical :: paw
    real(dp) :: tolerance,ecut
    real(dp) :: oracle_factor,oracle_min_occ
    type(chebfi_t) :: chebfi
    type(xg_t) :: Eig0_XW, Res0_XW
    type(xgBlock_t) :: Eig0_paral, Res0_paral
    ! arrays
    real(dp) :: tsec(2)
 
    ! *********************************************************************

    call timab(tim_sliceAll_dos,1,tsec)
    
    space          = sliceAll%space
    spacedim       = sliceAll%spacedim
    neigenpairs    = sliceAll%neigenpairs
    space_res      = sliceAll%space_res
    bandpp         = sliceAll%bandpp
    mdeg_filter    = sliceAll%mdeg_filter
    tolerance      = sliceAll%tolerance
    ecut           = sliceAll%ecut
    paral_kgb      = sliceAll%paral_kgb
    spacecom       = sliceAll%spacecom
    me_g0          = sliceAll%me_g0
    me_g0_fft      = sliceAll%me_g0_fft
    eigenproblem   = sliceAll%eigenproblem
    paw            = sliceAll%paw
    comm_rows      = sliceAll%comm_rows
    comm_cols      = sliceAll%comm_cols
    nbdbuf         = sliceAll%nbdbuf
    oracle         = sliceAll%oracle
    oracle_factor  = sliceAll%oracle_factor
    oracle_min_occ = sliceAll%oracle_min_occ

    !gpu_option       = sliceAll%gpu_option
    ! FIXME forced CPU, to be done in GPU
    gpu_option = ABI_GPU_DISABLED
    gpu_kokkos_nthrd = sliceAll%gpu_kokkos_nthrd

    ! ==================== INITIALIZATION ===========================
    
    ! Initialize DOS workspace
    write(std_out,'(a,i0)') '-----> allocate chebfi for DOS, bandpp=', bandpp
    call chebfi_init(chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp,&
&                    mdeg_filter,nbdbuf,space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,&
&                    comm_cols,oracle,oracle_factor,oracle_min_occ,gpu_option,&
&                    gpu_kokkos_nthrd=gpu_kokkos_nthrd)

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

    ! Store to sliceAll object without distribution across MPI processes
    ! All processors have all eigenvalues and residuals
    ! FIXME remove if condition after partialcopy is debugged on GPU
    !if (gpu_option==ABI_GPU_OPENMP) then
    !    write(std_out,'(a)') 'Copy (Eig0,Res0) from gpu'
    !    call xgBlock_copy_from_gpu(Eig0_paral)
    !    call xgBlock_copy_from_gpu(Res0_paral)
    !end if
    if (xmpi_comm_size(comm_cols)>1) then

        ! Initialize entire array with zeros, every MPI contains this array
        call xgBlock_zero(sliceAll%Eig0%self)
        call xgBlock_zero(sliceAll%Res0%self)

        ! Recover rank of current MPI process
        my_rank = xmpi_comm_rank(comm_cols)
        shift_row = my_rank * bandpp

        ! Reshape in order to use blockCopy on columns (requires same number of rows)
        call xgBlock_reshape(Eig0_paral,(/1,bandpp/))     
        call xgBlock_reshape(Res0_paral,(/1,bandpp/))     
 
        ! Fill entire array with entry=(current MPI part | zero otherwise)
        call slice_blockCopy(Eig0_paral,sliceAll%Eig0%self,1,shift_row+1,bandpp,shift_row+bandpp)
        call slice_blockCopy(Res0_paral,sliceAll%Res0%self,1,shift_row+1,bandpp,shift_row+bandpp)
        
        ! Sum entire object across MPI, every MPI contains the same entire array
        call xgBlock_mpi_sum(sliceAll%Eig0%self,comm=comm_cols)
        call xgBlock_mpi_sum(sliceAll%Res0%self,comm=comm_cols)

        ! Undo reshape
        call xgBlock_reshape(sliceAll%Eig0%self,(/neigenpairs,1/))     
        call xgBlock_reshape(sliceAll%Res0%self,(/neigenpairs,1/))

    else
        call xgBlock_copy(Eig0_paral,sliceAll%Eig0%self)
        call xgBlock_copy(Res0_paral,sliceAll%Res0%self)
    end if
    
    ! Debug; problem is -7 is too low
    !write(std_out,*) 'Rayleigh values (after MPI row coms)'
    !call xgBlock_print(sliceAll%Eig0%self, std_out)

    ! Free workspace
    write(std_out,'(a)') 'free chebfi for DOS <-----'
    call chebfi_free(chebfi)
    call xg_free(Eig0_XW)
    call xg_free(Res0_XW)
 
    call timab(tim_sliceAll_dos,2,tsec)

end subroutine sliceAll_dos
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_split
!! NAME
!! sliceAll_split
!! 
!! FUNCTION
!! Split spectrum to slices. Always performed on CPU
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

subroutine sliceAll_split(sliceAll,idxAll,ndegAll,sboundAll,pband,npbandSlice)

    implicit none

    !Arguments ------------------------------------
    type(sliceAll_t), intent(inout) :: sliceAll
    integer, pointer, intent(inout) :: idxAll(:,:)
    integer, pointer, intent(inout) :: ndegAll(:)
    integer, pointer, intent(inout) :: pband(:)
    integer, pointer, intent(inout) :: npbandSlice(:)
    real(dp), pointer, intent(inout) :: sboundAll(:,:)
    
    !Local variables-------------------------------
    integer :: j,k,jmax,spos,nvec_ovlp,nvec,k1,k2
    integer :: nv_pad,ndeg,npband,comm_cols
    integer :: nslice,neigenpairs,nline
    integer :: ndeg_max = 200
    integer :: balance_option
    integer :: ipt,npt,iptL,iptR
    real(dp) :: ramp
    real(dp) :: tol12 = 1.0e-12
    real(dp) :: ecut,low,upp,glb,gub,c,r
    real(dp) :: lj,uj,wj,finL,finR,foutL,foutR
    real(dp) :: wlj,wuj
    real(dp) :: fptL,fptR,ptL,ptR,fptIn
    real(dp) :: a_,b_ ! target interval scaled in -1,1
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    real(dp), allocatable :: slice_cut(:)
    real(dp), allocatable :: resid_cut(:)
    real(dp), pointer :: resid_(:,:)
    real(dp), pointer :: theta_(:,:)
    real(dp), allocatable, target :: resid(:)
    real(dp), allocatable, target :: theta(:)
    !real(dp), pointer :: resid_ptr(:) => NULL()
    !real(dp), pointer :: theta_ptr(:) => NULL()
    real(dp) :: tsec(2)

! *********************************************************************

    call timab(tim_sliceAll_split,1,tsec)
   
    ! Interval slicing parameters
    npband         = sliceAll%npband
    comm_cols      = sliceAll%comm_cols
    neigenpairs    = sliceAll%neigenpairs
    nslice         = sliceAll%nslice
    ecut           = sliceAll%ecut
    nline          = sliceAll%mdeg_filter
    balance_option = sliceAll%balfilter
    ramp           = sliceAll%ramp

    ! Set pointers to eigenvalue and residual memory
    Eig0_all = sliceAll%Eig0%self
    Res0_all = sliceAll%Res0%self

    ! Results could be complex (with null imaginary part), so neigenpairs has to be in cols, not rows
    call xgBlock_reverseMap(Eig0_all,theta_,rows=1,cols=neigenpairs)
    call xgBlock_reverseMap(Res0_all,resid_,rows=1,cols=neigenpairs)

    ! Save theta_,resid_ first row in theta,resid
    ABI_MALLOC(theta,(neigenpairs))
    ABI_MALLOC(resid,(neigenpairs))
    theta(1:neigenpairs) = theta_(1,1:neigenpairs)
    resid(1:neigenpairs) = resid_(1,1:neigenpairs)
    !theta_ptr => theta
    !resid_ptr => resid

    ! Sort thetas in increasing order and store permutation
    call sort_dp(neigenpairs,theta,pband,tol12)
    
    ! TODO sort resid wrt pband as well

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
    sliceAll%gub = gub
    sliceAll%glb = glb
   
    ! Initialize spectral cuts
    ABI_MALLOC(slice_cut,(nslice+1))
    slice_cut(:) = 0.d0
    slice_cut(1) = low
    slice_cut(nslice+1) = upp
    
    ABI_MALLOC(resid_cut,(nslice+1))
    resid_cut(:) = 0.d0
    resid_cut(1) = sqrt(resid(pband(1)))
    resid_cut(nslice+1) = sqrt(resid(pband(neigenpairs)))

    ! Compute spectral cuts on interior slices
    select case(balance_option)
    case(1)
        ! Balance interval widths
        slice_cut(2:nslice) = (/ (low+(upp-low)/nslice*j, j=1,nslice-1) /)
        ! resid cut not implemented in that case
    case(2)
        ! Balance number of vectors
        slice_cut(2:nslice) = (/ (theta(neigenpairs/nslice*j), j=1,nslice-1) /)
        resid_cut(2:nslice) = (/ (sqrt(resid(pband(neigenpairs/nslice*j))), j=1,nslice-1) /)
    end select

    ! Operation count estimation for each slice
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

        ! Uncomment following two lines to use overlap width that depends on residuals
        !wlj = resid_cut(j)
        !wuj = resid_cut(j+1)
        ! IML 28/01/2025 this does not guarantee a posteriori ratio < ramp. 
        ! Also results in very high ndeg (like 139). 

        write(std_out,*) '------------ /Divide/ Slice ',j
        
        ! Uncomment following two lines to use wj fixed overlap width
        wlj = wj
        wuj = wj
        write(std_out,*) '       overlap widths=', wlj, wuj
        ! *******************************************************************
        
        ndeg = nline

        ! Filter support [l-w,u+w) scaled to [-1,1)
        a_ = (lj-wlj-c)/r
        b_ = (uj+wuj-c)/r
       
        ! Count eigenvalues in slice cut plus overlap 
        call count_values(lj-wlj,uj+wuj,theta,spos,neigenpairs,k1,k2,nvec)
        nvec_ovlp = nvec
        write(std_out,*) '    after overlap', k1,k2

        ! Add padding for MPI efficiency
        if (npband > 1) then
            nv_pad = pad_size(nvec,npband)
            call shift_indices(spos,neigenpairs,nv_pad,k1,k2,nvec)
        end if
        write(std_out,*) '    after     pad', k1,k2

        ! Oscillations **********************************************
        !                                         IML 29/01/2025
        ! The convergence ratio r0/rN depends on the ratio
        ! [min f(lambda_in)]/[max f(lambda_out)] approximated as
        ! f(l)/f(l-w) and f(u)/f(u+w). However this approximation
        ! might not be accurate, due to oscillations.
        ! ***********************************************************

        ! Compute optimal degree in slice with overlap
        ! For Chebyshev, must compute the ratio of convergence
        if (j > 1) then
            ndeg = 4
            finL = 0.d0; finR = 0.d0; foutL = 1.d0; foutR = 1.d0
            do while ( (finL/foutL < ramp) .and. (finR/foutR < ramp) .and. (ndeg<ndeg_max) )
                finL  = bandpassIndicator_sca((lj-c)/r,a_,b_,ndeg)
                foutL = bandpassIndicator_sca(a_      ,a_,b_,ndeg)
                finR  = bandpassIndicator_sca((uj-c)/r,a_,b_,ndeg)
                foutR = bandpassIndicator_sca(b_      ,a_,b_,ndeg)
                ndeg = ndeg + 1
            end do
        else
            !write(std_out,*) 'first slice apriori=', 1.d0/cheb_poly1(lj,12,uj+wuj,ecut)
        end if
        ! FIXME impose fixed ndeg to test speed
        !ndeg = 62

        ! Compare max f(lambda_out) (oscillations) vs. f(l-w),f(u+w) (approximation)
        ! For first slice situation is quite different, ndeg not tuned for the moment
        !write(std_out,*) ' '
        if (j>1) then
            !write(std_out,*) 'Filter oscillations for slice=', lj,uj
            !write(std_out,*) '                        width=', uj-lj
            !write(std_out,*) '                    scaled to=', (lj-c)/r,(uj-c)/r
            npt = 1000 ! number of points to test for oscillations
            ! Maximum value in [glb,l) and [u,ulb), outside slice [u,l)
            foutL = bandpassIndicator_sca(a_,a_,b_,ndeg)
            foutR = bandpassIndicator_sca(b_,a_,b_,ndeg)
            iptL = maxloc((/ (bandpassIndicator_sca((glb+(ipt-1)*(lj-glb)/npt-c)/r,a_,b_,ndeg),ipt=1,npt) /),dim=1)
            iptR = maxloc((/ (bandpassIndicator_sca((uj+(ipt-1)*(gub-uj)/npt-c)/r,a_,b_,ndeg),ipt=1,npt) /),dim=1)
            ptL = (glb+(iptL-1)*(lj-glb)/npt-c)/r
            ptR = (uj+(iptR-1)*(gub-uj)/npt-c)/r
            fptL = bandpassIndicator_sca(ptL,a_,b_,ndeg)
            fptR = bandpassIndicator_sca(ptR,a_,b_,ndeg)
            !write(std_out,*) '***** Outer: max   f(t)=', fptL,fptR 
            !write(std_out,*) '             argmax   t=', ptL,ptR
            !write(std_out,*) '                  estim=', foutL,foutR
            ! Maximum value inside slice [u,l)
            finL = bandpassIndicator_sca((lj-c)/r,a_,b_,ndeg)
            finR = bandpassIndicator_sca((uj-c)/r,a_,b_,ndeg)
            iptL = minloc((/ (bandpassIndicator_sca((lj+(ipt-1)*(uj-lj)/npt-c)/r,a_,b_,ndeg),ipt=1,npt) /),dim=1)
            ptL = (lj+(iptL-1)*(uj-lj)/npt-c)/r
            fptIn = bandpassIndicator_sca(ptL,a_,b_,ndeg)
            !write(std_out,*) '***** Inner: min   f(t)=', fptIn
            !write(std_out,*) '             argmin   t=', ptL
            !write(std_out,*) '                  estim=', finL,finR
            !write(std_out,*) '======= Convergence rate inner/outer=',fptIn/max(fptL,fptR)
            ! This convergence rate should be compared to the a posteriori convergence rate
        else
            !write(std_out,*) 'First slice=', lj,uj
            !write(std_out,*) '      width=', uj-lj
            !write(std_out,*) '  scaled to=', (lj-c)/r,(uj-c)/r
        end if
        !write(std_out,*) ' '

        ! Plot filter
        npt = 100
        if (j>1) then
            !write(std_out,*) ' '
            !write(std_out,*) 'Plot filter ==== x | f(x)'
            
            ! This is for convergence rate
            !fptR = max(bandpassIndicator_sca((lj-wlj-c)/r,a_,b_,ndeg),bandpassIndicator_sca((uj+wuj-c)/r,a_,b_,ndeg))
            do ipt=1,npt
                ptL = (lj + (ipt-1)*(uj-lj)/npt - c)/r
                fptL = bandpassIndicator_sca(ptL,a_,b_,ndeg)
                !write(std_out,*) ptL, fptL
                !write(std_out,*) ptL, fptL, 'conv rate approx=', fptL/fptR ! current inn
            end do

            !write(std_out,*) ' '
        else
            !write(std_out,*) ' '
            fptR = cheb_poly1(uj+wuj,12,uj+wuj,gub) ! max out
            do ipt=1,npt
                ! uj,gub is the interval mapped to -1,1
                ! in this interval Chebyshev poly is bounded by 1
                ! so uj,gub is the interval to ignore
                ptL = lj + (ipt-1)*(uj-lj)/npt
                fptL = cheb_poly1(ptL,12,uj+wuj,gub)
                ! FIXME need to take absolute value: T_n<0 for n odd
                !write(std_out,*) ptL, fptL, 'conv rate approx=', fptL/fptR ! current inn
            end do
            !write(std_out,*) ' '
        end if

        ! Print slice interval info
        write(std_out,*) 'Without overlap=', lj, uj
        write(std_out,*) '          width=', uj-lj
        write(std_out,*) 'With    overlap=', lj-wlj, uj+wuj
        write(std_out,*) '          width=', uj+wuj-lj+wlj
        write(std_out,*) '  nvec(balance)=', neigenpairs/nslice
        write(std_out,*) '      nvec_ovlp=', nvec_ovlp
        write(std_out,*) '           nvec=', nvec
        write(std_out,*) '           ndeg=', ndeg
        write(std_out,*) ' '
        ! end print

        ! Store slice parameters
        idxAll(j,1:2) = (/k1,k2/)
        ndegAll(j) = ndeg
        sboundAll(j,1:4) = (/lj,uj,lj-wlj,uj+wuj/) 
           
    end do

    ! Update pointers for all slices
    sliceAll%pband => pband
    sliceAll%idx => idxAll 
    sliceAll%ndeg => ndegAll
    sliceAll%sbound => sboundAll
    sliceAll%npbandSlice => npbandSlice ! TODO distribute

    ! Free workspace not needed
    if (allocated(slice_cut)) ABI_FREE(slice_cut)
    if (allocated(resid_cut)) ABI_FREE(resid_cut)
    if (allocated(theta)) ABI_FREE(theta)
    if (allocated(resid)) ABI_FREE(resid)
   
    call timab(tim_sliceAll_split,2,tsec)

end subroutine sliceAll_split
!!***

!----------------------------------------------------------------------

!!****f* m_slice/sliceAll_merge
!! NAME
!! sliceAll_merge
!! 
!! FUNCTION
!! Find converged eigenvalues in each slice using a criterion.
!! Return range of first and last index to merge per slice, both 
!! computed from a residual criterion on eigenvalues.
!! 
!! SOURCE

subroutine sliceAll_merge(sliceAll,idx_ovlp,idx_merge,merge_option)

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
    real(dp) :: convr_L,convr_R ! convergence rates
    type(xgBlock_t) :: Eig0_all, Res0_all
    ! arrays
    integer, pointer :: pband(:)
    real(dp), pointer :: theta0(:,:), resid0(:,:)
    real(dp), pointer :: thetaN(:,:), residN(:,:)
    !real(dp), pointer :: theta0_(:), resid0_(:)
    !real(dp), pointer :: thetaN_(:), residN_(:)
    real(dp), allocatable :: convr_slice(:)
    real(dp), allocatable :: resid0_slice(:)
    real(dp), allocatable :: residN_slice(:)
    real(dp), allocatable :: residNapprox_slice(:)
    real(dp), allocatable :: theta_slice(:)
    real(dp) :: tsec(2)
 
    ! *********************************************************************

    call timab(tim_sliceAll_merge,1,tsec)

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

    ! Workaround for dimensions 
    !theta0_ => theta0(1,1:neigenpairs)
    !resid0_ => resid0(1,1:neigenpairs)
    !thetaN_ => thetaN(1,1:nband_ovlp)
    !residN_ => residN(1,1:nband_ovlp)

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
        ABI_MALLOC(convr_slice,(nband_slice))
        
        ABI_MALLOC(residN_slice,(nband_slice))
        ABI_MALLOC(resid0_slice,(nband_slice))
        ABI_MALLOC(residNapprox_slice,(nband_slice))

        !theta_slice(1:nband_slice) = 1.d0/thetaN(1,j1:j2)+(supp-slow)/2.d0 ! after slicing
        theta_slice(1:nband_slice) = thetaN(1, j1:j2) ! after slicing
        !write(std_out,*) 'debug', i1,i2,j1,j2
        residN_slice(1:nband_slice) = sqrt(residN(1,j1:j2))
        resid0_slice(1:nband_slice) = sqrt(resid0(1,pband(i1:i2)))

        convr_slice(1:nband_slice) = sqrt(resid0(1,pband(i1:i2)))/sqrt(residN(1,j1:j2))

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
        write(std_out,*) 'max residual rN       in Partition =', maxval(sqrt(residN(1,j1:j2)))
        write(std_out,*) 'min rate r0/rN        in Partition =', minval(convr_slice(k1:k2))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Partition (kept) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Partition (kept) =', maxval(theta_slice(k1:k2))
 
        ! Compute a posteriori estimation r0/rN ~ min(f(lambda_in))/max(f(lambda_out))
        !k_out = max(k1-1,1)
        !convr_L = bandpassIndicator_sca((theta_slice(k1)-c)/r,(flow-c)/r,(fupp-c)/r,ndeg) / &
        !&          bandpassIndicator_sca((theta_slice(k_out)-c)/r,(flow-c)/r,(fupp-c)/r,ndeg)
        !k_out = min(k2+1,nband_slice)
        !convr_R = bandpassIndicator_sca((theta_slice(k2)-c)/r,(flow-c)/r,(fupp-c)/r,ndeg) / &
        !          bandpassIndicator_sca((theta_slice(k_out)-c)/r,(flow-c)/r,(fupp-c)/r,ndeg)
        !write(std_out,*) '        L and R approximation in slice=', convr_L, convr_R

        !if (nvec_count < nband_slice) write(std_out,*) 'slice is not full' 

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
        write(std_out,*) 'worst convergence rate Support        =', minval(convr_slice(k1:k2))
        ! Print boundaries of this kept range
        write(std_out,*) 'min eigenvalue in Support (converged) =', minval(theta_slice(k1:k2))
        write(std_out,*) 'max eigenvalue in Support (converged) =', maxval(theta_slice(k1:k2))

        ! Number of eigenvalues with convergence ratio smaller than target
        ! normally this number should be larger than number of kept vals in slice
        !call count_values(1.1d0,10000000.d0,convr_slice,spos,nband_slice,k1,k2,nvec_count)
        !write(std_out,*) 'converged eigvals nb=', nvec_count, minval(convr_slice)

        ! Compute approximation rn approx r0*max(fout)/fin for every filter support eigenvalue [flow,fupp)
!        if (islice==1) then
!            residNapprox_slice(k1:k2) = (/(resid0_slice(k)/cheb_poly1(theta_slice(k),12,fupp,gub), k=k1,k2)/)
!        else
!            ! convr_L=max f(outer eigenvalues)
!            convr_L = max(bandpassIndicator_sca((flow-c)/r,(flow-c)/r,(fupp-c)/r,ndeg),&
!&                         bandpassIndicator_sca((fupp-c)/r,(flow-c)/r,(fupp-c)/r,ndeg))
!            residNapprox_slice(k1:k2) = (/(resid0_slice(k)*convr_L/&
!&                         bandpassIndicator_sca((theta_slice(k)-c)/r,(flow-c)/r,(fupp-c)/r,ndeg),k=k1,k2)/) 
!        end if 
        
        ! Print results eigenvalue, residualN and residual0 (useful to get convergence rate)
        !write(std_out,*) ''
        !write(std_out,*) 'min computed eigenvalue in slice=', minval(theta_slice(1:nband_slice))
        !write(std_out,*) 'max computed eigenvalue in slice=', maxval(theta_slice(1:nband_slice))
        !write(std_out,*) ''
        !write(std_out,*) 'eig=          | resN=         | res0=        ' 
        !do k=k1,k2
!            write(std_out,*) 'eig=', theta_slice(k), 'res_ex=', residN_slice(k), 'res_app=', residNapprox_slice(k)!, 
!            write(std_out,*) theta_slice(k), residN_slice(k), resid0_slice(k) !, 
!&            'rerr=', abs(residN_slice(k) - residNapprox_slice(k))
!        end do
!        write(std_out,*) ''


        if (allocated(theta_slice)) ABI_FREE(theta_slice)
        if (allocated(convr_slice)) ABI_FREE(convr_slice)
        if (allocated(residN_slice)) ABI_FREE(residN_slice)
        if (allocated(resid0_slice)) ABI_FREE(resid0_slice)
        if (allocated(residNapprox_slice)) ABI_FREE(residNapprox_slice)

    end do
    if (nband_conv<neigenpairs) ABI_ERROR("not enough eigenvalues converged")
    if (nband_conv>neigenpairs) ABI_ERROR("too many eigenvalues converged")
   
    !write(std_out,*) 'theta0=', theta0(1,pband(:))
    !write(std_out,*) 'thetaN=', thetaN(1,:)

    !write(std_out,*) 'resid0=', resid0(1,pband(:))
    !write(std_out,*) 'residN=', residN(1,:)

    ! TODO 
    ! * count how many thetaN_ are in low,upp for every slice
    ! * count how many thetaN_ are outside current slice, and if they converged
    ! * count how many are in overlap region
    ! This will help diagnostic convergence "slice full"

    !if (merge_option==0) then
    !    ! Present loop can be parallelized
    !    slow = theta0(1,1)
    !    supp = theta0(1,neigenpairs)
    !    do islice=1,nslice
    !        i1 = sliceAll%idx(islice,1)
    !        i2 = sliceAll%idx(islice,2)
    !        nband_slice = i2 - i1 + 1
    !        ! FIXME array range
    !        !do iband=i1,i2
    !        !    write(std_out,*) 'resid0/residN=', resid0(1,iband)/residN(1,iband)
    !        !end do
    !        !upp = sliceAll%sbound(islice)
    !        ! Kept indices
    !        j1 = maxloc(thetaN_, dim=1, mask=(thetaN_ < slow)) + 1
    !        j2 = j1 + nband_slice - 1
    !        !i2 = maxloc(theta, dim=1, mask=(theta < upp))
    !        idx_merge(islice,1) = j1
    !        idx_merge(islice,2) = j2
    !        slow = thetaN(1,j2)
    !    end do
    !else if (merge_option==1) then
    !    ! Present loop is not parallel
    !end if
 
    call timab(tim_sliceAll_merge,2,tsec)

end subroutine sliceAll_merge
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run
!! NAME
!! slice_run
!! 
!! FUNCTION
!! Run Spectrum Slicing on a given slice.
!! First slice has Chebyshev filter next ones bandpass Jackson-Chebyshev expansion.
!! 
!! SOURCE

subroutine slice_run(slice,getAX_BX,getBm1X,nspinor)

    implicit none

    ! Arguments ------------------------------------
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

    ! Local variables-------------------------------    
    integer :: islice
    ! typed objects
    type(xgBlock_t) :: X0
    type(xgBlock_t) :: eigen
    type(xgBlock_t) :: residu
    type(xgBlock_t) :: occ ! IL TODO not used in slicing but passed as arg
                           ! idea is to use it to store residuals before slice
 
    ! *********************************************************************

    islice = slice%islice
    X0     = slice%xgx0
    eigen  = slice%xgeigen
    residu = slice%xgresidu
    occ    = slice%xgocc

    if (islice==1) then
        write(std_out,'(a)') 'chebfi_run= ...'
        call chebfi_run(slice%chebfi,X0,getAX_BX,getBm1X,eigen,occ,residu,nspinor)
    else
        write(std_out,'(a)') 'slice_run= ...'
        call slice_run_bandpass(slice%chebfi,slice,X0,getAX_BX,getBm1X,eigen,occ,residu,nspinor)
    end if
    
    ! IML TODO add count
    ! In early SCF iterations, count the number of Ritz values that fall within the
    ! perturbed spectral interval of each slice, where the size of the perturbation
    ! is related to the residual norm of each Ritz pair. We could therefore terminate
    ! the subspace iterations when the counts no longer change.
 
end subroutine slice_run
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_run_bandpass
!! NAME
!! slice_run_bandpass
!!
!! FUNCTION
!! Apply the Spectrum Slicing algorithm on a set of vectors.
!!
!! INPUTS
!!  chebfi    = memory workspace for eigenpairs, residuals
!!  slice     = spectral slice parameters
!!  mpi_enreg = information about MPI parallelization
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X= pointer to the function giving B^-1|X>
!!           B is typically the overlap operator S
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!  slice <type(chebfi_t)>=all data used to apply Polynomial Filtering algorithm 
!!  on a single spectral slice
!!  eigen= Full eigenvalues (initial values on entry)
!!  residu= residuals, i.e. norm of (A-lambdaB)|X>
!!  X0= Full set of vectors (initial values on entry)
!!
!! SOURCE

subroutine slice_run_bandpass(chebfi,slice,X0,getAX_BX,getBm1X,eigen,occ,residu,nspinor)

 implicit none

!Arguments ------------------------------------
 type(chebfi_t) , intent(inout) :: chebfi
 type(slice_t)  , intent(inout) :: slice
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
 integer :: space
 integer :: spacedim
 integer :: neigenpairs
 integer :: nline
 integer :: gpu_option
 integer :: nrows, ncols
 integer :: iline, ilinep1, iband, ierr
 real(dp) :: sigma ! shift for harmonic RR
 real(dp) :: tolerance
 real(dp) :: one_over_r
 real(dp) :: two_over_r
 real(dp) :: center, radius
 real(dp) :: ck, mu, damp, tau  ! bandpass filter parameters
 real(dp) :: alow,bupp,low,upp  ! slice interval 
 !type(xg_t) :: DivResults
 type(xg_t) :: ChebyExpansion   ! Chebyshev expansion for vectors
!arrays
 real(dp) :: tsec(2)

! *********************************************************************

 write(std_out,*) 'TRACE initializing'
 
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

 ! Filter parameters passed from slice object
 ! Global spectral interval
 radius = (slice%gub - slice%glb)/2.d0   ! entire spectrum radius
 center = (slice%gub + slice%glb)/2.d0   ! entire spectrum center
 one_over_r = 1/radius
 two_over_r = 2/radius
 ! Target local to amplify scaled in [-1,1)
 nline = slice%degree                    ! polynomial filter degree
 low = slice%low                         ! filter support low bound
 upp = slice%upp                         ! filter support upper bound
 alow = (low - center) / radius          ! scaled filter support low
 bupp = (upp - center) / radius          ! scaled filter support upp

 ! Transpose
 if (chebfi%paral_kgb == 1) then

   call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
     STATE_LINALG,TRANS_ALL2ALL,chebfi%comm_rows,chebfi%comm_cols,0,0,chebfi%me_g0_fft,gpu_option=gpu_option)

   call xgTransposer_copyConstructor(chebfi%xgTransposerAX,chebfi%xgTransposerX,chebfi%AX%self,chebfi%xAXColsRows,STATE_LINALG)
   call xgTransposer_copyConstructor(chebfi%xgTransposerBX,chebfi%xgTransposerX,chebfi%BX%self,chebfi%xBXColsRows,STATE_LINALG)

   chebfi%xgTransposerX%gpu_kokkos_nthrd  = chebfi%gpu_kokkos_nthrd
   chebfi%xgTransposerAX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd
   chebfi%xgTransposerBX%gpu_kokkos_nthrd = chebfi%gpu_kokkos_nthrd

   write(std_out,*) 'TRACE start MPI Transpose'
   call timab(tim_slice2_transpose,1,tsec)
   ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
   call xgTransposer_transpose(chebfi%xgTransposerX,STATE_COLSROWS)
   chebfi%xgTransposerAX%state = STATE_COLSROWS
   chebfi%xgTransposerBX%state = STATE_COLSROWS
   ABI_NVTX_END_RANGE()
   call timab(tim_slice2_transpose,2,tsec)
   write(std_out,*) 'TRACE finished MPI Transpose'

 else
   call xgBlock_setBlock(chebfi%X, chebfi%xXColsRows, spacedim, neigenpairs)         !use xXColsRows instead of X notion
   call xgBlock_setBlock(chebfi%AX%self, chebfi%xAXColsRows, spacedim, neigenpairs)  !use xAXColsRows instead of AX notion
   call xgBlock_setBlock(chebfi%BX%self, chebfi%xBXColsRows, spacedim, neigenpairs)
 end if

 write(std_out,*) 'TRACE start Hamiltonian application (hopefuly on GPU)'
 ! AX_next=A*X -> 1 Hamiltonian application
 call timab(tim_slice2_getAX_BX,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
 call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
 call xgBlock_zero_im_g0(chebfi%xAXColsRows)
 call xgBlock_zero_im_g0(chebfi%xBXColsRows)
 ABI_NVTX_END_RANGE()
 call timab(tim_slice2_getAX_BX,2,tsec)
 write(std_out,*) 'TRACE finished Hamiltonian application'

 ! B-orthonormalize X, BX and AX
 !call xg_Borthonormalize(chebfi%xXColsRows,chebfi%xBxColsRows,ierr,1,gpu_option,AX=chebfi%xAXColsRows)
 ! IL TODO Deflate vectors 10/03/2025

 if (chebfi%paral_kgb == 1) then
   call timab(tim_slice2_barrier,1,tsec)
   call xmpi_barrier(chebfi%spacecom)
   call timab(tim_slice2_barrier,2,tsec)
 end if

 write(std_out,*) 'TRACE initialize Chebyshev expansion (hopefuly on GPU)'
 ! Compute Chebyshev polynomial expansion on X iteratively on iline=0,nline
 ! Initialize Xsum = 0 (bands are distributed)
 call timab(tim_slice2_expansion,1,tsec)
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
 call timab(tim_slice2_expansion,2,tsec)

 write(std_out,*) 'TRACE start Slice core'
 ABI_NVTX_START_RANGE(NVTX_SLICE_CORE)
 do iline = 0, nline - 1  

    !write(std_out,*) 'TRACE start next order'
    ! X_next=2/r*(AX_next-c*X_next)-X_prev, -> iline+1 Hamiltonian applications
    ABI_NVTX_START_RANGE(NVTX_SLICE_NEXT_ORDER)
    call slice_computeNextOrderChebfiPolynom(chebfi, iline, center, one_over_r, two_over_r, getBm1X)
    ABI_NVTX_END_RANGE()
    !write(std_out,*) 'TRACE finished next order'

    ! xXColsRows=X_next
    call timab(tim_slice2_swap,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_SWAP_BUF)
    call chebfi_swapInnerBuffers(chebfi, nrows, ncols)
    ABI_NVTX_END_RANGE()
    call timab(tim_slice2_swap,2,tsec)

    ! Add term into expansion
    !Xsum = damp(i+1)*mu(i+1)*X_next + Xsum
    call timab(tim_slice2_expansion,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_EXPANSION)
    ilinep1 = iline + 1
    mu = 2/Pi * (SIN(ilinep1*ACOS(alow)) - SIN(ilinep1*ACOS(bupp)))/ilinep1
    damp = ((1 - ilinep1/(nline+2))*SIN(ck)*COS(ilinep1*ck) + 1/(nline+2)*COS(ck)*SIN(ilinep1*ck))/SIN(ck)
    call xgBlock_saxpy(ChebyExpansion%self, mu*damp, chebfi%xXColsRows)
   
    ! Store term before exit
    ! AX_next=A*X_next -> iline+2 Hamiltonian applications
    if (iline==nline-1) then
        ! X_next=Xsum (copy Xsum to X_next)
        call timab(tim_slice2_copy, 1, tsec)
        call xgBlock_copy(ChebyExpansion%self, chebfi%xXColsRows)
        call timab(tim_slice2_copy, 2, tsec)
    end if
    ABI_NVTX_END_RANGE()
    call timab(tim_slice2_expansion,2,tsec)

    call timab(tim_slice2_getAX_BX,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_SLICE_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_slice2_getAX_BX,2,tsec)

 end do ! end iline
 ABI_NVTX_END_RANGE()
 write(std_out,*) 'TRACE finished Slice core'

 if (chebfi%paral_kgb == 1) then
   call timab(tim_slice2_barrier,1,tsec)
   call xmpi_barrier(chebfi%spacecom)
   call timab(tim_slice2_barrier,2,tsec)
 end if

 ! Free work space
 call xg_free(ChebyExpansion)

 ! Transpose back (MPI)
 call timab(tim_slice2_transpose,1,tsec)
 ABI_NVTX_START_RANGE(NVTX_SLICE_TRANSPOSE)
 if (chebfi%paral_kgb == 1) then
   call xmpi_barrier(chebfi%spacecom)

   call xgTransposer_transpose(chebfi%xgTransposerX, STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerAX,STATE_LINALG)
   call xgTransposer_transpose(chebfi%xgTransposerBX,STATE_LINALG)

   if (xmpi_comm_size(chebfi%spacecom) == 1) then 
     !only one MPI proc reset buffers to right addresses (because of X-Xcolwise swaps)
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
 call timab(tim_slice2_transpose,2,tsec)

 ! Get eigenvectors from X
 ABI_NVTX_START_RANGE(NVTX_SLICE_RR)
 !call xg_Borthonormalize(chebfi%X,chebfi%BX%self,ierr,1,gpu_option,AX=chebfi%AX%self)
 call xg_RayleighRitz(chebfi%X,chebfi%AX%self,chebfi%BX%self,chebfi%eigenvalues,ierr,0,&
&                     tim_slice2_RR,gpu_option,solve_ax_bx=.true.)
 ABI_NVTX_END_RANGE()
 if ( ierr /= 0 ) then
    ABI_WARNING("RayleighRitz did not work, but continue anyway.")
 end if

 ! chebfi%AX=AX-eig*BX
 call timab(tim_slice2_residu, 1, tsec)
 if (chebfi%paw) then
   call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%BX%self,chebfi%AX%self)
 else
   call xgBlock_colwiseCymax(chebfi%AX%self,chebfi%eigenvalues,chebfi%X,chebfi%AX%self)
 end if
 ! residu=|AX-eig*BX|^2 (norm squared!)
 call xgBlock_colwiseNorm2(chebfi%AX%self, residu)
 call timab(tim_slice2_residu, 2, tsec)

 ! Store modified chebfi workspace X result to slice
 call timab(tim_slice2_copy, 1, tsec)
 call xgBlock_copy(chebfi%X, X0)
 call timab(tim_slice2_copy, 2, tsec)

#if defined(HAVE_GPU_CUDA) && defined(HAVE_YAKL)
   if (gpu_option==ABI_GPU_KOKKOS) then
     call gpu_device_synchronize()
   end if
#endif

 if (chebfi%paral_kgb == 1) then
   call xgTransposer_free(chebfi%xgTransposerX)
   call xgTransposer_free(chebfi%xgTransposerAX)
   call xgTransposer_free(chebfi%xgTransposerBX)
 end if

end subroutine slice_run_bandpass
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_computeNextOrderChebfiPolynom
!! NAME
!! slice_computeNextOrderChebfiPolynom
!!
!! FUNCTION
!! From P_n(B-^1.A)|X> (where P_n is the Chebyshev polynom of order n),
!!   computes P_n+1(B-^1.A)|X>
!! NOTE; this is the same as the function in m_chebfi2 just with different
!! timers+nvtx vars given as module variables
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

subroutine slice_computeNextOrderChebfiPolynom(chebfi,ideg,center,one_over_r,two_over_r,getBm1X)

 implicit none

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
   !write(std_out,*) 'TRACE start Bm1X (invovl)'
   call timab(tim_slice2_invovl, 1, tsec)
   ABI_NVTX_START_RANGE(NVTX_SLICE_GET_BM1X)
   call getBm1X(chebfi%xAXColsRows, chebfi%X_next)
   ABI_NVTX_END_RANGE()
   call timab(tim_slice2_invovl, 2, tsec)
   !write(std_out,*) 'TRACE finished Bm1X (invovl)'
 else
   call timab(tim_slice2_copy, 1, tsec)
   call xgBlock_copy(chebfi%xAXColsRows,chebfi%X_next)
   call timab(tim_slice2_copy, 2, tsec)
 end if

 !write(std_out,*) 'TRACE start postinvovl'
 call timab(tim_slice2_postinvovl, 1, tsec)
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
 call timab(tim_slice2_postinvovl, 2, tsec)
! write(std_out,*) 'TRACE finished postinvovl'

end subroutine slice_computeNextOrderChebfiPolynom
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
    integer :: gpu_option,gpu_kokkos_nthrd
    integer :: nbdbuf, oracle
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
    nbdbuf           = sliceAll%nbdbuf
    oracle           = sliceAll%oracle
    oracle_factor    = sliceAll%oracle_factor
    oracle_min_occ   = sliceAll%oracle_min_occ

    ! Various parameters of current slice
    nband  = slice%nband
    bandpp = slice%bandpp
    ndeg   = slice%degree
  
    call slice_free(slice)

    ! Every slice MPI proc has entire space
    ! FIXME EW,RW devrait être alloués exactement comme xgx0 et xgresidu in chebfiwf 
    write(std_out,'(a,i0)') '-----> Allocating slice result memory nband_slice=', nband
    call xg_init(slice%XW,space,spacedim,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%EW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%RW,SPACE_R,1,nband,spacecom,me_g0=me_g0,gpu_option=gpu_option)
    call xg_init(slice%OCCW,SPACE_R,nband,1,gpu_option=gpu_option)! fill with zero

    ! Every slice MPI proc has a part (bandpp) of this space
    write(std_out,'(a,i0)') '-----> Allocating slice working memory bandpp=', bandpp
    call chebfi_init(slice%chebfi,nband,spacedim,tolerance,ecut,paral_kgb,bandpp,&
&                    ndeg,nbdbuf,space,1,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols,&
&                    oracle,oracle_factor,oracle_min_occ,gpu_option,&
&                    gpu_kokkos_nthrd=gpu_kokkos_nthrd)

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

!!****f* m_slice/slice_split
!! NAME
!! slice_split
!! 
!! FUNCTION
!! Split spectrum to slices. 
!!
!! SOURCE

! IML 24/01/2025 not used for the moment
subroutine slice_split(slices,nslice,nband,npband,theta,resid,iperm,ecut,ramp,balance,nline)

    implicit none

    !Arguments ------------------------------------
    type(slice_t), pointer, intent(inout) :: slices(:)
    integer, pointer, intent(inout) :: iperm(:)
    integer, intent(in) :: nband,balance,nline,npband,nslice
    real(dp), intent(in) :: ecut,ramp
    real(dp), target, intent(inout) :: theta(nband)
    real(dp), target, intent(in) :: resid(nband)
    
    !Local variables-------------------------------
    integer :: i,j,k,jmax,spos,nvec_ovlp,nvec,k1,k2
    integer :: nv_offset,nv_extra,nv_pad,ndeg
    integer :: jperm(nband-1)
    real(dp) :: tol12 = 1.0e-12
    real(dp) :: low,upp,glb,gub,gam,center,radius
    real(dp) :: lj,uj,lopj,uopj,wj
    real(dp) :: consdiff(nband-1)
    real(dp) :: slice_cut(nslice+1)
    real(dp) :: tsec(2)

! *********************************************************************

    call timab(tim_sliceAll_split,1,tsec)
   
    ! Sort thetas in increasing order
    call sort_dp(nband,theta,iperm,tol12)
    
    ! Spectrum bounds
    low = theta(1)
    upp = theta(nband)
    ! Bounds that contain entire spectrum (guaranteed)
    glb = low - sqrt(resid(1))
    gub = ecut
   
    ! Initialize spectral cuts
    slice_cut(:) = 0.d0
    slice_cut(1) = low
    slice_cut(nslice+1) = upp

    ! Compute spectral cuts on interior slices
    select case(balance)
    case(0) ! Balance spectral gaps = where eigenvalues are less concentrated

        ! Sort consecutive differences by increasing order
        jperm = (/ (k, k=1,nband-1) /)
        consdiff = (/ (theta(k+1) - theta(k), k=1,nband-1) /)
        call sort_dp(nband-1,consdiff,jperm,tol12)

        ! Take median of largest gaps
        do i=1,nslice
            jmax = jperm(nband-i)
            slice_cut(i+1) = (theta(jmax) + theta(jmax+1)) / 2.d0
        end do

    case(1)

        ! Balance interval widths
        slice_cut(2:nslice) = (/ (low+(upp-low)/nslice*j, j=1,nslice-1) /)

    case(2)

        ! Balance number of vectors
        slice_cut(2:nslice) = (/ (theta(nband/nslice*j), j=1,nslice-1) /)

    end select

    ! Operation count estimation for each slice
    do j=1,nslice
            
        ! Slice position, first (1), interior (2), last (3)
        spos = 2
        if (j==1) spos = 1
        if (j==nslice) spos = 3

        ! Count eigenvalues in slice cut plus overlap 
        lj = slice_cut(j)
        uj = slice_cut(j+1)
        wj = (uj - lj) / 10.d0
        call count_values(lj-wj,uj+wj,theta,spos,nband,k1,k2,nvec)
        nvec_ovlp = nvec

        ! Add extra fraction of eigenvalues
        ! This prevents eigenvalues from moving between slices
        write(std_out,*) 'before', k1,k2
        nv_offset = 15 ! add this as ABINIT input var
        nv_extra = ceiling(nv_offset * nvec / 100.d0)
        call shift_indices(spos,nband,nv_extra,k1,k2,nvec)
        write(std_out,*) 'after +15%', k1,k2
         
        ! Add padding for MPI efficiency
        if (npband > 1) then
            nv_pad = pad_size(nvec,npband)
            call shift_indices(spos,nband,nv_pad,k1,k2,nvec)
        end if
        write(std_out,*) 'after pad', k1,k2

        ! Compute optimal degree in enlarged slice (Overlap+Pad)
        ndeg = nline
        lopj = theta(k1)
        uopj = theta(k2)
        center = (glb + gub) / 2.d0
        radius = (gub - glb) / 2.d0
        if (j > 1) then
            if (uopj - lopj < 0.d1) then
                write(std_out,*) 'Slice ',j,' of ',i,'is too thin.'
                continue
            end if
            call slice_optimize_bandpass(lopj,uopj,center,radius,ramp,gam,ndeg)
        end if

        write(std_out,*) 'slice ',j,' of ',i,':nvec_ovlp=',nvec_ovlp,',nvec=',nvec,',ndeg=',ndeg
            
    
    end do
    
    call timab(tim_sliceAll_split,2,tsec)

end subroutine slice_split
!!***

!----------------------------------------------------------------------

!!****f* m_slice/slice_optimize_bandpass
!! NAME
!! slice_optimize_bandpass
!! 
!! FUNCTION
!! Optimize center and degree of Schofield (2012) filter for
!! achieving target amplification tau>0 on interval [a,b).
!! Center and radius are of the largest interval.
!! 
!! SOURCE

subroutine slice_optimize_bandpass(a,b,center,radius,tau,gam,ndeg)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in ) :: a, b, center, radius, tau
    real(dp), intent(out) :: gam
    integer , intent(out) :: ndeg

    !Local variables-------------------------------
    integer :: k
    integer :: ndeg_max = 200
    integer :: root_maxiter = 50
    real(dp) :: ared, bred, ared0, bred0
    real(dp) :: eps = 1e-7

! *********************************************************************
    
    ared = (a-center)/radius
    bred = (b-center)/radius
    gam = (ared + bred) * 0.5d0 ! initialize center on midpoint
    
    do k=2,ndeg_max

        ! Balance filter center s.t. f(a)=f(b)
        ared0 = ared; bred0 = bred
        gam = find_root(ared0,bred0,root_maxiter,eps,ared,bred,k)
        
        ! Interval is very thin, impose median
        if (gam > bred .or. gam < ared) then
            gam = (ared + bred) * 0.5d0 
        end if
        
        ! Optimize filter degree s.t. amplification < tau outside [a,b)
        if (ABS(bandpass_sca(bred,k,gam)) < tau .and. ABS(bandpass_sca(ared,k,gam)) < tau) then
            ndeg = k
            exit
        end if
        
    end do

end subroutine slice_optimize_bandpass
!!***

!----------------------------------------------------------------------

!!****f* m_slice/find_root
!! NAME
!! find_root
!!
!! FUNCTION
!! Solves equation f(x) = 0 using bisection method on [x0,x1]
!!
!! INPUTS
!!  x0=       initial lower bound of interval
!!  x1=       initial upper bound of interval
!!  fun=      pointer to the function giving f(x) for any x 
!!  maxiter=  maximum bisections to perform
!!  esp=      error tolerance used as termination criterion
!!
!! OUTPUT
!!  root
!!
!! SIDE EFFECTS
!!
!! SOURCE

function find_root(x0, x1, maxiter, eps, ared, bred, k) result(root)

    implicit none

    !Arguments ------------------------------------
    real(dp), intent(in   ) :: eps, ared, bred
    integer , intent(in   ) :: maxiter, k
    real(dp), intent(inout) :: x0, x1

    real(dp) :: root
    
    !Local variables-------------------------------
    real(dp) :: y0, y1, y2, x2
    integer  :: i
    
! *********************************************************************

    root = 1d10
    do i = 1, maxiter
        y0 = bandpass_sca(bred, k, x0) - bandpass_sca(ared, k, x0) ! fun(x0)
        y1 = bandpass_sca(bred, k, x1) - bandpass_sca(ared, k, x1) ! fun(x1)
        if (ABS(y0) < eps) then
            root = x0
            exit
        else if (ABS(y1) < eps) then
            root = x1
            exit
        else
            x2 = (x0 + x1) * 0.5
            y2 = bandpass_sca(bred, k, x2) - bandpass_sca(ared, k, x2) ! fun(x2)
            if (y0 * y2 < 0) then
                x1 = x2
            else
                x0 = x2
            end if
        end if
    end do

end function find_root
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
!! SIDE EFFECTS
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
!! SIDE EFFECTS
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

!!****f* m_slice/slice_merge
!! NAME
!! slice_merge
!! 
!! FUNCTION
!! Returns indices of eigenpairs whose eigenvalues lie in [a,b]
!! 
!! SOURCE

subroutine slice_merge(slice,eig,lowb,uppb,imin,imax)

    implicit none

    type(slice_t), intent(inout) :: slice
    real(dp), intent(in) :: lowb, uppb, eig
    integer, intent(inout) :: imin, imax

    !imin = slice%idx(i), eig < lowb   + 1 
    !imax = slice%idx(i), eig < uppb

end subroutine slice_merge
!!***

!----------------------------------------------------------------------

!!****f* m_slice/count_values
!! NAME
!! count_values
!! 
!! FUNCTION
!! Return number of consecutive theta values in interval [a,b)
!! and their index range (first and last indices).
!! Slice position (spos) takes into account limiting indices.
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

!!****f* m_slice/shift_indices
!! NAME
!! shift_indices
!! 
!! FUNCTION
!! Shift first and last slice indices by some number
!! 
!! SOURCE

subroutine shift_indices(spos,nband,nv_extra,i1,i2,nvec)

    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: spos,nband,nv_extra
    integer, intent(inout) :: i1,i2,nvec
    !Local variables-------------------------------
    integer :: nv_extra_half

! *********************************************************************

    nv_extra_half = Int(ceiling(nv_extra / 2.d0))
    select case(spos)
    case(1) ! first slice 
        i1 = 1
        i2 = i2 + nv_extra
    case(3) ! last slice
        i1 = i1 - nv_extra
        i2 = nband
    case(2) ! interior slice
        i1 = i1 - nv_extra_half
        i2 = i2 + nv_extra - nv_extra_half
    end select
    nvec = i2 - i1 + 1

end subroutine shift_indices
!!***

!----------------------------------------------------------------------

!!****f* m_slice/pad_size
!! NAME
!! pad_size
!! 
!! FUNCTION
!! Insert padding to size for efficient parallelisation over nproc 
!! processes. Result is minimum increment of size_a divisible by nproc.
!!
!! SOURCE

function pad_size(size_a, nproc) result(incr_a)

    implicit none

    !Arguments ------------------------------------
    integer, intent(in) :: size_a
    integer, intent(in) :: nproc
    
    !Local variables-------------------------------
    integer :: incr_a

! *********************************************************************

    incr_a = nproc - modulo(size_a, nproc)

end function pad_size
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
    on_host = xomp_is_initial_device()
    write(std_out,'(a,i0,a,i0,a,i0)') 'Memcopy: gpu_optionA ',gpua,' gpu_optionB ',gpub, ' on_host ', on_host
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
    ncolsA_block = a2_ - a1 + 1
    ncolsB_block = b2_ - b1 + 1

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

!!****f* m_slice/slice_findOptimalNumMpiProcs
!! NAME
!! slice_findOptimalNumMpiProcs
!! 
!! FUNCTION
!! Distributes 'nproc' number of MPI processes across 'nslice' slices from 
!! given parameters: number of vectors 'nvec' and filter degree 'ndeg' per slice.
!! The distribution is optimal in theory, in the sense that is it obtained as 
!! the solution to the discrete optimisation problem under constraint:
!! 
!!     min_{x(1),..,x(s)} max_{1,..,s}   nvec(i) * ndeg(i) / x(i)
!!     x(1) + .. + x(s) = nproc
!! 
!! where x(i) is the unknown number of MPI processes used by i-th slice.
!! The cost function to be minimised represents the operation count of 
!! the polynomial filtering, which is 'ndeg' Hamiltonian applications 
!! applied to 'nvec/x' vectors using the MPI band distribution.
!!
!! SOURCE

subroutine slice_findOptimalNumMpiProcs(nproc,nslice,idx,ndeg,nproc_opt)

    implicit none

    !Arguments ------------------------------------
    integer,          intent(in   ) :: nproc
    integer,          intent(in   ) :: nslice
    integer, pointer, intent(in   ) :: idx(:,:)
    integer, pointer, intent(in   ) :: ndeg(:)
    integer, pointer, intent(inout) :: nproc_opt(:)
    
    !Local variables-------------------------------
    integer :: i
    integer :: i_most_charged
    real(dp) :: sum_tot
    integer, allocatable :: nvec(:)

! *********************************************************************

    ABI_MALLOC(nvec, (nslice))
    nvec(:) = (/ (idx(i,2) - idx(i,1) + 1, i=1,nslice) /)
    sum_tot = dot_product(nvec,ndeg)

    write(std_out,*) 'DEBUG find optimal num of mpi'
    write(std_out,*) 'nvec=', nvec(:)
    write(std_out,*) 'ndeg=', ndeg(:)
    write(std_out,*) 'sum_tot=', sum_tot
    
    ! for all slices, initialize to floor
    nproc_opt(1:nslice) = (/ (max(floor(nproc*nvec(i)*ndeg(i)/sum_tot),1), i=1,nslice) /)
    
    ! find max charged slice
    i_most_charged = maxloc( (/ (nvec(i)*ndeg(i), i=1,nslice) /), dim=1)
    write(std_out,*) 'most charged slice is', i_most_charged
    ! add remainder in max charge slice to sum to nproc
    nproc_opt(i_most_charged) = nproc - sum(nproc_opt(1:nslice)) + nproc_opt(i_most_charged)

    if (allocated(nvec)) ABI_FREE(nvec)
   
    ! must also assure that nvec divides nproc_opt!
    ! possibly add pad ?
 
end subroutine slice_findOptimalNumMpiProcs
!!***

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



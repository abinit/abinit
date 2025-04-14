


    type, public :: polyfi_t 

        type(chebfi_t) :: chebfi
        real(dp) :: mineig_global
        real(dp) :: maxeig_global
        real(dp) :: lambda_minus
        real(dp) :: lambda_plus

    end type polyfi_t


    subroutine polyfi_init(

chebfi,neigenpairs,spacedim,tolerance,ecut,paral_kgb,bandpp, &
                       ndeg_filter,nbdbuf,space,eigenProblem,spacecom,me_g0,me_g0_fft,paw,comm_rows,comm_cols, &
                       oracle,oracle_factor,oracle_min_occ,gpu_option,gpu_kokkos_nthrd,gpu_thread_limit,from_linalg

polyfi,nband_sub,dtset%tolwfr_diago,dtset%ecut,&
     dtset%paral_kgb,space,1,spacecom_sub,&
     me_g0,me_g0_fft,l_paw,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
     mineig_global,maxeig_global,lambda_minus,lambda_plus,nrowsLinalg,!!!
     l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd,&
     gpu_thread_limit=dtset%gpu_thread_limi
        )


        

 call chebfi_init(polyfi%chebfi,nband_sub,dtset%tolwfr_diago,dtset%ecut,&
     dtset%paral_kgb,space,1,spacecom_sub,&
     me_g0,me_g0_fft,l_paw,l_mpi_enreg%comm_spinorfft,l_mpi_enreg%comm_band,&
     mineig_global,maxeig_global,lambda_minus,lambda_plus,nrowsLinalg,!!!
     l_gs_hamk%gpu_option,gpu_kokkos_nthrd=dtset%gpu_kokkos_nthrd,&
     gpu_thread_limit=dtset%gpu_thread_limit)


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

    ! Oracle hard-coded deactivation
    oracle = 0
    oracle_factor = 1.d0
    oracle_min_occ = 1

 end subroutine polyfi_init
  

!----------------------------------------------------------------------

!!****f* m_polyfi/polyfi_run
!! NAME
!! polyfi_run
!! 
!! FUNCTION
!! Diagonalize individual slice in parallel. Remember this is in distributed 
!! parallel region and there is no inter-slice communication.
!!
!! INPUT
!! X0     =eigenvector guess distributed for slice in colsrows representation
!! chebfi =data structure used to compute Chebyshev polynomial 
!! 
!! SOURCE

subroutine polyfi_run(chebfi,X0,getAX_BX,getBm1X,eigen,residu,nspinor)

    implicit none

    !Arguments ------------------------------------    
    type(chebfi_t), intent(inout) :: chebfi
    type(xgBlock_t), intent(inout) :: X0
    integer       , intent(in   ) :: nspinor
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
    integer :: color,my_rank,ierr
    integer, target, allocatable :: nrowsLinalg(:)
    integer, pointer :: nrowsLinalg_ptr(:) => null()
    type(xgBlock_t) :: X0
    type(chebfi_t) :: chebfi

    ! *********************************************************************
 
    ! Memory workspace Is Chebfi
    chebfi = polyfi%chebfi

    ! Parameters
    total_spacedim = chebfi%total_spacedim
    neigenpairs = chebfi%neigenpairs
    spacecom = chebfi%spacecom
    comm_size = xmpi_comm_size(spacecom)

    ! Arrays
    if(.not.allocated(nrowsLinalg)) ABI_MALLOC(nrowsLinalg,(comm_size))
    nrowsLinalg_ptr => nrowsLinalg
    nrowsLinalg(:) = polyfi%nrows_blockrows

    ! Note: we want to use the memory space of X0 but not the
    ! same pointer because it is common for all slices. For this
    ! reason we create a new xgBlock independent of X0 for slice.
    ! Do not do that: chebfi%xXColsRows = X0!!
    call xgBlock_setBlock(X0, chebfi%xXColsRows, total_spacedim, neigenpairs)

    if (polyfi%gpu_option==ABI_GPU_OFFLOAD) then
        ! Because X0 is on CPU but chebfi%xXColsRows on GPU
        ! FIXME either after or before setBlock
        ! call xgBlock_copy_from_gpu(chebfi%xXColsRows)
        ! call xgBlock_set_gpu_option(chebfi%xXColsRows,ABI_GPU_OFFLOAD)
    end if

    ! Restrict all communications to current subcommunicator
    call xgBlock_setComm(chebfi%xXColsRows,spacecom)
    call xgBlock_setComm(chebfi%X,spacecom)

    ! Create slice transposer using subcommunicator and allocate chebfi%X
    call xgTransposer_constructor(chebfi%xgTransposerX,chebfi%X,chebfi%xXColsRows,nspinor,&
        STATE_COLSROWS,TRANS_ALL2ALL,chebfi%comm_rows,spacecom,0,0,chebfi%me_g0_fft,&
        gpu_option=chebfi%gpu_option,gpu_thread_limit=chebfi%gpu_thread_limit,&
        custom_ncolsColsRows=.true.,nrowsLinalg_sub=nrowsLinalg_ptr)
        ! true to allow different bandpp (avoid pad)

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
        !call chebfi_applyLowpassFilter(slice,getAX_BX,getBm1X,nspinor)
        lambda_minus = slice%maxeig_global ! computed during compute_spectrum
        call polyfit_applyLowpassFilter(chebfi,getAX_BX,getBm1X,lambda_minus,nspinor)
    else
        !call chebfi_applyBandpassFilter(slice,getAX_BX,getBm1X,nspinor)

        call chebfi_applyBandpassFilter(chebfi,getAX_BX,getBm1X,glb,gub,lb,ub,nspinor)
    end if

    ! Transpose to linalg state
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
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

    ! Unitary test
    ABI_CHECK(rows(slice%X_linalg)==nrowsLinalg(xmpi_comm_rank(spacecom)),'wrong linalg representation')
    write(*,'(a,i6,i6)') '# proc has # rows of slice X ', xmpi_comm_rank(spacecom), rows(slice%X_linalg)

    ! Apply Rayleigh-Ritz for each MPI row
    ABI_NVTX_START_RANGE(NVTX_POLYFI_RR)
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

    
    ! Transpose to colsrows state (X only)
    call timab(tim_transpose,1,tsec)
    ABI_NVTX_START_RANGE(NVTX_POLYFI_TRANSPOSE)
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

    ! if gpu: This is important! Because X0 is on CPU
    ! call xgBlock_copy_from_gpu(slice%X)
    ! call xgBlock_set_gpu_option(slice%X,ABI_GPU_DISABLED)

    ! Copy slice solution to the extended buffer (requires colsrows state)
    call xgBlock_copy(slice%X,X0)
    ! FIXME same for eigen, residu?

    ! TODO we can also compute the merge on individual slices
    ! otherwise we do it outside?

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

end subroutine chebfi_run_slice
!!***

! ----------------------------------------------------------------------

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
    ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
    call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
    call xgBlock_zero_im_g0(chebfi%xAXColsRows)
    call xgBlock_zero_im_g0(chebfi%xBXColsRows)
    ABI_NVTX_END_RANGE()
    call timab(tim_getAX_BX,2,tsec)

    ABI_NVTX_START_RANGE(NVTX_POLYFI_CORE)
    do ideg = 0, ndeg_filter - 1

        ABI_NVTX_START_RANGE(NVTX_POLYFI_NEXT_ORDER)
        call chebfi_computeNextOrderChebfiPolynom(chebfi, ideg, center, one_over_r, two_over_r, getBm1X)
        ABI_NVTX_END_RANGE()

        ABI_NVTX_START_RANGE(NVTX_POLYFI_SWAP_BUF)
        if (chebfi%paral_kgb == 0) then
            call chebfi_swapInnerBuffers(chebfi, spacedim, neigenpairs)
        else
            call chebfi_swapInnerBuffers(chebfi, chebfi%total_spacedim, bandpp)
        end if
        ABI_NVTX_END_RANGE()

        !A * Psi
        call timab(tim_getAX_BX,1,tsec)
        ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
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

!!****f* m_polyfi/polyfi_bandpass
!! NAME
!! polyfi_bandpass
!!
!! FUNCTION
!! Apply Bandpass filter using Chebyshev-Jackson polynomial on a set of vectors.
!!
!! INPUTS
!!  polyfi  = polynomial filtering datastructure
!!  getAX_BX= pointer to the function giving A|X> and B|X>
!!            A is typically the Hamiltonian H, and B the overlap operator S
!!  getBm1X = pointer to the function giving B^-1|X>
!!            B is typically the overlap operator S
!!
!! SIDE EFFECTS
!!  polyfi <type(polyfi_t)>=all data used to apply Polynomial Filtering algorithm 
!!  polyfi%chebfi%xXColsRows= Filtered vectors
!!
!! SOURCE

subroutine polyfi_bandpass(polyfi,getAX_BX,getBm1X,nspinor)

    implicit none

    ! Arguments ------------------------------------
    type(polyfi_t), intent(inout) :: polyfi
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

    ! Memory buffers of chebfi datastructure
    chebfi = polyfi%chebfi

    ! Parameters
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
    mineig_global = polyfi%mineig_globa
    maxeig_global = polyfi%maxeig_global
    lambda_minus = polyfi%lambda_minus
    lambda_plus = polyfi%lambda_plus
   
    call xg_init(ChebyExpansion, chebfi%space, chebfi%, ncols, chebfi%spacecom, gpu_option=chebfi%gpu_option)
    call xgBlock_zero(ChebyExpansion%self)

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
        ABI_NVTX_START_RANGE(NVTX_POLYFI_NEXT_ORDER)
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
        ABI_NVTX_START_RANGE(NVTX_POLYFI_GET_AX_BX)
        call getAX_BX(chebfi%xXColsRows,chebfi%xAXColsRows,chebfi%xBXColsRows)
        call xgBlock_zero_im_g0(chebfi%xAXColsRows)
        call xgBlock_zero_im_g0(chebfi%xBXColsRows)
        ABI_NVTX_END_RANGE()
    
    end do ! end iline
    ABI_NVTX_END_RANGE()

    ! All slice processes wait for filter done
    ! FIXME Is this necessary? No communication happens actually
    if (chebfi%paral_kgb == 1) then
        call xmpi_barrier(chebfi%spacecom)
    end if

    ! Free Chebyshev expansion workspace
    call xg_free(ChebyExpansion)

end subroutine polyfi_bandpass
!!***



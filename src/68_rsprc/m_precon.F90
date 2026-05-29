!!****m* ABINIT/m_precon
!! NAME
!!  m_precon
!!
!! FUNCTION
!!  Object used for chi0-based preconditioning of the SCF.
!!
!! SOURCE

#include "abi_common.h"

module m_precon

    use iso_c_binding
    use defs_abitypes,          only : MPI_type
    use defs_basis
    use m_dtset
    use m_dtfil
    use m_xmpi

    use defs_datatypes,         only : pseudopotential_type
    use defs_wvltypes
    use m_atomdata,             only : atom_length
    use m_bandfft_kpt,          only : bandfft_kpt, bandfft_kpt_get_ikpt, bandfft_kpt_set_ikpt
    use m_cgprj,                only : ctocprj
    use m_cgtools
    use m_dfpt_mkvxc,           only : dfpt_mkvxc, dfpt_mkvxc_noncoll
    use m_fft,                  only : fourdp, fourwf, fftpac, zerosym
    use m_fftcore,              only : sphereboundary
    use m_fourier_interpol,     only : transgrid
    use m_iterative_solvers,    only : cg_linear_solver, gmres_linear_solver
    use m_kg,                   only : ph1d3d
    use m_mkrho
    use m_mpinfo,               only : proc_distrb_cycle, proc_distrb_band
    use m_occ,                  only : getnel
    use m_paw_dmft
    use m_pawrhoij,             only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free
    use m_pawcprj,              only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_mpi_allgather, pawcprj_free
    use m_pawang,               only : pawang_type
    use m_pawfgr,               only : pawfgr_type
    use m_pawtab,               only : pawtab_type
    use m_pawfgrtab,            only : pawfgrtab_type
    use m_paw_finegrid,         only : pawgylmg
    use m_paw_occupancies,      only : pawmkrhoij
    use m_paw_mkrho,            only : pawmkrho
    use m_paw_nhat,             only : pawsushat
    use m_prep_kgb,             only : prep_getghc, prep_index_wavef_bandpp, prep_fourwf
    use m_spacepar,             only : symrhg

    implicit none
    private

    type, public :: precon_object
        integer  :: iprcel
        real(dp) :: dielng, diemix
        !Geometry :
        real(dp) :: gprimd(3, 3), rprimd(3, 3), gmet(3, 3), rmet(3, 3)
        real(dp) :: ucvol, dvol
        !PAW :
        type(pseudopotential_type), pointer :: psps
        integer :: unpaw
        integer, pointer :: dimcprj(:), mcprj, usecprj
        type(pawcprj_type), pointer :: cprj(:, :)
        type(pawang_type), pointer :: pawang
        type(pawfgr_type), pointer :: pawfgr
        type(pawfgrtab_type), pointer :: pawfgrtab(:)
        type(pawtab_type), pointer :: pawtab(:)
        real(dp), pointer :: ylm(:, :)
        real(dp), pointer :: ylmgr(:, :, :)
        !To compute (weighted) densities and other quantities :
        real(dp), pointer :: fermie
        real(dp), pointer :: cg(:, :), eigen(:), occ(:), ph1d(:, :), phnons(:, :, :)
        integer, pointer  :: kg(:, :), npwarr(:), irrzon(:, :, :)
        integer, pointer :: atindx(:), atindx1(:), nattyp(:)
        integer, pointer :: symrec(:, :, :), indsym(:, :, :)
        real(dp), pointer :: xred(:, :)
        !For ffts :
        integer :: nfftprc      ! Number of fft grid points for preconditioned quantities (densities and/or potentials).
        integer :: ngfftprc(18) ! All needed information about the 3D FFT for preconditioned quantities.
        !For Kxc :
        integer :: nkxc
        real(dp), pointer :: kxc(:, :)
        real(dp), pointer :: rhor(:, :)
        real(dp), pointer :: vxc(:, :)

        !Logical variables :
        logical :: use_precon
        logical :: use_ldos
        logical :: use_paw_rhoij
        logical :: use_dos
        logical :: use_kxc
        logical :: use_ridgereg
        logical :: use_indices_arrays
        logical :: use_precomputed_rhoi
        logical :: use_precomputed_psii

        !For LDOS preconditioner :
        real(dp) :: tdos
        real(dp), allocatable :: ldos(:, :)
        real(dp), allocatable :: dos(:)

        !Usefull
        integer, allocatable :: cg_indices(:, :, :, :)
        integer, allocatable :: kg_indices(:, :)

        !Array to store precomputed ffts of cg :
        real(dp), allocatable :: precomputed_rhoi(:, :, :)
        integer, allocatable :: precomputed_rhoi_indices(:, :, :)
        real(dp), allocatable :: precomputed_psii(:, :, :, :)
        integer, allocatable :: precomputed_psii_indices(:, :, :)

        !Preconditioner parameters
        integer :: precon_verbose
        !Linear solver parameters
        integer :: linsolve_maxiter
        real(dp) :: linsolve_rtol, ridge_param
        !chi0_diag/quasidiag parameters
        real(dp) :: deigvals_tol_fp
        real(dp) :: precon_tsmear

    contains
        procedure :: init => precon_init                    ! Initialize the precon_object.
        procedure :: init_kxc => precon_init_kxc            ! Initialize kxc in the precon_object.
        procedure :: update => precon_update                ! Update the precon_object according to iprcel.
        procedure :: free => precon_free                    ! Dealocate arrays that are allocated in precon_init.

        procedure :: apply_dielmat => apply_dielmat         ! Apply the dielectric matrix to an input vector.
        procedure :: apply_adjdielmat => apply_adjdielmat   ! Apply the adjoint dielectric matrix to an input vector.
        procedure :: apply_precon => apply_precon           ! Apply the preconditioner to an input vector.
        
    end type precon_object

contains 

    ! TODO :  
    ! - debug quasidiag
    ! - non-collinear magnetism with band paral (Hyprid precon)
    ! - non coll : what spin representations are use for : kxc (dfpt_mkvxc_noncoll), prcref.
    ! - TODO : check places where it is assumed that nspinor=2 => nspden=4 and nsppol=2 => nspden=2
    !                                       on peut avoir nspden=1 dans les deux cas.
    ! - Linear solver : the tol is the absolute tol -> change to relative tol (in gmresm?)

    !****f* m_precon/precon_init
    !! NAME
    !!  precon_init
    !!
    !! FUNCTION
    !!  Initialize the precon_object.
    !!
    !! INPUTS
    !!  dtset   = all input variables for this dataset
    !!  atindx  = index table for atoms (see gstate.f)
    !!  atindx1 = index table for atoms, inverse of atindx (see gstate.f)
    !!  cg      = wf in G space
    !!  cprj    =
    !!  dimcprj = dimension of the cprj array
    !!  eigen   = array of eigenvalues
    !!  fermie  = fermi energie
    !!  gprimd  = dimensional reciprocal space primitive translations
    !!  irrzon  = irreducible zone data
    !!  kg      = reduced planewave coordinates
    !!  nattyp  = number of atoms of each type.
    !!  nfftmix = number of planewaves in the mixing/preconditioning grid
    !!  ngfftmix = 
    !!  npwarr  = number of planewaves and boundary planewaves at each k
    !!  pawang  = paw angular mesh and related data
    !!  pawfgr  = 
    !!  pawfgrtab = 
    !!  pawtab  = 
    !!  phnons  = nonsymmorphic translation phases
    !!  psps    = pseudopotential data
    !!  rprimd  = dimensional real space primitive translations
    !!  ucvol   = unit cell volume
    !!  xred    = reduced dimensionless atomic coordinates
    !!
    !! SOURCE
    subroutine precon_init(this, dtset, atindx, atindx1, cg, cprj, dimcprj, dtfil, eigen, fermie, &
        &   gmet, gprimd, indsym, irrzon, kg, mcprj, nattyp, nfftmix, ngfftmix, npwarr, occ, pawang, pawfgr, pawfgrtab, &
        &   pawtab, ph1d, phnons, psps, rhor, rmet, rprimd, symrec, ucvol, usecprj, vxc, xred, ylm)

        !Arguments ------------------------------------
        class(precon_object), intent(out) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        real(dp), intent(in) :: ucvol
        real(dp), intent(in), target :: fermie
        integer, intent(in) :: nfftmix
        integer, intent(in), target :: mcprj
        integer, intent(in), target :: usecprj
        !arrays
        real(dp), intent(in) :: gprimd(:, :), rprimd(:, :), gmet(:, :), rmet(:, :)
        integer, intent(in) :: ngfftmix(:)
        integer, intent(in), target  :: irrzon(:, :, :), kg(:, :), npwarr(:)
        integer, intent(in), target :: atindx(:), atindx1(:), nattyp(:)
        integer, intent(in), target :: symrec(:, :, :), indsym(:, :, :)
        real(dp), intent(in), target :: cg(:, :), eigen(:), occ(:), phnons(:, :, :), ph1d(:, :)
        real(dp), intent(in), target :: rhor(:, :), vxc(:, :)
        real(dp), intent(in), target :: xred(:, :)
        type(datafiles_type),intent(in) :: dtfil
        type(pseudopotential_type), intent(in), target :: psps
        integer, intent(in), target :: dimcprj(:)
        type(pawcprj_type), intent(in), target :: cprj(:, :)
        type(pawang_type), intent(in), target :: pawang
        type(pawfgr_type), intent(in), target :: pawfgr
        type(pawfgrtab_type), intent(in), target :: pawfgrtab(:)
        type(pawtab_type), intent(in), target :: pawtab(:)
        real(dp), intent(in), target :: ylm(:, :)

        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'precon%init')

        this%iprcel = dtset%iprcel
        this%use_precon = .false.

        if (this%iprcel >= 200 .and. this%iprcel < 300) then

            this%use_precon = .true.
            
            !Logical variables that describe the preconditioner :
            !   200 -> LDOS with RPA
            !   201 -> LDOS without RPA
            !   202 -> Hybrid (LDOS+diag)
            !   203 -> Hybrid for antiferro (LDOS+quasidiag) - WIP - Not documented
            !   210 -> Kerker with dielng   - Not documented
            !   211 -> Kerker with DOS      - Not documented
            !   212 -> Kerker with DOS+ and DOS - Not documented
            !   299 -> No preconditioning   - Not documented

            ! this%use_ldos = .true. activates the computation of the ldos.
            this%use_ldos = .false.                     
            if (this%iprcel == 200) this%use_ldos = .true.
            if (this%iprcel == 201) this%use_ldos = .true.
            if (this%iprcel == 202) this%use_ldos = .true.
            if (this%iprcel == 211) this%use_ldos = .true.
            if (this%iprcel == 212) this%use_ldos = .true.
            
            this%use_paw_rhoij = .false.

            ! If this%use_dos = .true. we will use this%dos.
            this%use_dos = .false.                     
            if (this%iprcel == 211) this%use_dos = .true.
            if (this%iprcel == 212) this%use_dos = .true.
            
            ! this%use_kxc = .true. activates the use of the exchange and correlation kernel.
            ! If this%use_kxc = .false. the RPA will be used.
            this%use_kxc = .false.                      
            if (this%iprcel == 201) this%use_kxc = .true.
            if (this%iprcel == 202) this%use_kxc = .true.

            ! this%use_ridgereg = .true. activates the use of an adapted linear solver.
            this%use_ridgereg = .false.                 
            
            ! this%use_indices_arrays = .true. indicates that we will use the arrays this%cg_indices and this%kg_indices.
            this%use_indices_arrays = .false.
            if (this%iprcel == 202) this%use_indices_arrays = .true.
            if (this%iprcel == 203) this%use_indices_arrays = .true.

            this%use_precomputed_rhoi = .false.
            if (this%iprcel == 202 .and. dtset%precon_in_memory==1) this%use_precomputed_rhoi = .true.

            this%use_precomputed_psii = .false.
            if (this%iprcel == 203 .and. dtset%precon_in_memory==1) this%use_precomputed_psii = .true.
            
            ! Other than here, iprcel is only used in apply_chi0, apply_dielmat and apply_adjdielmat.

            !Constant data from dtset
            this%dielng = dtset%dielng
            this%diemix = dtset%diemix
            this%nfftprc = nfftmix              ! FFT grid for preconditioned densities and/or potentials :
            this%ngfftprc = ngfftmix            ! same grid as the one used for mixing.
            !Other constants
            this%dvol   = ucvol/this%nfftprc    ! factor for integrals in real space (on the preconditioning FFT grid) : sum(f) * dvol ~ integral f
            this%gprimd = gprimd
            this%rprimd = rprimd
            this%gmet   = gmet
            this%rmet   = rmet
            this%ucvol  = ucvol
            !Pointers
            this%atindx => atindx 
            this%atindx1 => atindx1 
            this%cg     => cg
            this%eigen  => eigen
            this%fermie => fermie
            this%indsym => indsym
            this%irrzon => irrzon
            this%kg     => kg
            this%nattyp => nattyp
            this%npwarr => npwarr
            this%occ    => occ
            this%ph1d   => ph1d
            this%phnons => phnons
            this%psps   => psps
            this%rhor   => rhor
            this%symrec => symrec
            this%vxc    => vxc
            this%xred   => xred
            
            !PAW :
            if (psps%usepaw==1) then
                this%unpaw      = dtfil%unpaw
                this%cprj       => cprj
                this%usecprj    => usecprj
                this%dimcprj    => dimcprj
                this%mcprj      => mcprj
                this%pawang     => pawang
                this%pawfgr     => pawfgr
                this%pawfgrtab  => pawfgrtab
                this%pawtab     => pawtab
                this%ylm        => ylm
            end if
            
            !Initializing LDOS specific variables
            if (this%use_ldos) then
                !Allocating the array containing ldos
                ABI_MALLOC(this%ldos, (this%nfftprc, dtset%nspden))
            end if

            if (this%use_dos) then
                !Allocating the array containing the dos
                ABI_MALLOC(this%dos, (dtset%nspden))
            end if
            
            !Initializing variables needed for Kxc
            if (this%use_kxc) then
                !Preparing the allocation of Kxc
                if (dtset%xclevel==1) then  !LDA
                    this%nkxc = 2*min(dtset%nspden,2)-1
                else if (dtset%xclevel==2)then  !GGA+...
                    if (dtset%nspden==1) then
                        this%nkxc = 7
                    else if (dtset%nspden==2) then
                        this%nkxc = 19
                    else
                        ABI_BUG("chi0-based preconditioner (iprcel=2**): kxc not implemented for nspden > 2 (non-coll magn) in GGA")
                        ! TODO : update this if it gets implemented.
                    end if
                end if
            end if

            !Preconditioner parameters
            this%precon_verbose = dtset%precon_verbose
            !Linear solver parameters
            this%linsolve_maxiter = dtset%precon_ls_maxite
            this%linsolve_rtol = dtset%precon_ls_rtol
                ! For inversion of non positive definite (adjointe) dielectric matrix
            this%ridge_param = 0.01
            
            !chi0_diag/quasidiag parameters
            this%deigvals_tol_fp = tol10
            this%precon_tsmear = dtset%tsmear
            if (dtset%precon_tsmear /= zero) then   ! dtset%precon_tsmear = 0.0 means we use dtset%tsmear 
                                                    ! (the smearing temperature is not adjusted for the preonditioner).
                this%precon_tsmear = dtset%precon_tsmear    ! Only used for chi0_diag and chi0_quasidiag for now.
            end if

            !Usefull : indices mapping arrays
            if (this%use_indices_arrays)then
                ABI_MALLOC(this%cg_indices, (2*dtset%nspinor, dtset%mband, dtset%nkpt, dtset%nsppol))
                ABI_MALLOC(this%kg_indices, (2, dtset%nkpt))
            end if
            if (this%use_precomputed_rhoi) then
                ABI_MALLOC(this%precomputed_rhoi_indices, (dtset%mband, dtset%nkpt, dtset%nsppol))
            end if
            if (this%use_precomputed_psii) then
                ABI_MALLOC(this%precomputed_psii_indices, (dtset%mband, dtset%nkpt, dtset%nsppol))
            end if
            
        end if
        
    end subroutine precon_init

    !****f* m_precon/precon_init_kxc
    !! NAME
    !!  precon_init_kxc
    !!
    !! FUNCTION
    !!  Initialize the exchange and correlation kernel (kxc) in the precon_object.
    !!  This is needed because kxc needs to be initialized at a specific time.
    !!
    !! INPUTS
    !!  kxc = exchange and correlation kernel.
    !!
    !! SOURCE
    subroutine precon_init_kxc(this, kxc)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        real(dp), intent(in), target :: kxc(:, :)

        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'precon%init_kxc')
        if (this%use_precon) then
            if (this%use_kxc) then
                this%kxc => kxc
            end if
        end if

    end subroutine precon_init_kxc

    !****f* m_precon/precon_update
    !! NAME
    !!  precon_update
    !!
    !! FUNCTION
    !!  Update the precon_object :
    !!      For preconditioners using the LDOS (iprcel = 200 or 201) : 
    !!          Compute the new ldos (local density of state) with current wavefunctions 
    !!          and the new tdos (total density of state = integral of ldos).
    !!
    !! INPUTS
    !!  dtset     = All input variables for this dataset.
    !!  mpi_enreg = Information about MPI parallelization.
    !!
    !! SOURCE
    subroutine precon_update(this, dtset, mpi_enreg)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type), intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        
        !Local variables-------------------------------
        integer :: ispden

        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'precon%update')

        if (this%use_precon) then
            
            ! Indices in cg array
            if (this%use_indices_arrays) then
                call compute_cg_indices(dtset, mpi_enreg, this%npwarr, this%cg_indices)
                call compute_kg_indices(dtset, mpi_enreg, this%npwarr, this%kg_indices)
            end if

            !LDOS
            if (this%use_ldos) then 
                !update ldos
                call compute_ldos(this, dtset, mpi_enreg, this%ldos)
                !update tdos
                this%tdos = sum(this%ldos(:, 1)) * this%dvol
                ! TODO : More options to control when the ldos is updated

                ! Extremely Inefficient way to compute the DOS ...
                if (this%use_dos) then
                    do ispden = 1, dtset%nspden
                        this%dos(ispden) = sum(this%ldos(:, ispden)) * this%dvol
                    end do
                end if

            end if

            !Diag/Quasidiag chi0
            if (this%use_precomputed_psii) then
                call precompute_psii(this, dtset, mpi_enreg)
                ! Will also precompute rhoi.
            elseif (this%use_precomputed_rhoi) then
                call precompute_rhoi(this, dtset, mpi_enreg)
            end if

        end if
    end subroutine precon_update

    !****f* m_precon/precon_free_update
    !! NAME
    !!  precon_update
    !!
    !! FUNCTION
    !!  Deallocates arrays that might have been allocated in 'precon_update'.
    !!
    !! INPUTS
    !!  dtset     = All input variables for this dataset.
    !!  mpi_enreg = Information about MPI parallelization.
    !!
    !! SOURCE
    subroutine precon_free_update(this)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        
        ! *************************************************************************
        if (this%use_precon) then
            
            if (this%use_precomputed_rhoi) then
                ABI_FREE(this%precomputed_rhoi)
            end if
            if (this%use_precomputed_psii) then
                ABI_FREE(this%precomputed_psii)
            end if

        end if
    end subroutine precon_free_update

    !****f* m_precon/precon_free
    !! NAME
    !!  precon_free
    !!
    !! FUNCTION
    !!  Dealocate arrays that are allocated in precon_init.
    !!
    !! SOURCE
    subroutine precon_free(this)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'precon%free')
        if (this%use_precon) then

            if (this%use_indices_arrays) then
                ABI_FREE(this%cg_indices)
                ABI_FREE(this%kg_indices)
            end if
            
            if (this%use_precomputed_rhoi) then
                ABI_FREE(this%precomputed_rhoi_indices)
            end if
            if (this%use_precomputed_psii) then
                ABI_FREE(this%precomputed_psii_indices)
            end if

            if (this%use_ldos) then
                !Deallocating the array containing ldos and tdos
                ABI_FREE(this%ldos)
            end if

            if (this%use_dos) then
                ABI_FREE(this%dos)
            end if

        end if
    end subroutine precon_free

    !****f* m_precon/compute_r
    !! NAME
    !!  compute_r
    !!
    !! FUNCTION
    !!  Computes the array of r-vectors (in REDUCED coordinates).
    !!
    !! INPUTS
    !!  ngfft   = All needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!
    !! OUTPUTS
    !!  r_vectors(3, :) = 3 coordinates of the r_vectors.
    !!
    !! TODO : this is not used, maybe delete it ?
    !!
    !! SOURCE
    subroutine compute_r(ngfft, r_vectors)

        !Arguments ------------------------------------
        real(dp), intent(out) :: r_vectors(:, :)
        integer, intent(in) :: ngfft(:)

        !Local variables-------------------------------
        integer :: n1, n2, n3, i1, i2, i3, i_r

        ! *************************************************************************
        
        n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
        do i3=1,n3
            do i2=1,n2
                do i1=1,n1
                    i_r = 1 + (i1-1) + (i2-1)*n1 + (i3-1)*n1*n2
                    r_vectors(1, i_r) = real(i1-1)/n1
                    r_vectors(2, i_r) = real(i2-1)/n2
                    r_vectors(3, i_r) = real(i3-1)/n3
                end do
            end do
        end do

    end subroutine compute_r

    !****f* m_precon/get_r_vector
    !! NAME
    !!  get_r_vector
    !!
    !! FUNCTION
    !!  Get the vector r (in REDUCED coordinates) of index ifft
    !!
    !! INPUTS
    !!  ifft    = Index of sought vector r.
    !!  ngfft   = All needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!
    !! OUTPUT
    !!  r(3)    = sought vector r
    !!
    !! SOURCE
    function get_r_vector(ifft, ngfft) result(r)
        
        !Arguments ------------------------------------
        integer, intent(in) :: ifft
        integer, intent(in) :: ngfft(:)
        
        !Local variables-------------------------------
        integer :: n1, n2, n3, i1, i2, i3
        
        !Returned variable-------------------------------
        real(dp) :: r(3)
                
        ! *************************************************************************
        
        n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
        i1 = modulo((ifft-1), n1) + 1
        i2 = modulo((ifft-1)/n1, n2) + 1
        i3 = ((ifft-1)/n1)/n2 + 1
        r(1) = real(i1-1)/n1
        r(2) = real(i2-1)/n2
        r(3) = real(i3-1)/n3

    end function get_r_vector

    !****f* m_precon/get_g_vector
    !! NAME
    !!  get_g_vector
    !!
    !! FUNCTION
    !!  Get the vector g (in REDUCED coordinates) of index ifft
    !!
    !! INPUTS
    !!  ifft    = Index of sought vector g.
    !!  ngfft   = All needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!
    !! OUTPUT
    !!  g(3)    = sought vector g
    !!
    !! SOURCE
    function get_g_vector(ifft, ngfft) result(g)
        
        !Arguments ------------------------------------
        integer, intent(in) :: ifft
        integer, intent(in) :: ngfft(:)

        !Local variables-------------------------------
        integer :: n1, n2, n3, i1, i2, i3
        
        !Returned variable-------------------------------
        integer(dp) :: g(3)
                
        ! *************************************************************************
        
        n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
        ! (ifft-1) = (i1-1) + (i2-1)*n1 + (i3-1)*n1*n2
        i1 = modulo(ifft-1, n1) +1
        i2 = modulo((ifft-1)/n1, n2) +1
        i3 = ((ifft-1)/n1)/n2 +1
        g(1) = i1 - (i1/(n1/2+2))*n1-1
        g(2) = i2 - (i2/(n2/2+2))*n2-1
        g(3) = i3 - (i3/(n3/2+2))*n3-1

    end function get_g_vector

    !!****f* ABINIT/to_pauli
    !! NAME
    !!  to_pauli
    !!
    !! FUNCTION
    !!  Basis change from the default spin-basis to the Pauli basis for potentials and densities
    !!  in the direct (r) space.
    !!
    !! INPUT/OUTPUT
    !!  opt                = 0 : v is a potential
    !!                       1 : v is a density
    !!  v(nfft, nspden) = On input : Potential/density in the default spin-basis.
    !!                    On output : Potential/density in the Pauli basis.
    !!
    !! SOURCE
    subroutine to_pauli(this, opt, v)
        class(precon_object), intent(in) :: this
        !Arguments ------------------------------------
        real(dp), intent(inout) ::  v(:, :)
        integer :: opt
        !Local variables-------------------------------
        integer :: nspden
        real(dp), allocatable :: temp(:)
        
        ! *************************************************************************
        nspden = size(v, 2)

        !sigma_0, ... , sigma_3 are the Pauli matrices.
        if (opt == 0) then      !v is a potential
            if (nspden == 2) then
                !On input v(:, 1) is the spin-up potential and v(:, 2) is the spin-down potential.
                !On output the entire potential is v(:, 1)*sigma_0 + v(:, 2)*sigma_3.
                v(:, 1) = 0.5_dp*(v(:, 1) + v(:, 2))
                v(:, 2) = v(:, 1) - v(:, 2)
            else if (nspden == 4) then
                !                                   v(:, 1)             |  v(:, 3) + i*v(:, 4)
                !On input the entire potential is   --------------------|----------------------
                !                                   v(:, 3) - i*v(:, 4) |  v(:, 2)
                ! (see dotprod_vn in m_cgtools)
                !On output the entire potential is 
                !   v(:, 1)*sigma_0 + v(:, 2)*sigma_1 + v(:, 3)*sigma_2  + v(:, 4)*sigma_3
                !       v(:, 1) + v(:, 4)   |   v(:, 2) - i*v(:, 3)
                !   =   --------------------|----------------------
                !       v(:, 2) + i*v(:, 3) |   v(:, 1) - v(:, 4)
                
                ABI_MALLOC(temp, (size(v, 1)))

                v(:, 1) = 0.5_dp*(v(:, 1) + v(:, 2))
                temp = v(:, 4)
                v(:, 4) = v(:, 1) - v(:, 2) ! = 0.5_dp*(v(:, 1)_old  v(:, 2))
                v(:, 2) = v(:, 3)
                v(:, 3) = -temp

                ABI_FREE(temp)
            end if

        else if (opt == 1) then !v is a density
            if (nspden == 2) then
                !On input v(:, 1) is the total density and v(:, 2) is the spin-up density.
                !On output v(:, 1) is the total density and v(:, 2) is the spin density.
                v(:, 2) = 2*v(:, 2) - v(:, 1)
            end if
            !If nspden=4, the density is already given in the Pauli basis.
        end if
    end subroutine to_pauli

    !!****f* ABINIT/from_pauli
    !! NAME
    !!  from_pauli
    !!
    !! FUNCTION
    !!  Basis change from the Pauli basis to the Abinit default spin-basis for potentials and densities
    !!  in the real space.
    !!
    !! INPUT/OUTPUT
    !!  opt             = 0 : v is a potential
    !!                    1 : v is a density
    !!  v(nfft, nspden) = On input : Potential/density in the Pauli basis
    !!                    On output : Potential/density in the default spin-basis.
    !!
    !! SOURCE
    subroutine from_pauli(this, opt, v)
        class(precon_object), intent(in) :: this
        !Arguments ------------------------------------
        real(dp), intent(inout) ::  v(:, :)
        integer :: opt
        !Local variables-------------------------------
        integer :: nspden
        real(dp), allocatable :: temp(:)
        
        ! *************************************************************************
        nspden = size(v, 2)

        !sigma_0, ... , sigma_3 are the Pauli matrices.
        if (opt == 0) then      !v is a potential
            if (nspden == 2) then
                !On input the entire potential is v(:, 1)*sigma_0 + v(:, 2)*sigma_3.
                !On output v(:, 1) is the spin-up potential and v(:, 2) is the spin-down potential.
                v(:, 1) = v(:, 1) + v(:, 2)
                v(:, 2) = v(:, 1) - 2*v(:, 2)
            else if (nspden == 4) then

                !On input the entire potential is 
                !   v(:, 1)*sigma_0 + v(:, 2)*sigma_1 + v(:, 3)*sigma_2  + v(:, 4)*sigma_3
                !       v(:, 1) + v(:, 4)   |   v(:, 2) - i*v(:, 3)
                !   =   --------------------|----------------------
                !       v(:, 2) + i*v(:, 3) |   v(:, 1) - v(:, 4)
                !                                   v(:, 1)             |  v(:, 3) + i*v(:, 4)
                !On output the entire potential is  --------------------|----------------------
                !                                   v(:, 3) - i*v(:, 4) |  v(:, 2)
                ! (see dotprod_vn in m_cgtools)

                ABI_MALLOC(temp, (size(v, 1)))

                temp = v(:, 2)
                v(:, 2) = v(:, 1) - v(:, 4)
                v(:, 1) = v(:, 1) + v(:, 4)
                v(:, 3) = temp
                v(:, 4) = -v(:, 3)

                ABI_FREE(temp)

            end if

        else if (opt == 1) then !v is a density
            if (nspden == 2) then
                !On input v(:, 1) is the total density and v(:, 2) is the spin density.
                !On output v(:, 1) is the total density and v(:, 2) is the spin-up density.
                v(:, 2) = 0.5_dp*(v(:, 1) + v(:, 2))
            end if
            !If nspden=4, the density is already given in the Pauli basis.
        end if
    end subroutine from_pauli

    !****f* m_precon/apply_vc
    !! NAME
    !!  apply_vc
    !!
    !! FUNCTION
    !!  Apply the Coulomb kernel vc to a vector (in place) in the Pauli basis.
    !!
    !! INPUTS
    !!
    !! SIDE EFFECTS
    !!  vec_r (nfftprc, nspden) = Vector (in direct space) to which the Coulomb kernel vc is applied (in place).
    !!                            When nspden > 1 vec_r is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_vc(this, dtset, mpi_enreg, vec_r)
        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        integer :: ifft, ispden, cplex
        real(dp) :: g_cart_2
        integer :: n1, n2, n3
        real(dp) :: vec_g(2, this%nfftprc, 1)
        
        ! *************************************************************************
        
        !In the sigma_0, 1, 2, 3 (pauli) basis :
        !   The sigma_0 component of the density is multiplied by 4pi/G^2
        !   and the rest is 0.

        cplex = 1   ! vec is REAL
        n1=this%ngfftprc(1) ; n2=this%ngfftprc(2) ; n3=this%ngfftprc(3)

        ! FFT
        call fourdp(cplex, vec_g, vec_r(:, 1), -1, mpi_enreg, this%nfftprc, 1, this%ngfftprc, 0)
       
        ispden = 1
        do ifft = 2, this%nfftprc
            g_cart_2 = norm2(two_pi * matmul(this%gprimd, get_g_vector(ifft, this%ngfftprc)))**2
            vec_g(:, ifft, ispden) = (2*two_pi/g_cart_2) * vec_g(:, ifft, ispden)
        end do

        ! Set contribution of unbalanced components to zero.
        call zerosym(vec_g(:, :, 1), 2, n1, n2, n3)

        ! iFFT
        call fourdp(cplex, vec_g, vec_r(:, 1), 1, mpi_enreg, this%nfftprc, 1, this%ngfftprc, 0)
        
        do ispden = 2, dtset%nspden
            vec_r(:, ispden) = 0
        end do

    end subroutine apply_vc

    !****f* m_precon/apply_kxc
    !! NAME
    !!  apply_kxc
    !!
    !! FUNCTION
    !!  Apply the exchange and correlation kernel Kxc to a vector (in real space).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Informations about MPI parallelization.
    !!  vec_r (nfftprc, nspden) = Vector (in real space) to which the exchange and correlation kernel Kxc is applied.
    !!                            When nspden > 1 vec_r is in the default Abinit spin-basis.
    !!
    !! OUTPUTS
    !!  Kxc_vec_r (nfftprc, nspden) = Resulting vector (in real space) containing the application of Kxc to vec_r.
    !!                            When nspden > 1 Kxc_vec_r is in the Pauli-basis.
    !!
    !! SOURCE
    subroutine apply_kxc(this, dtset, mpi_enreg, vec_r, Kxc_vec_r)
        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: vec_r(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: Kxc_vec_r(this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        !scalars
        integer :: cplex, n3xccc, nhatdim, nhat1dim, nhat1grdim, nkxc, option, optnc, usexcnhat
        logical :: non_magnetic_xc
        !arrays
        real(dp), allocatable :: vec_r_default(:, :)
        real(dp), allocatable :: nhat(:, :), nhat1(:, :), nhat1gr(:, :, :)
        real(dp) :: dummy_xccc3d1(0), qphon(3)
        
        ! *************************************************************************

        if (size(this%kxc, 1) /= this%nfftprc) then
            ABI_BUG("chi0-based preconditioner : size(kxc, 1) /= nfftprc")
        end if

        !Applying Kxc : 
        cplex = 1   ! Input vector is real in real (direct) space.
        non_magnetic_xc = .false.
        nkxc = size(this%kxc, 2)

        usexcnhat = 0                                                           !
        nhat1dim = 0                                                            ! 
        ABI_MALLOC(nhat1, (cplex*this%nfftprc, dtset%nspden*nhat1dim))          ! PAW (TODO?)
        nhat1grdim = 0                                                          !
        ABI_MALLOC(nhat1gr, (cplex*this%nfftprc, dtset%nspden, 3*nhat1grdim))   !

        option = 2  ! Treats only density change (no core_correction)
        n3xccc = 0  !   -> Core-correction set to 0.
        qphon = 0.0_dp ! phonon vector

        if (dtset%nspden==1) then
            call dfpt_mkvxc(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfftprc, this%ngfftprc, nhat1, nhat1dim, &
            &               nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden, n3xccc, option, &
            &               qphon, vec_r, this%rprimd, usexcnhat, Kxc_vec_r, dummy_xccc3d1)
        
        else if (dtset%nspden==2) then
            ! Basis change to the default Abinit spin-basis for densities in case 
            ABI_MALLOC(vec_r_default, ((this%nfftprc), dtset%nspden))
            vec_r_default = vec_r
            call from_pauli(this, 1, vec_r_default)
            call dfpt_mkvxc(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfftprc, this%ngfftprc, nhat1, nhat1dim, &
            &               nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden, n3xccc, option, &
            &               qphon, vec_r_default, this%rprimd, usexcnhat, Kxc_vec_r, dummy_xccc3d1)
            ABI_FREE(vec_r_default)
        
        else if (dtset%nspden==4) then
            nhatdim = 0
            ABI_MALLOC(nhat, (this%nfftprc, dtset%nspden*nhatdim))              !
            optnc = 1   ! Compute the whole 2x2 Vres matrix
            call dfpt_mkvxc_noncoll(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfftprc, this%ngfftprc, nhat, nhatdim, &
            &               nhat1, nhat1dim, nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden,      &
            &               n3xccc, optnc, option, qphon, this%rhor, vec_r, this%rprimd, usexcnhat,        &
            &               this%vxc, Kxc_vec_r, dummy_xccc3d1)
            ! TODO : add an option for the precon to only compute V_11 V_22
            ABI_FREE(nhat)
            ! TODO noncoll : check in what spin-basis Kxc_vec_r is returned and adapt 'to_pauli'.
        end if

        ABI_FREE(nhat1)
        ABI_FREE(nhat1gr)

        call to_pauli(this, 0, Kxc_vec_r)

    end subroutine apply_kxc

    !****f* m_precon/apply_kernel
    !! NAME
    !! apply_kernel
    !!
    !! FUNCTION
    !!  Apply the kernel vc or (vc + Kxc) depending in iprcel to a vector (in place).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Informations about MPI parallelization.
    !!
    !! SIDE EFFECTS
    !!  vec_r (nfftprc, nspden) = Vector (in direct space) to which the kernel is applied (in place).
    !!                            When nspden > 1 vec_r is in the default Abinit spin-basis.
    !!
    !! SOURCE
    subroutine apply_kernel(this, dtset, mpi_enreg, vec_r)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        real(dp), allocatable :: Kxc_vec_r(:, :)

        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_kernel')
        
        ! RPA : LDOS/Kerker model - only vc
        if (.not. this%use_kxc) then
            
            call apply_vc(this, dtset, mpi_enreg, vec_r)    ! Apply vc in place

        ! No RPA : vc and Kxc
        else
            
            ! Apply Kxc
            ABI_MALLOC(Kxc_vec_r, (this%nfftprc, dtset%nspden))
            call apply_Kxc(this, dtset, mpi_enreg, vec_r, Kxc_vec_r)
            
            ! Apply vc in place
            call apply_vc(this, dtset, mpi_enreg, vec_r)    ! Apply vc in place
            
            ! Add Kxc_vec_r to vec_r
            vec_r = vec_r + Kxc_vec_r
            ABI_FREE(Kxc_vec_r)

        end if

    end subroutine apply_kernel

    !****f* m_precon/precon_increased_tsmear
    !! NAME
    !!  precon_increased_tsmear
    !!
    !! FUNCTION
    !!
    !! INPUTS
    !!
    !! OUTPUT
    !!  increased_tsmear = 
    !!
    !! SOURCE
    function precon_increased_tsmear(dtset) result(increased_tsmear)

        !Arguments ------------------------------------
        !scalars
        type(dataset_type),intent(in) :: dtset
        
        !Local variables-------------------------------
        !scalars
        real(dp) :: kpt_density, kpt_density_ref, tsmear_max

        !Returned variable-------------------------------
        real(dp) :: increased_tsmear
        
        ! *************************************************************************

        tsmear_max = 0.1
        kpt_density_ref = 6e5
        kpt_density = 6e5 ! TODO
        increased_tsmear = max(min(kpt_density/kpt_density_ref*dtset%tsmear, tsmear_max), dtset%tsmear)
        
    end function precon_increased_tsmear

    !****f* m_precon/derivative_occ
    !! NAME
    !!  derivative_occ
    !!
    !! FUNCTION
    !!  Compute the derivative of the occupation (f) of the band corresponding to eigenval 
    !!  with respect to the fermie temperature.
    !!      f = integral_((eigenval - fermie)/tsmear)^infty delta(t) dt
    !!      f' = -1/tsmear * delta((eigenval - fermie)/tsmear)
    !!
    !! INPUTS
    !!  occopt   = option for occupancies, determines delta
    !!  eigenval = eigenvalue
    !!  fermie   = fermi energie
    !!  tsmear   = smearing temperature
    !!
    !! OUTPUT
    !!  fprim    = occupation derivative
    !!
    !! SOURCE
    function derivative_occ(occopt, eigenval, fermie, tsmear) result(fprim)

        !Arguments ------------------------------------
        !scalars
        real(dp), intent(in) :: eigenval, fermie, tsmear
        integer, intent(in) :: occopt
        
        !Local variables-------------------------------
        !scalars
        real(dp) :: x, delta, a
        
        !Returned variable-------------------------------
        real(dp) :: fprim
        
        ! *************************************************************************
        
        x = (eigenval - fermie)/tsmear
    
        if (occopt<=2) then
            ABI_BUG("chi0-based preconditioner : Non-metallic occupation")
        else if (occopt==3) then
        !Fermi-Dirac smearing
            delta = exp(-abs(x))/(1+exp(-abs(x)))**2    !To avoid overflow of exp.
        else if (occopt==4) then
        !Cold Smearing
            a = -0.5634
            delta = (1.5+x*(-1.5*a+x*(-1.0+a*x)))*exp(-x**2)/sqrt(pi)
        else if (occopt==5) then
        !Cold Smearing
            a = -0.8165
            delta = (1.5+x*(-1.5*a+x*(-1.0+a*x)))*exp(-x**2)/sqrt(pi)
        else if (occopt==6) then
        !Smering of Methfessel and Paxton
            a = 0.0
            delta = (1.5+x*(-1.5*a+x*(-1.0+a*x)))*exp(-x**2)/sqrt(pi)
        else if (occopt==7) then
        !Gaussian smearing
            delta = exp(-x**2)/sqrt(pi)
        else if (occopt==8) then
        !Uniform smearing
            ABI_BUG("iprcel=2** : preconditioner needs a smooth smearing function")
        else if (occopt==9) then
        !Fermi-Dirac occupation is enforced with two distinct quasi-Fermi levels
            ABI_BUG("iprcel=2** : preconditioner not implemented for this smearing function")
        end if
    
        fprim = -1/tsmear * delta
        
    end function derivative_occ

    !****f* m_precon/compute_weighted_density
    !! NAME
    !!  compute_weighted_density
    !!
    !! FUNCTION
    !!  Wrapper for mkrho, symrhg and PAW : 
    !!  Compute a density-like quantity where the occupation are replaced by some weights
    !!      w_rho = sum_i weight_i |psi_i|^2 .
    !!  w_rho that has the same size as the preconditioned density/potential.
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  weights     = Weights that replace the occupations in the computation of density.
    !!
    !! OUTPUT
    !!  w_rhor      = "Weighted density" in real space, in the Pauli basis.
    !!
    !! SOURCE
    subroutine compute_weighted_density(this, dtset, mpi_enreg, weights, w_rhor)

        !Arguments ------------------------------------
        !scalars
        class(precon_object), intent(in) :: this
        type(dataset_type), intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: weights(:)
        real(dp), intent(out) :: w_rhor(this%nfftprc, dtset%nspden)

        !Local variables-------------------------------
        !scalars
        type(pawrhoij_type) :: pawrhoij(mpi_enreg%my_natom*this%psps%usepaw)
        integer :: dummy_int, mband_cprj, my_nspinor, mcprj_tmp
        integer :: mcg, cplex, cplex_rhoij
        real(dp) :: compch_fft
        integer :: optin, optout, optgrid
        !arrays
        real(dp) :: qphon(3)
        real(dp), allocatable :: w_rhowfg(:, :), w_rhowfr(:, :)
        type(pawcprj_type), allocatable :: cprj_tmp(:,:)
        type(paw_dmft_type)     :: dummy_paw_dmft
        type(wvl_wf_type)       :: dummy_wvl_wfs
        type(wvl_denspot_type)  :: dummy_wvl_den
        real(dp)                :: dummy_ylmgr(0, 0, 0)
        real(dp), allocatable   :: dummy_rhog(:, :), dummy_rhogf(:, :)

        ! *************************************************************************

        if (this%psps%usepaw==0) then
            ABI_MALLOC(w_rhowfr, (dtset%nfft, dtset%nspden))
            ABI_MALLOC(w_rhowfg, (2, dtset%nfft))
        else
            ABI_MALLOC(w_rhowfr, (this%pawfgr%nfftc, dtset%nspden))
            ABI_MALLOC(w_rhowfg, (2, this%pawfgr%nfftc))
        end if
        !w_rhowfr = zero
        !w_rhowfg = zero

        ! Compute the weighted density (w_rhor) using mkrho with weights in place of the occupations.
        mcg = size(this%cg)
        dummy_paw_dmft%use_dmft = 0
        dummy_paw_dmft%use_sc_dmft = 0
        call mkrho(this%cg, dtset, this%gprimd, this%irrzon, this%kg, mcg, mpi_enreg, this%npwarr, weights, &
        &   dummy_paw_dmft, this%phnons, w_rhowfg, w_rhowfr, this%rprimd, 0, this%ucvol, dummy_wvl_den, dummy_wvl_wfs, option=0, printout=.false.)
        ! symrhg already called in mkrho

        if (this%psps%usepaw==0) then
        ! In NC : the weighted density is directly w_rhowfr.
            if (this%nfftprc == dtset%nfft) then
                w_rhor = w_rhowfr
            else
                ABI_BUG("iprcel=2** : nfftprc /= nfft in Norm-conserving not implemented.")
            end if
        else
        ! In PAW : 
            ! First option : Add rhoij terms to w_rhowfr and transfer to fine grid. -UNUSED-
            if (this%use_paw_rhoij .and. this%nfftprc == this%pawfgr%nfft) then
                
                !Compute the rhoij equivalent for the weighted density.
                !   Sum_{n,k} {weight(n,k)*<Cnk|p_i><p_j|Cnk>}.
                
                my_nspinor = max(1, dtset%nspinor/mpi_enreg%nproc_spinor)
                mband_cprj = dtset%mband / mpi_enreg%nproc_band
                
                !Initialize pawrhoij
                cplex_rhoij = 1
                call pawrhoij_alloc(pawrhoij, cplex_rhoij, dtset%nspden, dtset%nspinor, &
                &       dtset%nsppol, dtset%typat, pawtab=this%pawtab)
                
                !Compute pawrhoij
                if (this%usecprj == 1) then ! cprj is saved in memory
                    call pawmkrhoij(this%atindx, this%atindx1, this%cprj, this%dimcprj, dtset%istwfk, dtset%kptopt, dtset%mband,&
                    &       mband_cprj, this%mcprj, dtset%mkmem, mpi_enreg, dtset%natom, dtset%nband, dtset%nkpt, dtset%nspden, &
                    &       dtset%nspinor, dtset%nsppol, weights, dtset%paral_kgb, dummy_paw_dmft, pawrhoij, this%unpaw,        &
                    &       dtset%usewvl, dtset%wtk)
                else                        ! cprj is computed on the fly
                    mcprj_tmp = my_nspinor * mband_cprj * dtset%mkmem * dtset%nsppol
                    ABI_MALLOC(cprj_tmp, (dtset%natom, mcprj_tmp))
                    call pawcprj_alloc(cprj_tmp, 0, this%dimcprj)
                    call ctocprj(this%atindx, this%cg, 1, cprj_tmp, this%gmet, this%gprimd, 0, 0, 0, dtset%istwfk, this%kg,     &
                    &       dtset%kptns, mcg, mcprj_tmp, dtset%mgfft, dtset%mkmem, mpi_enreg, this%psps%mpsang, dtset%mpw,      &
                    &       dtset%natom, this%nattyp, dtset%nband, dtset%natom, dtset%ngfft, dtset%nkpt, dtset%nloalg,          &
                    &       this%npwarr, dtset%nspinor, dtset%nsppol, dtset%nsppol, dtset%ntypat, dtset%paral_kgb, this%ph1d,   &
                    &       this%psps, this%rmet, dtset%typat, this%ucvol, this%unpaw, this%xred, this%ylm, dummy_ylmgr)
                    call pawmkrhoij(this%atindx, this%atindx1, cprj_tmp, this%dimcprj, dtset%istwfk, dtset%kptopt,              &
                    &       dtset%mband, mband_cprj, mcprj_tmp, dtset%mkmem, mpi_enreg, dtset%natom, dtset%nband, dtset%nkpt,   &
                    &       dtset%nspden, dtset%nspinor, dtset%nsppol, weights, dtset%paral_kgb, dummy_paw_dmft, pawrhoij,      &
                    &       this%unpaw, dtset%usewvl, dtset%wtk)
                    call pawcprj_free(cprj_tmp)
                    ABI_FREE(cprj_tmp)
                end if

                !Compute the total weighted density (adding PAW-correction).
                cplex = 1
                dummy_int=0
                qphon = 0
                call pawmkrho(1, compch_fft, cplex, this%gprimd, dummy_int, this%indsym, dummy_int, mpi_enreg,                  &
                &       mpi_enreg%my_natom, dtset%natom, dtset%nspden, dtset%nsym, dtset%ntypat, dtset%paral_kgb, this%pawang,  & 
                &       this%pawfgr, this%pawfgrtab, dtset%pawprtvol, pawrhoij, pawrhoij, this%pawtab, qphon, w_rhowfg,         &
                &       w_rhowfr, w_rhor, this%rprimd, dtset%symafm, this%symrec, dtset%typat, this%ucvol, dtset%usewvl, this%xred)

                call pawrhoij_free(pawrhoij)
            
            ! Second option : Transfer density from coarse to fine grid without rhoij corrections.
            elseif (this%nfftprc == this%pawfgr%nfft) then
                !Transfering the weighted density to the fine (PAW) grid, no rhoij correction added.
                cplex = 1
                optgrid = 1 ! coarse to fine
                optin = 0   ! real space
                optout = 0  !
                ABI_MALLOC(dummy_rhog, (2, this%pawfgr%nfftc))
                ABI_MALLOC(dummy_rhogf, (2, this%pawfgr%nfft))
                call transgrid(cplex, mpi_enreg, dtset%nspden, optgrid, optin, optout, dtset%paral_kgb, this%pawfgr, &
                &              dummy_rhog, dummy_rhogf, w_rhowfr, w_rhor)
                ABI_FREE(dummy_rhog)
                ABI_FREE(dummy_rhogf)
                !call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, this%nfftprc, ?nfftot, dtset%ngfft, 1, dtset%nsppol, dtset%nsym, &
                !this%phnons, rhog, rhor, this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)
            else
                ABI_BUG("iprcel=2** : nfftprc /= pawfgr%nfft in PAW not implemented.")
            end if

        end if

        !With collinear spins the weighted density is not returned in the Pauli (tot/spin) basis by mkrho.
        if (dtset%nspden == 2) then
            !spin = 2up - tot
            w_rhor(:, 2) = 2*w_rhor(:, 2) - w_rhor(:, 1)
        end if

        ABI_FREE(w_rhowfr)
        ABI_FREE(w_rhowfg)
    
    end subroutine compute_weighted_density

    !****f* m_precon/compute_ldos
    !! NAME
    !!  compute_ldos
    !!
    !! FUNCTION
    !!  Compute the local density of states defined as
    !!      ldos = sum_nk f'_nk |u_nk|^2 .
    !!  where f'_nk is the derivative of the occupation (nk) with respect to the fermi energie.
    !!  When 'nspden'>1, the ldos is returned in the Pauli-basis.
    !!
    !! INPUTS
    !!  dtset       = all input variables for this dataset
    !!  mpi_enreg   = informations about MPI parallelization
    !!
    !! OUTPUT
    !!  ldos        = local density of state
    !!
    !! SOURCE
    subroutine compute_ldos(this, dtset, mpi_enreg, ldos)

        !Arguments ------------------------------------
        !scalars
        !scalars
        class(precon_object), intent(in) :: this
        type(dataset_type), intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(out) :: ldos(:, :)

        !Local variables-------------------------------
        !scalars
        integer :: maxocc, i_eigen
        !integer :: ikpt, iband, isppol, nband_k, i_eigen
        !arrays
        real(dp), allocatable :: ldos_weights(:)

        ! *************************************************************************

        !compute weights
        ABI_MALLOC(ldos_weights, (dtset%mband*dtset%nkpt*dtset%nsppol))
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        ldos_weights = 0

        do i_eigen = 1, dtset%mband*dtset%nkpt*dtset%nsppol
            ldos_weights(i_eigen) = -derivative_occ(dtset%occopt, this%eigen(i_eigen), this%fermie, dtset%tsmear) * maxocc
        end do

        !Compute ldos using mkrho with ldos_weights in place of the occupations
        call compute_weighted_density(this, dtset, mpi_enreg, ldos_weights, ldos)
        ABI_FREE(ldos_weights)

    end subroutine compute_ldos

    !****f* m_precon/apply_chi0_dfermie
    !! NAME
    !!  apply_chi0_dfermie
    !!
    !! FUNCTION
    !!  
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!
    !! SIDE EFFECTS
    !!  vec_r (nfftprc, nspden) = Vector (in real space) to which the model chi0 operator is applied (in place).
    !!                            When nspden > 1 vec_r is in the Pauli spin-basis.
    !!
    !! SOURCE
    subroutine apply_chi0_dfermie(this, dtset, vec_r)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        integer :: ispden, jspden
        real(dp) :: delta_fermie

        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_chi0_dfermie')

        ! Precompute the dot product between the ldos and vec for each spin coordinate
        delta_fermie = zero
        do jspden = 1, dtset%nspden
            delta_fermie = delta_fermie + 1/this%tdos * dot_product(this%ldos(:, jspden), vec_r(:, jspden)) * this%dvol
        end do
        do ispden = 1, dtset%nspden
                vec_r(:, ispden) = delta_fermie * this%ldos(:, ispden)
        end do

    end subroutine apply_chi0_dfermie

    !****f* m_precon/apply_chi0_ldos
    !! NAME
    !!  apply_chi0_ldos
    !!
    !! FUNCTION
    !!  Apply the ldos model chi0 operator to the vector vec_r (in place) in the Pauli basis.
    !!
    !! INPUTS
    !!  mpi_enreg    = Information about MPI parallelization.
    !!
    !! SIDE EFFECTS
    !!  vec_r (nfftprc, nspden) = Vector (in direct space) to which the model chi0 operator is applied (in place).
    !!                            When nspden > 1 vec_r is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0_ldos(this, dtset, mpi_enreg, vec_r)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        type(dataset_type),intent(in) :: dtset
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: ispden
        !arrays
        real(dp), allocatable :: work_r(:, :)
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_chi0_ldos')
       
        if (abs(this%tdos) > epsilon(this%tdos)) then   !Checking that tdos is not 0.
            ABI_MALLOC(work_r, (this%nfftprc, dtset%nspden))
            
            !1) chi0(v)(r)_1 = -sum_ispden ldos_ispden(r)*v_ispden(r).
            work_r(:, 1) = 0
            do ispden = 1, dtset%nspden
                work_r(:, 1) = work_r(:, 1) - this%ldos(:, ispden)*vec_r(:, ispden)
            end do
            !2) chi0(v)(r)_ispden = -ldos_ispden(r)*v_1(r) for ispden>1.
            do ispden = 2, dtset%nspden
                work_r(:, ispden) = -this%ldos(:, ispden)*vec_r(:, 1)
            end do
            !3) Apply the part comming from the variations of the Fermi-level.
            call apply_chi0_dfermie(this, dtset, vec_r)
            vec_r = work_r + vec_r
            
            ABI_FREE(work_r)
        else 
            vec_r = 0
        end if

    end subroutine apply_chi0_ldos

    !!****f* ABINIT/apply_adjdielmat_ldos
    !! NAME
    !!  apply_adjdielmat
    !!
    !! FUNCTION
    !!  Apply the ldos- adjoint dielectric matrix I-chi0_ldos*vc to the density rho_r (given in the direct space).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  rho_r       = Density vector (in direct space).
    !!
    !! OUTPUT
    !!  adjdielmat_rho_r = adjdielmat * rho_r
    !!
    !! NOTES
    !!
    !! SOURCE
    subroutine apply_adjdielmat_ldos(this, dtset, mpi_enreg, rho_r, adjdielmat_rho_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: rho_r(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: adjdielmat_rho_r(this%nfftprc, dtset%nspden)
        
        ! *************************************************************************

        adjdielmat_rho_r = rho_r
        !1) Apply vc (in the Pauli basis)
        call apply_vc(this, dtset, mpi_enreg, adjdielmat_rho_r)
        !2) Apply chi0_ldos (in the Pauli basis)
        call apply_chi0_ldos(this, dtset, mpi_enreg, adjdielmat_rho_r)
        !3) adjdielmat_rho_r = rho_r - vc * chi0 * rho_r = adjdielmat * rho_r
        adjdielmat_rho_r = rho_r - adjdielmat_rho_r

    end subroutine apply_adjdielmat_ldos

    !!****f* ABINIT/apply_dielmat_ldos
    !! NAME
    !!  apply_adjdielmat
    !!
    !! FUNCTION
    !!  Apply the ldos- dielectric matrix I-vc*chi0_ldos to the potential v_r (given in the direct space).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  v_r       = Density vector (in direct space).
    !!
    !! OUTPUT
    !!  dielmat_v_r = dielmat * v_r
    !!
    !! NOTES
    !!
    !! SOURCE
    subroutine apply_dielmat_ldos(this, dtset, mpi_enreg, v_r, dielmat_v_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: v_r(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: dielmat_v_r(this%nfftprc, dtset%nspden)
        
        ! *************************************************************************

        dielmat_v_r = v_r
        !1) Apply chi0_ldos (in the Pauli basis)
        call apply_chi0_ldos(this, dtset, mpi_enreg, dielmat_v_r)
        !2) Apply vc (in the Pauli basis)
        call apply_vc(this, dtset, mpi_enreg, dielmat_v_r)
        !3) dielmat_v_r = v_r - vc * chi0 * v_r = dielmat * v_r
        dielmat_v_r = v_r - dielmat_v_r

    end subroutine apply_dielmat_ldos

    ! TODO : document 3 next subroutine

    ! Return the index of the eigenvalue corresponding to iband, ikpt, isppol in the 'eigen' array
    function get_eigen_index(dtset, iband, ikpt, isppol) result(i_eigen)
        
        !Arguments ------------------------------------
        !scalars
        type(dataset_type),intent(in) :: dtset
        integer, intent(in) :: iband, ikpt, isppol
        
        !Returned variable-------------------------------
        integer :: i_eigen

        ! *************************************************************************
        
        i_eigen = iband + (ikpt-1)*dtset%mband + (isppol-1)*dtset%mband*dtset%nkpt
 
    end function get_eigen_index

    subroutine compute_cg_indices(dtset, mpi_enreg, npwarr, cg_indices)
        
        !Arguments ------------------------------------
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: npwarr(:)
        integer :: cg_indices(2*dtset%nspinor, dtset%mband, dtset%nkpt, dtset%nsppol)
            ! If nspinor=2 : cg_indices(1, iband, ikpt, isppol):cg_indices(2, iband, ikpt, isppol) is the range of the spin up
            !                cg_indices(3, iband, ikpt, isppol):cg_indices(4, iband, ikpt, isppol) is the range of the spin down
            !                for the band of indices (iband, ikpt, isppol).

        !Local variables ------------------------------
        integer :: iband, ikpt, isppol, i_cg
        
        ! *************************************************************************

        cg_indices = zero
        i_cg = 1
        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt
                do iband = 1, dtset%nband(ikpt)
                    if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, iband, iband, isppol, mpi_enreg%me_kpt)) then
                        cycle
                    end if
                    cg_indices(1, iband, ikpt, isppol) = i_cg
                    cg_indices(2, iband, ikpt, isppol) = cg_indices(1, iband, ikpt, isppol) + npwarr(ikpt) - 1
                    i_cg = cg_indices(2, iband, ikpt, isppol) + 1
                    if (dtset%nspinor==2) then
                        cg_indices(3, iband, ikpt, isppol) = i_cg
                        cg_indices(4, iband, ikpt, isppol) = cg_indices(3, iband, ikpt, isppol) + npwarr(ikpt) - 1
                        i_cg = cg_indices(4, iband, ikpt, isppol) + 1
                end if
                end do
            end do
        end do

    end subroutine compute_cg_indices

    subroutine compute_kg_indices(dtset, mpi_enreg, npwarr, kg_indices)
        
        !Arguments ------------------------------------
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: npwarr(:)
        integer, intent(inout) :: kg_indices(2, dtset%nkpt)
        !Returned variable ----------------------------
        integer :: i_kg, ikpt, isppol

        ! *************************************************************************

        kg_indices = zero
        do isppol =1, dtset%nsppol
            i_kg = 1
            do ikpt = 1, dtset%nkpt
                    if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, dtset%nband(ikpt+(isppol-1)*dtset%nkpt), isppol, mpi_enreg%me_kpt)) then
                        cycle
                    end if
                    kg_indices(1, ikpt) = i_kg
                    kg_indices(2, ikpt) = kg_indices(1, ikpt) + npwarr(ikpt) - 1
                    i_kg = kg_indices(2, ikpt) + 1
            end do
        end do

    end subroutine compute_kg_indices

    ! Wrapper for fftpac and transgrid (PAW) :
    subroutine transfer_grid(this, dtset, mpi_enreg, ispden, rho_aug_r, rho_r)
        
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer :: ispden
        !arrays
        real(dp), intent(inout) :: rho_r(:, :)
        real(dp), intent(inout) :: rho_aug_r(:, :, :)

        !Local variables-------------------------------
        integer :: n1, n2, n3, n4, n5, n6
        integer :: cplex, optgrid, optin, optout
        real(dp), allocatable :: rho_coarse_r(:, :)
        real(dp), allocatable :: dummy_rhog(:, :), dummy_rhogf(:, :)

        ! *************************************************************************

        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)

        if (this%psps%usepaw==0) then
            ! In NC, the preconditioning grid should be the density/potential grid.
            if (this%nfftprc == n1*n2*n3) then
                call fftpac(ispden, mpi_enreg, 1, n1, n2, n3, n4, n5, n6, dtset%ngfft, rho_r, rho_aug_r, 1)
            else
                ABI_BUG("chi0-based preconditioner : nfftprc /= nfft in norm-conserving not implemented.")
            end if
        else
            ! In PAW, the preconditioning grid should be the fine grid.
            if (this%nfftprc == this%pawfgr%nfft) then
                ABI_MALLOC(rho_coarse_r, (dtset%nfft, 1))
                ! Augmented grid to coarse grid :
                call fftpac(1, mpi_enreg, 1, n1, n2, n3, n4, n5, n6, dtset%ngfft, rho_coarse_r, rho_aug_r, 1)
                ! Coarse grid to fine grid :
                cplex = 1
                optgrid = 1 ! coarse to fine
                optin = 0   ! real space
                optout = 0  !
                ABI_MALLOC(dummy_rhog, (2, this%pawfgr%nfftc))
                ABI_MALLOC(dummy_rhogf, (2, this%pawfgr%nfft))
                call transgrid(cplex, mpi_enreg, 1, optgrid, optin, optout, dtset%paral_kgb, this%pawfgr, dummy_rhog, dummy_rhogf, rho_coarse_r, rho_r(:, ispden))
                ABI_FREE(dummy_rhog)
                ABI_FREE(dummy_rhogf)
                ABI_FREE(rho_coarse_r)
            else
                ABI_BUG("chi0-based preconditioner : nfftprc /= pawfgr%nfft in PAW not implemented.")
            end if

        end if

    end subroutine transfer_grid

    ! Compute rho_i in collinear case without PAW corrections
    subroutine compute_rhoi_coll(this, dtset, mpi_enreg, iband, ikpt, isppol, rhoi_r)
        
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: isppol, ikpt, iband 
        real(dp), intent(inout) :: rhoi_r(this%nfftprc, 1)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex, optin, optout, optgrid
        integer :: ndat, option, tim_fourwf
        integer :: i_cg(2), i_kg(2)
        integer :: n1, n2, n3, n4, n5, n6
        integer :: istwf_k, npw_k
        !arrays
        integer :: gbound(2*dtset%mgfft+8,2)
        integer, allocatable :: kg_k(:, :)
        real(dp), allocatable :: rhoi_aug_r(:, :, :)
        !dummy arguments
        integer :: dummy_int
        real(dp) :: dummy_fofgout(2, 0), dummy_fofrout(2, dtset%ngfft(4), dtset%ngfft(5), dtset%ngfft(6))
        
        ! *************************************************************************

        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)
        
        ! TODO : check iband, ikpt, isppol belong to proc and return error if not
        !if () then
        !    ABI_BUG("chi0-based preconditioner (iprcel=2**): nfftprc /= nfft in norm-conserving not implemented.")
        !end if

        ! No spin or collinear spins - Wafefunctions have one spin component.
        if (dtset%nspinor == 1) then

            ! 1) Compute orbital density in real space (augmented basis) :
            
            ! Input parameters for fourwf :
            option = 1          ! Computes the density.
            ndat = 1            ! Only one FFT.
            tim_fourwf = 0
            ABI_MALLOC(rhoi_aug_r, (n4, n5, n6))
            rhoi_aug_r = zero   ! Initialization for fourwf (accumulation).
            istwf_k = dtset%istwfk(ikpt)    ! Option parameter that describes the storage of wfs at this kpt.

            if (.not. (dtset%paral_kgb == 1 .and. dtset%npband > 1)) then
                npw_k = this%npwarr(ikpt)               ! Number of plane-wave at this kpt.
                ABI_MALLOC(kg_k, (3, npw_k))
                i_kg = this%kg_indices(:, ikpt)
                kg_k = this%kg(:, i_kg(1):i_kg(2))      ! Reduced plane-wave coordinate (k+G) of this kpt.
                call sphereboundary(gbound, istwf_k, kg_k, dtset%mgfft, npw_k)
                i_cg = this%cg_indices(:, iband, ikpt, isppol)
                call fourwf(1, rhoi_aug_r, this%cg(:, i_cg(1):i_cg(2)), dummy_fofgout, dummy_fofrout,  &
                &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                &           dummy_int, n4, n5, n6, option, tim_fourwf, one, one)
            else
                ABI_BUG("'compute_rhoi_coll' should not be called with band parallelization.")
            end if
            ABI_FREE(kg_k)

            ! 2) Transfer rhoi_aug_r defined on the augmented (wavefunction) fft-grid to the preconditioning fft-grid.
            call transfer_grid(this, dtset, mpi_enreg, 1, rhoi_aug_r, rhoi_r)
            ABI_FREE(rhoi_aug_r)

            !4) Normalize rhoi_r.
            rhoi_r(:, 1) = rhoi_r(:, 1) / (sum(rhoi_r(:, 1)) * this%dvol) !Normalizing rho_ii_r.

        else
            ABI_BUG("chi0-based preconditioner : compute_rhoi_coll called with non-collinear magnetism")
        end if

    end subroutine compute_rhoi_coll

    ! Compute rho_i in non-collinear case
    subroutine compute_rhoi_noncoll(this, dtset,  mpi_enreg, iband, ikpt, isppol, rhoi_r)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: isppol, ikpt, iband 
        real(dp), intent(inout) :: rhoi_r(this%nfftprc, 1)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex, optin, optout, optgrid
        integer :: ndat, option, tim_fourwf, ier
        integer :: i_cg(4), i_kg(2)
        integer :: n1, n2, n3, n4, n5, n6
        integer :: istwf_k, npw_k
        integer :: ispden
        !arrays
        integer :: gbound(2*dtset%mgfft+8,2)
        integer, allocatable :: kg_k(:, :)
        real(dp), allocatable :: rhoi_aug_r(:, :, :, :)
        real(dp), allocatable :: rhoi_coarse_r(:, :)
        real(dp), allocatable :: psi_r_up(:, :, :, :), psi_r_down(:, :, :, :)
        !dummy arguments
        integer :: dummy_int
        real(dp) :: dummy_fofgout(2, 0), dummy_fofrout(2, dtset%ngfft(4), dtset%ngfft(5), dtset%ngfft(6))
        real(dp) :: dummy_denpot(dtset%ngfft(4), dtset%ngfft(5), dtset%ngfft(6))
        real(dp), allocatable :: dummy_rhog(:, :), dummy_rhogf(:, :)
        
        ! *************************************************************************
        
        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)
        
        ! Non collinear spins - Wafefunctions have two spin component.
        if (dtset%nspinor == 2) then

            !1) Compute psi_up and psi_down in real space :
            ABI_MALLOC(psi_r_up, (2, n4, n5, n6))
            ABI_MALLOC(psi_r_down, (2, n4, n5, n6))
            ! Input parameters for fourwf :
            option = 0          ! Only do the FFT.
            ndat = 1
            tim_fourwf = 0
            i_cg = this%cg_indices(:, iband, ikpt, isppol)
            istwf_k = dtset%istwfk(ikpt)    ! Option parameter that describes the storage of wfs at this kpt.
            !FFT for psi_up
            call fourwf(dummy_int, dummy_denpot, this%cg(:, i_cg(1):i_cg(2)), dummy_fofgout, psi_r_up,  &
            &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
            &           dummy_int, n4, n5, n6, option, tim_fourwf, one, one)
            !FFT for psi_down
            call fourwf(dummy_int, dummy_denpot, this%cg(:, i_cg(3):i_cg(4)), dummy_fofgout, psi_r_down,  &
            &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
            &           dummy_int, n4, n5, n6, option, tim_fourwf, one, one)
            !   (done separately for convenience & readability)
            
            !2) Compute the 4 components of the orbital density rhoi_aug_r :
            !   (in the Pauli basis for spins and real augmented basis for space)
            ABI_MALLOC(rhoi_aug_r, (n4, n5, n6, 4))
            ispden = 1  ! rho_sigma0(r) = |psi_up(r)|^2 + |psi_down(r)|^2
            rhoi_aug_r(:, :, :, ispden) = psi_r_up(1, :, :, :)**2 + psi_r_up(2, :, :, :)**2 + psi_r_down(1, :, :, :)**2 + psi_r_down(2, :, :, :)**2
            ispden = 2  ! rho_sigma1(r) = psi_up(r)* . psi_down(r) + psi_down(r)* . psi_up(r) = 2 Re(psi_up(r)*psi_down(r))
                        ! (* = conjugate)
            rhoi_aug_r(:, :, :, ispden) = 2*( psi_r_up(1, :, :, :)*psi_r_down(1, :, :, :) + psi_r_up(2, :, :, :)*psi_r_down(2, :, :, :) )
            ispden = 3  ! rho_sigma2(r) = i*(psi_down(r)* psi_up(r) - psi_up(r)* psi_down(r)) = 2 Im(psi_up(r)*psi_down(r))
            rhoi_aug_r(:, :, :, ispden) = 2*( psi_r_up(2, :, :, :)*psi_r_down(1, :, :, :) - psi_r_up(1, :, :, :)*psi_r_down(2, :, :, :) )
            ispden = 4  ! rho_sigma3(r) = |psi_up|^2 - |psi_down|^2
            rhoi_aug_r(:, :, :, ispden) = psi_r_up(1, :, :, :)**2 + psi_r_up(2, :, :, :)**2 - psi_r_down(1, :, :, :)**2 + psi_r_down(2, :, :, :)**2

            ABI_FREE(psi_r_up)
            ABI_FREE(psi_r_down)

            ! Change grid
            do ispden = 1, dtset%nspden
                call transfer_grid(this, dtset, mpi_enreg, ispden, rhoi_aug_r(:, :, :, ispden), rhoi_r)
            end do

            ABI_FREE(rhoi_aug_r)

        else
            ABI_BUG("chi0-based preconditioner : compute_rhoi_noncoll called with collinear magnetism")
        end if

    end subroutine compute_rhoi_noncoll

    subroutine precompute_psii(this, dtset, mpi_enreg)
         
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
       
        !Local variables-------------------------------
        !scalars
        integer :: n1, n2, n3, n4, n5, n6
        integer :: nspin, i_rhoi, isppol, ikpt, i_kpt_sppol, nband_k, rank, iband, iband1, iband2
        integer :: i_psii, option, ndat, istwf_k, npw_k, idat, ispinor, icplex, tim_fourwf
        integer :: i_kg(2), i_cg_iband1(2*dtset%nspinor), i_cg_iband2(2*dtset%nspinor)
        real(dp) :: sum_rhoi_r
        integer :: ifft
        logical :: band_paral
        !arrays
        integer, allocatable :: needed_bands_bounds(:, :)
        integer, allocatable :: needed_bands_number(:)
        integer, allocatable :: kg_k(:, :)
        integer :: gbound_k(2*dtset%mgfft+8,2)
        real(dp), allocatable :: psii_aug(:, :, :, :)
        real(dp), allocatable :: rhoi_aug(:, :, :)
        real(dp), allocatable :: delta_rho_coarse_r(:, :)
        !for band parall
        integer :: option_fourwf, blocksize, iblock, ibandblock1, ibandblock2, nbdblock, nfft_blocks
        integer :: i_cg_ibandblock1(2*dtset%nspinor), i_cg_ibandblock2(2*dtset%nspinor)
        real(dp), allocatable :: dummy_occ_k(:)
        !dummy arguments
        real(dp) ::  dummy_denpot(0, dtset%ngfft(5), dtset%ngfft(6)), dummy_fofgout(2, 0)
        
        ! *************************************************************************

        band_paral = (dtset%paral_kgb == 1 .and. dtset%npband > 1)

        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)

        ABI_MALLOC(needed_bands_bounds, (2, dtset%nkpt*dtset%nsppol))
        ABI_MALLOC(needed_bands_number, (dtset%nkpt*dtset%nsppol))
        call get_needed_bands_chi0diag(this, dtset, mpi_enreg, needed_bands_bounds, needed_bands_number)

        ! Allocate the array containing the precomputed psii
        ABI_MALLOC(this%precomputed_psii, (2, this%nfftprc, dtset%nspinor, sum(needed_bands_number)))
        this%precomputed_psii_indices = zero
        ABI_MALLOC(this%precomputed_rhoi, (this%nfftprc, dtset%nspinor, sum(needed_bands_number)))
        this%precomputed_rhoi_indices = zero
        i_psii = 1
        
        ABI_MALLOC(psii_aug, (2, n4, n5, n6*dtset%mband))
        ABI_MALLOC(rhoi_aug, (n4, n5, n6))

        !Loop over spins and kpoints
        do isppol = 1, dtset%nsppol
            do ikpt = 1, dtset%nkpt
                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt

                ! MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                nband_k = dtset%nband(i_kpt_sppol)
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if

                iband1 = needed_bands_bounds(1, i_kpt_sppol)
                iband2 = needed_bands_bounds(2, i_kpt_sppol)

                if (band_paral) then

                    if (dtset%nspinor==1) then
                        option_fourwf = 0
                        ndat = mpi_enreg%bandpp
                        
                        nbdblock=nband_k / (mpi_enreg%nproc_band * mpi_enreg%bandpp)
                        blocksize=nband_k / nbdblock
                        do iblock = 1, nbdblock     ! Loop over (LOBPCG) blocks

                            ibandblock1 = blocksize*(iblock-1) + 1  
                            ibandblock2 = blocksize*(iblock)
                            i_cg_ibandblock1 = this%cg_indices(:, ibandblock1, ikpt, isppol)
                            i_cg_ibandblock2 = this%cg_indices(:, ibandblock2, ikpt, isppol)    ! Changer

                            ABI_MALLOC(psii_aug, (2, n4, n5, n6*ndat))
                            ABI_MALLOC(dummy_occ_k, (nband_k))
                            !ABI_MALLOC(dummy_denpot, (n4, n5, n6))

                            call bandfft_kpt_set_ikpt(ikpt, mpi_enreg)
                            nfft_blocks = 1     ! TODO : what is that ??
                            call prep_fourwf(dummy_denpot, blocksize, this%cg(:, i_cg_ibandblock1(1):i_cg_ibandblock2(2)),    &
                            &           psii_aug(:, :, :, 1:n6*ndat), iblock, dtset%istwfk(ikpt), dtset%mgfft, mpi_enreg, nband_k,      &
                            &           ndat, dtset%ngfft, this%npwarr(ikpt),                                       &
                            &           n4, n5, n6, dummy_occ_k, option_fourwf, this%ucvol, dtset%wtk(ikpt), nfft_blocks)

                            ABI_FREE(dummy_occ_k)
                            !ABI_FREE(dummy_denpot)

                        end do
                    else
                        ABI_BUG("TODO non-collinear magnetisme in precon")
                    end if

                else
                    
                    if (dtset%nspinor==1) then
                        
                        option = 0
                        ndat = needed_bands_number(i_kpt_sppol)
                        istwf_k = dtset%istwfk(ikpt)        ! Option parameter that describes the storage of wfs at this kpt.
                        npw_k = this%npwarr(ikpt)           ! Number of plane-wave at this kpt.
                        ABI_MALLOC(kg_k, (3, npw_k))
                        i_kg = this%kg_indices(:, ikpt)
                        kg_k = this%kg(:, i_kg(1):i_kg(2))
                        call sphereboundary(gbound_k, istwf_k, kg_k, dtset%mgfft, npw_k)    ! Computes gbound.
                        tim_fourwf = 0
                        i_cg_iband1 = this%cg_indices(:, iband1, ikpt, isppol)
                        i_cg_iband2 = this%cg_indices(:, iband2, ikpt, isppol)
                        
                        call fourwf(0, dummy_denpot, this%cg(:, i_cg_iband1(1):i_cg_iband2(2)), dummy_fofgout,  &
                        &           psii_aug(:, :, :, 1:n6*ndat), gbound_k, gbound_k, istwf_k, kg_k, kg_k,      &
                        &           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, npw_k,                    &
                        &           n4, n5, n6, option, tim_fourwf, one, one)
                        
                        ABI_FREE(kg_k)

                    else
                        ABI_BUG("TODO non-collinear magnetisme in precon")
                    end if
                end if
                
                rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)

                do iband = iband1, iband2
                    
                    ! Check if this band belong to the current processor (assuming bands are distributed in order).
                    if (.not. (1 + mpi_enreg%bandpp*rank <= iband .and. iband <= mpi_enreg%bandpp*(rank+1))) then
                        cycle
                    end if

                    if (band_paral) then
                        idat = iband - mpi_enreg%bandpp*rank
                    else
                        idat = iband
                    end if
                        
                    rhoi_aug = zero
                    call cg_addtorho(n1, n2, n3, n4, n5, n6, 1, one, one, psii_aug(:, :, :, (idat-1)*n6+1:idat*n6), rhoi_aug)

                    ! Grid transfers:
                    ispinor = 1
                    do icplex = 1, 2
                        call transfer_grid(this, dtset, mpi_enreg, 1, psii_aug(icplex, :, :, idat:idat+n6-1), this%precomputed_psii(icplex, :, ispinor:ispinor, i_psii))
                        ! TODO : this will create a copy ... Maybe change dim order ?
                    end do
                    call transfer_grid(this, dtset, mpi_enreg, 1, rhoi_aug, this%precomputed_rhoi(:, :, i_psii))

                    ! Normalize psii:
                    sum_rhoi_r = 0
                    do ispinor = 1, dtset%nspinor
                        do ifft=1, this%nfftprc
                            sum_rhoi_r = sum_rhoi_r + this%precomputed_rhoi(ifft, ispinor, i_psii)
                        end do
                        this%precomputed_psii(:, :, ispinor, i_psii) = this%precomputed_psii(:, :, ispinor, i_psii) / sqrt(sum_rhoi_r * this%dvol)
                        this%precomputed_rhoi(:, ispinor, i_psii) = this%precomputed_rhoi(:, ispinor, i_psii) / (sum_rhoi_r * this%dvol)
                    end do

                    ! Save the index for iband, ikpt, isppol in precomputed_psii and precomputed_rhoi
                    this%precomputed_psii_indices(iband, ikpt, isppol) = i_psii
                    this%precomputed_rhoi_indices(iband, ikpt, isppol) = i_psii
                    i_psii = i_psii + 1

                end do  !iband

                
            end do  !ikpt
        end do  !isppol

        ABI_FREE(rhoi_aug)
        ABI_FREE(psii_aug)
        ABI_FREE(needed_bands_number)
        ABI_FREE(needed_bands_bounds)

    end subroutine precompute_psii

    subroutine precompute_rhoi(this, dtset, mpi_enreg)
         
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
       
        !Local variables-------------------------------
        !scalars
        integer :: nspin, i_rhoi, isppol, ikpt, i_kpt_sppol, nband_k, rank, iband, iband1, iband2
        !arrays
        integer, allocatable :: needed_bands_bounds(:, :)
        integer, allocatable :: needed_bands_number(:)
        !for band parall
        integer :: option_fourwf, ndat, blocksize, iblock, option, ibandblock1, ibandblock2, nbdblock, nfft_blocks
        integer :: n1, n2, n3, n4, n5, n6
        integer :: idat
        integer :: i_cg_ibandblock1(2*dtset%nspinor), i_cg_ibandblock2(2*dtset%nspinor)
        real(dp), allocatable :: dummy_occ_k(:)
        real(dp), allocatable :: dummy_denpot(:, :, :)
        real(dp), allocatable :: rhoi_aug(:, :, :)
        real(dp), allocatable :: psii_aug(:, :, :, :)
        
        ! *************************************************************************

        ABI_MALLOC(needed_bands_bounds, (2, dtset%nkpt*dtset%nsppol))
        ABI_MALLOC(needed_bands_number, (dtset%nkpt*dtset%nsppol))
        call get_needed_bands_chi0diag(this, dtset, mpi_enreg, needed_bands_bounds, needed_bands_number)

        ! Allocate the array containing the precomputed rhoi
        if (dtset%nspinor==1) then
            nspin = 1   ! Number of spin components in the orbital densities (rhoi).
        else if (dtset%nspinor==2) then
            nspin = 4
        else
            ABI_BUG("nspinor /= 1 or 2")
        end if
        ABI_MALLOC(this%precomputed_rhoi, (this%nfftprc, nspin, sum(needed_bands_number)))
        this%precomputed_rhoi_indices = zero
        i_rhoi = 1

        !Loop over spins and kpoints
        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt
                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt

                ! MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                nband_k = dtset%nband(i_kpt_sppol)
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if

                iband1 = needed_bands_bounds(1, i_kpt_sppol)
                iband2 = needed_bands_bounds(2, i_kpt_sppol)
                
                if (.not. (dtset%paral_kgb == 1 .and. dtset%npband > 1)) then  
                ! No parallelization over band : We loop over needed bands to compute and save the orbital densities.

                    do iband = iband1, iband2

                        ! No spin or collinear spins - Wafefunctions have one spin component.
                        if (dtset%nspinor == 1) then    
                            call compute_rhoi_coll(this, dtset, mpi_enreg, iband, ikpt, isppol, this%precomputed_rhoi(:, :, i_rhoi))
                            this%precomputed_rhoi_indices(iband, ikpt, isppol) = i_rhoi
                            i_rhoi = i_rhoi + 1
                        end if

                        ! Non collinear spins - Wavefunctions have two spins components.
                        if (dtset%nspinor == 2) then
                            call compute_rhoi_noncoll(this, dtset, mpi_enreg, iband, ikpt, isppol, this%precomputed_rhoi(:, :, i_rhoi))
                            this%precomputed_rhoi_indices(iband, ikpt, isppol) = i_rhoi
                            i_rhoi = i_rhoi + 1
                        end if
                        
                    end do  !iband

                else
                ! Parallelization over band : We first do all the needed fft, with 'prep_fourwf' that will take care of 
                ! the transpose from the linalg representation to the fft representation. 
                ! Then we fill the 'precomputed_rhoi' array band per band.
                    
                    n1 = dtset%ngfft(1)
                    n2 = dtset%ngfft(2)
                    n3 = dtset%ngfft(3)
                    n4 = dtset%ngfft(4)
                    n5 = dtset%ngfft(5)
                    n6 = dtset%ngfft(6)

                    if (dtset%nspinor==1) then

                        option_fourwf = 0
                        ndat = mpi_enreg%bandpp
                        
                        nbdblock=nband_k / (mpi_enreg%nproc_band * mpi_enreg%bandpp)
                        blocksize=nband_k / nbdblock
                        do iblock = 1, nbdblock     ! Loop over (LOBPCG) blocks

                            ibandblock1 = blocksize*(iblock-1) + 1  
                            ibandblock2 = blocksize*(iblock)
                            i_cg_ibandblock1 = this%cg_indices(:, ibandblock1, ikpt, isppol)
                            i_cg_ibandblock2 = this%cg_indices(:, ibandblock2, ikpt, isppol)

                            ABI_MALLOC(psii_aug, (2, n4, n5, n6*ndat))
                            ABI_MALLOC(dummy_occ_k, (nband_k))
                            ABI_MALLOC(dummy_denpot, (n4, n5, n6))

                            call bandfft_kpt_set_ikpt(ikpt, mpi_enreg)
                            nfft_blocks = 1     ! TODO : what is that ??
                            call prep_fourwf(dummy_denpot, blocksize, this%cg(:, i_cg_ibandblock1(1):i_cg_ibandblock2(2)),    &
                            &           psii_aug, iblock, dtset%istwfk(ikpt), dtset%mgfft, mpi_enreg, nband_k,      &
                            &           ndat, dtset%ngfft, this%npwarr(ikpt),                                       &
                            &           n4, n5, n6, dummy_occ_k, option_fourwf, this%ucvol, dtset%wtk(ikpt), nfft_blocks)

                            ABI_FREE(dummy_occ_k)
                            ABI_FREE(dummy_denpot)

                        end do

                    else
                        ABI_BUG("non-collinear magnetisme in precon TODO")
                        ! TODO : not copy the above but adapt it by multiplying ndat by nspinor i.e. remove if.
                    end if

                    idat = 0
                    ABI_MALLOC(rhoi_aug, (n4, n5, n6))

                    ! Loop over the bands of this proc (assuming bands are distributed in order).
                    rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)
                    do iband = 1 + mpi_enreg%bandpp*rank, mpi_enreg%bandpp*(rank+1)

                        idat = idat + 1

                        ! Check if this band is needed
                        if (.not. (iband1 <= iband .and. iband <= iband2)) then
                            cycle
                        end if

                        if (dtset%nspinor==1) then

                            !Compute rhoi_aug from psii_aug
                            rhoi_aug = zero
                            call cg_addtorho(n1, n2, n3, n4, n5, n6, 1, one, one, psii_aug(:, :, :, (idat-1)*n6+1:idat*n6), rhoi_aug)
                        
                            ! Grid transfer
                            call transfer_grid(this, dtset, mpi_enreg, 1, rhoi_aug, this%precomputed_rhoi(:, :, i_rhoi))
                            
                            ! Normalize
                            this%precomputed_rhoi(:, 1, i_rhoi) = this%precomputed_rhoi(:, 1, i_rhoi) / &
                            &                                     (sum(this%precomputed_rhoi(:, 1, i_rhoi)) * this%dvol)

                            ! Save the index
                            this%precomputed_rhoi_indices(iband, ikpt, isppol) = i_rhoi
                            i_rhoi = i_rhoi + 1
 
                        else
                            ABI_BUG("non-collinear magnetisme in precon TODO")
                        end if

                    end do
                    
                    ABI_FREE(rhoi_aug)
                    ABI_FREE(psii_aug)

                end if

            end do  !ikpt
        end do  !isppol

        ABI_FREE(needed_bands_number)
        ABI_FREE(needed_bands_bounds)

    end subroutine precompute_rhoi

    subroutine get_needed_bands_chi0diag(this, dtset, mpi_enreg, needed_bands_bounds, needed_bands_number)
        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        integer :: needed_bands_number(dtset%nsppol*dtset%nkpt)
        integer :: needed_bands_bounds(2, dtset%nsppol*dtset%nkpt)
       
        !Local variables-------------------------------
        !scalars
        integer :: i_eigen, i_kpt_sppol, ikpt, isppol, iband, nband_k, rank
        real(dp) :: fp, maxocc

        ! *************************************************************************
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        needed_bands_number = zero

        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt
                
                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt
                nband_k = dtset%nband(i_kpt_sppol)
                needed_bands_bounds(1, i_kpt_sppol) = nband_k + 1
                needed_bands_bounds(2, i_kpt_sppol) = 0

                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if
                
                rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)

                do iband = 1 + mpi_enreg%bandpp*rank, mpi_enreg%bandpp*(rank+1)
                    
                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)  ! Index of (iband, ikpt, isppol) in eigen array.
                    fp = derivative_occ(dtset%occopt, this%eigen(i_eigen), this%fermie, this%precon_tsmear) * maxocc

                    if (abs(fp) > this%deigvals_tol_fp) then
                        needed_bands_bounds(1, i_kpt_sppol) = min(needed_bands_bounds(1, i_kpt_sppol), iband)   ! iband_min
                        needed_bands_bounds(2, i_kpt_sppol) = max(needed_bands_bounds(2, i_kpt_sppol), iband)   ! iband_max
                    end if
                    
                end do

                needed_bands_number(i_kpt_sppol) = max(needed_bands_bounds(2, i_kpt_sppol) - needed_bands_bounds(1, i_kpt_sppol) + 1, 0)

            end do  !ikpt
        end do  !isppol

    end subroutine get_needed_bands_chi0diag

    subroutine compute_delta_occ(this, dtset, mpi_enreg, delta_V, delta_occ)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: delta_V(this%nfftprc, dtset%nspden)     ! In Pauli basis
        real(dp), intent(inout) :: delta_occ(:)
        
        !Local variables-------------------------------
        !scalars
        integer :: nband_k, nspin, rank
        integer :: i_eigen, ikpt, iband, isppol, ier, ispden
        integer :: iband1, iband2, i_kpt_sppol
        real(dp) :: fp, eigenval, maxocc
        real(dp) :: dos_fermie, delta_occ_tot, delta_fermie
        !arrays
        integer, allocatable :: needed_bands_bounds(:, :)
        integer, allocatable :: needed_bands_number(:)
        real(dp), allocatable :: rhoi_r(:, :)
        !real(dp), allocatable :: doccde(:), occ(:)

        ! *************************************************************************

        ABI_MALLOC(needed_bands_bounds, (2, dtset%nkpt*dtset%nsppol))
        ABI_MALLOC(needed_bands_number, (dtset%nkpt*dtset%nsppol))
        call get_needed_bands_chi0diag(this, dtset, mpi_enreg, needed_bands_bounds, needed_bands_number)

        ! Compute the delta_occ = fi' * <rhoii, vec> 
        delta_occ = zero
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        
        ! Allocate the arrays that will contain rhoi
        if (dtset%nspinor==1) then
            nspin = 1   ! Number of spin components in the orbital densities (rhoi).
        else if (dtset%nspinor==2) then
            nspin = 4
        else
            ABI_BUG("nspinor /= 1 or 2")
        end if
        ABI_MALLOC(rhoi_r, (this%nfftprc, nspin))
        dos_fermie = zero
        delta_occ_tot = zero

        ! Occupation derivatives could also be computetd with 'getnel' ... What is best ? TODO
            !ABI_MALLOC(doccde, (size(this%occ)))
            !ABI_MALLOC(occ, (size(this%occ)))
            !option=1
            !call getnel(doccde, dummy_real, this%eigen, entropy, this%fermie, this%fermie, maxocc, &
            !&       dtset%mband, dtset%nband, nelect, dtset%nkpt, dtset%nsppol, occ, dtset%occopt, &
            !&       option, dtset%tphysel, dtset%tsmear, dummy_int, dtset%wtk)
            !ABI_FREE(occ)

        ! 1) Eigenvalue variations
        !Loop over spins and kpoints
        do isppol =1, dtset%nsppol
            
            do ikpt = 1, dtset%nkpt
                
                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt
                nband_k = dtset%nband(i_kpt_sppol)

                ! MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if
                
                iband1 = needed_bands_bounds(1, i_kpt_sppol)
                iband2 = needed_bands_bounds(2, i_kpt_sppol)
                
                do iband = iband1, iband2   ! Loop over needed bands

                    rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)
                    if (.not.(1 + mpi_enreg%bandpp*rank <= iband .and. iband <=mpi_enreg%bandpp*(rank+1))) then
                        cycle
                    end if  ! TODO : create a cycle function, that is also correct in LOBPCG

                    !Indices
                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)  ! Index of (iband, ikpt, isppol) in eigen array.
                    
                    !2.1) Computing f'(eig_i - fermie).
                    eigenval = this%eigen(i_eigen)
                    fp = derivative_occ(dtset%occopt, eigenval, this%fermie, this%precon_tsmear) * maxocc
                    !fp = doccde(i_eigen)   ! Same as derivative_occ (probably more robust)
                    dos_fermie = dos_fermie + fp * dtset%wtk(ikpt)

                    ! No spin or collinear spins - Wafefunctions have one spin component.
                    if (dtset%nspinor == 1) then    

                        !2.2) Computing rho_i = |psi_i|^2 using fourwf (if fp is not 0).
                        if (this%use_precomputed_rhoi) then
                            rhoi_r = this%precomputed_rhoi(:, :, this%precomputed_rhoi_indices(iband, ikpt, isppol))    ! TODO : useless copy here, use a pointer ?
                        else
                            call compute_rhoi_coll(this, dtset, mpi_enreg, iband, ikpt, isppol, rhoi_r)
                        end if
                        
                        !2.3) delta_occ(i) = fp_i * dot(rho_i, delta_V)
                        ! dot-product in the up/down basis (equal to the dot product in the Pauli basis) :
                        ! rhoi_r has only one spin-component corresponding to isppol (up=down if nsppol=1, up or down if nsppol=2).
                        if (dtset%nspden == 1) then
                            delta_occ(i_eigen) = fp * dot_product(rhoi_r(:, 1), delta_V(:, isppol)) * this%dvol
                        elseif (dtset%nspden == 2) then
                            delta_occ(i_eigen) = fp * dot_product(rhoi_r(:, 1), delta_V(:, 1) + (1-2*(isppol-1)) * delta_V(:, 2)) * this%dvol
                            ! delta_V(:, 1) + (1-2*(isppol-1)) * delta_V(:, 2) is delta_V in the (up/down) coordinate 'isppol'.
                        end if

                    ! Non collinear spins - Wavefunctions have two spins components.
                    elseif (dtset%nspinor == 2 .and. dtset%nspden==4) then

                        !2.2) Computing rho_i (4-dim, in pauli basis).
                        if (this%use_precomputed_rhoi) then
                            rhoi_r = this%precomputed_rhoi(:, :, this%precomputed_rhoi_indices(iband, ikpt, isppol))    ! TODO : useless copy here, use a pointer ?
                        else
                            call compute_rhoi_noncoll(this, dtset, mpi_enreg, iband, ikpt, isppol, rhoi_r)
                        end if
                        
                        ! dot product in Pauli basis :
                        delta_occ(i_eigen) = zero
                        do ispden = 1, 4
                            ! rhoi_r has 4 spin-components in the pauli basis that all needs to be multiplied to the corresponding component in delta_V
                            delta_occ(i_eigen) = delta_occ(i_eigen) + fp * dot_product(rhoi_r(:, ispden), delta_V(:, ispden)) * this%dvol
                        end do

                    else
                        ABI_BUG("TODO non-collinear magnetism with nspden=/4")
                        ! The calculation above is probably true even when nspden=1, TODO: check
                    end if

                    delta_occ_tot = delta_occ_tot + delta_occ(i_eigen) * dtset%wtk(ikpt)

                end do

            end do  !ikpt
        end do  !isppol
        
        ABI_FREE(rhoi_r)

        !MPI parallelization over kpoints : sum delta_occ on all processors.
        ier = 0
        call xmpi_sum(delta_occ, mpi_enreg%comm_kptband, ier)
        call xmpi_sum(delta_occ_tot, mpi_enreg%comm_kptband, ier)
        call xmpi_sum(dos_fermie, mpi_enreg%comm_kptband, ier)
        
        ! 2) Fermi-level variation

        delta_fermie = delta_occ_tot / dos_fermie

        !Loop over spins and kpoints to do delta_occ(eigenvalue) -= f'(eigenvalue) * delta_fermie
        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt
                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt
                nband_k = dtset%nband(i_kpt_sppol)

                iband1 = needed_bands_bounds(1, i_kpt_sppol)
                iband2 = needed_bands_bounds(2, i_kpt_sppol)
                do iband = iband1, iband2   ! Loop over needed bands

                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)
                    eigenval = this%eigen(i_eigen)
                    fp = derivative_occ(dtset%occopt, eigenval, this%fermie, this%precon_tsmear) * maxocc
                    !fp = doccde(i_eigen)
                    delta_occ(i_eigen) = delta_occ(i_eigen) - fp * delta_fermie

                end do
            end do  !ikpt
        end do  !isppol
        !ABI_FREE(doccde)

        ABI_FREE(needed_bands_bounds)
        ABI_FREE(needed_bands_number)

    end subroutine compute_delta_occ

    subroutine compute_delta_wf(this, dtset, mpi_enreg, delta_V, delta_wf)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type), intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        real(dp), intent(inout) :: delta_V(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: delta_wf(:, :)
    
        !Local variables-------------------------------
        !scalars
        integer :: n1, n2, n3, n4, n5, n6
        integer :: isppol, ispden, ikpt, nband_k, npw_k, istwf_k, iband, jband
        integer :: i_eigen, j_eigen, ier, ndat, option
        integer :: i_kg(2), j_cg(2), i_cg(2)
        real(dp) :: fi, fj, ddiff, coeff
        integer :: tim_fourwf
        integer :: gbound(2*dtset%mgfft+8,2)
        real(dp) :: dotr, doti
        !arrays
        real(dp), allocatable :: delta_V_wf_i_r(:, :,:,:)
        real(dp), allocatable :: delta_V_wf_i(:,:)
        real(dp), allocatable :: delta_V_aug(:,:,:)
        !dummy
        integer :: dummy_int
        real(dp) :: dummy_real
        real(dp) :: dummy_denpot(0, dtset%ngfft(5), dtset%ngfft(6)), dummy_fofg(2, 0)
    
        ! ***************************************************************************
        write(6,*)'chi0diel compute_delta_wf'; flush(6) !DEBUG
        ! NOT WORKING
    
        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)
    
        delta_wf = zero
    
        ABI_MALLOC(delta_V_wf_i_r, (2, n4, n5, n6))
        ABI_MALLOC(delta_V_wf_i, (2, dtset%mpw))
        ABI_MALLOC(delta_V_aug, (n4, n5, n6))
    
        !Loop over spins and kpoints
        do isppol = 1, dtset%nsppol
            
            ispden = isppol     !TODO noncoll
            ! Transfer delta_V to the augmented (wavefunction) fft-grid
            call fftpac(ispden, mpi_enreg, dtset%nspden, n1, n2, n3, n4, n5, n6, dtset%ngfft, delta_V, delta_V_aug, 2)
    
            do ikpt = 1, dtset%nkpt
                nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
                !MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if
                npw_k = this%npwarr(ikpt)       ! Number of plane-wave at this kpt.
                istwf_k = dtset%istwfk(ikpt)    ! Option parameter that describes the storage of wfs at this kpt.
                i_kg = this%kg_indices(:, ikpt)
                call sphereboundary(gbound, istwf_k, this%kg(:, i_kg(1):i_kg(2)), dtset%mgfft, npw_k)    ! Computes gbound.
    
                do iband = 1, nband_k

                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)  ! Index of (iband, ikpt, isppol) in eigen array.
                    if (abs(this%eigen(i_eigen) - this%fermie) > tol1) then ! TODO : condition as input
                        cycle
                    end if
                    fi = this%occ(i_eigen)

                    i_cg = this%cg_indices(:, iband, ikpt, isppol)

                    !Input parameters for fourwf :
                    ndat = 1            ! Only one FFT.
                    tim_fourwf = 0
                    ! Multiply the wave function by the potential variation delta_V (ifft -> real space multiplication -> fft)
                    option = 0      ! ifft
                    call fourwf(dummy_int, dummy_denpot, this%cg(:, i_cg(1):i_cg(2)), dummy_fofg, delta_V_wf_i_r,  &
                    &           gbound, gbound, istwf_k, this%kg(:, i_kg(1):i_kg(2)), this%kg(:, i_kg(1):i_kg(2)), &
                    &           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                    &           dummy_int, n4, n5, n6, option, tim_fourwf, dummy_real, dummy_real)
                    call cg_vlocpsi(n4, n5, n6, n4, n5, n6, 1, 1, delta_V_aug, delta_V_wf_i_r)
                    option = 3      ! fft
                    call fourwf(dummy_int, dummy_denpot, dummy_fofg, delta_V_wf_i(:, 1:npw_k), delta_V_wf_i_r,  &
                    &           gbound, gbound, istwf_k, this%kg(:, i_kg(1):i_kg(2)), this%kg(:, i_kg(1):i_kg(2)), &
                    &           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                    &           dummy_int, n4, n5, n6, option, tim_fourwf, dummy_real, dummy_real)
                    write(6,*)'chi0diel compute_delta_wf - ok2 '; flush(6) !DEBUG
                    ! result stored in delta_V_wf_i(:, i_cg(1):i_cg(2))
    
                    do jband = 1, nband_k
                        write(6,*)'chi0diel compute_delta_wf - if yes '; flush(6) !DEBUG

                        if (abs(this%eigen(j_eigen) - this%fermie) > tol1) then ! TODO : condition as input
                            cycle
                        end if

                        if (iband == jband) then
                            cycle
                        end if  ! i=j contribution computed in compute_delta_occ

                        j_eigen = get_eigen_index(dtset, jband, ikpt, isppol)  ! Index of (jband, ikpt, isppol) in eigen array.
                        j_cg = this%cg_indices(:, jband, ikpt, isppol)
                        fj = this%occ(j_eigen)

                        ! Compute an equivalent of coeff=1/(eigen_i - eigen_j) that ensure correct compensation of the terms in compute_delta_occ
                        if (abs(this%eigen(i_eigen) - this%eigen(j_eigen)) < dtset%tsmear * tol10) then
                            ddiff = derivative_occ(dtset%occopt, (this%eigen(i_eigen)+this%eigen(j_eigen))/2, this%fermie, dtset%tsmear)
                        else
                            ddiff = (fi - fj)/(this%eigen(i_eigen) - this%eigen(j_eigen))
                        end if
                        coeff = ddiff * fi/(fi**2+fj**2)  ! From DFTK
                        write(6,*)'chi0diel compute_delta_wf - ok3 '; flush(6) !DEBUG

                        ! Compute dot product between wavefunction (j) and delta_V
                        call dotprod_g(dotr, doti, istwf_k, npw_k, 2, this%cg(:, j_cg(1):j_cg(2)), delta_V_wf_i(:, 1:npw_k), 0, mpi_enreg%comm_spinorfft) ! TODO : check dotprof(psi_i, psi_i) = 1
                        write(6,*)'chi0diel compute_delta_wf - ok4 '; flush(6) !DEBUG

                        delta_wf(1, i_cg(1):i_cg(2)) = delta_wf(1, i_cg(1):i_cg(2)) + &
                        &                              coeff * ( dotr * this%cg(1, j_cg(1):j_cg(2)) - doti * this%cg(2, j_cg(1):j_cg(2)) )
                        delta_wf(2, i_cg(1):i_cg(2)) = delta_wf(2, i_cg(1):i_cg(2)) + &
                        &                              coeff * ( dotr * this%cg(2, j_cg(1):j_cg(2)) + doti * this%cg(1, j_cg(1):j_cg(2)) )
                        write(6,*)'chi0diel compute_delta_wf - ok5 '; flush(6) !DEBUG
                    
                    end do
                end do
            end do
        end do
    
        ABI_FREE(delta_V_wf_i_r)
        ABI_FREE(delta_V_wf_i)
        ABI_FREE(delta_V_aug)
    
        !MPI parallelization over kpoints : sum delta_occ on all processors.
        ier = 0
        call xmpi_sum(delta_wf, mpi_enreg%comm_kpt, ier)

    end subroutine compute_delta_wf

    subroutine compute_delta_rho_from_delta_occ_only(this, dtset, mpi_enreg, delta_occ, delta_rho)
        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: delta_occ(size(this%eigen))
        real(dp), intent(inout) :: delta_rho(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        integer :: iband, isppol, ispden, ikpt, i_eigen, nband_k
        integer :: ier, rank
        integer :: maxocc
        integer :: iband1, iband2, i_kpt_sppol
        real(dp) :: fp
        !arrays
        integer, allocatable :: needed_bands_bounds(:, :)
        integer, allocatable :: needed_bands_number(:)
        real(dp), allocatable :: rhoi_r(:, :), delta_rho_g(:, :)
        
        ! *************************************************************************

        ABI_MALLOC(needed_bands_bounds, (2, dtset%nkpt*dtset%nsppol))
        ABI_MALLOC(needed_bands_number, (dtset%nkpt*dtset%nsppol))
        call get_needed_bands_chi0diag(this, dtset, mpi_enreg, needed_bands_bounds, needed_bands_number)

        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        delta_rho = zero

        !Loop over spins and kpoints
        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt

                i_kpt_sppol = ikpt+(isppol-1)*dtset%nkpt
                nband_k = dtset%nband(i_kpt_sppol)

                ! MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if
                
                iband1 = needed_bands_bounds(1, i_kpt_sppol)
                iband2 = needed_bands_bounds(2, i_kpt_sppol)
                
                do iband = iband1, iband2   ! Loop over needed bands

                    rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)
                    if (.not.(1 + mpi_enreg%bandpp*rank <= iband .and. iband <=mpi_enreg%bandpp*(rank+1))) then
                        cycle
                    end if  ! TODO : create a cycle function, that is also correct in LOBPCG

                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)  ! Index of (iband, ikpt, isppol) in eigen array.
                    fp = derivative_occ(dtset%occopt, this%eigen(i_eigen), this%fermie, this%precon_tsmear) * maxocc

                    ! Compute and add the contribution to delta_rhol.
                    if (dtset%nspinor == 1) then
                        ! Collinear magnetism or no magnetism
                    
                        ispden = isppol
                        if (this%use_precomputed_rhoi) then
                        ! We use the precomputed orbital density.
                            delta_rho(:, ispden) = delta_rho(:, ispden) + &
                            &                      dtset%wtk(ikpt) * delta_occ(i_eigen) * &
                            &                      this%precomputed_rhoi(:, 1, this%precomputed_rhoi_indices(iband, ikpt, isppol))
                        else
                        ! We need to recompute the orbital density.
                            ABI_MALLOC(rhoi_r, (this%nfftprc, 1))
                            call compute_rhoi_coll(this, dtset, mpi_enreg, iband, ikpt, isppol, rhoi_r)
                            delta_rho(:, ispden) = delta_rho(:, ispden) + &
                            &                      dtset%wtk(ikpt) * delta_occ(i_eigen) * rhoi_r(:, 1)
                            ABI_FREE(rhoi_r)
                        end if
                        ! delta_rho in up/down representation = expected representation by symrhg

                        ! TODO : here it is assumed that nsppol=2 => nspden=2 ... 

                    elseif (dtset%nspden == 4) then
                        ! Non collinear magnetism

                        if (this%use_precomputed_rhoi) then
                        ! We use the precomputed orbital density.
                            do ispden = 1, 4
                                delta_rho(:, ispden) = delta_rho(:, ispden) + &
                                &                      dtset%wtk(ikpt) * delta_occ(i_eigen) * &
                                &                      this%precomputed_rhoi(:, ispden, this%precomputed_rhoi_indices(iband, ikpt, isppol))
                            end do
                        else
                        ! We need to recompute the orbital density.
                            ABI_MALLOC(rhoi_r, (this%nfftprc, 4))
                            call compute_rhoi_noncoll(this, dtset, mpi_enreg, iband, ikpt, isppol, rhoi_r)
                            do ispden = 1, 4
                                delta_rho(:, ispden) = delta_rho(:, ispden) + &
                                &                      dtset%wtk(ikpt) * delta_occ(i_eigen) * rhoi_r(:, ispden)
                            end do
                            ABI_FREE(rhoi_r)
                        end if
                        ! delta_rho in Pauli representation = NOT the expected representation by symrhg

                    else
                        ABI_BUG("TODO non-collinear magnetism with nspden=/4")
                    end if
                    
                end do  ! iband

            end do  ! ikpt
        end do  ! isppol

        ABI_FREE(needed_bands_bounds)
        ABI_FREE(needed_bands_number)

        ! MPI parallelization over kpoints and bands : sum delta_rho on all processors.
        ier = 0
        call xmpi_sum(delta_rho, mpi_enreg%comm_kptband, ier)

        ! Symmetrization
        ABI_MALLOC(delta_rho_g, (2, this%nfftprc))
        if (dtset%nspinor == 1) then
        ! In collinear magnetism, we can readily apply symrhg to delta_rho (correct spin representation).
            call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, dtset%nfft, dtset%nfft, dtset%ngfft, dtset%nspden, dtset%nsppol, &
            &   dtset%nsym, this%phnons, delta_rho_g, delta_rho, this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)
        else
        ! In non-collinear magnetism, we apply symrhg independantly to each spin component, 
        ! to avoid having to change the spin representation of delta_rho.
            do ispden = 1, dtset%nspden
                call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, dtset%nfft, dtset%nfft, dtset%ngfft, 1, 1, &
                &   dtset%nsym, this%phnons, delta_rho_g, delta_rho(:, ispden:ispden), this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)
            end do
        end if
        ABI_FREE(delta_rho_g)

        call to_pauli(this, 1, delta_rho)

    end subroutine compute_delta_rho_from_delta_occ_only

    subroutine apply_chi0_diag(this, dtset, mpi_enreg, vec_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        !arrays
        real(dp), allocatable :: delta_occ(:)
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_chi0_diag')
       
        ABI_MALLOC(delta_occ, (size(this%eigen)))
        call compute_delta_occ(this, dtset, mpi_enreg, vec_r, delta_occ)
        call compute_delta_rho_from_delta_occ_only(this, dtset, mpi_enreg, delta_occ, vec_r)

        ABI_FREE(delta_occ)

    end subroutine apply_chi0_diag

    subroutine compute_delta_rho(this, dtset, mpi_enreg, delta_occ, delta_wf, delta_rho)
        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: delta_occ(size(this%eigen))
        real(dp), intent(inout) :: delta_wf(2, size(this%cg, 2))
        real(dp), intent(inout) :: delta_rho(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        integer :: iband, isppol, ispden, ikpt, istwf_k, i_eigen, nband_k, npw_k
        integer :: ier, rank
        integer :: n1, n2, n3, n4, n5, n6
        integer :: tim_fourwf, ndat, option
        integer :: maxocc
        real(dp) :: fp
        !arrays
        integer :: i_cg(2), i_kg(2)
        integer :: gbound(2*dtset%mgfft+8,2)
        integer, allocatable :: kg_k(:, :)
        real(dp), allocatable :: rho_aug_r_i(:, :, :), wf_aug_r_i(:, :, :, :), delta_wf_aug_r_i(:, :, :, :), delta_rho_aug_r(:, :, :, :)
        real(dp), allocatable :: delta_rho_coarse_r(:, :)
        real(dp), allocatable :: delta_rho_g(:, :)
        !dummy arguments
        integer :: dummy_int
        real(dp) :: dummy_real
        real(dp) ::  dummy_denpot(0, dtset%ngfft(5), dtset%ngfft(6)), dummy_fofgout(2, 0), dummy_fofrout(2, dtset%ngfft(4), dtset%ngfft(5), dtset%ngfft(6))
        
        ! *************************************************************************
        !  TODO : NOT WORKING
        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)

        ABI_MALLOC(delta_rho_aug_r, (n4, n5, n6, dtset%nspden))
        delta_rho_aug_r = zero
        ABI_MALLOC(rho_aug_r_i, (n4, n5, n6))
        ABI_MALLOC(wf_aug_r_i, (2, n4, n5, n6))
        ABI_MALLOC(delta_wf_aug_r_i, (2, n4, n5, n6))

        !Loop over spins and kpoints
        do isppol =1, dtset%nsppol
            do ikpt = 1, dtset%nkpt

                nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)

                !MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if

                istwf_k = dtset%istwfk(ikpt)        ! Option parameter that describes the storage of wfs at this kpt.
                if (.not. (dtset%paral_kgb == 1 .and. dtset%npband > 1)) then
                    npw_k = this%npwarr(ikpt)           ! Number of plane-wave at this kpt.
                    ABI_MALLOC(kg_k, (3, npw_k))
                    i_kg = this%kg_indices(:, ikpt)
                    kg_k = this%kg(:, i_kg(1):i_kg(2))
                    call sphereboundary(gbound, istwf_k, kg_k, dtset%mgfft, npw_k)    ! Computes gbound.
                else
                    npw_k = bandfft_kpt(ikpt)%npw_tot           ! Number of plane-wave (after transpose) at this kpt.
                    ABI_MALLOC(kg_k, (3, npw_k))
                    kg_k = bandfft_kpt(ikpt)%kg_k_gather        ! Reduced plane-wave coordinate (k+G) of this kpt (after transpose).
                    call sphereboundary(gbound, istwf_k, kg_k, dtset%mgfft, npw_k)
                end if      ! TODO : Not very nice to have a if here ....

                rank = xmpi_comm_rank(mpi_enreg%comm_bandfft)
                do iband = 1 + mpi_enreg%bandpp*rank, mpi_enreg%bandpp*(rank+1)

                    !Indices
                    i_eigen = get_eigen_index(dtset, iband, ikpt, isppol)       ! Index of (iband, ikpt, isppol) in eigen array.
                    fp = derivative_occ(dtset%occopt, this%eigen(i_eigen), this%fermie, this%precon_tsmear) * maxocc

                    if (abs(fp) > this%deigvals_tol_fp) then
                    
                        !Input parameters for fourwf :
                        ndat = 1            ! Only one FFT.
                        tim_fourwf = 0
                        rho_aug_r_i = zero
                        option = 1

                        !1) Contribution of the occupation variation
                        ! IFFT for wavefunction
                        if (.not. (dtset%paral_kgb == 1 .and. dtset%npband > 1)) then
                            i_cg = this%cg_indices(:, iband, ikpt, isppol)  ! Indices range of (iband, ikpt, isppol) in cg array
                            call fourwf(1, rho_aug_r_i, this%cg(:, i_cg(1):i_cg(2)), dummy_fofgout, wf_aug_r_i,  &
                            &           gbound, gbound, istwf_k, kg_k, kg_k, &
                            &           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            &           dummy_int, n4, n5, n6, option, tim_fourwf, one, one)
                        else
                            !i_cg = this%cg_fft_indices(:, iband, ikpt, isppol)  ! Indices range of (iband, ikpt, isppol) in cg array
                            !call fourwf(1, rho_aug_r_i, this%cg_fft(:, i_cg(1):i_cg(2)), dummy_fofgout, wf_aug_r_i,  &
                            !&           gbound, gbound, istwf_k, kg_k, kg_k, &
                            !&           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            !&           dummy_int, n4, n5, n6, option, tim_fourwf, one, one)
                        end if  ! TODO : change this to remove if ! (use pointers ?)

                        ! Compute and add the contribution to delta_rho
                        if (dtset%nspinor == 1) then
                            ispden = isppol
                            delta_rho_aug_r(:, :, :, ispden) = delta_rho_aug_r(:, :, :, ispden) + &
                            &                                  dtset%wtk(ikpt)/this%ucvol * delta_occ(i_eigen) * rho_aug_r_i
                        else
                            ABI_BUG("Non collinear not implemented")
                        end if
                        
                        !2) Contribution of wavefunction variation (if applicable)
                        if (.false.) then 
                        !if (norm2(delta_wf(:, i_cg(1):i_cg(2))) > tol14) then ! TODO : what tol ?
                            ! TODO : his part is probably incorrect ! 
                            ! IFFT for delta_wf (wavefunction variation)
                            call fourwf(1, dummy_denpot, delta_wf(:, i_cg(1):i_cg(2)), dummy_fofgout, delta_wf_aug_r_i,  &
                            &           gbound, gbound, istwf_k, this%kg(:, i_kg(1):i_kg(2)), this%kg(:, i_kg(1):i_kg(2)), &
                            &           dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            &           dummy_int, n4, n5, n6, 0, tim_fourwf, dummy_real, dummy_real)
                            
                            ! Sum over MPI processes (for parallelization over bands).
                            ier = 0
                            call xmpi_sum(delta_wf_aug_r_i, mpi_enreg%comm_bandfft, ier)
                            
                            ! Compute and add the contribution to delta_rho
                            if (dtset%nspinor == 1) then
                                ispden = isppol
                                delta_rho_aug_r(:, :, :, ispden) = delta_rho_aug_r(:, :, :, ispden) + &
                                &                                  dtset%wtk(ikpt) * this%occ(i_eigen) * &
                                &                                  ((delta_wf_aug_r_i(1, :, : , :)*wf_aug_r_i(1, :, :, :)) + &
                                &                                   (delta_wf_aug_r_i(2, :, : , :)*wf_aug_r_i(2, :, :, :))) ! TODO  : false ?
                            else
                                ABI_BUG("Non collinear not implemented")
                            end if
                        end if
                    end if

                end do  ! iband

                ABI_FREE(kg_k)
            end do  ! ikpt
        end do  ! isppol

        ABI_FREE(rho_aug_r_i)
        ABI_FREE(wf_aug_r_i)
        ABI_FREE(delta_wf_aug_r_i)
                  
        do ispden=1, dtset%nspden
            call transfer_grid(this, dtset, mpi_enreg, ispden, delta_rho_aug_r(:, :, :, ispden), delta_rho)
        end do
        
        ABI_FREE(delta_rho_aug_r)

        ! MPI parallelization over kpoints and bands : sum delta_rho on all processors.
        ier = 0
        call xmpi_sum(delta_rho, mpi_enreg%comm_kptband, ier)

        ABI_MALLOC(delta_rho_g, (2, this%nfftprc))
        call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, dtset%nfft, dtset%nfft, dtset%ngfft, dtset%nspden, dtset%nsppol, &
        &   dtset%nsym, this%phnons, delta_rho_g, delta_rho, this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)  ! TODO : carefull with tnons and symrel (do they get updated during structure optimization ?)
        ABI_FREE(delta_rho_g)
        ! TODO : deal with symmetries when spin (nsppol in non coll)

        call to_pauli(this, 1, delta_rho)

    end subroutine compute_delta_rho

    subroutine apply_chi0_quasidiag(this, dtset, mpi_enreg, vec_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        !arrays
        real(dp), allocatable :: delta_occ(:), delta_wf(:, :)
        
        ! *************************************************************************
       
        ABI_MALLOC(delta_occ, (size(this%eigen)))
        ABI_MALLOC(delta_wf, (2, size(this%cg, 2)))
        call compute_delta_occ(this, dtset, mpi_enreg, vec_r, delta_occ)
        delta_wf = zero
        !call compute_delta_wf(this, dtset, mpi_enreg, vec_r, delta_wf)
        vec_r = zero
        call compute_delta_rho(this, dtset, mpi_enreg, delta_occ, delta_wf, vec_r)

        ABI_FREE(delta_occ)
        ABI_FREE(delta_wf)

    end subroutine apply_chi0_quasidiag

    !****f* m_precon/apply_chi0
    !! NAME
    !!  apply_chi0
    !!
    !! FUNCTION
    !!  Apply the model (determined by iprcel) chi0 operator to the vector vec_r (in place).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  ispden      = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_r (nfftprc, nspden) = Vector (in direct space) to which the model chi0 operator is applied (in place).
    !!                            When nspden > 1 vec_r is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0(this, dtset, mpi_enreg, vec_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_r(this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: ispden
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_chi0')
       
        !Kerker with user_defined parameter dielng
        if (this%iprcel == 210) then
            ispden = 1
            vec_r(:, ispden) = (-1/(4*pi*(this%dielng)**2)) * vec_r(:, ispden)
            do ispden = 2, dtset%nspden
                vec_r(:, ispden) = 0
            end do
        end if

        !Kerker based on the DOS.
        if (this%iprcel == 211) then
            ispden = 1
            vec_r(:, ispden) = -this%dos(1) * vec_r(:, ispden)
            do ispden = 2, dtset%nspden
                vec_r(:, ispden) = 0
            end do
        end if

        !Kerker based on the DOS, that can be different between the up and down chanels in spin-polarized cases.
        if (this%iprcel == 212) then
            vec_r(:, 1) = 0
            do ispden = 1, dtset%nspden
                vec_r(:, 1) = vec_r(:, 1) - this%dos(ispden) * vec_r(:, ispden)
            end do
            do ispden = 2, dtset%nspden
                vec_r(:, ispden) = - this%dos(ispden) * vec_r(:, 1)
            end do
        end if
        
        !LDOS model
        if (this%iprcel == 200 .or. this%iprcel == 201) then
            call apply_chi0_ldos(this, dtset, mpi_enreg, vec_r)
        end if

    end subroutine apply_chi0

    !!****f* ABINIT/apply_adjdielmat
    !! NAME
    !!  apply_adjdielmat
    !!
    !! FUNCTION
    !!  Apply the adjoint dielectric matrix I-chi0*vc to the density rho_r (given in the Fourier space)
    !!  with a model chi0 operator.
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  rho_r       = Density vector (in direct space, in Pauli basis).
    !!
    !! OUTPUT
    !!  adjdielmat_rho_r = adjdielmat * rho_r
    !!
    !! NOTES
    !!
    !! SOURCE
    subroutine apply_adjdielmat(this, dtset, mpi_enreg, rho_r, adjdielmat_rho_r)

        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) :: rho_r(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: adjdielmat_rho_r(this%nfftprc, dtset%nspden)
        !Local variables-------------------------------
        real(dp), allocatable :: chi0_kxc_rho_r(:, :)
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_adjdielmat')
        if (this%use_precon) then

            if (this%iprcel == 299) then
            ! P=I : No preconditioning
                adjdielmat_rho_r = rho_r

            elseif (this%iprcel == 200) then
            ! More efficient implementation for the LDOS preconditioner.
                call apply_adjdielmat_ldos(this, dtset, mpi_enreg, rho_r, adjdielmat_rho_r)

            elseif (this%iprcel == 202) then
            ! When iprcel = 205, P = (I - chi0_ldos*vc - chi0_diag*Kxc)
                
                !1) Compute adjdielmat_rho_r = rho_r - chi0_ldos * vc *rho_r
                call apply_adjdielmat_ldos(this, dtset, mpi_enreg, rho_r, adjdielmat_rho_r)

                !2) Add -(chi0_diag * Kxc * rho_r) to adjdielmat_rho_r
                ABI_MALLOC(chi0_kxc_rho_r, (this%nfftprc, dtset%nspden))
                
                !2.2) Apply Kxc to vec_r
                call apply_kxc(this, dtset, mpi_enreg, rho_r, chi0_kxc_rho_r)
                !2.3) Apply chi0_diag to Kxc*vec_r (in place)
                call apply_chi0_diag(this, dtset, mpi_enreg, chi0_kxc_rho_r)
                !2.5) Add this contribution to adjdielmat_rho_r
                adjdielmat_rho_r = adjdielmat_rho_r - chi0_kxc_rho_r
                
                ABI_FREE(chi0_kxc_rho_r)

            elseif (this%iprcel == 203) then
            ! When iprcel = 203, P = (I - chi0_ldos*vc - chi0_quasidiag*Kxc)
                
                !1) Compute adjdielmat_rho_r = rho_r - chi0_ldos * vc *rho_r
                call apply_adjdielmat_ldos(this, dtset, mpi_enreg, rho_r, adjdielmat_rho_r)

                !2) Add -(chi0_diag * Kxc * rho_r) to adjdielmat_rho_r
                ABI_MALLOC(chi0_kxc_rho_r, (this%nfftprc, dtset%nspden))
                
                !2.2) Apply Kxc to vec_r
                call apply_kxc(this, dtset, mpi_enreg, rho_r, chi0_kxc_rho_r)
                !2.3) Apply chi0_quasidiag to Kxc*vec_r (in place)
                call apply_chi0_quasidiag(this, dtset, mpi_enreg, chi0_kxc_rho_r)
                !2.5) Add this contribution to adjdielmat_rho_r
                adjdielmat_rho_r = adjdielmat_rho_r - chi0_kxc_rho_r
                
                ABI_FREE(chi0_kxc_rho_r)

            else
            ! In the general case, P = (I - K*chi0_model) where K and chi0_model are defined 
            ! in the subroutine apply_kernel and apply_chi0 (depending on iprcel).
                
                adjdielmat_rho_r = rho_r
                !1) Apply the Kernel (vc or vc + Kxc depending on iprcel)
                call apply_kernel(this, dtset, mpi_enreg, adjdielmat_rho_r)
                !2) Applythe model chi0 operator
                call apply_chi0(this, dtset, mpi_enreg, adjdielmat_rho_r)
                !3) adjdielmat_rho_r = rho_r - K * chi0 * rho_r = adjdielmat * rho_r
                adjdielmat_rho_r = rho_r - adjdielmat_rho_r
                
            end if

        end if
    end subroutine apply_adjdielmat
    !!***
    
    !!****f* ABINIT/apply_dielmat
    !! NAME
    !!  apply_dielmat
    !!
    !! FUNCTION
    !!  Apply the model dielectric matrix (I - K*chi0) to the potential v_r (given in the direct/real space), 
    !!  with a model chi0 operator (contained in precon).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  v_r         = Potential vector (in real space, in Pauli basis)
    !!
    !! OUTPUT
    !!  dielmat_v_r = dielmat * v_r
    !!
    !! NOTES
    !!
    !! SOURCE
    subroutine apply_dielmat(this, dtset, mpi_enreg, v_r, dielmat_v_r)
    
        !Arguments ------------------------------------
        class(precon_object) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(in) ::  v_r(this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: dielmat_v_r(this%nfftprc, dtset%nspden)
        !Local variables-------------------------------
        real(dp), allocatable :: kxc_chi0_v_r(:, :), chi0_v_r(:, :)
        
        ! *************************************************************************
        if (this%precon_verbose>1) call wrtout(std_out, 'apply_dielmat')
        if (this%use_precon) then

            if (this%iprcel == 299) then
            ! P=I : No preconditioning
                dielmat_v_r = v_r
            
            elseif (this%iprcel == 200) then
            ! More efficient implementation for the LDOS preconditioner.
                call apply_dielmat_ldos(this, dtset, mpi_enreg, v_r, dielmat_v_r) 

            elseif (this%iprcel == 202) then
            ! When iprcel = 202, P = (I - vc*chi0_ldos - Kxc*chi0_diag)
                
                !1) Compute dielmat_v_r = v_r - vc * chi0_ldos *v_r
                call apply_dielmat_ldos(this, dtset, mpi_enreg, v_r, dielmat_v_r) 

                !2) Add -(Kxc * chi0_quasidiag * v_r) to dielmat_v_r
                
                !2.1) Apply chi0_quasidiag to v_r
                ABI_MALLOC(chi0_v_r, (this%nfftprc, dtset%nspden))
                chi0_v_r = v_r
                call apply_chi0_diag(this, dtset, mpi_enreg, chi0_v_r)
                !2.2) Apply Kxc to chi0*v_r
                ABI_MALLOC(kxc_chi0_v_r, (this%nfftprc, dtset%nspden))
                call apply_kxc(this, dtset, mpi_enreg, chi0_v_r, kxc_chi0_v_r)
                !2.3) Add this contribution to dielmat_v_r
                dielmat_v_r = dielmat_v_r - kxc_chi0_v_r
                ABI_FREE(kxc_chi0_v_r)
                ABI_FREE(chi0_v_r)

            elseif (this%iprcel == 203) then
            ! When iprcel = 203 , P = (I - vc*chi0_ldos - Kxc*chi0_quasidiag)
                
                !1) Compute dielmat_v_r = v_r - vc * chi0_ldos *v_r
                call apply_dielmat_ldos(this, dtset, mpi_enreg, v_r, dielmat_v_r) 

                !2) Add -(Kxc * chi0_quasidiag * v_r) to dielmat_v_r
                
                !2.1) Apply chi0_quasidiag to v_r
                ABI_MALLOC(chi0_v_r, (this%nfftprc, dtset%nspden))
                chi0_v_r = v_r
                call apply_chi0_quasidiag(this, dtset, mpi_enreg, chi0_v_r)
                !2.2) Apply Kxc to chi0*v_r
                ABI_MALLOC(kxc_chi0_v_r, (this%nfftprc, dtset%nspden))
                call apply_kxc(this, dtset, mpi_enreg, chi0_v_r, kxc_chi0_v_r)
                !2.3) Add this contribution to dielmat_v_r
                dielmat_v_r = dielmat_v_r - kxc_chi0_v_r
                ABI_FREE(kxc_chi0_v_r)
                ABI_FREE(chi0_v_r)

            else
            ! In the general case, P = (I - K*chi0_model) where K and chi0_model are defined 
            ! in the subroutine apply_kernel and apply_chi0 (depending on iprcel).
                
                dielmat_v_r = v_r
                !1) Apply the model chi0 operator
                call apply_chi0(this, dtset, mpi_enreg, dielmat_v_r)
                !2) Apply the Kernel (vc or vc + Kxc depending on iprcel)
                call apply_kernel(this, dtset, mpi_enreg, dielmat_v_r)
                !3) dielmat_v_r = v_r - K * chi0 * v_r = dielmat * v_r
                dielmat_v_r = v_r - dielmat_v_r
                
            end if

        end if
    end subroutine apply_dielmat
    !!***

    !****f* m_precon/apply_precon
    !! NAME
    !!  apply_precon
    !!
    !! FUNCTION
    !!  Apply a the preconditioner P^-1 to a given input vector 'vresid' where P is a model for the 
    !!  dielectric matrix or its adjoint based of a model of the non interacting susceptibility chi0.
    !!  More precisely, 
    !!      - if we are preconditioning potentials ('optres'=0), P is a model of the dielectric matrix (I-K*chi0)
    !!      - if we are preconditioning densities ('optres'=1), P is a model of the adjoint dielectric matrix (I-chi0*K).
    !!  The approximation are defined in 'apply_dielmat' and 'apply_adjdielmat' by Abinit input 'iprcel'.
    !!  
    !!  The preconditioner is applied by solving the linear equation P * 'vrespc' = 'vresid' iteratively, 
    !!  either using the GMRES method if P is well conditioned (positive definite) or using Ridge/Tikhonov regularization 
    !!  with the conjugate gradient if P might be ill-conditioned.
    !!
    !! INPUTS
    !!  dtset      = All input variables for this dataset.
    !!  mpi_enreg  = Information about MPI parallelization.
    !!  optreal    = Integer flag indicating whether the input is in real space (1) or reciprocal space (2).
    !!  optres     = Integer flag indicating whether we are preconditioning densities (1) or potentials (0).
    !!  vresid     = Residual vector to which the preconditioner is applied.
    !!
    !! OUTPUTS
    !!  vrespc     = Preconditioned residual vector.
    !!
    !! NOTES
    !!  'vresid' and 'vrespc' have a different shape than the typical density/potential vectors 
    !!  in the rest of this file, to match the shape needed in 'm_prcref'.
    !!
    !! SOURCE
    subroutine apply_precon(this, dtset, mpi_enreg, optreal, optres, vresid, vrespc)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        integer :: optreal, optres
        !arrays
        real(dp), intent(in) :: vresid(optreal*this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: vrespc(optreal*this%nfftprc, dtset%nspden)
        !Local variables-------------------------------
        !scalars
        integer :: ispden, start_ispden, end_ispden, n
        !arrays
        real(dp), allocatable :: rhs(:), est(:), P_rhs(:)
        real(dp), allocatable :: work_g(:, :, :)

        ! *************************************************************************
        if (this%use_precon) then

            !0.0) Update precon object
            call precon_update(this, dtset, mpi_enreg)

            ! The preconditioned density/potential residual vrespc = P^-1 * vresid is computed 
            ! by sovling the linear equation P * vrespc = vresid approximately with GMRES.

            n = dtset%nspden*1*this%nfftprc
            
            !0.1) Convert the input to direct/real space if needed.
            if (optreal==0) then
                ! vresid is given in the Fourier space : We need to do a ifft.
                ABI_MALLOC(work_g, (2, this%nfftprc, dtset%nspden))
                work_g = reshape(vresid, (/2,this%nfftprc, dtset%nspden/))
                call fourdp(1, work_g, vrespc(1:this%nfftprc, :), 1, mpi_enreg, this%nfftprc, dtset%nspden, this%ngfftprc, 0)
            else
                vrespc = vresid
            end if

            !0.2) Convert the input to the Pauli basis
            call to_pauli(this, optres, vrespc)

            !1) Right-hand side : rhs is vresid (flattened) in the direct/real space.
            ABI_MALLOC(rhs, (n))
            do ispden = 1, dtset%nspden
                ! Indices of the ispden component in the flattened (this%nfftprc, dtset%nspden)-array 'rhs'.
                start_ispden = 1+(ispden-1)*this%nfftprc
                end_ispden = ispden*this%nfftprc
                rhs(start_ispden:end_ispden) = vrespc(1:this%nfftprc, ispden)
            end do

            !2) Initial guess :
            ABI_MALLOC(est, (n))
            est = 0
            ! Is est = rhs a better starting point ?
            !est = rhs

            !3) Resolution of the linear system :
            if (this%precon_verbose>0) call wrtout(std_out, '| chi0diel preconditioning - Linear Solver:')
            if (this%use_ridgereg) then
                ! P is ill-conditionned :
                ! Ridge/Tikhonov regularization and CG : 
                ! We solve (P^*P + ridge_param*I) * est = P * rhs
                ! (P^*P + ridge_param*I) is self-adjoint and can be solved with CG.
                ABI_MALLOC(P_rhs, (n))
                call matvec(n, rhs, P_rhs)  ! TODO apply P_adj !!!
                call cg_linear_solver(n, ridge_matvec, P_rhs, est, (this%linsolve_maxiter-1)/2+1, this%linsolve_rtol, this%precon_verbose>0)
                ABI_FREE(P_rhs)
            else
                ! P is well conditionned :
                ! GMRES (P is not self-adjoint)
                call gmres_linear_solver(n, matvec, rhs, est, this%linsolve_maxiter, this%linsolve_rtol, this%precon_verbose>0)

            end if

            !4.0) Free arrays that might have been allocated by precon_update
            call precon_free_update(this)
        
            !4) Reshaping the final result :
            do ispden = 1, dtset%nspden
                ! Indices of the ispden component in the flattened (this%nfftprc, dtset%nspden)-array 'est'.
                start_ispden = 1+(ispden-1)*this%nfftprc
                end_ispden = ispden*this%nfftprc
                if (optreal==1) then
                    ! vrespc must be returned in the direct/real space.
                    vrespc(:, ispden) = est(start_ispden:end_ispden)
                else
                    ! vrespc must be returned in the fourier space, we need to do a fft.
                    call fourdp(1, work_g, est(start_ispden:end_ispden), -1, mpi_enreg, this%nfftprc, 1, this%ngfftprc, 0)
                    vrespc = reshape(work_g, (/2*this%nfftprc, dtset%nspden/))
                end if
            end do
            ! vrespc must be returned in the default Abinit spin-basis.
            call from_pauli(this, optres, vrespc)
            ! TODO noncoll : Check this basis change for potentials in Fourier space with non-collinear magnetism.
        
            ABI_FREE(rhs)
            ABI_FREE(est)
        
        end if
       
        contains

        ! Subroutine matvec that applies the model (adjoint-) dielectric matrix. -----------
        subroutine matvec(n_, x, y)
            integer, intent(in) :: n_
            real(dp), intent(inout), target :: x(n_), y(n_)
            type(c_ptr) :: x_c, y_c
            real(dp), pointer :: x_2d(:, :), y_2d(:, :)
        
        ! **********************************************************************************
        
            ! C-pointers to match the flattened arrays x and y to their 3D versions needed by
            ! 'apply_adjdielmat' and 'apply_dielmat'.
            x_c = c_loc(x)
            call c_f_pointer(x_c, x_2d, shape=[this%nfftprc, dtset%nspden])
            y_c = c_loc(y)
            call c_f_pointer(y_c, y_2d, shape=[this%nfftprc, dtset%nspden])
        
            if (optres==1) then
                ! We are preconditioning density residual so P models the adjoint dielectric matrix.
                call this%apply_adjdielmat(dtset, mpi_enreg, x_2d, y_2d)
            else if (optres==0) then
                ! We are preconditioning potential residual so P models the dielectric matrix.
                call this%apply_dielmat(dtset, mpi_enreg, x_2d, y_2d)
            end if

        end subroutine matvec ! ------------------------------------------------------------

        ! Subroutine ridge_matvec that applies the operator (P^*P + ridge_param*I) needed for ridge regularization.
        subroutine ridge_matvec(n_, x, y)
            integer, intent(in) :: n_
            real(dp), intent(inout), target :: x(n_), y(n_)
            type(c_ptr) :: x_c, y_c
            real(dp), pointer :: x_2d(:, :), y_2d(:, :), temp_2d(:, :)
        
        ! **********************************************************************************
        
            ! C pointers to match the flattened arrays x and y to their 3D versions needed by
            ! 'apply_adjdielmat' and 'apply_dielmat'.
            x_c = c_loc(x)
            call c_f_pointer(x_c, x_2d, shape=[this%nfftprc, dtset%nspden])
            y_c = c_loc(y)
            call c_f_pointer(y_c, y_2d, shape=[this%nfftprc, dtset%nspden])
            ABI_MALLOC(temp_2d, (this%nfftprc, dtset%nspden))  ! Temporary array for intermediate result
        
            if (optres==1) then
                ! We are preconditioning density residual so P models the adjoint dielectric matrix.
                call this%apply_adjdielmat(dtset, mpi_enreg, x_2d, temp_2d)
                call this%apply_dielmat(dtset, mpi_enreg, temp_2d, y_2d)
            else if (optres==0) then
                ! We are preconditioning potential residual so P models the dielectric matrix.
                call this%apply_dielmat(dtset, mpi_enreg, x_2d, temp_2d)
                call this%apply_adjdielmat(dtset, mpi_enreg, temp_2d, y_2d)
            end if
            y_2d = y_2d + this%ridge_param*x_2d
            ABI_FREE(temp_2d)

        end subroutine ridge_matvec ! ---------------------------------------------------------
        
    end subroutine apply_precon

end module m_precon
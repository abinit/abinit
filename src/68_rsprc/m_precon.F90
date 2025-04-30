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

    use defs_abitypes,  only : MPI_type
    use defs_basis
    use m_dtset
    use m_xmpi

    use defs_wvltypes
    use m_atomdata,     only : atom_length
    use m_dfpt_mkvxc,   only : dfpt_mkvxc
    use m_fft,          only : fourdp, fourwf, fftpac
    use m_fftcore,      only : sphereboundary
    use m_mkrho
    use m_mpinfo,       only : proc_distrb_cycle
    use m_paw_dmft
    use m_spacepar,     only : symrhg

#if defined HAVE_LINALG_MKL_OMATCOPY
    use mkl_rci, only : dfgmres, dfgmres_check, dfgmres_get, dfgmres_init
#endif
    
    implicit none
    private
    public :: precon_object
    public :: linsolve

    type, public :: precon_object
        integer  :: iprcel, nspden
        real(dp) :: dielng, diemix
        !Geometry :
        real(dp) :: gprimd(3, 3), rprimd(3, 3)
        real(dp) :: ucvol, dvol
        !For LDOS preconditioner :
        real(dp), pointer :: fermie
        real(dp), pointer :: cg(:, :), eigen(:), phnons(:, :, :)
        integer, pointer  :: kg(:, :), npwarr(:), irrzon(:, :, :)
        real(dp) :: tdos
        real(dp), allocatable :: ldos(:, :)
        !For local polarizability preconditioner :
        real(dp) :: gc
        real(dp), allocatable :: loc_pola(:, :)
            !To compute loc_pola :
        integer, pointer :: atindx1(:), nattyp(:)
        real(dp), pointer :: xred(:, :)
        real(dp), pointer :: rhor(:, :)
        !For ffts :
        integer  :: nfft
        !For Kxc :
        integer :: nkxc
        logical :: need_kxc
        real(dp), pointer :: kxc(:, :)

    contains
        procedure :: init => precon_init            ! Initializes the precon_object.
        procedure :: init_kxc => precon_init_kxc    ! Initializes kxc in the precon_object.
        procedure :: update => precon_update        ! Updates the precon_object according to iprcel.
        procedure :: free => precon_free            ! Dealocate arrays that are allocated in precon_init.
        procedure :: save => precon_save            ! Saves the LDOS or local polarizability contained in the precon_object in a file.
        procedure :: to_pauli => to_pauli           ! Basis change from the default abinit spin basis to the Pauli basis 
                                                    ! for density and potentials.
        procedure :: from_pauli => from_pauli       ! Basis change from the Pauli basis to the default abinit spin basis 
                                                    ! for density and potentials.
        procedure :: apply_kernel => apply_kernel   ! Applies the coulomb kernel vc to an input vector.
        procedure :: apply_chi0 => apply_chi0       ! Applies the model chi0 operator to an input vector.
        procedure :: save_applied_op_g => save_applied_op_g ! Saves the application of an operator in reciprocal space (for code validation).
        procedure :: save_applied_op_r => save_applied_op_r ! Saves the application of an operator in direct space (for code validation).

    end type precon_object

contains 

    !****f* m_precon/precon_init
    !! NAME
    !! precon_init
    !!
    !! FUNCTION
    !! Initializes the precon_object.
    !!
    !! INPUTS
    !!  dtset    = all input variables for this dataset
    !!  gprimd   = dimensional reciprocal space primitive translations
    !!  rprimd   = dimensional real space primitive translations
    !!  ucvol    = unit cell volume
    !!  cg       = wf in G space
    !!  eigen    = array of eigenvalues
    !!  fermie   = fermi energie
    !!  irrzon   = irreducible zone data
    !!  kg       = reduced planewave coordinates
    !!  npwarr   = number of planewaves and boundary planewaves at each k
    !!  phnons   = nonsymmorphic translation phases
    !!
    !! SOURCE
    subroutine precon_init(this, dtset, atindx1, cg, eigen, fermie, gprimd, &
        &   irrzon, kg, nattyp, npwarr, phnons, rhor, rprimd, ucvol, xred)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        real(dp), intent(in) :: ucvol
        real(dp), intent(in), target :: fermie

        !arrays
        real(dp), intent(in) :: gprimd(:, :), rprimd(:, :)
        integer, intent(in), target  :: irrzon(:, :, :), kg(:, :), npwarr(:)
        integer, intent(in), target :: atindx1(:), nattyp(:)
        real(dp), intent(in), target :: cg(:, :), eigen(:), phnons(:, :, :)
        real(dp), intent(in), target :: rhor(:, :)
        real(dp), intent(in), target :: xred(:, :)

        ! *************************************************************************
        !Constant data from dtset
        this%dielng = dtset%dielng
        this%diemix = dtset%diemix
        this%iprcel = dtset%iprcel
        this%nfft   = dtset%nfft
        this%nspden = dtset%nspden
        !Other constants
        this%dvol   = ucvol/this%nfft ! factor for integrals in real space: sum(f) * dvol ~ integral f
        this%gprimd = gprimd
        this%rprimd = rprimd
        this%ucvol  = ucvol
        this%need_kxc = .false.
        !Pointers
        this%atindx1 => atindx1 
        this%cg     => cg
        this%eigen  => eigen
        this%fermie => fermie
        this%irrzon => irrzon
        this%kg     => kg
        this%nattyp => nattyp
        this%npwarr => npwarr
        this%phnons => phnons
        this%rhor   => rhor
        this%xred   => xred
        !Initializing LDOS specific variables
        if (this%iprcel == 202) then
            !Allocating the array containing ldos
            ABI_MALLOC(this%ldos, (this%nfft, this%nspden))
        end if
        !Initializing ... specific variables
        if (this%iprcel == 203) then
            !Preparing the allocation of Kxc
            this%need_kxc = .true.
            if (dtset%xclevel==1) then  !LDA
                this%nkxc = 2*min(dtset%nspden,2)-1
            else if (dtset%xclevel==2)then  !GGA
                if (dtset%nspden==1) then
                    this%nkxc = 7
                else if (dtset%nspden==2) then
                    this%nkxc = 19
                end if
            end if
        end if
        !Initializing loc_pola specific variables
        if (this%iprcel == 204) then
            !this%gc = 1.0 ! TODO
            !Allocating the array containing the local polarizability
            !ABI_MALLOC(this%loc_pola, (this%nfft, this%nspden))
        end if

    end subroutine precon_init

    !****f* m_precon/precon_init_kxc
    !! NAME
    !! precon_init_kxc
    !!
    !! FUNCTION
    !! Initializes the exchange and correlation kernel (kxc) in the precon_object.
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

        if (this%need_kxc) then
            this%kxc => kxc
        end if

    end subroutine precon_init_kxc

    !****f* m_precon/precon_update
    !! NAME
    !! update
    !!
    !! FUNCTION
    !! Updates the precon_object :
    !!      For LDOS preconditioning (iprcel=202) : 
    !!          Computes the new ldos (local density of state) with current wavefunctions 
    !!          and the new tdos (total density of state = integral of ldos).
    !!      For Pola preconditioning :
    !!          Computes the local polarizability if istep = 1 (first SCF iteration).
    !!
    !! INPUTS
    !!  dtset     = All input variables for this dataset.
    !!  istep     = SCF step.
    !!  mpi_enreg = Information about MPI parallelization.
    !!
    !! SOURCE
    subroutine precon_update(this, dtset, istep, mpi_enreg)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type), intent(in) :: dtset
        integer, intent(in) :: istep
        type(MPI_type), intent(inout) :: mpi_enreg
        
        !Local variables-------------------------------
        integer :: ispden

        ! *************************************************************************

        !LDOS
        if (this%iprcel == 202) then 
            !update ldos
            call compute_ldos(dtset, this%cg, this%eigen, this%fermie, this%gprimd, this%irrzon, &
            &   this%kg, mpi_enreg, this%nfft, this%npwarr, this%phnons, this%rprimd, this%ucvol, &
            &   this%ldos)
            !update tdos
            this%tdos = sum(this%ldos(:, 1)) * this%dvol
            ! TODO : More options to control when the ldos is updated
        write(6,*)'chi0diel update : tdos, ldos(1:5, 1)', this%tdos, this%ldos(1:5, 1); flush(6) !DEBUG
            
        end if
        !write(6,*)'chi0diel update : kxc', this%kxc(1:5, :); flush(6) !DEBUG

        !Local polarizability
        if (this%iprcel == 204) then
            if (istep == 1) then    
                call compute_loc_pola(dtset, this%atindx1, this%gprimd, this%nattyp, this%nfft, this%nspden, &
                &   mpi_enreg, this%rhor, this%rprimd, this%xred, &
                &   this%loc_pola)
            end if
        end if
       
    end subroutine precon_update

    !****f* m_precon/precon_free
    !! NAME
    !! precon_free
    !!
    !! FUNCTION
    !! Dealocate arrays that are allocated in precon_init.
    !!
    !! SOURCE
    subroutine precon_free(this)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        
        ! *************************************************************************
       
        if (this%iprcel == 202) then
            !Deallocating the array containing ldos and tdos
            ABI_FREE(this%ldos)
        end if

        if (this%iprcel == 204) then
            ABI_FREE(this%loc_pola)
        end if

    end subroutine precon_free

    !****f* m_precon/compute_r
    !! NAME
    !! compute_r
    !!
    !! FUNCTION
    !! Computes the array of r-vectors (in REDUCED coordinates).
    !!
    !! INPUTS
    !!  ngfft   = All needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!
    !! OUTPUTS
    !!  r_vectors(3, :) = 3 coordinates of the r_vectors.
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
    !! get_r_vector
    !!
    !! FUNCTION
    !! Get the vector r (in REDUCED coordinates) of index ifft
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
    !! get_g_vector
    !!
    !! FUNCTION
    !! Get the vector g (in REDUCED coordinates) of index ifft
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
        i1 = modulo((ifft-1), n1) + 1
        i2 = modulo((ifft-1)/n1, n2) + 1
        i3 = ((ifft-1)/n1)/n2 + 1
        g(1) = modulo(i1-1 + n1/2, n1) - n1/2
        g(2) = modulo(i2-1 + n2/2, n2) - n2/2
        g(3) = modulo(i3-1 + n3/2, n3) - n3/2

    end function get_g_vector

    !****f* m_precon/precon_save
    !! NAME
    !! precon_save
    !!
    !! FUNCTION
    !! Saves the LDOS contained in the precon_object in a file named ldos.txt.
    !!     "     locale polarizability                "              loc_pola.txt. 
    !! (For code validation)
    !!
    !! INPUTS
    !!  ngfft   = All needed informations about the 3D FFT.
    !!  ispden  = Index of spin-density component.
    !!
    !! SOURCE
    subroutine precon_save(this, ngfft, ispden)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        integer :: ispden
        integer, intent(in) :: ngfft(:)

        !Local variables-------------------------------
        integer :: io, n, i

        ! *************************************************************************
        
        n = size(this%ldos, 1)

        if (this%iprcel==202) then
            n = size(this%ldos(:, ispden))
            ! Writing the file
            open(newunit=io, file="ldos.txt", status="replace", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%ldos(i, ispden)
                end do
            close(io)
        end if 

        if (this%iprcel==204) then
            n = size(this%loc_pola(:, ispden))
            ! Writing the file
            open(newunit=io, file="loc_pola.txt", status="replace", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%loc_pola(i, ispden)
                end do
            close(io)
        end if 

    end subroutine precon_save

    !!****f* ABINIT/to_pauli
    !! NAME
    !!  from_pauli
    !!
    !! FUNCTION
    !!  Basis change from the default spin-basis to the Pauli basis for potentials and densities
    !!  in the reciprocal (G-) space.
    !!
    !! INPUT/OUTPUT
    !!  opt                = 0 : v is a potential
    !!                       1 : v is a density
    !!  v(2, nfft, nspden) = On input : Potential/density in the default spin-basis.
    !!                       On output : Potential/density in the Pauli basis.
    !!
    !! SOURCE
    subroutine to_pauli(this, opt, v)
        class(precon_object), intent(in) :: this
        !Arguments ------------------------------------
        real(dp), intent(inout) ::  v(:, :, :)
        integer :: opt
        !Local variables-------------------------------
        integer :: nspden, nfft
        real(dp), allocatable :: temp(:)
        
        ! *************************************************************************
        nspden = size(v, 3)
        nfft = size(v, 2)

        !sigma_0, ... , sigma_3 are the Pauli matrices.
        if (opt == 0) then      !v is a potential
            if (nspden == 2) then
                !On input v(:, :, 1) is the spin-up potential and v(:, :, 2) is the spin-down potential.
                !On output the entire potential is v(:, :, 1)*sigma_0 + v(:, :, 2)*sigma_3.
                v(:, :, 1) = 0.5_dp*(v(:, :, 1) + v(:, :, 2))
                v(:, :, 2) = v(:, :, 1) - v(:, :, 2)
            else if (nspden == 4) then
                !                                   v(:, :, 1)  |   v(:, :, 2)
                !On input the entire potential is   ------------|---------------
                !                                   v(:, :, 3)  |   v(:, :, 4)
                !On output the entire potential is 
                !   v(:, :, 1)*sigma_0 + v(:, :, 2)*sigma_1 + v(:, :, 3)*sigma_2  + v(:, :, 4)*sigma_3
                v(:, :, 1) = 0.5_dp*(v(:, :, 1) + v(:, :, 4))
                v(:, :, 4) = v(:, :, 1) - v(:, :, 4)
                v(:, :, 2) = 0.5_dp*(v(:, :, 2) + v(:, :, 3))
                !v(:, :, 3) = i*(v(:, :, 2) - v(:, :, 3)) :
                ABI_MALLOC(temp, (nfft))
                temp = (v(1, :, 2) - v(1, :, 3))
                v(1, :, 3) = -1*(v(2, :, 2) - v(2, :, 3))  !  Re(i*) = -Im()
                v(2, :, 3) = temp                          !  Im(i*) = Re()
                ABI_FREE(temp)
            end if

        else if (opt == 1) then !v is a density
            if (nspden == 2) then
                !On input v(:, :, 1) is the total density and v(:, :, 2) is the spin-up density.
                !On output v(:, :, 1) is the total density and v(:, :, 2) is the spin density.
                v(:, :, 2) = 2*v(:, :, 2) - v(:, :, 1)
            end if
            !If nspden=4, the density is already given in the Pauli basis.
        end if
    end subroutine to_pauli

        !!****f* ABINIT/to_pauli
    !! NAME
    !!  from_pauli
    !!
    !! FUNCTION
    !!  Basis change from the default spin-basis to the Pauli basis for potentials and densities
    !!  in the direct (r-) space.
    !!
    !! INPUT/OUTPUT
    !!  opt                = 0 : v is a potential (TODO)
    !!                       1 : v is a density
    !!  v(nfft, nspden) = On input : Potential/density in the default spin-basis.
    !!                       On output : Potential/density in the Pauli basis.
    !!
    !! SOURCE
    subroutine to_pauli_r(this, opt, v)
        class(precon_object), intent(in) :: this
        !Arguments ------------------------------------
        real(dp), intent(inout) ::  v(:, :)
        integer :: opt
        !Local variables-------------------------------
        integer :: nspden, nfft
        real(dp), allocatable :: temp(:)
        
        ! *************************************************************************
        nspden = size(v, 2)
        nfft = size(v, 1)

        !sigma_0, ... , sigma_3 are the Pauli matrices.
        if (opt == 0) then      !v is a potential
            ABI_BUG("to_pauli in realspace TODO")

        else if (opt == 1) then !v is a density
            if (nspden == 2) then
                !On input v(:, 1) is the total density and v(:, 2) is the spin-up density.
                !On output v(:, 1) is the total density and v(:, 2) is the spin density.
                v(:, 2) = 2*v(:, 2) - v(:, 1)
            end if
            !If nspden=4, the density is already given in the Pauli basis.
        end if
    end subroutine to_pauli_r

    !!****f* ABINIT/from_pauli
    !! NAME
    !!  from_pauli
    !!
    !! FUNCTION
    !!  Basis change from the Pauli basis to the default spin-basis for potentials and densities
    !!  in the reciprocal (G-) space.
    !!
    !! INPUT/OUTPUT
    !!  opt                = 0 : v is a potential
    !!                       1 : v is a density
    !!  v(2, nfft, nspden) = On input : Potential/density in the Pauli basis
    !!                       On output : Potential/density in the default spin-basis.
    !!
    !! SOURCE
    subroutine from_pauli(this, opt, v)
        class(precon_object), intent(in) :: this
        !Arguments ------------------------------------
        real(dp), intent(inout) ::  v(:, :, :)
        integer :: opt
        !Local variables-------------------------------
        integer :: nspden, nfft
        real(dp), allocatable :: temp(:)
        
        ! *************************************************************************
        nspden = size(v, 3)
        nfft = size(v, 2)

        !sigma_0, ... , sigma_3 are the Pauli matrices.
        if (opt == 0) then      !v is a potential
            if (nspden == 2) then
                !On input the entire potential is v(:, :, 1)*sigma_0 + v(:, :, 2)*sigma_3.
                !On output v(:, :, 1) is the spin-up potential and v(:, :, 2) is the spin-down potential.
                v(:, :, 1) = v(:, :, 1) + v(:, :, 2)
                v(:, :, 2) = v(:, :, 1) - 2*v(:, :, 2)
            else if (nspden == 4) then
                !On input the entire potential is 
                !   v(:, :, 1)*sigma_0 + v(:, :, 2)*sigma_1 + v(:, :, 3)*sigma_2  + v(:, :, 4)*sigma_3
                !                                    v(:, :, 1)  |   v(:, :, 2)
                !On output the entire potential is   ------------|---------------
                !                                    v(:, :, 3)  |   v(:, :, 4)
                v(:, :, 1) = v(:, :, 1) + v(:, :, 4)
                v(:, :, 4) = v(:, :, 1) - 2*v(:, :, 4)
                !v(:, :, 2) = v(:, :, 2) - i*v(:, :, 3)
                !v(:, :, 3) = v(:, :, 2) + i* v(:, :, 3) :
                ABI_MALLOC(temp, (nfft))
                temp = v(3, 1, :)
                v(1, :, 3) = v(1, :, 2) - v(2, :, 3) !  Re(i*) = -Im()
                v(2, :, 2) = v(2, :, 2) + temp       !  Im(i*) = Re()
                ABI_FREE(temp)
                v(:, :, 2) = 2*v(:, :, 2) - v(:, :, 3)
            end if

        else if (opt == 1) then !v is a density
            if (nspden == 2) then
                !On input v(:, :, 1) is the total density and v(:, :, 2) is the spin density.
                !On output v(:, :, 1) is the total density and v(:, :, 2) is the spin-up density.
                v(:, :, 2) = 0.5_dp*(v(:, :, 1) + v(:, :, 2))
            end if
            !If nspden=4, the density is already given in the Pauli basis.
        end if
    end subroutine from_pauli

    !****f* m_precon/apply_vc
    !! NAME
    !! apply_vc
    !!
    !! FUNCTION
    !!
    !! INPUTS
    !!
    !! SOURCE
    subroutine apply_vc(this, ngfft, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, this%nspden)
        
        !Local variables-------------------------------
        integer :: ifft, ispden
        real(dp) :: g_cart_2
        
        ! *************************************************************************
        
        !In the sigma_0, 1, 2, 3 basis :
        !   The sigma_0 component of the density is multiplied by 4pi/G^2
        !   and the rest is 0.
        ispden = 1
        do ifft = 2, this%nfft
            g_cart_2 = norm2(two_pi * matmul(this%gprimd, get_g_vector(ifft, ngfft)))**2
            vec_g(:, ifft, ispden) = (2*two_pi/g_cart_2) * vec_g(:, ifft, ispden)
        end do

        do ispden = 2, this%nspden
            vec_g(:, :, ispden) = 0
        end do

    end subroutine apply_vc

    !****f* m_precon/apply_kxc
    !! NAME
    !! apply_kxc
    !!
    !! FUNCTION
    !!
    !! INPUTS
    !!
    !! SOURCE
    subroutine apply_kxc(this, dtset, mpi_enreg, ngfft, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, dtset%nspden)
        
        !Local variables-------------------------------
        !scalars
        integer :: cplex, n3xccc, nhat1dim, nhat1grdim, nkxc, option, usexcnhat
        logical :: non_magnetic_xc
        !arrays
        real(dp), allocatable :: Kxc_vec_r(:, :), vec_r(:, :)
        real(dp), allocatable :: vec_g_old(:, :, :) !DEBUG
        real(dp), allocatable :: nhat1(:, :), nhat1gr(:, :, :)
        real(dp) :: dummy_xccc3d1(0), dummy_qphon(3)
        
        ! *************************************************************************
        !write(6,*)'chi0diel apply_kxc : vec_g tot ', vec_g(1, 1:10, 1); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : vec_g spin ', vec_g(1, 1:10, 2); flush(6) !DEBUG
        ! TODO : remove all the write

        !Basis change to get vec_g in the default Abinit spin-basis to apply existing routines
        call from_pauli(this, 1, vec_g)
        ABI_MALLOC(vec_g_old, (2, this%nfft, dtset%nspden))    !DEBUG
        vec_g_old = vec_g   !DEBUG
        !write(6,*)'chi0diel apply_kxc : vec_g default1 ', vec_g(1, 1:10, 1); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : vec_g default2 ', vec_g(1, 1:10, 2); flush(6) !DEBUG

        !ifft to get vec_g in the real-space basis required by dfpt_mkvxc
        ABI_MALLOC(vec_r, (this%nfft, dtset%nspden))
        call fourdp(1, vec_g, vec_r, 1, mpi_enreg, this%nfft, dtset%nspden, ngfft, 0)
        !write(6,*)'chi0diel apply_kxc : vec_r up ', vec_r(1:10, 2); flush(6) !DEBUG

        !Applying Kxc : 
        cplex = 1   ! Input vector is real in real (direct) space.
        non_magnetic_xc = .false.
        nkxc = size(this%kxc, 2)
        !write(6,*)'chi0diel apply_kxc : nkxc ', nkxc; flush(6) !DEBUG

        usexcnhat = 0                                               !
        nhat1dim = 0                                                ! 
        ABI_MALLOC(nhat1, (cplex*this%nfft, dtset%nspden*nhat1dim)) ! PAW - TODO ?
        nhat1grdim = 0                                              !
        ABI_MALLOC(nhat1gr, (cplex*this%nfft, dtset%nspden, 3*nhat1grdim))  !

        option = 2  ! Treats only density change (no core_correction)
        n3xccc = 0  !   -> Core-correction set to 0.

        ABI_MALLOC(Kxc_vec_r, (this%nfft, dtset%nspden))
        !write(6,*)'chi0diel apply_kxc : this%kxc(1:5, 1)', this%kxc(1:5, 1); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : this%kxc(1:5, 2)', this%kxc(1:5, 2); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : this%kxc(1:5, 3)', this%kxc(1:5, 3); flush(6) !DEBUG

        call dfpt_mkvxc(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfft, ngfft, nhat1, nhat1dim, &
        &               nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden, n3xccc, option, &
        &               dummy_qphon, vec_r, this%rprimd, usexcnhat, Kxc_vec_r, dummy_xccc3d1)
        ABI_FREE(nhat1)
        ABI_FREE(nhat1gr)
        !write(6,*)'chi0diel apply_kxc : Kxc_vec_r up', Kxc_vec_r(1:10, 2); flush(6) !DEBUG

        !fft to return the result Kxc_vec_r in the plane-wave basis.
        call fourdp(1, vec_g, Kxc_vec_r, -1, mpi_enreg, this%nfft, dtset%nspden, ngfft, 0)
        !write(6,*)'chi0diel apply_kxc : Kxc_vec_g default1 ', vec_g(1, 1:10, 1); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : Kxc_vec_g default2 ', vec_g(1, 1:10, 2); flush(6) !DEBUG
        call this%save_applied_op_r(dtset, ngfft, vec_r, Kxc_vec_r, "applied_kxc_r.txt")  !DEBUG

        ABI_FREE(Kxc_vec_r)
        ABI_FREE(vec_r)

        call this%save_applied_op_g(dtset, ngfft, vec_g_old, vec_g, "applied_kxc.txt")  !DEBUG
        ABI_FREE(vec_g_old) !DEBUG

        !Basis change to get the resulting vector in the Pauli basis used in the rest of the preconditioning routines.
        call to_pauli(this, 0, vec_g)
        !write(6,*)'chi0diel apply_kxc : Kxc_vec_g tot ', vec_g(1, 1:10, 1); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_kxc : Kxc_vec_g spin ', vec_g(1, 1:10, 2); flush(6) !DEBUG

    end subroutine apply_kxc

    !****f* m_precon/apply_kernel
    !! NAME
    !! apply_vc
    !!
    !! FUNCTION
    !!
    !! INPUTS
    !!
    !! SOURCE
    subroutine apply_kernel(this, dtset, mpi_enreg, ngfft, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, dtset%nspden)
        
        !Local variables-------------------------------
        real(dp), allocatable :: Kxc_vec_g(:, :, :)

        ! *************************************************************************
        write(6,*)'chi0diel apply_kernel'; flush(6) !DEBUG
        
        !LDOS/Kerker model - only vc (RPA)
        if (this%iprcel == 201 .or. this%iprcel == 202) then
            call apply_vc(this, ngfft, vec_g)
        end if

        !Deigvals (name?) model - K * vec_g = vc * vec_tot + Kxc * vec_spin
        !   RPA on the total density
        !   Kxc restricted to the spin density
        if (this%iprcel == 203 .and. dtset%nspden==2) then
            write(6,*)'chi0diel apply_kernel iprcel 203'; flush(6) !DEBUG
            ABI_MALLOC(Kxc_vec_g, (2, this%nfft, dtset%nspden))
            Kxc_vec_g = vec_g   !Copy ?
            Kxc_vec_g(:, :, 1) = 0  ! Kxc only applied to the spin component
            call apply_vc(this, ngfft, vec_g)
            call apply_Kxc(this, dtset, mpi_enreg, ngfft, Kxc_vec_g)
            Kxc_vec_g(:, :, 1) = 0
            vec_g = vec_g + Kxc_vec_g
            ABI_FREE(Kxc_vec_g)
        end if

    end subroutine apply_kernel

    !****f* m_precon/derivative_occ
    !! NAME
    !! derivative_occ
    !!
    !! FUNCTION
    !! Computes the derivative of the occupation (f) of the band corresponding to eigenval 
    !! with respect to the fermie temperature.
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
            ABI_BUG("Non-metallic occupation")
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
            ABI_BUG("LDOS preconditioning needs a smooth smearing function")
        else if (occopt==9) then
        !Fermi-Dirac occupation is enforced with two distinct quasi-Fermi levels
            ABI_BUG("LDOS preconditioning not implemented for this smearing function")
        end if
    
        fprim = -1/tsmear * delta
        
    end function derivative_occ

    !****f* m_precon/compute_ldos
    !! NAME
    !! compute_ldos
    !!
    !! FUNCTION
    !!  Computes the local density of states defined as
    !!      ldos = sum_nk f'_nk |u_nk|^2 .
    !!  where f'_nk is the derivative of the occupation (nk) with respect to the fermi energie.
    !!
    !! INPUTS
    !!  occopt   = option for occupancies, determines delta
    !!  eigenval = eigenvalue
    !!  fermie   = fermi energie
    !!  tsmear   = smearing temperature
    !!  cg       = wf in G space
    !!  dtset    = all input variables for this dataset
    !!  ucvol    = unit cell volume
    !!  rprimd   = dimensional real space primitive translations
    !!  gprimd   = dimensional reciprocal space primitive translations
    !!  irrzon   = irreducible zone data
    !!  kg       = reduced planewave coordinates
    !!  npwarr   = number of planewaves and boundary planewaves at each k
    !!  phnons   = nonsymmorphic translation phases
    !!  nfft     = (effective) number of FFT grid points (for this processor)
    !!
    !! OUTPUT
    !!  ldos     = local density of state
    !!
    !! SOURCE
    subroutine compute_ldos(dtset, cg, eigen, fermie, gprimd, irrzon, &
        &   kg, mpi_enreg, nfft, npwarr, phnons, rprimd, ucvol, &
        &   ldos)

        !Arguments ------------------------------------
        !scalars
        integer, intent(in) :: nfft
        real(dp), intent(in) :: fermie
        !arrays
        real(dp), intent(in) :: eigen(:)
        real(dp), intent(out) :: ldos(:, :)
        !integer, intent(in) :: atindx(:), atindx1(:)
        !cprj, dimcprj, mband_cprj TODO 

        !mkrho arguments
        real(dp), intent(in) :: ucvol
        type(MPI_type), intent(inout) :: mpi_enreg
        type(dataset_type), intent(in) :: dtset
        integer, intent(in) :: irrzon(:, :, :)
        integer, intent(in) :: kg(:,:), npwarr(:)
        real(dp), intent(in) :: gprimd(:,:)
        real(dp), intent(in) :: cg(:,:)
        real(dp), intent(in) :: phnons(:, :, :)
        real(dp), intent(in) :: rprimd(3,3)
        
        !Local variables-------------------------------
        !scalars
        integer :: mband, i, mcg, maxocc, nfftot
        !arrays
        real(dp), allocatable :: ldos_weights(:)

        !dummy mkrho arguments
        real(dp)                :: rhog(2, nfft)
        type(paw_dmft_type)     :: paw_dmft
        type(wvl_wf_type)       :: wvl_wfs
        type(wvl_denspot_type)  :: wvl_den

        ! *************************************************************************

        !compute weights
        mband = size(eigen)
        ABI_MALLOC(ldos_weights, (mband))
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        do i=1, mband
            ldos_weights(i) = -derivative_occ(dtset%occopt, eigen(i), fermie, dtset%tsmear) * maxocc
        end do
        write(6,*)'chi0diel compute_ldos : mband, dtset%occopt, fermie, dtset%tsmear, maxocc', mband, dtset%occopt, fermie, dtset%tsmear, maxocc; flush(6) !DEBUG

        write(6,*)'chi0diel compute_ldos : ldos_weights', ldos_weights(1:10); flush(6) !DEBUG

        !Compute ldos using mkrho with ldos_weights in place of the occupations
        mcg = size(cg)
        paw_dmft%use_dmft = 0
        paw_dmft%use_sc_dmft = 0
        call mkrho(cg, dtset, gprimd, irrzon, kg, mcg, mpi_enreg, npwarr, ldos_weights, &
        &   paw_dmft, phnons, rhog, ldos, rprimd, 0, ucvol, wvl_den, wvl_wfs, option=0)
        ABI_FREE(ldos_weights)
        write(6,*)'chi0diel compute_ldos : ldos before sym', ldos(1:10, 1); flush(6) !DEBUG
        
        !TODO nfftmix != dtset%nfft en PAW grille
        nfftot = dtset%ngfft(1) * dtset%ngfft(2) * dtset%ngfft(3)
        call symrhg(1, gprimd, irrzon, mpi_enreg, nfft, nfftot, dtset%ngfft, dtset%nspden, dtset%nsppol, &
        &   dtset%nsym, phnons, rhog, ldos, rprimd, dtset%symafm, dtset%symrel, dtset%tnons)
        write(6,*)'chi0diel compute_ldos : ldos after sym', ldos(1:10, 1); flush(6) !DEBUG

        !if (psps%usepaw==1) then
        !    mcprj = size(cprj)
            !Computing the rhoij equivalent for the local density of states (ldos)
            !   Sum_{n,k} {ldos_weight(n,k)*<Cnk|p_i><p_j|Cnk>}.

            !1) Allocating pawrhoij object
        !    pawrhoij_ldos = ..
            !2) Computing
        !    call pawmkrhoij(atindx, atindx1, cprj, dimcprj, dtset%istwfk, dtset%kptopt, dtset%mband, mband_cprj, &
        !    &       mcprj, dtset%mkmem, mpi_enreg, dtset%natom, dtset%nband, dtset%nkpt, dtset%nspden, dtset%nspinor, &
        !    &       dtset%nsppol, ldos_weights, dtset%paral_kgb, paw_dmft, pawrhoij_ldos, dtfil%unpaw, dtset%usewvl, dtset%wtk)

            !Computing the total pseudo (compensated) ldos
        !    call pawmkrho(1, compch_fft, cplex, gprimd, idir, indsym, ipert, mpi_enreg, &
        !    &       my_natom, natom, dtset%nspden, dtset%nsym, dtset%ntypat, dtset%paral_kgb, pawang, pawfgr, pawfgrtab, &
        !    &       dtset%pawprtvol, pawrhoij_ldos, pawrhoij_ldos, pawtab, qpt, rhowfg, rhowfr, rhor, rprimd, dtset%symafm, &
        !    &       symrec, dtset%typat, ucvol, dtset%usewvl, xred, pawnhat=nhat)
        !end if

        !With collinear spins the ldos is not returned in the Pauli (tot/spin) basis by mkrho.
        if (dtset%nspden == 2) then
            !spin = 2*up - tot
            ldos(:, 2) = 2*ldos(:, 2) - ldos(:, 1)
        end if

    end subroutine compute_ldos

    !****f* m_precon/complex_mult
    !! NAME
    !! complex_mult
    !!
    !! FUNCTION
    !! Multiply to complex numbers given as size two arrays.
    !!
    !! SOURCE
    function complex_mult(z1, z2) result(z3)
        !Arguments ------------------------------------
        real(dp), intent(in) :: z1(2), z2(2)
        
        !Returned variable-------------------------------
        real(dp) :: z3(2)
                
        ! *************************************************************************
        
        !Performs z3 = z1*z2, the complex multiplication
        z3(1) = z1(1)*z2(1)-z1(2)*z2(2)
        z3(2) = z1(1)*z2(2)+z1(2)*z2(1)

    end function complex_mult

    !****f* m_precon/compute_loc_pola
    !! NAME
    !! compute_loc_pola
    !!
    !! FUNCTION
    !! computes the local polarizability estimate
    !!
    !! INPUTS
    !!
    !! SIDE EFFECTS
    !!
    !! SOURCE
    subroutine compute_loc_pola(dtset, atindx1, gprimd, nattyp, nfft, &
        &   nspden, mpi_enreg, rhor, rprimd, xred, &
        &   loc_pola)

        !Arguments ------------------------------------
        !scalars
        type(dataset_type),intent(in) :: dtset
        integer, intent(in) :: nfft, nspden
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: atindx1(:), nattyp(:)
        real(dp), intent(in) :: rhor(:, :), rprimd(3, 3), gprimd(3, 3), xred(:, :)
        real(dp), intent(out) :: loc_pola(:, :)
       
        !Local variables-------------------------------
        !scalars
        integer :: itypat, iattyp, iatom, ig, ir, ispden
        integer :: re, im
        real(dp) :: l_atom
        !arrays
        integer :: g(3)
        real(dp) :: form_factor(2), structure_factor(2)
        real(dp) :: r_atom(3)
        real(dp), allocatable:: rhor0(:), rhor0_atom(:)
        real(dp), allocatable :: r2(:)
        real(dp), allocatable:: rhog0_atom(:, :)
        
        ! *************************************************************************
        
        loc_pola = 0.0

        ! Arrays allocation
        ABI_MALLOC(rhor0, (nfft))
        ABI_MALLOC(rhor0_atom, (nfft))
        ABI_MALLOC(r2, (nfft))
        ABI_MALLOC(rhog0_atom, (2, nfft))
        re = 1
        im = 2

        do itypat = 1, dtset%ntypat !Loop over the types of atom
            l_atom = atom_length(dtset%densty(itypat, 1), dtset%ziontypat(itypat), dtset%znucl(itypat))   ! Atomic decay length
            
            do iattyp = 1, nattyp(itypat) !Loop over the atom (of this type)
                iatom = atindx1(iattyp)
                r_atom = xred(:, iatom) !in reduced coordinates
                
                !1) Computing rhor0_atom (with Gaussians)
                do ig = 1, nfft !Loop over the fft grid
                    g = get_g_vector(ig, dtset%ngfft) !in reduced coordinates
                    ! structure_factor(g) = exp(-i*2pi*dot(g, r_atom))
                    structure_factor(re) = cos(-two_pi*dot_product(g, r_atom))
                    structure_factor(im) = sin(-two_pi*dot_product(g, r_atom))
                    ! form_factor(g) = exp(-(2pi*l_atom*g)^2) (Gaussian)
                    form_factor(re) = exp(-(two_pi*l_atom*norm2(matmul(gprimd, g)))**2)
                    form_factor(im) = 0 
                    ! rhog0_atom = structure_factor * form_factor (multiplication in g-space)
                    rhog0_atom(:, ig) = complex_mult(form_factor, structure_factor)
                end do
                !ifft to get rhor0_atom in real space
                call fourdp(1, rhog0_atom, rhor0_atom, 1, mpi_enreg, nfft, 1, dtset%ngfft, 0) ! irfft (cplex=1)

                !2) Computing norm(r-r_atom)^2 (periodized with sin)
                !TODO : better + mistake (pic dans loc_pola)
                do ir = 1, nfft
                    r2(ir) = norm2(matmul(rprimd, 1/two_pi*sin(two_pi * (get_r_vector(ir, dtset%ngfft) - r_atom)))) ** 2
                end do

                !3) loc_pola = sum_atom rho_atom * |r-r_atom|^2
                do ispden = 1, nspden
                    ! TODO : add atom-specific coefficients
                    loc_pola(:, ispden) = loc_pola(:, ispden) * rhor0/(rhor0+rhor0_atom) + rhor(:, ispden) * rhor0_atom/(rhor0+rhor0_atom) * r2
                end do
                rhor0 = rhor0 + rhor0_atom
            
            end do 
        end do

        ABI_FREE(rhor0)
        ABI_FREE(rhor0_atom)
        ABI_FREE(r2)
        ABI_FREE(rhog0_atom)

    end subroutine compute_loc_pola

    !****f* m_precon/save_applied_op_g
    !! NAME
    !! save_applied_op_g
    !!
    !! FUNCTION
    !! Save the applied operator vec/op_vec in the reciprocal (G) space in a file.
    !!
    !! SOURCE
    subroutine save_applied_op_g(this, dtset, ngfft, vec, op_vec, filename)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        character(len=*), intent(in) :: filename ! Filename for the output file
        !arrays
        real(dp), intent(in) :: vec(:, :, :), op_vec(:, :, :)
        integer, intent(in) :: ngfft(:)
       
        !Local variables-------------------------------
        !scalars
        integer :: io, n, i, ispden
        
        ! *************************************************************************
        n = size(vec, 2)
        ! Writing the file
        open(newunit=io, file=filename, status="replace", action="write")
            do ispden=1, dtset%nspden
                do i=1, n
                    write (io, '(*(G0.6,:,","))') two_pi*matmul(this%gprimd, get_g_vector(i, ngfft)),   &
                    &                             vec(1, i, ispden), vec(2, i, ispden),                 &
                    &                             op_vec(1, i, ispden), op_vec(2, i, ispden)
                end do
            end do
        close(io)

    end subroutine save_applied_op_g

    !****f* m_precon/save_applied_op_r
    !! NAME
    !! save_applied_op_r
    !!
    !! FUNCTION
    !! Save the applied operator vec/op_vec in the direct (real) space in a file.
    !!
    !! SOURCE
    subroutine save_applied_op_r(this, dtset, ngfft, vec, op_vec, filename)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        character(len=*), intent(in) :: filename ! Filename for the output file
        !arrays
        real(dp), intent(in) :: vec(:, :), op_vec(:, :)
        integer, intent(in) :: ngfft(:)
       
        !Local variables-------------------------------
        !scalars
        integer :: io, n, i, ispden
        
        ! *************************************************************************
        n = size(vec, 1)
        ! Writing the file
        open(newunit=io, file=filename, status="replace", action="write")
            do ispden=1, dtset%nspden
                do i=1, n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)),      &
                    &                             vec(i, ispden), op_vec(i, ispden)                
                end do
            end do
        close(io)

    end subroutine save_applied_op_r

    !****f* m_precon/apply_chi0_ldos
    !! NAME
    !! apply_chi0_ldos
    !!
    !! FUNCTION
    !! Applies the ldos-model chi0 operator to the vector vec_g (in place).
    !!
    !! INPUTS
    !!  mpi_enreg    = Information about MPI parallelization.
    !!  ngfft        = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  ispden       = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_g (nspden, 2, :) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                         When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0_ldos(this, mpi_enreg, ngfft, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, this%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden
        !arrays
        real(dp), allocatable :: work_r(:), vec_r_1(:)
        real(dp), allocatable :: work_g(:, :)
        
        ! *************************************************************************
       
        if (abs(this%tdos) > epsilon(this%tdos)) then   !Checking that tdos is not 0.
            cplex = 1
            ABI_MALLOC(vec_r_1, (cplex*this%nfft))    
            ABI_MALLOC(work_r, (cplex*this%nfft))
            
            !1) ifft to get (the sigma_0 component of) vec in the real space
            ispden = 1
            call fourdp(cplex, vec_g(:, :, ispden), vec_r_1, 1, mpi_enreg, this%nfft, 1, ngfft, 0)

            do ispden = 1, this%nspden
                !2) chi0(v)(r)_ispden = -ldos_ispden(r)*v_1(r) + 1/dos * ldos_ispden(r)*integral(ldos_1(r')*v_1(r')*dr')
                work_r = -this%ldos(:, ispden)*vec_r_1 &
                                + 1/this%tdos * dot_product(this%ldos(:, 1), vec_r_1)*this%dvol * this%ldos(:, ispden)

                !3) fft to get vec back in the reciprocal space
                call fourdp(cplex, vec_g(:, :, ispden), work_r, -1, mpi_enreg, this%nfft, 1, ngfft, 0)
            end do

            ABI_FREE(work_r)
            ABI_FREE(vec_r_1)
        else 
            vec_g = 0
        end if

    end subroutine apply_chi0_ldos

    !****f* m_precon/apply_chi0_locpola
    !! NAME
    !! apply_chi0_ldos
    !!
    !! FUNCTION
    !! Applies the local polarizability model chi0 operator to the vector vec_g (in place).
    !! TODO
    !!
    !! INPUTS
    !!  mpi_enreg    = Information about MPI parallelization.
    !!  ngfft        = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  ispden       = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_g (nspden, 2, :) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                         When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0_locpola(this, mpi_enreg, ngfft, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, this%nspden)
       
        !Local variables-------------------------------
        !scalars
        !arrays
        
        ! *************************************************************************
       
        ! TODO : spin
        !    ABI_MALLOC(work_g, (2, this%nfft))
        !    ABI_MALLOC(vec_g_saved, (2, this%nfft))
        !    vec_g_saved = vec_g
        !    vec_g = 0.0
        !    do i=1, 3
        !        work_g = 0.0
        !        !1) Multiplication by d_i(g) in reciprocal space
        !        do i_g = 1, this%nfft
        !            !g = two_pi * matmul(this%gprimd, g_vectors(:, i_g))
        !            g = two_pi * matmul(this%gprimd, get_g_vector(i_g, ngfft))
        !            work_g(:, i_g) = complex_mult( [0.0_dp, g(i)/sqrt(1+(norm2(g)/this%gc)**2)], vec_g_saved(:, i_g) )
        !        end do

        !        !2) Multiplication by the local polarizability in real space
        !        call fourdp(cplex, work_g, work_r, 1, mpi_enreg, this%nfft, 1, ngfft, 0) !ifft
        !        work_r = this%loc_pola(:, ispden) * work_r                                         !local multiplication
        !        call fourdp(cplex, work_g, work_r, -1, mpi_enreg, this%nfft, 1, ngfft, 0) !fft

        !        !3) Multiplication by d_i(g) in reciprocal space
        !        do i_g = 1, this%nfft
        !            g = two_pi * matmul(this%gprimd, get_g_vector(i_g, ngfft))
        !            work_g(:, i_g) = complex_mult( [0.0_dp, g(i)/sqrt(1+(norm2(g)/this%gc)**2)], work_g(:, i_g) )
        !        end do

        !        vec_g = vec_g + work_g

        !    end do

    end subroutine apply_chi0_locpola

!****f* m_precon/compute_weights_chi0_diag
    !! NAME
    !! compute_weights_chi0_diag
    !!
    !! FUNCTION
    !! Computes the weights needed to compute the application od the model chi0 to a vector (vec_g)
    !! such that chi0 * vec = sum_i weight_i * rho_ii
    !! that is 
    !!      weight_i = f'i * dot(vec, rho_ii).
    !!
    !! INPUTS
    !!  dtset       =
    !!  mgfft       =
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  ngfft       = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  vec_g (2, nfft, nspden)      = Vector (in G-space) to which the model chi0 operator is to be applied.
    !!                              When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! OuTPUTS
    !!  weights(:) = ...
    !!                              When nspden > 1 chi0_vec_g is in the Pauli basis.
    !!
    !! SOURCE
    ! TODO : Name 
    subroutine compute_weights_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfft, vec_g, weights)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer,intent(in) :: mgfft
        !arrays
        integer, intent(in) :: ngfft(:) 
        real(dp), intent(inout) :: vec_g(2, this%nfft, dtset%nspden)
        real(dp), intent(inout) :: weights(:)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden, istwf_k, maxocc, mcg, my_nspinor, mband_mem, nband_k, ndat, nfftot, npw_k, option, tim_fourwf
        integer :: i_eigen, i_cg, j_cg, ikg, ikpt, iband, isppol, ier
        integer :: n1, n2, n3, n4, n5, n6
        real(dp) :: fp, eigenval
        integer :: dummy_int
        real(dp) :: dummy_real
        !arrays
        integer :: gbound(2*dtset%mgfft+8,2)
        integer, allocatable :: kg_k(:, :)
        real(dp), allocatable :: vec_r(:, :)
        real(dp), allocatable :: rhoaug_r_ii(:, :, :, :), rho_r_ii(:, :)
        real(dp) :: vec_g_old(2, this%nfft, dtset%nspden)   !DEBUG

        !dummy fourwfk arguments
        real(dp), allocatable ::  dummy_denpot(:, :, :), dummy_fofgout(:, :)
        
        ! *************************************************************************
        write(6,*)'chi0diel apply_chi0_mag'; flush(6) !DEBUG
        vec_g_old = vec_g   !DEBUG
       
        n1 = ngfft(1)
        n2 = ngfft(2)
        n3 = ngfft(3)
        n4 = ngfft(4)
        n5 = ngfft(5)
        n6 = ngfft(6)
        if (this%nfft /= n1*n2*n3) then
            ABI_BUG("Mismatch nfft != n1*n2*n3")
        end if

        !1) Computing vec_r = input vector in real space
        !   - in the Pauli basis if nspinor = 2
        !   - in the up/down basis if nspinor = 1
        if (dtset%nspinor == 1) then
            call from_pauli(this, 0, vec_g)
        end if
        ABI_MALLOC(vec_r, (this%nfft, dtset%nspden))
        write(6,*)'chi0diel apply_chi0_mag before fourdp vec_g(1, 1:10, 2): ', vec_g(1, 1:10, 2); flush(6) !DEBUG
        call fourdp(1, vec_g, vec_r, 1, mpi_enreg, this%nfft, 2, ngfft, 0)
        
        !2) Computing the weights = fi' * <rhoii, vec_spin> TODO
        weights = 0
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        write(6,*)'chi0diel apply_chi0_mag maxocc, this%dvol: ', maxocc, this%dvol; flush(6) !DEBUG
        
        !Allocating the array that will contain rho_ii
        ABI_MALLOC(rhoaug_r_ii, (2, n4, n5, n6))
        if (dtset%nspinor==1) then
            ABI_MALLOC(rho_r_ii, (n1*n2*n3, 1))
        else if (dtset%nspinor==2) then
            ABI_MALLOC(rho_r_ii, (n1*n2*n3, 4))
        else
            ABI_BUG("nspinor != 1 or 2")
        end if

        !Allocating dummy arrays
        ABI_MALLOC(dummy_denpot, (0, n5, n6))
        ABI_MALLOC(dummy_fofgout, (2, 0))

        i_cg = 0    ! Starting index for (ikpt, isppol) in cg array.
        my_nspinor = max(1, dtset%nspinor/mpi_enreg%nproc_spinor)

        do isppol =1, dtset%nsppol
            ikg = 0 ! Starting index for ikpt in kg.
            
            do ikpt = 1, dtset%nkpt

                nband_k = dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
                mband_mem = nband_k !
                if (dtset%paral_kgb==0) then
                    mband_mem = nband_k/mpi_enreg%nproc_band
                end if

                !MPI parallelization over kpoints : cycle if kpt does not belong to current processor.
                if (proc_distrb_cycle(mpi_enreg%proc_distrb, ikpt, 1, nband_k, isppol, mpi_enreg%me_kpt)) then
                    cycle
                end if

                npw_k = this%npwarr(ikpt)       ! Number of plane-wave at this kpt.
                istwf_k = dtset%istwfk(ikpt)    ! Option parameter that describes the storage of wfs at this kpt.
                ABI_MALLOC(kg_k, (3, npw_k))
                kg_k = this%kg(:, 1+ikg:npw_k+ikg)     ! Reduced plane-wave coordinate (k+G) of this kpt.
                write(6,*)'chi0diel apply_chi0_mag sphereboundary'; flush(6) !DEBUG

                call sphereboundary(gbound, istwf_k, kg_k, dtset%mgfft, npw_k)    ! Computes gbound.
                
                write(6,*)'chi0diel apply_chi0_mag sphereboundary done'; flush(6) !DEBUG

                do iband = 1, nband_k

                    !Indices
                    i_eigen = iband + (ikpt-1)*dtset%mband + (isppol-1)*dtset%mband*dtset%nkpt  ! Index of (iband, ikpt, isppol) in eigen array.
                    j_cg = i_cg + (iband-1) * npw_k * my_nspinor    ! Index of (iband, ikpt, isppol) in cg array.                    
                    
                    !2.1) Computing f'(eig_i - fermie).
                    eigenval = this%eigen(i_eigen)
                    fp = derivative_occ(dtset%occopt, eigenval, this%fermie, dtset%tsmear) * maxocc

                    if (abs(fp) > tol14) then   !TODO : choose tol ?
                        
                        !2.2) Computing rho_ii = |psi_i|^2 using fourwf (if fp is not 0).
                        
                        !Input parameters for fourwf :
                        option = 0      ! FFT from reciprocal to direct space.
                        ndat = 1        ! Only one FFT. TODO : Do them all in once ? 
                        tim_fourwf = 0

                        call fourwf(dummy_int, dummy_denpot, this%cg(:, j_cg+1:j_cg+npw_k), dummy_fofgout, rhoaug_r_ii,  &
                        &           gbound, gbound, istwf_k, kg_k, kg_k, mgfft, mpi_enreg, ndat, ngfft, npw_k, &
                        &           dummy_int, n4, n5, n6, option, tim_fourwf, dummy_real, dummy_real)

                        rhoaug_r_ii(1, :, :, :) = rhoaug_r_ii(1, :, :, :)**2 + rhoaug_r_ii(2, :, :, :)**2  !|.|^2
                        
                        if (dtset%nspinor == 2) then ! Non collinear spins.
                            do ispden=1, 4
                                !Transfer the rhoaug_r_ii defined on the large fft-grid to the smaller density/potential fft-grid.
                                call fftpac(ispden, mpi_enreg, 4, n1, n2, n3, n4, n5, n6, ngfft, rho_r_ii, rhoaug_r_ii(1, :, :, :), 1)
                                rho_r_ii(:, ispden) = rho_r_ii(:, ispden) / (sum(rho_r_ii(:, ispden)) * this%dvol) !Normalizing rho_ii_r.
                                ! dot product in Pauli basis :
                                weights(i_eigen) = weights(i_eigen) + fp * maxocc * dot_product(rho_r_ii(:, ispden), vec_r(:, ispden)) * this%dvol
                            end do
                        end if

                        if (dtset%nspinor == 2) then ! Collinear spins or no spin.
                            call fftpac(1, mpi_enreg, 1, n1, n2, n3, n4, n5, n6, ngfft, rho_r_ii, rhoaug_r_ii(1, :, :, :), 1)
                            rho_r_ii(:, 1) = rho_r_ii(:, 1) / (sum(rho_r_ii(:, 1)) * this%dvol) !Normalizing rho_ii_r.
                            weights(i_eigen) = fp * maxocc * dot_product(rho_r_ii(:, 1), vec_r(:, isppol)) * this%dvol
                        end if

                        if (ikpt==1) then
                        write(6,*)'--- chi0diel apply_chi0_mag:  i_eigen, iband, ikpt, isppol', i_eigen, iband, ikpt, isppol; flush(6) !DEBUG
                        write(6,*)'chi0diel apply_chi0_mag: eigenval', eigenval; flush(6) !DEBUG
                        write(6,*)'chi0diel apply_chi0_mag: fp', fp; flush(6) !DEBUG
                        write(6,*)'chi0diel apply_chi0_mag: rho_r_ii(1:5, 1)', rho_r_ii(1:5, 1); flush(6) !DEBUG
                        write(6,*)'chi0diel apply_chi0_mag: sum(rho_r_ii(:, 1))', sum(rho_r_ii(:, 1)); flush(6) !DEBUG
                        write(6,*)'chi0diel apply_chi0_mag: dot_product(rho_r_ii(:, 1), vec_r(:, 2))*this%dvol', dot_product(rho_r_ii(:, 1), vec_r(:, 2))*this%dvol; flush(6) !DEBUG
                        end if

                    end if
                    
                end do
                write(6,*)'chi0diel apply_chi0_mag:  i_cg, isppol', i_cg, isppol; flush(6) !DEBUG

                if (dtset%mkmem /= 0) then
                    i_cg = i_cg + npw_k * my_nspinor * mband_mem
                    ikg = ikg + npw_k
                !else ???
                end if
                ABI_FREE(kg_k)

            end do  !ikpt
        end do  !isppol

        ABI_FREE(rhoaug_r_ii)
        ABI_FREE(rho_r_ii)
        ABI_FREE(vec_r)
        ABI_FREE(dummy_denpot)
        ABI_FREE(dummy_fofgout)

        !MPI parallelization over kpoints : sum weights on all processors.
        ier = 0
        call xmpi_sum(weights, mpi_enreg%comm_kpt, ier)
        
        write(6,*)'chi0diel apply_chi0_mag weights(16:17)', weights(16:17); flush(6) !DEBUG
        write(6,*)'chi0diel apply_chi0_mag weights(48:49)', weights(48:49); flush(6) !DEBUG
        write(6,*)'chi0diel apply_chi0_mag weights', weights; flush(6) !DEBUG

    end subroutine compute_weights_chi0_diag

    !****f* m_precon/apply_chi0_diag
    !! NAME
    !! apply_chi0_mag
    !!
    !! FUNCTION
    !! Applies the model ... of the chi0 operator to the vector vec_g
    !!
    !! INPUTS
    !!  dtset
    !!  mgfft
    !!  mpi_enreg    = Information about MPI parallelization.
    !!  ngfft        = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  vec_g (2, nfft, nspden)      = Vector (in G-space) to which the model chi0 operator is applied.
    !!                              When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! OuTPUTS
    !!  chi0_vec_g (2, nfft, nspden) = ...
    !!                              When nspden > 1 chi0_vec_g is in the Pauli basis.
    !!
    !! SOURCE
    ! TODO : Name 
    subroutine apply_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfft, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer,intent(in) :: mgfft
        !arrays
        integer, intent(in) :: ngfft(:) 
        real(dp), intent(inout) :: vec_g(2, this%nfft, dtset%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden, istwf_k, maxocc, mcg, my_nspinor, mband_mem, nband_k, ndat, nfftot, npw_k, option, tim_fourwf
        integer :: i_eigen, i_cg, j_cg, ikg, ikpt, iband, isppol, ier, sign_isppol
        integer :: n1, n2, n3, n4, n5, n6
        real(dp) :: fp, eigenval
        integer :: dummy_int
        real(dp) :: dummy_real
        !arrays
        real(dp), allocatable :: chi0_vec_r(:, :)
        real(dp), allocatable :: weights(:)
        real(dp) :: vec_g_old(2, this%nfft, dtset%nspden)   !DEBUG

        !dummy mkrho arguments
        type(paw_dmft_type)     :: paw_dmft
        type(wvl_wf_type)       :: wvl_wfs
        type(wvl_denspot_type)  :: wvl_den
        
        ! *************************************************************************
        write(6,*)'chi0diel apply_chi0_mag'; flush(6) !DEBUG
        vec_g_old = vec_g   !DEBUG

        !1) Computing the weights = sum_i fi' * <rhoii, vec> TODO
        ABI_MALLOC(weights, (dtset%mband*dtset%nkpt*dtset%nsppol))
        call compute_weights_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfft, vec_g, weights)
        
        !2) Computing chi0 * vec using mkrho with custom the weights in place of the occupations.
        ! In place : vec_g = chi0 * vec_g
        ABI_MALLOC(chi0_vec_r, (this%nfft, dtset%nspden))
        mcg = size(this%cg, 2)
        paw_dmft%use_dmft = 0
        paw_dmft%use_sc_dmft = 0
        call mkrho(this%cg, dtset, this%gprimd, this%irrzon, this%kg, mcg, mpi_enreg, this%npwarr, weights, &
        &   paw_dmft, this%phnons, vec_g, chi0_vec_r, this%rprimd, 0, this%ucvol, wvl_den, wvl_wfs, option=0)
        ABI_FREE(weights)
        
        ! TODO :  nfftmix != this%nfft en PAW grille
        nfftot = ngfft(1) * ngfft(2) * ngfft(3)
        call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, this%nfft, nfftot, ngfft, dtset%nspden, dtset%nsppol, &
        &   dtset%nsym, this%phnons, vec_g, chi0_vec_r, this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)
        ! TODO : PAW
        !call this%save_applied_op_r(dtset, ngfft, vec_r, chi0_vec_r, "applied_chi0_mag_r.txt")  !DEBUG

        ABI_FREE(chi0_vec_r)
        write(6,*)'chi0diel apply_chi0_mag vec_g default1', vec_g(1, 1:10,  1); flush(6) !DEBUG
        write(6,*)'chi0diel apply_chi0_mag vec_g default2', vec_g(1, 1:10,  2); flush(6) !DEBUG

        !Basis change : With collinear spins vec_g is not returned in the Pauli (tot/spin) basis by mkrho.
        call to_pauli(this, 1, vec_g)
        write(6,*)'chi0diel apply_chi0_mag vec_g tot', vec_g(1, 1:10,  1); flush(6) !DEBUG
        write(6,*)'chi0diel apply_chi0_mag vec_g spin', vec_g(1, 1:10,  2); flush(6) !DEBUG

        call this%save_applied_op_g(dtset, ngfft, vec_g_old, vec_g, "applied_chi0_mag.txt")  !DEBUG

    end subroutine apply_chi0_diag

    !****f* m_precon/apply_chi0
    !! NAME
    !! apply_chi0
    !!
    !! FUNCTION
    !! Applies the model chi0 operator to the vector vec_g (in place)
    !!
    !! INPUTS
    !!  mpi_enreg    = Information about MPI parallelization.
    !!  ngfft        = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  ispden       = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_g (nspden, 2, :) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                         When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0(this, dtset, mgfft, mpi_enreg, ngfft, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: mgfft
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft, this%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden
        integer :: i, i_g
        !arrays
        real(dp), allocatable :: work_r(:), vec_r_1(:)
        real(dp), allocatable :: work_g(:, :), vec_g_saved(:, :), vec_g_old(:, :, :) !DEBUG
        real(dp) :: g(3)
        
        ! *************************************************************************
        write(6,*)'chi0diel apply_chi0'; flush(6) !DEBUG
       
        !Kerker
        if (this%iprcel == 201) then
            ispden = 1
            vec_g(:, :, ispden) = (-1/(4*pi*(this%dielng)**2)) * vec_g(:, :, ispden)
            do ispden = 2, this%nspden
                vec_g(:, :, ispden) = 0
            end do
        end if
        
        !LDOS model
        if (this%iprcel == 202) then
            call apply_chi0_ldos(this, mpi_enreg, ngfft, vec_g)
        end if

        !Deigvals (name?) model
        if (this%iprcel == 203) then
            call apply_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfft, vec_g)
        end if

        !Local polarizability model - Not implemented
        if (this%iprcel == 204) then
            call apply_chi0_locpola(this, mpi_enreg, ngfft, vec_g)
        end if

    end subroutine apply_chi0

!-------------------------------------------------------------------------------------------
!This should be moved to a specific file
!-------------------------------------------------------------------------------------------

    subroutine call_FGMRES(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp),intent(in) :: rhs(:)
        real(dp),intent(inout) :: est(:)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        !Local variables-------------------------------
        !MKL FGMRES
        integer :: RCI_request, itercount, size_vres
        integer :: ipar(128)
        real(dp) :: dpar(128)
        real(dp), allocatable :: tmp(:)

        ! *************************************************************************
        
        !FGMRES initialization

        ABI_MALLOC(tmp, ((2*gmres_maxiter+1)*n + gmres_maxiter*(gmres_maxiter+9)/2 + 1))
        call dfgmres_init(n, est, rhs, RCI_request, ipar, dpar, tmp)
        !setting FGMRES parameters
        ipar(7) = 0              ! control verbosity : no warning message
        ipar(5) = gmres_maxiter  ! maximum number of iterations
        ipar(8) = 1              ! dfgmres routine performs the stopping test for the maximum number of iterations ipar(4)≤ipar(5)
        ipar(9) = 1              ! dfgmres routine performs the residual stopping test dpar(5)≤dpar(4)=dpar(1)*dpar(3)+dpar(2)
        ipar(10) = 0             ! no user defined stopping tests
        ipar(11) = 0             ! non-preconditioned GMRES
        ipar(12) = 1             ! dfgmres routine performs the automatic test dpar(7)≤dpar(8)
        ipar(15) = gmres_maxiter ! number of the non-restarted FGMRES iterations (no restart here)
        dpar(1) = gmres_rtol     ! relative tolerance
        ! dpar(2) = 0.01          ! absolute tolerance  (DFTK default=0.01)
        
        !FGMRES iterations
        
        call dfgmres_check(n, est, rhs, RCI_request, ipar, dpar, tmp)
        call dfgmres(n, est, rhs, RCI_request, ipar, dpar, tmp)
        
        do
            if (RCI_request==-1) then
            !    maximum number of iterations is reached
                call dfgmres_get(n, est, rhs, RCI_request, ipar, dpar, tmp, itercount)
                exit
            else if (RCI_request==0) then
            !    successful completion of the task
                call dfgmres_get(n, est, rhs, RCI_request, ipar, dpar, tmp, itercount)
                exit
            else  if (RCI_request==1) then
            !    multiply the matrix P by tmp(ipar(22)) and put the result in tmp(ipar(23))
                call matvec(n, tmp(ipar(22):ipar(22)+2*size_vres-1), tmp(ipar(23):ipar(23)+2*size_vres-1))
            !    proceed with FGMRES iterations
                call dfgmres(2*size_vres, est, rhs, RCI_request, ipar, dpar, tmp)
        !---------------------------------------------------------------------
        !  FGMRES Errors
            else if (RCI_request==-10) then
                ABI_BUG('FGMRES : attempt to divide by zero')
                exit
            else if (RCI_request==-11) then
                ABI_BUG('FGMRES : infinite cycle')
                exit
            else if (RCI_request==-12) then
                ABI_BUG('FGMRES : errors were found in the method parameters')
                exit
            ! RCI_request = 2, 3, 4 should not happen with this choice of parameters
            else
                ABI_BUG('FGMRES : RCI_request has unexpected value')
            end if
        !---------------------------------------------------------------------
        end do
    end subroutine call_FGMRES

    subroutine call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp),intent(in) :: rhs(n)
        real(dp),intent(inout) :: est(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
        !Local variables-------------------------------
        integer :: its, info, m
        real(dp) :: res, del
        real(dp), allocatable :: h(:, :), v(:, :)

        ! *************************************************************************

        m = gmres_maxiter
        ABI_MALLOC(h, (m+1, m))
        ABI_MALLOC(v, (n, m+1))
        res = gmres_rtol
        del = 0
        its = gmres_maxiter  ! No restart
        info = 1
        call gmresm(m, n, est, rhs, matvec, psolve, dotprd, h, v, res, del, its, info)
        ABI_FREE(h)
        ABI_FREE(v)

        contains
        ! Dummy :  No preconditioning
        subroutine psolve(n_, x)
            integer, intent(in) :: n_
            real(dp), intent(inout) :: x
        end subroutine psolve
        ! Dot product
        function dotprd(n_, a, b) result(c)
            integer, intent(in) :: n_
            real(dp), intent(inout) :: a(n_), b(n_)
            real(dp) :: c
            ! ***********************
            c = dot_product(a, b)
        end function dotprd

    end subroutine call_gmresm

    subroutine linsolve(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
        !Arguments ------------------------------------
        integer, intent(in) :: n, gmres_maxiter
        real(dp), intent(in) :: gmres_rtol
        real(dp), intent(in) :: rhs(n)
        real(dp), intent(inout) :: est(n)
        interface
            subroutine matvec(n_, x, y)
                integer, intent(in) :: n_
                double precision, intent(inout), target :: x(n_), y(n_)
            end subroutine matvec
        end interface
      
        ! *************************************************************************
        
        !TODO : dirty check of MKL availability
#if defined HAVE_LINALG_MKL_OMATCOPY
        write(6,*)'using FGMRES'; flush(6) !DEBUG
        call call_FGMRES(n, matvec, rhs, est, gmres_maxiter, gmres_rtol)
#else
        write(6,*)'using gmresm'; flush(6) !DEBUG
        call call_gmresm(n, matvec, est, rhs, gmres_maxiter, gmres_rtol)
#endif
        write(6,*)'chi0diel : gmres done'; flush(6) !DEBUG
      
    end subroutine linsolve

!----------------------------------------------------------------------
! Openpipeflow.org.  If used in your work, please cite
! Willis, A. (2017) SoftwareX 6, 124-127.
! https://doi.org/10.1016/j.softx.2017.05.003 (open access)
!                                      Thanks in advance! Ashley 2019.
!----------------------------------------------------------------------
! solve A x = b for x ;  
! minimise |Ax-b| subject to constraint |x| < delta .
! requires lapack routines dgelsy, dgesvd.
!----------------------------------------------------------------------
! m	  gmres dimension
! n 	  dimension of x
! x	  on input:  guess for x, can be 0
!         on exit:  solution x, subject to constraint if del>0
! b	  input b
! matvec  performs y := A x, call matvec(N,x, y)
! psolve  preconditioner, solve M x_out = x_in, call psolve(N,x)
! dotprd  dot product, d = dotprd(n,a,b)
! h       Hessian matrix,  size (m+1)*m
! v       Krylov subspace, size n*(m+1)
! res	  on input: |Ax-b|/|b|<res; 
!         on exit:  residual reached
! del     on input: if(del>0) then the x returned is the hookstep
!         on exit:  norm of next b predicted by hook
! its	  on input: max num its; 
!         on exit:  number of its taken
! info	  on input: if(info==1) print* residuals
!                   if(info==2) recalc hookstep with new del
! 	  on exit:  0 sucessful, 1 method breakdown, 2 max its
!							A.P.Willis 2008
!----------------------------------------------------------------------

 subroutine gmresm(m,n,x,b,matvec,psolve,dotprd,h,v,res,del,its,info)
   implicit none
   integer,          intent(in)    :: m
   integer,          intent(in)    :: n
   double precision, intent(inout) :: x(n)
   double precision, intent(in)    :: b(n)
   external                        :: matvec,psolve
   double precision, external      :: dotprd
   double precision, intent(inout) :: h(m+1,m)
   double precision, intent(inout) :: v(n,m+1)
   double precision, intent(inout) :: res
   double precision, intent(inout) :: del
   integer,          intent(inout) :: its
   integer,          intent(inout) :: info
   double precision :: tol,res_,stgn, w(n), z(n)
   double precision :: h_(m+1,m), y(m+1), p(m+1), work(4*m+1)
   integer :: imx, piv(m), rank, i
   double precision, save :: beta
   integer, save :: j
   logical :: done   

   if(info==2) then
      call hookstep(j,h,m,beta,del, y)
      z = matmul(v(:,1:j),y(1:j))
      call psolve(n, z)
      x = z
      info = 0
      return
   end if	 

   tol = res
   imx = its
   its = 0
   v   = 0d0

 1 continue
   res_ = 1d99
   stgn = 1d0 - 1d-14
 
   beta = dsqrt(dotprd(n,x,x)) 
   if(beta==0d0)  w = 0d0
   if(beta/=0d0)  call matvec(n,x, w)
   w = b - w
   beta = dsqrt(dotprd(n,w,w)) 
   v(:,1) = w / beta
     
   h = 0d0
   do j = 1, m
      its = its + 1
      z = v(:,j)      
      call psolve(n, z)
      call matvec(n, z, w)
      do i = 1, j
         h(i,j) = dotprd(n,w,v(1,i))
         w = w - h(i,j)*v(:,i)
      end do
      h(j+1,j) = dsqrt(dotprd(n,w,w))
      v(:,j+1) = w / h(j+1,j)
          
      p(1) = beta
      p(2:j+1) = 0d0
      h_(1:j+1,1:j) = h(1:j+1,1:j)
      call dgelsy(j+1,j,1,h_,m+1,p,m+1,piv,m,rank,work,4*m+1,i)
      if(i/=0) stop 'gmresm: dgelsy'
      y = p

      p(1:j+1) = - matmul(h(1:j+1,1:j),y(1:j))
      p(1) = p(1) + beta
      res = dsqrt(dot_product(p(1:j+1),p(1:j+1)))
      if(info==1) print*, 'gmresm: it=', its,' res=', real(res)
      
      done = (res<=tol .or. its==imx .or. res>res_)
      if(done .or. j==m) then
        if(del>0d0)  call hookstep(j,h,m,beta,del, y)
         z = matmul(v(:,1:j),y(1:j))
         call psolve(n, z)
         x = x + z
        if(its==imx) info = 2
        if(res>res_) info = 1
        if(res<=tol) info = 0
         if(done)     return
        if(del>0d0)  print*, 'gmres: warning! restart affects hookstep'
         goto 1       ! (j==m) restart
      end if
      res_ = res*stgn

   end do   
 
 end subroutine gmresm
 
 
!-----------------------------------------------------------------
! replace y with a vector that generates a hookstep
! c.f. Viswanath (2008) arXiv:0809.1498
!-----------------------------------------------------------------
 subroutine hookstep(j,h,m,beta,del, y)
   implicit none
   integer,          intent(in)    :: j, m
   double precision, intent(in)    :: h(m+1,j), beta
   double precision, intent(inout) :: del
   double precision, intent(out)   :: y(j)
   double precision :: a(j+1,j), s(j), u(j+1,j+1), vt(j,j), work(5*(j+1))
   double precision :: p(j+1), q(j), mu, qn
   integer :: info
   
   a = h(1:j+1,1:j)
   
   call dgesvd('A','A',j+1,j,a,j+1,s,u,j+1,vt,j,work,5*(j+1),info)
   if(info/=0) stop 'hookstep: dgesvd'
   
   p(1:j) = beta * u(1,1:j)   

   mu = max(s(j)*s(j)*1d-6,1d-99)
   qn = 1d99
   do while(qn>del)
      mu = mu * 1.1d0
      q = p(1:j)*s/(mu+s*s)
      qn = dsqrt(dot_product(q,q))
   end do

   y = matmul(q,vt)

   p = - matmul(h(1:j+1,1:j),y(1:j))
   p(1) = p(1) + beta
   del = dsqrt(dot_product(p,p))
 
 end subroutine hookstep

!-------------------------------------------------------------------------------------------


end module m_precon
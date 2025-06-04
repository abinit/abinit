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

    use defs_abitypes,      only : MPI_type
    use defs_basis
    use m_dtset
    use m_dtfil
    use m_xmpi

    use defs_datatypes,     only : pseudopotential_type
    use defs_wvltypes
    use m_atomdata,         only : atom_length
    use m_dfpt_mkvxc,       only : dfpt_mkvxc, dfpt_mkvxc_noncoll
    use m_fft,              only : fourdp, fourwf, fftpac
    use m_fftcore,          only : sphereboundary
    use m_mkrho
    use m_mpinfo,           only : proc_distrb_cycle
    use m_paw_dmft
    use m_pawrhoij,         only : pawrhoij_type
    use m_pawcprj,          only : pawcprj_type
    use m_pawang,           only : pawang_type
    use m_pawfgr,           only : pawfgr_type
    use m_pawtab,           only : pawtab_type
    use m_pawfgrtab,        only : pawfgrtab_type
    use m_paw_occupancies,  only : pawmkrhoij
    use m_paw_mkrho,        only : pawmkrho
    use m_spacepar,         only : symrhg

#if defined HAVE_LINALG_MKL_OMATCOPY
    use mkl_rci, only : dfgmres, dfgmres_check, dfgmres_get, dfgmres_init
#endif
    
    implicit none
    private
    public :: precon_object
    public :: linsolve

    type, public :: precon_object
        integer  :: iprcel
        !integer  :: nspden
        real(dp) :: dielng, diemix
        !Geometry :
        real(dp) :: gprimd(3, 3), rprimd(3, 3)
        real(dp) :: ucvol, dvol
        !PAW :
        type(pseudopotential_type), pointer :: psps
        integer :: mband_cprj
        integer, pointer :: dimcprj(:)
        integer :: unpaw
        type(pawcprj_type), pointer :: cprj(:, :)
        type(pawang_type), pointer :: pawang
        type(pawfgr_type), pointer :: pawfgr
        type(pawfgrtab_type), pointer :: pawfgrtab(:)
        type(pawtab_type), pointer :: pawtab(:)
        !To compute (weighted) densities and other quantities :
        real(dp), pointer :: fermie
        real(dp), pointer :: cg(:, :), eigen(:), phnons(:, :, :)
        integer, pointer  :: kg(:, :), npwarr(:), irrzon(:, :, :)
        integer, pointer :: atindx(:), atindx1(:), nattyp(:)
        integer, pointer :: symrec(:, :, :), indsym(:, :, :)
        real(dp), pointer :: xred(:, :)
        !For ffts :
        integer :: nfftprc  ! Number of fft grid points for preconditioned densities and/or potentials.
        !For Kxc :
        integer :: nkxc
        logical :: need_kxc
        real(dp), pointer :: kxc(:, :)
        real(dp), pointer :: rhor(:, :)
        real(dp), pointer :: vxc(:, :)

        !For LDOS preconditioner :
        real(dp) :: tdos
        real(dp), allocatable :: ldos(:, :)
        !For local polarizability preconditioner :
        real(dp) :: gc
        real(dp), allocatable :: loc_pola(:, :)

    contains
        procedure :: init => precon_init            ! Initialize the precon_object.
        procedure :: init_kxc => precon_init_kxc    ! Initialize kxc in the precon_object.
        procedure :: update => precon_update        ! Update the precon_object according to iprcel.
        procedure :: free => precon_free            ! Dealocate arrays that are allocated in precon_init.
        procedure :: save => precon_save            ! Save the LDOS or local polarizability contained in the precon_object in a file.
        procedure :: to_pauli => to_pauli           ! Basis change from the default abinit spin basis to the Pauli basis 
                                                    ! for density and potentials.
        procedure :: from_pauli => from_pauli       ! Basis change from the Pauli basis to the default abinit spin basis 
                                                    ! for density and potentials.
        procedure :: apply_kernel => apply_kernel   ! Apply the Coulomb kernel (and exchange and correlation kernel depending on iprcel) 
                                                    ! to an input vector.
        procedure :: apply_chi0 => apply_chi0       ! Apply the model chi0 operator to an input vector.
        procedure :: save_applied_op_g => save_applied_op_g ! Save the application of an operator in reciprocal space (for code validation).
        procedure :: save_applied_op_r => save_applied_op_r ! Save the application of an operator in direct space (for code validation).

    end type precon_object

contains 

    !****f* m_precon/precon_init
    !! NAME
    !!  precon_init
    !!
    !! FUNCTION
    !!  Initialize the precon_object.
    !!
    !! INPUTS
    !!  dtset   = all input variables for this dataset
    !!  atindx  =
    !!  atindx1 = index table for atoms, inverse of atindx (see gstate.f)
    !!  cg      = wf in G space
    !!  cprj    =
    !!  dimcprj =
    !!  eigen   = array of eigenvalues
    !!  fermie  = fermi energie
    !!  gprimd  = dimensional reciprocal space primitive translations
    !!  irrzon  = irreducible zone data
    !!  kg      = reduced planewave coordinates
    !!  nattyp  = number of atoms of each type.
    !!  nfftmix = 
    !!  npwarr  = number of planewaves and boundary planewaves at each k
    !!  pawang  =
    !!  pawfgr  =
    !!  pawfgrtab =
    !!  pawtab  =
    !!  phnons  = nonsymmorphic translation phases
    !!  psps    =
    !!  rprimd  = dimensional real space primitive translations
    !!  ucvol   = unit cell volume
    !!  xred    = reduced dimensionless atomic coordinates
    !!
    !! SOURCE
    subroutine precon_init(this, dtset, atindx, atindx1, cg, cprj, dimcprj, dtfil, eigen, fermie, &
        &   gprimd, indsym, irrzon, kg, mband_cprj, nattyp, nfftmix, npwarr, pawang, pawfgr, pawfgrtab, &
        &   pawtab, phnons, psps, rhor, rprimd, symrec, ucvol, vxc, xred)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        real(dp), intent(in) :: ucvol
        real(dp), intent(in), target :: fermie
        integer, intent(in) :: nfftmix, mband_cprj
        !arrays
        real(dp), intent(in) :: gprimd(:, :), rprimd(:, :)
        integer, intent(in), target  :: irrzon(:, :, :), kg(:, :), npwarr(:)
        integer, intent(in), target :: atindx(:), atindx1(:), nattyp(:)
        integer, intent(in), target :: symrec(:, :, :), indsym(:, :, :)
        real(dp), intent(in), target :: cg(:, :), eigen(:), phnons(:, :, :)
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

        ! *************************************************************************
        !Constant data from dtset
        this%dielng = dtset%dielng
        this%diemix = dtset%diemix
        this%iprcel = dtset%iprcel
        this%nfftprc = nfftmix 
        !dtset%nspden = dtset%nspden
        !Other constants
        this%dvol   = ucvol/dtset%nfft ! factor for integrals in real space: sum(f) * dvol ~ integral f
        this%gprimd = gprimd
        this%rprimd = rprimd
        this%ucvol  = ucvol
        this%need_kxc = .false.
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
        this%phnons => phnons
        this%psps   => psps
        this%rhor   => rhor
        this%symrec => symrec
        this%vxc    => vxc
        this%xred   => xred
        
        !PAW :
        if (psps%usepaw==1) then
            this%unpaw      = dtfil%unpaw
            this%mband_cprj = mband_cprj
            this%dimcprj    => dimcprj
            this%cprj       => cprj
            this%pawang     => pawang
            this%pawfgr     => pawfgr
            this%pawfgrtab  => pawfgrtab
            this%pawtab     => pawtab
        end if
        
        !Initializing LDOS specific variables
        if (this%iprcel == 202) then
            !Allocating the array containing ldos
            ABI_MALLOC(this%ldos, (this%nfftprc, dtset%nspden)) ! TODO
        end if
        
        !Initializing chi0_diag specific variables
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

        if (this%need_kxc) then
            this%kxc => kxc
        end if

    end subroutine precon_init_kxc

    !****f* m_precon/precon_update
    !! NAME
    !!  precon_update
    !!
    !! FUNCTION
    !!  Update the precon_object :
    !!      For LDOS preconditioning (iprcel=202) : 
    !!          Compute the new ldos (local density of state) with current wavefunctions 
    !!          and the new tdos (total density of state = integral of ldos).
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
            call compute_ldos(this, dtset, mpi_enreg, this%ldos)
            !update tdos
            this%tdos = sum(this%ldos(:, 1)) * this%dvol
            ! TODO : More options to control when the ldos is updated
        end if

        !Local polarizability
        if (this%iprcel == 204) then
        end if
       
    end subroutine precon_update

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
       
        if (this%iprcel == 202) then
            !Deallocating the array containing ldos and tdos
            ABI_FREE(this%ldos)
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
        i1 = modulo((ifft-1), n1) + 1
        i2 = modulo((ifft-1)/n1, n2) + 1
        i3 = ((ifft-1)/n1)/n2 + 1
        g(1) = modulo(i1-1 + n1/2, n1) - n1/2
        g(2) = modulo(i2-1 + n2/2, n2) - n2/2
        g(3) = modulo(i3-1 + n3/2, n3) - n3/2

    end function get_g_vector

    !****f* m_precon/precon_save
    !! NAME
    !!  precon_save
    !!
    !! FUNCTION
    !!  Saves the LDOS contained in the precon_object in a file named ldos.txt.
    !!      "     locale polarizability                "              loc_pola.txt. 
    !!  (For code validation)
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
    !!  Basis change from the Abinit default spin-basis to the Pauli basis for potentials and densities
    !!  in the reciprocal (G) space.
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
    !!  in the direct (r) space.
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
    !!  Basis change from the Pauli basis to the Abinit default spin-basis for potentials and densities
    !!  in the reciprocal (G) space.
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
    !!  apply_vc
    !!
    !! FUNCTION
    !!  Apply the Coulomb kernel vc to a vector (in place).
    !!
    !! INPUTS
    !!  ngfftprc   = All needed informations about the 3D FFT for preconditioned quantities.
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the Coulomb kernel vc is applied (in place).
    !!                               When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_vc(this, dtset, ngfftprc, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type),intent(in) :: dtset
        !arrays
        integer, intent(in) :: ngfftprc(:)
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        integer :: ifft, ispden
        real(dp) :: g_cart_2
        
        ! *************************************************************************
        
        !In the sigma_0, 1, 2, 3 basis :
        !   The sigma_0 component of the density is multiplied by 4pi/G^2
        !   and the rest is 0.
        ispden = 1
        do ifft = 2, this%nfftprc
            g_cart_2 = norm2(two_pi * matmul(this%gprimd, get_g_vector(ifft, ngfftprc)))**2
            vec_g(:, ifft, ispden) = (2*two_pi/g_cart_2) * vec_g(:, ifft, ispden)
        end do

        do ispden = 2, dtset%nspden
            vec_g(:, :, ispden) = 0
        end do

    end subroutine apply_vc

    !****f* m_precon/apply_kxc
    !! NAME
    !!  apply_kxc
    !!
    !! FUNCTION
    !!  Apply the exchange and correlation kernel Kxc to a vector (in place).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Informations about MPI parallelization.
    !!  ngfftprc    = All needed informations about the 3D FFT for preconditioned quantities.
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the exchange and correlation kernel Kxc is applied (in place).
    !!                               When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_kxc(this, dtset, mpi_enreg, ngfftprc, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfftprc(:)
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        !scalars
        integer :: cplex, n3xccc, nhatdim, nhat1dim, nhat1grdim, nkxc, option, optnc, usexcnhat
        logical :: non_magnetic_xc
        !arrays
        real(dp), allocatable :: Kxc_vec_r(:, :), vec_r(:, :)
        !real(dp), allocatable :: vec_g_old(:, :, :)            !DEBUG
        real(dp), allocatable :: nhat(:, :), nhat1(:, :), nhat1gr(:, :, :)
        real(dp) :: dummy_xccc3d1(0), dummy_qphon(3)
        
        ! *************************************************************************

        if (size(this%kxc, 1) /= this%nfftprc) then
            ABI_BUG("iprcel=203 : size(kxc, 1) /= nfftprc")
        end if

        !Basis change to get vec_g in the default Abinit spin-basis to apply existing routines
        call from_pauli(this, 1, vec_g)
        !ABI_MALLOC(vec_g_old, (2, this%nfft, dtset%nspden))    !DEBUG
        !vec_g_old = vec_g                                      !DEBUG

        !ifft to get vec_g in the real-space basis required by dfpt_mkvxc
        ABI_MALLOC(vec_r, (this%nfftprc, dtset%nspden))
        call fourdp(1, vec_g, vec_r, 1, mpi_enreg, this%nfftprc, dtset%nspden, ngfftprc, 0)

        !Applying Kxc : 
        cplex = 1   ! Input vector is real in real (direct) space.
        non_magnetic_xc = .false.
        nkxc = size(this%kxc, 2)

        usexcnhat = 0                                               !
        nhat1dim = 0                                                ! 
        ABI_MALLOC(nhat1, (cplex*this%nfftprc, dtset%nspden*nhat1dim)) ! PAW - TODO ?
        nhat1grdim = 0                                              !
        ABI_MALLOC(nhat1gr, (cplex*this%nfftprc, dtset%nspden, 3*nhat1grdim))  !

        option = 2  ! Treats only density change (no core_correction)
        n3xccc = 0  !   -> Core-correction set to 0.

        ABI_MALLOC(Kxc_vec_r, (this%nfftprc, dtset%nspden))

        if (dtset%nspden==1 .or. dtset%nspden==2) then
            call dfpt_mkvxc(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfftprc, ngfftprc, nhat1, nhat1dim, &
            &               nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden, n3xccc, option, &
            &               dummy_qphon, vec_r, this%rprimd, usexcnhat, Kxc_vec_r, dummy_xccc3d1)
        end if
        
        if (dtset%nspden==4) then
            nhatdim = 0
            ABI_MALLOC(nhat, (this%nfftprc, dtset%nspden*nhatdim)) ! TODO : PAW
            optnc = 1   ! Compute the whole 2x2 Vres matrix
            call dfpt_mkvxc_noncoll(cplex, dtset%ixc ,this%kxc, mpi_enreg, this%nfftprc, ngfftprc, nhat, nhatdim, &
            &               nhat1, nhat1dim, nhat1gr, nhat1grdim, nkxc, non_magnetic_xc, dtset%nspden,      &
            &               n3xccc, optnc, option, dummy_qphon, this%rhor, vec_r, this%rprimd, usexcnhat,        &
            &               this%vxc, Kxc_vec_r, dummy_xccc3d1)
           ABI_FREE(nhat)
        end if

        ABI_FREE(nhat1)
        ABI_FREE(nhat1gr)

        !fft to return the result Kxc_vec_r in the plane-wave basis.
        call fourdp(1, vec_g, Kxc_vec_r, -1, mpi_enreg, this%nfftprc, dtset%nspden, ngfftprc, 0)
        !call this%save_applied_op_r(dtset, ngfftprc, vec_r, Kxc_vec_r, "applied_kxc_r.txt")  !DEBUG

        ABI_FREE(Kxc_vec_r)
        ABI_FREE(vec_r)

        !call this%save_applied_op_g(dtset, ngfftprc, vec_g_old, vec_g, "applied_kxc.txt")  !DEBUG
        !ABI_FREE(vec_g_old)                                    !DEBUG

        !Basis change to get the resulting vector in the Pauli basis used in the rest of the preconditioning routines.
        call to_pauli(this, 0, vec_g)

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
    !!  ngfftprc    = All needed informations about the 3D FFT for preconditioned quantities.
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the kernel is applied (in place).
    !!                               When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_kernel(this, dtset, mpi_enreg, ngfftprc, vec_g)
        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type),intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfftprc(:)
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
        
        !Local variables-------------------------------
        real(dp), allocatable :: Kxc_vec_g(:, :, :)

        ! *************************************************************************
        write(6,*)'chi0diel apply_kernel'; flush(6) !DEBUG
        
        ! RPA : LDOS/Kerker model - only vc
        if (this%iprcel == 201 .or. this%iprcel == 202) then
            call apply_vc(this, dtset, ngfftprc, vec_g)
        end if

        ! No RPA : vc and Kxc
        if (this%iprcel == 203) then
            ABI_MALLOC(Kxc_vec_g, (2, this%nfftprc, dtset%nspden))
            Kxc_vec_g = vec_g   !Copy ?
            call apply_vc(this, dtset, ngfftprc, vec_g)    ! Apply vc in place
            call apply_Kxc(this, dtset, mpi_enreg, ngfftprc, Kxc_vec_g)
            vec_g = vec_g + Kxc_vec_g
            ABI_FREE(Kxc_vec_g)
        end if

    end subroutine apply_kernel

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
    !!  w_rhor      = "Weighted density" in real space.
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
        integer :: dummy_int
        integer :: mcg, cplex, nfftot, mcprj
        real(dp) :: compch_fft
        !arrays
        real(dp) :: qphon(3)
        real(dp), allocatable :: w_rhowfg(:, :), w_rhowfr(:, :)
        type(paw_dmft_type)     :: dummy_paw_dmft
        type(wvl_wf_type)       :: dummy_wvl_wfs
        type(wvl_denspot_type)  :: dummy_wvl_den

        ! *************************************************************************

        if (this%psps%usepaw==0) then
            ABI_MALLOC(w_rhowfr, (dtset%nfft, dtset%nspden))
            ABI_MALLOC(w_rhowfg, (2, dtset%nfft))
        else
            ABI_MALLOC(w_rhowfr, (this%pawfgr%nfftc, dtset%nspden))
            ABI_MALLOC(w_rhowfg, (2, this%pawfgr%nfftc))
        end if

        ! Compute the weighted density (w_rhor) using mkrho with weights in place of the occupations.
        mcg = size(this%cg)
        dummy_paw_dmft%use_dmft = 0
        dummy_paw_dmft%use_sc_dmft = 0
        call mkrho(this%cg, dtset, this%gprimd, this%irrzon, this%kg, mcg, mpi_enreg, this%npwarr, weights, &
        &   dummy_paw_dmft, this%phnons, w_rhowfg, w_rhowfr, this%rprimd, 0, this%ucvol, dummy_wvl_den, dummy_wvl_wfs, option=0)
        
        ! Symmetrize the weighted density (w_rhor)
        nfftot = dtset%ngfft(1) * dtset%ngfft(2) * dtset%ngfft(3)
        call symrhg(1, this%gprimd, this%irrzon, mpi_enreg, dtset%nfft, nfftot, dtset%ngfft, dtset%nspden, dtset%nsppol, &
        &   dtset%nsym, this%phnons, w_rhowfg, w_rhowfr, this%rprimd, dtset%symafm, dtset%symrel, dtset%tnons)

        if (this%psps%usepaw==0) then
        ! In NC : the weighted density is directly w_rhowfr.
            if (this%nfftprc == dtset%nfft) then
                w_rhor = w_rhowfr
            else
                ABI_BUG("iprcel=2** : nfftprc /= nfft not implemented. TODO")
            end if
        else
        ! In PAW : Add rhoij terms to  w_rhowfr.
            if (this%nfftprc == this%pawfgr%nfft) then
                !Compute the rhoij equivalent for the weighted density.
                !   Sum_{n,k} {weight(n,k)*<Cnk|p_i><p_j|Cnk>}.
                mcprj = size(this%cprj)
                call pawmkrhoij(this%atindx, this%atindx1, this%cprj, this%dimcprj, dtset%istwfk, dtset%kptopt, dtset%mband, this%mband_cprj, &
                &       mcprj, dtset%mkmem, mpi_enreg, dtset%natom, dtset%nband, dtset%nkpt, dtset%nspden, dtset%nspinor, &
                &       dtset%nsppol, weights, dtset%paral_kgb, dummy_paw_dmft, pawrhoij, this%unpaw, dtset%usewvl, dtset%wtk)

                !Compute the total weighted density.
                cplex = 1
                dummy_int=0
                qphon = 0
                call pawmkrho(1, compch_fft, cplex, this%gprimd, dummy_int, this%indsym, dummy_int, mpi_enreg, &
                &       mpi_enreg%my_natom, dtset%natom, dtset%nspden, dtset%nsym, dtset%ntypat, dtset%paral_kgb, this%pawang, this%pawfgr, this%pawfgrtab, &
                &       dtset%pawprtvol, pawrhoij, pawrhoij, this%pawtab, qphon, w_rhowfg, w_rhowfr, w_rhor, this%rprimd, dtset%symafm, &
                &       this%symrec, dtset%typat, this%ucvol, dtset%usewvl, this%xred)
            else
                ABI_BUG("iprcel=2** : nfftprc /= pawfgr%nfft not implemented. TODO")
            end if
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
        integer :: mband, i, maxocc
        !arrays
        real(dp), allocatable :: ldos_weights(:)

        ! *************************************************************************

        !compute weights
        mband = size(this%eigen)
        ABI_MALLOC(ldos_weights, (mband))
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        do i=1, mband
            ldos_weights(i) = -derivative_occ(dtset%occopt, this%eigen(i), this%fermie, dtset%tsmear) * maxocc
        end do

        !Compute ldos using mkrho with ldos_weights in place of the occupations
        call compute_weighted_density(this, dtset, mpi_enreg, ldos_weights, ldos)
        ABI_FREE(ldos_weights)

        !With collinear spins the ldos is not returned in the Pauli (tot/spin) basis by mkrho.
        if (dtset%nspden == 2) then
            !spin = 2*up - tot
            ldos(:, 2) = 2*ldos(:, 2) - ldos(:, 1)
        end if

    end subroutine compute_ldos

    !****f* m_precon/complex_mult
    !! NAME
    !!  complex_mult
    !!
    !! FUNCTION
    !!  Multiply two complex numbers given as size two arrays.
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

    !****f* m_precon/save_applied_op_g
    !! NAME
    !!  save_applied_op_g
    !!
    !! FUNCTION
    !!  Save the applied operator vec/op_vec in the reciprocal (G) space in a file.
    !!  For code validation only.
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
    !!  save_applied_op_r
    !!
    !! FUNCTION
    !!  Save the applied operator vec/op_vec in the direct (real) space in a file.
    !!  For code validation only.
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
    !!  apply_chi0_ldos
    !!
    !! FUNCTION
    !!  Apply the ldos model chi0 operator to the vector vec_g (in place).
    !!
    !! INPUTS
    !!  mpi_enreg    = Information about MPI parallelization.
    !!  ngfftprc     = Contain all needed information about 3D FFT, see ~abinit/doc/variables/gstate/#ngfft.
    !!  ispden       = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                             When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0_ldos(this, dtset, mpi_enreg, ngfftprc, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        type(dataset_type),intent(in) :: dtset
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        integer, intent(in) :: ngfftprc(:)
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
       
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
            ABI_MALLOC(vec_r_1, (cplex*this%nfftprc))    
            ABI_MALLOC(work_r, (cplex*this%nfftprc))
            
            !1) ifft to get (the sigma_0 component of) vec in the real space
            ispden = 1
            call fourdp(cplex, vec_g(:, :, ispden), vec_r_1, 1, mpi_enreg, this%nfftprc, 1, ngfftprc, 0)

            do ispden = 1, dtset%nspden
                !2) chi0(v)(r)_ispden = -ldos_ispden(r)*v_1(r) + 1/dos * ldos_ispden(r)*integral(ldos_1(r')*v_1(r')*dr')
                work_r = -this%ldos(:, ispden)*vec_r_1 &
                                + 1/this%tdos * dot_product(this%ldos(:, 1), vec_r_1)*this%dvol * this%ldos(:, ispden)

                !3) fft to get vec back in the reciprocal space
                call fourdp(cplex, vec_g(:, :, ispden), work_r, -1, mpi_enreg, this%nfftprc, 1, ngfftprc, 0)
            end do

            ABI_FREE(work_r)
            ABI_FREE(vec_r_1)
        else 
            vec_g = 0
        end if

    end subroutine apply_chi0_ldos

    !****f* m_precon/compute_weights_chi0_diag
    !! NAME
    !!  compute_weights_chi0_diag
    !!
    !! FUNCTION
    !!  Compute the weights needed to compute the application of the diagonal model chi0 to a vector (vec_g)
    !!  such that chi0 * vec = sum_i weight_i * rho_ii
    !!  that is 
    !!      weight_i = f'i * dot(vec, rho_ii).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the model chi0 operator is to be applied.
    !!                               When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! OuTPUTS
    !!  weights(:) = Weights needed for the application of the model chi0 to vec_g. To be used in apply_chi0_diag.
    !!
    !! SOURCE
    subroutine compute_weights_chi0_diag(this, dtset, mpi_enreg, vec_g, weights)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        !arrays
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
        real(dp), intent(inout) :: weights(:)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden, istwf_k, maxocc, mcg, my_nspinor, mband_mem, nband_k, ndat, nfftot, npw_k, nspin, option, tim_fourwf
        integer :: i_eigen, i_cg, j_cg, ikg, ikpt, iband, isppol, ier
        integer :: n1, n2, n3, n4, n5, n6
        real(dp) :: fp, eigenval
        integer :: dummy_int
        real(dp) :: dummy_real
        real(dp) :: weight_r, weight_i
        !arrays
        integer :: gbound(2*dtset%mgfft+8,2)
        integer, allocatable :: kg_k(:, :)
        real(dp), allocatable :: vec_r(:, :)
        real(dp), allocatable :: rhoaug_r_ii(:, :, :, :), rho_r_ii(:, :)
        real(dp), allocatable :: psi_r_i_up(:, :, :, :), psi_r_i_down(:, :, :, :)
        real(dp) :: vec_g_old(2, this%nfftprc, dtset%nspden)   !DEBUG

        !dummy fourwfk arguments
        real(dp), allocatable ::  dummy_denpot(:, :, :), dummy_fofgout(:, :), dummy_fofrout(:, :, :, :)
        
        ! *************************************************************************
        write(6,*)'chi0diel apply_chi0_mag'; flush(6) !DEBUG
        vec_g_old = vec_g   !DEBUG
       
        n1 = dtset%ngfft(1)
        n2 = dtset%ngfft(2)
        n3 = dtset%ngfft(3)
        n4 = dtset%ngfft(4)
        n5 = dtset%ngfft(5)
        n6 = dtset%ngfft(6)
        if (this%nfftprc /= n1*n2*n3) then      ! TODO :  Incorrect for PAW
            ABI_BUG("Mismatch nfft /= n1*n2*n3")
        end if

        !1) Computing vec_r = input vector in real space
        !   - in the up/down basis if nspinor = 1
        !   - in the Pauli basis if nspinor = 2
        if (dtset%nspinor == 1) then
            call from_pauli(this, 0, vec_g)
        end if
        ABI_MALLOC(vec_r, (this%nfftprc, dtset%nspden))
        write(6,*)'chi0diel apply_chi0_mag before fourdp vec_g(1, 1:10, 2): ', vec_g(1, 1:10, 2); flush(6) !DEBUG
        call fourdp(1, vec_g, vec_r, 1, mpi_enreg, this%nfftprc, 2, dtset%ngfft, 0)
        
        !2) Computing the weights = fi' * <rhoii, vec_spin> TODO
        !   vec_spin is the spin component (in basis {tot, spin}) of the input vector vec_g
        weights = 0
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        write(6,*)'chi0diel apply_chi0_mag maxocc, this%dvol: ', maxocc, this%dvol; flush(6) !DEBUG
        
        !Allocating the arrays that will contain rho_ii
        if (dtset%nspinor==1) then
            nspin = 1   ! Number of spin components in the orbital densities.
        else if (dtset%nspinor==2) then
            nspin = 4
            ABI_MALLOC(psi_r_i_up, (2, n4, n5, n6))
            ABI_MALLOC(psi_r_i_down, (2, n4, n5, n6))
        else
            ABI_BUG("nspinor /= 1 or 2")
        end if
        ABI_MALLOC(rhoaug_r_ii, (n4, n5, n6, nspin))
        ABI_MALLOC(rho_r_ii, (n1*n2*n3, nspin))

        !Allocating dummy arrays
        ABI_MALLOC(dummy_denpot, (0, n5, n6))
        ABI_MALLOC(dummy_fofgout, (2, 0))
        ABI_MALLOC(dummy_fofrout, (2, n4, n5, n6))

        i_cg = 0    ! Starting index for (ikpt, isppol) in cg array.
        my_nspinor = max(1, dtset%nspinor/mpi_enreg%nproc_spinor)
        if (my_nspinor /= dtset%nspinor) then
            ABI_BUG("SCF preconditioner 'chi0_diag' incompatible with spinor parallelisation")  ! TODO ?
        end if

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

                call sphereboundary(gbound, istwf_k, kg_k, dtset%mgfft, npw_k)    ! Computes gbound.
                
                do iband = 1, nband_k

                    !Indices
                    !TODO : verifier mband = nbandk ? + verifier si eigen est distribué pour paral
                    i_eigen = iband + (ikpt-1)*dtset%mband + (isppol-1)*dtset%mband*dtset%nkpt  ! Index of (iband, ikpt, isppol) in eigen array.
                    j_cg = i_cg + (iband-1) * npw_k * my_nspinor    ! Index of (iband, ikpt, isppol) in cg array.                    
                    
                    !2.1) Computing f'(eig_i - fermie).
                    eigenval = this%eigen(i_eigen)
                    fp = derivative_occ(dtset%occopt, eigenval, this%fermie, dtset%tsmear) * maxocc

                    if (abs(fp) > tol14) then   !TODO : choose tol ?
                        
                        !2.2) Computing rho_ii = |psi_i|^2 using fourwf (if fp is not 0).
                        
                        ! No spin or collinear spins - Wafefunctions have one spin component.
                        if (dtset%nspinor == 1) then
                            
                            !a) Compute orbital density in real space (augmented basis) :
                            !Input parameters for fourwf :
                            option = 1          ! Computes the density.
                            ndat = 1            ! Only one FFT.
                            tim_fourwf = 0
                            rhoaug_r_ii = zero  ! Initialization for fourwf (accumulation).
                            weight_r = 1
                            weight_i = 1

                            call fourwf(1, rhoaug_r_ii, this%cg(:, j_cg+1:j_cg+npw_k), dummy_fofgout, dummy_fofrout,  &
                            &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            &           dummy_int, n4, n5, n6, option, tim_fourwf, weight_r, weight_i)

                            !b) Compute the dot product between the input vector and the orbital density and update the weights :
                            ! Transfer rhoaug_r_ii defined on the large fft-grid to the smaller density/potential fft-grid.
                            call fftpac(1, mpi_enreg, 1, n1, n2, n3, n4, n5, n6, dtset%ngfft, rho_r_ii, rhoaug_r_ii(:, :, :, 1), 1)   ! DEBUG : weights for nspinor=1, nsppol=2 ok (factor from fft?)
                            rho_r_ii(:, 1) = rho_r_ii(:, 1) / (sum(rho_r_ii(:, 1)) * this%dvol) !Normalizing rho_ii_r.   [ sqrnorm_v, meanfft_r, dotprod_vn in cg_tools/ cgtools? ]
                            ! dot product:
                            weights(i_eigen) = fp * maxocc * dot_product(rho_r_ii(:, 1), vec_r(:, isppol)) * this%dvol  ! dot-product in the up/down basis
                        end if

                        ! Non collinear spins - Wavefunctions have two spins components.
                        if (dtset%nspinor == 2) Then

                            !a) Compute psi_up and psi_down in real space :
                            ! Input parameters for fourwf :
                            option = 0          ! Only do the FFT.
                            ndat = 1
                            tim_fourwf = 0
                            !FFT for psi_up
                            call fourwf(dummy_int, dummy_denpot, this%cg(:, j_cg+1:j_cg+npw_k), dummy_fofgout, psi_r_i_up,  &
                            &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            &           dummy_int, n4, n5, n6, option, tim_fourwf, dummy_real, dummy_real)
                            !FFT for psi_down
                            call fourwf(dummy_int, dummy_denpot, this%cg(:, j_cg+npw_k+1:j_cg+2*npw_k), dummy_fofgout, psi_r_i_down,  &
                            &           gbound, gbound, istwf_k, kg_k, kg_k, dtset%mgfft, mpi_enreg, ndat, dtset%ngfft, npw_k, &
                            &           dummy_int, n4, n5, n6, option, tim_fourwf, dummy_real, dummy_real)
                            !   (done separately for convenience & readability)
                            
                            !b) Compute the 4 components of the orbital density rhoaug_r_ii :
                            !   (in the Pauli basis for spins and real augmented basis for space)
                            ispden = 1  ! rho_sigma0 = |psi_up|^2 + |psi_down|^2
                            rhoaug_r_ii(:, :, :, ispden) = psi_r_i_up(1, :, :, :)**2 + psi_r_i_up(2, :, :, :)**2 + psi_r_i_down(1, :, :, :)**2 + psi_r_i_down(2, :, :, :)**2
                            ispden = 2  ! rho_sigma1 = psi_up* psi_down + psi_down* psi_up
                            rhoaug_r_ii(:, :, :, ispden) = 2*( psi_r_i_up(1, :, :, :)*psi_r_i_down(1, :, :, :) + psi_r_i_up(2, :, :, :)*psi_r_i_down(2, :, :, :) )
                            ispden = 3  ! rho_sigma2 = i*(psi_down* psi_up - psi_up* psi_down)
                            rhoaug_r_ii(:, :, :, ispden) = 2*( psi_r_i_up(2, :, :, :)*psi_r_i_down(1, :, :, :) - psi_r_i_up(1, :, :, :)*psi_r_i_down(2, :, :, :) )
                            ispden = 4  ! rho_sigma3 = |psi_up|^2 - |psi_down|^2
                            rhoaug_r_ii(:, :, :, ispden) = psi_r_i_up(1, :, :, :)**2 + psi_r_i_up(2, :, :, :)**2 - psi_r_i_down(1, :, :, :)**2 + psi_r_i_down(2, :, :, :)**2
                       
                            !c) Compute the dot product between the input vector and the orbital density and update the weights :
                            do ispden=1, 4
                                ! Transfer rhoaug_r_ii defined on the large fft-grid to the smaller density/potential fft-grid.
                                call fftpac(ispden, mpi_enreg, 4, n1, n2, n3, n4, n5, n6, dtset%ngfft, rho_r_ii(:, ispden), rhoaug_r_ii(:, :, :, ispden), 1)
                                rho_r_ii(:, ispden) = rho_r_ii(:, ispden) / (sum(rho_r_ii(:, ispden)) * this%dvol) !Normalizing rho_ii_r. TODO : check that normalizing each component make sense.
                                ! dot product in Pauli basis :
                                weights(i_eigen) = weights(i_eigen) + fp * maxocc * dot_product(rho_r_ii(:, ispden), vec_r(:, ispden)) * this%dvol  ! dotprod_vn
                            end do
                        end if

                        ! TODO use accumulation option of fourdp to make the computation more efficient ? Impossible IMO

                        !if (ikpt==1) then
                        !write(6,*)'--- chi0diel apply_chi0_mag:  i_eigen, iband, ikpt, isppol', i_eigen, iband, ikpt, isppol; flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: eigenval', eigenval; flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: fp', fp; flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: rhoaug_r_ii( 1:5, 5, 5)', rhoaug_r_ii(1:5, 5, 5); flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: rho_r_ii(1:5, 1)', rho_r_ii(1:5, 1); flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: sum(rho_r_ii(:, 1))', sum(rho_r_ii(:, 1)); flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: dot_product(rho_r_ii(:, 1), vec_r(:, 2))*this%dvol', dot_product(rho_r_ii(:, 1), vec_r(:, 2))*this%dvol; flush(6) !DEBUG
                        !write(6,*)'chi0diel apply_chi0_mag: weights(i_eigen)', weights(i_eigen); flush(6) !DEBUG
                        !end if

                    end if
                    
                end do

                if (dtset%mkmem /= 0) then
                    i_cg = i_cg + npw_k * dtset%nspinor * mband_mem
                    !i_cg = i_cg + npw_k * my_nspinor * mband_mem
                    ikg = ikg + npw_k
                end if
                ABI_FREE(kg_k)

            end do  !ikpt
        end do  !isppol

        ABI_FREE(rhoaug_r_ii)
        if (dtset%nspinor==2) then
            ABI_FREE(psi_r_i_up)
            ABI_FREE(psi_r_i_down)
        end if
        ABI_FREE(rho_r_ii)
        ABI_FREE(vec_r)
        ABI_FREE(dummy_denpot)
        ABI_FREE(dummy_fofgout)
        ABI_FREE(dummy_fofrout)

        !MPI parallelization over kpoints : sum weights on all processors.
        ier = 0
        call xmpi_sum(weights, mpi_enreg%comm_kpt, ier)
        
        !write(6,*)'chi0diel apply_chi0_mag weights(16:17)', weights(16:17); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_chi0_mag weights(48:49)', weights(48:49); flush(6) !DEBUG
        !write(6,*)'chi0diel apply_chi0_mag weights', weights; flush(6) !DEBUG

    end subroutine compute_weights_chi0_diag

    !****f* m_precon/apply_chi0_diag
    !! NAME
    !!  apply_chi0_mag
    !!
    !! FUNCTION
    !!  Apply the diagonal model chi0 operator to the vector vec_g.
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mgfft       = Maximum size of 1D FFTs for the "coarse" grid.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  ngfftprc    = Contain all needed information about 3D FFTfor preconditioned quantities.
    !!                See ~abinit/doc/variables/gstate/#ngfft
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                               When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    ! TODO : Name 
    subroutine apply_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfftprc, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer,intent(in) :: mgfft
        !arrays
        integer, intent(in) :: ngfftprc(:) 
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
       
        !Local variables-------------------------------
        !scalars
        integer :: cplex
        integer :: ispden, istwf_k, maxocc, mcg, my_nspinor, mband_mem, nband_k, ndat, nfftot, npw_k, option, tim_fourwf
        integer :: i_eigen, i_cg, j_cg, ikg, ikpt, iband, isppol, ier, sign_isppol
        integer :: n1, n2, n3, n4, n5, n6
        real(dp) :: fp, eigenval
        real(dp) :: dummy_real
        !arrays
        real(dp), allocatable :: chi0_vec_r(:, :)
        real(dp), allocatable :: weights(:)
        real(dp) :: vec_g_old(2, this%nfftprc, dtset%nspden)   !DEBUG

        !dummy mkrho arguments
        type(paw_dmft_type)     :: paw_dmft
        type(wvl_wf_type)       :: wvl_wfs
        type(wvl_denspot_type)  :: wvl_den
        
        ! *************************************************************************
        write(6,*)'chi0diel apply_chi0_mag'; flush(6) !DEBUG
        vec_g_old = vec_g   !DEBUG

        !1) Computing the weights : weight_i = sum_i fi' * <rhoii, vec>
        ABI_MALLOC(weights, (dtset%mband*dtset%nkpt*dtset%nsppol))
        call compute_weights_chi0_diag(this, dtset, mpi_enreg, vec_g, weights)
        
        !2) Computing chi0 * vec using mkrho with custom weights in place of the occupations.
        ABI_MALLOC(chi0_vec_r, (this%nfftprc, dtset%nspden))
        call compute_weighted_density(this, dtset, mpi_enreg, weights, chi0_vec_r) ! TODO : change compute_weighted_density to return the vector directly in G-space.
        call fourdp(1, vec_g, chi0_vec_r, -1, mpi_enreg, this%nfftprc, dtset%nspden, ngfftprc, 0)
        ABI_FREE(chi0_vec_r)
        ABI_FREE(weights)

        !Basis change : With collinear spins vec_g is not returned in the Pauli (tot/spin) basis by symrhg.
        call to_pauli(this, 1, vec_g)

        call this%save_applied_op_g(dtset, ngfftprc, vec_g_old, vec_g, "applied_chi0_mag.txt")  !DEBUG

    end subroutine apply_chi0_diag

    !****f* m_precon/apply_chi0
    !! NAME
    !!  apply_chi0
    !!
    !! FUNCTION
    !!  Apply the model (determined by iprcel) chi0 operator to the vector vec_g (in place).
    !!
    !! INPUTS
    !!  dtset       = All input variables for this dataset.
    !!  mgfft       = Maximum size of 1D FFTs for the "coarse" grid.
    !!  mpi_enreg   = Information about MPI parallelization.
    !!  ngfftprc    = Contain all needed information about 3D FFTfor preconditioned quantities.
    !!                See ~abinit/doc/variables/gstate/#ngfft 
    !!  ispden      = Index of spin-density component.
    !!
    !! SIDE EFFECTS
    !!  vec_g (2, nfftprc, nspden) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!                         When nspden > 1 vec_g is in the Pauli basis.
    !!
    !! SOURCE
    subroutine apply_chi0(this, dtset, mgfft, mpi_enreg, ngfftprc, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: mgfft
        !arrays
        integer, intent(in) :: ngfftprc(:)
        real(dp), intent(inout) :: vec_g(2, this%nfftprc, dtset%nspden)
       
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
            do ispden = 2, dtset%nspden
                vec_g(:, :, ispden) = 0
            end do
        end if
        
        !LDOS model
        if (this%iprcel == 202) then
            call apply_chi0_ldos(this, dtset, mpi_enreg, ngfftprc, vec_g)
        end if

        !Diagonal chi0 model
        if (this%iprcel == 203) then
            call apply_chi0_diag(this, dtset, mgfft, mpi_enreg, ngfftprc, vec_g)
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
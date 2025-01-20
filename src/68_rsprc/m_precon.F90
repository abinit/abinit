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
    
    use defs_abitypes, only : MPI_type
    use defs_basis
    use defs_wvltypes

    use m_atomdata,  only : atom_length
    use m_dtset
    use m_fft,      only : fourdp
    use m_mkrho
    use m_paw_dmft
    
    implicit none
    private
    public :: precon_object

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
        real(dp), allocatable :: ldos(:, :), tdos(:)
        !For local polarizability preconditioner :
        real(dp) :: gc
        real(dp), allocatable :: loc_pola(:, :)
            !To compute loc_pola :
        integer, pointer :: atindx1(:), nattyp(:)
        real(dp), pointer :: xred(:, :)
        real(dp), pointer :: rhor(:, :)
        !For ffts :
        integer  :: nfft

    contains
        procedure :: init => precon_init        ! Initializes the precon_object.
        procedure :: update => precon_update    ! Updates the precon_object according to iprcel.
        procedure :: free => precon_free        ! Dealocate arrays that are allocated in precon_init.
        procedure :: save => precon_save        ! Saves the LDOS or local polarizability contained in the precon_object in a file.
        procedure :: apply_chi0 => apply_chi0   ! Applies the model chi0 operator to an imput vector.
        procedure :: save_applied_op => save_applied_op ! Saves the application of an operator.

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
        !Pointers
        this%atindx1=> atindx1 
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
            !Allocating the arrays containing ldos and tdos
            ABI_MALLOC(this%ldos, (this%nfft, this%nspden))
            ABI_MALLOC(this%tdos, (this%nspden))
        end if
        !Initializing loc_pola specific variables
        if (this%iprcel == 203) then
            this%gc = 1.0 ! TODO
            !Allocating the array containing the local polarizability
            ABI_MALLOC(this%loc_pola, (this%nfft, this%nspden))
        end if

    end subroutine precon_init

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
            do ispden=1,this%nspden
                this%tdos = sum(this%ldos(:, ispden)) * this%dvol
            end do
        end if

        !Local polarizability
        if (this%iprcel == 203) then
            if (istep == 1) then    ! TODO : More option to control when the ldos is updated
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
            ABI_FREE(this%tdos)
        end if

        if (this%iprcel == 203) then
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
    !! Get the vector r (in reduced coordinates) of index ifft
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
    !! Get the vector g (in reduced coordinates) of index ifft
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
        g(1) = i1-1
        g(2) = i2-1
        g(3) = i3-1

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
        class(precon_object), intent(inout) :: this
        integer :: ispden
        integer, intent(in) :: ngfft(:)

        !Local variables-------------------------------
        logical :: exist
        integer :: io, n, i

        ! *************************************************************************
       
        if (this%iprcel==202) then
            n = size(this%ldos(:, ispden))
            ! Writing the file
            inquire(file="ldos.txt", exist=exist)
            if (exist) then
                open(newunit=io, file="ldos.txt", status="replace", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%ldos(i, ispden)
                end do
                close(io)
            else
                open(newunit=io, file="ldos.txt", status="new", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%ldos(i, ispden)
                end do
                close(io)
            end if 
        end if 

        if (this%iprcel==203) then
            n = size(this%loc_pola(:, ispden))
            ! Writing the file
            inquire(file="loc_pola.txt", exist=exist)
            if (exist) then
                open(newunit=io, file="loc_pola.txt", status="replace", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%loc_pola(i, ispden)
                end do
                close(io)
            else
                open(newunit=io, file="loc_pola.txt", status="new", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), this%loc_pola(i, ispden)
                end do
                close(io)
            end if 
        end if 

    end subroutine precon_save

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
           delta = 1/(exp(x/2)+exp(-x/2))**2
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
        
         fprim = 1/tsmear * delta
        
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
        integer :: mband, i, mcg, maxocc
        !arrays
        real(dp), allocatable :: ldos_wheights(:)

        !dummy mkrho arguments
        real(dp)                :: rhog(2, nfft)
        type(paw_dmft_type)     :: paw_dmft
        type(wvl_wf_type)       :: wvl_wfs
        type(wvl_denspot_type)  :: wvl_den

        ! *************************************************************************

        !compute wheights
        mband = size(eigen)
        ABI_MALLOC(ldos_wheights, (mband))
        maxocc = two / (dtset%nsppol * dtset%nspinor)   !Maximum number of occupations (1 or 2)
        do i=1, mband
            ldos_wheights(i) = derivative_occ(dtset%occopt, eigen(i), fermie, dtset%tsmear) * maxocc
        end do

        !Compute ldos using mkrho with ldos_wheights in place of the occupations
        mcg = size(cg)
        paw_dmft%use_dmft = 0
        paw_dmft%use_sc_dmft = 0
        call mkrho(cg, dtset, gprimd, irrzon, kg, mcg, mpi_enreg, npwarr, ldos_wheights, &
        &   paw_dmft, phnons, rhog, ldos, rprimd, 0, ucvol, wvl_den, wvl_wfs, option=0)

        ABI_FREE(ldos_wheights)

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
                !write(6,*)'    compute_loc_pola : rhor0_atom', rhor0_atom; flush(6) !DEBUG

                !2) Computing norm(r-r_atom)^2 (periodized with sin)
                !TODO : better + mistake (pic dans loc_pola)
                do ir = 1, nfft
                    r2(ir) = norm2(matmul(rprimd, 1/two_pi*sin(two_pi * (get_r_vector(ir, dtset%ngfft) - r_atom)))) ** 2
                end do
                !write(6,*)'    compute_loc_pola : r2', r2; flush(6) !DEBUG

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

    subroutine save_applied_op(this, ngfft, optspace, vec, op_vec)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        integer, intent(in) :: optspace ! 1 for real space, 2 for Fourier space
        !arrays
        real(dp), intent(in) :: vec(:), op_vec(:)
        integer, intent(in) :: ngfft(:)
       
        !Local variables-------------------------------
        !scalars
        logical :: exist
        integer :: io, n, i
        
        ! *************************************************************************
        n = size(vec)
        ! Writing the file
        inquire(file="applied_op.txt", exist=exist)
        if (exist) then
            open(newunit=io, file="applied_op.txt", status="replace", action="write")
            do i=1,n
                if(optspace==1) then
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), vec(i), op_vec(i)
                else if(optspace==2) then
                    write (io, '(*(G0.6,:,","))') matmul(this%gprimd, get_g_vector(i, ngfft)), vec(i), op_vec(i)
                end if
            end do
            close(io)
        else
            open(newunit=io, file="applied_op.txt", status="new", action="write")
            do i=1,n
                if(optspace==1) then
                    write (io, '(*(G0.6,:,","))') matmul(this%rprimd, get_r_vector(i, ngfft)), vec(i), op_vec(i)
                else if(optspace==2) then
                    write (io, '(*(G0.6,:,","))') matmul(this%gprimd, get_g_vector(i, ngfft)), vec(i), op_vec(i)
                end if
            end do
            close(io)
        end if 
    end subroutine save_applied_op

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
    !!  vec_g (2, :) = Vector (in G-space) to which the model chi0 operator is applied (in place).
    !!
    !! SOURCE
    subroutine apply_chi0(this, mpi_enreg, ngfft, ispden, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: ispden
        !arrays
        integer, intent(in) :: ngfft(:)
        real(dp), intent(inout) :: vec_g(2, this%nfft)
       
        !Local variables-------------------------------
        !scalars
        integer :: size_vec, cplex
        integer :: i, i_g
        !arrays
        real(dp), allocatable :: work_r(:)
        real(dp), allocatable :: work_g(:, :), vec_g_saved(:, :)
        real(dp) :: g(3)
        
        ! *************************************************************************
       

        if (this%iprcel == 201) then
        !Kerker
            vec_g = (-1/(4*pi*(this%dielng)**2)) * vec_g
        end if
        
        if (this%iprcel == 202) then
        !LDOS model
            !1) ifft to get vec in the real space
            size_vec = size(vec_g, 2)
            cplex = 1                           ! TODO : cplex as argument ?
            ABI_MALLOC(work_r, (cplex*size_vec))
            call fourdp(cplex, vec_g, work_r, 1, mpi_enreg, size_vec, 1, ngfft, 0)

            !2) chi0(v)(r) = -ldos(r)*v(r) + 1/dos * ldos(r)*integral(ldos(r')*v(r')*dr')
            work_r = -this%ldos(:, ispden)*work_r &
                         + 1/this%tdos(ispden) * dot_product(this%ldos(:, ispden), work_r)*this%dvol * this%ldos(:, ispden)

            !3) fft to get vec back in the reciprocal space
            call fourdp(cplex, vec_g, work_r, -1, mpi_enreg, size_vec, 1, ngfft, 0)
            ABI_FREE(work_r)
        end if

        if (this%iprcel == 203) then
        !Local polarizability model
            ABI_MALLOC(work_g, (2, this%nfft))
            ABI_MALLOC(vec_g_saved, (2, this%nfft))
            vec_g_saved = vec_g
            vec_g = 0.0
            do i=1, 3
                work_g = 0.0
                !1) Multiplication by d_i(g) in reciprocal space
                do i_g = 1, this%nfft
                    !g = two_pi * matmul(this%gprimd, g_vectors(:, i_g))
                    g = two_pi * matmul(this%gprimd, get_g_vector(i_g, ngfft))
                    work_g(:, i_g) = complex_mult( [0.0_dp, g(i)/sqrt(1+(norm2(g)/this%gc)**2)], vec_g_saved(:, i_g) )
                end do

                !2) Multiplication by the local polarizability in real space
                call fourdp(cplex, work_g, work_r, 1, mpi_enreg, size_vec, 1, ngfft, 0) !ifft
                work_r = this%loc_pola(:, ispden) * work_r                                         !local multiplication
                call fourdp(cplex, work_g, work_r, -1, mpi_enreg, size_vec, 1, ngfft, 0) !fft

                !3) Multiplication by d_i(g) in reciprocal space
                do i_g = 1, this%nfft
                    g = two_pi * matmul(this%gprimd, get_g_vector(i_g, ngfft))
                    work_g(:, i_g) = complex_mult( [0.0_dp, g(i)/sqrt(1+(norm2(g)/this%gc)**2)], work_g(:, i_g) )
                end do

                vec_g = vec_g + work_g

            end do
        end if

    end subroutine apply_chi0

end module m_precon
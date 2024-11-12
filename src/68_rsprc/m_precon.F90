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
    
    use defs_basis
    use m_dtset
    use m_fft,      only : fourdp
    use defs_abitypes, only : MPI_type

    use m_mkrho
    use m_paw_dmft
    use defs_wvltypes
    use defs_wvltypes
    
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
        integer  :: occopt
        real(dp) :: tsmear
        real(dp), pointer :: fermie
        real(dp), pointer :: cg(:, :), eigen(:), phnons(:, :, :)
        integer, pointer  :: kg(:, :), npwarr(:), irrzon(:, :, :)
        real(dp), allocatable :: ldos(:, :), tdos(:)
        !For ffts :
        integer  :: nfft

    contains
        procedure :: init => precon_init
        procedure :: update => precon_update
        procedure :: free => precon_free
        procedure :: save_ldos => save_ldos
        procedure :: apply_chi0 => apply_chi0

    end type precon_object

contains 

    !****f* m_precon/update
    !! NAME
    !! update
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
    subroutine precon_init(this, dtset, gprimd, rprimd, ucvol, cg, eigen, fermie, irrzon, kg, npwarr, phnons)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        !scalars
        type(dataset_type),intent(in) :: dtset
        real(dp), intent(in) :: ucvol
        real(dp), intent(in), target :: fermie

        !arrays
        real(dp), intent(in) :: gprimd(:, :), rprimd(:, :)
        real(dp), intent(in), target :: cg(:, :), eigen(:), phnons(:, :, :)
        integer, intent(in), target  :: irrzon(:, :, :), kg(:, :), npwarr(:)

        ! *************************************************************************
        !Constant data from dtset
        this%dielng = dtset%dielng
        this%diemix = dtset%diemix
        this%iprcel = dtset%iprcel
        this%nfft   = dtset%nfft
        this%nspden = dtset%nspden
        this%occopt = dtset%occopt
        this%tsmear = dtset%tsmear
        !Other constants
        this%gprimd = gprimd
        this%rprimd = rprimd
        this%ucvol  = ucvol
        this%dvol   = ucvol/this%nfft ! factor for integrals in real space: sum(f) * dvol ~ integral f
        !Pointers
        this%cg     => cg
        this%eigen  => eigen
        this%fermie => fermie
        this%irrzon => irrzon
        this%kg     => kg
        this%npwarr => npwarr
        this%phnons => phnons
        !Initializing LDOS specific variables
        if (this%iprcel == 202) then
            !Allocating the array containing ldos and tdos
            ABI_MALLOC(this%ldos, (this%nfft, this%nspden))
            ABI_MALLOC(this%tdos, (this%nspden))
        end if

    end subroutine precon_init

    !****f* m_precon/update
    !! NAME
    !! update
    !!
    !! FUNCTION
    !! Updates the precon_object :
    !!      For LDOS preconditioning (iprcel=...) : 
    !!          computes the new ldos (local density of state) with current wavefunctions 
    !!          and the new tdos (total density of state = integral of ldos) 
    !!
    !! INPUTS
    !!  dtset     = all input variables for this dataset
    !!  mpi_enreg = information about MPI parallelization
    !!
    !! SOURCE
    subroutine precon_update(this, dtset, mpi_enreg)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        type(dataset_type), intent(in) :: dtset
        type(MPI_type), intent(inout) :: mpi_enreg
        
        !Local variables-------------------------------
        integer :: ispden

        ! *************************************************************************

        if (this%iprcel == 202) then 
            !update ldos
            call compute_ldos(this%occopt, this%eigen, this%fermie, this%tsmear, this%cg,   &
            &   dtset, this%ucvol, this%rprimd, this%gprimd,                                &
            &   this%irrzon, this%kg, this%npwarr, this%phnons,                             &
            &   this%nfft, mpi_enreg,                                                       &
            &   this%ldos)
            !update tdos
            do ispden=1,this%nspden
                this%tdos = sum(this%ldos(:, ispden)) * this%dvol
            end do
        end if
       
    end subroutine precon_update

    subroutine precon_free(this)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        ! *************************************************************************
       
        if (this%iprcel == 202) then
            !Deallocating the array containing ldos and tdos
            ABI_FREE(this%ldos)
            ABI_FREE(this%tdos)
        end if
    end subroutine precon_free

    subroutine save_ldos(this, ispden)

        !Arguments ------------------------------------
        class(precon_object), intent(inout) :: this
        integer :: ispden

        !Local variables-------------------------------
        logical :: exist
        integer :: io, n, i

        ! *************************************************************************
       
        if (this%iprcel==202) then
            n = size(this%ldos)
            inquire(file="ldos.txt", exist=exist)
            if (exist) then
                open(newunit=io, file="ldos.txt", status="replace", action="write")
                do i=1,n
                    ! TODO : compute r 
                    write (io, '(*(G0.6,:,","))') this%ldos(i, ispden)
                end do
                close(io)
            else
                open(newunit=io, file="ldos.txt", status="new", action="write")
                do i=1,n
                    write (io, '(*(G0.6,:,","))') this%ldos(i, ispden)
                end do
                close(io)
            end if 
        end if 

    end subroutine save_ldos

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
    subroutine compute_ldos(occopt, eigen, fermie, tsmear, cg,  &
        &   dtset, ucvol, rprimd, gprimd,                       &
        &   irrzon, kg, npwarr, phnons,                         &
        &   nfft, mpi_enreg,                                    &
        &   ldos)

        !Arguments ------------------------------------
        !scalars
        integer, intent(in) :: occopt, nfft
        real(dp), intent(in) :: tsmear, fermie
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
        integer :: mband, i, mcg
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
        do i=1, mband
            ldos_wheights(i) = derivative_occ(occopt, eigen(i), fermie, tsmear)
        end do

        mcg = size(cg)
        paw_dmft%use_dmft = 0
        paw_dmft%use_sc_dmft = 0
        !compute ldos using mkrho with ldos_wheights in place of the occupations
        call mkrho(cg, dtset, gprimd, irrzon, kg, mcg, mpi_enreg, npwarr, ldos_wheights, &
        &   paw_dmft, phnons, rhog, ldos, rprimd, 0, ucvol, wvl_den, wvl_wfs, option=0)

        ABI_FREE(ldos_wheights)

    end subroutine compute_ldos

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
    subroutine apply_chi0(this, mpi_enreg, ngfft, ispden, g_vectors, vec_g)

        !Arguments ------------------------------------
        class(precon_object), intent(in) :: this
        !scalars
        type(MPI_type), intent(in) :: mpi_enreg
        integer, intent(in) :: ispden
        !arrays
        integer, intent(in) :: ngfft(:)
        integer, intent(in) :: g_vectors(:, :)
        real(dp), intent(inout) :: vec_g(2, this%nfft)
       
        !Local variables-------------------------------
        !scalars
        integer :: size_vec, cplex
        !real(dp) :: 
        real(dp), allocatable :: work_r(:)
        
        ! *************************************************************************
       

        if (this%iprcel == 201) then
        !Kerker
            vec_g = (-1/(4*pi*(this%dielng)**2)) * vec_g
        else if (this%iprcel == 202) then
        !LDOS
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


    end subroutine apply_chi0

end module m_precon
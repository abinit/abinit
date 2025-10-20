!!****m* ABINIT/m_frohlich
!! NAME
!!  m_frohlich
!!
!! FUNCTION
!!  Description
!!
!! COPYRIGHT
!!  Copyright (C) 2018-2025 ABINIT group (VV, XG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_frohlich

 use defs_basis
 use m_abicore
 use m_errors
 use m_crystal
 use m_ebands
 use m_efmas_defs
 use m_ifc
 use m_dtset

 use m_fstrings,            only : sjoin, itoa
 use m_gaussian_quadrature, only : cgqf

 implicit none

 private
!!***

!----------------------------------------------------------------------

!!****t* m_frohlich/frohlich_t
!! NAME
!!  frohlich_t
!!
!! FUNCTION
!!  Description
!!
!! SOURCE

 type,public :: frohlich_t

  integer :: kind = 0
   ! Type of the Fr\"ohlich model
   ! used to access the calculations of various properties
   ! 0 -> ndeg != 1 && ndeg != 3, generalized Fr\"ohlich model
   ! 1 -> ndeg = 1, standard (possibly anisotropic) Fr\"ohlich model
   ! 2 -> ndeg = 3, cubic generalized Fr\"ohlich model with 3-fold degeneracy

 ! Geometry ------------------------------------

  real(dp) :: ucvol
   ! Real space unit cell volume

  real(dp) :: gmet(3,3)
    ! Reciprocal space metric


 ! Electroncic subspace -----------------------

  logical :: isinitel = .false.
   ! Flag indicating that the electronic subspace has been initialized

  real(dp) :: kpt(3)
   ! k-point characterizing the electronic subspace, i.e. the k-point at which
   ! the effective mass tensor is obtained (usually, CBM or VBM)

  integer :: ndeg
   ! Number of degenerate bands taken into account

  complex(dp), allocatable :: eig2_diag_cart(:,:,:,:)
   ! Band curvature double tensor in Cartesian coordinates
   ! (3, 3, ndeg, ndeg)

  real(dp) :: band_params(3)
   ! Parameters describing electronic bands
   ! The meaning of this arrayd depends on the value of the kind variable
   ! kind = 0 -> undefined
   ! kind = 1 -> inverse effetive masses along the 100, 010 and 001 directions
   ! kind = 2 -> Luttinger-Kohn parameters A, B, C

  logical :: lutt_warn(3) = .false.

  real(dp), allocatable :: sqrt_efmas_avg(:)
   ! Square root of effective mass averaged over q-sphere for each band
   ! (ndeg)

  logical, allocatable :: saddle_warn(:)
   ! Signals if the k-point characterizing the electronic subspace is a
   ! saddle-point for each of the degenerate bands
   ! (ndeg)

  real(dp) :: sqrt_efmas_tot
   ! Total square root of effective mass average

  real(dp), allocatable :: invefmas(:,:)
   ! Inverse electronic effective masses along
   ! kind = 1 -> 100, 010, 001 directions
   ! kind = 2 -> 100, 110, 111 directions
   ! (ndeg, 3)


 ! Vibrational (phonon) subspace ---------------

  logical :: isinitph = .false.
   ! Flag indicating that the phonon subspace has been initialized

  integer :: natom
   ! Number of atoms in the cell

  integer :: nqdir
   ! Number of points for spherical integration used to compute ZPR and other
   ! quantities

  real(dp), allocatable :: unit_qdir(:,:)
   ! Unit q-vectors representing reciprocal space directions used to compute
   ! ZPR and other quantities
   ! (3, nqdir)

  real(dp), allocatable :: weights_qdir(:)
   ! Gaussian quadrature weights used for spherical intragraion over q-vectors
   ! (nqdir)

  real(dp), allocatable :: dielt_qdir(:)
   ! High-frequency dielectric constant for each q-vector direction
   ! (nqdir)

  real(dp), allocatable :: phfreq_qdir(:,:)
   ! Phonon frequencies for each mode and q-vector direction
   ! (3*natom, nqdir)

  real(dp), allocatable :: polarity_qdir(:,:,:)
   ! Mode polarity vectors for each mode and q-vector direction
   ! (3, 3*natom, nqdir)

  real(dp), allocatable :: proj_polarity_qdir(:,:)
   ! Projections of mode polarity vectors for each mode and q-vector direction
   ! (3*natom, nqdir)

  real(dp), allocatable :: investar(:,:)
   ! Inverse effective dielectric constant for each mode and q-vector
   ! (3*natom, nqdir)

  real(dp), allocatable :: dielavg(:)
   ! Dielectric average over q-vectors for each mode (Eq. (26) of [deMelo2023])
   ! (3*natom)

  logical, allocatable :: isiractive(:)
   ! Flags to detect the infrared-active phonon modes
   ! (3*natom)

  real(dp) :: dielt_eff
   ! Effective dielectric constant in the strong-coupling regime

  real(dp) :: phfreq_eff
   ! Effective LO phonon frequency in the strong-coupling regime


 ! Weak-coupling parameters --------------------

  real(dp), allocatable :: zpr_band(:)
   ! Zero-point renormalization energy for each band (Eq. (17) of [deMelo2023])
   ! (ndeg)

  logical :: sign_warn = .false.
   ! Sginals an error if a saddle-point is encountered or bands contribute to
   ! the ZPR with different signs

  real(dp) :: zpr_gamma
   ! Correction to the ZPR taking into account the infrared divergence of the
   ! electron-phonon coupling in the Fr\"ohlich model

  real(dp) :: zpr
   ! Total ZPR for these bands and k-point

  real(dp) :: enpol_wc
   ! Fr\"ohlich polaron formation energy in the weak-coupling regime
   ! (Eq. (26) of [deMelo2023])

  real(dp), allocatable :: zpr_k(:,:,:)
   ! Direction dependent ZPR (Eq. (86) of [Guster2021])
   ! (3, ndeg, 3)

  real(dp), allocatable :: invpolmas(:,:)
   ! Inverse Fr\"ohlich polaron effective masses in the weak-coupling regime
   ! (Sec. III A, B of [Guster2021])
   ! kind = 1 -> 100, 010, 001 directions
   ! kind = 2 -> 100, 110, 111 directions
   ! (ndeg, 3)


  contains

    procedure :: init_ph => frohlich_init_ph
     ! Initialization of a vibrational (phonon) subspace parameters

    procedure :: init_el => frohlich_init_el
     ! Initialization of an electronic subspace parameters

    procedure :: free_ph => frohlich_free_ph
     ! Free memory allocated to the vibrational subspace

    procedure :: free_el => frohlich_free_el
     ! Free memory allocated to the electronic subspace
     ! and other related quantities

    procedure :: calc_zpr => frohlich_calc_zpr
     ! Calculate the zero-point renormalization energy corresponding to the
     ! weak-coupling treatment of the Fr\"ohlich model

    procedure :: calc_polaronmass => frohlich_calc_polaronmass
     ! Calculate the polaron effective mass corresponding to the weak-coupling
     ! treatment of the Fr\"ohlich model; available for kind = 1 or 2

    ! procedure :: ncwrite => frohlich_ncwrite
     ! Write main dimensions and header on a netcdf file


 end type frohlich_t
!!***

public :: frohlichmodel_zpr         ! Main routine to compute ZPR
public :: frohlichmodel_polaronmass ! Main routine to compute polaron effective
                                    ! masses

contains !=====================================================================
!!***

!!****f* m_frohlich/frohlichmodel_polaronmass
!! NAME
!!  frohlichmodel_polaronmass
!!
!! FUNCTION
!!  Main routine to compute the polaron effective masses of the generalized
!! Fr\"ohlich model and other related quantities
!!
!! INPUTS
!!  cryst<crystal_t>=Structure defining the unit cell
!!  dtset<dataset_type>=All input variables for this dataset.
!!  efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band
!! degeneracy at each k point
!!  efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!   efmasval(:,:)%eig2_diag band curvature double tensor
!!  ifc<ifc_type>=contains the dynamical matrix and the IFCs.
!!
!! NOTES
!!  This routine has to be merged with the frohlichmodel_zpr routine, and their
!!  text output needs to be refined. For now, they are being kept to be
!!  compatible with the legacy unit tests.
!!
!! SOURCE

subroutine frohlichmodel_polaronmass(frohlich, cryst, dtset, efmasdeg, efmasval, ifc)

!Arguments ------------------------------------
!scalars
 class(frohlich_t),intent(inout) :: frohlich
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ifc_type),intent(in) :: ifc
!arrays
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables-------------------------------
!scalar
 integer :: nu, ikpt, ideg, ndeg
 integer :: iqdir
!arrays
 real(dp) :: kpt(3)

! *************************************************************************

 ! Initalize phonon and dielectric subspace
 call frohlich%init_ph(cryst, dtset%efmas_ntheta, ifc)

 ! For each k-point (with possible degeneracy), initialize an inistance of the
 ! generalize Fr\"ohlich model and calculate related quantities
 do ikpt=1,dtset%nkpt
   kpt(:) = dtset%kptns(:, ikpt)

   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)
     ndeg = efmasdeg(ikpt)%degs_bounds(2, ideg) - &
       efmasdeg(ikpt)%degs_bounds(1, ideg) + 1

     call frohlich%init_el(cryst, kpt, ndeg, efmasval(ideg, ikpt)%eig2_diag)
     call frohlich%calc_zpr()
     call frohlich%calc_polaronmass()

     ! Luttinger parameters if applicable
     if (frohlich%ndeg == 3) then

       if (.not. (any(frohlich%saddle_warn))) then
         if (any(frohlich%lutt_warn)) then

           write(ab_out, '(2a)') ch10, &
             ' Luttinger parameters could not be determined:'

           if (frohlich%lutt_warn(1)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (100) direction.'
           endif

           if (frohlich%lutt_warn(2)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (111) direction.'
           endif

           if (frohlich%lutt_warn(3)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (110) direction.'
           endif

           write(ab_out, '(a)') ch10
         else
           write(ab_out, '(a,3f14.6)') &
             '   Luttinger parameters (A, B, C) (a.u.): ', &
             frohlich%band_params(:)
         endif
       endif

       ! Print inverse electronic effective masses in the output
       write(ab_out, '(a)') repeat('-', 80)
       write(ab_out, '(a)') '   Polaron properties from the generalized &
         &Froehlich model'
       write(ab_out, '(a)') repeat('-', 80)
       write(ab_out, '(a)') '   Polar modes'
       write(ab_out, '(a)') '   ##      Frequency(meV)            Epsilon*'

       ! For a cubic material, dielectric tensor and phonon frequencies at Gamma
       ! do not depend on the q-vector direction
       iqdir = 1
       do nu = 4,3*cryst%natom
         if (frohlich%isiractive(nu)) then
          write(ab_out,'(2x,i3,5x,f15.6,5x,f15.6)') nu, &
            frohlich%phfreq_qdir(nu, iqdir)*Ha_eV*1000.0_dp, &
            one/frohlich%investar(nu, iqdir)
         endif
       enddo

       write(ab_out, '(a)') ' '
       write(ab_out, '(a,f10.2)') '   ZPR (meV): ', &
         frohlich%zpr_k(1, 1, 1)*Ha_eV*1000.0_dp
       write(ab_out, '(a)') ' '
       write(ab_out, '(a)') '   Electronic effective mass (a.u.) &
         &along 3 directions'
       write(ab_out, '(a, 3f15.6)')'    Direction 100:         ', &
         one/frohlich%invefmas(:, 1)
       write(ab_out, '(a, 3f15.6)')'    Direction 110:         ', &
         one/frohlich%invefmas(:, 2)
       write(ab_out, '(a, 3f15.6)')'    Direction 111:         ', &
         one/frohlich%invefmas(:, 3)

     ! Print inverse polaron effective masses in the output
       write(ab_out, '(a)') ' '
       write(ab_out, '(a)') '   Polaron effective mass (a.u.) along 3 directions'
       write(ab_out, '(a, 3f15.6)') '    Direction 100:         ', &
         one/frohlich%invpolmas(:, 1)
       write(ab_out, '(a, 3f15.6)') '    Direction 110:         ', &
         one/frohlich%invpolmas(:, 2)
       write(ab_out, '(a, 3f15.6)') '    Direction 111:         ', &
         one/frohlich%invpolmas(:, 3)
       write(ab_out, '(a)')' '
       write(ab_out, '(a)')'   Sum rule of inverse polaron masses check-up &
         &(for convergence purposes):'
       write(ab_out,'(a, 3f15.6)')'    Direction 100:         ', &
         sum(frohlich%invpolmas(:, 1))
       write(ab_out,'(a, 3f15.6)')'    Direction 110:         ', &
         sum(frohlich%invpolmas(:, 2))
       write(ab_out,'(a, 3f15.6)')'    Direction 111:         ', &
         sum(frohlich%invpolmas(:, 3))
     endif

     call frohlich%free_el()
   enddo
 enddo

 call frohlich%free_ph()

end subroutine frohlichmodel_polaronmass
!!***


!!****f* m_frohlich/frohlichmodel_zpr
!! NAME
!!  frohlichmodel_zpr
!!
!! FUNCTION
!!  Main routine to compute the ZPR of the generalized Fr\"ohlich model and
!! other related quantities
!!
!! INPUTS
!!  cryst<crystal_t>=Structure defining the unit cell
!!  dtset<dataset_type>=All input variables for this dataset.
!!  efmasdeg(nkpt_rbz) <type(efmasdeg_type)>= information about the band
!! degeneracy at each k point
!!  efmasval(mband,nkpt_rbz) <type(efmasdeg_type)>= double tensor datastructure
!!   efmasval(:,:)%eig2_diag band curvature double tensor
!!  ifc<ifc_type>=contains the dynamical matrix and the IFCs.
!!
!! SOURCE

subroutine frohlichmodel_zpr(frohlich, cryst, dtset, efmasdeg, efmasval, ifc)

!Arguments ------------------------------------
!scalars
 class(frohlich_t),intent(inout) :: frohlich
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ifc_type),intent(in) :: ifc
!arrays
 type(efmasdeg_type), intent(in) :: efmasdeg(:)
 type(efmasval_type), intent(in) :: efmasval(:,:)

!Local variables-------------------------------
!scalar
 integer :: nu, ikpt, ideg, ndeg, iband
!arrays
 real(dp) :: kpt(3)

! *************************************************************************

 ! Initalize phonon and dielectric subspace
 call frohlich%init_ph(cryst, dtset%efmas_ntheta, ifc)

 ! Dielectric average
 write(ab_out, '(a)') repeat('-', 80)
 write(ab_out, '(a)') ' Dielectric average (EQ. 25 Melo2022)'
 write(ab_out, '(a)') repeat('-', 80)
 write(ab_out, '(a)') &
   ' Mode    <1/epsilon*SQRT(w_LO/2)>              Cumulative sum'

 do nu=4,3*cryst%natom
   write(ab_out, '(i5,f28.12,f28.12)') nu, &
     frohlich%dielavg(nu), sum(frohlich%dielavg(1:nu))
 enddo
 write(ab_out, '(a)') repeat('-', 80)

 ! Infrared ZPR correction (does not depend on the electronic subspace)
 write(ab_out, '(6a,f14.6,a,f14.6,a)') ch10, &
   ' Rough correction to the ZPR, to take into account the missing q=0 piece &
   &using Frohlich model:', ch10, &
   ' (+ for occupied states, - for unoccupied states) * zpr_q0_fact / &
   &(Nqpt_full_bz)**(1/3) ', ch10, &
   ' where Nqpt_full_bz=number of q wavevectors in full BZ, and zpr_q0_fact=',&
   frohlich%zpr_gamma, ' Ha=', frohlich%zpr_gamma*Ha_eV, ' eV'

 ! For each k-point (with possible degeneracy), initialize an inistance of the
 ! generalize Fr\"ohlich model and calculate ZPR and other related quantities
 do ikpt=1,dtset%nkpt
   kpt(:) = dtset%kptns(:, ikpt)

   do ideg=efmasdeg(ikpt)%deg_range(1),efmasdeg(ikpt)%deg_range(2)
     ndeg = efmasdeg(ikpt)%degs_bounds(2, ideg) - &
       efmasdeg(ikpt)%degs_bounds(1, ideg) + 1

     call frohlich%init_el(cryst, kpt, ndeg, efmasval(ideg, ikpt)%eig2_diag)
     call frohlich%calc_zpr()

     ! Print ZPR results for each k-point
     if (ndeg == 1) then
       write(ab_out, '(2a,3(f6.3,a),i5)') ch10, &
         ' - At k-point (', kpt(1), ',', kpt(2), ',', kpt(3), '), band ', &
         efmasdeg(ikpt)%degs_bounds(1,ideg)
     else
       write(ab_out, '(2a,3(f6.3,a),i5,a,i5)') ch10, &
         ' - At k-point (', kpt(1), ',', kpt(2), ',', kpt(3), '), bands ', &
         efmasdeg(ikpt)%degs_bounds(1,ideg), ' through ', &
         efmasdeg(ikpt)%degs_bounds(2,ideg)
     endif

     ! Luttinger parameters if applicable
     if (frohlich%kind == 2) then

       if (.not. (any(frohlich%saddle_warn))) then
         if (any(frohlich%lutt_warn)) then

           write(ab_out, '(2a)') ch10, &
             ' Luttinger parameters could not be determined:'

           if (frohlich%lutt_warn(1)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (100) direction.'
           endif

           if (frohlich%lutt_warn(2)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (111) direction.'
           endif

           if (frohlich%lutt_warn(3)) then
             write(ab_out, '(a)')  '     Predicted degeneracies for &
               &deg_dim = 3 are not met for (110) direction.'
           endif

           write(ab_out, '(a)') ch10
         else
           write(ab_out, '(a,3f14.6)') &
             ' Luttinger parameters (A, B, C) [at. units]: ', &
             frohlich%band_params(:)
         endif
       endif
     endif

     ! Effective mass average and ZPR
     do iband=1,ndeg
       if (frohlich%saddle_warn(iband)) then
         write(ab_out, '(a,i5,a)') ' Band ', &
           efmasdeg(ikpt)%degs_bounds(1, ideg) + iband - 1, ' SADDLE POINT - &
           &Frohlich effective mass and ZPR cannot be defined. '
       else
         write(ab_out, '(a,i5,a,f14.10)') ' Band ', &
           efmasdeg(ikpt)%degs_bounds(1, ideg) + iband - 1, ' Angular average &
           &effective mass for Frohlich model (<m**0.5>)**2= ', &
           sign(frohlich%sqrt_efmas_avg(iband)**2, &
           frohlich%sqrt_efmas_avg(iband))
       endif
     enddo

     if (.not. frohlich%sign_warn) then
       write(ab_out, '(2a)') &
         ' Angular and band average effective mass and ZPR for Frohlich model.'

       write(ab_out, '(a,es16.6)') ' Value of     (<<m**0.5>>)**2 = ', &
         (sum(abs(frohlich%sqrt_efmas_avg(:))) / ndeg)**2

       write(ab_out, '(a,es16.6)') ' Absolute Value of <<m**0.5>> = ', &
         (sum(abs(frohlich%sqrt_efmas_avg(:))) / ndeg)

       write(ab_out, '(a,es16.6,a,es16.6,a)') &
         ' ZPR from Frohlich model      = ', frohlich%zpr, ' Ha=', &
         frohlich%zpr*Ha_eV,' eV'
     else
       write(ab_out, '(a)') ' Angular and band average effective mass for &
         &Frohlich model cannot be defined because of a sign problem.'
     endif

     call frohlich%free_el()
   enddo
 enddo

 call frohlich%free_ph()

end subroutine frohlichmodel_zpr
!!***


!!****f* m_frohlich/frohlich_free_el
!! NAME
!!  frohlich_free_el
!!
!! FUNCTION
!!  Deallocate dynamic memory related to the electronic subspace and other
!! related quantities
!!
!! INPUTS
!!
!! SOURCE

subroutine frohlich_free_el(self)

!Arguments ------------------------------------
 class(frohlich_t),intent(inout) :: self
! *************************************************************************

 self%isinitel = .false.

 ! real
 ABI_SFREE(self%sqrt_efmas_avg)
 ABI_SFREE(self%zpr_band)
 ABI_SFREE(self%invpolmas)
 ABI_SFREE(self%zpr_k)
 ABI_SFREE(self%invefmas)

 ! complex
 ABI_SFREE(self%eig2_diag_cart)

 ! logical
 ABI_SFREE(self%saddle_warn)

end subroutine frohlich_free_el
!!***


!!****f* m_frohlich/frohlich_free_ph
!! NAME
!!  frohlich_free_ph
!!
!! FUNCTION
!!  Deallocate dynamic memory related to the vibrational subspace
!!
!! INPUTS
!!
!! SOURCE

subroutine frohlich_free_ph(self)

!Arguments ------------------------------------
 class(frohlich_t),intent(inout) :: self
! *************************************************************************

 self%isinitph = .false.

 ! real
 ABI_SFREE(self%unit_qdir)
 ABI_SFREE(self%weights_qdir)
 ABI_SFREE(self%dielt_qdir)
 ABI_SFREE(self%phfreq_qdir)
 ABI_SFREE(self%polarity_qdir)
 ABI_SFREE(self%proj_polarity_qdir)
 ABI_SFREE(self%investar)
 ABI_SFREE(self%dielavg)

 ! logical
 ABI_SFREE(self%isiractive)

end subroutine frohlich_free_ph
!!***


!!****f* m_frohlich/frohlich_calc_polaronmass
!! NAME
!!  frohlich_calc_polaronmass
!!
!! FUNCTION
!!  Description
!!
!! INPUTS
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! OUTPUT
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! SOURCE

subroutine frohlich_calc_polaronmass(self)

!Arguments ------------------------------------
!scalars
 class(frohlich_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: ii
 integer :: ikdir, nkdir
 integer :: ik, nkgrid
 integer :: lwork, info
 integer :: nu, iband
 integer :: ixi, nxi
 integer :: iqdir
 integer :: sigma
 real(dp) :: deltak, minefmas
 real(dp) :: xi, qlen
 real(dp) :: tmp
!arrays
 real(dp) :: id33(3, 3)
 real(dp) :: kpt(3), unit_kdir(3, 3)
 real(dp) :: k_plus_q(3)
 real(dp), allocatable :: eigenvec(:,:), eigenval(:)
 real(dp), allocatable :: invefmas(:,:)
 real(dp), allocatable :: lutt_eigenvec(:,:,:,:), lutt_eigenval(:,:,:)
 real(dp), allocatable :: lk_h(:,:)
 real(dp), allocatable :: phfreq(:), investar(:)
 real(dp), allocatable :: intsum(:,:,:,:)
 real(dp), allocatable :: zpr(:,:,:), zpr_ddk(:,:)
 real(dp), allocatable :: intsuminv(:,:)
 complex(dp), allocatable :: work(:)
! *************************************************************************

 ! TODO: check if the electronic and phonon parts are initialized
 ABI_MALLOC(self%invpolmas, (self%ndeg, 3))

 if (self%kind == 2) then
   ! Polaron effective masses in the triply-degenerate case
   ! (Eq. (86) of [Guster2021])

   ! Initialization of the diagonalization routine
   ! TODO: this part is common for many routines -> initialize once and store
   ! necessary variables as local state constants in the datatype
   ABI_MALLOC(eigenvec, (self%ndeg, self%ndeg))
   ABI_MALLOC(eigenval, (self%ndeg))
   ABI_MALLOC(work, (3*self%ndeg - 2))
   lwork = -1
   call dsyev('V', 'U', self%ndeg, eigenvec(:,:), self%ndeg, eigenval(:), &
     work, lwork, info)
   lwork = int(work(1))
   ABI_FREE(work)
   ABI_MALLOC(work, (lwork))

   ! Inverse effectiv mass tensor for 100, 110 and 111 directions
   ! used to obtain the material-dependent characteristic wavevector length
   nkdir = 3
   unit_kdir(:, 1) = (/1, 0, 0/)
   unit_kdir(:, 2) = (/1, 1, 0/)/sqrt(2.0)
   unit_kdir(:, 3) = (/1, 1, 1/)/sqrt(3.0)
   ABI_MALLOC(invefmas, (self%ndeg, nkdir))
   ABI_MALLOC(self%invefmas, (self%ndeg, nkdir))

   ! 100: 2A, 2B (two-fold)
   invefmas(1, 1) = two*self%band_params(1)
   invefmas(2, 1) = two*self%band_params(2)
   invefmas(3, 1) = two*self%band_params(2)
   ! 110: A + B + C, 2B, A + B - C
   invefmas(1, 2) = self%band_params(1) + self%band_params(2) + &
     self%band_params(3)
   invefmas(2, 2) = min(two*self%band_params(2), self%band_params(1) + &
     self%band_params(2) - self%band_params(3))
   invefmas(3, 2) = max(two*self%band_params(2), self%band_params(1) + &
     self%band_params(2) - self%band_params(3))
   ! 111: 2/3*(A + 2B + 2C), 2/3*(A + 2B - C) (two-fold)
   invefmas(1, 3) = two*(self%band_params(1) + two*self%band_params(2) + &
     two*self%band_params(3)) / three
   invefmas(2, 3) = two*(self%band_params(1) + two*self%band_params(2) - &
     self%band_params(3)) / three
   invefmas(3, 3) = invefmas(2, 3)

   self%invefmas(:,:) = invefmas(:,:)

   ! Setting up the parameter to obtain the second derivative of the ZPR with
   ! the finite differences (FD) method (Eqs. (80), (86) of [Guster2021])
   ! FD: number of points
   nkgrid = 3
   ! FD: material-dependent length scale <--> lowest optical frequency
   iqdir = 1
   minefmas = one / maxval(abs(invefmas(:,:)))
   deltak = sqrt(two*minefmas*self%phfreq_qdir(4, iqdir) / 1000.0)

   ! FD: Luttinger eigenvalues and eigenvectors in the symmetry inequivalent
   ! cubic directions to be used in the ZPR calculations
   ABI_MALLOC(lutt_eigenvec, (self%ndeg, self%ndeg, nkgrid, nkdir))
   ABI_MALLOC(lutt_eigenval, (self%ndeg, nkgrid, nkdir))

   lutt_eigenval(:,:,:) = zero
   do ikdir=1,nkdir
     do ik=1,nkgrid
       kpt(:) = (ik - one)*deltak*unit_kdir(:, ikdir)
       call lk_hamiltonian(self%band_params(:), kpt(:), &
         lutt_eigenvec(:,:, ik, ikdir))
       call dsyev('V', 'U', self%ndeg, lutt_eigenvec(:,:, ik, ikdir), &
         self%ndeg, lutt_eigenval(:, ik, ikdir), work(:), lwork, info)
     enddo
   enddo

   ! Cubic system is assumed: phonon frequencies and permitivity at Gamma
   ! do not depend on the q-vector direction
   iqdir = 1
   ABI_MALLOC(phfreq, (3*self%natom))
   ABI_MALLOC(investar, (3*self%natom))
   phfreq(:) = self%phfreq_qdir(:, iqdir)
   investar(:) = self%investar(:, iqdir)

   ! Effective mass sign
   sigma = 1
   if (self%band_params(1) < 0) sigma = -1

   ! 3x3 Identity matrix
   id33(:,:) = zero
   do ii=1,3
     id33(ii, ii) = one
   enddo

   ! FD: main loop
   ! nqidr = 2*ntheta**2 (ntheta = nxi)
   nxi = sqrt(one*self%nqdir/2)
   ABI_MALLOC(intsum, (self%ndeg, self%ndeg, nkgrid, nkdir))
   ABI_MALLOC(zpr, (nkgrid, self%ndeg, nkdir))
   ABI_MALLOC(self%zpr_k, (nkgrid, self%ndeg, nkdir))
   ABI_MALLOC(lk_h, (self%ndeg, self%ndeg))
   ABI_MALLOC(intsuminv, (self%ndeg, self%ndeg))

   ! Values of ZPR (Eq. (86) of [Guster2021]) to be used in FD
   zpr(:,:,:) = zero
   do nu=4,3*self%natom

     ! Summation over the infrared-active phonon modes
     if (self%isiractive(nu)) then

       do ikdir=1,nkdir
         do iband=1,self%ndeg
           do ik=1,nkgrid
             kpt(:)  = (ik - one)*deltak*unit_kdir(:, ikdir)

             do ixi=0,nxi
               xi = ixi*pi / (two*nxi)
               if (ixi == nxi) xi = xi - tol8

               ! q-vector length = (omega/A)^{1/2}*tan(xi)
               qlen = sqrt(phfreq(nu) / abs(self%band_params(1)))*tan(xi)
               do iqdir=1,self%nqdir
                 k_plus_q(:) = kpt(:) + qlen*self%unit_qdir(:, iqdir)
                 call lk_hamiltonian(self%band_params(:), k_plus_q(:), &
                   lk_h(:,:))

                 intsum(:,:, ik, ikdir) = &
                   abs(lutt_eigenval(iband, ik, ikdir)*id33(:,:)) - (sigma * &
                   lk_h(:,:) + phfreq(nu)*id33(:,:))

                 call mat3inv(intsum(:,:, ik, ikdir), intsuminv(:,:))

                 tmp = dot_product(lutt_eigenvec(:, iband, 2, ikdir), &
                   matmul(intsuminv(:,:), lutt_eigenvec(:, iband, 2, ikdir)))

                 zpr(ik, iband, ikdir) = zpr(ik, iband, ikdir) + &
                   investar(nu)*phfreq(nu)*tmp*sqrt(phfreq(nu) / &
                   abs(self%band_params(1))) / cos(xi)**2 * &
                   self%weights_qdir(iqdir)
               enddo
             enddo
           enddo
         enddo
       enddo

     endif
   enddo
   zpr(:,:,:) = half*quarter*piinv*sigma/nxi * zpr(:,:,:)
   self%zpr_k(:,:,:) = zpr(:,:,:)

   ! FD: actual finite differences to obtain the second derivative of ZPR
   ! (Eq. (80) of [Guster2021])
   ABI_MALLOC(zpr_ddk, (self%ndeg, nkdir))

   zpr_ddk(:,:) = zero
   do ikdir=1,nkdir
     do iband=1,self%ndeg
       ! Copied from the previous implementation
       ! Does not look like a usual FD expression...
       zpr_ddk(iband, ikdir) = four/three *(zpr(2, iband, ikdir) - &
         zpr(1, iband, ikdir) - (zpr(3, iband, ikdir) - &
         zpr(1, iband, ikdir)) / 16.0_dp) * two/deltak**2
     enddo
   enddo

   do ikdir=1,nkdir
     do iband=1,self%ndeg
       self%invpolmas(iband, ikdir) = &
         two/deltak**2 * lutt_eigenval(iband, 2, ikdir) + zpr_ddk(iband, ikdir)
     enddo
   enddo

   ABI_FREE(eigenvec)
   ABI_FREE(eigenval)
   ABI_FREE(work)
   ABI_FREE(invefmas)
   ABI_FREE(lutt_eigenvec)
   ABI_FREE(lutt_eigenval)
   ABI_FREE(phfreq)
   ABI_FREE(investar)
   ABI_FREE(intsum)
   ABI_FREE(zpr)
   ABI_FREE(lk_h)
   ABI_FREE(zpr_ddk)
   ABI_FREE(intsuminv)
 endif

end subroutine frohlich_calc_polaronmass
!!***


!!****f* m_frohlich/frohlich_calc_zpr
!! NAME
!!  frohlich_calc_zpr
!!
!! FUNCTION
!!  Description
!!
!! INPUTS
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! OUTPUT
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! SOURCE

subroutine frohlich_calc_zpr(self)

!Arguments ------------------------------------
!scalars
 class(frohlich_t), intent(inout) :: self

!Local variables-------------------------------
!scalars
 integer :: iqdir, nu
 integer :: iband, jband
 integer :: lwork, info
!arrays
 real(dp), allocatable :: efmas_qdir(:,:)
 real(dp), allocatable :: invefmas_avg(:)
 real(dp), allocatable :: eigenval(:), rwork(:)
 complex(dp), allocatable :: eigenvec(:,:), work(:)
 complex(dp), allocatable :: f3d(:,:)
 logical, allocatable :: efmas_pos(:)
! *************************************************************************

 ! TODO: check if the electronic and phonon parts are initialized

 ! Initialization of the diagonalization routine for degenerate case
 ABI_MALLOC(eigenval, (self%ndeg))
 if (self%ndeg > 1) then
   ABI_MALLOC(eigenvec, (self%ndeg, self%ndeg))
   lwork = -1
   ABI_MALLOC(rwork, (3*self%ndeg - 2))
   ABI_MALLOC(work, (1))
   call zheev('V', 'U', self%ndeg, eigenvec(:,:), self%ndeg, eigenval(:), &
     work(:), lwork, rwork(:), info)
   lwork = int(work(1))
   ABI_FREE(work)
   ABI_MALLOC(work, (lwork))
 endif

 ABI_MALLOC(f3d, (self%ndeg, self%ndeg))
 ABI_MALLOC(invefmas_avg, (self%ndeg))
 ABI_MALLOC(efmas_pos, (self%ndeg))
 ABI_MALLOC(efmas_qdir, (self%ndeg, self%nqdir))
 ABI_MALLOC(self%sqrt_efmas_avg, (self%ndeg))
 ABI_MALLOC(self%saddle_warn, (self%ndeg))
 ABI_MALLOC(self%zpr_band, (self%ndeg))

 ! Effective mass tensor in q-vector directions
 do iqdir=1,self%nqdir

   ! Band curvature tensor (inverse effective mass)
   do iband=1,self%ndeg
     do jband=1,self%ndeg
       f3d(iband, jband) = dot_product(self%unit_qdir(:, iqdir), matmul( &
         self%eig2_diag_cart(:,:, iband, jband), self%unit_qdir(:, iqdir)))
     enddo
   enddo

   ! Band curvature tensor diagonalization
   if (self%ndeg == 1) then
     eigenval(1) = f3d(1, 1)
   else
     eigenvec(:,:) = f3d(:,:)
     eigenval(:) = zero
     work(:) = zero
     rwork(:) = zero
     call zheev('V', 'U', self%ndeg, eigenvec(:,:), self%ndeg, eigenval(:), &
       work(:), lwork, rwork(:), info)
     ABI_CHECK(info == 0, sjoin("zheev returned info: ", itoa(info)))
   endif

   efmas_qdir(:, iqdir) = one/eigenval(:)
 enddo

 ! Integration over the q-sphere: square root of the effective mass average and
 ! ZPR for each band (Eq. (17) of [deMelo2023])
 self%saddle_warn(:) = .false.
 self%sqrt_efmas_avg(:) = zero
 self%zpr_band(:) = zero
 invefmas_avg(:) = zero
 do iqdir=1,self%nqdir
   ! square root of the effective mass average
   self%sqrt_efmas_avg(:) = self%sqrt_efmas_avg(:) + &
     self%weights_qdir(iqdir)*sqrt(abs(efmas_qdir(:, iqdir)))

   ! ZPR
   do nu=4,3*self%natom
     self%zpr_band(:) = self%zpr_band(:) + (self%weights_qdir(iqdir) * &
       self%investar(nu, iqdir)*sqrt(self%phfreq_qdir(nu, iqdir) * &
       abs(efmas_qdir(:, iqdir))))
   enddo

   ! Inverse effective mass average: used to obtain the sign of the ZPR
   invefmas_avg(:) = invefmas_avg(:) + &
     self%weights_qdir(iqdir)/efmas_qdir(:, iqdir)

   ! Check for saddle-points
   if (iqdir == 1) efmas_pos(:) = (efmas_qdir(:, iqdir) > 0)
   do iband=1,self%ndeg
     if (efmas_pos(iband) .neqv. (efmas_qdir(iband, iqdir) > 0)) then
       self%saddle_warn(iband) = .true.
     endif
   enddo
 enddo
 self%sqrt_efmas_avg(:) = quarter*piinv*self%sqrt_efmas_avg(:)
 self%zpr_band(:) = sqrthalf*quarter*piinv*self%zpr_band(:)
 invefmas_avg(:) = quarter*piinv*invefmas_avg(:)

 ! Check for the sign problem and caclulate total ZPR and polaron formation
 ! energy in the weak-coupling regime
 self%sign_warn = .false.
 do iband=1,self%ndeg
   if (self%saddle_warn(iband) .or. &
     (efmas_pos(iband) .neqv. efmas_pos(1))) then
     self%sign_warn = .true.
   else
     self%sqrt_efmas_avg(iband) = &
       sign(self%sqrt_efmas_avg(iband), invefmas_avg(iband))

     self%zpr_band(iband) = sign(self%zpr_band(iband), invefmas_avg(iband))
   endif
 enddo

 if (.not. self%sign_warn) then
   self%zpr = -sum(self%zpr_band(:)) / self%ndeg
   self%sqrt_efmas_tot = sum(self%sqrt_efmas_avg(:)) / self%ndeg
   self%enpol_wc = -self%sqrt_efmas_tot*sum(self%dielavg(:))
 endif

 ABI_FREE(eigenval)
 ABI_FREE(f3d)
 ABI_FREE(invefmas_avg)
 ABI_FREE(efmas_pos)
 ABI_FREE(efmas_qdir)
 if (self%ndeg > 1) then
   ABI_FREE(eigenvec)
   ABI_FREE(rwork)
   ABI_FREE(work)
 endif

end subroutine frohlich_calc_zpr
!!***


!!****f* m_frohlich/frohlich_init_el
!! NAME
!!  frohlich_init_el
!!
!! FUNCTION
!!  Description
!!
!! INPUTS
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!  cryst<crystal_t> = Structure defining the unit cell
!!    %rprimd(3, 3) = real space primitvie vectors
!!  kpt(3) = k-point characterizing the electronic subspace
!!  ndeg = Number of degenerate bands
!!  eig2_diag(3, 3, ndeg, ndeg) = Band curvature double tensor
!!
!! OUTPUT
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! SOURCE

subroutine frohlich_init_el(self, cryst, kpt, ndeg, eig2_diag)

!Arguments ------------------------------------
!scalars
 class(frohlich_t), intent(inout) :: self
 type(crystal_t), intent(in) :: cryst
 integer,intent(in) :: ndeg
!arrays
 real(dp), intent(in) :: kpt(3)
 complex(dp), intent(in) :: eig2_diag(3, 3, ndeg, ndeg)

!Local variables-------------------------------
!scalars
 integer :: lwork, info
 integer :: iband, jband
 integer :: idir, ipar
!arrays
 real(dp) :: unit_kdir(3, 3)
 real(dp) :: eigenval(ndeg), lutt_eigenval(ndeg, ndeg)
 real(dp), allocatable :: rwork(:)
 complex(dp) :: eigenvec(ndeg, ndeg), lutt_eigenvec(ndeg, ndeg)
 complex(dp), allocatable :: work(:)
 logical :: lutt_found(3)
! *************************************************************************

 self%kpt(:) = kpt(:)
 self%ndeg = ndeg

 ! Initialize the band curvatutre double tensor in Cartesian coordiantes
 ABI_MALLOC(self%eig2_diag_cart, (3, 3, ndeg, ndeg))

 do iband=1,ndeg
   do jband=1,ndeg
     self%eig2_diag_cart(:,:, iband, jband) = one/two_pi**2 * &
       matmul(matmul(cryst%rprimd(:,:),eig2_diag(:,:, iband, jband)), &
       transpose(cryst%rprimd(:,:)))
   enddo
 enddo

 ! Determine the type of the Fr\"ohlich model
 ! TODO: check if a material actually has cubic symmetry for kind = 3?
 if (ndeg == 1) then
   ! Standard or Anisotropic Fr\"ohlich model
   self%kind = 1

   ! Assuming the single band is parabolic in the three Cartesian directions
   unit_kdir(:, 1) = (/1, 0, 0/)
   unit_kdir(:, 2) = (/0, 1, 0/)
   unit_kdir(:, 3) = (/0, 0, 1/)

   do idir=1,3
     self%band_params(idir) = dot_product(unit_kdir(:, idir), &
       matmul(self%eig2_diag_cart(:,:, 1, 1), unit_kdir(:, idir)))
   enddo

 else if (ndeg == 3) then
   ! Cubic generalized Fr\"ohlich model with triply degenerate bands
   self%kind = 2

   ! Symmety inequivalent cubic directions to obtain the Luttinger parameters
   unit_kdir(:, 1) = (/1, 0, 0/)
   unit_kdir(:, 2) = (/1, 1, 0/)/sqrt(2.0)
   unit_kdir(:, 3) = (/1, 1, 1/)/sqrt(3.0)

   ! Initialize the diagonalization routine
   lwork = -1
   ABI_MALLOC(rwork, (3*ndeg - 2))
   ABI_MALLOC(work, (1))
   call zheev('V', 'U', ndeg, eigenvec(:,:), ndeg, eigenval, work(:), lwork, &
     rwork(:), info)
   lwork=int(work(1))
   ABI_FREE(work)
   ABI_MALLOC(work, (lwork))

   lutt_eigenval(:, :) = zero
   ! Inverse effective mass tensor in the symmetry inequivalent directions
   do idir=1,3

     do iband=1,ndeg
       do jband=1,ndeg
         lutt_eigenvec(iband, jband) = dot_product(unit_kdir(:, idir), &
           matmul(self%eig2_diag_cart(:, :, iband, jband), unit_kdir(:, idir)))
       enddo
     enddo

     work(:) = zero
     rwork(:) = zero
     call zheev('V', 'U', ndeg, lutt_eigenvec(:,:), ndeg, &
       lutt_eigenval(idir, :), work(:), lwork, rwork(:), info)
     ABI_CHECK(info == 0, sjoin("zheev returned info: ", itoa(info)))
  enddo

  ABI_FREE(work)
  ABI_FREE(rwork)

  ! TODO: this looks ugly
  ! is there any better wa to analyze and set the Luttinger parameters?

  ! Check degeneracies in the (100) direction, get A and B luttinger params
  ! Inverse effective masses are: 2*A and 2*B (twofold)
  if (abs(lutt_eigenval(1, 2) - lutt_eigenval(1, 3)) < tol5) then
    self%band_params(1)=half*lutt_eigenval(1, 1)
    self%band_params(2)=half*half*(lutt_eigenval(1, 2) + lutt_eigenval(1, 3))
  else if (abs(lutt_eigenval(1, 2) - lutt_eigenval(1, 1)) < tol5) then
    self%band_params(1) = half*half*(lutt_eigenval(1, 2) + lutt_eigenval(1, 1))
    self%band_params(2) = half*lutt_eigenval(1, 3)
  else
    self%lutt_warn(1) = .true.
  endif

  ! Check degeneracies in the (111) direction, get C luttinger parameter
  ! Inverse effective masses are 2/3*(A + 2B - C) (twofold) and
  ! 2/3*(A + 2B + 2C)
  if (abs(lutt_eigenval(3, 2) - lutt_eigenval(3, 3)) < tol5) then
    self%band_params(3) = self%band_params(1) + 2*self%band_params(2) - &
      onehalf*half*(lutt_eigenval(3, 2) + lutt_eigenval(3, 3))
  else if (abs(lutt_eigenval(3, 2) - lutt_eigenval(3, 1)) < tol5) then
    self%band_params(3) = self%band_params(1) + 2*self%band_params(2) - &
      onehalf*half*(lutt_eigenval(3, 2) + lutt_eigenval(3, 1))
  else
    self%lutt_warn(2) = .true.
  endif

  ! Verifty that the (110) direction inverse effective masses are coherent
  ! with the Luttinger parameters: 2*B, A + B + C, A + B - C
  lutt_found(:) = .false.
  do ipar=1,ndeg
    if (abs(lutt_eigenval(2, ipar) - 2*self%band_params(2)) < tol4) then
      lutt_found(1) = .true.
    else if (abs(lutt_eigenval(2, ipar) - (self%band_params(1) + &
      self%band_params(2) - self%band_params(3))) < tol4) then
      lutt_found(2) = .true.
    else if (abs(lutt_eigenval(2, ipar) - (self%band_params(1) + &
      self%band_params(2) + self%band_params(3))) < tol4) then
      lutt_found(3) = .true.
    endif
  enddo

  if (.not. (all(lutt_found))) self%lutt_warn(3) = .true.

 else
   ! Generalized Fr\"ohlich model
   self%kind = 0
 endif

 ! The electronic subspace has been initialized
 self%isinitel = .true.

end subroutine frohlich_init_el
!!***


!!****f* m_frohlich/frohlich_init_ph
!! NAME
!!  frohlich_init_ph
!!
!! FUNCTION
!!  Description
!!
!! INPUTS
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!  cryst<crystal_t> = Structure defining the unit cell
!!    %natom = number of atoms in the unit cell
!!    %ucvol = real space unit cell volume
!!  efmas_ntheta = number of points required for spherical integration used
!!  to obtain the effective mass tensor [Laflamme2016]
!!  ifc<ifc_type> = Dynamical matrix and interatomic force constants
!!
!! OUTPUT
!!  self<frohlich_t> = Datatype gathering information on the Fr\"ohlich model
!!
!! SOURCE

subroutine frohlich_init_ph(self, cryst, efmas_ntheta, ifc)

!Arguments ------------------------------------
!scalars
 class(frohlich_t), intent(inout) :: self
 type(crystal_t), intent(in) :: cryst
 type(ifc_type), intent(in) :: ifc
 integer, intent(in) :: efmas_ntheta
!arrays

!Local variables-------------------------------
!scalars
 integer :: ntheta, nphi, nqdir
 integer :: iphi, itheta, iqdir, nu
! integer :: ikpt, ideg, iband, jband
! integer :: ndeg
 real(dp) :: weight
 real(dp) :: weight_phi, phi_radians
 real(dp) :: costheta, sintheta, cosphi, sinphi
!arrays
! real(dp) :: kpt(3)
 real(dp), allocatable :: gq_points_theta(:), gq_weights_theta(:)
 real(dp), allocatable :: gq_points_cosphi(:), gq_points_sinphi(:)

! *************************************************************************

 self%natom = cryst%natom
 self%ucvol = cryst%ucvol
 self%gmet(:,:) = cryst%gmet(:,:)

 ! Initialization of integrals
 ! This part allocates and initializes arrays unit_qdir(3,nqdir) and
 ! weights_qdir(nqdir) used for spherical integration over q-points directions
 ntheta = efmas_ntheta
 nphi = 2*ntheta
 nqdir = nphi*ntheta

 ABI_MALLOC(gq_points_theta, (ntheta))
 ABI_MALLOC(gq_weights_theta, (ntheta))
 ABI_MALLOC(gq_points_cosphi, (nphi))
 ABI_MALLOC(gq_points_sinphi, (nphi))

 self%nqdir = nqdir
 ABI_MALLOC(self%unit_qdir, (3, nqdir))
 ABI_MALLOC(self%weights_qdir, (nqdir))

 call cgqf(ntheta, 1, zero, zero, zero, pi, gq_points_theta, gq_weights_theta)
 weight_phi = two*pi/real(nphi, dp)

 do iphi=1,nphi
   phi_radians = weight_phi * (iphi-1)
   gq_points_cosphi(iphi) = cos(phi_radians)
   gq_points_sinphi(iphi) = sin(phi_radians)
 enddo

 nqdir = 0
 do itheta=1,ntheta
   costheta = cos(gq_points_theta(itheta))
   sintheta = sin(gq_points_theta(itheta))
   weight = gq_weights_theta(itheta)*weight_phi*sintheta

   do iphi=1,nphi
     cosphi = gq_points_cosphi(iphi)
     sinphi = gq_points_sinphi(iphi)
     nqdir = nqdir + 1

     self%unit_qdir(1, nqdir) = sintheta*cosphi
     self%unit_qdir(2, nqdir) = sintheta*sinphi
     self%unit_qdir(3, nqdir) = costheta
     self%weights_qdir(nqdir) = weight
   enddo
 enddo

 ABI_FREE(gq_points_theta)
 ABI_FREE(gq_weights_theta)
 ABI_FREE(gq_points_cosphi)
 ABI_FREE(gq_points_sinphi)

 ! Initialization of the generalized Fr\"ohlich parameters

 ABI_MALLOC(self%dielt_qdir, (nqdir))
 ABI_MALLOC(self%phfreq_qdir, (3*self%natom, nqdir))
 ABI_MALLOC(self%polarity_qdir, (3, 3*self%natom, nqdir))
 ABI_MALLOC(self%proj_polarity_qdir, (3*self%natom, nqdir))
 ABI_MALLOC(self%investar, (3*self%natom, nqdir))
 ABI_MALLOC(self%dielavg, (3*self%natom))
 ABI_MALLOC(self%isiractive, (3*self%natom))

 ! Phonon frequencies and mode polarity vectors for each q-vector direction
 call ifc%calcnwrite_nana_terms(cryst, nqdir, self%unit_qdir(:,:), &
   phfrq2l=self%phfreq_qdir, polarity2l=self%polarity_qdir)

 ! High-frequency dielectric constant for each q-vector direction
 do iqdir=1,nqdir
   self%dielt_qdir(iqdir) = dot_product(self%unit_qdir(:, iqdir), &
     matmul(ifc%dielt(:,:), self%unit_qdir(:, iqdir)))
 enddo

 self%investar(:,:) = zero
 self%dielavg(:) = zero
 ! Prjections of mode polarities on unit q-vectors,
 ! inverse effective dielectric constants for each mode and q-vector directions
 ! and dielectric average over q-vectors (Eqs. (22), (26) of [deMelo2023])
 do iqdir=1,nqdir
   do nu=4,3*self%natom
     self%proj_polarity_qdir(nu, iqdir) = &
       dot_product(self%unit_qdir(:, iqdir), self%polarity_qdir(:, nu, iqdir))

     self%investar(nu, iqdir) = (self%proj_polarity_qdir(nu, iqdir) / &
       self%dielt_qdir(iqdir) / self%phfreq_qdir(nu, iqdir))**2 * &
       four*pi/self%ucvol

     self%dielavg(nu) = self%dielavg(nu) + (self%weights_qdir(iqdir) * &
       self%investar(nu, iqdir)*sqrt(self%phfreq_qdir(nu, iqdir)))
   enddo
 enddo
 self%dielavg(:) = sqrthalf*quarter*piinv * self%dielavg(:)

 self%isiractive(:) = .false.
 ! Check for the infrared-active phonon modes: dielectric average is non-zero
 do nu=4,3*self%natom
   if (self%dielavg(nu) > tol10) self%isiractive(nu) = .true.
 enddo

 ! Rough correction for the ZPR in the infrared limit of the Fr\"ohlich model
 self%zpr_gamma = zero
 do nu=4,3*self%natom
   do iqdir=1,self%nqdir
     self%zpr_gamma = self%zpr_gamma + &
       self%investar(nu, iqdir)*self%weights_qdir(iqdir)
   enddo
 enddo
!Do not remove this useless line: a bug in the ifx2025 compiler is avoided thanks to it ...
 write(std_out,'(a,es16.6)')' frohlich_init_ph : after the loop, self%zpr_gamma=',self%zpr_gamma
 self%zpr_gamma = quarter*piinv*self%zpr_gamma
 self%zpr_gamma = two*(three*quarter*piinv/self%ucvol)**third * self%zpr_gamma

 ! The phonon subspace has been initialized
 self%isinitph = .true.

end subroutine frohlich_init_ph
!!***


!!****f* m_frohlich/mat3inv
!! NAME
!!  mat3inv
!!
!! FUNCTION
!!  Inverts a 3x3 matrix of real elements
!!
!! INPUTS
!!  mm(3, 3) = real matrix to be inverted
!!
!! OUTPUT
!!  mit(3, 3) = inverse of the input matrix
!!
!! SOURCE

subroutine mat3inv(mm, mit)

!Arguments ------------------------------------
!arrays
 real(dp), intent(in) :: mm(3, 3)
 real(dp), intent(out) :: mit(3, 3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd
 character(len=500) :: msg
!arrays
 real(dp) :: tt(3, 3)

! *************************************************************************

 ! Minors and determinant
 tt(1, 1) = mm(2, 2)*mm(3, 3) - mm(3, 2)*mm(2, 3)
 tt(1, 2) = mm(1, 3)*mm(3, 2) - mm(1, 2)*mm(3, 3)
 tt(1, 3) = mm(1, 2)*mm(2, 3) - mm(1, 3)*mm(2, 2)
 tt(2, 1) = mm(2, 3)*mm(3, 1) - mm(2, 1)*mm(3, 3)
 tt(2, 2) = mm(1, 1)*mm(3, 3) - mm(3, 1)*mm(1, 3)
 tt(2, 3) = mm(1, 3)*mm(2, 1) - mm(1, 1)*mm(2, 3)
 tt(3, 1) = mm(2, 1)*mm(3, 2) - mm(2, 2)*mm(3, 1)
 tt(3, 2) = mm(1, 2)*mm(3, 1) - mm(1, 1)*mm(3, 2)
 tt(3, 3) = mm(1, 1)*mm(2, 2) - mm(2, 1)*mm(1, 2)

 dd = mm(1, 1)*tt(1, 1) + mm(2, 1)*tt(2, 1) + mm(3, 1)*tt(3, 1)

 ! Make sure the matrix is not singular
 if (dd /=0) then
   mit(:,:) = tt(:,:) / dd
 else
   write(msg, '(2a,2x,9(i0,1x),a)') 'Attempting to invert real array',ch10, &
     mm,' ==> determinant is zero.'
   ABI_ERROR(msg)
 endif

end subroutine mat3inv
!!***


!!****f* m_frohlich/lk_hamiltonian
!! NAME
!!  lk_hamiltonian
!!
!! FUNCTION
!!  Construct a Luttinger-Kohn Hamiltonian at a given k-point for a set of the
!! Luttinger-Kohn parameters A, B, C
!!
!! INPUTS
!!  band_params(3) = Luttinger-Kohn parameters A, B, C
!!  kpt(3) = wavevector at which the Hamiltonian is obtained (dimensional)
!!
!! OUTPUT
!!  h_lk(3, 3) = Luttingher-Kohn Hamiltonian, H(k)
!!
!! SOURCE

subroutine lk_hamiltonian(band_params, kpt, h_lk)

!Arguments ------------------------------------
!arrays
 real(dp), intent(in) :: band_params(3)
 real(dp), intent(in) :: kpt(3)
 real(dp), intent(out) :: h_lk(3, 3)

! *************************************************************************

 h_lk(1, 1) = band_params(1)*kpt(1)**2 + band_params(2)*(kpt(2)**2 + kpt(3)**2)
 h_lk(2, 2) = band_params(1)*kpt(2)**2 + band_params(2)*(kpt(1)**2 + kpt(3)**2)
 h_lk(3, 3) = band_params(1)*kpt(3)**2 + band_params(2)*(kpt(1)**2 + kpt(2)**2)
 h_lk(1, 2) = band_params(3)*kpt(1)*kpt(2)
 h_lk(1, 3) = band_params(3)*kpt(1)*kpt(3)
 h_lk(2, 3) = band_params(3)*kpt(2)*kpt(3)
 ! Symmetric matrix
 h_lk(2, 1) = h_lk(1, 2)
 h_lk(3, 1) = h_lk(1, 3)
 h_lk(3, 2) = h_lk(2, 3)

end subroutine lk_hamiltonian
!!***

end module m_frohlich
!!***

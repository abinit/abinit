!!****m* ABINIT/m_rta
!! NAME
!!  m_rta
!!
!! FUNCTION
!!  This module provides objects and procedures to compute transport properties by solving the
!!  linearized Boltzmann transport equation (BTE) in the relaxation-time approximation (RTA).
!!  or within the Iterative Bolztmann Equation (IBTE).
!!  The RTA has two different flavors: Self-energy Relaxation Time Approximation (SERTA)
!!  in which the back-scattering term is completely ignored and the Momentum-Relaxation Time Approximation (MRTA)
!!  in which backscattering is partly accounter for by multiplying the e-ph self-energy integrand function
!!  by the efficiency factor alpha that depends on the incoming/outgoing electron group velocity.
!!  The implementation assumes e-ph scattering although additional scattering mechanisms (e.g. ionized impurities)
!!  can be easily included once an appropriate model is added to the ab-initio e-ph scattering rates.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2026 ABINIT group (HM, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!!  Dimensional analysis for conductivity (sigma), mobility (mu).
!!
!!   [sigma] = Siemens/m  with S = Ampere/Volt = Ohm^-1
!!   [mu] = S L^2 Q
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_rta

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_copy
 use m_nctk
 use m_wfk
 use m_ephtk
 use m_sigmaph
 use m_dtset
 use m_dtfil
 use m_krank
 use netcdf

 use defs_datatypes,   only : pseudopotential_type
 use m_ebands,         only : ebands_t, gaps_t, edos_t, klinterp_t, klinterp_new
 use m_io_tools,       only : flush_unit, open_file
 use m_time,           only : cwtime, cwtime_report
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : bisect, simpson_int, safe_div, arth
 use m_fstrings,       only : strcat, sjoin, itoa, ltoa, stoa, ftoa, yesno
 use m_kpts,           only : kpts_timrev_from_kptopt, kpts_map
 use m_occ,            only : occ_fd, occ_dfde
 use m_pawtab,         only : pawtab_type
 use m_ddk,            only : ddkstore_t

 implicit none

 private
!!****

 public :: rta_driver      ! Compute transport properties within the RTA (SERTA and MRTA)
 public :: ibte_driver     ! Compute transport properties within the IBTE.
!!****

!----------------------------------------------------------------------

!!****t* m_rta/rta_t
!! NAME
!! rta_t
!!
!! FUNCTION
!! Container for transport quantities in the RTA
!!
!! SOURCE

type,public :: rta_t

   integer :: nsppol
   ! Number of independent spin polarizations.

   integer :: nspinor
   ! Number of spinor components.

   integer :: nkcalc
   ! Number of computed k-points i.e. k-points inside the sigma_erange energy window.

   integer :: ntemp
   ! Number of temperatures.

   integer :: nw
   ! Number of energies (chemical potentials) at which transport quantities are computed
   ! Same number of energies used in DOS.

   integer :: bmin, bmax, bsize
   ! Only bands between bmin and bmax are considered in the integrals
   ! as we don't compute linewidths for all states.
   ! bmin = minval(%bstart_ks); bmax = maxval(%bstop_ks)
   ! bisze = bmax - bmin + 1

   integer :: nrta
   ! Number of relaxation-time approximations used (1 for SERTA, 2 for MRTA)

   real(dp) :: eph_extrael
   ! Extra electrons per unit cell used to compute SERTA lifetimes in sigmaph.

   real(dp) :: eph_fermie
   ! Fermi level specified in the input file when computing the SIGEPH file.

   real(dp) :: transport_extrael
   ! Extra electrons per unit cell specified in the input file when computing the SIGEPH file.

   real(dp) :: transport_fermie
   ! Fermi level specified in the input file when computing the SIGEPH file.

   logical :: assume_gap
   ! True if we are dealing with a semiconductor.
   ! This parameter is initialized from the value of sigma_erange provided by the user.

   integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc, nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all degenerate states should be included when symmetries are used.

   integer,allocatable :: bstop_ks(:,:)
   ! bstop_ks(nkcalc, nsppol)

   integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc, nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all degenerate states should be included when symmetries are used.

  integer,allocatable :: kcalc2ibz(:,:)
   !kcalc2ibz(nkcalc, 6))
   ! Mapping ikcalc --> IBZ as reported by listkk.

  integer,allocatable :: kcalc2ebands(:,:)
   ! (6, nkcalc)
   ! Mapping ikcalc --> ebands IBZ
   ! Note that this array is not necessarily equal to kcalc2ibz computed in sigmaph
   ! because we may have used sigma_ngkpt to downsample the initial nkpt mesh.
   ! This array is computed in get_ebands and is equal to kcalc2ibz if sigma_nkpt == ngkpt

   real(dp),allocatable :: kTmesh(:)
   ! (%ntemp)
   ! List of k * T temperatures at which to compute the transport

   real(dp),allocatable :: eph_mu_e(:)
   ! (%ntemp)
   ! Chemical potential at this carrier concentration and temperature from sigeph (lifetime)

   real(dp),allocatable :: transport_mu_e(:)
   ! (%ntemp)
   ! Chemical potential at this carrier concentration and temperature

   real(dp),allocatable :: eminmax_spin(:,:)
   ! (2, %nsppol))
   ! min/Max energy of the original ebands object

   real(dp),allocatable :: linewidths(:,:,:,:,:)
   ! (ntemp, bmin:bmax, nkpt, nsppol, nrta)
   ! Linewidth in the IBZ computed in the SERTA/MRTA.
   ! Non-zero only for the kcalc k-points.

   real(dp),allocatable :: vbks(:,:,:,:)
   ! (3, bmin:bmax, nkpt, nsppol))
   ! band velocity in Cartesian coordinates in the IBZ
   ! Non-zero only for the kcalc k-points.

   type(gaps_t) :: gaps
   ! gaps of original ebands object. Only if assume_gap

   type(ebands_t) :: ebands
   ! bandstructure object used to compute the transport properties
   ! Allocate using only the relevant bands for transport
   ! including valence states to allow to compute different doping

   type(edos_t) :: edos
   ! electronic density of states
   ! edos%mesh is the mesh used for vv_dos, vvtau_dos and tau_dos
   ! (%nw)

   real(dp),allocatable :: tau_dos(:,:,:,:)
   ! tau(e) (isotropic average for tau_nk for SERTA and MRTA.
   ! (nw, ntemp, nsppol, nrta)

   real(dp),allocatable :: vv_dos(:,:,:,:)
   ! (v x v)  DOS
   ! (nw, 3, 3, nsppol)

   real(dp),allocatable :: vvtau_dos(:,:,:,:,:,:)
   ! (v x v * tau) DOS
   ! (nw, 3, 3, ntemp, nsppol, nrta)

   real(dp),allocatable :: n_ehst(:,:,:)
    ! (2, %nsppol, %ntemp)
    ! Number of electrons (e) and holes (h) per unit cell
    ! The first dimension is for electrons/holes.
    ! If nsppol == 2, the second dimension is the number of e/h for spin else the total number of e/h summed over spins.

   real(dp),allocatable :: l0(:,:,:,:,:,:), l1(:,:,:,:,:,:), l2(:,:,:,:,:,:)
   ! (3, 3, nw, nsppol, ntemp, nrta)
   ! Onsager coefficients in Cartesian coordinates

   real(dp),allocatable :: l11_mu(:,:,:,:,:), l12_mu(:,:,:,:,:), l22_mu(:,:,:,:,:)
   ! (3, 3, nsppol, ntemp, nrta)
   ! Onsager coefficients in Cartesian coordinates, calculated at the correct mu (at the exact temperature)

   real(dp),allocatable :: sigma(:,:,:,:,:,:)
   real(dp),allocatable :: seebeck(:,:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:,:)
   real(dp),allocatable :: pi(:,:,:,:,:,:)
   real(dp),allocatable :: zte(:,:,:,:,:,:)
   ! (3, 3, nw, nsppol, ntemp, nrta)
   ! Transport coefficients in Cartesian coordinates

   real(dp),allocatable :: mobility(:,:,:,:,:,:,:)
   ! Mobility
   ! (3, 3, nw, %ntemp, 2, %nsppol, nrta)
   ! 5-th index is for e-h

   real(dp),allocatable :: conductivity(:,:,:,:,:)
   ! Conductivity at the Fermi level
   ! (3, 3, %ntemp, %nsppol, nrta)

   real(dp),allocatable :: resistivity(:,:,:,:)
   ! (3, 3, %ntemp, nrta)

   real(dp),allocatable :: n(:,:,:)
   ! (nw, ntemp, 2) carrier density for e/h (n/cm^3)

   real(dp),allocatable :: mobility_mu(:,:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp, nrta)
   ! mobility for electrons and holes (third dimension) at transport_mu_e(ntemp)
   ! Third dimension is for electron/hole

   real(dp),allocatable :: conductivity_mu(:,:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp, nrta)
   ! Conductivity in Siemens * cm-1
   ! computed by summing over k-points rather that by performing an energy integration).

   real(dp),allocatable :: resistivity_mu(:,:,:,:)
   ! (3, 3, ntemp, nrta)
   ! Resistivity obtained by inverting conductivity after summing of e- and holes and spins
   ! computed by summing over k-points rather that by performing an energy integration).



 contains

   procedure :: compute_rta
   procedure :: compute_rta_mobility
   procedure :: print_rta_txt_files
   procedure :: write_tensor
   procedure :: free => rta_free
   procedure :: rta_ncwrite

 end type rta_t
!!***

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_driver
!! NAME
!! rta_driver
!!
!! FUNCTION
!! Driver to compute transport properties within the RTA.
!!
!! INPUTS
!! dtfil<datafiles_type>=variables related to files.
!! ngfftc(18)=Coarse FFT mesh.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! cryst<crystal_t>=Crystalline structure
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine rta_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfftc(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
 type(rta_t) :: rta
!arrays
 integer :: units(2)

! *************************************************************************

 units = [std_out, ab_out]
 call wrtout(units, ch10//' Entering transport RTA computation driver.')
 call wrtout(units, sjoin("- Reading carrier lifetimes from:", dtfil%filsigephin), newlines=1, do_flush=.True.)

 ! Initialize RTA object
 rta = rta_new(dtset, dtfil, ngfftc, cryst, ebands, pawtab, psps, comm)

 ! Compute RTA transport quantities
 call rta%compute_rta(cryst, dtset, dtfil, comm)

 call rta%free()

end subroutine rta_driver
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_new
!! NAME
!! rta_new
!!
!! FUNCTION
!! Build object to compute RTA transport quantities.
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  dtfil<datafiles_type>=variables related to files.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  comm=MPI communicator.
!!
!! SOURCE

type(rta_t) function rta_new(dtset, dtfil, ngfftc, cryst, ebands, pawtab, psps, comm) result (new)

!Arguments -------------------------------------
 integer, intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfftc(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
 integer,parameter :: sppoldbl1 = 1, master = 0
 integer :: ierr, spin, nprocs, my_rank, ik_ibz, ib, irta, itemp, ndat, nsppol, idat, mband, ikpt
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname_dense
 type(ebands_t) :: tmp_ebands, ebands_dense
 type(klinterp_t) :: klinterp
 type(ddkstore_t) :: ds
 type(sigmaph_t) :: sigmaph
 type(krank_t) :: krank
!arrays
 integer :: kptrlatt(3,3), units(2), sigma_ngkpt(3)
 integer,allocatable :: indkk(:,:)
 real(dp) :: extrael_fermie(2), sigma_erange(2)
 real(dp),allocatable :: values_bksd(:,:,:,:), vals_bsd(:,:,:), tmp_array4(:,:,:,:), tmp_array5(:,:,:,:,:)
!************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 units = [std_out, ab_out]

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! Use sigma_erange to understand if we are dealing with a metal or a semiconductor.
 new%assume_gap = (.not. all(dtset%sigma_erange < zero) .or. dtset%gw_qprange /= 0)

 ! Read data from SIGEPH file.
 sigmaph = sigmaph_read(dtfil%filsigephin, dtset, xmpi_comm_self, msg, ierr, keep_open=.true., &
                        extrael_fermie=extrael_fermie, sigma_ngkpt=sigma_ngkpt, sigma_erange=sigma_erange)
 ABI_CHECK(ierr == 0, msg)

 !if (any(sigma_erange /= zero)) then
 !  ABI_CHECK(all(dtset%sigma_erange /= zero), "sigma_erange is required in input with a value compatible with SIGEPH.nc")
 !  ! Make sure that the two values are consistent
 !  ! Cannot switch to metallic case if SigmaPH file was produced assuming gapped-system
 !  if (.not. (all(dtset%sigma_erange < zero) .eqv. all(sigma_erange < zero))) then
 !    ABI_ERROR("The values of sigma_erange from input and SIGEPH are not compatible")
 !  end if
 !end if

 ! How many RTA approximations have we computed in sigmaph? (SERTA, MRTA ...?)
 new%nrta = 2; if (sigmaph%mrta == 0) new%nrta = 1

 ! Copy important arrays from sigmaph file.
 ! Allocate temperature arrays (use same values as the ones used in the SIGEPH calculation).
 new%ntemp = sigmaph%ntemp
 call alloc_copy(sigmaph%kTmesh, new%kTmesh)

 new%nkcalc = sigmaph%nkcalc
 call alloc_copy(sigmaph%bstart_ks, new%bstart_ks)
 call alloc_copy(sigmaph%bstop_ks, new%bstop_ks)
 call alloc_copy(sigmaph%nbcalc_ks, new%nbcalc_ks)
 call alloc_copy(sigmaph%kcalc2ibz, new%kcalc2ibz)

 new%bmin = minval(sigmaph%bstart_ks); new%bmax = maxval(sigmaph%bstop_ks)
 !new%bmin = 1; new%bmax = ebands%mband ! This for debugging purposes, results should not change
 new%bsize = new%bmax - new%bmin + 1

 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor
 nsppol = new%nsppol

 ABI_MALLOC(new%eminmax_spin, (2, nsppol))
 new%eminmax_spin = ebands%get_minmax("eig")

 if (new%assume_gap) then
   ! Get gaps
   new%gaps = ebands%get_gaps(ierr)
   if (ierr /= 0) then
     do spin=1, nsppol
       ABI_WARNING(trim(new%gaps%errmsg_spin(spin)))
       new%gaps%vb_max(spin) = ebands%fermie - 1 * eV_Ha
       new%gaps%cb_min(spin) = ebands%fermie + 1 * eV_Ha
     end do
     !ABI_ERROR("ebands_get_gaps returned non-zero exit status. See above warning messages...")
     ABI_WARNING("ebands_get_gaps returned non-zero exit status. See above warning messages...")
   end if
   if (my_rank == master) call new%gaps%print(units)
 end if

 ! =================================================
 ! Read lifetimes and new%ebands from SIGEPH.nc file
 ! =================================================
 ! After this point we have:
 !
 !      vbks(3, bmin:bmax, nkpt, nsppol)
 !      linewidths(self%ntemp, bmin:bmax, nkpt, nsppol, 2)
 !
 if (any(dtset%sigma_ngkpt /= 0)) then
   ! If integrals are computed with the sigma_ngkpt k-mesh, we need to downsample ebands.
   !call wrtout(units, sjoin(" SIGMAPH file used sigma_ngkpt:", ltoa(sigma_ngkpt)))
   call wrtout(units, sjoin(" Computing integrals with downsampled sigma_ngkpt:", ltoa(dtset%sigma_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)

   tmp_ebands = ebands%downsample(cryst, kptrlatt, dtset%sigma_nshiftk, dtset%sigma_shiftk)
   new%ebands = sigmaph%get_ebands(cryst, tmp_ebands, [new%bmin, new%bmax], &
                                   new%kcalc2ebands, new%linewidths, new%vbks, xmpi_comm_self)
   call tmp_ebands%free()
 else
   !call wrtout(units, sjoin(" Computing integrals with SIGEPH k-mesh:", ebands_kmesh2str(ebands))
   new%ebands = sigmaph%get_ebands(cryst, ebands, [new%bmin, new%bmax], &
                                   new%kcalc2ebands, new%linewidths, new%vbks, xmpi_comm_self)
   kptrlatt = new%ebands%kptrlatt
 end if

 !print *, "linewidth_serta", maxval(abs(new%linewidths(:,:,:,:,1)))
 !print *, "linewidth_mrta", maxval(abs(new%linewidths(:,:,:,:,2)))
 !print *, "max velocities", maxval(abs(new%vbks))

 if ( &
     dtset%useria == 888 .and. &
     (dtset%getwfkfine /= 0 .or. dtset%irdwfkfine /= 0 .or. dtset%getwfkfine_filepath /= ABI_NOFILE)) then

   ! In principle only getwfkfine_filepath is used here
   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)
   ABI_CHECK(nctk_try_fort_or_ncfile(wfk_fname_dense, msg) == 0, msg)

   call wrtout(units, " EPH double grid interpolation: will read energies from: "//trim(wfk_fname_dense), newlines=1)
   mband = new%ebands%mband

   !ebands_dense = wfk_read_ebands(wfk_fname_dense, comm)

   tmp_ebands = wfk_read_ebands(wfk_fname_dense, comm)
   ebands_dense = tmp_ebands%chop(1, mband)
   call tmp_ebands%free()
   if (my_rank == master) then
     write(std_out, *)" Using kptrlatt: ", ebands_dense%kptrlatt
     write(std_out, *)"       shiftk: ", ebands_dense%shiftk
   end if
   ABI_CHECK_IEQ(mband, ebands_dense%mband, "Inconsistent number of bands for the fine and dense grid:")

   ! Compute v_{nk} on the dense grid in Cartesian coordinates.
   ! vdiago(3, bmin:bmax, nkpt, nsppol)
   ! NB: We select bands in [bmin:bmax] but all k-points in the IBZ are computed!
   ds%only_diago = .True.; ds%bmin = new%bmin; ds%bmax = new%bmax; ds%mode = "cart"
   call ds%compute_ddk(wfk_fname_dense, "", dtset, psps, pawtab, ngfftc, comm)

   ! Transfer data to new%vbks
   ABI_MOVE_ALLOC(ds%vdiago, new%vbks)
   call ds%free()
   !print *, "vbks:", new%vbks

   ! Linear interpolation in k-space of the linewidths from input SIGEPH to the dense IBZ provided by fine WFK file.
   ! First of all transfer linewidths to values_bksd to prepare call to klinterp_new.
   ndat = new%ntemp * new%nrta
   ABI_MALLOC(values_bksd, (new%bmin:new%bmax, new%ebands%nkpt, nsppol, ndat))

   do irta=1,new%nrta
     do spin=1,nsppol
       do ik_ibz=1,new%ebands%nkpt
         do ib=new%bmin,new%bmax
           do itemp=1,new%ntemp
             idat = itemp + new%ntemp * (irta - 1)
             values_bksd(ib, ik_ibz, spin, idat) = new%linewidths(itemp, ib, ik_ibz, spin, irta) ! - t(e)
           end do
         end do
       end do
     end do
   end do

   ! Build linear interpolator for linewidths (use only bsize bands)
   klinterp = klinterp_new(cryst, new%ebands%kptrlatt, new%ebands%nshiftk, new%ebands%shiftk, new%ebands%kptopt, &
                           new%ebands%kptns, new%bsize, new%ebands%nkpt, nsppol, ndat, values_bksd, comm)
   ABI_FREE(values_bksd)

   ! HERE we re-malloc new%ebands and %linewidths on the fine k-mesh.
   ! The call must be executed here, once klinterp has been built.
   ! After this point we can use new%ebands to allocate stuff.
   call ebands_dense%move_alloc(new%ebands)

   ! Unlinke the ebands stored in SIGEPH, the eigens read from WFK_FINE have not been
   ! shifted with the scissors operator or updated according to extrael_fermie so do it now.
   call ephtk_update_ebands(dtset, new%ebands, "GS energies read from WFK_FINE")

   ! And now interpolate linewidths on the fine k-mesh
   ! Note: k-points that close to the edge of the pocket may get zero linewidths
   ! One may fix the problem by using lw(e).
   ABI_REMALLOC(new%linewidths, (new%ntemp, new%bmin:new%bmax, new%ebands%nkpt, nsppol, new%nrta))
   ABI_MALLOC(vals_bsd, (new%bmin:new%bmax, nsppol, ndat))

   ierr = 0
   do ik_ibz=1,new%ebands%nkpt

     call klinterp%eval_bsd(new%ebands%kptns(:, ik_ibz), vals_bsd)
     !vals_bsd = vals_bsd + lw(e)

     if (any(vals_bsd < zero)) then
       ierr = ierr + 1
       where (vals_bsd < zero) vals_bsd = zero
     end if

     ! Transfer data.
     do spin=1,nsppol
       do irta=1,new%nrta
         do itemp=1,new%ntemp
           idat = itemp + new%ntemp * (irta - 1)
           do ib=new%bmin,new%bmax
             new%linewidths(itemp, ib, ik_ibz, spin, irta) = vals_bsd(ib, spin, idat)
           end do
         end do
       end do
     end do

   end do ! ik_ibz

   if (ierr /= 0) then
     ! This should never happen for linear interpolation.
     ABI_WARNING(sjoin("Linear interpolation produced:", itoa(ierr), " k-points with negative linewidths"))
   end if

   ABI_FREE(vals_bsd)
   call klinterp%free()
 end if

 ! FIXME: I think transport_ngkpt is buggy, wrong ne(T), weird zeros if MRTA ...
 ! Do we really need this option? Can't we replace it with sigma_ngkpt and eph_task 7?

 if (any(dtset%transport_ngkpt /= 0)) then
   ! Perform further downsampling (useful for debugging purposes)
   call wrtout(units, " Downsampling the k-mesh before computing transport:")
   call wrtout(units, sjoin(" Using transport_ngkpt: ", ltoa(dtset%transport_ngkpt)))
   kptrlatt = 0
   kptrlatt(1, 1) = dtset%transport_ngkpt(1)
   kptrlatt(2, 2) = dtset%transport_ngkpt(2)
   kptrlatt(3, 3) = dtset%transport_ngkpt(3)
   tmp_ebands = new%ebands%downsample(cryst, kptrlatt, 1, [zero, zero, zero])

   ! Map the points of the downsampled bands to dense ebands
   ABI_MALLOC(indkk, (6, tmp_ebands%nkpt))

   call krank%from_kptrlatt(new%ebands%nkpt, new%ebands%kptns, new%ebands%kptrlatt, compute_invrank=.False.)

   if (kpts_map("symrec", ebands%kptopt, cryst, krank, tmp_ebands%nkpt, tmp_ebands%kptns, indkk) /= 0) then
     write(msg, '(3a)' ) &
       "Error while downsampling ebands in the transport driver",ch10, &
       "The k-point could not be generated from a symmetrical one."
     ABI_ERROR(msg)
   end if

   call krank%free()

   ! Downsampling linewidths and velocities.
   ABI_MOVE_ALLOC(new%linewidths, tmp_array5)
   ABI_REMALLOC(new%linewidths, (new%ntemp, new%bmin:new%bmax, tmp_ebands%nkpt, nsppol, new%nrta))
   do ikpt=1,tmp_ebands%nkpt
     new%linewidths(:,:,ikpt,:,:) = tmp_array5(:,:,indkk(1, ikpt),:,:)
   end do
   ABI_FREE(tmp_array5)

   ABI_MOVE_ALLOC(new%vbks, tmp_array4)
   ABI_REMALLOC(new%vbks, (3, new%bmin:new%bmax, tmp_ebands%nkpt, nsppol))
   do ikpt=1,tmp_ebands%nkpt
     new%vbks(:,:,ikpt,:) = tmp_array4(:,:,indkk(1, ikpt),:)
   end do
   ABI_FREE(tmp_array4)

   !print *, "after downsampling linewidths"
   !print *, "linewidth_serta", maxval(abs(new%linewidths(:,:,:,:,1)))
   !print *, "linewidth_mrta", maxval(abs(new%linewidths(:,:,:,:,2)))

   ABI_FREE(indkk)
   call tmp_ebands%move_alloc(new%ebands)
 end if

 ! Same doping case as in sigmaph file.
 ABI_MALLOC(new%eph_mu_e, (new%ntemp))
 ABI_MALLOC(new%transport_mu_e, (new%ntemp))

 new%eph_extrael = extrael_fermie(1)
 new%eph_fermie = extrael_fermie(2)
 new%transport_fermie = dtset%eph_fermie
 new%transport_extrael = dtset%eph_extrael
 new%eph_mu_e = sigmaph%mu_e
 new%transport_mu_e = sigmaph%mu_e

 if (new%transport_fermie /= zero) new%transport_mu_e = new%transport_fermie

 if (new%transport_fermie == zero .and. new%transport_extrael /= new%eph_extrael) then

   if (new%transport_extrael /= new%eph_extrael) then
     write(msg,'(2(a,e18.8),3a)') &
       ' extrael from SIGEPH: ',new%transport_extrael, ' and input file: ',new%eph_extrael, "differ", ch10, &
       ' Will recompute the chemical potential'
     call wrtout(std_out, msg)
   end if

   ! Compute Fermi level for different T values.
   call ebands%get_muT_with_fd(new%ntemp, new%kTmesh, dtset%spinmagntarget, dtset%prtvol, new%transport_mu_e, comm)
 end if

 ! TODO: Implement possible change of sigma_erange, useful for convergence studies
 !   1) Run sigmaph with relatively large sigma_erange.
 !   2) Decrease energy window in the transport part to analyze the behaviour of transport tensors.

 ! sigmaph is not needed anymore. Free it.
 sigmaph%ncid = nctk_noid
 call sigmaph%free()

 call cwtime_report(" rta_new", cpu, wall, gflops)

end function rta_new
!!***

!----------------------------------------------------------------------

!!****f* m_rta/compute_rta
!! NAME
!! compute_rta
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! dtfil<datafiles_type>=variables related to files.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine compute_rta(self, cryst, dtset, dtfil, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 class(rta_t),intent(inout) :: self
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst

!Local variables ------------------------------
 integer,parameter :: nvecs0 = 0, master = 0
 integer :: nsppol, nkibz, ib, ik_ibz, iw, spin, ii, jj, itemp, irta, itens_, iscal, cnt
 integer :: ntens, edos_intmeth, ifermi, iel, nvals, my_rank
 integer :: ncid
 !character(len=500) :: msg
 character(len=fnlen) :: path
 real(dp) :: emin, emax, edos_broad, edos_step, max_occ, kT, Tkelv, linewidth, fact0, cpu, wall, gflops
!arrays
 integer :: units(2)
 real(dp) :: vr(3), dummy_vecs(1,1,1,1,1), work_33(3,3), S_33(3,3), mat33(3,3)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:,:), out_valsdos(:,:,:,:), dummy_dosvecs(:,:,:,:,:)
 real(dp),allocatable :: out_tensdos(:,:,:,:,:,:), tau_vals(:,:,:,:,:), l0inv_33nw(:,:,:)
!************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm)
 units = [std_out, ab_out]

 ! Basic dimensions
 nsppol = self%ebands%nsppol; nkibz = self%ebands%nkpt

 ! Allocate v x v tensors with and without the lifetimes. Eq 8 of [[cite:Madsen2018]]
 ! The total number of tensorial entries is ntens and accounts for nrta
 ! Remember that we haven't computed all the k-points in the IBZ hence we can have zero linewidths
 ! or very small values when the states are at the band edge so we use safe_dif to avoid SIGFPE.
 ! Also, note how we store only the states in the energy window.

 nvals = self%ntemp * self%nrta
 ABI_CALLOC(tau_vals, (self%ntemp, self%nrta, self%bmin:self%bmax, nkibz, nsppol))

 ntens = (1 + self%ntemp) * self%nrta
 ABI_CALLOC(vv_tens, (3, 3, 1 + self%ntemp, self%nrta, self%bmin:self%bmax, nkibz, nsppol))

 cnt = 0
 do spin=1,nsppol
   do ik_ibz=1,nkibz
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     do ib=self%bmin,self%bmax

       vr(:) = self%vbks(:, ib, ik_ibz, spin)
       ! Store outer product (v_bks x v_bks) in vv_tens. This part does not depend on T and irta.
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, 1:self%nrta, ib, ik_ibz, spin) = vr(ii) * vr(jj)
         end do
       end do

       ! Multiply by the lifetime (SERTA and MRTA)
       do irta=1,self%nrta
         do itemp=1,self%ntemp
           linewidth = self%linewidths(itemp, ib, ik_ibz, spin, irta)
           call safe_div(vv_tens(:,:, 1, irta, ib, ik_ibz, spin), two * linewidth, zero, &
                         vv_tens(:,:, 1 + itemp, irta, ib, ik_ibz, spin))
           call safe_div(one, two * linewidth, zero, tau_vals(itemp, irta, ib, ik_ibz, spin))
         end do
       end do

     end do
   end do
 end do

 !call xmpi_sum(vv_tens, comm, ierr)
 !call xmpi_sum(tau_vals, comm, ierr)
 call cwtime_report(" compute_rta_loop1", cpu, wall, gflops)

 ! Compute DOS and VV_DOS and VV_TAU_DOS
 ! Define integration method and mesh step.
 edos_intmeth = 2; if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae
 if (edos_step == zero) edos_step = 0.001
 !if (edos_step == zero) edos_step = ten / Ha_meV
 edos_broad = dtset%tsmear

 ! Set default energy range for DOS
 ! If sigma_erange is set, get emin and emax from this variable
 ! MG: TODO This value should be read from SIGEPH
 ! Recheck metals
 if (self%assume_gap) then
   emin = huge(one); emax = -huge(one)
   do spin=1,self%ebands%nsppol
     if (dtset%sigma_erange(1) >= zero) emin = min(emin, self%gaps%vb_max(spin) + tol2 * eV_Ha - dtset%sigma_erange(1))
     if (dtset%sigma_erange(2) >= zero) emax = max(emax, self%gaps%cb_min(spin) - tol2 * eV_Ha + dtset%sigma_erange(2))
   end do
   ABI_CHECK(emin /=  huge(one), "Cannot initialize emin")
   ABI_CHECK(emax /= -huge(one), "Cannot initialize emax")
 else
   emin = minval(self%eminmax_spin(1, :)); emin = emin - tol1 * abs(emin)
   emax = maxval(self%eminmax_spin(2, :)); emax = emax + tol1 * abs(emax)
 end if

 ! Compute DOS, vv_dos and vvtau_DOS (v x v tau)
 !
 !    out_valsdos: (nw, 2, nvals, nsppol) array with DOS for scalar quantities if nvals > 0
 !    out_tensdos: (nw, 2, 3, 3, ntens,  nsppol) array with DOS weighted by tensorial terms if ntens > 0
 !
 !  Vectors and tensors are in Cartesian coordinates.
 !  Note how we compute the DOS only between [emin, emax] to save time and memory
 !  this implies that IDOS and edos%ifermi are ill-defined

 self%edos = self%ebands%get_edos_matrix_elements(cryst, self%bsize, &
                                             nvals, tau_vals, nvecs0, dummy_vecs, ntens, vv_tens, &
                                             edos_intmeth, edos_step, edos_broad, &
                                             out_valsdos, dummy_dosvecs, out_tensdos, comm, &
                                             brange=[self%bmin, self%bmax], erange=[emin, emax])

 if (my_rank == master) then
   call self%edos%print([std_out, ab_out], header="Computation of DOS, VV_DOS and VVTAU_DOS")
 end if

 call cwtime_report(" compute_rta_edos", cpu, wall, gflops)

 ! Unpack data stored in out_tensdos with shape (nw, 2, 3, 3, ntens, nsppol)
 self%nw = self%edos%nw
 ABI_MALLOC(self%tau_dos, (self%nw, self%ntemp, nsppol, self%nrta))
 ! TODO: Exchange dims?
 ABI_MALLOC(self%vv_dos, (self%nw, 3, 3, nsppol))
 ABI_MALLOC(self%vvtau_dos, (self%nw, 3, 3, self%ntemp, nsppol, self%nrta))

 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp+1

       itens_ = itemp + (irta - 1) * (self%ntemp + 1)
       if (itemp == 1) then
         self%vv_dos(:,:,:,spin) = out_tensdos(:, 1, :, :, itens_, spin)
       else
         self%vvtau_dos(:,:,:, itemp-1, spin, irta) = out_tensdos(:, 1, :, :, itens_, spin)
       end if

     end do
   end do
 end do

 ! Transfer data for tau(e)
 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp
       iscal = itemp + (irta - 1) * self%ntemp
       self%tau_dos(:, itemp, spin, irta) = out_valsdos(:, 1, iscal, spin)
     end do
   end do
 end do

 ! Free memory
 ABI_SFREE(out_tensdos)
 ABI_SFREE(tau_vals)
 ABI_SFREE(out_valsdos)
 ABI_SFREE(dummy_dosvecs)
 ABI_SFREE(vv_tens)

 ! Compute Onsager coefficients. Eq 9 of [[cite:Madsen2018]]
 ! See also Eqs 41, page 11 of https://arxiv.org/pdf/1402.6979.pdf
 !
 !      L^\alpha(\mu, T) = \int de \sigma(e, T) (e - mu)^\alpha (-df/de)
 !
 ! with \sigma(e, T) stored in vvtau_dos

 ABI_MALLOC(self%l0, (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_MALLOC(self%l1, (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_MALLOC(self%l2, (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))

 call onsager(0, self%l0)
 call onsager(1, self%l1)
 call onsager(2, self%l2)

 call cwtime_report(" compute_rta_onsanger", cpu, wall, gflops)

 ! Compute transport tensors, Eqs 12-15 of [[cite:Madsen2018]] and convert to SI units.
 ABI_CALLOC(self%sigma,   (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%seebeck, (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%kappa,   (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%pi,      (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%zte,     (3, 3, self%nw, self%nsppol, self%ntemp, self%nrta))

 ! Sigma = L0
 !TODO : missing maxocc??
 fact0 = (siemens_SI / Bohr_meter / cryst%ucvol)
 self%sigma = fact0 * self%l0

 ! Used to stored L0^-1
 ABI_MALLOC(l0inv_33nw, (3, 3, self%nw))

 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp

       TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one

       ! S = -1/T L0^-1 L1 = -1/T sigma L1
       do iw=1,self%nw
         call inv33(self%l0(:, :, iw, spin, itemp, irta), work_33)
         l0inv_33nw(:,:,iw) = work_33
         self%seebeck(:,:,iw,spin,itemp,irta) = - (volt_SI / TKelv) * matmul(work_33, self%l1(:,:,iw,spin,itemp,irta))
       end do

       ! kappa = 1/T [L2 - L1 L0^-1 L1]
       ! HM: Check why do we need minus sign here to get consistent results with Boltztrap!
       ! MG: Likely because of a different definition of kappa.
       do iw=1,self%nw
         work_33 = self%l1(:, :, iw, spin, itemp, irta)
         work_33 = self%l2(:, :, iw, spin, itemp, irta) - matmul(work_33, matmul(l0inv_33nw(:, :, iw), work_33))
         !self%kappa(:,:,iw,spin, itemp,spin,irta) = - (volt_SI**2 * fact0 / TKelv) * work_33
         self%kappa(:,:,iw,spin,itemp,irta) = + (volt_SI**2 * fact0 / TKelv) * work_33
       end do

       ! Peltier pi = -L1 L0^-1
       do iw=1,self%nw
         work_33 = self%l1(:, :, iw, spin, itemp, irta)
         self%pi(:,:,iw,spin,itemp,irta) = - volt_SI * matmul(work_33, l0inv_33nw(:, :, iw))
       end do

       ! ZT: S^T sigma S k^-1 T (tensor form with k=k_electronic only):
       do iw=1,self%nw
         S_33 = self%seebeck(:,:,iw,spin,itemp,irta)
         S_33 = matmul(matmul(transpose(S_33), self%sigma(:,:,iw,spin,itemp,irta)), S_33)
         call inv33(self%kappa(:,:,iw,spin,itemp,irta), work_33)
         self%zte(:,:,iw,spin,itemp,irta) = matmul(S_33, work_33) * TKelv
       end do

     end do
   end do
 end do

 ABI_FREE(l0inv_33nw)

 !Here are computed the transport tensors at exact mu (exact temperature and not based on a grid)

!  do irta=1,self%nrta
!   do spin=1,nsppol
!     do itemp=1,self%ntemp
!
!       TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one
!
!       ! S = -1/T L11^-1 L12 = -1/T sigma^-1 L12
!         call inv33(self%l11_mu(:, :, spin, itemp, irta), work_33)
!        ! l0inv_33nw(:,:) = work_33
!         self%seebeck_mu(:,:,spin,itemp,irta) = - (volt_SI / TKelv) * matmul(work_33, self%l12_mu(:,:,spin,itemp,irta))
!
!
!       ! kappa = 1/T [L22 - L21 L11^-1 L12] (in RTA L12=L21)
!         call inv33(self%l11_mu(:, :, spin, itemp, irta), work_33)
!         work_33 = matmul(work_33,self%l12_mu(:, :, spin, itemp, irta))
!         work_33 = self%l22_mu(:, :, spin, itemp, irta) - matmul(self%l12_mu(:, :, spin, itemp, irta), work_33)
!         self%kappa_mu(:,:,spin,itemp,irta) = + (volt_SI**2 * fact0 / TKelv) * work_33
!
!
!         ! Peltier pi = -L12 L11^-1 (in RTA and only in RTA, if TR sym. is not broken)
!
!         call inv33(self%l11_mu(:, :, spin, itemp, irta), work_33)
!         self%pi_mu(:,:,spin,itemp,irta) = - volt_SI * matmul(self%l12_mu(:, :, spin, itemp, irta), work_33)
!
!     end do ! itemp
!   end do !spin
! end do  !irta



 ! Compute the index of the Fermi level and handle possible out of range condition.
 ifermi = bisect(self%edos%mesh, self%ebands%fermie)
 if (ifermi == 0 .or. ifermi == self%nw) then
   ABI_ERROR("Bisection could not find the index of the Fermi level in edos%mesh!")
 end if

 max_occ = two / (self%nspinor * self%nsppol)

 ! Conductivity
 ABI_MALLOC(self%conductivity, (3, 3, self%ntemp, self%nsppol, self%nrta))
 do spin=1,self%nsppol
   do itemp=1,self%ntemp
     do irta=1,self%nrta
       do jj=1,3
         do ii=1,3
           self%conductivity(ii,jj,itemp,spin,irta) = self%sigma(ii, jj, ifermi, spin, itemp, irta) * 0.01 !m^-1 to cm^-1
         end do
       end do
     end do
   end do ! itemp
 end do ! spin

 ABI_MALLOC(self%resistivity, (3, 3, self%ntemp, self%nrta))
 do irta=1,self%nrta
  do itemp=1,self%ntemp
    work_33 = sum(self%conductivity(:,:,itemp,:,irta), dim=3)
    call inv33(work_33, mat33); mat33 = 1e+6_dp * mat33
    self%resistivity(:, :, itemp, irta) = mat33
  end do
 end do



 ! Mobility
 ABI_MALLOC(self%n, (self%nw, self%ntemp, 2))
 ABI_MALLOC(self%mobility, (3, 3, self%nw, self%ntemp, 2, self%nsppol, self%nrta))

 do spin=1,self%nsppol
   do itemp=1,self%ntemp
     ! Compute carrier density
     kT = self%kTmesh(itemp)

     ! MG TODO: I think that here we should use mu_e instead of ifermi.
     ! Compute carrier density of electrons (ifermi:self%nw)
     do iw=1,self%nw ! doping
       self%n(iw,itemp,1) = carriers(self%edos%mesh, self%edos%dos(:,spin) * max_occ, ifermi, self%nw, &
                                     kT, self%edos%mesh(iw)) / cryst%ucvol / Bohr_meter**3
     end do

     ! Compute carrier density of holes (1:ifermi)
     do iw=1,self%nw ! doping
       self%n(iw,itemp,2) = carriers(self%edos%mesh, self%edos%dos(:,spin) * max_occ, 1, ifermi, &
                                     kT, self%edos%mesh(iw)) / cryst%ucvol / Bohr_meter**3
     end do

     self%n(:,itemp,2) = self%n(self%nw,itemp,2) - self%n(:,itemp,2)

     ! Compute mobility
     do irta=1,self%nrta
       do iel=1,2
         do iw=1,self%nw
           do jj=1,3
             do ii=1,3
               call safe_div(self%sigma(ii, jj, iw, spin, itemp, irta) * 100**2, &
                             e_Cb * self%n(iw, itemp, iel), &
                             zero, self%mobility(ii, jj, iw, itemp, iel, spin, irta))
             end do
           end do
         end do
       end do
     end do
   end do ! itemp
 end do ! spin

 ! Compute RTA mobility
 call self%compute_rta_mobility(cryst, comm)

 if (my_rank == master) then
   ! Print RTA results to stdout and other external txt files (for the test suite)
   call self%print_rta_txt_files(cryst, dtset, dtfil)

   ! Creates the netcdf file used to store the results of the calculation.
   path = strcat(dtfil%filnam_ds(4), "_RTA.nc")
   call wrtout(units, ch10//sjoin("- Writing RTA transport results to:", path))
   NCF_CHECK(nctk_open_create(ncid, path , xmpi_comm_self))
   call self%rta_ncwrite(cryst, dtset, ncid)
   NCF_CHECK(nf90_close(ncid))
 end if

 call cwtime_report(" compute_rta", cpu, wall, gflops)

contains

 real(dp) function carriers(wmesh, dos, istart, istop, kT, mu)

 !Arguments -------------------------------------------
 real(dp),intent(in) :: kT, mu
 real(dp),intent(in) :: wmesh(self%nw), dos(self%nw)
 integer,intent(in) :: istart, istop

 !Local variables -------------------------------------
 integer :: iw
 real(dp) :: kernel(self%nw), integral(self%nw)

 kernel = zero
 do iw=istart,istop
   kernel(iw) = dos(iw) * occ_fd(wmesh(iw), kT, mu)
 end do
 call simpson_int(self%nw, edos_step, kernel, integral)
 carriers = integral(self%nw)

 end function carriers

 ! Compute L^\alpha(\mu, T) = \int de \sigma(e, T) (e - mu)^\alpha (-df/de)
 subroutine onsager(order, lorder)

 !Arguments -------------------------------------------
 integer,intent(in) :: order
 real(dp),intent(out) :: lorder(3, 3, self%nw, self%nsppol, self%ntemp, self%nrta)

 !Local variables -------------------------------------
 integer :: spin, iw, imu, irta
 real(dp) :: mu, ee, kT
 real(dp) :: kernel(self%nw,3,3,self%nsppol), integral(self%nw)

 ! Get spin degeneracy
 max_occ = two / (self%nspinor * self%nsppol)

 do irta=1,self%nrta
   do itemp=1,self%ntemp
     kT = self%kTmesh(itemp)
     ! Loop over chemical potentials mu
     do imu=1,self%nw
       mu = self%edos%mesh(imu)

       ! Build integrand for given mu
       do iw=1,self%nw
         ee = self%edos%mesh(iw)
         if (order > 0) then
           kernel(iw,:,:,:) = - max_occ * self%vvtau_dos(iw,:,:,itemp,:,irta) * (ee - mu)** order * occ_dfde(ee, kT, mu)
         else
           kernel(iw,:,:,:) = - max_occ * self%vvtau_dos(iw,:,:,itemp,:,irta) * occ_dfde(ee, kT, mu)
         end if
       end do

       ! Integrate with simpson_int
       do spin=1,self%nsppol
         do jj=1,3
           do ii=1,3
             call simpson_int(self%nw, edos_step, kernel(:,ii,jj, spin), integral)
             lorder(ii, jj, imu, spin, itemp, irta) = integral(self%nw)
           end do
         end do
       end do

     end do ! imu
   end do ! itemp
 end do ! irta

 end subroutine onsager

end subroutine compute_rta
!!***

!----------------------------------------------------------------------

!!****f* m_rta/compute_rta_mobility
!! NAME
!! compute_rta_mobility
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! comm=MPI communicator.
!!
!! SOURCE

subroutine compute_rta_mobility(self, cryst, comm)

!Arguments ------------------------------------
 class(rta_t),intent(inout) :: self
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: comm

!Local variables ------------------------------
 integer :: nsppol, nkibz, ib, ik_ibz, spin, ii, jj, itemp, ieh, cnt, nprocs, irta, time_opt
 real(dp) :: eig_nk, mu_e, linewidth, fact, fact0, max_occ, kT, wtk, cpu, wall, gflops
 real(dp) :: vr(3), vv_tens(3,3), vv_tenslw(3,3), work_33(3,3), mat33(3,3) !, tmp_tens(3,3)
!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 nprocs = xmpi_comm_size(comm)
 nkibz = self%ebands%nkpt; nsppol = self%ebands%nsppol

 time_opt = 0 ! This to preserve the previous behaviour in which TR was not used.
 !time_opt = -1 ! This to preserve the previous behaviour in which TR was not used.

 ABI_CALLOC(self%mobility_mu, (3, 3, 2, nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%conductivity_mu, (3, 3, 2, nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%l11_mu, (3, 3, nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%l12_mu, (3, 3, nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%l22_mu, (3, 3, nsppol, self%ntemp, self%nrta))
 ABI_CALLOC(self%resistivity_mu, (3, 3, self%ntemp, self%nrta))
 self%conductivity_mu=zero
 self%l11_mu=zero
 self%l12_mu=zero
 self%l22_mu=zero
 self%resistivity_mu=zero
! ABI_CALLOC(self%seebeck_mu, (3, 3, 2, nsppol, self%ntemp, self%nrta))
! ABI_CALLOC(self%kappa_mu, (3, 3, 2, nsppol, self%ntemp, self%nrta))
! ABI_CALLOC(self%peltier_mu, (3, 3, 2, nsppol, self%ntemp, self%nrta))
 ! Compute (e/h) carriers per unit cell at the different temperatures.
 ABI_CALLOC(self%n_ehst, (2, self%nsppol, self%ntemp))
 call self%ebands%get_carriers(self%ntemp, self%kTmesh, self%transport_mu_e, self%n_ehst)

 ! Compute conductivity_mu i.e. results in which lifetimes have been computed in a consistent way
 ! with the same the Fermi level. In all the other cases, indeed, we assume that tau does not depend on ef.
 !
 ! sigma_RTA = -S e^2 / (N_k Omega) sum_\nk (v_\nk \otimes v_\nk) \tau_\nk (df^0/de_\nk)
 !
 ! with S the spin degeneracy factor.
 !
 ! TODO: Implement other tensors. Compare these results with the ones obtained with spectral sigma
 ! In principle, they should be the same, in practice the integration of sigma requires enough resolution
 ! around the band edge.
 !print *, "in RTA max velocities", maxval(abs(self%vbks))
 cnt = 0
 do spin=1,nsppol
   do ik_ibz=1,nkibz
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     wtk = self%ebands%wtk(ik_ibz)

     do ib=self%bmin,self%bmax
       eig_nk = self%ebands%eig(ib, ik_ibz, spin)

       ! Store outer product in vv_tens
       vr(:) = self%vbks(:, ib, ik_ibz, spin)
       ! Don't remove this if: it makes the loop a bit faster and, most importantly,
       ! it prevents intel from miscompiling the code.
       if (all(abs(vr) == zero)) cycle

       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj) = vr(ii) * vr(jj)
         end do
       end do

       ! Symmetrize tensor.
       !print *, "intens", vv_tens
       vv_tens = cryst%symmetrize_cart_tens33(vv_tens, time_opt)
       !print *, "out_tens", vv_tens

       ! Multiply by the lifetime (SERTA or MRTA)
       do irta=1,self%nrta
         do itemp=1,self%ntemp
           kT = self%kTmesh(itemp)
           mu_e = self%transport_mu_e(itemp)
           ieh = 2; if (eig_nk >= mu_e) ieh = 1
           linewidth = self%linewidths(itemp, ib, ik_ibz, spin, irta)
           !print *, linewidth, wtk, occ_dfde(eig_nk, kT, mu_e), "tens", vv_tens
           call safe_div( - wtk * vv_tens * occ_dfde(eig_nk, kT, mu_e), two * linewidth, zero, vv_tenslw)
           self%conductivity_mu(:, :, ieh, spin, itemp, irta) = self%conductivity_mu(:, :, ieh, spin, itemp, irta) &
             + vv_tenslw(:, :)
           self%l11_mu(:, :, spin, itemp, irta) = sum(self%conductivity_mu(:,:,:, spin, itemp, irta), dim=3) !sum over ieh (holes and e-)
        !Here I implement the other onsager coefficients at the correct mu
           call safe_div( - wtk * vv_tens * occ_dfde(eig_nk, kT, mu_e)*(eig_nk-mu_e), two * linewidth, zero, vv_tenslw)
           self%l12_mu(:, :, spin, itemp, irta)=self%l12_mu(:, :, spin, itemp, irta)+vv_tenslw(:, :)
           call safe_div( - wtk * vv_tens * occ_dfde(eig_nk, kT, mu_e)*(eig_nk-mu_e)**2, two * linewidth, zero, vv_tenslw)
           self%l22_mu(:, :, spin, itemp, irta)=self%l22_mu(:, :, spin, itemp, irta)+vv_tenslw(:, :)
         end do
       end do
     end do

   end do ! ik_ibz
 end do ! spin

 !call xmpi_sum(self%conductivity_mu, comm, ierr)

 ! Get units conversion factor including spin degeneracy.
 max_occ = two / (self%nspinor * self%nsppol)
 fact0 = max_occ * (siemens_SI / Bohr_meter / cryst%ucvol) / 100
 self%conductivity_mu = fact0 * self%conductivity_mu  ! siemens cm^-1
 self%l11_mu = max_occ * self%l11_mu
 self%l12_mu = max_occ * self%l12_mu
 self%l22_mu = max_occ * self%l22_mu

!Same for resistivity_mu (at correct mu)

! ABI_MALLOC(self%resistivity_mu, (3, 3, self%ntemp, self%nrta))
 do irta=1,self%nrta
  do itemp=1,self%ntemp

    work_33 = sum(self%l11_mu(:,:,:,itemp,irta), dim=3) !sum over spins
    work_33 = work_33*fact0
    call inv33(work_33, mat33); mat33 = 1e+6_dp * mat33
    self%resistivity_mu(:, :, itemp, irta) = mat33
  end do
 end do


 !TODO: implement for mobility as well (verify)
 ! Scale by the carrier concentration
 fact = 100**3 / e_Cb
 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp
       do ieh=1,2 ! e/h
         call safe_div(fact * self%conductivity_mu(:,:,ieh,spin,itemp, irta), &
                       self%n_ehst(ieh, spin, itemp) / cryst%ucvol / Bohr_meter**3, zero, &
                       self%mobility_mu(:,:,ieh,spin,itemp,irta))
       end do
     end do
   end do
 end do

 call cwtime_report(" compute_rta_mobility", cpu, wall, gflops)

end subroutine compute_rta_mobility
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_ncwrite
!! NAME
!! rta_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! ncid=Netcdf file handle.
!!
!! SOURCE

subroutine rta_ncwrite(self, cryst, dtset, ncid)

!Arguments --------------------------------------
 class(rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: ncid

!Local variables --------------------------------
 integer :: ncerr, ii
 real(dp) :: cpu, wall, gflops
 real(dp) :: work(dtset%nsppol)
!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ! Write to netcdf file
 NCF_CHECK(cryst%ncwrite(ncid))
 NCF_CHECK(self%ebands%ncwrite(ncid))
 NCF_CHECK(self%edos%ncwrite(ncid))

 !nctk_copy from sigeph?
 !    nctkarr_t("eph_ngqpt_fine", "int", "three"), &

 ncerr = nctk_def_dims(ncid, [ &
    nctkdim_t("ntemp", self%ntemp), nctkdim_t("nrta", self%nrta), nctkdim_t("nsppol", self%nsppol), &
    nctkdim_t("nkcalc", self%nkcalc), nctkdim_t("nkibz", self%ebands%nkpt) &
 ], defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_arrays(ncid, [ &
    nctkarr_t('transport_ngkpt', "int", "three"), &
    nctkarr_t('sigma_erange', "dp", "two"), &
    nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
    nctkarr_t("kcalc2ebands", "int", "six, nkcalc"), &
    nctkarr_t("kibz", "dp", "three, nkibz"), &
    nctkarr_t('kTmesh', "dp", "ntemp"), &
    nctkarr_t('transport_mu_e', "dp", "ntemp"), &
    nctkarr_t('n_ehst', "dp", "two, nsppol, ntemp"), &
    nctkarr_t('eph_mu_e', "dp", "ntemp"), &
    nctkarr_t('vb_max', "dp", "nsppol"), &
    nctkarr_t('cb_min', "dp", "nsppol"), &
    nctkarr_t('vv_dos', "dp", "edos_nw, three, three, nsppol"), &
    nctkarr_t('vvtau_dos', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('tau_dos', "dp", "edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('L0', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('L1', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('L2', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('sigma', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('kappa', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('zte', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('seebeck', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('pi', "dp", "three, three, edos_nw, nsppol, ntemp, nrta"), &
    nctkarr_t('mobility', "dp", "three, three, edos_nw, ntemp, two, nsppol, nrta"), &
    nctkarr_t('conductivity', "dp", "three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('resistivity', "dp", "three, three, ntemp, nrta"), &
    nctkarr_t('N', "dp", "edos_nw, ntemp, two"), &
    !nctkarr_t('conductivity_mu',"dp", "three, three, two, nsppol, ntemp, nrta")], &
    nctkarr_t('mobility_mu', "dp", "three, three, two, nsppol, ntemp, nrta")], &
 defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: "assume_gap"])
 NCF_CHECK(ncerr)
 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
    "eph_extrael", "eph_fermie", "transport_extrael", "transport_fermie"])
 NCF_CHECK(ncerr)

 ! Write data.
 ii = 0; if (self%assume_gap) ii = 1
 ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: "assume_gap"], [ii], datamode=.True.)

 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   "eph_extrael", "eph_fermie", "transport_extrael", "transport_fermie"], &
   [self%eph_extrael, self%eph_fermie, self%transport_extrael, self%transport_fermie])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "transport_ngkpt"), dtset%transport_ngkpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_erange"), dtset%sigma_erange))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ibz"), self%kcalc2ibz))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ebands"), self%kcalc2ebands))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kibz"), self%ebands%kptns))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
 if (self%assume_gap) then
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vb_max"), self%gaps%vb_max))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cb_min"), self%gaps%cb_min))
 else
   ! Set vbm and cbm to fermie if metal.
   work(:) = self%ebands%fermie
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vb_max"), work))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cb_min"), work))
 end if
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eph_mu_e"), self%eph_mu_e))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "transport_mu_e"), self%transport_mu_e))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "n_ehst"), self%n_ehst))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vv_dos"), self%vv_dos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vvtau_dos"),  self%vvtau_dos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "tau_dos"),  self%tau_dos))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L0"), self%l0))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L1"), self%l1))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "L2"), self%l2))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma"),   self%sigma))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kappa"),   self%kappa))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "zte"),   self%zte))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "seebeck"), self%seebeck))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "pi"),      self%pi))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "N"), self%n))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility"), self%mobility))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "conductivity"), self%conductivity))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "resistivity"), self%resistivity))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility_mu"), self%mobility_mu))
 !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "conductivity_mu"), self%conductivity_mu))
 !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "resistivity_mu"), self%resistivity_mu))

 call cwtime_report(" rta_ncwrite", cpu, wall, gflops)

end subroutine rta_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_rta/print_rta_txt_files
!! NAME
!! print_rta_txt_files
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! dtfil<datafiles_type>=variables related to files.
!!
!! SOURCE

subroutine print_rta_txt_files(self, cryst, dtset, dtfil)

!Arguments --------------------------------------
 class(rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil

!Local variables --------------------------------
 integer :: itemp, spin, irta, ii, nsp
 real(dp) :: TKelv
 character(len=500) :: msg, pre, rta_type
 integer :: units(2)
 character(len=2) :: components(3)
 real(dp) :: mat33(3,3),  work33(3,3)
!************************************************************************

 units = [std_out, ab_out]
 call wrtout(units, ch10//' Transport (RTA) calculation results, chemical potential adjusted with T:', newlines=1)
 components = ["xx", "yy", "zz"]

 do irta=1,self%nrta
   if (irta == 1) rta_type = "SERTA"
   if (irta == 2) rta_type = "MRTA"

   if (self%assume_gap) then
     ! SemiConductor
     do ii=1,3
       call wrtout(units, sjoin(" Cartesian component of", rta_type, "mobility tensor:", components(ii)))
       write(msg, "(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
       call wrtout(units, msg)

       do spin=1,self%nsppol
         if (self%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)

         do itemp=1,self%ntemp
           write(msg,"(f16.2,2e16.2,2f16.2)") &
             self%kTmesh(itemp) / kb_HaK, &
             self%n_ehst(1, spin, itemp) / cryst%ucvol / Bohr_cm**3, &
             self%n_ehst(2, spin, itemp) / cryst%ucvol / Bohr_cm**3, &
             self%mobility_mu(ii, ii, 1, spin, itemp, irta), &
             self%mobility_mu(ii, ii, 2, spin, itemp, irta)
           call wrtout(units, msg)
         end do ! itemp
       end do ! spin
       call wrtout(units, ch10)
     end do ! ii

   else
     ! Metals. Print conductivity (spin resolved) and resistivity (no spin resolved)
!    do ii=1,2
!      if (ii == 1) msg = sjoin(" Conductivity [Siemens cm^-1] using ", rta_type, "approximation")
!      if (ii == 2) msg = sjoin(" Resistivity [micro-Ohm cm] using ", rta_type, "approximation")
!      call wrtout(units, msg)
!
!      nsp = self%nsppol; if (ii == 2) nsp = 1
!      do spin=1,nsp
!        if (nsp == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
!        write(msg, "(4a16)") 'Temperature (K)', 'xx', 'yy', 'zz'
!        call wrtout(units, msg)
!        do itemp=1,self%ntemp
!          if (ii == 1) then
!            mat33 = self%conductivity(:,:,itemp,spin,irta)
!          else
!            mat33 = self%resistivity(:,:,itemp,irta)
!          end if
!          write(msg,"(f16.2,3e16.6)") self%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3)
!          call wrtout(units, msg)
!        end do !itemp
!      end do !spin
!      call wrtout(units, ch10)
!    end do

     !For coefficients calculated at exact mu (i.e for mu calculated at exact temp.):
     !TODO: do the same for SM and thus mobility

     do ii=1,2
       if (ii == 1) msg = sjoin(" Conductivity [Siemens cm^-1] using ", rta_type, "approximation")
       if (ii == 2) msg = sjoin(" Resistivity [micro-Ohm cm] using ", rta_type, "approximation")
       call wrtout(units, msg)

       do itemp=1, self%ntemp
         nsp = self%nsppol; if (ii == 2) nsp = 1
         do spin=1,nsp
           if (nsp == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
           write(msg, "(4a16)") 'Temperature (K)', 'xx', 'yy', 'zz'
           call wrtout(units, msg)

           if (ii == 1) then
             mat33 = sum(self%conductivity_mu(:,:,:,spin, itemp, irta), dim=3)
           else
             mat33 = self%resistivity_mu(:,:,itemp,irta)
           end if
           write(msg,"(f16.2,3e16.6)") self%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3)
           call wrtout(units, msg)
         end do !spin
       end do !itemp
       call wrtout(units, ch10)
     end do
   end if

   msg = sjoin(" Seebeck [Volts / Kelvin] using ", rta_type, "approximation")
   call wrtout(units, msg)
   do itemp=1, self%ntemp
     !TODO : what is the following line? Resets 0 Kelvin to 1 Kelvin minimum T??? a bit arbitrary. Should handle T=0 more cleanly
     TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one
     nsp = self%nsppol;
     do spin=1,nsp
       if (nsp == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(4a16)") 'Temperature (K)', 'xx', 'yy', 'zz'
       call wrtout(units, msg)
       call inv33(self%l11_mu(:, :, spin, itemp, irta), work33)
       mat33 = - (volt_SI / TKelv) * matmul(work33, self%l12_mu(:,:,spin,itemp,irta))
       write(msg,"(f16.2,3e16.6)") TKelv, mat33(1,1), mat33(2,2), mat33(3,3)
       call wrtout(units, msg)
     end do !spin
   end do !itemp
   call wrtout(units, ch10)

   msg = sjoin(" Kappa [W/m*K] using ", rta_type, "approximation")
   call wrtout(units, msg)
   do itemp=1, self%ntemp
     TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one
     nsp = self%nsppol;
     do spin=1,nsp
       if (nsp == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(4a16)") 'Temperature (K)', 'xx', 'yy', 'zz'
       call wrtout(units, msg)
         call inv33(self%l11_mu(:, :, spin, itemp, irta), work33)
         mat33 = + (volt_SI**2 * (siemens_SI / Bohr_meter / cryst%ucvol) / TKelv) * (self%l22_mu(:,:,spin,itemp,irta) - matmul(self%l12_mu(:,:,spin,itemp,irta),matmul(work33,self%l12_mu(:,:,spin,itemp,irta))))
       write(msg,"(f16.2,3e16.6)") TKelv, mat33(1,1), mat33(2,2), mat33(3,3)
       call wrtout(units, msg)
     end do !spin
   end do !itemp
   call wrtout(units, ch10)

   msg = sjoin(" Peltier [Volts] using ", rta_type, "approximation")
   call wrtout(units, msg)
   do itemp=1, self%ntemp
     TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one
     nsp = self%nsppol;
     do spin=1,nsp
       if (nsp == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(4a16)") 'Temperature (K)', 'xx', 'yy', 'zz'
       call wrtout(units, msg)
         call inv33(self%l11_mu(:, :, spin, itemp, irta), work33)
         mat33 = - volt_SI * matmul(self%l12_mu(:,:,spin,itemp,irta), work33)
       write(msg,"(f16.2,3e16.6)") TKelv, mat33(1,1), mat33(2,2), mat33(3,3)
       call wrtout(units, msg)
     end do !spin
   end do !itemp
   call wrtout(units, ch10)



 end do ! irta

 do irta=1,self%nrta
   select case (irta)
   case (1)
     pre = "_SERTA"
   case (2)
     pre = "_MRTA"
   case default
     ABI_ERROR(sjoin("Don't know how to handle irta:", itoa(irta)))
   end select
   call self%write_tensor(dtset, irta, "sigma", self%sigma(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SIGMA"))
   call self%write_tensor(dtset, irta, "seebeck", self%seebeck(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SBK"))
   call self%write_tensor(dtset, irta, "kappa", self%kappa(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_KAPPA"))
   call self%write_tensor(dtset, irta, "zte", self%zte(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_ZTE"))
   call self%write_tensor(dtset, irta, "pi", self%pi(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_PI"))
 end do

end subroutine print_rta_txt_files
!!***

!----------------------------------------------------------------------

!!****f* m_rta/write_tensor
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine write_tensor(self, dtset, irta, header, values, path)

!Arguments --------------------------------------
 class(rta_t),intent(in) :: self
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: irta
 character(len=*),intent(in) :: header
 real(dp),intent(in) :: values(:,:,:,:,:)
 character(len=*),intent(in) :: path

!Local variables --------------------------------
 integer :: itemp, iw, ount
 character(len=500) :: msg, rta_type
 real(dp),allocatable :: tmp_values(:,:,:,:,:)
!************************************************************************

 if (open_file(trim(path), msg, newunit=ount, form="formatted", action="write", status='unknown') /= 0) then
   ABI_ERROR(msg)
 end if

 if (irta == 1) rta_type = "RTA type: Self-energy relaxation time approximation (SERTA)"
 if (irta == 2) rta_type = "RTA type: Momentum relaxation time approximation (MRTA)"

 ! write header
 write(ount, "(2a)")"# ", trim(header)
 write(ount, "(2a)")"# ", trim(rta_type)
 ! TODO: Units ?
 write(ount, "(a, 3(i0, 1x))")"#", dtset%transport_ngkpt
 write(ount, "(a)")"#"

 ! This to improve portability of the unit tests.
 call alloc_copy(values, tmp_values)
 where (abs(values) > tol30)
   tmp_values = values
 else where
   tmp_values = zero
 end where

 ! (nw, 3, 3, nsppol, ntemp)
 if (self%nsppol == 1) then
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)")"# Energy [Ha], (xx, yx, zx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor."
     do iw=1,self%nw
       write(ount, "(10(es16.6))")self%edos%mesh(iw), tmp_values(:, :, iw, 1, itemp)
     end do
   end do
  write(ount, "(a)")""
 else
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)") &
       "# Energy [Ha], (xx, yx, zx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor for spin up followed by spin down."
     do iw=1,self%nw
       write(ount, "(19(es16.6))")self%edos%mesh(iw), tmp_values(:, :, iw, 1, itemp), tmp_values(:, :, iw, 2, itemp)
     end do
   end do
  write(ount, "(a)")""
 end if

 close(ount)

 ABI_FREE(tmp_values)

end subroutine write_tensor
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_free
!! NAME
!! rta_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! INPUTS
!!
!! SOURCE

subroutine rta_free(self)

!Arguments --------------------------------------
 class(rta_t),intent(inout) :: self

 ABI_SFREE(self%n)
 ABI_SFREE(self%vv_dos)
 ABI_SFREE(self%vvtau_dos)
 ABI_SFREE(self%tau_dos)
 ABI_SFREE(self%bstart_ks)
 ABI_SFREE(self%bstop_ks)
 ABI_SFREE(self%nbcalc_ks)
 ABI_SFREE(self%kcalc2ibz)
 ABI_SFREE(self%kcalc2ebands)
 ABI_SFREE(self%kTmesh)
 ABI_SFREE(self%eminmax_spin)
 ABI_SFREE(self%eph_mu_e)
 ABI_SFREE(self%transport_mu_e)
 ABI_SFREE(self%vbks)
 ABI_SFREE(self%linewidths)
 ABI_SFREE(self%l0)
 ABI_SFREE(self%l1)
 ABI_SFREE(self%l2)
 ABI_SFREE(self%sigma)
 ABI_SFREE(self%mobility)
 ABI_SFREE(self%conductivity)
 ABI_SFREE(self%resistivity)
 ABI_SFREE(self%seebeck)
 ABI_SFREE(self%kappa)
 ABI_SFREE(self%zte)
 ABI_SFREE(self%pi)
 ABI_SFREE(self%mobility_mu)
 ABI_SFREE(self%conductivity_mu)
 ABI_SFREE(self%resistivity_mu)
 ABI_SFREE(self%l11_mu)
 ABI_SFREE(self%l12_mu)
 ABI_SFREE(self%l22_mu)
 ABI_SFREE(self%n_ehst)

 call self%ebands%free()
 call self%gaps%free()
 call self%edos%free()

end subroutine rta_free
!!***

!----------------------------------------------------------------------

!!****f* m_rta/ibte_driver
!! NAME
!! ibte_driver
!!
!! FUNCTION
!! Driver to compute transport properties within the IBTE.
!!
!! INPUTS
!! dtfil<datafiles_type>=variables related to files.
!! ngfftc(18)=Coarse FFT mesh
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! cryst<crystal_t>=Crystalline structure
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine ibte_driver(dtfil, ngfftc, dtset, ebands, cryst, pawtab, psps, comm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfftc(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
 integer,parameter :: master = 0
 integer :: iet, max_et
 integer :: spin, ikcalc, nkcalc, nbsum, nbcalc, itemp, iter, ierr, bsize
 integer :: nkibz, nsppol, band_k, ik_ibz, bmin, bmax, band_sum, ntemp, ii, jj, iq_sum, nsp
 integer :: ikq_ibz, isym_kq, trev_kq, cnt, tag, nprocs, receiver, my_rank, isym, itime, isym_lgk
 integer :: ncid, grp_ncid, ncerr
 real(dp) :: kT, mu_e, e_nk, dfde_nk, tau_nk, lw_nk, max_adiff, cpu, wall, gflops, abs_tol, rtmp
 logical :: send_data
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(rta_t) :: ibte
!arrays
 integer :: units(2), dims(4)
 logical,allocatable :: converged(:)
 real(dp) :: vec3(3), sym_vec(3), mat33(3,3), f_kq(3), work33(3,3)
 real(dp) :: inv_sig_p(3,3)
 real(dp) :: fsum_eh(3,2,ebands%nsppol), max_adiff_spin(ebands%nsppol)
 real(dp) :: onsager(3,3,3,ebands%nsppol)
 real(dp),pointer :: sig_p(:,:,:,:), mob_p(:,:,:,:), sbk_p(:,:,:), kappa_p(:,:,:), pi_p(:,:,:)
 real(dp),target,allocatable :: ibte_sigma(:,:,:,:,:), ibte_mob(:,:,:,:,:), ibte_rho(:,:,:)
 real(dp),target,allocatable :: ibte_seebeck(:,:,:,:), ibte_kappa(:,:,:,:), ibte_pi(:,:,:,:)
 real(dp),allocatable :: grp_srate(:,:,:,:), fkn_in(:,:,:,:), fkn_out(:,:,:,:), fkn_efield(:,:,:,:), fkn_serta(:,:,:,:), taukn_serta(:,:,:,:)
 real(dp) :: fact_sbk
 character(len=2) :: components(3)
 real(dp), allocatable :: sig_gen(:,:,:,:), mob_gen(:,:,:,:), sig_l21(:,:,:,:), sig_l22(:,:,:,:), mob_21(:,:,:,:), mob_22(:,:,:,:)
 type :: scatk_t

   integer :: rank = xmpi_undefined_rank

   integer :: nq_ibzk_eff
   ! Number of effective q-points in the IBZ(k)

   integer :: lgk_nsym
   ! Number of symmetry operations in the little group of k.

   integer,allocatable :: lgk_sym2glob(:,:)
   ! lgk_sym2glob(2, lgk_nsym)
   ! Mapping isym_lg --> [isym, itime]
   ! where isym is the index of the operation in the global array **crystal%symrec**
   ! and itim is 2 if time-reversal T must be included else 1. Depends on ikcalc

   integer,allocatable :: kq_symtab(:,:)
   ! kq_symtab(6, nq_ibzk_eff)

   real(dp),allocatable :: vals(:,:,:,:)
   ! (nq_ibzk_eff, nbsum, nbcalc, ntemp)
 end type scatk_t

 type(scatk_t),target,allocatable :: sr(:,:)
 type(scatk_t),pointer :: sr_p

! *************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 units = [std_out, ab_out]

 call wrtout(units, ch10//' Entering IBTE driver.')
 call wrtout(units, sjoin("- Reading SERTA lifetimes and e-ph scattering operator from:", &
             dtfil%filsigephin), newlines=1, do_flush=.True.)

 ! Initialize IBTE object
 ibte = rta_new(dtset, dtfil, ngfftc, cryst, ebands, pawtab, psps, comm)

 ! Compute RTA transport quantities
 call ibte%compute_rta(cryst, dtset, dtfil, comm)

 nkcalc = ibte%nkcalc
 nkibz = ibte%ebands%nkpt; nsppol = ibte%nsppol; ntemp = ibte%ntemp
 bmin = ibte%bmin; bmax = ibte%bmax
 !call wrtout(std_out, sjoin(" nkcalc", itoa(nkcalc), "bmin:", itoa(bmin), "bmax:", itoa(bmax)))

 ! Get units conversion factor including spin degeneracy.
!max_occ = two / (self%nspinor * self%nsppol)
!fact_sigma = max_occ * (siemens_SI / Bohr_meter) / 100.  ! Siemens / cm
! fact_mob   = fact_sigma * 100.**3 / e_Cb * Bohr_meter**3 ! cm^2 / V / s
 fact_sbk   = -volt_SI * kb_HaK                           ! Volt / Kelvin
!call wrtout(std_out,"factors: ")
! write (std_out,*) fact_sigma,  fact_sbk

 !call ibte%read_scattering()
 ! Loops and memory are distributed over k-points and collinear spins
 ABI_MALLOC(sr, (nkcalc, nsppol))
 cnt = 0
 do spin=1,nsppol
   do ikcalc=1,nkcalc
     cnt = cnt + 1
     sr(ikcalc, spin)%rank = mod(cnt, nprocs)
   end do
 end do

 call cwtime(cpu, wall, gflops, "start")
 ! Master reads and sends data to the rank treating (ikcalc, spin).
 if (my_rank == master) then
   NCF_CHECK(nctk_open_read(ncid, dtfil%filsigephin, xmpi_comm_self))
 end if

 do spin=1,nsppol
   do ikcalc=1,nkcalc
     sr_p => sr(ikcalc, spin)
     receiver = sr_p%rank
     send_data = master /= receiver
     if (.not. any(my_rank == [master, receiver])) cycle
     !call wrtout(std_out, sjoin(" Sending data from my_rank:", itoa(my_rank), " to:", itoa(receiver)))

     if (my_rank == master) then
       ! Get ncid of group used to store scattering rate for this k-point.
       ncerr = nf90_inq_ncid(ncid, strcat("srate_k", itoa(ikcalc), "_s", itoa(spin)), grp_ncid)
       if (ncerr /= NF90_NOERR) then
         ABI_ERROR("Cannot find collision terms in SIGEPH file. Rerun eph_task -4 step with ibte_prep 1.")
       end if
       NCF_CHECK(nctk_get_dim(grp_ncid, "nq_ibzk_eff", sr_p%nq_ibzk_eff))
       NCF_CHECK(nctk_get_dim(grp_ncid, "nbsum", nbsum))
       NCF_CHECK(nctk_get_dim(grp_ncid, "nbcalc", nbcalc))
       NCF_CHECK(nctk_get_dim(grp_ncid, "lgk_nsym", sr_p%lgk_nsym))
       dims = [sr_p%nq_ibzk_eff, nbsum, nbcalc, sr_p%lgk_nsym]
     end if

     if (send_data) then
       tag = size(dims)
       if (my_rank == master) call xmpi_send(dims, receiver, tag, comm, ierr)
       if (my_rank == receiver) then
         call xmpi_recv(dims, master, tag, comm, ierr)
         sr_p%nq_ibzk_eff = dims(1); nbsum = dims(2); nbcalc = dims(3); sr_p%lgk_nsym = dims(4)
       end if
     end if

     ! Note that the size along the (n, m) axis does not depend on the kcalc index.
     ABI_CALLOC(sr_p%vals, (sr_p%nq_ibzk_eff, bmin:bmax, bmin:bmax, ntemp))
     ABI_MALLOC(sr_p%kq_symtab, (6, sr_p%nq_ibzk_eff))
     ABI_MALLOC(sr_p%lgk_sym2glob, (2, sr_p%lgk_nsym))

     if (my_rank == master) then
       NCF_CHECK(nf90_get_var(grp_ncid, nctk_idname(grp_ncid, "kq_symtab"), sr_p%kq_symtab))
       NCF_CHECK(nf90_get_var(grp_ncid, nctk_idname(grp_ncid, "lgk_sym2glob"), sr_p%lgk_sym2glob))

       ! Note that on file, we have:
       !
       !      nctkarr_t("srate", "dp", "nq_ibzk_eff, nbsum, nbcalc, ntemp")
       !
       ! but in terms of n, m indices we have that the:
       !
       !   n index: bstart_ks bstop_ks
       !   m index: bsum_start up to bsum_stop to account for phonon emission/absorption.
       !
       ! so we have to insert the values in bmin:bmax slice.
       ! TODO: Recheck this part.
       ABI_MALLOC(grp_srate, (sr_p%nq_ibzk_eff, nbsum, nbcalc, ntemp))
       NCF_CHECK(nf90_get_var(grp_ncid, nctk_idname(grp_ncid, "srate"), grp_srate))
       ii = ibte%bstart_ks(ikcalc, spin)
       jj = ibte%bstop_ks(ikcalc, spin)
       sr_p%vals(:, bmin:bmax, ii:jj, :) = grp_srate(:, 1:bmax-bmin+1, 1:nbcalc, :)
       ABI_SFREE(grp_srate)
     end if

     if (send_data) then
       tag = size(sr_p%vals)
       if (my_rank == master) then
         tag = tag + 1; call xmpi_send(sr_p%vals, receiver, tag, comm, ierr)
         tag = tag + 1; call xmpi_send(sr_p%kq_symtab, receiver, tag, comm, ierr)
         tag = tag + 1; call xmpi_send(sr_p%lgk_sym2glob, receiver, tag, comm, ierr)
       end if
       if (my_rank == receiver) then
         tag = tag + 1; call xmpi_recv(sr_p%vals, master, tag, comm, ierr)
         tag = tag + 1; call xmpi_recv(sr_p%kq_symtab, master, tag, comm, ierr)
         tag = tag + 1; call xmpi_recv(sr_p%lgk_sym2glob, master, tag, comm, ierr)
       end if
     end if

     if (send_data .and. my_rank /= receiver) call free_sr_ks(ikcalc, spin)

   end do ! spin
 end do ! ikcalc

 if (my_rank == master) then
   NCF_CHECK(nf90_close(ncid))
 end if

 call cwtime_report(" sigeph IO", cpu, wall, gflops)

 ! Solve the linearized BTE with B = 0.
 !
 !   F_\nk = e df/de_\nk v_\nk \tau^0 + \tau^0 \sum_{mq} Srate_{nk,mq} F_{m,k+q}
 !
 ! where F is a vector in Cartesian coordinates and tau^0 is the SERTA relaxation time.
 !
 ! Take advantage of the following symmetry properties:
 !
 ! 1. F_k = F_Sk.
 ! 2. F_{-k} = -F_k if TR symmetry.
 ! 3. The q-space integration is reduced to the IBZ(k) using the symmetries of the little group of k.

 !call ibte%solve_ibte(solver_type=1)

 bsize = bmax - bmin + 1
 ABI_CALLOC(fkn_in, (3, nkibz, bmin:bmax, nsppol))
 ABI_CALLOC(fkn_out, (3, nkibz, bmin:bmax, nsppol))
 ABI_CALLOC(fkn_serta, (3, nkibz, bmin:bmax, nsppol))
 ABI_CALLOC(taukn_serta, (3, nkibz, bmin:bmax, nsppol))
 ABI_MALLOC(ibte_sigma, (3, 3, 2, nsppol, ntemp))
 ABI_MALLOC(ibte_seebeck, (3, 3, nsppol, ntemp))
 ABI_MALLOC(ibte_pi, (3, 3, nsppol, ntemp))
 ABI_MALLOC(ibte_mob, (3, 3, 2, nsppol, ntemp))
 ABI_MALLOC(converged, (ntemp))
 ABI_MALLOC(sig_gen, (3, 3, 2, nsppol))
 ABI_MALLOC(sig_l21, (3, 3, 2, nsppol))
 ABI_MALLOC(sig_l22, (3, 3, 2, nsppol))
 ABI_MALLOC(mob_gen, (3, 3, 2, nsppol))
 ABI_MALLOC(mob_21, (3, 3, 2, nsppol))
 ABI_MALLOC(mob_22, (3, 3, 2, nsppol))
 ABI_MALLOC(ibte_kappa, (3, 3, nsppol, ntemp))
 abs_tol = dtset%ibte_abs_tol

 ! If the fermi level is inside the gap, F_k is gonna be very small
 ! hence once should use a much smaller tolerance to converge.
 ! According to numerical tests, a reasonable value of abs_tol can be estimated from the free carrier density using:
 !
 !  1e-20 * e_density (in cm**-3)
 !
 ! These are the numerical values used to derive the fit:

 !  max_adiff = np.array([9e-12, 5.6e-6, 2.7e-6, 8.2e-60, 2.9e-46, 6.9e-6, 2.9e-8, 9.5E+00, 7.3E-04, 9.2E-04, 8.4E-04])
 !  e_density = np.array([0.97e8, 0.86e13, 0.26e14, 0.1e-40, 0.89e-26, 0.45e14, 0.28e12, 0.10E+19, 0.12E+16, 0.10E+15, 0.90E+14])

 if (abs_tol <= zero) then
   rtmp = minval(ibte%n_ehst, mask=ibte%n_ehst > zero) * ibte%nsppol
   abs_tol = 1e-20 * rtmp / cryst%ucvol / Bohr_cm**3
   call wrtout(std_out, " Input ibte_abs_tol <= zero ==> computing abs tolerance from minimal carrier density over all T")
   call wrtout(std_out, " using: abs_tol = 1e-20 * e_density (in cm**-3)")
   call wrtout(std_out, sjoin(" abs_tol:", ftoa(abs_tol), " from carrier_density:", ftoa(rtmp / cryst%ucvol / Bohr_cm**3)))
 end if

 if (my_rank == master) then
   path = strcat(dtfil%filnam_ds(4), "_RTA.nc")
   call wrtout(units, ch10//sjoin("- Writing IBTE transport results to:", path))
   NCF_CHECK(nctk_open_modify(ncid, path , xmpi_comm_self))

   ncerr = nctk_def_dims(ncid, [ &
      nctkdim_t("nkibz", nkibz), nctkdim_t("bsize", bsize), nctkdim_t("nkcalc", ibte%nkcalc) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t('fkn_out_sigma', "dp", "three, nkibz, bsize, nsppol, ntemp"), &
     nctkarr_t('ibte_sigma', "dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('ibte_mob', "dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t("kcalc2ibz", "int", "nkcalc, six"), &
     nctkarr_t("kcalc2ebands", "int", "six, nkcalc"), &
     nctkarr_t("kibz", "dp", "three, nkibz"), &
     nctkarr_t('ibte_rho', "dp", "three, three, ntemp") &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   NCF_CHECK(nctk_set_datamode(ncid))
 end if

 cnt = 0
 do itemp=1,ntemp
 !do itemp=ntemp, 1, -1
   cnt = cnt + 1
   kT = max(ibte%kTmesh(itemp), one * kb_HaK)
   mu_e = ibte%eph_mu_e(itemp)
   sig_p => ibte_sigma(:,:,:,:,itemp)
   mob_p => ibte_mob(:,:,:,:,itemp)
   sbk_p => ibte_seebeck(:,:,:,itemp)
   kappa_p => ibte_kappa(:,:,:,itemp)
   pi_p => ibte_pi(:,:,:,itemp)

   call wrtout(std_out," Value of kT: ", pre_newlines=1, newlines=1)
        write (std_out,*) kT


   ! Precompute tau_serta and fkn_serta for this T: f^'_nk v_\nk * \tau^0
   do spin=1,nsppol
     do ikcalc=1,nkcalc
       !ik_ibz = ibte%kcalc2ibz(ikcalc, 1)
       ik_ibz = ibte%kcalc2ebands(1, ikcalc)
       do band_k=ibte%bstart_ks(ikcalc, spin), ibte%bstop_ks(ikcalc, spin)
         lw_nk = ibte%linewidths(itemp, band_k, ik_ibz, spin, 1)  ! 1 --> SERTA linewidths.
         call safe_div(one, two * lw_nk, zero, tau_nk)
         taukn_serta(:, ik_ibz, band_k, spin) = tau_nk
         e_nk = ebands%eig(band_k, ik_ibz, spin)
         dfde_nk = occ_dfde(e_nk, kT, mu_e)
         fkn_serta(:, ik_ibz, band_k, spin) = tau_nk * dfde_nk * ibte%vbks(:, band_k, ik_ibz, spin)
       end do
     end do
   end do

   call wrtout(std_out, sjoin(" Begin IBTE loop for itemp:", itoa(itemp), ", KT:", ftoa(kT / kb_HaK), "[K]"), &
               pre_newlines=1, newlines=1)

   ! iter = 0 --> Compute SERTA transport tensors just for initial reference.
   call ibte_calc_tensors(ibte, cryst, itemp, kT, mu_e, fkn_serta, onsager, sig_p, mob_p, fsum_eh, comm, iet)

   !TODO fix dtset here ans use maxocc, CAREFUL, it has to be included in all L coeff. I think
!max_occ = two / (self%nspinor * self%nsppol)
! fact_sigma = max_occ * (siemens_SI / Bohr_meter) / 100.  ! Siemens / cm
! fact_mob   = fact_sigma * 100.**3 / e_Cb * Bohr_meter**3 ! cm^2 / V / s
 !fact_sbk   = volt_SI * kb_HaK                           ! Volt / Kelvin
  ! sig_p = fact_sigma * sig_p  ! siemens cm^-1
 !  mob_p = fact_mob   * mob_p  !

   ! Print mobility for semiconductors, conductivity for metals.
   if (ibte%assume_gap) then
     do spin=1,nsppol
       mat33 = sum(mob_p(:,:,:,spin), dim=3)
       write(msg, "(i5,1x,es9.1, *(1x, f16.2))") &
         0, zero, mat33(1,1), mat33(2,2), mat33(3,3), sum(fsum_eh(:,:,spin))
     end do
   else
     do spin=1,nsppol
       mat33 = sum(sig_p(:,:,:,spin), dim=3)
       write(msg, "(i5,1x,es9.1, *(1x, es16.6))") &
         0, zero, mat33(1,1), mat33(2,2), mat33(3,3), sum(fsum_eh(:,:,spin))
     end do
   end if
   call wrtout(std_out, msg)

   fkn_in = fkn_serta
   ! TODO: B-field
   ! Initialize fkn_in either from SERTA or from previous T.
   !if (cnt == 1) fkn_in = fkn_serta
   !if (cnt > 1 ) fkn_in = fkn_out
   fkn_out = zero

   max_et=2 !max_et =1 pour F^E et max_et=2 pour F^E et F^T, iet = 3 pour calculer kappa

! start here with loop over E and T perturbations: once F^Eps is done, init F^T and converge it as well
   electherm_loop: do iet = 1, max_et

     if (iet == 2) then
     ! initialize F^T from F^E
       do spin=1,nsppol
          do ikcalc=1,nkcalc
            ik_ibz = ibte%kcalc2ebands(1, ikcalc)
            do band_k=ibte%bstart_ks(ikcalc, spin), ibte%bstop_ks(ikcalc, spin)
               e_nk = ebands%eig(band_k, ik_ibz, spin)
              !Which one of the following is correct ?
              ! fkn_in(:, ik_ibz, band_k, spin)  = fkn_serta(:, ik_ibz, band_k, spin) * (e_nk-mu_e)/kT
               !fkn_in(:, ik_ibz, band_k, spin)  = fkn_efield(:, ik_ibz, band_k, spin) * (e_nk-mu_e)/kT
               fkn_in(:, ik_ibz, band_k, spin)  = fkn_efield(:, ik_ibz, band_k, spin) * (e_nk-mu_e)
               fkn_serta(:, ik_ibz, band_k, spin)  = fkn_serta(:, ik_ibz, band_k, spin) * (e_nk-mu_e)
            end do !band_k
          end do !ikcalc
        end do !spin
       ABI_CHECK(converged(itemp), "Warning: E field BTE not converged, I will not try to converge the T gradient")
       if (.not. converged(itemp)) exit electherm_loop
     end if ! iet == 2  ?????


     ! Begin iterative solver.
     iter_loop: do iter=1,dtset%ibte_niter

       ! call wrtout(std_out," check beginning loop, iter, iet:", pre_newlines=1, newlines=1)
       ! write (std_out,*) iter, iet

        ! Loop over the nk index in F_nk.
       do spin=1,nsppol
          do ikcalc=1,nkcalc
            sr_p => sr(ikcalc, spin)
            if (sr_p%rank /= my_rank) cycle ! MPI parallelism
            ik_ibz = ibte%kcalc2ebands(1, ikcalc)
            do band_k=ibte%bstart_ks(ikcalc, spin), ibte%bstop_ks(ikcalc, spin)

              ! Summing over the q-points in the effective IBZ(k) and the m band index.
              ! Results stored in vec3. Integration weights are already included.
              vec3 = zero
              do band_sum=ibte%bmin, ibte%bmax
                do iq_sum=1, sr_p%nq_ibzk_eff
                  ikq_ibz = sr_p%kq_symtab(1, iq_sum); isym_kq = sr_p%kq_symtab(2, iq_sum)
                  trev_kq = sr_p%kq_symtab(6, iq_sum) !; g0_kq = sr_p%kq_symtab(3:5, iq_sum)
                  ! Build F_{m,k+q} in the effective IBZ(k) from fkn_in using symmetries (need k+q --> IBZ map)
                  ! Use transpose(R) because we are using the tables for the wavefunctions
                  ! In this case listkk has been called with symrec and use_symrec=False
                  ! so q_bz = S^T q_ibz where S is the isym_kq symmetry
                  ! vkq = matmul(transpose(cryst%symrel_cart(:,:,isym_kq)), vkq)
                  mat33 = transpose(cryst%symrel_cart(:,:,isym_kq))
                  f_kq = matmul(mat33, fkn_in(:, ikq_ibz, band_sum, spin))
                  if (trev_kq == 1) f_kq = -f_kq
                  vec3 = vec3 + sr_p%vals(iq_sum, band_sum, band_k, itemp) * f_kq(:)
                end do ! iq_sum
              end do ! band_k

              ! Symmetrize intermediate results using the operations of the little group of k.
              sym_vec = zero
              do isym_lgk=1,sr_p%lgk_nsym
                isym = sr_p%lgk_sym2glob(1, isym_lgk)
                itime = sr_p%lgk_sym2glob(2, isym_lgk)
                mat33 = transpose(cryst%symrel_cart(:,:,isym))
                !if(itime == 1) mat33 = -mat33 ! FIXME: here there's a different convention for TR used in m_lgroup
                if (itime == 2) mat33 = -mat33
                sym_vec = sym_vec + matmul(mat33, vec3)
              end do
              sym_vec = taukn_serta(:, ik_ibz, band_k, spin) * sym_vec / sr_p%lgk_nsym
              ! if iet 2 serta and fkn_in have been multiplied by e-mu / T
              fkn_out(:, ik_ibz, band_k, spin) = fkn_serta(:, ik_ibz, band_k, spin) + sym_vec
            end do ! band_k
          end do ! ikcalc
       end do ! spin


       call xmpi_sum(fkn_out, comm, ierr)

       ! Write fkn_out_sigma to disk.
       if (my_rank == master) then
           NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "fkn_out_sigma"), fkn_out, start=[1,1,1,1,itemp]))
       end if

       do spin=1,nsppol
         max_adiff_spin(spin) = maxval(abs(fkn_out(:,:,:,spin) - fkn_in(:,:,:,spin)))
       end do
       max_adiff = maxval(max_adiff_spin)


       mob_gen=sig_gen
       mob_21=sig_gen
       mob_22=sig_gen

!Here I test if F^T=(eps-mu)F^E
!        if (iet==2) then
!         do spin=1,nsppol
!          do ikcalc=1,nkcalc
!            ik_ibz = ibte%kcalc2ebands(1, ikcalc)
!            do band_k=ibte%bstart_ks(ikcalc, spin), ibte%bstop_ks(ikcalc, spin)
!               e_nk = ebands%eig(band_k, ik_ibz, spin)
              !Which one of the following is correct ?
              ! fkn_in(:, ik_ibz, band_k, spin)  = fkn_serta(:, ik_ibz, band_k, spin) * (e_nk-mu_e)/kT
               !fkn_in(:, ik_ibz, band_k, spin)  = fkn_efield(:, ik_ibz, band_k, spin) * (e_nk-mu_e)/kT
!               fkn_out(:, ik_ibz, band_k, spin)  = fkn_efield(:, ik_ibz, band_k, spin) * (e_nk-mu_e)
               !fkn_serta(:, ik_ibz, band_k, spin)  = fkn_serta(:, ik_ibz, band_k, spin) * (e_nk-mu_e)
!            end do
!          end do
!        end do
!       end if


       ! Compute transport tensors from fkn_out (= F_eps or F_T)
       if(iet==1) call ibte_calc_tensors(ibte, cryst, itemp, kT, mu_e, fkn_out, onsager, sig_p, mob_p, fsum_eh, comm, iet)
! call flush_unit(std_out)

      !For seebeck, use another array(sig_gen) because sig_p points towards ibte_sigma
      if(iet==2) call ibte_calc_tensors(ibte, cryst, itemp, kT, mu_e, fkn_out, onsager, sig_gen, mob_gen, fsum_eh, comm, iet)

       !check to find the bug
! call wrtout(std_out," check 3", pre_newlines=1, newlines=1)
! call flush_unit(std_out)

       ! Print mobility for semiconductors or conductivity for metals.
       if (iet==1) then! for F_eps and charge transport
         !sig_p = fact_sigma * sig_gen  ! siemens cm^-1
         !mob_p = fact_mob   * mob_gen  ! cm2/V/s

!call wrtout(std_out, "check 3bis: sig_gen, fact_sigma, sig_p")
 !               write(std_out,*) , sig_gen, fact_sigma, sig_p
          !check to find the bug
! call wrtout(std_out," check 4", pre_newlines=1, newlines=1)

         if (ibte%assume_gap) then
           do spin=1,nsppol
             mat33 = sum(mob_p(:,:,:,spin), dim=3)
             write(msg, "(i5,1x,es9.1,*(1x, f16.2))") &
               iter, max_adiff_spin(spin), mat33(1,1), mat33(2,2), mat33(3,3), sum(fsum_eh(:,:,spin))
           end do
         else
           do spin=1,nsppol
             mat33 = sum(sig_p(:,:,:,spin), dim=3)
             write(msg, "(i5,1x,es9.1,*(1x, es16.6))") &
               iter, max_adiff_spin(spin), mat33(1,1), mat33(2,2), mat33(3,3), sum(fsum_eh(:,:,spin))
           end do
         end if
         call wrtout(std_out, msg)

       else if (iet == 2) then! for Seebeck, Seebeck=1/T sig_p^-1 * sig_gen
         ! output sig_p^-1 * sig_gen
         !pay attention when implementing if sigma=zero, bug ! impossible to get S.
         ! use matr3inv to inverse the matrix sigma but verify also that the det is not equal to zero !
         ! maybe type the command use "..." to be able to use the functions, check in the code

         !call wrtout(std_out, "check : sig_gen, sig_p")
         !write(std_out,*) sig_gen, sig_p

         !TODO: implement for inv_sig_p with two spins
         do spin=1,nsppol
           call inv33(sum(sig_p(:,:,:,spin),dim=3), inv_sig_p)
           !We divide by 100 because sig_p is in Siemens cm^-1 and we want to retrieve meters
          ! work33 = matmul (inv_sig_p, sum(sig_gen(:,:,:,spin),dim=4))
           mat33 = matmul (inv_sig_p, sum(sig_gen(:,:,:,spin),dim=3)) / (kT)
          ! mat33 = (volt_SI / (ibte%kTmesh(itemp) / kb_HaK)) * matmul (inv_sig_p, sum(sig_gen(:,:,:,spin),dim=3))
           write(msg, "(i5,1x,es9.1,*(1x, es16.6))") &
               iter, max_adiff_spin(spin), mat33(1,1), mat33(2,2), mat33(3,3), sum(fsum_eh(:,:,spin))
           sbk_p (:,:,spin) = mat33
         end do
         call wrtout(std_out, msg)
       end if ! iet == 1, 2

       if (ibte%assume_gap) then
         write(msg, "(a5,1x,a9,*(1x, a16))")" ITER", "max_adiff", "mobility_e+h", "sum_k(df_k)"
       else
         write(msg, "(a5,1x,a9,*(1x, a16))")" ITER", "max_adiff", "conductivity", "sum_k(df_k)"
       end if
       call wrtout(std_out, msg)


       ! Check for convergence by testing max_k |F_k^i - F_k^{i-1}|.
       call wrtout (std_out, "Print check convergence: max_adiff, abs_tol")
       write(std_out,*) max_adiff, abs_tol
       converged(itemp) = max_adiff < abs_tol
       if (converged(itemp)) then
         call wrtout(std_out, sjoin(" IBTE solver converged after:", itoa(iter), &
                     "iterations within ibte_abs_tol:", ftoa(abs_tol)), pre_newlines=1)

         exit iter_loop
       else
         ! Linear mixing of fkn_in and fkn_out.
         fkn_in = (one - dtset%ibte_alpha_mix) * fkn_in + dtset%ibte_alpha_mix * fkn_out
         fkn_out = zero
       end if
        !check where is the problem with iter loop

        !call wrtout(std_out," check iet:", pre_newlines=1, newlines=1)
        !write (std_out,*) iet

     end do iter_loop

     if (.not. converged(itemp)) then
       msg = sjoin("Not converged after:", itoa(dtset%ibte_niter), "max iterations")
       call wrtout(ab_out, msg, pre_newlines=1, newlines=1)
       ABI_WARNING(msg)
     end if

     !Change of unit here because sig_p is used in to calculate sbk in SI units
        if (iet == 1) fkn_efield = fkn_out
        if (iet==2) then
                !iet =3
                call ibte_calc_tensors(ibte, cryst, itemp, kT, mu_e, fkn_efield, onsager, sig_l21, mob_21, fsum_eh, comm, 3)
                call ibte_calc_tensors(ibte, cryst, itemp, kT, mu_e, fkn_out, onsager, sig_l22, mob_22, fsum_eh, comm, 3)
                do spin=1, nsppol
                call inv33(sum(sig_p(:,:,:,spin),dim=3), inv_sig_p)
                kappa_p(:,:,spin)=(1/kT)*(sum(sig_l22(:,:,:,spin), dim=3)-matmul(sum(sig_l21(:,:,:,spin),dim=3),sbk_p(:,:,spin)*kT))
                pi_p(:,:,spin)= matmul (sum(sig_l21(:,:,:,spin),dim=3), inv_sig_p) ! PI=L21/L11 in IBTE
                end do
        end if
   end do electherm_loop
! end E/T loop

!Units conversion for sigma, at this step in order to keep atomic units for the calculation of S
!Truncated in purpose for nspinor=2 and nsppol=1 compared to how its done in the calc_tensor routine because dont have access to self in this routine
sig_p=sig_p*(siemens_SI / Bohr_meter / cryst%ucvol) / 100
!check if sig_p and sbk_p are empty
!call wrtout(std_out, "sig_p, sbk_p")
!write(std_out,*)  sig_p, sbk_p

 end do ! itemp

 !TODO: do it for spins as well
 ABI_MALLOC(ibte_rho, (3, 3, ntemp))
 do itemp=1,ntemp
   work33 = sum(ibte_sigma(:,:,:,1,itemp), dim=3)
   if (ibte%nsppol == 2) work33 = work33 + sum(ibte_sigma(:,:,:,2,itemp), dim=3)
   call inv33(work33, mat33)
   ibte_rho(:, :, itemp) = 1e+6_dp * mat33
 end do

!ibte_sigma rho and seebeck seem empty, let's check


!write(*,*) "sigma: ", ibte_sigma
!write(*,*) "rho: ", ibte_rho
!write(*,*) "seebeck: ", ibte_seebeck


 if (my_rank == master) then
   ! Write final results to main output.
   components = ["xx", "yy", "zz"]
   if (ibte%assume_gap) then
     ! SemiConductor
     do ii=1,3
       call wrtout(units, sjoin(" Cartesian component of IBTE mobility tensor:", components(ii)))
       write(msg, "(a16,2(a32),a16)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]', "Converged"
       call wrtout(units, msg)

       do spin=1,ibte%nsppol
         if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)

         do itemp=1,ibte%ntemp
           write(msg,"(f16.2,2e16.2,2f16.2,a16)") &
             ibte%kTmesh(itemp) / kb_HaK, &
             ibte%n_ehst(1, spin, itemp) / cryst%ucvol / Bohr_cm **3, &
             ibte%n_ehst(2, spin, itemp) / cryst%ucvol / Bohr_cm **3, &
             ibte_mob(ii, ii, 1, spin, itemp), ibte_mob(ii, ii, 2, spin, itemp), &
             yesno(converged(itemp))
           call wrtout(units, msg)
         end do ! itemp
       end do ! spin
       ! TODO:HERE look into adding seebeck output for semiconductors as well, check what happens with e and h components in tensors
       ! routine
       call wrtout(units, ch10)
     end do ! ii

   else
     ! Metals. Print conductivity (spin resolved) and resistivity (no spin resolved)
     do ii=1,2
       if (ii == 1) msg = " Conductivity [Siemens cm^-1] using IBTE"
       if (ii == 2) msg = " Resistivity [micro-Ohm cm] using IBTE"
       call wrtout(units, msg)

       nsp = ibte%nsppol; if (ii == 2) nsp = 1
       do spin=1,nsp
         if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
         write(msg, "(5a16)") 'Temperature (K)', 'xx', 'yy', 'zz', "Converged"
         call wrtout(units, msg)
         do itemp=1,ibte%ntemp
           if (ii == 1) then
             mat33 = sum(ibte_sigma(:,:,:,spin,itemp), dim=3)
           else
             mat33 = ibte_rho(:,:,itemp)
           end if
           write(msg,"(f16.2,3e16.6,a16)") &
             ibte%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3), yesno(converged(itemp))
           call wrtout(units, msg)
         end do ! itemp
       end do ! spin
       call wrtout(units, ch10)
     end do ! ii

     ! HERE add output of Seebeck and kappa_el coefficients
     ! TODO: add off diagonal Seebeck and sigma coefficients
     msg = " Seebeck [Volts / Kelvin] using IBTE"
     call wrtout(units, msg)
     do spin=1,ibte%nsppol
       if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(5a16)") 'Temperature (K)', 'xx', 'yy', 'zz', "Converged"
       call wrtout(units, msg)
       do itemp=1,ibte%ntemp
         mat33 = fact_sbk * ibte_seebeck(:,:,spin,itemp)
         write(msg,"(f16.2,3e16.6,a16)") &
           ibte%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3), yesno(converged(itemp))
         call wrtout(units, msg)
       end do ! itemp
     end do ! spin
     call wrtout(units, ch10)


     msg = " Kappa [W/m*K] using IBTE"
     call wrtout(units, msg)
     do spin=1,ibte%nsppol
       if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(5a16)") 'Temperature (K)', 'xx', 'yy', 'zz', "Converged"
       call wrtout(units, msg)
       do itemp=1,ibte%ntemp
         mat33 = volt_SI**2 * kb_HaK * (siemens_SI / Bohr_meter / cryst%ucvol)* ibte_kappa(:,:,spin,itemp)
         write(msg,"(f16.2,3e16.6,a16)") &
           ibte%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3), yesno(converged(itemp))
         call wrtout(units, msg)
       end do ! itemp
     end do ! spin
     call wrtout(units, ch10)

     msg = " Peltier [Volts] using IBTE"
     call wrtout(units, msg)
     do spin=1,ibte%nsppol
       if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)
       write(msg, "(5a16)") 'Temperature (K)', 'xx', 'yy', 'zz', "Converged"
       call wrtout(units, msg)
       do itemp=1,ibte%ntemp
         mat33 = -volt_SI * ibte_pi(:,:,spin,itemp)
         write(msg,"(f16.2,3e16.6,a16)") &
           ibte%kTmesh(itemp) / kb_HaK, mat33(1,1), mat33(2,2), mat33(3,3), yesno(converged(itemp))
         call wrtout(units, msg)
       end do ! itemp
     end do ! spin
     call wrtout(units, ch10)

     msg = "Carrier density: "
     call wrtout(units, msg)
       write(msg, "(a16,a16,a32)") 'Temperature [K]', ' chem pot [eV]  ', 'e/h density [cm^-3]'
       call wrtout(units, msg)

       do spin=1,ibte%nsppol
         if (ibte%nsppol == 2) call wrtout(units, sjoin(" For spin:", stoa(spin)), newlines=1)

         do itemp=1,ibte%ntemp
           write(msg,"(f16.2,3e16.2)") &
             ibte%kTmesh(itemp) / kb_HaK, &
             ibte%eph_mu_e(itemp) * eV_Ha, &
             ibte%n_ehst(1, spin, itemp) / cryst%ucvol / Bohr_cm **3, &
             ibte%n_ehst(2, spin, itemp) / cryst%ucvol / Bohr_cm **3
           call wrtout(units, msg)
         end do ! itemp
       end do ! spin
       call wrtout(units, ch10)




   end if

   !pre = "_IBTE"
   !call ibte%write_tensor(dtset, irta, "sigma", ibte%sigma(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SIGMA"))
   !call ibte%write_tensor(dtset, irta, "seebeck", ibte%seebeck(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SBK"))
   !call ibte%write_tensor(dtset, irta, "kappa", ibte%kappa(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_KAPPA"))
   !call ibte%write_tensor(dtset, irta, "zte", ibte%zte(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_ZTE"))
   !call ibte%write_tensor(dtset, irta, "pi", ibte%pi(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_PI"))

   ! Print IBTE results to stdout and other external txt files (for the test suite)
   !call ibte%print_rta_txt_files(cryst, dtset, dtfil)
   ! Creates the netcdf file used to store the results of the calculation.
   path = strcat(dtfil%filnam_ds(4), "_RTA.nc")
   call wrtout(units, ch10//sjoin("- Writing IBTE transport results to:", path))
   NCF_CHECK(nctk_open_modify(ncid, path , xmpi_comm_self))

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t('ibte_sigma', "dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('ibte_mob', "dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('ibte_rho', "dp", "three, three, ntemp"), &
     nctkarr_t('ibte_kappa', "dp", "three, three, nsppol, ntemp"), &
     nctkarr_t('ibte_seebeck', "dp", "three, three, nsppol, ntemp"), &
     nctkarr_t('ibte_pi', "dp", "three, three, nsppol, ntemp") &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_sigma"), ibte_sigma))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_mob"), ibte_mob))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_rho"), ibte_rho))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_seebeck"), ibte_seebeck))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_pi"), ibte_pi))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ibte_kappa"), ibte_kappa))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ibz"), ibte%kcalc2ibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ebands"), ibte%kcalc2ebands))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kibz"), ibte%ebands%kptns))
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Free memory
 ABI_FREE(fkn_serta)
 ABI_FREE(taukn_serta)
 ABI_FREE(fkn_in)
 ABI_FREE(fkn_out)
 ABI_FREE(ibte_sigma)
 ABI_FREE(ibte_seebeck)
 ABI_FREE(ibte_pi)
 ABI_FREE(ibte_mob)
 ABI_FREE(ibte_rho)
 ABI_FREE(ibte_kappa)
 ABI_FREE(converged)
 ABI_FREE(sig_gen)
 ABI_FREE(sig_l21)
 ABI_FREE(sig_l22)
 ABI_FREE(mob_gen)
 ABI_FREE(mob_21)
 ABI_FREE(mob_22)

 do spin=1,nsppol
   do ikcalc=1,nkcalc
     call free_sr_ks(ikcalc, spin)
   end do
 end do
 ABI_FREE(sr)

 call ibte%free()

contains

subroutine free_sr_ks(ikc, isp)
  integer,intent(in) :: ikc, isp
  ABI_SFREE(sr(ikc, isp)%vals)
  ABI_SFREE(sr(ikc, isp)%lgk_sym2glob)
  ABI_SFREE(sr(ikc, isp)%kq_symtab)
end subroutine free_sr_ks

end subroutine ibte_driver
!!***

!----------------------------------------------------------------------

!!****f* m_rta/ibte_calc_tensors
!! NAME
!! ibte_calc_tensors
!!
!! FUNCTION
!!   calculate transport tensors within iBTE method, from v x F expressions for current and sigma etc
!!   this routine now stays in atomic units to accommodate calculation of sigma and seebeck
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! comm=MPI communicator.
!!
!! SOURCE

subroutine ibte_calc_tensors(self, cryst, itemp, kT, mu_e, fk, onsager, sigma_eh, mob_eh, fsum_eh, comm, iet)

!Arguments ------------------------------------
 class(rta_t),intent(inout) :: self
 type(crystal_t),intent(in) :: cryst
 integer,intent(in) :: itemp
 real(dp),intent(in) :: kT, mu_e
 real(dp),intent(in) :: fk(3, self%ebands%nkpt, self%bmin:self%bmax, self%nsppol)
 real(dp),intent(out) :: sigma_eh(3,3,2,self%nsppol), mob_eh(3,3,2,self%nsppol)
 real(dp),intent(out) :: fsum_eh(3,2,self%nsppol), onsager(3,3,3,self%nsppol)
 integer,intent(in) :: comm, iet !iet to know which tensor we calculate, 1 for L11, 2 for L12 and 3 for L21 and L22

!Local variables ------------------------------
!scalars
 integer :: nsppol, nkibz, ib, ik_ibz, spin, ii, jj, ieh, cnt, nprocs, ia, time_opt
 real(dp) :: eig_nk, max_occ, wtk, emu_alpha, fact, fact0
 !real(dp) :: fact_sigma, fact_mob
 !arrays
 real(dp) :: vr(3), vv_tens(3,3)
!************************************************************************

 ABI_UNUSED(kt)

 nprocs = xmpi_comm_size(comm)

 ! Copy important dimensions
 nkibz = self%ebands%nkpt; nsppol = self%ebands%nsppol
 time_opt = 0 ! This to preserve the previous behaviour in which TR was not used.

 ! sigma_IBTE = (-S e^ / omega sum_\nk) (v_\nk \otimes F_\nk)
 ! with S the spin degeneracy factor.
 sigma_eh = zero; fsum_eh = zero; onsager = zero

 ! Compute mobility_mu i.e. results in which lifetimes have been computed in a consistent way
 ! with the same the Fermi level. In all the other cases, indeed, we assume that tau does not depend on ef.
 !

 !TODO: rewrite this because every L (onsager coeff) is given with a minus sign and it's just the case for L11 and L22 normally
 ! Fortunately these minus signs compensate each other in the code or are suppressed by putting a -1 factor in front of sbk and PI
 ! It is not because of the sign of e !
 cnt = 0
 do spin=1,nsppol
   do ik_ibz=1,nkibz
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     wtk = self%ebands%wtk(ik_ibz)

     do ib=self%bmin,self%bmax
       eig_nk = self%ebands%eig(ib, ik_ibz, spin)

       ! Compute outer product in vv_tens and symmetrize tensor.
       vr(:) = self%vbks(:, ib, ik_ibz, spin)
       do ia=1,3
         if (ia == 1) then
           emu_alpha = one
         else
           emu_alpha = (eig_nk - mu_e) ** (ia - 1)
         end if

         do ii=1,3
           do jj=1,3
             vv_tens(ii, jj) = vr(ii) * fk(jj, ik_ibz, ib, spin)
           end do
         end do
         vv_tens = vv_tens * emu_alpha
         vv_tens = cryst%symmetrize_cart_tens33(vv_tens, time_opt)
        ! TODO: write the expression below in another way, without the useless sigma_eh - and without the sign - (adapt the
        ! calculations of coefficients in the main routine then)
         if (ia == 1 .and. iet /= 3) then
           ieh = 2; if (eig_nk >= mu_e) ieh = 1
           sigma_eh(:,:,ieh,spin) = sigma_eh(:,:,ieh,spin) - wtk * vv_tens
           fsum_eh(:,ieh,spin) = fsum_eh(:,ieh,spin) + wtk * cryst%symmetrize_cart_vec3(fk(:, ik_ibz, ib, spin), time_opt)

         !Here the tensorial product in order to get L21 and L22 to calculate kappa
         else if (ia==2 .and. iet==3) then
                 ieh = 2; if (eig_nk >= mu_e) ieh = 1
                 sigma_eh(:,:,ieh,spin) = sigma_eh(:,:,ieh,spin) - wtk * vv_tens
           fsum_eh(:,ieh,spin) = fsum_eh(:,ieh,spin) + wtk * cryst%symmetrize_cart_vec3(fk(:, ik_ibz, ib, spin), time_opt)
         end if



         onsager(:,:,ia,spin) = onsager(:,:,ia,spin) - wtk * vv_tens
       end do ! ia

     end do ! ib
   end do ! ik_ibz
 end do ! spin

 if (iet==1) then ! In order to output the mobility (only for L11, the conductivity, thus iet==1)
   max_occ = two / (self%nspinor * self%nsppol)
   fact0 = max_occ * (siemens_SI / Bohr_meter / cryst%ucvol) / 100
   fact = 100**3 / e_Cb

   sigma_eh = fact0 * sigma_eh  ! siemens cm^-1
   fsum_eh = fsum_eh / cryst%ucvol

 ! Scale by the carrier concentration.
   do spin=1,nsppol
     do ieh=1,2
       call safe_div(sigma_eh(:,:,ieh,spin) * fact, &
                     self%n_ehst(ieh, spin, itemp) / cryst%ucvol / Bohr_meter**3, zero, mob_eh(:,:,ieh,spin))
     end do
   end do

   !Here I rescale sigma_eh to output correctly the Onsager coeff L11
   sigma_eh = sigma_eh / fact0
   end if
 !call xmpi_sum(sigma_eh, comm, ierr)
 !call xmpi_sum(onsager, comm, ierr)

! max_occ = two / (self%nspinor * self%nsppol)
! sigma_eh = max_occ * sigma_eh / cryst%ucvol
 !fsum_eh = fsum_eh / cryst%ucvol

 max_occ = two / (self%nspinor * self%nsppol)
 !Take into account spin degeneracy for all Onsager coefficient:
 sigma_eh= max_occ * sigma_eh



! fact0 = max_occ * (siemens_SI / Bohr_meter / cryst%ucvol) / 100
! fact = 100**3 / e_Cb
! sigma_eh = fact0 * sigma_eh



 ! Scale by the carrier concentration.
! do spin=1,nsppol
!   do ieh=1,2
!     call safe_div(sigma_eh(:,:,ieh,spin) * fact, &
!                   self%n_ehst(ieh, spin, itemp) / cryst%ucvol / Bohr_meter**3, zero, mob_eh(:,:,ieh,spin))

!check to find the bug
!call wrtout(std_out, " ieh, shape(sigma_eh) vol shape(mob_eh) ", pre_newlines=1, newlines=1)
!write (std_out,*)  ieh, shape(sigma_eh), cryst%ucvol, shape(mob_eh)
!call wrtout(std_out," check safe_div, in order: sigma_eh, self%n_ehst, cryst%ucvol, zero", pre_newlines=1, newlines=1)
!        write (std_out,*) sigma_eh(:,:,ieh,spin), self%n_ehst(ieh, spin, itemp), cryst%ucvol, zero
! call flush_unit(std_out)

!        call safe_div(sigma_eh(:,:,ieh,spin), &
!                   self%n_ehst(ieh, spin, itemp) / cryst%ucvol, zero, mob_eh(:,:,ieh,spin))
!mob_eh(:,:,ieh,spin) = sigma_eh(:,:,ieh,spin)   / (self%n_ehst(ieh, spin, itemp) / cryst%ucvol)
           !check to find the bug
!           call wrtout(std_out," check safe_div, mob_eh: ", pre_newlines=1, newlines=1)
!        write (std_out,*) mob_eh(:,:,ieh,spin)
! call flush_unit(std_out)

!   end do
! end do
!In order to be able to compile, to have a value for the dummy argument mob_eh

!call wrtout(std_out, "shape(sigma_eh), sigma_eh,  shape(mob_eh), mob_eh ", pre_newlines=1, newlines=1)
!write (std_out,*) shape(sigma_eh),sigma_eh, shape(mob_eh), mob_eh

!mob_eh=zero

!mob_eh=sigma_eh


!max_occ = two / (self%nspinor * self%nsppol)
! fact_sigma = max_occ * (siemens_SI / Bohr_meter) / 100.  ! Siemens / cm
! fact_mob   = fact_sigma * 100.**3 / e_Cb * Bohr_meter**3 ! cm^2 / V / s
! fact_sbk   = volt_SI * kb_HaK                           ! Volt / Kelvin
!   sigma_eh = fact_sigma * sigma_eh  ! siemens cm^-1
 !  mob_eh = fact_mob   * mob_eh  !
!call wrtout(std_out, "factors: fact_sigma, fact_mob", pre_newlines=1, newlines=1)
!   write (std_out,*) fact_sigma, fact_mob

end subroutine ibte_calc_tensors
!!***

! Invert 3x3 matrix, copied from matr3inv
pure subroutine inv33(aa, ait)

!Arguments ------------------------------------
!arrays
 real(dp),intent(in) :: aa(3,3)
 real(dp),intent(out) :: ait(3,3)

!Local variables-------------------------------
!scalars
 real(dp) :: dd,det,t1,t2,t3

! *************************************************************************

 t1 = aa(2,2) * aa(3,3) - aa(3,2) * aa(2,3)
 t2 = aa(3,2) * aa(1,3) - aa(1,2) * aa(3,3)
 t3 = aa(1,2) * aa(2,3) - aa(2,2) * aa(1,3)
 det = aa(1,1) * t1 + aa(2,1) * t2 + aa(3,1) * t3

 ! Make sure matrix is not singular
 if (abs(det) > 100 * tiny(one)) then
   dd = one / det
   ait(1,1) = t1 * dd
   ait(2,1) = t2 * dd
   ait(3,1) = t3 * dd
   ait(1,2) = (aa(3,1)*aa(2,3)-aa(2,1)*aa(3,3)) * dd
   ait(2,2) = (aa(1,1)*aa(3,3)-aa(3,1)*aa(1,3)) * dd
   ait(3,2) = (aa(2,1)*aa(1,3)-aa(1,1)*aa(2,3)) * dd
   ait(1,3) = (aa(2,1)*aa(3,2)-aa(3,1)*aa(2,2)) * dd
   ait(2,3) = (aa(3,1)*aa(1,2)-aa(1,1)*aa(3,2)) * dd
   ait(3,3) = (aa(1,1)*aa(2,2)-aa(2,1)*aa(1,2)) * dd
   ait = transpose(ait)
 else
   ait = zero
 end if

 end subroutine inv33

end module m_rta
!!***

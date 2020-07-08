!!****m* ABINIT/m_rta
!! NAME
!!  m_rta
!!
!! FUNCTION
!!  Module to compute transport properties using the
!!  Boltzmann transport equation (BTE) in the relaxation-time approximation (RTA).
!!  Initially for phonon-limited carrier mobility.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (HM, MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
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
 use m_ebands
 use m_nctk
 use m_wfk
 use m_ephtk
 use m_sigmaph
 use m_dtset
 use m_dtfil
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : pseudopotential_type, ebands_t
 use m_io_tools,       only : open_file
 use m_time,           only : cwtime, cwtime_report
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : bisect, simpson_int, safe_div
 use m_fstrings,       only : strcat, sjoin, itoa, ltoa, stoa, ftoa
 use m_kpts,           only : listkk, kpts_timrev_from_kptopt
 use m_occ,            only : occ_fd, occ_dfde
 use m_pawtab,         only : pawtab_type
 use m_ddk,            only : ddkstore_t

 implicit none

 private
!!****

 public :: rta_driver ! Main entry point for transport calculations withing RTA
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
   ! number of independent spin polarizations.

   integer :: nspinor
   ! number of spinorial components.

   integer :: ntemp
   ! number of temperatures.

   integer :: nw
   ! number of energies (chemical potentials) at which transport quantities are computed
   ! Same number of energies used in DOS.

   integer :: bmin, bmax, bsize

   integer :: nrta
   ! Number of relaxation-time approximations used (1 for SERTA, 2 for MRTA)

   real(dp) :: eph_extrael
   ! extra electrons per unit cell from sigeph (lifetimes)

   real(dp) :: eph_fermie
   ! Fermi level from input file from sigeph (lifetimes)

   real(dp) :: transport_extrael
   ! extra electrons per unit cell

   real(dp) :: transport_fermie
   ! Fermi level from input file

   logical :: assume_gap

   real(dp),allocatable :: kTmesh(:)
   ! (%ntemp)
   ! a list of temperatures at which to compute the transport

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
    ! Linewidth computed in the SERTA/MRTA

   real(dp),allocatable :: velocity(:,:,:,:)
   ! (3, bmin:bmax, nkpt, nsppol))
   ! band velocity in Cartesian coordinates.

   type(gaps_t) :: gaps
   ! gaps of original ebands object

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

   real(dp),allocatable :: ne(:)
   ! (%ntemp) number of electrons at transport_mu_e(ntemp)

   real(dp),allocatable :: nh(:)
   ! (%ntemp) number of holes at transport_mu_e(ntemp)

   real(dp),allocatable :: l0(:,:,:,:,:,:)
   real(dp),allocatable :: l1(:,:,:,:,:,:)
   real(dp),allocatable :: l2(:,:,:,:,:,:)
   ! (3, 3, nw, ntemp, nsppol, nrta)
   ! Onsager coeficients

   real(dp),allocatable :: sigma(:,:,:,:,:,:)
   real(dp),allocatable :: seebeck(:,:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:,:)
   real(dp),allocatable :: pi(:,:,:,:,:,:)
   real(dp),allocatable :: zte(:,:,:,:,:,:)
   ! (3, 3, nw, ntemp, nsppol, nrta)
   ! Transport coefficients

   real(dp),allocatable :: mobility(:,:,:,:,:,:,:)
   ! Mobility
   ! (3, 3, nw, %ntemp, 2, %nsppol, nrta)
   ! 5-th index is for e-h

   real(dp),allocatable :: n(:,:,:)
   ! (nw, ntemp, 2) carrier density for e/h (n/cm^3)

   real(dp),allocatable :: mobility_mu(:,:,:,:,:,:)
   ! (3, 3, 2, ntemp, nsppol, nrta)
   ! mobility for electrons and holes (third dimension) at transport_mu_e(ntemp)
   ! Third index is for electron/hole

 contains

    procedure :: compute => rta_compute
    procedure :: compute_mobility => rta_compute_mobility
    procedure :: print_txt_files => rta_print_txt_files
    procedure :: write_tensor => rta_write_tensor
    procedure :: free => rta_free

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
!! General driver to compute transport properties
!!
!! INPUTS
!! dtfil<datafiles_type>=variables related to files.
!! ngfftc(18)=Coarse FFT meshe
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
 integer,parameter :: master = 0
 integer :: ierr, my_rank
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=500) :: msg
 character(len=fnlen) :: path
 type(sigmaph_t) :: sigmaph
 type(rta_t) :: rta
!arrays
 integer :: unts(2)
 real(dp) :: extrael_fermie(2)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)
 unts = [std_out, ab_out]

 call wrtout(unts, ch10//' Entering transport RTA computation driver.')
 call wrtout(unts, sjoin("- Reading carrier lifetimes from:", dtfil%filsigephin), newlines=1, do_flush=.True.)

 sigmaph = sigmaph_read(dtfil%filsigephin, dtset, xmpi_comm_self, msg, ierr, keep_open=.true., extrael_fermie=extrael_fermie)
 ABI_CHECK(ierr == 0, msg)

 ! Initialize RTA object
 ! TODO: Should store more metadata: energy window, nkcalc ....
 rta = rta_new(dtset, dtfil, ngfftc, sigmaph, cryst, ebands, pawtab, psps, extrael_fermie, comm)

 ! sigmaph is not needed anymore. Free it.
 sigmaph%ncid = nctk_noid
 call sigmaph%free()

 ! Compute RTA transport quantities
 call rta%compute(cryst, dtset, comm)

 ! Compute RTA mobility
 call rta%compute_mobility(cryst, dtset, comm)

 if (my_rank == master) then
   ! Print RTA results to stdout and other external txt files (for the test suite)
   call rta%print_txt_files(cryst, dtset, dtfil)

   ! Creates the netcdf file used to store the results of the calculation.
#ifdef HAVE_NETCDF
   path = strcat(dtfil%filnam_ds(4), "_RTA.nc")
   call wrtout(unts, ch10//sjoin("- Writing RTA transport results to:", path))
   NCF_CHECK(nctk_open_create(ncid, path , xmpi_comm_self))
   call rta_ncwrite(rta, cryst, dtset, ncid)
   NCF_CHECK(nf90_close(ncid))
 end if
#endif

 ! Free memory
 call rta%free()

end subroutine rta_driver
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_new
!! NAME
!! rta_new
!!
!! FUNCTION
!! Compute transport quantities in the relaxation time approximation
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  dtfil<datafiles_type>=variables related to files.
!!  sigmaph<sigmaph_t>=Object with e-ph self-energy results.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  extrael_fermie
!!  comm=MPI communicator.
!!
!! SOURCE

type(rta_t) function rta_new(dtset, dtfil, ngfftc, sigmaph, cryst, ebands, pawtab, psps, extrael_fermie, comm) result (new)

!Arguments -------------------------------------
 integer, intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(sigmaph_t),intent(in) :: sigmaph
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfftc(18)
 real(dp),intent(in) :: extrael_fermie(2)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
 integer,parameter :: sppoldbl1 = 1, master = 0
 integer :: ierr, spin, nprocs, my_rank, timrev, ik_ibz, ib, irta, itemp, ndat, nsppol, idat, mband, ikpt
 real(dp) :: dksqmax, cpu, wall, gflops
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname_dense
 type(ebands_t) :: tmp_ebands, ebands_dense
 type(klinterp_t) :: klinterp
 type(ddkstore_t) :: ds
!arrays
 integer :: kptrlatt(3,3), unts(2)
 integer,allocatable :: indkk(:,:)
 real(dp),allocatable :: values_bksd(:,:,:,:), vals_bsd(:,:,:)
 real(dp),allocatable :: tmp_array4(:,:,:,:), tmp_array5(:,:,:,:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 unts = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")

 new%assume_gap = (.not. all(dtset%sigma_erange < zero) .or. dtset%gw_qprange /= 0)

 ! Allocate temperature arrays (same as the ones used in the SIGEPH calculation).
 new%ntemp = sigmaph%ntemp
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 new%kTmesh = sigmaph%kTmesh

 ! How many RTA approximations have we computed in sigmaph? (SERTA, MRTA)
 new%nrta = 2; if (sigmaph%mrta == 0) new%nrta = 1

 new%bmin = minval(sigmaph%bstart_ks); new%bmax = maxval(sigmaph%bstop_ks)
 !new%bmin = 1; new%bmax = ebands%mband  ! This for debugging purposes
 new%bsize = new%bmax - new%bmin + 1

 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor
 nsppol = new%nsppol

 ABI_MALLOC(new%eminmax_spin, (2, nsppol))
 new%eminmax_spin = ebands_get_minmax(ebands, "eig")

 if (new%assume_gap) then
   ! Information about the gaps
   new%gaps = ebands_get_gaps(ebands, ierr)
   if (ierr /= 0) then
     do spin=1, nsppol
       MSG_WARNING(trim(new%gaps%errmsg_spin(spin)))
       new%gaps%vb_max(spin) = ebands%fermie - 1 * eV_Ha
       new%gaps%cb_min(spin) = ebands%fermie + 1 * eV_Ha
     end do
     MSG_WARNING("ebands_get_gaps returned non-zero exit status. See above warning messages...")
   end if

   if (my_rank == master) then
     call new%gaps%print(unit=std_out); call new%gaps%print(unit=ab_out)
   end if
 end if

 ! =================================================
 ! Read lifetimes and new%ebands from SIGEPH.nc file
 ! =================================================
 ! After this point we have:
 !
 !      velocity(3, bmin:bmax, nkpt, nsppol)
 !      linewidths(self%ntemp, bmin:bmax, nkpt, nsppol, 2)
 !
 if (any(dtset%sigma_ngkpt /= 0)) then
   ! If integrals are computed with the sigma_ngkpt k-mesh, we need to downsample ebands.
   call wrtout(unts, sjoin(" Computing integrals with downsampled sigma_ngkpt:", ltoa(dtset%sigma_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)

   tmp_ebands = ebands_downsample(ebands, cryst, kptrlatt, dtset%sigma_nshiftk, dtset%sigma_shiftk)
   new%ebands = sigmaph%get_ebands(cryst, tmp_ebands, [new%bmin, new%bmax], new%linewidths, new%velocity, xmpi_comm_self)
   call ebands_free(tmp_ebands)
 else
   new%ebands = sigmaph%get_ebands(cryst, ebands, [new%bmin, new%bmax], new%linewidths, new%velocity, xmpi_comm_self)
   kptrlatt = ebands%kptrlatt
 end if

 !print *, "linewidth_serta", maxval(abs(new%linewidths(:,:,:,:,1)))
 !print *, "linewidth_mrta", maxval(abs(new%linewidths(:,:,:,:,2)))

 if ( &
     dtset%useria == 888 .and. &
     (dtset%getwfkfine /= 0 .or. dtset%irdwfkfine /= 0 .or. dtset%getwfkfine_filepath /= ABI_NOFILE)) then

   ! In principle only getwfkfine_filepath is used here
   wfk_fname_dense = trim(dtfil%fnameabi_wfkfine)
   ABI_CHECK(nctk_try_fort_or_ncfile(wfk_fname_dense, msg) == 0, msg)

   call wrtout(unts, " EPH double grid interpolation: will read energies from: "//trim(wfk_fname_dense), newlines=1)
   mband = new%ebands%mband

   !ebands_dense = wfk_read_ebands(wfk_fname_dense, comm)

   tmp_ebands = wfk_read_ebands(wfk_fname_dense, comm)
   ebands_dense = ebands_chop(tmp_ebands, 1, mband)
   call ebands_free(tmp_ebands)
   if (my_rank == master) then
     write(std_out, *)" Using kptrlatt: ", ebands_dense%kptrlatt
     write(std_out, *)"       shiftk: ", ebands_dense%shiftk
   end if
   ABI_CHECK_IEQ(mband, ebands_dense%mband, "Inconsistent number of bands for the fine and dense grid:")

   ! Compute v_{nk} on the dense grid in Cartesian coordinates.
   ! vdiago(3, bmin:bmax, nkpt, nsppol)
   ! NB: We select states in [bmin:bmax] but all k-points in the IBZ are computed!
   ds%only_diago = .True.; ds%bmin = new%bmin; ds%bmax = new%bmax; ds%mode = "cart"
   call ds%compute_ddk(wfk_fname_dense, "", dtset, psps, pawtab, ngfftc, comm)

   ! Transfer data to new%velocity
   ABI_MOVE_ALLOC(ds%vdiago, new%velocity)
   call ds%free()
   !print *, "velocity:", new%velocity

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
   call ebands_move_alloc(ebands_dense, new%ebands)

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
     MSG_WARNING(sjoin("Linear interpolation produced:", itoa(ierr), "negative linewidths"))
   end if

   ABI_FREE(vals_bsd)
   call klinterp%free()
 end if

 ! FIXME: I think transport_ngkpt is buggy, wrong ne(T), weird zeros if MRTA ...
 ! Do we really need this option? Can't we replace it with sigma_ngkpt and eph_task 7?

 if (any(dtset%transport_ngkpt /= 0)) then
   ! Perform further downsampling (usefull for debugging purposes)
   call wrtout(unts, " Downsampling the k-mesh before computing transport:")
   call wrtout(unts, sjoin(" Using transport_ngkpt: ", ltoa(dtset%transport_ngkpt)))
   kptrlatt = 0
   kptrlatt(1, 1) = dtset%transport_ngkpt(1)
   kptrlatt(2, 2) = dtset%transport_ngkpt(2)
   kptrlatt(3, 3) = dtset%transport_ngkpt(3)
   tmp_ebands = ebands_downsample(new%ebands, cryst, kptrlatt, 1, [zero, zero, zero])

   ! Map the points of the downsampled bands to dense ebands
   timrev = kpts_timrev_from_kptopt(ebands%kptopt)
   ABI_MALLOC(indkk, (tmp_ebands%nkpt, 6))

   call listkk(dksqmax, cryst%gmet, indkk, new%ebands%kptns, tmp_ebands%kptns, &
               new%ebands%nkpt, tmp_ebands%nkpt, cryst%nsym, &
               sppoldbl1, cryst%symafm, cryst%symrec, timrev, comm, exit_loop=.True., use_symrec=.True.)

   if (dksqmax > tol12) then
      write(msg, '(3a,es16.6,a)' ) &
       "Error while downsampling ebands in the transport driver",ch10, &
       "The k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10
      MSG_ERROR(msg)
   end if

   ! Downsampling linewidths and velocities.
   ABI_MOVE_ALLOC(new%linewidths, tmp_array5)
   ABI_REMALLOC(new%linewidths, (new%ntemp, new%bmin:new%bmax, tmp_ebands%nkpt, nsppol, new%nrta))
   do ikpt=1,tmp_ebands%nkpt
     new%linewidths(:,:,ikpt,:,:) = tmp_array5(:,:,indkk(ikpt, 1),:,:)
   end do
   ABI_FREE(tmp_array5)

   ABI_MOVE_ALLOC(new%velocity, tmp_array4)
   ABI_REMALLOC(new%velocity, (3, new%bmin:new%bmax, tmp_ebands%nkpt, nsppol))
   do ikpt=1,tmp_ebands%nkpt
     new%velocity(:,:,ikpt,:) = tmp_array4(:,:,indkk(ikpt, 1),:)
   end do
   ABI_FREE(tmp_array4)

   !call downsample_array5(new%linewidths, indkk, tmp_ebands%nkpt)
   !call downsample_array4(new%velocity, indkk, tmp_ebands%nkpt)

   !print *, "after downsampling linewidths"
   !print *, "linewidth_serta", maxval(abs(new%linewidths(:,:,:,:,1)))
   !print *, "linewidth_mrta", maxval(abs(new%linewidths(:,:,:,:,2)))

   ABI_FREE(indkk)
   call ebands_move_alloc(tmp_ebands, new%ebands)
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

   ! Compute Fermi level at various T
   call ebands_get_muT_with_fd(ebands, new%ntemp, new%kTmesh, dtset%spinmagntarget, dtset%prtvol, new%transport_mu_e, comm)
 end if

 call cwtime_report(" rta_new", cpu, wall, gflops)

! contains
!
! subroutine downsample_array4(array, indkk, nkpt)
!
!   ! (ntemp, max_nbcalc, nkcalc, nsppol)
!   real(dp),allocatable,intent(inout) :: array(:,:,:,:)
!   integer,intent(in) :: indkk(:,:)
!
!   integer :: ikpt, nkpt
!   integer :: tmp_shape(4)
!   real(dp),allocatable :: tmp_array(:,:,:,:)
!
!   ABI_MOVE_ALLOC(array, tmp_array)
!   tmp_shape = shape(array)
!   tmp_shape(3) = nkpt
!   ABI_MALLOC(array, (tmp_shape(1), tmp_shape(2), tmp_shape(3), tmp_shape(4)))
!   do ikpt=1,nkpt
!     array(:,:,ikpt,:) = tmp_array(:,:,indkk(ikpt, 1),:)
!   end do
!   ABI_FREE(tmp_array)
!
! end subroutine downsample_array4
!
! subroutine downsample_array5(array, indkk, nkpt)
!
!   ! (ntemp, max_nbcalc, nkcalc, nsppol,:)
!   real(dp),allocatable,intent(inout) :: array(:,:,:,:,:)
!   integer,intent(in) :: indkk(:,:)
!
!   integer :: ikpt, nkpt
!   integer :: tmp_shape(5)
!   real(dp),allocatable :: tmp_array(:,:,:,:,:)
!
!   ABI_MOVE_ALLOC(array, tmp_array)
!   tmp_shape = shape(array)
!   tmp_shape(3) = nkpt
!   ABI_MALLOC(array, (tmp_shape(1), tmp_shape(2), tmp_shape(3), tmp_shape(4), tmp_shape(5)))
!   do ikpt=1,nkpt
!     array(:,:,ikpt,:,:) = tmp_array(:,:,indkk(ikpt, 1),:,:)
!   end do
!   ABI_FREE(tmp_array)
!
! end subroutine downsample_array5

end function rta_new
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_compute
!! NAME
!! rta_compute
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine rta_compute(self, cryst, dtset, comm)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 class(rta_t),intent(inout) :: self
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst

!Local variables ------------------------------
 integer,parameter :: nvecs0 = 0, master = 0
 integer :: nsppol, nkpt, ib, ik_ibz, iw, spin, ii, jj, itemp, irta, itens, iscal, cnt
 integer :: ntens, edos_intmeth, ifermi, iel, nvals, my_rank
 real(dp) :: emin, emax, edos_broad, edos_step, max_occ, kT, Tkelv, linewidth, fact0
 real(dp) :: cpu, wall, gflops
 real(dp) :: vr(3), dummy_vecs(1,1,1,1,1), work_33(3,3), S_33(3,3)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:,:), out_valsdos(:,:,:,:), dummy_dosvecs(:,:,:,:,:)
 real(dp),allocatable :: out_tensdos(:,:,:,:,:,:), tau_vals(:,:,:,:,:), l0inv_33nw(:,:,:)
 !character(len=500) :: msg

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm)

 ! Basic dimensions
 nsppol = self%ebands%nsppol; nkpt = self%ebands%nkpt

 ! Allocate vv tensors with and without the lifetimes. Eq 8 of [[cite:Madsen2018]]
 ! The total number of tensorial entries is ntens and accounts for nrta
 ! Remember that we haven't computed all the k-points in the IBZ hence we can have zero linewidths
 ! or very small values when the states is at the band edge so use safe_dif to avoid SIGFPE.
 ! Also, note how we store only the states in the energy window

 nvals = self%ntemp * self%nrta
 ABI_CALLOC(tau_vals, (self%ntemp, self%nrta, self%bmin:self%bmax, nkpt, nsppol))

 ntens = (1 + self%ntemp) * self%nrta
 ABI_CALLOC(vv_tens, (3, 3, 1 + self%ntemp, self%nrta, self%bmin:self%bmax, nkpt, nsppol))

 cnt = 0
 do spin=1,nsppol
   do ik_ibz=1,nkpt
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     do ib=self%bmin,self%bmax

       vr(:) = self%velocity(:, ib, ik_ibz, spin)
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
 call cwtime_report(" rta_compute_loop1", cpu, wall, gflops)

 ! Compute DOS and VV_DOS and VV_TAU_DOS
 ! Define integration method and mesh step.
 edos_intmeth = 2; if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae
 if (edos_step == zero) edos_step = 0.001
 !if (edos_step == zero) edos_step = ten / Ha_meV
 edos_broad = dtset%tsmear

 ! Set default energy range for DOS
 ! If sigma_erange is set, get emin and emax from this variable
 ! MG: TODO This value should be read from SIGPEH
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
 !    out_valsdos: (nw, 2, nvals, nsppol) array with DOS for scalar quantities if nvals > 0
 !    out_tensdos: (nw, 2, 3, 3, ntens,  nsppol) array with DOS weighted by tensorial terms if ntens > 0
 !
 !  Vectors and tensors are in Cartesian coordinates.
 !  Note how we compute the DOS only between [emin, emax] to save time and memory
 !  this implies that IDOS and edos%ifermi are ill-defined

 self%edos = ebands_get_edos_matrix_elements(self%ebands, cryst, self%bsize, &
                                             nvals, tau_vals, nvecs0, dummy_vecs, ntens, vv_tens, &
                                             edos_intmeth, edos_step, edos_broad, &
                                             out_valsdos, dummy_dosvecs, out_tensdos, comm, &
                                             brange=[self%bmin, self%bmax], erange=[emin, emax])

 if (my_rank == master) then
   call self%edos%print(unit=std_out, header="Computation of DOS, VV_DOS and VVTAU_DOS")
   call self%edos%print(unit=ab_out,  header="Computation of DOS, VV_DOS and VVTAU_DOS")
 end if

 call cwtime_report(" rta_compute_edos", cpu, wall, gflops)

 ! Unpack data stored in out_tensdos with shape (nw, 2, 3, 3, ntens, nsppol)
 self%nw = self%edos%nw
 ABI_MALLOC(self%tau_dos, (self%nw, self%ntemp, nsppol, self%nrta))
 ! TODO: Exchange dims?
 ABI_MALLOC(self%vv_dos, (self%nw, 3, 3, nsppol))
 ABI_MALLOC(self%vvtau_dos, (self%nw, 3, 3, self%ntemp, nsppol, self%nrta))

 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp+1

       itens = itemp + (irta - 1) * (self%ntemp + 1)
       if (itemp == 1) then
         self%vv_dos(:,:,:,spin) = out_tensdos(:, 1, :, :, itens, spin)
       else
         self%vvtau_dos(:,:,:, itemp-1, spin, irta) = out_tensdos(:, 1, :, :, itens, spin)
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

 ABI_MALLOC(self%l0, (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_MALLOC(self%l1, (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_MALLOC(self%l2, (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))

 call onsager(0, self%l0)
 call onsager(1, self%l1)
 call onsager(2, self%l2)

 call cwtime_report(" rta_compute_onsanger", cpu, wall, gflops)

 ! Compute transport tensors, Eqs 12-15 of [[cite:Madsen2018]] and convert to SI units.
 ABI_CALLOC(self%sigma,   (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%seebeck, (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%kappa,   (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%pi,      (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%zte,     (3, 3, self%nw, self%ntemp, self%nsppol, self%nrta))

 ! Sigma = L0
 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)
 self%sigma = fact0 * self%l0

 ! Used to stored L0^-1
 ABI_MALLOC(l0inv_33nw, (3, 3, self%nw))

 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp

       TKelv = self%kTmesh(itemp) / kb_HaK; if (TKelv < one) Tkelv = one

       ! S = -1/T L0^-1 L1 = -1/T sigma L1
       do iw=1,self%nw
         call inv33(self%l0(:, :, iw, itemp, spin, irta), work_33)
         l0inv_33nw(:,:,iw) = work_33
         self%seebeck(:,:,iw,itemp,spin,irta) = - (volt_SI / TKelv) * matmul(work_33, self%l1(:,:,iw,itemp,spin,irta))
       end do

       ! kappa = 1/T [L2 - L1 L0^-1 L1]
       ! HM: Check why do we need minus sign here to get consistent results with Boltztrap!
       ! MG: Likely because of a different definition of kappa.
       do iw=1,self%nw
         work_33 = self%l1(:, :, iw, itemp, spin, irta)
         work_33 = self%l2(:, :, iw, itemp, spin, irta) - matmul(work_33, matmul(l0inv_33nw(:, :, iw), work_33))
         !self%kappa(:,:,iw,itemp,spin,irta) = - (volt_SI**2 * fact0 / TKelv) * work_33
         self%kappa(:,:,iw,itemp,spin,irta) = + (volt_SI**2 * fact0 / TKelv) * work_33
       end do

       ! Peltier pi = -L1 L0^-1
       do iw=1,self%nw
         work_33 = self%l1(:, :, iw, itemp, spin, irta)
         self%pi(:,:,iw,itemp,spin,irta) = - volt_SI * matmul(work_33, l0inv_33nw(:, :, iw))
       end do

       ! ZT:  S^T sigma S k^-1 T (tensor form with k=k_electronic only):
       do iw=1,self%nw
         S_33 = self%seebeck(:,:,iw,itemp,spin,irta)
         S_33 = matmul(matmul(transpose(S_33), self%sigma(:,:,iw,itemp,spin,irta)), S_33)
         call inv33(self%kappa(:,:,iw,itemp,spin,irta), work_33)
         self%zte(:,:,iw,itemp,spin,irta) = matmul(S_33, work_33) * TKelv
       end do

     end do
   end do
 end do

 ABI_FREE(l0inv_33nw)

 ! Compute the index of the Fermi level and handle possible out of range condition.
 ifermi = bisect(self%edos%mesh, self%ebands%fermie)
 if (ifermi == 0 .or. ifermi == self%nw) then
   MSG_ERROR("Bisection could not find the index of the Fermi level in edos%mesh!")
 end if

 ! Mobility
 max_occ = two / (self%nspinor * self%nsppol)
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
               call safe_div(self%sigma(ii, jj, iw, itemp, spin, irta) * 100**2, &
                             e_Cb * self%n(iw, itemp, iel), &
                             zero, self%mobility(ii, jj, iw, itemp, iel, spin, irta))
             end do
           end do
         end do
       end do
     end do
   end do ! itemp
 end do ! spin

 call cwtime_report(" rta_compute", cpu, wall, gflops)

contains

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
 real(dp),intent(out) :: lorder(3, 3, self%nw, self%ntemp,self%nsppol,self%nrta)

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
           kernel(iw,:,:,:) = - max_occ * self%vvtau_dos(iw,:,:,itemp,:,irta) * (ee - mu)**order * occ_dfde(ee, kT, mu)
         else
           kernel(iw,:,:,:) = - max_occ * self%vvtau_dos(iw,:,:,itemp,:,irta) * occ_dfde(ee, kT, mu)
         end if
       end do

       ! Integrate with simpson_int
       do spin=1,self%nsppol
         do jj=1,3
           do ii=1,3
             call simpson_int(self%nw, edos_step, kernel(:,ii,jj, spin), integral)
             lorder(ii, jj, imu, itemp, spin, irta) = integral(self%nw)
           end do
         end do
       end do

     end do ! imu
   end do ! itemp
 end do ! irta

 end subroutine onsager

end subroutine rta_compute
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_compute_mobility
!! NAME
!! rta_compute_mobility
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! comm=MPI communicator.
!!
!! SOURCE

subroutine rta_compute_mobility(self, cryst, dtset, comm)

!Arguments ------------------------------------
 class(rta_t),intent(inout) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: comm

!Local variables ------------------------------
 integer :: nsppol, nkpt, ib, ik_ibz, spin, ii, jj, itemp, ieh, cnt, nprocs, irta
 real(dp) :: eig_nk, mu_e, linewidth, fact, fact0, max_occ, kT, wtk
 real(dp) :: cpu, wall, gflops
 real(dp) :: vr(3), vv_tens(3,3), vv_tenslw(3,3)

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ABI_UNUSED(dtset%natom)
 nprocs = xmpi_comm_size(comm)

 ! Copy important dimensions
 nkpt = self%ebands%nkpt; nsppol = self%ebands%nsppol
 max_occ = two / (self%nspinor * self%nsppol)

 ABI_MALLOC(self%mobility_mu, (3, 3, 2, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%ne, (self%ntemp))
 ABI_CALLOC(self%nh, (self%ntemp))

 ! Compute carrier concentration
 do spin=1,nsppol
   do ik_ibz=1,nkpt
     wtk = self%ebands%wtk(ik_ibz)
     do ib=1, self%ebands%nband(ik_ibz + (spin-1)* self%ebands%nkpt)
       eig_nk = self%ebands%eig(ib, ik_ibz, spin)

       do itemp=1,self%ntemp
         kT = self%kTmesh(itemp)
         mu_e = self%transport_mu_e(itemp)
         if (eig_nk >= mu_e) then
           self%ne(itemp) = self%ne(itemp) + wtk * occ_fd(eig_nk, kT, mu_e) * max_occ
         else
           self%nh(itemp) = self%nh(itemp) + wtk * (one - occ_fd(eig_nk, kT, mu_e)) * max_occ
         end if
       end do

     end do
   end do
 end do

 ! Get units conversion factor and spin degeneracy
 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)
 fact = max_occ * fact0 / e_Cb * 100**2

 ! Compute mobility_mu i.e. results in which lifetimes have been computed in a consistent way
 ! with the same the Fermi level. In all the other cases, indeed, we assume that tau does not depend on ef.
 !
 ! TODO: Implement other tensors. Compare these results with the ones obtained with spectral sigma
 ! In principle, they should be the same, in practice the integration of sigma requires enough resolution
 ! around the band edge.
 self%mobility_mu = zero
 cnt = 0
 do spin=1,nsppol
   do ik_ibz=1,nkpt
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     wtk = self%ebands%wtk(ik_ibz)

     do ib=self%bmin,self%bmax
       eig_nk = self%ebands%eig(ib, ik_ibz, spin)

       ! Store outer product in vv_tens
       vr(:) = self%velocity(:, ib, ik_ibz, spin)
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj) = vr(ii) * vr(jj)
         end do
       end do
       ! Symmetrize tensor
       vv_tens = cryst%symmetrize_cart_tens33(vv_tens)

       ! Multiply by the lifetime (SERTA or MRTA)
       do irta=1,self%nrta
         do itemp=1,self%ntemp
           kT = self%kTmesh(itemp)
           mu_e = self%transport_mu_e(itemp)
           ieh = 2; if (eig_nk >= mu_e) ieh = 1
           linewidth = self%linewidths(itemp, ib, ik_ibz, spin, irta)
           call safe_div( - wtk * vv_tens(:, :) * occ_dfde(eig_nk, kT, mu_e), two * linewidth, zero, vv_tenslw(:, :))
           self%mobility_mu(:, :, ieh, itemp, spin, irta) = self%mobility_mu(:, :, ieh, itemp, spin, irta) &
             + vv_tenslw(:, :)
         end do
       end do
     end do

   end do ! ik_ibz
 end do ! spin

 !call xmpi_sum(self%mobility_mu, comm, ierr)

 ! Scale by the carrier concentration
 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp
       ! electrons
       call safe_div(fact * self%mobility_mu(:,:,1,itemp,spin,irta), &
                     self%ne(itemp) / cryst%ucvol / Bohr_meter**3, zero, self%mobility_mu(:,:,1,itemp,spin,irta) )
       ! holes
       call safe_div(fact * self%mobility_mu(:,:,2,itemp,spin,irta), &
                     self%nh(itemp) / cryst%ucvol / Bohr_meter**3, zero, self%mobility_mu(:,:,2,itemp,spin,irta) )
     end do
   end do
 end do

 call cwtime_report(" rta_compute_mobility", cpu, wall, gflops)

end subroutine rta_compute_mobility
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
 type(rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 integer,intent(in) :: ncid

!Local variables --------------------------------
 integer :: ncerr
 real(dp) :: cpu, wall, gflops
 real(dp) :: work(dtset%nsppol)

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

#ifdef HAVE_NETCDF
 ! Write to netcdf file
 ncerr = nctk_def_dims(ncid, [ &
    nctkdim_t("ntemp", self%ntemp), nctkdim_t("nrta", self%nrta), nctkdim_t("nsppol", self%nsppol)], defmode=.True.)
 NCF_CHECK(ncerr)

 NCF_CHECK(cryst%ncwrite(ncid))
 NCF_CHECK(ebands_ncwrite(self%ebands, ncid))
 NCF_CHECK(self%edos%ncwrite(ncid))

 !nctk_copy from sigeph?
 !    nctkarr_t("eph_ngqpt_fine", "int", "three"), &

 ncerr = nctk_def_arrays(ncid, [ &
    nctkarr_t('transport_ngkpt', "int", "three"), &
    nctkarr_t('sigma_erange', "dp", "two"), &
    nctkarr_t('kTmesh', "dp", "ntemp"), &
    nctkarr_t('transport_mu_e', "dp", "ntemp"), &
    nctkarr_t('eph_mu_e', "dp", "ntemp"), &
    nctkarr_t('vb_max', "dp", "nsppol"), &
    nctkarr_t('cb_min', "dp", "nsppol"), &
    nctkarr_t('vv_dos', "dp", "edos_nw, three, three, nsppol"), &
    nctkarr_t('vvtau_dos', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('tau_dos', "dp", "edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('L0', "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('L1', "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('L2', "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('sigma',   "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('kappa',   "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('zte',   "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('seebeck', "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('pi',      "dp", "three, three, edos_nw, ntemp, nsppol, nrta"), &
    nctkarr_t('mobility',"dp", "three, three, edos_nw, ntemp, two, nsppol, nrta"), &
    nctkarr_t('N',  "dp", "edos_nw, ntemp, two"), &
    nctkarr_t('mobility_mu',"dp", "three, three, two, ntemp, nsppol, nrta")], &
 defmode=.True.)
 NCF_CHECK(ncerr)

 ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
    "eph_extrael", "eph_fermie", "transport_extrael", "transport_fermie"])
 NCF_CHECK(ncerr)

 NCF_CHECK(nctk_set_datamode(ncid))
 ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   "eph_extrael", "eph_fermie", "transport_extrael", "transport_fermie"], &
   [self%eph_extrael, self%eph_fermie, self%transport_extrael, self%transport_fermie])
 NCF_CHECK(ncerr)

 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "transport_ngkpt"), dtset%transport_ngkpt))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "sigma_erange"), dtset%sigma_erange))
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
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility_mu"), self%mobility_mu))
#endif

 call cwtime_report(" rta_ncwrite", cpu, wall, gflops)

end subroutine rta_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_print_txt_files
!! NAME
!! rta_print_txt_files
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! dtfil<datafiles_type>=variables related to files.
!!
!! SOURCE

subroutine rta_print_txt_files(self, cryst, dtset, dtfil)

!Arguments --------------------------------------
 class(rta_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil

!Local variables --------------------------------
 integer :: itemp, spin, irta, ii
 integer :: unts(2)
 character(len=500) :: msg, pre, rta_type
 character(len=2) :: components(3)

!************************************************************************

 unts = [std_out, ab_out]

 call wrtout(unts, ch10//' Transport (RTA) calculation results:', newlines=1)
 components = ["xx", "yy", "zz"]

 do irta=1,self%nrta
   if (irta == 1) rta_type = "SERTA"
   if (irta == 2) rta_type = "MRTA"

   do ii=1,3
     call wrtout(unts, sjoin(" Cartesian component of", rta_type, "mobility tensor:", components(ii)))
     write(msg, "(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
     call wrtout(unts, msg)

     do spin=1,self%nsppol
       if (self%nsppol == 2) call wrtout(unts, sjoin(" For spin:", stoa(spin)), newlines=1)
       do itemp=1,self%ntemp
         write(msg,"(f16.2,2e16.2,2f16.2)") &
           self%kTmesh(itemp) / kb_HaK, &
           self%ne(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
           self%nh(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
           self%mobility_mu(ii, ii, 1, itemp, spin, irta), self%mobility_mu(ii, ii, 2, itemp, spin, irta)
           !self%transport_mu_e(itemp) * Ha_eV
         call wrtout(unts, msg)
       end do ! itemp
     end do ! spin
     call wrtout(unts, ch10)
   end do

   call wrtout(unts, ch10)
 end do ! irta

 do irta=1,self%nrta
   select case (irta)
   case (1)
     pre = "_SERTA"
   case (2)
     pre = "_MRTA"
   case default
     MSG_ERROR(sjoin("Don't know how to handle irta:", itoa(irta)))
   end select
   call self%write_tensor(dtset, irta, "sigma", self%sigma(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SIGMA"))
   call self%write_tensor(dtset, irta, "seebeck", self%seebeck(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_SBK"))
   call self%write_tensor(dtset, irta, "kappa", self%kappa(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_KAPPA"))
   call self%write_tensor(dtset, irta, "zte", self%zte(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_ZTE"))
   call self%write_tensor(dtset, irta, "pi", self%pi(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_PI"))
 end do

end subroutine rta_print_txt_files
!!***

!----------------------------------------------------------------------

!!****f* m_rta/rta_write_tensor
!! NAME
!!
!! FUNCTION
!!
!! INPUTS
!!
!! SOURCE

subroutine rta_write_tensor(self, dtset, irta, header, values, path)

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
   MSG_ERROR(msg)
 end if

 if (irta == 1) rta_type = "RTA type: Self-energy relaxation time approximation (SERTA)"
 if (irta == 2) rta_type = "RTA type: Momentum relaxation time approximation (MRTA)"

 ! write header
 write(ount, "(2a)")"# ", trim(header)
 write(ount, "(a)")"# ", trim(rta_type)
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

 ! (nw, 3, 3, ntemp, nsppol)
 if (self%nsppol == 1) then
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)")"# Energy [Ha], (xx, yx, yx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor."
     do iw=1,self%nw
       write(ount, "(10(es16.6))")self%edos%mesh(iw), tmp_values(:, :, iw, itemp, 1)
     end do
   end do
  write(ount, "(a)")""
 else
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)") &
       "# Energy [Ha], (xx, yx, yx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor for spin up followed by spin down."
     do iw=1,self%nw
       write(ount, "(19(es16.6))")self%edos%mesh(iw), tmp_values(:, :, iw, itemp, 1), tmp_values(:, :, iw, itemp, 2)
     end do
   end do
  write(ount, "(a)")""
 end if

 close(ount)

 ABI_FREE(tmp_values)

end subroutine rta_write_tensor
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
 ABI_SFREE(self%kTmesh)
 ABI_SFREE(self%eminmax_spin)
 ABI_SFREE(self%eph_mu_e)
 ABI_SFREE(self%transport_mu_e)
 ABI_SFREE(self%velocity)
 ABI_SFREE(self%linewidths)
 ABI_SFREE(self%l0)
 ABI_SFREE(self%l1)
 ABI_SFREE(self%l2)
 ABI_SFREE(self%sigma)
 ABI_SFREE(self%mobility)
 ABI_SFREE(self%seebeck)
 ABI_SFREE(self%kappa)
 ABI_SFREE(self%zte)
 ABI_SFREE(self%pi)
 ABI_SFREE(self%mobility_mu)
 ABI_SFREE(self%nh)
 ABI_SFREE(self%ne)
 ABI_SFREE(self%n)

 call ebands_free(self%ebands)
 call self%gaps%free()
 call self%edos%free()

end subroutine rta_free
!!***

end module m_rta
!!***

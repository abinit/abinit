!!****m* ABINIT/m_rta
!! NAME
!!  m_rta
!!
!! FUNCTION
!!  Module to compute transport properties using the
!!  Boltzmann transport equation (BTE) in the relaxation-time approximation.
!!  Initially for electron mobility limited by electron-phonon scattering.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (HM)
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
 use iso_c_binding
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hide_blas
 use m_copy
 use m_ebands
 use m_nctk
 use m_sigmaph
 use m_dtset
 use m_dtfil
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t
 use m_io_tools,       only : open_file
 use m_time,           only : cwtime, cwtime_report
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : bisect, simpson_int, safe_div
 use m_fstrings,       only : strcat, sjoin, itoa, ltoa, stoa, ftoa
 use m_kpts,           only : listkk, kpts_timrev_from_kptopt
 use m_occ,            only : occ_fd, occ_dfd

 implicit none

 private
!!****

 public :: rta_driver ! main entry point for transport calculations withing RTA
!!****

!----------------------------------------------------------------------

!!****t* m_rta/rta_t
!! NAME
!! rta_t
!!
!! FUNCTION
!! Container for transport quantities in the relaxation time approximation
!!
!! SOURCE

type,public :: rta_t

   integer :: nsppol
   ! number of independent spin polarizations

   integer :: nspinor
   ! number of spinorial components

   integer :: ntemp
   ! number of temperatures

   integer :: ndop
   ! number of carrier concentrations at which to evaluate chemical potential energy

   integer :: nw
   ! number of frequencies at which transport quantities are computed

   integer :: nrta
   ! Number of relaxation-time approximations used (SERTA, MRTA)

   real(dp) :: eph_extrael
   ! extra electrons per unit cell from sigeph (lifetimes)

   real(dp) :: eph_fermie
   ! Fermi level from input file from sigeph (lifetimes)

   real(dp) :: transport_extrael
   ! extra electrons per unit cell

   real(dp) :: transport_fermie
   ! Fermi level from input file

   real(dp),allocatable :: kTmesh(:)
   ! a list of temperatures at which to compute the transport

   real(dp),allocatable :: eph_mu_e(:)
   ! (%ntemp, %ndop)
   ! Chemical potential at this carrier concentration and temperature from sigeph (lifetime)

   real(dp),allocatable :: transport_mu_e(:)
   ! (%ntemp, %ndop)
   ! Chemical potential at this carrier concentration and temperature

   real(dp),allocatable :: eminmax_spin(:,:)
   ! min/max energy of the original ebands object

   real(dp),allocatable :: linewidth_serta(:,:,:,:)
   ! Linewidth computed in the self-energy relaxation time aproximation

   real(dp),allocatable :: linewidth_mrta(:,:,:,:)
   ! Linewidth computed in the momentum relaxation time approximation

   real(dp),allocatable :: velocity(:,:,:,:)
   ! band velocity in Cartesian coordinates.

   type(gaps_t) :: gaps
   ! get gaps of original ebands object

   !integer :: nmu
   ! number of dopings

   !real(dp) :: nmesh(:)
   ! a list of carrier concentrations at which to compute transport

   type(ebands_t) :: ebands
   ! bandstructure object used to compute the transport properties
   ! Allocate using only the relevant bands for transport
   ! including valence states to allow to compute different doping

   type(edos_t) :: edos
   ! electronic density of states
   ! edos%mesh is the mesh used for vv_dos, vvtau_dos and tau_dos
   ! (%nw)

   real(dp),allocatable :: tau_dos(:,:,:,:)
   ! tau(e)  DOS
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
   ! (nw, 3, 3, ntemp, nsppol, nrta)
   ! Onsager coeficients

   real(dp),allocatable :: sigma(:,:,:,:,:,:)
   real(dp),allocatable :: seebeck(:,:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:,:)
   real(dp),allocatable :: pi(:,:,:,:,:,:)
   ! (nw, 3, 3, ntemp, nsppol, nrta)
   ! Transport coefficients

   real(dp),allocatable :: mobility(:,:,:,:,:,:,:)
   ! Mobility
   ! (%nw, 3, 3, %ntemp, 2, %nsppol, nrta)
   ! 5-th index is for e-h

   real(dp),allocatable :: n(:,:,:)
   ! (nw, ntemp,2) carrier density for e/h (n/cm^3)

   real(dp),allocatable :: mobility_mu(:,:,:,:,:,:)
   ! (2, 3, 3, ntemp, nsppol, nrta)
   ! mobility for electrons and holes (first dimension) at transport_mu_e(ntemp)
   ! First index is for e-h

 contains

    procedure :: compute => rta_compute
    procedure :: compute_mobility => rta_compute_mobility
    procedure :: print => rta_print
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
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! cryst<crystal_t>=Crystalline structure
!! comm=MPI communicator.
!!
!! SOURCE

subroutine rta_driver(dtfil, dtset, ebands, cryst, comm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands

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
 real(dp) :: extrael_fermie(2)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)

 call wrtout([std_out, ab_out], ch10//' Entering transport RTA computation driver.')
 call wrtout([std_out, ab_out], sjoin("- Reading carrier lifetimes from:", dtfil%filsigephin), newlines=1)

 sigmaph = sigmaph_read(dtfil%filsigephin, dtset, xmpi_comm_self, msg, ierr, &
                        keep_open=.true., extrael_fermie=extrael_fermie)
 ABI_CHECK(ierr == 0, msg)

 ! Initialize RTA object
 ! TODO: Should store more metadata: energy window, nkcalc ....
 rta = rta_new(dtset, sigmaph, cryst, ebands, extrael_fermie, comm)

 ! sigmaph is not needed anymore. Free it.
 sigmaph%ncid = nctk_noid
 call sigmaph%free()

 ! Compute RTA transport quantities
 call rta%compute(cryst, dtset, comm)

 ! Compute RTA mobility
 call rta%compute_mobility(cryst, dtset, comm)

 if (my_rank == master) then
   ! Print RTA results to stdout and other external txt files (for the test suite)
   call rta%print(cryst, dtset, dtfil)

   ! Creates the netcdf file used to store the results of the calculation.
#ifdef HAVE_NETCDF
   path = strcat(dtfil%filnam_ds(4), "_RTA.nc")
   call wrtout([std_out, ab_out], ch10//sjoin("- Writing RTA transport results to:", path))
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
!!  sigmaph<sigmaph_t>=Object with e-ph self-energy results.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  extrael_fermie
!!  comm=MPI communicator.
!!
!! SOURCE

type(rta_t) function rta_new(dtset, sigmaph, cryst, ebands, extrael_fermie, comm) result (new)

!Arguments -------------------------------------
 integer, intent(in) :: comm
 type(sigmaph_t),intent(in) :: sigmaph
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 real(dp),intent(in) :: extrael_fermie(2)

!Local variables ------------------------------
 integer,parameter :: occopt3 = 3, sppoldbl1 = 1, master = 0
 integer :: ierr, spin, nprocs, my_rank, timrev
 real(dp) :: dksqmax
 real(dp) :: cpu, wall, gflops
 character(len=500) :: msg
 type(ebands_t) :: tmp_ebands
!arrays
 integer :: kptrlatt(3,3)
 integer,allocatable :: indkk(:,:)

!************************************************************************

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call cwtime(cpu, wall, gflops, "start")

 ! Allocate temperature arrays (same as the ones used in the SIGEPH calculation).
 new%ntemp = sigmaph%ntemp
 ABI_MALLOC(new%kTmesh, (new%ntemp))
 new%kTmesh = sigmaph%kTmesh

 ! How many RTA approximations have we computed in sigmaph? (SERTA, MRTA)
 new%nrta = 2; if (sigmaph%mrta == 0) new%nrta = 1

 ! Information about the gaps
 new%nsppol = ebands%nsppol; new%nspinor = ebands%nspinor
 ABI_MALLOC(new%eminmax_spin, (2, ebands%nsppol))
 new%eminmax_spin = get_minmax(ebands, "eig")

 ierr = get_gaps(ebands, new%gaps)
 if (ierr /= 0) then
   do spin=1, ebands%nsppol
     MSG_WARNING(trim(new%gaps%errmsg_spin(spin)))
     new%gaps%vb_max(spin) = ebands%fermie - 1 * eV_Ha
     new%gaps%cb_min(spin) = ebands%fermie + 1 * eV_Ha
   end do
   MSG_WARNING("get_gaps returned non-zero exit status. See above warning messages...")
 end if

 if (my_rank == master) then
   call new%gaps%print(unit=std_out)
   call new%gaps%print(unit=ab_out)
 end if

 ! Read lifetimes to ebands object
 if (any(dtset%sigma_ngkpt /= 0)) then
   ! If integrals are computed with sigma_ngkpt k-mesh, we need to downsample ebands.
   call wrtout([std_out, ab_out], sjoin(" Computing integrals with downsampled sigma_ngkpt:", ltoa(dtset%sigma_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%sigma_ngkpt(1); kptrlatt(2,2) = dtset%sigma_ngkpt(2); kptrlatt(3,3) = dtset%sigma_ngkpt(3)

   tmp_ebands = ebands_downsample(ebands, cryst, kptrlatt, 1, [zero, zero, zero])
   new%ebands = sigmaph%get_ebands(cryst, tmp_ebands, new%linewidth_serta, new%linewidth_mrta, &
                                   new%velocity, xmpi_comm_self, ierr)
   call ebands_free(tmp_ebands)
 else
   new%ebands = sigmaph%get_ebands(cryst, ebands, new%linewidth_serta, new%linewidth_mrta, &
                                   new%velocity, xmpi_comm_self, ierr)
 end if

 !print *, "linewidth_serta", maxval(abs(new%linewidth_serta))
 !print *, "linewidth_mrta", maxval(abs(new%linewidth_mrta))

 ! FIXME: I think transport_ngkpt is buggy, wrong ne(T), weird zeros if MRTA ...
 ! Do we really need this option?

 if (any(dtset%transport_ngkpt /= 0)) then
   ! Perform further downsampling (usefull for debugging purposes)
   call wrtout([std_out, ab_out], " Downsampling the k-mesh before computing transport:")
   call wrtout([std_out, ab_out], sjoin(" Using transport_ngkpt: ", ltoa(dtset%transport_ngkpt)))
   kptrlatt = 0
   kptrlatt(1,1) = dtset%transport_ngkpt(1)
   kptrlatt(2,2) = dtset%transport_ngkpt(2)
   kptrlatt(3,3) = dtset%transport_ngkpt(3)
   tmp_ebands = ebands_downsample(new%ebands, cryst, kptrlatt, 1, [zero, zero, zero])

   ! Map the points of the downsampled bands to dense ebands
   timrev = kpts_timrev_from_kptopt(ebands%kptopt)
   ABI_MALLOC(indkk, (tmp_ebands%nkpt, 6))
   call listkk(dksqmax, cryst%gmet, indkk, new%ebands%kptns, tmp_ebands%kptns, &
               new%ebands%nkpt, tmp_ebands%nkpt, cryst%nsym, &
               sppoldbl1, cryst%symafm, cryst%symrec, timrev, comm, exit_loop=.True., use_symrec=.True.)

   if (dksqmax > tol12) then
      write(msg, '(3a,es16.6,a)' ) &
       "Error while downsampling ebands in transport driver",ch10, &
       "The k-point could not be generated from a symmetrical one. dksqmax: ",dksqmax, ch10
      MSG_ERROR(msg)
   end if

   ! Downsampling linewidths and velocities.
   call downsample_array(new%linewidth_serta, indkk, tmp_ebands%nkpt)
   call downsample_array(new%linewidth_mrta,  indkk, tmp_ebands%nkpt)
   call downsample_array(new%velocity,        indkk, tmp_ebands%nkpt)

   !print *, "after downsamplinglinewidth_serta", maxval(abs(new%linewidth_serta))
   !print *, "fter downsamplinglinewidth_mrta", maxval(abs(new%linewidth_mrta))

   ABI_SFREE(indkk)
   call ebands_free(new%ebands)
   call ebands_copy(tmp_ebands, new%ebands)
   call ebands_free(tmp_ebands)
   !call ebands_move_alloc(tmp_ebands, new%ebands)
 end if

 ! Same doping case as sigmaph
 new%ndop = 1
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

 contains
 subroutine downsample_array(array, indkk, nkpt)

   real(dp),allocatable,intent(inout) :: array(:,:,:,:)
   integer,intent(in) :: indkk(:,:)

   integer :: ikpt, nkpt
   integer :: tmp_shape(4)
   real(dp),allocatable :: tmp_array(:,:,:,:)

   ABI_MOVE_ALLOC(array, tmp_array)
   tmp_shape = shape(array)
   tmp_shape(3) = nkpt
   ABI_MALLOC(array, (tmp_shape(1), tmp_shape(2), tmp_shape(3), tmp_shape(4)))
   do ikpt=1,nkpt
     array(:,:,ikpt,:) = tmp_array(:,:,indkk(ikpt,1),:)
   end do
   ABI_FREE(tmp_array)

 end subroutine downsample_array

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
 integer :: nsppol, nkpt, mband, ib, ik, iw, spin, ii, jj, itemp, irta, itens, iscal, cnt
 integer :: ntens, edos_intmeth, ifermi, iel, nvals, my_rank
 real(dp) :: emin, emax, edos_broad, edos_step, max_occ, kT, linewidth, fact0
 real(dp) :: cpu, wall, gflops
 real(dp) :: vr(3), dummy_vecs(1,1,1,1,1)
 real(dp),allocatable :: vv_tens(:,:,:,:,:,:,:), out_valsdos(:,:,:,:), dummy_dosvecs(:,:,:,:,:)
 real(dp),allocatable :: out_tensdos(:,:,:,:,:,:), tau_vals(:,:,:,:,:)
 !character(len=500) :: msg

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")
 my_rank = xmpi_comm_rank(comm)

 ! Basic dimensions
 nsppol = self%ebands%nsppol; nkpt = self%ebands%nkpt; mband = self%ebands%mband

 ! Allocate vv tensors with and without the lifetimes. Eq 8 of [[cite:Madsen2018]]
 ! The total number of tensorial entries is ntens and accounts for nrta
 ! Remember that we haven't computed all the k-points in the IBZ hence we can have zero linewidths
 ! or very small values when the states is at the band edge so use safe_dif to avoid SIGFPE.
 nvals = self%ntemp * self%nrta
 ABI_CALLOC(tau_vals, (self%ntemp, self%nrta, mband, nkpt, nsppol))

 ntens = (1 + self%ntemp) * self%nrta
 ABI_CALLOC(vv_tens, (3, 3, 1 + self%ntemp, self%nrta, mband, nkpt, nsppol))

 cnt = 0
 do spin=1,nsppol
   do ik=1,nkpt
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     do ib=1,mband

       vr(:) = self%velocity(:, ib, ik, spin)
       ! Store outer product (v_bks x v_bks) in vv_tens. This part does not depend on T and irta.
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj, 1, 1:self%nrta, ib, ik, spin) = vr(ii) * vr(jj)
         end do
       end do

       ! Multiply by the lifetime (SERTA and MRTA)
       do irta=1,self%nrta
         do itemp=1,self%ntemp
           if (irta == 1) linewidth = abs(self%linewidth_serta(itemp, ib, ik, spin))
           if (irta == 2) linewidth = abs(self%linewidth_mrta(itemp, ib, ik, spin))
           call safe_div(vv_tens(:, :, 1, irta, ib, ik, spin), linewidth, zero, vv_tens(:, :, 1+itemp, irta, ib, ik, spin))
           call safe_div(one / two, linewidth, zero, tau_vals(itemp, irta, ib, ik, spin))
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
 edos_intmeth = 2
 if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae
 if (edos_step == zero) edos_step = 0.001
 if (edos_step == zero) edos_step = ten / Ha_meV
 edos_broad = dtset%tsmear

 ! Set default energy range for DOS
 emin = minval(self%eminmax_spin(1,:)); emin = emin - 0.1_dp * abs(emin)
 emax = maxval(self%eminmax_spin(2,:)); emax = emax + 0.1_dp * abs(emax)

 ! If sigma_erange is set, get emin and emax from this variable
 do spin=1,self%ebands%nsppol
   if (dtset%sigma_erange(1) >= zero) emin = self%gaps%vb_max(spin) + tol2 * eV_Ha - dtset%sigma_erange(1)
   if (dtset%sigma_erange(2) >= zero) emax = self%gaps%cb_min(spin) - tol2 * eV_Ha + dtset%sigma_erange(2)
 end do

 ! Compute DOS, vv_dos and vvtau_DOS (v x v tau)
 !    out_valsdos: (nw, 2, nvals, nsppol) array with DOS for scalar quantities if nvals > 0
 !    out_tensdos: (nw, 2, 3, 3, ntens,  nsppol) array with DOS weighted by tensorial terms if ntens > 0
 !  Vectors and tensors are in Cartesian coordinates.
 !  Note how we compute the DOS only between [emin, emax] to save time and memory
 !  this implies that IDOS and edos%ifermi are ill-defined

 self%edos = ebands_get_edos_matrix_elements(self%ebands, cryst, &
                                             nvals, tau_vals, nvecs0, dummy_vecs, ntens, vv_tens, &
                                             edos_intmeth, edos_step, edos_broad, comm, &
                                             out_valsdos, dummy_dosvecs, out_tensdos, emin=emin, emax=emax)

 if (my_rank == master) then
   call self%edos%print(unit=std_out, header="Computation of DOS, VV_DOS and VVTAU_DOS")
   call self%edos%print(unit=ab_out,  header="Computation of DOS, VV_DOS and VVTAU_DOS")
 end if

 call cwtime_report(" rta_compute_edos", cpu, wall, gflops)

 ! Unpack data stored in out_tensdos with shape (nw, 2, 3, 3, ntens, nsppol)
 self%nw = self%edos%nw
 ABI_MALLOC(self%tau_dos, (self%nw, self%ntemp, nsppol, self%nrta))
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

 ! Transfer tau(e)
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
 ABI_MALLOC(self%l0, (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_MALLOC(self%l1, (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_MALLOC(self%l2, (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))

 call onsager(0, self%l0)
 call onsager(1, self%l1)
 call onsager(2, self%l2)

 call cwtime_report(" rta_compute_onsanger", cpu, wall, gflops)

 ! Compute transport quantities, Eqs 12-15 of [[cite:Madsen2018]]
 ABI_MALLOC(self%sigma,   (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%seebeck, (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%kappa,   (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_MALLOC(self%pi,      (self%nw, 3, 3, self%ntemp, self%nsppol, self%nrta))

 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)

 self%sigma = fact0 * self%l0
 call safe_div(volt_SI * self%l1, self%l0, zero, self%pi)

 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp
       kT = self%kTmesh(itemp) / kb_HaK
       call safe_div(volt_SI * self%l1(:,:,:,itemp,spin,irta), kT * self%l0(:,:,:,itemp,spin,irta), zero, &
                     self%seebeck(:,:,:,itemp,spin,irta))

       ! HM: to write it as a single division I do: kappa = L1^2/L0 + L2 = (L1^2 + L2*L0)/L0
       ! Check why do we need minus sign here to get consistent results with Boltztrap!
       ! MG: Very likely because we don't use the same conventions in the defintion of the Onsager coefficients.
       call safe_div( -volt_SI**2 * fact0 * &
         (self%l1(:,:,:,itemp,spin,irta)**2 - self%l2(:,:,:,itemp,spin,irta) * self%l0(:,:,:,itemp,spin,irta)), &
         kT * self%l0(:,:,:,itemp,spin,irta), zero, self%kappa(:,:,:,itemp,spin,irta))
     end do
   end do
 end do

 ! Compute the index of the Fermi level and handle possible out of range condition.
 ifermi = bisect(self%edos%mesh, self%ebands%fermie)
 if (ifermi == 0 .or. ifermi == self%nw) then
   MSG_ERROR("Bisection could not find the index of the Fermi level in edos%mesh!")
 end if

 ! Mobility
 ABI_MALLOC(self%n, (self%nw, self%ntemp, 2))
 ABI_MALLOC(self%mobility, (self%nw, 3, 3, self%ntemp, 2, self%nsppol, self%nrta))

 self%mobility = zero
 max_occ = two / (self%nspinor * self%nsppol)

 do spin=1,self%nsppol
   do itemp=1,self%ntemp
     ! Compute carrier density
     kT = self%kTmesh(itemp)

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
         do ii=1,3
           do jj=1,3
             do iw=1,self%nw
               call safe_div(self%sigma(iw, ii, jj, itemp, spin, irta) * 100**2, &
                             e_Cb * self%n(iw, itemp, iel), &
                             zero, self%mobility(iw, ii, jj, itemp, iel, spin, irta))
             end do
           end do
         end do
       end do
     end do
   end do ! itemp
 end do ! spin

 call cwtime_report(" rta_compute", cpu, wall, gflops)

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

 subroutine onsager(order, lorder)

 !Arguments -------------------------------------------
 integer,intent(in) :: order
 real(dp),intent(out) :: lorder(self%nw,3,3,self%ntemp,self%nsppol,self%nrta)

 !Local variables -------------------------------------
 integer :: spin, iw, imu, irta
 real(dp) :: fact, mu, ee, kT
 real(dp) :: kernel(self%nw,3,3,self%nsppol), integral(self%nw)

 ! Get spin degeneracy
 max_occ = two / (self%nspinor*self%nsppol)
 ! 2 comes from linewidth-lifetime relation because we divided by the linewidth and now by (2 * linewidth)
 fact = max_occ / two

 do irta=1,self%nrta
   do itemp=1,self%ntemp
     kT = self%kTmesh(itemp)
     do imu=1,self%nw
       mu = self%edos%mesh(imu)

       ! Build integrand for mu
       do iw=1,self%nw
         ee = self%edos%mesh(iw)
         if (order > 0) then
           kernel(iw,:,:,:) = fact * self%vvtau_dos(iw,:,:,itemp,:,irta) * (mu - ee)**order * occ_dfd(ee, kT, mu)
         else
           kernel(iw,:,:,:) = fact * self%vvtau_dos(iw,:,:,itemp,:,irta) * occ_dfd(ee, kT, mu)
         end if
       end do

       ! Integrate with simpson_int
       do spin=1,self%nsppol
         do jj=1,3
           do ii=1,3
             call simpson_int(self%nw, edos_step, kernel(:,ii,jj, spin), integral)
             lorder(imu, ii, jj, itemp, spin, irta) = integral(self%nw)
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
 integer :: nsppol, nkpt, mband, ib, ik, spin, ii, jj, itemp, ielhol, nvalence, cnt, nprocs, irta
 real(dp) :: eig_nk, mu_e, linewidth, fact, fact0, max_occ, kT, wtk
 real(dp) :: cpu, wall, gflops
 real(dp) :: vr(3), vv_tens(3,3), vv_tenslw(3,3)

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ABI_UNUSED(dtset%natom)
 nprocs = xmpi_comm_size(comm)

 ! create alias for dimensions
 mband = self%ebands%mband; nkpt = self%ebands%nkpt; nsppol = self%ebands%nsppol

 ! Compute index of valence band
 ! TODO:
 !  1) Should add nelect0 to ebands to keep track of intrinsic
 !  2) Generalize expression to metals using mu_e
 max_occ = two / (self%nspinor * self%nsppol)
 nvalence = nint((self%ebands%nelect - self%eph_extrael) / max_occ)

 ABI_MALLOC(self%mobility_mu, (2, 3, 3, self%ntemp, self%nsppol, self%nrta))
 ABI_CALLOC(self%ne, (self%ntemp))
 ABI_CALLOC(self%nh, (self%ntemp))

 ! Compute carrier concentration
 do spin=1,nsppol
   do ik=1,nkpt
     wtk = self%ebands%wtk(ik)

     ! number of holes
     do ib=1,nvalence
       eig_nk = self%ebands%eig(ib, ik, spin)
       do itemp=1,self%ntemp
         kT = self%kTmesh(itemp)
         mu_e = self%transport_mu_e(itemp)
         self%nh(itemp) = self%nh(itemp) + wtk * (one - occ_fd(eig_nk, kT, mu_e)) * max_occ
       end do
     end do

     ! number of electrons
     do ib=nvalence+1,mband
       eig_nk = self%ebands%eig(ib, ik, spin)
       do itemp=1,self%ntemp
         kT = self%kTmesh(itemp)
         mu_e = self%transport_mu_e(itemp)
         self%ne(itemp) = self%ne(itemp) + wtk * occ_fd(eig_nk, kT, mu_e) * max_occ
       end do
     end do

   end do
 end do

 ! Get units conversion factor and spin degeneracy
 fact0 = (Time_Sec * siemens_SI / Bohr_meter / cryst%ucvol)
 fact = max_occ * fact0 / e_Cb / two * 100**2

 ! Compute mobility_mu i.e. results in which lifetimes have been computed in a consistent way
 ! with the same the Fermi level. In all the other cases, indeed we assume tau does not depend on ef.
 self%mobility_mu = zero
 cnt = 0
 do spin=1,nsppol
   do ik=1,nkpt
     !cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
     wtk = self%ebands%wtk(ik)

     do ib=1,mband
       ielhol = 2; if (ib > nvalence) ielhol = 1
       eig_nk = self%ebands%eig(ib, ik, spin)

       ! Store outer product in vv_tens
       vr(:) = self%velocity(:,ib,ik,spin)
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
           mu_e = self%transport_mu_e(itemp)
           kT = self%kTmesh(itemp)
           if (irta == 1) linewidth = abs(self%linewidth_serta(itemp, ib, ik, spin))
           if (irta == 2) linewidth = abs(self%linewidth_mrta(itemp, ib, ik, spin))
           call safe_div( wtk * vv_tens(:, :) * occ_dfd(eig_nk, kT, mu_e), linewidth, zero, vv_tenslw(:, :))
           self%mobility_mu(ielhol, :, :, itemp, spin, irta) = self%mobility_mu(ielhol, :, :, itemp, spin, irta) &
             + vv_tenslw(:, :)
         end do
       end do
     end do

   end do ! kpt
 end do ! spin

 !call xmpi_sum(self%mobility_mu, comm, ierr)

 ! Scale by the carrier concentration
 do irta=1,self%nrta
   do spin=1,nsppol
     do itemp=1,self%ntemp
       ! electrons
       call safe_div(fact * self%mobility_mu(1,:,:,itemp,spin,irta), &
                     self%ne(itemp) / cryst%ucvol / Bohr_meter**3, zero, self%mobility_mu(1,:,:,itemp,spin,irta) )
       ! holes
       call safe_div(fact * self%mobility_mu(2,:,:,itemp,spin,irta), &
                     self%nh(itemp) / cryst%ucvol / Bohr_meter**3, zero, self%mobility_mu(2,:,:,itemp,spin,irta) )
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
    nctkarr_t('L0', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('L1', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('L2', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('sigma',   "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('kappa',   "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('seebeck', "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('pi',      "dp", "edos_nw, three, three, ntemp, nsppol, nrta"), &
    nctkarr_t('mobility',"dp", "edos_nw, three, three, ntemp, two, nsppol, nrta"), &
    nctkarr_t('N',  "dp", "edos_nw, ntemp, two"), &
    nctkarr_t('mobility_mu',"dp", "two, three, three, ntemp, nsppol, nrta")], &
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
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vb_max"), self%gaps%vb_max))
 NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "cb_min"), self%gaps%cb_min))
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

!!****f* m_rta/rta_print
!! NAME
!! rta_print
!!
!! FUNCTION
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!!
!! SOURCE

subroutine rta_print(self, cryst, dtset, dtfil)

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
     call wrtout(unts, sjoin("Cartesian component of", rta_type, "mobility tensor:", components(ii)))
     write(msg, "(a16,a32,a32)") 'Temperature [K]', 'e/h density [cm^-3]', 'e/h mobility [cm^2/Vs]'
     call wrtout(unts, msg)

     do spin=1,self%nsppol
       if (self%nsppol == 2) call wrtout([std_out, ab_out], sjoin(" For spin:", stoa(spin)), newlines=1)
       do itemp=1,self%ntemp
         write(msg,"(f16.2,2e16.2,2f16.2)") &
           self%kTmesh(itemp) / kb_HaK, &
           self%ne(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
           self%nh(itemp) / cryst%ucvol / (Bohr_meter * 100)**3, &
           self%mobility_mu(1, ii, ii, itemp, spin, irta), self%mobility_mu(2, ii, ii, itemp, spin, irta)
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
   call self%write_tensor(dtset, irta, "pi", self%pi(:,:,:,:,:,irta), strcat(dtfil%filnam_ds(4), pre, "_PI"))
 end do

end subroutine rta_print
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

 ! (nw, 3, 3, ntemp, nsppol)
 if (self%nsppol == 1) then
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)")"# Energy [Ha], (xx, yx, yx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor."
     do iw=1,self%nw
       write(ount, "(10(es16.6))")self%edos%mesh(iw), values(iw, :, :, itemp, 1)
     end do
   end do
  write(ount, "(a)")""
 else
   do itemp=1, self%ntemp
     write(ount, "(/, a, 1x, f16.2)")"# T = ", self%kTmesh(itemp) / kb_HaK
     write(ount, "(a)") &
       "# Energy [Ha], (xx, yx, yx, xy, yy, zy, xz, yz, zz) Cartesian components of tensor for spin up followed by spin down."
     do iw=1,self%nw
       write(ount, "(19(es16.6))")self%edos%mesh(iw), values(iw, :, :, itemp, 1), values(iw, :, :, itemp, 2)
     end do
   end do
  write(ount, "(a)")""
 end if

 close(ount)

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
 ABI_SFREE(self%linewidth_mrta)
 ABI_SFREE(self%linewidth_serta)
 ABI_SFREE(self%l0)
 ABI_SFREE(self%l1)
 ABI_SFREE(self%l2)
 ABI_SFREE(self%sigma)
 ABI_SFREE(self%mobility)
 ABI_SFREE(self%seebeck)
 ABI_SFREE(self%kappa)
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

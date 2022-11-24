!!****m* ABINIT/m_cumulant
!! NAME
!!  m_cumulant
!!
!! FUNCTION
!!  Module to compute the e-ph spectral functions within the cumulant formalism
!!  and, optionally, transport properties within the Kubo formalism.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (JCA, MG)
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

module m_cumulant


 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_ebands
 use m_nctk
 use m_sigmaph
 use m_dtset
 use m_dtfil
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes,   only : ebands_t
 use defs_abitypes,     only : MPI_type
 use m_io_tools,       only : open_file, file_exists, is_open
 use m_time,           only : cwtime, cwtime_report
 use m_crystal,        only : crystal_t
 use m_numeric_tools,  only : simpson_cplx, arth, c2r, simpson, safe_div, simpson_int, ctrap, linfit, linspace
 use m_fstrings,       only : strcat, sjoin, itoa, ltoa, stoa, ftoa
 use m_distribfft,     only : init_distribfft_seq
 !use m_kpts,           only : kpts_timrev_from_kptopt
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_fft,             only : fourdp
 use m_fftcore,       only : ngfft_seq,fftalg_isavailable
 use m_occ,            only : occ_fd, occ_dfde
 !use m_pawtab,         only : pawtab_type

 implicit none

 private
!!****

 public :: cumulant_driver                 ! Main entry point for cumulant computations
!!****

!----------------------------------------------------------------------

!!****t* m_cumulant/cumulant_t
!! NAME
!! cumulant_t
!!
!! FUNCTION
!! Container for cumulant and transport quantities
!!
!! SOURCE

 type,public :: cumulant_t

  integer :: my_nspins
   ! Number of spins treated by this MPI rank

  integer :: nkcalc
   ! Total number of k-points computed. global variable

  integer :: my_nkcalc
   ! Number of k-points treated by this MPI rank i.e. local variable.

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin). global variable

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: nqbz
   ! Number of q-points in the (dense) BZ for sigma integration

  integer :: nqibz
   ! Number of q-points in the (dense) IBZ for sigma integration


  integer :: nbsum
   ! Total number of bands used in sum over states without taking into account MPI distribution.

  integer :: natom3
   ! 3 * natom.

  integer :: nwr
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Read from SIGEPH.nc
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 0 (taken from dtset%nfreqsp)

  integer :: nwr_ce
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Read from SIGEPH.nc
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 0 (taken from dtset%nfreqsp)


  integer :: ntemp
   ! Number of temperatures.

  integer :: comm
   !

  !type(xcomm_t) :: all_comm

  type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins (high-level)

  type(xcomm_t) :: kcalc_comm
   ! MPI communicator for parallelism over k-points

  type(xcomm_t) :: wt_comm
    ! MPI communicator for time/frequency domain (low-level)
    ! NB: memory is not distributed

  type(xcomm_t) :: ncwrite_comm
   ! MPI communicator for parallel netcdf IO used to write results for the different k-points/spins

  !integer :: coords(5)
   ! Cartesian coordinates of this processor in the Cartesian grid.

  ! real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  integer :: ngqpt(3)
   ! Number of divisions in the Q mesh in the BZ.


  integer,allocatable :: kcalc2ebands(:,:)
   ! Mapping ikcalc --> ebands IBZ
   ! Note that this array is not necessarily equation to kcalc2ibz computed in sigmaph
   ! because we may have used sigma_nkpt to downsample the initial nkpt mesh.
   ! This array is computed in get_ebands and is equal to kcalc2ibz if sigma_nkpt == ngkpt



   real(dp),allocatable :: linewidths(:,:,:,:,:)
   ! (ntemp, bmin:bmax, nkpt, nsppol, nrta)
   ! Linewidth in the IBZ computed in the SERTA/MRTA.
   ! Non-zero only for the kcalc k-points.


   real(dp),allocatable :: vbks(:,:,:,:)
   ! (3, bmin:bmax, nkpt, nsppol))
   ! band velocity in Cartesian coordinates in the IBZ
   ! Non-zero only for the kcalc k-points.



   integer :: bmin, bmax, bsize
   ! Only bands between bmin and bmax are considered in the integrals
   ! as we don't compute linewidths for all states.
   ! bmin = minval(%bstart_ks); bmax = maxval(%bstop_ks)
   ! bisze = bmax - bmin + 1


   type(ebands_t) :: ebands
   ! bandstructure object used to compute the transport properties
   ! Allocate using only the relevant bands for transport
   ! including valence states to allow to compute different doping


  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.


  real(dp) :: tolcum

  integer :: debug

  integer,allocatable :: bstart_ks(:,:)
   ! bstart_ks(nkcalc, nsppol)
   ! Initial KS band index included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  !integer,allocatable :: bstop_ks(:,:)
   ! bstop_ks(nkcalc, nsppol)

  integer,allocatable :: nbcalc_ks(:,:)
   ! nbcalc_ks(nkcalc, nsppol)
   ! Number of bands included in self-energy matrix elements for each k-point in kcalc.
   ! Depends on spin because all denerate states should be included when symmetries are used.

  integer,allocatable :: coords_kws(:,:)
   ! (2, self%nsppol))
   ! Cartesian coordinates of this processor in the Cartesian grid.

  integer,allocatable :: my_spins(:)
   ! my_spins(my_nspins)
   ! Indirect table giving the spin indices treated by this rank.
   ! Used only the collinear case with nspinor == 1

  integer,allocatable :: my_ikcalc(:)
   ! my_ikcalc(my_nkcalc)
   ! List of ikcalc indices treated by this pool if k-point parallelism is activated.

  real(dp),allocatable :: time_mesh(:,:,:,:,:)
  ! time_mesh(nwr, max_nbcalc, my_nkcalc, nsppol)
  ! Time mesh in atomic units

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.
   ! This array is not MPI-distributed.

  real(dp),allocatable :: kTmesh(:)
   ! kTmesh(ntemp)
   ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
   ! mu_e(ntemp)
   ! chemical potential of electrons for the different temperatures.


   real(dp),allocatable :: l0(:,:,:,:,:), l1(:,:,:,:,:), l2(:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp)
   ! Onsager coeficients in Cartesian coordinates



  integer,allocatable :: kcalc2ibz(:,:)
   !kcalc2ibz(nkcalc, 6))
   ! Mapping ikcalc --> IBZ as reported by listkk.


  real(dp),allocatable :: dw_vals(:,:,:,:)
   !  dw_vals(ntemp, max_nbcalc, my_nkcalc, nsppol)
   !  Debye-Waller term (static).


  real(dp),allocatable :: e0vals(:,:,:)
   ! (max_nbcalc, my_nkcalc, nsppol)
   ! KS energies at the calculated states retrieved from EPH output file

  real(dp),allocatable :: wrmesh_b(:,:,:,:)
  real(dp),allocatable :: wrmesh_ce(:,:,:,:)
  ! wrmesh_b(nwr, max_nbcalc, my_nkcalc, nsppol))
  ! Frequency mesh along the real axis (Ha units) used for the different bands
  ! Each mesh is **centered** on the corresponding KS energy.

  complex(dpc),allocatable :: vals_e0ks(:,:,:,:)
   ! vals_e0ks(ntemp, max_nbcalc, my_nkcalc, nsppol))
   ! Sigma_eph(omega=eKS, kT, band, ikcalc, spin).
   ! Fan-Migdal + Debye-Waller

  complex(dpc),allocatable :: vals_wr(:,:,:,:,:)
   ! vals_wr(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Sigma_eph(omega, kT, band, ikcalc, spin).
   ! enk_KS corresponds to nwr/2 + 1.

     complex(dpc),allocatable :: ct_vals(:,:,:,:,:)
   ! ct_vals(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Cumulant function (time, kT, band, ikcalc, spin).

     complex(dpc),allocatable :: c1(:,:,:,:,:)
   ! ct_vals(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Cumulant function (time, kT, band, ikcalc, spin).

     complex(dpc),allocatable :: c2(:,:,:,:,:)
   ! ct_vals(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Cumulant function (time, kT, band, ikcalc, spin).

     complex(dpc),allocatable :: c3(:,:,:,:,:)
   ! ct_vals(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Cumulant function (time, kT, band, ikcalc, spin).


     complex(dpc),allocatable :: gt_vals(:,:,:,:,:)
   ! vals_wr(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Green's function in time domain (time, kT, band, ikcalc, spin).

     complex(dpc),allocatable :: gw_vals(:,:,:,:,:)
   ! gw_vals(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Green's function in frequency domain(omega, kT, band) for given (ikcalc, spin).

   ! real(dp),allocatable :: ce_spfunc_wr(:,:,:,:,:)
   ! ce_spfunc_wr(nwr, ntemp, max_nbcalc, my_nkcalc, nsppol)
   ! Absorption spectrum (omega, kT, band, ikcalc, spin).


   real(dp),allocatable :: seebeck(:,:,:,:,:)
   real(dp),allocatable :: kappa(:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp)
   ! Transport coefficients in Cartesian coordinates


   real(dp),allocatable :: conductivity_mu(:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp)
   ! Conductivity in Siemens * cm-1
   ! computed by summing over k-points rather that by performing an energy integration).

   real(dp),allocatable :: print_dfdw(:,:)

   real(dp),allocatable :: mobility_mu(:,:,:,:,:)
   ! (3, 3, 2, nsppol, ntemp)
   ! mobility for electrons and holes (third dimension) at transport_mu_e(ntemp)
   ! Third dimension is for electron/hole


   real(dp),allocatable :: transport_mu_e(:)
   ! (%ntemp)
   ! Chemical potential at this carrier concentration and temperature


   real(dp) :: transport_fermie
   ! Fermi level specified in the input file when computing the SIGEPH file.


   real(dp) :: transport_extrael
   ! Extra electrons per unit cell specified in the input file when computing the SIGEPH file.


   real(dp),allocatable :: n_ehst(:,:,:)
    ! (2, %nsppol, %ntemp)
    ! Number of electrons (e) and holes (h) per unit cell
    ! The first dimension is for electrons/holes.
    ! If nsppol == 2, the second dimension is the number of e/h for spin else the total number of e/h summed over spins.


   real(dp) :: eph_extrael
   ! Extra electrons per unit cell used to compute SERTA lifetimes in sigmaph.

   real(dp) :: eph_fermie
   ! Fermi level specified in the input file when computing the SIGEPH file.

   integer         :: ce_ngfft(18)
   integer         :: ce_ngfft_g(18)
   type(mpi_type) :: ce_mpi_enreg



 contains

    procedure :: init => cumulant_init
    procedure :: compute => cumulant_compute
    procedure :: kubo_transport => cumulant_kubo_transport
    procedure :: sigmaph_ncread => cumulant_sigmaph_ncread
    procedure :: ncwrite => cumulant_ncwrite
    !procedure :: print_txt_files => cumulant_print_txt_files
    procedure :: free => cumulant_free

 end type cumulant_t

!!***

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_driver
!! NAME
!! cumulant_driver
!!
!! FUNCTION
!! General driver to compute transport properties
!!
!! INPUTS
!! dtfil<datafiles_type>=variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! cryst<crystal_t>=Crystalline structure
!! comm=Initial MPI communicator with all procs provided by user.
!!
!! PARENTS
!!      m_eph_driver
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_driver(dtfil, dtset, ebands, cryst, comm)

!Arguments ------------------------------------
!scalars
 integer, intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands

!Local variables ------------------------------
 integer,parameter :: master = 0
 integer :: my_rank, ierr
#ifdef HAVE_NETCDF
 integer :: ncid
#endif
 character(len=1000) :: msg
 character(len=fnlen) :: path
 type(cumulant_t) :: cumulant
 type(sigmaph_t) :: sigmaph
!arrays
 !real(dp) :: extrael_fermie(2), sigma_erange(2)
 integer :: unts(2) !, sigma_ngkpt(3)

! *************************************************************************

 my_rank = xmpi_comm_rank(comm)
 unts = [std_out, ab_out]

 !sigeph_filepath = dtfil%filsigephin
 call wrtout(unts, ch10//' Entering cumulant expansion computation driver.')

 ! Initialize cumulant object
 call cumulant%init(dtset, dtfil, cryst, ebands, comm, sigmaph)
 !if (cumulant%all_comm%me == -1) goto 100

 ! Compute C(t), G(t), G(w) and A(w) using cumulant expansion
 call cumulant%compute()

 ! Compute transport properties within the Kubo formalism.
 if (any(abs(dtset%sigma_erange) > zero)) call cumulant%kubo_transport(dtset, cryst)!path, ncid)

 ! Output netcdf file with cumulant results e.g. A_nk(w).
 ! Use EPH prefix because we may also have EE and PHE cumulants.
 path = strcat(dtfil%filnam_ds(4), "_EPH_CUMULANT.nc")
 call cumulant%ncwrite(path, cryst, ebands)

 !if (my_rank == master) then
 !  ! Print cumulant expansion results to stdout and other external txt files (for the test suite)
 !  !call cumulant%print_txt_files(cryst, dtset, dtfil)
 !end if

 ! Free memory
100 call cumulant%free()

end subroutine cumulant_driver
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_init
!! NAME
!! cumulant_init
!!
!! FUNCTION
!!   Initialization of the cumulant, ex: time mesh
!!
!! INPUTS
!! dtset<dataset_type>=All input variables for this dataset.
!! dtfil<datafiles_type>=variables related to files.
!! comm=Initial MPI communicator with all procs provided by user.
!!
!! PARENTS
!!      m_cumulant
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_init(self, dtset, dtfil, cryst, ebands, comm, sigmaph )

!Arguments --------------------------------------
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(sigmaph_t), intent(out) :: sigmaph
 class(cumulant_t),intent(inout) :: self

!Local variables --------------------------------
 integer, parameter :: master = 0
 integer :: ncerr, ncid, my_rank, nprocs, ierr, spin, ikcalc, ib, cnt, nbands, color !my_spin, my_kcalc,
 real(dp) :: cpu, wall, gflops, rsize
 logical :: is_prime
 character(len=500) :: msg
 character(len=fnlen) :: sigeph_filepath
 integer :: unts(2), facts(2)
 real(dp), allocatable :: rtmp_vals_wr(:,:,:,:,:,:), rtmp_vals_e0ks(:,:,:,:,:)
 type(ebands_t) :: tmp_ebands
 real(dp) :: extrael_fermie(2),sigma_erange(2)
 integer :: sigma_ngkpt(3)
#ifdef HAVE_MPI
 integer :: ndims, comm_cart, me_cart
 logical :: reorder
 integer,allocatable :: dims(:)
 logical,allocatable :: periods(:), keepdim(:)
#endif

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 self%comm = comm
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)


 unts = [std_out, ab_out]

 ! Read data from out_SIGEPH.out file
 sigeph_filepath = dtfil%filsigephin
 call wrtout(unts, sjoin("- Reading Sigma results from:", sigeph_filepath), newlines=1, do_flush=.True.)
 call self%sigmaph_ncread(sigeph_filepath, ncid, comm)
 sigmaph = sigmaph_read(sigeph_filepath, dtset, comm, msg, ierr, keep_open = .false., &
         extrael_fermie=extrael_fermie, sigma_ngkpt=sigma_ngkpt, sigma_erange=sigma_erange)

 self%bmin = minval(sigmaph%bstart_ks); self%bmax = maxval(sigmaph%bstop_ks)
 ABI_CALLOC(self%vbks, (3, self%bmin:self%bmax, ebands%nkpt, sigmaph%nsppol))

if (any(abs(dtset%sigma_erange) > zero)) tmp_ebands = sigmaph%get_ebands(cryst, ebands, [self%bmin, self%bmax], &
                                   self%kcalc2ebands, self%linewidths, self%vbks, xmpi_comm_self)
 self%eph_extrael = extrael_fermie(1)
 self%eph_fermie = extrael_fermie(2)

 ! Possibility to increase nwr in case of interpolation in cumulant to add extra points
 ! TODO: not working yet
 self%nwr_ce = self%nwr!*4 -1 ! Odd

 self%ieta = j_dpc * sigmaph%ieta
 self%ebands = tmp_ebands
 self%mu_e = sigmaph%mu_e
  
 ! Setting variables to launch FFT calculations later on
 call ngfft_seq(self%ce_ngfft, [self%nwr, 1, 1])
 self%ce_ngfft(4:6) = self%ce_ngfft(1:3)

 call initmpi_seq(self%ce_mpi_enreg)
 call init_distribfft_seq(self%ce_mpi_enreg%distribfft, 'c', self%ce_ngfft(2), self%ce_ngfft(3), 'all')
 call init_distribfft_seq(self%ce_mpi_enreg%distribfft, 'f', self%ce_ngfft(2), self%ce_ngfft(3), 'all')


 call ngfft_seq(self%ce_ngfft_g, [self%nwr_ce, 1, 1])
 self%ce_ngfft_g(4:6) = self%ce_ngfft_g(1:3)


 ! Setting debugging ( higher verbosity )
 self%tolcum = dtset%tolcum
 if (self%tolcum > zero) then
   self%debug = 0
 else
   self%debug = 1
 endif
 self%tolcum = abs(self%tolcum)

 call wrtout(unts, " Cumulant parameters:")
 call wrtout(unts, sjoin(" Number of spins:                   ", itoa(self%nsppol)))
 call wrtout(unts, sjoin(" Number of spinor components:       ", itoa(self%nspinor)))
 call wrtout(unts, sjoin(" Number of k-points computed:       ", itoa(self%nkcalc)))
 call wrtout(unts, sjoin(" Maximum number of bands computed:  ", itoa(self%max_nbcalc)))
 call wrtout(unts, sjoin(" Number of frequencies in Sigma(w): ", itoa(self%nwr)))
 call wrtout(unts, sjoin(" Number of frequencies in G(w): ", itoa(self%nwr_ce)))
 call wrtout(unts, sjoin(" Number of Temperatures:            ", itoa(self%ntemp)))

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! At present, there are two levels of MPI parallelism: kcalc and time/frequencies.
 ! Only the kcalc parallelism leads to memory distribution (my_nkcalc) whereas the secon level is only
 ! used to distribute the loop over frequencies/time points.
 ! Note the following:
 !
 ! 1) We distribute contigous blocks of kcalc k-points so that we can read/write data in mykcalc chuncks.

 ! TODO:
 !  *) Spin and nkcalc(nsppol)
 !  *) Decide whther it makes sense to create a new comm if nprocs > nkcalc  and mod(nprocs, nkcalc) /= 0
 !  *) No need anymore for parallelism in time/frequency


 self%kcalc_comm%nproc = nprocs
! self%wt_comm%nproc = nprocs / self%kcalc_comm%nproc

 if (nprocs > self%nkcalc) then
   call ifact2(nprocs, facts, is_prime)
   if (is_prime) then
     ! TODO
     ABI_ERROR("prime numbers are evil!")
   end if
   if (facts(1) > self%nkcalc) then
      ! Avoid procs with no k-points.
      facts(1:2) = facts(2:1:-1)
   end if
   self%kcalc_comm%nproc = facts(1)
 !  self%wt_comm%nproc = facts(2)
   !do cnt=self%nkcalc,1,-1
   !  if (mod(nprocs, cnt) == 0 .and. mod(self%nkcalc, cnt) == 0) then
   !    self%kcalc_comm%nproc = cnt; self%my_nkcalc = self%nkcalc / cnt; exit
   !  end if
   !end do
 end if

 ABI_ICALLOC(self%coords_kws, (2, self%nsppol))
#ifdef HAVE_MPI
 ! Create 2d cartesian communicator: kpoints in Sigma_k, (time/frequency)
 ! FIXME: Fix spin
 spin = 1
 ndims = 1
 ABI_MALLOC(dims, (ndims))
 ABI_MALLOC(periods, (ndims))
 ABI_MALLOC(keepdim, (ndims))
 periods(:) = .False.; reorder = .False.
! dims = [self%kcalc_comm%nproc, self%wt_comm%nproc]
 dims = [self%kcalc_comm%nproc]

 call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm_cart, ierr)
 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, self%coords_kws(:, spin), ierr)

 ! Create communicator for kpoints.
 keepdim = .False.; keepdim(1) = .True.
 call MPI_CART_SUB(comm_cart, keepdim, self%kcalc_comm%value, ierr); self%kcalc_comm%me = xmpi_comm_rank(self%kcalc_comm%value)

 ! Create communicator for w/t.
! keepdim = .False.; keepdim(2) = .True.
! call MPI_CART_SUB(comm_cart, keepdim, self%wt_comm%value, ierr); self%wt_comm%me = xmpi_comm_rank(self%wt_comm%value)

 ! Create communicator for spins.
 !keepdim = .False.; keepdim(5) = .True.
 !call MPI_CART_SUB(comm_cart, keepdim, self%spin_comm%value, ierr); self%spin_comm%me = xmpi_comm_rank(self%spin_comm%value)

 ! Create communicator for parallel IO
 ! If parallelism over time/frequencies is ON, we select only a column of procs (kcalc dimension)
 ! to avoid writing the same data multiple times as omega/time are not MPI distributed
 ! Obviously I'm assuming HDF5 + MPI-IO
 color = xmpi_undefined; if (self%coords_kws(2, spin) == 0) color = 1
 call xmpi_comm_split(comm, color, my_rank, self%ncwrite_comm%value, ierr)
 if (color == 1) then
   self%ncwrite_comm%me = xmpi_comm_rank(self%ncwrite_comm%value)
   self%ncwrite_comm%nproc = xmpi_comm_size(self%ncwrite_comm%value)
 else
   call self%ncwrite_comm%set_to_null()
 end if

 ABI_FREE(dims)
 ABI_FREE(periods)
 ABI_FREE(keepdim)
 call xmpi_comm_free(comm_cart)
#endif

 call wrtout(unts, &
   sjoin("- Using: ", itoa(self%kcalc_comm%nproc), "MPI procs to distribute:", itoa(self%nkcalc), "k-points."))
 !call wrtout(unts, &
!   sjoin("- Using: ", itoa(self%wt_comm%nproc), "MPI procs to parallelize:", itoa(self%nwr_ce), "frequency/time loops."))
 call wrtout(std_out, sjoin("- Performing parallel IO with:", itoa(self%ncwrite_comm%nproc), "procs"))

 if (self%ncwrite_comm%nproc /= 1 .and. .not. nctk_has_mpiio) then
   ABI_WARNING("Cumulant requires netcdf with MPI-IO support! Continuing anyway but runtime error is expected!")
 end if

 ! Consistency check.
 !if (self%kcalc_comm%nproc * self%wt_comm%nproc /= nprocs) then
 if (self%kcalc_comm%nproc /= nprocs) then
   write(msg, "(a,i0,3a, 3(a,1x,i0))") &
     "Cannot create 2d Cartesian grid with total nprocs: ", nprocs, ch10, &
     "Idle processes are not supported. The product of comm%nproc should be equal to nprocs.", ch10, &
     "kcalc_nproc (", self%kcalc_comm%nproc, ") != ", & ! x wt_nproc (", self%wt_comm%nproc, ") != ", &
     self%kcalc_comm%nproc ! * self%wt_comm%nproc
   ABI_ERROR(msg)
 end if

 ! Here we distribute the k-points inside kcalc_comm.
 call xmpi_split_block(self%nkcalc, self%kcalc_comm%value, self%my_nkcalc, self%my_ikcalc)
 ABI_CHECK(self%my_nkcalc > 0, sjoin("nkcalc (", itoa(self%nkcalc), ") > nprocs (", itoa(nprocs), ")"))

 ib = mod(self%nkcalc, self%my_nkcalc)
 call xmpi_sum(ib, self%comm, ierr)
 if (ib /= 0) then
   ABI_COMMENT("The number of MPI procs should be divisible by nkcalc to distribute memory equally!")
 end if

 ! Distribute spins and create mapping to spin index.
 if (self%nsppol == 2) then
   ABI_ERROR("cumulant with nsppol 2 is not supported.")
   call xmpi_split_block(self%nsppol, comm, self%my_nspins, self%my_spins)
   !ABI_CHECK(self%my_nspins > 0, sjoin("nsppol (", itoa(self%nsppol), ") > spin_comm_nproc (", itoa(comm), ")"))
 else
   ! No nsppol parallelism DOH!
   self%my_nspins = 1
   ABI_MALLOC(self%my_spins, (self%my_nspins))
   self%my_spins = 1
 end if

 !ABI_CALLOC(self%ce_spfunc_wr, (self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_CALLOC(self%gw_vals, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))

 if (self%debug == 1) then
   ABI_CALLOC(self%time_mesh, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
   ABI_CALLOC(self%ct_vals, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
   ABI_CALLOC(self%c1, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
   ABI_CALLOC(self%c2, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
   ABI_CALLOC(self%c3, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
   ABI_CALLOC(self%gt_vals, (self%nwr_ce, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 endif

 ABI_MALLOC(rtmp_vals_wr, (2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(rtmp_vals_e0ks, (2, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%vals_wr, (self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%vals_e0ks, (self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%wrmesh_b, (self%nwr, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%wrmesh_ce, (self%nwr_ce, self%max_nbcalc, self%my_nkcalc, self%nsppol))
! ABI_MALLOC(self%wrmesh_ce, (self%nwr_ce, self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%e0vals, (self%max_nbcalc, self%my_nkcalc, self%nsppol))
 ABI_MALLOC(self%dw_vals, (self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol))

 ! Estimate memory (only nwr arrays)
 rsize = five; if (self%debug == 1) rsize = rsize + 3
 rsize = rsize * self%nwr_ce * self%ntemp * self%max_nbcalc * self%my_nkcalc * self%nsppol
 call xmpi_max(rsize, self%comm, ierr)
 write(msg,'(a,f8.1,a)')' Memory needed for cumulant arrays per MPI proc: ',dp * rsize * b2Mb, ' [Mb] <<< MEM'
 call wrtout(std_out, msg)

 ! Each MPI proc reads a chunk of my_nkcalc entries starting at (ikcalc, spin) from the SIGEPH.nc file.
 ! so that memory is MPI-distributed.
 ! On disk, we have the following netcdf arrays:
 !
 !         nctkarr_t("wrmesh_b", "dp", "nwr, max_nbcalc, nkcalc, nsppol"), &
 !         nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
 !         nctkarr_t("vals_e0ks", "dp", "two, ntemp, max_nbcalc, nkcalc, nsppol"), &
 !         nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
 !
 ! Note the nkcalc global dimension instead of my_nkcalc.
 ! Results are in a.u.

 ! Offse of this proc
 spin = self%my_spins(1)
 ikcalc = self%my_ikcalc(1)

#ifdef HAVE_NETCDF
 ncerr = nf90_get_var(ncid, vid("vals_wr"), rtmp_vals_wr, start=[1,1,1,1,ikcalc,spin], &
                      count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%my_nspins])
 NCF_CHECK(ncerr)

 ncerr = nf90_get_var(ncid, vid("wrmesh_b"), self%wrmesh_b, start=[1,1,ikcalc,spin], &
                      count=[self%nwr, self%max_nbcalc, self%my_nkcalc, self%my_nspins])
 NCF_CHECK(ncerr)

 ncerr = nf90_get_var(ncid, vid("vals_e0ks"), rtmp_vals_e0ks, start=[1,1,1,ikcalc,spin],  &
                      count=[2, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%my_nspins])
 NCF_CHECK(ncerr)

 ncerr = nf90_get_var(ncid, vid("ks_enes"), self%e0vals, start=[1,ikcalc,spin], &
                      count=[self%max_nbcalc, self%my_nkcalc, self%my_nspins])
 NCF_CHECK(ncerr)

 ncerr = nf90_get_var(ncid, vid("dw_vals"), self%dw_vals, start=[1,1,ikcalc,spin], &
                      count=[self%ntemp, self%max_nbcalc, self%my_nkcalc, self%my_nspins])
 NCF_CHECK(ncerr)

 ! Close the file here but the Kubo equation needs to read v_nk from file.
 ! We will have to rationalize this part.
 NCF_CHECK(nf90_close(ncid))
#endif

 ! Convert from real arrays to Fortran complex
 self%vals_wr =  (rtmp_vals_wr(1,:,:,:,:,:) + j_dpc * rtmp_vals_wr(2,:,:,:,:,:))!*Ha_eV
 self%vals_e0ks = (rtmp_vals_e0ks(1,:,:,:,:) + j_dpc * rtmp_vals_e0ks(2,:,:,:,:))!*Ha_eV

! self%wrmesh_b = self%wrmesh_b*Ha_eV
! self%e0vals = self%e0vals*Ha_eV
! self%max_nbcalc=self%max_nbcalc*Ha_eV

! self%e0vals = self%e0vals * Ha_eV
! self%wrmesh_b = self%wrmesh_b * Ha_eV

 ABI_SFREE(rtmp_vals_wr)
 ABI_SFREE(rtmp_vals_e0ks)

 call cwtime_report(" cumulant init", cpu, wall, gflops)

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
end function vid

end subroutine cumulant_init
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_compute
!! NAME
!! cumulant_compute
!!
!! FUNCTION
!!  Compute cumulant.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_compute(self)

!Arguments ------------------------------------
 class(cumulant_t),intent(inout) :: self
!Local variables ------------------------------
 integer,parameter :: master = 0
 integer :: nbands, ib, ikcalc, it, iw, spin, itemp, comm
 integer :: nwr, nwr_ce
 integer :: my_rank, nprocs, ierr, my_spin, my_ik
 real(dp) :: time, omega, init_t, time_step, time_step_init, time_max, init_w, end_w, wmesh_step_init
 real(dp) :: cpu, wall, gflops, cpu_kloop, wall_kloop, gflops_kloop
 real(dp) :: wr_step, wr_step_ce, Ha_fs
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: temp_g(:,:,:), temp_r(:,:), temp_r_cplx(:,:), temp_g_ce(:,:,:)
 real(dp),allocatable :: betaoverw2(:)
 complex(dpc),allocatable :: betaoverw2c(:), temp_reflex(:)
 real(dp),allocatable :: wrmesh_shifted(:), wrmesh_shifted_ce(:), beta(:), inv_wrmesh_shifted_sq(:), c3(:)
 real(dp),allocatable :: time_mesh(:), time_mesh_temp(:)
 real(dp) :: output_value, output_c3, output_test2r, output_test2i
 real(dp) :: m_fit_re, b_fit_re, m_fit_im, b_fit_im, res_re, res_im
 complex(dpc),allocatable :: c1(:), ct_temp(:), c_temp(:)
 complex(dpc),allocatable :: c2(:), ct(:), gt(:), gw(:), g1(:)
 complex(dpc) :: output_test, output_test2
 Integer :: sts, itest
 integer :: fftalg
 logical :: use_fft

!************************************************************************

! Initialization of parallel variables
 comm = self%comm
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call wrtout(std_out, " Computing cumulant. This may take some time depending on the dims of the problem...")
 call cwtime(cpu_kloop, wall_kloop, gflops_kloop, "start")

 ! Initialization of variables and allocation of arrays
 nwr = self%nwr
 nwr_ce = self%nwr_ce
 ABI_CALLOC(c1, (nwr))
 ABI_CALLOC(c_temp, (nwr_ce))
 ABI_CALLOC(temp_reflex, (nwr/2))
 ABI_CALLOC(c2, (nwr))
 ABI_CALLOC(c3, (nwr))
 ABI_CALLOC(g1, (nwr_ce))
 ABI_CALLOC(gw, (nwr_ce))
 ABI_CALLOC(time_mesh, (nwr_ce))
 ABI_CALLOC(time_mesh_temp, (nwr))
 ABI_CALLOC(ct, (nwr_ce))
 ABI_CALLOC(ct_temp, (nwr))
 ABI_CALLOC(gt, (nwr_ce))
 ABI_CALLOC(temp_g, (2,nwr,1))
 ABI_CALLOC(temp_g_ce, (2,nwr_ce,1))

 ABI_CALLOC(temp_r, (nwr,1))
 ABI_CALLOC(temp_r_cplx, (2*nwr_ce,1))
 !ABI_CALLOC(inv_wrmesh_shifted_sq, (nwr_ce))
 ABI_CALLOC(wrmesh_shifted, (nwr))
 ABI_CALLOC(wrmesh_shifted_ce, (nwr_ce))
 ABI_CALLOC(beta, (nwr))
 ABI_CALLOC(betaoverw2, (nwr))
 ABI_CALLOC(betaoverw2c, (nwr))
 Ha_fs = 8.955433106


 fftalg = self%ce_ngfft(7)
 ! Loops are MPI-parallelized over k-points and spin.
 do my_spin=1,self%my_nspins
   spin = self%my_spins(my_spin)
   do my_ik=1,self%my_nkcalc
     ikcalc= self%my_ikcalc(my_ik)
     !if (ikcalc /= 1) cycle ! For debugging purpose
     nbands = self%nbcalc_ks(ikcalc, spin)
     call cwtime(cpu, wall, gflops, "start")

     do ib=1,nbands
       !if (ib > 1) cycle ! MG HACK To be able to run tests quickly.

       ! Shifting w-mesh to center at KS-energy and remove the value w=0.0 from integrations ( Principal Value )
       wr_step = self%wrmesh_b(2,ib,my_ik,spin) - self%wrmesh_b(1,ib,my_ik,spin)
       wrmesh_shifted(:) = self%wrmesh_b(:,ib,my_ik,spin) - self%e0vals(ib,my_ik,spin) &
               - 0.5 * wr_step


       do itemp=1,self%ntemp
        ! if (itemp > 1) cycle ! MG HACK To be able to run tests quickly.

         ! Defining time mesh
         ! all variables with nwr_ce size are to use after interpolation of the
         ! cumulant function in case we want to add more points
         time_max = (log(self%tolcum) / ( - abs(aimag(self%vals_e0ks(itemp, ib, my_ik,spin)) )) )

         init_t = 0.0
         time_step_init = 1.0/(nwr - 1.0)
         time_mesh_temp(:) = arth(init_t, time_step_init, nwr) 
         time_mesh_temp = time_mesh_temp * 2*PI/wr_step
         time_step = time_mesh_temp(2) - time_mesh_temp(1)

         time_mesh(:) = time_mesh_temp(:) !! arth(init_t, time_step_init, nwr_ce) 
         time_mesh = time_mesh * 2*PI/wr_step

         ! Defining frequency mesh for after interpolation
         
         !! TODO: Fix problems from when adding more points to time mesh by
         !interpolation

         !!wmesh_step_init = 1.0/(nwr_ce/2.0 -1 )
         !!init_w = -1.0
         !!end_w = 0.0
         !!wrmesh_shifted_ce(:nwr_ce/2+1) = linspace(init_w,end_w, nwr_ce/2+1)
         !!init_w = 0.0
         !!end_w = 1.0
         !!wrmesh_shifted_ce(nwr_ce/2+1:) = linspace(init_w,end_w, nwr_ce/2+1)
 
         !!wr_step_ce = wrmesh_shifted_ce(2) - wrmesh_shifted_ce(1)
         !!wrmesh_shifted_ce(:) = wrmesh_shifted_ce(:) * PI/time_step - 0.5 *wr_step_ce
         !!self%wrmesh_ce(:,ib,my_ik,spin) = wrmesh_shifted_ce(:) + self%e0vals(ib,my_ik,spin) + 0.5 *wr_step_ce
         !!wrmesh_shifted_ce = wrmesh_shifted

         if (self%debug == 1) self%time_mesh(:,itemp,ib,my_ik,spin) = time_mesh_temp(:)! * 0.24188845385 *10E-4 ! picoseconds ( I think )

         if (time_mesh_temp(nwr) < time_max) ABI_WARNING(sjoin("KBT",itoa(my_ik),itoa(ib),itoa(itemp)))
         if (time_mesh_temp(nwr) < time_max) ABI_WARNING(sjoin("Time mesh smaller than needed to reach tolerance value. Actual value of", ftoa(time_mesh_temp(nwr)), ", desirable ",ftoa(time_max),". Increase nfreqsp."))
         beta(:) = abs( aimag( self%vals_wr(:, itemp, ib, my_ik, spin ) ) ) / pi

         ! Calculation of the different terms of the cumulant ( c1, c2 and c3 ),
         betaoverw2(:) = beta(:) / (wrmesh_shifted(:) ** 2)

         c2(:) = -j_dpc * real( self%vals_e0ks( itemp, ib, my_ik, spin ) ) * time_mesh_temp(:)

         output_c3 = simpson(wr_step, betaoverw2)
         c3(:) = -1.0 * output_c3


         temp_r(:,1) = betaoverw2(:)

         use_fft = .True.       
         if (use_fft) then
           call fourdp(1, temp_g, temp_r, -1, self%ce_mpi_enreg, nwr, 1, self%ce_ngfft , 0 )
           c1(:) = temp_g(1, :, 1) + j_dpc* temp_g(2, :, 1)
           c1 = c1 * wr_step * nwr
           c1(2::2) = -1.0 * c1(2::2)
         else
           ABI_COMMENT("FFT not available, using DFT instead. Slower but reliable. Parallelism over frequency for integration of C(t)")
           do it=1, nwr
!            if (mod(it, self%wt_comm%nproc) /= self%wt_comm%me) cycle  ! MPI parallelism over time
             time = time_mesh_temp(it)
               
             ! Discrete Fourier Transform of C1 ( cumulant function without the +iwt and -1 parts
             c_temp(:) = betaoverw2(:) *  exp( - j_dpc * time * wrmesh_shifted(:) )
             c1(it) = simpson_cplx( nwr, wr_step, c_temp)
           enddo
             
         endif


         ! The cumulant function as sum of the different components
         ct_temp(:) = c1 + c2 + c3
         ! Fitting the cumulant function
    	 res_re = linfit(nwr,time_mesh_temp(:), real(ct_temp(:)), m_fit_re, b_fit_re)
    	 res_im = linfit(nwr,time_mesh_temp(:), aimag(ct_temp(:)), m_fit_im, b_fit_im)
!         call xmpi_sum(ct, self%wt_comm%value, ierr)
         if (self%debug == 1) then ! .and. self%wt_comm%nproc > 1) then
                self%ct_vals(:, itemp, ib, my_ik, spin) = ct(:)
                self%c1(:, itemp, ib, my_ik, spin) = c1(:)
                self%c2(:, itemp, ib, my_ik, spin) = c2(:)
                self%c3(:, itemp, ib, my_ik, spin) = c3(:)
         endif
	    
         ! Adding extra interpolated points to the cumulant function
    	 !!ct(:nwr/2) = ct_temp(:nwr/2)
    	 !!ct(nwr/2:) = time_mesh(nwr/2:) * m_fit_re + b_fit_re + j_dpc * ( time_mesh(nwr/2:) * m_fit_im + b_fit_im )

    	 !!do iw=1, nwr/2
	     !!    temp_reflex(nwr/2-iw) = - ct(iw) + time_mesh(iw) * m_fit_re + b_fit_re + j_dpc * ( time_mesh(iw) * m_fit_im + b_fit_im )
    	 !!enddo
	     !!ct(nwr_ce - nwr/2:) =  ct(nwr_ce - nwr/2:) + temp_reflex(:)
    	 ct(:) = ct_temp(:)
         ! Retarded Green's function in time domain
         gt = - j_dpc * exp( ct )

         ! Collect data if wt parallelism.
!         call xmpi_sum(gt, self%wt_comm%value, ierr)

         if (self%debug == 1) then ! .and. self%wt_comm%nproc > 1) then
           self%gt_vals(:, itemp, ib, my_ik, spin) = gt(:)
         end if

         use_fft = .True.
         if (use_fft) then
           
           ! Fast Fourier Transform to obtain the Green's function in frequency
           ! domain
           temp_g_ce(1,:,1) = real(gt(:))
           temp_g_ce(2,:,1) = aimag(gt(:))
           call fourdp(2, temp_g_ce, temp_r_cplx, 1, self%ce_mpi_enreg, nwr_ce, 1, self%ce_ngfft_g , 0 )
           gw(:) = temp_r_cplx(1::2,1) + j_dpc* temp_r_cplx(2::2,1)

           ! TODO use ig2gfft from 52_fft_mpi_noabirule/m_fftcore.F90 instead of the two following lines
           self%gw_vals(:int(nwr_ce/2.0), itemp, ib, my_ik, spin) = gw(int(nwr_ce/2.0):)
           self%gw_vals(int(nwr_ce/2.0):, itemp, ib, my_ik, spin) = gw(:int(nwr_ce/2.0))

 
           self%gw_vals(:, itemp, ib, my_ik, spin) = self%gw_vals(:, itemp, ib, my_ik, spin) * time_step

           ! FFT is different from integration methods and the two extreme points
           ! need to be compensated
           self%gw_vals(:, itemp, ib, my_ik, spin) = self%gw_vals(:, itemp, ib, my_ik, spin)  & 
                 - 0.5* gt(1)*time_step - 0.5 * gt(nwr_ce)*exp(j_dpc*wrmesh_shifted_ce(:)*time_mesh(nwr_ce)) * time_step


         else
           ABI_COMMENT("FFT not available, using DFT instead. Slower but reliable. Parallelism over time for integration of G(t)")

           do iw=1, nwr_ce
!             if (mod(iw, self%wt_comm%nproc) /= self%wt_comm%me) cycle  ! MPI parallelism over freqs
             omega = wrmesh_shifted_ce(iw)
         
             ! Discrete Fourier Transform of the Green's function
             g1(:) = exp( j_dpc * omega * time_mesh(:) ) * gt(:)
         
             ! Retarded Green's function in frequency domain
             self%gw_vals(iw, itemp, ib, my_ik, spin) = simpson_cplx(nwr_ce, time_step, g1)
         
             !if (my_rank == 0) write(ab_out, *)"gw_vals",  self%gw_vals(iw, itemp, ib, my_ik, spin)
      

           end do ! iw
         end if
         
         

         ! Collect data if wt parallelism.
!         call xmpi_sum(self%gw_vals(:, itemp, ib, my_ik, spin) , self%wt_comm%value, ierr)

       end do  ! itemp
     end do ! ib

     write(msg,'(4(a,i0),a,f8.2)') " k-point [", ikcalc, "/", self%nkcalc, "]"
     call cwtime_report(msg, cpu, wall, gflops)
   end do ! my_ik
 end do ! my_spin

 
 ABI_SFREE(c1)
 ABI_SFREE(c2)
 ABI_SFREE(c3)
 ABI_SFREE(ct)
 ABI_SFREE(gt)
 ABI_SFREE(g1)
 ABI_SFREE(gw)
 ABI_SFREE(beta)
 ABI_SFREE(temp_g)
 ABI_SFREE(temp_g_ce)
 ABI_SFREE(temp_r)
 ABI_SFREE(temp_r_cplx)
 ABI_SFREE(time_mesh)
 ABI_SFREE(wrmesh_shifted)
 ABI_SFREE(wrmesh_shifted_ce)
 !ABI_SFREE(inv_wrmesh_shifted_sq)
 ABI_SFREE(betaoverw2)

 call cwtime_report(" cumulant_compute", cpu_kloop, wall_kloop, gflops_kloop)

! contains

! real function lreg(x,y, n, nsteps) result(m,b)
! !
! ! Determines a linear regression from x and y values
! ! What are the best m and b values to produce y= m*x + b ?
! !
! !
! ! Cost function ( Root Mean Squared Error ): J = 1/n sum_i^n (pred_i - y_i)^2
! ! where pred is the predicted value and y the true value
! ! 
! ! Goal: minimize J
! ! How? Using Gradient Descent
! ! - Learning rate is the step that the new value will be ( too small, takes longer; too large, it can be instable )
! ! - Initial m or b are chosen randomly
! ! - n is the number of points
! !
! ! new b = old b - 2*(learning rate)/n sum_i^n (pred(x_i) - y) * x_i
! ! new m = old m  - 2*(learning rate)/n sum_i^n (pred(x_i) - y) 
! !
!
! integer, intent(in) :: n, nsteps
! real(dp), intent(in) :: x(n), y(n)
! !real(dp) :: m,b
! integer :: istep, i
! real(dp) :: lrate, cost, acc
! real(dp) :: pred(n)
!
! m = 0 ! Initial guesses
! b = 0
! do istep=1,nsteps
!   
!  pred(:) = m * x(:) + b ! Prediction with the new coefficients m, b
!
!  ! Check accuracy comparing the linear regression and the data
!  do i=1, n
!    if (abs(y(i) ) < 1e-4 ) cycle
!    acc = acc + sum( abs(pred(i) - y(i))/y(i) )
!  end do i
!  acc = 1 - acc
!  print *, "Accuracy: ",istep, acc
!
!  cost = 1.0/n * sum(pred(:) - y(:))**2
!
!  print *, "Cost: ", istep, cost
!
!  ! Update coefficients
!
!  m = m - 2.0*lrate/n * sum(pred(:) - y(:))
!  b = b - 2.0*lrate/n * sum((pred(:) - y(:))*x(:))
! enddo ! istep
!
! end function lreg

!  complex function trapz(f_size, f_step, f)
!
!   integer,intent(in) :: f_size
!   real(dp),intent(in) :: f_step
!   complex(dpc),intent(in) :: f(f_size)
!
!   trapz = ( sum(f) - 0.5* f(1) - 0.5* f(f_size) )* f_step
!
!  end function trapz

end subroutine cumulant_compute
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_kubo_transport
!! NAME
!! cumulant_kubo_transport
!!
!! FUNCTION
!!  Compute conductivity/mobility
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_kubo_transport(self, dtset, cryst)

!Arguments ------------------------------------
 class(cumulant_t),intent(inout) :: self
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
! type(ebands_t),intent(in) :: ebands

!Local variables ------------------------------
 integer,parameter :: master = 0
 integer :: cnt
 integer :: nbands, ib, ikcalc, it, iw, spin, itemp, comm, ik_ibz
 integer :: ieh, ib_eph
 integer :: ii, jj, time_opt, isym_k, trev_k, dfdw
 integer :: my_rank, nprocs, ierr, my_spin, my_ik
 real(dp) :: time, omega, init_t, time_step, time_max
 real(dp) :: cpu, wall, gflops, cpu_kloop, wall_kloop, gflops_kloop
 real(dp) :: eig_nk, sp_func, wr_step
 real(dp) :: integration, mu_e, max_occ, fact0, fact
 character(len=500) :: msg
!arrays
 real(dp),allocatable :: kernel(:), test_Aw(:), test_dfdw(:),dfdw_acc(:),Aw(:),Aw_l0(:),dfdw_l0(:)
 real(dp) :: int_Aw, int_dfdw
 real(dp) :: vr(3), vv_tens(3,3), S(3,3), wtk, spfunc
 real(dp), allocatable :: onsager_coeff(:,:)
 real(dp) :: work_33(3,3),l0inv_33nw(3,3,2)
 real(dp) :: Tkelv


!************************************************************************

 comm = self%comm
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 
 call wrtout(std_out, " Computing conductivity using Kubo-Greenwood method.")
 call cwtime(cpu_kloop, wall_kloop, gflops_kloop, "start")


 ABI_MALLOC(kernel, (self%nwr))
! ABI_MALLOC(test_Aw, (self%ntemp))
! ABI_MALLOC(test_dfdw, (self%ntemp))
 ABI_MALLOC(Aw, (self%nwr))
 ABI_MALLOC(dfdw_acc, (self%nwr))
 ABI_MALLOC(Aw_l0, (self%ntemp))
 ABI_MALLOC(dfdw_l0, (self%ntemp))

 ABI_CALLOC(self%l0, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%l1, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%l2, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%seebeck, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%kappa, (3, 3, 2, self%nsppol, self%ntemp))



 ABI_CALLOC(self%mobility_mu, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%conductivity_mu, (3, 3, 2, self%nsppol, self%ntemp))
 ABI_CALLOC(self%print_dfdw, (self%nwr, self%ntemp))

 ABI_MALLOC(self%transport_mu_e, (self%ntemp))

 ABI_CALLOC(self%n_ehst, (2, self%nsppol, self%ntemp))

 self%transport_fermie = dtset%eph_fermie
 self%transport_extrael = dtset%eph_extrael
 self%transport_mu_e = self%mu_e
 if (self%transport_fermie /= zero) self%transport_mu_e = self%transport_fermie

 if (self%transport_fermie == zero .and. self%transport_extrael /= self%eph_extrael) then

   if (self%transport_extrael /= self%eph_extrael) then
     write(msg,'(2(a,e18.8),3a)') &
       ' extrael from SIGEPH: ',self%transport_extrael, ' and input file: ',self%eph_extrael, "differ", ch10, &
       ' Will recompute the chemical potential'
     call wrtout(std_out, msg)
   end if

   ! Compute Fermi level for different T values.
   call ebands_get_muT_with_fd(self%ebands, self%ntemp, self%kTmesh, dtset%spinmagntarget, dtset%prtvol, self%transport_mu_e, comm)
 end if


 call ebands_get_carriers(self%ebands, self%ntemp, self%kTmesh, self%transport_mu_e, self%n_ehst)


 time_opt =0
 cnt = 0
 do my_spin=1,self%my_nspins
   spin = self%my_spins(my_spin)
   do my_ik=1,self%my_nkcalc
     ikcalc= self%my_ikcalc(my_ik)
     !if (ikcalc > 1) cycle
     ik_ibz = self%kcalc2ibz(ikcalc, 1) 
     isym_k = self%kcalc2ibz(ikcalc, 2)
     trev_k = self%kcalc2ibz(ikcalc, 6)

     wtk = self%ebands%wtk(ik_ibz)
     S = transpose(cryst%symrel_cart(:,:,isym_k))

     nbands = self%nbcalc_ks(ikcalc, spin)
     call cwtime(cpu, wall, gflops, "start")

     do ib=self%bmin,self%bmax
       !if (ib > self%bmin) cycle ! MG HACK To be able to run tests quickly.
       ib_eph = ib - self%bmin + 1
       eig_nk = self%ebands%eig(ib, ik_ibz, spin)
       wr_step = self%wrmesh_b(2,ib_eph,my_ik,spin) - self%wrmesh_b(1,ib_eph,my_ik,spin)
       vr(:) = self%vbks(:, ib, ik_ibz, spin)
       ! Store outer product (v_bks x v_bks) in vv_tens. This part does not depend on T and irta.
       do ii=1,3
         do jj=1,3
           vv_tens(ii, jj) = vr(ii) * vr(jj)
         end do
       end do
       ! Calculation of the velocity tensor
       vv_tens = cryst%symmetrize_cart_tens33(vv_tens, time_opt)
         do itemp=1,self%ntemp
	 !if (itemp > 1) cycle
         
         Tkelv = self%kTmesh(itemp) / kb_HaK; if (Tkelv < one) Tkelv = one
           do iw=1, self%nwr
!             if (mod(iw, self%wt_comm%nproc) /= self%wt_comm%me) cycle  ! MPI parallelism over freqs

                !  Preparing all elements needed for conductivity
                omega = self%wrmesh_b(iw,ib_eph,my_ik,my_spin)
                sp_func = -aimag (self%gw_vals(iw, itemp, ib_eph, my_ik, my_spin) ) / pi
!                test_Aw(itemp) = test_Aw(itemp) + sp_func
                dfdw = occ_dfde(omega, self%kTmesh(itemp), self%mu_e(itemp))
                self%print_dfdw(iw,itemp) = dfdw
!                test_dfdw(itemp) = test_dfdw(itemp) + dfdw
                kernel(iw) = - dfdw * sp_func**2 
                Aw(iw) = sp_func**2
                dfdw_acc(iw) = dfdw

           end do !iw
         mu_e = self%transport_mu_e(itemp)
         ieh = 2; if (eig_nk >= mu_e) ieh = 1
         integration = simpson( wr_step, kernel)
         int_Aw = simpson(wr_step,Aw)
         int_dfdw = simpson(wr_step,dfdw_acc)
         ! Calculation of the conductivity
         self%l0( :, :, ieh, spin, itemp ) = self%l0( :, :, ieh, spin, itemp ) + integration*vv_tens(:,:)*wtk
         Aw_l0(itemp) = Aw_l0(itemp) + int_Aw*wtk*vv_tens(1,1) 
         dfdw_l0(itemp) = dfdw_l0(itemp) + int_dfdw*wtk*vv_tens(1,1)
         call inv33(self%l0(:, :, ieh, spin, itemp), work_33)

         l0inv_33nw(:,:,ieh) = work_33
         self%l1( :, :, ieh, spin, itemp ) = self%l0( :, :, ieh, spin, itemp )*(eig_nk - self%mu_e(itemp)) 
         self%l2( :, :, ieh, spin, itemp ) = self%l1( :, :, ieh, spin, itemp )*(eig_nk - self%mu_e(itemp)) 

         self%seebeck(:,:,ieh,spin,itemp) = matmul(work_33, self%l1(:,:,ieh,spin,itemp)) / Tkelv

         work_33 = self%l1(:, :, ieh, spin, itemp)
         work_33 = self%l2(:, :, ieh, spin, itemp) - matmul(work_33, matmul(l0inv_33nw(:, :, ieh), work_33))
         self%kappa(:,:,ieh,spin,itemp) = work_33 / Tkelv
         !self%conductivity_mu( :, :, ieh, spin, itemp ) = self%conductivity_mu( :, :, ieh, spin, itemp ) + integration*vv_tens(:,:)*wtk
!         call xmpi_sum(self%conductivity_mu(:, :, ieh, spin, itemp) , self%wt_comm%value, ierr)
         !print *,"cond",ikcalc, itemp, self%conductivity_mu(1, 1, 1,1, itemp)
         end do ! itemp

     end do !ib

   end do ! my_ik
   ! Collect data if k-points parallelism.
   !call xmpi_sum(self%conductivity_mu , self%kcalc_comm%value, ierr)
 end do !my_spin
 max_occ = two / (self%nspinor * self%nsppol)
 fact0 = max_occ * (siemens_SI / Bohr_meter / cryst%ucvol) / 100
 self%conductivity_mu = fact0 * self%l0  ! siemens cm^-1
 self%seebeck = - volt_SI  * max_occ * self%seebeck
 self%kappa = + volt_SI**2 * fact0 * self%kappa
 do itemp=1, self%ntemp
         Tkelv = self%kTmesh(itemp) / kb_HaK; if (Tkelv < one) Tkelv = one
        print *, Tkelv, fact0*Aw_l0(itemp), fact0*dfdw_l0(itemp)
 end do

 ! Scale by the carrier concentration
 fact = 100**3 / e_Cb
 do my_spin=1,self%my_nspins
   spin = self%my_spins(my_spin)
     do itemp=1,self%ntemp
       do ieh=1,2 ! e/h
         call safe_div(fact * self%conductivity_mu(:,:,ieh,spin,itemp), &
                       self%n_ehst(ieh, spin, itemp) / cryst%ucvol / Bohr_meter**3, zero, &
                       self%mobility_mu(:,:,ieh,spin,itemp))
       end do
     end do
   end do


 call cwtime_report(" cumulant_kubo_transport", cpu_kloop, wall_kloop, gflops_kloop)

 ABI_SFREE(kernel)
! ABI_SFREE(test_Aw)
! ABI_SFREE(test_dfdw)
 ABI_SFREE(Aw)
 ABI_SFREE(dfdw_acc)
 ABI_SFREE(Aw_l0)
 ABI_SFREE(dfdw_l0)

end subroutine cumulant_kubo_transport

!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_sigmaph_ncread
!! NAME
!! cumulant_sigmaph_ncread
!!
!! FUNCTION
!!   read out_SIGPH.nc file
!!
!! INPUTS
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! ncid=Netcdf file handle.
!!
!! PARENTS
!!      m_cumulant
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_sigmaph_ncread(self, path, ncid, comm)

!Arguments --------------------------------------
 class(cumulant_t),intent(inout) :: self
 character(len=fnlen),intent(in) :: path
 integer,intent(in) :: comm
 integer,intent(out) :: ncid

!Local variables --------------------------------
 integer :: ncerr,ierr
 real(dp) :: cpu, wall, gflops
 character(len=1000) :: msg
 !real(dp), allocatable :: vals_wr(:,:,:,:,:,:),vals_e0ks(:,:,:,:,:)

!************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 ! Open netcdf file
 msg = "Netcdf not activated at configure time!"
 ierr = 1
#ifdef HAVE_NETCDF
 ierr = 0

 if (.not. file_exists(path)) then
   msg = sjoin("Cannot find file", path)
   ierr = 1; return
 end if

 call cwtime(cpu, wall, gflops, "start")
 NCF_CHECK(nctk_open_read(ncid, path, comm))

 NCF_CHECK(nctk_get_dim(ncid, "nkcalc", self%nkcalc))
 NCF_CHECK(nctk_get_dim(ncid, "max_nbcalc", self%max_nbcalc))
 NCF_CHECK(nctk_get_dim(ncid, "nsppol", self%nsppol))
 NCF_CHECK(nctk_get_dim(ncid, "number_of_spinor_components", self%nspinor))
 NCF_CHECK(nctk_get_dim(ncid, "ntemp", self%ntemp))
 NCF_CHECK(nctk_get_dim(ncid, "nwr", self%nwr))
 NCF_CHECK(nctk_get_dim(ncid, "nqbz", self%nqbz))
 NCF_CHECK(nctk_get_dim(ncid, "nqibz", self%nqibz))
 NCF_CHECK(nctk_get_dim(ncid, "natom3", self%natom3))

 ABI_MALLOC(self%kcalc, (3, self%nkcalc))
 ABI_MALLOC(self%nbcalc_ks, (self%nkcalc, self%nsppol))
 ABI_MALLOC(self%bstart_ks, (self%nkcalc, self%nsppol))
 ABI_MALLOC(self%kcalc2ibz, (self%nkcalc, 6))
 ABI_MALLOC(self%kTmesh, (self%ntemp))
 !ABI_MALLOC(self%wrmesh_b, (self%nwr, self%max_nbcalc, self%nkcalc, self%nsppol))
 !ABI_MALLOC(self%vals_wr, ( self%nwr, self%ntemp, self%max_nbcalc, self%nkcalc, self%nsppol))
 !ABI_MALLOC(self%vals_e0ks, ( self%ntemp, self%max_nbcalc, self%nkcalc, self%nsppol))
 !ABI_MALLOC(self%e0vals, (self%max_nbcalc, self%nkcalc, self%nsppol))
 NCF_CHECK(nf90_get_var(ncid, vid("ngqpt"), self%ngqpt))
 NCF_CHECK(nf90_get_var(ncid, vid("kcalc"), self%kcalc))
 NCF_CHECK(nf90_get_var(ncid, vid("kTmesh"), self%kTmesh))
 NCF_CHECK(nf90_get_var(ncid, vid("nbcalc_ks"), self%nbcalc_ks))
 NCF_CHECK(nf90_get_var(ncid, vid("bstart_ks"), self%bstart_ks))
 NCF_CHECK(nf90_get_var(ncid, vid("kcalc2ibz"), self%kcalc2ibz))

 call cwtime_report(" sigmaph_ncread", cpu, wall, gflops)
#endif

contains
 integer function vid(vname)
   character(len=*),intent(in) :: vname
   vid = nctk_idname(ncid, vname)
end function vid

end subroutine cumulant_sigmaph_ncread
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_ncwrite
!! NAME
!! cumulant_ncwrite
!!
!! FUNCTION
!!
!! INPUTS
!! path=Filenae of output netcdf file.
!! cryst<crystal_t>=Crystalline structure
!! dtset<dataset_type>=All input variables for this dataset.
!! ncid=Netcdf file handle.
!!
!! PARENTS
!!      m_cumulant
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_ncwrite(self, path, cryst, ebands)

!Arguments --------------------------------------
 class(cumulant_t),intent(in) :: self
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 character(len=*),intent(in) :: path

!Local variables --------------------------------
 integer,parameter :: master = 0
 integer :: ncerr, ncid, my_rank, ikcalc, spin, comm, my_ik, ib, itemp
 real(dp) :: cpu, wall, gflops

!************************************************************************

 !comm = self%comm my_rank =
 call wrtout([std_out, ab_out], ch10//sjoin("- Writing cumulant results to:", path))
 call cwtime(cpu, wall, gflops, "start")

#ifdef HAVE_NETCDF
 ! Only one proc create the file, write structure and define basic dimensions.
 ! Then we reopen the file in MPI-IO mode.

 if (xmpi_comm_rank(self%comm) == master) then

   NCF_CHECK(nctk_open_create(ncid, path, xmpi_comm_self))

   ! Write to netcdf file
   NCF_CHECK(cryst%ncwrite(ncid))
   ! FIXME: Cannot write ebands because it crashes in
   !   k_dependent = "no"; if (any(ebands%nband(1) /= ebands%nband)) k_dependent = "yes"
   ! Should understand why!
   !NCF_CHECK(ebands_ncwrite(ebands, ncid))

   ! Add cumulant dimensions.
   ncerr = nctk_def_dims(ncid, [ &
     nctkdim_t("nkcalc", self%nkcalc), nctkdim_t("max_nbcalc", self%max_nbcalc), &
     nctkdim_t("nsppol", self%nsppol), nctkdim_t("ntemp", self%ntemp), &
     nctkdim_t("nqbz", self%nsppol), nctkdim_t("nqibz", self%ntemp), &
     nctkdim_t("nwr", self%nwr)], &
     defmode=.True.)
   NCF_CHECK(ncerr)
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "eph_task",  "nbsum", "bsum_start", "bsum_stop", "symdynmat", &
     "ph_intmeth", "eph_intmeth", "qint_method", "eph_transport", &
     "imag_only", "symv1scf", "dvdb_add_lr", "mrta", "ibte_prep"])
   NCF_CHECK(ncerr)
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "eta", "wr_step", "eph_fsewin", "eph_fsmear", "eph_extrael", "eph_fermie", "ph_wstep", "ph_smear", "eph_phwinfact"])
   NCF_CHECK(ncerr)


   ! Define arrays. Note nkcalc instead of my_nkcalc
   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("bstart_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("nbcalc_ks", "int", "nkcalc, nsppol"), &
     nctkarr_t("ngqpt", "dp", "three"), &
     nctkarr_t("kcalc", "dp", "three, nkcalc"), &
     nctkarr_t("kcalc2ibz", "dp", " nkcalc, six"), &
     nctkarr_t("kTmesh", "dp", "ntemp"), &
     !nctkarr_t("mu_e", "dp", "ntemp"), &
     nctkarr_t("wrmesh_b", "dp", "nwr, max_nbcalc, nkcalc, nsppol"), &
     !nctkarr_t("vals_wr", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("gw_vals", "dp", "two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ks_enes", "dp", "max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("dw_vals", "dp", "ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t('dfdw',"dp", "nwr, ntemp"), &
     nctkarr_t('conductivity_mu',"dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('mobility_mu', "dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('seebeck',"dp", "three, three, two, nsppol, ntemp"), &
     nctkarr_t('kappa',"dp", "three, three, two, nsppol, ntemp") &
   ])
   NCF_CHECK(ncerr)

   if (self%debug == 1) then
     ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("time_mesh", "dp", "nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("ct_vals", "dp","two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("c1", "dp","two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("c2", "dp","two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("c3", "dp","two, nwr, ntemp, max_nbcalc, nkcalc, nsppol"), &
     nctkarr_t("gt_vals", "dp","two, nwr, ntemp, max_nbcalc, nkcalc, nsppol") ] )
     NCF_CHECK(ncerr)
   endif

   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "eta"), self%ieta))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nbsum"), self%nbsum))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ngqpt"), self%ngqpt))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc"), self%kcalc))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kTmesh"), self%kTmesh))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "nbcalc_ks"), self%nbcalc_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "bstart_ks"), self%bstart_ks))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kcalc2ibz"), self%kcalc2ibz))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "dfdw"), self%print_dfdw))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "conductivity_mu"), self%conductivity_mu))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "mobility_mu"), self%mobility_mu))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "seebeck"), self%seebeck))
   NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "kappa"), self%kappa))
   ! FIXME This part is wrong since these arrays are MPI distributed
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "wrmesh_b"), self%wrmesh_b))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "vals_wr"), c2r(self%vals_wr)))
   !NCF_CHECK(nf90_put_var(ncid, nctk_idname(ncid, "ks_enes"), self%e0vals))

   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Barrier to avoid race conditions.
 !call wrtout(std_out, "before barrier")
 call xmpi_barrier(self%comm)

 ! Only the procs in the ncwrite_comm communicator write to disk.
 if (self%ncwrite_comm%value == xmpi_comm_null) goto 100

 ! open file for parallel-IO mode inside comm. All procs in ncwrite_comm enter this part.
 call wrtout(std_out, sjoin(" Performing parallel IO with:", itoa(self%ncwrite_comm%nproc), "procs"))
 NCF_CHECK(nctk_open_modify(ncid, path, self%ncwrite_comm%value))

 NCF_CHECK(nctk_set_datamode(ncid))

 ! Activate collectve IO.
 ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "gw_vals"), nf90_collective)
 NCF_CHECK(ncerr)

 ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "wrmesh_b"), nf90_collective)
 NCF_CHECK(ncerr)

 ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "ks_enes"), nf90_collective)
 NCF_CHECK(ncerr)

 ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "dw_vals"), nf90_collective)
 NCF_CHECK(ncerr)



 if (self%debug == 1) then
   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "time_mesh"), nf90_collective)
   NCF_CHECK(ncerr)

   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "ct_vals"), nf90_collective)
   NCF_CHECK(ncerr)

   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "c1"), nf90_collective)
   NCF_CHECK(ncerr)

   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "c2"), nf90_collective)
   NCF_CHECK(ncerr)

   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "c3"), nf90_collective)
   NCF_CHECK(ncerr)

   ncerr = nf90_var_par_access(ncid, nctk_idname(ncid, "gt_vals"), nf90_collective)
   NCF_CHECK(ncerr)
 endif
 spin = self%my_spins(1)
 ikcalc = self%my_ikcalc(1) ! index of the first kcalc treated by this rank.

 ! Start to write my **contiguous block** of kpoints from this **global** location
 ! Each MPI proc writes my_nkcalc entries.

 ncerr = nf90_put_var(ncid, nctk_idname(ncid, "gw_vals"), c2r(self%gw_vals), &
                      start=[1,1,1,1,ikcalc,spin], &
                      count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
 NCF_CHECK(ncerr)

 ncerr = nf90_put_var(ncid, nctk_idname(ncid, "wrmesh_b"), self%wrmesh_b, &
                      start=[1,1,ikcalc,spin], &
                      count=[self%nwr, self%max_nbcalc, self%my_nkcalc, self%nsppol])
 NCF_CHECK(ncerr)

 ncerr = nf90_put_var(ncid, nctk_idname(ncid, "ks_enes"), self%e0vals, &
                      start=[1,ikcalc,spin], &
                      count=[self%max_nbcalc, self%my_nkcalc, self%nsppol])
 NCF_CHECK(ncerr)

 ncerr = nf90_put_var(ncid, nctk_idname(ncid, "dw_vals"), self%dw_vals, &
                      start=[1,1,ikcalc,spin], &
                      count=[self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
 NCF_CHECK(ncerr)




 if (self%debug == 1) then

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "time_mesh"), self%time_mesh, &
                        start=[1,1,1,ikcalc,spin], &
                        count=[self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "ct_vals"), c2r(self%ct_vals), &
                        start=[1,1,1,1,ikcalc,spin], &
                        count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "c1"), c2r(self%c1), &
                        start=[1,1,1,1,ikcalc,spin], &
                        count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "c2"), c2r(self%c2), &
                        start=[1,1,1,1,ikcalc,spin], &
                        count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "c3"), c2r(self%c3), &
                        start=[1,1,1,1,ikcalc,spin], &
                        count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)

   ncerr = nf90_put_var(ncid, nctk_idname(ncid, "gt_vals"), c2r(self%gt_vals), &
                        start=[1,1,1,1,ikcalc,spin], &
                        count=[2, self%nwr, self%ntemp, self%max_nbcalc, self%my_nkcalc, self%nsppol])
   NCF_CHECK(ncerr)
 endif

 NCF_CHECK(nf90_close(ncid))

 ! Write to ab_out for automatic testing.
 if (xmpi_comm_rank(self%comm) == master .and. is_open(ab_out)) then
   write(ab_out, "(/,a)")" Print first 10 frequenciew in gw_vals array (re-im) for testing purposes:"
   write(ab_out, "(2(a, i0))")" spin: ", spin, ", ikcalc: ", ikcalc
   my_ik = 1
   do ib=1,self%nbcalc_ks(ikcalc, spin)
     do itemp=1,self%ntemp
       write(ab_out, "(2(a,i0))")" gw_vals for itemp:", itemp, "ib: ", ib
       write(ab_out, "(*(es13.5))")dble(self%gw_vals(1:min(10, self%nwr), itemp, ib, my_ik, spin))
       write(ab_out, "(*(es13.5))")aimag(self%gw_vals(1:min(10, self%nwr), itemp, ib, my_ik, spin))
     end do
  end do
 end if

 100 call cwtime_report(" cumulant_ncwrite", cpu, wall, gflops)
#endif

end subroutine cumulant_ncwrite
!!***

!----------------------------------------------------------------------

!!****f* m_cumulant/cumulant_free
!! NAME
!! cumulant_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! INPUTS
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

subroutine cumulant_free(self)

!Arguments --------------------------------------
 class(cumulant_t),intent(inout) :: self

!************************************************************************

 ABI_SFREE(self%nbcalc_ks)
 ABI_SFREE(self%bstart_ks)
 ABI_SFREE(self%kcalc2ibz)
 ABI_SFREE(self%coords_kws)
 ABI_SFREE(self%my_spins)
 ABI_SFREE(self%my_ikcalc)
 ABI_SFREE(self%e0vals)
 ABI_SFREE(self%dw_vals)
 ABI_SFREE(self%gw_vals)
 !ABI_SFREE(self%ce_spfunc_wr)
 ABI_SFREE(self%conductivity_mu)
 ABI_SFREE(self%mobility_mu)
 ABI_SFREE(self%print_dfdw)
 ABI_SFREE(self%seebeck)
  ABI_SFREE(self%kappa)
 ABI_SFREE(self%l0)
 ABI_SFREE(self%l1)
 ABI_SFREE(self%l2)
 ABI_SFREE(self%time_mesh)
 ABI_SFREE(self%ct_vals)
 ABI_SFREE(self%c1)
 ABI_SFREE(self%c2)
 ABI_SFREE(self%c3)
 ABI_SFREE(self%gt_vals)
 ABI_SFREE(self%wrmesh_b)
 ABI_SFREE(self%vals_e0ks)
 ABI_SFREE(self%vals_wr)
 ABI_SFREE(self%kcalc)
 ABI_SFREE(self%kTmesh)
 call destroy_mpi_enreg(self%ce_mpi_enreg)

 call self%kcalc_comm%free()
 call self%spin_comm%free()
! call self%wt_comm%free()
 call self%ncwrite_comm%free()

end subroutine cumulant_free
!!***


integer function gcd(m, n) result(answer)
 ! The greatest common divisor (GCD) of two nonzero integers
 ! i.e. the largest positive integer that divides each of the integers.
 ! gcd(a, 0) = gcd(0, a) = |a|
 ! gcd(0, 0) is commonly defined as 0.

 integer,intent(in)  :: m, n
 integer :: irest,ifirst

 ifirst = iabs(m)
 answer = iabs(n)
 if (answer ==  0) then
    answer = ifirst
 else
    do
       irest = mod(ifirst,answer)
       if (irest == 0)  exit
       ifirst = answer
       answer = irest
    end do
    answer= iabs(answer)
 end if

end function gcd

integer function lcm(a, b)
 integer,intent(in) :: a, b
 lcm = iabs(a * b) / gcd(a,b)
end function lcm
!
!integer function gcd(a,b)
! integer :: a,b,t
! do while (b/=0)
!     t = b
!     b = mod(a,b)
!     a = t
! end do
! gcd = abs(a)
!end function gcd

! Two Factor decomposotionns of positive integer
! The first entry in facts is the largest factor unless is_print is True.
! 1 facts (1, 1) is_prime: True
! 2 facts (2, 1) is_prime: False
! 12 facts (6, 2) is_prime: False
! 11 facts (1, 11) is_prime: True
! 17 facts (1, 17) is_prime: True
! 33 facts (11, 3) is_prime: False
! 35 facts (7, 5) is_prime: False

subroutine ifact2(nn, facts, is_prime)

 integer,intent(in) :: nn
 integer,intent(out) :: facts(2)
 logical,intent(out) :: is_prime

!Local variables ------------------------------
 integer :: start, ii

! *************************************************************************

 ABI_CHECK(nn > 0, sjoin("invalid nn:", itoa(nn)))

 start = nn / 2 + 1
 do ii=start, 1, -1
   if (mod(nn, ii) == 0) exit
 end do

 facts = [ii, nn / ii]
 is_prime = facts(1) == 1

end subroutine ifact2
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


end module m_cumulant
!!***

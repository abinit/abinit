!!****m* ABINIT/m_gstore_sigeph
!! NAME
!! m_gstore_sigeph
!!
!! FUNCTION
!!  Compute (diagonal) matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  in the KS basis using precomputed e-ph matrix elements.
!!  See also m_sigmaph, for a version in which the g-matrix elements are computed on-the-fly.
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gstore_sigeph

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use netcdf
 use m_nctk
 use m_dvdb,           only : dvdb_t
 use m_crystal,        only : crystal_t
 use m_hamiltonian,    only : gs_hamiltonian_type, rf_hamiltonian_type
 use m_dtset,          only : dataset_type
 use m_dtfil,          only : datafiles_type
 use m_wfd,            only : wfd_t
 use m_ephtk
 !use m_mkffnl
 use m_sigtk
 use m_cgtools

 !use m_time,           only : cwtime, cwtime_report, sec2str
 use m_io_tools,       only : iomode_from_fname !, file_exists, is_open, open_file, flush_unit
 use m_numeric_tools,  only : arth !, c2r, get_diag, linfit, iseven, simpson_cplx, simpson, print_arr, inrange
 use m_fstrings,       only : tolower, itoa, ftoa, sjoin, ktoa, ltoa, strcat, replace_ch0, yesno, string_in
 use m_cgtools,        only : cg_zgemm !, cg_zdotc
 use m_kg,             only : getph !, mkkpg, mkkin
 use defs_datatypes,   only : pseudopotential_type
 use defs_abitypes,    only : mpi_type
 use m_hdr,            only : hdr_type, fform_from_ext
 use m_geometry,       only : phdispl_cart2red_nmodes
 use m_ebands,         only : ebands_t, gaps_t
 use m_kpts,           only : kpts_timrev_from_kptopt, kpts_map !, kpts_sort, kpts_pack_in_stars, kptrlatt_from_ngkpt
 use m_ioarr,          only : read_rhor
 use m_fftcore,        only : ngfft_seq !, sphereboundary, get_kg, kgindex
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_ifc,            only : ifc_type
 use m_dfpt_cgwf,      only : stern_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawcprj,        only : pawcprj_type, pawcprj_free
 use m_pstat,          only : pstat_proc
 use m_occ,            only : occ_be, occ_fd
 use m_lgroup,         only : lgroup_t
 use m_gstore,         only : gstore_t, gqk_t

 implicit none

 private
 public :: gstore_sigeph

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha
!!***

!----------------------------------------------------------------------

!!****t* m_gstore_sigeph/sep_t
!! NAME
!! sep_t
!!
!! FUNCTION
!! Container for the (diagonal) matrix elements of the electron-phonon self-energy
!! in the KS representation i.e. Sigma_eph(omega, T, band, k, spin).
!! Provides methods to compute QP corrections, spectral functions, QP linewidths and
!! save the results to netcdf file.
!!
!! TODO
!!  Fix problem with spin parallelism and output of results.
!!
!! SOURCE

 type,public :: sep_t

  integer :: nwr = 0
   ! Number of frequency points along the real axis for Sigma(w) and spectral function A(w)
   ! Odd number so that the mesh is centered on the KS energy.
   ! The spectral function is computed only if nwr > 0 (taken from dtset%nfreqsp)

  integer :: ntemp = 0
  ! Number of temperatures.

  logical :: imag_only

  real(dp) :: wr_step
   ! Step of the linear mesh along the real axis (Ha units).

  complex(dp) :: ieta = zero
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  real(dp),allocatable :: kTmesh(:)
  ! kTmesh(ntemp)
  ! List of temperatures (kT units).

  real(dp),allocatable :: mu_e(:)
  ! mu_e(ntemp)
  ! chemical potential of electrons for the different temperatures.

  complex(dp),allocatable :: vals_e0ks(:,:,:)
  ! Sigma_eph(omega=eKS, kT, band) for given (ikcalc, spin).
  ! Fan-Migdal + Debye-Waller

  complex(dp),allocatable :: fan_vals(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! Fan-Migdal

  complex(dp),allocatable :: fan_stern_vals(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! Fan-Migdal adiabatic Sternheimer part

  complex(dp),allocatable :: dvals_de0ks(:,:,:)
  ! (ntemp, nb_k, nkcalc)
  ! d Re Sigma_eph(omega, kT, band, kcalc) / d omega (omega=eKS)

  real(dp),allocatable :: dw_vals(:,:,:)
  !  dw_vals(ntemp, nb_k, nkcalc) for given (ikcalc, spin)
  !  Debye-Waller term (static).

  real(dp),allocatable :: dw_stern_vals(:,:,:)
   !  dw_stern_vals(ntemp, nb_k, nkcalc)
   !  Debye-Waller Sternheimer term (static) .

  complex(dp),allocatable :: vals_wr(:,:,:,:)
   ! vals_wr(nwr, ntemp, nb_k, nkcalc)
   ! Sigma_eph(omega, kT, band)
   ! enk_KS corresponds to nwr/2 + 1.

  real(dp),allocatable :: wrmesh_b(:,:,:)
   ! wrmesh_b(nwr, nb_k, nkcalc)
   ! Frequency mesh along the real axis (Ha units) used for the different bands
   ! Each mesh is **centered** on the corresponding KS energy.

  !integer :: phmesh_size
   ! Number of phonon frequencies in phonon mesh used for Eliashberg functions and
   ! and other omega-resolved quantities.

  !real(dp),allocatable :: phmesh(:)
   ! phmesh(phmesh_size)
   ! phonon mesh in Ha.

 contains

   procedure :: gather_and_write_results => sep_gather_and_write_results
   ! Write main dimensions and header of sigmaph on a netcdf file.

   procedure :: free => sep_free
   ! Free dynamic memory
 end type sep_t
!!***

contains
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/gstore_sigeph
!! NAME
!!  gstore_sigeph
!!
!! FUNCTION
!!  Compute matrix elements of the e-ph self-energy (Fan Migdal + Debye Waller).
!!  using precomputed e-ph matrix elements.
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! dtfil<datafiles_type>=Variables related to files.
!! cryst: Crystalline structure
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!!
!! OUTPUT
!!
!! SOURCE

subroutine gstore_sigeph(wfk0_path, ngfft, ngfftf, dtset, dtfil, cryst, ebands, dvdb, ifc, &
                         pawfgr, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb
 type(ifc_type),target,intent(in) :: ifc
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: comm
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 !@type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
 integer,parameter :: master = 0, with_cplex1 = 1, cplex1 = 1, pawread0 = 0, ndat1 = 1, berryopt0 = 0, istwfk_1 = 1
 integer :: n1, n2, n3, n4, n5, n6, nb_k, nb_kq, glob_nk, ntemp, cplex
 integer :: spin, my_is, my_ik, my_iq, my_ip, in_k, im_kq, ierr, gap_err, my_rank, ip1, ip2, nu, ipc, idir, ipert
 integer :: it, ik_ibz, ikq_ibz, band_k, band_kq, timrev_k, ii, ikcalc, natom, natom3, nsppol, nspden, nspinor, nkpt !,ik_bz
 integer :: isym_k,isym_kq,trev_k,trev_kq
 integer :: istwf_k, istwf_kq, npw_k, npw_kq, nkpg_kq
 integer :: nfft, nfftf, mgfft, mgfftf, nkpg !,nkpg1,cnt, enough_stern
 integer :: usecprj, mpw, nbsum, ibsum_kq, band_me, u1_band ! usevnl,optlocal,optnl,opt_gvnlx1,
 real(dp) :: wqnu, gkq2, weight_q, eig0nk, eig0mk, eig0mkq, ediff, gmod2, hmod2, gdw2, rfact, gdw2_stern, rtmp !,nqnu,gkq2,gkq2_pf,
 !real(dp) :: cpu, wall, gflops
 logical :: q_is_gamma, intra_band, same_band, isirr_k, isirr_kq, stern_use_cache
 complex(dp) :: cfact !, sig_cplx
 character(len=5000) :: msg
 type(gaps_t) :: gaps
 type(lgroup_t) :: lg_myk
 type(gstore_t) :: gstore
 type(sep_t) :: sigma
 type(hdr_type) :: pot_hdr
 type(crystal_t) :: pot_cryst
 type(wfd_t) :: wfd
 !type(u1_cache_t) :: u1c
 type(stern_t) :: stern
 type(gs_hamiltonian_type) :: gs_ham_kq
 type(rf_hamiltonian_type) :: rf_ham_kq
!arrays
 integer :: gmax(3), g0_k(3), g0_kq(3), work_ngfft(18) !, g0_q(3),
 integer :: units(2), my_kqmap(6)
 integer,allocatable :: phmodes_skip(:)
 real(dp) :: kk(3), kk_ibz(3), kq_ibz(3), qpt(3), kq(3), fermie1_idir_ipert(3,cryst%natom), dotri(2)
 real(dp),allocatable :: vtrial(:,:), work(:,:,:,:)
 real(dp),allocatable :: kinpw_k(:), kinpw_kq(:),kpg_kq(:,:),kpg_k(:,:)
 real(dp),allocatable :: ffnl_k(:,:,:,:),ffnl_kq(:,:,:,:),ph3d_k(:,:,:),ph3d_kq(:,:,:),v1scf(:,:,:,:)
 real(dp) :: displ_red_nu(2, 3, cryst%natom)
 !real(dp),allocatable :: gkq_atm(:,:,:), gkq_nu(:,:,:) gaussw_qnu(:)
 real(dp),allocatable :: cg1s_kq(:,:,:,:), h1kets_kq_allperts(:,:,:,:)
 integer,allocatable :: nband(:,:), wfd_istwfk(:), kg_kq(:,:) !, , qselect(:), !, kg_k(:,:),
 integer,allocatable :: gbound_kq(:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:) !, ihave_ikibz_spin(:,:)
 real(dp) :: vec_natom3(2, 3*cryst%natom) ! zpr_frohl_sphcorr(3*cryst%natom),
 real(dp),allocatable :: bra_kq(:,:), kets_k(:,:,:)
 real(dp),allocatable :: stern_ppb(:,:,:,:), stern_fan_t(:), stern_dw(:,:,:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: rfact_t(:), nqnu_t(:), f_mkq(:) !, f_nk(:)
 complex(dp),allocatable :: cfact_t(:), cfact2_t(:), cfact_wr(:), tpp_red(:,:) !,fmw_frohl_sphcorr(:,:,:,:),
 type(pawrhoij_type),allocatable :: pot_pawrhoij(:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:), cwaveprj(:,:)
!----------------------------------------------------------------------

 units = [std_out, ab_out]
 my_rank = xmpi_comm_rank(comm)

 call wrtout(std_out, " Computing Fan-Migdal + DW self-energy from GSTORE.nc", pre_newlines=1)

 ! Init gstore and MPI grid from file and dtset.
 ! The Fan-Migdal SE requires |g(k,q)|^2 in the phonon representation but
 ! we also need complex g(k,q=Gamma) in the atom representation for the DW term in the RIA.
 call gstore%from_ncpath(dtfil%filgstorein, with_cplex1, dtset, cryst, ebands, ifc, &
                         "phonon", dtset%gstore_gname, .True., comm)

 natom = cryst%natom; natom3 = 3 * cryst%natom; nkpt = ebands%nkpt
 nsppol = dtset%nsppol; nspden = dtset%nspden; nspinor = dtset%nspinor

 ! Consistency check.
 ierr = 0
 if (gstore%qzone /= "bz") then
   ABI_ERROR_NOSTOP("gstore_sigeph assumes qzone == `bz`", ierr)
 end if
 if (gstore%has_used_lgk /= 0) then
   ABI_ERROR_NOSTOP("gstore_sigeph does not support use_lgk /=0.", ierr)
 end if
 if (gstore%has_used_lgq /= 0) then
   ABI_ERROR_NOSTOP("gstore_sigeph does not support use_lgq /=0.", ierr)
 end if
 if (ierr /= 0) then
   write(msg,'(a,i0,5a)')&
     'Checking consistency of input data against itself gave ',ierr,' inconsistencies.',ch10,&
     'The details of the problems can be FOUND ABOVE (or in output or log file), in an earlier WARNING.',ch10,&
     'In parallel, the details might not even be printed there. Then, try running in sequential to see the details.'
   ABI_ERROR(msg)
 end if

 ! TODO: For off-diagonal terms, see https://journals.aps.org/prb/abstract/10.1103/PhysRevB.101.121102

 ! Check consistency of little group options
 ABI_CHECK(gstore%check_little_group(dtset, msg) == 0, msg)

 ! FFT meshes from input file, not necessarily equal to the ones found in the external files.
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Initialize parameters in sigma object.
 sigma%imag_only = .False.
 sigma%ieta = + j_dpc * dtset%zcut

 ! Build (linear) mesh of K * temperatures. tsmesh(1:3) = [start, step, num]
 call dtset%get_ktmesh(sigma%ntemp, sigma%kTmesh)

 ! Compute the chemical potential at the different physical temperatures with Fermi-Dirac.
 ! TODO: One should check that nband is > nbocc to avoid inaccuracies in mu_e.
 ABI_MALLOC(sigma%mu_e, (sigma%ntemp))
 sigma%mu_e(:) = ebands%fermie
 if (dtset%eph_fermie == zero) then
   call ebands%get_muT_with_fd(sigma%ntemp, sigma%ktmesh, dtset%spinmagntarget, dtset%prtvol, sigma%mu_e, gstore%comm)
 end if

 ! Compute gaps.
 gaps = ebands%get_gaps(gap_err)
 if (gap_err /= 0) then
   ABI_ERROR("Cannot compute fundamental and direct gap (likely metal)")
 end if

 if (my_rank == master) then
   call gaps%print(units, kTmesh=sigma%ktmesh, mu_e=sigma%mu_e, &
                   header="Gaps, band edges and relative position wrt Fermi level")
 end if
 call gaps%free()

 ! Frequency mesh for sigma(w) and spectral functions.
 call dtset%get_wrmesh_for_sigeph(sigma%nwr, sigma%wr_step)

 ! Setup a mask to skip accumulating the contribution of certain phonon modes.
 call ephtk_set_phmodes_skip(dtset%natom, dtset%eph_phrange, phmodes_skip)

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 ! 2) Perform the setup needed for the non-local factors:
 !
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 usecprj = 0
 call gs_ham_kq%init(psps, pawtab, nspinor, nsppol, nspden, natom,&
  dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg,&
  comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab,&
  usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 if (dtset%eph_stern /= 0 .and. .not. sigma%imag_only) then
   ! Prepare call to Sternheimer solver.
   ! The static correction to FM_nk is:
   !    \sum_{qnu} (2n_qnu + 1) <H^1_{qnu} psi_nk| psi^1_{nk; qnu}>

   ! Allocate work space arrays.
   ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
   ! Note nvloc in vlocal (we will select one/four spin components afterwards)
   ABI_CALLOC(vtrial, (nfftf, nspden))
   ABI_CALLOC(vlocal, (n4, n5, n6, gs_ham_kq%nvloc))

   ! Read the GS potential (vtrial) from input POT file.
   ! In principle one may store vtrial in the DVDB but getpot_filepath is simpler to implement.
   call wrtout(units, sjoin(" Reading GS KS potential for Sternheimer from: ", dtfil%filpotin))
   call read_rhor(dtfil%filpotin, cplex1, dtset%nspden, nfftf, ngfftf, pawread0, mpi_enreg, vtrial, pot_hdr, pot_pawrhoij, comm, &
                  allow_interp=.True., want_varname="vtrial")
   pot_cryst = pot_hdr%get_crystal()
   if (gstore%cryst%compare(pot_cryst, header=" Comparing input crystal with POT crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and POT do not agree! Check messages above!")
   end if
   call pot_cryst%free(); call pot_hdr%free()

   ! Initialize the wave function descriptor.
   ! Only wavefunctions for the symmetrical image of the k/k+q wavevectors treated by this MPI rank are stored.
   ABI_MALLOC(nband, (nkpt, nsppol))
   ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
   ABI_MALLOC(keep_ur, (dtset%mband, nkpt ,nsppol))

   nband = dtset%mband; bks_mask = .False.; keep_ur = .False.

   ! Initialize bks_mask
   call gstore%fill_bks_mask(dtset%mband, nkpt, nsppol, bks_mask)

   !if (dtset%userie == 124) then
   !  ! Debugging section have all states on each MPI rank.
   !  bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
   !end if

   ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
   ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
   ! that will be used to symmetrize the wavefunctions in G-space.
   call gstore%get_mpw_gmax(dtset%ecut, mpw, gmax)

   ! Init work_ngfft
   gmax = gmax + 4 ! FIXME: this is to account for umklapp, should also consider Gamma-only and istwfk
   gmax = 2*gmax + 1

   call ngfft_seq(work_ngfft, gmax)
   !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
   ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

   ! Allocate PW-arrays. Note mpw in kg_kq
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
   ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
   ABI_MALLOC(wfd_istwfk, (nkpt))
   wfd_istwfk = 1

   call wfd%init(cryst, pawtab, psps, keep_ur, dtset%mband, nband, nkpt, nsppol, bks_mask,&
                 dtset%nspden, nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
                 dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

   call wfd%print([std_out], header="Wavefunctions for Sternheimer.")
   call pstat_proc%print(_PSTAT_ARGS_)

   ABI_FREE(nband)
   ABI_FREE(bks_mask)
   ABI_FREE(keep_ur)
   ABI_FREE(wfd_istwfk)

   ! Read wavefunctions.
   call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

   ! Open the DVDB file
   call dvdb%open_read(ngfftf, xmpi_comm_self)
   ABI_CHECK(dvdb%has_fields("pot1", msg), msg)

   !if (sigma%pert_comm%nproc > 1) then
   !  ! Activate parallelism over perturbations
   !  call dvdb%set_pert_distrib(sigma%my_npert, natom3, sigma%my_pinfo, sigma%pert_table, sigma%pert_comm%value)
   !end if

   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, gstore%qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, xmpi_comm_self)

   !enough_stern = 0
   ABI_MALLOC(tpp_red, (natom3, natom3))
   ABI_MALLOC(gbound_kq, (2*wfd%mgfft+8, 2))
 end if ! eph_stern /= 0

 ! Allocate work space arrays used inside the loops. Then we are ready to go!
 ntemp = sigma%ntemp
 ABI_MALLOC(nqnu_t, (ntemp))
 ABI_MALLOC(f_mkq, (ntemp))
 ABI_MALLOC(cfact_t, (ntemp))
 ABI_MALLOC(cfact2_t, (ntemp))
 ABI_MALLOC(rfact_t, (ntemp))
 ABI_MALLOC(stern_fan_t, (ntemp))

 call pstat_proc%print(_PSTAT_ARGS_)

 ! Loop over collinear spins.
 do my_is=1,gstore%my_nspins
   associate (gqk => gstore%gqk(my_is), cryst => gstore%cryst)
   spin = gstore%my_spins(my_is); nb_k = gqk%nb_k; nb_kq = gqk%nb_kq; glob_nk = gqk%glob_nk

   ABI_CHECK(allocated(gqk%my_g2), "my_g2 is not allocated")
   ABI_CHECK(allocated(gqk%my_wnuq), "my_wnuq is not allocated")

   ! Allocate arrays for self-energy matrix elements.
   ABI_CALLOC(sigma%vals_e0ks, (ntemp, nb_k, glob_nk))
   ABI_CALLOC(sigma%dvals_de0ks, (ntemp, nb_k, glob_nk))
   ABI_CALLOC(sigma%fan_vals, (ntemp, nb_k, glob_nk))
   ABI_CALLOC(sigma%fan_stern_vals, (ntemp, nb_k, glob_nk))
   ABI_CALLOC(sigma%dw_vals, (ntemp, nb_k, glob_nk))
   ABI_CALLOC(sigma%dw_stern_vals, (ntemp, nb_k, glob_nk))

   ! Prepare computation of Sigma_{nk}(w) and spectral function.
   if (sigma%nwr > 0) then
     ABI_CALLOC(sigma%vals_wr, (sigma%nwr, ntemp, nb_k, glob_nk))
     ABI_CALLOC(sigma%wrmesh_b, (sigma%nwr, nb_k, glob_nk))
   end if

   ABI_CALLOC(stern_dw, (2, natom3, natom3, nb_k))

   ! Loop over my k-points in |n,k>.
   do my_ik=1,gqk%my_nk
     kk = gqk%my_kpts(:, my_ik)

     ik_ibz = gqk%my_k2ibz(1, my_ik); isym_k = gqk%my_k2ibz(2, my_ik)
     trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5, my_ik)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     ABI_CHECK(isirr_k, "For the time being the k-point in Sigma_{nk} must be in the IBZ")
     kk_ibz = ebands%kptns(:,ik_ibz)

     ! Will store results in sigma% using glob_ik index.
     ikcalc = gqk%my_k2glob(my_ik)

     ! Compute the little group of the k-point so that we can sum g(k,q) only for q in the IBZ_k.
     if (dtset%gstore_use_lgk /= 0) then
       timrev_k = kpts_timrev_from_kptopt(ebands%kptopt)
       call lg_myk%init(cryst, kk, timrev_k, gstore%nqbz, gstore%qbz, gstore%nqibz, gstore%qibz, xmpi_comm_self)
     end if

     if (sigma%nwr > 0) then
       ! Build linear mesh **centered** around the KS energy.
       do in_k=1,nb_k
         band_k = in_k + gqk%bstart_k - 1
         eig0nk = ebands%eig(band_k, ik_ibz, spin) - sigma%wr_step * (sigma%nwr / 2)
         sigma%wrmesh_b(:,in_k, ikcalc) = arth(eig0nk, sigma%wr_step, sigma%nwr)
      end do
     end if

     if (dtset%eph_stern /= 0) then
       npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
       call gs_ham_kq%eph_setup_k("k", kk, istwfk_1, npw_k, wfd%kdata(ik_ibz)%kg_k, dtset, cryst, psps, & ! in
                                  nkpg, kpg_k, ffnl_k, kinpw_k, ph3d_k, xmpi_comm_self)                   ! out
     end if

     ! Sum over my q-points.
     do my_iq=1,gqk%my_nq
       call gqk%myqpt(my_iq, gstore, weight_q, qpt); q_is_gamma = sum(qpt**2) < tol14

       ! weight_q is computed here. It depends whether we are summing over the full BZ or IBZ_k.
       weight_q = one / gstore%nqbz
       if (dtset%gstore_use_lgk /= 0) then
         ii = lg_myk%findq_ibzk(qpt); if (ii == -1) cycle; weight_q = lg_myk%weights(ii)
       end if

       ! Find the image of k+q in the IBZ.
       kq = kk + qpt
       if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq, my_kqmap) /= 0) then
         ABI_ERROR(sjoin("Cannot map k+q to IBZ with k+q:", ktoa(kq)))
       end if
       ikq_ibz = my_kqmap(1)

       ikq_ibz = my_kqmap(1); isym_kq = my_kqmap(2)
       trev_kq = my_kqmap(6); g0_kq = my_kqmap(3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)

       if (dtset%eph_stern /= 0 .and. .not. sigma%imag_only) then
         ! Activate Sternheimer.
         ! NB: Assume adiabatic AHC expression to compute the contribution of states above gqk%nb

         ! Get istwf_kq, npw_kq, kg_kq for k+q.
         call wfd%get_gvec_gbound(cryst%gmet, dtset%ecut, kq, ikq_ibz, isirr_kq, dtset%nloalg, & ! in
                                  istwf_kq, npw_kq, kg_kq, nkpg_kq, kpg_kq, gbound_kq)           ! out

         call gs_ham_kq%eph_setup_k("kq", kq, istwfk_1, npw_kq, kg_kq, dtset, cryst, psps, &   ! in
                                    nkpg, kpg_kq, ffnl_kq, kinpw_kq, ph3d_kq, xmpi_comm_self)  ! out

         ! Fourier interpolation of the DFPT potentials.
         call dvdb%get_ftqbz(qpt, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)

         ! Build global array with GS wavefunctions cg_kq at k+q to prepare call to dfpt_cgwf.
         ! NB: bsum_range is not compatible with Sternheimer.
         ! There's a check at the level of the parser in chkinp.

         ! NOTE: in the present version, we need to gather all nbsum bands on each core before calling dfpt_cgwf.
         ! In principle one can call dfpt_cgwf in band-para mode but then
         ! we are obliged to call the sternheimer solver with one psi1 and all procs in bsum_comm
         ! just to to be able to apply the projector operator.
         ! The present version is not memory efficient and leads to a big load imbalance if
         ! bsum%comm%nproc > nband_calc_ks

         stern_use_cache = merge(.True., .False., dtset%eph_stern == 1)
         fermie1_idir_ipert = zero ! FIXME: This is needed for metals.
         nbsum = gqk%nb_kq
         call stern%init(dtset, npw_k, npw_kq, nspinor, nbsum, nbsum, fermie1_idir_ipert, &
                         stern_use_cache, work_ngfft, mpi_enreg, xmpi_comm_self)

         ABI_MALLOC(bra_kq, (2, npw_kq*nspinor))
         do ibsum_kq=1, nbsum
           ! Reconstruct u_kq(G) from the IBZ image.
           call wfd%rotate_cg(ibsum_kq, ndat1, spin, kq_ibz, npw_kq, kg_kq, istwf_kq, &
                              cryst, my_kqmap, gbound_kq, work_ngfft, work, bra_kq)
           stern%cgq(:, :, ibsum_kq) = bra_kq
         end do ! ibsum_kq
         ABI_FREE(bra_kq)

         ! Loop over all 3*natom perturbations (Each core prepares its own potentials)
         ! In the inner loop, we calculate H1 * psi_k, stored in h1kets_kq on the k+q sphere.
         ! Allocate vlocal1 with correct cplex. Note nvloc
         ABI_MALLOC_OR_DIE(vlocal1, (cplex*n4, n5, n6, gs_ham_kq%nvloc, gqk%my_npert), ierr)
         ABI_CALLOC(cg1s_kq, (2, npw_kq*nspinor, natom3, nb_k))

         ! h1kets_kq are MPI distributed inside pert_comm but we need off-diagonal pp' terms --> collect results.
         ABI_CALLOC(h1kets_kq_allperts, (2, npw_kq*nspinor, natom3, nb_k))
         ! Compute S_pp' = <D_{qp} vscf u_nk|u'_{nk+q p'}>
         ABI_CALLOC(stern_ppb, (2, natom3, natom3, nb_k))

         do my_ip=1, gqk%my_npert
           ipc = gqk%my_pertcases(my_ip); idir = mod(ipc-1, 3) + 1; ipert = (ipc - idir) / 3 + 1
           !print *, "idir, ipert:", idir, ipert

           ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
           ! Each CPU prepares its own potentials.
           call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_ham_kq%nvloc, &
             pawfgr, mpi_enreg, vtrial, v1scf(:,:,:,my_ip), vlocal, vlocal1(:,:,:,:,my_ip))

           ! Continue to initialize the Hamiltonian (call it here to support dfpt_cgwf Sternheimer).
           call gs_ham_kq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

           ! Prepare application of the NL part.
           call rf_ham_kq%init(cplex, gs_ham_kq, ipert, has_e1kbsc=.true.)
           call rf_ham_kq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

           ! Wait for gatherv operation
           !if (.not. stern%has_band_para .and. cgq_request /= xmpi_request_null) call xmpi_wait(cgq_request, ierr)
           ABI_MALLOC(kets_k, (2, npw_k*nspinor, nb_k))

           do in_k=1,nb_k
             band_k = in_k + gqk%bstart_k - 1
             stern%bands_treated_now(:) = 0; stern%bands_treated_now(band_k) = 1
             stern%rank_band = 0; u1_band = band_k; band_me = band_k

             ! Init entry in cg1s_kq, either from cache or with zeros.
             cg1s_kq(:,:,ipc,in_k) = zero
             call wfd%copy_cg(band_k, ik_ibz, spin, kets_k(1, 1, in_k))

             !print *, "Stern for band_k", band_k, " with nb_kq:", gqk%nb_kq
             call stern%solve(u1_band, band_me, idir, ipert, qpt, gs_ham_kq, rf_ham_kq, &
                              ebands%eig(:,ik_ibz,spin), ebands%eig(:,ikq_ibz,spin), &
                              kets_k(:,:,in_k), cwaveprj0, cg1s_kq(:,:,ipc,in_k), cwaveprj, msg, ierr)
             ABI_CHECK(ierr == 0, msg)

             ! Store H(1) applied to GS wavefunction Psi_nk(0)
             ! stern%gh1c_n stores <G|H1|C0 band,k>
             h1kets_kq_allperts(:,:,ipc,in_k) = stern%gh1c_n
           end do ! in_k

           ABI_FREE(kets_k)
           call rf_ham_kq%free()
         end do ! my_ip  (loop over my perturbations)

         ! Compute <D_p H psi_nk | D_p' psi_nk> and store it in stern_ppb
         do in_k=1,nb_k
           call cg_zgemm("C", "N", npw_kq*nspinor, natom3, natom3, &
             h1kets_kq_allperts(:,:,:,in_k), cg1s_kq(:,:,:,in_k), stern_ppb(:,:,:,in_k))

            ! Save data for Debye-Waller that is performed outside the q-loop.
            if (q_is_gamma) stern_dw(:,:,:,in_k) = stern_ppb(:,:,:,in_k)
         end do

         ABI_FREE(cg1s_kq)
         ABI_FREE(v1scf)
         ABI_FREE(vlocal1)
         ABI_FREE(h1kets_kq_allperts)
         call stern%free()
       end if ! eph_stern

       ! Sum over my phonon modes.
       do my_ip=1,gqk%my_npert
         nu = my_ip + gqk%my_pert_start - 1
         wqnu = gqk%my_wnuq(my_ip, my_iq)
         ! Ignore unstable modes or modes that should be skipped.
         if (ephtk_skip_phmode(nu, wqnu, phmodes_skip, dtset%eph_phrange_w)) cycle

         nqnu_t(:) = occ_be(wqnu, sigma%kTmesh, zero)

         ! Compute T_pp'(q,nu) matrix in reduced coordinates.
         if (dtset%eph_stern /= 0) then
           call phdispl_cart2red_nmodes(natom, 1, cryst%gprimd, gqk%my_displ_cart(:,:,:,my_ip,my_iq), displ_red_nu)
           call sigtk_dw_tpp_red(natom, displ_red_nu, tpp_red)
         end if

         ! Sum over bands in |m,k+q>
         do im_kq=1,gqk%nb_kq
           band_kq = im_kq + gqk%bstart_kq - 1
           eig0mkq = ebands%eig(band_kq, ikq_ibz, spin)

           ! Compute electronic occ for all Temps (note mu_e(it) Fermi level)
           do it=1,ntemp
             f_mkq(it) = occ_fd(eig0mkq, sigma%kTmesh(it), sigma%mu_e(it))
           end do

           ! Loop over the n index in |n,k>.
           do in_k=1,nb_k
             band_k = in_k + gqk%bstart_k - 1
             eig0nk = ebands%eig(band_k, ik_ibz, spin)
             ediff = eig0nk - eig0mk
             intra_band = q_is_gamma .and. ediff <= TOL_EDIFF
             same_band = band_k == band_kq

             ! The frequency dependent part evaluated at eig0nk for all T.
             if (dtset%eph_ahc_type == 1) then
               cfact_t(:) =  (nqnu_t + f_mkq      ) / (eig0nk - eig0mkq + wqnu + sigma%ieta) + &
                             (nqnu_t - f_mkq + one) / (eig0nk - eig0mkq - wqnu + sigma%ieta)
             else
               cfact_t(:) =  (two * nqnu_t + one) / (eig0nk - eig0mkq + sigma%ieta)
             end if

             gkq2 = weight_q * gqk%my_g2(my_ip, im_kq, my_iq, in_k, my_ik) !; print *, "gkq2", gkq2
             cfact_t = cfact_t * gkq2

             ! Compute contribution to Fan-Migdal for M > gqk%nb_kq
             !if (dtset%eph_stern /= 0) then
             if (dtset%eph_stern /= 0 .and. im_kq == 1 .and. gqk%qpt_kpt_comm%me == 0) then
               ! sum_{pp'} d_p* Stern_{pp'} d_p' with d = displ_red_nu and S = stern_ppb(:,:,:,in_k)
               vec_natom3 = zero
               call cg_zgemm("N", "N", natom3, natom3, 1, stern_ppb(:,:,:,in_k), displ_red_nu, vec_natom3)
               dotri = cg_zdotc(natom3, displ_red_nu, vec_natom3)
               !write(std_out, *)"dotri:", dotri
               rfact = dotri(1)
               rfact = rfact * weight_q / (two * wqnu)
               stern_fan_t = (two * nqnu_t(:) + one) * rfact
               sigma%fan_stern_vals(:, in_k, ikcalc) = sigma%fan_stern_vals(:, in_k, ikcalc) + stern_fan_t
               cfact_t = cfact_t + stern_fan_t
             end if

             sigma%vals_e0ks(:, in_k, ikcalc) = sigma%vals_e0ks(:, in_k, ikcalc) + cfact_t
             sigma%fan_vals(:, in_k, ikcalc) = sigma%fan_vals(:, in_k, ikcalc) + cfact_t

             ! Derivative of FM sigma at eig0nk for all T.
             ! Accumulate d(Re Sigma) / dw(w=eKS) for state in_k
             !cfact(x) =  (nqnu_t + f_mkq      ) / (x - eig0mkq + wqnu + sigma%ieta) + &
             !            (nqnu_t - f_mkq + one) / (x - eig0mkq - wqnu + sigma%ieta)
             gmod2 = (eig0nk - eig0mkq + wqnu) ** 2
             hmod2 = (eig0nk - eig0mkq - wqnu) ** 2
             rfact_t(:) = (nqnu_t + f_mkq      ) * (-gmod2 + aimag(sigma%ieta)**2) / (gmod2 + aimag(sigma%ieta)**2) ** 2 + &
                          (nqnu_t - f_mkq + one) * (-hmod2 + aimag(sigma%ieta)**2) / (hmod2 + aimag(sigma%ieta)**2) ** 2

             sigma%dvals_de0ks(:, in_k, ikcalc) = sigma%dvals_de0ks(:, in_k, ikcalc) + gkq2 * rfact_t

             ! Accumulate Sigma(w) for state in_k if spectral function is wanted.
             if (sigma%nwr > 0) then
               ! Zcut version
               do it=1,ntemp
                 cfact_wr(:) = (nqnu_t(it) + f_mkq(it)      ) / (sigma%wrmesh_b(:,in_k, ikcalc) - eig0mkq + wqnu + sigma%ieta) + &
                               (nqnu_t(it) - f_mkq(it) + one) / (sigma%wrmesh_b(:,in_k, ikcalc) - eig0mkq - wqnu + sigma%ieta)
                 cfact_wr(:) = gkq2 * cfact_wr(:)

                 !if (intra_band .and. sigma%frohl_model == 1)  then
                 !  ! Add Frohlich correction to Sigma_nk(w)
                 !  cfact_wr(:) = zero; if (same_band) cfact_wr(:) = fmw_frohl_sphcorr(:,nu,it,in_k)
                 !end if

                 sigma%vals_wr(:,it,in_k,ikcalc) = sigma%vals_wr(:,it,in_k,ikcalc) + cfact_wr(:)

                 ! Add static term from Sternheimer to Sigma(w) as well.
                 !if (dtset%eph_stern /= 0) then
                 !  sigma%vals_wr(:, it, in_k, ikcalc) = sigma%vals_wr(:, it, in_k, ikcalc) + rtmp
                 !end if
               end do
             end if ! nwr > 0

             gdw2 = gqk%my_gdw2(my_ip, im_kq, my_iq, in_k, my_ik)
             !print *, "gdw2", gdw2

             ! Accumulate DW for each T, add it to Sigma(e0) and Sigma(w) as well
             ! - (2 n_{q\nu} + 1) * gdw2 / (e_nk - e_mk)
             if (abs(ediff) > EPHTK_WTOL) then
               cfact_t(:) = - weight_q * gdw2 * (two * nqnu_t + one)  / (ediff + sigma%ieta)
             else
               cfact_t(:) = zero
             endif

             if (dtset%eph_stern /= 0 .and. im_kq == 1 .and. gqk%qpt_kpt_comm%me == 0) then

               ! Compute DW term for M > gqk%nb_kq
               cfact = zero
               do ip2=1,natom3
                 do ip1=1,natom3
                   cfact = cfact + tpp_red(ip1, ip2) * cmplx(stern_dw(1,ip1,ip2,in_k), stern_dw(2,ip1,ip2,in_k), kind=dp)
                 end do
               end do
               ! There's no 1/two here because I don't symmetrize the expression.
               ! TODO: Test symmetrization, real quantity? add support for the different Eliashberg functions with Stern
               gdw2_stern = real(cfact) / (four * wqnu)

               ! Add contribution due to the Sternheimer. ediff is absorbed in Sternheimer.
               cfact2_t = - weight_q * gdw2_stern * (two * nqnu_t(it) + one)
               cfact_t = cfact_t + cfact2_t
               sigma%dw_stern_vals(:, in_k, ikcalc) = sigma%dw_stern_vals(:, in_k, ikcalc) + real(cfact2_t)
             end if

             sigma%dw_vals(:, in_k, ikcalc) = sigma%dw_vals(:, in_k, ikcalc) + real(cfact_t)
             sigma%vals_e0ks(:, in_k, ikcalc) = sigma%vals_e0ks(:, in_k, ikcalc) + real(cfact_t)

             if (sigma%nwr > 0) then
               ! Add static DW term to Sigma(w).
               do it=1,ntemp
                 sigma%vals_wr(:, it, in_k, ikcalc) = sigma%vals_wr(:, it, in_k, ikcalc) + real(cfact_t(it))
               end do
             end if
           end do ! in_k

         end do ! im_kq
       end do ! my_ip

       ABI_SFREE(kpg_kq)
       ABI_SFREE(ffnl_kq)
       ABI_SFREE(kinpw_kq)
       ABI_SFREE(ph3d_kq)
       ABI_SFREE(stern_ppb)
     end do ! my_iq

     call lg_myk%free()

     ABI_SFREE(kpg_k)
     ABI_SFREE(ffnl_k)
     ABI_SFREE(kinpw_k)
     ABI_SFREE(ph3d_k)
   end do ! my_ik

   ABI_SFREE(stern_dw)

   call sigma%gather_and_write_results(gstore, gqk, dtset, ebands)
   end associate
 end do ! my_is

 ABI_FREE(nqnu_t)
 ABI_FREE(f_mkq)
 ABI_FREE(cfact_t)
 ABI_FREE(cfact2_t)
 ABI_FREE(rfact_t)
 ABI_FREE(stern_fan_t)
 ABI_FREE(phmodes_skip)
 ABI_FREE(ph1d)
 ABI_SFREE(vtrial)
 ABI_SFREE(cfact_wr)
 ABI_SFREE(vtrial)
 ABI_SFREE(vlocal)
 ABI_SFREE(kg_kq)
 ABI_SFREE(gbound_kq)
 ABI_SFREE(tpp_red)
 ABI_SFREE(work)


 call wfd%free(); call gstore%free(); call sigma%free(); call gs_ham_kq%free()

end subroutine gstore_sigeph
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/sep_gather_and_write_results
!! NAME
!!  sep_gather_and_write_results
!!
!! FUNCTION
!!  Collect results for a given spin, average results in the degenerate subspace.
!!  Finally, write results to ab_out and netcdf file.
!!
!! INPUTS
!!
!! OUTPUT
!!
!! SOURCE

subroutine sep_gather_and_write_results(sigma, gstore, gqk, dtset, ebands)

!Local variables-------------------------------
 class(sep_t),intent(inout) :: sigma
 type(gstore_t),intent(in) :: gstore
 type(gqk_t),intent(in) :: gqk
 type(ebands_t),intent(in) :: ebands
 type(dataset_type),intent(in) :: dtset
 integer,parameter :: max_ntemp = 50
 integer :: it, in_k, ikcalc, ik_bz, spin, ierr, bstart_k, bstop_k, cnt, ndeg
 integer :: band_k,ik_ibz,ib_val,ib_cond,jj,ideg,ii,iw,nstates !, nb_k
 !integer :: nq_ibzk_eff, nelem, imyq, iq_ibz_k, sr_ncid
 logical :: changed !, iwrite
 real(dp) :: ravg,kse,kse_prev,dw,fan0,ks_gap,kse_val,kse_cond,qpe_oms,qpe_oms_val,qpe_oms_cond
 real(dp) :: ravg2 ! invsig2fmts, tau
 complex(dp) :: sig0c,zc,qpe,qpe_prev,qpe_val,qpe_cond,cavg1,cavg2,cavg3,cavg4
 character(len=500) :: this_gtype ! msg
 !integer :: grp_ncid, ncerr
 type(degtab_t) :: degtab
!arrays
 real(dp) :: kcalc(3)
 !integer, allocatable :: recvcounts(:), displs(:), nq_rank(:), kq_symtab(:,:), my_kq_symtab(:,:)
 integer,allocatable :: degblock(:,:)
 real(dp) :: qp_gaps(sigma%ntemp),qpoms_gaps(sigma%ntemp)
 !real(dp),allocatable :: aw(:,:,:), a2few_avg(:,:), gather_srate(:,:,:,:), grp_srate(:,:,:,:)
 real(dp) :: ks_enes(gqk%nb_k), ze0_vals(sigma%ntemp, gqk%nb_k)
 !real(dp) :: gfw_avg(sigma%phmesh_size, 3)
 complex(dp) :: qpoms_enes(sigma%ntemp, gqk%nb_k),qp_enes(sigma%ntemp, gqk%nb_k)
!! *************************************************************************

 spin = gqk%spin

 ! Sum partial terms inside qgk%comm.
 call xmpi_sum(sigma%vals_e0ks, gqk%comm%value, ierr)
 call xmpi_sum(sigma%dvals_de0ks, gqk%comm%value, ierr)
 call xmpi_sum(sigma%fan_vals, gqk%comm%value, ierr)
 call xmpi_sum(sigma%fan_stern_vals, gqk%comm%value, ierr)
 call xmpi_sum(sigma%dw_vals, gqk%comm%value, ierr)
 call xmpi_sum(sigma%dw_stern_vals, gqk%comm%value, ierr)
 if (sigma%nwr > 0) call xmpi_sum(sigma%vals_wr, gqk%comm%value, ierr)

 ! Only procs inside ncwrite_comm perform IO (ab_out and ncid)
 !iwrite = self%ncwrite_comm%value /= xmpi_comm_null
 !if (.not. iwrite) return

 this_gtype = "KS"
 if (gstore%gtype == "gwpt" .and. dtset%gstore_gname == "gvals") this_gtype = "GWPT"
 if (gstore%gtype == "gwpt" .and. dtset%gstore_gname == "gvals_ks") this_gtype = "KS"

 ! Write legend.
 if (spin == 1) then
   write(ab_out,"(a)")repeat("=", 80)
   write(ab_out,"(a)")" Final results in eV."
   write(ab_out,"(a)")" Notations:"
   write(ab_out,"(a)")"     eKS: Kohn-Sham energy. eQP: quasi-particle energy."
   write(ab_out,"(a)")"     eQP - eKS: Difference between the QP and the KS energy."
   write(ab_out,"(a)")"     SE1(eKS): Real part of the self-energy computed at the KS energy, SE2 for imaginary part."
   write(ab_out,"(a)")"     Z(eKS): Renormalization factor."
   write(ab_out,"(a)")"     FAN: Real part of the Fan term at eKS. DW: Debye-Waller term."
   write(ab_out,"(a)")"     DeKS: KS energy difference between this band and band-1, DeQP same meaning but for eQP."
   write(ab_out,"(a)")"     OTMS: On-the-mass-shell approximation with eQP ~= eKS + Sigma(omega=eKS)"
   write(ab_out,"(a)")"     TAU(eKS): Lifetime in femtoseconds computed at the KS energy."
   write(ab_out,"(a)")"     mu_e: Fermi level for given (T, nelect)"
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
   write(ab_out,"(2a)")" Using g(k,q) of type: ", trim(this_gtype)
   write(ab_out,"(a)")" "
   write(ab_out,"(a)")" "
 end if

 ! Compute QP energies and Gaps (Note that I'm assuming a non-magnetic semiconductor!)
 ib_val = nint(ebands%nelect / (two / ebands%nspinor)); ib_cond = ib_val + 1

 do ikcalc=1,gqk%glob_nk
   ik_bz = gstore%kglob2bz(ikcalc, spin)
   ik_ibz = gstore%kbz2ibz(1, ik_bz)
   kcalc = gstore%kbz(:, ik_bz)

   if (abs(dtset%symsigma) == 1) then
     ! Average self-energy matrix elements in the degenerate subspace.
     ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
     ! Here we make sure that all the degenerate states are included.
     ! Store also band indices of the degenerate sets, used to average final results.
     bstart_k = gqk%bstart_k; bstop_k = gqk%bstop_k
     call ebands%enclose_degbands(ik_ibz, spin, bstart_k, bstop_k, changed, TOL_EDIFF, degblock=degblock)
     !if (changed) then
     !  ABI_WARNING("Changed")
     !end if
     bstart_k = gqk%bstart_k; bstop_k = gqk%bstop_k

     ! Store band indices used for averaging (shifted by bstart_k)
     ndeg = size(degblock, dim=2)
     ABI_MALLOC(degtab%bids, (ndeg))

     do ii=1,ndeg
       ! Make sure boundaries are within the input nk states.
       ! In principle the nk states should be initialized so that all degenerate states are included.
       degblock(1, ii) = max(degblock(1, ii), bstart_k)
       degblock(2, ii) = min(degblock(2, ii), bstop_k)
       cnt = degblock(2, ii) - degblock(1, ii) + 1
       ABI_MALLOC(degtab%bids(ii)%vals, (cnt))
       degtab%bids(ii)%vals = [(jj, jj= &
         degblock(1, ii) - bstart_k + 1, &
         degblock(2, ii) - bstart_k + 1)]
     end do

     ! Average self-energy matrix elements in the degenerate subspace.
     do ideg=1,size(degtab%bids)
       associate (bids => degtab%bids(ideg)%vals)
       nstates = size(bids)

       do it=1,sigma%ntemp
         ! Average QP(T) and Z(T).
         cavg1 = sum(sigma%vals_e0ks(it, bids(:), ikcalc)) / nstates
         cavg2 = sum(sigma%dvals_de0ks(it, bids(:), ikcalc)) / nstates
         cavg3 = sum(sigma%fan_vals(it, bids(:), ikcalc)) / nstates
         cavg4 = sum(sigma%fan_stern_vals(it, bids(:), ikcalc)) / nstates
         ravg = sum(sigma%dw_vals(it, bids(:), ikcalc)) / nstates
         ravg2 = sum(sigma%dw_stern_vals(it, bids(:), ikcalc)) / nstates

         do ii=1,nstates
           sigma%vals_e0ks(it, bids(ii), ikcalc) = cavg1
           sigma%dvals_de0ks(it, bids(ii), ikcalc) = cavg2
           sigma%fan_vals(it, bids(ii), ikcalc) = cavg3
           sigma%fan_stern_vals(it, bids(ii), ikcalc) = cavg4
           sigma%dw_vals(it, bids(ii), ikcalc) = ravg
           sigma%dw_stern_vals(it, bids(ii), ikcalc) = ravg2
         end do ! ii

         if (sigma%nwr > 0) then
           ! Average Sigma(omega, T)
           do iw=1,sigma%nwr
             cavg1 = sum(sigma%vals_wr(iw, it, bids(:), ikcalc)) / nstates
             do ii=1,nstates
               sigma%vals_wr(iw, it, bids(ii), ikcalc) = cavg1
             end do
           end do
         end if

       end do ! it
       end associate
     end do ! ideg

     call degtab%free()
     ABI_FREE(degblock)
   end if ! symsigma == +1

   kse_val = huge(one) * tol6; kse_cond = huge(one) * tol6
   qp_enes = huge(one) * tol6; qpoms_enes = huge(one) * tol6
   ks_enes = huge(one) * tol6; ze0_vals = huge(one) * tol6
   ks_gap = -one; qpoms_gaps = -one; qp_gaps = -one

   ! Loop over temperatures.
   do it=1,sigma%ntemp
     ! Write header.
     if (it <= max_ntemp) then
       if (ebands%nsppol == 1) then
         write(ab_out,"(3a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", T: ", sigma%kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", sigma%mu_e(it) * Ha_eV
       else
         write(ab_out,"(3a,i1,a,f6.1,a,f8.3)") &
           "K-point: ", trim(ktoa(kcalc)), ", spin: ", spin, ", T: ",sigma%kTmesh(it) / kb_HaK, &
           " [K], mu_e: ", sigma%mu_e(it) * Ha_eV
       end if
       if (sigma%imag_only) then
         write(ab_out,"(a)")"   B    eKS    SE2(eKS)  TAU(eKS)  DeKS"
       else
         write(ab_out,"(a)")"   B    eKS     eQP    eQP-eKS   SE1(eKS)  SE2(eKS)  Z(eKS)  FAN(eKS)   DW      DeKS     DeQP"
       end if
     end if

     ! Loop over band n_k for this k-point and spin.
     do in_k=1,gqk%nb_k
       band_k = in_k + bstart_k - 1
       kse = ebands%eig(band_k, ik_ibz, spin)
       ks_enes(in_k) = kse
       sig0c = sigma%vals_e0ks(it, in_k, ikcalc)
       dw = sigma%dw_vals(it, in_k, ikcalc)
       fan0 = real(sig0c) - dw
       ! Compute QP energies with On-the-Mass-Shell approximation and first renormalization i.e. Z(eKS)
       ! TODO: Note that here I use the full Sigma including the imaginary part
       !zc = one / (one - sigma%dvals_de0ks(it, in_k))
       zc = one / (one - real(sigma%dvals_de0ks(it, in_k, ikcalc)))
       ze0_vals(it, in_k) = real(zc)
       qpe = kse + real(zc) * real(sig0c)
       qpe_oms = kse + real(sig0c)
       if (in_k == 1) then
         kse_prev = kse; qpe_prev = qpe
       end if
       if (band_k == ib_val) then
         kse_val = kse; qpe_val = qpe; qpe_oms_val = qpe_oms
       end if
       if (band_k == ib_cond) then
         kse_cond = kse; qpe_cond = qpe; qpe_oms_cond = qpe_oms
       end if

       if (it <= max_ntemp) then
         if (sigma%imag_only) then
           ! 1/tau  = 2 Imag(Sigma)
           !invsig2fmts = Time_Sec * 1e+15 / two
           !tau = 999999.0_dp
           !if (abs(aimag(sig0c)) > tol16) tau = invsig2fmts / abs(aimag(sig0c))
           !tau = min(tau, 999999.0_dp)
           !write(ab_out, "(i4,2(f8.3,1x),f8.1,1x,f8.3)") &
           !    band_k, kse * Ha_eV, aimag(sig0c) * Ha_eV, tau, (kse - kse_prev) * Ha_eV
         else
           write(ab_out, "(i4, 10(f8.3,1x))") &
             band_k, kse * Ha_eV, real(qpe) * Ha_eV, (real(qpe) - kse) * Ha_eV, &
             real(sig0c) * Ha_eV, aimag(sig0c) * Ha_eV, real(zc), &
             fan0 * Ha_eV, dw * Ha_eV, (kse - kse_prev) * Ha_eV, real(qpe - qpe_prev) * Ha_eV
         end if
       end if

       if (in_k > 1) then
         kse_prev = kse; qpe_prev = qpe
       end if
       qpoms_enes(it, in_k) = qpe_oms
       qp_enes(it, in_k) = qpe
       if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
         ! We have enough states to compute the gap.
         if (it == 1) ks_gap = kse_cond - kse_val
         qpoms_gaps(it) = qpe_oms_cond - qpe_oms_val
         qp_gaps(it) = real(qpe_cond - qpe_val)
       end if
     end do ! in_k

     ! Print KS and QP gaps.
     if (it <= max_ntemp) then
       if (.not. sigma%imag_only) then
         if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
           write(ab_out, "(a)")" "
           write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, &
             "(assuming bval:", ib_val, " ==> bcond:", ib_cond, ")"
           write(ab_out, "(2(a,f8.3),a)")" QP gap: ",qp_gaps(it) * Ha_eV," (OTMS: ",qpoms_gaps(it) * Ha_eV, ")"
           write(ab_out, "(2(a,f8.3),a)")" QP_gap - KS_gap: ",(qp_gaps(it) - ks_gap) * Ha_eV,&
               " (OTMS: ",(qpoms_gaps(it) - ks_gap) * Ha_eV, ")"
           write(ab_out, "(a)")" "
         end if
       else
         if (kse_val /= huge(one) * tol6 .and. kse_cond /= huge(one) * tol6) then
           write(ab_out, "(a)")" "
           write(ab_out, "(a,f8.3,1x,2(a,i0),a)")" KS gap: ",ks_gap * Ha_eV, "(assuming bval:",ib_val," ==> bcond:",ib_cond,")"
           write(ab_out, "(a)")" "
         end if
       end if
       write(ab_out, "(a)")repeat("=", 92)
     end if

   end do ! it
 end do ! ikcalc

 if (sigma%ntemp > max_ntemp) then
   write(ab_out, "(a,i0,a)")" No more than ", max_ntemp, " temperatures are written to the main output file."
   write(ab_out, "(2a)")" Please use SIGEPH.nc file and AbiPy to analyze the results.",ch10
 end if

end subroutine sep_gather_and_write_results
!!***

!----------------------------------------------------------------------

!!****f* m_gstore_sigeph/sep_free
!! NAME
!!  sep_free
!!
!! FUNCTION
!!  Free dynamic memory.
!!
!! SOURCE

subroutine sep_free(sigma)

!Arguments ------------------------------------
 class(sep_t),intent(inout) :: sigma
! *********************************************************************

 ABI_SFREE(sigma%kTmesh)
 ABI_SFREE(sigma%mu_e)
 ABI_SFREE(sigma%vals_e0ks)
 ABI_SFREE(sigma%dvals_de0ks)
 ABI_SFREE(sigma%fan_vals)
 ABI_SFREE(sigma%fan_stern_vals)
 ABI_SFREE(sigma%dw_vals)
 ABI_SFREE(sigma%dw_stern_vals)
 ABI_SFREE(sigma%vals_wr)
 ABI_SFREE(sigma%wrmesh_b)

end subroutine sep_free
!!***

end module m_gstore_sigeph
!!***

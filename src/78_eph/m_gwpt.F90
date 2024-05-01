!!****m* ABINIT/m_gwpt
!! NAME
!!  m_gwpt
!!
!! FUNCTION
!!  Compute the electron-phonon matrix elements with the GWPT formalism
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2024 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwpt

 use, intrinsic :: iso_c_binding
#ifdef HAVE_MPI2
 use mpi
#endif
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_hide_blas
 use m_copy
 use m_ifc
 use m_ebands
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_skw
 use m_krank
 use m_lgroup
 use m_ephwg
 use m_sort
 use m_hdr
 use m_sigtk
 use m_ephtk
 use netcdf
 use m_nctk
 use m_rf2
 use m_dtset
 use m_dtfil
 use m_clib
 use m_mkffnl
 use m_screen

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_gwdefs,         only : GW_Q0_DEFAULT
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, print_arr, inrange
 use m_io_tools,       only : iomode_from_fname, file_exists, is_open, open_file, flush_unit
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex
 use m_cgtk,           only : cgtk_rotate, cgtk_change_gsphere
 use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map
 use m_kg,             only : getph, mkkpg
 use m_bz_mesh,        only : isamek, kmesh_t, find_qmesh
 use m_gsphere,        only : gsphere_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_ioarr,          only : read_rhor
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
 use m_dfpt_cgwf,      only : stern_t
 use m_phonons,        only : phstore_t, phstore_new
 use m_io_screening,   only : hscr_t, get_hscr_qmesh_gsph
 use m_vcoul,          only : vcoul_t
 use m_gstore,         only : gstore_t, gqk_t

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

!----------------------------------------------------------------------

!!****t* m_gwpt/gwpt_t
!! NAME
!! gwpt_t
!!
!! FUNCTION
!! Container for the (diagonal) matrix elements of the electron-phonon self-energy
!! in the KS representation i.e. gwpt_eph(omega, T, band, k, spin).
!! Provides methods to compute QP corrections, spectral functions, QP linewidths and
!! save the results to netcdf file.
!!
!! SOURCE

 type,public :: gwpt_t

  integer :: nkcalc
   ! Number of k-points computed (inside energy window)

  integer :: max_nbcalc
   ! Maximum number of bands computed (max over nkcalc and spin).

  integer :: nsppol
   ! Number of independent spin polarizations.

  integer :: nspinor
   ! Number of spinor components.

  integer :: symsigma
   ! 1 if matrix elements should be symmetrized.
   ! Required when the sum over q in the BZ is replaced by IBZ(k).

  integer :: timrev
   ! timrev = 1 if the use of time-reversal is allowed; 0 otherwise

  integer :: nbsum
   ! Total number of bands used in sum over states without taking into account MPI distribution.

  integer :: bsum_start, bsum_stop
   ! First and last band included in self-energy sum without taking into account MPI distribution inside bsum_comm
   ! nbsum = bsum_stop - bsum_start + 1

  integer :: my_bsum_start, my_bsum_stop
   ! Initial and final band index included in self-energy sum
   ! Processor-dependent if Re-Im calculation.
   ! Processor-independent and computed at runtime on the basis of the nk states in gwpt_{nk} if imag_only

  integer :: my_npert
   ! Number of atomic perturbations or phonon modes treated by this MPI rank.
   ! Note that natom3 are equally distributed. This allows us to use allgather instead of allgatherv

  type(xcomm_t) :: pert_comm
   ! MPI communicator for parallelism over atomic perturbations.

  !type(xcomm_t) :: qb_comm
   ! MPI communicator used to distribute (band_sum, q-points)

  type(xcomm_t) :: qpt_comm
   ! MPI communicator for q-points

  type(xcomm_t) :: bsum_comm
   ! MPI communicator for bands in self-energy sum

  type(xcomm_t) :: kcalc_comm
   ! MPI communicator for parallelism over k-points (high-level)

  type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins (high-level)

  !type(xcomm_t) :: pqb_comm
    ! MPI communicator for the (perturbation, band_sum, qpoint_sum)

  type(xcomm_t) :: ncwrite_comm
   ! MPI communicator for parallel netcdf IO used to write results for the different k-points/spins

  integer :: coords_pqbks(5)
   ! Cartesian coordinates of this processor in the Cartesian grid.

  integer :: nqbz
   ! Number of q-points in the (dense) BZ for gwpt integration

  integer :: nqibz
   ! Number of q-points in the (dense) IBZ for gwpt integration

  integer :: ncid = nctk_noid
   ! Netcdf file handle used to save results.

  !integer :: mpw
   ! Maximum number of PWs for all possible k+q

  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  integer :: ngqpt(3)
   ! Number of divisions in the Q mesh in the BZ.

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

  integer :: my_nspins
   ! Number of spins treated by this MPI rank

  integer,allocatable :: my_spins(:)
   ! my_spins(my_nspins)
   ! Indirect table giving the spin indices treated by this MPI rank.
   ! Used only in the collinear case with nsppol = 2 and nspinor == 1

  integer :: my_nkcalc
   ! Number of k-points treated by this MPI rank

  integer,allocatable :: my_ikcalc(:)
   ! my_ikcalc(my_nkcalc)
   ! List of ikcalc indices treated by this pool if k-point parallelism is activated.

  integer(i1b),allocatable :: itreat_qibz(:)
   ! itreat_qibz(nqibz)
   ! Table used to distribute potentials over q-points in the IBZ.
   ! The loop over qpts in the IBZ(k) is MPI distributed inside qpt_comm according to this table.
   ! 0 if this IBZ point is not treated by this proc.
   ! 1 if this IBZ is treated.

  integer,allocatable :: my_pinfo(:,:)
   ! my_pinfo(3, my_npert)
   ! my_pinfo(1, ip) gives the `idir` index of the ip-th perturbation.
   ! my_pinfo(2, ip) gives the `ipert` index of the ip-th perturbation.
   ! my_pinfo(3, ip) gives `pertcase`=idir + (ipert-1)*3

  integer,allocatable :: pert_table(:,:)
   ! pert_table(2, natom3)
   ! pert_table(1, npert): rank of the processor treating this atomic perturbation.
   ! pert_table(2, npert): imyp index in my_pinfo table, -1 if this rank is not treating ipert.

  integer,allocatable:: ind_qbz2ibz(:,:)
   ! (6, %nqibz)
   ! Mapping qBZ to IBZ

  integer,allocatable :: nbsum_rank(:,:)
   ! (%bsum_comm%nproc, 2)
   ! (rank+1, 1): Number of bands treated by rank in %bsum_comm.
   ! (rank+1, 2): bsum_start of MPI rank
   ! Available only if .not. imag_only

  real(dp),allocatable :: kcalc(:,:)
   ! kcalc(3, nkcalc)
   ! List of k-points where the self-energy is computed.

  real(dp),allocatable :: qbz(:,:)
   ! qbz(3, nqbz)
   ! Reduced coordinates of the q-points in the full BZ.

  real(dp),allocatable :: qibz(:,:)
   ! qibz(3, nqibz)
   ! Reduced coordinates of the q-points in the IBZ (full simmetry of the system).

  real(dp),allocatable :: wtq(:)
   ! wtq(nqibz)
   ! Weights of the q-points in the IBZ (normalized to one).

  !real(dp),allocatable :: vcar_calc(:,:,:,:)
   ! (3, max_nbcalc, nkcalc, nsppol))
   ! Diagonal elements of velocity operator in cartesian coordinates for all states in gwpt_nk.

  type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

  contains

    procedure :: print => gwpt_print
     ! Print results to main output file.

    procedure :: free => gwpt_free
     ! Free gwpt object

 end type gwpt_t
!!***

 public :: gwpt_run        ! Main entry point to compute GWpt matrix elements

 real(dp),private,parameter :: TOL_EDIFF = 0.001_dp * eV_Ha

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_gwpt/gwpt_run
!! NAME
!!  gwpt_run
!!
!! FUNCTION
!!  Compute e-ph matrix elements with the GWPT formalism.
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! wfk_hdr=Header of the WFK file.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! NOTES
!!
!!  1) Conventions used for g-spheres and the periodic part of the KS states:
!!
!!   _kq   --> k + q
!!   _kmp  --> k - p
!!   _kqmp  --> k + q - p
!!
!!  2) The routines used to symmetrize wavefunctions and DFPT scattering potentials
!!     expect in input symmetry tables generated using different conventions.
!!     For the wavefunctions we use the symrel convention while for the scattering potentials we use the symrec convention.
!!     We encode this in the name of the variable using e.g. mapc_qq for the symrec convention (C) and mapl_k convention (L)
!!
!!  3) The DFPT routines operate on double-precision wavefunctions stored in arrays with real/imag part e.g. cg(1:2,npw_k)
!!     while the GW routines operate on complex arrays of kind=gwpc where gwpc is defined at configure-time.
!!     The default value of gwpc is single-precision
!!     We use the following conventions:
!!
!!       cg_kq, cr_kq
!!       cg1_kqmp, cr1_kqmp
!!       ug_gk, ur_kq
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwpt_run(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, ifc, wfk_hdr, &
                    pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb !, drhodb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(hdr_type),intent(in) :: wfk_hdr
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18),ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_getgh1c1 = 1, berryopt0 = 0, istw1 = 1, ider0 = 0, idir0 = 0, istwfk1 = 1
 integer,parameter :: useylmgr = 0, useylmgr1 =0, master = 0, ndat1 = 1
 integer,parameter :: igscq0 = 0, icgq0 = 0, usedcwavef0 = 0, nbdbuf0 = 0, quit0 = 0, cplex1 = 1, pawread0 = 0
 integer :: band_me, nband_me, stern_comm
 integer :: my_rank,nsppol,nkpt,iq_ibz,iq_bz,my_npert
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs
 integer :: ibsum_kq, ib_k, u1c_ib_k, band_ks, u1_band, ib_sum, ii !, jj, iw !ib_kq
 !integer :: u1_master, ip
 integer :: mcgq, mgscq !, ig, ispinor !, ifft !nband_kq,
 integer :: my_spin, spin, idir,ipert, npw_pp !,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: isym_q, trev_q
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_k_ibz, istwf_k_ibz
 integer :: ikq_ibz, isym_kq, trev_kq, npw_kq, istwf_kq,  npw_kq_ibz, istwf_kq_ibz
 integer :: ikmp_ibz, isym_kmp, trev_kmp, npw_kmp, istwf_kmp, npw_kmp_ibz, istwf_kmp_ibz
 integer :: ikqmp_ibz, isym_kqmp, trev_kqmp, npw_kqmp, istwf_kqmp, npw_kqmp_ibz, istwf_kqmp_ibz
 integer :: mpw,ierr,imyq,ignore_kq, ignore_ibsum_kq ! band,
 integer :: n1,n2,n3,n4,n5,n6,nspden,nu, mqmem, mm_kq, nn_k
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,cnt,imyp !, restart
 integer :: tot_nlines_done, nlines_done, nline_in, grad_berry_size_mpw1, enough_stern
 integer :: nbcalc_ks,nbsum,bsum_start, bsum_stop, bstart_ks,ikcalc,bstart,bstop, sendcount !iatom,
 integer :: ipp_bz, comm_rpt, nqlwl, ebands_timrev ! osc_npw,
 integer :: ffnl_k_request, ffnl_kq_request, nelem, cgq_request
 integer :: nkibz, nkbz, qptopt, qtimrev, my_iq, my_ik
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all,cpu_ks,wall_ks,gflops_ks
 real(dp) :: cpu_setk, wall_setk, gflops_setk, cpu_qloop, wall_qloop, gflops_qloop
 real(dp) :: ecut,eshift,weight_q,rfact,ediff,q0rad, out_resid, bz_vol
 logical :: isirr_k, isirr_kq, isirr_kmp, isirr_kqmp, isirr_q, gen_eigenpb, qq_is_gamma, pp_is_gamma
 logical :: stern_use_cache, intra_band, same_band, stern_has_band_para, use_ftinterp
 type(krank_t) :: qrank_ibz ! krank_ibz,
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_ham_kq
 type(rf_hamiltonian_type) :: rf_ham_kq
 type(gwpt_t) :: gwpt
 type(ddkop_t) :: ddkop
 type(rf2_t) :: rf2
 type(crystal_t) :: pot_cryst, den_cryst
 type(hdr_type) :: pot_hdr, den_hdr
 type(stern_t) :: stern_kmp !, stern_kqmp
 type(kmesh_t) :: pp_mesh, kmesh
 type(gsphere_t) :: gsph_c
 type(hscr_t) :: hscr
 type(vcoul_t) :: vcp
 type(screen_t),target :: W
 type(screen_info_t) :: W_info
 type(gstore_t),target :: gstore
 type(gqk_t),pointer :: gqk
 character(len=fnlen) :: w_fname, gstore_filepath
 character(len=5000) :: msg
!arrays
 integer :: g0_k(3), g0_kq(3), g0_kmp(3), g0_kqmp(3)
 integer :: mapl_k(6), mapl_kq(6), mapl_kqmp(6), mapl_kmp(6), mapc_qq(6), mapc_qq2dvdb(6)
 integer :: units(2), work_ngfft(18), gmax(3)
 integer,allocatable :: bands_treated_now(:)
 integer(i1b),allocatable :: itreatq_dvdb(:)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), kg_kmp(:,:), kg_kqmp(:,:)
 integer,allocatable :: gbound_k(:,:), gbound_kq(:,:), gbound_kmp(:,:), gbound_kqmp(:,:), gbound_pp(:,:)
 integer,allocatable :: rank_band(:), root_bcalc(:) ! osc_indpw(:), osc_gvecq(:,:),
 integer,allocatable :: nband(:,:), qselect(:), wfd_istwfk(:)
 integer,allocatable :: ibzspin_2ikcalc(:,:)
 integer, allocatable :: displs(:), recvcounts(:)
 integer,allocatable :: qibz2dvdb(:)
 integer,pointer :: kg_pp(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3), kqmp(3), kmp(3), pp(3), kmp_ibz(3), kqmp_ibz(3)
 real(dp) :: phfrq(3*cryst%natom), dotri(2), qq_ibz(3), qpt(3) !qpt_cart(3),
 real(dp) :: tsec(2) ! vk(3), vkq(3),
 real(dp) :: vec_natom3(2, 3*cryst%natom)
 real(dp) :: wqnu,gkq2,eig0nk,eig0mkq !eig0mk,
 real(dp),allocatable,target :: cgq(:,:,:)
 real(dp),allocatable :: qlwl(:,:),vcar_calc(:,:,:,:)
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:) ! grad_berry(:,:),
 real(dp),allocatable :: ffnl_k(:,:,:,:),ffnl_kq(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:),gkq0_atm(:,:,:,:)
 real(dp),allocatable :: gscq(:,:,:), out_eig1_k(:), cg1s_kq(:,:,:,:), h1kets_kq_allperts(:,:,:,:)
 real(dp),allocatable :: dcwavef(:, :), gh1c_n(:, :), ghc(:,:), gsc(:,:), stern_ppb(:,:,:,:), stern_dw(:,:,:,:)
 logical,allocatable :: ihave_ikibz_spin(:,:), bks_mask(:,:,:),keep_ur(:,:,:)
 real(dp),allocatable :: ug_k(:,:), ug_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cg_work(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: vtrial(:,:),gvnlx1(:,:),gvnlxc(:,:),work(:,:,:,:), rhor(:,:)
 real(dp),allocatable :: gs1c(:,:), gkq_allgather(:,:,:)
 !real(dp),allocatable :: wtk(:), kibz(:,:), kbz(:,:)
 !real(dp),allocatable :: phfreqs_qibz(:,:), pheigvec_qibz(:,:,:,:), eigvec_qpt(:,:,:)
 real(dp) :: ylmgr_dum(1,1,1)
 complex(gwpc),allocatable :: ur_k(:), ur_kq(:), ur_kmp(:), ur_kqmp(:), cwork_ur(:) !, workq_ug(:)
 complex(gwpc),allocatable :: crhog_bkmp_nk(:), crhog_mkq_bkqmp(:)
 !complex(dp),allocatable :: g_xc(:,:), g_ks(:,:), g_sigma(:,:)
 real(dp),allocatable :: cg_kmp(:,:), cg_kqmp(:,:), cg1_kqmp(:,:) !, cg1_kqmp(:,:)
 !real(dp),allocatable :: full_cg1_kmp(:,:), full_cg1_kqmp(:,:), full_cg1_kqmp(:,:) !, full_cg1_kqmp(:,:)
 !complex(gwpc),allocatable :: full_ur1_kmp(:), full_ur1_kqmp(:,:), full_ur1_kqmp(:,:) !, full_ur1_kqmp(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:), cwaveprj(:,:)
 type(pawrhoij_type),allocatable :: pot_pawrhoij(:), den_pawrhoij(:)
#if defined HAVE_MPI && !defined HAVE_MPI2_INPLACE
 integer :: me
 real(dp),allocatable :: cgq_buf(:)
 real(dp),pointer :: cgq_ptr(:)
#endif

!************************************************************************

 ! TODO LIST:
 ! 1) Restart capabilities:
 !      Check if a GSTORE.nc file is already available and read mask telling us the (qpt, spin, kpt)
 !      that have been computed. Use this mask to cycle loops.
 ! 2) Use filtering techniques in (qpt, kpt) and (band1, band2) space. Can use gstore input vars to implement that
 ! 3)

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 units = [std_out, ab_out]

 ! Check if a previous GSTORE file is present to restart the calculation
 ! Here we try to read an existing GSTORE file if eph_restart == 1.
 ! and we compare the variables with the state of the code
 !restart = 0; ierr = 1; gstore_filepath = strcat(dtfil%filnam_ds(4), "_GSTORE.nc")
 !if (my_rank == master .and. dtset%eph_restart == 1) then
 !  Read mask
 !end if

 gstore_filepath = strcat(dtfil%filnam_ds(4), "_GSTORE.nc")
 call gstore%init(gstore_filepath, dtset, wfk_hdr, cryst, ebands, ifc, comm)

 !call gstore%from_ncpath(gstore_filepath, with_cplex0, dtset, cryst, ebands, ifc, comm)

 !if (my_rank == master .and. dtset%eph_restart == 1) then
 !  if (ierr == 0) then
 !    if (any(gwpt_restart%qp_done /= 1)) then
 !      !call gwpt%compare(gwpt_restart)
 !      ! Get list of QP states that have been computed.
 !      gwpt%qp_done = gwpt_restart%qp_done
 !      restart = 1
 !      call wrtout(units, "- Restarting from previous GSTORE.nc file")
 !      call wrtout(units, sjoin("- Number of k-points completed:", itoa(count(gwpt%qp_done == 1)), "/", itoa(gwpt%nkcalc)))
 !    else
 !      restart = 0; gwpt%qp_done = 0
 !      msg = sjoin("Found GSTORE.nc file with all QP entries already computed.", ch10, &
 !                  "Will overwrite:", gstore_filepath, ch10, &
 !                  "Keeping backup copy in:", strcat(gstore_filepath, ".bkp"))
 !      call wrtout(ab_out, sjoin("WARNING: ", msg))
 !      ABI_WARNING(msg)
 !      ! Keep backup copy
 !      ABI_CHECK(clib_rename(gstore_filepath, strcat(gstore_filepath, ".bkp")) == 0, "Failed to rename SIGPEPH file.")
 !    end if
 !  end if
 !  call gwpt_restart%free()
 !end if

 !call xmpi_bcast(restart, master, comm, ierr)
 !call xmpi_bcast(gwpt%qp_done, master, comm, ierr)

 !if (restart == 0) then
 !  call gwpt%write(dtset, cryst, ebands, wfk_hdr, dtfil, comm)
 !else
 !  ! Open file inside ncwrite_comm to perform parallel IO if kpt parallelism.
 !  if (gwpt%ncwrite_comm%value /= xmpi_comm_null) then
 !    NCF_CHECK(nctk_open_modify(gwpt%ncid, gstore_filepath, gwpt%ncwrite_comm%value))
 !    NCF_CHECK(nctk_set_datamode(gwpt%ncid))
 !  end if
 !end if

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor
 nspden = dtset%nspden; nkpt = ebands%nkpt

 ! FFT meshes from input file, not necessarly equal to the ones found in the external files.
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 ecut = dtset%ecut ! dtset%dilatmx

 ! Construct object to store final results.
 gwpt = gwpt_new(dtset, ecut, cryst, ebands, ifc, dtfil, qrank_ibz, comm)
 call qrank_ibz%free()

 !if (my_rank == master) then
 !  call gwpt%print(dtset, ab_out)
 !  call gwpt%print(dtset, std_out)
 !end if

 ! Compute mpw and gmax
 ! TODO: here we have the additional pp to be considered !
 ! This is the maximum number of PWs for all possible k+q treated.
 call gstore%get_mpw_gmax(ecut, mpw, gmax)

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, shouls also consider Gamma-only and istwfk
 gmax = 2*gmax + 1
 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! Initialize the wave function descriptor.
 ! Each node has all k-points and spins and bands between my_bsum_start and my_bsum_stop
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkpt ,nsppol))

 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.

 ! Mapping gwpt_{k,s} states to IBZ. -1 if not computed
 !ABI_MALLOC(ibzspin_2ikcalc, (nkpt, nsppol))
 !ibzspin_2ikcalc = -1

 !! Each node needs the wavefunctions for gwpt_{nk}
 !! TODO: kcalc should depend on the spin!
 !do spin=1,gwpt%nsppol
 !  do ikcalc=1,gwpt%nkcalc
 !    ik_ibz = gwpt%kcalc2ibz(ikcalc, 1)
 !    bstart = gwpt%bstart_ks(ikcalc, spin)
 !    bstop = bstart + gwpt%nbcalc_ks(ikcalc, spin) - 1
 !    bks_mask(bstart:bstop, ik_ibz, spin) = .True.
 !    ibzspin_2ikcalc(ik_ibz, spin) = ikcalc
 !  end do
 !end do
 !ABI_FREE(ibzspin_2ikcalc)

 bks_mask(gwpt%my_bsum_start:gwpt%my_bsum_stop,:,:) = .True.

 !if (dtset%userie == 124) then
 !  ! Uncomment this line to have all states on each MPI rank.
 !  bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
 !end if

 ! This table is needed when computing the imaginary part:
 ! k+q states outside the energy window are not read hence their contribution won't be included.
 ! Error is small provided calculation is close to convergence.
 ! To reduce the error one should increase the value of phwinfact
 !ABI_MALLOC(ihave_ikibz_spin, (nkpt, nsppol))
 !ihave_ikibz_spin = .False.
 !do spin=1,gwpt%nsppol
 !  do ik_ibz=1,ebands%nkpt
 !    if (any(bks_mask(:, ik_ibz, spin))) ihave_ikibz_spin(ik_ibz, spin) = .True.
 !  end do
 !end do

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd_init(wfd, cryst, pawtab, psps, keep_ur, dtset%mband, nband, nkpt, nsppol, bks_mask,&
               nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print(header="Wavefunctions for GWPT calculation.", mode_paral='PERS')

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 !call wfd%change_ngfft(cryst, psps, ngfft)

 ! if PAW, one has to solve a generalized eigenproblem
 ! Be careful here because I will need sij_opt == -1
 usecprj = 0; gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1

 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))
 ABI_MALLOC(cwaveprj, (natom, nspinor*usecprj))
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red, (2, 3, cryst%natom, natom3))

 ABI_MALLOC(gbound_k, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(gbound_kq, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(gbound_kmp, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(gbound_kqmp, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(gbound_pp, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(cg_work, (2, mpw*wfd%nspinor))

 ! ============================
 ! Compute vnk matrix elements
 ! ============================
 ! Diagonal elements of velocity operator in cartesian coordinates for all states in gwpt_nk.

 ABI_CALLOC(vcar_calc, (3, gwpt%max_nbcalc, gwpt%nkcalc, nsppol))
 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

 ! Consider only the nk states in gwpt_nk
 ! All gwpt_nk states are available on each node so MPI parallelization is easy.
 call cwtime(cpu_ks, wall_ks, gflops_ks, "start", msg=" Computing v_nk matrix elements for all states in gwpt_nk...")
 cnt = 0
 do spin=1,nsppol
   do ikcalc=1,gwpt%nkcalc
     kk = gwpt%kcalc(:, ikcalc)
     bstart_ks = gwpt%bstart_ks(ikcalc, spin)
     ik_ibz = gwpt%kcalc2ibz(ikcalc, 1)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do ib_k=1,gwpt%nbcalc_ks(ikcalc, spin)
       cnt = cnt + 1; if (mod(cnt, nprocs) /= my_rank) cycle ! MPI parallelism.
       band_ks = ib_k + bstart_ks - 1
       call wfd%copy_cg(band_ks, ik_ibz, spin, cg_work)
       eig0nk = ebands%eig(band_ks, ik_ibz, spin)
       vcar_calc(:, ib_k, ikcalc, spin) = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cg_work, cwaveprj0)
     end do

   end do
 end do
 call xmpi_sum(vcar_calc, comm, ierr)
 ! Write v_nk to disk.
 !if (my_rank == master) then
 !  NCF_CHECK(nf90_put_var(gwpt%ncid, nctk_idname(gwpt%ncid, "vcar_calc"), vcar_calc))
 !end if
 ABI_FREE(vcar_calc)

 call ddkop%free()
 call cwtime_report(" Velocities", cpu_ks, wall_ks, gflops_ks)

 ! Precompute phonon frequencies and eigenvectors in the IBZ.
 ! These quantities are then used to symmetrize quantities for q in the IBZ(k) in order
 ! to reduce the number of calls to ifc%fourq (expensive if dipdip == 1)

 ! Radius of sphere with volume equivalent to the micro zone.
 q0rad = two_pi * (three / (four_pi * cryst%ucvol * gstore%nqbz)) ** third
 bz_vol = two_pi**3 / cryst%ucvol

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1   ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2      ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output

 !ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 ! 2) Perform the setup needed for the non-local factors:
 !
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_ham_kq, psps, pawtab, nspinor, nsppol, nspden, natom,&
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg,&
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab,&
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 ! Allocate work space arrays.
 ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, gs_ham_kq%nvloc))

 if (dtset%eph_stern /= 0) then
   ! Read the GS potential (vtrial) from input POT file
   ! In principle one may store vtrial in the DVDB but getpot_filepath is simpler to implement.
   call wrtout(units, sjoin(" Reading GS KS potential for Sternheimer from: ", dtfil%filpotin))
   call read_rhor(dtfil%filpotin, cplex1, nspden, nfftf, ngfftf, pawread0, mpi_enreg, vtrial, pot_hdr, pot_pawrhoij, comm, &
                  allow_interp=.True.)
   pot_cryst = pot_hdr%get_crystal()
   if (cryst%compare(pot_cryst, header=" Comparing input crystal with POT crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and POT do not agree! Check messages above!")
   end if
   call pot_cryst%free(); call pot_hdr%free()
 end if

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 !call drhodb%open_read(ngfftf, xmpi_comm_self)

 my_npert = gwpt%my_npert
 if (gwpt%pert_comm%nproc > 1) then
   !  Activate parallelism over perturbations
   call dvdb%set_pert_distrib(gwpt%my_npert, natom3, gwpt%my_pinfo, gwpt%pert_table, gwpt%pert_comm%value)
   !call drhodb%set_pert_distrib(gwpt%my_npert, natom3, gwpt%my_pinfo, gwpt%pert_table, gwpt%pert_comm%value)
 end if

 ! Find correspondence IBZ --> set of q-points in DVDB.
 ! use_ftinterp selects whether DFPT potentials should be read from the DVDB or Fourier-interpolated on the fly.
 ! Activate FT interpolation automatically if required q-points in the IBZ are not found in the DVDB.

 ! qibz2dvdb gives the mapping gstore%ibz --> dvdb%ibz
 ! TODO: Make sure we have the same ibz in rho1%ibz
 use_ftinterp = .False.
 ABI_MALLOC(qibz2dvdb, (gstore%nqibz))
 if (dvdb%find_qpts(gstore%nqibz, gstore%qibz, qibz2dvdb, comm) /= 0) then
   call wrtout(units, " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
   use_ftinterp = .True.
 else
   call wrtout(units, " DVDB file contains all q-points in the IBZ --> Reading DFPT potentials from file.")
   use_ftinterp = .False.
 end if

 if (use_ftinterp) then
   ! Use ddb_ngqpt q-mesh to compute the real-space represention of DFPT v1scf potentials to prepare Fourier interpolation.
   ! R-points are distributed inside comm_rpt
   ! Note that when R-points are distributed inside qpt_comm we cannot interpolate potentials on-the-fly
   ! inside the loop over q-points.
   ! In this case, indeed, the interpolation must be done in gwpt_setup_qloop once we know the q-points contributing
   ! to the integral and the potentials must be cached.
   !FIXME: qpt_comm is buggy.
   !if (gwpt%imag_only) comm_rpt = xmpi_comm_self
   !comm_rpt = gwpt%bsum_comm%value
   comm_rpt = xmpi_comm_self
   qptopt = ebands%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
   !call drhodb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)

   ! Build q-cache in the *dense* IBZ using the global mask qselect and itreat_qibz.
   ABI_MALLOC(qselect, (gstore%nqibz))
   qselect = 1
   call dvdb%ftqcache_build(nfftf, ngfftf, gstore%nqibz, gstore%qibz, dtset%dvdb_qcache_mb, qselect, gwpt%itreat_qibz, comm)

 else
   ABI_MALLOC(qselect, (dvdb%nqpt))
   qselect = 1
 end if

 call dvdb%print(prtvol=dtset%prtvol)
 !call drhodb%print(prtvol=dtset%prtvol)

 if (.not. use_ftinterp) then
   ! Need to translate itreat_qibz into itreatq_dvdb.
   ABI_ICALLOC(itreatq_dvdb, (dvdb%nqpt))
   do iq_ibz=1,gstore%nqibz
     if (gwpt%itreat_qibz(iq_ibz) == 0) cycle
     db_iqpt = qibz2dvdb(iq_ibz)
     ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(gstore%qibz(:, iq_ibz)), "in the DVDB file."))
     itreatq_dvdb(db_iqpt) = 1
   end do
   call dvdb%qcache_read(nfftf, ngfftf, dtset%dvdb_qcache_mb, qselect, itreatq_dvdb, comm)
   ABI_FREE(itreatq_dvdb)
 end if

 ABI_FREE(qibz2dvdb)
 ABI_FREE(qselect)

 ! ================
 ! HANDLE SCREENING
 ! ================
 ABI_CALLOC(rhor, (nfftf, nspden))

 nqlwl = 0; w_fname = ABI_NOFILE
 if (dtset%getscr /= 0 .or. dtset%irdscr /= 0 .or. dtset%getscr_filepath /= ABI_NOFILE) then
   w_fname = dtfil%fnameabi_scr
 end if

 call kmesh%init(cryst, wfk_hdr%nkpt, wfk_hdr%kptns, dtset%kptopt)

 !call wrtout(std_out, "Begin SCR part")
 if (w_fname /= ABI_NOFILE) then
   ! Read g-sphere and pp_mesh from SCR file.
   call get_hscr_qmesh_gsph(w_fname, dtset, cryst, hscr, pp_mesh, gsph_c, qlwl, comm)
   call hscr%free()
   nqlwl = size(qlwl, dim=2)
   w_info%use_mdf = MDL_NONE
 else
   ! Init pp_mesh from the K-mesh reported in the WFK file.
   call find_qmesh(pp_mesh, cryst, kmesh)
   ! The G-sphere for W and Sigma_c is initialized from ecuteps.
   call gsph_c%init(cryst, 0, ecut=dtset%ecuteps)
   dtset%npweps = gsph_c%ng

   ! We also need the density for the model dielectric function
   call read_rhor(dtfil%fildensin, cplex1, nspden, nfftf, ngfftf, pawread0, mpi_enreg, rhor, den_hdr, den_pawrhoij, comm, &
                  allow_interp=.True.)
   den_cryst = den_hdr%get_crystal()
   if (cryst%compare(den_cryst, header=" Comparing input crystal with POT crystal") /= 0) then
     ABI_ERROR("Crystal structure from WFK and DEN do not agree! Check messages above!")
   end if
   call den_cryst%free(); call pot_hdr%free()
   w_info%use_mdf = MDL_BECHSTEDT
   w_info%eps_inf = dtset%mdf_epsinf
   ABI_CHECK(w_info%eps_inf > zero, "Model dielectric function requires the specification of mdf_epsinf")
 end if

 if (nqlwl == 0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i0,a,3f9.6)')&
     "The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
     "Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)
   ABI_COMMENT(msg)
 end if

 ! TODO: In principle one should have gpsh_x as well.
 call vcp%init(gsph_c, cryst, pp_mesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_c%ng, &
               nqlwl, qlwl, comm)
 call kmesh%free()

 ! Init Wc.
 ! Incore or out-of-core solution?
 mqmem = 0; if (dtset%gwmem /10 == 1) mqmem = pp_mesh%nibz
 w_info%invalid_freq = dtset%gw_invalid_freq
 w_info%mat_type = MAT_INV_EPSILON
 call W%init(w_info, cryst, pp_mesh, gsph_c, vcp, w_fname, mqmem, dtset%npweps, &
             dtset%iomode, ngfftf, nfftf, nsppol, nspden, rhor, dtset%prtvol, comm)

 ABI_FREE(rhor)
 ABI_FREE(qlwl)

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ! Allocate g-vectors centered on k, k+q, k-p, and k+q-p
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))
 ABI_MALLOC(kg_kmp, (3, mpw))
 ABI_MALLOC(kg_kqmp, (3, mpw))

 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

 ABI_MALLOC(ur_k, (wfd%nfft*nspinor))
 ABI_MALLOC(ur_kq, (wfd%nfft*nspinor))
 ABI_MALLOC(ur_kmp, (wfd%nfft*nspinor))
 ABI_MALLOC(ur_kqmp, (wfd%nfft*nspinor))
 ABI_MALLOC(cwork_ur, (wfd%nfft*nspinor))
 ABI_MALLOC(cg_kmp, (2, mpw*nspinor))
 ABI_MALLOC(cg_kqmp, (2, mpw*nspinor))
 ! First order change (full term including active space)
 ABI_MALLOC(cg1_kqmp, (2, mpw*nspinor))

 ! Find correspondence IBZ_k --> IBZ
 ! Assume qptopt == kptopt unless value is specified in input
 !qptopt = ebands%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
 !qtimrev = kpts_timrev_from_kptopt(qptopt)

 nbsum = dtset%mband
 nbsum = 1 ! DEBUG

 stern_use_cache = merge(.True., .False., dtset%eph_stern == 1)

 ! =================================================
 ! Loop over spins and q-points (usually in the IBZ)
 ! =================================================
 do my_spin=1,gstore%my_nspins
   spin = gstore%my_spins(my_spin)
   gqk => gstore%gqk(my_spin)
   do my_iq=1,gqk%my_nq
     if (my_iq > 3) cycle ! DEBUG
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)
     qq_is_gamma = sum(qpt**2) < tol14

     ! Note symrec conventions here.
     iq_bz = gqk%my_q2bz(my_iq)
     iq_ibz = gqk%my_q2ibz(1, my_iq) !; isym_q = gqk%my_q2ibz(2, my_iq)
     !trev_q = gqk%my_q2ibz(6, my_iq); g0_q = gqk%my_q2ibz(3:5, my_iq)
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     !tsign_q = 1; if (trev_q == 1) tsign_q = -1
     qq_ibz = gstore%qibz(:, iq_ibz)
     mapc_qq = gqk%my_q2ibz(:, my_iq)
     !ABI_CHECK(isirr_q, "The q-point is usually in the IBZ!")

     !if (done_qbz_spin(iq_bz, spin) == 1) cycle
     call cwtime(cpu_ks, wall_ks, gflops_ks, "start")
     call wrtout(std_out, sjoin(" Computing g^GW(k,q) for q:", ktoa(qpt), "spin:", itoa(spin)))
     !print *, "iq_ibz:", iq_ibz, "qpt:", qpt, "qq_ibz:", qq_ibz

     ! Compute phonons for this qpt.
     call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red)

     ! ====================================
     ! Get DFPT potentials for this q-point
     ! ====================================
     ! After this branch we have allocated:
     !   v1scf(cplex, nfftf, nspden, my_npert))
     !   vxc1(cplex, nfft, nspden, my_npert)

     if (use_ftinterp) then
       ! Use Fourier interpolation to get DFPT potentials and DFPT densities for this qpt.
       call dvdb%get_ftqbz(cryst, qpt, qq_ibz, mapc_qq, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)
       !call drhodb%get_rho1_ftqbz(cryst, qpt, qq_ibz, mapc_qq, drho_cplex, nfftf, ngfftf, vxc1, gqk%pert_comm%value)

     else
       ! Read and reconstruct the dvscf potentials and the densities for this qpt and my_npert perturbations.
       db_iqpt = dvdb%findq(qq_ibz)
       ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file."))
       ! The first entry in mapc_qq2dvdb gives the index in dvdb%qpts.
       ! The other entries in mapc_qq are OK as they refer to symmetries
       mapc_qq2dvdb = mapc_qq; mapc_qq2dvdb(1) = db_iqpt
       call dvdb%readsym_qbz(cryst, qpt, mapc_qq2dvdb, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)

       !db_iqpt = drhodb%findq(qq_ibz)
       !ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DRHODB file."))
       !mapc_qq2dvdb = mapc_qq; mapc_qq2dvdb(1) = db_iqpt
       !call drhodb%readsym_rho1_qbz(cryst, qpt, mapc_qq2dvdb, drho_cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)
     end if
     !ABI_CHECK_IEQ(cplex, drho_cplex, "Different values of cplex for v1 and rho1!")

     ! Allocate vlocal1 with correct cplex. Note nvloc
     ABI_MALLOC_OR_DIE(vlocal1, (cplex*n4, n5, n6, gs_ham_kq%nvloc, my_npert), ierr)

     ! ========================================================================
     ! Loop over k-points in the e-ph matrix elements (usually in the full BZ).
     ! ========================================================================
     do my_ik=1,gqk%my_nk
       if (my_ik > 12) cycle ! DEBUG
       ! Symmetry indices for kk.
       kk = gqk%my_kpts(:, my_ik)
       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       ik_ibz = gqk%my_k2ibz(1, my_ik) ; isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5,my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       mapl_k = gqk%my_k2ibz(:, my_ik)
       print *, "my_ik", my_ik, " of my_nk:", gqk%my_nk

       kk_ibz = ebands%kptns(:,ik_ibz)
       istwf_k_ibz = wfd%istwfk(ik_ibz); npw_k_ibz = wfd%npwarr(ik_ibz)
       !print *, "ik_ibz:", ik_ibz, "kk:", kk, "kk_ibz:", kk_ibz

       ! Get npw_k, kg_k for k
       call wfd%get_gvec_gbound(cryst%gmet, ecut, kk, ik_ibz, isirr_k, istwf_k, npw_k, kg_k, gbound_k)

       ! Compute k+G vectors
       nkpg = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk, nkpg, npw_k)

       ! Compute nonlocal form factors ffnl_k at (k+G)
       ABI_MALLOC(ffnl_k, (npw_k, 1, psps%lmnmax, psps%ntypat))

       call mkffnl_objs(cryst, psps, 1, ffnl_k, ider0, idir0, kg_k, kpg_k, kk, nkpg, npw_k, ylm_k, ylmgr_dum, &
                        comm=gqk%pert_comm%value, request=ffnl_k_request)

       ! Find k + q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered as:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qq_ibz

       if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, 1, kq, mapl_kq) /= 0) then
         write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k+q could not be generated from a symmetrical one.",trim(ltoa(kq))
         ABI_ERROR(msg)
       end if
       ikq_ibz = mapl_kq(1); isym_kq = mapl_kq(2); trev_kq = mapl_kq(6); g0_kq = mapl_kq(3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)
       istwf_kq_ibz = wfd%istwfk(ikq_ibz); npw_kq_ibz = wfd%npwarr(ikq_ibz)

       ! Get npw_kq, kg_kq for k+q.
       call wfd%get_gvec_gbound(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, istwf_kq, npw_kq, kg_kq, gbound_kq)

       ! Compute k+q+G vectors
       nkpg1 = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl_kq at (k+q+G)
       ABI_MALLOC(ffnl_kq, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl_objs(cryst, psps, 1, ffnl_kq, ider0, idir0, kg_kq, kpg1_k, kq, nkpg1, npw_kq, ylm_kq, ylmgr_kq, &
                        comm=gqk%pert_comm%value, request=ffnl_kq_request)

       ! Double loop over the (m, n) bands indices in the <mm_kq,kq|Delta_{q,nu}Sigma|k,nn_k> elements.
       ABI_MALLOC(ug_k, (2, npw_k*nspinor))
       ABI_MALLOC(ug_kq, (2, npw_kq*nspinor))

       do mm_kq=1,1
         call wfd%rotate_cg(mm_kq, ndat1, spin, kq_ibz, npw_kq, kg_kq, istwf_kq, &
                            cryst, mapl_kq, gbound_kq, work_ngfft, work, ug_kq, urs_kbz=ur_kq)
       end do ! mm_kq
       do nn_k=1,1
         call wfd%rotate_cg(nn_k, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                            cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_k)
       end do ! nn_k
       !if (cplex == 1)
       !  g_xc(mm_kq, nn_k) = sum(GWPC_CONJG(ur_kq) * ur_k * vxc1(1,:,ispin, imyip)) / nfftf
       !else

       ! =======================================
       ! Sum over the momenta pp in the full BZ.
       ! =======================================
       ! Be careful here because pp should run over the list of wavevectors in the screening!
       ! as pp_mesh%bz is not necessarily equivalent to the k-mesh for the wavefunctions.
       ! and we have to use ipp_bz to get/symmetrize W(pp_bz).
       do ipp_bz=1,pp_mesh%nbz
         pp = pp_mesh%bz(:,ipp_bz)
         pp_is_gamma = sum(pp**2) < tol14
         !print *, "Begin sum over pp: ", trim(ktoa(pp))
         if (ipp_bz > 12) cycle ! DEBUG

         ! Symmetry tables and g-sphere centered on k-p
         kmp = kk - pp
         if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, 1, kmp, mapl_kmp) /= 0) then
           write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k-p could not be generated from a symmetrical one.",trim(ltoa(kmp))
           ABI_ERROR(msg)
         end if
         ikmp_ibz = mapl_kmp(1); isym_kmp = mapl_kmp(2); trev_kmp = mapl_kmp(6); g0_kmp = mapl_kmp(3:5)
         isirr_kmp = (isym_kmp == 1 .and. trev_kmp == 0 .and. all(g0_kmp == 0))
         kmp_ibz = ebands%kptns(:, ikmp_ibz)
         istwf_kmp_ibz = wfd%istwfk(ikmp_ibz); npw_kmp_ibz = wfd%npwarr(ikmp_ibz)
         ! Get npw_kmp, kg_kmp for k-p
         call wfd%get_gvec_gbound(cryst%gmet, ecut, kmp, ikmp_ibz, isirr_kmp, istwf_kmp, npw_kmp, kg_kmp, gbound_kmp)

         ! Symmetry tables and g-sphere centered on k+q-p
         kqmp = kq - pp
         if (kpts_map("symrel", ebands_timrev, cryst, gstore%krank_ibz, 1, kqmp, mapl_kqmp) /= 0) then
           write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k+q-p could not be generated from a symmetrical one.",trim(ltoa(kqmp))
           ABI_ERROR(msg)
         end if
         ikqmp_ibz = mapl_kqmp(1); isym_kqmp = mapl_kqmp(2); trev_kqmp = mapl_kqmp(6); g0_kqmp = mapl_kqmp(3:5)
         isirr_kqmp = (isym_kqmp == 1 .and. trev_kqmp == 0 .and. all(g0_kqmp == 0))
         kqmp_ibz = ebands%kptns(:, ikqmp_ibz)
         istwf_kqmp_ibz = wfd%istwfk(ikqmp_ibz); npw_kqmp_ibz = wfd%npwarr(ikqmp_ibz)
         ! Get npw_kqmp, kg_kqmp for k+q-p
         call wfd%get_gvec_gbound(cryst%gmet, ecut, kqmp, ikqmp_ibz, isirr_kqmp, istwf_kqmp, npw_kqmp, kg_kqmp, gbound_kqmp)

         ! This is the g-sphere for W_{gg'}(pp).
         ! Note that in this case, the sphere is Gamma-centered i.e. it does not depend on the pp wavevector
         npw_pp = W%npw; kg_pp => W%gvec
         call sphereboundary(gbound_pp, istwfk1, kg_pp, wfd%mgfft, npw_pp)

         ABI_MALLOC(crhog_bkmp_nk, (npw_pp))
         ABI_MALLOC(crhog_mkq_bkqmp, (npw_pp))

         ! Prepare the object for applying W_qbz.
         ! FIXME: Sq = q+G0 with non-zero G0 is not supported.
         call W%symmetrizer(ipp_bz, cryst, gsph_c, pp_mesh, vcp)
         !call W%w0gemv("N", in_npw, nspinor, only_diago, cone_gw, czero_gw, in_ket, out_ket)

         ! We need two stern_t objects to compute the first order change of the wavefunctions at k-p and k+q-p.
         ! Clearly we should not duplicate the work when pp = 0.
         ! When pp == 0, we also get g_ks via stern_solve
         ! Also, one should handle more carefully the integration in g_sigma around pp = Gamma in the case of semiconductors.

         stern_has_band_para = .False.
         !if (stern_has_band_para) then
         !  nband_me = sigma%my_bsum_stop - sigma%my_bsum_start + 1
         !  stern_comm = sigma%bsum_comm%value
         !else
           nband_me = nbsum
           stern_comm = xmpi_comm_self
         !end if

#if 0
         call stern_kmp%init(dtset, npw_kmp, npw_kqmp, nspinor, nbsum, nband_me, stern_use_cache, work_ngfft, mpi_enreg, stern_comm)

         ! =========================================================
         ! Get GS wavefunctions at k+q-p and store them in stern_kmp.
         ! =========================================================

         !do ib_sum=sigma%my_bsum_start, sigma%my_bsum_stop
         do ib_sum=1, nbsum
           call wfd%rotate_cg(ib_sum, ndat1, spin, kqmp_ibz, npw_kqmp, kg_kqmp, istwf_kqmp, &
                              cryst, mapl_kqmp, gbound_kqmp, work_ngfft, work, cg_kqmp)

           ! NB: cg_kqmp is dimensioned with mpw
           if (stern_kmp%has_band_para) then
             NOT_IMPLEMENTED_ERROR()
             ii = -1
             !ii = ib_sum - sigma%my_bsum_start + 1
             stern_kmp%cgq(:,:,ii) = cg_kqmp(:,1:npw_kqmp*nspinor)
           else
             stern_kmp%cgq(:,:,ib_sum) = cg_kqmp(:,1:npw_kqmp*nspinor)
             !print *, "stern_kmp%cgq:", stern_kmp%cgq(:, 1:3, ib_sum)
           end if
         end do
#endif

         ! =========================================================
         ! Get GS wavefunctions at ?? and store them in stern_kqmp.
         ! =========================================================
         !call stern_kqmp%init(dtset, npw_kmp, npw_kqmp, nspinor, nbsum, nband_me, stern_use_cache, work_ngfft, mpi_enreg, stern_comm)
         !call stern_kqmp%free()

         ! ==========================
         ! Sum over bands up to nbsum
         ! ==========================
         !do ib_sum=sigma%my_bsum_start, sigma%my_bsum_stop
         do ib_sum=1, nbsum
           ! <bsum,k-p| e^{-ip+G'}|n,k>
           call wfd%rotate_cg(ib_sum, ndat1, spin, kmp_ibz, npw_kmp, kg_kmp, istwf_kmp, &
                             cryst, mapl_kmp, gbound_kmp, work_ngfft, work, cg_kmp, urs_kbz=ur_kmp)

           cwork_ur = GWPC_CONJG(ur_kmp) * ur_k
           call fft_ur(npw_pp, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_pp, gbound_pp, cwork_ur, crhog_bkmp_nk)

           ! <m,k+q| e^{ip+G}|bsum,k+q-p> --> compute <k| e^{-i(q+G)}|k+q> with FFT and take the CC.
           call wfd%rotate_cg(ib_sum, ndat1, spin, kqmp_ibz, npw_kqmp, kg_kqmp, istwf_kqmp, &
                              cryst, mapl_kqmp, gbound_kqmp, work_ngfft, work, cg_kqmp, urs_kbz=ur_kqmp)

           cwork_ur = GWPC_CONJG(ur_kqmp) * ur_kqmp
           call fft_ur(npw_pp, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_pp, gbound_pp, cwork_ur, crhog_mkq_bkqmp)
           crhog_mkq_bkqmp = GWPC_CONJG(crhog_mkq_bkqmp)
         end do ! ib_sum

         ABI_FREE(crhog_bkmp_nk)
         ABI_FREE(crhog_mkq_bkqmp)

         ! ========================================
         ! Loop over my set of atomic perturbations
         ! ========================================
         ! For each perturbation:
         !    - setup H1 from vlocal1.
         !    For each band in band_sum:
         !        - Solve the Sternheimer non-self-consistently and get the KS e-ph matrix elements
         !        - Build the full first-order wavefunction including the active subspace.

         if (ffnl_k_request /= xmpi_request_null) call xmpi_wait(ffnl_k_request, ierr)
         if (ffnl_kq_request /= xmpi_request_null) call xmpi_wait(ffnl_kq_request, ierr)

         do imyp=1,gqk%my_npert
           idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)

           ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
           ! Each CPU prepares its own potentials.
           call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_ham_kq%nvloc, &
                                      pawfgr, mpi_enreg, vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))

           ! Continue to initialize the Hamiltonian (call it here to support dfpt_cgwf Sternheimer).
           call gs_ham_kq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

           ! Prepare application of the NL part.
           call init_rf_hamiltonian(cplex, gs_ham_kq, ipert, rf_ham_kq, has_e1kbsc=.true.)
           call rf_ham_kq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

           ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
           call getgh1c_setup(gs_ham_kq, rf_ham_kq, dtset, psps, kk, kq, idir, ipert, &  ! In
             cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
             npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
             dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnl_k, ffnl_kq, ph3d, ph3d1, &  ! Out
             reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

           ! Compute H(1) applied to GS wavefunction Psi_nk(0)
           !do ib_k=1,nbcalc_ks
           !  if (gwpt%bsum_comm%skip(ib_k, root=root_bcalc(ib_k))) cycle ! MPI parallelism inside bsum_comm
           !                                                              ! Store rank treating ib_k in root_bcalc
           !  band_ks = ib_k + bstart_ks - 1
           !  eig0nk = ebands%eig(band_ks, ik_ibz, spin)
           !  ! Use scissor shift on 0-order eigenvalue
           !  eshift = eig0nk - dtset%dfpt_sciss
           !
           !  call getgh1c(berryopt0, kets_k(:,:,ib_k), cwaveprj0, h1kets_kq(:,:,imyp, ib_k), &
           !    grad_berry, gs1c, gs_ham_kq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
           !    optnl, opt_gvnlx1, rf_ham_kq, sij_opt, tim_getgh1c1, usevnl)
           !end do

           ! ================
           ! NSCF Sternheimer
           ! ================
           ! This is just to test stern_kmp%solve for pp == 0, other values of pp won't work.
#if 0
           if (pp_is_gamma .and. nbsum /= 1) then
           do ib_sum=1, nbsum
             stern_kmp%bands_treated_now(:) = 0; stern_kmp%bands_treated_now(ib_sum) = 1

             call wfd%rotate_cg(ib_sum, ndat1, spin, kmp_ibz, npw_kmp, kg_kmp, istwf_kmp, &
                                cryst, mapl_kmp, gbound_kmp, work_ngfft, work, cg_kmp)

             stern_kmp%rank_band = 0; u1_band = ib_sum; band_me = ib_sum
             if (stern_kmp%has_band_para) then
               NOT_IMPLEMENTED_ERROR()
             end if
             call stern_kmp%solve(u1_band, band_me, idir, ipert, qpt, gs_ham_kq, rf_ham_kq, &
                                  ebands%eig(:,ikmp_ibz,spin), ebands%eig(:,ikqmp_ibz,spin), &
                                  cg_kmp, cwaveprj0, cg1_kqmp, cwaveprj, msg, ierr) ! , &
                                  !full_cg1=full_cg1_kmpq, full_ur1=full_ur1_kmpq)
             ABI_CHECK(ierr == 0, msg)
           end do ! ibsum
           end if
#endif

           call rf_ham_kq%free()
           ABI_FREE(kinpw1)
           ABI_FREE(dkinpw)
           ABI_FREE(ph3d)
           ABI_SFREE(ph3d1)
         end do ! imyp

         !call stern_kmp%free()
       end do ! ipp_bz

       ABI_FREE(kpg_k)
       ABI_FREE(ffnl_k)
       ABI_FREE(kpg1_k)
       ABI_FREE(ffnl_kq)
       !
       ABI_FREE(ug_k)
       ABI_FREE(ug_kq)
     end do ! my_ik

     ABI_FREE(vlocal1)
     ABI_FREE(v1scf)
     call cwtime_report(" One q-point", cpu_ks, wall_ks, gflops_ks)
   end do ! iq_ibz
 end do ! my_spin

 call cwtime_report(" gwpt_eph full calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Free memory
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(kg_kmp)
 ABI_FREE(kg_kqmp)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(cg_work)
 ABI_FREE(ur_k)
 ABI_FREE(ur_kq)
 ABI_FREE(ur_kmp)
 ABI_FREE(ur_kqmp)
 ABI_FREE(cwork_ur)
 ABI_FREE(cg_kmp)
 ABI_FREE(cg_kqmp)
 ABI_FREE(cg1_kqmp)
 !ABI_FREE(ihave_ikibz_spin)
 !ABI_FREE(grad_berry)
 ABI_FREE(vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_FREE(gbound_k)
 ABI_FREE(gbound_kq)
 ABI_FREE(gbound_kmp)
 ABI_FREE(gbound_kqmp)
 ABI_FREE(gbound_pp)

 call gs_ham_kq%free()
 call wfd%free()
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 call vcp%free()
 call w%free()
 call pp_mesh%free()
 call gsph_c%free()
 call gstore%free()
 call gwpt%free()

 ! This to make sure that the parallel output of GSTORE is completed
 call xmpi_barrier(comm)
 call cwtime_report(" gwpt_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

end subroutine gwpt_run
!!***

!----------------------------------------------------------------------

!!****f* m_gwpt/gwpt_new
!! NAME
!!  gwpt_new
!!
!! FUNCTION
!!  Creation method (allocates memory, initialize data from input vars).
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  ecut=Cutoff energy for wavefunctions.
!!  cryst<crystal_t>=Crystalline structure
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!!  dtfil<datafiles_type>=variables related to files.
!!  comm=MPI communicator
!!
!! SOURCE

type(gwpt_t) function gwpt_new(dtset, ecut, cryst, ebands, ifc, dtfil, qrank_ibz, comm) result(gwpt)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
 type(datafiles_type),intent(in) :: dtfil
 type(krank_t),intent(out) :: qrank_ibz

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, istwfk1 = 1
 integer :: my_rank,my_nshiftq,cnt,nprocs,ik_ibz,ndeg, iq_ibz, qptopt, qtimrev
 integer :: ii, ierr, spin, gap_err, ikcalc, qprange_, bstop !it, ipw, i1,i2,i3,
 integer :: jj, bstart, natom, natom3
 integer :: isym_k, trev_k, mband, nrest, color
 character(len=5000) :: msg
 real(dp) :: cpu_all, wall_all, gflops_all, cpu, wall, gflops
 logical :: changed, isirr_k
 type(gaps_t) :: gaps
 type(krank_t) :: krank_ibz
!arrays
 integer :: g0_k(3), units(2), indkk_k(6,1), qptrlatt(3,3)  !band_block(2),
 integer,allocatable :: temp(:,:), degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:), iperm(:)
 real(dp):: my_shiftq(3,1), kk(3) !, kq(3)
#ifdef HAVE_MPI
 integer,parameter :: ndims = 5
 integer :: comm_cart, me_cart
 logical :: reorder
 integer :: dims(ndims)
 logical :: periods(ndims), keepdim(ndims)
#endif

! *************************************************************************

 ABI_UNUSED(dtfil%ireadden)
 ABI_UNUSED(ifc%natom)


 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 call cwtime(cpu, wall, gflops, "start")

 units = [std_out, ab_out]

 ! Copy important dimensions.
 gwpt%nsppol = ebands%nsppol; gwpt%nspinor = ebands%nspinor; mband = dtset%mband
 natom = cryst%natom; natom3 = cryst%natom * 3

 ! Broadening parameter from zcut
 gwpt%ieta = + j_dpc * dtset%zcut

 ! Define q-mesh for integration of the self-energy.
 ! Either q-mesh from DVDB (no interpolation) or eph_ngqpt_fine (Fourier interpolation if q not in DDB)
 gwpt%ngqpt = dtset%ddb_ngqpt; my_nshiftq = 1; my_shiftq(:,1) = dtset%ddb_shiftq
 if (all(dtset%eph_ngqpt_fine /= 0)) then
   gwpt%ngqpt = dtset%eph_ngqpt_fine; my_shiftq = 0
 end if

 ! Setup IBZ, weights and BZ.
 ! Assume qptopt == kptopt unless value is specified in input
 qptrlatt = 0; qptrlatt(1, 1) = gwpt%ngqpt(1); qptrlatt(2, 2) = gwpt%ngqpt(2); qptrlatt(3, 3) = gwpt%ngqpt(3)
 qptopt = ebands%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
 qtimrev = kpts_timrev_from_kptopt(qptopt)
 call kpts_ibz_from_kptrlatt(cryst, qptrlatt, qptopt, my_nshiftq, my_shiftq, &
                             gwpt%nqibz, gwpt%qibz, gwpt%wtq, gwpt%nqbz, gwpt%qbz, bz2ibz=gwpt%ind_qbz2ibz)

 ! HM: the bz2ibz produced above is incomplete, I do it here using listkk
 ABI_MALLOC(temp, (6, gwpt%nqbz))

 qrank_ibz = krank_from_kptrlatt(gwpt%nqibz, gwpt%qibz, qptrlatt, compute_invrank=.False.)
 if (kpts_map("symrec", qtimrev, cryst, qrank_ibz, gwpt%nqbz, gwpt%qbz, temp) /= 0) then
   ABI_ERROR("Cannot map qBZ to qIBZ!")
 end if

 gwpt%ind_qbz2ibz(1,:) = temp(1,:)
 gwpt%ind_qbz2ibz(2,:) = temp(2,:)
 gwpt%ind_qbz2ibz(3,:) = temp(6,:)
 gwpt%ind_qbz2ibz(4,:) = temp(3,:)
 gwpt%ind_qbz2ibz(5,:) = temp(4,:)
 gwpt%ind_qbz2ibz(6,:) = temp(5,:)
 ABI_FREE(temp)
!END DEBUG

 gaps = ebands_get_gaps(ebands, gap_err)

 ! ======================================================
 ! Select k-point and bands where corrections are wanted
 ! ======================================================
 !
 ! if symsigma == +1, we have to include all degenerate states in the set
 ! because the final QP corrections will be obtained by averaging the results in the degenerate subspace.
 ! We initialize IBZ(k) here so that we have all the basic dimensions of the run and it's possible
 ! to distribuite the calculations among processors.
 gwpt%symsigma = dtset%symsigma; gwpt%timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 call cwtime_report(" gwpt_new: k-points", cpu, wall, gflops)

 ! TODO: nkcalc should be spin dependent (similar piece of code in m_gwr).
 if (dtset%nkptgw /= 0) then

   ! Treat the k-points and bands specified in the input file via kptgw and bdgw.
   call sigtk_kcalc_from_nkptgw(dtset, mband, gwpt%nkcalc, gwpt%kcalc, gwpt%bstart_ks, gwpt%nbcalc_ks)

 else

   if (any(abs(dtset%sigma_erange) > zero)) then
     ! Use sigma_erange and (optionally) gwpt_ngkpt
     call sigtk_kcalc_from_erange(dtset, cryst, ebands, gaps, gwpt%nkcalc, gwpt%kcalc, gwpt%bstart_ks, gwpt%nbcalc_ks, comm)

   else
     ! Use qp_range to select the interesting k-points and the corresponding bands.
     !
     !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
     ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
     !          bands above and below the Fermi level.
     ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
     !          Include all occupied states and `num` empty states.

     qprange_ = dtset%gw_qprange
     if (gap_err /= 0 .and. qprange_ == 0) then
       ABI_WARNING("Cannot compute fundamental and direct gap (likely metal). Will replace qprange 0 with qprange 1")
       qprange_ = 1
     end if

     if (qprange_ /= 0) then
       call sigtk_kcalc_from_qprange(dtset, cryst, ebands, qprange_, gwpt%nkcalc, gwpt%kcalc, gwpt%bstart_ks, gwpt%nbcalc_ks)
     else
       ! qprange is not specified in the input.
       ! Include direct and fundamental KS gap or include states depending on the position wrt band edges.
       call sigtk_kcalc_from_gaps(dtset, ebands, gaps, gwpt%nkcalc, gwpt%kcalc, gwpt%bstart_ks, gwpt%nbcalc_ks)
     end if
   end if

 end if ! nkptgw /= 0

 ! The k-point and the symmetries connecting the BZ k-point to the IBZ.
 ABI_MALLOC(gwpt%kcalc2ibz, (gwpt%nkcalc, 6))
 if (abs(gwpt%symsigma) == 1) then
   ABI_MALLOC(gwpt%degtab, (gwpt%nkcalc, gwpt%nsppol))
 end if

 ! Workspace arrays used to compute degeneracy tables.
 ABI_ICALLOC(degblock_all, (2, mband, gwpt%nkcalc, gwpt%nsppol))
 ABI_ICALLOC(ndeg_all, (gwpt%nkcalc, gwpt%nsppol))

 krank_ibz = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)
 ierr = 0

 do ikcalc=1,gwpt%nkcalc
   if (mod(ikcalc, nprocs) /= my_rank) then
     gwpt%kcalc2ibz(ikcalc, :) = 0
     gwpt%bstart_ks(ikcalc, :) = 0
     gwpt%nbcalc_ks(ikcalc, :) = 0
     cycle ! MPI parallelism inside comm
   end if

   ! Note symrel and use_symrel.
   ! These are the conventions for the symmetrization of the wavefunctions used in cgtk_rotate.
   kk = gwpt%kcalc(:, ikcalc)

   if (kpts_map("symrel", gwpt%timrev, cryst, krank_ibz, 1, kk, indkk_k) /= 0) then
      write(msg, '(11a)' )&
       "The WFK file cannot be used to compute self-energy corrections at k-point: ",trim(ktoa(kk)),ch10,&
       "The k-point cannot be generated from a symmetrical one.", ch10,&
       "q-mesh: ",trim(ltoa(gwpt%ngqpt)),", k-mesh (from kptrlatt): ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
       'Action: check your WFK file and the (k, q) point input variables.'
      ABI_ERROR(msg)
   end if

   ! TODO: Invert dims and update abipy
   gwpt%kcalc2ibz(ikcalc, :) = indkk_k(:, 1)

   ik_ibz = indkk_k(1,1); isym_k = indkk_k(2,1)
   trev_k = indkk_k(6, 1); g0_k = indkk_k(3:5,1)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   !kk_ibz = ebands%kptns(:,ik_ibz)
   if (.not. isirr_k) then
     ABI_WARNING(sjoin("The k-point in gwpt_{nk} must be in the IBZ but got:", ktoa(kk)))
     ierr = ierr + 1
   end if

   ! We will have to average the QP corrections over degenerate states if symsigma=1 is used.
   ! Here we make sure that all the degenerate states are included.
   ! Store also band indices of the degenerate sets, used to average final results.
   if (abs(gwpt%symsigma) == 1) then
     cnt = 0
     do spin=1,gwpt%nsppol
       bstop = gwpt%bstart_ks(ikcalc, spin) + gwpt%nbcalc_ks(ikcalc, spin) - 1
       call ebands_enclose_degbands(ebands, ik_ibz, spin, gwpt%bstart_ks(ikcalc, spin), bstop, changed, TOL_EDIFF, &
                                    degblock=degblock)
       if (changed) then
         gwpt%nbcalc_ks(ikcalc, spin) = bstop - gwpt%bstart_ks(ikcalc, spin) + 1
         cnt = cnt + 1
         if (cnt < 5) then
           write(msg,'(2(a,i0),2a,2(1x,i0))') &
             "Not all the degenerate states for ikcalc: ",ikcalc,", spin: ",spin,ch10, &
             "were included in the bdgw set. bdgw has been automatically changed to: ",gwpt%bstart_ks(ikcalc, spin), bstop
           ABI_COMMENT(msg)
         end if
         write(msg,'(2(a,i0),2a)') &
           "The number of included states: ", bstop, &
           " is larger than the number of bands in the input ",dtset%nband(ik_ibz + (spin-1)*ebands%nkpt),ch10,&
           "Action: Increase nband."
         ABI_CHECK(bstop <= dtset%nband(ik_ibz + (spin-1)*ebands%nkpt), msg)
       end if

       ! Store band indices used for averaging (shifted by bstart_ks)
       ndeg = size(degblock, dim=2)
       ndeg_all(ikcalc, spin) = ndeg
       degblock_all(:, 1:ndeg, ikcalc, spin) = degblock(:, 1:ndeg)

       ABI_FREE(degblock)
     end do
   end if ! symsigma
 end do ! ikcalc

 call krank_ibz%free()
 ABI_CHECK(ierr == 0, "kptgw wavevectors must be in the IBZ read from the WFK file.")

 ! Collect data
 call xmpi_sum(gwpt%kcalc2ibz, comm, ierr)
 call xmpi_sum(gwpt%bstart_ks, comm, ierr)
 call xmpi_sum(gwpt%nbcalc_ks, comm, ierr)

 ! Build degtab tables.
 if (abs(gwpt%symsigma) == 1) then
   call xmpi_sum(ndeg_all, comm, ierr)
   call xmpi_sum(degblock_all, comm, ierr)
   do ikcalc=1,gwpt%nkcalc
     do spin=1,gwpt%nsppol
       ndeg = ndeg_all(ikcalc, spin)
       ABI_MALLOC(gwpt%degtab(ikcalc, spin)%bids, (ndeg))
       do ii=1,ndeg
         cnt = degblock_all(2, ii, ikcalc, spin) - degblock_all(1, ii, ikcalc, spin) + 1
         ABI_MALLOC(gwpt%degtab(ikcalc, spin)%bids(ii)%vals, (cnt))
         gwpt%degtab(ikcalc, spin)%bids(ii)%vals = [(jj, jj= &
           degblock_all(1, ii, ikcalc, spin) - gwpt%bstart_ks(ikcalc, spin) + 1, &
           degblock_all(2, ii, ikcalc, spin) - gwpt%bstart_ks(ikcalc, spin) + 1)]
       end do
     end do
   end do
 end if
 ABI_FREE(degblock_all)
 ABI_FREE(ndeg_all)

 call cwtime_report(" gwpt_new: kptgw", cpu, wall, gflops)

 ! Now we can finally compute max_nbcalc
 gwpt%max_nbcalc = maxval(gwpt%nbcalc_ks)

 ABI_MALLOC(gwpt%bstop_ks, (gwpt%nkcalc, gwpt%nsppol))
 gwpt%bstop_ks = gwpt%bstart_ks + gwpt%nbcalc_ks - 1

 ! Define number of bands included in self-energy summation as well as the band range.
 ! This value depends on the kind of calculation as imag_only can take advantage of
 ! the energy window aroud the band edges.
 !
 ! Notes about MPI version.
 !
 ! Loops are MPI parallelized over q-points
 ! wavefunctions are NOT distributed but only states between my_bsum_start and my_bsum_stop
 ! are allocated and read from file.
 ! perturbations and q-points in the IBZ can also be distributed.

 gwpt%bsum_start = 1; gwpt%bsum_stop = mband
 if (all(dtset%sigma_bsum_range /= 0)) then
   gwpt%bsum_start = max(dtset%sigma_bsum_range(1), 1)
   gwpt%bsum_stop = min(dtset%sigma_bsum_range(2), mband)
 end if
 gwpt%nbsum = gwpt%bsum_stop - gwpt%bsum_start + 1

 ! ========================
 ! === MPI DISTRIBUTION ===
 ! ========================
 ! Init for sequential execution.
 gwpt%my_npert = natom3

 if (any(dtset%eph_np_pqbks /= 0)) then
   ! Use parameters from input file.
   gwpt%pert_comm%nproc = dtset%eph_np_pqbks(1)
   gwpt%qpt_comm%nproc  = dtset%eph_np_pqbks(2)
   gwpt%bsum_comm%nproc = dtset%eph_np_pqbks(3)
   gwpt%kcalc_comm%nproc = dtset%eph_np_pqbks(4)
   gwpt%spin_comm%nproc = dtset%eph_np_pqbks(5)
   gwpt%my_npert = natom3 / gwpt%pert_comm%nproc
   ABI_CHECK(gwpt%my_npert > 0, "pert_comm_nproc cannot be greater than 3 * natom.")
   ABI_CHECK(mod(natom3, gwpt%pert_comm%nproc) == 0, "pert_comm_nproc must divide 3 * natom.")
 else
   ! Automatic grid generation.
   if (gwpt%pert_comm%nproc == 1) then
     ! Try again with more procs.
     do cnt=natom3,2,-1
       if (mod(nprocs, cnt) == 0 .and. mod(natom3, cnt) == 0) then
         gwpt%pert_comm%nproc = cnt; gwpt%my_npert = natom3 / cnt; exit
       end if
     end do
   end if

   if (gwpt%my_npert == natom3 .and. nprocs > 1) then
     ABI_WARNING("The number of MPI procs should be divisible by 3*natom to reduce memory requirements!")
   end if

   ! Define number of procs for q-points and bands. nprocs is divisible by pert_comm%nproc.
   ! Try to distribute equally nbsum first.
   nrest = nprocs / gwpt%pert_comm%nproc
   do bstop=nrest,1,-1
     if (mod(gwpt%nbsum, bstop) == 0 .and. mod(nprocs, gwpt%pert_comm%nproc * bstop) == 0) then
       gwpt%bsum_comm%nproc = bstop; gwpt%qpt_comm%nproc = nrest / gwpt%bsum_comm%nproc
       exit
     end if
   end do
 end if

 ! Consistency check.
 if (gwpt%pert_comm%nproc * gwpt%qpt_comm%nproc * gwpt%bsum_comm%nproc * gwpt%kcalc_comm%nproc * gwpt%spin_comm%nproc /= nprocs) then
   write(msg, "(a,i0,3a, 6(a,1x,i0))") &
     "Cannot create 5d Cartesian grid with total nprocs: ", nprocs, ch10, &
     "Idle processes are not supported. The product of the `nprocs_*` vars should be equal to nprocs.", ch10, &
     "pert_nproc (", gwpt%pert_comm%nproc, ") x qpt_nproc (", gwpt%qpt_comm%nproc, ") x bsum_nproc (", gwpt%bsum_comm%nproc, &
     ") x kcalc_nproc (", gwpt%kcalc_comm%nproc, ") x spin_nproc (", gwpt%spin_comm%nproc, ") != ", nprocs
   ABI_ERROR(msg)
 end if

 gwpt%coords_pqbks = 0
#ifdef HAVE_MPI
 ! Create 5d cartesian communicator: 3*natom perturbations, q-points in IBZ, bands in gwpt sum, kpoints in gwpt_k, spins
 ! FIXME: Fix spin
 periods(:) = .False.; reorder = .False.
 dims = [gwpt%pert_comm%nproc, gwpt%qpt_comm%nproc, gwpt%bsum_comm%nproc, gwpt%kcalc_comm%nproc, gwpt%spin_comm%nproc]
 ! Try New distrib ?
 !dims = [gwpt%pert_comm%nproc, gwpt%bsum_comm%nproc, gwpt%qpt_comm%nproc, gwpt%kcalc_comm%nproc, gwpt%spin_comm%nproc]

 call MPI_CART_CREATE(comm, ndims, dims, periods, reorder, comm_cart, ierr)
 ! Find the index and coordinates of the current processor
 call MPI_COMM_RANK(comm_cart, me_cart, ierr)
 call MPI_CART_COORDS(comm_cart, me_cart, ndims, gwpt%coords_pqbks, ierr)

 ! Create communicator to distribute natom3 perturbations.
 keepdim = .False.; keepdim(1) = .True.; call gwpt%pert_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for qpoints in self-energy integration.
 keepdim = .False.; keepdim(2) = .True.; call gwpt%qpt_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for bands for self-energy summation
 keepdim = .False.; keepdim(3) = .True.; call gwpt%bsum_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for kpoints.
 keepdim = .False.; keepdim(4) = .True.; call gwpt%kcalc_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for spins.
 keepdim = .False.; keepdim(5) = .True.; call gwpt%spin_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for the (band_sum, qpoint_sum) loops
 !keepdim = .False.; keepdim(2:3) = .True.; call gwpt%qb_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for the (perturbation, band_sum, qpoint_sum)
 !keepdim = .False.; keepdim(1:3) = .True.; call gwpt%pqb_comm%from_cart_sub(comm_cart, keepdim)

 call xmpi_comm_free(comm_cart)
#endif

 ! Distribute k-points and create mapping to ikcalc index.
 call xmpi_split_cyclic(gwpt%nkcalc, gwpt%kcalc_comm%value, gwpt%my_nkcalc, gwpt%my_ikcalc)
 ABI_CHECK(gwpt%my_nkcalc > 0, sjoin("nkcalc (", itoa(gwpt%nkcalc), ") < kcalc_comm_nproc (", itoa(gwpt%kcalc_comm%nproc), ")"))

 ! Distribute spins and create mapping to spin index.
 if (gwpt%nsppol == 2) then
   call xmpi_split_block(gwpt%nsppol, gwpt%spin_comm%value, gwpt%my_nspins, gwpt%my_spins)
   ABI_CHECK(gwpt%my_nspins > 0, sjoin("nsppol (", itoa(gwpt%nsppol), ") < spin_comm_nproc (", itoa(gwpt%spin_comm%nproc), ")"))
 else
   ! No nsppol parallelism DOH!
   gwpt%my_nspins = 1
   ABI_MALLOC(gwpt%my_spins, (gwpt%my_nspins))
   gwpt%my_spins = 1
 end if

 ! Create MPI communicator for parallel netcdf IO used to write results for the different k-points.
 ! This communicator is defined only on the processes that will perform IO.
 call gwpt%ncwrite_comm%set_to_null()

 if (gwpt%kcalc_comm%nproc == 1 .and. gwpt%spin_comm%nproc == 1) then
   ! Easy-peasy: only master in comm_world performs IO.
   if (my_rank == master) call gwpt%ncwrite_comm%set_to_self()
 else
    ! Create subcommunicator by selecting one proc per kpoint-spin subgrid.
    ! Since we write to ab_out in gwpt_gather_and_write, make sure that ab_out is connected!
    ! This means gwpt_nk resuls will be spread among multiple ab_out files.
    ! Only SIGPEPH.nc will contain all the results.
    ! Remember that now all nc define operations must be done inside ncwrite_comm
    ! Obviously I'm assuming HDF5 + MPI-IO
    !
    ! NB: If MPI_UNDEFINED is passed as the colour value, the subgroup in which the calling
    ! MPI process will be placed is MPI_COMM_NULL

    color = xmpi_undefined; if (all(gwpt%coords_pqbks(1:3) == 0)) color = 1
    call xmpi_comm_split(comm, color, my_rank, gwpt%ncwrite_comm%value, ierr)
    if (color == 1) then
      gwpt%ncwrite_comm%me = xmpi_comm_rank(gwpt%ncwrite_comm%value)
      gwpt%ncwrite_comm%nproc = xmpi_comm_size(gwpt%ncwrite_comm%value)
      if (my_rank == master) then
        call wrtout(units, &
          sjoin("- Using parallelism over k-points/spins. Cannot write full results to main output", ch10, &
                "- All procs except master will write to dev_null. Use SIGEPH.nc to analyze results."))
        !write(std_out, *)"ncwrite_comm_me:", gwpt%ncwrite_comm%me, "ncwrite_comm%nproc:", gwpt%ncwrite_comm%nproc
      end if
      if (.not. is_open(ab_out)) then
        !if (open_file(strcat(dtfil%filnam_ds(2), "_rank_", itoa(gwpt%ncwrite_comm%me)), msg, unit=ab_out, &
        if (open_file(NULL_FILE, msg, unit=ab_out, form="formatted", action="write", status='unknown') /= 0) then
          ABI_ERROR(msg)
        end if
      end if
    else
      call gwpt%ncwrite_comm%set_to_null()
    end if
 end if

 ! Build table with list of perturbations treated by this CPU inside pert_comm
 call ephtk_set_pertables(cryst%natom, gwpt%my_npert, gwpt%pert_table, gwpt%my_pinfo, gwpt%pert_comm%value)

 ! Split bands among the procs inside bsum_comm using block distribution.
 call xmpi_split_work(gwpt%nbsum, gwpt%bsum_comm%value, gwpt%my_bsum_start, gwpt%my_bsum_stop)
 if (gwpt%my_bsum_start == gwpt%nbsum + 1) then
   ABI_ERROR("gwpt code does not support idle processes! Decrease ncpus or increase nband or use eph_np_pqbks input var.")
 end if
 gwpt%my_bsum_start = gwpt%bsum_start + gwpt%my_bsum_start - 1
 gwpt%my_bsum_stop = gwpt%bsum_start + gwpt%my_bsum_stop - 1
 ABI_MALLOC(gwpt%nbsum_rank, (gwpt%bsum_comm%nproc, 3))
 ii = gwpt%my_bsum_stop - gwpt%my_bsum_start + 1
 call xmpi_allgather(ii, gwpt%nbsum_rank(:,1), gwpt%bsum_comm%value, ierr)
 ii = gwpt%my_bsum_start
 call xmpi_allgather(ii, gwpt%nbsum_rank(:,2), gwpt%bsum_comm%value, ierr)

 call wrtout(std_out, sjoin(" Global bands for self-energy sum, bsum_start: ", itoa(gwpt%bsum_start), &
   " bsum_bstop:", itoa(gwpt%bsum_stop)))
 call wrtout(std_out, sjoin(" Allocating and treating bands from my_bsum_start: ", itoa(gwpt%my_bsum_start), &
   " up to my_bsum_stop:", itoa(gwpt%my_bsum_stop)))

 ! Distribute DFPT potentials (IBZ q-points) inside qpt_comm.
 ! Note that we distribute IBZ instead of the full BZ or the IBZ_k inside the loop over ikcalc.
 ! This means that the load won't be equally distributed but memory will scale with qpt_comm%nproc.
 ! To reduce load imbalance, we sort the qibz points by norm and use cyclic distribution inside qpt_comm
 ABI_ICALLOC(gwpt%itreat_qibz, (gwpt%nqibz))
 call sort_rpts(gwpt%nqibz, gwpt%qibz, cryst%gmet, iperm)
 do ii=1,gwpt%nqibz
   iq_ibz = iperm(ii)
   if (mod(ii, gwpt%qpt_comm%nproc) == gwpt%qpt_comm%me) gwpt%itreat_qibz(iq_ibz) = 1
 end do
 ABI_FREE(iperm)

 call wrtout(std_out, sjoin("P Number of q-points in the IBZ treated by this proc: " ,itoa(count(gwpt%itreat_qibz == 1))))

 ! ================================================================
 ! Allocate arrays used to store final results and set them to zero
 ! ================================================================
 call cwtime_report(" MPI setup", cpu, wall, gflops)

 bstart = gwpt%bsum_start

 !if (my_rank == master) then
 !  msg = "Gaps, band edges and relative position wrt Fermi level"
 !  call gaps%print(unit=std_out, kTmesh=gwpt%ktmesh, header=msg)
 !  call gaps%print(unit=ab_out, kTmesh=gwpt%ktmesh, header=msg)
 !end if
 call gaps%free()

 call cwtime_report(" gwpt_new: all", cpu_all, wall_all, gflops_all)

end function gwpt_new
!!***

subroutine gwpt_free(gwpt)

!Arguments ------------------------------------
 class(gwpt_t),intent(inout) :: gwpt

! *************************************************************************

 ! integer
 ABI_SFREE(gwpt%bstart_ks)
 ABI_SFREE(gwpt%bstop_ks)
 ABI_SFREE(gwpt%nbcalc_ks)
 ABI_SFREE(gwpt%kcalc2ibz)
 ABI_SFREE(gwpt%my_ikcalc)
 ABI_SFREE(gwpt%my_spins)
 ABI_SFREE(gwpt%itreat_qibz)
 ABI_SFREE(gwpt%my_pinfo)
 ABI_SFREE(gwpt%pert_table)
 ABI_SFREE(gwpt%ind_qbz2ibz)
 ABI_SFREE(gwpt%nbsum_rank)

 ! real
 ABI_SFREE(gwpt%kcalc)

 ABI_SFREE(gwpt%qbz)
 ABI_SFREE(gwpt%qibz)
 ABI_SFREE(gwpt%wtq)

 ! complex

 ! datatypes.
 if (allocated(gwpt%degtab)) then
   call degtab_array_free(gwpt%degtab)
   ABI_FREE(gwpt%degtab)
 end if

 ! Deallocate MPI communicators
 call gwpt%pert_comm%free()
 call gwpt%qpt_comm%free()
 call gwpt%bsum_comm%free()
 !call gwpt%qb_comm%free()
 call gwpt%kcalc_comm%free()
 call gwpt%spin_comm%free()
 !call gwpt%pqb_comm%free()
 call gwpt%ncwrite_comm%free()

 ! Close netcdf file.
 if (gwpt%ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(gwpt%ncid))
 end if

end subroutine gwpt_free
!!***

!!****f* m_gwpt/gwpt_print
!! NAME
!!  gwpt_print
!!
!! FUNCTION
!!  Print info on gwpt datatype.
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  unt=Fortran unit number
!!
!! SOURCE

subroutine gwpt_print(gwpt, dtset, unt)

!Arguments ------------------------------------
 integer,intent(in) :: unt
 type(dataset_type),intent(in) :: dtset
 class(gwpt_t),intent(in) :: gwpt

!Local variables-------------------------------
 integer :: ikc, is
 !character(len=5000) :: msg

! *************************************************************************

 if (unt == dev_null) return

 ! Write dimensions
 write(unt,"(/,a)")sjoin(" Number of bands in gwpt sum:", itoa(gwpt%nbsum))
 write(unt,"(a)")sjoin(" From bsum_start:", itoa(gwpt%bsum_start), "to bsum_stop:", itoa(gwpt%bsum_stop))
 if (dtset%eph_stern /= 0) then
   write(unt, "(a)")" Treating high-energy bands with Sternheimer and static gwpt-energy."
   write(unt, "(a, es16.6, a, i0)")" Tolwfr:", dtset%tolwfr, ", nline: ", dtset%nline
 end if
 write(unt,"(a)")sjoin(" Symsigma: ",itoa(gwpt%symsigma), "Timrev:", itoa(gwpt%timrev))
 write(unt,"(a)")sjoin(" Imaginary shift in the denominator (zcut): ", ftoa(aimag(gwpt%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
 write(unt,"(a)")sjoin(" Ab-initio q-mesh from DDB file:", ltoa(dtset%ddb_ngqpt))
 write(unt,"(a)")sjoin(" Q-mesh used for gwpt-energy integration [ngqpt]:", ltoa(gwpt%ngqpt))
 write(unt,"(a)")sjoin(" Number of q-points in the IBZ:", itoa(gwpt%nqibz))
 write(unt,"(a)")sjoin(" asr:", itoa(dtset%asr), "chneut:", itoa(dtset%chneut))
 write(unt,"(a)")sjoin(" dipdip:", itoa(dtset%dipdip), "symdynmat:", itoa(dtset%symdynmat))

 write(unt,"(a, i0)")" Number of k-points for gwpt-energy corrections: ", gwpt%nkcalc
 if (any(abs(dtset%sigma_erange) /= zero)) then
   write(unt, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
 end if
 write(unt,"(a)")" List of k-points for gwpt-energy corrections:"
 do ikc=1,gwpt%nkcalc
   if (ikc > 10) then
     write(unt, "(2a)")" nkcalc > 10. Stop printing more k-point information.",ch10
     exit
   end if
   do is=1,gwpt%nsppol
     if (gwpt%nsppol == 2) write(unt,"(a,i1,a)")" For spin: ",is, ", ikcalc, spin, kpt, bstart, bstop"
     write(unt, "(2(i4,2x),a,2(i4,1x))") &
       ikc, is, trim(ktoa(gwpt%kcalc(:,ikc))), gwpt%bstart_ks(ikc,is), gwpt%bstart_ks(ikc,is) + gwpt%nbcalc_ks(ikc,is) - 1
     end do
 end do

 write(unt, "(/,a)")" === MPI parallelism ==="
 write(unt, "(2(a,i0))")"P Allocating and summing bands from my_bsum_start: ", gwpt%my_bsum_start, &
     " up to my_bsum_stop: ", gwpt%my_bsum_stop
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", gwpt%pert_comm%nproc
 write(unt, "(a,i0)")"P Number of perturbations treated by this CPU: ", gwpt%my_npert
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", gwpt%qpt_comm%nproc
 write(unt, "(2(a,i0))")"P Number of q-points in the IBZ treated by this proc: " , &
     count(gwpt%itreat_qibz == 1), " of ", gwpt%nqibz
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over bands: ", gwpt%bsum_comm%nproc
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over spins: ", gwpt%spin_comm%nproc
 write(unt, "(a,i0)")"P Number of CPUs for parallelism over k-points: ", gwpt%kcalc_comm%nproc
 write(unt, "(2(a,i0),/)")"P Number of k-point in gwpt_nk treated by this proc: ", gwpt%my_nkcalc, " of ", gwpt%nkcalc

end subroutine gwpt_print
!!***

end module m_gwpt
!!***

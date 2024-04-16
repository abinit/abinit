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
 use m_dfpt_cgwf,      only : dfpt_cgwf
 use m_phonons,        only : phstore_t, phstore_new
 use m_io_screening,   only : hscr_t, get_hscr_qmesh_gsph
 use m_vcoul,          only : vcoul_t

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

  type(xcomm_t) :: qb_comm
   ! MPI communicator used to distribute (band_sum, q-points)

  type(xcomm_t) :: qpt_comm
   ! MPI communicator for q-points

  type(xcomm_t) :: bsum_comm
   ! MPI communicator for bands in self-energy sum

  type(xcomm_t) :: kcalc_comm
   ! MPI communicator for parallelism over k-points (high-level)

  type(xcomm_t) :: spin_comm
   ! MPI communicator for parallelism over spins (high-level)

  type(xcomm_t) :: pqb_comm
    ! MPI communicator for the (perturbation, band_sum, qpoint_sum)

  type(xcomm_t) :: ncwrite_comm
   ! MPI communicator for parallel netcdf IO used to write results for the different k-points/spins

  integer :: coords_pqbks(5)
   ! Cartesian coordinates of this processor in the Cartesian grid.

  integer :: nqbz
   ! Number of q-points in the (dense) BZ for gwpt integration

  integer :: nqibz
   ! Number of q-points in the (dense) IBZ for gwpt integration

  integer :: nqibz_k
   ! Number of q-points in the IBZ(k). Depends on ikcalc.

  integer :: my_nqibz_k
   ! Number of q-points in the IBZ(k) treated by this MPI proc. Depends on ikcalc.
   ! Differs from nqibz_k only if imag with tetra because in this case we can introduce a cutoff on the weights

  integer :: lgk_nsym
   ! Number of symmetries in the little group of k. Depends on ikcalc.

  integer :: ncid = nctk_noid
   ! Netcdf file handle used to save results.

  integer :: mpw
   ! Maximum number of PWs for all possible k+q

  complex(dpc) :: ieta
   ! Used to shift the poles in the complex plane (Ha units)
   ! Corresponds to `i eta` term in equations.

  logical :: use_ftinterp = .False.
   ! whether DFPT potentials should be read from the DVDB or Fourier-interpolated on the fly.

  integer :: gmax(3)

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

  integer,allocatable :: myq2ibz_k(:)
   ! myq2ibz_k(my_nqibz_k)
   ! Mapping my q-point index --> index in nqibz_k arrays (IBZ_k)
   ! Differs from nqibz_k only if imag with tetra because in this case we can introduce a cutoff.

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

  integer,allocatable:: indkk_kq(:, :)
   ! (6, %nqibz_k))
   ! Mapping k+q --> initial IBZ. Depends on ikcalc.
   ! These table used the conventions for the symmetrization of the wavefunctions expected by cgtk_rotate.
   ! In this case listkk has been called with symrel and use_symrec=False

  integer,allocatable :: ind_q2dvdb_k(:,:)
   ! (6, %nqibz_k))
   ! Mapping qibz_k --> IBZ found in DVDB file.
   ! Used when DFPT potentials are read from DVDB file so that we know how to access/symmetrize v1scf
   ! Depends on ikcalc.

  integer,allocatable :: ind_ibzk2ibz(:,:)
   ! (6, %nqibz_k))
   ! Mapping qibz_k --> IBZ defined by eph_ngqpt_fine.
   ! Depends on ikcalc.

  integer,allocatable :: qibz2dvdb(:)
   ! (%nqibz))
   ! Mapping dvdb%ibz --> %ibz

  integer, allocatable :: lgk_sym2glob(:, :)
   ! lgk_sym2glob(2, lgk_nsym)
   ! Mapping isym_lg --> [isym, itime]
   ! where isym is the index of the operation in the global array **crystal%symrec**
   ! and itim is 2 if time-reversal T must be included else 1. Depends on ikcalc

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

  real(dp),allocatable :: qibz_k(:,:)
   ! qibz(3, nqibz_k)
   ! Reduced coordinates of the q-points in the IBZ(k). Depends on ikcalc.

  real(dp),allocatable :: wtq_k(:)
   ! wtq(nqibz_k)
   ! Weights of the q-points in the IBZ(k) (normalized to one). Depends on ikcalc.

  real(dp),allocatable :: vcar_calc(:,:,:,:)
   ! (3, max_nbcalc, nkcalc, nsppol))
   ! Diagonal elements of velocity operator in cartesian coordinates for all states in gwpt_nk.

  integer, allocatable :: qp_done(:,:)
   ! qp_done(kcalc, spin)
   ! Keep track of the QP states already computed for restart of the calculation

  type(degtab_t),allocatable :: degtab(:,:)
   ! (nkcalc, nsppol)
   ! Table used to average QP results in the degenerate subspace if symsigma == 1

  contains

    procedure :: setup_kcalc => gwpt_setup_kcalc
     ! Return tables used to perform the sum over q-points for given k-point.

    procedure :: print => gwpt_print
     ! Print results to main output file.

    procedure :: free => gwpt_free
      ! Free gwpt object

 end type gwpt_t
!!***

 public :: gwpt_run        ! Main entry point to compute GWpt matrix elements
 !private :: gwpt_new   ! Creation method (allocates memory, initialize data from input vars).

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
 type(dvdb_t),intent(inout) :: dvdb
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
 integer :: band_me, nband_me
 integer :: my_rank,nsppol,nkpt,iq_ibz,iq_ibz_k,my_npert ! iq_ibz_frohl,iq_bz_frohl,
 integer :: cplex,db_iqpt,natom,natom3,ipc,nspinor,nprocs, qptopt ! = 1
 integer :: ibsum_kq, ib_k, u1c_ib_k, band_ks, u1_band !, ii, jj, iw !ib_kq, ibsum,
 !integer :: u1_master, ip
 integer :: mcgq, mgscq !, ig, ispinor !, ifft !nband_kq,
 integer :: idir,ipert !,ip1,ip2,idir1,ipert1,idir2,ipert2
 integer :: ik_ibz,ikq_ibz,isym_k,isym_kq,trev_k,trev_kq, isym_q, trev_q
 integer :: my_spin, spin, istwf_k, istwf_kq, istwf_kqirr, npw_k, npw_kq, npw_kqirr
 integer :: mpw,ierr,imyq,ignore_kq, ignore_ibsum_kq ! band,
 integer :: n1,n2,n3,n4,n5,n6,nspden,nu, mqmem, jm, jn
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,cnt,imyp !, restart
 integer :: tot_nlines_done, nlines_done, nline_in, grad_berry_size_mpw1, enough_stern
 integer :: nbcalc_ks,nbsum,bsum_start, bsum_stop, bstart_ks,my_ikcalc,ikcalc,bstart,bstop, sendcount !iatom,
 integer :: comm_rpt, nqlwl, ebands_timrev ! osc_npw,
 integer :: ffnlk_request, ffnl1_request, nelem, cgq_request
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all,cpu_ks,wall_ks,gflops_ks
 real(dp) :: cpu_setk, wall_setk, gflops_setk, cpu_qloop, wall_qloop, gflops_qloop
 real(dp) :: ecut,eshift,weight_q,rfact,ediff,q0rad, out_resid
 real(dp) :: bz_vol
 logical :: isirr_k, isirr_kq, gen_eigenpb, q_is_gamma, isirr_q, use_ifc_fourq, use_u1c_cache, intra_band, same_band
 type(krank_t) :: krank
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(gwpt_t) :: gwpt !, gwpt_restart
 type(ddkop_t) :: ddkop
 type(rf2_t) :: rf2
 type(crystal_t) :: pot_cryst, den_cryst
 type(hdr_type) :: pot_hdr, den_hdr
 type(phstore_t) :: phstore
 type(u1cache_t) :: u1c
 type(kmesh_t) :: qmesh, kmesh
 type(gsphere_t) :: gsph_c
 type(hscr_t) :: hscr
 type(vcoul_t) :: vcp
 type(screen_t) :: W
 type(screen_info_t) :: W_info
 character(len=fnlen) :: w_fname
 character(len=5000) :: msg
!arrays
 integer :: g0_k(3),g0_kq(3), units(2), work_ngfft(18), gmax(3), indkk_kq(6,1)
 integer,allocatable :: bands_treated_now(:)
 integer(i1b),allocatable :: itreatq_dvdb(:)
 integer,allocatable :: kg_k(:,:),kg_kq(:,:),nband(:,:), qselect(:), wfd_istwfk(:)
 integer,allocatable :: gbound_kq(:,:), osc_gbound_q(:,:), rank_band(:), root_bcalc(:) ! osc_indpw(:), osc_gvecq(:,:),
 integer,allocatable :: ibzspin_2ikcalc(:,:)
 integer, allocatable :: displs(:), recvcounts(:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3),qpt(3),phfrq(3*cryst%natom), dotri(2),qq_ibz(3) !qpt_cart(3),
 real(dp) :: tsec(2) ! vk(3), vkq(3),
 real(dp) :: vec_natom3(2, 3*cryst%natom)
 real(dp) :: wqnu,gkq2,eig0nk,eig0mkq !eig0mk,
 !real(dp) :: gdw2, gdw2_stern, rtmp
 real(dp),allocatable,target :: cgq(:,:,:)
 real(dp),allocatable :: qlwl(:,:) ! doccde(:),eigen(:),occfact(:),
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: grad_berry(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:)
 real(dp),allocatable :: ffnlk(:,:,:,:),ffnl1(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:),gkq0_atm(:,:,:,:)
 real(dp),allocatable :: gscq(:,:,:), out_eig1_k(:), cg1s_kq(:,:,:,:), h1kets_kq_allperts(:,:,:,:)
 real(dp),allocatable :: dcwavef(:, :), gh1c_n(:, :), ghc(:,:), gsc(:,:), stern_ppb(:,:,:,:), stern_dw(:,:,:,:)
 logical,allocatable :: ihave_ikibz_spin(:,:), bks_mask(:,:,:),keep_ur(:,:,:)
 real(dp),allocatable :: bra_kq(:,:),kets_k(:,:,:),h1kets_kq(:,:,:,:),cgwork(:,:)
 real(dp),allocatable :: ph1d(:,:),vlocal(:,:,:,:),vlocal1(:,:,:,:,:)
 real(dp),allocatable :: ylm_kq(:,:),ylm_k(:,:),ylmgr_kq(:,:,:)
 real(dp),allocatable :: vtrial(:,:),gvnlx1(:,:),gvnlxc(:,:),work(:,:,:,:), vcar_ibz(:,:,:,:), rhor(:,:)
 real(dp),allocatable :: gs1c(:,:)
 real(dp),allocatable :: gkq_allgather(:,:,:)
 !real(dp),allocatable :: phfreqs_qibz(:,:), pheigvec_qibz(:,:,:,:), eigvec_qpt(:,:,:)
 real(dp) :: ylmgr_dum(1,1,1)
 !logical,allocatable :: osc_mask(:)
 !real(dp),allocatable :: gkq2_lr(:,:,:)
 !complex(dpc),allocatable :: osc_ks(:,:)
 complex(gwpc),allocatable :: ur_k(:), ur_kq(:) !, work_ur(:), workq_ug(:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:), cwaveprj(:,:)
 type(pawrhoij_type),allocatable :: pot_pawrhoij(:), den_pawrhoij(:)
#if defined HAVE_MPI && !defined HAVE_MPI2_INPLACE
 integer :: me
 real(dp),allocatable :: cgq_buf(:)
 real(dp),pointer :: cgq_ptr(:)
#endif

!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 units = [std_out, ab_out]

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
 gwpt = gwpt_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm)

 ! Check if a previous netcdf file is present and restart the calculation
 ! Here we try to read an existing SIGEPH file if eph_restart == 1.
 ! and we compare the variables with the state of the code (i.e. new gwpt generated in gwpt_new)
 !restart = 0; ierr = 1; sigeph_filepath = strcat(dtfil%filnam_ds(4), "_SIGEPH.nc")
 !if (my_rank == master .and. dtset%eph_restart == 1) then
 !  gwpt_restart = gwpt_read(sigeph_filepath, dtset, xmpi_comm_self, msg, ierr)
 !end if

 !if (my_rank == master .and. dtset%eph_restart == 1) then
 !  if (ierr == 0) then
 !    if (any(gwpt_restart%qp_done /= 1)) then
 !      !call gwpt%compare(gwpt_restart)
 !      ! Get list of QP states that have been computed.
 !      gwpt%qp_done = gwpt_restart%qp_done
 !      restart = 1
 !      call wrtout(units, "- Restarting from previous SIGEPH.nc file")
 !      call wrtout(units, sjoin("- Number of k-points completed:", itoa(count(gwpt%qp_done == 1)), "/", itoa(gwpt%nkcalc)))
 !    else
 !      restart = 0; gwpt%qp_done = 0
 !      msg = sjoin("Found SIGEPH.nc file with all QP entries already computed.", ch10, &
 !                  "Will overwrite:", sigeph_filepath, ch10, &
 !                  "Keeping backup copy in:", strcat(sigeph_filepath, ".bkp"))
 !      call wrtout(ab_out, sjoin("WARNING: ", msg))
 !      ABI_WARNING(msg)
 !      ! Keep backup copy
 !      ABI_CHECK(clib_rename(sigeph_filepath, strcat(sigeph_filepath, ".bkp")) == 0, "Failed to rename SIGPEPH file.")
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
 !    NCF_CHECK(nctk_open_modify(gwpt%ncid, sigeph_filepath, gwpt%ncwrite_comm%value))
 !    NCF_CHECK(nctk_set_datamode(gwpt%ncid))
 !  end if
 !end if

 if (my_rank == master) then
   call gwpt%print(dtset, ab_out)
   call gwpt%print(dtset, std_out)
 end if
 my_npert = gwpt%my_npert


 ! This is the maximum number of PWs for all possible k+q treated.
 mpw = gwpt%mpw; gmax = gwpt%gmax

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
 ABI_MALLOC(ibzspin_2ikcalc, (nkpt, nsppol))
 ibzspin_2ikcalc = -1

 ! Each node needs the wavefunctions for gwpt_{nk}
 ! TODO: kcalc should depend on the spin!
 do spin=1,gwpt%nsppol
   do ikcalc=1,gwpt%nkcalc
     ik_ibz = gwpt%kcalc2ibz(ikcalc, 1)
     bstart = gwpt%bstart_ks(ikcalc, spin)
     bstop = bstart + gwpt%nbcalc_ks(ikcalc, spin) - 1
     bks_mask(bstart:bstop, ik_ibz, spin) = .True.
     ibzspin_2ikcalc(ik_ibz, spin) = ikcalc
   end do
 end do

 bks_mask(gwpt%my_bsum_start:gwpt%my_bsum_stop,:,:) = .True.

 !if (dtset%userie == 124) then
 !  ! Uncomment this line to have all states on each MPI rank.
 !  bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
 !end if

 ! This table is needed when computing the imaginary part:
 ! k+q states outside the energy window are not read hence their contribution won't be included.
 ! Error is small provided calculation is close to convergence.
 ! To reduce the error one should increase the value of phwinfact
 ABI_MALLOC(ihave_ikibz_spin, (nkpt, nsppol))
 ihave_ikibz_spin = .False.
 do spin=1,gwpt%nsppol
   do ik_ibz=1,ebands%nkpt
     if (any(bks_mask(:, ik_ibz, spin))) ihave_ikibz_spin(ik_ibz, spin) = .True.
   end do
 end do

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

 ! if PAW, one has to solve a generalized eigenproblem
 ! Be careful here because I will need sij_opt == -1
 usecprj = 0; gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1

 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))
 ABI_MALLOC(cwaveprj, (natom, nspinor*usecprj))
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(gbound_kq, (2*wfd%mgfft+8, 2))
 ABI_MALLOC(osc_gbound_q, (2*wfd%mgfft+8, 2))

 ! ============================
 ! Compute vnk matrix elements
 ! ============================
 ABI_MALLOC(cgwork, (2, mpw*wfd%nspinor))
 ABI_CALLOC(gwpt%vcar_calc, (3, gwpt%max_nbcalc, gwpt%nkcalc, nsppol))

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
       call wfd%copy_cg(band_ks, ik_ibz, spin, cgwork)
       eig0nk = ebands%eig(band_ks, ik_ibz, spin)
       gwpt%vcar_calc(:, ib_k, ikcalc, spin) = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cgwork, cwaveprj0)
     end do

   end do
 end do
 call xmpi_sum(gwpt%vcar_calc, comm, ierr)

 ! Write v_nk to disk.
 !if (my_rank == master) then
 !  NCF_CHECK(nf90_put_var(gwpt%ncid, nctk_idname(gwpt%ncid, "vcar_calc"), gwpt%vcar_calc))
 !end if

 ABI_FREE(cgwork)
 call ddkop%free()
 call cwtime_report(" Velocities", cpu_ks, wall_ks, gflops_ks)

 ! Precompute phonon frequencies and eigenvectors in the IBZ.
 ! These quantities are then used to symmetrize quantities for q in the IBZ(k) in order
 ! to reduce the number of calls to ifc%fourq (expensive if dipdip == 1)

 use_ifc_fourq = .False. !use_ifc_fourq = .True. !use_ifc_fourq = dtset%userib == 123
 phstore = phstore_new(cryst, ifc, gwpt%nqibz, gwpt%qibz, use_ifc_fourq, gwpt%pert_comm%value)
 call cwtime_report(" phonons in the IBZ", cpu_ks, wall_ks, gflops_ks)

 ! Radius of sphere with volume equivalent to the micro zone.
 q0rad = two_pi * (three / (four_pi * cryst%ucvol * gwpt%nqbz)) ** third
 bz_vol = two_pi**3 / cryst%ucvol

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1   ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2      ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0 ! gvnlx1 is output

 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 ! 2) Perform the setup needed for the non-local factors:
 !
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 call init_hamiltonian(gs_hamkq, psps, pawtab, nspinor, nsppol, nspden, natom,&
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg,&
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab,&
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 ! Allocate work space arrays.
 ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, gs_hamkq%nvloc))

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

 if (gwpt%pert_comm%nproc > 1) then
   !  Activate parallelism over perturbations
   call dvdb%set_pert_distrib(gwpt%my_npert, natom3, gwpt%my_pinfo, gwpt%pert_table, gwpt%pert_comm%value)
 end if

 ! Find correspondence IBZ --> set of q-points in DVDB.
 ! Activate FT interpolation automatically if required q-points in the IBZ are not found in the DVDB.
 gwpt%use_ftinterp = .False.
 ABI_MALLOC(gwpt%qibz2dvdb, (gwpt%nqibz))
 if (dvdb%find_qpts(gwpt%nqibz, gwpt%qibz, gwpt%qibz2dvdb, comm) /= 0) then
   call wrtout(units, " Cannot find eph_ngqpt_fine q-points in DVDB --> Activating Fourier interpolation.")
   gwpt%use_ftinterp = .True.
 else
   call wrtout(units, " DVDB file contains all q-points in the IBZ --> Reading DFPT potentials from file.")
   gwpt%use_ftinterp = .False.
 end if

 if (gwpt%use_ftinterp) then
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

   ! Build q-cache in the *dense* IBZ using the global mask qselect and itreat_qibz.
   ABI_MALLOC(qselect, (gwpt%nqibz))
   qselect = 1
   call dvdb%ftqcache_build(nfftf, ngfftf, gwpt%nqibz, gwpt%qibz, dtset%dvdb_qcache_mb, qselect, gwpt%itreat_qibz, comm)

 else
   ABI_MALLOC(qselect, (dvdb%nqpt))
   qselect = 1
 end if

 call dvdb%print(prtvol=dtset%prtvol)

 if (.not. gwpt%use_ftinterp) then
   ! Need to translate itreat_qibz into itreatq_dvdb.
   ABI_ICALLOC(itreatq_dvdb, (dvdb%nqpt))
   do iq_ibz=1,gwpt%nqibz
     if (gwpt%itreat_qibz(iq_ibz) == 0) cycle
     db_iqpt = gwpt%qibz2dvdb(iq_ibz)
     ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(gwpt%qibz(:, iq_ibz)), "in the DVDB file."))
     itreatq_dvdb(db_iqpt) = 1
   end do
   call dvdb%qcache_read(nfftf, ngfftf, dtset%dvdb_qcache_mb, qselect, itreatq_dvdb, comm)
   ABI_FREE(itreatq_dvdb)
 end if

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

 if (w_fname /= ABI_NOFILE) then
   ! Read gsphere and qmesh from SCR file.
   call get_hscr_qmesh_gsph(w_fname, dtset, cryst, hscr, qmesh, gsph_c, qlwl, comm)
   call hscr%free()
   nqlwl = size(qlwl, dim=2)
   w_info%use_mdf = MDL_NONE
 else
   ! Init Qmesh from the K-mesh reported in the WFK file.
   call find_qmesh(qmesh, cryst, kmesh)
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
 call vcp%init(gsph_c, cryst, qmesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_c%ng, &
               nqlwl, qlwl, comm)

 ! Init Wc.
 ! Incore or out-of-core solution?
 mqmem = 0; if (dtset%gwmem /10 == 1) mqmem = qmesh%nibz
 w_info%invalid_freq = dtset%gw_invalid_freq
 w_info%mat_type = MAT_INV_EPSILON
 call W%init(w_info, cryst, qmesh, gsph_c, vcp, w_fname, mqmem, dtset%npweps, &
             dtset%iomode, ngfftf, nfftf, nsppol, nspden, rhor, dtset%prtvol, comm)

 ABI_FREE(rhor)
 ABI_FREE(qlwl)

 ! Allocate g-vectors for k and k+q
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 ! Spherical Harmonics for useylm == 1.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

 ebands_timrev = kpts_timrev_from_kptopt(ebands%kptopt)

 ABI_MALLOC(ur_k, (wfd%nfft*nspinor))
 ABI_MALLOC(ur_kq, (wfd%nfft*nspinor))
 !ABI_MALLOC(work_ur, (wfd%nfft*nspinor))
 !ABI_MALLOC(gkq2_lr, (gwpt%eph_doublegrid%ndiv, nbcalc_ks, gwpt%my_npert))

 ! Build krank object to find k-points
 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 ! Find correspondence IBZ_k --> IBZ
 !ABI_MALLOC(iqk2dvdb, (6, self%nqibz_k))

 ! Assume qptopt == kptopt unless value is specified in input
 !qptopt = ebands%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
 !qtimrev = kpts_timrev_from_kptopt(qptopt)
 !qptrlatt = 0; qptrlatt(1,1) = self%ngqpt(1); qptrlatt(2,2) = self%ngqpt(2); qptrlatt(3,3) = self%ngqpt(3)
 !qrank = krank_from_kptrlatt(self%nqibz, self%qibz, qptrlatt, compute_invrank=.False.)

 !if (kpts_map("symrec", qtimrev, cryst, qrank, self%nqibz_k, self%qibz_k, iqk2dvdb) /= 0) then
 !  write(msg, '(3a)' )&
 !    "At least one of the q points in the IBZ_k could not be generated from one in the IBZ.", ch10,&
 !    "Action: check your DVDB file and use eph_task to interpolate the potentials on a denser q-mesh."
 !  ABI_ERROR(msg)
 !end if
 !call qrank%free()

 ! Loop over (spin, qbiz, atom_pert) in Sigma^{spin}_{qibz, ipert)
 do my_spin=1,gwpt%my_nspins
   spin = gwpt%my_spins(my_spin)
   do iq_ibz=1,gwpt%nqibz
     !iq_ibz = gwpt%my_iqbz(my_iq)
     qq_ibz = gwpt%qibz(:, iq_ibz)
     q_is_gamma = sum(qq_ibz**2) < tol14
     ! Compute phonons for this qq_ibz
     call ifc%fourq(cryst, qq_ibz, phfrq, displ_cart, out_displ_red=displ_red)

     ! ====================================
     ! Get DFPT potentials for this q-point
     ! ====================================
     cplex = 2 ! FIXME
     ABI_MALLOC(v1scf, (cplex, nfftf, nspden, my_npert))
     !if (gwpt%use_ftinterp) then
     !  ! Use Fourier interpolation to get DFPT potentials for this qpt (hopefully in cache).
     !  db_iqpt = gwpt%ind_ibzk2ibz(1, iq_ibz_k)
     !  call dvdb%get_ftqbz(cryst, qpt, qq_ibz, gwpt%ind_ibzk2ibz(:, iq_ibz_k), cplex, nfftf, ngfftf, v1scf, &
     !                      gwpt%pert_comm%value)
     !else
     !  ! Read and reconstruct the dvscf potentials for qpt and my_npert perturbations.
     !  ! This call allocates v1scf(cplex, nfftf, nspden, my_npert))
     !  db_iqpt = gwpt%ind_q2dvdb_k(1, iq_ibz_k)
     !  ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file."))
     !  call dvdb%readsym_qbz(cryst, qpt, gwpt%ind_q2dvdb_k(:,iq_ibz_k), cplex, nfftf, ngfftf, v1scf, gwpt%pert_comm%value)
     !end if

     ! Allocate vlocal1 with correct cplex. Note nvloc
     ABI_MALLOC_OR_DIE(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert), ierr)

     do my_ikcalc=1,gwpt%my_nkcalc
       ikcalc = gwpt%my_ikcalc(my_ikcalc)

       ! Symmetry indices for kk.
       kk = gwpt%kcalc(:, ikcalc)
       ik_ibz = gwpt%kcalc2ibz(ikcalc, 1); isym_k = gwpt%kcalc2ibz(ikcalc, 2)
       trev_k = gwpt%kcalc2ibz(ikcalc, 6); g0_k = gwpt%kcalc2ibz(ikcalc, 3:5)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       ABI_CHECK(isirr_k, "For the time being the k-point in gwpt_{nk} must be in the IBZ")
       kk_ibz = ebands%kptns(:,ik_ibz)
       ! Get npw_k, kg_k for k
       call wfd%get_gvec_kq(cryst%gmet, ecut, kk, ik_ibz, isirr_k, istwf_k, npw_k, kg_k)

       ! Compute k+G vectors
       nkpg = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg_k, (npw_k, nkpg))
       if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk, nkpg, npw_k)

       ! Compute nonlocal form factors ffnlk at (k+G)
       ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))

       call mkffnl_objs(cryst, psps, 1, ffnlk, ider0, idir0, kg_k, kpg_k, kk, nkpg, npw_k, ylm_k, ylmgr_dum, &
                        comm=gwpt%pert_comm%value, request=ffnlk_request)

       ! Find k + q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered as:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qq_ibz
       if (kpts_map("symrel", ebands_timrev, cryst, krank, 1, kq, indkk_kq) /= 0) then
         write(msg, '(2a)' ) &
          "The WFK file cannot be used to compute phonon linewidths.",ch10
          !"At least one of the k-points on the FS could not be generated from a symmetrical one.", ch10, &
          !"q-mesh: ", trim(ltoa(gamma_ngqpt)), ", k-mesh (from kptrlatt) ", trim(ltoa(get_diag(ebands%kptrlatt))), &
          !'Action: check your WFK file and the (k, q) point input variables.'
          ABI_ERROR(msg)
       end if
       ikq_ibz = indkk_kq(1, 1); isym_kq = indkk_kq(2, 1)
       trev_kq = indkk_kq(6, 1); g0_kq = indkk_kq(3:5, 1)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)

       ! Get npw_kq, kg_kq for k+q.
       call wfd%get_gvec_kq(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, istwf_kq, npw_kq, kg_kq)

       ! Compute k+q+G vectors
       nkpg1 = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))
       call mkffnl_objs(cryst, psps, 1, ffnl1, ider0, idir0, kg_kq, kpg1_k, kq, nkpg1, npw_kq, ylm_kq, ylmgr_kq, &
                        comm=gwpt%pert_comm%value, request=ffnl1_request)

       ! Double loop over the bands indices in the e-ph matrix elements.
       do jm=1,1
         if (isirr_kq) then
           call wfd%get_ur(jm, ikq_ibz, spin, ur_kq)
         else
         end if

         do jn=1,1
           if (isirr_k) then
             ! Copy u_kq(G)
             !call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
             call wfd%get_ur(jn, ik_ibz, spin, ur_k)
           else
             !! Reconstruct u_kq(G) from the IBZ image.
             !call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)
             !call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1, &
             !                 npw_kqirr, wfd%kdata(ikq_ibz)%kg_k, &
             !                 npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
           end if
         end do ! jn
       end do ! jm

       do imyp=1,gwpt%my_npert
         idir = gwpt%my_pinfo(1, imyp); ipert = gwpt%my_pinfo(2, imyp); ipc = gwpt%my_pinfo(3, imyp)

         ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
         ! Each CPU prepares its own potentials.
         call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc, &
           pawfgr, mpi_enreg, vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))

         ! Continue to initialize the Hamiltonian (call it here to support dfpt_cgwf Sternheimer).
         call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)
         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         if (ffnlk_request /= xmpi_request_null) call xmpi_wait(ffnlk_request, ierr)
         if (ffnl1_request /= xmpi_request_null) call xmpi_wait(ffnl1_request, ierr)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &  ! In
           cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
           npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
           dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
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
         !    grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
         !    optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c1, usevnl)
         !end do

         ! NSCF solution of the Sternheimer equation for all the bands included in the sum over states.
         !do iq_sum=1,nqbz
         !end do

         ABI_FREE(dkinpw)
         ABI_FREE(kinpw1)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)
       end do ! imyp

     end do ! my_ikcalc

     ABI_FREE(kpg_k)
     ABI_FREE(kpg1_k)
     ABI_FREE(ffnlk)
     ABI_FREE(ffnl1)
     ABI_FREE(vlocal1)
     ABI_FREE(v1scf)
   end do ! iq_ibz
 end do ! spin

 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(ur_k)
 ABI_FREE(ur_kq)


 call krank%free()
 call vcp%free()
 call w%free()
 call qmesh%free()
 call gsph_c%free()
 call kmesh%free()

 RETURN

 ! Loop over k-points in gwpt_nk. Loop over spin is internal as we operate on nspden components at once.
 do my_ikcalc=1,gwpt%my_nkcalc
   ikcalc = gwpt%my_ikcalc(my_ikcalc)

   ! Check if this (kpoint, spin) was already calculated
   if (all(gwpt%qp_done(ikcalc, :) == 1)) cycle
   call cwtime(cpu_ks, wall_ks, gflops_ks, "start")

   ! Find IBZ(k) for q-point integration.
   call cwtime(cpu_setk, wall_setk, gflops_setk, "start")
   ! FIXME invert spin but checks shape of the different arrays!
   call gwpt%setup_kcalc(dtset, cryst, ebands, ikcalc, dtset%prtvol, gwpt%pqb_comm%value)

   ! Symmetry indices for kk.
   kk = gwpt%kcalc(:, ikcalc)
   ik_ibz = gwpt%kcalc2ibz(ikcalc, 1); isym_k = gwpt%kcalc2ibz(ikcalc, 2)
   trev_k = gwpt%kcalc2ibz(ikcalc, 6); g0_k = gwpt%kcalc2ibz(ikcalc, 3:5)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   ABI_CHECK(isirr_k, "For the time being the k-point in gwpt_{nk} must be in the IBZ")
   kk_ibz = ebands%kptns(:,ik_ibz)
   npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)

   ! Allocate PW-arrays. Note mpw in kg_kq
   ABI_MALLOC(kg_k, (3, npw_k))
   kg_k = wfd%kdata(ik_ibz)%kg_k
   ABI_MALLOC(kg_kq, (3, mpw))

   ! Spherical Harmonics for useylm == 1.
   ABI_MALLOC(ylm_k, (mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylm_kq, (mpw, psps%mpsang**2 * psps%useylm))
   ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))

   ! Compute k+G vectors
   nkpg = 3*dtset%nloalg(3)
   ABI_MALLOC(kpg_k, (npw_k, nkpg))
   if (nkpg > 0) call mkkpg(kg_k, kpg_k, kk, nkpg, npw_k)

   ! Compute nonlocal form factors ffnlk at (k+G)
   ABI_MALLOC(ffnlk, (npw_k, 1, psps%lmnmax, psps%ntypat))

   call mkffnl_objs(cryst, psps, 1, ffnlk, ider0, idir0, kg_k, kpg_k, kk, nkpg, npw_k, ylm_k, ylmgr_dum, &
                    comm=gwpt%pert_comm%value, request=ffnlk_request)

   call cwtime_report(" Setup kcalc", cpu_setk, wall_setk, gflops_setk)

   ! TODO: Spin should be treated in a more flexible and scalable way --> kcalc and bdgw should depend on spin.
   ! Introduce other comm and cartesian dimension for spin
   do my_spin=1,gwpt%my_nspins
     spin = gwpt%my_spins(my_spin)

     ! Check if this kpoint and spin was already calculated
     if (gwpt%qp_done(ikcalc, spin) == 1) cycle

     !call timab(1900, 1, tsec)
     ! Bands in gwpt_nk to compute and number of bands in sum over states.
     bstart_ks = gwpt%bstart_ks(ikcalc, spin)
     nbcalc_ks = gwpt%nbcalc_ks(ikcalc, spin)
     bsum_start = gwpt%bsum_start; bsum_stop = gwpt%bsum_stop
     nbsum = gwpt%nbsum
     ABI_MALLOC(root_bcalc, (nbcalc_ks))

     ! Allocate eph matrix elements.
     ABI_MALLOC(gkq_atm, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkq_nu, (2, nbcalc_ks, natom3))
     ABI_MALLOC(gkq_allgather, (2, nbcalc_ks * natom3, 2))

     ! Allocate arrays for Debye-Waller
     ABI_CALLOC_OR_DIE(gkq0_atm, (2, nbcalc_ks, gwpt%my_bsum_start:gwpt%my_bsum_stop, natom3), ierr)
     if (dtset%eph_stern /= 0) then
       ABI_CALLOC(stern_dw, (2, natom3, natom3, nbcalc_ks))
       enough_stern = 0
       use_u1c_cache = merge(.True., .False., dtset%eph_stern == 1)
       tot_nlines_done = 0
     end if

     ! Load ground-state wavefunctions for which corrections are wanted (available on each node)
     ! Note: One should rotate the wavefunctions if kk is not in the IBZ (not implemented)
     ABI_MALLOC(kets_k, (2, npw_k*nspinor, nbcalc_ks))

     !if (osc_ecut /= zero) then
     !  ABI_MALLOC(ur_k, (wfd%nfft*nspinor, nbcalc_ks))
     !  ABI_MALLOC(ur_kq, (wfd%nfft*nspinor))
     !  ABI_MALLOC(work_ur, (wfd%nfft*nspinor))
     !  ABI_MALLOC(gkq2_lr, (gwpt%eph_doublegrid%ndiv, nbcalc_ks, gwpt%my_npert))
     !end if

     do ib_k=1,nbcalc_ks
       band_ks = ib_k + bstart_ks - 1
       call wfd%copy_cg(band_ks, ik_ibz, spin, kets_k(1, 1, ib_k))
       !if (osc_ecut > zero) call wfd%get_ur(band_ks, ik_ibz, spin, ur_k(1, ib_k))
     end do

     ! Distribute q-points, compute tetra weigths.
     !call gwpt_setup_qloop(gwpt, dtset, cryst, ebands, dvdb, spin, ikcalc, nfftf, ngfftf, gwpt%pqb_comm%value)
     !call timab(1900, 2, tsec)

     ! ==========================================
     ! Integration over my q-points in the IBZ(k)
     ! ==========================================
     call cwtime(cpu_qloop, wall_qloop, gflops_qloop, "start")
     ignore_kq = 0; ignore_ibsum_kq = 0

     do imyq=1,gwpt%my_nqibz_k
       call cwtime(cpu, wall, gflops, "start")
       iq_ibz_k = gwpt%myq2ibz_k(imyq)
       qpt = gwpt%qibz_k(:, iq_ibz_k)
       q_is_gamma = sum(qpt**2) < tol14

       iq_ibz = gwpt%ind_ibzk2ibz(1, iq_ibz_k)
       isym_q = gwpt%ind_ibzk2ibz(2, iq_ibz_k)
       trev_q = gwpt%ind_ibzk2ibz(6, iq_ibz_k)
       ! Don't test if umklapp == 0 because we use the periodic gauge: phfreq(q+G) = phfreq(q) and eigvec(q) = eigvec(q+G)
       isirr_q = (isym_q == 1 .and. trev_q == 0)
       !qq_ibz = gwpt%qibz(:, iq_ibz)

       ! Find k + q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered as:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qpt
       ikq_ibz = gwpt%indkk_kq(1, iq_ibz_k); isym_kq = gwpt%indkk_kq(2, iq_ibz_k)
       trev_kq = gwpt%indkk_kq(6, iq_ibz_k); g0_kq = gwpt%indkk_kq(3:5, iq_ibz_k)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)
       !nband_kq = ebands%nband(ikq_ibz + (spin-1) * ebands%nkpt)

       ! This can happen if we have loaded the wavefunctions inside the energy range.
       !if (gwpt%imag_only .and. .not. ihave_ikibz_spin(ikq_ibz, spin)) then
       !  ignore_kq = ignore_kq + 1; cycle
       !end if

       ! ====================================
       ! Get DFPT potentials for this q-point
       ! ====================================
       if (gwpt%use_ftinterp) then
         ! Use Fourier interpolation to get DFPT potentials for this qpt (hopefully in cache).
         db_iqpt = gwpt%ind_ibzk2ibz(1, iq_ibz_k)
         qq_ibz = gwpt%qibz(:, db_iqpt)
         call dvdb%get_ftqbz(cryst, qpt, qq_ibz, gwpt%ind_ibzk2ibz(:, iq_ibz_k), cplex, nfftf, ngfftf, v1scf, &
                             gwpt%pert_comm%value)
       else
         ! Read and reconstruct the dvscf potentials for qpt and my_npert perturbations.
         ! This call allocates v1scf(cplex, nfftf, nspden, my_npert))
         db_iqpt = gwpt%ind_q2dvdb_k(1, iq_ibz_k)
         ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file."))
         call dvdb%readsym_qbz(cryst, qpt, gwpt%ind_q2dvdb_k(:,iq_ibz_k), cplex, nfftf, ngfftf, v1scf, gwpt%pert_comm%value)
       end if

       ! Rotate phonon frequencies and displacements for q in BZ. Non-blocking operation inside pert_comm
       !call timab(1901, 1, tsec)

       call phstore%async_rotate(cryst, ifc, iq_ibz, gwpt%qibz(:, iq_ibz), qpt, isym_q, trev_q)
       !call ifc%fourq(cryst, qpt, phfrq, displ_cart, out_displ_red=displ_red, comm=gwpt%pert_comm%value)

       ! Get npw_kq, kg_kq for k+q.
       call wfd%get_gvec_kq(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, istwf_kq, npw_kq, kg_kq)
       !call timab(1901, 2, tsec)
       !call timab(1902, 1, tsec)

       istwf_kqirr = wfd%istwfk(ikq_ibz); npw_kqirr = wfd%npwarr(ikq_ibz)
       ABI_MALLOC(bra_kq, (2, npw_kq*nspinor))
       ABI_MALLOC(cgwork, (2, npw_kqirr*nspinor))

       !if (osc_ecut /= zero) then
       !  ! Finds the boundary of the basis sphere of G vectors (for this kq point)
       !  ! for use in improved zero padding of ffts in 3 dimensions.
       !  call sphereboundary(gbound_kq, istwf_kq, kg_kq, wfd%mgfft, npw_kq)

       !  ! Compute "small" G-sphere centered on qpt and gbound for zero-padded FFT for oscillators.
       !  call get_kg(qpt, istw1, abs(osc_ecut), cryst%gmet, osc_npw, osc_gvecq)
       !  call sphereboundary(osc_gbound_q, istw1, osc_gvecq, wfd%mgfft, osc_npw)

       !  ! Compute correspondence G-sphere --> FFT mesh.
       !  ABI_MALLOC(osc_indpw, (osc_npw))
       !  ABI_MALLOC(osc_mask, (osc_npw))
       !  call kgindex(osc_indpw, osc_gvecq, osc_mask, wfd%mpi_enreg, ngfft, osc_npw)
       !  ABI_FREE(osc_mask)

       !  ABI_MALLOC(workq_ug, (npw_kq*nspinor))
       !  ABI_MALLOC(osc_ks, (osc_npw*nspinor, nbcalc_ks))
       !end if

       ! Allocate array to store H1 |psi_nk> for all 3*natom perturbations
       ABI_MALLOC_OR_DIE(h1kets_kq, (2, npw_kq*nspinor, my_npert, nbcalc_ks), ierr)

       ! Allocate vlocal1 with correct cplex. Note nvloc
       ABI_MALLOC_OR_DIE(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert), ierr)

       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))
       ABI_MALLOC(gvnlx1, (2, npw_kq*nspinor))

       ! Set up the spherical harmonics (Ylm) at k and k+q. See also dfpt_looppert
       !if (psps%useylm == 1) then
       !   optder = 0; if (useylmgr == 1) optder = 1
       !   call initylmg(cryst%gprimd, kg_k, kk, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1, &
       !     [npw_k], dtset%nsppol, optder, cryst%rprimd, ylm_k, ylmgr)
       !   call initylmg(cryst%gprimd, kg_kq, kq, mkmem1, mpi_enreg, psps%mpsang, mpw, nband, mkmem1, &
       !     [npw_kq], dtset%nsppol, optder, cryst%rprimd, ylm_kq, ylmgr_kq)
       !end if

       ! Compute k+q+G vectors
       nkpg1 = 3*dtset%nloalg(3)
       ABI_MALLOC(kpg1_k, (npw_kq, nkpg1))
       if (nkpg1 > 0) call mkkpg(kg_kq, kpg1_k, kq, nkpg1, npw_kq)

       ! Compute nonlocal form factors ffnl1 at (k+q+G)
       ABI_MALLOC(ffnl1, (npw_kq, 1, psps%lmnmax, psps%ntypat))

       call mkffnl_objs(cryst, psps, 1, ffnl1, ider0, idir0, kg_kq, kpg1_k, kq, nkpg1, npw_kq, ylm_kq, ylmgr_kq, &
                        comm=gwpt%pert_comm%value, request=ffnl1_request)

       if (dtset%eph_stern /= 0) then
         ! Build global array with GS wavefunctions cg_kq at k+q to prepare call to dfpt_cgwf.
         ! NB: bsum_range is not compatible with Sternheimer.
         ! There's a check at the level of the parser in chkinp.

         call timab(1908, 1, tsec)
         ABI_CALLOC(cg1s_kq, (2, npw_kq*nspinor, natom3, nbcalc_ks))

         ! NOTE that in the present version we need to gather all nbsum bands
         ! on each core before calling dfpt_cgwf.
         ! In principle one can call dfpt_cgwf in band-para mode but then
         ! we are obliged to call the sternheimer solver with one psi1 and all procs in bsum_comm
         ! just to to be able to apply the projector operator.
         ! The present version is not memory efficient and leads to a big load imbalance if
         ! bsum%comm%nproc > nband_calc_ks

!#define DEV_BAND_PARA

#ifdef DEV_BAND_PARA
         nband_me = gwpt%my_bsum_stop - gwpt%my_bsum_start + 1
#else
         nband_me = nbsum
#endif

         ABI_MALLOC(cgq, (2, npw_kq * nspinor, nband_me))
         ABI_MALLOC(gscq, (2, npw_kq * nspinor, nband_me*psps%usepaw))

         do ibsum_kq=gwpt%my_bsum_start, gwpt%my_bsum_stop

           if (isirr_kq) then
              call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
            else
              ! Reconstruct u_kq(G) from the IBZ image.
              call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)
              call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1, &
                               npw_kqirr, wfd%kdata(ikq_ibz)%kg_k, &
                               npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
            end if

#ifdef DEV_BAND_PARA
            ii = ibsum_kq - gwpt%my_bsum_start + 1
            cgq(:,:,ii) = bra_kq
#else
            cgq(:, :, ibsum_kq) = bra_kq
#endif
         end do

         cgq_request = xmpi_request_null

#ifndef DEV_BAND_PARA
         if (gwpt%bsum_comm%nproc > 1) then
           ! If band parallelism, need to gather all bands nbsum bands.
           ! FIXME: This part is network intensive, one can avoid it by calling dfpt_cgwf in band-para mode.
           !call xmpi_sum(cgq, gwpt%bsum_comm%value, ierr)
           !call xmpi_isum_ip(cgq, gwpt%bsum_comm%value, cgq_request, ierr)

           nelem = 2 * npw_kq * nspinor
           call gwpt%bsum_comm%prep_gatherv(nelem, gwpt%nbsum_rank(:,1), sendcount, recvcounts, displs)
#ifdef HAVE_MPI
           !call MPI_ALLGATHERV(MPI_IN_PLACE, sendcount, MPI_DOUBLE_PRECISION, cgq, recvcounts, displs, &
           !                    MPI_DOUBLE_PRECISION, gwpt%bsum_comm%value, ierr)

#if defined HAVE_MPI2_INPLACE
           call MPI_IALLGATHERV(MPI_IN_PLACE, sendcount, MPI_DOUBLE_PRECISION, cgq, recvcounts, displs, &
                                MPI_DOUBLE_PRECISION, gwpt%bsum_comm%value, cgq_request, ierr)
#else
           ABI_MALLOC(cgq_buf,(sendcount))
           me=1+xmpi_comm_rank(gwpt%bsum_comm%value)
           cgq_buf(1:sendcount)=cgq_ptr(displs(me)+1:displs(me)+sendcount)
           call c_f_pointer(c_loc(cgq),cgq_ptr,[2*npw_kq*nspinor*nband_me])
           call MPI_IALLGATHERV(cgq_buf, sendcount, MPI_DOUBLE_PRECISION, cgq_ptr, recvcounts, displs, &
                                MPI_DOUBLE_PRECISION, gwpt%bsum_comm%value, cgq_request, ierr)
           ABI_FREE(cgq_buf)
#endif
           call xmpi_requests_add(+1)
#endif

           ABI_FREE(recvcounts)
           ABI_FREE(displs)
         end if
#endif
         call timab(1908, 2, tsec)
       end if  ! eph_stern

       ! Loop over all 3*natom perturbations (Each core prepares its own potentials)
       ! In the inner loop, we calculate H1 * psi_k, stored in h1kets_kq on the k+q sphere.
       do imyp=1,my_npert
         idir = gwpt%my_pinfo(1, imyp); ipert = gwpt%my_pinfo(2, imyp); ipc = gwpt%my_pinfo(3, imyp)

         ! Set up local potential vlocal1 with proper dimensioning, from vtrial1 taking into account the spin.
         ! Each CPU prepares its own potentials.
         call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc, &
           pawfgr, mpi_enreg, vtrial, v1scf(:,:,:,imyp), vlocal, vlocal1(:,:,:,:,imyp))

         ! Continue to initialize the Hamiltonian (call it here to support dfpt_cgwf Sternheimer).
         call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)

         ! Prepare application of the NL part.
         call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)
         call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,imyp), with_nonlocal=.true.)

         if (ffnlk_request /= xmpi_request_null) call xmpi_wait(ffnlk_request, ierr)
         if (ffnl1_request /= xmpi_request_null) call xmpi_wait(ffnl1_request, ierr)

         ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
         call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk, kq, idir, ipert, &  ! In
           cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
           npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
           dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
           reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

         ! Compute H(1) applied to GS wavefunction Psi_nk(0)
         do ib_k=1,nbcalc_ks
           if (gwpt%bsum_comm%skip(ib_k, root=root_bcalc(ib_k))) cycle ! MPI parallelism inside bsum_comm
                                                                       ! Store rank treating ib_k in root_bcalc
           band_ks = ib_k + bstart_ks - 1
           eig0nk = ebands%eig(band_ks, ik_ibz, spin)
           ! Use scissor shift on 0-order eigenvalue
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, kets_k(:,:,ib_k), cwaveprj0, h1kets_kq(:,:,imyp, ib_k), &
             grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
             optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c1, usevnl)
         end do

         do ib_k=1,nbcalc_ks
           call xmpi_bcast(h1kets_kq(:,:,imyp,ib_k), root_bcalc(ib_k), gwpt%bsum_comm%value, ierr)
         end do

         if (dtset%eph_stern /= 0) then
           call timab(1909, 1, tsec)
           ! Activate Sternheimer. Note that we are still inside the MPI loop over my_npert.
           ! NB: Assume adiabatic AHC expression to compute the contribution of states above nbsum.

#ifdef DEV_BAND_PARA
           ! Prepare band parallelism in dfpt_cgwf via mpi_enreg.
           mpi_enreg%comm_band = gwpt%bsum_comm%value
           mpi_enreg%me_band = gwpt%bsum_comm%me
           mpi_enreg%nproc_band = gwpt%bsum_comm%nproc
#endif

           ABI_CALLOC(out_eig1_k, (2*nbsum**2))
           ABI_MALLOC(dcwavef, (2, npw_kq*nspinor*usedcwavef0))
           ABI_MALLOC(gh1c_n, (2, npw_kq*nspinor))
           ABI_MALLOC(ghc, (2, npw_kq*nspinor))
           ABI_MALLOC(gsc, (2, npw_kq*nspinor))
           ABI_MALLOC(gvnlxc, (2, npw_kq*nspinor))

           ! TODO: grad_berry is problematic because in dfpt_cgwf, the array is declared with
           !
           !  real(dp),intent(in) :: grad_berry(2,mpw1*nspinor,nband)
           !
           ! and
           !
           !  npw1_k = number of plane waves at this k+q point
           !
           ! So in principle we should allocate lot of memory to avoid bound checking error!
           ! For the time being use mpw1 = 0 because mpw1 is not used in this call to dfpt_cgwf
           ! still it's clear that the treatment of this array must be completely refactored in the DFPT code.
           !
           grad_berry_size_mpw1 = 0

           !TODO: to distribute cgq and kets memory, use mband_mem per core in band comm, but coordinate everyone with
           ! the following array (as opposed to the distribution of cg1 which is done in the normal dfpt calls
           ABI_MALLOC(bands_treated_now, (nbsum))
           ABI_MALLOC (rank_band, (nbsum))
           rank_band = 0

           nline_in = min(100, npw_kq); if (dtset%nline > nline_in) nline_in = min(dtset%nline, npw_kq)

#ifndef DEV_BAND_PARA
           ! Wait for gatherv operation
           if (cgq_request /= xmpi_request_null) call xmpi_wait(cgq_request, ierr)
#endif

           do ib_k=1,nbcalc_ks
             band_ks = ib_k + bstart_ks - 1
             bands_treated_now(:) = 0; bands_treated_now(band_ks) = 1

#ifdef DEV_BAND_PARA
             ! Init rank_band and band_me from nbsum_rank.
             rank_band = -1; band_me = 1
             do ip=1,gwpt%bsum_comm%nproc
               ii = gwpt%nbsum_rank(ip,2)
               jj = gwpt%nbsum_rank(ip,2) + gwpt%nbsum_rank(ip,1) -1
               rank_band(ii:jj) = ip - 1
               if (inrange(band_ks, [ii, jj])) u1_master = ip - 1
             end do
             if (inrange(band_ks, [gwpt%my_bsum_start, gwpt%my_bsum_stop])) then
               band_me = band_ks - gwpt%my_bsum_start + 1
               u1_band = band_ks
             else
               band_me = 1
               u1_band = -band_ks
             end if
#else
             rank_band = 0
             band_me = band_ks
             u1_band  = band_ks
             if (gwpt%bsum_comm%skip(ib_k)) cycle ! MPI parallelism inside bsum_comm
#endif

             ! Init entry in cg1s_kq, either from cache or with zeros.
             if (use_u1c_cache) then
               u1c_ib_k = u1c%find_band(band_ks)
               if (u1c_ib_k /= -1) then
                 call cgtk_change_gsphere(nspinor, &
                                          u1c%prev_npw_kq, istwfk1, u1c%prev_kg_kq, u1c%prev_cg1s_kq(1,1,ipc,u1c_ib_k), &
                                          npw_kq, istwfk1, kg_kq, cg1s_kq(1,1,ipc,ib_k), work_ngfft, work)
               else
                 cg1s_kq(:,:,ipc,ib_k) = zero
               end if

             else
               cg1s_kq(:,:,ipc,ib_k) = zero
             end if

             mcgq = npw_kq * nspinor * nband_me
             mgscq = npw_kq * nspinor * nband_me * psps%usepaw
             nlines_done = 0
             call timab(1909, 2, tsec)

             call dfpt_cgwf(u1_band, band_me, rank_band, bands_treated_now, berryopt0, &
               cgq, cg1s_kq(:,:,ipc,ib_k), kets_k(:,:,ib_k), &  ! Important stuff
               cwaveprj, cwaveprj0, rf2, dcwavef, &
               ebands%eig(:, ik_ibz, spin), ebands%eig(:, ikq_ibz, spin), out_eig1_k, &
               ghc, gh1c_n, grad_berry, gsc, gscq, &
               gs_hamkq, gvnlxc, gvnlx1, icgq0, idir, ipert, igscq0, &
               mcgq, mgscq, mpi_enreg, grad_berry_size_mpw1, cryst%natom, nbsum, nband_me, &
               nbdbuf0, nline_in, npw_k, npw_kq, nspinor, &
               opt_gvnlx1, dtset%prtvol, quit0, out_resid, rf_hamkq, dtset%dfpt_sciss, -one, dtset%tolwfr, &
               usedcwavef0, dtset%wfoptalg, nlines_done)

             tot_nlines_done = tot_nlines_done + nlines_done

#ifdef DEV_BAND_PARA
             call xmpi_bcast(cg1s_kq(:,:,ipc,ib_k), u1_master, gwpt%bsum_comm%value, ierr)
#endif

             ! Handle possible convergence error.
             if (u1_band > 0) then
               if (out_resid > dtset%tolwfr) then
                 write(msg, "(a,i0,a, 2(a,es13.5), 2a,i0,a)") &
                   " Sternheimer didn't convergence for band: ", band_ks, ch10, &
                   " resid:", out_resid, " >= tolwfr: ", dtset%tolwfr, ch10, &
                   " after nline: ", nlines_done, " iterations. Increase nline and/or tolwfr."
                 ABI_ERROR(msg)
               else if (out_resid < zero) then
                 ABI_ERROR(sjoin(" resid: ", ftoa(out_resid), ", nlines_done:", itoa(nlines_done)))
               end if

               if (my_rank == master .and. (enough_stern <= 5 .or. dtset%prtvol > 10)) then
                 write(std_out, "(2(a,es13.5),a,i0)") &
                   " Sternheimer converged with resid: ", out_resid, " <= tolwfr: ", dtset%tolwfr, &
                   " after nlines_done: ", nlines_done
                 enough_stern = enough_stern + 1
               end if
             end if
           end do ! ib_k

#ifdef DEV_BAND_PARA
           ! Revert changes in mpi_enreg.
           mpi_enreg%comm_band = xmpi_comm_self
           mpi_enreg%me_band = 0
           mpi_enreg%nproc_band = 1
#endif

           ABI_FREE(bands_treated_now)
           ABI_FREE(rank_band)
           ABI_FREE(out_eig1_k)
           ABI_FREE(dcwavef)
           ABI_FREE(gh1c_n)
           ABI_FREE(ghc)
           ABI_FREE(gsc)
           ABI_FREE(gvnlxc)
           if (imyp == my_npert) then
             ABI_FREE(cgq)
             ABI_FREE(gscq)
           end if
           !call timab(1909, 2, tsec)
         end if ! sternheimer

         call rf_hamkq%free()
         ABI_FREE(kinpw1)
         ABI_FREE(dkinpw)
         ABI_FREE(ph3d)
         ABI_SFREE(ph3d1)
       end do ! imyp  (loop over perturbations)

       !call timab(1902, 2, tsec)
       ABI_FREE(gs1c)
       ABI_FREE(gvnlx1)
       ABI_FREE(vlocal1)
       ABI_FREE(v1scf)

       ! Wait from phonon frequencies and displacements inside pert_comm
       call phstore%wait(cryst, phfrq, displ_cart, displ_red)

       if (dtset%eph_stern /= 0) then
         call timab(1910, 1, tsec)
         ! Add contribution to Fan-Migdal self-energy coming from Sternheimer.
         ! NB: All procs inside (bsum_comm x pert_comm) enter here!

         ! Store |Psi_1> to init Sternheimer solver for the next q-point.
         call u1c%store(qpt, npw_kq, nspinor, natom3, bstart_ks, nbcalc_ks, kg_kq, cg1s_kq)

         ! h1kets_kq are MPI distributed inside pert_comm but we need off-diagonal pp' terms --> collect results.
         ABI_CALLOC(h1kets_kq_allperts, (2, npw_kq*nspinor, natom3, nbcalc_ks))

         ! Compute S_pp' = <D_{qp} vscf u_nk|u'_{nk+q p'}>
         ABI_CALLOC(stern_ppb, (2, natom3, natom3, nbcalc_ks))

         do ib_k=1,nbcalc_ks
           if (gwpt%bsum_comm%skip(ib_k)) cycle ! MPI parallelism inside bsum_comm

           call xmpi_sum(cg1s_kq(:,:,:,ib_k), gwpt%pert_comm%value, ierr)

           ! TODO
           !nelem = 2*npw_kq*nspinor*gwpt%my_npert
           !call MPI_ALLGATHER(MPI_IN_PLACE, nelem, MPI_DOUBLE_PRECISION, cg1s_kq(:,:,:,ib_k), nelem, &
           !                   MPI_DOUBLE_PRECISION, gwpt%pert_comm%value, ierr)

           call xmpi_allgather(h1kets_kq(:,:,:,ib_k), 2*npw_kq*nspinor*gwpt%my_npert, &
                               h1kets_kq_allperts(:,:,:,ib_k), gwpt%pert_comm%value, ierr)

           call cg_zgemm("C", "N", npw_kq*nspinor, natom3, natom3, &
             h1kets_kq_allperts(:,:,:,ib_k), cg1s_kq(:,:,:,ib_k), stern_ppb(:,:,:,ib_k))

           ! Save data for Debye-Waller that is performed outside the q-loop.
           if (q_is_gamma) stern_dw(:,:,:,ib_k) = stern_ppb(:,:,:,ib_k)
         end do

         ABI_FREE(cg1s_kq)
         ABI_FREE(h1kets_kq_allperts)

         if (q_is_gamma) call xmpi_sum(stern_dw, gwpt%bsum_comm%value, ierr)

         ! Compute contribution to Fan-Migdal for M > gwpt%nbsum
         do imyp=1,my_npert
           nu = gwpt%my_pinfo(3, imyp)
           wqnu = phfrq(nu)

           do ib_k=1,nbcalc_ks
             if (gwpt%bsum_comm%skip(ib_k)) cycle ! MPI parallelism inside bsum_comm
             ! sum_{pp'} d_p* Stern_{pp'} d_p' with d = displ_red(:,:,:,nu) and S = stern_ppb(:,:,:,ib_k)
             vec_natom3 = zero
             call cg_zgemm("N", "N", natom3, natom3, 1, stern_ppb(:,:,:,ib_k), displ_red(:,:,:,nu), vec_natom3)
             dotri = cg_zdotc(natom3, displ_red(:,:,:,nu), vec_natom3)
             !write(std_out, *)"dotri:", dotri
             rfact = dotri(1)
             !rfact = cg_real_zdotc(natom3, displ_red(:,:,:,nu), vec_natom3)
             rfact = rfact * gwpt%wtq_k(iq_ibz_k) / (two * wqnu)
           end do
         end do

         ABI_FREE(stern_ppb)
         call timab(1910, 2, tsec)
       end if ! eph_stern /= 0

       ! ==============================================
       ! Sum over m bands parallelized inside bsum_comm
       ! ==============================================
       call timab(1903, 1, tsec)

       do ibsum_kq=gwpt%my_bsum_start, gwpt%my_bsum_stop
         call timab(1904, 1, tsec)
         ! This can happen if we have loaded the wavefunctions inside the energy range.
         !if (.not. wfd%ihave_ug(ibsum_kq, ikq_ibz, spin)) then
         !    ignore_ibsum_kq = ignore_ibsum_kq + 1; cycle
         !end if

         ! Symmetrize k+q wavefunctions in the BZ from IBZ (if needed).
         if (isirr_kq) then
           ! Copy u_kq(G)
           call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, bra_kq)
         else
           ! Reconstruct u_kq(G) from the IBZ image.
           ! Use cgwork as workspace array, results stored in bra_kq
           ! g0_kq = g0ibz_kq + g0bz_kq
           call wfd%copy_cg(ibsum_kq, ikq_ibz, spin, cgwork)
           call cgtk_rotate(cryst, kq_ibz, isym_kq, trev_kq, g0_kq, nspinor, ndat1, &
                            npw_kqirr, wfd%kdata(ikq_ibz)%kg_k, &
                            npw_kq, kg_kq, istwf_kqirr, istwf_kq, cgwork, bra_kq, work_ngfft, work)
         end if

         ! Get gkk(kcalc, q, idir_ipert) in the atomic representation.
         ! No need to handle istwf_kq because it's always 1.
         gkq_atm = zero; cnt = 0
         do imyp=1,my_npert
           ipc = gwpt%my_pinfo(3, imyp)
           ! Calculate <u_(band,k+q)^(0)|H_(k+q,k)^(1)|u_(band,k)^(0)> for this pert (NC psps) istwf_k always 1
           do ib_k=1,nbcalc_ks
             gkq_atm(:, ib_k, ipc) = cg_zdotc(npw_kq*nspinor, bra_kq, h1kets_kq(:,:,imyp,ib_k))
             cnt = cnt + 1
             gkq_allgather(:,cnt, 1) = gkq_atm(:, ib_k, ipc)
           end do
           !call cg_zgemv("C", npw_kq*nspinor, nbcalc_ks, h1kets_kq(:,:,:,imyp), bra_kq, gkq_atm(:,:,ipc))
         end do
         call timab(1904, 2, tsec)
         call timab(1905, 1, tsec)
         !ii = nbcalc_ks * my_npert
         !call cg_zgemm("H", "N", npw_kq*nspinor, ii, ii, h1kets_kq, bra_kq, gkq_atm)
         !call cg_zgemm("H", "N", npw_kq*nspinor, ii, ii, bra_kq, h1kets_kq, gkq_atm)

         ! Get gkk(kcalc, q, nu) in the phonon representation.
         ! Need to gather all perts distributed in pert_comm
         if (gwpt%pert_comm%nproc > 1) then
           call xmpi_allgather(gkq_allgather(:,:,1), 2 * nbcalc_ks * my_npert, gkq_allgather(:,:,2), &
                               gwpt%pert_comm%value, ierr)
           do cnt=1,nbcalc_ks*natom3
             ipc = 1 + (cnt - 1) / nbcalc_ks
             ib_k = 1 + mod(cnt - 1, nbcalc_ks)
             gkq_atm(:, ib_k, ipc) = gkq_allgather(:, cnt, 2)
           end  do
         end if

         call ephtk_gkknu_from_atm(1, nbcalc_ks, 1, natom, gkq_atm, phfrq, displ_red, gkq_nu)

         ! bsum_2 and bsum_3 are hotspots.
         call timab(1905, 2, tsec)
         call timab(1906, 1, tsec)

         !if (osc_ecut > zero) then
         !  workq_ug = cmplx(bra_kq(1, :), bra_kq(2, :), kind=gwpc)
         !  call fft_ug(npw_kq, wfd%nfft, nspinor, ndat1, wfd%mgfft, wfd%ngfft, &
         !              istwf_kq, kg_kq, gbound_kq, workq_ug, ur_kq)

         !  ! We need <k+q| e^{iq+G}|k> --> compute <k| e^{-i(q+G)}|k+q> with FFT and take CC.
         !  do ib_k=1,nbcalc_ks
         !    work_ur = ur_kq * conjg(ur_k(:, ib_k))
         !    ! Call zero-padded FFT routine.
         !    call fftpad(work_ur, ngfft, n1, n2, n3, n1, n2, n3, nspinor, wfd%mgfft, -1, osc_gbound_q)

         !    ! Need results on the G-sphere --> Transfer data from FFT to G-sphere.
         !    do ispinor=1,nspinor
         !      do ig=1,osc_npw
         !        ifft = osc_indpw(ig) + (ispinor-1) * wfd%nfft
         !        osc_ks(ig + (ispinor -1) * osc_npw, ib_k) = conjg(work_ur(ifft))
         !      end do
         !    end do

         !    !band_ks = ib_k + bstart_ks - 1
         !    !if (ibsum_kq == band_ks) then
         !    !if (ibsum_kq == band_ks .and. all(abs(qpt) < tol12)) then
         !    !  write(std_out,"(a,i0,2a)")" Ene and Oscillator for band: ", band_ks, ", and q-point: ", trim(ktoa(qpt))
         !    !  write(std_out,*)ebands%eig(band_ks, ik_ibz, spin) * Ha_eV, osc_ks(:2,ib_k)
         !    !end if
         !  end do
         !end if

         eig0mkq = ebands%eig(ibsum_kq, ikq_ibz, spin)

         ! q-weight for naive integration
         weight_q = gwpt%wtq_k(iq_ibz_k)
         call timab(1906, 2, tsec)
         call timab(1907, 1, tsec)

         ! Accumulate contribution to the FM self-energy
         do imyp=1,my_npert
           nu = gwpt%my_pinfo(3, imyp)
           ! Ignore unstable modes or modes that should be skipped.
           wqnu = phfrq(nu)

           ! For each band in gwpt_{nk}
           do ib_k=1,nbcalc_ks
             band_ks = ib_k + bstart_ks - 1
             eig0nk = ebands%eig(band_ks, ik_ibz, spin)
             gkq2 = weight_q * (gkq_nu(1,ib_k,nu) ** 2 + gkq_nu(2,ib_k,nu) ** 2)
             ediff = eig0nk - eig0mkq
             intra_band = q_is_gamma .and. ediff <= TOL_EDIFF
             same_band = ibsum_kq == band_ks

           end do ! ib_k
         end do ! imyp
         call timab(1907, 2, tsec)

       end do ! ibsum_kq (sum over bands at k+q)
       call timab(1903, 2, tsec)

       ABI_FREE(bra_kq)
       ABI_FREE(cgwork)
       ABI_FREE(h1kets_kq)
       ABI_FREE(kpg1_k)
       ABI_FREE(ffnl1)

       !if (osc_ecut /= zero) then
       !  ABI_FREE(osc_gvecq)
       !  ABI_FREE(osc_indpw)
       !  ABI_FREE(osc_ks)
       !  ABI_FREE(workq_ug)
       !end if

       if (imyq <= 10 .or. mod(imyq, 100) == 0) then
         write(msg,'(4(a,i0),a)') " k-point [",my_ikcalc,"/",gwpt%my_nkcalc, "] q-point [",imyq,"/",gwpt%my_nqibz_k,"]"
         call cwtime_report(msg, cpu, wall, gflops)
       end if
     end do ! iq_ibz_k (sum over q-points in IBZ_k)

     call cwtime_report(" Fan-Migdal q-loop", cpu_qloop, wall_qloop, gflops_qloop)

     ! Print cache stats.
     if (gwpt%use_ftinterp) then
       call dvdb%ft_qcache%report_stats()
       if (dvdb%ft_qcache%v1scf_3natom_request /= xmpi_request_null) call xmpi_wait(dvdb%ft_qcache%v1scf_3natom_request, ierr)
     else
       call dvdb%qcache%report_stats()
     end if

     ABI_FREE(kets_k)
     ABI_FREE(gkq_atm)
     ABI_FREE(gkq_nu)
     ABI_FREE(gkq_allgather)

     !if (osc_ecut /= zero) then
     !  ABI_FREE(ur_k)
     !  ABI_FREE(ur_kq)
     !  ABI_FREE(work_ur)
     !  ABI_FREE(gkq2_lr)
     !end if

     if (my_rank == master) then
       if (ignore_kq /= 0) write(std_out, "(a, 1x, i0)")" Number of ignored k+q points:", ignore_kq
       if (ignore_ibsum_kq /= 0) write(std_out, "(a, 1x, i0)")" Number of ignored (k+q, m) states:", ignore_ibsum_kq
       !if (dtset%eph_stern /= 0) then
       !  call wrtout(std_out, sjoin(" Total number of NSCF Sternheimer iterations:", itoa(tot_nlines_done)))
       !end if
     end if

     ! Collect results inside pqb_comm and write results for this (k-point, spin) to NETCDF file.
     !call gwpt%gather_and_write(dtset, ebands, ikcalc, spin, gwpt%pqb_comm%value)

     ABI_SFREE(root_bcalc)
   end do ! spin

   ABI_FREE(kg_k)
   ABI_FREE(kg_kq)
   ABI_FREE(ylm_k)
   ABI_FREE(ylm_kq)
   ABI_FREE(ylmgr_kq)
   ABI_FREE(kpg_k)
   ABI_FREE(ffnlk)

   call cwtime_report(" One ikcalc k-point", cpu_ks, wall_ks, gflops_ks)
 end do ! ikcalc

 call cwtime_report(" gwpt_eph full calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Free memory
 ABI_FREE(ihave_ikibz_spin)
 ABI_FREE(grad_berry)
 ABI_FREE(vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red)
 ABI_FREE(gbound_kq)
 ABI_FREE(osc_gbound_q)
 ABI_FREE(ibzspin_2ikcalc)
 ABI_SFREE(vcar_ibz)

 call gs_hamkq%free()
 call wfd%free()
 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj)
 call phstore%free()
 call u1c%free()
 call gwpt%free()

 ! This to make sure that the parallel output of SIGEPH is completed
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

type(gwpt_t) function gwpt_new(dtset, ecut, cryst, ebands, ifc, dtfil, comm) result(gwpt)

!Arguments ------------------------------------
 integer,intent(in) :: comm
 real(dp),intent(in) :: ecut
 type(crystal_t),intent(in) :: cryst
 type(dataset_type),intent(in) :: dtset
 type(ebands_t),intent(in) :: ebands
 type(ifc_type),intent(in) :: ifc
 type(datafiles_type),intent(in) :: dtfil

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, istwfk1 = 1
 integer :: my_rank,ik,my_nshiftq,my_mpw,cnt,nprocs,ik_ibz,ndeg, iq_ibz, qptopt, qtimrev
 integer :: onpw, ii, ipw, ierr, spin, gap_err, ikcalc, qprange_, bstop !it,
 integer :: jj, bstart, natom, natom3
 integer :: isym_k, trev_k, mband, i1,i2,i3, nrest, color
 character(len=5000) :: msg
 real(dp) :: estep, cpu_all, wall_all, gflops_all, cpu, wall, gflops
 logical :: changed, isirr_k
 type(gaps_t) :: gaps
 type(krank_t) :: krank, qrank
!arrays
 integer :: intp_nshiftk
 integer :: intp_kptrlatt(3,3), g0_k(3), units(2), indkk_k(6,1), my_gmax(3), band_block(2), qptrlatt(3,3)
 integer,allocatable :: temp(:,:), degblock(:,:), degblock_all(:,:,:,:), ndeg_all(:,:), iperm(:)
 real(dp):: params(4), my_shiftq(3,1), kk(3), kq(3), intp_shiftk(3)
#ifdef HAVE_MPI
 integer,parameter :: ndims = 5
 integer :: comm_cart, me_cart
 logical :: reorder
 integer :: dims(ndims)
 logical :: periods(ndims), keepdim(ndims)
#endif

! *************************************************************************

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

 qrank = krank_from_kptrlatt(gwpt%nqibz, gwpt%qibz, qptrlatt, compute_invrank=.False.)
 if (kpts_map("symrec", qtimrev, cryst, qrank, gwpt%nqbz, gwpt%qbz, temp) /= 0) then
   ABI_ERROR("Cannot map qBZ to qIBZ!")
 end if

 call qrank%free()

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

 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)
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

   if (kpts_map("symrel", gwpt%timrev, cryst, krank, 1, kk, indkk_k) /= 0) then
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

 call krank%free()
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

 ! Compute mpw and gmax
 call ephtk_get_mpw_gmax(gwpt%nkcalc, gwpt%kcalc, ecut, cryst%gmet, gwpt%mpw, gwpt%gmax, comm)
 call wrtout(std_out, sjoin(' Optimal value of mpw:', itoa(gwpt%mpw), "gmax:", ltoa(gwpt%gmax)))
 call cwtime_report(" gwpt_new: mpw", cpu, wall, gflops)

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

   ! TODO: Spin
   ! Automatic grid generation over q-points and spins.
   !if (gwpt%nsppol == 2 .and. mod(nprocs, 2) == 0) then
   !  spin_comm%nproc = 2
   !  gwpt%qpt_comm%nproc = nprocs / 2
   !else
   !  gwpt%qpt_comm%nproc = nprocs
   !end if

   ! Handle parallelism over perturbations first.
   ! Use MPI communicator to distribute the 3 * natom perturbations to reduce memory requirements for DFPT potentials.
   ! Ideally, perturbations are equally distributed --> total number of CPUs should be divisible by 3 * natom.
   ! or at least, divisible by one integer i for i in [2, 3 * natom - 1].

   ! Try to have 3 perts per proc first because the q-point parallelism is more efficient.
   ! The memory for W(R,r,ipert) will increase though.
   !do cnt=natom,2,-1
   !  if (mod(nprocs, cnt) == 0 .and. mod(natom3, cnt) == 0) then
   !    gwpt%pert_comm%nproc = cnt; gwpt%my_npert = natom3 / cnt; exit
   !  end if
   !end do

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
 keepdim = .False.; keepdim(2:3) = .True.; call gwpt%qb_comm%from_cart_sub(comm_cart, keepdim)
 ! Create communicator for the (perturbation, band_sum, qpoint_sum)
 keepdim = .False.; keepdim(1:3) = .True.; call gwpt%pqb_comm%from_cart_sub(comm_cart, keepdim)

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
 ABI_ICALLOC(gwpt%qp_done, (gwpt%nkcalc, gwpt%nsppol))

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
 ABI_SFREE(gwpt%myq2ibz_k)
 ABI_SFREE(gwpt%itreat_qibz)
 ABI_SFREE(gwpt%my_pinfo)
 ABI_SFREE(gwpt%pert_table)
 ABI_SFREE(gwpt%ind_qbz2ibz)
 ABI_SFREE(gwpt%indkk_kq)
 ABI_SFREE(gwpt%ind_q2dvdb_k)
 ABI_SFREE(gwpt%ind_ibzk2ibz)
 ABI_SFREE(gwpt%qibz2dvdb)
 ABI_SFREE(gwpt%lgk_sym2glob)
 ABI_SFREE(gwpt%nbsum_rank)

 ! real
 ABI_SFREE(gwpt%kcalc)
 ABI_SFREE(gwpt%vcar_calc)
 ABI_SFREE(gwpt%qp_done)
 ABI_SFREE(gwpt%qbz)
 ABI_SFREE(gwpt%qibz)
 ABI_SFREE(gwpt%wtq)
 ABI_SFREE(gwpt%qibz_k)
 ABI_SFREE(gwpt%wtq_k)

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
 call gwpt%qb_comm%free()
 call gwpt%kcalc_comm%free()
 call gwpt%spin_comm%free()
 call gwpt%pqb_comm%free()
 call gwpt%ncwrite_comm%free()

 ! Close netcdf file.
 if (gwpt%ncid /= nctk_noid) then
   NCF_CHECK(nf90_close(gwpt%ncid))
 end if

end subroutine gwpt_free
!!***

!!****f* m_gwpt/gwpt_setup_kcalc
!! NAME
!!  gwpt_setup_kcalc
!!
!! FUNCTION
!!  Prepare calculations of self-energy matrix elements for ikcalc index.
!!
!! INPUTS
!!  dtset<dataset_type>=All input variables for this dataset.
!!  cryst<crystal_t> = Crystal structure.
!!  dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!!  ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!!  ikcalc=Index of the k-point to compute.
!!  prtvol= Verbosity level
!!  comm= MPI communicator
!!
!! SOURCE

subroutine gwpt_setup_kcalc(self, dtset, cryst, ebands, ikcalc, prtvol, comm)

!Arguments ------------------------------------
 integer,intent(in) :: ikcalc, prtvol, comm
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 class(gwpt_t),target,intent(inout) :: self
 type(ebands_t),intent(in) :: ebands

!Local variables-------------------------------
 integer,parameter :: master = 0
 integer :: spin, my_rank, iq_ibz, nprocs, qtimrev, qptopt  !, nbcalc_ks !, bstart_ks
 integer :: ikpt, ibz_k, isym_k, itim_k !isym_lgk,
 real(dp) :: cpu, wall, gflops
 character(len=5000) :: msg
 type(lgroup_t),target :: lgk
 type(lgroup_t),pointer :: lgk_ptr
 type(krank_t) :: krank, qrank
!arrays
 integer :: qptrlatt(3,3)
 integer,allocatable :: iqk2dvdb(:,:)
 real(dp) :: kk(3)
 real(dp),allocatable :: kq_list(:,:)

! *************************************************************************

 ABI_SFREE(self%qibz_k)
 ABI_SFREE(self%wtq_k)

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 kk = self%kcalc(:, ikcalc)

 call wrtout(std_out, sjoin(ch10, repeat("=", 92)))
 msg = sjoin("[", itoa(ikcalc), "/", itoa(self%nkcalc), "]")
 call wrtout(std_out, sjoin(" Computing self-energy matrix elements for k-point:", ktoa(kk), msg))
 ! TODO Integrate with spin parallelism.
 spin = 1
 write(msg, "(3(a, i0))")" Treating ", self%nbcalc_ks(ikcalc, spin), " band(s) in gwpt_nk between: ", &
   self%bstart_ks(ikcalc, spin)," and: ", self%bstart_ks(ikcalc, spin) + self%nbcalc_ks(ikcalc, spin) - 1
 call wrtout(std_out, msg)
 write(msg, "(2(a,i0))")"P Allocating and summing bands from my_bsum_start: ", self%my_bsum_start, &
     " up to my_bsum_stop: ", self%my_bsum_stop
 call wrtout(std_out, msg)
 if (dtset%eph_stern /= 0) then
   if (dtset%eph_stern ==  1) call wrtout(std_out, " Sternheimer method activated with cache for u1_nk")
   if (dtset%eph_stern == -1) call wrtout(std_out, " Sternheimer method activated WITHOUT cache for u1_nk!")
 end if

 call cwtime(cpu, wall, gflops, "start")

 if (self%symsigma == 0) then
   ! Do not use symmetries in BZ sum_q --> nqibz_k == nqbz
   self%nqibz_k = self%nqbz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = self%qbz; self%wtq_k = one / self%nqbz
   call wrtout(std_out, sjoin(" symgsigma = 0 --> Integration done over full BZ with nqbz:", itoa(self%nqibz_k)))

   ! Store little group symmetries (well, just 1)
   self%lgk_nsym = 1
   ABI_REMALLOC(self%lgk_sym2glob, (2, self%lgk_nsym))
   self%lgk_sym2glob(:, 1) = [1, 1]

 else if (abs(self%symsigma) == 1) then
   ! Use the symmetries of the little group of the k-point
   ! Pack points in *shells* to minimise cache misses.
   lgk = lgroup_new(cryst, kk, self%timrev, self%nqbz, self%qbz, self%nqibz, self%qibz, comm)
   lgk_ptr => lgk

   ! Store little group symmetries.
   self%lgk_nsym = lgk_ptr%nsym_lg
   ABI_REMALLOC(self%lgk_sym2glob, (2, self%lgk_nsym))
   self%lgk_sym2glob = lgk_ptr%lgsym2glob

   call wrtout(std_out, sjoin(" Number of operations in little group(k):", itoa(lgk_ptr%nsym_lg), &
     "(including time-reversal symmetry)"))
   call wrtout(std_out, sjoin(" Number of q-points in the IBZ(k):", itoa(lgk_ptr%nibz)))

   if (dtset%prtvol > 0) call lgk_ptr%print(unit=std_out, prtvol=dtset%prtvol)

   ! TODO: Pointers instead of copies to save space?
   self%nqibz_k = lgk_ptr%nibz
   ABI_MALLOC(self%qibz_k, (3, self%nqibz_k))
   ABI_MALLOC(self%wtq_k, (self%nqibz_k))
   self%qibz_k = lgk_ptr%ibz; self%wtq_k = lgk_ptr%weights
 else
   ABI_ERROR(sjoin("Wrong symsigma:", itoa(self%symsigma)))
 end if

 call cwtime_report(" lgroup_symgwpt", cpu, wall, gflops)

 ! TODO: Cleanup

 if (self%symsigma == 0) then
   ! Find correspondence IBZ_k --> IBZ
   ABI_MALLOC(iqk2dvdb, (6, self%nqibz_k))

   ! Assume qptopt == kptopt unless value is specified in input
   qptopt = ebands%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
   qtimrev = kpts_timrev_from_kptopt(qptopt)
   qptrlatt = 0; qptrlatt(1,1) = self%ngqpt(1); qptrlatt(2,2) = self%ngqpt(2); qptrlatt(3,3) = self%ngqpt(3)
   qrank = krank_from_kptrlatt(self%nqibz, self%qibz, qptrlatt, compute_invrank=.False.)

   if (kpts_map("symrec", qtimrev, cryst, qrank, self%nqibz_k, self%qibz_k, iqk2dvdb) /= 0) then
     write(msg, '(3a)' )&
       "At least one of the q points in the IBZ_k could not be generated from one in the IBZ.", ch10,&
       "Action: check your DVDB file and use eph_task to interpolate the potentials on a denser q-mesh."
     ABI_ERROR(msg)
   end if
   call qrank%free()

   ABI_REMALLOC(self%ind_ibzk2ibz, (6, self%nqibz_k))
   do iq_ibz=1,self%nqibz_k
     self%ind_ibzk2ibz(:, iq_ibz) = iqk2dvdb(:, iq_ibz)
   end do
   ABI_FREE(iqk2dvdb)

 else if (abs(self%symsigma) == 1) then

   ! IBZ_k --> BZ --> IBZ
   ABI_REMALLOC(self%ind_ibzk2ibz, (6, self%nqibz_k))
   self%ind_ibzk2ibz = 0
   do ikpt=1,self%nqbz
     ibz_k    = lgk_ptr%bz2ibz_smap(1,ikpt)
     !isym_lgk = lgk_ptr%bz2ibz_smap(2,ikpt)
     !isym_k   = lgk_ptr%lgsym2glob(1,isym_lgk)
     !itim_k   = lgk_ptr%lgsym2glob(2,isym_lgk)
     isym_k   = lgk_ptr%bz2ibz_smap(2,ikpt)
     itim_k   = lgk_ptr%bz2ibz_smap(3,ikpt)
     ! I assume that isym=1 and itim_k=0 is identity but still verify the kpoint
     if (isym_k /= 1 .or. itim_k /= 1 .or. any(lgk_ptr%bz2ibz_smap(4:,ikpt) /= 0)) cycle
     ! check IBZ_k --> BZ
     ABI_CHECK(sum(abs(self%qbz(:,ikpt) - self%qibz_k(:,ibz_k))) < tol8, 'Wrong mapping')
     ! IBZ_k --> IBZ
     !self%ind_ibzk2ibz(:, ibz_k) = self%ind_qbz2ibz(:,ikpt)
     self%ind_ibzk2ibz(1, ibz_k) = self%ind_qbz2ibz(1, ikpt)
     self%ind_ibzk2ibz(2, ibz_k) = self%ind_qbz2ibz(2, ikpt)
     self%ind_ibzk2ibz(6, ibz_k) = self%ind_qbz2ibz(3, ikpt)
     self%ind_ibzk2ibz(3:5, ibz_k) = self%ind_qbz2ibz(4:6, ikpt)
   end do
   do ikpt=1,self%nqibz_k
     ABI_CHECK(self%ind_ibzk2ibz(1, ikpt) /= 0, 'Did not find mapping')
   end do
   call lgk%free()
 else
   ABI_ERROR(sjoin("Wrong symsigma:", itoa(self%symsigma)))
 endif

 call cwtime_report(" IBZ_k --> IBZ", cpu, wall, gflops)

 if (.not. self%use_ftinterp) then
   ! Find correspondence IBZ_k --> set of q-points in DVDB.
   ! Need to handle q_bz = S q_ibz by symmetrizing the potentials already available in the DVDB.
   !
   ! Note:
   !   q --> -q symmetry is always used for phonons.
   !   we use symrec instead of symrel (see also m_dvdb)
   ! IBZ_K -> BZ -> IBZ -> DVDB
   ABI_REMALLOC(self%ind_q2dvdb_k, (6, self%nqibz_k))
   self%ind_q2dvdb_k = self%ind_ibzk2ibz
   do ikpt=1,self%nqibz_k
     self%ind_q2dvdb_k(1, ikpt) = self%qibz2dvdb(self%ind_ibzk2ibz(1, ikpt))
   end do
   call cwtime_report(" IBZ_k --> DVDB", cpu, wall, gflops)
 end if

 ! Find k+q in the extended zone and extract symmetry info.
 ! Be careful here because there are two umklapp vectors to be considered:
 !
 !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
 !
 ! Note symrel and use_symrec=.False. in get_mapping.
 ! This means that this table can be used to symmetrize wavefunctions in cgtk_rotate.
 !
 ABI_MALLOC(kq_list, (3, self%nqibz_k))
 do iq_ibz=1,self%nqibz_k
   kq_list(:, iq_ibz) = kk + self%qibz_k(:,iq_ibz)
 end do

 ! Use iqk2dvdb as workspace array.
 ABI_MALLOC(iqk2dvdb, (6, self%nqibz_k))

 krank = krank_from_kptrlatt(ebands%nkpt, ebands%kptns, ebands%kptrlatt, compute_invrank=.False.)

 if (kpts_map("symrel", self%timrev, cryst, krank, self%nqibz_k, kq_list, iqk2dvdb) /= 0) then
   write(msg, '(11a)' )&
    "The WFK file cannot be used to compute self-energy corrections at k: ", trim(ktoa(kk)), ch10,&
    "At least one of the k+q points could not be generated from a symmetrical one.", ch10,&
    "Q-mesh: ",trim(ltoa(self%ngqpt)),", K-mesh (from kptrlatt) ",trim(ltoa(get_diag(dtset%kptrlatt))),ch10, &
    "Action: check your WFK file and the k/q point input variables."
   ABI_ERROR(msg)
 end if

 call krank%free()

 ABI_FREE(kq_list)

 ABI_REMALLOC(self%indkk_kq, (6, self%nqibz_k))
 do iq_ibz=1,self%nqibz_k
   self%indkk_kq(:, iq_ibz) = iqk2dvdb(:,iq_ibz)
 end do
 ABI_FREE(iqk2dvdb)

 call cwtime_report(" k+q --> ebands", cpu, wall, gflops)

end subroutine gwpt_setup_kcalc
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

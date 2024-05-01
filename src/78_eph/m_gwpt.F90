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

 public :: gwpt_run        ! Main entry point to compute GWpt matrix elements

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
 integer,parameter :: useylmgr = 0, useylmgr1 =0, master = 0, ndat1 = 1, with_cplex0 = 0
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
 integer :: n1,n2,n3,n4,n5,n6,nspden,nu, mqmem, mm_kq, nn_k, restart
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg,nkpg1,cnt,imyp !, restart
 integer :: tot_nlines_done, nlines_done, nline_in, grad_berry_size_mpw1, enough_stern
 integer :: nbcalc_ks,nbsum,bsum_start, bsum_stop, bstart_ks,ikcalc,bstart,bstop, sendcount !iatom,
 integer :: ipp_bz, comm_rpt, nqlwl, ebands_timrev ! osc_npw,
 integer :: ffnl_k_request, ffnl_kq_request, nelem, cgq_request, root_ncid
 integer :: nkibz, nkbz, qptopt, my_iq, my_ik
 real(dp) :: cpu,wall,gflops,cpu_all,wall_all,gflops_all,cpu_ks,wall_ks,gflops_ks
 real(dp) :: cpu_setk, wall_setk, gflops_setk, cpu_qloop, wall_qloop, gflops_qloop
 real(dp) :: ecut,eshift,weight_q,rfact,ediff,q0rad, out_resid, bz_vol
 logical :: isirr_k, isirr_kq, isirr_kmp, isirr_kqmp, isirr_q, gen_eigenpb, qq_is_gamma, pp_is_gamma
 logical :: stern_use_cache, intra_band, same_band, stern_has_band_para, use_ftinterp
 complex(dpc) :: ieta
 type(krank_t) :: qrank_ibz ! krank_ibz,
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_ham_kq
 type(rf_hamiltonian_type) :: rf_ham_kq
 !type(gwpt_t) :: gwpt
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
 integer,allocatable :: displs(:), recvcounts(:), qibz2dvdb(:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:), my_iqibz_inds(:), iperm(:)
 integer(i1b),allocatable :: itreat_qibz(:)
 integer,pointer :: kg_pp(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3), kqmp(3), kmp(3), pp(3), kmp_ibz(3), kqmp_ibz(3)
 real(dp) :: phfrq(3*cryst%natom), dotri(2), qq_ibz(3), qpt(3) !qpt_cart(3),
 real(dp) :: tsec(2) ! vk(3), vkq(3),
 real(dp) :: vec_natom3(2, 3*cryst%natom)
 real(dp) :: wqnu,gkq2,eig0nk,eig0mkq !eig0mk,
 real(dp),allocatable,target :: cgq(:,:,:)
 real(dp),allocatable :: qlwl(:,:),vk_cart_ibz(:,:,:,:)
 real(dp),allocatable :: displ_cart(:,:,:,:),displ_red(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg1_k(:,:),kpg_k(:,:),dkinpw(:) ! grad_berry(:,:),
 real(dp),allocatable :: ffnl_k(:,:,:,:),ffnl_kq(:,:,:,:),ph3d(:,:,:),ph3d1(:,:,:),v1scf(:,:,:,:)
 real(dp),allocatable :: gkq_atm(:,:,:),gkq_nu(:,:,:),gkq0_atm(:,:,:,:)
 real(dp),allocatable :: gscq(:,:,:), out_eig1_k(:), cg1s_kq(:,:,:,:), h1kets_kq_allperts(:,:,:,:)
 real(dp),allocatable :: dcwavef(:, :), gh1c_n(:, :), ghc(:,:), gsc(:,:), stern_ppb(:,:,:,:), stern_dw(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
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

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor
 nspden = dtset%nspden; nkpt = ebands%nkpt

 ieta = +j_dpc * dtset%zcut

 ! Check if a previous GSTORE file is present to restart the calculation
 ! Here we try to read an existing GSTORE file if eph_restart == 1.
 ! and we compare the variables with the state of the code

 restart = 0; ierr = 1; gstore_filepath = strcat(dtfil%filnam_ds(4), "_GSTORE.nc")
 if (my_rank == master .and. dtset%eph_restart == 1) then
    if (file_exists(gstore_filepath)) then
       !NCF_CHECK(nctk_open_read(root_ncid, gstore_filepath, xmpi_comm_self))
       !NCF_CHECK(nctk_get_dim(ncid, "gstore_nqbz", gstore%nqbz))
       !ABI_MALLOC(done_qbz_spin, (gstore%nqbz, nsppol))
       !NCF_CHECK(nf90_get_var(root_ncid, nctk_idname(root_ncid, "gstore_done_qbz_spin"), done_qbz_spin))
       !NCF_CHECK(nf90_close(root_ncid))

       !restart = 1
       !if (any(sigma_restart%qp_done /= 1)) then
       !  restart = 1
       !  call wrtout(units, "- Restarting from previous SIGEPH.nc file")
       !  call wrtout(units, sjoin("- Number of k-points completed:", itoa(count(sigma%qp_done == 1)), "/", itoa(sigma%nkcalc)))
       !else
       !  restart = 0; sigma%qp_done = 0
       !  msg = sjoin("Found SIGEPH.nc file with all QP entries already computed.", ch10, &
       !              "Will overwrite:", gstore_filepath, ch10, &
       !              "Keeping backup copy in:", strcat(gstore_filepath, ".bkp"))
       !  call wrtout(ab_out, sjoin("WARNING: ", msg))
       !  ABI_WARNING(msg)
       !  ! Keep backup copy
       !  ABI_CHECK(clib_rename(gstore_filepath, strcat(gstore_filepath, ".bkp")) == 0, "Failed to rename SIGPEPH file.")
       !end if
    end if
 end if

 call xmpi_bcast(restart, master, comm, ierr)
 !call xmpi_bcast(gwpt%qp_done, master, comm, ierr)

 if (restart == 0) then
   call gstore%init(gstore_filepath, dtset, wfk_hdr, cryst, ebands, ifc, comm)
   ABI_ICALLOC(done_qbz_spin, (gstore%nqbz, nsppol))
 else
   call gstore%from_ncpath(gstore_filepath, with_cplex0, dtset, cryst, ebands, ifc, comm)
   NOT_IMPLEMENTED_ERROR()
   ! TODO: Read mask
 end if

 !ndone = count(done_qbz_spin == 1)

 !if (restart == 0) then
 !  call gwpt%write(dtset, cryst, ebands, wfk_hdr, dtfil, comm)
 !else
 !  ! Open file inside ncwrite_comm to perform parallel IO if kpt parallelism.
 !  if (gwpt%ncwrite_comm%value /= xmpi_comm_null) then
 !    NCF_CHECK(nctk_open_modify(gwpt%ncid, gstore_filepath, gwpt%ncwrite_comm%value))
 !    NCF_CHECK(nctk_set_datamode(gwpt%ncid))
 !  end if
 !end if

 ! FFT meshes from input file, not necessarly equal to the ones found in the external files.
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3)
 n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 ecut = dtset%ecut ! dtset%dilatmx

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

 ! TODO: Distribute wavefunctions according to k, q and pp
 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.
 !bks_mask(gwpt%my_bsum_start:gwpt%my_bsum_stop,:,:) = .True.
 bks_mask = .True.
 !call gstore%fill_bks_mask_with_pp(mband, nkibz, nsppol, bks_mask)

 !if (dtset%userie == 124) then
 !  ! Uncomment this line to have all states on each MPI rank.
 !  bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
 !end if

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

#if 0
 ! TODO: This part will be reintegrated afterwards if needed.
 ! ============================
 ! Compute vnk matrix elements
 ! ============================
 ! Diagonal elements of velocity operator in cartesian coordinates for all states in gwpt_nk.

 ABI_CALLOC(vk_cart_ibz, (3, gwpt%max_nbcalc, gwpt%nkcalc, nsppol))
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
       vk_cart_ibz(:, ib_k, ikcalc, spin) = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cg_work, cwaveprj0)
     end do

   end do
 end do
 call xmpi_sum(vk_cart_ibz, comm, ierr)

 ! Write v_nk to disk.
 !if (my_rank == master) then
 !  NCF_CHECK(nf90_put_var(gwpt%ncid, nctk_idname(gwpt%ncid, "vk_cart_ibz"), vk_cart_ibz))
 !end if
 ABI_FREE(vk_cart_ibz)

 call ddkop%free()
 call cwtime_report(" Velocities", cpu_ks, wall_ks, gflops_ks)
#endif

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

 ! Activate parallelism over perturbations at the level of the DVDB
 call gstore%set_perts_distrib(cryst, dvdb, my_npert)
 !call gstore%set_perts_distrib(crys, drhodb, my_npert)

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

 ! Distribute DFPT potentials (IBZ q-points) inside qpt_comm.
 ! Note that we distribute IBZ instead of the full BZ or the IBZ_k inside the loop over ikcalc.
 ! This means that the load won't be equally distributed but memory will scale with qpt_comm%nproc.
 ! To reduce load imbalance, we sort the qibz points by norm and use cyclic distribution inside qpt_comm

 ! itreat_qibz(nqibz)
 ! Table used to distribute potentials over q-points in the IBZ.
 ! The loop over qpts in the IBZ(k) is MPI distributed inside qpt_comm according to this table.
 ! 0 if this IBZ point is not treated by this proc.
 ! 1 if this IBZ is treated.

 ! TODO: Recheck this part
 ABI_ICALLOC(itreat_qibz, (gstore%nqibz))
 itreat_qibz = 1
 !call sort_rpts(gstore%nqibz, gstore%qibz, cryst%gmet, iperm)
 !do ii=1,gstore%nqibz
 !  iq_ibz = iperm(ii)
 !  do my_spin=1,gstore%my_nspins
 !    gqk => gstore%gqk(my_spin)
 !    if (mod(ii, gqk%qpt_comm%nproc) == gqk%qpt_comm%me) itreat_qibz(iq_ibz) = 1
 !  end do
 !end do
 !ABI_FREE(iperm)

 call wrtout(std_out, sjoin("P Number of q-points in the IBZ treated by this proc: " ,itoa(count(itreat_qibz == 1))))

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
   call dvdb%ftqcache_build(nfftf, ngfftf, gstore%nqibz, gstore%qibz, dtset%dvdb_qcache_mb, qselect, itreat_qibz, comm)

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
     if (itreat_qibz(iq_ibz) == 0) cycle
     db_iqpt = qibz2dvdb(iq_ibz)
     ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(gstore%qibz(:, iq_ibz)), "in the DVDB file."))
     itreatq_dvdb(db_iqpt) = 1
   end do
   call dvdb%qcache_read(nfftf, ngfftf, dtset%dvdb_qcache_mb, qselect, itreatq_dvdb, comm)
   ABI_FREE(itreatq_dvdb)
 end if

 ABI_FREE(itreat_qibz)
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

     if (done_qbz_spin(iq_bz, spin) == 1) cycle
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
     !print *, "my_npert", my_npert
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
       !print *, "my_ik", my_ik, " of my_nk:", gqk%my_nk

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
       ! TODO: Should oredere nbz in shells so that one can reduce the memory required to store W(pp)
       ! if we activate pp-parallelism.
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
                                  cg_kmp, cwaveprj0, cg1_kqmp, cwaveprj, msg, ierr) !, &
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

 ! This to make sure that the parallel output of GSTORE is completed
 call xmpi_barrier(comm)
 call cwtime_report(" gwpt_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

end subroutine gwpt_run
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

!!  subroutine gwpt_print(gwpt, dtset, unt)
!!
!!  !Arguments ------------------------------------
!!   integer,intent(in) :: unt
!!   type(dataset_type),intent(in) :: dtset
!!   class(gwpt_t),intent(in) :: gwpt
!!
!!  !Local variables-------------------------------
!!   integer :: ikc, is
!!   !character(len=5000) :: msg
!!
!!  ! *************************************************************************
!!
!!   if (unt == dev_null) return
!!
!!   ! Write dimensions
!!   write(unt,"(/,a)")sjoin(" Number of bands in gwpt sum:", itoa(gwpt%nbsum))
!!   write(unt,"(a)")sjoin(" From bsum_start:", itoa(gwpt%bsum_start), "to bsum_stop:", itoa(gwpt%bsum_stop))
!!   if (dtset%eph_stern /= 0) then
!!     write(unt, "(a)")" Treating high-energy bands with Sternheimer and static gwpt-energy."
!!     write(unt, "(a, es16.6, a, i0)")" Tolwfr:", dtset%tolwfr, ", nline: ", dtset%nline
!!   end if
!!   write(unt,"(a)")sjoin(" Symsigma: ",itoa(gwpt%symsigma), "Timrev:", itoa(gwpt%timrev))
!!   write(unt,"(a)")sjoin(" Imaginary shift in the denominator (zcut): ", ftoa(aimag(gwpt%ieta) * Ha_eV, fmt="f5.3"), "[eV]")
!!   write(unt,"(a)")sjoin(" Ab-initio q-mesh from DDB file:", ltoa(dtset%ddb_ngqpt))
!!   write(unt,"(a)")sjoin(" Q-mesh used for gwpt-energy integration [ngqpt]:", ltoa(gwpt%ngqpt))
!!   write(unt,"(a)")sjoin(" Number of q-points in the IBZ:", itoa(gwpt%nqibz))
!!   write(unt,"(a)")sjoin(" asr:", itoa(dtset%asr), "chneut:", itoa(dtset%chneut))
!!   write(unt,"(a)")sjoin(" dipdip:", itoa(dtset%dipdip), "symdynmat:", itoa(dtset%symdynmat))
!!
!!   write(unt,"(a, i0)")" Number of k-points for gwpt-energy corrections: ", gwpt%nkcalc
!!   if (any(abs(dtset%sigma_erange) /= zero)) then
!!     write(unt, "(a, 2(f6.3, 1x), a)")" sigma_erange: ", dtset%sigma_erange(:) * Ha_eV, " (eV)"
!!   end if
!!   write(unt,"(a)")" List of k-points for gwpt-energy corrections:"
!!   do ikc=1,gwpt%nkcalc
!!     if (ikc > 10) then
!!       write(unt, "(2a)")" nkcalc > 10. Stop printing more k-point information.",ch10
!!       exit
!!     end if
!!     do is=1,gwpt%nsppol
!!       if (gwpt%nsppol == 2) write(unt,"(a,i1,a)")" For spin: ",is, ", ikcalc, spin, kpt, bstart, bstop"
!!       write(unt, "(2(i4,2x),a,2(i4,1x))") &
!!         ikc, is, trim(ktoa(gwpt%kcalc(:,ikc))), gwpt%bstart_ks(ikc,is), gwpt%bstart_ks(ikc,is) + gwpt%nbcalc_ks(ikc,is) - 1
!!       end do
!!   end do
!!
!!   write(unt, "(/,a)")" === MPI parallelism ==="
!!   write(unt, "(2(a,i0))")"P Allocating and summing bands from my_bsum_start: ", gwpt%my_bsum_start, &
!!       " up to my_bsum_stop: ", gwpt%my_bsum_stop
!!   write(unt, "(a,i0)")"P Number of CPUs for parallelism over perturbations: ", gwpt%pert_comm%nproc
!!   write(unt, "(a,i0)")"P Number of perturbations treated by this CPU: ", gwpt%my_npert
!!   write(unt, "(a,i0)")"P Number of CPUs for parallelism over q-points: ", gwpt%qpt_comm%nproc
!!   !write(unt, "(2(a,i0))")"P Number of q-points in the IBZ treated by this proc: " , &
!!   !    count(gwpt%itreat_qibz == 1), " of ", gwpt%nqibz
!!   write(unt, "(a,i0)")"P Number of CPUs for parallelism over bands: ", gwpt%bsum_comm%nproc
!!   write(unt, "(a,i0)")"P Number of CPUs for parallelism over spins: ", gwpt%spin_comm%nproc
!!   write(unt, "(a,i0)")"P Number of CPUs for parallelism over k-points: ", gwpt%kcalc_comm%nproc
!!   write(unt, "(2(a,i0),/)")"P Number of k-point in gwpt_nk treated by this proc: ", gwpt%my_nkcalc, " of ", gwpt%nkcalc
!!
!!  end subroutine gwpt_print
!!***

end module m_gwpt
!!***

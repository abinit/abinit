!!****m* ABINIT/m_gwpt
!! NAME
!!  m_gwpt
!!
!! FUNCTION
!!  Compute the electron-phonon matrix elements with the GWPT formalism
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2025 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt.
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwpt

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 USE_MPI
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_hide_blas
 use m_copy
 use m_ifc
 use m_wfk
 use m_ddb
 use m_ddk
 use m_dvdb
 use m_fft
 use m_hamiltonian
 use m_pawcprj
 use m_wfd
 use m_krank
 use m_sort
 use m_hdr
 use m_sigtk
 use m_ephtk
 use netcdf
 use m_nctk
 use m_dtset
 use m_dtfil
 use m_clib
 use m_mkffnl
 use m_screen
 use m_xcdata

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : pseudopotential_type
 use m_gwdefs,         only : GW_Q0_DEFAULT
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, print_arr, inrange
 use m_io_tools,       only : iomode_from_fname, file_exists, is_open, open_file, flush_unit
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex, print_ngfft
 use m_cgtk,           only : cgtk_rotate, cgtk_change_gsphere
 use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map
 use m_kg,             only : getph
 use m_bz_mesh,        only : isamek, kmesh_t
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
 use m_gstore,         only : gstore_t, gqk_t, gstore_check_restart, GSTORE_GMODE_ATOM, GSTORE_GMODE_PHONON
 use m_rhotoxc,        only : rhotoxc
 use m_drivexc,        only : check_kxc
 use m_occ,            only : get_fact_spin_tol_empty
 use m_ppmodel,        only : PPM_HYBERTSEN_LOUIE, PPM_GODBY_NEEDS
 use m_ebands,         only : ebands_t
 use m_pstat,          only : pstat_proc
 use m_ppmodel,       only : ppmodel_t

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 public :: gwpt_run  ! Main entry point to compute GWPT e-ph matrix elements

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
!! dvdb=Database with the DFPT SCF potentials.
!! drhovdb=Database with the DFPT SCF densities.
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
!!     The default value of gwpc is single-precision.
!!     We use the following conventions for the buffers used to store the wavefunctions:
!!
!!       cg_kq, cr_kq
!!       cg1_kqmp, cr1_kqmp
!!
!! OUTPUT
!!
!! SOURCE

subroutine gwpt_run(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, dvdb, drhodb, ifc, wfk_hdr, &
                    pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),target,intent(in) :: ebands
 type(dvdb_t),intent(inout) :: dvdb, drhodb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(hdr_type),intent(in) :: wfk_hdr
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18), ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: LOG_MODQ = 1, LOG_MODK = 4, LOG_MODP = 4
 integer,parameter :: tim_getgh1c1 = 1, berryopt0 = 0, istw1 = 1, ider0 = 0, idir0 = 0, istwfk1 = 1
 integer,parameter :: useylmgr = 0, useylmgr1 = 0, master = 0, ndat1 = 1, with_cplex0 = 0, n3xccc0 = 0
 integer,parameter :: igscq0 = 0, icgq0 = 0, usedcwavef0 = 0, nbdbuf0 = 0, quit0 = 0, cplex1 = 1, pawread0 = 0
 integer :: band, band_me, nband_me, stern_comm, nkpt, my_rank, nsppol, iq_ibz, iq_bz, my_npert
 integer :: cplex,drho_cplex,nkxc,nk3xc,option,usexcnhat,db_iqpt,natom,natom3,ipc,nspinor,nprocs !, gsum_master
 integer :: ib_sum, ii, ib, u1_band !,u1c_ib_k,  jj, iw !ib_kq, band_ks, ib_k, ibsum_kq, u1_master, ip
 integer :: my_is, spin, idir,ipert, ig, max_npw_xc, min_npw_xc, npw_x, npw_c, nw_nk, nw_mkq
 integer :: my_pp_start_spin(dtset%nsppol), my_pp_stop_spin(dtset%nsppol), my_npp(dtset%nsppol)
 integer :: isym_q, trev_q
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_k_ibz, istwf_k_ibz
 integer :: ikq_ibz, isym_kq, trev_kq, npw_kq, istwf_kq,  npw_kq_ibz, istwf_kq_ibz
 integer :: ikmp_ibz, isym_kmp, trev_kmp, npw_kmp, istwf_kmp, npw_kmp_ibz, istwf_kmp_ibz
 integer :: ikqmp_ibz, isym_kqmp, trev_kqmp, npw_kqmp, istwf_kqmp, npw_kqmp_ibz, istwf_kqmp_ibz, mpw,ierr,nqbz,ncerr !,spad
 integer :: n1,n2,n3,n4,n5,n6,nspden, mqmem, m_kq, n_k, restart, root_ncid, spin_ncid, usecprj !,sij_opt
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg_k,nkpg_kq,nkpg_kqmp,nkpg_kmp,imyp, cnt, nvloc, iw_nk, iw_mkq, ndone, nmiss
 integer :: my_ipp, ipp_bz, ipp_ibz, isym_pp, itim_pp, comm_rpt, nqlwl
 integer :: qptopt, my_iq, my_ik, qbuf_size, iqbuf_cnt, nb ! nelem,
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_qq, wall_qq, gflops_qq, cpu_kk, wall_kk, gflops_kk, cpu_pp, wall_pp, gflops_pp
 !real(dp) :: cpu_all, wall_all, gflops_all, cpu_qq, wall_qq, gflops_qq, cpu_kk, wall_kk, gflops_kk
 !real(dp) :: drude_plsmf,my_plsmf
 real(dp) :: fact_spin, theta_mu_minus_e0i, tol_empty, tol_empty_in !, e_mkq, e_nk ! e0i,
 real(dp),ABI_CONTIGUOUS pointer :: qp_ene(:,:,:), qp_occ(:,:,:)
 real(dp) :: ecut,weight_q,bigexc,bigsxc,vxcavg ! ediff, eshift, q0rad, bz_vol,
 logical :: isirr_k, isirr_kq, isirr_kmp, isirr_kqmp, qq_is_gamma, pp_is_gamma, isirr_q
 logical :: stern_use_cache, stern_has_band_para, use_ftinterp
 logical :: print_time_qq, print_time_kk, print_time_pp, non_magnetic_xc, need_x_kmp, need_x_kqmp
 complex(dpc) :: ieta
 type(wfd_t) :: wfd
 type(gs_hamiltonian_type) :: gs_ham_kqmp, gs_ham_kmp
 type(rf_hamiltonian_type) :: rf_ham_kqmp, rf_ham_kmp
 type(ddkop_t) :: ddkop
 type(crystal_t) :: pot_cryst, den_cryst
 type(hdr_type) :: pot_hdr, den_hdr
 type(stern_t) :: stern_kmp, stern_kqmp
 type(kmesh_t) :: pp_mesh, kmesh
 type(gsphere_t) :: gsph_c
 type(gsphere_t),target :: gsph_x
 type(hscr_t) :: hscr
 type(vcoul_t) :: vcp
 type(screen_t),target :: screen
 type(screen_info_t) :: w_info
 type(gstore_t),target :: gstore
 type(gqk_t),pointer :: gqk
 type(xcdata_type) :: xcdata
 !type(ppmodel_t) :: PPm
 character(len=fnlen) :: screen_filepath, gstore_filepath
 character(len=5000) :: msg, qkp_string
!arrays
 integer :: nbsum, my_bsum_start(dtset%nsppol), my_bsum_stop(dtset%nsppol), my_nbsum(dtset%nsppol)
 integer :: g0_k(3), g0_kq(3), g0_kmp(3), g0_kqmp(3)
 integer :: mapl_k(6), mapl_kq(6), mapl_kqmp(6), mapl_kmp(6), mapc_qq(6), mapc_qq2dvdb(6)
 integer :: units(2), work_ngfft(18), gmax(3)
 integer(i1b),allocatable :: itreatq_dvdb(:)
 integer,allocatable :: my_pp_inds(:) ! bands_treated_now(:),
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), kg_kmp(:,:), kg_kqmp(:,:)
 integer,allocatable :: gbound_k(:,:), gbound_kq(:,:), gbound_kmp(:,:), gbound_kqmp(:,:), gbound_c(:,:), gbound_x(:,:)
 integer,allocatable :: nband(:,:), wfd_istwfk(:), count_bk(:,:)
 integer,allocatable :: qibz2dvdb(:) !, displs(:), recvcounts(:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:) !, my_iqibz_inds(:) !, iperm(:)
 integer(i1b),allocatable :: itreat_qibz(:)
 integer, ABI_CONTIGUOUS pointer :: kg_c(:,:), kg_x(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3), kqmp(3), kmp(3), pp(3), kmp_ibz(3), kqmp_ibz(3)
 real(dp) :: phfr_qq(3*cryst%natom), qq_ibz(3), qpt(3)
 !real(dp) :: tsec(2) ! vk(3), vkq(3),
 real(dp) :: eig0nk !,eig0mkq !eig0mk,
 complex(gwpc) :: ctmp_gwpc, xdot_tmp
!arrays
 real(dp) :: fermie1_idir_ipert(3,cryst%natom), ylmgr_dum(1,1,1), dum_nhat(0), dum_xccc3d(0)
 real(dp),allocatable :: qlwl(:,:), vk_cart_ibz(:,:,:), displ_cart_qq(:,:,:,:),displ_red_qq(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg_k(:,:),kpg_kq(:,:),kpg_kmp(:,:),kpg_kqmp(:,:), dkinpw(:)
 real(dp),allocatable :: ffnl_kmp(:,:,:,:),ffnl_kqmp(:,:,:,:)
 real(dp),allocatable :: ph3d_kmp(:,:,:), ph3d1_kqmp(:,:,:), ph3d_kqmp(:,:,:), ph3d1_kmp(:,:,:)
 real(dp),allocatable, target :: vxc1_qq(:,:,:,:)
 real(dp),allocatable :: gsig_atm(:,:,:,:),gsig_nu(:,:,:,:), gxc_atm(:,:,:,:), gxc_nu(:,:,:,:), gks_atm(:,:,:,:), gks_nu(:,:,:,:)
 real(dp),allocatable :: cg_work(:,:), ug_k(:,:), ug_kq(:,:) !,kets_k(:,:,:),h1kets_kq(:,:,:,:)
 real(dp),allocatable :: ph1d(:,:), vlocal(:,:,:,:), vlocal1_qq(:,:,:,:,:), v1scf_qq(:,:,:,:), vlocal1_mqq(:,:,:,:,:), v1scf_mq(:,:,:,:)
 real(dp),allocatable :: ylm_k(:,:), ylm_kq(:,:), ylm_kmp(:,:), ylm_kqmp(:,:)
 real(dp),allocatable :: ylmgr_kq(:,:,:), ylmgr_kmp(:,:,:), ylmgr_kqmp(:,:,:)
 real(dp),allocatable :: vtrial(:,:), work(:,:,:,:), rhor(:,:), vxc(:,:), kxc(:,:)
 real(dp),allocatable :: omegame0i_nk(:), omegame0i_mkq(:), omegas_nk(:), omegas_mkq(:)
 real(dp),allocatable :: my_gbuf(:,:,:,:,:,:) !, buf_wqnu(:,:), buf_eigvec_cart(:,:,:,:,:)
 real(dp),allocatable :: cg_kmp(:,:), cg_kqmp(:,:), cg1_kqmp(:,:), cg1_kmp(:,:), full_cg1_kqmp(:,:), full_cg1_kmp(:,:)
 complex(dp), contiguous, pointer :: cvxc1_qq_ptr(:,:,:)
 complex(gwpc),allocatable :: ur_kmp(:), ur_kqmp(:), cwork_ur(:), rhotwg_c(:), rhotwg_x(:), vc_sqrt_gx(:) !, rhotwgp(:)
 complex(gwpc),allocatable :: full_ur1_kqmp(:), full_ur1_kmp(:), sigcme_nk(:), sigcme_mkq(:), ur_nk(:,:), ur_mkq(:,:)
 complex(gwpc),allocatable :: vec_gwc_nk(:,:,:), vec_gwc_mkq(:,:,:), vec_gx_nk(:,:), vec_gx_mkq(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:), cwaveprj(:,:)
 type(pawrhoij_type),allocatable :: pot_pawrhoij(:), den_pawrhoij(:)
!************************************************************************

 ! Problems to be addressed:
 !
 ! 1)
 ! Sigma is usually split into Sigma_c(w) and Sigma_x where Sigma_x is the static Fock operator
 ! evaluated with KS orbitals. The advantage of such partitioning is that Sigma_x = iGv
 ! can be computed by summing over occupied states only. Sigma_x requires more G-vectors to converge
 ! as the bare Coulomb interaction goes as 1/|q+G|^2 that is not integrable in 3D but this "expensive"
 ! operations are needed only inside a sum over bands that is restricted to occupied states.
 ! On the other hand, Sigma_c(w) is way more expensive as we have to sum a large number of empty states
 ! while taking the w-dependence of the screening into account.
 ! Fortunately, all the operations can be restricted to a small G-sphere of kinetic energy ecuteps that can be handled
 ! with a coarser FFT mesh.
 ! Another distinct advantage of such splitting is that one can handle the divergence in vc(q,g) for |q+g| --> 0
 ! using well know techniques from GW and the anisotropic behavior of eps-1(q) for q --> 0 in low-dimensional systems.
 ! The disavantage is that one needs to compute the GWPT e-ph matrix in two steps, first Sigma_c and then Sigma_x,
 ! so certain operations such as the k-point mapping, and the computation of the form factors are performed twice
 ! Note, however, that MG believes that Sigma_x is a much better approximation than v_xc when one is interested
 ! in the e-ph matrix elements connecting low-energy states such as band edges to high-energy states.
 !
 ! 2)
 ! We need to solve the NSCF Sternheimer for q and -q. In principle one can solve the equation only at q
 ! and then use spatial inversion or TR to get the solution at -q but this requires solving the Sternheimer
 ! for all the pp wavevectors in the BZ (or better in the IBZ_{q,k,alpha). The use of symmetries is rendered complicated
 ! by the parallelism over pp but perhaps one can precompute \Delta psi with all MPI procs and write the results to temporary files.

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm); units = [std_out, ab_out]
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 ! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 ecut = dtset%ecut ! dtset%dilatmx
 ieta = +j_dpc * dtset%zcut

 ! Set tolerance used to decide if a band is empty
 ! and normalization of theta_mu_minus_esum. If nsppol == 2, qp_occ $\in [0,1]$
 tol_empty_in = 0.01
 call get_fact_spin_tol_empty(nsppol, nspinor, tol_empty_in, fact_spin, tol_empty)

 qp_ene => ebands%eig; qp_occ => ebands%occ

 ! Check if a previous GSTORE.nc file is present to restart the calculation if dtset%eph_restart == 1,
 ! and use done_qbz_spin mask to cycle the loops below if restart /= 0.
 gstore_filepath = strcat(dtfil%filnam_ds(4), "_GSTORE.nc")
 call gstore_check_restart(gstore_filepath, dtset, nqbz, done_qbz_spin, restart, comm)

 if (restart == 0) then
   ! Build new gstore object from dtset input variables.
   call gstore%init(gstore_filepath, dtset, dtfil, wfk_hdr, cryst, ebands, ifc, comm)
   ABI_REMALLOC(done_qbz_spin, (gstore%nqbz, nsppol))
   done_qbz_spin = 0
 else
   ! Init gstore from pre-existent file.
   call gstore%from_ncpath(gstore_filepath, with_cplex0, dtset, cryst, ebands, ifc, comm)
 end if

 ! Open GSTORE.nc file and go to data mode.
 NCF_CHECK(nctk_open_modify(root_ncid, gstore%path, comm))
 NCF_CHECK(nctk_set_datamode(root_ncid))

 !if (my_rank == master) call gwpt%print(dtset, std_out)

 call gstore%get_missing_qbz_spin(done_qbz_spin, ndone, nmiss)
 !call wrtout(units, sjoin("- Number of q-points/spin completed:", itoa(count(done_qbz_spin == 1)), "/", itoa(sigma%nkcalc)))

 ! ================
 ! HANDLE SCREENING
 ! ================
 ! Init gpsh_c for the correlated part.
 nqlwl = 0; screen_filepath = ABI_NOFILE
 if (dtset%getscr /= 0 .or. dtset%irdscr /= 0 .or. dtset%getscr_filepath /= ABI_NOFILE) screen_filepath = dtfil%fnameabi_scr

 call kmesh%init(cryst, wfk_hdr%nkpt, wfk_hdr%kptns, dtset%kptopt)

 w_info%use_ppm = dtset%ppmodel
 if (screen_filepath /= ABI_NOFILE) then
   ! Read g-sphere and pp_mesh from SCR file.
   call get_hscr_qmesh_gsph(screen_filepath, dtset, cryst, hscr, pp_mesh, gsph_c, qlwl, comm)
   call hscr%free()
   nqlwl = size(qlwl, dim=2)
   w_info%use_mdf = MDL_NONE
 else
   ! Init pp_mesh from the K-mesh reported in the WFK file.
   call pp_mesh%find_qmesh(cryst, kmesh)
   ! The G-sphere for W and Sigma_c is initialized from ecuteps.
   call gsph_c%init(cryst, 0, ecut=dtset%ecuteps)
   dtset%npweps = gsph_c%ng
   w_info%use_mdf = MDL_BECHSTEDT
   w_info%eps_inf = dtset%mdf_epsinf
   !w_info%use_ppm = PPM_HYBERTSEN_LOUIE
   ABI_CHECK(w_info%eps_inf > zero, "Model dielectric function requires the specification of mdf_epsinf")
   ABI_CHECK(w_info%use_ppm /= PPM_GODBY_NEEDS, "Godby needs PPM is not compatible with model dielectric function")
 end if

#if 0
 call test_charge(nfftf,ks_ebands%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
                  Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)

 my_plsmf = drude_plsmf; if (Dtset%ppmfrq > tol6) my_plsmf = Dtset%ppmfrq
 call PPm%init(epsm1%mqmem, epsm1%nqibz, epsm1%npwe, Sigp%ppmodel, my_plsmf, Dtset%gw_invalid_freq)
 call ppm%new_setup(iq_ibz, Cryst, Qmesh, npwe, nomega, omega, epsm1_ggw, nfftf, gvec, ngfftf, rhor_tot)
 capp ppm%free()
#endif

 call pstat_proc%print(_PSTAT_ARGS_)

 if (nqlwl == 0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i0,a,3f9.6)')&
     "The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
     "Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)
   ABI_COMMENT(msg)
 end if

 ! Init g-sphere for the exchange part from ecutsigx
 call gsph_c%extend(cryst, dtset%ecutsigx, gsph_x)

 ! TODO:
 ! Here we sort the pp_mesh by stars so that we can split the pp wavevectors in blocks and therefore
 ! reduce the number of wavevectors in the IBZ that must be stored in memory.
 !call pp_mesh%pack_by_stars()

 ! Distribute the sum over pp wavevectors inside pp_sum_comm using block distribution.
 my_pp_start_spin = -1; my_pp_stop_spin = 0
 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)
   call xmpi_split_block(pp_mesh%nbz, gqk%pp_sum_comm%value, my_npp(spin), my_pp_inds)
   if (my_npp(spin) > 0) then
     my_pp_start_spin(spin) = my_pp_inds(1); my_pp_stop_spin(spin) = my_pp_inds(my_npp(spin))
   end if
   ABI_SFREE(my_pp_inds)
 end do ! my_is

 ! Initialize Coulomb term on the IBZ of the pp_mesh. Use the largest G-sphere.
 npw_x = gsph_x%ng; npw_c = gsph_c%ng
 max_npw_xc = max(npw_x, npw_c)
 min_npw_xc = min(npw_x, npw_c)
 if (gsph_x%ng >= gsph_c%ng) then
   call vcp%init(gsph_x, cryst, pp_mesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_x%ng, &
                 nqlwl, qlwl, comm)
 else
   call vcp%init(gsph_c, cryst, pp_mesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_c%ng, &
                 nqlwl, qlwl, comm)
 end if
 if (my_rank == master) call vcp%print([std_out], prtvol=dtset%prtvol)
 call kmesh%free()

 call gsph_x%print(unit=std_out, prtvol=dtset%prtvol)
 call gsph_c%print(unit=std_out, prtvol=dtset%prtvol)
 ABI_CHECK_IGE(npw_x, 1, "npw_x <= 1")
 ABI_CHECK_IGE(npw_c, 1, "npw_c <= 1")

 ! Initialize the wave function descriptor.
 ! Each node has all k-points and spins and bands between my_bsum_start and my_bsum_stop
 ! TODO: One can exploit qq, kk and pp parallelism to find the wavevectors in the IBZ
 ! that will be needed in the loops below and allocate only these wavevectors so that memory scales.

 nbsum = dtset%mband
 my_bsum_start = 1; my_bsum_stop = nbsum; my_nbsum = my_bsum_stop - my_bsum_start + 1
 stern_has_band_para = .False.
 ! FIXME: This term is needed in metals.
 fermie1_idir_ipert = zero

 nkpt = wfk_hdr%nkpt
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkpt, nsppol))

 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.
 !bks_mask(my_bsum_start:my_bsum_stop,:,:) = .True.
 ! Distribute wavefunctions according to the set of kk, qq and pp wavevectors treated by this MPI proc.

 ! Also, compute mpw and gmax including the additional pp
 ! This is the maximum number of PWs for all possible k+q treated.
 call gstore%fill_bks_mask_pp_mesh(ecut, dtset%mband, nkpt, nsppol, my_pp_start_spin, my_pp_stop_spin, pp_mesh, &
                                   my_bsum_start, my_bsum_stop, bks_mask, mpw, gmax)

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, should also consider Gamma-only and istwfk
 gmax = 2*gmax + 1

 if (dtset%userie == 124) then
   ! Uncomment this line to have all states on each MPI rank.
   bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
 end if

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd%init(cryst, pawtab, psps, keep_ur, dtset%mband, nband, nkpt, nsppol, bks_mask,&
               nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call wfd%print([std_out], header="Wavefunctions for GWPT calculation.")
 call pstat_proc%print(_PSTAT_ARGS_)

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! FFT meshes from input file, not necessarily equal to the ones found in the external files.
 ! NB: ur arrays are always allocated with nfft and not with product(ngfft(4:6)).
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 call print_ngfft([std_out], ngfft, header="FFT mesh")

 ! Set the FFT mesh
 call wfd%change_ngfft(cryst, psps, ngfft)

 usecprj = dtset%usepaw
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))
 ABI_MALLOC(cwaveprj, (natom, nspinor*usecprj))
 ABI_MALLOC(displ_cart_qq, (2, 3, natom, natom3))
 ABI_MALLOC(displ_red_qq, (2, 3, natom, natom3))
 ABI_MALLOC(gbound_k, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_kq, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_kmp, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_kqmp, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_c, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_x, (2*mgfft+8, 2))
 ABI_MALLOC(cg_work, (2, mpw*nspinor))
 ABI_MALLOC(full_ur1_kqmp, (nfft*nspinor))
 ABI_MALLOC(full_ur1_kmp, (nfft*nspinor))

 ! ============================
 ! Compute vnk matrix elements
 ! ============================
 ! Diagonal elements of velocity operator in cartesian coordinates for all kk in the IBZ.
 ! Use ndone to understand if velocities have been already compured in a previous run.

 if (gstore%with_vk /= 0 .and. ndone == 0) then
   call wrtout(std_out, " computing and writing velocity operator matrix elements in the ibz")
   call wrtout(std_out, " note that not all the k-points in the ibz are computed when kfilter is activated!")
   call cwtime(cpu_kk, wall_kk, gflops_kk, "start")

   ! On disk, we have:
   !    nctkarr_t("vk_cart_ibz", "dp", "three, nb, gstore_nkibz"))
   !    nctkarr_t("vkmat_cart_ibz", "dp", "two, three, nb, nb, gstore_nkibz")))

   ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

   do my_is=1,gstore%my_nspins
     spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is); nb = gqk%nb

     ! Be careful as wavefunctions might be replicated.
     ! Use count_bk to count how many states has been computed
     ! in parallel in order to rescale the results.
     if (gstore%with_vk == 1) then
       ABI_CALLOC(vk_cart_ibz, (3, gqk%nb, gstore%nkibz))
       ABI_ICALLOC(count_bk, (nb, gstore%nkibz))
     else
       ABI_ERROR("gstore%with_vk 2 not implemented")
     end if

     do my_ik=1,gqk%my_nk
       ! The k-point and the symmetries relating the BZ k-point to the IBZ.
       kk = gqk%my_kpts(:, my_ik)
       ik_ibz = gqk%my_k2ibz(1, my_ik) ; isym_k = gqk%my_k2ibz(2, my_ik)
       trev_k = gqk%my_k2ibz(6, my_ik); g0_k = gqk%my_k2ibz(3:5,my_ik)
       isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
       if (.not. isirr_k) cycle

       npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
       call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

       do band=gqk%bstart, gqk%bstop
         call wfd%copy_cg(band, ik_ibz, spin, cg_work)
         eig0nk = ebands%eig(band, ik_ibz, spin)
         ib = band - gqk%bstart + 1
         vk_cart_ibz(:, ib, ik_ibz) = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cg_work, cwaveprj0)
         count_bk(ib, ik_ibz) = count_bk(ib, ik_ibz) + 1
       end do
     end do ! my_ik

     call xmpi_sum(count_bk, gqk%comm%value, ierr)
     call xmpi_sum(vk_cart_ibz, gqk%comm%value, ierr)

     do ik_ibz=1, gstore%nkibz
       do band=gqk%bstart, gqk%bstop
         ib = band - gqk%bstart + 1
         if (count_bk(ib, ik_ibz) == 0) cycle
         do ii=1,3
           vk_cart_ibz(ii,ib,ik_ibz) = vk_cart_ibz(ii,ib,ik_ibz) / count_bk(ib, ik_ibz)
         end do
       end do
     end do

     ! Write v_nk to disk.
     !if (gqk%comm%me == master) then
       NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))
       NCF_CHECK(nf90_put_var(spin_ncid, spin_vid("vk_cart_ibz"), vk_cart_ibz))
     !end if

     ABI_FREE(vk_cart_ibz)
     ABI_FREE(count_bk)
   end do ! my_is

   call ddkop%free()
   call cwtime_report(sjoin(" Computation of v_k group velocities with with_vk:", itoa(gstore%with_vk)), cpu_kk, wall_kk, gflops_kk)
 end if ! ndone /= 0

 ! Radius of sphere with volume equivalent to the micro zone.
 !q0rad = two_pi * (three / (four_pi * cryst%ucvol * gstore%nqbz)) ** third
 !bz_vol = two_pi**3 / cryst%ucvol

 ! Open the DVDB file with first-order potentials and drhodb with the first-order densities.
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 call drhodb%open_read(ngfftf, xmpi_comm_self)
 ABI_CHECK(dvdb%has_fields("pot1", msg), msg)
 ABI_CHECK(drhodb%has_fields("den1", msg), msg)

 ! Make sure that dvdb and drhodb have the same q-points
 ABI_CHECK_IEQ(dvdb%nqpt, drhodb%nqpt, "Different number of q-points in DVDB and DRHODB")
 ierr = 0
 do ii=1,dvdb%nqpt
   if (any(dvdb%qpts(:, ii) /= drhodb%qpts(:, ii))) then
     ierr = ierr + 1; call wrtout(std_out, sjoin(ktoa(dvdb%qpts(:, ii)), " /= ", ktoa(drhodb%qpts(:, ii))))
   end if
 end do
 ABI_CHECK(ierr == 0, "Found different q-points in DVDB and DRHODB. See messages above!")

 ! Check if the q-points are present in the DVDB
 !qptopt = dtset%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
 !call dvdb%need_ftinterp(nq_path, qpath%points, qptopt, qmap_symrec, use_ftinterp)

 !if (.not. use_ftinterp .and. dtset%eph_use_ftinterp /= 0) then
 !  ABI_WARNING("Enforcing FT interpolation for q-points even if it's not strictly needed.")
 !  use_ftinterp = .True.
 !end if

 !if (use_ftinterp) then
 !  call wrtout(units, " Cannot find all q-points in the DVDB --> Activating Fourier interpolation.")
 !  call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, xmpi_comm_self)
 !else
 !  call wrtout(units, " DVDB file contains all q-points along the path --> Reading DFPT potentials from file.")
 !end if

 ! Activate parallelism over perturbations at the level of the DVDB, my_npert is output
 call gstore%set_perts_distrib(cryst, dvdb, my_npert)
 call gstore%set_perts_distrib(cryst, drhodb, my_npert)
 !print *, "Treating my_npert", my_npert

 ! This part is taken from dfpt_vtorho
 !==== Initialize most of the Hamiltonian (and derivative) ====
 ! 1) Allocate all arrays and initialize quantities that do not depend on k and spin.
 ! 2) Perform the setup needed for the non-local factors:
 !
 ! Norm-conserving: Constant kleimann-Bylander energies are copied from psps to gs_ham_kqmp.
 ! PAW: Initialize the overlap coefficients and allocate the Dij coefficients.

 ! Get one-dimensional structure factor information on the coarse grid.
 ABI_MALLOC(ph1d, (2,3*(2*mgfft+1)*natom))
 call getph(cryst%atindx, natom, n1, n2, n3, ph1d, cryst%xred)

 call gs_ham_kqmp%init(psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 call gs_ham_kmp%init(psps, pawtab, nspinor, nsppol, nspden, natom, &
   dtset%typat, cryst%xred, nfft, mgfft, ngfft, cryst%rprimd, dtset%nloalg, &
   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab, mpi_spintab=mpi_enreg%my_isppoltab, &
   usecprj=usecprj, ph1d=ph1d, nucdipmom=dtset%nucdipmom, gpu_option=dtset%gpu_option)

 ! Allocate work space arrays.
 ! vtrial and vlocal are required for Sternheimer (H0). DFPT routines do not need it.
 ! Note nvloc in vlocal (we will select one/four spin components afterwards)
 nvloc = gs_ham_kqmp%nvloc
 ABI_CALLOC(vtrial, (nfftf, nspden))
 ABI_CALLOC(vlocal, (n4, n5, n6, nvloc))

 ! Read the GS potential (vtrial) from input POT file
 ! In principle one may store GS vtrial in the DVDB but getpot_filepath is simpler to implement.
 call wrtout(units, sjoin(" Reading KS GS potential for Sternheimer from: ", dtfil%filpotin))
 call read_rhor(dtfil%filpotin, cplex1, nspden, nfftf, ngfftf, pawread0, mpi_enreg, vtrial, pot_hdr, pot_pawrhoij, comm, &
                allow_interp=.True., want_varname="vtrial")
 pot_cryst = pot_hdr%get_crystal()
 if (cryst%compare(pot_cryst, header=" Comparing input crystal with POT crystal") /= 0) then
   ABI_ERROR("Crystal structure from WFK and POT do not agree! Check messages above!")
 end if
 call pot_cryst%free(); call pot_hdr%free()

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
   ! Use ddb_ngqpt q-mesh to compute the real-space representation of DFPT v1scf_qq potentials to prepare Fourier interpolation.
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
   call drhodb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, comm_rpt)
 else
 end if

 call dvdb%print(prtvol=dtset%prtvol)
 call drhodb%print(prtvol=dtset%prtvol)

 if (.not. use_ftinterp) then
   ! Need to translate itreat_qibz into itreatq_dvdb.
   ABI_ICALLOC(itreatq_dvdb, (dvdb%nqpt))
   do iq_ibz=1,gstore%nqibz
     if (itreat_qibz(iq_ibz) == 0) cycle
     db_iqpt = qibz2dvdb(iq_ibz)
     ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(gstore%qibz(:, iq_ibz)), "in the DVDB file."))
     itreatq_dvdb(db_iqpt) = 1
   end do
   ABI_FREE(itreatq_dvdb)
 end if

 ABI_FREE(itreat_qibz)
 ABI_FREE(qibz2dvdb)

 ! The GS density is needed for vxc1_qq and the model dielectric function
 ABI_CALLOC(rhor, (nfftf, nspden))
 call read_rhor(dtfil%fildensin, cplex1, nspden, nfftf, ngfftf, pawread0, mpi_enreg, rhor, den_hdr, den_pawrhoij, comm, &
                allow_interp=.True., want_varname="density")
 den_cryst = den_hdr%get_crystal()
 if (cryst%compare(den_cryst, header=" Comparing input crystal with DEN crystal") /= 0) then
   ABI_ERROR("Crystal structure from WFK and DEN do not agree! Check messages above!")
 end if
 call den_cryst%free(); call den_hdr%free()

 ! Init Wc. In-core or out-of-core solution?
 !mqmem = 0; if (dtset%gwmem /10 == 1) mqmem = pp_mesh%nibz
 mqmem = pp_mesh%nibz
 w_info%invalid_freq = dtset%gw_invalid_freq
 w_info%mat_type = MAT_INV_EPSILON
 w_info%wint_method = WINT_PPMODEL

 ! possible memory leak here
 !   === STACK OF LEN 3 ===
 !   {'A:botsq@m_ppmodel.F90:255': [<var=botsq, A@m_ppmodel.F90:255, addr=0x60000370d3c0, size_mb=0.000>],
 !    'A:eig@m_ppmodel.F90:257': [<var=eig, A@m_ppmodel.F90:257, addr=0x60000370d410, size_mb=0.000>],
 !    'A:otq@m_ppmodel.F90:256': [<var=otq, A@m_ppmodel.F90:256, addr=0x60000370d490, size_mb=0.000>]}

 call screen%init(w_info, cryst, pp_mesh, gsph_c, vcp, screen_filepath, mqmem, dtset%npweps, &
                  dtset%iomode, ngfftf, nfftf, nsppol, nspden, rhor, dtset%prtvol, comm)
 ABI_FREE(qlwl)
 !call screen%print(units)
 !call screen%ppm%print(units)

 ! Allocate g-vectors centered on k, k+q, k-p, and k+q-p
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))
 ABI_MALLOC(kg_kmp, (3, mpw))
 ABI_MALLOC(kg_kqmp, (3, mpw))
 ! Spherical Harmonics for useylm == 1.
 ! FIXME: These arrays should allocated with npw_k, npw_kq
 ! but should recheck the API used to symmetrized wavefunctions.
 ABI_MALLOC(ylm_k, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kq, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kmp, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylm_kqmp, (mpw, psps%mpsang**2 * psps%useylm))
 ABI_MALLOC(ylmgr_kq, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))
 ABI_MALLOC(ylmgr_kmp, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))
 ABI_MALLOC(ylmgr_kqmp, (mpw, 3, psps%mpsang**2 * psps%useylm * useylmgr1))
 ! GS wavefunctions
 ABI_MALLOC(ur_kmp, (nfft*nspinor))
 ABI_MALLOC(ur_kqmp, (nfft*nspinor))
 ABI_MALLOC(cwork_ur, (nfft*nspinor))
 ABI_MALLOC(cg_kmp, (2, mpw*nspinor))
 ABI_MALLOC(cg_kqmp, (2, mpw*nspinor))
 ! First order change (full term including active space)
 ABI_MALLOC(cg1_kqmp, (2, mpw*nspinor))
 ABI_MALLOC(cg1_kmp, (2, mpw*nspinor))

 !stern_use_cache = merge(.True., .False., dtset%eph_stern == 1)
 stern_use_cache = .False.

 if (my_rank == master) call gstore%print([std_out])
 call pstat_proc%print(_PSTAT_ARGS_)

 ! This parameter defines the size of the q-buffer used to store the g(k, q) e-ph matrix elements
 ! for all the k-point treated by this MPI rank.
 ! Increasing the buffer size increases the memory requirements
 ! but it leads to better performance as the number of IO operations is decreased.
 ! TODO: Should compute it on the basis of my_nkpt and my_nqpt
 qbuf_size = 1
 !qbuf_size = 16
 call wrtout(std_out, sjoin(" Begin computation of GWPT e-ph matrix elements with qbuf_size:", itoa(qbuf_size)), pre_newlines=1)

 ! A similar piece of code is used in m_respfn_driver.
 ! option 2 for xc and kxc (no paramagnetic part if xcdata%nspden=1). Note dum_xccc3d to ignore NLCC
 nkxc = 2*min(dtset%nspden,2)-1; if (dtset%xclevel==2) nkxc = 12*min(dtset%nspden,2)-5
 call xcdata_init(xcdata, dtset=dtset)
 non_magnetic_xc = (dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

 call check_kxc(dtset%ixc, dtset%optdriver)
 ABI_MALLOC(kxc, (nfft, nkxc))
 ABI_MALLOC(vxc, (nfft, dtset%nspden))

 nk3xc=1; option=2; usexcnhat=0
 call rhotoxc(bigexc, bigsxc, kxc, mpi_enreg, nfft, ngfft, &
              dum_nhat, 0, dum_nhat, 0, nkxc, nk3xc, non_magnetic_xc, n3xccc0, option, rhor, &
              cryst%rprimd, usexcnhat, vxc, vxcavg, dum_xccc3d, xcdata)
 call pstat_proc%print(_PSTAT_ARGS_)

 ! ===================================================
 ! Loop over MPI distributed spins in Sigma (gqk%comm)
 ! ===================================================

 do my_is=1,gstore%my_nspins
   spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is); my_npert = gqk%my_npert
   ABI_CHECK_IEQ(my_npert, gqk%my_npert, "my_npert")
   ABI_CHECK_IGEQ(nbsum, gqk%bstop, "nband must be greater than the max band in the e-ph matrix elements")

   NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

   ! TODO: Should introduce the possibility of specifying different nb states
   ! for the incoming and the intermediates states.
   nb = gqk%nb
   ABI_MALLOC(iq_buf, (2, qbuf_size))
   ABI_MALLOC(gsig_atm, (2, nb, nb, natom3))
   ABI_MALLOC(gsig_nu, (2, nb, nb, natom3))
   ABI_MALLOC(gks_atm, (2, nb, nb, natom3))
   ABI_MALLOC(gks_nu, (2, nb, nb, natom3))
   ABI_MALLOC(gxc_atm, (2, nb, nb, natom3))
   ABI_MALLOC(gxc_nu, (2, nb, nb, natom3))

   ABI_MALLOC_OR_DIE(ur_nk,  (nfft*nspinor, gqk%bstart:gqk%bstop), ierr)
   ABI_MALLOC_OR_DIE(ur_mkq, (nfft*nspinor, gqk%bstart:gqk%bstop), ierr)

   ! Inside the loops we compute gsig_nu(2, nb, nb, natom3)
   ABI_MALLOC_OR_DIE(my_gbuf, (gqk%cplex, nb, nb, natom3, gqk%my_nk, qbuf_size), ierr)

   ! Allocate memory to deal with frequencies in Sigma(w).
   nw_nk = 1 + (gqk%bstop - gqk%bstart + 1)
   nw_mkq = 1 + (gqk%bstop - gqk%bstart + 1)
   ABI_MALLOC(omegame0i_nk, (nw_nk))
   ABI_MALLOC(omegame0i_mkq, (nw_mkq))
   ABI_MALLOC(omegas_nk, (nw_nk))
   ABI_MALLOC(omegas_mkq, (nw_mkq))
   ABI_MALLOC(sigcme_nk, (nw_nk))
   ABI_MALLOC(sigcme_mkq, (nw_mkq))

   ! Correlated part
   ABI_MALLOC_OR_DIE(vec_gwc_nk, (npw_c*nspinor, nw_nk, gqk%bstart:gqk%bstop), ierr)
   ABI_MALLOC_OR_DIE(vec_gwc_mkq, (npw_c*nspinor, nw_mkq, gqk%bstart:gqk%bstop), ierr)
   ! Exchange part
   ABI_MALLOC_OR_DIE(vec_gx_nk, (npw_x*nspinor,  gqk%bstart:gqk%bstop), ierr)
   ABI_MALLOC_OR_DIE(vec_gx_mkq, (npw_x*nspinor, gqk%bstart:gqk%bstop), ierr)

   ! ============================================================
   ! Loop over MPI distributed q-points in Sigma_q (gqk%qpt_comm)
   ! ============================================================
   do my_iq=1,gqk%my_nq
     !if (my_iq <= 4) cycle
     call gqk%myqpt(my_iq, gstore, weight_q, qpt)
     qq_is_gamma = sum(qpt**2) < tol14

     ! Note symrec conventions here.
     iq_bz = gqk%my_q2bz(my_iq)

     ! Handle possible restart.
     if (done_qbz_spin(iq_bz, spin) == 1) then
       call wrtout(std_out, sjoin(" iq_bz:", itoa(iq_bz), ", spin: ", itoa(spin), " already computed --> skipping iteration"))
       cycle
     end if

     iq_ibz = gqk%my_q2ibz(1, my_iq) ; isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq) !; g0_q = gqk%my_q2ibz(3:5, my_iq)
     isirr_q = (isym_q == 1 .and. trev_q == 0)
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     !tsign_q = 1; if (trev_q == 1) tsign_q = -1
     qq_ibz = gstore%qibz(:, iq_ibz)
     mapc_qq = gqk%my_q2ibz(:, my_iq)

     print_time_qq = my_rank == 0 .and. (my_iq <= LOG_MODQ .or. mod(my_iq, LOG_MODQ) == 0)
     if (print_time_qq) then
       call cwtime(cpu_qq, wall_qq, gflops_qq, "start")
       call inds2str(0, sjoin(" Computing g^Sigma(k, q) for qpt:", ktoa(qpt)), my_iq, gqk%my_nq, gqk%glob_nq, msg)
       call wrtout(std_out, sjoin(msg, ", for spin:", itoa(spin)), pre_newlines=1)
     end if
     !print *, "iq_ibz:", iq_ibz, "qpt:", qpt, "qq_ibz:", qq_ibz

     ! FIXME: Here I should compute stuff in the IBZ and then rotate in order to fix the gauge in g
     ! Get phonon frequencies and eigenvectors for the corresponding q-point in the IBZ.
     call ifc%fourq(cryst, qpt, phfr_qq, displ_cart_qq, out_displ_red=displ_red_qq)

     if (isirr_q) then
       !displ_cart_qbz = displ_cart_qibz; displ_red_qbz = displ_red_qibz; pheigvec_qbz = pheigvec_qibz
     else
       !! Rotate phonon eigenvectors from q_ibz to q_bz.
       !! This part is needed to enforce the gauge in the ph eigenvectors, including e(-q) = e(q)^*
       !call pheigvec_rotate(cryst, qq_ibz, isym_q, trev_q, pheigvec_qibz, pheigvec_qbz, displ_cart_qbz, &
       !                     displ_red_qbz=displ_red_qbz)
     end if

     ! ==================================================
     ! Get DFPT potentials and densities for this q-point
     ! ==================================================
     ! After this branch we know `cplex` and we have allocated:
     !
     !   v1scf_qq(cplex, nfftf, nspden, my_npert))
     !   vxc1_qq(cplex, nfft, nspden, my_npert)
     !
     ! Important: vxc1_qq does not include the contribution due to the model core charge (if any).

     if (use_ftinterp) then
       ! Use Fourier interpolation to get DFPT potentials and DFPT densities for this qpt.
       call dvdb%get_ftqbz(qpt, cplex, nfftf, ngfftf, v1scf_qq, gqk%pert_comm%value)

       call drhodb%get_vxc1_ftqbz(dtset, cryst, qpt, drho_cplex, nfftf, ngfftf, nkxc, kxc, &
                                  vxc1_qq, non_magnetic_xc, usexcnhat, gqk%pert_comm%value)
     else
       ! Read and reconstruct the dvscf potentials and the densities for this qpt and my_npert perturbations.
       db_iqpt = dvdb%findq(qq_ibz)
       ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DVDB file."))
       ! The first entry in mapc_qq2dvdb gives the index in dvdb%qpts. The other entries in mapc_qq are OK as they refer to symmetries.
       mapc_qq2dvdb = mapc_qq; mapc_qq2dvdb(1) = db_iqpt
       call dvdb%readsym_qbz(cryst, qpt, mapc_qq2dvdb, cplex, nfftf, ngfftf, v1scf_qq, gqk%pert_comm%value)

       db_iqpt = drhodb%findq(qq_ibz)
       ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qpt), "in DRHODB file."))
       mapc_qq2dvdb = mapc_qq; mapc_qq2dvdb(1) = db_iqpt
       call drhodb%read_vxc1_qbz(dtset, cryst, qpt, mapc_qq2dvdb, drho_cplex, nfftf, ngfftf, nkxc, kxc, &
                                 vxc1_qq, non_magnetic_xc, usexcnhat, gqk%pert_comm%value)
     end if

     ABI_CHECK_IEQ(cplex, drho_cplex, "Different values of cplex for v1 and rho1!")

     cvxc1_qq_ptr => null(); if (cplex == 2) call c_f_pointer(c_loc(vxc1_qq), cvxc1_qq_ptr, [nfft, nspden, my_npert])

     ! Allocate vlocal1_qq with correct cplex. Note nvloc
     ABI_MALLOC_OR_DIE(vlocal1_qq, (cplex*n4, n5, n6, nvloc, my_npert), ierr)
     ABI_MALLOC_OR_DIE(vlocal1_mqq, (cplex*n4, n5, n6, nvloc, my_npert), ierr)

     ! Build DFPT potential at -qq
     v1scf_mq = v1scf_qq; if (cplex == 2) v1scf_mq(2,:,:,:) = -v1scf_mq(2,:,:,:)

     ! =============================================================
     ! Loop over k-points in the e-ph matrix elements (gqk%kpt_comm)
     ! =============================================================

     do my_ik=1,gqk%my_nk
       ! NB: All procs in gqk%pert_comm and gqk%bsum_com and gqk%pp_sum_comm enter this section.

       iqbuf_cnt = 1 + mod(my_iq - 1, qbuf_size)
       iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]

       ! Set entry to zero. Important as there are cycle instructions inside these loops
       ! and we don't want random numbers written to disk.
       my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = zero
       gks_atm = zero

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

       print_time_kk = my_rank == 0 .and. (my_ik <= LOG_MODK .or. mod(my_ik, LOG_MODK) == 0)
       if (print_time_kk) call cwtime(cpu_kk, wall_kk, gflops_kk, "start")

       ! Get npw_k, kg_k for kk
       call wfd%get_gvec_gbound(cryst%gmet, ecut, kk, ik_ibz, isirr_k, dtset%nloalg, & ! in
                                istwf_k, npw_k, kg_k, nkpg_k, kpg_k, gbound_k)         ! out

       ! Find k + q in the extended zone and extract symmetry info.
       ! Be careful here because there are two umklapp vectors to be considered as:
       !
       !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
       !
       kq = kk + qq_ibz

       if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq, mapl_kq) /= 0) then
         write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k+q could not be generated from a symmetrical one.",trim(ltoa(kq))
         ABI_ERROR(msg)
       end if
       ikq_ibz = mapl_kq(1); isym_kq = mapl_kq(2); trev_kq = mapl_kq(6); g0_kq = mapl_kq(3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)
       istwf_kq_ibz = wfd%istwfk(ikq_ibz); npw_kq_ibz = wfd%npwarr(ikq_ibz)

       ! Get npw_kq, kg_kq for k+q.
       call wfd%get_gvec_gbound(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, dtset%nloalg, &  ! in
                                istwf_kq, npw_kq, kg_kq, nkpg_kq, kpg_kq, gbound_kq)      ! out

       ABI_MALLOC(ug_k, (2, npw_k*nspinor))
       ABI_MALLOC(ug_kq, (2, npw_kq*nspinor))

       ! Precompute ur_nk and ur_mkq for all m and n indices treated by me
       ! TODO: Can distribute operations inside gqk%pert_comm
       do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
         call wfd%rotate_cg(n_k, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                            cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_nk(:,n_k))
       end do

       do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
         call wfd%rotate_cg(m_kq, ndat1, spin, kq_ibz, npw_kq, kg_kq, istwf_kq, &
                            cryst, mapl_kq, gbound_kq, work_ngfft, work, ug_kq, urs_kbz=ur_mkq(:,m_kq))
       end do

       ! ===========================
       ! Compute <m,k+q|vxc1_qq|n,k>
       ! ===========================
       ! Parallelized in gqk%pert_ppsum_comm
       gxc_atm = czero; cnt = 0
       do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
         do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
            cnt = cnt + 1
            if (gqk%pp_sum_comm%skip(cnt)) cycle ! MPI parallelism inside pp_sum_comm
            do imyp=1,gqk%my_npert
              if (cplex == 1) then
                ctmp_gwpc = sum(GWPC_CONJG(ur_mkq(:,m_kq)) * ur_nk(:,n_k) * vxc1_qq(1,:,spin,imyp)) / nfftf
              else
                ctmp_gwpc = sum(GWPC_CONJG(ur_mkq(:,m_kq)) * ur_nk(:,n_k) * cvxc1_qq_ptr(:,spin,imyp)) / nfftf
              end if
              ipc = gqk%my_pertcases(imyp)
              gxc_atm(:, m_kq, n_k, ipc) = [real(ctmp_gwpc), aimag(ctmp_gwpc)]
            end do ! imyp
         end do ! n_k
       end do ! m_kq

       call xmpi_sum(gxc_atm, gqk%pert_ppsum_comm%value, ierr)
       ! TODO: this is an all_gatherv but oh well.
       !call xmpi_sum(gxc_atm, gqk%pert_comm%value, ierr)

       ! Get KS g_xc in the phonon representation.
       call ephtk_gkknu_from_atm(nb, nb, 1, natom, gxc_atm, phfr_qq, displ_red_qq, gxc_nu)

       ! ===========================================================
       ! MPI sum over the pp momenta in the full BZ gqk%pp_sum_comm
       ! ===========================================================
       !
       ! Be careful here because pp should run over the list of wavevectors in the screening matrix!
       ! as pp_mesh%bz is not necessarily equivalent to the k-mesh for the wavefunctions,
       ! and we have to use the ipp_bz index to symmetrize W(pp_bz) from W(pp_ibz).
       !
       ! TODO: Should order nbz in shells so that one can reduce the memory required
       ! to store W(pp) if pp_parallelism is activated.

       gsig_atm = zero

       do ipp_bz=my_pp_start_spin(spin), my_pp_stop_spin(spin)
         ! NB: All procs in gqk%pert_comm and gqk%bsum_com enter this section.

         my_ipp = ipp_bz - my_pp_start_spin(spin) + 1
         print_time_pp = my_rank == 0 .and. (my_ipp <= LOG_MODP .or. mod(my_ipp, LOG_MODP) == 0)
         if (print_time_pp) call cwtime(cpu_pp, wall_pp, gflops_pp, "start")

         pp = pp_mesh%bz(:,ipp_bz)
         pp_is_gamma = sum(pp**2) < tol14
         qkp_string = sjoin("While treating qpt: ", ktoa(qpt), "kpt:", ktoa(kk), "pp:", ktoa(pp), ch10)

         ! Symmetry tables and g-sphere centered on k-p.
         kmp = kk - pp
         if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kmp, mapl_kmp) /= 0) then
           write(msg, '(4a)' )"k-mesh is not closed!",ch10, &
             "k-p could not be generated from a symmetrical one.",trim(ltoa(kmp))
           ABI_ERROR(msg)
         end if
         ikmp_ibz = mapl_kmp(1); isym_kmp = mapl_kmp(2); trev_kmp = mapl_kmp(6); g0_kmp = mapl_kmp(3:5)
         isirr_kmp = (isym_kmp == 1 .and. trev_kmp == 0 .and. all(g0_kmp == 0))
         kmp_ibz = ebands%kptns(:, ikmp_ibz)
         istwf_kmp_ibz = wfd%istwfk(ikmp_ibz); npw_kmp_ibz = wfd%npwarr(ikmp_ibz)

         ! Get npw_kmp, kg_kmp for k-p.
         call wfd%get_gvec_gbound(cryst%gmet, ecut, kmp, ikmp_ibz, isirr_kmp, dtset%nloalg, &  ! in
                                  istwf_kmp, npw_kmp, kg_kmp, nkpg_kmp, kpg_kmp, gbound_kmp)   ! out

         ! Compute nonlocal form factors ffnl_kmp at (k-p+G).
         ABI_MALLOC(ffnl_kmp, (npw_kmp, 1, psps%lmnmax, psps%ntypat))

         call mkffnl_objs(cryst, psps, 1, ffnl_kmp, ider0, idir0, kg_kmp, kpg_kmp, kmp, nkpg_kmp, &
                          npw_kmp, ylm_kmp, ylmgr_dum) !, comm=gqk%pert_comm%value, request=ffnl_kmp_request)


         ! Symmetry tables and g-sphere centered on k+q-p.
         kqmp = kq - pp
         if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kqmp, mapl_kqmp) /= 0) then
           write(msg,'(4a)')"k-mesh is not closed!",ch10, "k+q-p could not be generated from a symmetrical one.",trim(ltoa(kqmp))
           ABI_ERROR(msg)
         end if
         ikqmp_ibz = mapl_kqmp(1); isym_kqmp = mapl_kqmp(2); trev_kqmp = mapl_kqmp(6); g0_kqmp = mapl_kqmp(3:5)
         isirr_kqmp = (isym_kqmp == 1 .and. trev_kqmp == 0 .and. all(g0_kqmp == 0))
         kqmp_ibz = ebands%kptns(:, ikqmp_ibz)
         istwf_kqmp_ibz = wfd%istwfk(ikqmp_ibz); npw_kqmp_ibz = wfd%npwarr(ikqmp_ibz)

         ! Get npw_kqmp, kg_kqmp for k+q-p
         call wfd%get_gvec_gbound(cryst%gmet, ecut, kqmp, ikqmp_ibz, isirr_kqmp, dtset%nloalg, &    ! in
                                  istwf_kqmp, npw_kqmp, kg_kqmp, nkpg_kqmp, kpg_kqmp, gbound_kqmp)  ! out

         ! Compute nonlocal form factors ffnl_kmp at (k-p+G).
         ABI_MALLOC(ffnl_kqmp, (npw_kqmp, 1, psps%lmnmax, psps%ntypat))
         ABI_MALLOC(full_cg1_kqmp, (2, npw_kqmp*nspinor))
         ABI_MALLOC(full_cg1_kmp, (2, npw_kmp*nspinor))

         call mkffnl_objs(cryst, psps, 1, ffnl_kqmp, ider0, idir0, kg_kqmp, kpg_kqmp, kqmp, nkpg_kqmp, &
                          npw_kqmp, ylm_kqmp, ylmgr_dum) ! , comm=gqk%pert_comm%value, request=ffnl_kqmp_request)

         ! ====================================
         ! This is the g-sphere for W_{gg'}(pp)
         ! ====================================
         ! Note that in this case, the sphere is always Gamma-centered i.e. it does not depend on the pp wavevector
         ABI_CHECK_IEQ(npw_c, screen%npw, "npw_c == screen%npw")
         kg_c => screen%gvec(:, 1:npw_c)
         kg_x => gsph_x%gvec(:, 1:npw_x)
         call sphereboundary(gbound_c, istwfk1, kg_c, mgfft, npw_c)
         call sphereboundary(gbound_x, istwfk1, kg_x, mgfft, npw_x)

         ABI_MALLOC(rhotwg_c, (npw_c*nspinor))
         ABI_MALLOC(rhotwg_x, (npw_x*nspinor))
         ABI_MALLOC(vc_sqrt_gx, (npw_x))

         ! We need two stern_t objects to compute the first order change of the wavefunctions at k-p and k+q-p.
         ! Clearly, we should not duplicate the work when pp = 0.
         ! When pp == 0, we also get the gks matrix elements after stern_solve.
         ! Alternatively, one can solve the Sternheimer in the IBZ(kappa, alpha), store the results on disk
         ! and then use symmetries to reconstruct delta_u in the full BZ on the fly assuming spatial inversion or TR.
         ! Also, one should handle more carefully the integration in g_sigma around pp = Gamma in the case of semiconductors.

         nband_me = nbsum; stern_comm = xmpi_comm_self
         if (stern_has_band_para) then
           !nband_me = sigma%my_bsum_stop - sigma%my_bsum_start + 1; stern_comm = sigma%bsum_comm%value
           NOT_IMPLEMENTED_ERROR()
         end if

         ! =========================================================
         ! Get GS wavefunctions at k+q-p and store them in stern_kmp.
         ! =========================================================
         call stern_kmp%init(dtset, npw_kmp, npw_kqmp, nspinor, nbsum, nband_me, fermie1_idir_ipert, &
                             stern_use_cache, work_ngfft, mpi_enreg, stern_comm)

         do ib_sum=my_bsum_start(spin), my_bsum_stop(spin)
           call wfd%rotate_cg(ib_sum, ndat1, spin, kqmp_ibz, npw_kqmp, kg_kqmp, istwf_kqmp, &
                              cryst, mapl_kqmp, gbound_kqmp, work_ngfft, work, cg_kqmp)

           ! NB: cg_kqmp is dimensioned with mpw --> have to slice cg_kqmp
           ii = ib_sum; if (stern_kmp%has_band_para) ii = ib_sum - my_bsum_start(spin) + 1
           stern_kmp%cgq(:,:,ii) = cg_kqmp(:,1:npw_kqmp*nspinor)
         end do ! ib_sum

         ! ========================================================
         ! Get GS wavefunctions at k-p and store them in stern_kqmp
         ! ========================================================
         call stern_kqmp%init(dtset, npw_kqmp, npw_kmp, nspinor, nbsum, nband_me, fermie1_idir_ipert, &
                              stern_use_cache, work_ngfft, mpi_enreg, stern_comm)

         do ib_sum=my_bsum_start(spin), my_bsum_stop(spin)
           call wfd%rotate_cg(ib_sum, ndat1, spin, kmp_ibz, npw_kmp, kg_kmp, istwf_kmp, &
                              cryst, mapl_kmp, gbound_kmp, work_ngfft, work, cg_kmp)

           ! NB: cg_kmp is dimensioned with mpw --> have to slice cg_kmp
           ii = ib_sum; if (stern_kqmp%has_band_para) ii = ib_sum - my_bsum_start(spin) + 1
           stern_kqmp%cgq(:,:,ii) = cg_kmp(:,1:npw_kmp*nspinor)
         end do ! ib_sum

#if 1
         ! Prepare the object for applying W(pp_bz).
         ! FIXME: Sq = q+G0 with non-zero G0 is not supported.
         call screen%rotate_iqbz(ipp_bz, cryst, gsph_c, pp_mesh, vcp)

         ! Find the corresponding irred pp-point in the pp_mesh.
         call pp_mesh%get_bz_item(ipp_bz, pp, ipp_ibz, isym_pp, itim_pp)

         ! Get Fourier components of the Coulomb interaction in the BZ
         ! In 3D systems, neglecting umklapp: vc(Sq,sG) = vc(q,G) = 4pi/|q+G|**2
         ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.

         ! NOTE: vc_sqrt_gx is dimensioned with npw_x --> use rottb table from gsph_x.
         do ig=1,npw_x
           vc_sqrt_gx(gsph_x%rottb(ig, itim_pp, isym_pp)) = vcp%vc_sqrt(ig, ipp_ibz)
         end do

         !e_nk = qp_ene(n_k, ik_ibz, spin)
         !e_mkq = qp_ene(m_kq, ikq_ibz, spin)

         ! =====================================
         ! Sum over bands (n' in the equations)
         ! =====================================
         do ib_sum=my_bsum_start(spin), my_bsum_stop(spin)
           !if (ib_sum > 8) cycle
           ! NB: All procs in gqk%pert_comm enter this part.

           ! Get u_{n',k-p}(r), stored in ur_kmp
           call wfd%rotate_cg(ib_sum, ndat1, spin, kmp_ibz, npw_kmp, kg_kmp, istwf_kmp, &
                              cryst, mapl_kmp, gbound_kmp, work_ngfft, work, cg_kmp, urs_kbz=ur_kmp)

           ur_kmp = GWPC_CONJG(ur_kmp)

           ! =====================================
           ! Precompute oscillator matrix elements
           ! =====================================
           ! These terms do not depend on (idir, ipert) and can be reused in the loop over perturbations below.
           ! If the bands in the sum are distributed, one has to transmit the (m, n) indices.

           theta_mu_minus_e0i = fact_spin * qp_occ(ib_sum, ikmp_ibz, spin)
           need_x_kmp = (abs(theta_mu_minus_e0i / fact_spin) >= tol_empty) ! allow negative occ numbers
           !need_x_kmp = .True.

           ! Contract immediately over g' with the frequency convolution: \int de' Wc_{gg'}(pp, e') / (omega - e_{bsum, kmp) - e')
           ! Store results in vec_gwc_nk(:,:,n_k)
           if (gqk%pert_comm%nproc > 1) vec_gwc_nk = zero

           do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
             !if gqk%pert_comm%skip(n_k) cycle ! MPI parallelism inside pert_comm

             ! <bsum,k-p|e^{-i(p+G')}r|n,k> * vc_sqrt(p,G')
             cwork_ur = ur_kmp * ur_nk(:,n_k)

             if (need_x_kmp) then
               call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_x)
               call multiply_by_vc_sqrt("N", npw_x, nspinor, vc_sqrt_gx, rhotwg_x)
               vec_gx_nk(:,n_k) = rhotwg_x(1:npw_c*nspinor)
               ! FIXME: This is wrong if nspinor == 2
               rhotwg_c(:) = rhotwg_x(1:npw_c*nspinor)

             else
               call fft_ur(npw_c, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_c, gbound_c, cwork_ur, rhotwg_c)
               call multiply_by_vc_sqrt("N", npw_c, nspinor, vc_sqrt_gx, rhotwg_c)
             end if

             ! FIXME: npw_x should be npw_c here
             ! Prepare list of omegas: first e_nk then e_mkq for all m indices.
             omegas_nk(1) = qp_ene(n_k, ik_ibz, spin); cnt = 1
             do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
               cnt = cnt + 1; omegas_nk(cnt) = qp_ene(m_kq, ikq_ibz, spin)
             end do
             omegame0i_nk = omegas_nk - qp_ene(ib_sum, ikmp_ibz, spin)

             call screen%calc_sigc("N", nw_nk, omegame0i_nk, theta_mu_minus_e0i, dtset%zcut, &
                                   nspinor, npw_c, npw_c, rhotwg_c, vec_gwc_nk(:,:,n_k), sigcme_nk)
           end do ! n_k

           ! TODO: this is an all_gatherv but oh well.
           !call xmpi_sum(vec_gwc_nk, gqk%pert_comm, ierr)

           ! Get u_{n',k+q-p}(r), stored in ur_kqmp
           call wfd%rotate_cg(ib_sum, ndat1, spin, kqmp_ibz, npw_kqmp, kg_kqmp, istwf_kqmp, &
                              cryst, mapl_kqmp, gbound_kqmp, work_ngfft, work, cg_kqmp, urs_kbz=ur_kqmp)

           ur_kqmp = GWPC_CONJG(ur_kqmp)
           theta_mu_minus_e0i = fact_spin * qp_occ(ib_sum, ikqmp_ibz, spin)
           omegame0i_mkq = omegas_mkq - qp_ene(ib_sum, ikqmp_ibz, spin)

           need_x_kqmp = (abs(theta_mu_minus_e0i / fact_spin) >= tol_empty) ! allow negative occ numbers
           !need_x_kqmp = .True.

           ! Contract immediately over g with the frequency convolution: \int de' Wc_{gg'}(pp, e') / (omega - e_{bsum, kqmp) - e')
           ! Store results in vec_gwc_mkq(:,:,m_kq),
           if (gqk%pert_comm%nproc > 1) vec_gwc_mkq = zero

           do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
             !if (gqk%pert_sumcomm%skip(m_kq)) cycle ! MPI parallelism inside pert_comm

             ! <m,k+q|e^{i(p+G)}r|bsum,k+q-p> * vc_sqrt(p,G) -> Exchange bra and ket and take CC of the FFT in multiply_by_vc_sqrt
             cwork_ur = ur_kqmp * ur_mkq(:,m_kq)

             if (need_x_kqmp) then
               call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_x)
               call multiply_by_vc_sqrt("C", npw_x, nspinor, vc_sqrt_gx, rhotwg_x)
               vec_gx_mkq(:,m_kq) = rhotwg_x(1:npw_c*nspinor)
               rhotwg_c(:) = rhotwg_x(1:npw_c*nspinor)

             else
               call fft_ur(npw_c, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_c, gbound_c, cwork_ur, rhotwg_c)
               call multiply_by_vc_sqrt("C", npw_c, nspinor, vc_sqrt_gx, rhotwg_c)
             end if

             ! Prepare list of omegas: first e_mkq then e_nk for all n indices.
             omegas_mkq(1) = qp_ene(m_kq, ikq_ibz, spin); cnt = 1
             do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
               cnt = cnt + 1; omegas_mkq(cnt) = qp_ene(n_k, ik_ibz, spin)
             end do

             call screen%calc_sigc("T", nw_mkq, omegame0i_mkq, theta_mu_minus_e0i, dtset%zcut, &
                                   nspinor, npw_c, npw_c, rhotwg_c, vec_gwc_mkq(:,:,m_kq), sigcme_mkq)
           end do ! m_kq
           ! TODO: this is an all_gatherv but oh well.
           !call xmpi_sum(vec_gwc_mkq, gqk%pert_comm, ierr)

           ! ========================================
           ! Loop over my set of atomic perturbations
           ! ========================================
           ! For each perturbation:
           !    - setup H1 from vlocal1_qq or vlocal1_mqq.
           !    For each band in band_sum:
           !        - Solve the Sternheimer non-self-consistently and get the KS e-ph matrix elements.
           !        - Build the full first-order wavefunction including the active subspace.
           ! TODO: Should create array of gs_ham(my_npert) and rf_ham(my_npert) but I'm not sure the GPU version supports
           !       multiple instances.

           do imyp=1,gqk%my_npert
             ! NB: Only one proc enters this section. No MPI parallelism allowed here.

             idir = dvdb%my_pinfo(1, imyp); ipert = dvdb%my_pinfo(2, imyp); ipc = dvdb%my_pinfo(3, imyp)
             !print *, "For kk, ", kk, "pp:", pp, "idir, ipert", idir, ipert

             ! Set up local potential vlocal1_qq with proper dimensioning, from vtrial1 taking into account the spin
             ! and prepare application of the NL part. Each CPU prepares its own potentials.
             call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, nvloc, &
                                        pawfgr, mpi_enreg, vtrial, v1scf_qq(:,:,:,imyp), vlocal, vlocal1_qq(:,:,:,:,imyp))

             call gs_ham_kqmp%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)
             call rf_ham_kqmp%init(cplex, gs_ham_kqmp, ipert, has_e1kbsc=.true.)
             call rf_ham_kqmp%load_spin(spin, vlocal1=vlocal1_qq(:,:,:,:,imyp), with_nonlocal=.true.)

             ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
             ! TODO: Replace this call with low-level operations
             call getgh1c_setup(gs_ham_kqmp, rf_ham_kqmp, dtset, psps, kmp, kqmp, idir, ipert, &           ! In
               natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_kmp, &                                   ! In
               npw_kmp, npw_kqmp, useylmgr1, kg_kmp, ylm_kmp, kg_kqmp, ylm_kqmp, ylmgr_kqmp, &             ! In
               dkinpw, nkpg_kmp, nkpg_kqmp, kpg_kmp, kpg_kqmp, kinpw1, ffnl_kmp, ffnl_kqmp, ph3d_kmp, ph3d1_kqmp , & ! InOut
               reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)                                ! Reuse some arrays

             !call gs_hamkq%eph_setup_k("k" , kk_bz, istwf_k, npw_k, kg_k, dtset, cryst, psps, &       ! in
             !                          nkpg_k, kpg_k, ffnl_k, kinpw_k, ph3d_k, gqk%pert_comm%value)   ! out

             ! =====================
             ! NSCF Sternheimer at q
             ! =====================
             ! Compute Delta_{q,idir,ipert} \psi_{bsum, k-p}
             call wfd%rotate_cg(ib_sum, ndat1, spin, kmp_ibz, npw_kmp, kg_kmp, istwf_kmp, &
                                cryst, mapl_kmp, gbound_kmp, work_ngfft, work, cg_kmp)

             stern_kmp%bands_treated_now(:) = 0; stern_kmp%bands_treated_now(ib_sum) = 1
             stern_kmp%rank_band = 0; u1_band = ib_sum; band_me = ib_sum
             if (stern_kmp%has_band_para) then
               NOT_IMPLEMENTED_ERROR()
             end if

             call stern_kmp%solve(u1_band, band_me, idir, ipert, qpt, gs_ham_kqmp, rf_ham_kqmp, &
                                  ebands%eig(:,ikmp_ibz,spin), ebands%eig(:,ikqmp_ibz,spin), &
                                  cg_kmp, cwaveprj0, cg1_kqmp, cwaveprj, msg, ierr, &
                                  full_cg1=full_cg1_kqmp, full_ur1=full_ur1_kqmp)

             ! TODO: Last states may fail to converge
             ABI_WARNING_IF(ierr /= 0, sjoin("Stern at +q", qkp_string, msg))

             ! Store KS e-ph matrix elements for this perturbation.
             if (pp_is_gamma) then
               ib = ib_sum - gqk%bstart + 1
               gks_atm(:,:,ib,ipc) = stern_kmp%eig1_k(:, gqk%bstart:gqk%bstop, ib_sum)
             end if

             ABI_FREE(kinpw1)
             ABI_FREE(dkinpw)
             ABI_FREE(ph3d_kmp)
             ABI_SFREE(ph3d1_kqmp)

             ! <m,k+q|e^{i(p+G)r}|Delta_q psi_{bsum,k-p}> --> exchange bra and ket and take the CC of the FFT.
             full_ur1_kqmp = GWPC_CONJG(full_ur1_kqmp)

             do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
               cwork_ur = full_ur1_kqmp * ur_mkq(:,m_kq)

               if (need_x_kqmp) then
                 call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_x)
                 rhotwg_x = GWPC_CONJG(rhotwg_x)
                 rhotwg_c(:) = rhotwg_x(1:npw_c*nspinor)
               else
                 call fft_ur(npw_c, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_c, gbound_c, cwork_ur, rhotwg_c)
                 rhotwg_c(:) = GWPC_CONJG(rhotwg_c)
               end if

               do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
                 if (m_kq == n_k) then
                   ctmp_gwpc = sum(rhotwg_c(:) * vec_gwc_nk(:,1,n_k))
                 else
                   ! Take the average
                   iw_mkq = m_kq - gqk%bstart + 2
                   ctmp_gwpc = half * sum(rhotwg_c(:) * (vec_gwc_nk(:,1,n_k) + vec_gwc_nk(:,iw_mkq,n_k)))
                 end if

                 if (need_x_kqmp) then
                   ! TODO recheck
                   xdot_tmp = - xdotc(npw_x*nspinor, rhotwg_x, 1, vec_gx_nk(:,n_k), 1)
                   !ctmp_gwpc = ctmp_gwpc + xdot_tmp ! * theta_mu_minus_e0i
                 end if

                 gsig_atm(1, m_kq, n_k, ipc) = gsig_atm(1, m_kq, n_k, ipc) + real(ctmp_gwpc)
                 gsig_atm(2, m_kq, n_k, ipc) = gsig_atm(2, m_kq, n_k, ipc) + aimag(ctmp_gwpc)
               end do ! n_k
             end do ! m_kq

             ! FIXME
!if (.False.) then
if (.not. qq_is_gamma) then
             ! ===========================
             ! Same operations but for -qq
             ! ===========================
             call gs_ham_kmp%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)
             call rf_ham_kmp%init(cplex, gs_ham_kmp, ipert, has_e1kbsc=.true.)
             call rf_ham_kmp%load_spin(spin, vlocal1=vlocal1_mqq(:,:,:,:,imyp), with_nonlocal=.true.)

             ! TODO: Replace this call with low-level operations
             call getgh1c_setup(gs_ham_kmp, rf_ham_kmp, dtset, psps, kqmp, kmp, idir, ipert, &             ! In
               natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_kmp, &                                   ! In
               npw_kqmp, npw_kmp, useylmgr1, kg_kqmp, ylm_kqmp, kg_kmp, ylm_kmp, ylmgr_kmp, &              ! In
               dkinpw, nkpg_kqmp, nkpg_kmp, kpg_kqmp, kpg_kmp, kinpw1, ffnl_kqmp, ffnl_kmp, ph3d_kqmp, ph3d1_kmp, & ! InOut
               reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)                                ! Reuse some arrays

             !call gs_hamkq%eph_setup_k("k" , kk_bz, istwf_k, npw_k, kg_k, dtset, cryst, psps, &       ! in
             !                          nkpg_k, kpg_k, ffnl_k, kinpw_k, ph3d_k, gqk%pert_comm%value)   ! out

             ! ======================
             ! NSCF Sternheimer at -q
             ! ======================
             ! Compute Delta_{-q,idir,ipert} \psi_{bsum, k+q-p}
             call wfd%rotate_cg(ib_sum, ndat1, spin, kqmp_ibz, npw_kqmp, kg_kqmp, istwf_kqmp, &
                                cryst, mapl_kqmp, gbound_kqmp, work_ngfft, work, cg_kqmp)

             stern_kqmp%bands_treated_now(:) = 0; stern_kqmp%bands_treated_now(ib_sum) = 1
             stern_kqmp%rank_band = 0; u1_band = ib_sum; band_me = ib_sum
             if (stern_kqmp%has_band_para) then
               NOT_IMPLEMENTED_ERROR()
             end if

             call stern_kqmp%solve(u1_band, band_me, idir, ipert, -qpt, gs_ham_kmp, rf_ham_kmp, &
                                   ebands%eig(:,ikqmp_ibz,spin), ebands%eig(:,ikmp_ibz,spin), &
                                   cg_kqmp, cwaveprj0, cg1_kmp, cwaveprj, msg, ierr, &
                                   full_cg1=full_cg1_kmp, full_ur1=full_ur1_kmp)

             ! TODO: Last states may fail to converge
             ABI_WARNING_IF(ierr /= 0, sjoin("Stern at -q:", qkp_string, msg))

             ABI_FREE(kinpw1)
             ABI_FREE(dkinpw)
             ABI_FREE(ph3d_kqmp)
             ABI_SFREE(ph3d1_kmp)

             full_ur1_kmp = GWPC_CONJG(full_ur1_kmp)

             do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%m_start, gqk%m_stop

               ! <Delta_{-q} psi_{bsum,k+q-p}|e^{-i(p+G')r}|n,k>
               cwork_ur = full_ur1_kqmp * ur_nk(:,n_k)

               if (need_x_kmp) then
                 call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_x)
                 rhotwg_c(:) = rhotwg_x(1:npw_c*nspinor)
               else
                 call fft_ur(npw_c, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_c, gbound_c, cwork_ur, rhotwg_c)
               end if

               do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%n_start, gqk%n_stop
                 if (m_kq == n_k) then
                   ctmp_gwpc = sum(rhotwg_c(:) * vec_gwc_mkq(:,1,m_kq))
                 else
                   ! Take the average
                   iw_nk = n_k - gqk%bstart + 2
                   ctmp_gwpc = half * sum(rhotwg_c(:) * (vec_gwc_mkq(:,1,m_kq) + vec_gwc_mkq(:,iw_nk,m_kq)))
                 end if

                 if (need_x_kmp) then
                   xdot_tmp = - xdotc(npw_x*nspinor, vec_gx_mkq(:,m_kq), 1, rhotwg_x, 1)
                   ! TODO recheck
                   !ctmp_gwpc = ctmp_gwpc + xdot_tmp ! * theta_mu_minus_e0i
                 end if

                 gsig_atm(1, m_kq, n_k, ipc) = gsig_atm(1, m_kq, n_k, ipc) +  real(ctmp_gwpc)
                 gsig_atm(2, m_kq, n_k, ipc) = gsig_atm(2, m_kq, n_k, ipc) + aimag(ctmp_gwpc)
               end do ! n_k
             end do ! m_kq
end if ! .not qq_is_gamma.

           end do  ! imyp (my perturbations)
           call rf_ham_kqmp%free(); call rf_ham_kmp%free()
         end do ! ib_sum (sum over bands)
#endif

         if (print_time_pp) then
           call inds2str(1, " My pp-point:", my_ipp, my_npp(spin), pp_mesh%nbz, msg)
           call cwtime_report(msg, cpu_pp, wall_pp, gflops_pp); if (my_ipp == LOG_MODP) call wrtout(std_out, "...", do_flush=.True.)
         end if

         ABI_FREE(kpg_kmp)
         ABI_FREE(kpg_kqmp)
         ABI_FREE(ffnl_kmp)
         ABI_FREE(ffnl_kqmp)
         ABI_FREE(full_cg1_kqmp)
         ABI_FREE(full_cg1_kmp)
         ABI_FREE(rhotwg_c)
         ABI_FREE(rhotwg_x)
         ABI_FREE(vc_sqrt_gx)
         call stern_kmp%free(); call stern_kqmp%free()
       end do ! ipp_bz

       ABI_FREE(kpg_k)
       ABI_FREE(kpg_kq)
       ABI_FREE(ug_k)
       ABI_FREE(ug_kq)

       ! Here we are outside of the loop over pp_sum, band_sum and perturbations.
       ! Collect gsig_atm and gks_atm inside pert_ppsum_comm so that all procs can operate on the data.
       call xmpi_sum_master(gsig_atm, master, gqk%pert_ppsum_bsum_comm%value, ierr)
       call xmpi_sum_master(gks_atm , master, gqk%pert_ppsum_bsum_comm%value, ierr)

       ! TODO gks_atm and gks_nu
       !gsig_atm = gsig_atm / (cryst%ucvol * pp_mesh%nbz)
       gsig_atm = gsig_atm + gks_atm - gxc_atm

       select case (gstore%gmode)
       case (GSTORE_GMODE_PHONON)
         ! FIXME Perhaps it's gonna be easier if we only support GMODE_ATOM
         ! Get g^{Sigma} in the phonon representation.
         call ephtk_gkknu_from_atm(nb, nb, 1, natom, gsig_atm, phfr_qq, displ_red_qq, gsig_nu)

         ! Save e-ph matrix elements in the buffer.
         select case (gqk%cplex)
         case (1)
           my_gbuf(1,:,:,:, my_ik, iqbuf_cnt) = gsig_nu(1,:,:,:) ** 2 + gsig_nu(2,:,:,:) ** 2
         case (2)
           my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = gsig_nu
         end select

       case (GSTORE_GMODE_ATOM)
         ! Save e-ph matrix elements in the buffer.
         select case (gqk%cplex)
         case (1)
           my_gbuf(1,:,:,:, my_ik, iqbuf_cnt) = gsig_atm(1,:,:,:)**2 + gsig_atm(2,:,:,:)** 2
         case (2)
           my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = gsig_atm
         end select

       case default
         ABI_ERROR(sjoin("Invalid gstore%gmode:", gstore%gmode))
       end select

       ! Dump buffer
       if (iqbuf_cnt == qbuf_size) call dump_my_gbuf()

       if (print_time_kk) then
         call inds2str(3, "My k-point", my_ik, gqk%my_nk, gqk%glob_nk, msg)
         call cwtime_report(msg, cpu_kk, wall_kk, gflops_kk); if (my_ik == LOG_MODK) call wrtout(std_out, "...", do_flush=.True.)
       end if
     end do ! my_ik

     ABI_FREE(v1scf_qq)
     ABI_FREE(vlocal1_qq)
     ABI_FREE(vlocal1_mqq)
     ABI_FREE(vxc1_qq)

     if (print_time_qq) then
       call inds2str(2, "My q-point", my_iq, gqk%my_nq, gqk%glob_nq, msg)
       call cwtime_report(msg, cpu_qq, wall_qq, gflops_qq); if (my_iq == LOG_MODQ) call wrtout(std_out, "...", do_flush=.True.)
     end if
   end do ! iq_ibz

   ! Dump the remainder.
   if (iqbuf_cnt /= 0) call dump_my_gbuf()

   ABI_FREE(ur_nk)
   ABI_FREE(ur_mkq)
   ABI_FREE(iq_buf)
   ABI_FREE(my_gbuf)
   ABI_FREE(gsig_atm)
   ABI_FREE(gsig_nu)
   ABI_FREE(gks_atm)
   ABI_FREE(gks_nu)
   ABI_FREE(gxc_atm)
   ABI_FREE(gxc_nu)
   ABI_FREE(omegas_nk)
   ABI_FREE(omegas_mkq)
   ABI_FREE(omegame0i_nk)
   ABI_FREE(omegame0i_mkq)
   ABI_FREE(vec_gwc_nk)
   ABI_FREE(vec_gwc_mkq)
   ABI_FREE(vec_gx_nk)
   ABI_FREE(vec_gx_mkq)
   ABI_FREE(sigcme_nk)
   ABI_FREE(sigcme_mkq)
 end do ! my_is

 call cwtime_report(" gwpt_eph full calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Set gstore_completed to 1 so that we can easily check if restarted is needed.
 !if (my_rank == master) then
   NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_completed"), 1))
 !end if
 NCF_CHECK(nf90_close(root_ncid))
 call xmpi_barrier(comm)

 ! Output some of the results to ab_out for testing purposes
 call gstore%print_for_abitests(dtset)

 ! Free memory
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(kg_kmp)
 ABI_FREE(kg_kqmp)
 ABI_FREE(ylm_k)
 ABI_FREE(ylm_kq)
 ABI_FREE(ylm_kmp)
 ABI_FREE(ylm_kqmp)
 ABI_FREE(ylmgr_kq)
 ABI_FREE(ylmgr_kmp)
 ABI_FREE(ylmgr_kqmp)
 ABI_FREE(cg_work)
 ABI_FREE(ur_kmp)
 ABI_FREE(ur_kqmp)
 ABI_FREE(full_ur1_kqmp)
 ABI_FREE(full_ur1_kmp)
 ABI_FREE(cwork_ur)
 ABI_FREE(cg_kmp)
 ABI_FREE(cg_kqmp)
 ABI_FREE(cg1_kqmp)
 ABI_FREE(cg1_kmp)
 ABI_FREE(vtrial)
 ABI_FREE(work)
 ABI_FREE(ph1d)
 ABI_FREE(vlocal)
 ABI_FREE(displ_cart_qq)
 ABI_FREE(displ_red_qq)
 ABI_FREE(gbound_k)
 ABI_FREE(gbound_kq)
 ABI_FREE(gbound_kmp)
 ABI_FREE(gbound_kqmp)
 ABI_FREE(gbound_c)
 ABI_FREE(gbound_x)
 ABI_FREE(done_qbz_spin)
 ABI_FREE(rhor)
 ABI_FREE(kxc)
 ABI_FREE(vxc)

 call gs_ham_kqmp%free(); call gs_ham_kmp%free(); call wfd%free(); call vcp%free()
 call screen%free(); call pp_mesh%free(); call gsph_c%free(); call gsph_x%free(); call gstore%free()
 call pawcprj_free(cwaveprj0); call pawcprj_free(cwaveprj)
 ABI_FREE(cwaveprj0)
 ABI_FREE(cwaveprj)

 call xmpi_barrier(comm) ! This to make sure that the parallel output of GSTORE is completed
 call cwtime_report(" gwpt_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

contains

subroutine inds2str(level, prefix, my_ik, my_nk, nk_tot, out_str)
 character(len=*),intent(in) :: prefix
 integer,intent(in) :: level, my_ik, my_nk, nk_tot
 character(len=*),intent(out) :: out_str

 out_str = sjoin(prefix, itoa(my_ik), "/", itoa(my_nk), "[", itoa(nk_tot), "]")
 out_str = repeat(' ', 4 * level) // trim(out_str)
end subroutine  inds2str

subroutine dump_my_gbuf()

 ! This function is called inside the double loop over (my_is, my_iq) or when we exit
 ! from the my_iq loop to dump the remainder that is still in the q-buffer,
 ! All the MPI procs in the (kpt_comm x pert_comm) grid shall call this contained routine
 ! as we have side-effects i.e. iqbuf_cnt set to 0.

 ! On disk we have the global arrays:
 !
 !      nctkarr_t("gvals", "dp", "gstore_cplex, nb, nb, natom3, glob_nk, glob_nq")
 !
 ! while the local MPI buffers are dimensioned as follows:
 !
 !      my_gbuf(2, nb, nb, natom3, gqk%my_nk, qbuf_size)

 ! If parallelism over perturbation is activated, only the procs treating the first perturbation
 ! i.e. the procs treating different k-points for this q are involved in IO
 ! as all the local buffers store results for all natom3 perturbations.

 integer :: ii, iq_bz, iq_glob, my_iq
 !return

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! FIXME: Recheck this part as we have way more levels of parallelism in GWPT
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 if (gqk%coords_qkpb_sumbp(3) /= 0) goto 10 ! Yes, I'm very proud of this GOTO.
 if (gqk%pert_ppsum_bsum_comm%me /= 0) goto 10 ! Yes, I'm very proud of this GOTO.

 !iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]
 my_iq = iq_buf(1, 1)
 iq_glob = my_iq + gqk%my_qstart - 1

 !print *, "in dump_my_gbuf with start: ", [1, 1, 1, 1, gqk%my_kstart, iq_glob]
 !print *, "                  count; ", [gqk%cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%my_nk, iqbuf_cnt]
 !print *, "my_gbuf", my_gbuf(:,:,:,natom3,1,1)

 ! NB: this is an individual IO operation
 ncerr = nf90_put_var(spin_ncid, spin_vid("gvals"), my_gbuf, &
                      start=[1, 1, 1, 1, gqk%my_kstart, iq_glob], &
                      count=[gqk%cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%my_nk, iqbuf_cnt])
 NCF_CHECK(ncerr)

 ! Only one proc sets the entry in done_qbz_spin to 1 for all the q-points in the buffer.
 !if (all(gqk%coords_qkpb_sumbp(2:3) == [0, 0]))  then
   do ii=1,iqbuf_cnt
     iq_bz = iq_buf(2, ii)
     NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_done_qbz_spin"), 1, start=[iq_bz, spin]))
   end do
 !end if

 ! Zero the counter before returning
10 iqbuf_cnt = 0

 !NCF_CHECK(nf90_sync(spin_ncid))
 !NCF_CHECK(nf90_sync(root_ncid))

end subroutine dump_my_gbuf

integer function root_vid(var_name)
  character(len=*),intent(in) :: var_name
  root_vid = nctk_idname(root_ncid, var_name)
end function root_vid

integer function spin_vid(var_name)
  character(len=*),intent(in) :: var_name
  spin_vid = nctk_idname(spin_ncid, var_name)
end function spin_vid

end subroutine gwpt_run
!!***

!!****f* m_gwpt/multiply_by_vc_sqrt
!! NAME
!!  multiply_by_vc_sqrt
!!
!! FUNCTION
!! Multiply rhotwg by the square root of the Coulomb term taking into account nspinor.
!!
!! INPUTS
!!
!! SOURCE

subroutine multiply_by_vc_sqrt(trans, npw, nspinor, vc_sqrt_gx, rhotwg)

 character(len=1),intent(in) :: trans
 integer,intent(in) :: npw, nspinor
 complex(gwpc),intent(in) :: vc_sqrt_gx(npw)
 complex(gwpc),intent(inout) :: rhotwg(npw*nspinor)

!Local variables ------------------------------
 integer :: ii, spad
!************************************************************************

 select case (trans)
 case ("N")
   do ii=1,nspinor
     spad = (ii-1) * npw
     rhotwg(spad+1:spad+npw) = rhotwg(spad+1:spad+npw) * vc_sqrt_gx(1:npw)
   end do

 case ("C")
   ! Take the complex conjugate of rhotwg.
   do ii=1,nspinor
     spad = (ii-1) * npw
     rhotwg(spad+1:spad+npw) = GWPC_CONJG(rhotwg(spad+1:spad+npw)) * vc_sqrt_gx(1:npw)
   end do

 case default
   ABI_ERROR(sjoin("Invalid trans", trans))
 end select

end subroutine multiply_by_vc_sqrt
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
!!   write(unt, "(2(a,i0),/)")"P Number of k-points in gwpt_nk treated by this proc: ", gwpt%my_nkcalc, " of ", gwpt%nkcalc
!!
!!  end subroutine gwpt_print
!!***

end module m_gwpt
!!***

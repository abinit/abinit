!!****m* ABINIT/m_eph_path
!! NAME
!!  m_eph_path
!!
!! FUNCTION
!!  Compute e-ph matrix elements g(k,q) along an arbitrary path either in k- or q-space.
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

module m_eph_path

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_ifc
 use m_ddb
 use m_dvdb
 use m_copy
 use m_hamiltonian
 use m_pawcprj
 use m_ephtk
 use netcdf
 use m_nctk
 use m_dtset
 use m_dtfil

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : pseudopotential_type
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_cgtools,        only : cg_zdotc
 use m_crystal,        only : crystal_t
 use m_ebands,         only : ebands_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawfgr,         only : pawfgr_type
 use m_cgwf,           only : nscf_t
 use m_bz_mesh,        only : kpath_t
 use m_wfd,            only : u0_cache_t

 implicit none

 private
!!***

 public :: eph_path_run

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_eph_path/eph_path_run
!! NAME
!!  eph_path_run
!!
!! FUNCTION
!!  Compute e-ph matrix elements g(k,q) along an arbitrary path either in k- or q-space.
!!  Wavefunctions at k and k+q are computed non-self-consistently by invoking the CG eigensolver
!!  starting from the GS potential read from file.
!!  The DFPT potentials at q are usually obtained via Fourier interpolation, but it is also possible
!!  to use fully ab-initio potentials provided the DVDB file contains all the q-points along the path.
!!  This requires performing DFPT calculations for all the q-points, and then merging
!!  all the POT1 files with the mrgdv utility.
!!
!! INPUTS
!! dtfil<datafiles_type>=Variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! cryst<crystal_t>=crystal structure parameters
!! wfk_ebands: electron bands from the input WFK file (used to propagate nelect and fermie to GPATH.nc)
!! dvdb<dbdb_type>=Database with the DFPT SCF potentials.
!! ifc<ifc_type>=interatomic force constants and corresponding real space grid info.
!! pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!! Only writing.
!!
!! SOURCE

subroutine eph_path_run(dtfil, dtset, cryst, wfk_ebands, dvdb, ifc, pawfgr, pawang, pawrad, pawtab, psps, comm)

!Arguments ------------------------------------
!scalars
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: wfk_ebands
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 integer,intent(in) :: comm
!arrays
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwfk_1 = 1, tim_getgh1c = 1, berryopt0 = 0, useylmgr1 = 0, master = 0, ndims=3, paral_kgb0 = 0, ndat1 = 1
 integer :: sij_opt,usecprj,usevnl,optlocal,optnl,opt_gvnlx1, nu
 integer :: spin, iq, ik, nk_path, nq_path, ierr, npw_k, npw_kq, my_rank, nprocs, n1, n2, n3, n4, n5, n6, cplex
 integer :: natom, natom3, nsppol, nspden, nspinor, qptopt, comm_cart, me_cart
 integer :: nfft,nfftf,mgfft,mgfftf, my_npert, my_ip, idir, ipert, ipc, ncerr, ncid, my_nkpath, my_nqpath
 integer :: in_k, im_kq, my_is, my_ik, my_iq, nband, nb_in_g, ii, band_n, band_m, bstart, bstop, my_nspins, np, tot_nscf_ierr
 real(dp) :: cpu_all,wall_all,gflops_all, eig0nk, eshift
 logical :: qq_is_gamma, need_ftinterp, gen_eigenpb, use_cg_k, use_cg_kq, use_cache
 type(gs_hamiltonian_type) :: gs_ham_k, gs_ham_kq
 type(rf_hamiltonian_type) :: rf_ham_kq
 type(nscf_t) :: nscf
 type(kpath_t) :: qpath, kpath
 type(xcomm_t) :: kpt_comm, qpt_comm, pert_comm
 type(u0_cache_t) :: ucache_kq, ucache_k
 character(len=fnlen) :: gpath_path
 character(len=5000) :: msg
 character(len=10) :: priority
!!arrays
 integer :: units(2), ngfft(18),ngfftf(18), coords_spin(ndims), dims(ndims)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), qmap_symrec(:,:), my_ik_inds(:), my_iq_inds(:), my_spins(:), my_iperts(:)
 integer,allocatable :: pert_table(:,:), my_pinfo(:,:)
 real(dp) :: kk(3), qq(3), kq(3), phfreqs(3*cryst%natom), phfreqs_ev(3*cryst%natom), fake_path(3,2)
 real(dp),allocatable :: grad_berry(:,:), kinpw_k(:), kinpw_kq(:)
 real(dp),allocatable :: cg_k(:,:,:), cg_kq(:,:,:), gsc_k(:,:,:), gsc_kq(:,:,:),eig_k(:), eig_kq(:)
 real(dp),allocatable :: v1scf(:,:,:,:), vlocal1(:,:,:,:), vlocal(:,:,:,:), gkq_atm(:,:,:,:), gkq_nu(:,:,:,:), gkq2_nu(:,:,:)
 real(dp),allocatable :: gvnlx1(:,:), gs1c(:,:), h1_kets_kq(:,:,:), displ_cart(:,:,:,:),displ_red_qq(:,:,:,:)
 real(dp),allocatable :: kpg_k(:,:), ph3d_k(:,:,:), ffnl_k(:,:,:,:), vlocal_k(:,:,:,:)
 real(dp),allocatable :: kpg_kq(:,:), ph3d_kq(:,:,:), ffnl_kq(:,:,:,:), vlocal_kq(:,:,:,:), real_vec(:)
 logical :: reorder, periods(ndims), keepdim(ndims)
 type(pawcprj_type),allocatable  :: cwaveprj0(:,:)
 type(xcomm_t),allocatable :: comm_my_is(:)
!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 call cwtime(cpu_all, wall_all, gflops_all, "start")
 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 units = [std_out, ab_out]

 ! Copy important dimensions.
 natom = cryst%natom; natom3 = 3 * natom; nsppol = dtset%nsppol; nspinor = dtset%nspinor; nspden = dtset%nspden

 ! Build (k/q)-path. NB: Input variables have been already checked for consistency in chkinp.
 select case (dtset%eph_fix_korq)
 case ("k")
   call qpath%init(dtset%ph_qpath(:,1:dtset%ph_nqpath), cryst%gprimd, dtset%ph_ndivsm)
   nq_path = qpath%npts
   call qpath%print(units, header=sjoin("q-point path for g(k,q) with fixed k:", ktoa(dtset%eph_fix_wavevec)), prtvol=dtset%prtvol)
   fake_path(:,1) = dtset%eph_fix_wavevec; fake_path(:,2) = dtset%eph_fix_wavevec + one
   call kpath%init(fake_path, cryst%gprimd, 0)
   nk_path = 1

 case ("q")
   call kpath%init(dtset%kptbounds(:,1:dtset%nkpath), cryst%gprimd, dtset%ndivsm)
   nk_path = kpath%npts
   call kpath%print(units, header=sjoin("k-point path for g(k,q) with fixed q:", ktoa(dtset%eph_fix_wavevec)), prtvol=dtset%prtvol)
   fake_path(:,1) = dtset%eph_fix_wavevec; fake_path(:,2) = dtset%eph_fix_wavevec + one
   call qpath%init(fake_path, cryst%gprimd, 0)
   nq_path = 1

 case default
   ABI_ERROR(sjoin("Invalid value of eph_fix_korq:", dtset%eph_fix_korq))
 end select

 ! Define band range and nb_in_g from eph_path_brange.
 nband = dtset%mband; bstart = dtset%eph_path_brange(1); bstop = dtset%eph_path_brange(2)
 if (bstart <= 0) bstart = 1
 if (bstop <= 0) bstop = nband
 nb_in_g = bstop - bstart + 1

 ! The values of eph_path_brange must be validated at this level!
 ABI_CHECK_IRANGE(bstart, 1, nband, "Wrong eph_path_brange(1)")
 ABI_CHECK_IRANGE(bstop, 1, nband, "Wrong eph_path_brange(2)")
 ABI_CHECK_IGEQ(bstop, bstart, "eph_path_brange(2) < eph_path_brange(1)")
 call wrtout(units, sjoin("Computing g with eph_path_brange:", ltoa([bstart, bstop])))

 ! Distribute spins inside input comm.
 call xmpi_split_nsppol(comm, nsppol, my_nspins, my_spins, comm_my_is)

 ! ==================
 ! MPI cartesian grid
 ! ==================
 my_npert = natom3
 do my_is=1,my_nspins
   spin = my_spins(my_is)
   np = comm_my_is(my_is)%nproc; priority = "12"

   if (any(dtset%eph_np_pqbks /= 0)) then
     ! Take MPI grid from input.
     pert_comm%nproc = dtset%eph_np_pqbks(1)
     qpt_comm%nproc  = dtset%eph_np_pqbks(2)
     ABI_CHECK_IEQ(dtset%eph_np_pqbks(3), 1, "Band parallelism not implemented in eph_path")
     kpt_comm%nproc = dtset%eph_np_pqbks(4)

   else
     ! Automatic MPI grid generation.
     if (nk_path == 1) then
       call xmpi_distrib_2d(np, priority, nq_path, natom3, qpt_comm%nproc, pert_comm%nproc, ierr)
     else
       call xmpi_distrib_2d(np, priority, nk_path, natom3, kpt_comm%nproc, pert_comm%nproc, ierr)
     end if
     ABI_CHECK(ierr == 0, sjoin("Cannot distribute nprocs:", itoa(np), " with priority: ", priority, ". Decrease MPI nprocs"))
   end if

   ! Consistency check
   write(msg, "(a,i2,a,3(i0,1x))")"P Cartesian grid for spin", spin, ": (pert_comm%nproc, qpt_comm%nproc, kpt_comm%nproc) = ", &
                                  pert_comm%nproc, qpt_comm%nproc, kpt_comm%nproc
   call wrtout(units, msg)

   if (pert_comm%nproc * qpt_comm%nproc * kpt_comm%nproc /= np) then
     write(msg, "(a,i0,3a, 4(a,1x,i0))") &
       "Cannot create Cartesian grid with nproc(spin): ", np, ch10, &
       "Idle processes are not supported. The product of the `nproc_*` vars should be equal to nproc.", ch10, &
       "qpt_nproc (", qpt_comm%nproc, ") x kpt_nproc (", kpt_comm%nproc, ") x pert_nproc", pert_comm%nproc, &
       ") != ", np
     ABI_ERROR(msg)
   end if

   ! For each spin treated by this rank, create MPI cartesian communicator of rank ndims.
   periods(:) = .False.; reorder = .False.
   dims = [pert_comm%nproc, qpt_comm%nproc, kpt_comm%nproc]

#ifdef HAVE_MPI
   call MPI_CART_CREATE(comm_my_is(my_is)%value, ndims, dims, periods, reorder, comm_cart, ierr)
   ! Find the index and coordinates of the current processor
   call MPI_COMM_RANK(comm_cart, me_cart, ierr)
   call MPI_CART_COORDS(comm_cart, me_cart, ndims, coords_spin, ierr)

   ! Communicator for q-points in g(k,q)
   keepdim = .False.; keepdim(1) = .True.; call pert_comm%from_cart_sub(comm_cart, keepdim)
   keepdim = .False.; keepdim(2) = .True.; call qpt_comm%from_cart_sub(comm_cart, keepdim)
   keepdim = .False.; keepdim(3) = .True.; call kpt_comm%from_cart_sub(comm_cart, keepdim)
   call xmpi_comm_free(comm_cart)
#endif
 end do ! my_is

 ! Distribute k-points (q-points) inside kpt_comm (qpt_comm) using block distribution.
 call xmpi_split_block(nk_path, kpt_comm%value, my_nkpath, my_ik_inds)
 call xmpi_split_block(nq_path, qpt_comm%value, my_nqpath, my_iq_inds)
 call xmpi_split_block(natom3, pert_comm%value, my_npert, my_iperts)
 ABI_FREE(my_iperts)

 ! Idle processors are not supported (tested).
 ABI_CHECK_IGEQ(my_nkpath, 1, "Too many procs for k-point parallelism.")
 ABI_CHECK_IGEQ(my_nqpath, 1, "Too many procs for q-point parallelism.")
 ABI_CHECK_IGEQ(my_npert, 1, "Too many procs for perturbation parallelism.")

 if (pert_comm%nproc > 1) then
   ! Build table with list of perturbations treated by this CPU inside pert_comm
   call ephtk_set_pertables(cryst%natom, my_npert, pert_table, my_pinfo, pert_comm%value)
   ! Activate parallelism over perturbations
   call dvdb%set_pert_distrib(my_npert, natom3, my_pinfo, pert_table, pert_comm%value)
   ABI_FREE(my_pinfo)
   ABI_FREE(pert_table)
   ABI_WARNING("Parallelism over perturbations should be tested!")
 end if

 ! Load KS potential from file.
 call nscf%init(dtset, dtfil, cryst, comm)

 ! FFT meshes (taken from the GS POT file)
 ngfft = nscf%ngfft; ngfftf = nscf%ngfftf
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Open the DVDB file and make sure we have POT1 files.
 call dvdb%open_read(ngfftf, xmpi_comm_self)
 ABI_CHECK(dvdb%has_fields("pot1", msg), msg)

 ! Check if all the q-points are present in the DVDB.
 qptopt = dtset%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt
 call dvdb%need_ftinterp(nq_path, qpath%points, qptopt, qmap_symrec, need_ftinterp)

 if (.not. need_ftinterp .and. dtset%eph_use_ftinterp /= 0) then
   ABI_COMMENT("Enforcing FT interpolation for q-points even if it's not strictly needed.")
   need_ftinterp = .True.
 end if

 if (need_ftinterp) then
   call wrtout(units, " Cannot find all q-points in the DVDB --> Activating Fourier interpolation.")
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, xmpi_comm_self)
 else
   call wrtout(units, " DVDB file contains all q-points along the path --> Reading DFPT potentials from file.")
 end if

 ! Prepare call to getgh1c
 usevnl = 0
 optlocal = 1    ! local part of H^(1) is computed in gh1c=<G|H^(1)|C>
 optnl = 2       ! non-local part of H^(1) is totally computed in gh1c=<G|H^(1)|C>
 opt_gvnlx1 = 0  ! gvnlx1 is output
 usecprj = 0

 ABI_MALLOC(gvnlx1, (2, usevnl))
 ABI_MALLOC(grad_berry, (2, nspinor*(berryopt0/4)))
 ABI_MALLOC(cwaveprj0, (natom, nspinor*usecprj))
 ABI_MALLOC(displ_cart, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(displ_red_qq, (2, 3, cryst%natom, natom3))
 ABI_MALLOC(gkq_atm, (2, nb_in_g, nb_in_g, natom3))
 ABI_MALLOC(gkq_nu, (2, nb_in_g, nb_in_g, natom3))
 ABI_MALLOC(gkq2_nu, (nb_in_g, nb_in_g, natom3))

 ! Master writes metadata to GPATH file.
 gpath_path = strcat(dtfil%filnam_ds(4), "_GPATH.nc")

 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(ncid, gpath_path, xmpi_comm_self))
   ! Add crystalline structure.
   NCF_CHECK(cryst%ncwrite(ncid))

   ! Write dimensions.
   ncerr = nctk_def_dims(ncid, [ &
      nctkdim_t("nspinor", nspinor), nctkdim_t("nspden", nspden), nctkdim_t("nsppol", nsppol), &
      nctkdim_t("nband", nband), nctkdim_t("nb_in_g", nb_in_g), &
      nctkdim_t("nq_path", nq_path), nctkdim_t("nk_path", nk_path), &
      nctkdim_t("natom", cryst%natom), nctkdim_t("natom3", natom3), nctkdim_t("number_of_phonon_modes", natom3) &
   ], defmode=.True.)
   NCF_CHECK(ncerr)

   ! integer scalars
   ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
     "bstart", "bstop", "dvdb_add_lr", "used_ftinterp" &
   ])
   NCF_CHECK(ncerr)

   ! double precision scalars
   ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
     "nelect", "fermie"  &
   ])
   NCF_CHECK(ncerr)

   ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t("kpoints", "dp", "three, nk_path"), &
     nctkarr_t("qpoints", "dp", "three, nq_path"), &
     nctkarr_t("qweights", "dp", "nq_path"), &
     nctkarr_t("all_eigens_k", "dp", "nband, nk_path, nsppol"), &
     nctkarr_t("all_eigens_kq", "dp", "nband, nq_path, nsppol"), &
     nctkarr_t("gkq2_nu", "dp", "nb_in_g, nb_in_g, natom3, nq_path, nk_path, nsppol"), &
     nctkarr_t("phfreqs", "dp", "natom3, nq_path"), &
     nctkarr_t("phdispl_cart", "dp", "two, three, natom, natom3, nq_path"), &
     nctkarr_t("eph_fix_korq", "c", "one"), &
     nctkarr_t("eph_fix_wavevec", "dp", "three") &
   ])
   NCF_CHECK(ncerr)

   NCF_CHECK(nf90_def_var_fill(ncid, vid("gkq2_nu"), NF90_FILL, -one))
   NCF_CHECK(nf90_def_var_fill(ncid, vid("phfreqs"), NF90_FILL, -one))

   ! Write data.
   NCF_CHECK(nctk_set_datamode(ncid))
   ii = merge(1, 0, need_ftinterp)
   ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
     "bstart", "bstop", "dvdb_add_lr", "used_ftinterp"], &
     [bstart, bstop, dtset%dvdb_add_lr, ii  &
   ])
   NCF_CHECK(ncerr)

   ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
     "nelect", "fermie"], &
     [wfk_ebands%nelect, wfk_ebands%fermie &
   ])
   NCF_CHECK(ncerr)

   ! arrays
   NCF_CHECK(nf90_put_var(ncid, vid("kpoints"), kpath%points(:,1:nk_path)))
   NCF_CHECK(nf90_put_var(ncid, vid("qpoints"), qpath%points(:,1:nq_path)))
   ABI_MALLOC(real_vec, (nq_path))
   real_vec = one
   NCF_CHECK(nf90_put_var(ncid, vid("qweights"), real_vec))
   ABI_FREE(real_vec)
   NCF_CHECK(nf90_put_var(ncid, vid("eph_fix_korq"), dtset%eph_fix_korq))
   NCF_CHECK(nf90_put_var(ncid, vid("eph_fix_wavevec"), dtset%eph_fix_wavevec))

   ! Compute non-analytic phonons for q--> 0 in polar materials.
   if (nq_path > 1 .and. (any(ifc%zeff /= zero))) call ifc%calcnwrite_nana_terms_qpath(qpath, cryst, ncid, units)
   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Make sure the netcdf file has been written by master before continuing.
 call xmpi_barrier(comm)

 ! All procs open the GPATH file here.
 ! FIXME
 !NCF_CHECK(nctk_open_modify(ncid, gpath_path, comm))
 NCF_CHECK(nctk_open_modify(ncid, gpath_path, xmpi_comm_self))

 ! The cache allows one the reuse the wavefunctions of the previous k/q to init the NSCF cycle
 ! It usually reduces the number of iterations by 3-4 but it requires more memory.
 tot_nscf_ierr = 0
 use_cache = .True. !; use_cache = .False.
 call ucache_k%init(use_cache .and. my_nkpath > 1, ngfft)
 call ucache_kq%init(use_cache .and. my_nqpath > 1, ngfft)

 ! Loop over spins (MPI parallelized)
 do my_is=1,my_nspins
   spin = my_spins(my_is)

   ! Loop over k-points in k-path (MPI parallelized).
   do my_ik=1,my_nkpath
     ik = my_ik_inds(my_ik); kk = kpath%points(:, ik)
     !print *, "ik, kk", ik, kk

     ! Prepare NSCF run at k.
     ! gs_ham_k has pointers to the *_k arrays in output so we cannot deallocate them till the end.
     ! This is the reason why we use vlocal_k although this term does not depend on k
     call nscf%setup_kpt(spin, kk, istwfk_1, nband, cryst, dtset, psps, pawtab, pawfgr, &              ! in
                         npw_k, kg_k, kpg_k, ph3d_k, kinpw_k, ffnl_k, vlocal_k, cg_k, gsc_k, gs_ham_k) ! out

     ! Cache to initialize u_{nk}(g).
     use_cg_k = (my_ik > 1 .and. ucache_k%use_cache)
     if (use_cg_k) call ucache_k%get_kpt(kk, istwfk_1, npw_k, nspinor, nband, kg_k, cg_k)

     ! Compute u_{nk}(g)
     call nscf%solve_kpt(spin, kk, istwfk_1, nband, cryst, dtset, dtfil, gs_ham_k, &
                         use_cg_k, npw_k, cg_k, gsc_k, eig_k, msg, ierr)

     ABI_WARNING_IF(ierr /= 0, msg)
     tot_nscf_ierr = tot_nscf_ierr + ierr

     call ucache_k%store_kpt(kk, istwfk_1, npw_k, nspinor, nband, kg_k, cg_k)

     !if (pert_comm%me == master) then
     NCF_CHECK(nf90_put_var(ncid, vid("all_eigens_k"), eig_k, start=[1,ik,spin]))
     !end if

     ! Make sure all procs in pert_comm have the same wavefunctions at k.
     call xmpi_bcast(cg_k, master, pert_comm%value, ierr)
     if (psps%usepaw == 1) call xmpi_bcast(gsc_k, master, pert_comm%value, ierr)

     ! Allocate vlocal. Note nvloc
     ABI_MALLOC(vlocal, (n4, n5, n6, gs_ham_k%nvloc))

     ! Loop over q-points in q-path (MPI parallelized).
     ! All procs in pert_comm enter this loop with the same ik/iq indices.
     do my_iq=1,my_nqpath
       iq = my_iq_inds(my_iq); qq = qpath%points(:,iq); qq_is_gamma = sum(qq**2) < tol14
       kq = kk + qq
       !print *, "iq, kq", iq, kq

       ! Prepare NSCF run at k+q.
       ! gs_ham_kq has pointers to the *_kq arrays in output so we cannot deallocate them till the end.
       ! This is the reason why we use vlocal_kq although this term does not depend on k+q.
       call nscf%setup_kpt(spin, kq, istwfk_1, nband, cryst, dtset, psps, pawtab, pawfgr, &                         ! in
                           npw_kq, kg_kq, kpg_kq, ph3d_kq, kinpw_kq, ffnl_kq, vlocal_kq, cg_kq, gsc_kq, gs_ham_kq)  ! out

       ! Cache for u_{m k+q}(g).
       use_cg_kq = (my_iq > 1 .and. ucache_kq%use_cache)
       if (use_cg_kq) call ucache_kq%get_kpt(kq, istwfk_1, npw_kq, nspinor, nband, kg_kq, cg_kq)

       ! We can use cg_k as input for the NSCF for a very quick return.
       if (qq_is_gamma) then
         use_cg_kq = .True.; cg_kq = cg_k
       end if

       ! Compute u_{m k+q}(g)
       call nscf%solve_kpt(spin, kq, istwfk_1, nband, cryst, dtset, dtfil, gs_ham_kq, &
                           use_cg_kq, npw_kq, cg_kq, gsc_kq, eig_kq, msg, ierr)

       ABI_WARNING_IF(ierr /= 0, msg)
       tot_nscf_ierr = tot_nscf_ierr + ierr

       call ucache_kq%store_kpt(kq, istwfk_1, npw_kq, nspinor, nband, kg_kq, cg_kq)

       ! This to have the same gauge when qq = 0.
       if (qq_is_gamma) cg_kq = cg_k

       ! Make sure all procs in pert_comm have the same wavefunctions at k+q.
       call xmpi_bcast(cg_kq, master, pert_comm%value, ierr)
       if (psps%usepaw == 1) call xmpi_bcast(gsc_kq, master, pert_comm%value, ierr)

       ! Get phonons for this q-point.
       call ifc%fourq(cryst, qq, phfreqs, displ_cart, out_displ_red=displ_red_qq)
       phfreqs_eV = phfreqs * Ha_eV

       !if (my_ik == 1 .and. pert_comm%me == master) then
       if (my_ik == 1) then
         NCF_CHECK(nf90_put_var(ncid, vid("all_eigens_kq"), eig_kq, start=[1,iq,spin]))
         ! Write phonons for this q.
         if (spin == 1) then
           NCF_CHECK(nf90_put_var(ncid, vid("phfreqs"), phfreqs_ev, start=[1,iq]))
           NCF_CHECK(nf90_put_var(ncid, vid("phdispl_cart"), displ_cart, start=[1,1,1,1,iq]))
         end if
       end if

       ! For PAW, one has to solve a generalized eigenproblem.
       gen_eigenpb = psps%usepaw == 1; sij_opt = 0; if (gen_eigenpb) sij_opt = 1
       ABI_MALLOC(gs1c, (2, npw_kq*nspinor*((sij_opt+1)/2)))
       ABI_MALLOC(h1_kets_kq, (2, npw_kq*nspinor, nb_in_g))

       ! ====================================
       ! Get DFPT potentials for this q-point
       ! ====================================
       ! After this branch we have allocated v1scf(cplex, nfftf, nspden, my_npert)).
       if (need_ftinterp) then
         call dvdb%get_ftqbz(qq, cplex, nfftf, ngfftf, v1scf, pert_comm%value)
       else
         ! Read and reconstruct the dvscf potentials for qq and my_npert perturbations.
         call dvdb%readsym_qbz(cryst, qq, qmap_symrec(:,iq), cplex, nfftf, ngfftf, v1scf, pert_comm%value)
       end if

       ! Allocate vlocal1 with correct cplex. Note nvloc.
       ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_ham_kq%nvloc))

       ! Load the k/k+q dependent parts of the Hamiltonian.
       ! NB: In this routine we have to use gs_ham_k to have {k+q}_H0_k.
       ! Using gs_ham_kq would be wrong as it would lead to {k+q}_H0_{k+q}.

       call gs_ham_k%load_kprime(kpt_kp=kq, npw_kp=npw_kq, istwf_kp=istwfk_1, kg_kp=kg_kq, kpg_kp=kpg_kq, kinpw_kp=kinpw_kq, &
                                 ph3d_kp=ph3d_kq, ffnl_kp=ffnl_kq, compute_ph3d=.true., compute_gbound=.true.)
       !call gs_ham_k%print([std_out], "gs_ham_k after load", dtset%prtvol)

       ! Loop over my atomic perturbations: apply H1_{kappa, alpha} and compute gkq_atm.
       gkq_atm = zero
       do my_ip=1,my_npert
         idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

         ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
         call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_ham_kq%nvloc,&
                                    pawfgr, nscf%mpi_enreg, nscf%vtrial, v1scf(:,:,:,my_ip), vlocal, vlocal1)

         ! Prepare application of the NL part.
         call rf_ham_kq%init(cplex, gs_ham_k, ipert, has_e1kbsc=.true.)
         call rf_ham_kq%load_spin(spin, vlocal1=vlocal1, with_nonlocal=.true.)
         ! Load k-dependent part in the 1st-order Hamiltonian datastructure
         !call rf_ham_kq%load_k(npw_k=npw_k)

         ! Calculate dvscf * psi_k, results stored in h1_kets_kq on the k+q sphere.
         ! Compute H(1) applied to GS wavefunction Psi(0).
         do in_k=1,nb_in_g
           band_n = in_k + bstart - 1
           eig0nk = eig_k(band_n)
           ! Use scissor shift on 0-order eigenvalue.
           eshift = eig0nk - dtset%dfpt_sciss

           call getgh1c(berryopt0, cg_k(:,:,band_n), cwaveprj0, h1_kets_kq(:,:,in_k), &
                        grad_berry, gs1c, gs_ham_k, gvnlx1, idir, ipert, [eshift], nscf%mpi_enreg, ndat1, optlocal, &
                        optnl, opt_gvnlx1, rf_ham_kq, sij_opt, tim_getgh1c, usevnl)
         end do ! in_k

         ! Calculate <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation. No need to handle istwf_kq because it's always 1.
         !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(band_m)
         do in_k=1,nb_in_g
           !print *, "maxval(abs(h1_kets_kq(:,:,in_k))): ", maxval(abs(h1_kets_kq(:,:,in_k)))
           do im_kq=1,nb_in_g
             band_m = im_kq + bstart - 1
             gkq_atm(:, im_kq, in_k, ipc) = cg_zdotc(npw_kq*nspinor, cg_kq(:,:,band_m), h1_kets_kq(:,:,in_k))
           end do
         end do

       end do ! my_ip

       ! Collect gkq_atm inside pert_comm so that all procs can operate on the data.
       if (pert_comm%nproc > 1) call xmpi_sum(gkq_atm, pert_comm%value, ierr)

       ! From atom to phonon mode representation. Results stored in gkq_nu.
       call ephtk_gkknu_from_atm(nb_in_g, nb_in_g, 1, natom, gkq_atm, phfreqs, displ_red_qq, gkq_nu)
       !print *, "gkq_atm:", gkq_atm; print *, "displ_red_qq:", displ_red_qq; print *, "gkq_nu:", gkq_nu

       ! Write |g|^2 for this q.
       !if (pert_comm%me == master) then
       gkq2_nu = gkq_nu(1,:,:,:)**2 + gkq_nu(2,:,:,:)**2
       NCF_CHECK(nf90_put_var(ncid, vid("gkq2_nu"), gkq2_nu, start=[1,1,1,iq,ik,spin]))
       !end if

       ABI_FREE(gs1c)
       ABI_FREE(vlocal1)
       ABI_FREE(v1scf)
       ABI_FREE(vlocal_kq)
       ABI_FREE(ph3d_kq)
       ABI_FREE(kpg_kq)
       ABI_FREE(kinpw_kq)
       ABI_FREE(ffnl_kq)
       ABI_FREE(kg_kq)
       ABI_FREE(eig_kq)
       ABI_FREE(cg_kq)
       ABI_FREE(gsc_kq)
       ABI_FREE(h1_kets_kq)
       call gs_ham_kq%free(); call rf_ham_kq%free()
     end do ! my_iq

     ABI_FREE(vlocal)
     ABI_FREE(vlocal_k)
     ABI_FREE(ph3d_k)
     ABI_FREE(kpg_k)
     ABI_FREE(kinpw_k)
     ABI_FREE(ffnl_k)
     ABI_FREE(kg_k)
     ABI_FREE(eig_k)
     ABI_FREE(cg_k)
     ABI_FREE(gsc_k)
     call gs_ham_k%free()
   end do ! my_ik
 end do ! my_is

 NCF_CHECK(nf90_close(ncid))
 call xmpi_barrier(comm)

 ! ===========================================
 ! Write results to ab_out for automatic tests
 ! ===========================================
 call xmpi_sum(tot_nscf_ierr, comm, ierr)
 tot_nscf_ierr = int(tot_nscf_ierr / dble(pert_comm%nproc))

 if (my_rank == master) then
   if (tot_nscf_ierr == 0) then
     call wrtout(units, &
       sjoin("Computation of g(k,q) completed. All NSCF runs converged within tolwfr: ", ftoa(dtset%tolwfr)), pre_newlines=1)
   else
     msg = sjoin("WARNING:", itoa(tot_nscf_ierr), "NSCF runs did not converge within tolwfr: ", ftoa(dtset%tolwfr), ". Increase nstep!")
     call wrtout(ab_out, msg)
     ABI_WARNING(msg)
   end if

   NCF_CHECK(nctk_open_read(ncid, gpath_path, xmpi_comm_self))

   ! Write k/q wavevectors.
   call wrtout(units, "kpoints:")
   do ik=1, nk_path
     call wrtout(units, sjoin(char(9), itoa(ik), ktoa(kpath%points(:,ik))))
   end do
   call wrtout(units, "qpoints:")
   do iq=1, nq_path
     call wrtout(units, sjoin(char(9), itoa(iq), ktoa(qpath%points(:,iq))))
   end do

   ! Write KS eigenvalues.
   if (nq_path > 1) then
     ABI_MALLOC(eig_kq, (nband))
     do spin=1,nsppol
       call wrtout(units, sjoin(" Energies_kq in eV for spin:", itoa(spin)))
       do iq=1, nq_path
         if (all(iq /= [1, 2, nq_path-1, nq_path])) cycle
         NCF_CHECK(nf90_get_var(ncid, vid("all_eigens_kq"), eig_kq, start=[1,iq,spin]))
         do ii=0,(nband-1)/8
           write(msg, '(a, 8es16.6)' )' ene:',(eig_kq(band_m) * Ha_eV, band_m=1+ii*8,min(nband,8+ii*8))
           call wrtout(units, msg)
         end do
       end do
     end do
     ABI_FREE(eig_kq)

   else
     ! nk_path > 1
     ABI_MALLOC(eig_k, (nband))
     do spin=1,nsppol
       call wrtout(units, sjoin(" Energies_k in eV for spin:", itoa(spin)))
       do ik=1, nk_path
         NCF_CHECK(nf90_get_var(ncid, vid("all_eigens_k"), eig_k, start=[1,ik,spin]))
         if (all(ik /= [1, 2, nk_path-1, nk_path])) cycle
         do ii=0,(nband-1)/8
           write(msg, '(a, 8es16.6)' )' ene:',(eig_k(band_n) * Ha_eV, band_n=1+ii*8,min(nband,8+ii*8))
           call wrtout(units, msg)
         end do
       end do
     end do
     ABI_FREE(eig_k)
   end if

   ! Save g^2 in nc format without any average. This operation will be performed by AbiPy (need ph freqs and eigenergies)
   call wrtout(units, " Writing sqrt(1/N_b^2 \sum_{mn} |g_{mn,nu}(k, q)|^2) in meV for testing purpose.", pre_newlines=2)
   write(msg, "(1x,4(a5,1x),a16)") "nu","iq", "ik", "spin", "|g| in meV"
   call wrtout(units, msg)

   do spin=1,nsppol
     do ik=1, nk_path
       if (all(ik /= [1, 2, nk_path-1, nk_path])) cycle
       do iq=1, nq_path
         if (all(iq /= [1, 2, nq_path-1, nq_path])) cycle
         NCF_CHECK(nf90_get_var(ncid, vid("gkq2_nu"), gkq2_nu, start=[1,1,1,iq,ik,spin]))
         do nu=1,natom3
           write(msg, "(1x,4(i5,1x),es16.6)") nu, iq, ik, spin, sqrt(sum(gkq2_nu(:,:, nu)) / nb_in_g**2) * Ha_meV
           call wrtout(units, msg)
         end do
       end do ! iq
     end do ! ik
   end do ! spin

   NCF_CHECK(nf90_close(ncid))
 end if ! master

 ! Free memory.
 ABI_FREE(my_ik_inds)
 ABI_FREE(my_iq_inds)
 ABI_FREE(gvnlx1)
 ABI_FREE(grad_berry)
 ABI_FREE(qmap_symrec)
 ABI_FREE(gkq_atm)
 ABI_FREE(gkq_nu)
 ABI_FREE(gkq2_nu)
 ABI_FREE(displ_cart)
 ABI_FREE(displ_red_qq)
 ABI_FREE(my_spins)

 call pawcprj_free(cwaveprj0)
 ABI_FREE(cwaveprj0)
 do my_is=1,my_nspins
   call comm_my_is(my_is)%free()
 end do
 ABI_FREE(comm_my_is)

 call qpath%free(); call kpath%free(); call ucache_k%free(); call ucache_kq%free()
 call qpt_comm%free(); call kpt_comm%free(); call pert_comm%free(); call nscf%free()

 call cwtime_report(" eph_path: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)
 !stop

contains
 integer function vid(var_name)
   character(len=*),intent(in) :: var_name
   vid = nctk_idname(ncid, var_name)
 end function vid

end subroutine eph_path_run
!!***

end module m_eph_path
!!***

!!****m* ABINIT/m_wqk
!! NAME
!!  m_wqk
!!
!! FUNCTION
!!  Compute
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

module m_wqk

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 USE_MPI
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_wfk
 use m_fft
 use m_wfd
 use m_krank
 use m_sort
 use m_hdr
 use netcdf
 use m_nctk
 use m_dtset
 use m_dtfil

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : pseudopotential_type
 use m_gwdefs,         only : GW_Q0_DEFAULT
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, print_arr
 use m_io_tools,       only : iomode_from_fname, file_exists, is_open, open_file, flush_unit
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex, print_ngfft
 use m_cgtk,           only : cgtk_rotate
 use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map
 use m_bz_mesh,        only : isamek, kmesh_t
 use m_gsphere,        only : gsphere_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
 use m_io_screening,   only : hscr_t, get_hscr_qmesh_gsph, read_screening
 use m_vcoul,          only : vcoul_t
 use m_occ,            only : get_fact_spin_tol_empty
 use m_ebands,         only : ebands_t
 use m_pstat,          only : pstat_proc

 implicit none

 private
!!***

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif

 public :: wqk_run  ! Main entry point to compute wqk matrix elements.

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_wqk/wqk_run
!! NAME
!!  wqk_run
!!
!! FUNCTION
!!  Compute W_qk
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
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

subroutine wqk_run(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, wfk_hdr, &
                   pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),target,intent(in) :: ebands
 type(hdr_type),intent(in) :: wfk_hdr
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18), ngfftf(18)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: LOG_MODQ = 1, LOG_MODK = 4, LOG_MODP = 4
 integer,parameter :: istwfk1 = 1, master = 0, ndat1 = 1, cplex1 = 1, pawread0 = 0
 integer :: band, nkpt, my_rank, nsppol, iq_ibz, iq_bz
 integer :: cplex,nspinor,nprocs, ii, ib, my_is, spin, npw_c
 integer :: isym_q, trev_q
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_k_ibz, istwf_k_ibz
 integer :: ikq_ibz, isym_kq, trev_kq, npw_kq, istwf_kq,  npw_kq_ibz, istwf_kq_ibz
 integer :: mpw,ierr,nqbz,ncerr !,spad
 integer :: n1,n2,n3,n4,n5,n6,nspden, mqmem, m_kq, n_k, restart, root_ncid, spin_ncid
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg_k,nkpg_kq,cnt, ndone, nmiss
 integer :: my_ipp, ipp_bz, nqlwl, scr_iomode
 integer :: qptopt, my_iq, my_ik, qbuf_size, iqbuf_cnt, nb
 real(dp) :: cpu_all, wall_all, gflops_all, cpu_qq, wall_qq, gflops_qq, cpu_kk, wall_kk, gflops_kk
 real(dp) :: ecut
 logical :: isirr_k, isirr_kq, qq_is_gamma, isirr_q
 logical :: print_time_qq, print_time_kk
 type(wfd_t) :: wfd
 type(crystal_t) :: den_cryst
 type(hdr_type) :: den_hdr
 type(kmesh_t) :: qmesh, kmesh
 type(gsphere_t),target :: gsph_c
 type(hscr_t),target :: hscr
 type(vcoul_t) :: vcp
 character(len=fnlen) :: screen_filepath
 character(len=5000) :: msg, qkp_string
!arrays
 integer :: nbsum, my_bsum_start(dtset%nsppol), my_bsum_stop(dtset%nsppol), my_nbsum(dtset%nsppol)
 integer :: g0_k(3), g0_kq(3), units(2), work_ngfft(18), gmax(3)
 integer :: mapl_k(6), mapl_kq(6), mapc_qq(6)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), my_pp_inds(:)
 integer,allocatable :: gbound_k(:,:), gbound_kq(:,:), gbound_c(:,:), nband(:,:), wfd_istwfk(:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:)
 !integer(i1b),allocatable :: itreat_qibz(:)
 !integer, ABI_CONTIGUOUS pointer :: kg_c(:,:)
 real(dp) :: kk(3),kq(3),kk_ibz(3),kq_ibz(3), kqmp(3), kmp(3), pp(3), kmp_ibz(3), kqmp_ibz(3), qq_ibz(3), qpt(3)
 !complex(gwpc) :: ctmp_gwpc, xdot_tmp
!arrays
 real(dp),allocatable :: qlwl(:,:), kpg_k(:,:), kpg_kq(:,:), cg_work(:,:), ug_k(:,:), ug_kq(:,:)
 real(dp),allocatable :: work(:,:,:,:), my_gbuf(:,:,:,:,:,:)
 complex(gwpc),allocatable :: cwork_ur(:), rhotwg_c(:), vc_sqrt_gc(:), ur_nk(:,:), ur_mkq(:,:), epsm1_ggw(:,:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm); units = [std_out, ab_out]
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 nsppol = ebands%nsppol; nspinor = ebands%nspinor; nspden = dtset%nspden
 ecut = dtset%ecut ! dtset%dilatmx

 ! Open WQK.nc file and go to data mode.
 NCF_CHECK(nctk_open_modify(root_ncid, strcat(dtfil%filnam_ds(4), "_WQK.nc") , comm))
 NCF_CHECK(nctk_set_datamode(root_ncid))

 !call wrtout(units, sjoin("- Number of q-points/spin completed:", itoa(count(done_qbz_spin == 1)), "/", itoa(sigma%nkcalc)))

 ! ================
 ! HANDLE SCREENING
 ! ================
 ! Init gsph_c for the correlated part.
 screen_filepath = dtfil%fnameabi_scr
 ABI_CHECK(dtfil%fnameabi_scr /= ABI_NOFILE, "SCR file must be specified")

 call kmesh%init(cryst, wfk_hdr%nkpt, wfk_hdr%kptns, dtset%kptopt)

 ! Read g-sphere for correlation and qmesh from SCR file.
 call get_hscr_qmesh_gsph(screen_filepath , dtset, cryst, hscr, qmesh, gsph_c, qlwl, comm)
 call hscr%print(units, dtset%prtvol, header="Header of the SCR file")

 nqlwl = size(qlwl, dim=2)
 if (nqlwl == 0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i0,a,3f9.6)')&
     "The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
     "Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)
   ABI_COMMENT(msg)
 end if

 ! TODO:
 ! Here we sort the qmesh by stars so that we can split the pp wavevectors in blocks and therefore
 ! reduce the number of wavevectors in the IBZ that must be stored in memory.
 !call qmesh%pack_by_stars()

 ! Distribute the sum over pp wavevectors inside pp_sum_comm using block distribution.
 !my_pp_start_spin = -1; my_pp_stop_spin = 0
 !do my_is=1,gstore%my_nspins
 !  spin = gstore%my_spins(my_is); gqk => gstore%gqk(my_is)
 !  call xmpi_split_block(qmesh%nbz, gqk%pp_sum_comm%value, my_npp(spin), my_pp_inds)
 !  if (my_npp(spin) > 0) then
 !    my_pp_start_spin(spin) = my_pp_inds(1); my_pp_stop_spin(spin) = my_pp_inds(my_npp(spin))
 !  end if
 !  ABI_SFREE(my_pp_inds)
 !end do ! my_is

 ! Initialize Coulomb term on the IBZ of the qmesh. Use largest G-sphere.
 npw_c = gsph_c%ng
 call vcp%init(gsph_c, cryst, qmesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_c%ng, &
                 nqlwl, qlwl, comm)
 ABI_FREE(qlwl)

 if (my_rank == master)  then
   call kmesh%print(units, header="K-mesh for wavefunctions", prtvol=dtset%prtvol)
   call gsph_c%print(units, dtset%prtvol, header="G-sphere for correlation")
   call vcp%print(units, prtvol=dtset%prtvol)
 end if
 call kmesh%free()

 ABI_CHECK_IGE(npw_c, 1, "npw_c <= 1")

 ! Initialize the wave function descriptor.
 ! Each node has all k-points and spins and bands between my_bsum_start and my_bsum_stop
 ! TODO: One can exploit qq, kk and pp parallelism to find the wavevectors in the IBZ
 ! that will be needed in the loops below and allocate only these wavevectors so that memory scales.

 !nbsum = dtset%mband
 !my_bsum_start = 1; my_bsum_stop = nbsum; my_nbsum = my_bsum_stop - my_bsum_start + 1

 nkpt = wfk_hdr%nkpt
 ABI_MALLOC(nband, (nkpt, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkpt, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkpt, nsppol))

 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.
 !bks_mask(my_bsum_start:my_bsum_stop,:,:) = .True.
 ! Distribute wavefunctions according to the set of kk, qq and pp wavevectors treated by this MPI proc.

 ! Also, compute mpw and gmax including the additional pp
 ! This is the maximum number of PWs for all possible k+q treated.
 !call gstore%fill_bks_mask_qmesh(ecut, dtset%mband, nkpt, nsppol, my_pp_start_spin, my_pp_stop_spin, qmesh, &
 !                                  my_bsum_start, my_bsum_stop, bks_mask, mpw, gmax)

 !mpw ??

 ! Init work_ngfft
 gmax = gmax + 4 ! FIXME: this is to account for umklapp, should also consider Gamma-only and istwfk
 gmax = 2*gmax + 1

 !if (dtset%userie == 124) then
 !  ! Debugging section have all states on each MPI rank.
 !  bks_mask = .True.; call wrtout(std_out, " Storing all bands for debugging purposes.")
 !end if

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkpt))
 wfd_istwfk = 1

 call wfd%init(cryst, pawtab, psps, keep_ur, dtset%mband, nband, nkpt, nsppol, bks_mask,&
               nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

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
 call wfd%print(units, header="Wavefunctions for GWPT calculation.")

 ABI_MALLOC(gbound_k, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_kq, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_c, (2*mgfft+8, 2))
 ABI_MALLOC(cg_work, (2, mpw*nspinor))

 ! Allocate workspace arrays.
 ! Find correspondence IBZ --> set of q-points in DVDB.

 !ABI_ICALLOC(itreat_qibz, (gstore%nqibz))
 !itreat_qibz = 1
 !call wrtout(std_out, sjoin("P Number of q-points in the IBZ treated by this proc: " ,itoa(count(itreat_qibz == 1))))
 !ABI_FREE(itreat_qibz)

 ! Initialize plasmon-pole object.
 mqmem = qmesh%nibz

 ! Read symmetrized em1 from file and build ppmodel parameters.
 ! TODO: MPI-shared memory + compute my set of pp-vectors on ppm%new_setup
 scr_iomode = iomode_from_fname(screen_filepath)
 ABI_MALLOC(epsm1_ggw, (npw_c, npw_c, hscr%nomega))

 do iq_ibz=1,qmesh%nibz
   call read_screening("inverse_dielectric_function", screen_filepath, &
                       npw_c, 1, hscr%nomega, epsm1_ggw, scr_iomode, comm, iqiA=iq_ibz)
 end do ! iq_ibz

 ABI_FREE(epsm1_ggw)
 call hscr%free()

 ! Allocate g-vectors centered on k, k+q, k-p, and k+q-p.
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))
 ! Spherical Harmonics for useylm == 1.
 ! FIXME: These arrays should be allocated with npw_k, npw_kq
 ! but should recheck the API used to symmetrized wavefunctions.
 ! GS wavefunctions
 ABI_MALLOC(cwork_ur, (nfft*nspinor))

 call pstat_proc%print(_PSTAT_ARGS_)

 ! This parameter defines the size of the q-buffer used to store the g(k, q) e-ph matrix elements
 ! for all the k-point treated by this MPI rank.
 ! Increasing the buffer size increases the memory requirements
 ! but it leads to better performance as the number of IO operations is decreased.
 ! TODO: Should compute it on the basis of my_nkpt and my_nqpt
 qbuf_size = 1
 !qbuf_size = 16
 call wrtout(std_out, sjoin(" Begin computation of GWPT e-ph matrix elements with qbuf_size:", itoa(qbuf_size)), pre_newlines=1)

 ! ===============================
 ! Loop over MPI distributed spins
 ! ===============================

#if 0
 do my_is=1,gstore%my_nspins
   !NCF_CHECK(nf90_inq_ncid(root_ncid, strcat("gqk", "_spin", itoa(spin)), spin_ncid))

   ! TODO: Should introduce the possibility of specifying different nb states
   ! for the incoming and the intermediate states.
   nb = gqk%nb
   ABI_MALLOC(iq_buf, (2, qbuf_size))

   ABI_MALLOC_OR_DIE(ur_nk,  (nfft*nspinor, gqk%bstart:gqk%bstop), ierr)
   ABI_MALLOC_OR_DIE(ur_mkq, (nfft*nspinor, gqk%bstart:gqk%bstop), ierr)
   ABI_MALLOC_OR_DIE(my_gbuf, (gqk%cplex, nb, nb, natom3, gqk%my_nk, qbuf_size), ierr)

   ! ============================================================
   ! Loop over MPI distributed q-points in Sigma_q (gqk%qpt_comm)
   ! ============================================================
   do my_iq=1,gqk%my_nq
     qq_is_gamma = sum(qpt**2) < tol14
     iq_bz = gqk%my_q2bz(my_iq)

     ! Handle possible restart.
     !if (done_qbz_spin(iq_bz, spin) == 1) then
     !  call wrtout(std_out, sjoin(" iq_bz:", itoa(iq_bz), ", spin: ", itoa(spin), " already computed --> skipping iteration"))
     !  cycle
     !end if

     ! Note symrec conventions here as needed to symmetrize DFPT potentials.
     iq_ibz = gqk%my_q2ibz(1, my_iq) ; isym_q = gqk%my_q2ibz(2, my_iq)
     trev_q = gqk%my_q2ibz(6, my_iq) !; g0_q = gqk%my_q2ibz(3:5, my_iq)
     isirr_q = (isym_q == 1 .and. trev_q == 0)
     !isirr_q = (isym_q == 1 .and. trev_q == 0 .and. all(g0_q == 0))
     !tsign_q = 1; if (trev_q == 1) tsign_q = -1
     !qq_ibz = gstore%qibz(:, iq_ibz)
     mapc_qq = gqk%my_q2ibz(:, my_iq)

     print_time_qq = my_rank == 0 .and. (my_iq <= LOG_MODQ .or. mod(my_iq, LOG_MODQ) == 0)
     if (print_time_qq) then
       call cwtime(cpu_qq, wall_qq, gflops_qq, "start")
       call inds2str(0, sjoin(" Computing g^Sigma(k, q) for qpt:", ktoa(qpt)), my_iq, gqk%my_nq, gqk%glob_nq, msg)
       call wrtout(std_out, sjoin(msg, ", for spin:", itoa(spin)), pre_newlines=1)
     end if

     !if (isirr_q) then
     !else
     !end if

     ! ==================
     ! Loop over k-points
     ! ==================

     do my_ik=1,gqk%my_nk
       ! NB: All procs in gqk%pert_comm and gqk%bsum_com and gqk%pp_sum_comm enter this section.
       iqbuf_cnt = 1 + mod(my_iq - 1, qbuf_size)
       iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]

       ! Set entry to zero. Important as there are cycle instructions inside these loops
       ! and we don't want random numbers written to disk.
       my_gbuf(:,:,:,:, my_ik, iqbuf_cnt) = zero

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

       !if (kpts_map("symrel", ebands%kptopt, cryst, gstore%krank_ibz, 1, kq, mapl_kq) /= 0) then
       !  write(msg, '(4a)' )"k-mesh is not closed!",ch10, "k+q could not be generated from a symmetrical one.",trim(ltoa(kq))
       !  ABI_ERROR(msg)
       !end if
       ikq_ibz = mapl_kq(1); isym_kq = mapl_kq(2); trev_kq = mapl_kq(6); g0_kq = mapl_kq(3:5)
       isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
       kq_ibz = ebands%kptns(:, ikq_ibz)
       istwf_kq_ibz = wfd%istwfk(ikq_ibz); npw_kq_ibz = wfd%npwarr(ikq_ibz)

       ! Get npw_kq, kg_kq for k+q.
       call wfd%get_gvec_gbound(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, dtset%nloalg, &  ! in
                                istwf_kq, npw_kq, kg_kq, nkpg_kq, kpg_kq, gbound_kq)      ! out

       ABI_MALLOC(ug_k, (2, npw_k*nspinor))
       ABI_MALLOC(ug_kq, (2, npw_kq*nspinor))

       ! Precompute ur_nk and ur_mkq for all m and n indices treated by me.
       ! TODO: Can distribute operations inside gqk%pert_comm
       do n_k=gqk%bstart, gqk%bstop ! do n_k=gqk%n_start, gqk%n_stop
         call wfd%rotate_cg(n_k, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                            cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_nk(:,n_k))
       end do

       do m_kq=gqk%bstart, gqk%bstop ! do m_kq=gqk%m_start, gqk%m_stop
         call wfd%rotate_cg(m_kq, ndat1, spin, kq_ibz, npw_kq, kg_kq, istwf_kq, &
                            cryst, mapl_kq, gbound_kq, work_ngfft, work, ug_kq, urs_kbz=ur_mkq(:,m_kq))
       end do

       ! Be careful here because pp should run over the list of wavevectors in the screening matrix!
       ! as qmesh%bz is not necessarily equivalent to the k-mesh for the wavefunctions,
       ! and we have to use the ipp_bz index to symmetrize W(pp_bz) from W(pp_ibz).
       !
       ! TODO: Should order nbz in shells so that one can reduce the memory required
       ! to store W(pp) if pp_parallelism is activated.

       ABI_FREE(kpg_k)
       ABI_FREE(kpg_kq)
       ABI_FREE(ug_k)
       ABI_FREE(ug_kq)

       ! Dump buffer
       !if (iqbuf_cnt == qbuf_size) call dump_my_gbuf()

       if (print_time_kk) then
         call inds2str(3, "My k-point", my_ik, gqk%my_nk, gqk%glob_nk, msg)
         call cwtime_report(msg, cpu_kk, wall_kk, gflops_kk); if (my_ik == LOG_MODK) call wrtout(std_out, "...", do_flush=.True.)
       end if
     end do ! my_ik

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
 end do ! my_is

 call cwtime_report(" vkk calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Set vqk_completed to 1 so that we can easily check if restarted is needed.
 !if (my_rank == master) then
 !  NCF_CHECK(nf90_put_var(root_ncid, root_vid("vqk_completed"), 1))
 !end if
 NCF_CHECK(nf90_close(root_ncid))
 call xmpi_barrier(comm)

 ! Output some of the results to ab_out for testing purposes
 !call gstore%print_for_abitests(dtset)

 ! Free memory
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(cg_work)
 ABI_FREE(cwork_ur)
 ABI_FREE(work)
 ABI_FREE(gbound_k)
 ABI_FREE(gbound_kq)
 ABI_FREE(gbound_c)
 ABI_FREE(done_qbz_spin)

 call wfd%free(); call vcp%free()
 call qmesh%free(); call gsph_c%free()

 call xmpi_barrier(comm) ! This to make sure that the parallel output of GSTORE is completed
 call cwtime_report(" wqk_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)
#endif

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

 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! FIXME: Recheck this part as we have way more levels of parallelism in GWPT
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !if (gqk%coords_qkpb_sumbp(3) /= 0) goto 10 ! Yes, I'm very proud of this GOTO.
 !if (gqk%pert_ppsum_bsum_comm%me /= 0) goto 10 ! Yes, I'm very proud of this GOTO.

 !iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]
 !my_iq = iq_buf(1, 1)
 !iq_glob = my_iq + gqk%my_qstart - 1

 !print *, "in dump_my_gbuf with start: ", [1, 1, 1, 1, gqk%my_kstart, iq_glob]
 !print *, "                  count; ", [gqk%cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%my_nk, iqbuf_cnt]
 !print *, "my_gbuf", my_gbuf(:,:,:,natom3,1,1)

 ! NB: this is an individual IO operation
 !ncerr = nf90_put_var(spin_ncid, spin_vid("gvals"), my_gbuf, &
 !                     start=[1, 1, 1, 1, gqk%my_kstart, iq_glob], &
 !                     count=[gqk%cplex, gqk%nb, gqk%nb, gqk%natom3, gqk%my_nk, iqbuf_cnt])
 !NCF_CHECK(ncerr)

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

end subroutine wqk_run
!!***

end module m_wqk
!!***

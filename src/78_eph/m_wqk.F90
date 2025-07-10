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
 use m_sort
 use m_hdr
 use netcdf
 use m_nctk
 use m_dtset
 use m_dtfil
 use m_fstab
 use m_ephtk

 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : pseudopotential_type
 use m_gwdefs,         only : GW_Q0_DEFAULT, cone_gw, czero_gw
 use m_hide_blas,      only : xgemm, xdotc
 use m_special_funcs,  only : gaussian
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
 use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
 use m_numeric_tools,  only : arth, linspace !, get_diag, print_arr
 use m_io_tools,       only : iomode_from_fname
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, print_ngfft
 !use m_krank,          only : krank_t, krank_from_kptrlatt
 !use m_cgtk,           only : cgtk_rotate
 !use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm
 use m_crystal,        only : crystal_t
 use m_kpts,           only : kpts_ibz_from_kptrlatt, kpts_timrev_from_kptopt, kpts_map
 use m_bz_mesh,        only : isamek, kmesh_t
 use m_gsphere,        only : gsphere_t
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
 use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
 use m_io_screening,   only : hscr_t, get_hscr_qmesh_gsph !, read_screening
 use m_vcoul,          only : vcoul_t
 !use m_occ,            only : get_fact_spin_tol_empty
 use m_ebands,         only : ebands_t, edos_t
 use m_pstat,          only : pstat_proc
 use m_sigtk,          only : sigtk_multiply_by_vc_sqrt
 use m_screening,      only : epsm1_t
 use m_ddk,            only : ddkop_t, ddkop_new

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
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!
!! SOURCE

subroutine wqk_run(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, wfk_hdr, &
                   pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 character(len=*),intent(in) :: wfk0_path
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),target,intent(in) :: ebands
 type(hdr_type),intent(in) :: wfk_hdr
 type(pseudopotential_type),intent(in) :: psps
 type(mpi_type),intent(inout) :: mpi_enreg
 integer,intent(in) :: comm
!arrays
 integer,intent(in) :: ngfft(18), ngfftf(18)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: LOG_MODQ = 1, LOG_MODK = 4, istwfk1 = 1, master = 0, ndat1 = 1
 integer :: id_required, approx_type, ikxc, option_test, nkxc, ig
 integer :: nkibz, my_rank, nsppol, iq_ibz, iq_bz, isym_qq, itim_qq
 integer :: cplex,nspinor,nprocs, ii, ib, my_is, spin, max_npw_xc, min_npw_xc,  npw_x, npw_c, bstart_k, bstop_k, nband_k, ik_bz
 integer :: isym_q, trev_q, bstart_kq, bstop_kq, nband_kq, ncols, mpw, ierr, ncerr !,spad
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_k_ibz, istwf_k_ibz
 integer :: ikq_ibz, isym_kq, trev_kq, npw_kq, istwf_kq,  npw_kq_ibz, istwf_kq_ibz
 integer :: m_kq, n_k, root_ncid, ib_k, ib_kq, prtwqk, ne ! n1,n2,n3,n4,n5,n6,
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg_k,nkpg_kq,cnt, edos_intmeth, nqlwl, scr_iomode
 integer :: ikfs_bz, ikq_fs, my_ik, my_nk, qbuf_size, iqbuf_cnt !, my_iq,
 real(dp) :: cpu_all, wall_all, gflops_all, cpu, wall, gflops, cpu_qq, wall_qq, gflops_qq, cpu_kk, wall_kk, gflops_kk
 real(dp) :: edos_step, edos_broad, ecut, e_nk, e_mkq, e_min, e_max, de
 logical :: isirr_k, isirr_kq, qq_is_gamma, isirr_q, remove_exchange, print_time_qq, print_time_kk
 type(wfd_t) :: wfd
 type(kmesh_t) :: qmesh, kmesh
 type(gsphere_t),target :: gsph_x, gsph_c
 type(hscr_t),target :: hscr
 type(vcoul_t) :: vcp
 type(edos_t) :: edos
 type(epsm1_t) :: epsm1
 type(ddkop_t) :: ddkop
 !type(krank_t) :: krank_ibz
 character(len=fnlen) :: screen_filepath
 character(len=5000) :: msg, qkp_string
!arrays
 integer :: g0_k(3), g0_kq(3), units(2), work_ngfft(18), gmax(3), mapl_k(6), mapl_kq(6)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), gbound_k(:,:), gbound_kq(:,:), gbound_c(:,:), gbound_x(:,:)
 integer,allocatable :: iq_buf(:,:), done_qbz_spin(:,:), kbz2ibz(:,:), nband(:,:), wfd_istwfk(:)
 !integer(i1b),allocatable :: itreat_qibz(:)
 integer, ABI_CONTIGUOUS pointer :: kg_c(:,:), kg_x(:,:)
 real(dp) :: kk(3), kq(3), kk_ibz(3), kq_ibz(3), qq_bz(3) ! qq_ibz(3),
 complex(gwpc) :: ctmp_gwpc, mu
!arrays
 real(dp) :: n0(ebands%nsppol)
 real(dp),allocatable :: qlwl(:,:), kpg_k(:,:), kpg_kq(:,:), ug_k(:,:,:), ug_kq(:,:), cg_work(:,:)
 real(dp),allocatable :: work(:,:,:,:), e_mesh(:), vcart_ibz(:,:,:,:) !, my_gbuf(:,:,:,:,:,:)
 complex(gwpc),allocatable :: cwork_ur(:), rhotwg_x(:,:,:),rhotwg_c(:,:,:), w_rhotwg_c(:,:,:),vc_sqrt_gx(:), ur_nk(:,:), ur_mkq(:,:)
 complex(gwpc),allocatable :: kxcg(:,:)
 complex(dp),allocatable :: w_ee(:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
 type(fstab_t),target,allocatable :: fstab(:)
!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
 end if
 if (dtset%nsppol == 2) then
   ABI_ERROR("nsppol implemented")
 end if

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm); units = [std_out, ab_out]
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 nsppol = ebands%nsppol; nspinor = ebands%nspinor; ecut = dtset%ecut ! dtset%dilatmx

 ! Compute electron DOS.
 edos_intmeth = 2; if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 edos_step = 0.01 * eV_Ha; edos_broad = 0.3 * eV_Ha
 edos = ebands%get_edos(cryst, edos_intmeth, edos_step, edos_broad, comm)

 ! Get DOS per spin channel
 n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos%print(unit=ab_out)
   call edos%print(unit=std_out)
   !path = strcat(dtfil%filnam_ds(4), "_EDOS")
   !call wrtout(ab_out, sjoin("- Writing electron DOS to file:", path, ch10))
   !call edos%write(path)
 end if

 ! Find Fermi surface k-points
 ! TODO: support kptopt, change setup of k-points if tetra: fist tetra weights then k-points on the Fermi surface!
 ABI_MALLOC(fstab, (nsppol))
 call fstab_init(fstab, ebands, cryst, dtset, comm)
 if (my_rank == master) then
   call fstab_print(fstab, unit=std_out)
   call fstab_print(fstab, unit=ab_out)
 end if

 spin = 1
 associate (fs => fstab(spin))
 ! Build linear mesh for W(e,e').
 e_min = minval(ebands%eig(fs%bmin,:,spin)) - 0.1_dp * eV_Ha
 e_max = maxval(ebands%eig(fs%bmax,:,spin)) + 0.1_dp * eV_Ha
 de = 0.002 * eV_Ha
 ne = 1 + (e_max - e_min) / de
 ABI_CALLOC(w_ee, (ne, ne))
 ABI_MALLOC(e_mesh, (ne))
 e_mesh = linspace(e_min, e_max, ne)
 end associate

 ! ================
 ! HANDLE SCREENING
 ! ================
 ! Init gsph_c for the correlated part.
 screen_filepath = dtfil%fnameabi_scr
 ABI_CHECK(dtfil%fnameabi_scr /= ABI_NOFILE, "SCR file must be specified")

 ! Read g-sphere for correlation and qmesh from SCR file.
 call get_hscr_qmesh_gsph(screen_filepath, dtset, cryst, hscr, qmesh, gsph_c, qlwl, comm)
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

 ! Init g-sphere for the exchange part from ecutsigx.
 call gsph_c%extend(cryst, dtset%ecutsigx, gsph_x)

 ! TODO:
 ! Here we sort the pp_mesh by stars so that we can split the pp wavevectors in blocks and therefore
 ! reduce the number of wavevectors in the IBZ that must be stored in memory.
 !call qmesh%pack_by_stars()

 !ABI_ICALLOC(itreat_qibz, (gstore%nqibz))
 !itreat_qibz = 1
 !call wrtout(std_out, sjoin("P Number of q-points in the IBZ treated by this proc: " ,itoa(count(itreat_qibz == 1))))
 !ABI_FREE(itreat_qibz)

 ! Initialize Coulomb term on the IBZ of the qmesh. Use largest G-sphere.
 call kmesh%init(cryst, wfk_hdr%nkpt, wfk_hdr%kptns, dtset%kptopt)

 npw_x = gsph_x%ng; npw_c = gsph_c%ng
 min_npw_xc = min(npw_x, npw_c)
 max_npw_xc = max(npw_x, npw_c)
 if (gsph_x%ng >= gsph_c%ng) then
   call vcp%init(gsph_x, cryst, qmesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_x%ng, &
                 nqlwl, qlwl, comm)
 else
   call vcp%init(gsph_c, cryst, qmesh, kmesh, dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, gsph_c%ng, &
                 nqlwl, qlwl, comm)
 end if

 ABI_FREE(qlwl)

 if (my_rank == master)  then
   call kmesh%print(units, header="K-mesh for wavefunctions", prtvol=dtset%prtvol)
   call gsph_x%print(units, dtset%prtvol, header="G-sphere for exchange")
   call gsph_c%print(units, dtset%prtvol, header="G-sphere for correlation")
   call vcp%print(units, prtvol=dtset%prtvol)
 end if
 call kmesh%free()

 ABI_CHECK_IGE(npw_x, 1, "npw_x <= 1")
 ABI_CHECK_IGE(npw_c, 1, "npw_c <= 1")

 ! Note that in this case, the sphere is always Gamma-centered i.e. it does not depend on the qq wavevector
 kg_c => gsph_c%gvec(:, 1:npw_c)
 kg_x => gsph_x%gvec(:, 1:npw_x)

 ! Initialize the wave function descriptor.
 ! Only wavefunctions on the FS are stored in wfd.
 ! Need all k-points on the FS because of k+q, spin is not distributed for the time being.
 ! It would be possible to reduce the memory allocated per MPI-rank via OpenMP.

 nkibz = wfk_hdr%nkpt
 ABI_MALLOC(nband, (nkibz, nsppol))
 ABI_MALLOC(bks_mask, (dtset%mband, nkibz, nsppol))
 ABI_MALLOC(keep_ur, (dtset%mband, nkibz, nsppol))
 nband = dtset%mband; bks_mask = .False.; keep_ur = .False.

 do spin=1,ebands%nsppol
   associate (fs => fstab(spin))
   do ik_bz=1,fs%nkfs
     ik_ibz = fs%indkk_fs(1, ik_bz)
     bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz)
     bks_mask(bstart_k:bstart_k+nband_k-1, ik_ibz, spin) = .True.
   end do
   end associate
 end do

 ! Impose istwfk=1 for all k points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 ABI_MALLOC(wfd_istwfk, (nkibz))
 wfd_istwfk = 1

 call wfd%init(cryst, pawtab, psps, keep_ur, dtset%mband, nband, nkibz, nsppol, bks_mask,&
               dtset%nspden, nspinor, ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
               dtset%nloalg, dtset%prtvol, dtset%pawprtvol, comm)

 call pstat_proc%print(_PSTAT_ARGS_)

 ABI_FREE(nband)
 ABI_FREE(bks_mask)
 ABI_FREE(keep_ur)
 ABI_FREE(wfd_istwfk)

 ! Read wavefunctions.
 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

 ! TODO
 ! mpw is the maximum number of plane-waves over k and k+q where k and k+q are in the BZ.
 ! we also need the max components of the G-spheres (k, k+q) in order to allocate the workspace array work
 ! that will be used to symmetrize the wavefunctions in G-space.
 ! This is the maximum number of PWs for all possible k+q treated.

 call ephtk_get_mpw_gmax(nkibz, wfk_hdr%kptns, ecut, cryst%gmet, mpw, gmax, comm)

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)
#if 0
 call cwtime(cpu, wall, gflops, "start", msg=" Computing v_nk matrix elements for all states on the FS...")
 ii = huge(1); jj = -1
 do spin=1,nsppol
   ii = min(ii, fstab(spin)%bmin)
   jj = max(jj, fstab(spin)%bmax)
 end do
 ABI_CALLOC(vcart_ibz, (3, ii:jj, nkibz, nsppol))
 ABI_MALLOC(cg_work, (2, mpw * nspinor))

 cnt = 0
 do spin=1,nsppol
   associate (fs => fstab(spin))
   do ik_ibz=1,ebands%nkpt
     kk = ebands%kptns(:, ik_ibz)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     ! NB: The two checks below are global --> all procs will cycle.
     if (all(bks_mask(:, ik_ibz, spin) .eqv. .False.)) cycle
     if (npw_k == 1) cycle
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do band_k=fs%bmin,fs%bmax
       if (.not. bks_mask(band_k, ik_ibz, spin)) cycle
       call wfd%copy_cg(band_k, ik_ibz, spin, cg_work)
       eig0nk = ebands%eig(band_k, ik_ibz, spin)
       vk = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cg_work, cwaveprj0)
       vcart_ibz(:, band_k, ik_ibz, spin) = vk

       ! TODO: Use ebands_get_edos_matrix_elements
       ! reald(dp) :: vv_fs(3,3,nsppol)
       !vv_fs = zero
       !sigma = fs%eph_fsmear
       !if (fs%eph_fsmear < zero) then
       !  sigma = max(maxval([(abs(dot_product(fs%vk(:, ib2), fs%kmesh_cartvec(:,ii))), ii=1,3)]), fs%min_smear)
       !end if
       !fs_wtk = gaussian(emesh - ebands%eig(band_k, ik_ibz, spin), sigma) / (one * fs%nktot)
       !do ii=1,3
       !  do jj=1,3
       !    vv_fs(ii,jj,spin) = vv_fs(ii,jj,spin) + vk(ii) * vk(jj) * fs_wtk
       !  end do
       !end do
       !call xmpi_sum(vv_fs, comm, ierr)

     end do
   end do
   end associate
 end do ! spin

 call xmpi_sum(vcart_ibz, comm, ierr)
 ABI_FREE(cg_work)
 ABI_FREE(vcart_ibz)
 call cwtime_report(" Velocities", cpu, wall, gflops)
#endif
 call ddkop%free()

 call ngfft_seq(work_ngfft, gmax)
 !write(std_out,*)"work_ngfft(1:3): ",work_ngfft(1:3)
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! FFT meshes from input file, not necessarily equal to the ones found in the external files.
 ! NB: ur arrays are always allocated with nfft and not with product(ngfft(4:6)).
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 !n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 call print_ngfft([std_out], ngfft, header="FFT mesh")

 ! Set the FFT mesh
 call wfd%change_ngfft(cryst, psps, ngfft)
 call wfd%print(units, header="Wavefunctions for GWPT calculation.")

 ABI_MALLOC(gbound_k, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_kq, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_c, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_x, (2*mgfft+8, 2))

 ! ====================================
 ! This is the g-sphere for W_{gg'}(pp)
 ! ====================================
 ! Init g-sphere for the exchange part from ecutsigx.
 call sphereboundary(gbound_c, istwfk1, kg_c, mgfft, npw_c)
 call sphereboundary(gbound_x, istwfk1, kg_x, mgfft, npw_x)

 ABI_MALLOC(vc_sqrt_gx, (npw_x))

 ! Allocate g-vectors centered on k, k+q, k-p, and k+q-p.
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kq, (3, mpw))

 spin = 1
 associate (fs => fstab(spin))
 ABI_MALLOC_OR_DIE(ur_nk,  (nfft*nspinor, fs%bmax-fs%bmin+1), ierr)
 ABI_MALLOC_OR_DIE(ur_mkq, (nfft*nspinor, fs%bmax-fs%bmin+1), ierr)
 !ABI_MALLOC_OR_DIE(my_gbuf, (gqk%cplex, nb, nb, natom3, gqk%my_nk, qbuf_size), ierr)
 ABI_MALLOC(cwork_ur, (nfft*nspinor))

 ! Read symmetrized em1 from file and build ppmodel parameters.
 ! TODO: MPI-shared memory + compute my set of pp-vectors on ppm%new_setup
 id_required = 4; ikxc = 0; approx_type = 0; option_test = 0; nkxc = 0
 ABI_MALLOC(kxcg, (nfftf, nkxc))
 scr_iomode = iomode_from_fname(screen_filepath)

 call epsm1%from_file(screen_filepath, qmesh%nibz, npw_c, comm)
 call epsm1%mkdump(vcp, npw_c, kg_c, nkxc, kxcg, id_required, approx_type, &
                   ikxc, option_test, dtfil%fnameabo_scr, scr_iomode, nfftf, ngfft, comm)

 call epsm1%malloc_epsm1_qbz(npw_c, epsm1%nomega)
 ABI_FREE(kxcg)

 ! This parameter defines the size of the q-buffer used to store the g(k, q) e-ph matrix elements
 ! for all the k-point treated by this MPI rank.
 ! Increasing the buffer size increases the memory requirements
 ! but it leads to better performance as the number of IO operations is decreased.
 ! TODO: Should compute it on the basis of my_nkpt and my_nqpt
 qbuf_size = 1
 !qbuf_size = 16
 call wrtout(std_out, sjoin(" Begin computation of Wqk matrix elements with qbuf_size:", itoa(qbuf_size)), pre_newlines=1)
 call pstat_proc%print(_PSTAT_ARGS_)

 ! Open WQK.nc file and go to data mode.
 prtwqk = 1
 if (prtwqk /= 0) then
   NCF_CHECK(nctk_open_create(root_ncid, strcat(dtfil%filnam_ds(4), "_WQK.nc") , comm))
   NCF_CHECK(nctk_set_datamode(root_ncid))
 end if

 !ABI_MALLOC(iq_buf, (2, qbuf_size))

 ! ==================
 ! Loop over k-points
 ! ==================
 mu = zero
 do ikfs_bz=1,fs%nkfs
   my_ik = ikfs_bz
   my_nk = fs%nkfs
   ! NB: All procs in gqk%pert_comm and gqk%bsum_com and gqk%pp_sum_comm enter this section.
   !iqbuf_cnt = 1 + mod(my_iq - 1, qbuf_size)
   !iq_buf(:, iqbuf_cnt) = [my_iq, iq_bz]

   ! Set entry to zero. Important as there are cycle instructions inside these loops
   ! and we don't want random numbers written to disk.
   !my_gbuf(:,:,:,:, ikfs_bz, iqbuf_cnt) = zero
   print_time_kk = my_rank == 0 .and. (ikfs_bz <= LOG_MODK .or. mod(ikfs_bz, LOG_MODK) == 0)
   if (print_time_kk) call cwtime(cpu_kk, wall_kk, gflops_kk, "start")

   ! Symmetry indices for kk.
   kk = fs%kpts(:, ikfs_bz)
   ! The k-point and the symmetries relating the BZ k-point to the IBZ.
   ik_ibz = fs%indkk_fs(1, ikfs_bz) ; isym_k = fs%indkk_fs(2, ikfs_bz)
   trev_k = fs%indkk_fs(6, ikfs_bz); g0_k = fs%indkk_fs(3:5,ikfs_bz)
   isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
   mapl_k = fs%indkk_fs(:, ikfs_bz)
   !print *, "ikfs_bz", ikfs_bz, " of my_nk:", gqk%my_nk

   kk_ibz = ebands%kptns(:,ik_ibz)
   istwf_k_ibz = wfd%istwfk(ik_ibz); npw_k_ibz = wfd%npwarr(ik_ibz)
   !print *, "ik_ibz:", ik_ibz, "kk:", kk, "kk_ibz:", kk_ibz

   ! Get npw_k, kg_k for kk
   call wfd%get_gvec_gbound(cryst%gmet, ecut, kk, ik_ibz, isirr_k, dtset%nloalg, & ! in
                            istwf_k, npw_k, kg_k, nkpg_k, kpg_k, gbound_k)         ! out

   ! Precompute ur_nk and for all n_k bands.
   bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz); bstop_k = bstart_k + nband_k - 1
   ABI_MALLOC(ug_k, (2, npw_k*nspinor, nband_k))
   call wfd%rotate_cg(bstart_k, nband_k, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                      cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_nk)

   ! ==================
   ! Loop over q-points
   ! ==================
   do iq_bz=1,qmesh%nbz
     qq_bz = qmesh%bz(:, iq_bz)
     qq_is_gamma = sum(qq_bz**2) < tol14

     ! Find k + q in the extended zone and extract symmetry info.
     ! Be careful here because there are two umklapp vectors to be considered as:
     !
     !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
     !
     kq = kk + qq_bz
     ikq_fs = fs%findkg0(kq, g0_kq)
     if (ikq_fs == -1) cycle

     mapl_kq = fs%indkk_fs(:, ikq_fs)
     ikq_ibz = mapl_kq(1); isym_kq = mapl_kq(2); trev_kq = mapl_kq(6); g0_kq = mapl_kq(3:5)
     isirr_kq = (isym_kq == 1 .and. trev_kq == 0 .and. all(g0_kq == 0))
     kq_ibz = ebands%kptns(:, ikq_ibz)
     istwf_kq_ibz = wfd%istwfk(ikq_ibz); npw_kq_ibz = wfd%npwarr(ikq_ibz)

     ! Get npw_kq, kg_kq for k+q.
     call wfd%get_gvec_gbound(cryst%gmet, ecut, kq, ikq_ibz, isirr_kq, dtset%nloalg, &  ! in
                              istwf_kq, npw_kq, kg_kq, nkpg_kq, kpg_kq, gbound_kq)      ! out
     ABI_FREE(kpg_kq)

     ! Precompute ur_mkq for all m_kq bands.
     bstart_kq = fs%bstart_cnt_ibz(1, ikq_ibz); nband_kq = fs%bstart_cnt_ibz(2, ikq_ibz); bstop_kq = bstart_kq + nband_kq - 1
     ABI_MALLOC(ug_kq, (2, npw_kq*nspinor))
     ABI_MALLOC(rhotwg_x, (npw_x*nspinor, nband_k, nband_kq))

     do m_kq=bstart_kq, bstop_kq
       ib_kq = m_kq - bstart_kq + 1
       call wfd%rotate_cg(m_kq, ndat1, spin, kq_ibz, npw_kq, kg_kq, istwf_kq, &
                          cryst, mapl_kq, gbound_kq, work_ngfft, work, ug_kq, urs_kbz=ur_mkq(:,ib_kq))

       do n_k=bstart_k, bstop_k
         ib_k = n_k - bstart_k + 1
         cwork_ur = conjg(ur_nk(:,ib_k)) * ur_mkq(:,ib_kq)
         call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_x(:, ib_k, ib_kq))
         call sigtk_multiply_by_vc_sqrt("N", npw_x, nspinor, 1, vc_sqrt_gx, rhotwg_x(:, ib_k, ib_kq))
       end do
     end do

     ABI_FREE(ug_kq)

     ! Find the corresponding irred pp-point in the pp_mesh.
     call qmesh%get_bz_item(iq_bz, qq_bz, iq_ibz, isym_qq, itim_qq)

     ! Get Fourier components of the Coulomb interaction in the BZ
     ! In 3D systems, neglecting umklapp: vc(Sq,sG) = vc(q,G) = 4pi/|q+G|**2
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.

     ! NOTE: vc_sqrt_gx is dimensioned with npw_x --> use rottb table from gsph_x.
     do ig=1,npw_x
       vc_sqrt_gx(gsph_x%rottb(ig, itim_qq, isym_qq)) = vcp%vc_sqrt(ig, iq_ibz)
     end do

     remove_exchange = .True.
     call epsm1%rotate_iqbz(iq_bz, hscr%nomega, npw_c, gsph_c, qmesh, remove_exchange)

     ! Compute \sum_{gg'} (M_g^{kn,k'm})^* W_{gg'}(q) M_g'^{kn,k'm} with k' = k + q.

     ! TODO: nspinor
     ABI_MALLOC(rhotwg_c, (npw_c*nspinor, nband_k, nband_kq))
     rhotwg_c(:,:,:) = rhotwg_x(1:npw_c,:,:)

     ABI_MALLOC(w_rhotwg_c, (npw_c*nspinor, nband_k, nband_kq))
     ncols = nband_k * nband_kq

     !call xgemm("N", "N", npw_c, ncols, npw_c, cone_gw, epsm1%epsm1_qbz(:,:,1), npw_c, rhotwg_c(:,:,1), npw_c*nspinor, &
     !           czero_gw, w_rhotwg_c(:,:,1), npw_c)

     do ib_k=1, nband_k
       do ib_kq=1, nband_kq
         w_rhotwg_c(:,ib_k,ib_kq) = matmul(epsm1%epsm1_qbz(:,:,1), rhotwg_c(:,ib_k,ib_kq))
       end do
     end do
     !print *, "epsm1_qibz:", maxval(abs(epsm1%epsm1_qbz(:,:,1)))
     !if (.not. qq_is_gamma) then
     !if (qq_is_gamma) then
     !print *,  " rhotwg_c",  maxval(abs(rhotwg_c)), maxloc(abs(rhotwg_c))
     !end if
     !print *, "w_rhotwg_c:", maxval(abs(w_rhotwg_c))

     do n_k=bstart_k, bstop_k
       ib_k = n_k - bstart_k + 1
       e_nk = ebands%eig(n_k, ik_ibz, spin)
       !gaussian(e_nk - ef, smear)
       do m_kq=bstart_kq, bstop_kq
         ib_kq = m_kq - bstart_kq + 1
         e_mkq = ebands%eig(m_kq, ikq_ibz, spin)
         !gaussian(e_mkq - ef, smear)
         !rhotwg_x(npw_x*nspinor, ib_k, ib_kq)); w_rhotwg_c(npw_c*nspinor, ib_k, ib_kq)
         ctmp_gwpc = xdotc(npw_c*nspinor, rhotwg_c(:,ib_k, ib_kq), 1, w_rhotwg_c(:,ib_k, ib_kq), 1)
         if (remove_exchange) then
           ctmp_gwpc = ctmp_gwpc + xdotc(npw_x*nspinor, rhotwg_x(:,ib_k, ib_kq), 1, rhotwg_x(:,ib_k, ib_kq), 1)
         end if
         !print *, "ctmp_gwpc:", ctmp_gwpc
         !mu = mu + ctmp_gwpc *
         !w_ee =
       end do
     end do

     ABI_FREE(rhotwg_x)
     ABI_FREE(rhotwg_c)
     ABI_FREE(w_rhotwg_c)

     !print_time_qq = my_rank == 0 .and. (iq_bz <= LOG_MODQ .or. mod(iq_bz, LOG_MODQ) == 0)
     !if (print_time_qq) then
     !  call cwtime(cpu_qq, wall_qq, gflops_qq, "start")
     !  call inds2str(0, sjoin(" Computing g^Sigma(k, q) for qq_bz:", ktoa(qq_bz)), iq_bz, gqk%my_nq, gqk%glob_nq, msg)
     !  call wrtout(std_out, sjoin(msg, ", for spin:", itoa(spin)), pre_newlines=1)
     !end if

     ! Dump buffer
     !if (iqbuf_cnt == qbuf_size) call dump_my_gbuf()

     !if (print_time_qq) then
     !  call inds2str(2, "My q-point", my_iq, my_nq, qmesh%nbz, msg)
     !  call cwtime_report(msg, cpu_qq, wall_qq, gflops_qq); if (iq_bz == LOG_MODQ) call wrtout(std_out, "...", do_flush=.True.)
     !end if
   end do ! iq_ibz

    if (print_time_kk) then
      call inds2str(1, "My k-point", my_ik, my_nk, fs%nkfs, msg)
      call cwtime_report(msg, cpu_kk, wall_kk, gflops_kk); if (ikfs_bz == LOG_MODK) call wrtout(std_out, "...", do_flush=.True.)
    end if

    ABI_FREE(ug_k)
    ABI_FREE(kpg_k)
 end do ! ikfs_bz
 end associate

 ! Dump the remainder.
 !if (iqbuf_cnt /= 0) call dump_my_gbuf()

 ABI_FREE(ur_nk)
 ABI_FREE(ur_mkq)
 !ABI_FREE(iq_buf)
 !ABI_FREE(my_gbuf)
 !end do ! my_is

 !call xmpi_sum_master(w_ee, master, comm, ierr)

 call cwtime_report(" wqk calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Set vqk_completed to 1 so that we can easily check if restarted is needed.
 if (prtwqk /= 0) then
   !if (my_rank == master) then
   !  NCF_CHECK(nf90_put_var(root_ncid, root_vid("vqk_completed"), 1))
   !end if
   NCF_CHECK(nf90_close(root_ncid))
   call xmpi_barrier(comm)
 end if

 ! Output some of the results to ab_out for testing purposes
 !call gstore%print_for_abitests(dtset)

 ! Free memory
 ABI_FREE(kg_k)
 ABI_FREE(kg_kq)
 ABI_FREE(cwork_ur)
 ABI_FREE(work)
 ABI_FREE(gbound_k)
 ABI_FREE(gbound_kq)
 ABI_FREE(gbound_c)
 ABI_FREE(gbound_x)
 !ABI_FREE(done_qbz_spin)
 ABI_FREE(vc_sqrt_gx)
 ABI_FREE(w_ee)
 ABI_FREE(e_mesh)

 call wfd%free(); call vcp%free(); call qmesh%free(); call gsph_x%free(); call gsph_c%free(); call hscr%free(); call edos%free()
 call epsm1%free()

 do spin=1,ebands%nsppol
   call fstab(spin)%free()
 end do
 ABI_FREE(fstab)

 call xmpi_barrier(comm) ! This to make sure that the parallel output of GSTORE is completed
 call cwtime_report(" wqk_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

contains

subroutine inds2str(ilevel, prefix, my_ik, my_nk, nk_tot, out_str)
 character(len=*),intent(in) :: prefix
 integer,intent(in) :: ilevel, my_ik, my_nk, nk_tot
 character(len=*),intent(out) :: out_str

 out_str = sjoin(prefix, itoa(my_ik), "/", itoa(my_nk), "[", itoa(nk_tot), "]")
 out_str = repeat(' ', 4 * ilevel) // trim(out_str)
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

 !integer :: ii, iq_bz, iq_glob, my_iq

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
 !  do ii=1,iqbuf_cnt
 !    iq_bz = iq_buf(2, ii)
 !    NCF_CHECK(nf90_put_var(root_ncid, root_vid("gstore_done_qbz_spin"), 1, start=[iq_bz, spin]))
 !  end do
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

!integer function spin_vid(var_name)
!  character(len=*),intent(in) :: var_name
!  spin_vid = nctk_idname(spin_ncid, var_name)
!end function spin_vid

end subroutine wqk_run
!!***

end module m_wqk
!!***

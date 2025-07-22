!!****m* ABINIT/m_wkk
!! NAME
!!  m_wkk
!!
!! FUNCTION
!!  Computation of the matrix elements of the screened interaction W_kk'.
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

module m_wkk

 !use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_mpinfo
 use m_errors
 use m_fft
 use m_wfd
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
 use m_numeric_tools,  only : arth, linspace ! print_arr
 use m_io_tools,       only : iomode_from_fname
 use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, print_ngfft
 use m_crystal,        only : crystal_t
 use m_bz_mesh,        only : kmesh_t, findqg0
 use m_htetra,         only : htetra_t
 use m_gsphere,        only : gsphere_t
 use m_pawtab,         only : pawtab_type
 use m_io_screening,   only : hscr_t, get_hscr_qmesh_gsph
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

 public :: wkk_run  ! Main entry point to compute Wkk' matrix elements.

!----------------------------------------------------------------------

contains  !=====================================================
!!***

!----------------------------------------------------------------------

!!****f* m_wkk/wkk_run
!! NAME
!!  wkk_run
!!
!! FUNCTION
!!  Computation of the matrix elements of the screened interaction W_kk'
!!  between two Cooper pairs.
!!
!! INPUTS
!! wfk0_path=String with the path to the GS unperturbed WFK file.
!! dtfil<datafiles_type>=Variables related to files.
!! ngfft(18),ngfftf(18)=Coarse and Fine FFT meshes.
!! dtset<dataset_type>=All input variables for this dataset.
!! cryst: Crystalline structure
!! ebands<ebands_t>=The GS KS band structure (energies, occupancies, k-weights...)
!! wfk_hdr=Header of the WFK file.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!! comm=MPI communicator.
!!
!! OUTPUT
!!  Results are written to file in netcdf format.
!!
!! SOURCE

subroutine wkk_run(wfk0_path, dtfil, ngfft, ngfftf, dtset, cryst, ebands, wfk_hdr, &
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
 integer :: cplex,nspinor,nprocs, ii, ib, my_is, spin, npw_x, npw_c ! max_npw_xc, min_npw_xc,
 integer :: bstart_k, bstop_k, nband_k, ncols, mpw, ierr, ncerr
 integer :: bstart_kp, bstop_kp, nband_kp, ik_bz, bmin, bmax, max_nb
 integer :: ikp_ibz, isym_kp, trev_kp, npw_kp, istwf_kp, npw_kp_ibz, istwf_kp_ibz
 integer :: ik_ibz, isym_k, trev_k, npw_k, istwf_k, npw_k_ibz, istwf_k_ibz
 integer :: m_k, im_k, n_kp, in_kp, root_ncid, wkk_mode, ne
 integer :: nfft,nfftf,mgfft,mgfftf,nkpg_kp,nkpg_k,cnt, edos_intmeth, nqlwl, scr_iomode
 integer :: ikp_bz, my_ikp, my_nkp !, my_iq,
 real(dp) :: cpu_all, wall_all, gflops_all, cpu, wall, gflops, cpu_kp, wall_kp, gflops_kp
 real(dp) :: edos_step, edos_broad, e_mk, e_nkp, e_min, e_max, e_step, smear_mk, smear_nkp, faq
 logical :: isirr_k, isirr_kp, qq_is_gamma, remove_exchange, print_time_kp
 type(wfd_t) :: wfd
 type(kmesh_t) :: qmesh, kmesh
 type(gsphere_t),target :: gsph_x, gsph_c
 type(hscr_t),target :: hscr
 type(vcoul_t) :: vcp
 type(edos_t) :: edos
 type(epsm1_t) :: epsm1
 type(ddkop_t) :: ddkop
 type(htetra_t) :: tetra
 character(len=fnlen) :: screen_filepath
 character(len=5000) :: msg, qkp_string
!arrays
 integer :: g0_k(3), g0_kp(3), g0_qq(3), units(2), work_ngfft(18), gmax(3), mapl_k(6), mapl_kp(6), mg0(3)
 integer,allocatable :: kg_kp(:,:), kg_k(:,:), gbound_kp(:,:), gbound_k(:,:), gbound_c(:,:), gbound_x(:,:)
 integer,allocatable :: nband(:,:), wfd_istwfk(:)
 integer, contiguous, pointer :: kg_c(:,:), kg_x(:,:)
 real(dp) :: kk_ibz(3), kk_bz(3), kp_ibz(3), kp_bz(3), qq_bz(3), kk_diff(3) ! qq_ibz(3)
 complex(gwpc) :: ctmp_gwpc, mu_c
 character(len=fnlen) :: path
!arrays
 real(dp) :: n0(ebands%nsppol)
 real(dp),allocatable :: qlwl(:,:), kpg_kp(:,:), kpg_k(:,:), ug_kp(:,:,:), ug_k(:,:), cg_work(:,:)
 real(dp),allocatable :: work(:,:,:,:), e_mesh(:), e_args(:), vcart_ibz(:,:,:,:), wgt_mk(:), wgt_nkp(:,:)
 complex(gwpc),allocatable :: cwork_ur(:), rhotwg_mn_x(:,:,:), rhotwg_mn_c(:,:,:), w_rhotwg_mn_c(:,:,:)
 complex(gwpc),allocatable :: vc_sqrt_gx(:), ur_nkp(:,:), ur_mk(:,:), kxcg(:,:)
 complex(dp),allocatable :: w_ee(:,:) ! w_kkp(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:), keep_ur(:,:,:)
 type(fstab_t),target,allocatable :: fstab(:)
!************************************************************************

 ABI_CHECK(psps%usepaw == 0, "PAW not implemented")
 ABI_CHECK(dtset%nsppol == 1, "nsppol 2 implemented!")
 ABI_CHECK(dtset%nspinor == 1, "nspinor 2 implemented!")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm); units = [std_out, ab_out]
 call cwtime(cpu_all, wall_all, gflops_all, "start")

 nsppol = ebands%nsppol; nspinor = ebands%nspinor

 ! Compute electron DOS.
 edos_intmeth = 2; if (dtset%prtdos /= 0) edos_intmeth = dtset%prtdos
 edos_step = dtset%dosdeltae; edos_broad = dtset%tsmear
 edos_step = 0.01 * eV_Ha; edos_broad = 0.3 * eV_Ha
 edos = ebands%get_edos(cryst, edos_intmeth, edos_step, edos_broad, comm)
 ! Get DOS per spin channel
 n0(:) = edos%gef(1:edos%nsppol)
 if (my_rank == master) then
   call edos%print(units)
   path = strcat(dtfil%filnam_ds(4), "_EDOS")
   call edos%write(path)
 end if

 ! Find Fermi surface k-points.
 ! TODO: support kptopt, change setup of k-points if tetra: fist tetra weights then k-points on the Fermi surface!
 ABI_MALLOC(fstab, (nsppol))
 call fstab_init(fstab, ebands, cryst, dtset, tetra, comm)
 call tetra%free()

 bmin = huge(1); bmax = -1
 do spin=1,nsppol
   bmin = min(bmin, fstab(spin)%bmin)
   bmax = max(bmax, fstab(spin)%bmax)
 end do
 max_nb = bmax - bmin + 1

 ! Build linear mesh for W(e,e').
 spin = 1
 associate (fs => fstab(spin))
 e_min = minval(ebands%eig(fs%bmin,:,spin)) - 0.1_dp * eV_Ha
 e_max = maxval(ebands%eig(fs%bmax,:,spin)) + 0.1_dp * eV_Ha
 e_step = 0.002 * eV_Ha
 ne = 1 + (e_max - e_min) / e_step
 ABI_CALLOC(w_ee, (ne, ne))
 ABI_MALLOC(e_mesh, (ne))
 ABI_MALLOC(e_args, (ne))
 ABI_MALLOC(wgt_mk, (ne))
 ABI_MALLOC(wgt_nkp, (ne, max_nb))
 e_mesh = linspace(e_min, e_max, ne)
 end associate

 ! ================
 ! HANDLE SCREENING
 ! ================
 screen_filepath = dtfil%fnameabi_scr
 ABI_CHECK(dtfil%fnameabi_scr /= ABI_NOFILE, "SCR file must be specified in input!")

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

 ! Initialize Coulomb term on the IBZ of the qmesh. Use largest G-sphere.
 call kmesh%init(cryst, wfk_hdr%nkpt, wfk_hdr%kptns, dtset%kptopt)
 faq = one / (cryst%ucvol*kmesh%nbz)

 ! TODO:
 ! Here we sort the pp_mesh by stars so that we can split the pp wavevectors in blocks and therefore
 ! reduce the number of wavevectors in the IBZ that must be stored in memory.
 !call qmesh%pack_by_stars()

 npw_x = gsph_x%ng; npw_c = gsph_c%ng !; min_npw_xc = min(npw_x, npw_c); max_npw_xc = max(npw_x, npw_c)
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
               dtset%nspden, nspinor, dtset%ecut, dtset%ecutsm, dtset%dilatmx, wfd_istwfk, ebands%kptns, ngfft,&
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

 call ephtk_get_mpw_gmax(nkibz, wfk_hdr%kptns, dtset%ecut, cryst%gmet, mpw, gmax, comm)
 ! FIXME
 !gmax = 2 * gmax

 ddkop = ddkop_new(dtset, cryst, pawtab, psps, wfd%mpi_enreg, mpw, wfd%ngfft)

#if 0
 call cwtime(cpu, wall, gflops, "start", msg=" Computing v_nk matrix elements for all states on the FS...")

 ABI_CALLOC(vcart_ibz, (3, bmin:bmax, nkibz, nsppol))
 ABI_MALLOC(cg_work, (2, mpw * nspinor))

 cnt = 0
 do spin=1,nsppol
   associate (fs => fstab(spin))
   do ik_ibz=1,ebands%nkpt
     kk_ibz = ebands%kptns(:, ik_ibz)
     npw_k = wfd%npwarr(ik_ibz); istwf_k = wfd%istwfk(ik_ibz)
     ! NB: The two checks below are global --> all procs will cycle.
     if (all(bks_mask(:, ik_ibz, spin) .eqv. .False.)) cycle
     if (npw_k == 1) cycle
     cnt = cnt + 1; if (mod(cnt, nproc) /= my_rank) cycle ! MPI parallelism.

     call ddkop%setup_spin_kpoint(dtset, cryst, psps, spin, kk_ibz, istwf_k, npw_k, wfd%kdata(ik_ibz)%kg_k)

     do band_k=fs%bmin,fs%bmax
       if (.not. bks_mask(band_k, ik_ibz, spin)) cycle
       call wfd%copy_cg(band_k, ik_ibz, spin, cg_work)
       eig0nk = ebands%eig(band_k, ik_ibz, spin)
       vk = ddkop%get_vdiag(eig0nk, istwf_k, npw_k, wfd%nspinor, cg_work, cwaveprj0)
       vcart_ibz(:, band_k, ik_ibz, spin) = vk
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
 ABI_MALLOC(work, (2, work_ngfft(4), work_ngfft(5), work_ngfft(6)))

 ! FFT meshes from input file, not necessarily equal to the ones found in the external files.
 ! NB: ur arrays are always allocated with nfft and not with product(ngfft(4:6)).
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))

 call print_ngfft([std_out], ngfft, header="FFT mesh")

 ! Set the FFT mesh
 call wfd%change_ngfft(cryst, psps, ngfft)
 call wfd%print(units, header="Wavefunctions for GWPT calculation.")

 ABI_MALLOC(gbound_kp, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_k, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_c, (2*mgfft+8, 2))
 ABI_MALLOC(gbound_x, (2*mgfft+8, 2))

 ! ====================================
 ! This is the g-sphere for W_{gg'}(qq)
 ! ====================================
 ! Init g-sphere for the exchange part from ecutsigx.
 call sphereboundary(gbound_c, istwfk1, kg_c, mgfft, npw_c)
 call sphereboundary(gbound_x, istwfk1, kg_x, mgfft, npw_x)

 ABI_MALLOC(vc_sqrt_gx, (npw_x))

 ! Allocate g-vectors centered on k, k'
 ABI_MALLOC(kg_k, (3, mpw))
 ABI_MALLOC(kg_kp, (3, mpw))

 spin = 1
 associate (fs => fstab(spin))
 ABI_MALLOC_OR_DIE(ur_nkp,  (nfft*nspinor, fs%bmax-fs%bmin+1), ierr)
 ABI_MALLOC_OR_DIE(ur_mk, (nfft*nspinor, fs%bmax-fs%bmin+1), ierr)
 ABI_MALLOC(cwork_ur, (nfft*nspinor))

 ! Read symmetrized em1 from file
 id_required = 4; ikxc = 0; approx_type = 0; option_test = 0; nkxc = 0
 ABI_MALLOC(kxcg, (nfftf, nkxc))
 scr_iomode = iomode_from_fname(screen_filepath)

 call epsm1%from_file(screen_filepath, qmesh%nibz, npw_c, comm)
 call epsm1%mkdump(vcp, npw_c, kg_c, nkxc, kxcg, id_required, approx_type, &
                   ikxc, option_test, dtfil%fnameabo_scr, scr_iomode, nfftf, ngfft, comm)

 call epsm1%malloc_epsm1_qbz(npw_c, epsm1%nomega)
 ABI_FREE(kxcg)

 call pstat_proc%print(_PSTAT_ARGS_)

 wkk_mode = 1
 ! Master creates the netcdf file used to store the results of the calculation.
 if (my_rank == master) then
   NCF_CHECK(nctk_open_create(root_ncid, strcat(dtfil%filnam_ds(4), "_WKK.nc") , xmpi_comm_self))
   NCF_CHECK(cryst%ncwrite(root_ncid))
   NCF_CHECK(ebands%ncwrite(root_ncid))
   NCF_CHECK(edos%ncwrite(root_ncid))

   !ncerr = nctk_def_dims(ncid, [ &
   !  nctkdim_t("smat_bsize1", smat_bsize1), nctkdim_t("smat_bsize2", smat_bsize2) &
   !  ], defmode=.True.)
   !NCF_CHECK(ncerr)
   !ncerr = nctk_def_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "gwr_completed", "scf_iteration" &
   !])
   !NCF_CHECK(ncerr)
   !ncerr = nctk_def_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "wr_step", "gwr_boxcutmin", "cosft_duality_error", "regterm" &
   !])
   !NCF_CHECK(ncerr)

   !ncerr = nctk_def_arrays(ncid, [ &
   !  nctkarr_t("gwr_task", "char", "character_string_length"), &
   !])
   !NCF_CHECK(ncerr)
   !NCF_CHECK(nctk_set_datamode(ncid))

   !ncerr = nctk_write_iscalars(ncid, [character(len=nctk_slen) :: &
   !  "gwr_completed", "sig_diago", "b1gw", "b2gw", "symsigma", "symchi", "scf_iteration"], &
   !  [0, merge(1, 0, gwr%sig_diago), gwr%b1gw, gwr%b2gw, gwr%dtset%symsigma, dtset%symchi, gwr%scf_iteration])
   !NCF_CHECK(ncerr)

   !ncerr = nctk_write_dpscalars(ncid, [character(len=nctk_slen) :: &
   !  "wr_step", "ecuteps", "ecut", "ecutwfn", "ecutsigx", "gwr_boxcutmin", &
   !  "cosft_duality_error", "regterm"], &
   !  [gwr%wr_step, dtset%ecuteps, dtset%ecut, dtset%ecutwfn, dtset%ecutsigx, dtset%gwr_boxcutmin, &
   !   gwr%ft_max_error(1), gwr%ft_max_error(2), gwr%ft_max_error(3), gwr%cosft_duality_error, regterm &
   !  ])
   !NCF_CHECK(ncerr)
   !NCF_CHECK(nf90_put_var(ncid, vid("gwr_task"), trim(dtset%gwr_task)))

   NCF_CHECK(nf90_close(root_ncid))
 end if
 call xmpi_barrier(comm)

 ! Open WKK.nc file and go to data mode.
 NCF_CHECK(nctk_open_modify(root_ncid, strcat(dtfil%filnam_ds(4), "_WKK.nc") , comm))
 NCF_CHECK(nctk_set_datamode(root_ncid))

 mg0 = [1, 1, 1]
 remove_exchange = .True.
 mu_c = zero

 ! Loop over k'-points in the energy window.
 do ikp_bz=1,fs%nkfs
   my_ikp = ikp_bz; my_nkp = fs%nkfs

   print_time_kp = my_rank == 0 .and. (ikp_bz <= LOG_MODK .or. mod(ikp_bz, LOG_MODK) == 0)
   if (print_time_kp) call cwtime(cpu_kp, wall_kp, gflops_kp, "start")

   ! The k-point and the symmetries relating the BZ k-point to the IBZ.
   kp_bz = fs%kpts(:, ikp_bz)
   ikp_ibz = fs%indkk_fs(1, ikp_bz) ; isym_kp = fs%indkk_fs(2, ikp_bz)
   trev_kp = fs%indkk_fs(6, ikp_bz); g0_kp = fs%indkk_fs(3:5,ikp_bz)
   isirr_kp = (isym_kp == 1 .and. trev_kp == 0 .and. all(g0_kp == 0))
   mapl_kp = fs%indkk_fs(:, ikp_bz)

   kp_ibz = ebands%kptns(:,ikp_ibz)
   istwf_kp_ibz = wfd%istwfk(ikp_ibz); npw_kp_ibz = wfd%npwarr(ikp_ibz)

   ! Get npw_kp, kg_kp for this kp_bz.
   call wfd%get_gvec_gbound(cryst%gmet, dtset%ecut, kp_bz, ikp_ibz, isirr_kp, dtset%nloalg, & ! in
                            istwf_kp, npw_kp, kg_kp, nkpg_kp, kpg_kp, gbound_kp)        ! out

   ! Rotate from IBZ to BZ and compute ur_nkp for all n_kp bands.
   bstart_kp = fs%bstart_cnt_ibz(1, ikp_ibz); nband_kp = fs%bstart_cnt_ibz(2, ikp_ibz); bstop_kp = bstart_kp + nband_kp - 1
   ABI_CALLOC(ug_kp, (2, npw_kp*nspinor, nband_kp))
   call wfd%rotate_cg(bstart_kp, nband_kp, spin, kp_ibz, npw_kp, kg_kp, istwf_kp, &
                      cryst, mapl_kp, gbound_kp, work_ngfft, work, ug_kp, urs_kbz=ur_nkp)
   do ii=1, nband_kp
    !print *, "isirr_kp:", isirr_kp
    !print *, "g:", sum(ug_kp(1,:,ii)**2 + ug_kp(2,:,ii)**2)
    !print *, "r:", sum(abs(ur_nkp(:,ii)) ** 2) / nfft ! * cryst%ucvol /
   end do

   ! Loop over k-points in the energy window.
   do ik_bz=1,fs%nkfs

     mapl_k = fs%indkk_fs(:, ik_bz)
     ik_ibz = mapl_k(1); isym_k = mapl_k(2); trev_k = mapl_k(6); g0_k = mapl_k(3:5)
     isirr_k = (isym_k == 1 .and. trev_k == 0 .and. all(g0_k == 0))
     kk_ibz = ebands%kptns(:, ik_ibz)
     istwf_k_ibz = wfd%istwfk(ik_ibz); npw_k_ibz = wfd%npwarr(ik_ibz)

     kk_bz = fs%kpts(:, ik_bz)
     kk_ibz = ebands%kptns(:,ik_ibz)

     ! Get npw_k, kg_k for this k.
     call wfd%get_gvec_gbound(cryst%gmet, dtset%ecut, kk_bz, ik_ibz, isirr_k, dtset%nloalg, &  ! in
                              istwf_k, npw_k, kg_k, nkpg_k, kpg_k, gbound_k)             ! out
     ABI_FREE(kpg_k)

     bstart_k = fs%bstart_cnt_ibz(1, ik_ibz); nband_k = fs%bstart_cnt_ibz(2, ik_ibz); bstop_k = bstart_k + nband_k - 1
     ABI_MALLOC(ug_k, (2, npw_k*nspinor))
     ABI_MALLOC(rhotwg_mn_x, (npw_x*nspinor, nband_k, nband_kp))

     do m_k=bstart_k, bstop_k
       ! Rotate from IBZ to BZ and compute ur_mk
       im_k = m_k - bstart_k + 1
       call wfd%rotate_cg(m_k, ndat1, spin, kk_ibz, npw_k, kg_k, istwf_k, &
                          cryst, mapl_k, gbound_k, work_ngfft, work, ug_k, urs_kbz=ur_mk(:,im_k))

       do n_kp=bstart_kp, bstop_kp
         in_kp = n_kp - bstart_kp + 1
         cwork_ur = conjg(ur_nkp(:,in_kp)) * ur_mk(:,im_k)
         call fft_ur(npw_x, nfft, nspinor, ndat1, mgfft, ngfft, istwfk1, kg_x, gbound_x, cwork_ur, rhotwg_mn_x(:, im_k, in_kp))
         call sigtk_multiply_by_vc_sqrt("N", npw_x, nspinor, ndat1, vc_sqrt_gx, rhotwg_mn_x(:, im_k, in_kp))
       end do
     end do

     ABI_FREE(ug_k)

     ! Identify q and G0 where q + G0 = k_GW - k_i
     kk_diff = kk_bz - kp_bz

     call findqg0(iq_bz, g0_qq, kk_diff, qmesh%nbz, qmesh%bz, mG0)

     ! Find the corresponding irred qq-point in qmesh.
     call qmesh%get_bz_item(iq_bz, qq_bz, iq_ibz, isym_qq, itim_qq)
     qq_is_gamma = sum(qq_bz**2) < tol14

     ! Find k + q in the extended zone and extract symmetry info.
     ! Be careful here because there are two umklapp vectors to be considered as:
     !
     !   k + q = k_bz + g0_bz = IS(k_ibz) + g0_ibz + g0_bz
     !
     ! Get Fourier components of the Coulomb interaction in the BZ
     ! In 3D systems, neglecting umklapp: vc(Sq,sG) = vc(q,G) = 4pi/|q+G|**2
     ! The same relation holds for 0-D systems, but not in 1-D or 2D systems. It depends on S.
     ! NOTE: vc_sqrt_gx is dimensioned with npw_x --> use rottb table from gsph_x.
     do ig=1,npw_x
       vc_sqrt_gx(gsph_x%rottb(ig, itim_qq, isym_qq)) = vcp%vc_sqrt(ig, iq_ibz)
     end do

     call epsm1%rotate_iqbz(iq_bz, hscr%nomega, npw_c, gsph_c, qmesh, remove_exchange)

     ! Compute \sum_{gg'} (M_g^{kn,k'm})^* W_{gg'}(q) M_g'^{kn,k'm} with k' = k + q.
     ! TODO: nspinor
     ABI_MALLOC(rhotwg_mn_c, (npw_c*nspinor, nband_k, nband_kp))
     rhotwg_mn_c(:,:,:) = rhotwg_mn_x(1:npw_c,:,:)

     ncols = nband_k * nband_kp
     ABI_MALLOC(w_rhotwg_mn_c, (npw_c*nspinor, nband_k, nband_kp))

     call xgemm("N", "N", npw_c, ncols, npw_c, cone_gw, epsm1%epsm1_qbz(:,:,1), npw_c, rhotwg_mn_c(:,:,1), npw_c*nspinor, &
                czero_gw, w_rhotwg_mn_c(:,:,1), npw_c)

     !do in_kp=1, nband_k
     !  do im_k=1, nband_kq
     !    w_rhotwg_mn_c(:,in_kp,im_k) = matmul(epsm1%epsm1_qbz(:,:,1), rhotwg_mn_c(:,im_k,in_kp))
     !  end do
     !end do
     !print *,  " rhotwg_mn_c",  maxval(abs(rhotwg_mn_c)), maxloc(abs(rhotwg_mn_c))
     !print *, "w_rhotwg_mn_c:", maxval(abs(w_rhotwg_mn_c))

     ! Precompute delta(e - e_nkp - ef).
     do n_kp=bstart_kp, bstop_kp
       in_kp = n_kp - bstart_kp + 1
       e_nkp = ebands%eig(n_kp, ik_ibz, spin)
       e_args = e_mesh - e_nkp
       smear_nkp = 0.1_dp * eV_Ha
       wgt_nkp(:, im_k) = gaussian(e_args, smear_nkp)
     end do

     do n_kp=bstart_kp, bstop_kp
       in_kp = n_kp - bstart_kp + 1
       e_nkp = ebands%eig(n_kp, ikp_ibz, spin)
       !fs%tetra_wtk(n_kp - fs%bmin + 1, ikp_ibz)
       do m_k=bstart_k, bstop_k
         im_k = m_k - bstart_k + 1
         e_mk = ebands%eig(m_k, ik_ibz, spin)
         !fs%tetra_wtk(m_k - fs%bmin + 1, ik_ibz)
         smear_mk = 0.1_dp * eV_Ha
         ctmp_gwpc = xdotc(npw_c*nspinor, rhotwg_mn_c(:,im_k,in_kp), 1, w_rhotwg_mn_c(:,im_k,in_kp), 1)
         if (remove_exchange) then
           ctmp_gwpc = ctmp_gwpc + xdotc(npw_x*nspinor, rhotwg_mn_x(:,im_k,in_kp), 1, rhotwg_mn_x(:,im_k,in_kp), 1)
         end if
         ctmp_gwpc = faq * ctmp_gwpc
         !print *, "ctmp_gwpc:", ctmp_gwpc
         !mu_c = mu_c + ctmp_gwpc * gaussian(e_mk - ebands%fermie, smear_mk) * gaussian(e_mp, smear_mk)
         ! Computes: A := A + alpha * x * y^T
         !w_ee =
         !call dger(ne, ne, alpha, x, incx, y, incy, w_ee, lda)
       end do
     end do

     ABI_FREE(rhotwg_mn_x)
     ABI_FREE(rhotwg_mn_c)
     ABI_FREE(w_rhotwg_mn_c)
   end do ! ik_bz

    if (print_time_kp) then
      call inds2str(1, "My kp-point", my_ikp, my_nkp, fs%nkfs, msg)
      call cwtime_report(msg, cpu_kp, wall_kp, gflops_kp); if (ikp_bz == LOG_MODK) call wrtout(std_out, "...", do_flush=.True.)
    end if

    ABI_FREE(ug_kp)
    ABI_FREE(kpg_kp)
 end do ! ikp_bz
 end associate

 ABI_FREE(ur_nkp)
 ABI_FREE(ur_mk)
 !end do ! my_is

 call xmpi_sum(w_ee, comm, ierr)

 call cwtime_report(" wkk calculation", cpu_all, wall_all, gflops_all, end_str=ch10)

 ! Set vqk_completed to 1 so that we can easily check if restarted is needed.
 !if (wkk_mode /= 0) then
 !if (my_rank == master) then
 !  NCF_CHECK(nf90_put_var(root_ncid, root_vid("vqk_completed"), 1))
 !end if
 NCF_CHECK(nf90_close(root_ncid))
 call xmpi_barrier(comm)
 !end if

 ! Output some of the results to ab_out for testing purposes
 !call gstore%print_for_abitests(dtset)

 ! Free memory
 ABI_FREE(kg_kp)
 ABI_FREE(kg_k)
 ABI_FREE(cwork_ur)
 ABI_FREE(work)
 ABI_FREE(gbound_kp)
 ABI_FREE(gbound_k)
 ABI_FREE(gbound_c)
 ABI_FREE(gbound_x)
 ABI_FREE(vc_sqrt_gx)
 ABI_FREE(w_ee)
 ABI_FREE(e_mesh)
 ABI_FREE(e_args)
 ABI_FREE(wgt_mk)
 ABI_FREE(wgt_nkp)

 call wfd%free(); call vcp%free(); call qmesh%free(); call gsph_x%free(); call gsph_c%free(); call hscr%free(); call edos%free()
 call epsm1%free()

 do spin=1,ebands%nsppol
   call fstab(spin)%free()
 end do
 ABI_FREE(fstab)

 call xmpi_barrier(comm) ! This to make sure that the parallel output of GSTORE is completed
 call cwtime_report(" wkk_run: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

contains

subroutine inds2str(ilevel, prefix, my_ik, my_nk, nk_tot, out_str)
 character(len=*),intent(in) :: prefix
 integer,intent(in) :: ilevel, my_ik, my_nk, nk_tot
 character(len=*),intent(out) :: out_str

 out_str = sjoin(prefix, itoa(my_ik), "/", itoa(my_nk), "[", itoa(nk_tot), "]")
 out_str = repeat(' ', 4 * ilevel) // trim(out_str)
end subroutine  inds2str

integer function root_vid(var_name)
  character(len=*),intent(in) :: var_name
  root_vid = nctk_idname(root_ncid, var_name)
end function root_vid

!integer function spin_vid(var_name)
!  character(len=*),intent(in) :: var_name
!  spin_vid = nctk_idname(spin_ncid, var_name)
!end function spin_vid

end subroutine wkk_run
!!***

end module m_wkk
!!***

!!****m* ABINIT/m_eph_path
!! NAME
!!  m_eph_path
!!
!! FUNCTION
!!  Compute e-ph matrix elements along a path
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
! use m_fft
 use m_hamiltonian
! use m_pawcprj
! use m_wfd
! use m_hdr
! use m_sigtk
! use m_ephtk
! use netcdf
! use m_nctk
 use m_dtset
 use m_dtfil
! use m_clib
! use m_mkffnl
!
 use defs_abitypes,    only : mpi_type
 use defs_datatypes,   only : ebands_t, pseudopotential_type
 use m_time,           only : cwtime, cwtime_report, timab, sec2str
! use m_fstrings,       only : itoa, ftoa, sjoin, ktoa, ltoa, strcat
! use m_numeric_tools,  only : arth, c2r, get_diag, linfit, iseven, simpson_cplx, simpson, print_arr, inrange
! use m_fftcore,        only : ngfft_seq, sphereboundary, get_kg, kgindex
! use m_cgtools,        only : cg_zdotc, cg_real_zdotc, cg_zgemm, fxphas_seq
 use m_crystal,        only : crystal_t
 use m_getgh1c,        only : getgh1c, rf_transgrid_and_pack, getgh1c_setup
 use m_pawang,         only : pawang_type
 use m_pawrad,         only : pawrad_type
 use m_pawtab,         only : pawtab_type
! use m_pawrhoij,       only : pawrhoij_type
 use m_pawfgr,         only : pawfgr_type
! use m_phonons,        only : phstore_t, phstore_new
 use m_pstat,          only : pstat_t
 use m_cgwf,           only : nscf_t
 use m_bz_mesh,        only : kpath_t, kpath_new  !isamek, make_path,

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
!!  Compute e-ph matrix elements along a path.
!!
!! INPUTS
!! dtset<dataset_type>=All input variables for this dataset.
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
!!
!! SOURCE

subroutine eph_path_run(dtfil, dtset, cryst, dvdb, ifc, &
                        pawfgr, pawang, pawrad, pawtab, psps, mpi_enreg, comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(dvdb_t),intent(inout) :: dvdb
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(pawfgr_type),intent(in) :: pawfgr
 type(ifc_type),intent(in) :: ifc
 type(mpi_type),intent(inout) :: mpi_enreg
!arrays
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: istwfk_1 = 1 ! ,tim_getgh1c1 = 1, berryopt0 = 0, ider0 = 0, idir0 = 0,
! integer,parameter :: useylmgr = 0, useylmgr1 =0, master = 0, ndat1 = 1, cplex1 = 1, pawread0 = 0
 integer :: spin, iqpt, nband, ierr, npw_k, npw_kq, my_rank, nprocs, n1, n2, n3, n4, n5, n6, cplex
 integer :: natom, natom3, nsppol, nspden, nspinor, cnt, qptopt
 integer :: nfft,nfftf,mgfft,mgfftf !, nkpg, nkpg1, qbuf_size, iqbuf_cnt, root_ncid, spin_ncid, ncerr
 real(dp) :: cpu_all,wall_all,gflops_all
 logical :: qq_is_gamma, use_ftinterp
! complex(dp) :: cfact,dka,dkap,dkpa,dkpap, cnum, sig_cplx, cfact2
 type(gs_hamiltonian_type) :: gs_hamk, gs_hamkq
 type(nscf_t) :: nscf
 type(kpath_t) :: qpath
 character(len=5000) :: msg
!!arrays
 integer :: units(2), ngfft(18),ngfftf(18)
 integer,allocatable :: kg_k(:,:), kg_kq(:,:), qmap_symrec(:,:)
 real(dp) :: kpt(3), qpt(3), kq(3)
 real(dp),allocatable :: cg_k(:,:,:), cg_kq(:,:,:), gsc_k(:,:,:), gsc_kq(:,:,:),eig_k(:), eig_kq(:)
 real(dp),allocatable :: v1scf(:,:,:,:), vlocal1(:,:,:,:,:) ! gkq_atm(:,:,:,:),gkq_nu(:,:,:,:)
! !real(dp),allocatable :: phfreqs_qibz(:,:), pheigvec_qibz(:,:,:,:), eigvec_qpt(:,:,:)
! real(dp) :: ylmgr_dum(1,1,1)
!************************************************************************

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

 call cwtime(cpu_all, wall_all, gflops_all, "start")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 units = [std_out, ab_out]

 !call pstat%from_pid(); call pstat%print([std_out], header="Memory at the beginning of eph_path")

 !! Copy important dimensions
 natom = cryst%natom; natom3 = 3 * natom; nsppol = dtset%nsppol; nspinor = dtset%nspinor
 nspden = dtset%nspden

 ! Build q-path.
 qpath = kpath_new(dtset%ph_qpath(:,1:dtset%ph_nqpath), cryst%gprimd, dtset%ph_ndivsm)

 ! Load KS potential from file.
 call nscf%init(dtset, dtfil, cryst, comm)
 nband = dtset%mband

 ! FFT meshes (taken from the GS POT file)
 ngfft = nscf%ngfft; ngfftf = nscf%ngfftf
 nfftf = product(ngfftf(1:3)); mgfftf = maxval(ngfftf(1:3))
 nfft = product(ngfft(1:3)) ; mgfft = maxval(ngfft(1:3))
 n1 = ngfft(1); n2 = ngfft(2); n3 = ngfft(3); n4 = ngfft(4); n5 = ngfft(5); n6 = ngfft(6)

 ! Open the DVDB file
 call dvdb%open_read(ngfftf, xmpi_comm_self)

 qptopt = dtset%kptopt; if (dtset%qptopt /= 0) qptopt = dtset%qptopt

 call dvdb%need_ftinterp(qpath%npts, qpath%points, qptopt, qmap_symrec, use_ftinterp)
 if (.not. use_ftinterp .and. dtset%eph_use_ftinterp /= 0) then
   ABI_WARNING("Enforcing FT interpolation for q-points even if it's not strictly needed.")
   use_ftinterp = .True.
 end if
 ABI_FREE(qmap_symrec)

#if 0
 if (use_ftinterp) then
   call wrtout(units, " Cannot find all q-points in the DVDB --> Activating Fourier interpolation.")
   call dvdb%ftinterp_setup(dtset%ddb_ngqpt, qptopt, 1, dtset%ddb_shiftq, nfftf, ngfftf, xmpi_comm_self)

   ! Build q-cache in the *dense* IBZ using the global mask qselect and itreat_qibz.
   ABI_MALLOC(itreat_qibz, (gstore%nqibz))
   ABI_MALLOC(qselect, (gstore%nqibz))
   qselect = 0; itreat_qibz = 0
   call dvdb%ftqcache_build(nfftf, ngfftf, gstore%nqibz, gstore%qibz, dtset%dvdb_qcache_mb, qselect, itreat_qibz, gstore%comm)
   ABI_FREE(itreat_qibz)
   ABI_FREE(qselect)

 else
   call wrtout(units, " DVDB file contains all q-points in the IBZ --> Reading DFPT potentials from file.")
   ! Need to translate itreat_qibz into itreatq_dvdb.
   ! FIXME: Not used
   ABI_ICALLOC(qselect, (dvdb%nqpt))
   ABI_ICALLOC(itreatq_dvdb, (dvdb%nqpt))
   !do iq_ibz=1,gstore%nqibz
   !  if (itreat_qibz(iq_ibz) == 0) cycle
   !  db_iqpt = qibz2dvdb(iq_ibz)
   !  ABI_CHECK(db_iqpt /= -1, sjoin("Could not find IBZ q-point:", ktoa(gstore%qibz(:, iq_ibz)), "in the DVDB file."))
   !  itreatq_dvdb(db_iqpt) = 1
   !end do
   call dvdb%qcache_read(nfftf, ngfftf, dtset%dvdb_qcache_mb, qselect, itreatq_dvdb, comm)
   ABI_FREE(qselect)
   ABI_FREE(itreatq_dvdb)
 end if
#endif

 do spin=1,dtset%nsppol
   ! Compute psi_nk
   kpt = zero
   call nscf%solve(spin, kpt, istwfk_1, nband, cryst, dtset, dtfil, psps, pawtab, pawfgr, &
                   npw_k, kg_k, cg_k, gsc_k, eig_k, gs_hamk, msg, ierr) ! out
   ABI_CHECK(ierr == 0, msg)

   do iqpt=1,qpath%npts
     qpt = qpath%points(:,iqpt)
     qq_is_gamma = sum(qpt**2) < tol14
     kq = kpt + qpt

     ! Compute psi_mkq
     if (.not. qq_is_gamma) then
       call nscf%solve(spin, kq, istwfk_1, nband, cryst, dtset, dtfil, psps, pawtab, pawfgr, &
                       npw_kq, kg_kq, cg_kq, gsc_kq, eig_kq, gs_hamkq, msg, ierr) ! out
       ABI_CHECK(ierr == 0, msg)
     else
       call alloc_copy(kg_k, kg_kq)
       call alloc_copy(eig_k, eig_kq)
       call alloc_copy(cg_k, cg_kq)
       call alloc_copy(gsc_k, gsc_kq)
       call copy_hamiltonian(gs_hamkq, gs_hamk)
     end if

#if 0
     ! Version with qcache.
     if (use_ftinterp) then
       call dvdb%get_ftqbz(cryst, qq_bz, qq_ibz, gqk%my_q2ibz(:, my_iq), cplex, nfftf, ngfftf, v1scf, &
                           gqk%pert_comm%value)
     else
       ! Read and reconstruct the dvscf potentials for qpt and my_npert perturbations.
       db_iqpt = dvdb%findq(qq_ibz)
       ABI_CHECK(db_iqpt /= -1, sjoin("Could not find symmetric of q-point:", ktoa(qq_bz), "in DVDB file."))
       ! The first entry in mapc_qq2dvdb gives the index in dvdb%qpts.
       ! The other entries in mapc_qq are OK as they refer to symmetries.
       mapc_qq2dvdb = gqk%my_q2ibz(:, my_iq); mapc_qq2dvdb(1) = db_iqpt
       call dvdb%readsym_qbz(cryst, qq_bz, mapc_qq2dvdb, cplex, nfftf, ngfftf, v1scf, gqk%pert_comm%value)
     end if

     ABI_FREE(v1scf)
     ! Allocate vlocal1 with correct cplex. Note nvloc and my_npert.
     ABI_MALLOC(vlocal1, (cplex*n4, n5, n6, gs_hamkq%nvloc, my_npert))
     ABI_FREE(vlocal1)

     ! Set up local potential vlocal1 with proper dimensioning from vtrial1 taking into account the spin.
     do my_ip=1,my_npert
       call rf_transgrid_and_pack(spin, nspden, psps%usepaw, cplex, nfftf, nfft, ngfft, gs_hamkq%nvloc,&
                                  pawfgr, mpi_enreg, dummy_vtrial, v1scf(:,:,:,my_ip), vlocal, vlocal1(:,:,:,:,my_ip))
     end do

     ! Continue to initialize the GS Hamiltonian
     !call gs_hamkq%load_spin(spin, vlocal=vlocal, with_nonlocal=.true.)
#endif

     ! Loop over atomic perturbations, apply H1_{kappa, alpha} and compute e-ph matrix elements.
#if 0
     ! Loop over my atomic perturbations and compute gkq_atm.
     gkq_atm = zero
     do my_ip=1,my_npert
       idir = dvdb%my_pinfo(1, my_ip); ipert = dvdb%my_pinfo(2, my_ip); ipc = dvdb%my_pinfo(3, my_ip)

       ! Prepare application of the NL part.
       call init_rf_hamiltonian(cplex, gs_hamkq, ipert, rf_hamkq, has_e1kbsc=.true.)

       call rf_hamkq%load_spin(spin, vlocal1=vlocal1(:,:,:,:,my_ip), with_nonlocal=.true.)

       ! This call is not optimal because there are quantities in out that do not depend on idir,ipert
       call getgh1c_setup(gs_hamkq, rf_hamkq, dtset, psps, kk_bz, kq_bz, idir, ipert, &             ! In
                          cryst%natom, cryst%rmet, cryst%gprimd, cryst%gmet, istwf_k, &             ! In
                          npw_k, npw_kq, useylmgr1, kg_k, ylm_k, kg_kq, ylm_kq, ylmgr_kq, &         ! In
                          dkinpw, nkpg, nkpg1, kpg_k, kpg1_k, kinpw1, ffnlk, ffnl1, ph3d, ph3d1, &  ! Out
                          reuse_kpg_k=1, reuse_kpg1_k=1, reuse_ffnlk=1, reuse_ffnl1=1)              ! Reuse some arrays

       ! Calculate dvscf * psi_k, results stored in h1kets_kq on the k+q sphere.
       ! Compute H(1) applied to GS wavefunction Psi(0)
       do in_k=1,nband_k
         band_k = in_k + bstart_k - 1
         eig0nk = ebands%eig(band_k, ik_ibz, spin)
         ! Use scissor shift on 0-order eigenvalue
         eshift = eig0nk - dtset%dfpt_sciss

         call getgh1c(berryopt0, kets_k(:,:,in_k), cwaveprj0, h1kets_kq(:,:,in_k), &
                      grad_berry, gs1c, gs_hamkq, gvnlx1, idir, ipert, eshift, mpi_enreg, optlocal, &
                      optnl, opt_gvnlx1, rf_hamkq, sij_opt, tim_getgh1c, usevnl)
       end do

       call rf_hamkq%free()
       ABI_FREE(kinpw1)
       ABI_FREE(dkinpw)
       ABI_FREE(ph3d)
       ABI_SFREE(ph3d1)

       ! Calculate <psi_{k+q,j}|dvscf_q*psi_{k,i}> for this perturbation. No need to handle istwf_kq because it's always 1.
       do in_k=1,nband_k
         do im_kq=1,nband_kq
           gkq_atm(:, im_kq, in_k, ipc) = cg_zdotc(npw_kq*nspinor, bras_kq(1,1,im_kq), h1kets_kq(1,1,in_k))
         end do
       end do

     end do ! my_ip
#endif

     ABI_FREE(kg_kq)
     ABI_FREE(eig_kq)
     ABI_FREE(cg_kq)
     ABI_FREE(gsc_kq)
     call gs_hamkq%free()
   end do ! iqpt

   ABI_FREE(kg_k)
   ABI_FREE(eig_k)
   ABI_FREE(cg_k)
   ABI_FREE(gsc_k)
   call gs_hamk%free()
 end do ! spin

 call nscf%free()
 call qpath%free()

 ! Average matrix elements over degenerate states (electrons at k, k+q, and phonons

 call cwtime_report(" eph_path: MPI barrier before returning.", cpu_all, wall_all, gflops_all, end_str=ch10, comm=comm)

end subroutine eph_path_run
!!***

end module m_eph_path
!!***

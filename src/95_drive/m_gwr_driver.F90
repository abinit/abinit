!!****m* ABINIT/m_gwr_driver
!! NAME
!!  m_gwr_driver
!!
!! FUNCTION
!!  Driver for GWR calculations
!!
!! COPYRIGHT
!!  Copyright (C) 2021-2022 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  APACHE license version 2.0, see ~abinit/COPYING
!!  or https://www.apache.org/licenses/LICENSE-2.0 .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_gwr_driver

#ifdef HAVE_MPI2
 use mpi
#endif
 use defs_basis
 use defs_wvltypes
 use m_errors
 use m_abicore
 use m_xmpi
 use m_xomp
 use m_hdr
 use libxc_functionals
 use m_crystal
 use m_ebands
 use m_dtset
 use m_dtfil
 use m_wfk
 use m_distribfft
 use m_nctk
 use, intrinsic :: iso_c_binding

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : MPI_type
 use m_time,            only : timab
 use m_io_tools,        only : file_exists, open_file, get_unit, iomode_from_fname
 use m_time,            only : cwtime, cwtime_report, sec2str
 use m_fstrings,        only : strcat, sjoin, ftoa, itoa, string_in
 use m_slk,             only : matrix_scalapack
 use m_fftcore,         only : print_ngfft, get_kg
 use m_fft,             only : fourdp
 use m_ioarr,           only : read_rhor
 use m_energies,        only : energies_type, energies_init
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, &
                               pawrhoij_inquire_dim, pawrhoij_symrhoij, pawrhoij_unpack
 use m_pawdij,          only : pawdij, symdij_all
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_pwaves_lmn,  only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawpwij,         only : pawpwff_t, pawpwff_init, pawpwff_free, paw_rho_tw_g
 use m_kg,              only : getph !, getcut
 use m_wfd,             only : test_charge
 use m_pspini,          only : pspini
 use m_paw_correlations,only : pawpuxinit
 use m_paw_dmft,        only : paw_dmft_type
 use m_paw_sphharm,     only : setsym_ylm
 use m_paw_mkrho,       only : denfgr
 use m_paw_nhat,        only : nhatgrid, pawmknhat
 use m_paw_tools,       only : chkpawovlp, pawprt
 use m_paw_denpot,      only : pawdenpot
 use m_paw_init,        only : pawinit, paw_gencond
 use m_pawcprj,         only : pawcprj_type, pawcprj_free, pawcprj_alloc ! , paw_overlap
 use m_ksdiago,         only : ugb_t
 !use m_fock,            only : fock_type, fock_init, fock_destroy, fock_ACE_destroy, fock_common_destroy, &
 !                              fock_BZ_destroy, fock_update_exc, fock_updatecwaveocc
 use m_mkrho,           only : prtrhomxmn
 use m_melemts,         only : melflags_t
 use m_setvtr,          only : setvtr
 use m_vhxc_me,         only : calc_vhxc_me
 use m_gwr,             only : gwr_t
 !use m_ephtk,          only : ephtk_update_ebands

 implicit none

 private
!!***

 public :: gwr_driver

#ifdef HAVE_MPI1
 include 'mpif.h'
#endif
!!***

contains
!!***

!!****f* m_gwr_driver/gwr_driver
!! NAME
!!  gwr_driver
!!
!! FUNCTION
!! Main routine for GWR calculations.
!!
!! INPUTS
!! acell(3)=Length scales of primitive translations (bohr)
!! codvsn=Code version
!! dtfil<datafiles_type>=Variables related to files.
!! dtset<dataset_type>=All input variables for this dataset.
!! pawang<pawang_type)>=PAW angular mesh and related data.
!! pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! psps<pseudopotential_type>=Variables related to pseudopotentials.
!!   Before entering the first time in the routine, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! xred(3,natom)=Reduced atomic coordinates.
!!
!! NOTES
!!
!! ON THE USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ... Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! SOURCE

subroutine gwr_driver(codvsn, dtfil, dtset, pawang, pawrad, pawtab, psps, xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master = 0, cplex1 = 1, ipert0 = 0, idir0 = 0, optrhoij1 = 1
 integer :: ii, comm, nprocs, my_rank, mgfftf, nfftf, omp_ncpus, work_size, nks_per_proc
 integer :: ierr, spin, ik_ibz, nband_k, iomode__, color, io_comm
 real(dp) :: eff, mempercpu_mb, max_wfsmem_mb, nonscal_mem
 real(dp) :: ecore, ecut_eff, ecutdg_eff, cpu, wall, gflops, diago_cpu, diago_wall, diago_gflops
 logical, parameter :: is_dfpt = .false.
 logical :: read_wfk, write_wfk, cc4s_task, rectangular
 character(len=500) :: msg
 character(len=fnlen) :: wfk_path, den_path, kden_path, out_path
 type(hdr_type) :: wfk_hdr, den_hdr, kden_hdr, owfk_hdr
 type(crystal_t) :: cryst, den_cryst, wfk_cryst
 type(ebands_t) :: ks_ebands, owfk_ebands
 type(pawfgr_type) :: pawfgr
 type(wvl_data) :: wvl
 type(mpi_type) :: mpi_enreg_seq
 type(gwr_t) :: gwr
 type(wfk_t) :: owfk
!arrays
 real(dp), parameter :: k0(3) = zero
 integer :: cplex, cplex_dij, cplex_rhoij
 integer :: gnt_option,has_dijU,has_dijso,ider,izero
 integer :: istep, moved_atm_inside, moved_rhor, n3xccc, sc_mode
 !integer :: ngrvdw,nhatgrdim,nkxc,nkxc1,nprocs,nscf,nspden_rhoij,nzlmopt,optene
 integer :: ndij !,ndim,nfftf,nfftf_tot,nkcalc,gwc_nfft,gwc_nfftot,gwx_nfft,gwx_nfftot
 integer :: ngrvdw, nhatgrdim, nkxc, nkxc1, nspden_rhoij, optene, nzlmopt
 integer :: optcut, optgr0, optgr1, optgr2, optrad, psp_gencond, option
 integer :: rhoxsp_method, usexcnhat !, use_umklp
 real(dp) :: compch_fft, compch_sph !,r_s,rhoav,alpha
 !real(dp) :: drude_plsmf !,my_plsmf,ecut_eff,ecutdg_eff,ehartree
 real(dp) :: gsqcutc_eff, gsqcutf_eff, gsqcut_shp
 real(dp) :: vxcavg !,vxcavg_qp ucvol,
 real(dp) :: gw_gsq !, gsqcut, gwc_gsq, gwx_gsq,
 logical :: call_pawinit
 type(energies_type) :: KS_energies
 type(melflags_t) :: KS_mflags
 type(paw_dmft_type) :: Paw_dmft
 type(ugb_t) :: ugb
 type(xmpi_pool2d_t) :: diago_pool
 !type(fock_type),pointer :: fock => null()
!arrays
 integer :: ngfftc(18),ngfftf(18),units(2) ! gwc_ngfft(18),gwx_ngfft(18),
 integer,allocatable :: nq_spl(:), l_size_atm(:) !,qp_vbik(:,:),tmp_gfft(:,:)
 integer,allocatable :: tmp_kstab(:,:,:), npwarr_ik(:), gvec_(:,:), istwfk_ik(:), nband_iks(:,:)
 real(dp) :: strsxc(6), diago_info(3, dtset%nkpt, dtset%nsppol),tsec(2)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grvdw(:,:),qmax(:)
 real(dp),allocatable :: ks_nhat(:,:),ks_nhatgr(:,:,:),ks_rhog(:,:)
 real(dp),allocatable :: ks_rhor(:,:),ks_vhartr(:), ks_vtrial(:,:), ks_vxc(:,:)
 real(dp),allocatable :: ks_taur(:,:) !, ks_vxctau(:,:), xcctau3d(:)
 real(dp),allocatable :: kxc(:,:), ph1d(:,:), ph1df(:,:) !qp_kxc(:,:),
 real(dp),allocatable :: vpsp(:), xccc3d(:), dijexc_core(:,:,:) !, dij_hf(:,:,:)
 real(dp),allocatable :: eig_k(:), occ_k(:)
 real(dp),contiguous,pointer :: cg_k_ptr(:,:)
 type(Paw_an_type),allocatable :: KS_paw_an(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:)
 type(pawpwff_t),allocatable :: Paw_pwff(:)
 !type(pawcprj_type),allocatable :: cprj_k(:,:)

!************************************************************************

 ! This part performs the initialization of the basic objects used to perform e-ph calculations:
 !
 !     1) Crystal structure `cryst`
 !     2) Ground state band energies: `ks_ebands`
 !     5) Pseudos and PAW basic objects.
 !
 ! Once we have these objects, we can call specialized routines for e-ph calculations.
 ! Notes:
 !
 !   * Any modification to the basic objects mentioned above should be done here (e.g. change of efermi)
 !   * This routines shall not allocate big chunks of memory. The CPU-demanding sections should be
 !     performed in the subdriver that will employ different MPI distribution schemes optimized for that particular task.

 ! abirules!
 if (.False.) write(std_out,*)xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 units(:) = [std_out, ab_out]

 call cwtime(cpu, wall, gflops, "start")

! write(msg,'(7a)')&
! ' SIGMA: Calculation of the GW corrections ',ch10,ch10,&
! ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
! ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
! call wrtout(units, msg)
!
#if defined HAVE_GW_DPC
 write(msg,'(a,i2,a)')'.Using double precision arithmetic; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic; gwpc = ',gwpc,ch10
#endif
 call wrtout(units, msg)

 ! autoparal section
 ! TODO: This just to activate autoparal in AbiPy. Lot of things should be improved.
 if (dtset%max_ncpus /= 0) then
   write(ab_out,'(a)')"--- !Autoparal"
   write(ab_out,"(a)")"# Autoparal section for GWR runs"
   write(ab_out,"(a)")   "info:"
   write(ab_out,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ab_out,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ab_out,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ab_out,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ab_out,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ab_out,"(a,i0)")"    mband: ",dtset%mband
   write(ab_out,"(a,i0)")"    gwr_task: ",trim(dtset%gwr_task)

   work_size = dtset%nkpt * dtset%nsppol
   ! Non-scalable memory in Mb i.e. memory that is not distributed with MPI.
   nonscal_mem = zero
   max_wfsmem_mb = (two * dp * dtset%mpw * dtset%mband * dtset%nkpt * dtset%nsppol * dtset%nspinor * b2Mb) * 1.1_dp

   ! List of configurations.
   ! Assuming an OpenMP implementation with perfect speedup!
   write(ab_out,"(a)")"configurations:"

   do ii=1,dtset%max_ncpus
     nks_per_proc = work_size / ii
     nks_per_proc = nks_per_proc + mod(work_size, ii)
     eff = (one * work_size) / (ii * nks_per_proc)
     ! Add the non-scalable part and increase by 10% to account for other datastructures.
     mempercpu_mb = (max_wfsmem_mb + nonscal_mem) * 1.1_dp

     do omp_ncpus=1,1 !xomp_get_max_threads()
       write(ab_out,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
       write(ab_out,"(a,i0)")"      mpi_ncpus: ",ii
       write(ab_out,"(a,i0)")"      omp_ncpus: ",omp_ncpus
       write(ab_out,"(a,f12.9)")"      efficiency: ",eff
       write(ab_out,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
     end do
   end do
   write(ab_out,'(a)')"..."
   ABI_STOP("Stopping now!")
 end if

 cryst = dtset%get_crystal(img=1)

 ! Some variables need to be initialized/nullify at start
 usexcnhat = 0
 call energies_init(KS_energies)

 den_path = dtfil%fildensin; wfk_path = dtfil%fnamewffk; kden_path = dtfil%filkdensin

 if (my_rank == master) then
   ! Initialize filenames. Accept files in Fortran or in netcdf format.
   if (nctk_try_fort_or_ncfile(den_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find DEN file:", den_path, ". Error:", msg))
   end if
   call wrtout(units, sjoin("- Reading GS density from: ", den_path))

   if (dtset%usekden == 1) then
     if (nctk_try_fort_or_ncfile(kden_path, msg) /= 0) then
       ABI_ERROR(sjoin("Cannot find KDEN file:", kden_path, ". Error:", msg))
     end if
     call wrtout(units, sjoin("- Reading KDEN kinetic energy density from: ", kden_path))
   end if
   call wrtout(ab_out, ch10//ch10)
 end if ! master

 ! Broadcast filenames (needed if we are using netcdf files)
 call xmpi_bcast(den_path, master, comm, ierr)
 call xmpi_bcast(kden_path, master, comm, ierr)

 ! TODO: FFT meshes for DEN/POT should be initialized from the DEN file instead of the dtset.
 ! Interpolating the DEN indeed breaks degeneracies in the vxc matrix elements.

 call pawfgr_init(pawfgr, dtset, mgfftf, nfftf, ecut_eff, ecutdg_eff, ngfftc, ngfftf, &
                  gsqcutc_eff=gsqcutc_eff, gsqcutf_eff=gsqcutf_eff, gmet=cryst%gmet, k0=k0)

 call print_ngfft(ngfftc, header='Coarse FFT mesh for the wavefunctions')
 call print_ngfft(ngfftf, header='Dense FFT mesh for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg_seq)
 call init_distribfft_seq(mpi_enreg_seq%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 call init_distribfft_seq(mpi_enreg_seq%distribfft, 'f', ngfftf(2), ngfftf(3), 'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset, dtfil, ecore, psp_gencond, gsqcutc_eff, gsqcutf_eff, pawrad, pawtab, psps, cryst%rprimd, comm_mpi=comm)

 ! ============================
 ! ==== PAW initialization ====
 ! ============================
 if (dtset%usepaw == 1) then
   call chkpawovlp(cryst%natom, cryst%ntypat, dtset%pawovlp, pawtab, cryst%rmet, cryst%typat, cryst%xred)

   cplex_dij = dtset%nspinor; cplex = 1; ndij = 1

   ABI_MALLOC(ks_pawrhoij, (cryst%natom))
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij, nspden_rhoij=nspden_rhoij, &
                             nspden=dtset%nspden, spnorb=dtset%pawspnorb, cpxocc=dtset%pawcpxocc)
   call pawrhoij_alloc(ks_pawrhoij, cplex_rhoij, nspden_rhoij, dtset%nspinor, dtset%nsppol, cryst%typat, pawtab=pawtab)

   ! Test if we have to call pawinit
   gnt_option = 1; if (dtset%pawxcdev == 2 .or. (dtset%pawxcdev == 1 .and. dtset%positron /= 0)) gnt_option = 2
   call paw_gencond(dtset, gnt_option, "test", call_pawinit)
   !call_pawinit = .True.

   if (psp_gencond == 1 .or. call_pawinit) then
     call timab(553, 1, tsec)
     gsqcut_shp = two * abs(dtset%diecut) * dtset%dilatmx**2 / pi**2
     call pawinit(dtset%effmass_free, gnt_option, gsqcut_shp, zero, dtset%pawlcutd, dtset%pawlmix, &
                  psps%mpsang, dtset%pawnphi, cryst%nsym, dtset%pawntheta, pawang, pawrad, &
                  dtset%pawspnorb, pawtab, dtset%pawxcdev, dtset%ixc, dtset%usepotzero)
     call timab(553,2,tsec)

     ! Update internal values
     call paw_gencond(dtset, gnt_option, "save", call_pawinit)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:cryst%ntypat)%has_kij   = 2
     if (pawtab(1)%has_nabla==1) pawtab(1:cryst%ntypat)%has_nabla = 2
   end if

   psps%n1xccc = maxval(pawtab(1:cryst%ntypat)%usetcore)

   ! Initialize optional flags in Pawtab to zero
   ! Cannot be done in Pawinit since the routine is called only if some parts. are changed
   pawtab(:)%has_nabla = 0
   pawtab(:)%lamb_shielding = zero

   call setsym_ylm(cryst%gprimd, pawang%l_max-1, cryst%nsym, dtset%pawprtvol, cryst%rprimd, cryst%symrec, pawang%zarot)

   ! Initialize and compute data for DFT+U
   Paw_dmft%use_dmft = Dtset%usedmft
   call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla, &
     is_dfpt,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,dtset%nspinor,Cryst%ntypat,dtset%optdcmagpawu,Pawang,Dtset%pawprtvol, &
     Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu,dtset%ucrpa)

   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   !call pawrhoij_copy(wfk_hdr%pawrhoij, KS_Pawrhoij)

   !  Evaluate form factor of radial part of phi.phj-tphi.tphj.
   gw_gsq = max(Dtset%ecutsigx, Dtset%ecuteps) / (two*pi**2)

   ! Set up q-grid, make qmax 20% larger than largest expected.
   ABI_MALLOC(nq_spl, (Psps%ntypat))
   ABI_MALLOC(qmax, (Psps%ntypat))
   qmax = SQRT(gw_gsq)*1.2d0
   nq_spl = Psps%mqgrid_ff
   ! write(std_out,*)"using nq_spl",nq_spl,"qmax=",qmax

   rhoxsp_method = 1  ! Arnaud-Alouani (default in sigma)
   !rhoxsp_method = 2 ! Shiskin-Kresse
   if (Dtset%pawoptosc /= 0) rhoxsp_method = Dtset%pawoptosc

   ABI_MALLOC(paw_pwff, (psps%ntypat))
   call pawpwff_init(Paw_pwff, rhoxsp_method, nq_spl, qmax, cryst%gmet, pawrad, pawtab, psps)

   ABI_FREE(nq_spl)
   ABI_FREE(qmax)

   ! Variables/arrays related to the fine FFT grid
   ABI_CALLOC(ks_nhat, (nfftf, Dtset%nspden))

   ABI_MALLOC(pawfgrtab, (cryst%natom))
   call pawtab_get_lsize(pawtab, l_size_atm, cryst%natom, cryst%typat)

   cplex = 1
   call pawfgrtab_init(pawfgrtab, cplex, l_size_atm, dtset%nspden, dtset%typat)
   ABI_FREE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat = maxval(Pawtab(:)%usexcnhat)
   ! * 0 if Vloc in atomic data is Vbare    (Blochl's formulation)
   ! * 1 if Vloc in atomic data is VH(tnzc) (Kresse's formulation)
   call wrtout(std_out, sjoin(' using usexcnhat: ', itoa(usexcnhat)))
   !
   ! Identify parts of the rectangular grid where the density has to be calculated
   optcut = 0; optgr0 = Dtset%pawstgylm; optgr1 = 0; optgr2 = 0; optrad = 1 - Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(cryst%atindx1, cryst%gmet, cryst%natom, cryst%natom, cryst%nattyp, ngfftf, cryst%ntypat,&
    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   call pawfgrtab_print(Pawfgrtab,Cryst%natom,unit=std_out,prtvol=Dtset%pawprtvol)

 else
   ABI_MALLOC(Paw_pwff, (0))
   ABI_MALLOC(Pawfgrtab, (0))
 end if ! End of PAW Initialization

 ! Allocate these arrays anyway, since they are passed to subroutines.
 ABI_MALLOC_IFNOT(ks_nhat, (nfftf, 0))
 ABI_MALLOC_IFNOT(dijexc_core, (1, 1, 0))

 !=============================================
 ! Read density and compare crystal structures
 ! ============================================
 ABI_MALLOC(ks_rhor, (nfftf, dtset%nspden))

 call read_rhor(den_path, cplex1, dtset%nspden, nfftf, ngfftf, dtset%usepaw, mpi_enreg_seq, ks_rhor, &
                den_hdr, ks_pawrhoij, comm, allow_interp=.False.)

 den_cryst = den_hdr%get_crystal()
 if (cryst%compare(den_cryst, header=" Comparing input crystal with DEN crystal") /= 0) then
   ABI_ERROR("Crystal structure from input and from DEN file do not agree! Check messages above!")
 end if
 call den_cryst%free(); call den_hdr%free()

 ABI_MALLOC(ks_taur, (nfftf, dtset%nspden * dtset%usekden))
 if (dtset%usekden == 1) then
   call read_rhor(kden_path, cplex1, dtset%nspden, nfftf, ngfftf, 0, mpi_enreg_seq, ks_taur, &
                  kden_hdr, ks_pawrhoij, comm, allow_interp=.False.)
   call kden_hdr%free()
   call prtrhomxmn(std_out, mpi_enreg_seq, nfftf, ngfftf, dtset%nspden, 1, ks_taur, optrhor=1, ucvol=cryst%ucvol)
 end if

 !========================================
 !==== Additional computation for PAW ====
 !========================================
 nhatgrdim = 0
 if (dtset%usepaw == 1) then
   ! Calculate the compensation charge nhat.
   if (Dtset%xclevel==2) nhatgrdim = usexcnhat * Dtset%pawnhatxc
   cplex = 1; ider = 2 * nhatgrdim; izero = 0
   if (nhatgrdim > 0) then
     ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,3*nhatgrdim))
   end if
   if (nhatgrdim == 0) then
     ABI_MALLOC(ks_nhatgr,(0,0,0))
   end if

   call pawmknhat(compch_fft,cplex,ider,idir0,ipert0,izero,Cryst%gprimd,&
                  Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
                  Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,k0,Cryst%rprimd,&
                  Cryst%ucvol,dtset%usewvl,Cryst%xred)

   ! === Evaluate onsite energies, potentials, densities ===
   !  * Initialize variables/arrays related to the PAW spheres.
   !  * Initialize also lmselect (index of non-zero LM-moments of densities).
   ABI_MALLOC(KS_paw_ij, (Cryst%natom))
   has_dijso = Dtset%pawspnorb; has_dijU = merge(0, 1, Dtset%usepawu == 0)

   call paw_ij_nullify(KS_paw_ij)
   call paw_ij_init(KS_paw_ij,cplex,Dtset%nspinor,Dtset%nsppol,&
     Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
     has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=1,&
     has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1, &
     has_dijfock=dtset%usefock)

   nkxc1 = 0
   ABI_MALLOC(KS_paw_an, (Cryst%natom))
   call paw_an_nullify(KS_paw_an)
   call paw_an_init(KS_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,0,Dtset%nspden,&
     cplex,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=1,has_vxctau=dtset%usekden)

   !  Calculate onsite vxc with and without core charge.
   nzlmopt=-1; option=0; compch_sph=greatest_real
   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert0,&
     Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
     Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
     Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
     Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)

 else
   ABI_MALLOC(ks_nhatgr, (0, 0, 0))
   ABI_MALLOC(ks_paw_ij, (0))
   ABI_MALLOC(ks_paw_an, (0))
 end if ! PAW

 !call test_charge(nfftf,ks_ebands%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
 !                 Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)

 ! For PAW, add the compensation charge on the FFT mesh, then get rho(G).
 ! NB: ks_nhat is already included in the density stored on file so we don't need to add it.to ks_rhor
 !if (dtset%usepaw==1) ks_rhor = ks_rhor + ks_nhat

 ! TODO: Overloaded interface with units or just change the API to accept units
 call prtrhomxmn(std_out, mpi_enreg_seq, nfftf, ngfftf, dtset%nspden, 1, ks_rhor, ucvol=cryst%ucvol)
 call prtrhomxmn(ab_out, mpi_enreg_seq, nfftf, ngfftf, dtset%nspden, 1, ks_rhor, ucvol=cryst%ucvol)

 if (dtset%usekden==1) then
   call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_taur,optrhor=1,ucvol=cryst%ucvol)
   call prtrhomxmn(ab_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_taur,optrhor=1,ucvol=cryst%ucvol)
 end if

 ! FFT n(r) --> n(g)
 ABI_MALLOC(ks_rhog, (2, nfftf))
 call fourdp(cplex1, ks_rhog, ks_rhor(:, 1), -1, mpi_enreg_seq, nfftf, 1, ngfftf, 0)

 ! Compute structure factor phases and large sphere cutoff
 ABI_MALLOC(ph1d, (2, 3 * (2 * Dtset%mgfft + 1) * Cryst%natom))
 ABI_MALLOC(ph1df, (2, 3 * (2 * mgfftf + 1) * Cryst%natom))

 call getph(cryst%atindx, cryst%natom, ngfftc(1), ngfftc(2), ngfftc(3), ph1d, cryst%xred)

 if (psps%usepaw == 1 .and. pawfgr%usefinegrid == 1) then
   call getph(cryst%atindx, cryst%natom, ngfftf(1), ngfftf(2), ngfftf(3), ph1df, cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

 ! The following steps have been gathered in the setvtr routine:
 !  - get Ewald energy and Ewald forces
 !  - compute local ionic pseudopotential vpsp
 !  - eventually compute 3D core electron density xccc3d
 !  - eventually compute vxc and vhartr
 !  - set up ks_vtrial
 !
 !*******************************************************************
 !**** NOTE THAT HERE Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
 !*******************************************************************

 ngrvdw = 0
 ABI_MALLOC(grvdw, (3, ngrvdw))
 ABI_MALLOC(grchempottn, (3, cryst%natom))
 ABI_MALLOC(grewtn, (3, cryst%natom))
 nkxc = 0
 if (dtset%nspden == 1) nkxc = 2
 if (dtset%nspden >= 2) nkxc = 3 ! check GGA and spinor, quite a messy part!!!
 ! In case of MGGA, fxc and kxc are not available and we dont need them (for now ...)
 if (dtset%ixc < 0 .and. libxc_functionals_ismgga()) nkxc = 0
 if (nkxc /= 0) then
   ABI_MALLOC(kxc, (nfftf, nkxc))
 end if

 n3xccc = 0; if (psps%n1xccc /= 0) n3xccc = nfftf
 ABI_MALLOC(xccc3d, (n3xccc))
 ABI_MALLOC(ks_vhartr, (nfftf))
 ABI_MALLOC(ks_vtrial, (nfftf, dtset%nspden))
 ABI_MALLOC(vpsp, (nfftf))
 ABI_MALLOC(ks_vxc, (nfftf, dtset%nspden))

 ! TODO: I don't think direct diago can be used with mega-GGA due to the functional derivative wrt KS states.
 ! TB-BK should be OK though.

 !ABI_MALLOC(ks_vxctau, (nfftf, dtset%nspden * dtset%usekden))
 !ABI_MALLOC(xcctau3d, (n3xccc * dtset%usekden))
 !ABI_FREE(ks_vxctau)
 !ABI_FREE(xcctau3d)

 optene = 4; moved_atm_inside = 0; moved_rhor = 0; istep = 1

 call setvtr(Cryst%atindx1,Dtset,KS_energies,cryst%gmet,cryst%gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,&
             istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
             Cryst%nattyp,nfftf,ngfftf,ngrvdw,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
             optene,pawrad,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,cryst%rmet,cryst%rprimd,strsxc,&
             Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,Wvl,xccc3d,Cryst%xred, &
             taur=ks_taur) !xcctau3d=xcctau3d, vxctau=ks_vxctau)

 ABI_FREE(grvdw)
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)

 !============================
 !==== Compute KS PAW Dij ====
 !============================
 if (dtset%usepaw == 1) then
   call timab(561,1,tsec)

   ! Calculate the unsymmetrized Dij.
   call pawdij(cplex,Dtset%enunit,Cryst%gprimd,ipert0,&
               Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
               Dtset%nspden,Cryst%ntypat,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
               Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
               k0,Dtset%spnorbscl,Cryst%ucvol,dtset%cellcharge(1),ks_vtrial,ks_vxc,Cryst%xred,&
               nucdipmom=Dtset%nucdipmom)

   ! Symmetrize KS Dij
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert0,&
                   Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,KS_paw_ij,Pawang,&
                   Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)

   ! Output the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Cryst%natom,KS_paw_ij,KS_Pawrhoij,Pawtab)
   call timab(561,2,tsec)
 end if

 call cwtime_report(" prepare gwr_driver_init", cpu, wall, gflops)

 if (string_in(dtset%gwr_task, "HDIAGO, HDIAGO_FULL, CC4S, CC4S_FULL")) then
   ! ==========================================
   ! Direct diagonalization of the Hamiltonian
   ! ==========================================
   ABI_MALLOC(nband_iks, (dtset%nkpt, dtset%nsppol))
   ABI_MALLOC(npwarr_ik, (dtset%nkpt))
   ABI_MALLOC(istwfk_ik, (dtset%nkpt))
   istwfk_ik = 1

   ! Compute npw_k from ecut so that we can update the header and redefine %mpw
   do ik_ibz=1,dtset%nkpt
     !if (dtset%istwfk(ik_ibz) == 2) istwfk_ik(ik_ibz) = 2  ! TODO: istwkf 2 is not yet supported.
     call get_kg(dtset%kptns(:,ik_ibz), istwfk_ik(ik_ibz), dtset%ecut, cryst%gmet, npwarr_ik(ik_ibz), gvec_)
     ABI_FREE(gvec_)
   end do
   dtset%mpw = maxval(npwarr_ik)

   ! CC4S does not need to output the WFK file.
   write_wfk = string_in(dtset%gwr_task, "HDIAGO, HDIAGO_FULL")

   ! Use input nband or min of npwarr_ik to set the number of bands.
   if (string_in(dtset%gwr_task, "HDIAGO, CC4S")) nband_iks(:,:) = maxval(dtset%nband)
   if (string_in(dtset%gwr_task, "HDIAGO_FULL, CC4S_FULL")) nband_iks(:,:) = minval(npwarr_ik)
   cc4s_task = string_in(dtset%gwr_task, "CC4S, CC4S_FULL")
   if (cc4s_task) then
     ABI_CHECK_IEQ(dtset%nkpt, 1, "CC4S interface does not support more than one k-point.")
     !ABI_CHECK(dtset%nkpt == 1 .and. all(abs(dtset%kptns(:,1)) < tol12), "CC4S requires Gamma-only sampling")
   end if

   ! Build header with new npwarr and nband.
   owfk_ebands = ebands_from_dtset(dtset, npwarr_ik, nband=nband_iks)
   owfk_ebands%eig = zero
   call hdr_init(owfk_ebands, codvsn, dtset, owfk_hdr, pawtab, 0, psps, wvl%descr)

   ! Change the value of istwfk taken from dtset.
   ABI_REMALLOC(owfk_hdr%istwfk, (dtset%nkpt))
   owfk_hdr%istwfk(:) = istwfk_ik

   ! Build pools to distribute (kpt, spin). Try to get rectangular grids in each pool to improve efficiency in slk diago.
   rectangular = .True.; if (dtset%nkpt == 1) rectangular = .False.
   call diago_pool%from_dims(dtset%nkpt, dtset%nsppol, comm, rectangular=rectangular)
   diago_info = zero

   !if (dtset%usefock == 1) then
   !  ! Initialize data_type fock for the calculation. See also m_scfcv_core.
   !  ABI_CHECK(dtset%usepaw == 0, "FOCK with PAW not coded")
   !  !call fock_from_wfk(fock, dtset, cryst, pawang, pawfgr, pawtab)
   !end if

   if (write_wfk) then
     ! Master writes header and Fortran record markers
     out_path = dtfil%fnameabo_wfk; if (dtset%iomode == IO_MODE_ETSF) out_path = nctk_ncify(out_path)
     iomode__ = iomode_from_fname(out_path)
     call wrtout(std_out, sjoin(" Writing wavefunctions to file:", out_path))
     if (my_rank == master) then
       call owfk%open_write(owfk_hdr, out_path, 0, iomode__, get_unit(), xmpi_comm_self, write_hdr=.True., write_frm=.True.)
       call owfk%close()
     end if
     call xmpi_barrier(comm)
   end if

   do spin=1,dtset%nsppol
     do ik_ibz=1,dtset%nkpt
       if (.not. diago_pool%treats(ik_ibz, spin)) cycle
       call cwtime(diago_cpu, diago_wall, diago_gflops, "start")

       nband_k = nband_iks(ik_ibz, spin)
       call ugb%from_diago(spin, istwfk_ik(ik_ibz), dtset%kptns(:,ik_ibz), dtset%ecut, nband_k, ngfftc, nfftf, &
                           dtset, pawtab, pawfgr, ks_paw_ij, cryst, psps, ks_vtrial, eig_k, diago_pool%comm%value)

       call cwtime(diago_cpu, diago_wall, diago_gflops, "stop")
       if (diago_pool%comm%me == 0) diago_info(1, ik_ibz, spin) = diago_wall
       call cwtime(diago_cpu, diago_wall, diago_gflops, "start")

       owfk_ebands%eig(1:nband_k, ik_ibz, spin) = eig_k(1:nband_k)

       if (write_wfk) then
         ! occupancies are set to zero. Client code is responsible for recomputing occ and fermie when reading this WFK.
         ABI_CALLOC(occ_k, (nband_k))
         color = merge(1, 0, ugb%my_nband > 0)
         call xmpi_comm_split(diago_pool%comm%value, color, diago_pool%comm%me, io_comm, ierr)

         if (ugb%my_nband > 0) then
           ABI_CHECK(all(shape(ugb%cg_k) == [2, ugb%npwsp, ugb%my_nband]), "Wrong shape")
           ABI_CHECK_IEQ(ugb%npw_k, owfk_hdr%npwarr(ik_ibz), "Wronk npw_k")
           call c_f_pointer(c_loc(ugb%cg_k), cg_k_ptr, shape=[2, ugb%npwsp * ugb%my_nband])

           ! Reopen file inside io_comm.
           call owfk%open_write(owfk_hdr, out_path, 0, iomode__, get_unit(), io_comm, &
                                write_hdr=.False., write_frm=.False.)

           !sc_mode = merge(xmpio_single, xmpio_collective, ugb%has_idle_procs)
           sc_mode = xmpio_collective
           call owfk%write_band_block([ugb%my_bstart, ugb%my_bstop], ik_ibz, spin, sc_mode, &
                                       kg_k=ugb%kg_k, cg_k=cg_k_ptr, &
                                       eig_k=owfk_ebands%eig(:, ik_ibz, spin), occ_k=occ_k)
           call owfk%close()
         end if
         call xmpi_comm_free(io_comm)
         ABI_FREE(occ_k)
       end if

       call cwtime(diago_cpu, diago_wall, diago_gflops, "stop")
       if (diago_pool%comm%me == 0) diago_info(2:3, ik_ibz, spin) = [diago_wall, dble(diago_pool%comm%nproc)]

       if (cc4s_task) call cc4s_gamma(spin, ik_ibz, dtset, dtfil, cryst, owfk_ebands, psps, pawtab, paw_pwff, ugb)

       ABI_FREE(eig_k)
       call ugb%free()
     end do ! ik_ibz
   end do ! spin

   call wrtout(std_out, " Direct diago completed by this MPI pool. Other pools might take more time if k != 0")

   call xmpi_sum_master(diago_info, master, comm, ierr)
   if (my_rank == master) then
     do spin=1,dtset%nsppol
       do ik_ibz=1,dtset%nkpt
         associate (info => diago_info(:, ik_ibz, spin))
         write(std_out, "(2(a,i0),5a,i0)") " ik_ibz: ", ik_ibz, ", spin: ", spin, &
           ", diago_wall: ", trim(sec2str(info(1))), ", io_wall: ", trim(sec2str(info(2))), ", nprocs: ", int(info(3))
         end associate
       end do
     end do
   end if

   ! Collect eigenvalues for the different k-points/spins
   do spin=1,dtset%nsppol
     do ik_ibz=1,dtset%nkpt
       if (diago_pool%treats(ik_ibz, spin) .and. diago_pool%comm%me /= 0) owfk_ebands%eig(:, ik_ibz, spin) = zero
     end do
   end do
   call xmpi_sum(owfk_ebands%eig, comm, ierr)

   call ebands_update_occ(owfk_ebands, dtset%spinmagntarget, prtvol=dtset%prtvol, fermie_to_zero=.False.)

   if (my_rank == master) then
     call ebands_print_gaps(owfk_ebands, ab_out, header="KS gaps after direct diagonalization")
     call ebands_print_gaps(owfk_ebands, std_out, header="KS gaps after direct diagonalization")
     if (cc4s_task) call cc4s_write_eigens(owfk_ebands, dtfil)
   end if

   ABI_FREE(npwarr_ik)
   ABI_FREE(istwfk_ik)
   ABI_FREE(nband_iks)
   call owfk_hdr%free(); call ebands_free(owfk_ebands); call diago_pool%free()

   ! Deallocate exact exchange data at the end of the calculation
   !if (dtset%usefock == 1) then
   !  if (fock%fock_common%use_ACE/=0) call fock_ACE_destroy(fock%fockACE)
   !  !call fock_common_destroy(fock%fock_common)
   !  call fock_BZ_destroy(fock%fock_BZ)
   !  call fock_destroy(fock)
   !  nullify(fock)
   !end if

 else
   ! ====================================================
   ! === This is the real GWR stuff once all is ready ===
   ! ====================================================
   ABI_CHECK(dtset%usepaw == 0, "PAW in GWR not yet implemented.")
   read_wfk = .True.
   if (read_wfk) then
     if (my_rank == master) then
       if (nctk_try_fort_or_ncfile(wfk_path, msg) /= 0) then
          ABI_ERROR(sjoin("Cannot find GS WFK file:", wfk_path, ". Error:", msg))
       end if
       call wrtout(units, sjoin("- Reading GS states from WFK file:", wfk_path))
     end if

     ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
     call xmpi_bcast(wfk_path, master, comm, ierr)

     ! Construct crystal and ks_ebands from the GS WFK file.
     ks_ebands = wfk_read_ebands(wfk_path, comm, out_hdr=wfk_hdr)
     call wfk_hdr%vs_dtset(dtset)

     wfk_cryst = wfk_hdr%get_crystal()
     if (cryst%compare(wfk_cryst, header=" Comparing input crystal with WFK crystal") /= 0) then
       ABI_ERROR("Crystal structure from input and from WFK file do not agree! Check messages above!")
     end if
     !call wfk_cryst%print(header="crystal structure from WFK file")
     call wfk_cryst%free()

     ! Make sure that ef is inside the gap if semiconductor.
     call ebands_update_occ(ks_ebands, dtset%spinmagntarget, prtvol=dtset%prtvol, fermie_to_zero=.True.)

     ! Here we change the GS bands (Fermi level, scissors operator ...)
     ! All the modifications to ebands should be done here.
     !call ephtk_update_ebands(dtset, ks_ebands, "Ground state energies")
   end if

   call gwr%init(dtset, dtfil, cryst, psps, pawtab, ks_ebands, mpi_enreg_seq, comm)
   if (gwr%idle_proc) goto 100

   !=== Calculate Vxc(b1,b2,k,s)=<b1,k,s|v_{xc}|b2,k,s> for all the states included in GW ===
   !  * This part is parallelized within wfd%comm since each node has all GW wavefunctions.
   !  * Note that vH matrix elements are calculated using the true uncutted interaction.

   call KS_mflags%reset()
   !if (rdm_update) then
   !KS_mflags%has_hbare=1
   !KS_mflags%has_kinetic=1
   !end if
   KS_mflags%has_vhartree=1
   KS_mflags%has_vxc     =1
   KS_mflags%has_vxcval  =1
   if (Dtset%usepawu /= 0  ) KS_mflags%has_vu      = 1
   if (Dtset%useexexch /= 0) KS_mflags%has_lexexch = 1
   if (Dtset%usepaw==1 .and. Dtset%gw_sigxcore == 1) KS_mflags%has_sxcore = 1
   !if (gwcalctyp<10        )  KS_mflags%only_diago = 1 ! off-diagonal elements only for SC on wavefunctions.
   KS_mflags%only_diago = 1 ! off-diagonal elements only for SC on wavefunctions.

   ! Load wavefunctions for Sigma_nk in gwr%kcalc_wfd.
   call gwr%load_kcalc_wfd(wfk_path, tmp_kstab)

   ! Compute gwr%ks_me matrix elements.
   if (.not. string_in(dtset%gwr_task, "RPA_ENERGY")) then
     ! FIXME: This routine allocates (nband, nband) matrices and should be rewritten!
     call calc_vhxc_me(gwr%kcalc_wfd, ks_mflags, gwr%ks_me, cryst, dtset, nfftf, ngfftf, &
                       ks_vtrial, ks_vhartr, ks_vxc, psps, pawtab, ks_paw_an, pawang, pawfgrtab, ks_paw_ij, dijexc_core, &
                       ks_rhor, usexcnhat, ks_nhat, ks_nhatgr, nhatgrdim, tmp_kstab, taur=ks_taur)
     if (my_rank == master) call gwr%ks_me%print(header="KS matrix elements", unit=std_out)
   end if

   ABI_FREE(tmp_kstab)

   if (read_wfk) then
     ! Build Green's function in imaginary-time from WFK file
     call gwr%read_ugb_from_wfk(wfk_path)
   else
     !call gwr%get_ugb_from_vtrial(ngfftf, ks_vtrial)
   end if

   ! Now call high-level routines depending on gwr_task.
   select case (dtset%gwr_task)
   case ("RPA_ENERGY")
     call gwr%rpa_energy()
   case ("G0W0")
     call gwr%run_g0w0()
   case ("G0V")
     call gwr%build_sigxme()
   case ("EGEW", "EGW0", "G0EW")
     call gwr%run_energy_scf()
   case default
     ABI_ERROR(sjoin("Invalid gwr_task:", dtset%gwr_task))
   end select
 end if

 !=====================
 !==== Free memory ====
 !=====================
100 call xmpi_barrier(comm)
 ABI_FREE(ks_nhat)
 ABI_FREE(ks_nhatgr)
 ABI_FREE(dijexc_core)
 call pawfgr_destroy(pawfgr)

 if (dtset%usepaw == 1) then
   ! Deallocation for PAW.
   call pawrhoij_free(ks_pawrhoij)
   ABI_FREE(ks_pawrhoij)
   call pawfgrtab_free(pawfgrtab)
   !ABI_FREE(pawfgrtab)
   call paw_ij_free(ks_paw_ij)
   ABI_FREE(ks_paw_ij)
   call paw_an_free(ks_paw_an)
   !ABI_FREE(ks_paw_an)
   call pawpwff_free(Paw_pwff)
 end if

 ABI_FREE(ph1d)
 ABI_FREE(ph1df)
 ABI_FREE(ks_rhor)
 ABI_FREE(ks_rhog)
 ABI_FREE(ks_taur)
 ABI_FREE(kxc)
 ABI_FREE(xccc3d)
 ABI_FREE(ks_vhartr)
 ABI_FREE(ks_vtrial)
 ABI_FREE(vpsp)
 ABI_FREE(ks_vxc)
 ! PAW
 ABI_SFREE(paw_pwff)
 ABI_SFREE(pawfgrtab)
 ABI_SFREE(ks_paw_an)

 call cryst%free(); call wfk_hdr%free(); call ebands_free(ks_ebands); call destroy_mpi_enreg(mpi_enreg_seq)
 call gwr%free()

end subroutine gwr_driver
!!***

!!****f* m_gwr_driver/cc4s_write_eigens
!! NAME
!!  cc4s_write_eigens
!!
!! FUNCTION
!!  Write eigenvalues in CC4S format. Only master proc should call this routine.
!!
!! INPUTS

subroutine cc4s_write_eigens(ebands, dtfil)

!Arguments ------------------------------------
 type(ebands_t),intent(in) :: ebands
 type(datafiles_type),intent(in) :: dtfil

!Local variables-------------------------------
 integer :: unt, ik_ibz, spin, band
 character(len=500) :: msg
 character(len=fnlen) :: filepath
! *************************************************************************

 ! See https://manuals.cc4s.org/user-manual/objects/EigenEnergies.html
 filepath = trim(dtfil%filnam_ds(4))//'_EigenEnergies.yaml'
 write(ab_out, "(3a)")ch10," Writing Eigenenergies metadata to file: ", trim(filepath)
 if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
   ABI_ERROR(msg)
 end if
 write(unt,'(A)')    'version: 100'
 write(unt,'(A)')    'type: Tensor'
 write(unt,'(A)')    'scalarType: Real64'
 write(unt,'(A)')    'dimensions:'
 write(unt,'(A,I0)') '- length: ',ebands%mband * ebands%nkpt * ebands%nsppol
 write(unt,'(A)')    '  type: State'
 write(unt,'(A)')    'elements:'
 write(unt,'(A)')    '  type: TextFile'
 write(unt,'(A)')    'unit: 1.0     # Hartree units'
 write(unt,'(A)')    'metaData:'
 write(unt,'(A,E22.15)')    '  fermiEnergy: ',ebands%fermie
 write(unt,'(A)')    '  energies:'

 do spin=1,ebands%nsppol
   do ik_ibz=1,ebands%nkpt
     do band=1,ebands%nband(ik_ibz + (spin-1)*ebands%nkpt)
       write(unt,"(a,e22.15)") '  - ',ebands%eig(band,ik_ibz,spin)
     end do
   end do
 end do
 close(unt)

 filepath = trim(dtfil%filnam_ds(4))//'_EigenEnergies.elements'
 write(ab_out, "(3a)")ch10," Writing Eigenenergies to file: ", trim(filepath)
 if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
   ABI_ERROR(msg)
 end if

 do spin=1,ebands%nsppol
   do ik_ibz=1,ebands%nkpt
     do band=1,ebands%nband(ik_ibz + (spin-1)*ebands%nkpt)
       write(unt,"(e22.15)") ebands%eig(band,ik_ibz,spin)
     end do
   end do
 end do

 close(unt)

end subroutine cc4s_write_eigens
!!***

!!****f* m_gwr_driver/cc4s_gamma
!! NAME
!!  cc4s_gamma
!!
!! FUNCTION
!! Interface with CC4S code.
!! Compute <b1,k|e^{-iGr}|b2,k> matrix elements and store them to disk
!!
!! INPUTS

subroutine cc4s_gamma(spin, ik_ibz, dtset, dtfil, cryst, ebands, psps, pawtab, paw_pwff, ugb)

 use m_numeric_tools, only : blocked_loop
 use m_gwdefs,        only : GW_Q0_DEFAULT
 use m_fftcore,       only : sphereboundary !, getng
 use m_fft_mesh,      only : setmesh
 use m_fft,           only : uplan_t ! fft_ug, fft_ur,
 use m_vcoul,         only : vcgen_t
 use m_pstat,         only : pstat_t
 use m_pawpwij,       only : pawpwij_t, pawpwij_init, pawpwij_free

!Arguments ------------------------------------
 integer,intent(in) :: spin, ik_ibz
 type(dataset_type),intent(in) :: dtset
 type(datafiles_type),intent(in) :: dtfil
 type(crystal_t),intent(in) :: cryst
 type(ebands_t),intent(in) :: ebands
 type(Pseudopotential_type),intent(in) :: psps
 type(pawpwff_t),intent(in) :: paw_pwff(dtset%ntypat*dtset%usepaw)
 type(pawtab_type),intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(ugb_t),target,intent(in) :: ugb

!Local variables-------------------------------
!scalars
 integer,parameter :: mG0(3) = 0, master = 0
 integer :: nproc, my_rank, my_ib2st, npw_k, nspinor, m_npw, npwvec, ig, mpierr, fh, comm, buf_size, ierr
 integer :: band1, band1_start, batch1_size, n1dat, idat1, m_istwfk, iatom, dim_rtwg
 integer :: band2, band2_start, batch2_size, n2dat, idat2, units(2), ii, unt, nqibz_, nqbz_, nkbz_, test_unt, M_
 integer(XMPI_OFFSET_KIND) :: offset
 real(dp) :: cpu, wall, gflops, qpt(3), qbz_(3,1), gcart(3), kpt(3), max_abs_err, abs_err, my_gw_qlwl(3), mem_mb, bz_vol
 character(len=500) :: msg
 character(len=fnlen) :: filepath, cvx_filepath
 logical :: k_is_gamma,  debug_this
 type(uplan_t) :: uplan_1, uplan_2, uplan_m
 type(vcgen_t) :: vcgen
 !type(pstat_t) :: pstat
 integer :: u_ngfft(18), u_nfft, u_mgfft, enforce_sym, method, nlmn_atm(cryst%natom)
 integer,pointer :: gvec_max(:,:)
 integer,allocatable,target :: m_gvec(:,:)
 complex(dp),allocatable :: ug1_batch(:,:), ur1_batch(:,:), ur2_batch(:,:), ur12_batch(:,:), ug12_batch(:,:), work(:)
 complex(gwpc),allocatable :: sqrt_vc(:), paw_rhotwg(:)
 type(pawpwij_t),allocatable :: pwij(:)
 type(pawcprj_type),allocatable :: cprj1(:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 comm = ugb%comm; nproc = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 units = [std_out, ab_out]
 npw_k = ugb%npw_k; nspinor = ugb%nspinor
 bz_vol = two_pi**3/cryst%ucvol

 debug_this = merge(.False., .True., nproc > 1)
 debug_this = .False.

 if (dtset%prtvol > 10) call ugb%print([std_out], dtset%prtvol, header="ugb for CC4S")

 ! m_gvec is the g-sphere for the oscillators M computed from ecuteps (half-sphere if wavefunctions have TR).
 kpt = dtset%kptns(:,ik_ibz); k_is_gamma = all(abs(kpt) < tol12)
 m_istwfk = 1; if (ugb%istwf_k == 2) m_istwfk = 2
 call get_kg(kpt, m_istwfk, dtset%ecuteps, cryst%gmet, m_npw, m_gvec, kin_sorted=.False.)

 ! Setup FFT mesh
 u_ngfft = dtset%ngfft
 method = 2
 if (dtset%fftgw==00 .or. dtset%fftgw==01) method=0
 if (dtset%fftgw==10 .or. dtset%fftgw==11) method=1
 if (dtset%fftgw==20 .or. dtset%fftgw==21) method=2
 if (dtset%fftgw==30 .or. dtset%fftgw==31) method=3
 enforce_sym = mod(dtset%fftgw, 10)
 ! Gamma only --> we don't need to rotate wavefunctions in the BZ
 if (k_is_gamma) enforce_sym = 0

 npwvec = npw_k; gvec_max => ugb%kg_k
 if (m_npw > npw_k) then
   npwvec = m_npw; gvec_max => m_gvec
 end if
 call setmesh(cryst%gmet, gvec_max, u_ngfft, npwvec, m_npw, npw_k, u_nfft, method, mG0, cryst, enforce_sym, unit=std_out)
 u_mgfft = maxval(u_ngfft(1:3))
 qpt = zero

 ! Get full BZ associated to ebands
 !call kpts_ibz_from_kptrlatt(cryst, ebands%kptrlatt, ebands%kptopt, ebands%nshiftk, ebands%shiftk, &
 !  nkibz, kibz, wtk, nkbz, kbz, bz2ibz=bz2ibz)

 nqibz_ = 1; nqbz_ = 1; qbz_ = zero; nkbz_ = 1
 ! TODO: MC technique does not seem to work as expected, even in the legacy code.
 call vcgen%init(cryst, ebands%kptrlatt, nkbz_, nqibz_, nqbz_, qbz_, &
                 dtset%rcut, dtset%gw_icutcoul, dtset%vcutgeo, dtset%ecuteps, comm)

 ! NB: npweps = m_npw
 ABI_MALLOC(sqrt_vc, (m_npw))
 my_gw_qlwl(:) = GW_Q0_DEFAULT; if (dtset%gw_nqlwl > 0) my_gw_qlwl = dtset%gw_qlwl(:,1)
 call vcgen%get_vc_sqrt(qpt, m_npw, m_gvec, my_gw_qlwl, cryst, sqrt_vc, comm)
 sqrt_vc(1) = zero
 call vcgen%free()

 if (my_rank == master) then
   call wrtout(units, " Computing oscilator matrix elements for CC4S.")
   do ii=1,size(units)
     call print_ngfft(u_ngfft, header='FFT mesh for wavefunctions', unit=units(ii))
   end do

   ! =====================
   ! Write files for CC4S
   ! =====================
   if (ik_ibz == 1 .and. spin == 1) then
     ! Write g-vector files, see https://manuals.cc4s.org/user-manual/objects/GridVectors.html
     filepath = trim(dtfil%filnam_ds(4))//'_GridVectors.yaml'
     write(ab_out, "(3a)")ch10," Writing Gridvectors metadata to file: ", trim(filepath)
     if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     write(unt,'(a)')    'version: 100'
     write(unt,'(a)')    'type: Tensor'
     write(unt,'(a)')    'scalarType: Real64'
     write(unt,'(a)')    'dimensions:'
     write(unt,'(a)')    '  - length: 3'
     write(unt,'(a)')    '    type: Vector'
     write(unt,'(a,i0)') '  - length: ',m_npw
     write(unt,'(a)')    '    type: Momentum'
     write(unt,'(a)')    'elements:'
     write(unt,'(a)')    '  type: TextFile'
     write(unt,'(a)')    'unit: 1.0  # Bohr^-1'
     ! The last three lines correspond to the reciprocal lattice vectors (including the factor 2pi)
     write(unt,'(a)')    'metaData:'
     write(unt,'("  Gi: [",E22.15,",",E22.15,",",E22.15,"]")') &
        two_pi*cryst%gprimd(1,1), two_pi*cryst%gprimd(2,1), two_pi*cryst%gprimd(3,1)
     write(unt,'("  Gj: [",E22.15,",",E22.15,",",E22.15,"]")') &
        two_pi*cryst%gprimd(1,2), two_pi*cryst%gprimd(2,2), two_pi*cryst%gprimd(3,2)
     write(unt,'("  Gk: [",E22.15,",",E22.15,",",E22.15,"]")') &
       two_pi*cryst%gprimd(1,3), two_pi*cryst%gprimd(2,3), two_pi*cryst%gprimd(3,3)
     close(unt)

     ! Write g-vectors
     filepath = trim(dtfil%filnam_ds(4))//'_GridVectors.elements'
     if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
       ABI_ERROR(msg)
     end if

     !TODO: k-points to be implemented.
     ! MG-TODO: I believe these are cart coords
     do ig=1,m_npw
       gcart = two_pi * matmul(cryst%gprimd, m_gvec(:,ig))
       do ii=1,3
         !write(unt,*) two_pi*GVEC_FULL(ii,ig,KQ)
         write(unt, *) gcart(ii)
       end do
     end do
     close(unt)

     ! https://manuals.cc4s.org/user-manual/objects/CoulombVertex.html
     filepath = trim(dtfil%filnam_ds(4))//'_CoulombVertex.yaml'
     write(ab_out, "(3a)")ch10, ' Writing CoulombVertex metadata to file: ', trim(filepath)
     if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     write(unt,'(a)')    'version: 100'
     write(unt,'(a)')    'type: Tensor'
     write(unt,'(a)')    'scalarType: Complex64'
     write(unt,'(a)')    'dimensions:'
     write(unt,'(a,i0)') '- length: ',m_npw
     write(unt,'(a)')    '  type: AuxiliaryField'
     !write(unt,'(a,i0)') '- length: ',(NBANDSDUMP)*WDES%ISPIN
     write(unt,'(a,i0)') '- length: ',ugb%nband_k
     write(unt,'(a)')    '  type: State'
     write(unt,'(a,i0)') '- length: ',ugb%nband_k
     write(unt,'(a)')    '  type: State'
     write(unt,'(a)')    'elements:'
     write(unt,'(a)')    '  type: IeeeBinaryFile'
     write(unt,'(a)')    'unit: 1.0   # Atomic units'
     !write(unt,'(a)')    'unit: 0.1917011272153577       # = sqrt(Eh/eV)'
     write(unt,'(a)')    'metaData:'
     if (m_istwfk == 2) then
       write(unt,'(a)')    '  halfGrid: 1'
     else
       write(unt,'(a)')    '  halfGrid: 0'
     end if
     close(unt)

     ! https://manuals.cc4s.org/user-manual/objects/CoulombPotential.html
     filepath = trim(dtfil%filnam_ds(4))//'_CoulombPotential.yaml'
     write(ab_out, "(3a)")ch10, ' Writing CoulombPotential metadata to file: ', trim(filepath)
     if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     write(unt,'(a)')    'version: 100'
     write(unt,'(a)')    'type: Tensor'
     write(unt,'(a)')    'scalarType: Real64'
     write(unt,'(a)')    'dimensions:'
     write(unt,'(a,i0)') '  - length: ',m_npw
     write(unt,'(a)')    '    type: Momentum'
     write(unt,'(a)')    'elements:'
     write(unt,'(a)')    '  type: TextFile'
     write(unt,'(a)')    'unit: 1.0     # Atomic units '
     !rite(unt,'(a)')    'unit: 0.2479966649373453       # =(Eh/eV*Bohr^3/Angstrom^3)'
     close(unt)

     filepath = trim(dtfil%filnam_ds(4))//'_CoulombPotential.elements'
     write(ab_out, "(3a)")ch10, ' Writing CoulombPotential data to file: ', trim(filepath)
     if (open_file(filepath, msg, newunit=unt, access="stream", form="formatted", status="replace", action="write") /= 0) then
       ABI_ERROR(msg)
     end if
     do ig=1,m_npw
       write(unt,*)real(sqrt_vc(ig) * conjg(sqrt_vc(ig)), kind=dp)
       !write(unt,*)real(sqrt_vc(ig)**2)
     end do
     !KQ=1 !k-points to be implemented
     !DO NG=1,NGVECTOR
     !  write(unt,*) REAL(POTFAK_FULL(NG,KQ)*CONJG(POTFAK_FULL(NG,KQ)),kind=q)
     !ENDDO
     close(unt)
   end if ! ik_ibz == 1 .and. spin == 1
 end if ! my_rank == master

 ! Open binary file to store CoulombVertex
 cvx_filepath = trim(dtfil%filnam_ds(4))//'_CoulombVertex.elements'
 if (my_rank == master) write(ab_out, "(3a)")ch10, ' Writing CoulombVertex data to file: ', trim(cvx_filepath)
#ifdef HAVE_MPI_IO
 call MPI_FILE_OPEN(comm, cvx_filepath, MPI_MODE_CREATE + MPI_MODE_WRONLY, xmpio_info, fh, mpierr)
 ABI_CHECK_MPI(mpierr, "MPI_FILE_OPEN")
#else
 ABI_ERROR("CC4S interface requires MPI-IO!")
#endif

 if (debug_this .and. my_rank == 0) then
   if (open_file("test_mg", msg, newunit=test_unt, form="formatted", status="replace", action="write") /= 0) then
     ABI_ERROR(msg)
   end if
 end if

 ! Define batch sizes and allocate workspace arrays.
 ! Increasing this value improves efficiency (less communication) at the price of more memory.
 !call pstat%from_pid(); call pstat%print([std_out], reload=.True.)

 batch1_size = min(48, ugb%nband_k); batch2_size = min(48, ugb%nband_k)
 !batch1_size = 1; batch2_size = 1
 call wrtout(std_out, sjoin(" Using batch1_size:", itoa(batch1_size), ", batch2_size:",  itoa(batch2_size)))

 mem_mb = (two * u_nfft * nspinor * batch1_size + &
           two * u_nfft * nspinor * batch2_size * two + &
           two * m_npw * nspinor * batch2_size) *  dp * b2Mb
 call wrtout(std_out, sjoin(" Memory for workspace arrays ", ftoa(mem_mb, fmt="f8.1"), "[Mb] <<< MEM"))

 ABI_MALLOC(ur1_batch, (u_nfft * nspinor, batch1_size))
 ABI_MALLOC(ur2_batch, (u_nfft * nspinor, batch2_size))
 ABI_MALLOC(ur12_batch, (u_nfft * nspinor, batch2_size))
 ABI_MALLOC(ug12_batch, (m_npw * nspinor, batch2_size))

 if (psps%usepaw == 1) then
   ! Evaluate oscillator matrix elements btw partial waves. Note q=Gamma
   ABI_MALLOC(pwij, (psps%ntypat))
   call pawpwij_init(pwij, m_npw, qpt, m_gvec, cryst%rprimd, psps, pawtab, paw_pwff)
   do iatom=1,cryst%natom
     nlmn_atm(iatom) = pawtab(cryst%typat(iatom))%lmn_size
   end do
   ABI_MALLOC(cprj1,  (cryst%natom, nspinor*batch1_size))
   call pawcprj_alloc(cprj1, 0, nlmn_atm)
   dim_rtwg = nspinor
   ABI_MALLOC(paw_rhotwg, (m_npw*dim_rtwg))
 end if

 ! TODO:
 ! 1) take advantage of M_{b1,b2}(g) = <b1|e^{-ig.r}|b2> => M_{b1,b2}(g) = M_{b2,b1}(-g)^*
 !    once I have a better understanding of the fileformat expected by CC4S.
 ! 2) Handle parallel IO if nsppol 2 (we are inside the spin loop that is already MPI distributed!)
 ! 3) Clarify ordering of CoulombVertex (b1,b2 vs b2,b1) and eigenvalues (spin?)
 ! 4) Treatment of q--> 0 in vc_coul
 ! 4) See other TODOs below.

 call uplan_1%init(npw_k, nspinor, batch1_size, u_ngfft, ugb%istwf_k, ugb%kg_k, dp, dtset%use_gpu_cuda)
 call uplan_2%init(npw_k, nspinor, batch2_size, u_ngfft, ugb%istwf_k, ugb%kg_k, dp, dtset%use_gpu_cuda)
 call uplan_m%init(m_npw, nspinor, batch2_size, u_ngfft, m_istwfk, m_gvec, dp, dtset%use_gpu_cuda)

 M_ = m_npw

 ! Blocked loop over group of b1 indices. NB: Assuming bands distributed in contiguous blocks.
 do band1_start=1, ugb%nband_k, batch1_size
   ! Collect n1dat bands starting from band1_start on each proc.
   n1dat = blocked_loop(band1_start, ugb%nband_k, batch1_size)

   call ugb%mat%collect_cplx(ugb%npwsp, n1dat, [1, band1_start], ug1_batch)
   if (psps%usepaw == 1) call ugb%collect_cprj(nspinor, n1dat, band1_start, cprj1)

   ! FFT: ug1_batch --> ur1_batch
   !call fft_ug(ugb%npw_k, u_nfft, nspinor, n1dat, u_mgfft, u_ngfft, ugb%istwf_k, ugb%kg_k, gbound_k, ug1_batch, ur1_batch)
   call uplan_1%execute_gr(n1dat, ug1_batch(:,1), ur1_batch(:,1))
   if (ugb%istwf_k /= 2) ur1_batch = conjg(ur1_batch)  ! Not needed if k == Gamma as ur1 is real.
   ABI_FREE(ug1_batch)

   ! Blocked loop over MY group of b2 indices (contiguous blocks)
   do band2_start=ugb%my_bstart, ugb%my_bstop, batch2_size
     n2dat = blocked_loop(band2_start, ugb%my_bstop, batch2_size)
     my_ib2st = band2_start - ugb%my_bstart + 1

     ! FFT: ugb%mat --> ur2_batch for n2dat states.
     !call fft_ug(ugb%npw_k, u_nfft, nspinor, n2dat, u_mgfft, u_ngfft, ugb%istwf_k, ugb%kg_k, gbound_k, &
     !            ugb%mat%buffer_cplx(:,my_ib2st), ur2_batch(:,1))

     call uplan_2%execute_gr(n2dat, ugb%mat%buffer_cplx(:,my_ib2st), ur2_batch(:,1))

     ! For each row of the submatrix, build n2dat products (band1, idat2) in r-space, then r --> g.
     do idat1=1,n1dat
       band1 = band1_start + idat1 - 1

       do idat2=1,n2dat
         ur12_batch(:,idat2) = ur1_batch(:,idat1) * ur2_batch(:,idat2)
       end do
       call uplan_m%execute_rg(n2dat, ur12_batch(:,1), ug12_batch(:,1))

       !call fft_ur(m_npw, u_nfft, nspinor, n2dat, u_mgfft, u_ngfft, m_istwfk, m_gvec, m_gbound, ur12_batch, ug12_batch)

       if (psps%usepaw == 1) then
         ! Add PAW on-site contributions
         do idat2=1,n2dat
           associate (cprj1_kmq => cprj1(:, 1 + (idat1-1)*nspinor), &
                      cprj2_k => ugb%cprj_k(:, 1 + (my_ib2st+idat2-2)*nspinor))  ! NB: ugb%cprj_k(2, nspinor*my_nband)
           paw_rhotwg = zero
           call paw_rho_tw_g(cryst, pwij, m_npw, dim_rtwg, nspinor, m_gvec, cprj1_kmq, cprj2_k, paw_rhotwg)
           ug12_batch(:,idat2) = ug12_batch(:,idat2) + paw_rhotwg
           end associate
         end do
       end if

       if (nspinor == 2) then
         ! Sum over spinors and repack data in the first n2dat positions to prepare IO operation.
         do idat2=1,n2dat
           ug12_batch(1:m_npw,idat2) = ug12_batch(1:m_npw,idat2) + ug12_batch(m_npw+1:,idat2)
         end do
         do idat2=2,n2dat,2
           ug12_batch(m_npw+1:,idat2-1) = ug12_batch(1:m_npw,idat2)
         end do
       end if

       do idat2=1,n2dat
         !if (band1 == band2_start + idat2 -1)  then
           write(std_out,*) " ug12_batch(g=0,band1,band2), band", ug12_batch(1,idat2), band1, band2_start + idat2 -1
         !end if
         ! Multiply by sqrt(vc(g))
         ug12_batch(:,idat2) = ug12_batch(:,idat2) * sqrt_vc(:) * sqrt(cryst%ucvol) ! * sqrt(bz_vol) ! * FIXME
       end do
       !write(std_out,*)" max(abs(ug12_batch)):", maxval(abs(ug12_batch(:,1:n2dat)))

#ifdef HAVE_MPI_IO
       ! Write ug12_batch using Stream-IO
       buf_size = m_npw * n2dat
       offset = ((band2_start-1) * m_npw + (band1-1) * m_npw * ugb%nband_k) * xmpi_bsize_dpc
       call MPI_FILE_WRITE_AT(fh, offset, ug12_batch, buf_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
       ABI_HANDLE_MPIERR(mpierr)

       if (my_rank == 0 .and. debug_this) then
         do idat2=1,n2dat
           band2 = band2_start + idat2 - 1
           write(test_unt,*)band1, band2, ug12_batch(1:M_,idat2)
         end do
      end if
#endif
     end do ! idat1

   end do ! band2_start
 end do ! band1_start

#ifdef HAVE_MPI_IO
 call MPI_FILE_CLOSE(fh, mpierr)
 ABI_CHECK_MPI(mpierr, "FILE_CLOSE!")
 call xmpi_barrier(comm)

 if (my_rank == 0) then
   buf_size = 4
   call wrtout(units, sjoin(" Writing Coulomb vertex for testing purposes with ng:", itoa(buf_size)), newlines=1, pre_newlines=1)
   ABI_MALLOC(work, (buf_size))
   call MPI_FILE_OPEN(xmpi_comm_self, cvx_filepath, MPI_MODE_RDONLY, xmpio_info, fh, mpierr)
   ABI_CHECK_MPI(mpierr, "MPI_FILE_OPEN")
   ierr = 0
   band1_loop: do band1=1, ugb%nband_k
   do band2=1, ugb%nband_k
     ierr = ierr + 1; if (ierr == 6) exit band1_loop
     !if (ig == 1) work(1) = zero
     offset = ((band2-1) * m_npw + (band1-1) * m_npw * ugb%nband_k) * xmpi_bsize_dpc
     call MPI_FILE_READ_AT(fh, offset, work, buf_size, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_HANDLE_MPIERR(mpierr)
     call wrtout(units, sjoin(" For band1:", itoa(band1), ", band2:", itoa(band2)))
     write(msg, "(*(1x, es12.5))")work(1:buf_size)
     call wrtout(units, msg)
   end do
   end do band1_loop
   ABI_FREE(work)
 end if
#endif

 ! =============
 ! DEBUG SECTION
 ! =============
#ifdef HAVE_MPI_IO
 if (my_rank == 0 .and. debug_this) then
   close(test_unt)
   if (open_file("test_mg", msg, newunit=test_unt, form="formatted", status="old", action="read") /= 0) then
     ABI_ERROR(msg)
   end if

   call MPI_FILE_OPEN(xmpi_comm_self, cvx_filepath, MPI_MODE_RDONLY, xmpio_info, fh, mpierr)
   ABI_CHECK_MPI(mpierr, "MPI_FILE_OPEN")

   ABI_MALLOC(work, (m_npw))
   max_abs_err = zero
   do ig=1, ugb%nband_k**2
     read(test_unt,*) band1, band2, ug12_batch(1:M_,1)
     offset = ((band2-1) * m_npw + (band1-1) * m_npw * ugb%nband_k) * xmpi_bsize_dpc
     call MPI_FILE_READ_AT(fh, offset, work, m_npw, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, mpierr)
     ABI_HANDLE_MPIERR(mpierr)

     abs_err = maxval(abs(ug12_batch(1:M_,1) - work(1:M_)))
     max_abs_err = max(max_abs_err, abs_err)
     if (abs_err > zero) write(std_out, *)" For ig:", ig, "/", ugb%nband_k**2, "abs_err", abs_err
     !write(std_out, *)"1:", ug12_batch(1:M_,1); write(std_out, *)"2:", work(1:M_)
   end do

   close(test_unt)
   call MPI_FILE_CLOSE(fh, mpierr)
   ABI_CHECK_MPI(mpierr, "FILE_CLOSE!")
   ABI_FREE(work)

   write(std_out,*)" max_abs_err:", max_abs_err
   ABI_CHECK(max_abs_err < tol16, sjoin("max_abs_err:", ftoa(max_abs_err)))
   call wrtout(std_out, " Debugging section OK!!!")
 end if
#endif

 call uplan_1%free(); call uplan_2%free(); call uplan_m%free()

 ABI_FREE(m_gvec)
 ABI_FREE(ur1_batch)
 ABI_FREE(ur2_batch)
 ABI_FREE(ur12_batch)
 ABI_FREE(ug12_batch)
 ABI_FREE(sqrt_vc)

 if (psps%usepaw == 1) then
   call pawpwij_free(pwij)
   ABI_FREE(pwij)
   call pawcprj_free(cprj1)
   ABI_FREE(cprj1)
   ABI_FREE(paw_rhotwg)
 end if

 call cwtime_report(" cc4s_gamma", cpu, wall, gflops)

end subroutine cc4s_gamma
!!***

end module m_gwr_driver
!!***

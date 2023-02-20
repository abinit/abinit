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
 use iso_c_binding

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : MPI_type
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
 use m_pawpwij,         only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_kg,              only : getph !, getcut
 use m_pspini,          only : pspini
 use m_paw_correlations,only : pawpuxinit
 use m_paw_dmft,        only : paw_dmft_type
 use m_paw_sphharm,     only : setsym_ylm
 use m_paw_mkrho,       only : denfgr
 use m_paw_nhat,        only : nhatgrid, pawmknhat
 use m_paw_tools,       only : chkpawovlp, pawprt
 use m_paw_denpot,      only : pawdenpot
 use m_paw_init,        only : pawinit, paw_gencond
 use m_pawcprj,         only : pawcprj_type, pawcprj_free !, pawcprj_alloc, paw_overlap
 use m_ksdiago,         only : ugb_t
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
!! rprim(3,3)=Dimensionless real space primitive translations.
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

subroutine gwr_driver(acell, codvsn, dtfil, dtset, pawang, pawrad, pawtab, psps, rprim, xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
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
 logical :: read_wfk, write_wfk, cc4s_task
 character(len=500) :: msg
 character(len=fnlen) :: wfk_path, den_path, kden_path, out_path
 type(hdr_type) :: wfk_hdr, den_hdr, kden_hdr, owfk_hdr
 type(crystal_t) :: cryst, den_cryst, wfk_cryst
 type(ebands_t) :: ks_ebands, owfk_ebands
 type(pawfgr_type) :: pawfgr
 type(wvl_data) :: wvl
 type(mpi_type) :: mpi_enreg
 type(gwr_t) :: gwr
 type(wfk_t) :: owfk
!arrays
 real(dp), parameter :: k0(3) = zero
 integer :: cplex, cplex_dij, cplex_rhoij
 integer :: gnt_option !,has_dijU,has_dijso,iab,bmin,bmax,irr_idx1,irr_idx2
 integer :: istep, moved_atm_inside, moved_rhor, n3xccc, sc_mode
 integer :: ndij !,ndim,nfftf,nfftf_tot,nkcalc,gwc_nfft,gwc_nfftot,gwx_nfft,gwx_nfftot
 integer :: ngrvdw, nhatgrdim, nkxc, nspden_rhoij, optene !nzlmopt,
 integer :: optcut, optgr0, optgr1, optgr2, optrad, psp_gencond !option,
 integer :: rhoxsp_method, usexcnhat !, use_aerhor,use_umklp
 !real(dp) :: compch_fft, compch_sph !,r_s,rhoav,alpha
 real(dp) :: gsqcutc_eff, gsqcutf_eff, gsqcut_shp
 real(dp) :: vxcavg !,vxcavg_qp ucvol,
 real(dp) :: gwc_gsq, gwx_gsq,gw_gsq !, gsqcut
 logical :: call_pawinit
 type(energies_type) :: KS_energies
 type(melflags_t) :: KS_mflags
 type(paw_dmft_type) :: Paw_dmft
 type(ugb_t), target :: ugb
 type(xmpi_pool2d_t) :: diago_pool
!arrays
 integer :: ngfftc(18),ngfftf(18),unts(2) ! gwc_ngfft(18),gwx_ngfft(18),
 integer,allocatable :: nq_spl(:), l_size_atm(:) !,qp_vbik(:,:),tmp_gfft(:,:)
 integer,allocatable :: tmp_kstab(:,:,:), npwarr_ik(:), gvec_(:,:), istwfk_ik(:), nband_iks(:,:)
 real(dp) :: strsxc(6), diago_info(3, dtset%nkpt, dtset%nsppol) !,tsec(2)
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

 call cwtime(cpu, wall, gflops, "start")

! write(msg,'(7a)')&
! ' SIGMA: Calculation of the GW corrections ',ch10,ch10,&
! ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
! ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
! call wrtout([std_out, ab_out], msg)
!
#if defined HAVE_GW_DPC
 write(msg,'(a,i2,a)')'.Using double precision arithmetic; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic; gwpc = ',gwpc,ch10
#endif
 call wrtout([std_out, ab_out], msg)

 if (psps%usepaw == 1) then
   ABI_ERROR("PAW not implemented")
   ABI_UNUSED((/pawang%nsym, pawrad(1)%mesh_size/))
 end if

#ifndef HAVE_LINALG_SCALAPACK
  ABI_ERROR("GWR code requires scalapack library")
#endif

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)
 unts(:) = [std_out, ab_out]

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

 den_path = dtfil%fildensin
 wfk_path = dtfil%fnamewffk
 kden_path = dtfil%filkdensin

 if (my_rank == master) then
   ! Initialize filenames. Accept files in Fortran or in netcdf format.
   ! Accept DEN file in Fortran or in netcdf format.
   if (nctk_try_fort_or_ncfile(den_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find DEN file:", den_path, ". Error:", msg))
   end if
   call wrtout(unts, sjoin("- Reading GS density from: ", den_path))

   if (dtset%usekden == 1) then
     if (nctk_try_fort_or_ncfile(kden_path, msg) /= 0) then
       ABI_ERROR(sjoin("Cannot find KDEN file:", kden_path, ". Error:", msg))
     end if
     call wrtout(unts, sjoin("- Reading KDEN kinetic energy density from: ", kden_path))
   end if
   call wrtout(ab_out, ch10//ch10)
 end if ! master

 ! Broadcast filenames (needed because they might have been changed if we are using netcdf files)
 call xmpi_bcast(den_path, master, comm, ierr)
 call xmpi_bcast(kden_path, master, comm, ierr)

 ! TODO: FFT meshes for DEN/POT should be initialized from the DEN file instead of the dtset.
 ! Interpolating the DEN indeed breaks degeneracies in the vxc matrix elements.

 call pawfgr_init(pawfgr, dtset, mgfftf, nfftf, ecut_eff, ecutdg_eff, ngfftc, ngfftf, &
                  gsqcutc_eff=gsqcutc_eff, gsqcutf_eff=gsqcutf_eff, gmet=cryst%gmet, k0=k0)

 call print_ngfft(ngfftc, header='Coarse FFT mesh for the wavefunctions')
 call print_ngfft(ngfftf, header='Dense FFT mesh for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft, 'c', ngfftc(2), ngfftc(3), 'all')
 call init_distribfft_seq(mpi_enreg%distribfft, 'f', ngfftf(2), ngfftf(3), 'all')

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

   if (psp_gencond == 1 .or. call_pawinit) then
     gsqcut_shp = two * abs(dtset%diecut) * dtset%dilatmx**2 / pi**2
     call pawinit(dtset%effmass_free, gnt_option, gsqcut_shp, zero, dtset%pawlcutd, dtset%pawlmix, &
                  psps%mpsang, dtset%pawnphi, cryst%nsym, dtset%pawntheta, pawang, pawrad, &
                  dtset%pawspnorb, pawtab, dtset%pawxcdev, dtset%ixc, dtset%usepotzero)

     ! Update internal values
     call paw_gencond(dtset, gnt_option, "save", call_pawinit)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:cryst%ntypat)%has_kij   = 2
     if (pawtab(1)%has_nabla==1) pawtab(1:cryst%ntypat)%has_nabla = 2
   end if

   psps%n1xccc = MAXVAL(pawtab(1:cryst%ntypat)%usetcore)

   ! Initialize optional flags in Pawtab to zero
   ! Cannot be done in Pawinit since the routine is called only if some pars. are changed
   pawtab(:)%has_nabla = 0

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
   rhoxsp_method = 1  ! Arnaud-Alouani (default in sigma)
   !rhoxsp_method = 2 ! Shiskin-Kresse
   if (Dtset%pawoptosc /= 0) rhoxsp_method = Dtset%pawoptosc

   ! The q-grid must contain the FFT mesh used for sigma_c and the G-sphere for the exchange part.
   ! We use the FFT mesh for sigma_c since COHSEX and the extrapolar method require oscillator
   ! strengths on the FFT mesh.
   !ABI_MALLOC(tmp_gfft,(3, gwc_nfftot))
   !call get_gftt(gwc_ngfft, k0, gmet, gwc_gsq, tmp_gfft)
   !ABI_FREE(tmp_gfft)

   gwx_gsq = Dtset%ecutsigx/(two*pi**2)
   gw_gsq = MAX(gwx_gsq, gwc_gsq)

   ! Set up q-grid, make qmax 20% larger than largest expected.
   ABI_MALLOC(nq_spl, (Psps%ntypat))
   ABI_MALLOC(qmax, (Psps%ntypat))
   qmax = SQRT(gw_gsq)*1.2d0
   nq_spl = Psps%mqgrid_ff
   ! write(std_out,*)"using nq_spl",nq_spl,"qmax=",qmax
   ABI_MALLOC(Paw_pwff, (Psps%ntypat))

   call pawpwff_init(Paw_pwff,rhoxsp_method,nq_spl,qmax,cryst%gmet,Pawrad,Pawtab,Psps)

   ABI_FREE(nq_spl)
   ABI_FREE(qmax)

   ! Variables/arrays related to the fine FFT grid
   ABI_MALLOC(ks_nhat, (nfftf, Dtset%nspden))
   ks_nhat = zero
   ABI_MALLOC(Pawfgrtab, (Cryst%natom))
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)

   cplex = 1
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Dtset%nspden,Dtset%typat)
   ABI_FREE(l_size_atm)
   !compch_fft=greatest_real
   usexcnhat = MAXVAL(Pawtab(:)%usexcnhat)
   ! * 0 if Vloc in atomic data is Vbare    (Blochl's formulation)
   ! * 1 if Vloc in atomic data is VH(tnzc) (Kresse's formulation)
   call wrtout(std_out, sjoin(' using usexcnhat: ', itoa(usexcnhat)))
   !
   ! Identify parts of the rectangular grid where the density has to be calculated
   optcut = 0; optgr0 = Dtset%pawstgylm; optgr1 = 0; optgr2 = 0; optrad = 1 - Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(cryst%atindx1, cryst%gmet, cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
    optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   call pawfgrtab_print(Pawfgrtab,Cryst%natom,unit=std_out,prtvol=Dtset%pawprtvol)

 else
   ABI_MALLOC(Paw_pwff, (0))
   ABI_MALLOC(Pawfgrtab, (0))
 end if ! End of PAW Initialization

 ! Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(ks_nhat)) then
   ABI_MALLOC(ks_nhat, (nfftf, 0))
 end if
 if (.not.allocated(dijexc_core)) then
   ABI_MALLOC(dijexc_core, (1, 1, 0))
 end if

 ! Read density and compare crystal structures.
 ABI_MALLOC(ks_rhor, (nfftf, dtset%nspden))

 call read_rhor(den_path, cplex1, dtset%nspden, nfftf, ngfftf, dtset%usepaw, mpi_enreg, ks_rhor, &
                den_hdr, ks_pawrhoij, comm) !, allow_interp=.True.)

 ! TODO: Overloaded interface with units or just change the API to accept units
 call prtrhomxmn(std_out, MPI_enreg, nfftf, ngfftf, dtset%nspden, 1, ks_rhor, ucvol=cryst%ucvol)
 call prtrhomxmn(ab_out, MPI_enreg, nfftf, ngfftf, dtset%nspden, 1, ks_rhor, ucvol=cryst%ucvol)

 den_cryst = den_hdr%get_crystal()
 if (cryst%compare(den_cryst, header=" Comparing input crystal with DEN crystal") /= 0) then
   ABI_ERROR("Crystal structure from input and from DEN file do not agree! Check messages above!")
 end if
 call den_cryst%free()
 call den_hdr%free()

 ! FFT n(r) --> n(g)
 ABI_MALLOC(ks_rhog, (2, nfftf))
 call fourdp(cplex1, ks_rhog, ks_rhor(:, 1), -1, mpi_enreg, nfftf, 1, ngfftf, 0)

 ABI_MALLOC(ks_taur, (nfftf, dtset%nspden * dtset%usekden))
 if (dtset%usekden == 1) then
   call read_rhor(kden_path, cplex1, dtset%nspden, nfftf, ngfftf, 0, mpi_enreg, ks_taur, &
                  kden_hdr, ks_pawrhoij, comm) !, allow_interp=.True.)
   call kden_hdr%free()
   call prtrhomxmn(std_out, MPI_enreg, nfftf, ngfftf, dtset%nspden, 1, ks_taur, optrhor=1, ucvol=cryst%ucvol)
 end if

 !========================================
 !==== Additional computation for PAW ====
 !========================================
 nhatgrdim = 0
 if (dtset%usepaw == 1) then
   ABI_ERROR("PAW in GWR not yet implemented.")
 else
   ABI_MALLOC(ks_nhatgr, (0, 0, 0))
   ABI_MALLOC(ks_paw_ij, (0))
   ABI_MALLOC(ks_paw_an, (0))
 end if ! PAW

 ! Compute structure factor phases and large sphere cutoff
 ABI_MALLOC(ph1d, (2, 3 * (2 * Dtset%mgfft + 1) * Cryst%natom))
 ABI_MALLOC(ph1df, (2, 3 * (2 * mgfftf + 1) * Cryst%natom))

 call getph(Cryst%atindx, Cryst%natom, ngfftc(1), ngfftc(2), ngfftc(3), ph1d, Cryst%xred)

 if (psps%usepaw == 1 .and. pawfgr%usefinegrid == 1) then
   call getph(Cryst%atindx, Cryst%natom, ngfftf(1), ngfftf(2), ngfftf(3), ph1df, Cryst%xred)
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
 ABI_MALLOC(grchempottn, (3, Cryst%natom))
 ABI_MALLOC(grewtn, (3, Cryst%natom))
 nkxc = 0
 if (Dtset%nspden == 1) nkxc = 2
 if (Dtset%nspden >= 2) nkxc = 3 ! check GGA and spinor, quite a messy part!!!
 ! In case of MGGA, fxc and kxc are not available and we dont need them (for now ...)
 if (Dtset%ixc < 0 .and. libxc_functionals_ismgga()) nkxc = 0
 if (nkxc /= 0) then
   ABI_MALLOC(kxc, (nfftf, nkxc))
 end if

 n3xccc = 0; if (Psps%n1xccc /= 0) n3xccc = nfftf
 ABI_MALLOC(xccc3d, (n3xccc))
 ABI_MALLOC(ks_vhartr, (nfftf))
 ABI_MALLOC(ks_vtrial, (nfftf, Dtset%nspden))
 ABI_MALLOC(vpsp, (nfftf))
 ABI_MALLOC(ks_vxc, (nfftf, Dtset%nspden))

 ! I don't think direct diago can be used with mega-GGA due to the functional derivative wrt KS states.
 ! TB-BK should be OK though.

 !ABI_MALLOC(ks_vxctau, (nfftf, dtset%nspden * dtset%usekden))
 !ABI_MALLOC(xcctau3d, (n3xccc * dtset%usekden))
 !ABI_FREE(ks_vxctau)
 !ABI_FREE(xcctau3d)

 optene = 4; moved_atm_inside = 0; moved_rhor = 0; istep = 1

 call setvtr(Cryst%atindx1,Dtset,KS_energies,cryst%gmet,cryst%gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,&
             istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg,&
             Cryst%nattyp,nfftf,ngfftf,ngrvdw,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
             optene,pawrad,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,cryst%rmet,cryst%rprimd,strsxc,&
             Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,Wvl,xccc3d,Cryst%xred, &
             taur=ks_taur) !xcctau3d=xcctau3d, vxctau=ks_vxctau)

 ABI_FREE(grvdw)
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)

 call cwtime_report(" prepare gwr_driver_init", cpu, wall, gflops)

 if (string_in(dtset%gwr_task, "HDIAGO, HDIAGO_FULL, CC4S, CC4S_FULL")) then
   ! ==========================================
   ! Direct diagonalization of the Hamiltonian
   ! ==========================================
   ABI_MALLOC(nband_iks, (dtset%nkpt, dtset%nsppol))
   ABI_MALLOC(npwarr_ik, (dtset%nkpt))
   ABI_MALLOC(istwfk_ik, (dtset%nkpt))
   istwfk_ik = 1 ! TODO: istwkf 2 is not yet supported.

   ! Compute npw_k from ecut so that we can update the header.
   do ik_ibz=1,dtset%nkpt
     call get_kg(dtset%kptns(:,ik_ibz), istwfk_ik(ik_ibz), dtset%ecut, cryst%gmet, npwarr_ik(ik_ibz), gvec_)
     ABI_FREE(gvec_)
   end do

   ! CC4S does not need to output the WFK file.
   write_wfk = string_in(dtset%gwr_task, "HDIAGO, HDIAGO_FULL")

   ! Use input nband or min of npwarr_ik to set the number of bands.
   if (string_in(dtset%gwr_task, "HDIAGO, CC4S")) nband_iks(:,:) = maxval(dtset%nband)
   if (string_in(dtset%gwr_task, "HDIAGO_FULL, CC4S_FULL")) nband_iks(:,:) = minval(npwarr_ik)
   cc4s_task = string_in(dtset%gwr_task, "CC4S, CC4S_FULL")
   if (cc4s_task) then
     ABI_CHECK(dtset%nkpt == 1 .and. all(abs(dtset%kptns(:,1)) < tol12), "CC4S requires Gamm-only sampling")
   end if

   ! Build header with new npwarr and nband.
   owfk_ebands = ebands_from_dtset(dtset, npwarr_ik, nband=nband_iks)
   owfk_ebands%eig = zero
   call hdr_init(owfk_ebands, codvsn, dtset, owfk_hdr, pawtab, 0, psps, wvl%descr)

   ! Change the value of istwfk taken from dtset.
   ABI_REMALLOC(owfk_hdr%istwfk, (dtset%nkpt))
   owfk_hdr%istwfk(:) = istwfk_ik

   ! Build pools to distribute (kpt, spin)
   ! Try to have rectangular grids in each pool to improve efficiency in scalapack diago.
   call diago_pool%from_dims(dtset%nkpt, dtset%nsppol, comm, rectangular=.True.)
   diago_info = zero

   if (write_wfk) then
     ! ===============================================
     ! Master writes header and Fortran record markers
     ! ===============================================
     out_path = dtfil%fnameabo_wfk; if (dtset%iomode == IO_MODE_ETSF) out_path = nctk_ncify(out_path)
     iomode__ = iomode_from_fname(out_path)
     call wrtout(std_out, sjoin(" Writing wavefunctions to file:", out_path))
     if (my_rank == master) then
       call owfk%open_write(owfk_hdr, out_path, 0, iomode__, get_unit(), xmpi_comm_self, &
                            write_hdr=.True., write_frm=.True.)
       call owfk%close()
     end if
     call xmpi_barrier(comm)

     ! Reopen file inside diago_pool%comm.
     !call owfk%open_write(owfk_hdr, out_path, 0, iomode__, get_unit(), diago_pool%comm%value, &
     !                     write_hdr=.False., write_frm=.False.)
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
        ! occupancies are set to zero.
        ! Client code is responsible for recomputing occ and fermie when reading this WFK.
        !ABI_CALLOC(occ_k, (owfk%mband))
        ABI_CALLOC(occ_k, (nband_k))
        !sc_mode = merge(xmpio_single, xmpio_collective, ugb%has_idle_procs)
        !print *, "About to write with ugb%my_nband", ugb%my_nband, "has_idle_procs:", ugb%has_idle_procs
        !call xmpi_barrier(comm)

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

       if (cc4s_task) call cc4s_gamma(spin, ik_ibz, dtset, cryst, ugb)

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

   ! Collect eigenvalues
   ! TODO: CHECK THIS PART WITH MANY PROCS AS I'VE SEEN WEIRD RESULTS ON LUMI
   ! Also, Si4x4x4 with gamma-only seems not to work
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
   end if

   ABI_FREE(npwarr_ik)
   ABI_FREE(istwfk_ik)
   ABI_FREE(nband_iks)

   call owfk_hdr%free()
   call ebands_free(owfk_ebands)
   call diago_pool%free()
   !call owfk%close()

 !else if (dtset%gwr_task == "CC4CS") then
   ! Diagonalize Hamiltonian at k = Gamma
   ! Compute oscillator matrix elements and save results to disk
   !call cc4cs()

 else
   ! ====================================================
   ! === This is the real GWR stuff once all is ready ===
   ! ====================================================
   read_wfk = .True.
   if (read_wfk) then
     if (my_rank == master) then
       if (nctk_try_fort_or_ncfile(wfk_path, msg) /= 0) then
          ABI_ERROR(sjoin("Cannot find GS WFK file:", wfk_path, ". Error:", msg))
       end if
       call wrtout(unts, sjoin("- Reading GS states from WFK file:", wfk_path))
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
     call wfk_cryst%free()
     !call wfk_cryst%print(header="crystal structure from WFK file")

     ! Make sure that ef is inside the gap if semiconductor.
     call ebands_update_occ(ks_ebands, dtset%spinmagntarget, prtvol=dtset%prtvol, fermie_to_zero=.True.)

     ! Here we change the GS bands (Fermi level, scissors operator ...)
     ! All the modifications to ebands should be done here.
     !call ephtk_update_ebands(dtset, ks_ebands, "Ground state energies")
   end if

   call gwr%init(dtset, dtfil, cryst, psps, pawtab, ks_ebands, mpi_enreg, comm)
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
   if (Dtset%usepawu /= 0  )  KS_mflags%has_vu      = 1
   if (Dtset%useexexch /= 0)  KS_mflags%has_lexexch = 1
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

   select case (dtset%gwr_task)
   case ("RPA_ENERGY")
     call gwr%rpa_energy()
   case ("G0W0")
     call gwr%run_g0w0()
   case ("G0V")
     call gwr%build_sigxme()
   case ("EGEW", "EGW0", "G0EW")
     call gwr%run_energy_scf()
   !case ("CHI0")
   !  call compute_chi0(wfk_path)
   !case ("CHI0_HEAD_WINGS")
   !  call compute_chi0_head_wings(wfk_path)
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
   ABI_FREE(pawfgrtab)
   call paw_ij_free(ks_paw_ij)
   ABI_FREE(ks_paw_ij)
   !call paw_an_free(ks_paw_an)
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

 call cryst%free()
 call wfk_hdr%free()
 call ebands_free(ks_ebands)
 call destroy_mpi_enreg(mpi_enreg)
 call gwr%free()

end subroutine gwr_driver
!!***

!!****f* m_gwr_driver/cc4s_gamma
!! NAME
!!  cc4s_gamma
!!
!! FUNCTION
!! Interface with CC4S code.
!! Compute <i,k|e^{-iGr}|j,k> matrix elements and store them to disk
!!
!! INPUTS

subroutine cc4s_gamma(spin, ik_ibz, dtset, cryst, ugb)

 use m_numeric_tools, only : blocked_loop
 use m_fftcore,       only : sphereboundary !, getng
 use m_fft_mesh,      only : setmesh
 use m_fft,           only : fft_ug, fft_ur, uplan_t

!Arguments ------------------------------------
 integer,intent(in) :: spin, ik_ibz
 type(dataset_type),intent(in) :: dtset
 type(crystal_t),intent(in) :: cryst
 type(ugb_t), target,intent(in) :: ugb

!Local variables-------------------------------
!scalars
 integer,parameter :: istwfk1 = 1, mG0(3) = 0
 integer :: nproc, my_rank, my_ib2, npw_k, nspinor, m_npw, npwvec !master,
 integer :: band1_start, band1_stop, batch1_size, n1dat, idat1
 integer :: band2_start, band2_stop, batch2_size, n2dat, idat2
 real(dp) :: cpu, wall, gflops
 type(uplan_t) :: uplan_1, uplan_2, uplan_m
 integer :: u_ngfft(18), u_nfft, u_mgfft, enforce_sym, method
 integer,pointer :: gvec_max(:,:)
 integer,allocatable :: gbound_k(:,:), m_gbound(:,:)
 integer,allocatable,target :: m_gvec(:,:)
 complex(dpc),allocatable :: ug1_batch(:,:), ur1_batch(:,:), ur2_batch(:,:), ur12_batch(:,:), ug12_batch(:,:)

! *************************************************************************

 call cwtime(cpu, wall, gflops, "start")

 nproc = xmpi_comm_size(ugb%comm); my_rank = xmpi_comm_rank(ugb%comm)
 npw_k = ugb%npw_k; nspinor = ugb%nspinor

 ! g-sphere for oscillators.
 call get_kg([zero, zero, zero], istwfk1, dtset%ecuteps, cryst%gmet, m_npw, m_gvec, kin_sorted=.False.)

 ! Setup FFT mesh
 u_ngfft = dtset%ngfft
 method = 2
 if (dtset%fftgw==00 .or. dtset%fftgw==01) method=0
 if (dtset%fftgw==10 .or. dtset%fftgw==11) method=1
 if (dtset%fftgw==20 .or. dtset%fftgw==21) method=2
 if (dtset%fftgw==30 .or. dtset%fftgw==31) method=3
 enforce_sym = MOD(dtset%fftgw, 10)
 enforce_sym = 0 ! Gamma only --> we don't need to rotate wavefunctions in the BZ

 npwvec = npw_k; gvec_max => ugb%kg_k
 if (m_npw > npw_k) then
   npwvec = m_npw; gvec_max => m_gvec
 end if
 call setmesh(cryst%gmet, gvec_max, u_ngfft, npwvec, m_npw, npw_k, u_nfft, method, mG0, cryst, enforce_sym, unit=std_out)
 !call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 u_mgfft = maxval(u_ngfft(1:3))

 ABI_MALLOC(gbound_k, (2 * u_mgfft + 8, 2))
 call sphereboundary(gbound_k, ugb%istwf_k, ugb%kg_k, u_mgfft, npw_k)

 ABI_MALLOC(m_gbound, (2 * u_mgfft + 8, 2))
 call sphereboundary(m_gbound, istwfk1, m_gvec, u_mgfft, m_npw)

 ! Allocate workspace arrays.
 batch1_size = min(24, ugb%nband_k)
 batch2_size = 4
 !batch2_size = 1

 ABI_MALLOC(ur1_batch, (u_nfft * nspinor, batch1_size))
 ABI_MALLOC(ur2_batch, (u_nfft * nspinor, batch2_size))
 ABI_MALLOC(ur12_batch, (u_nfft * nspinor**2, batch2_size))
 ABI_MALLOC(ug12_batch, (m_npw * nspinor**2, batch2_size))
 ABI_CHECK(nspinor == 1, "nspinor == 2 not implemented in CC4S")

 ! TODO: take advantage of
 ! M_{12}(g) = <1|e^{-ig.r}|2> => M_{12}(g) = M_{21}(-g)^*
 ! once I have a better understanding of the fileformat expected by CC4S

 call uplan_1%init(npw_k, nspinor, batch1_size, u_ngfft, ugb%istwf_k, ugb%kg_k, dp, dtset%use_gpu_cuda)
 call uplan_2%init(npw_k, nspinor, batch2_size, u_ngfft, ugb%istwf_k, ugb%kg_k, dp, dtset%use_gpu_cuda)
 call uplan_m%init(m_npw, nspinor, batch2_size, u_ngfft, istwfk1, m_gvec, dp, dtset%use_gpu_cuda)

 do band1_start=1, ugb%nband_k, batch1_size

   ! Collect n1dat bands starting from band1_start on each proc.
   n1dat = blocked_loop(band1_start, ugb%nband_k, batch1_size)
   band1_stop = band1_start + n1dat - 1
   call ugb%mat%collect_cplx(ugb%npwsp, n1dat, [1, band1_start], ug1_batch)

   ! ug1_batch --> ur1_batch
   !call fft_ug(ugb%npw_k, u_nfft, nspinor, n1dat, &
   !            u_mgfft, u_ngfft, ugb%istwf_k, ugb%kg_k, gbound_k, &
   !            ug1_batch, &      ! in
   !            ur1_batch)        ! out

   call uplan_1%execute_gr(n1dat, ug1_batch(:,1), ur1_batch(:,1))

   if (ugb%istwf_k /= 2) ur1_batch = conjg(ur1_batch)  ! Not needed if k == Gamma
   ABI_FREE(ug1_batch)

   ! NB: Assuming bands distributed in contiguous blocks.
   do band2_start=ugb%my_bstart, ugb%my_bstop, batch2_size
     my_ib2 = band2_start - ugb%my_bstart + 1
     n2dat = blocked_loop(band2_start, ugb%my_bstop, batch2_size)
     band2_stop = band2_start + n2dat - 1

     !call fft_ug(ugb%npw_k, u_nfft, nspinor, n2dat, &
     !            u_mgfft, u_ngfft, ugb%istwf_k, ugb%kg_k, gbound_k, &
     !            ugb%mat%buffer_cplx(:,my_ib2), &     ! in
     !            ur2_batch(:,1))                      ! out

     call uplan_2%execute_gr(n2dat, ugb%mat%buffer_cplx(:,my_ib2), ur2_batch(:,1))

     do idat1=1,n1dat
       do idat2=1,n2dat
         ur12_batch(:, idat2) = ur1_batch(:, idat1) * ur2_batch(:, idat2)
       end do

       ! FIXME: nspinor
       !call fft_ur(m_npw, u_nfft, nspinor, n2dat, &
       !            u_mgfft, u_ngfft, istwfk1, m_gvec, m_gbound, &
       !            ur12_batch, &    ! in
       !            ug12_batch)      ! out

       call uplan_m%execute_rg(n2dat, ur12_batch(:,1), ug12_batch(:,1))

       !do idat2=1,n2dat
       !  write(std_out, *)ug12_batch(1, idat2), idat2, "ug12_batch(1, idat2), idat"
       !end do

     end do ! idat1
   end do ! my_ib2

   ! TODO: Write ug12_batch.
   !  - Handle parallel IO if nsppol 2 (we are inside the spin loop that is already MPI distributed!
   !  - Format for nspinor == 2?

 end do ! band1_start

 call uplan_1%free()
 call uplan_2%free()
 call uplan_m%free()

 ABI_FREE(gbound_k)
 ABI_FREE(m_gbound)
 ABI_FREE(m_gvec)
 ABI_FREE(ur1_batch)
 ABI_FREE(ur2_batch)
 ABI_FREE(ur12_batch)
 ABI_FREE(ug12_batch)

 call cwtime_report(" cc4s_gamma", cpu, wall, gflops)

end subroutine cc4s_gamma
!!***

end module m_gwr_driver
!!***

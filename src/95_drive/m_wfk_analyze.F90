!!****m* ABINIT/m_wfk_analyze
!! NAME
!!  m_wfk_analyze
!!
!! FUNCTION
!!  Post-processing tools for WFK file
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2026 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfk_analyze

 use, intrinsic :: iso_c_binding
 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hdr
 use m_crystal
 use m_ebands
 use m_nctk
 use m_wfd
 use m_dtset
 use m_dtfil
 use m_distribfft

 use m_io_tools,        only : iomode_from_fname, get_unit
 use defs_datatypes,    only : pseudopotential_type
 use defs_abitypes,     only : mpi_type
 use m_time,            only : timab
 use m_fstrings,        only : strcat, sjoin, itoa, ftoa, ltoa
 use m_fftcore,         only : print_ngfft
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq, init_mpi_enreg
 use m_esymm,           only : esymm_t, esymm_free
 use m_ddk,             only : ddkstore_t
 use m_ksdiago,         only : psbands_t
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init, pawfgrtab_print
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, pawrhoij_inquire_dim
 use m_pawdij,          only : pawdij, symdij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_sphharm,     only : setsym_ylm
 use m_paw_init,        only : pawinit, paw_gencond
 use m_paw_nhat,        only : nhatgrid
 use m_paw_tools,       only : chkpawovlp
 use m_paw_correlations,only : pawpuxinit
 use m_paw_pwaves_lmn,  only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_classify_bands,  only : classify_bands
 use m_pspini,          only : pspini
 use m_sigtk,           only : sigtk_kpts_in_erange
 use m_iowf,            only : prtkbff
 use m_wfd_wannier,     only : wfd_run_wannier
 use m_wfk,             only : wfk_to_bz, wfk_t, wfk_read_eigenvalues, wfk_check_symtab

 implicit none

 private
!!***

 public :: wfk_analyze
!!***

contains
!!***

!!****f* ABINIT/wfk_analyze
!! NAME
!!  wfk_analyze
!!
!! FUNCTION
!! Main routine implementing postprocessing tools for the WFK file.
!! Main differences wrt cut3d:
!!
!!   - MPI support.
!!   - No interactive prompt.
!!   - Run the analysis once, store all the important results in netcdf files
!!     and use python tools to analyze data
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
!!   Before entering the first time in the routine, a significant part of psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, psps might be identical to the
!!   one of the previous dtset, in which case, no reinitialisation is scheduled in pspini.F90.
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

subroutine wfk_analyze(acell, codvsn, dtfil, dtset, pawang, pawrad, pawtab, psps, rprim, xred)

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
 integer,parameter :: master = 0, formeig0 = 0
 integer :: comm,nprocs,my_rank,mgfftf,nfftf !,nfftf_tot
 integer :: optcut,optgr0,optgr1,optgr2,optrad,psp_gencond,ii
 !integer :: option,option_test,option_dij,optrhoij
 integer :: band,ik_ibz,spin,first_band,last_band, nband_k, islice, ib, mpw, mcg, nb, npw_k
 integer :: ierr,usexcnhat, sc_mode, nspinor, nsto
 integer :: cplex,cplex_dij,cplex_rhoij,ndij,nspden_rhoij,gnt_option
 real(dp),parameter :: spinmagntarget=-99.99_dp
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff,gsqcut_shp, gs_fermie
 !real(dp) :: cpu,wall,gflops
 !real(dp) :: ex_energy,gsqcutc_eff,gsqcutf_eff,nelect,norm,oldefermi
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path, outwfk_path
 logical :: call_pawinit, use_paw_aeur
 type(hdr_type) :: wfk0_hdr, hdr_bz, out_hdr
 type(crystal_t) :: cryst, cryst_dtset
 type(ebands_t) :: ebands, ebands_bz
 type(pawfgr_type) :: pawfgr
 !type(paw_dmft_type) :: paw_dmft
 type(mpi_type) :: mpi_enreg
 type(wfd_t) :: wfd
 type(ddkstore_t) :: ds
 type(wfk_t) :: in_wfk, out_wfk
 !type(dataset_type) :: my_dtset
!arrays
 integer :: ngfftc(18),ngfftf(18), units(2), band_block(2), bstart
 integer,allocatable :: l_size_atm(:), kg_k(:,:)
 real(dp),parameter :: k0(3)=zero
 real(dp),pointer :: gs_eigen(:,:,:)
 real(dp),allocatable :: eig_k(:), occ_k(:), thetas(:) !, out_cg(:,:), work(:,:,:,:), allcg_k(:,:)
 real(dp),allocatable,target :: cg_k(:,:)
 complex(gwp),allocatable :: ur_ae(:)
 complex(dp),pointer :: cg_k_cplx(:,:)
 complex(dp),allocatable :: ps_ug(:,:)
 logical,allocatable :: keep_ur(:,:,:),bks_mask(:,:,:)
 real(dp) :: tsec(2)
 type(Pawrhoij_type),allocatable :: pawrhoij(:)
 type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)
 type(esymm_t),allocatable :: esymm(:,:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)
 type(psbands_t),allocatable :: psb_ks(:,:)
!************************************************************************

 DBG_ENTER('COLL')

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred
 units = [std_out, ab_out]

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 wfk0_path = dtfil%fnamewffk
 if (my_rank == master) then
   ! Accept WFK file in Fortran or netcdf format.
   if (nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     ABI_ERROR(sjoin("Cannot find GS WFK file:", ch10, msg))
   end if
 end if
 call xmpi_bcast(wfk0_path, master, comm, ierr)
 call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path))

 !call cwtime(cpu,wall,gflops,"start")

 ! Construct crystal and ebands from the GS WFK file.
 call wfk_read_eigenvalues(wfk0_path, gs_eigen, wfk0_hdr, comm) !,gs_occ)
 call wfk0_hdr%vs_dtset(dtset)
 nspinor = dtset%nspinor

 ! Get fermie from the GS calculation.
 ! NB: It might understimate the real Fermi level, especially if the den was computed on a shifted k-mesh
 ! at present it's only used to implement pseudobands
 gs_fermie = wfk0_hdr%fermie

 cryst = wfk0_hdr%get_crystal()
 call cryst%print(header="Crystal structure from WFK file")

 ! Compare structure with the one computed from input file.
 cryst_dtset = dtset%get_crystal(1)
 if (cryst%compare(cryst_dtset, header=" Comparing WFK crystal with crystal from dtset") /= 0) then
   ABI_ERROR("Crystal structure from WFK and dataser do not agree! Check messages above!")
 end if
 call cryst_dtset%free()

 call ebands%from_hdr(wfk0_hdr, maxval(wfk0_hdr%nband), gs_eigen)

 !call ebands%update_occ(spinmagntarget)
 call ebands%print([std_out], header="Ground state energies", prtvol=dtset%prtvol)
 ABI_FREE(gs_eigen)

 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
                  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)

 call print_ngfft([std_out], ngfftc, header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft([std_out], ngfftf, header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call mpi_enreg%distribfft%init_seq('c',ngfftc(2),ngfftc(3),'all')
 call mpi_enreg%distribfft%init_seq('f',ngfftf(2),ngfftf(3),'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,pawrad,pawtab,psps,cryst%rprimd,comm_mpi=comm)

 ! ============================
 ! ==== PAW initialization ====
 ! ============================
 if (dtset%usepaw == 1) then
   call chkpawovlp(cryst%natom,cryst%ntypat,dtset%pawovlp,pawtab,cryst%rmet,cryst%typat,cryst%xred)

   cplex_dij=nspinor; cplex=1; ndij=1

   ABI_MALLOC(pawrhoij,(cryst%natom))
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
                             nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
   call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,nspinor,dtset%nsppol,cryst%typat,pawtab=pawtab)

   ! Initialize values for several basic arrays
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1 .or. call_pawinit) then
     call timab(553,1,tsec)
     gsqcut_shp = two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     call pawinit(dtset%effmass_free,gnt_option,gsqcut_shp,zero,dtset%pawlcutd,dtset%pawlmix,&
                  psps%mpsang,dtset%pawnphi,cryst%nsym,dtset%pawntheta,pawang,Pawrad,&
                  dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%ixc,dtset%usepotzero)
     call timab(553,2,tsec)

     ! Update internal values
     call paw_gencond(dtset,gnt_option,"save",call_pawinit)

   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:cryst%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:cryst%ntypat)%has_nabla=2
   end if

   psps%n1xccc=MAXVAL(pawtab(1:cryst%ntypat)%usetcore)

   ! Initialize optional flags in pawtab to zero
   ! (Cannot be done in Pawinit since the routine is called only if some pars. are changed)
   pawtab(:)%has_nabla = 0
   pawtab(:)%usepawu   = 0
   pawtab(:)%useexexch = 0
   pawtab(:)%exchmix   =zero
   pawtab(:)%lamb_shielding   =zero

   call setsym_ylm(cryst%gprimd,pawang%l_max-1,cryst%nsym,dtset%pawprtvol,cryst%rprimd,cryst%symrec,pawang%zarot)

   ! Initialize and compute data for DFT+U
   !paw_dmft%use_dmft=dtset%usedmft
   !call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
   !    .false.,dtset%jpawu,dtset%lexexch,dtset%lpawu,cryst%ntypat,pawang,dtset%pawprtvol,&
   !    Pawrad,pawtab,dtset%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu)
   !ABI_CHECK(paw_dmft%use_dmft==0,"DMFT not available")
   !call destroy_sc_dmft(paw_dmft)

   if (my_rank == master) call pawtab_print(pawtab, unit=std_out)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(wfk0_hdr%pawrhoij,pawrhoij)

   ! Variables/arrays related to the fine FFT grid.
   ABI_MALLOC(pawfgrtab,(cryst%natom))
   call pawtab_get_lsize(pawtab,l_size_atm,cryst%natom,cryst%typat)
   cplex=1
   call pawfgrtab_init(pawfgrtab,cplex,l_size_atm,dtset%nspden,dtset%typat)
   ABI_FREE(l_size_atm)

   usexcnhat=maxval(pawtab(:)%usexcnhat)
   ! 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
   ! 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   call wrtout(std_out,sjoin("using usexcnhat= ",itoa(usexcnhat)))
   !
   ! Identify parts of the rectangular grid where the density has to be calculated ===
   !optcut=0; optgr0=dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-dtset%pawstgylm
   !if (dtset%xclevel==2 .and. usexcnhat>0) optgr1=dtset%pawstgylm
   optcut=1; optgr0=1; optgr1=1; optgr2=1; optrad=1

   call nhatgrid(cryst%atindx1,cryst%gmet,cryst%natom,cryst%natom,cryst%nattyp,ngfftf,cryst%ntypat,&
                optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,cryst%rprimd,cryst%typat,cryst%ucvol,cryst%xred)

   call pawfgrtab_print(pawfgrtab,cryst%natom,unit=std_out,prtvol=dtset%pawprtvol)

   !ABI_MALLOC(ks_nhat,(nfftf,dtset%nspden))
   !ks_nhat=zero
 else
   ABI_MALLOC(pawfgrtab,(0))
 end if !End of PAW Initialization

 select case (dtset%wfk_task)

 case (WFK_TASK_FULLBZ, WFK_TASK_OPTICS_FULLBZ)
   ! Read wfk0_path and build WFK in full BZ.
   if (my_rank == master) then
     outwfk_path = dtfil%fnameabo_wfk; if (dtset%iomode == IO_MODE_ETSF) outwfk_path = nctk_ncify(outwfk_path)
     call wfk_to_bz(wfk0_path, dtset, psps, pawtab, outwfk_path, hdr_bz, ebands_bz)
     call ebands_bz%free()

     ! Write KB form factors.
     if (dtset%prtkbff == 1 .and. dtset%iomode == IO_MODE_ETSF .and. dtset%usepaw == 0) then
       call prtkbff(outwfk_path, hdr_bz, psps, dtset%prtvol)
     end if
     call hdr_bz%free()
   end if
   call xmpi_barrier(comm)

   if (dtset%wfk_task == WFK_TASK_OPTICS_FULLBZ) then
     ! Calculate the DDK matrix elements from the WFK file in the full BZ.
     ! This is needed for computing non-linear properties in optics as symmetries are not
     ! implemented correctly.
     ds%only_diago = .False.
     call ds%compute_ddk(outwfk_path, dtfil%filnam_ds(4), dtset, psps, pawtab, ngfftc, comm)
     call ds%free()
   end if

 case (WFK_TASK_KPTS_ERANGE)
   call sigtk_kpts_in_erange(dtset, cryst, ebands, psps, pawtab, dtfil%filnam_ds(4), comm)

 case (WFK_TASK_DDK, WFK_TASK_DDK_DIAGO)
   ! Calculate the DDK matrix elements from the WFK file
   ds%only_diago = .False.; if (dtset%wfk_task == WFK_TASK_DDK_DIAGO) ds%only_diago = .True.
   call ds%compute_ddk(wfk0_path, dtfil%filnam_ds(4), dtset, psps, pawtab, ngfftc, comm)
   call ds%free()

 case (WFK_TASK_EINTERP)
   ! Band structure interpolation from eigenvalues computed on the k-mesh.
   call ebands%interpolate_kpath(dtset, cryst, [0, 0], dtfil%filnam_ds(4), comm)

 case (WFK_TASK_CHECK_SYMTAB)
   call wfk_check_symtab(wfk0_path, comm)

 case (WFK_TASK_CLASSIFY)
   ! Band classification.
   call read_wfd()

   ABI_MALLOC(esymm,(wfd%nkibz,wfd%nsppol))
   use_paw_aeur=.False. ! should pass ngfftf but the dense mesh is not forced to be symmetric

   do spin=1,wfd%nsppol
     do ik_ibz=1,wfd%nkibz
       first_band = 1
       last_band  = wfd%nband(ik_ibz,spin)
       call classify_bands(wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,wfd%ngfft,&
       cryst,ebands,pawtab,pawrad,pawang,psps,dtset%tolsym,esymm(ik_ibz,spin))
     end do
   end do

   call esymm_free(esymm)
   ABI_FREE(esymm)

 !case (WFK_TASK_UR)
 !  ! plot KSS wavefunctions. Change bks_mask to select particular states.
 !  ABI_MALLOC(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 !  bks_mask=.False.; bks_mask(1:4,1,1)=.True.
 !  call wfd%plot_ur(Cryst,Psps,Pawtab,Pawrad,ngfftf,bks_mask)
 !  ABI_FREE(bks_mask)

 case (WFK_TASK_PSEUDOBANDS)
   if (my_rank /= master) goto 100 ! NO MPI parallelism here

   ! out_hdr is the header the STO_WFK file.
   call wfk0_hdr%copy(out_hdr)

   ! Pre-compute slices for all k-points and spin so that we know the new number of bands in STO_WFK file.
   ABI_MALLOC(psb_ks, (ebands%nkpt, ebands%nsppol))
   do spin=1,ebands%nsppol
     do ik_ibz=1,ebands%nkpt
       associate (psb => psb_ks(ik_ibz, spin))
       nband_k = ebands%nband(ik_ibz + (spin-1)*ebands%nkpt)
       call psb%init(dtset, nband_k, ebands%eig(:, ik_ibz, spin), gs_fermie)
       ! Change the number of bands to account for pseudo bands.
       !print *, "nb_tot:", psb%nb_tot
       out_hdr%nband(ik_ibz + (spin-1)*ebands%nkpt) = psb%nb_tot
       end associate
     end do
   end do

   ! Compute new value of bantot
   ! TODO: Have to change all arrays in outhdr_hdr depending on nband_ks
   out_hdr%bantot = sum(out_hdr%nband)
   out_hdr%mband = maxval(out_hdr%nband)
   ABI_RECALLOC(out_hdr%occ, (out_hdr%bantot))

   outwfk_path = strcat(dtfil%filnam_ds(4), "_STO_WFK")
   if (dtset%iomode == IO_MODE_ETSF) outwfk_path = nctk_ncify(outwfk_path)
   call out_wfk%open_write(out_hdr, outwfk_path, formeig0, dtset%iomode, get_unit(), xmpi_comm_self)

   call in_wfk%open_read(wfk0_path, formeig0, iomode_from_fname(wfk0_path), get_unit(), xmpi_comm_self)

   ! The output arrays eig_k and occ_k contain the *full* set of eigenvalues and occupation
   ! factors stored in the file and are dimensioned with mband.
   ABI_MALLOC(eig_k, (wfk0_hdr%mband))
   ABI_MALLOC(occ_k, (wfk0_hdr%mband))
   mpw = maxval(wfk0_hdr%npwarr)
   ABI_MALLOC(kg_k, (3, mpw))
   !print *, "mpw:", mpw
   sc_mode = xmpio_single

   do spin=1,ebands%nsppol
     do ik_ibz=1,ebands%nkpt
       associate (psb => psb_ks(ik_ibz, spin))
       ! Read and write protected states.
       band_block = [1, psb%nb_protected]
       nb = band_block(2) - band_block(1) + 1
       npw_k = wfk0_hdr%npwarr(ik_ibz)
       mcg = npw_k * nspinor * nb
       ABI_MALLOC(cg_k, (2, mcg))

       call in_wfk%read_band_block(band_block, ik_ibz, spin, sc_mode, &
                                   kg_k=kg_k, cg_k=cg_k, eig_k=eig_k, occ_k=occ_k)

       !call wrtout(std_out, sjoin(" About to write islice:", itoa(0), "with band block:", ltoa(band_block)))
       eig_k(1:psb%nb_tot) = psb%ps_eig(:) ! Change eigenvalues to account for pseudo bands
       !occ_k = ???
       call out_wfk%write_band_block(band_block, ik_ibz, spin, sc_mode, &
                                     kg_k=kg_k, cg_k=cg_k, eig_k=eig_k, occ_k=occ_k)
       ABI_FREE(cg_k)

       ! ========================================
       ! Build pseudobands and write them to disk
       ! ========================================
       bstart = psb%nb_protected + 1
       do islice=1,psb%nslices
         band_block = psb%subspace(1:2, islice)
         nb = band_block(2) - band_block(1) + 1
         mcg = npw_k * nspinor * nb
         ABI_MALLOC(cg_k, (2, mcg))

         call in_wfk%read_band_block(band_block, ik_ibz, spin, sc_mode, cg_k=cg_k)
         call c_f_pointer(c_loc(cg_k), cg_k_cplx, [npw_k*nspinor, nb])

         ! Allocate pseudobands.
         nsto = psb%subspace(3, islice)
         ABI_CALLOC(ps_ug, (npw_k*nspinor, nsto))
         ABI_MALLOC(thetas, (nb))

         if (nsto == 1) then
           ! Use KS state.
           ps_ug = cg_k_cplx
         else
           ! Multiply by random phases.
           do ii=1,nsto
             call random_number(thetas)
             do ib=1,nb
               ps_ug(:,ii) = ps_ug(:,ii) + cg_k_cplx(:,ib) * exp(j_dpc*two_pi*thetas(ib)) / sqrt(one * nsto)
             end do
           end do
         end if

         band_block = [bstart, bstart + psb%subspace(3, islice) - 1]
         !call wrtout(std_out, sjoin(" About to write islice:", itoa(islice), "with band block:", ltoa(band_block)))
         call out_wfk%write_band_block(band_block, ik_ibz, spin, sc_mode, cg_k=cg_k)
                                       !kg_k=kg_k, cg_k=cg_k, eig_k=eig_k, occ_k=occ_k)
         bstart = bstart + psb%subspace(3, islice)

         ABI_FREE(thetas)
         ABI_FREE(cg_k)
         ABI_FREE(ps_ug)
       end do ! islice

       call psb%free()
       end associate
     end do ! ik_ibz
   end do ! spin

   ABI_FREE(eig_k)
   ABI_FREE(occ_k)
   ABI_FREE(kg_k)
   call in_wfk%close(); call out_wfk%close(); call out_hdr%free()

   ! DEBUG section. Try to read the output WFK file.
   !call out_wfk%open_read(outwfk_path, formeig0, iomode_from_fname(outwfk_path), get_unit(), xmpi_comm_self)
   !do spin=1,ebands%nsppol
   !  do ik_ibz=1,ebands%nkpt
   !    nband_k = out_wfk%hdr%nband(ik_ibz + (spin-1)*ebands%nkpt)
   !    band_block = [1, nband_k]
   !    call out_wfk%read_band_block(band_block, ik_ibz, spin, sc_mode, &
   !                                 kg_k=kg_k, cg_k=cg_k, eig_k=eig_k, occ_k=occ_k)
   !  end do
   !end do
   !call out_wfk%close()

 case (WFK_TASK_PAW_AEPSI)
   ! Compute AE PAW wavefunction in real space on the dense FFT mesh.
   call read_wfd()

   ABI_CHECK(wfd%usepaw == 1, "Not a PAW run")
   ABI_MALLOC(paw_onsite, (cryst%natom))
   call paw_pwaves_lmn_init(paw_onsite,cryst%natom,cryst%natom,cryst%ntypat, &
                            cryst%rprimd,cryst%xcart,pawtab,pawrad,pawfgrtab)

   ! Use dense FFT mesh
   call wfd%change_ngfft(cryst,psps,ngfftf)
   band = 1; spin = 1; ik_ibz = 1

   ABI_MALLOC(ur_ae, (wfd%nfft*wfd%nspinor))
   call wfd%paw_get_aeur(band,ik_ibz,spin,cryst,paw_onsite,psps,pawtab,pawfgrtab,ur_ae)
   ABI_FREE(ur_ae)

   call paw_pwaves_lmn_free(paw_onsite)
   ABI_FREE(paw_onsite)

 case (WFK_TASK_WANNIER)
   ! Construct Wannier functions.

   ! This part was implemented by gmatteo to debug GaAs with a 8x8x8 k-mesh.

   !if (wfk0_hdr%kptopt == 1) then
   !  ! Generate WFK in the full BZ (only master works here)
   !  outwfk_path = dtfil%fnameabo_wfk; if (dtset%iomode == IO_MODE_ETSF) outwfk_path = nctk_ncify(outwfk_path)
   !  if (my_rank == master) then
   !    call wrtout(units, sjoin("- Generating WFK file with kpoints in the full BZ and istwfk == 1", outwfk_path))
   !    call wfk_to_bz(wfk0_path, dtset, psps, pawtab, outwfk_path, hdr_bz, ebands_bz)
   !    call ebands_bz%free(); call hdr_bz%free()
   !  end if
   !  call xmpi_barrier(comm)
   !  my_dtset = dtset%copy()
   !  ebands_bz = wfk_read_ebands(outwfk_path, comm, hdr_bz)
   !  call hdr_transfer_nkpt_arrays(hdr_bz, my_dtset)
   !  my_dtset%kptopt = hdr_bz%kptopt
   !  call hdr_bz%vs_dtset(my_dtset)
   !  call wfd_run_wannier__(outwfk_path, my_dtset, ebands_bz, hdr_bz)
   !  call ebands_bz%free(); call hdr_bz%free(); call my_dtset%free()

   !else
     call wfk0_hdr%vs_dtset(dtset)
     call wfd_run_wannier__(wfk0_path, dtset, ebands, wfk0_hdr)
   !end if

 case default
   ABI_ERROR(sjoin("Wrong wfk_task:", itoa(dtset%wfk_task)))
 end select

100 continue

 ! Free memory
 call cryst%free(); call ebands%free(); call wfd%free(); call destroy_mpi_enreg(mpi_enreg); call wfk0_hdr%free()
 call pawfgr_destroy(pawfgr)

 ! Deallocation for PAW.
 if (dtset%usepaw==1) then
   call pawrhoij_free(pawrhoij)
   ABI_FREE(pawrhoij)
   call pawfgrtab_free(pawfgrtab)
   !call paw_ij_free(paw_ij)
   !ABI_FREE(paw_ij)
   !call paw_an_free(paw_an)
   !ABI_FREE(paw_an)
 end if
 ABI_FREE(pawfgrtab)

 DBG_EXIT('COLL')

 contains
!!***

!!****f* wfk_analyze/read_wfd
!! NAME
!!  read_wfd
!!
!! FUNCTION
!!  Initialize the wavefunction descriptor from file.
!!
!! SOURCE

subroutine read_wfd()

 ABI_MALLOC(keep_ur, (ebands%mband, ebands%nkpt, ebands%nsppol))
 ABI_MALLOC(bks_mask, (ebands%mband, ebands%nkpt, ebands%nsppol))
 keep_ur = .False.; bks_mask = .True.

 call wfd%init(cryst,pawtab,psps,keep_ur,ebands%mband,ebands%nband,ebands%nkpt,dtset%nsppol,bks_mask,&
   dtset%nspden,dtset%nspinor,ecut_eff,dtset%ecutsm,dtset%dilatmx,wfk0_hdr%istwfk,ebands%kptns,ngfftc,&
   dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

 ABI_FREE(keep_ur)
 ABI_FREE(bks_mask)

 call wfd%read_wfk(wfk0_path, iomode_from_fname(wfk0_path))

end subroutine read_wfd

subroutine wfd_run_wannier__(wfk_filepath, dtset_, ebands_, hdr_)

 type(dataset_type),intent(in) :: dtset_
 character(len=*),intent(in) :: wfk_filepath
 type(ebands_t),intent(in) :: ebands_
 type(hdr_type),intent(in) :: hdr_

 ABI_MALLOC(keep_ur, (ebands_%mband, ebands_%nkpt, ebands_%nsppol))
 ABI_MALLOC(bks_mask, (ebands_%mband, ebands_%nkpt, ebands_%nsppol))
 keep_ur = .False.; bks_mask = .True.

 ! Impose istwfk = 1 for all k-points. This is also done in respfn (see inkpts)
 ! wfd_read_wfk will handle a possible conversion if WFK contains istwfk /= 1.
 !wfk0_hdr%istwfk = 1; ebands%istwfk = 1; dtset%istwfk = 1

 call wfd%init(cryst, pawtab, psps, keep_ur, ebands_%mband, ebands_%nband, ebands_%nkpt, dtset_%nsppol, bks_mask, &
   dtset_%nspden,  dtset_%nspinor, ecut_eff, dtset_%ecutsm, dtset_%dilatmx, hdr_%istwfk, ebands_%kptns, ngfftc, &
   dtset_%nloalg, dtset_%prtvol, dtset_%pawprtvol, comm)

 ABI_FREE(keep_ur)
 ABI_FREE(bks_mask)
 call wfd%read_wfk(wfk_filepath, iomode_from_fname(wfk_filepath))

 call wfd_run_wannier(cryst=cryst, ebands=ebands_, hdr=hdr_, mpi_enreg=mpi_enreg, &
                      ngfftc=ngfftc, ngfftf=ngfftf, wfd=wfd, dtset=dtset_, dtfil=dtfil,  &
                      pawang=pawang, pawrad=pawrad, pawtab=pawtab, psps=psps)

end subroutine wfd_run_wannier__

end subroutine wfk_analyze
!!***

subroutine hdr_transfer_nkpt_arrays(hdr, dtset)

  use m_copy, only : alloc_copy

  class(hdr_type),intent(in) :: hdr
  type(dataset_type),intent(inout) :: dtset

  ABI_SFREE(dtset%istwfk)
  ABI_SFREE(dtset%kpt)
  ABI_SFREE(dtset%kptns)
  ABI_SFREE(dtset%occ_orig)
  ABI_SFREE(dtset%wtk)
  ABI_SFREE(dtset%kptns_hf)  ! Free HF k-points as well.
  ABI_SFREE(dtset%nband)

  dtset%nkpt = hdr%nkpt
  call alloc_copy(hdr%istwfk, dtset%istwfk)
  call alloc_copy(hdr%nband, dtset%nband)
  call alloc_copy(hdr%kptns, dtset%kpt)
  call alloc_copy(hdr%kptns, dtset%kptns)
  !call alloc_copy(hdr%occ, dtset%occ_orig(:,1)
  call alloc_copy(hdr%wtk, dtset%wtk)
  call alloc_copy(hdr%kptns, dtset%kptns_hf)

end subroutine hdr_transfer_nkpt_arrays

end module m_wfk_analyze
!!***

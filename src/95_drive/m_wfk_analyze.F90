!!****m* ABINIT/m_wfk_analyze
!! NAME
!!  m_wfk_analyze
!!
!! FUNCTION
!!  Post-processing tools for WFK file
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2020 ABINIT group (MG)
!!  This file is distributed under the terms of the
!!  GNU General Public License, see ~abinit/COPYING
!!  or http://www.gnu.org/copyleft/gpl.txt .
!!
!! PARENTS
!!
!! CHILDREN
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

module m_wfk_analyze

 use defs_basis
 use m_abicore
 use m_xmpi
 use m_errors
 use m_hdr
 use m_crystal
 use m_ebands
 use m_nctk
 use m_wfk
 use m_wfd
 use m_dtset
 use m_dtfil
 use m_distribfft

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : mpi_type
 use m_time,            only : timab
 use m_fstrings,        only : strcat, sjoin, itoa, ftoa
 use m_fftcore,         only : print_ngfft
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_esymm,           only : esymm_t, esymm_free
 use m_ddk,             only : ddk_compute
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
!! PARENTS
!!      driver
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
!! CHILDREN
!!
!! SOURCE

subroutine wfk_analyze(acell,codvsn,dtfil,dtset,pawang,pawrad,pawtab,psps,rprim,xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: master=0
 integer :: comm,nprocs,my_rank,mgfftf,nfftf !,nfftf_tot
 integer :: optcut,optgr0,optgr1,optgr2,optrad,psp_gencond !,ii
 !integer :: option,option_test,option_dij,optrhoij
 integer :: band,ik_ibz,spin,first_band,last_band
 integer :: ierr,usexcnhat
 integer :: cplex,cplex_dij,cplex_rhoij,ndij,nspden_rhoij,gnt_option
 real(dp),parameter :: spinmagntarget=-99.99_dp
 real(dp) :: ecore,ecut_eff,ecutdg_eff,gsqcutc_eff,gsqcutf_eff,gsqcut_shp
 !real(dp) :: cpu,wall,gflops
 !real(dp) :: ex_energy,gsqcutc_eff,gsqcutf_eff,nelect,norm,oldefermi
 character(len=500) :: msg
 character(len=fnlen) :: wfk0_path,wfkfull_path
 logical :: call_pawinit, use_paw_aeur
 type(hdr_type) :: wfk0_hdr
 type(crystal_t) :: cryst
 type(ebands_t) :: ebands
 type(pawfgr_type) :: pawfgr
 !type(paw_dmft_type) :: paw_dmft
 type(mpi_type) :: mpi_enreg
 type(wfd_t) :: wfd
!arrays
 integer :: ngfftc(18),ngfftf(18)
 integer,allocatable :: l_size_atm(:)
 real(dp),parameter :: k0(3)=zero
 !real(dp) :: nelect_per_spin(dtset%nsppol),n0(dtset%nsppol)
 real(dp),pointer :: gs_eigen(:,:,:)
 complex(gwpc),allocatable :: ur_ae(:)
 logical,allocatable :: keep_ur(:,:,:),bks_mask(:,:,:)
 real(dp) :: tsec(2)
 type(Pawrhoij_type),allocatable :: pawrhoij(:)
 type(pawfgrtab_type),allocatable :: pawfgrtab(:)
 !type(paw_ij_type),allocatable :: paw_ij(:)
 !type(paw_an_type),allocatable :: paw_an(:)
 type(esymm_t),allocatable :: esymm(:,:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 DBG_ENTER('COLL')

 ! abirules!
 if (.False.) write(std_out,*)acell,codvsn,rprim,xred

 comm = xmpi_world; nprocs = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 wfk0_path = dtfil%fnamewffk
 if (my_rank == master) then
   ! Accept WFK file in Fortran or netcdf format.
   if (nctk_try_fort_or_ncfile(wfk0_path, msg) /= 0) then
     MSG_ERROR(sjoin("Cannot find GS WFK file:", ch10, msg))
   end if
 end if
 call xmpi_bcast(wfk0_path, master, comm, ierr)
 call wrtout(ab_out, sjoin("- Reading GS states from WFK file:", wfk0_path) )

 !call cwtime(cpu,wall,gflops,"start")

 ! Costruct crystal and ebands from the GS WFK file.
 call wfk_read_eigenvalues(wfk0_path,gs_eigen,wfk0_hdr,comm) !,gs_occ)
 call wfk0_hdr%vs_dtset(dtset)

 cryst = wfk0_hdr%get_crystal()
 call cryst%print(header="crystal structure from WFK file")

 ebands = ebands_from_hdr(wfk0_hdr,maxval(wfk0_hdr%nband),gs_eigen)

 ! TODO:
 ! Make sure everything is OK if WFK comes from a NSCF run since occ are set to zero
 ! fermie is set to 0 if nscf!

 ! Here we change the GS bands (fermi level, scissors operator ...)
 ! All the modifications to ebands should be done here.

 !if (dtset%occopt /= ebands%occopt .or. abs(dtset%tsmear - ebands%tsmear) > tol12) then
 !  write(msg,"(2a,2(a,i0,a,f14.6,a))")&
 !    " Changing occupation scheme as input occopt and tsmear differ from those read from WFK file.",ch10,&
 !    "   From WFK file: occopt = ",ebands%occopt,", tsmear = ",ebands%tsmear,ch10,&
 !    "   From input:    occopt = ",dtset%occopt,", tsmear = ",dtset%tsmear,ch10
 !  call wrtout(ab_out,msg)
 !  call ebands_set_scheme(ebands,dtset%occopt,dtset%tsmear,spinmagntarget,dtset%prtvol)
 !end if

 !if (dtset%eph_fermie /= zero) then ! default value of eph_fermie is zero hence no tolerance is used!
 !  ABI_CHECK(abs(dtset%eph_extrael) <= tol12, "eph_fermie and eph_extrael are mutually exclusive")
 !  call wrtout(ab_out, sjoin(" Fermi level set by the user at:",ftoa(dtset%eph_fermie)))
 !  call ebands_set_fermie(ebands, dtset%eph_fermie, msg)
 !  call wrtout(ab_out,msg)

 !else if (abs(dtset%eph_extrael) > tol12) then
 !  NOT_IMPLEMENTED_ERROR()
 !  ! TODO: Be careful with the trick used in elphon for passing the concentration
 !  !call ebands_set_nelect(ebands, dtset%eph_extrael, spinmagntarget, msg)
 !  call wrtout(ab_out,msg)
 !end if

 !call ebands_update_occ(ebands, spinmagntarget)
 call ebands_print(ebands,header="Ground state energies",prtvol=dtset%prtvol)
 ABI_FREE(gs_eigen)

 ! TODO Recheck getng, should use same trick as that used in screening and sigma.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
                  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=cryst%gmet,k0=k0)

 call print_ngfft(ngfftc,header='Coarse FFT mesh used for the wavefunctions')
 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(mpi_enreg)
 call init_distribfft_seq(mpi_enreg%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(mpi_enreg%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,pawrad,pawtab,psps,cryst%rprimd,comm_mpi=comm)

 ! ============================
 ! ==== PAW initialization ====
 ! ============================
 if (dtset%usepaw == 1) then
   call chkpawovlp(cryst%natom,cryst%ntypat,dtset%pawovlp,pawtab,cryst%rmet,cryst%typat,cryst%xred)

   cplex_dij=dtset%nspinor; cplex=1; ndij=1

   ABI_MALLOC(pawrhoij,(cryst%natom))
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
                             nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
   call pawrhoij_alloc(pawrhoij,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,cryst%typat,pawtab=pawtab)

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

 case (WFK_TASK_FULLBZ)
   ! Read wfk0_path and build WFK in full BZ.
   if (my_rank == master) then
     wfkfull_path = dtfil%fnameabo_wfk; if (dtset%iomode == IO_MODE_ETSF) wfkfull_path = nctk_ncify(wfkfull_path)
     call wfk_tofullbz(wfk0_path, dtset, psps, pawtab, wfkfull_path)
   end if
   call xmpi_barrier(comm)

 case (WFK_TASK_KPTS_ERANGE)
   call sigtk_kpts_in_erange(dtset, cryst, ebands, psps, pawtab, dtfil%filnam_ds(4), comm)

 !case (WFK_TASK_LDOS)

 case (WFK_TASK_DDK, WFK_TASK_DDK_DIAGO)
   ! Calculate the DDK matrix elements from the WFK file
   call ddk_compute(wfk0_path, dtfil%filnam_ds(4), dtset, psps, pawtab, ngfftc, comm)

 case (WFK_TASK_EINTERP)
   ! Band structure interpolation from eigenvalues computed on the k-mesh.
   call ebands_interpolate_kpath(ebands, dtset, cryst, [0, 0], dtfil%filnam_ds(4), comm)

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

 case (WFK_TASK_PAW_AEPSI)
   ! Compute AE PAW wavefunction in real space on the dense FFT mesh.
   call read_wfd()

   ABI_CHECK(wfd%usepaw == 1, "Not a PAW run")
   ABI_MALLOC(paw_onsite, (cryst%natom))
   call paw_pwaves_lmn_init(paw_onsite,cryst%natom,cryst%natom,cryst%ntypat,&
   cryst%rprimd,cryst%xcart,pawtab,pawrad,pawfgrtab)

   ! Use dense FFT mesh
   call wfd%change_ngfft(cryst,psps,ngfftf)
   band = 1; spin = 1; ik_ibz = 1

   ABI_MALLOC(ur_ae, (wfd%nfft*wfd%nspinor))
   call wfd%paw_get_aeur(band,ik_ibz,spin,cryst,paw_onsite,psps,pawtab,pawfgrtab,ur_ae)
   ABI_FREE(ur_ae)

   call paw_pwaves_lmn_free(paw_onsite)
   ABI_FREE(paw_onsite)

 !case ("paw_aeden")

 case default
   MSG_ERROR(sjoin("Wrong task:", itoa(dtset%wfk_task)))
 end select

 ! Free memory
 call cryst%free()
 call ebands_free(ebands)
 call wfd%free()
 call pawfgr_destroy(pawfgr)
 call destroy_mpi_enreg(mpi_enreg)
 call wfk0_hdr%free()

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
!!  Initialize the wavefunction descriptor
!!
!! PARENTS
!!      wfk_analyze
!!
!! CHILDREN
!!
!! SOURCE

subroutine read_wfd()

! *************************************************************************

   ABI_MALLOC(keep_ur, (ebands%mband, ebands%nkpt, ebands%nsppol))
   ABI_MALLOC(bks_mask, (ebands%mband, ebands%nkpt, ebands%nsppol))
   keep_ur = .False.; bks_mask = .True.

   call wfd_init(wfd,cryst,pawtab,psps,keep_ur,ebands%mband,ebands%nband,ebands%nkpt,dtset%nsppol,bks_mask,&
     dtset%nspden,dtset%nspinor,ecut_eff,dtset%ecutsm,dtset%dilatmx,wfk0_hdr%istwfk,ebands%kptns,ngfftc,&
     dtset%nloalg,dtset%prtvol,dtset%pawprtvol,comm)

   ABI_FREE(keep_ur)
   ABI_FREE(bks_mask)

   call wfd%read_wfk(wfk0_path,IO_MODE_MPI)
   !call wfd%test_ortho(cryst, pawtab)

 end subroutine read_wfd

end subroutine wfk_analyze
!!***

end module m_wfk_analyze
!!***

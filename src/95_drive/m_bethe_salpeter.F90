!!****m* ABINIT/m_bethe_salpeter
!! NAME
!!  m_bethe_salpeter
!!
!! FUNCTION
!!  Main routine to calculate dielectric properties by solving the Bethe-Salpeter equation in
!!  Frequency-Reciprocal space on a transition (electron-hole) basis set.
!!
!! COPYRIGHT
!! Copyright (C) 1992-2009 EXC group (L.Reining, V.Olevano, F.Sottile, S.Albrecht, G.Onida)
!! Copyright (C) 2009-2020 ABINIT group (MG, YG)
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

module m_bethe_salpeter

 use defs_basis
 use defs_wvltypes
 use m_bs_defs
 use m_abicore
 use m_xmpi
 use m_errors
 use m_screen
 use m_nctk
 use m_distribfft
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_dtset
 use m_dtfil
 use m_crystal

 use defs_datatypes,    only : pseudopotential_type, ebands_t
 use defs_abitypes,     only : MPI_type
 use m_gwdefs,          only : GW_Q0_DEFAULT
 use m_time,            only : timab
 use m_fstrings,        only : strcat, sjoin, endswith, itoa
 use m_io_tools,        only : file_exists, iomode_from_fname
 use m_geometry,        only : mkrdim, metric, normv
 use m_hide_lapack,     only : matrginv
 use m_mpinfo,          only : destroy_mpi_enreg, initmpi_seq
 use m_fftcore,         only : print_ngfft
 use m_fft_mesh,        only : rotate_FFT_mesh, get_gftt, setmesh
 use m_fft,             only : fourdp
 use m_bz_mesh,         only : kmesh_t, kmesh_init, kmesh_free, get_ng0sh, kmesh_print, get_BZ_item, find_qmesh, make_mesh
 use m_double_grid,     only : double_grid_t, double_grid_init, double_grid_free
 use m_ebands,          only : ebands_init, ebands_print, ebands_copy, ebands_free, &
                               ebands_update_occ, ebands_get_valence_idx, ebands_apply_scissors, ebands_report_gap
 use m_kg,              only : getph
 use m_gsphere,         only : gsphere_t, gsph_free, gsph_init, print_gsphere, gsph_extend
 use m_vcoul,           only : vcoul_t, vcoul_init, vcoul_free
 use m_qparticles,      only : rdqps, rdgw  !, show_QP , rdgw
 use m_wfd,             only : wfd_init, wfd_t, test_charge
 use m_wfk,             only : wfk_read_eigenvalues
 use m_energies,        only : energies_type, energies_init
 use m_io_screening,    only : hscr_t, hscr_free, hscr_io, hscr_bcast, hscr_from_file, hscr_print
 use m_haydock,         only : exc_haydock_driver
 use m_exc_diago,       only : exc_diago_driver
 use m_exc_analyze,     only : exc_den
 use m_eprenorms,       only : eprenorms_t, eprenorms_free, eprenorms_from_epnc, eprenorms_bcast
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free,&
&                              pawrhoij_inquire_dim, pawrhoij_symrhoij
 use m_pawdij,          only : pawdij, symdij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_hr,          only : pawhur_t, pawhur_free, pawhur_init
 use m_pawpwij,         only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_paw_sphharm,     only : setsym_ylm
 use m_paw_denpot,      only : pawdenpot
 use m_paw_init,        only : pawinit,paw_gencond
 use m_paw_onsite,      only : pawnabla_init
 use m_paw_dmft,        only : paw_dmft_type
 use m_paw_mkrho,       only : denfgr
 use m_paw_nhat,        only : nhatgrid,pawmknhat
 use m_paw_tools,       only : chkpawovlp,pawprt
 use m_paw_correlations,only : pawpuxinit
 use m_exc_build,       only : exc_build_ham
 use m_setvtr,          only : setvtr
 use m_mkrho,           only : prtrhomxmn
 use m_pspini,          only : pspini
 use m_drivexc,         only : mkdenpos

 implicit none

 private
!!***

 public :: bethe_salpeter
!!***

contains
!!***

!!****f* m_bethe_salpeter/bethe_salpeter
!! NAME
!!  bethe_salpeter
!!
!! FUNCTION
!!  Main routine to calculate dielectric properties by solving the Bethe-Salpeter equation in
!!  Frequency-Reciprocal space on a transition (electron-hole) basis set.
!!
!! INPUTS
!! acell(3)=Length scales of primitive translations (bohr)
!! codvsn=Code version
!! Dtfil<datafiles_type>=Variables related to files.
!! Dtset<dataset_type>=All input variables for this dataset.
!! Pawang<pawang_type)>=PAW angular mesh and related data.
!! Pawrad(ntypat*usepaw)<pawrad_type>=Paw radial mesh and related data.
!! Pawtab(ntypat*usepaw)<pawtab_type>=Paw tabulated starting data.
!! Psps<pseudopotential_type>=Variables related to pseudopotentials.
!!   Before entering the first time in the routine, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters bethe_salpeter, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=Dimensionless real space primitive translations.
!! xred(3,natom)=Reduced atomic coordinates.
!!
!! Input files used during the calculation.
!!  KSS        : Kohn Sham electronic structure file.
!!  SCR (SUSC) : Files containing the symmetrized epsilon^-1 or the irreducible RPA polarizability,
!!               respectively. Used to construct the screening W.
!!  GW file    : Optional file with the GW QP corrections.
!!
!! OUTPUT
!!  Output is written on the main output file and on the following external files:
!!   * _RPA_NLF_MDF: macroscopic RPA dielectric function without non-local field effects.
!!   * _GW_NLF_MDF: macroscopic RPA dielectric function without non-local field effects calculated
!!                 with GW energies or the scissors operator.
!!   * _EXC_MDF: macroscopic dielectric function with excitonic effects obtained by solving the
!!              Bethe-Salpeter problem at different level of sophistication.
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
!!      bs_parameters_free,chkpawovlp,crystal_free,denfgr,destroy_mpi_enreg
!!      double_grid_free,ebands_free,ebands_update_occ,energies_init
!!      eprenorms_free,exc_build_ham,exc_den,exc_diago_driver
!!      exc_haydock_driver,fourdp,get_gftt,getph,gsph_free,hdr_free
!!      init_distribfft_seq,initmpi_seq,kmesh_free,metric,mkdenpos,mkrdim
!!      nhatgrid,paw_an_free,paw_an_init,paw_an_nullify,paw_gencond,paw_ij_free
!!      paw_ij_init,paw_ij_nullify,pawdenpot,pawdij,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawhur_free,pawhur_init,pawinit,pawmknhat
!!      pawnabla_init,pawprt,pawpuxinit,pawpwff_free,pawpwff_init
!!      pawrhoij_alloc,pawrhoij_copy,pawrhoij_free,pawtab_get_lsize
!!      pawtab_print,print_ngfft,prtrhomxmn,pspini,rdqps,rotate_fft_mesh
!!      screen_free,screen_init,screen_nullify,setsym_ylm,setup_bse
!!      setup_bse_interp,setvtr,symdij,test_charge,timab,vcoul_free,wfd_free
!!      wfd_init,wfd_mkrho,wfd_print,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate
!!      wfd_test_ortho,wfd_wave_free,wrtout,xmpi_bcast
!!
!! SOURCE

subroutine bethe_salpeter(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred)

!Arguments ------------------------------------
!scalars
 character(len=8),intent(in) :: codvsn
 type(datafiles_type),intent(inout) :: Dtfil
 type(dataset_type),intent(inout) :: Dtset
 type(pawang_type),intent(inout) :: Pawang
 type(pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3),xred(3,Dtset%natom)
 type(pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables ------------------------------
!scalars
 integer,parameter :: tim_fourdp0=0,level=40,ipert0=0,idir0=0,cplex1=1,master=0,option1=1
 integer :: band,cplex_rhoij,spin,ik_ibz,mqmem,iwarn
 integer :: has_dijU,has_dijso,gnt_option
 integer :: ik_bz,mband
 integer :: choice
 integer :: ider
 integer :: usexcnhat,nfft_osc,mgfft_osc
 integer :: isym,izero
 integer :: optcut,optgr0,optgr1,optgr2,option,optrad,optrhoij,psp_gencond
 integer :: ngrvdw,nhatgrdim,nkxc1,nprocs,nspden_rhoij,nzlmopt,ifft
 integer :: my_rank,rhoxsp_method,comm
 integer :: mgfftf,spin_opt,which_fixed
 integer :: nscf,nbsc,nkxc,n3xccc
 integer :: nfftf,nfftf_tot,nfftot_osc,my_minb,my_maxb
 integer :: optene,moved_atm_inside,moved_rhor,initialized,istep,ierr
 real(dp) :: ucvol,drude_plsmf,ecore,ecut_eff,ecutdg_eff,norm
 real(dp) :: gsqcutc_eff,gsqcutf_eff,gsqcut_shp
 real(dp) :: compch_fft,compch_sph,gsq_osc
 real(dp) :: vxcavg
 logical :: iscompatibleFFT,is_dfpt=.false.,paw_add_onsite,call_pawinit
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname,w_fname
 type(Pawfgr_type) :: Pawfgr
 type(excfiles) :: BS_files
 type(excparam) :: BSp
 type(paw_dmft_type) :: Paw_dmft
 type(MPI_type) :: MPI_enreg_seq
 type(crystal_t) :: Cryst
 type(kmesh_t) :: Kmesh,Qmesh
 type(gsphere_t) :: Gsph_x,Gsph_c,Gsph_x_dense,Gsph_c_dense
 type(Hdr_type) :: Hdr_wfk,Hdr_bse
 type(ebands_t) :: KS_BSt,QP_BSt,KS_BSt_dense,QP_BSt_dense
 type(Energies_type) :: KS_energies
 type(vcoul_t) :: Vcp
 type(wfd_t) :: Wfd,Wfd_dense
 type(screen_t) :: W
 type(screen_info_t) :: W_info
 type(wvl_data) :: wvl
 type(kmesh_t) :: Kmesh_dense,Qmesh_dense
 type(Hdr_type) :: Hdr_wfk_dense
 type(double_grid_t) :: grid
 type(vcoul_t) :: Vcp_dense
 type(eprenorms_t) :: Epren
!arrays
 integer :: ngfft_osc(18),ngfftc(18),ngfftf(18),nrcell(3)
 integer,allocatable :: ktabr(:,:),l_size_atm(:)
 integer,allocatable :: nband(:,:),nq_spl(:),irottb(:,:)
 integer,allocatable :: qp_vbik(:,:)
 integer,allocatable :: gfft_osc(:,:)
 real(dp),parameter :: k0(3)=zero
 real(dp) :: tsec(2),gmet(3,3),gprimd(3,3),qphon(3),rmet(3,3),rprimd(3,3),eh_rcoord(3),strsxc(6)
 real(dp),allocatable :: ph1df(:,:),prev_rhor(:,:),ph1d(:,:)
 real(dp),allocatable :: ks_nhat(:,:),ks_nhatgr(:,:,:),ks_rhog(:,:),ks_rhor(:,:),qp_aerhor(:,:)
 real(dp),allocatable :: qp_rhor(:,:),qp_rhog(:,:) !,qp_vhartr(:),qp_vtrial(:,:),qp_vxc(:,:)
 real(dp),allocatable :: qp_rhor_paw(:,:),qp_rhor_n_one(:,:),qp_rhor_nt_one(:,:),qp_nhat(:,:)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grvdw(:,:),qmax(:)
 real(dp),allocatable :: vpsp(:),xccc3d(:)
 real(dp),allocatable :: ks_vhartr(:),ks_vtrial(:,:),ks_vxc(:,:)
 real(dp),allocatable :: kxc(:,:) !,qp_kxc(:,:)
 complex(dpc),allocatable :: m_ks_to_qp(:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:)
 type(Pawrhoij_type),allocatable :: prev_Pawrhoij(:) !QP_pawrhoij(:),
 type(pawpwff_t),allocatable :: Paw_pwff(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(pawhur_t),allocatable :: Hur(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:)

!************************************************************************

 DBG_ENTER('COLL')

 call timab(650,1,tsec) ! bse(Total)
 call timab(651,1,tsec) ! bse(Init1)

 comm = xmpi_world; nprocs  = xmpi_comm_size(comm); my_rank = xmpi_comm_rank(comm)

 wfk_fname = dtfil%fnamewffk

 if (nctk_try_fort_or_ncfile(wfk_fname, msg) /= 0) then
   MSG_ERROR(msg)
 end if
 call xmpi_bcast(wfk_fname, master, comm, ierr)

 write(msg,'(8a)')&
  ' Exciton: Calculation of dielectric properties by solving the Bethe-Salpeter equation ',ch10,&
  ' in frequency domain and reciprocal space on a transitions basis set. ',ch10,&
  ' Based on a program developed by L. Reining, V. Olevano, F. Sottile, ',ch10,&
  ' S. Albrecht, and G. Onida. Incorporated in ABINIT by M. Giantomassi. ',ch10
 call wrtout([std_out, ab_out], msg)

#ifdef HAVE_GW_DPC
 if (gwpc/=8) then
   write(msg,'(6a)')ch10,&
    ' Number of bytes for double precision complex /=8 ',ch10,&
    ' Cannot continue due to kind mismatch in BLAS library ',ch10,&
    ' Some BLAS interfaces are not generated by abilint '
   MSG_ERROR(msg)
 end if
 write(msg,'(a,i2,a)')'.Using double precision arithmetic ; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic ; gwpc = ',gwpc,ch10
#endif
 call wrtout([std_out, ab_out], msg)

 !=== Some variables need to be initialized/nullify at start ===
 call energies_init(KS_energies)
 usexcnhat=0
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)
 !
 !=== Define FFT grid(s) sizes ===
 !* Be careful! This mesh is only used for densities, potentials and the matrix elements of v_Hxc. It is NOT the
 !(usually coarser) GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
 !See also NOTES in the comments at the beginning of this file.
 !NOTE: This mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

!TODO Recheck getng, should use same trick as that used in screening and sigma.
 call pawfgr_init(Pawfgr,Dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
  gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)

 ! === Initialization of basic objects including the BSp structure that defines the parameters of the run ===
 call setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
   Cryst,Kmesh,Qmesh,KS_BSt,QP_BSt,Hdr_wfk,Gsph_x,Gsph_c,Vcp,Hdr_bse,w_fname,Epren,comm,wvl%descr)

 if (BSp%use_interp) then
   call setup_bse_interp(Dtset,Dtfil,BSp,Cryst,Kmesh,Kmesh_dense,&
     Qmesh_dense,KS_BSt_dense,QP_BSt_dense,Gsph_x_dense,Gsph_c_dense,&
     Vcp_dense,Hdr_wfk_dense,ngfftf,grid,comm)
 end if

 !call timab(652,2,tsec) ! setup_bse

 nfftot_osc=PRODUCT(ngfft_osc(1:3))
 nfft_osc  =nfftot_osc  !no FFT //
 mgfft_osc =MAXVAL(ngfft_osc(1:3))

 call print_ngfft(ngfft_osc,header='FFT mesh used for oscillator strengths')

 !TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol

 !
 !============================
 !==== PAW initialization ====
 !============================
 if (Dtset%usepaw==1) then
   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,xred)

   ABI_MALLOC(KS_Pawrhoij,(Cryst%natom))
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
               nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
   call pawrhoij_alloc(KS_Pawrhoij,cplex_rhoij,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)

   ! Initialize values for several basic arrays ===
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or.call_pawinit) then
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     call pawinit(Dtset%effmass_free,gnt_option,gsqcut_shp,zero,Dtset%pawlcutd,Dtset%pawlmix,&
       Psps%mpsang,Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Pawang,Pawrad,&
       Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,Dtset%ixc,Dtset%usepotzero)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit)
   else
     if (Pawtab(1)%has_kij  ==1) Pawtab(1:Cryst%ntypat)%has_kij  =2
     if (Pawtab(1)%has_nabla==1) Pawtab(1:Cryst%ntypat)%has_nabla=2
   end if
   Psps%n1xccc=MAXVAL(Pawtab(1:Cryst%ntypat)%usetcore)

   ! Initialize optional flags in Pawtab to zero
   ! (Cannot be done in Pawinit since the routine is called only if some parameters are changed)
   Pawtab(:)%has_nabla = 0
   Pawtab(:)%usepawu   = 0
   Pawtab(:)%useexexch = 0
   Pawtab(:)%exchmix   =zero

   ! * Evaluate <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for the long wavelength limit.
   ! TODO solve problem with memory leak and clean this part as well as the associated flag
   call pawnabla_init(Psps%mpsang,Cryst%ntypat,Pawrad,Pawtab)

   call setsym_ylm(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)

   ! Initialize and compute data for DFT+U
   Paw_dmft%use_dmft=Dtset%usedmft
   call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla,&
      is_dfpt,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,Cryst%ntypat,Pawang,Dtset%pawprtvol,&
      Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu)
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     MSG_ERROR('BS equation with DFT+U not completely coded!')
   end if
   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(Hdr_wfk%pawrhoij,KS_Pawrhoij)

   ! Re-symmetrize rhoij ===
   ! this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1
!  call pawrhoij_symrhoij(KS_Pawrhoij,KS_Pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert0,&
!  &             Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,Pawtab,&
!  &             Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

   ! Evaluate form factor of radial part of phi.phj-tphi.tphj ===
   rhoxsp_method=1 ! Arnaud-Alouani
   if (Dtset%pawoptosc /= 0) rhoxsp_method = Dtset%pawoptosc

   ABI_MALLOC(gfft_osc,(3,nfftot_osc))
   call get_gftt(ngfft_osc,k0,gmet,gsq_osc,gfft_osc)
   ABI_FREE(gfft_osc)

   ! Set up q grids, make qmax 20% larger than largest expected:
   ABI_MALLOC(nq_spl,(Psps%ntypat))
   ABI_MALLOC(qmax,(Psps%ntypat))
   nq_spl = Psps%mqgrid_ff
   qmax = SQRT(gsq_osc)*1.2d0 ! qmax=Psps%qgrid_ff(Psps%mqgrid_ff)
   ABI_MALLOC(Paw_pwff,(Psps%ntypat))

   call pawpwff_init(Paw_pwff,rhoxsp_method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)

   ABI_FREE(nq_spl)
   ABI_FREE(qmax)
   !
   ! Variables/arrays related to the fine FFT grid ===
   ABI_MALLOC(ks_nhat,(nfftf,Dtset%nspden))
   ks_nhat=zero
   ABI_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   call pawfgrtab_init(Pawfgrtab,cplex1,l_size_atm,Dtset%nspden,Dtset%typat)
   ABI_FREE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
   ! * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
   ! * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i0)')' bethe_salpeter : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg)
   !
   ! Identify parts of the rectangular grid where the density has to be calculated ===
   optcut=0; optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)
 else
   ABI_MALLOC(Paw_pwff,(0))
 end if !End of PAW Initialization

 ! Consistency check and additional stuff done only for GW with PAW.
 if (Dtset%usepaw==1) then
   if (Dtset%ecutwfn < Dtset%ecut) then
     write(msg,"(5a)")&
      "WARNING - ",ch10,&
      "  It is highly recommended to use ecutwfn = ecut for GW calculations with PAW since ",ch10,&
      "  an excessive truncation of the planewave basis set can lead to unphysical results."
     call wrtout(ab_out,msg,'COLL')
   end if

   ABI_CHECK(Dtset%usedmft==0,"DMFT + BSE not allowed")
   ABI_CHECK(Dtset%useexexch==0,"LEXX + BSE not allowed")
 end if

 ! Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(ks_nhat))  then
   ABI_MALLOC(ks_nhat,(nfftf,0))
 end if

!==================================================
!==== Read KS band structure from the KSS file ====
!==================================================

 ! Initialize wave function handler, allocate wavefunctions.
 my_minb=1; my_maxb=BSp%nbnds; mband=BSp%nbnds
 ABI_MALLOC(nband,(Kmesh%nibz,Dtset%nsppol))
 nband=mband

!At present, no memory distribution, each node has the full set of states.
 ABI_MALLOC(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.TRUE.

 ABI_MALLOC(keep_ur,(mband,Kmesh%nibz,Dtset%nsppol))
 keep_ur=.FALSE.; if (MODULO(Dtset%gwmem,10)==1) keep_ur = .TRUE.

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,mband,nband,Kmesh%nibz,Dtset%nsppol,bks_mask,&
  Dtset%nspden,Dtset%nspinor,Dtset%ecutwfn,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,ngfft_osc,&
  Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm)

 ABI_FREE(bks_mask)
 ABI_FREE(nband)
 ABI_FREE(keep_ur)

 call wfd%print(header="Wavefunctions used to construct the e-h basis set",mode_paral='PERS')

 call timab(651,2,tsec) ! bse(Init1)
 call timab(653,1,tsec) ! bse(rdkss)

 call wfd%read_wfk(wfk_fname,iomode_from_fname(wfk_fname))

 ! This test has been disabled (too expensive!)
 if (.False.) call wfd%test_ortho(Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(653,2,tsec) ! bse(rdkss)
 call timab(655,1,tsec) ! bse(mkrho)

 !TODO : check the consistency of Wfd with Wfd_dense !!!
 if (BSp%use_interp) then
   ! Initialize wave function handler, allocate wavefunctions.
   my_minb=1; my_maxb=BSp%nbnds; mband=BSp%nbnds
   ABI_MALLOC(nband,(Kmesh_dense%nibz,Dtset%nsppol))
   nband=mband

   ! At present, no memory distribution, each node has the full set of states.
   ! albeit we allocate only the states that are used.
   ABI_MALLOC(bks_mask,(mband,Kmesh_dense%nibz,Dtset%nsppol))
   bks_mask=.False.
   do spin=1,Bsp%nsppol
     bks_mask(Bsp%lomo_spin(spin):Bsp%humo_spin(spin),:,spin) = .True.
   end do
   !bks_mask=.TRUE.

   ABI_MALLOC(keep_ur,(mband,Kmesh_dense%nibz,Dtset%nsppol))
   keep_ur=.FALSE.; if (MODULO(Dtset%gwmem,10)==1) keep_ur = .TRUE.

   call wfd_init(Wfd_dense,Cryst,Pawtab,Psps,keep_ur,mband,nband,Kmesh_dense%nibz,Dtset%nsppol,&
    bks_mask,Dtset%nspden,Dtset%nspinor,Dtset%ecutwfn,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk_dense%istwfk,Kmesh_dense%ibz,ngfft_osc,&
    Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm)

   ABI_FREE(bks_mask)
   ABI_FREE(nband)
   ABI_FREE(keep_ur)

   call wfd_dense%print(header="Wavefunctions on the dense K-mesh used for interpolation",mode_paral='PERS')

   call wfd_dense%read_wfk(Dtfil%fnameabi_wfkfine, iomode_from_fname(dtfil%fnameabi_wfkfine))

   ! This test has been disabled (too expensive!)
   if (.False.) call wfd_dense%test_ortho(Cryst,Pawtab,unit=std_out,mode_paral="COLL")
 end if

 !=== Calculate the FFT index of $(R^{-1}(r-\tau))$ ===
 !* S=\transpose R^{-1} and k_BZ = S k_IBZ
 !* irottb is the FFT index of $R^{-1} (r-\tau)$ used to symmetrize u_Sk.
 ABI_MALLOC(irottb,(nfftot_osc,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfft_osc,irottb,iscompatibleFFT)

 ABI_MALLOC(ktabr,(nfftot_osc,Kmesh%nbz))
 do ik_bz=1,Kmesh%nbz
   isym=Kmesh%tabo(ik_bz)
   do ifft=1,nfftot_osc
     ktabr(ifft,ik_bz)=irottb(ifft,isym)
   end do
 end do
 ABI_FREE(irottb)
 !
 !===========================
 !=== COMPUTE THE DENSITY ===
 !===========================
 !* Evaluate Planewave part (complete charge in case of NC pseudos).
 !
 ABI_MALLOC(ks_rhor,(nfftf,Wfd%nspden))
 call wfd%mkrho(Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)
 !
 !=== Additional computation for PAW ===
 nhatgrdim=0
 if (Dtset%usepaw==1) then
   !
   ! Calculate the compensation charge nhat.
   if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   ider=2*nhatgrdim; izero=0; qphon(:)=zero
   if (nhatgrdim>0)  then
     ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,3))
   end if

   call pawmknhat(compch_fft,cplex1,ider,idir0,ipert0,izero,Cryst%gprimd,&
     Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
     Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,qphon,Cryst%rprimd,&
     Cryst%ucvol,Dtset%usewvl,Cryst%xred)

   ! Evaluate onsite energies, potentials, densities ===
   !   * Initialize variables/arrays related to the PAW spheres.
   !   * Initialize also lmselect (index of non-zero LM-moments of densities).
   ABI_MALLOC(KS_paw_ij,(Cryst%natom))
   call paw_ij_nullify(KS_paw_ij)

   has_dijso=Dtset%pawspnorb
   has_dijU=merge(0,1,Dtset%usepawu==0)

   call paw_ij_init(KS_paw_ij,cplex1,Dtset%nspinor,Dtset%nsppol,&
     Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
     has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=0,has_dijxc_hat=0,has_dijxc_val=0,&
     has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   ABI_MALLOC(KS_paw_an,(Cryst%natom))
   call paw_an_nullify(KS_paw_an)

   nkxc1=0
   call paw_an_init(KS_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,0,Dtset%nspden,&
     cplex1,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=0)

   ! Calculate onsite vxc with and without core charge ===
   nzlmopt=-1; option=0; compch_sph=greatest_real

   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert0,&
     Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
     Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
     Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
     Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)
 end if !PAW

 if (.not.allocated(ks_nhatgr))  then
   ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
   Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)

 ! === For PAW, add the compensation charge on the FFT mesh, then get rho(G) ===
 if (Dtset%usepaw==1) ks_rhor(:,:)=ks_rhor(:,:)+ks_nhat(:,:)
 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_rhor,ucvol=ucvol)

 ABI_MALLOC(ks_rhog,(2,nfftf))

 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,tim_fourdp0)
 call timab(655,2,tsec) ! bse(mkrho)

 !
 ! The following steps have been gathered in the setvtr routine:
 ! - get Ewald energy and Ewald forces
 ! - compute local ionic pseudopotential vpsp
 ! - eventually compute 3D core electron density xccc3d
 ! - eventually compute vxc and vhartr
 ! - set up ks_vtrial
 !
 ! *******************************************************************
 ! **** NOTE THAT HERE Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
 ! *******************************************************************
 ngrvdw=0
 ABI_MALLOC(grvdw,(3,ngrvdw))
 ABI_MALLOC(grchempottn,(3,Cryst%natom))
 ABI_MALLOC(grewtn,(3,Cryst%natom))
 nkxc=0
!if (Wfd%nspden==1) nkxc=2
!if (Wfd%nspden>=2) nkxc=3 ! check GGA and spinor, quite a messy part!!!
 ABI_MALLOC(kxc,(nfftf,nkxc))

 n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
 ABI_MALLOC(xccc3d,(n3xccc))
 ABI_MALLOC(ks_vhartr,(nfftf))
 ABI_MALLOC(ks_vtrial,(nfftf,Wfd%nspden))
 ABI_MALLOC(vpsp,(nfftf))
 ABI_MALLOC(ks_vxc,(nfftf,Wfd%nspden))

 optene=4; moved_atm_inside=0; moved_rhor=0; initialized=1; istep=1
!
!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
 if (Dtset%mgfftdg/=mgfftf) then
!  write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
!  write(std_out,*)'HACKING Dtset%mgfftf'
!  Dtset%mgfftdg=mgfftf
 end if
 ABI_MALLOC(ph1d,(2,3*(2*Dtset%mgfft+1)*Cryst%natom))
 ABI_MALLOC(ph1df,(2,3*(2*mgfftf+1)*Cryst%natom))

 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)

 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
   call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

 ABI_FREE(ph1d)

 call setvtr(Cryst%atindx1,Dtset,KS_energies,Cryst%gmet,Cryst%gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,&
& istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
& Cryst%nattyp,nfftf,ngfftf,ngrvdw,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
& optene,Pawrad,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
& Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,wvl,xccc3d,Cryst%xred)

 ABI_FREE(ph1df)
 ABI_FREE(vpsp)

 ! ============================
 ! ==== Compute KS PAW Dij ====
 ! ============================
 if (Wfd%usepaw==1) then
   _IBM6("Another silly write for IBM6")
   call timab(561,1,tsec)
   !
   ! Calculate the unsymmetrized Dij.
   call pawdij(cplex1,Dtset%enunit,Cryst%gprimd,ipert0,&
     Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
     Dtset%nspden,Cryst%ntypat,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
     Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
     k0,Dtset%spnorbscl,Cryst%ucvol,dtset%charge,ks_vtrial,ks_vxc,Cryst%xred,&
     nucdipmom=Dtset%nucdipmom)

   ! Symmetrize KS Dij
   call symdij(Cryst%gprimd,Cryst%indsym,ipert0,&
     Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,KS_paw_ij,Pawang,&
     Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)

   ! Output the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Cryst%natom,KS_paw_ij,KS_Pawrhoij,Pawtab)
   call timab(561,2,tsec)
 end if

 ABI_FREE(kxc)
 ABI_FREE(xccc3d)
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)
 ABI_FREE(grvdw)

 !=== QP_BSt stores energies and occ. used for the calculation ===
 !* Initialize QP_BSt with KS values.
 !* In case of SC update QP_BSt using the QPS file.
 ABI_MALLOC(qp_rhor,(nfftf,Dtset%nspden))
 qp_rhor=ks_rhor

 ! AE density used for the model dielectric function.
 ABI_MALLOC(qp_aerhor,(nfftf,Dtset%nspden))
 qp_aerhor=ks_rhor

 ! PAW: Compute AE rhor. Under testing
 if (Wfd%usepaw==1 .and. BSp%mdlf_type/=0) then
   MSG_WARNING("Entering qp_aerhor with PAW")

   ABI_MALLOC(qp_rhor_paw   ,(nfftf,Wfd%nspden))
   ABI_MALLOC(qp_rhor_n_one ,(nfftf,Wfd%nspden))
   ABI_MALLOC(qp_rhor_nt_one,(nfftf,Wfd%nspden))

   ABI_MALLOC(qp_nhat,(nfftf,Wfd%nspden))
   qp_nhat = ks_nhat
   ! TODO: I pass KS_pawrhoij instead of QP_pawrhoij but in the present version there's no difference.

   call denfgr(Cryst%atindx1,Cryst%gmet,Wfd%comm,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,qp_nhat,&
&   Wfd%nspinor,Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,KS_pawrhoij,Pawtab,Dtset%prtvol,&
&   qp_rhor,qp_rhor_paw,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   norm = SUM(qp_rhor_paw(:,1))*Cryst%ucvol/PRODUCT(Pawfgr%ngfft(1:3))
   write(msg,'(a,F8.4)') '  QUASIPARTICLE DENSITY CALCULATED - NORM OF DENSITY: ',norm
   call wrtout(std_out,msg)
   write(std_out,*)"MAX", MAXVAL(qp_rhor_paw(:,1))
   write(std_out,*)"MIN", MINVAL(qp_rhor_paw(:,1))

   ABI_FREE(qp_nhat)
   ABI_FREE(qp_rhor_n_one)
   ABI_FREE(qp_rhor_nt_one)

   ! Use ae density for the model dielectric function.
   iwarn=0
   call mkdenpos(iwarn,nfftf,Wfd%nspden,option1,qp_rhor_paw,dtset%xc_denpos)
   qp_aerhor = qp_rhor_paw
   ABI_FREE(qp_rhor_paw)
 end if

!call copy_bandstructure(KS_BSt,QP_BSt)

 if (.FALSE.) then
!  $ m_ks_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}> $
   ABI_MALLOC(m_ks_to_qp,(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   m_ks_to_qp=czero
   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz
       do band=1,Wfd%nband(ik_ibz,spin)
         m_ks_to_qp(band,band,ik_ibz,spin)=cone ! Initialize the QP amplitudes with KS wavefunctions.
       end do
     end do
   end do
   !
   ! Now read m_ks_to_qp and update the energies in QP_BSt.
   ! TODO switch on the renormalization of n in sigma.
   ABI_MALLOC(prev_rhor,(nfftf,Wfd%nspden))
   ABI_MALLOC(prev_Pawrhoij,(Cryst%natom*Wfd%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Wfd%usepaw,Wfd%nspden,1,nscf,&
    nfftf,ngfftf,Cryst%ucvol,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_ks_to_qp,prev_rhor,prev_Pawrhoij)

   ABI_FREE(prev_rhor)
   if (Psps%usepaw==1.and.nscf>0) then
     call pawrhoij_free(prev_pawrhoij)
   end if
   ABI_FREE(prev_pawrhoij)
   !
   !  if (nscf>0.and.wfd_iam_master(Wfd)) then ! Print the unitary transformation on std_out.
   !  call show_QP(QP_BSt,m_ks_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
   !  end if
   !
   !  === Compute QP wfg as linear combination of KS states ===
   !  * Wfd%ug is modified inside calc_wf_qp
   !  * For PAW, update also the on-site projections.
   !  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
   !  TODO here we should use nbsc instead of nbnds

   call wfd%rotate(Cryst,m_ks_to_qp)

   ABI_FREE(m_ks_to_qp)
   !
   !  === Reinit the storage mode of Wfd as ug have been changed ===
   !  * Update also the wavefunctions for GW corrections on each processor
   call wfd%reset_ur_cprj()

   call wfd%test_ortho(Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

   ! Compute QP occupation numbers.
   call wrtout(std_out,'bethe_salpeter: calculating QP occupation numbers')

   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=0)
   ABI_MALLOC(qp_vbik,(QP_BSt%nkpt,QP_BSt%nsppol))
   qp_vbik(:,:) = ebands_get_valence_idx(QP_BSt)
   ABI_FREE(qp_vbik)

   call wfd%mkrho(Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
 end if

 ABI_MALLOC(qp_rhog,(2,nfftf))
 call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,0)

 ! States up to lomo-1 are useless now since only the states close to the gap are
 ! needed to construct the EXC Hamiltonian. Here we deallocate the wavefunctions
 ! to make room for the excitonic Hamiltonian that is allocated in exc_build_ham.
 ! and for the screening that is allocated below.
 ! Hereafter bands from 1 up to lomo-1 and all bands above humo+1 should not be accessed!
 ABI_MALLOC(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))

 bks_mask=.FALSE.
 do spin=1,Bsp%nsppol
   if (Bsp%lomo_spin(spin)>1) bks_mask(1:Bsp%lomo_spin(spin)-1,:,spin)=.TRUE.
   if (Bsp%humo_spin(spin)+1<=Wfd%mband) bks_mask(Bsp%humo_spin(spin)+1:,:,spin)=.TRUE.
 end do
 call wfd%wave_free(what="All", bks_mask=bks_mask)
 ABI_FREE(bks_mask)
 !
 ! ================================================================
 ! Build the screened interaction W in the irreducible wedge.
 ! * W(q,G1,G2) = vc^{1/2} (q,G1) e^{-1}(q,G1,G2) vc^{1/2) (q,G2)
 ! * Use Coulomb term for q-->0,
 ! * Only the first small Q is used, shall we average if nqlwl>1?
 ! ================================================================
 ! TODO clean this part and add an option to retrieve a single frequency to save memory.
 call timab(654,1,tsec) ! bse(rdmkeps^-1)

 call screen_nullify(W)
 if (BSp%use_coulomb_term) then !  Init W.
   ! Incore or out-of-core solution?
   mqmem=0; if (Dtset%gwmem/10==1) mqmem=Qmesh%nibz

   W_info%invalid_freq = Dtset%gw_invalid_freq
   W_info%mat_type = MAT_INV_EPSILON
   !W_info%mat_type = MAT_W
   !W_info%vtx_family
   !W_info%ixc
   !W_info%use_ada
   W_info%use_mdf = BSp%mdlf_type
   !W_info%use_ppm
   !W_info%vtx_test
   !W_info%wint_method
   !
   !W_info%ada_kappa
   W_info%eps_inf = BSp%eps_inf
   !W_info%drude_plsmf

   call screen_init(W,W_Info,Cryst,Qmesh,Gsph_c,Vcp,w_fname,mqmem,Dtset%npweps,&
    Dtset%iomode,ngfftf,nfftf_tot,Wfd%nsppol,Wfd%nspden,qp_aerhor,Wfd%prtvol,Wfd%comm)
 end if
 call timab(654,2,tsec) ! bse(rdmkeps^-1)
 !
 ! =============================================
 ! ==== Build of the excitonic Hamiltonian =====
 ! =============================================
 !
 call timab(656,1,tsec) ! bse(mkexcham)

 call exc_build_ham(BSp,BS_files,Cryst,Kmesh,Qmesh,ktabr,Gsph_x,Gsph_c,Vcp,&
& Wfd,W,Hdr_bse,nfftot_osc,ngfft_osc,Psps,Pawtab,Pawang,Paw_pwff)
 !
 ! Free Em1 to make room for the full excitonic Hamiltonian.
 call screen_free(W)

 call timab(656,2,tsec) ! bse(mkexcham)
 !
 ! =========================================
 ! ==== Macroscopic dielectric function ====
 ! =========================================
 call timab(657,1,tsec) ! bse(mkexceps)
 !
 ! First deallocate the internal %ur buffers to make room for the excitonic Hamiltonian.
 call timab(658,1,tsec) ! bse(wfd_wave_free)
 ABI_MALLOC(bks_mask,(Wfd%mband,Wfd%nkibz,Wfd%nsppol))
 bks_mask=.TRUE.
 call wfd%wave_free(what="Real_space", bks_mask=bks_mask)
 ABI_FREE(bks_mask)
 call timab(658,2,tsec) ! bse(wfd_wave_free)

 ! Compute the commutator [r,Vu] (PAW only).
 ABI_MALLOC(HUr,(Cryst%natom*Wfd%usepaw))

 call timab(659,1,tsec) ! bse(make_pawhur_t)
 if (Bsp%inclvkb/=0 .and. Wfd%usepaw==1 .and. Dtset%usepawu/=0) then !TODO here I need KS_Paw_ij
   MSG_WARNING("Commutator for DFT+U not tested")
   call pawhur_init(hur,Wfd%nsppol,Wfd%pawprtvol,Cryst,Pawtab,Pawang,Pawrad,KS_Paw_ij)
 end if
 call timab(659,2,tsec) ! bse(make_pawhur_t)

 select case (BSp%algorithm)
 case (BSE_ALGO_NONE)
   MSG_COMMENT("Skipping solution of the BSE equation")

 case (BSE_ALGO_DDIAGO, BSE_ALGO_CG)
   call timab(660,1,tsec) ! bse(exc_diago_driver)
   call exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps, Pawtab,Hur,Hdr_bse,drude_plsmf,Epren)
   call timab(660,2,tsec) ! bse(exc_diago_driver)

   if (.FALSE.) then ! Calculate electron-hole excited state density. Not tested at all.
     call exc_den(BSp,BS_files,ngfftf,nfftf,Kmesh,ktabr,Wfd)
   end if

   if (.FALSE.) then
     paw_add_onsite=.FALSE.; spin_opt=1; which_fixed=1; eh_rcoord=(/zero,zero,zero/); nrcell=(/2,2,2/)
!    call exc_plot(Bsp,Bs_files,Wfd,Kmesh,Cryst,Psps,Pawtab,Pawrad,paw_add_onsite,spin_opt,which_fixed,eh_rcoord,nrcell,ngfftf)
   end if

   if (BSp%use_interp) then
     MSG_ERROR("Interpolation technique not coded for diagonalization and CG")
   end if

 case (BSE_ALGO_Haydock)
   call timab(661,1,tsec) ! bse(exc_haydock_driver)

   if (BSp%use_interp) then

     call exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_BSt,Wfd,Psps,Pawtab,Hur,Epren,&
       Kmesh_dense,KS_BSt_dense,QP_BSt_dense,Wfd_dense,Vcp_dense,grid)
   else
     call exc_haydock_driver(BSp,BS_files,Cryst,Kmesh,Hdr_bse,KS_BSt,QP_BSt,Wfd,Psps,Pawtab,Hur,Epren)
   end if

   call timab(661,2,tsec) ! bse(exc_haydock_driver)

 case default
   write(msg,'(a,i0)')" Wrong BSE algorithm: ",BSp%algorithm
   MSG_ERROR(msg)
 end select

 call timab(657,2,tsec) ! bse(mkexceps)

 !=====================
 !==== Free memory ====
 !=====================
 ABI_FREE(ktabr)
 ABI_FREE(ks_vhartr)
 ABI_FREE(ks_vtrial)
 ABI_FREE(ks_vxc)
 ABI_FREE(ks_nhat)
 ABI_FREE(ks_nhatgr)
 ABI_FREE(ks_rhog)
 ABI_FREE(ks_rhor)
 ABI_FREE(qp_rhog)
 ABI_FREE(qp_rhor)
 ABI_FREE(qp_aerhor)
 !
 ! Free local data structures.
 call destroy_mpi_enreg(MPI_enreg_seq)
 call cryst%free()
 call gsph_free(Gsph_x)
 call gsph_free(Gsph_c)
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call Hdr_wfk%free()
 call Hdr_bse%free()
 call ebands_free(KS_BSt)
 call ebands_free(QP_BSt)
 call vcoul_free(Vcp)
 call pawhur_free(Hur)
 ABI_FREE(Hur)
 call bs_parameters_free(BSp)
 call wfd%free()
 call pawfgr_destroy(Pawfgr)
 call eprenorms_free(Epren)

 ! Free memory used for interpolation.
 if (BSp%use_interp) then
   call double_grid_free(grid)
   call wfd_dense%free()
   call gsph_free(Gsph_x_dense)
   call gsph_free(Gsph_c_dense)
   call kmesh_free(Kmesh_dense)
   call kmesh_free(Qmesh_dense)
   call ebands_free(KS_BSt_dense)
   call ebands_free(QP_BSt_dense)
   call vcoul_free(Vcp_dense)
   call Hdr_wfk_dense%free()
 end if

 ! Optional deallocation for PAW.
 if (Dtset%usepaw==1) then
   call pawrhoij_free(KS_Pawrhoij)
   ABI_FREE(KS_Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   ABI_FREE(Pawfgrtab)
   call paw_ij_free(KS_paw_ij)
   ABI_FREE(KS_paw_ij)
   call paw_an_free(KS_paw_an)
   ABI_FREE(KS_paw_an)
   call pawpwff_free(Paw_pwff)
 end if
 ABI_FREE(Paw_pwff)

 call timab(650,2,tsec) ! bse(Total)

 DBG_EXIT('COLL')

end subroutine bethe_salpeter
!!***

!!****f* m_bethe_salpeter/setup_bse
!! NAME
!!  setup_bse
!!
!! FUNCTION
!!  This routine performs the initialization of basic objects and quantities used for Bethe-Salpeter calculations.
!!  In particular the excparam data type that defines the parameters of the calculation is completely
!!  initialized starting from the content of Dtset and the parameters read from the external WFK and SCR (SUSC) file.
!!
!! INPUTS
!! codvsn=Code version
!! ngfft_gw(18)=Information about 3D FFT for density and potentials, see ~abinit/doc/variables/vargs.htm#ngfft
!! acell(3)=Length scales of primitive translations (bohr)
!! rprim(3,3)=Dimensionless real space primitive translations.
!! Dtset<dataset_type>=All input variables for this dataset.
!!  Some of them might be redefined here TODO
!! Dtfil=filenames and unit numbers used in abinit.
!! Psps <pseudopotential_type>=variables related to pseudopotentials
!! Pawtab(Psps%ntypat*Dtset%usepaw)<pawtab_type>=PAW tabulated starting data
!!
!! OUTPUT
!! Cryst<crystal_t>=Info on the crystalline Structure.
!! Kmesh<kmesh_t>=Structure defining the k-sampling for the wavefunctions.
!! Qmesh<kmesh_t>=Structure defining the q-sampling for the symmetrized inverse dielectric matrix.
!! Gsph_x<gsphere_t=Data type gathering info on the G-sphere for wave functions and e^{-1},
!! KS_BSt<ebands_t>=The KS band structure (energies, occupancies, k-weights...)
!! Vcp<vcoul_t>=Structure gathering information on the Coulomb interaction in reciprocal space,
!!   including a possible cutoff in real space.
!! ngfft_osc(18)=Contain all needed information about the 3D FFT for the oscillator matrix elements.
!!   See ~abinit/doc/variables/vargs.htm#ngfft
!! Bsp<excparam>=Basic parameters defining the Bethe-Salpeter run. Completely initialed in output.
!! Hdr_wfk<Hdr_type>=The header of the WFK file.
!! Hdr_bse<Hdr_type>=Local header initialized from the parameters used for the Bethe-Salpeter calculation.
!! BS_files<excfiles>=Files used in the calculation.
!! w_file=File name used to construct W. Set to ABI_NOFILE if no external file is used.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      ebands_apply_scissors,bsp_calctype2str,crystal_from_hdr,crystal_print
!!      ebands_copy,ebands_init,ebands_print,ebands_report_gap
!!      ebands_update_occ,eprenorms_bcast,eprenorms_from_epnc,find_qmesh
!!      get_bz_item,get_ng0sh,gsph_extend,gsph_init,hdr_init,hdr_update
!!      hdr_vs_dtset,hscr_bcast,hscr_free,hscr_from_file,hscr_print
!!      init_transitions,kmesh_init,kmesh_print,make_mesh,matrginv,metric
!!      mkrdim,pawrhoij_alloc,pawrhoij_copy,pawrhoij_free,print_bs_files
!!      print_bs_parameters,print_gsphere,print_ngfft,rdgw,setmesh,vcoul_init
!!      wfk_read_eigenvalues,wrtout,xmpi_bcast,xmpi_max,xmpi_split_work
!!
!! SOURCE

subroutine setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
& Cryst,Kmesh,Qmesh,KS_BSt,QP_bst,Hdr_wfk,Gsph_x,Gsph_c,Vcp,Hdr_bse,w_fname,Epren,comm,Wvl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=8),intent(in) :: codvsn
 character(len=fnlen),intent(out) :: w_fname
 type(dataset_type),intent(inout) :: Dtset
 type(datafiles_type),intent(in) :: Dtfil
 type(pseudopotential_type),intent(in) :: Psps
 type(excparam),intent(inout) :: Bsp
 type(hdr_type),intent(out) :: Hdr_wfk,Hdr_bse
 type(crystal_t),intent(out) :: Cryst
 type(kmesh_t),intent(out) :: Kmesh,Qmesh
 type(gsphere_t),intent(out) :: Gsph_x,Gsph_c
 type(ebands_t),intent(out) :: KS_BSt,QP_Bst
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(vcoul_t),intent(out) :: Vcp
 type(excfiles),intent(out) :: BS_files
 type(wvl_internal_type), intent(in) :: Wvl
 type(eprenorms_t),intent(out) :: Epren
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: ngfft_osc(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables ------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer(i8b) :: work_size,tot_nreh,neh_per_proc,il
 integer :: bantot,enforce_sym,ib,ibtot,ik_ibz,isppol,jj,method,iat,ount !ii,
 integer :: mband,io,nfftot_osc,spin,hexc_size,nqlwl,iq
 integer :: timrev,iq_bz,isym,iq_ibz,itim
 integer :: my_rank,nprocs,fform,npwe_file,ierr,my_k1, my_k2,my_nbks
 integer :: first_dig,second_dig,it
 real(dp) :: ucvol,qnorm
 real(dp):: eff,mempercpu_mb,wfsmem_mb,nonscal_mem,ug_mem,ur_mem,cprj_mem
 logical,parameter :: remove_inv=.FALSE.
 logical :: ltest,occ_from_dtset
 character(len=500) :: msg
 character(len=fnlen) :: gw_fname,test_file,wfk_fname
 character(len=fnlen) :: ep_nc_fname
 type(hscr_t) :: Hscr
!arrays
 integer :: ng0sh_opt(3),val_idx(Dtset%nsppol)
 integer,allocatable :: npwarr(:),val_indeces(:,:),nlmn_atm(:)
 real(dp) :: qpt_bz(3),minmax_tene(2)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3)
 real(dp) :: qred2cart(3,3),qcart2red(3,3)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 real(dp),allocatable :: igwene(:,:,:)
 real(dp),pointer :: energies_p(:,:,:)
 complex(dpc),allocatable :: gw_energy(:,:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)

!************************************************************************

 DBG_ENTER("COLL")

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 ! === Check for calculations that are not implemented ===
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband must be constant')
 ABI_CHECK(Dtset%nspinor==1,"nspinor==2 not coded")

 ! === Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume ===
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 ! Read energies and header from the WFK file.
 wfk_fname = dtfil%fnamewffk
 if (.not. file_exists(wfk_fname)) then
   wfk_fname = nctk_ncify(wfk_fname)
   MSG_COMMENT(sjoin("File not found. Will try netcdf file: ", wfk_fname))
 end if

 call wfk_read_eigenvalues(wfk_fname,energies_p,Hdr_wfk,comm)
 mband = MAXVAL(Hdr_wfk%nband)

 call hdr_wfk%vs_dtset(dtset)

 ! === Create crystal_t data type ===
 !remove_inv= .FALSE. !(nsym_kss/=Hdr_wfk%nsym)
 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry
            ! 2 => take advantage of time-reversal symmetry

 cryst = Hdr_wfk%get_crystal(gw_timrev=timrev, remove_inv=remove_inv)
 call cryst%print()
 !
 ! Setup of the k-point list and symmetry tables in the  BZ -----------------------------------
 if (Dtset%chksymbreak==0) then
   MSG_WARNING("Calling make_mesh")
   call make_mesh(Kmesh,Cryst,Dtset%kptopt,Dtset%kptrlatt,Dtset%nshiftk,Dtset%shiftk,break_symmetry=.TRUE.)
   ! TODO
   !Check if kibz from KSS file corresponds to the one returned by make_mesh.
 else
   call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt)
 end if
 BSp%nkibz = Kmesh%nibz  !We might allow for a smaller number of points....

 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 nqlwl = 0
 w_fname = ABI_NOFILE
 if (Dtset%getscr/=0 .or. Dtset%irdscr/=0 .or. dtset%getscr_filepath /= ABI_NOFILE) then
   w_fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   w_fname=Dtfil%fnameabi_sus
   MSG_ERROR("(get|ird)suscep not implemented")
 end if

 if (w_fname /= ABI_NOFILE) then ! Read dimensions from the external file

   if (.not. file_exists(w_fname)) then
     w_fname = nctk_ncify(w_fname)
     MSG_COMMENT(sjoin("File not found. Will try netcdf file: ", w_fname))
   end if

   if (my_rank==master) then
     ! Master reads npw and nqlwl from SCR file.
     call wrtout(std_out,sjoin('Testing file: ', w_fname))

     call hscr_from_file(hscr,w_fname,fform,xmpi_comm_self)
     ! Echo the header.
     if (Dtset%prtvol>0) call hscr_print(Hscr)

     npwe_file = Hscr%npwe ! Have to change %npweps if it was larger than dim on disk.
     nqlwl     = Hscr%nqlwl

     if (Dtset%npweps>npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
        "The number of G-vectors stored on file (",npwe_file,") is smaller than Dtset%npweps = ",Dtset%npweps,ch10,&
        "Calculation will proceed with the maximum available set, npwe_file = ",npwe_file
       MSG_WARNING(msg)
       Dtset%npweps = npwe_file
     else  if (Dtset%npweps<npwe_file) then
       write(msg,'(2(a,i0),2a,i0)')&
        "The number of G-vectors stored on file (",npwe_file,") is larger than Dtset%npweps = ",Dtset%npweps,ch10,&
        "Calculation will proceed with Dtset%npweps = ",Dtset%npweps
       MSG_COMMENT(msg)
     end if
   end if

   call hscr_bcast(Hscr,master,my_rank,comm)
   call xmpi_bcast(Dtset%npweps,master,comm,ierr)
   call xmpi_bcast(nqlwl,master,comm,ierr)

   if (nqlwl>0) then
     ABI_MALLOC(qlwl,(3,nqlwl))
     qlwl = Hscr%qlwl
   end if
   !
   ! Init Qmesh from the SCR file.
   call kmesh_init(Qmesh,Cryst,Hscr%nqibz,Hscr%qibz,Dtset%kptopt)

   ! The G-sphere for W and Sigma_c is initialized from the gvectors found in the SCR file.
   call gsph_init(Gsph_c,Cryst,Dtset%npweps,gvec=Hscr%gvec)

   call hscr_free(Hscr)
 else
   ! Init Qmesh from the K-mesh reported in the WFK file.
   call find_qmesh(Qmesh,Cryst,Kmesh)

   ! The G-sphere for W and Sigma_c is initialized from ecuteps.
   call gsph_init(Gsph_c,Cryst,0,ecut=Dtset%ecuteps)
   Dtset%npweps = Gsph_c%ng
 end if

 BSp%npweps = Dtset%npweps
 BSp%ecuteps = Dtset%ecuteps

 if (nqlwl==0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,nqlwl)= GW_Q0_DEFAULT
   write(msg,'(3a,i2,a,3f9.6)')&
     "The Header of the screening file does not contain the list of q-point for the optical limit ",ch10,&
     "Using nqlwl= ",nqlwl," and qlwl = ",qlwl(:,1)
   MSG_COMMENT(msg)
 end if
 write(std_out,*)"nqlwl and qlwl for Coulomb singularity and e^-1",nqlwl,qlwl

 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed,
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 call kmesh_print(Qmesh,"Q-mesh for the screening function",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Qmesh,"Q-mesh for the screening function",ab_out ,0           ,"COLL")

 do iq_bz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iq_bz,qpt_bz,iq_ibz,isym,itim)
   sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
   if (ANY(ABS(Qmesh%bz(:,iq_bz)-sq )>1.0d-4)) then
     write(std_out,*) sq,Qmesh%bz(:,iq_bz)
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
      'qpoint ',Qmesh%bz(:,iq_bz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
      'through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
      'however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do

 BSp%algorithm = Dtset%bs_algorithm
 BSp%nstates   = Dtset%bs_nstates
 Bsp%nsppol    = Dtset%nsppol
 Bsp%hayd_term = Dtset%bs_hayd_term

 ! Define the algorithm for solving the BSE.
 if (BSp%algorithm == BSE_ALGO_HAYDOCK) then
   BSp%niter       = Dtset%bs_haydock_niter
   BSp%haydock_tol = Dtset%bs_haydock_tol

 else if (BSp%algorithm == BSE_ALGO_CG) then
   ! FIXME For the time being use an hardcoded value.
   ! TODO change name in Dtset%
   BSp%niter       = Dtset%nstep !100
   BSp%cg_tolwfr   = Dtset%tolwfr
   BSp%nline       = Dtset%nline
   BSp%nbdbuf      = Dtset%nbdbuf
   BSp%nstates     = Dtset%bs_nstates
   MSG_WARNING("Check CG setup")
 else
   !BSp%niter       = 0
   !BSp%tol_iter    = HUGE(one)
 end if
 !
 ! Shall we include Local field effects?
 SELECT CASE (Dtset%bs_exchange_term)
 CASE (0,1)
   BSp%exchange_term = Dtset%bs_exchange_term
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong bs_exchange_term: ",Dtset%bs_exchange_term
   MSG_ERROR(msg)
 END SELECT
 !
 ! Treatment of the off-diagonal coupling block.
 SELECT CASE (Dtset%bs_coupling)
 CASE (0)
   BSp%use_coupling = 0
   msg = 'RESONANT ONLY CALCULATION'
 CASE (1)
   BSp%use_coupling = 1
   msg = ' RESONANT+COUPLING CALCULATION '
 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong bs_coupling: ",Dtset%bs_coupling
   MSG_ERROR(msg)
 END SELECT
 call wrtout(std_out,msg)

 BSp%use_diagonal_Wgg = .FALSE.
 Bsp%use_coulomb_term = .TRUE.
 BSp%eps_inf=zero
 Bsp%mdlf_type=0

 first_dig =MOD(Dtset%bs_coulomb_term,10)
 second_dig=Dtset%bs_coulomb_term/10

 Bsp%wtype = second_dig
 SELECT CASE (second_dig)
 CASE (BSE_WTYPE_NONE)
   call wrtout(std_out,"Coulomb term won't be calculated")
   Bsp%use_coulomb_term = .FALSE.

 CASE (BSE_WTYPE_FROM_SCR)
   call wrtout(std_out,"W is read from an external SCR file")
   Bsp%use_coulomb_term = .TRUE.

 CASE (BSE_WTYPE_FROM_MDL)
   call wrtout(std_out,"W is approximated with the model dielectric function")
   Bsp%use_coulomb_term = .TRUE.
   BSp%mdlf_type = MDL_BECHSTEDT
   BSp%eps_inf = Dtset%mdf_epsinf
   ABI_CHECK(Bsp%eps_inf > zero, "mdf_epsinf <= 0")

 CASE DEFAULT
   write(msg,'(a,i0)')" Wrong second digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
   MSG_ERROR(msg)
 END SELECT
 !
 ! Diagonal approximation or full matrix?
 BSp%use_diagonal_Wgg = .TRUE.
 if (Bsp%wtype /= BSE_WTYPE_NONE) then
   SELECT CASE (first_dig)
   CASE (0)
     call wrtout(std_out,"Using diagonal approximation W_GG")
     BSp%use_diagonal_Wgg = .TRUE.
   CASE (1)
     call wrtout(std_out,"Using full W_GG' matrix ")
     BSp%use_diagonal_Wgg = .FALSE.
   CASE DEFAULT
     write(msg,'(a,i0)')" Wrong first digit in bs_coulomb_term: ",Dtset%bs_coulomb_term
     MSG_ERROR(msg)
   END SELECT
 end if

 !TODO move the initialization of the parameters for the interpolation in setup_bse_interp

 BSp%use_interp = .FALSE.
 BSp%interp_mode = BSE_INTERP_YG
 BSp%interp_kmult(1:3) = 0
 BSp%prep_interp = .FALSE.
 BSp%sum_overlaps = .TRUE. ! Sum over the overlaps

 ! Printing ncham
 BSp%prt_ncham = .FALSE.

 ! Deactivate Interpolation Technique by default
! if (.FALSE.) then

 ! Reading parameters from the input file
 BSp%use_interp = (dtset%bs_interp_mode /= 0)
 BSp%prep_interp = (dtset%bs_interp_prep == 1)

 SELECT CASE (dtset%bs_interp_mode)
 CASE (0)
   ! No interpolation, do not print anything !
 CASE (1)
   call wrtout(std_out,"Using interpolation technique with energies and wavefunctions from dense WFK")
 CASE (2)
   call wrtout(std_out,"Interpolation technique with energies and wfn on dense WFK + treatment ABC of divergence")
 CASE (3)
   call wrtout(std_out,"Interpolation technique + divergence ABC along diagonal")
 CASE (4)
   call wrtout(std_out,"Using interpolation technique mode 1 with full computation of hamiltonian")
 CASE DEFAULT
   MSG_ERROR(sjoin("Wrong interpolation mode for bs_interp_mode:", itoa(dtset%bs_interp_mode)))
 END SELECT

 ! Read from dtset
 if(BSp%use_interp) then
   BSp%interp_method = dtset%bs_interp_method
   BSp%rl_nb = dtset%bs_interp_rl_nb
   BSp%interp_m3_width = dtset%bs_interp_m3_width
   BSp%interp_kmult(1:3) = dtset%bs_interp_kmult(1:3)
   BSp%interp_mode = dtset%bs_interp_mode
 end if

 ! Dimensions and parameters of the calculation.
 ! TODO one should add npwx as well
 !BSp%npweps=Dtset%npweps
 !BSp%npwwfn=Dtset%npwwfn

 ABI_MALLOC(Bsp%lomo_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%homo_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%lumo_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%humo_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%nbndv_spin, (Bsp%nsppol))
 ABI_MALLOC(Bsp%nbndc_spin, (Bsp%nsppol))

 ! FIXME use bs_loband(nsppol)
 Bsp%lomo_spin = Dtset%bs_loband
 write(std_out,*)"bs_loband",Dtset%bs_loband
 !if (Bsp%nsppol == 2) Bsp%lomo_spin(2) = Dtset%bs_loband

 ! Check lomo correct only for unpolarized semiconductors
 !if (Dtset%nsppol == 1 .and. Bsp%lomo > Dtset%nelect/2) then
 !  write(msg,'(a,i0,a,f8.3)') " Bsp%lomo = ",Bsp%lomo," cannot be greater than nelect/2 = ",Dtset%nelect/2
 !  MSG_ERROR(msg)
 !end if
 !
 ! ==============================================
 ! ==== Setup of the q for the optical limit ====
 ! ==============================================
 Bsp%inclvkb = Dtset%inclvkb

 qred2cart = two_pi*Cryst%gprimd
 qcart2red = qred2cart
 call matrginv(qcart2red,3,3)

 if (Dtset%gw_nqlwl==0) then
   BSp%nq = 6
   ABI_MALLOC(BSp%q,(3,BSp%nq))
   BSp%q(:,1) = (/one,zero,zero/)  ! (100)
   BSp%q(:,2) = (/zero,one,zero/)  ! (010)
   BSp%q(:,3) = (/zero,zero,one/)  ! (001)
   BSp%q(:,4) = MATMUL(qcart2red,(/one,zero,zero/)) ! (x)
   BSp%q(:,5) = MATMUL(qcart2red,(/zero,one,zero/)) ! (y)
   BSp%q(:,6) = MATMUL(qcart2red,(/zero,zero,one/)) ! (z)
 else
   BSp%nq = Dtset%gw_nqlwl
   ABI_MALLOC(BSp%q,(3,BSp%nq))
   BSp%q = Dtset%gw_qlwl
 end if

 do iq=1,BSp%nq ! normalization
   qnorm = normv(BSp%q(:,iq),Cryst%gmet,"G")
   BSp%q(:,iq) = BSp%q(:,iq)/qnorm
 end do
 !
 ! ======================================================
 ! === Define the flags defining the calculation type ===
 ! ======================================================
 Bsp%calc_type = Dtset%bs_calctype

 BSp%mbpt_sciss = zero ! Shall we use the scissors operator to open the gap?
 if (ABS(Dtset%mbpt_sciss)>tol6) BSp%mbpt_sciss = Dtset%mbpt_sciss

!now test input parameters from input and WFK file and assume some defaults
!
! TODO Add the possibility of using a randomly shifted k-mesh with nsym>1.
! so that densities and potentials are correctly symmetrized but
! the list of the k-point in the IBZ is not expanded.

 if (mband < Dtset%nband(1)) then
   write(msg,'(2(a,i0),3a,i0)')&
    'WFK file contains only ', mband,' levels instead of ',Dtset%nband(1),' required;',ch10,&
    'The calculation will be done with nbands= ',mband
   MSG_WARNING(msg)
   Dtset%nband(:) = mband
 end if

 BSp%nbnds = Dtset%nband(1) ! TODO Note the change in the meaning of input variables

 if (BSp%nbnds<=Dtset%nelect/2) then
   write(msg,'(2a,a,i0,a,f8.2)')&
    'BSp%nbnds cannot be smaller than homo ',ch10,&
    'while BSp%nbnds = ',BSp%nbnds,' and Dtset%nelect = ',Dtset%nelect
   MSG_ERROR(msg)
 end if

!TODO add new dim for exchange part and consider the possibility of having npwsigx > npwwfn (see setup_sigma).

 ! === Build enlarged G-sphere for the exchange part ===
 call gsph_extend(Gsph_c,Cryst,Dtset%ecutwfn,Gsph_x)
 call print_gsphere(Gsph_x,unit=std_out,prtvol=Dtset%prtvol)

 ! NPWVEC as the biggest between npweps and npwwfn. MG RECHECK this part.
 !BSp%npwwfn = Dtset%npwwfn
 Bsp%npwwfn = Gsph_x%ng  ! FIXME temporary hack
 BSp%npwvec=MAX(BSp%npwwfn,BSp%npweps)
 Bsp%ecutwfn = Dtset%ecutwfn

 ! Compute Coulomb term on the largest G-sphere.
 if (Gsph_x%ng > Gsph_c%ng ) then
   call vcoul_init(Vcp,Gsph_x,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Dtset%ecutsigx,Gsph_x%ng,&
     nqlwl,qlwl,ngfftf,comm)
 else
   call vcoul_init(Vcp,Gsph_c,Cryst,Qmesh,Kmesh,Dtset%rcut,Dtset%icutcoul,Dtset%vcutgeo,Dtset%ecutsigx,Gsph_c%ng,&
     nqlwl,qlwl,ngfftf,comm)
 end if

 ABI_FREE(qlwl)

 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))
 ABI_MALLOC(doccde,(bantot))
 ABI_MALLOC(eigen,(bantot))
 ABI_MALLOC(occfact,(bantot))
 doccde=zero; eigen=zero; occfact=zero

 ! Get occupation from input if occopt == 2
 occ_from_dtset = (Dtset%occopt == 2)

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ik_ibz=1,Dtset%nkpt
     do ib=1,Hdr_wfk%nband(ik_ibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=BSP%nbnds) then
         jj=jj+1
         eigen  (jj)=energies_p(ib,ik_ibz,isppol)
         if (occ_from_dtset) then
           !Not occupations must be the same for different images
           occfact(jj)=Dtset%occ_orig(ibtot,1)
         else
           occfact(jj)=Hdr_wfk%occ(ibtot)
         end if
       end if
     end do
   end do
 end do

 ABI_FREE(energies_p)
 !
 ! Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of
 ! symmetry operations in the old GW code (symmorphy and inversion)
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_MALLOC(npwarr,(Dtset%nkpt))
 npwarr=BSP%npwwfn

 call ebands_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt,&
&  dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
&  dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 !TODO Occupancies are zero if NSCF. One should calculate the occupancies from the energies when
 ! the occupation scheme for semiconductors is used.
 call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=Dtset%prtvol)

 call ebands_print(KS_BSt,"Band structure read from the WFK file",unit=std_out,prtvol=Dtset%prtvol)

 call ebands_report_gap(KS_BSt,header=" KS band structure",unit=std_out,mode_paral="COLL")

 ABI_MALLOC(val_indeces,(KS_BSt%nkpt,KS_BSt%nsppol))
 val_indeces = ebands_get_valence_idx(KS_BSt)

 do spin=1,KS_BSt%nsppol
   val_idx(spin) = val_indeces(1,spin)
   write(msg,'(a,i2,a,i0)')" For spin : ",spin," val_idx ",val_idx(spin)
   call wrtout(std_out,msg)
   if ( ANY(val_indeces(1,spin) /= val_indeces(:,spin)) ) then
     MSG_ERROR("BSE code does not support metals")
   end if
 end do

 ABI_FREE(val_indeces)
 !
 ! === Create the BSE header ===
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_bse,Pawtab,pertcase0,Psps,wvl)

 ! === Get Pawrhoij from the header of the WFK file ===
 ABI_MALLOC(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   call pawrhoij_alloc(Pawrhoij,1,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)
 end if

 call hdr_bse%update(bantot,1.0d20,1.0d20,1.0d20,Cryst%rprimd,occfact,Pawrhoij,Cryst%xred,dtset%amu_orig(:,1))

 ABI_FREE(occfact)

 if (Dtset%usepaw==1) call pawrhoij_free(Pawrhoij)
 ABI_FREE(Pawrhoij)

 ! Find optimal value for G-sphere enlargment due to oscillator matrix elements
 ! We will split k-points over processors
 call xmpi_split_work(Kmesh%nbz, comm, my_k1, my_k2)

 ! If there is no work to do, just skip the computation
 if (my_k2-my_k1+1 <= 0) then
   ng0sh_opt(:)=(/zero,zero,zero/)
 else
   ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true since bz is buggy
   ! * -one is used because we loop over all the possibile differences, unlike screening
   call get_ng0sh(my_k2-my_k1+1,Kmesh%bz(:,my_k1:my_k2),Kmesh%nbz,Kmesh%bz,&
&    Qmesh%nbz,Qmesh%bz,-one,ng0sh_opt)
 end if

 call xmpi_max(ng0sh_opt,BSp%mg0,comm,ierr)

 write(msg,'(a,3(i0,1x))') ' optimal value for ng0sh = ',BSp%mg0
 call wrtout(std_out,msg)

 ! === Setup of the FFT mesh for the oscilator strengths ===
 ! * ngfft_osc(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine ngfft_osc(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 ngfft_osc(1:18)=Dtset%ngfft(1:18); method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10)

 call setmesh(gmet,Gsph_x%gvec,ngfft_osc,BSp%npwvec,BSp%npweps,BSp%npwwfn,nfftot_osc,method,BSp%mg0,Cryst,enforce_sym)
 nfftot_osc=PRODUCT(ngfft_osc(1:3))

 call print_ngfft(ngfft_osc,"FFT mesh for oscillator matrix elements",std_out,"COLL",prtvol=Dtset%prtvol)
 !
 ! BSp%homo gives the
 !BSp%homo  = val_idx(1)
 ! highest occupied band for each spin
 BSp%homo_spin = val_idx

 ! TODO generalize the code to account for this unlikely case.
 !if (Dtset%nsppol==2) then
 !  ABI_CHECK(BSp%homo == val_idx(2),"Different valence indeces for spin up and down")
 !end if

 !BSp%lumo = BSp%homo + 1
 !BSp%humo = BSp%nbnds
 !BSp%nbndv = BSp%homo  - BSp%lomo + 1
 !BSp%nbndc = BSp%nbnds - BSp%homo

 BSp%lumo_spin = BSp%homo_spin + 1
 BSp%humo_spin = BSp%nbnds
 BSp%nbndv_spin = BSp%homo_spin  - BSp%lomo_spin + 1
 BSp%nbndc_spin = BSp%nbnds - BSp%homo_spin
 BSp%maxnbndv = MAXVAL(BSp%nbndv_spin(:))
 BSp%maxnbndc = MAXVAL(BSp%nbndc_spin(:))

 BSp%nkbz = Kmesh%nbz

 call ebands_copy(KS_BSt,QP_bst)
 ABI_MALLOC(igwene,(QP_bst%mband,QP_bst%nkpt,QP_bst%nsppol))
 igwene=zero

 call bsp_calctype2str(Bsp,msg)
 call wrtout(std_out,"Calculation type: "//TRIM(msg))

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator energy= ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the KS energies."
     call wrtout(std_out,msg)
     call ebands_apply_scissors(QP_BSt,BSp%mbpt_sciss)
   else
     write(msg,'(a,f8.2,a)')' Using KS energies since mbpt_sciss= ',BSp%mbpt_sciss*Ha_eV," [eV]."
     call wrtout(std_out,msg)
   end if

 CASE (BSE_HTYPE_RPA_QPENE) ! Read _GW files with the corrections TODO here I should introduce variable getgw
   gw_fname=TRIM(Dtfil%filnam_ds(4))//'_GW'
   gw_fname="__in.gw__"
   if (.not.file_exists(gw_fname)) then
     msg = " File "//TRIM(gw_fname)//" not found. Aborting now"
     MSG_ERROR(msg)
   end if

   call rdgw(QP_Bst,gw_fname,igwene,extrapolate=.FALSE.) ! here gwenergy is real

   do isppol=1,Dtset%nsppol
     write(std_out,*) ' k       GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(QP_bst%eig(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
     write(std_out,*) ' k       Im GW energies [eV]'
     do ik_ibz=1,Kmesh%nibz
       write(std_out,'(i3,7x,10f7.2/50(10x,10f7.2/))')ik_ibz,(igwene(ib,ik_ibz,isppol)*Ha_eV,ib=1,BSp%nbnds)
     end do
   end do
   !
   ! If required apply the scissors operator on top of the QP bands structure (!)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the QP energies!"
     MSG_COMMENT(msg)
     call ebands_apply_scissors(QP_BSt,BSp%mbpt_sciss)
   end if

 CASE (BSE_HTYPE_RPA_QP)
   MSG_ERROR("Not implemented error!")

 CASE DEFAULT
   write(msg,'(a,i0)')"Unknown value for Bsp%calc_type= ",Bsp%calc_type
   MSG_ERROR(msg)
 END SELECT

 call ebands_report_gap(QP_BSt,header=" QP band structure",unit=std_out,mode_paral="COLL")

 ! Transitions are ALWAYS ordered in c-v-k mode with k being the slowest index.
 ! FIXME: linewidths not coded.
 ABI_MALLOC(gw_energy,(BSp%nbnds,Kmesh%nibz,Dtset%nsppol))
 gw_energy = QP_BSt%eig

 BSp%have_complex_ene = ANY(igwene > tol16)

 ! Compute the number of resonant transitions, nreh, for the two spin channels and initialize BSp%Trans.
 ABI_MALLOC(Bsp%nreh,(Bsp%nsppol))

 ! Possible cutoff on the transitions.
 BSp%ircut = Dtset%bs_eh_cutoff(1)
 BSp%uvcut = Dtset%bs_eh_cutoff(2)

 call init_transitions(BSp%Trans,BSp%lomo_spin,BSp%humo_spin,BSp%ircut,Bsp%uvcut,BSp%nkbz,Bsp%nbnds,Bsp%nkibz,&
&  BSp%nsppol,Dtset%nspinor,gw_energy,QP_BSt%occ,Kmesh%tab,minmax_tene,Bsp%nreh)

 ! Setup of the frequency mesh for the absorption spectrum.
 ! If not specified, use the min-max resonant transition energy and make it 10% smaller|larger.

 !if (ABS(Dtset%bs_freq_mesh(1)) < tol6) then
 !   Dtset%bs_freq_mesh(1) = MAX(minmax_tene(1) - minmax_tene(1) * 0.1, zero)
 !end if

 if (ABS(Dtset%bs_freq_mesh(2)) < tol6) then
    Dtset%bs_freq_mesh(2) = minmax_tene(2) + minmax_tene(2) * 0.1
 end if

 Bsp%omegai = Dtset%bs_freq_mesh(1)
 Bsp%omegae = Dtset%bs_freq_mesh(2)
 Bsp%domega = Dtset%bs_freq_mesh(3)
 BSp%broad  = Dtset%zcut

 ! The frequency mesh (including the complex imaginary shift)
 BSp%nomega = (BSp%omegae - BSp%omegai)/BSp%domega + 1
 ABI_MALLOC(BSp%omega,(BSp%nomega))
 do io=1,BSp%nomega
   BSp%omega(io) = (BSp%omegai + (io-1)*BSp%domega)  + j_dpc*BSp%broad
 end do

 ABI_FREE(gw_energy)
 ABI_FREE(igwene)

 do spin=1,Bsp%nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',BSp%nreh(spin)
   call wrtout(std_out,msg)
 end do

 if (ANY(Bsp%nreh/=Bsp%nreh(1))) then
   write(msg,'(a,2(i0,1x))')" BSE code with different number of transitions for the two spin channels: ",Bsp%nreh
   MSG_WARNING(msg)
 end if
 !
 ! Create transition table vcks2t
 Bsp%lomo_min = MINVAL(BSp%lomo_spin)
 Bsp%homo_max = MAXVAL(BSp%homo_spin)
 Bsp%lumo_min = MINVAL(BSp%lumo_spin)
 Bsp%humo_max = MAXVAL(BSp%humo_spin)

 ABI_MALLOC(Bsp%vcks2t,(BSp%lomo_min:BSp%homo_max,BSp%lumo_min:BSp%humo_max,BSp%nkbz,Dtset%nsppol))
 Bsp%vcks2t = 0

 do spin=1,BSp%nsppol
   do it=1,BSp%nreh(spin)
     BSp%vcks2t(BSp%Trans(it,spin)%v,BSp%Trans(it,spin)%c,BSp%Trans(it,spin)%k,spin) = it
   end do
 end do

 hexc_size = SUM(Bsp%nreh); if (Bsp%use_coupling>0) hexc_size=2*hexc_size
 if (Bsp%nstates<=0) then
   Bsp%nstates=hexc_size
 else
   if (Bsp%nstates>hexc_size) then
      Bsp%nstates=hexc_size
      write(msg,'(2(a,i0),2a)')&
       "Since the total size of excitonic Hamiltonian ",hexc_size," is smaller than Bsp%nstates ",Bsp%nstates,ch10,&
       "the number of excitonic states nstates has been modified"
     MSG_WARNING(msg)
   end if
 end if

 msg=' Fundamental parameters for the solution of the Bethe-Salpeter equation:'
 call print_bs_parameters(BSp,unit=std_out,header=msg,mode_paral="COLL",prtvol=Dtset%prtvol)
 call print_bs_parameters(BSp,unit=ab_out, header=msg,mode_paral="COLL")

 if (ANY (Cryst%symrec(:,:,1) /= RESHAPE ( (/1,0,0,0,1,0,0,0,1/),(/3,3/) )) .or. ANY( ABS(Cryst%tnons(:,1)) > tol6) ) then
   write(msg,'(3a,9i2,2a,3f6.3,2a)')&
     "The first symmetry operation should be the Identity with zero tnons while ",ch10,&
     "symrec(:,:,1) = ",Cryst%symrec(:,:,1),ch10,&
     "tnons(:,1)    = ",Cryst%tnons(:,1),ch10,&
     "This is not allowed, sym_rhotwgq0 should be changed."
   MSG_ERROR(msg)
 end if
 !
 ! Prefix for generic output files.
 BS_files%out_basename = TRIM(Dtfil%filnam_ds(4))
 !
 ! Search for files to restart from.
 if (Dtset%gethaydock/=0 .or. Dtset%irdhaydock/=0) then
   BS_files%in_haydock_basename = TRIM(Dtfil%fnameabi_haydock)
 end if

 test_file = Dtfil%fnameabi_bsham_reso
 if (file_exists(test_file)) then
   BS_files%in_hreso = test_file
 else
   BS_files%out_hreso = TRIM(Dtfil%filnam_ds(4))//'_BSR'
 end if

 test_file = Dtfil%fnameabi_bsham_coup
 if (file_exists(test_file) ) then
   BS_files%in_hcoup = test_file
 else
   BS_files%out_hcoup = TRIM(Dtfil%filnam_ds(4))//'_BSC'
 end if
 !
 ! in_eig is the name of the input file with eigenvalues and eigenvectors
 ! constructed from getbseig or irdbseig. out_eig is the name of the output file
 ! produced by this dataset. in_eig_exists checks for the presence of the input file.
 !
 if (file_exists(Dtfil%fnameabi_bseig)) then
   BS_files%in_eig = Dtfil%fnameabi_bseig
 else
   BS_files%out_eig = TRIM(BS_files%out_basename)//"_BSEIG"
 end if

 call print_bs_files(BS_files,unit=std_out)
 !
 ! ==========================================================
 ! ==== Temperature dependence of the spectrum ==============
 ! ==========================================================
 BSp%do_ep_renorm = .FALSE.
 BSp%do_lifetime = .FALSE. ! Not yet implemented

 ep_nc_fname = 'test_EP.nc'
 if(file_exists(ep_nc_fname)) then
   BSp%do_ep_renorm = .TRUE.

   if(my_rank == master) then
     call eprenorms_from_epnc(Epren,ep_nc_fname)
   end if
   call eprenorms_bcast(Epren,master,comm)
 end if
 !
 ! ==========================================================
 ! ==== Final check on the parameters of the calculation ====
 ! ==========================================================
 if ( Bsp%use_coupling>0 .and. ALL(Bsp%algorithm /= [BSE_ALGO_DDIAGO, BSE_ALGO_HAYDOCK]) ) then
   MSG_ERROR("Resonant+Coupling is only available with the direct diagonalization or the haydock method.")
 end if

 ! autoparal section
 if (dtset%max_ncpus /=0 .and. dtset%autoparal /=0 ) then
   ount = ab_out
   ! TODO:
   ! nsppol and calculation with coupling!

   ! Temporary table needed to estimate memory
   ABI_MALLOC(nlmn_atm,(Cryst%natom))
   if (Dtset%usepaw==1) then
     do iat=1,Cryst%natom
       nlmn_atm(iat)=Pawtab(Cryst%typat(iat))%lmn_size
     end do
   end if

   tot_nreh = SUM(BSp%nreh)
   work_size = tot_nreh * (tot_nreh + 1) / 2

   write(ount,'(a)')"--- !Autoparal"
   write(ount,"(a)")'#Autoparal section for Bethe-Salpeter runs.'

   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ount,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ount,"(a,i0)")"    nkibz: ",Bsp%nkibz
   write(ount,"(a,i0)")"    nkbz: ",Bsp%nkbz
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    lomo_min: ",Bsp%lomo_min
   write(ount,"(a,i0)")"    humo_max: ",Bsp%humo_max
   write(ount,"(a,i0)")"    tot_nreh: ",tot_nreh
   !write(ount,"(a,i0)")"    nbnds: ",Ep%nbnds

   ! Wavefunctions are not distributed. We read all the bands
   ! from 1 up to Bsp%nbnds because we have to recompute rhor
   ! but then we deallocate all the states that are not used for the construction of the e-h
   ! before allocating the EXC hamiltonian. Hence we can safely use  (humo - lomo + 1) instead of Bsp%nbnds.
   !my_nbks = (Bsp%humo - Bsp%lomo +1) * Bsp%nkibz * Dtset%nsppol

   ! This one overestimates the memory but it seems to be safer.
   my_nbks = Bsp%nbnds * Dtset%nkpt * Dtset%nsppol

   ! Memory needed for Fourier components ug.
   ug_mem = two*gwpc*Dtset%nspinor*Bsp%npwwfn*my_nbks*b2Mb

   ! Memory needed for real space ur.
   ur_mem = zero
   if (MODULO(Dtset%gwmem,10)==1) then
     ur_mem = two*gwpc*Dtset%nspinor*nfftot_osc*my_nbks*b2Mb
   end if

   ! Memory needed for PAW projections Cprj
   cprj_mem = zero
   if (Dtset%usepaw==1) cprj_mem = dp*Dtset%nspinor*SUM(nlmn_atm)*my_nbks*b2Mb

   wfsmem_mb = ug_mem + ur_mem + cprj_mem

   ! Non-scalable memory in Mb i.e. memory that is not distributed with MPI:  wavefunctions + W
   nonscal_mem = (wfsmem_mb + two*gwpc*BSp%npweps**2*b2Mb) * 1.1_dp

   ! List of configurations.
   write(ount,"(a)")"configurations:"
   do il=1,dtset%max_ncpus
     if (il > work_size) cycle
     neh_per_proc = work_size / il
     neh_per_proc = neh_per_proc + MOD(work_size, il)
     eff = (one * work_size) / (il * neh_per_proc)

     ! EXC matrix is distributed.
     mempercpu_mb = nonscal_mem + two * dpc * neh_per_proc * b2Mb

     write(ount,"(a,i0)")"    - tot_ncpus: ",il
     write(ount,"(a,i0)")"      mpi_ncpus: ",il
     !write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
     write(ount,"(a,f12.9)")"      efficiency: ",eff
     write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
   end do

   write(ount,'(a)')"..."

   ABI_FREE(nlmn_atm)
   MSG_ERROR_NODUMP("aborting now")
 end if

 DBG_EXIT("COLL")

end subroutine setup_bse
!!***

!!****f* m_bethe_salpeter/setup_bse_interp
!! NAME
!!  setup_bse_interp
!!
!! FUNCTION
!!
!! INPUTS
!! ngfft_gw(18)=Information about 3D FFT for density and potentials, see ~abinit/doc/variables/vargs.htm#ngfft
!! acell(3)=Length scales of primitive translations (bohr)
!! rprim(3,3)=Dimensionless real space primitive translations.
!! Dtset<dataset_type>=All input variables for this dataset.
!!  Some of them might be redefined here TODO
!! Dtfil=filenames and unit numbers used in abinit. fnameabi_wfkfile is changed is Fortran file is not
!! found but a netcdf version with similar name is available.
!!
!! OUTPUT
!! Cryst<crystal_structure>=Info on the crystalline Structure.
!! Kmesh<BZ_mesh_type>=Structure defining the k-sampling for the wavefunctions.
!! Qmesh<BZ_mesh_type>=Structure defining the q-sampling for the symmetrized inverse dielectric matrix.
!! Gsph_x<gsphere_t=Data type gathering info on the G-sphere for wave functions and e^{-1},
!! KS_BSt<Bandstructure_type>=The KS band structure (energies, occupancies, k-weights...)
!! Vcp<vcoul_t>=Structure gathering information on the Coulomb interaction in reciprocal space,
!!   including a possible cutoff in real space.
!! ngfft_osc(18)=Contain all needed information about the 3D FFT for the oscillator matrix elements.
!!   See ~abinit/doc/variables/vargs.htm#ngfft
!! Bsp<excparam>=Basic parameters defining the Bethe-Salpeter run. Completely initialed in output.
!! Hdr_wfk<Hdr_type>=The header of the WFK file.
!! Hdr_bse<Hdr_type>=Local header initialized from the parameters used for the Bethe-Salpeter calculation.
!! w_file=File name used to construct W. Set to ABI_NOFILE if no external file is used.
!!
!! PARENTS
!!      bethe_salpeter
!!
!! CHILDREN
!!      ebands_apply_scissors,double_grid_init,ebands_copy,ebands_init,ebands_print
!!      ebands_report_gap,ebands_update_occ,find_qmesh,gsph_extend,gsph_init
!!      init_transitions,kmesh_init,kmesh_print,make_mesh,print_gsphere
!!      vcoul_init,wfk_read_eigenvalues,wrtout
!!
!! SOURCE

subroutine setup_bse_interp(Dtset,Dtfil,BSp,Cryst,Kmesh,&
& Kmesh_dense,Qmesh_dense,KS_BSt_dense,QP_bst_dense,Gsph_x,Gsph_c,Vcp_dense,Hdr_wfk_dense,ngfftf,grid,comm)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 type(dataset_type),intent(in) :: Dtset
 type(datafiles_type),intent(inout) :: Dtfil
 type(excparam),intent(inout) :: Bsp
 type(hdr_type),intent(out) :: Hdr_wfk_dense
 type(crystal_t),intent(in) :: Cryst
 type(kmesh_t),intent(in) :: Kmesh
 type(kmesh_t),intent(out) :: Kmesh_dense,Qmesh_dense
 type(ebands_t),intent(out) :: KS_BSt_dense,QP_Bst_dense
 type(double_grid_t),intent(out) :: grid
 type(vcoul_t),intent(out) :: Vcp_dense
 type(gsphere_t),intent(out) :: Gsph_x,Gsph_c
!arrays
 integer,intent(in) :: ngfftf(18)

!Local variables ------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer :: bantot_dense,ib,ibtot,ik_ibz,isppol,jj
 integer :: nbnds_kss_dense
 integer :: spin,hexc_size
 integer :: my_rank
 integer :: it
 integer :: nprocs
 integer :: is1,is2,is3,is4
 real(dp) :: nelect_hdr_dense
 logical,parameter :: remove_inv=.FALSE.
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname_dense
 integer :: nqlwl
!arrays
 integer,allocatable :: npwarr(:)
 real(dp),allocatable :: shiftk(:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:)
 real(dp),pointer :: energies_p_dense(:,:,:)
 complex(dpc),allocatable :: gw_energy(:,:,:)
 integer,allocatable :: nbands_temp(:)
 integer :: kptrlatt_dense(3,3)
 real(dp),allocatable :: qlwl(:,:)
 real(dp) :: minmax_tene(2)

!************************************************************************

 DBG_ENTER("COLL")

 kptrlatt_dense = zero

 my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)

 SELECT CASE(BSp%interp_mode)
 CASE (1,2,3,4)
   nbnds_kss_dense = -1
   wfk_fname_dense = Dtfil%fnameabi_wfkfine
   call wrtout(std_out,"BSE Interpolation: will read energies from: "//trim(wfk_fname_dense),"COLL")

   if (nctk_try_fort_or_ncfile(wfk_fname_dense, msg) /= 0) then
     MSG_ERROR(msg)
   end if

   Dtfil%fnameabi_wfkfine = wfk_fname_dense

   call wfk_read_eigenvalues(wfk_fname_dense,energies_p_dense,Hdr_wfk_dense,comm)
   nbnds_kss_dense = MAXVAL(Hdr_wfk_dense%nband)
 CASE DEFAULT
   MSG_ERROR("Not yet implemented")
 END SELECT

 nelect_hdr_dense = Hdr_wfk_dense%nelect

 if (ABS(Dtset%nelect-nelect_hdr_dense)>tol6) then
   write(msg,'(2(a,f8.2))')&
&   "File contains ", nelect_hdr_dense," electrons but nelect initialized from input is ",Dtset%nelect
   MSG_ERROR(msg)
 end if

 ! Setup of the k-point list and symmetry tables in the  BZ
 SELECT CASE(BSp%interp_mode)
 CASE (1,2,3,4)
   if(Dtset%chksymbreak == 0) then
     ABI_MALLOC(shiftk,(3,Dtset%nshiftk))
     kptrlatt_dense(:,1) = BSp%interp_kmult(1)*Dtset%kptrlatt(:,1)
     kptrlatt_dense(:,2) = BSp%interp_kmult(2)*Dtset%kptrlatt(:,2)
     kptrlatt_dense(:,3) = BSp%interp_kmult(3)*Dtset%kptrlatt(:,3)
     do jj = 1,Dtset%nshiftk
       shiftk(:,jj) = Bsp%interp_kmult(:)*Dtset%shiftk(:,jj)
     end do
     call make_mesh(Kmesh_dense,Cryst,Dtset%kptopt,kptrlatt_dense,Dtset%nshiftk,shiftk,break_symmetry=.TRUE.)
     ABI_FREE(shiftk)
   else
     !Initialize Kmesh with no wrapping inside ]-0.5;0.5]
     call kmesh_init(Kmesh_dense,Cryst,Hdr_wfk_dense%nkpt,Hdr_wfk_dense%kptns,Dtset%kptopt)
   end if
 CASE DEFAULT
   MSG_ERROR("Not yet implemented")
 END SELECT

 ! Init Qmesh
 call find_qmesh(Qmesh_dense,Cryst,Kmesh_dense)

 call gsph_init(Gsph_c,Cryst,0,ecut=Dtset%ecuteps)

 call double_grid_init(Kmesh,Kmesh_dense,Dtset%kptrlatt,BSp%interp_kmult,grid)

 BSp%nkibz_interp = Kmesh_dense%nibz  !We might allow for a smaller number of points....

 call kmesh_print(Kmesh_dense,"Interpolated K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh_dense,"Interpolated K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 if (nbnds_kss_dense < Dtset%nband(1)) then
   write(msg,'(2(a,i0),3a,i0)')&
    'Interpolated WFK file contains only ', nbnds_kss_dense,' levels instead of ',Dtset%nband(1),' required;',ch10,&
    'The calculation will be done with nbands= ',nbnds_kss_dense
   MSG_WARNING(msg)
   MSG_ERROR("Not supported yet !")
 end if

 ABI_MALLOC(nbands_temp,(Hdr_wfk_dense%nkpt*Hdr_wfk_dense%nsppol))
 do isppol=1,Hdr_wfk_dense%nsppol
   do ik_ibz=1,Hdr_wfk_dense%nkpt
     nbands_temp(ik_ibz+(isppol-1)*Hdr_wfk_dense%nkpt) = Dtset%nband(1)
   end do
 end do

 call gsph_extend(Gsph_c,Cryst,Dtset%ecutwfn,Gsph_x)
 call print_gsphere(Gsph_x,unit=std_out,prtvol=Dtset%prtvol)

 nqlwl=1
 ABI_MALLOC(qlwl,(3,nqlwl))
 qlwl(:,nqlwl)= GW_Q0_DEFAULT

 ! Compute Coulomb term on the largest G-sphere.
 if (Gsph_x%ng > Gsph_c%ng ) then
   call vcoul_init(Vcp_dense,Gsph_x,Cryst,Qmesh_dense,Kmesh_dense,Dtset%rcut,Dtset%icsing,&
&    Dtset%vcutgeo,Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
 else
   call vcoul_init(Vcp_dense,Gsph_c,Cryst,Qmesh_dense,Kmesh_dense,Dtset%rcut,Dtset%icsing,&
&    Dtset%vcutgeo,Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
 end if

 ABI_FREE(qlwl)

 bantot_dense=SUM(Hdr_wfk_dense%nband(1:Hdr_wfk_dense%nkpt*Hdr_wfk_dense%nsppol))
 ABI_MALLOC(doccde,(bantot_dense))
 ABI_MALLOC(eigen,(bantot_dense))
 ABI_MALLOC(occfact,(bantot_dense))
 doccde=zero; eigen=zero; occfact=zero

 jj=0; ibtot=0
 do isppol=1,Hdr_wfk_dense%nsppol
   do ik_ibz=1,Hdr_wfk_dense%nkpt
     do ib=1,Hdr_wfk_dense%nband(ik_ibz+(isppol-1)*Hdr_wfk_dense%nkpt)
       ibtot=ibtot+1
       if (ib<=BSP%nbnds) then
         jj=jj+1
         occfact(jj)=Hdr_wfk_dense%occ(ibtot)
         eigen  (jj)=energies_p_dense(ib,ik_ibz,isppol)
       end if
     end do
   end do
 end do

 ABI_FREE(energies_p_dense)

 ABI_MALLOC(npwarr,(kmesh_dense%nibz))
 npwarr=BSP%npwwfn

 call ebands_init(bantot_dense,KS_BSt_dense,Dtset%nelect,doccde,eigen,Hdr_wfk_dense%istwfk,Kmesh_dense%ibz,nbands_temp,&
&  Kmesh_dense%nibz,npwarr,Hdr_wfk_dense%nsppol,Hdr_wfk_dense%nspinor,Hdr_wfk_dense%tphysel,Hdr_wfk_dense%tsmear,&
&  Hdr_wfk_dense%occopt,occfact,Kmesh_dense%wt,&
&  hdr_wfk_dense%charge, hdr_wfk_dense%kptopt, hdr_wfk_dense%kptrlatt_orig, hdr_wfk_dense%nshiftk_orig, &
&  hdr_wfk_dense%shiftk_orig, hdr_wfk_dense%kptrlatt, hdr_wfk_dense%nshiftk, hdr_wfk_dense%shiftk)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 ABI_FREE(nbands_temp)

 ABI_FREE(occfact)

 !TODO Occupancies are zero if NSCF. One should calculate the occupancies from the energies when
 ! the occupation scheme for semiconductors is used.
 call ebands_update_occ(KS_BSt_dense,Dtset%spinmagntarget,prtvol=Dtset%prtvol)

 call ebands_print(KS_BSt_dense,"Interpolated band structure read from the WFK file",unit=std_out,prtvol=Dtset%prtvol)

 call ebands_report_gap(KS_BSt_dense,header="Interpolated KS band structure",unit=std_out,mode_paral="COLL")

 BSp%nkbz_interp = Kmesh_dense%nbz

 call ebands_copy(KS_BSt_dense,QP_bst_dense)

 SELECT CASE (Bsp%calc_type)
 CASE (BSE_HTYPE_RPA_KS)
   if (ABS(BSp%mbpt_sciss)>tol6) then
     write(msg,'(a,f8.2,a)')' Applying a scissors operator energy= ',BSp%mbpt_sciss*Ha_eV," [eV] on top of the KS energies."
     call wrtout(std_out,msg)
     call ebands_apply_scissors(QP_BSt_dense,BSp%mbpt_sciss)
   else
     write(msg,'(a,f8.2,a)')' Using KS energies since mbpt_sciss= ',BSp%mbpt_sciss*Ha_eV," [eV]."
     call wrtout(std_out,msg)
   end if
   !
 CASE (BSE_HTYPE_RPA_QPENE) ! Read _GW files with the corrections TODO here I should introduce variable getgw
   MSG_ERROR("Not yet implemented with interpolation !")
 CASE (BSE_HTYPE_RPA_QP)
   MSG_ERROR("Not implemented error!")
 CASE DEFAULT
   write(msg,'(a,i0)')"Unknown value for Bsp%calc_type= ",Bsp%calc_type
   MSG_ERROR(msg)
 END SELECT

 call ebands_report_gap(QP_BSt_dense,header=" Interpolated QP band structure",unit=std_out,mode_paral="COLL")

 ! Transitions are ALWAYS ordered in c-v-k mode with k being the slowest index.
 ! FIXME: linewidths not coded.
 ABI_MALLOC(gw_energy, (BSp%nbnds,Kmesh_dense%nibz,Dtset%nsppol))
 gw_energy = QP_BSt_dense%eig

 ABI_MALLOC(Bsp%nreh_interp,(Hdr_wfk_dense%nsppol))
 Bsp%nreh_interp=zero

 call init_transitions(BSp%Trans_interp,BSp%lomo_spin,BSp%humo_spin,BSp%ircut,Bsp%uvcut,BSp%nkbz_interp,Bsp%nbnds,&
&  Bsp%nkibz_interp,Hdr_wfk_dense%nsppol,Hdr_wfk_dense%nspinor,gw_energy,QP_BSt_dense%occ,Kmesh_dense%tab,minmax_tene,&
&  Bsp%nreh_interp)

 ABI_FREE(gw_energy)

 do spin=1,Dtset%nsppol
   write(msg,'(a,i2,a,i0)')" For spin: ",spin,' the number of resonant e-h transitions is: ',BSp%nreh_interp(spin)
   call wrtout(std_out,msg)
 end do

 if (ANY(Bsp%nreh_interp/=Bsp%nreh_interp(1))) then
   write(msg,'(a,(i0))')" BSE code does not support different number of transitions for the two spin channels",Bsp%nreh
   MSG_ERROR(msg)
 end if
 !
 ! Create transition table vcks2t
 is1=BSp%lomo_min;is2=BSp%homo_max;is3=BSp%lumo_min;is4=BSp%humo_max
 ABI_MALLOC(Bsp%vcks2t_interp, (is1:is2,is3:is4,BSp%nkbz_interp,Dtset%nsppol))
 Bsp%vcks2t_interp = 0

 do spin=1,Dtset%nsppol
   do it=1,BSp%nreh_interp(spin)
     BSp%vcks2t_interp(BSp%Trans_interp(it,spin)%v,BSp%Trans_interp(it,spin)%c, BSp%Trans_interp(it,spin)%k,spin) = it
   end do
 end do

 hexc_size = SUM(Bsp%nreh_interp); if (Bsp%use_coupling>0) hexc_size=2*hexc_size
 if (Bsp%nstates_interp<=0) then
   Bsp%nstates_interp=hexc_size
 else
   if (Bsp%nstates_interp>hexc_size) then
      Bsp%nstates_interp=hexc_size
      write(msg,'(2(a,i0),2a)')&
       "Since the total size of excitonic Hamiltonian ",hexc_size," is smaller than Bsp%nstates ",Bsp%nstates_interp,ch10,&
       "the number of excitonic states nstates has been modified"
     MSG_WARNING(msg)
   end if
 end if

 DBG_EXIT("COLL")

end subroutine setup_bse_interp
!!***

end module m_bethe_salpeter
!!***

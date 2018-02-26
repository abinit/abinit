!{\src2tex{textfont=tt}}
!!****f* ABINIT/bethe_salpeter
!! NAME
!!  bethe_salpeter
!!
!! FUNCTION
!!  Main routine to calculate dielectric properties by solving the Bethe-Salpeter equation in
!!  Frequency-Reciprocal space on a transition (electron-hole) basis set.
!!
!! COPYRIGHT
!! Copyright (C) 2009-2018 ABINIT group (M.Giantomassi, L. Reining, V. Olevano, F. Sottile, S. Albrecht, G. Onida)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      screen_free,screen_init,screen_nullify,setsymrhoij,setup_bse
!!      setup_bse_interp,setvtr,symdij,test_charge,timab,vcoul_free,wfd_free
!!      wfd_init,wfd_mkrho,wfd_print,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate
!!      wfd_test_ortho,wfd_wave_free,wrtout,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine bethe_salpeter(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_bs_defs
 use m_profiling_abi
 use m_xmpi
 use m_errors
 use m_screen
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr

 use m_fstrings,        only : strcat, sjoin, endswith
 use m_io_tools,        only : file_exists, iomode_from_fname
 use m_mpinfo,          only : destroy_mpi_enreg
 use m_fftcore,         only : print_ngfft
 use m_fft_mesh,        only : rotate_FFT_mesh, get_gftt
 use m_crystal,         only : crystal_t, crystal_free
 use m_crystal_io,      only : crystal_ncwrite
 use m_bz_mesh,         only : kmesh_t, kmesh_free, make_path
 use m_double_grid,     only : double_grid_t, double_grid_free
 use m_ebands,          only : ebands_update_occ, ebands_free, ebands_copy, get_valence_idx
 use m_gsphere,         only : gsphere_t, gsph_free
 use m_vcoul,           only : vcoul_t, vcoul_free
 use m_qparticles,      only : rdqps !, show_QP , rdgw
 use m_wfd,             only : wfd_init, wfd_free, wfd_print, wfd_t, wfd_test_ortho,&
&                              wfd_read_wfk, wfd_wave_free, wfd_rotate, wfd_reset_ur_cprj
 use m_energies,        only : energies_type, energies_init
 use m_haydock,         only : exc_haydock_driver
 use m_exc_diago,       only : exc_diago_driver
 use m_eprenorms,       only : eprenorms_t, eprenorms_free
 use m_pawang,          only : pawang_type
 use m_pawrad,          only : pawrad_type
 use m_pawtab,          only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,          only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,          only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,       only : pawfgrtab_type, pawfgrtab_free, pawfgrtab_init
 use m_pawrhoij,        only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy,&
&                              pawrhoij_free, pawrhoij_get_nspden, symrhoij
 use m_pawdij,          only : pawdij, symdij
 use m_pawfgr,          only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_pawhr,           only : pawhur_t, pawhur_free, pawhur_init
 use m_pawpwij,        only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_paw_dmft,        only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'bethe_salpeter'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_28_numeric_noabirule
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_64_psp
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_69_wfdesc
 use interfaces_71_bse
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
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
 integer :: band,spin,ik_ibz,mqmem,iwarn !ii,
 integer :: has_dijU,has_dijso,gnt_option,iomode
 integer :: ik_bz,mband
 integer :: choice
 integer :: ider !,ierr
 integer :: usexcnhat,nfft_osc,mgfft_osc
 integer :: isym,izero
 integer :: optcut,optgr0,optgr1,optgr2,option,optrad,optrhoij,psp_gencond
 integer :: ngrvdw,nhatgrdim,nkxc1,nprocs,nspden_rhoij,nzlmopt,ifft
 integer :: my_rank,rhoxsp_method,comm
 integer :: mgfftf,spin_opt,which_fixed
 integer :: nscf,nbsc,nkxc,n3xccc
 integer :: nfftf,nfftf_tot,nfftot_osc,my_minb,my_maxb
 integer :: optene,moved_atm_inside,moved_rhor,initialized,istep,ierr
 real(dp) :: ucvol,drude_plsmf,ecore,ecut_eff,ecutdg_eff,opt_ecut,norm
 real(dp) :: gsqcutc_eff,gsqcutf_eff,gsqcut_shp
 real(dp) :: compch_fft,compch_sph,gsq_osc
 real(dp) :: vxcavg 
 logical :: iscompatibleFFT,paw_add_onsite,call_pawinit
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname,w_fname
 type(Pawfgr_type) :: Pawfgr
 type(excfiles) :: BS_files
 type(excparam) :: BSp
 type(paw_dmft_type) :: Paw_dmft
 type(MPI_type) :: MPI_enreg_seq
 type(crystal_t) :: Cryst
 type(kmesh_t) :: Kmesh,Qmesh
 type(gsphere_t) :: Gsph_x,Gsph_c
 type(gsphere_t) :: Gsph_x_dense,Gsph_c_dense
 type(Hdr_type) :: Hdr_wfk,Hdr_bse
 type(ebands_t) :: KS_BSt,QP_BSt
 type(Energies_type) :: KS_energies
 type(vcoul_t) :: Vcp
 type(wfd_t) :: Wfd
 type(screen_t) :: W
 type(screen_info_t) :: W_info
 type(wvl_data) :: wvl
 type(kmesh_t) :: Kmesh_dense,Qmesh_dense
 type(Hdr_type) :: Hdr_wfk_dense
 type(ebands_t) :: KS_BSt_dense,QP_BSt_dense
 type(wfd_t) :: Wfd_dense
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
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
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
& ' Exciton: Calculation of dielectric properties by solving the Bethe-Salpeter equation ',ch10,&
& ' in frequency domain and reciprocal space on a transitions basis set. ',ch10,&
& ' Based on a program developed by L. Reining, V. Olevano, F. Sottile, ',ch10,&
& ' S. Albrecht, and G. Onida. Incorporated in ABINIT by M. Giantomassi. ',ch10
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

#ifdef HAVE_GW_DPC
 if (gwpc/=8) then
   write(msg,'(6a)')ch10,&
&   ' Number of bytes for double precision complex /=8 ',ch10,&
&   ' Cannot continue due to kind mismatch in BLAS library ',ch10,&
&   ' Some BLAS interfaces are not generated by abilint '
   MSG_ERROR(msg)
 end if
 write(msg,'(a,i2,a)')'.Using double precision arithmetic ; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic ; gwpc = ',gwpc,ch10
#endif
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

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
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))
 !
 ! * Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')
 !
 ! ===========================================
 ! === Open and read pseudopotential files ===
 ! ===========================================
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)

 ! === Initialization of basic objects including the BSp structure that defines the parameters of the run ===
 call setup_bse(codvsn,acell,rprim,ngfftf,ngfft_osc,Dtset,Dtfil,BS_files,Psps,Pawtab,BSp,&
& Cryst,Kmesh,Qmesh,KS_BSt,QP_BSt,Hdr_wfk,Gsph_x,Gsph_c,Vcp,Hdr_bse,w_fname,Epren,comm,wvl%descr)

 if (BSp%use_interp) then
   call setup_bse_interp(Dtset,Dtfil,BSp,Cryst,Kmesh,Kmesh_dense,&
&   Qmesh_dense,KS_BSt_dense,QP_BSt_dense,Gsph_x_dense,Gsph_c_dense,&
&   Vcp_dense,Hdr_wfk_dense,ngfftf,grid,comm)
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

   ABI_DT_MALLOC(KS_Pawrhoij,(Cryst%natom))
   nspden_rhoij=pawrhoij_get_nspden(Dtset%nspden,Dtset%nspinor,Dtset%pawspnorb)
   call pawrhoij_alloc(KS_Pawrhoij,Dtset%pawcpxocc,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)

   ! Initialize values for several basic arrays ===
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit) 

   if (psp_gencond==1.or.call_pawinit) then
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     call pawinit(gnt_option,gsqcut_shp,zero,Dtset%pawlcutd,Dtset%pawlmix,&
&     Psps%mpsang,Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Pawang,Pawrad,&
&     Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,Dtset%xclevel,Dtset%usepotzero)

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

   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)

   ! Initialize and compute data for LDA+U
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     Paw_dmft%use_dmft=Dtset%usedmft
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla,&
&     Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,Cryst%ntypat,Pawang,Dtset%pawprtvol,&
&     Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu)
     MSG_ERROR("BS equation with LDA+U not completely coded")
   end if
   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(Hdr_wfk%pawrhoij,KS_Pawrhoij)

   ! Re-symmetrize symrhoij ===
   ! this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1
!  call symrhoij(KS_Pawrhoij,KS_Pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert0,&
!  &  Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,Pawtab,&
!  &  Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

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
   ABI_DT_MALLOC(Paw_pwff,(Psps%ntypat))

   call pawpwff_init(Paw_pwff,rhoxsp_method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)

   ABI_FREE(nq_spl)
   ABI_FREE(qmax)
   !  
   ! Variables/arrays related to the fine FFT grid ===
   ABI_MALLOC(ks_nhat,(nfftf,Dtset%nspden))
   ks_nhat=zero
   ABI_DT_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   call pawfgrtab_init(Pawfgrtab,cplex1,l_size_atm,Dtset%nspden,Dtset%typat)
   ABI_FREE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
   ! * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
   ! * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i2)')' bethe_salpeter : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
   !  
   ! Identify parts of the rectangular grid where the density has to be calculated ===
   optcut=0; optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)
 else
   ABI_DT_MALLOC(Paw_pwff,(0))
 end if !End of PAW Initialization

 ! Consistency check and additional stuff done only for GW with PAW.
 if (Dtset%usepaw==1) then
   if (Dtset%ecutwfn < Dtset%ecut) then
     write(msg,"(5a)")&
&     "WARNING - ",ch10,&
&     "  It is highly recommended to use ecutwfn = ecut for GW calculations with PAW since ",ch10,&
&     "  an excessive truncation of the planewave basis set can lead to unphysical results."
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

!* Initialize wave function handler, allocate wavefunctions.
 my_minb=1; my_maxb=BSp%nbnds; mband=BSp%nbnds
 ABI_MALLOC(nband,(Kmesh%nibz,Dtset%nsppol))
 nband=mband

!At present, no memory distribution, each node has the full set of states.
 ABI_MALLOC(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.TRUE.

 ABI_MALLOC(keep_ur,(mband,Kmesh%nibz,Dtset%nsppol))
 keep_ur=.FALSE.
 if (MODULO(Dtset%gwmem,10)==1) keep_ur = .TRUE.

!opt_ecut=zero
 opt_ecut=Dtset%ecutwfn

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,BSp%npwwfn,mband,nband,Kmesh%nibz,Dtset%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,ngfft_osc,&
& Gsph_x%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 ABI_FREE(bks_mask)
 ABI_FREE(nband)
 ABI_FREE(keep_ur)

 call wfd_print(Wfd,header="Wavefunctions used to construct the e-h basis set",mode_paral='PERS')

 call timab(651,2,tsec) ! bse(Init1)
 call timab(653,1,tsec) ! bse(rdkss)
 
 call wfd_read_wfk(Wfd,wfk_fname,iomode_from_fname(wfk_fname))

 ! This test has been disabled (too expensive!)
 if (.False.) call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

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
   keep_ur=.FALSE.
   if (MODULO(Dtset%gwmem,10)==1) keep_ur = .TRUE.

!  opt_ecut=zero
   opt_ecut=Dtset%ecutwfn

   call wfd_init(Wfd_dense,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,BSp%npwwfn,mband,nband,Kmesh_dense%nibz,Dtset%nsppol,&
&   bks_mask,Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk_dense%istwfk,Kmesh_dense%ibz,ngfft_osc,&
&   Gsph_x_dense%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

   ABI_FREE(bks_mask)
   ABI_FREE(nband)
   ABI_FREE(keep_ur)

   call wfd_print(Wfd_dense,header="Wavefunctions on the dense K-mesh used for interpolation",mode_paral='PERS')

   iomode = iomode_from_fname(dtfil%fnameabi_wfkfine)
   call wfd_read_wfk(Wfd_dense,Dtfil%fnameabi_wfkfine,iomode)

   ! This test has been disabled (too expensive!)
   if (.False.) call wfd_test_ortho(Wfd_dense,Cryst,Pawtab,unit=std_out,mode_paral="COLL")
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
 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)
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
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
&   Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,qphon,Cryst%rprimd,&
&   Cryst%ucvol,Dtset%usewvl,Cryst%xred)

   ! Evaluate onsite energies, potentials, densities ===
   !   * Initialize variables/arrays related to the PAW spheres.
   !   * Initialize also lmselect (index of non-zero LM-moments of densities).
   ABI_DT_MALLOC(KS_paw_ij,(Cryst%natom))
   call paw_ij_nullify(KS_paw_ij)

   has_dijso=Dtset%pawspnorb
   has_dijU=Dtset%usepawu
   has_dijso=Dtset%pawspnorb
   has_dijU=Dtset%usepawu

   call paw_ij_init(KS_paw_ij,cplex1,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=0,has_dijxc_hat=0,has_dijxc_val=0,&
&   has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   ABI_DT_MALLOC(KS_paw_an,(Cryst%natom))
   call paw_an_nullify(KS_paw_an)

   nkxc1=0
   call paw_an_init(KS_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,0,Dtset%nspden,&
&   cplex1,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=0)

   ! Calculate onsite vxc with and without core charge ===
   nzlmopt=-1; option=0; compch_sph=greatest_real

   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert0,&
&   Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
&   Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
&   Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
&   Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)
 end if !PAW

 if (.not.allocated(ks_nhatgr))  then
   ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
 !
 ! === For PAW, add the compensation charge on the FFT mesh, then get rho(G) ===
 if (Dtset%usepaw==1) ks_rhor(:,:)=ks_rhor(:,:)+ks_nhat(:,:)
 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_rhor,ucvol=ucvol)

 ABI_MALLOC(ks_rhog,(2,nfftf))

 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp0)

 call timab(655,2,tsec) ! bse(mkrho)
 !
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
&   Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
&   Dtset%nspden,Cryst%ntypat,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
&   Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&   k0,Dtset%spnorbscl,Cryst%ucvol,dtset%charge,ks_vtrial,ks_vxc,Cryst%xred,&
&   nucdipmom=Dtset%nucdipmom)
   !  
   ! Symmetrize KS Dij
   call symdij(Cryst%gprimd,Cryst%indsym,ipert0,&
&   Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
   !  
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
   call wrtout(std_out,msg,'COLL')
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
!  $ m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}> $
   ABI_MALLOC(m_lda_to_qp,(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   m_lda_to_qp=czero
   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz
       do band=1,Wfd%nband(ik_ibz,spin)
         m_lda_to_qp(band,band,ik_ibz,spin)=cone ! Initialize the QP amplitudes with KS wavefunctions.
       end do
     end do
   end do
!  
!  * Now read m_lda_to_qp and update the energies in QP_BSt.
!  TODO switch on the renormalization of n in sigma.
   ABI_MALLOC(prev_rhor,(nfftf,Wfd%nspden))
   ABI_DT_MALLOC(prev_Pawrhoij,(Cryst%natom*Wfd%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Wfd%usepaw,Wfd%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Wfd%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,prev_rhor,prev_Pawrhoij)

   ABI_FREE(prev_rhor)
   if (Psps%usepaw==1.and.nscf>0) then
     call pawrhoij_free(prev_pawrhoij)
   end if
   ABI_DT_FREE(prev_pawrhoij)
   !  
   !  if (nscf>0.and.wfd_iam_master(Wfd)) then ! Print the unitary transformation on std_out.
   !  call show_QP(QP_BSt,m_lda_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
   !  end if
   !  
   !  === Compute QP wfg as linear combination of KS states ===
   !  * Wfd%ug is modified inside calc_wf_qp
   !  * For PAW, update also the on-site projections.
   !  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
   !  TODO here we should use nbsc instead of nbnds

   call wfd_rotate(Wfd,Cryst,m_lda_to_qp)

   ABI_FREE(m_lda_to_qp)
   !  
   !  === Reinit the storage mode of Wfd as ug have been changed ===
   !  * Update also the wavefunctions for GW corrections on each processor
   call wfd_reset_ur_cprj(Wfd)

   call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

   ! Compute QP occupation numbers.
   call wrtout(std_out,'bethe_salpeter: calculating QP occupation numbers','COLL')

   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=0)
   ABI_MALLOC(qp_vbik,(QP_BSt%nkpt,QP_BSt%nsppol))
   qp_vbik(:,:) = get_valence_idx(QP_BSt)
   ABI_FREE(qp_vbik)

   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
 end if

 ABI_MALLOC(qp_rhog,(2,nfftf))
 call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Wfd%paral_kgb,0)

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
 call wfd_wave_free(Wfd,"All",bks_mask)
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
&   Dtset%iomode,ngfftf,nfftf_tot,Wfd%nsppol,Wfd%nspden,qp_aerhor,Wfd%prtvol,Wfd%comm)
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
 call wfd_wave_free(Wfd,"Real_space",bks_mask)
 ABI_FREE(bks_mask)
 call timab(658,2,tsec) ! bse(wfd_wave_free)

 ! Compute the commutator [r,Vu] (PAW only).
 ABI_DT_MALLOC(HUr,(Cryst%natom*Wfd%usepaw))

 call timab(659,1,tsec) ! bse(make_pawhur_t)
 if (Bsp%inclvkb/=0 .and. Wfd%usepaw==1 .and. Dtset%usepawu/=0) then !TODO here I need KS_Paw_ij
   MSG_WARNING("Commutator for LDA+U not tested")
   call pawhur_init(hur,Wfd%nsppol,Wfd%pawprtvol,Cryst,Pawtab,Pawang,Pawrad,KS_Paw_ij)
 end if
 call timab(659,2,tsec) ! bse(make_pawhur_t)

 select case (BSp%algorithm)
 case (BSE_ALGO_NONE)
   MSG_COMMENT("Skipping solution of the BSE equation")

 case (BSE_ALGO_DDIAGO, BSE_ALGO_CG)
   call timab(660,1,tsec) ! bse(exc_diago_driver)
   call exc_diago_driver(Wfd,Bsp,BS_files,KS_BSt,QP_BSt,Cryst,Kmesh,Psps,&
&   Pawtab,Hur,Hdr_bse,drude_plsmf,Epren)
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
&     Kmesh_dense,KS_BSt_dense,QP_BSt_dense,Wfd_dense,Vcp_dense,grid)
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
 !* Destroy local data structures.
 call destroy_mpi_enreg(MPI_enreg_seq)
 call crystal_free(Cryst)
 call gsph_free(Gsph_x)
 call gsph_free(Gsph_c)
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call hdr_free(Hdr_wfk)
 call hdr_free(Hdr_bse)
 call ebands_free(KS_BSt)
 call ebands_free(QP_BSt)
 call vcoul_free(Vcp)
 call pawhur_free(Hur)
 ABI_DT_FREE(Hur)
 call bs_parameters_free(BSp)
 call wfd_free(Wfd)
 call pawfgr_destroy(Pawfgr)
 call eprenorms_free(Epren)

 ! Free memory used for interpolation.
 if (BSp%use_interp) then 
   call double_grid_free(grid)
   call wfd_free(Wfd_dense)
   call gsph_free(Gsph_x_dense)
   call gsph_free(Gsph_c_dense)
   call kmesh_free(Kmesh_dense)
   call kmesh_free(Qmesh_dense)
   call ebands_free(KS_BSt_dense)
   call ebands_free(QP_BSt_dense)
   call vcoul_free(Vcp_dense)
   call hdr_free(Hdr_wfk_dense)
 end if

 ! Optional deallocation for PAW.
 if (Dtset%usepaw==1) then 
   call pawrhoij_free(KS_Pawrhoij)
   ABI_DT_FREE(KS_Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   ABI_DT_FREE(Pawfgrtab)
   call paw_ij_free(KS_paw_ij)
   ABI_DT_FREE(KS_paw_ij)
   call paw_an_free(KS_paw_an)
   ABI_DT_FREE(KS_paw_an)
   call pawpwff_free(Paw_pwff)
 end if
 ABI_DT_FREE(Paw_pwff)

 call timab(650,2,tsec) ! bse(Total)

 DBG_EXIT('COLL')

end subroutine bethe_salpeter
!!***

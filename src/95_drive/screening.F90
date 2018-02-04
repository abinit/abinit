!{\src2tex{textfont=tt}}
!!****f* ABINIT/screening
!! NAME
!! screening
!!
!! FUNCTION
!! Calculate screening and dielectric functions
!!
!! COPYRIGHT
!! Copyright (C) 2001-2018 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<datafiles_type)>=variables related to file names and unit numbers.
!! Pawang<pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! Psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  Before entering the first time in screening, a significant part of Psps has been initialized:
!!  the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid, ntypat,n1xccc,usepaw,useylm,
!!  and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!  the call to pspini. The next time the code enters screening, Psps might be identical to the
!!  one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!! Output is written on the main output file.
!! The symmetrical inverse dielectric matrix is stored in the _SCR file
!!
!! SIDE EFFECTS
!!  Dtset<type(dataset_type)>=all input variables for this dataset
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut) for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ... are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg) for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...Total density, potentials, ... are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used. It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf) are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      apply_scissor,calc_rpa_functional,cchi0,cchi0q0,chi0_bksmask
!!      chi0q0_intraband,chi_free,chkpawovlp,coeffs_gausslegint,crystal_free
!!      destroy_mpi_enreg,ebands_copy,ebands_free,ebands_update_occ
!!      em1params_free,energies_init,fourdp,get_gftt,getph,gsph_free,hdr_free
!!      hscr_free,hscr_io,init_distribfft_seq,initmpi_seq,kmesh_free,kxc_ada
!!      kxc_driver,littlegroup_free,lwl_write,make_epsm1_driver,metric,mkrdim
!!      nhatgrid,output_chi0sumrule,paw_an_free,paw_an_init,paw_an_nullify
!!      paw_gencond,paw_ij_free,paw_ij_init,paw_ij_nullify,paw_pwaves_lmn_free
!!      paw_pwaves_lmn_init,pawdenpot,pawdij,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawinit,pawmknhat,pawnabla_init,pawprt
!!      pawpuxinit,pawpwff_free,pawpwff_init,pawrhoij_alloc,pawrhoij_copy
!!      pawrhoij_free,pawtab_get_lsize,pawtab_print,print_arr,print_ngfft
!!      prtrhomxmn,pspini,random_stopping_power,rdgw,rdqps,rotate_fft_mesh
!!      setsymrhoij,setup_screening,setvtr,spectra_free,spectra_repr
!!      spectra_write,symdij,symdij_all,test_charge,timab,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_free,wfd_init,wfd_mkrho,wfd_print
!!      wfd_read_wfk,wfd_rotate,wfd_test_ortho,write_screening,wrtout
!!      xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine screening(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_xmpi
 use m_xomp
 use m_errors
 use m_ab7_mixing
 use m_kxc
 use m_nctk
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use libxc_functionals
 use m_hdr

 use m_io_tools,      only : open_file, file_exists, iomode_from_fname
 use m_fstrings,      only : int2char10, sjoin, strcat, itoa
 use m_energies,      only : energies_type, energies_init
 use m_numeric_tools, only : print_arr, iseven, coeffs_gausslegint
 use m_geometry,      only : normv, vdotw
 use m_gwdefs,        only : GW_TOLQ0, GW_TOLQ, em1params_free, em1params_t, GW_Q0_DEFAULT
 use m_mpinfo,        only : destroy_mpi_enreg
 use m_crystal,       only : crystal_free, crystal_t
#ifdef HAVE_NETCDF
 use m_crystal_io,    only : crystal_ncwrite
#endif
 use m_ebands,        only : ebands_update_occ, ebands_copy, get_valence_idx, get_occupied, apply_scissor, &
&                            ebands_free, ebands_has_metal_scheme, ebands_ncwrite
 use m_bz_mesh,       only : kmesh_t, kmesh_free, littlegroup_t, littlegroup_free
 use m_gsphere,       only : gsph_free, gsphere_t
 use m_vcoul,         only : vcoul_t, vcoul_free
 use m_qparticles,    only : rdqps, rdgw, show_QP
 use m_screening,     only : make_epsm1_driver, lwl_write, chi_t, chi_free, chi_new
 use m_io_screening,  only : hscr_new, hscr_io, write_screening, hscr_free, hscr_t
 use m_spectra,       only : spectra_t, spectra_write, spectra_repr, spectra_free, W_EM_LF, W_EM_NLF, W_EELF
 use m_fftcore,       only : print_ngfft
 use m_fft_mesh,      only : rotate_FFT_mesh, cigfft, get_gftt
 use m_wfd,           only : wfd_init, wfd_free,  wfd_nullify, wfd_print, wfd_t, wfd_rotate, wfd_test_ortho,&
&                            wfd_read_wfk, wfd_test_ortho, wfd_copy, wfd_change_ngfft
 use m_chi0,          only : output_chi0sumrule
 use m_pawang,        only : pawang_type
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,        only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,        only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,     only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy,&
&                            pawrhoij_free, symrhoij, pawrhoij_get_nspden
 use m_pawdij,        only : pawdij, symdij_all
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawpwij,       only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_pawfgr,        only : pawfgr_type, pawfgr_init, pawfgr_destroy

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'screening'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_51_manage_mpi
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_64_psp
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_69_wfdesc
 use interfaces_70_gw
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Dtset%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Dtset%usepaw)

!Local variables ------------------------------
 character(len=4) :: ctype='RPA '
!scalars
 integer,parameter :: level30=30,tim_fourdp4=4,NOMEGA_PRINTED=15,master=0
 integer :: spin,ik_ibz,my_nbks
 integer :: choice,cplex,dim_kxcg,dim_wing,ount,omp_ncpus
 integer :: fform_chi0,fform_em1,gnt_option,iat,ider,idir,ierr,band
 integer :: ifft,ii,ikbz,ikxc,initialized,iomega,ios,ipert
 integer :: iqibz,iqcalc,is_qeq0,isym,izero,ifirst,ilast
 integer :: label,mgfftf,mgfftgw
 integer :: nt_per_proc,work_size
 integer :: moved_atm_inside,moved_rhor
 integer :: nbcw,nbsc,nbvw,nkxc,nkxc1,n3xccc,optene,istep
 integer :: nfftf,nfftf_tot,nfftgw,nfftgw_tot,ngrvdw,nhatgrdim,nprocs,nspden_rhoij
 integer :: nscf,nzlmopt,mband
 integer :: optcut,optgr0,optgr1,optgr2,option,approx_type,option_test,optgrad
 integer :: optrad,optrhoij,psp_gencond,my_rank
 integer :: rhoxsp_method,comm,test_type,tordering,unt_em1,unt_susc,usexcnhat
 real(dp) :: compch_fft,compch_sph,domegareal,e0,ecore,ecut_eff,ecutdg_eff
 real(dp) :: gsqcutc_eff,gsqcutf_eff,gsqcut_shp,omegaplasma,ucvol,vxcavg,gw_gsq,r_s
 real(dp) :: alpha,rhoav,opt_ecut,factor
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem,ug_mem,ur_mem,cprj_mem
 logical :: found,iscompatibleFFT,use_tr,is_first_qcalc
 logical :: add_chi0_intraband,update_energies,call_pawinit
 character(len=10) :: string
 character(len=500) :: msg
 character(len=80) :: bar
 type(ebands_t) :: KS_BSt,QP_BSt
 type(kmesh_t) :: Kmesh,Qmesh
 type(vcoul_t) :: Vcp
 type(crystal_t) :: Cryst
 type(em1params_t) :: Ep
 type(Energies_type) :: KS_energies
 type(gsphere_t) :: Gsph_epsG0,Gsph_wfn
 type(Hdr_type) :: Hdr_wfk,Hdr_local
 type(MPI_type) :: MPI_enreg_seq
 type(Pawfgr_type) :: Pawfgr
 type(hscr_t) :: Hem1,Hchi0
 type(wfd_t) :: Wfd,Wfdf
 type(spectra_t) :: spectra
 type(chi_t) :: chihw
 type(wvl_data) :: wvl_dummy
!arrays
 integer :: ibocc(Dtset%nsppol),ngfft_gw(18),ngfftc(18),ngfftf(18)
 integer,allocatable :: irottb(:,:),ktabr(:,:),ktabrf(:,:),l_size_atm(:)
 integer,allocatable :: ks_vbik(:,:),ks_occ_idx(:,:),qp_vbik(:,:),nband(:,:)
 integer,allocatable :: nq_spl(:),nlmn_atm(:),gw_gfft(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),k0(3),qtmp(3),rmet(3,3),rprimd(3,3),tsec(2),strsxc(6)
 real(dp),allocatable :: igwene(:,:,:),chi0_sumrule(:),ec_rpa(:),rspower(:)
 real(dp),allocatable :: nhat(:,:),nhatgr(:,:,:),ph1d(:,:),ph1df(:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhor_p(:,:),rhor_kernel(:,:),taug(:,:),taur(:,:)
 real(dp),allocatable :: z(:),zw(:),grchempottn(:,:),grewtn(:,:),grvdw(:,:),kxc(:,:),qmax(:)
 real(dp),allocatable :: ks_vhartr(:),vpsp(:),ks_vtrial(:,:),ks_vxc(:,:),xccc3d(:)
 complex(gwpc),allocatable :: arr_99(:,:),kxcg(:,:),fxc_ADA(:,:,:)
 complex(dpc),allocatable :: m_lda_to_qp(:,:,:,:)
 complex(dpc),allocatable :: chi0_lwing(:,:,:),chi0_uwing(:,:,:),chi0_head(:,:,:)
 complex(dpc),allocatable :: chi0intra_lwing(:,:,:),chi0intra_uwing(:,:,:),chi0intra_head(:,:,:)
 complex(gwpc),allocatable,target :: chi0(:,:,:),chi0intra(:,:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: epsm1(:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:)
 character(len=80) :: title(2)
 character(len=fnlen) :: gw_fname,wfk_fname,lwl_fname
 type(littlegroup_t),pointer :: Ltg_q(:)
 type(Paw_an_type),allocatable :: Paw_an(:)
 type(Paw_ij_type),allocatable :: Paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:),prev_Pawrhoij(:)
 type(pawpwff_t),allocatable :: Paw_pwff(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 DBG_ENTER("COLL")

 call timab(301,1,tsec) ! overall time
 call timab(302,1,tsec) ! screening(init

 write(msg,'(6a)')&
& ' SCREENING: Calculation of the susceptibility and dielectric matrices ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
 call wrtout(ab_out, msg,'COLL')
 call wrtout(std_out,msg,'COLL')

 if(dtset%ucrpa>0) then
   write(msg,'(6a)')ch10,&
&   ' cRPA Calculation: The calculation of the polarisability is constrained (ucrpa/=0)',ch10
   call wrtout(ab_out, msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if
#if defined HAVE_GW_DPC
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
 call wrtout(ab_out, msg,'COLL')

 ! === Initialize MPI variables, and parallelization level ===
 ! gwpara: 0--> sequential run, 1--> parallelism over k-points, 2--> parallelism over bands.
 ! gwpara==2, each node has both fully and partially occupied states while conduction bands are divided
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 if (my_rank == master) then
   wfk_fname = dtfil%fnamewffk
   if (nctk_try_fort_or_ncfile(wfk_fname, msg) /= 0) then
     MSG_ERROR(msg)
   end if
 end if
 call xmpi_bcast(wfk_fname, master, comm, ierr)

 ! Some variables need to be initialized/nullify at start
 call energies_init(KS_energies)
 usexcnhat=0

 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)

!=== Define FFT grid(s) sizes ===
! Be careful! This mesh is only used for densities and potentials. It is NOT the (usually coarser)
! GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
! See also NOTES in the comments at the beginning of this file.
! NOTE: The mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

 k0(:)=zero
 call pawfgr_init(Pawfgr,Dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 call print_ngfft(ngfftf,'Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))

 ! We can intialize MPI_enreg and fft distrib here, now ngfft are known
 call initmpi_seq(MPI_enreg_seq) ! Fake MPI_type for the sequential part.
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')

!=============================================
!==== Open and read pseudopotential files ====
!=============================================
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level30,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)

!=== Initialize dimensions and basic objects ===
 call setup_screening(codvsn,acell,rprim,ngfftf,wfk_fname,dtfil,Dtset,Psps,Pawtab,&
& ngfft_gw,Hdr_wfk,Hdr_local,Cryst,Kmesh,Qmesh,KS_BSt,Ltg_q,Gsph_epsG0,Gsph_wfn,Vcp,Ep,comm)

 call timab(302,2,tsec) ! screening(init)
 call print_ngfft(ngfft_gw,'FFT mesh used for oscillator strengths')

 nfftgw_tot=PRODUCT(ngfft_gw(1:3))
 mgfftgw   =MAXVAL (ngfft_gw(1:3))
 nfftgw    =nfftgw_tot ! no FFT //

!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol

!==========================
!=== PAW initialization ===
!==========================
 if (Dtset%usepaw==1) then
   call timab(315,1,tsec) ! screening(pawin

   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)

   ABI_DT_MALLOC(Pawrhoij,(Cryst%natom))
   nspden_rhoij=pawrhoij_get_nspden(Dtset%nspden,Dtset%nspinor,Dtset%pawspnorb)
   call pawrhoij_alloc(Pawrhoij,Dtset%pawcpxocc,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,&
&   Cryst%typat,pawtab=Pawtab)

   ! Initialize values for several basic arrays stored in Pawinit
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

!  Initialize optional flags in Pawtab to zero
!  (Cannot be done in Pawinit since the routine is called only if some pars. are changed)
   Pawtab(:)%has_nabla = 0
   Pawtab(:)%usepawu   = 0
   Pawtab(:)%useexexch = 0
   Pawtab(:)%exchmix   =zero
!
!  * Evaluate <phi_i|nabla|phi_j>-<tphi_i|nabla|tphi_j> for the long wavelength limit.
!  TODO solve problem with memory leak and clean this part as well as the associated flag
   call pawnabla_init(Psps%mpsang,Cryst%ntypat,Pawrad,Pawtab)

   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,rprimd,Cryst%symrec,Pawang%zarot)

!  * Initialize and compute data for LDA+U.
!  paw_dmft%use_dmft=dtset%usedmft
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla,&
&     Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,Cryst%ntypat,Pawang,Dtset%pawprtvol,&
&     Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu,dtset%ucrpa)
   end if

   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)

   ! Re-symmetrize symrhoij.
!  this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1; ipert=0; idir=0
!  call symrhoij(Pawrhoij,Pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,Cryst%natom,&
!  &             Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,Pawtab,&
!  &             Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
!
   ! Evaluate form factors for the radial part of phi.phj-tphi.tphj ===
   ! rhoxsp_method=1 ! Arnaud-Alouani
   ! rhoxsp_method=2 ! Shiskin-Kresse
   rhoxsp_method=2

   ! At least for ucrpa, the Arnaud Alouani is always a better choice but needs a larger cutoff
   if(dtset%ucrpa>0) rhoxsp_method=1
   if (Dtset%pawoptosc /= 0) rhoxsp_method = Dtset%pawoptosc

   ABI_MALLOC(gw_gfft,(3,nfftgw_tot))
   call get_gftt(ngfft_gw,(/zero,zero,zero/),gmet,gw_gsq,gw_gfft)
   ABI_FREE(gw_gfft)

!  Set up q grids, make qmax 20% larger than largest expected:
   ABI_MALLOC(nq_spl,(Psps%ntypat))
   ABI_MALLOC(qmax,(Psps%ntypat))
   nq_spl = Psps%mqgrid_ff
   qmax = SQRT(gw_gsq)*1.2d0  !qmax = Psps%qgrid_ff(Psps%mqgrid_ff)
   ABI_DT_MALLOC(Paw_pwff,(Psps%ntypat))

   call pawpwff_init(Paw_pwff,rhoxsp_method,nq_spl,qmax,gmet,Pawrad,Pawtab,Psps)
   ABI_FREE(nq_spl)
   ABI_FREE(qmax)

   ! Variables/arrays related to the fine FFT grid
   ABI_MALLOC(nhat,(nfftf,Dtset%nspden))
   nhat=zero; cplex=1
   ABI_DT_MALLOC(Pawfgrtab,(Cryst%natom))
   call pawtab_get_lsize(Pawtab,l_size_atm,Cryst%natom,Cryst%typat)
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Dtset%nspden,Dtset%typat)
   ABI_FREE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
!  * 0 --> Vloc in atomic data is Vbare    (Blochl s formulation)
!  * 1 --> Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i3)')' screening : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
!
!  Identify parts of the rectangular grid where the density has to be calculated.
   optcut=0;optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   call timab(315,2,tsec) ! screening(pawin
 else
   ! allocate empty structure for the sake of -fcheck=all...
   ABI_DT_MALLOC(Paw_pwff,(0))
   ABI_DT_MALLOC(Pawrhoij,(0))
   ABI_DT_MALLOC(Pawfgrtab,(0))
 end if ! End of PAW initialization.

 ! Consistency check and additional stuff done only for GW with PAW.
 ABI_DT_MALLOC(Paw_onsite,(Cryst%natom))
 if (Dtset%usepaw==1) then
   if (Dtset%ecutwfn < Dtset%ecut) then
     write(msg,"(5a)")&
&     "WARNING - ",ch10,&
&     "  It is highly recommended to use ecutwfn = ecut for GW calculations with PAW since ",ch10,&
&     "  an excessive truncation of the planewave basis set can lead to unphysical results."
     call wrtout(ab_out,msg,'COLL')
   end if

   ABI_CHECK(Dtset%useexexch==0,"LEXX not yet implemented in GW")
   ABI_CHECK(Dtset%usedmft==0,"DMFT + GW not available")

   if (Dtset%pawcross==1) then
     optgrad=1
     call paw_pwaves_lmn_init(Paw_onsite,Cryst%natom,Cryst%natom,Cryst%ntypat,Cryst%rprimd,&
&     Cryst%xcart,Pawtab,Pawrad,Pawfgrtab,optgrad)
   end if
 end if

!Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(nhat))  then
   ABI_MALLOC(nhat,(nfftf,0))
 end if

 call timab(316,1,tsec) ! screening(wfs

!=====================================================
!=== Prepare the distribution of the wavefunctions ===
!=====================================================
! valence and partially occupied are replicate on each node  while conduction bands are MPI distributed.
! This method is mandatory if gwpara==2 and/or we are using awtr==1 or the spectral method.
! If awtr==1, we evaluate chi0 taking advantage of time-reversal (speed-up~2)
! Useful indeces:
!       nbvw = Max. number of fully/partially occupied states over spin
!       nbcw = Max. number of unoccupied states considering the spin
!TODO:
!  Here for semiconducting systems we have to be sure that each processor has all the
!  states considered in the SCGW, moreover nbsc<nbvw
!  in case of SCGW vale and conduction has to be recalculated to avoid errors
!  if a metal becomes semiconductor or viceversa.
!  Ideally nbvw should include only the states v such that the transition
!  c-->v is taken into account in cchi0 (see GW_TOLDOCC). In the present implementation

 ABI_MALLOC(ks_occ_idx,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_MALLOC(ks_vbik   ,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_MALLOC(qp_vbik   ,(KS_BSt%nkpt,KS_BSt%nsppol))

 call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=0)
 ks_occ_idx = get_occupied(KS_BSt,tol8) ! tol8 to be consistent when the density
 ks_vbik    = get_valence_idx(KS_BSt)

 ibocc(:)=MAXVAL(ks_occ_idx(:,:),DIM=1) ! Max occupied band index for each spin.
 ABI_FREE(ks_occ_idx)

 use_tr=.FALSE.; nbvw=0
 if (Dtset%gwpara==2.or.Ep%awtr==1.or.Dtset%spmeth>0) then
   use_tr=.TRUE.
   nbvw=MAXVAL(ibocc)
   nbcw=Ep%nbnds-nbvw
   write(msg,'(4a,i0,2a,i0,2a,i0,a)')ch10,&
&   '- screening: taking advantage of time-reversal symmetry ',ch10,&
&   '- Maximum band index for partially occupied states nbvw = ',nbvw,ch10,&
&   '- Remaining bands to be divided among processors   nbcw = ',nbcw,ch10,&
&   '- Number of bands treated by each node ~',nbcw/nprocs,ch10
   call wrtout(ab_out,msg,'COLL')
   if (Cryst%timrev/=2) then
     MSG_ERROR('Time-reversal cannot be used since cryst%timrev/=2')
   end if
 end if

 mband=Ep%nbnds
 ABI_MALLOC(nband,(Kmesh%nibz,Dtset%nsppol))
 nband=mband
 ABI_MALLOC(bks_mask,(mband,Kmesh%nibz,Dtset%nsppol))
 ABI_MALLOC(keep_ur,(mband,Kmesh%nibz,Dtset%nsppol))
 bks_mask=.FALSE.; keep_ur=.FALSE.

 ! autoparal section
 if (dtset%max_ncpus /=0) then
   ount = ab_out
   ! Temporary table needed to estimate memory
   ABI_MALLOC(nlmn_atm,(Cryst%natom))
   if (Dtset%usepaw==1) then
     do iat=1,Cryst%natom
       nlmn_atm(iat)=Pawtab(Cryst%typat(iat))%lmn_size
     end do
   end if

   write(ount,'(a)')"--- !Autoparal"
   write(ount,"(a)")"# Autoparal section for Screening runs"
   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ount,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ount,"(a,i0)")"    gwpara: ",dtset%gwpara
   write(ount,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    nbnds: ",Ep%nbnds

   work_size = nbvw * nbcw * Kmesh%nibz**2 * Dtset%nsppol

   ! Non-scalable memory in Mb i.e. memory that is not distribute with MPI.
   nonscal_mem = (two*gwpc*Ep%npwe**2*(Ep%nomega*b2Mb)) * 1.1_dp

   ! List of configurations.
   ! Assuming an OpenMP implementation with perfect speedup!
   write(ount,"(a)")"configurations:"

   do ii=1,dtset%max_ncpus
     nt_per_proc = 0
     eff = HUGE(one)
     max_wfsmem_mb = zero

     do my_rank=0,ii-1
       call chi0_bksmask(Dtset,Ep,Kmesh,nbvw,nbcw,my_rank,ii,bks_mask,keep_ur,ierr)
       if (ierr /= 0) exit
       nt_per_proc = MAX(nt_per_proc, COUNT(bks_mask(1:nbvw,:,:)) * COUNT(bks_mask(nbvw+1:,:,:)))
       eff = MIN(eff, (one * work_size) / (ii * nt_per_proc))

       ! Memory needed for Fourier components ug.
       my_nbks = COUNT(bks_mask)
       ug_mem = two*gwpc*Dtset%nspinor*Ep%npwwfn*my_nbks*b2Mb

       ! Memory needed for real space ur.
       ur_mem = two*gwpc*Dtset%nspinor*nfftgw*COUNT(keep_ur)*b2Mb

       ! Memory needed for PAW projections Cprj
       cprj_mem = zero
       if (Dtset%usepaw==1) cprj_mem = dp*Dtset%nspinor*SUM(nlmn_atm)*my_nbks*b2Mb

       max_wfsmem_mb = MAX(max_wfsmem_mb, ug_mem + ur_mem + cprj_mem)
     end do
     if (ierr /= 0) cycle

     ! Add the non-scalable part and increase by 10% to account for other datastructures.
     mempercpu_mb = (max_wfsmem_mb + nonscal_mem) * 1.1_dp

     do omp_ncpus=1,xomp_get_max_threads()
       write(ount,"(a,i0)")"    - tot_ncpus: ",ii * omp_ncpus
       write(ount,"(a,i0)")"      mpi_ncpus: ",ii
       write(ount,"(a,i0)")"      omp_ncpus: ",omp_ncpus
       write(ount,"(a,f12.9)")"      efficiency: ",eff
       write(ount,"(a,f12.2)")"      mem_per_cpu: ",mempercpu_mb
     end do
   end do
   write(ount,'(a)')"..."

   ABI_FREE(nlmn_atm)
   MSG_ERROR_NODUMP("aborting now")
 else
   call chi0_bksmask(Dtset,Ep,Kmesh,nbvw,nbcw,my_rank,nprocs,bks_mask,keep_ur,ierr)
 end if

! Initialize the Wf_info object (allocate %ug and %ur if required).
 opt_ecut=Dtset%ecutwfn
!opt_ecut=zero

!call gsph_init(Gsph_wfn,Cryst,Ep%npwvec,gvec=gvec_kss)

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Ep%npwwfn,mband,nband,Ep%nkibz,Dtset%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,ngfft_gw,&
& Gsph_wfn%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 if (Dtset%pawcross==1) then
   call wfd_init(Wfdf,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Ep%npwwfn,mband,nband,Ep%nkibz,Dtset%nsppol,bks_mask,&
&   Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,ngfft_gw,&
&   Gsph_wfn%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)
 end if

 ABI_FREE(bks_mask)
 ABI_FREE(nband)
 ABI_FREE(keep_ur)

 call wfd_print(Wfd,mode_paral='PERS')
!FIXME: Rewrite the treatment of use_tr branches in cchi0 ...
!Use a different nbvw for each spin.
!Now use_tr means that one can use time-reversal symmetry.

!==================================================
!==== Read KS band structure from the KSS file ====
!==================================================
 call wfd_read_wfk(Wfd,wfk_fname,iomode_from_fname(wfk_fname))

 if (Dtset%pawcross==1) then
   call wfd_copy(Wfd,Wfdf)
   call wfd_change_ngfft(Wfdf,Cryst,Psps,ngfftf)
 end if

 ! This test has been disabled (too expensive!)
 if (.False.) call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(316,2,tsec) ! screening(wfs
 call timab(319,1,tsec) ! screening(1)

 if (Cryst%nsym/=Dtset%nsym .and. Dtset%usepaw==1) then
   MSG_ERROR('Cryst%nsym/=Dtset%nsym, check pawinit and symrhoij')
 end if

! Get the FFT index of $ (R^{-1}(r-\tau)) $
! S= $\transpose R^{-1}$ and k_BZ = S k_IBZ
! irottb is the FFT index of $ R^{-1} (r-\tau) $ used to symmetrize u_Sk.
 ABI_MALLOC(irottb,(nfftgw,Cryst%nsym))
 call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfft_gw,irottb,iscompatibleFFT)

 ABI_MALLOC(ktabr,(nfftgw,Kmesh%nbz))
 do ikbz=1,Kmesh%nbz
   isym=Kmesh%tabo(ikbz)
   do ifft=1,nfftgw
     ktabr(ifft,ikbz)=irottb(ifft,isym)
   end do
 end do
 ABI_FREE(irottb)

 if (Dtset%usepaw==1 .and. Dtset%pawcross==1) then
   ABI_MALLOC(irottb,(nfftf,Cryst%nsym))
   call rotate_FFT_mesh(Cryst%nsym,Cryst%symrel,Cryst%tnons,ngfftf,irottb,iscompatibleFFT)

   ABI_MALLOC(ktabrf,(nfftf,Kmesh%nbz))
   do ikbz=1,Kmesh%nbz
     isym=Kmesh%tabo(ikbz)
     do ifft=1,nfftf
       ktabrf(ifft,ikbz)=irottb(ifft,isym)
     end do
   end do
   ABI_FREE(irottb)
 else
   ABI_MALLOC(ktabrf,(0,0))
 end if

!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
 write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
 !if (Dtset%mgfftdg/=mgfftf) write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
 ABI_MALLOC(ph1d,(2,3*(2*Dtset%mgfft+1)*Cryst%natom))
 ABI_MALLOC(ph1df,(2,3*(2*mgfftf+1)*Cryst%natom))
 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)

 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
   call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

! Initialize QP_BSt using KS bands
! In case of SCGW, update QP_BSt using the QPS file.
 call ebands_copy(KS_BSt,QP_BSt)

 call timab(319,2,tsec) ! screening(1)

!============================
!==== Self-consistent GW ====
!============================
 if (Ep%gwcalctyp>=10) then
   call timab(304,1,tsec) ! KS => QP; [wfrg]

   ! Initialize with KS eigenvalues and eigenfunctions.
   ABI_MALLOC(m_lda_to_qp,(Wfd%mband,Wfd%mband,Wfd%nkibz,Wfd%nsppol))
   m_lda_to_qp = czero
   do spin=1,Wfd%nsppol
     do ik_ibz=1,Wfd%nkibz
       do band=1,Wfd%nband(ik_ibz,spin)
         m_lda_to_qp(band,band,:,:) = cone
       end do
     end do
   end do

   ! Read unitary transformation and QP energies.
   ! TODO switch on the renormalization of n in screening, QPS should report bdgw
   ABI_MALLOC(rhor_p,(nfftf,Dtset%nspden))
   ABI_DT_MALLOC(prev_Pawrhoij,(Cryst%natom*Psps%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Dtset%usepaw,Dtset%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Dtset%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,m_lda_to_qp,rhor_p,prev_Pawrhoij)

   ABI_FREE(rhor_p)
   ABI_DT_FREE(prev_Pawrhoij)

   ! FIXME this is to preserve the old implementation for the head and the wings in ccchi0q0
   ! But has to be rationalized
   KS_BSt%eig=QP_BSt%eig

   ! Calculate new occ. factors and fermi level.
   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget)
   qp_vbik(:,:) = get_valence_idx(QP_BSt)

   ! === Update only the wfg treated with GW ===
   ! For PAW update and re-symmetrize cprj in the full BZ, TODO add rotation in spinor space
   if (nscf/=0) call wfd_rotate(Wfd,Cryst,m_lda_to_qp)

   ABI_FREE(m_lda_to_qp)
   call timab(304,2,tsec)
 end if ! gwcalctyp>=10

 call timab(305,1,tsec) ! screening(densit
!
!=== In case update the eigenvalues ===
!* Either use a scissor operator or an external GW file.
 gw_fname = "__in.gw__"
 update_energies = file_exists(gw_fname)

 if (ABS(Ep%mbpt_sciss)>tol6) then
   write(msg,'(5a,f7.3,a)')&
&   ' screening : performing a first self-consistency',ch10,&
&   ' update of the energies in W by a scissor operator',ch10,&
&   ' applying a scissor operator of [eV] : ',Ep%mbpt_sciss*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   call apply_scissor(QP_BSt,Ep%mbpt_sciss)
 else if (update_energies) then
   write(msg,'(4a)')&
&   ' screening : performing a first self-consistency',ch10,&
&   ' update of the energies in W by a previous GW calculation via GW file: ',TRIM(gw_fname)
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   ABI_MALLOC(igwene,(QP_Bst%mband,QP_Bst%nkpt,QP_Bst%nsppol))
   call rdgw(QP_BSt,gw_fname,igwene,extrapolate=.FALSE.)
!  call rdgw(QP_BSt,gw_fname,igwene,extrapolate=.TRUE.)
   ABI_FREE(igwene)
   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget)
 end if

!========================
!=== COMPUTE DENSITY ====
!========================
!* Evaluate PW part (complete charge in case of NC pseudos)
!TODO this part has to be rewritten. If I decrease the tol on the occupations
!I have to code some MPI stuff also if use_tr==.TRUE.

 ABI_MALLOC(rhor,(nfftf,Dtset%nspden))
 ABI_MALLOC(taur,(nfftf,Dtset%nspden*Dtset%usekden))

 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,rhor)
 if (Dtset%usekden==1) then
   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,taur,optcalc=1)
 end if

 call timab(305,2,tsec) ! screening(densit

 nhatgrdim = 0
 if (Dtset%usepaw==1) then ! Additional computation for PAW.
   call timab(320,1,tsec) ! screening(paw

!  Add the compensation charge to the PW density.
   nhatgrdim=0; if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   cplex=1; ider=2*nhatgrdim; izero=0
   if (nhatgrdim>0)  then
     ABI_MALLOC(nhatgr,(nfftf,Dtset%nspden,3))
   else
     ABI_MALLOC(nhatgr,(nfftf,Dtset%nspden,0))
   end if
   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
&   Pawfgrtab,nhatgr,nhat,Pawrhoij,Pawrhoij,Pawtab,k0,Cryst%rprimd,Cryst%ucvol,dtset%usewvl,Cryst%xred)

!  === Evaluate onsite energies, potentials, densities ===
!  * Initialize variables/arrays related to the PAW spheres.
!  * Initialize also lmselect (index of non-zero LM-moments of densities).
   cplex=1
   ABI_DT_MALLOC(Paw_ij,(Cryst%natom))
   call paw_ij_nullify(Paw_ij)
   call paw_ij_init(Paw_ij,cplex,Dtset%nspinor,Wfd%nsppol,&
&   Wfd%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_exexch_pot=1,has_pawu_occ=1)

   nkxc1=0
   ABI_DT_MALLOC(Paw_an,(Cryst%natom))
   call paw_an_nullify(Paw_an)
   call paw_an_init(Paw_an,Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,&
&   cplex,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=0)

   nzlmopt=-1; option=0; compch_sph=greatest_real
   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert,Dtset%ixc,&
&   Cryst%natom,Cryst%natom,Dtset%nspden,Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,Paw_an,Paw_an,&
&   Paw_ij,Pawang,Dtset%pawprtvol,Pawrad,Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,&
&   Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)
   call timab(320,2,tsec) ! screening(paw
 else
   ABI_DT_MALLOC(Paw_ij,(0))
   ABI_DT_MALLOC(Paw_an,(0))
 end if ! usepaw

 call timab(321,1,tsec) ! screening(2)

 !JB : Should be remove : cf. l 839
 if (.not.allocated(nhatgr))  then
   ABI_MALLOC(nhatgr,(nfftf,Dtset%nspden,0))
 end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,rhor,ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,omegaplasma)

!For PAW, add the compensation charge the FFT mesh, then get rho(G).
 if (Dtset%usepaw==1) rhor(:,:)=rhor(:,:)+nhat(:,:)

 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,rhor,ucvol=ucvol)
 if(Dtset%usekden==1)then
   call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,taur,ucvol=ucvol,optrhor=1)
 end if

 if (dtset%gwgamma>0) then
   ABI_MALLOC(rhor_kernel,(nfftf,Dtset%nspden))
 end if

 ABI_MALLOC(rhog,(2,nfftf))
 ABI_MALLOC(taug,(2,nfftf*Dtset%usekden))
 call fourdp(1,rhog,rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp4)
 if(Dtset%usekden==1)then
   call fourdp(1,taug,taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp4)
 end if

!The following steps have been gathered in the setvtr routine:
!- get Ewald energy and Ewald forces
!- compute local ionic pseudopotential vpsp
!- eventually compute 3D core electron density xccc3d
!- eventually compute vxc and vhartr
!- set up ks_vtrial
!**************************************************************
!**** NOTE THAT Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
!**************************************************************

 ngrvdw=0
 ABI_MALLOC(grvdw,(3,ngrvdw))
 ABI_MALLOC(grchempottn,(3,Cryst%natom))
 ABI_MALLOC(grewtn,(3,Cryst%natom))
 nkxc=0
 if (Dtset%nspden==1) nkxc=2
 if (Dtset%nspden>=2) nkxc=3 ! check GGA and spinor that is messy !!!
 ! If MGGA, fxc and kxc are not available and we dont need them for the screening part (for now ...)
 if (Dtset%ixc<0 .and. libxc_functionals_ismgga()) nkxc=0
 if (nkxc/=0)  then
   ABI_MALLOC(kxc,(nfftf,nkxc))
 end if

 n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
 ABI_MALLOC(xccc3d,(n3xccc))
 ABI_MALLOC(ks_vhartr,(nfftf))
 ABI_MALLOC(ks_vtrial,(nfftf,Dtset%nspden))
 ABI_MALLOC(vpsp,(nfftf))
 ABI_MALLOC(ks_vxc,(nfftf,Dtset%nspden))

 optene=4; moved_atm_inside=0; moved_rhor=0; initialized=1; istep=1
 call setvtr(Cryst%atindx1,Dtset,KS_energies,Cryst%gmet,Cryst%gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,istep,kxc,mgfftf,&
& moved_atm_inside,moved_rhor,MPI_enreg_seq, &
& Cryst%nattyp,nfftf,ngfftf,ngrvdw,nhat,nhatgr,nhatgrdim,nkxc,Cryst%ntypat,&
& Psps%n1xccc,n3xccc,optene,pawrad,Pawtab,ph1df,Psps,rhog,rhor,Cryst%rmet,Cryst%rprimd,strsxc,Cryst%ucvol,usexcnhat,&
& ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,wvl_dummy,xccc3d,Cryst%xred,taug=taug,taur=taur)

 if (nkxc/=0)  then
   ABI_FREE(kxc)
 end if
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)
 ABI_FREE(grvdw)
 ABI_FREE(xccc3d)

!============================
!==== Compute KS PAW Dij ====
!============================
 if (Dtset%usepaw==1) then
   call timab(561,1,tsec)

!  Calculate unsymmetrized Dij.
   cplex=1; ipert=0; idir=0
   call pawdij(cplex,Dtset%enunit,Cryst%gprimd,ipert,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
&   Dtset%nspden,Cryst%ntypat,Paw_an,Paw_ij,Pawang,Pawfgrtab,Dtset%pawprtvol,&
&   Pawrad,Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,k0,Dtset%spnorbscl,&
&   Cryst%ucvol,dtset%charge,ks_vtrial,ks_vxc,Cryst%xred,&
&   nucdipmom=Dtset%nucdipmom)

!  Symmetrize KS Dij
#if 0
   call symdij(Cryst%gprimd,Cryst%indsym,ipert,Cryst%natom,Cryst%natom,&
&   Cryst%nsym,Cryst%ntypat,0,Paw_ij,Pawang,Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,&
&   Cryst%symrec)
#else
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert,Cryst%natom,Cryst%natom,&
&   Cryst%nsym,Cryst%ntypat,Paw_ij,Pawang,Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,&
&   Cryst%symrec)
#endif
!
!  Output of the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Cryst%natom,Paw_ij,Pawrhoij,Pawtab)
   call timab(561,2,tsec)
 end if
!
!=== Calculate the frequency mesh ===
!* First omega is always zero without broadening.
!FIXME what about metals? I think we should add eta, this means we need to know if the system is metallic, for example using occopt
!MS Modified to account for non-zero starting frequency (19-11-2010)
!MS Modified for tangent grid (07-01-2011)
 ABI_MALLOC(Ep%omega,(Ep%nomega))
 Ep%omega(1)=CMPLX(Ep%omegaermin,zero,kind=dpc)

 if (Ep%nomegaer>1) then ! Avoid division by zero.
   if (Dtset%gw_frqre_tangrid==0.and.Dtset%gw_frqre_inzgrid==0) then
     domegareal=(Ep%omegaermax-Ep%omegaermin)/(Ep%nomegaer-1)
     do iomega=2,Ep%nomegaer
       Ep%omega(iomega)=CMPLX(Ep%omegaermin+(iomega-1)*domegareal,zero,kind=dpc)
     end do
   else if (Dtset%gw_frqre_tangrid==1.and.Dtset%gw_frqre_inzgrid==0) then ! We have tangent transformed grid
     MSG_WARNING('EXPERIMENTAL - Using tangent transform grid for contour deformation.')
     Ep%omegaermax = Dtset%cd_max_freq
     Ep%omegaermin = zero
     ifirst=1; ilast=Ep%nomegaer
     if (Dtset%cd_subset_freq(1)/=0) then ! Only a subset of frequencies is being calculated
       ifirst=Dtset%cd_subset_freq(1); ilast=Dtset%cd_subset_freq(2)
     end if
     factor = Dtset%cd_halfway_freq/TAN(pi*quarter)
!    *Important*-here nfreqre is used because the step is set by the original grid
     domegareal=(ATAN(Ep%omegaermax/factor)*two*piinv)/(Dtset%nfreqre-1) ! Stepsize in transformed variable
     do iomega=1,Ep%nomegaer
       Ep%omega(iomega)=CMPLX(factor*TAN((iomega+ifirst-2)*domegareal*pi*half),zero,kind=dpc)
     end do
     Ep%omegaermin = REAL(Ep%omega(1))
     Ep%omegaermax = REAL(Ep%omega(Ep%nomegaer))
   else if (Dtset%gw_frqre_tangrid==0.and.Dtset%gw_frqre_inzgrid==1) then
     e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
     domegareal=one/(Ep%nomegaer)
     do iomega=1,Ep%nomegaer
       factor = (iomega-1)*domegareal
       Ep%omega(iomega)=CMPLX(e0*factor/(one-factor),zero,kind=dpc)
     end do
     Ep%omegaermin = REAL(Ep%omega(1))
     Ep%omegaermax = REAL(Ep%omega(Ep%nomegaer))
   else
     MSG_ERROR('Error in specification of real frequency grid')
   end if
 end if

 if (Ep%plasmon_pole_model.and.Ep%nomega==2) then
   e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
   Ep%omega(2)=CMPLX(zero,e0,kind=dpc)
 end if
!
!=== For AC, use Gauss-Legendre quadrature method ===
!* Replace $ \int_0^\infty dx f(x) $ with $ \int_0^1 dz f(1/z - 1)/z^2 $.
!* Note that the grid is not log as required by CD, so we cannot use the same SCR file.
 if (Ep%analytic_continuation) then
   ABI_MALLOC(z,(Ep%nomegaei))
   ABI_MALLOC(zw,(Ep%nomegaei))
   call coeffs_gausslegint(zero,one,z,zw,Ep%nomegaei)
   do iomega=1,Ep%nomegaei
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,one/z(iomega)-one,kind=dpc)
   end do
   ABI_FREE(z)
   ABI_FREE(zw)
 else if (Ep%contour_deformation.and.(Dtset%cd_customnimfrqs/=0)) then
   Ep%omega(Ep%nomegaer+1)=CMPLX(zero,Dtset%cd_imfrqs(1))
   do iomega=2,Ep%nomegaei
     if (Dtset%cd_imfrqs(iomega)<=Dtset%cd_imfrqs(iomega-1)) then
       MSG_ERROR(' Specified imaginary frequencies need to be strictly increasing!')
     end if
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,Dtset%cd_imfrqs(iomega))
   end do
 else if (Ep%contour_deformation.and.(Dtset%gw_frqim_inzgrid/=0)) then
   e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
   domegareal=one/(Ep%nomegaei+1)
   do iomega=1,Ep%nomegaei
     factor = iomega*domegareal
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,e0*factor/(one-factor),kind=dpc)
   end do
 else if (Ep%contour_deformation.and.(Ep%nomegaei/=0)) then
   e0=Dtset%ppmfrq; if (e0<0.1d-4) e0=omegaplasma
   do iomega=1,Ep%nomegaei
     Ep%omega(Ep%nomegaer+iomega)=CMPLX(zero,e0/(Dtset%freqim_alpha-two)&
&     *(EXP(two/(Ep%nomegaei+1)*LOG(Dtset%freqim_alpha-one)*iomega)-one),kind=dpc)
   end do
 end if

 if (Dtset%cd_full_grid/=0) then ! Full grid will be calculated
!  Grid values are added after the last imaginary freq.
   do ios=1,Ep%nomegaei
     do iomega=2,Ep%nomegaer
       Ep%omega(Ep%nomegaer+Ep%nomegaei+(ios-1)*(Ep%nomegaer-1)+(iomega-1)) = &
&       CMPLX(REAL(Ep%omega(iomega)),AIMAG(Ep%omega(Ep%nomegaer+ios)))
     end do
   end do
 end if

 ! Report frequency mesh for chi0.
 write(msg,'(2a)')ch10,' calculating chi0 at frequencies [eV] :'
 call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 do iomega=1,Ep%nomega
   write(msg,'(i3,2es16.6)')iomega,Ep%omega(iomega)*Ha_eV
   call wrtout(std_out,msg,'COLL'); call wrtout(ab_out,msg,'COLL')
 end do

 ! Allocate chi0, wings and array for chi0_sumrule check.
 ABI_MALLOC(chi0_sumrule,(Ep%npwe))

 write(msg,'(a,f12.1,a)')' Memory required for chi0 matrix= ',two*gwpc*Ep%npwe**2*Ep%nI*Ep%nJ*Ep%nomega*b2Mb," [Mb]."
 call wrtout(std_out,msg,'COLL')
 ABI_STAT_MALLOC(chi0,(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega), ierr)
 ABI_CHECK(ierr==0, "Out of memory in chi0")
!
!============================== END OF THE INITIALIZATION PART ===========================
!
!======================================================================
!==== Loop over q-points. Calculate \epsilon^{-1} and save on disc ====
!======================================================================
 call timab(321,2,tsec) ! screening(2)

 iqcalc = 0
 do iqibz=1,Qmesh%nibz
   call timab(306,1,tsec)
   is_first_qcalc=(iqibz==1)

   ! Selective q-point calculation.
   found=.FALSE.; label=iqibz
   if (Ep%nqcalc/=Ep%nqibz) then
     do ii=1,Ep%nqcalc
       qtmp(:)=Qmesh%ibz(:,iqibz)-Ep%qcalc(:,ii)
       found=(normv(qtmp,gmet,'G')<GW_TOLQ)
       if (found) then
         label=ii; EXIT !ii
       end if
     end do
     if (.not.found) CYCLE !iqibz
     qtmp(:)=Ep%qcalc(:,1)-Qmesh%ibz(:,iqibz)
     is_first_qcalc=(normv(qtmp,gmet,'G')<GW_TOLQ)
   end if
   iqcalc = iqcalc + 1

   bar=REPEAT('-',80)
   write(msg,'(4a,1x,a,i2,a,f9.6,2(",",f9.6),3a)')ch10,ch10,bar,ch10,&
&   ' q-point number ',label,'        q = (',(Qmesh%ibz(ii,iqibz),ii=1,3),') [r.l.u.]',ch10,bar
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   is_qeq0 = 0; if (normv(Qmesh%ibz(:,iqibz),gmet,'G')<GW_TOLQ0) is_qeq0=1

   call timab(306,2,tsec)

   if (is_qeq0==1) then
     ! Special treatment of the long wavelength limit.
     call timab(307,1,tsec)

     ABI_MALLOC(chi0_lwing,(Ep%npwe*Ep%nI,Ep%nomega,3))
     ABI_MALLOC(chi0_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,3))
     ABI_MALLOC(chi0_head,(3,3,Ep%nomega))

     call cchi0q0(use_tr,Dtset,Cryst,Ep,Psps,Kmesh,QP_BSt,KS_BSt,Gsph_epsG0,&
&     Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,Pawfgrtab,Paw_onsite,ktabr,ktabrf,nbvw,ngfft_gw,nfftgw,&
&     ngfftf,nfftf_tot,chi0,chi0_head,chi0_lwing,chi0_uwing,Ltg_q(iqibz),chi0_sumrule,Wfd,Wfdf)

     chihw = chi_new(ep%npwe, ep%nomega)
     chihw%head = chi0_head
     chihw%lwing = chi0_lwing
     chihw%uwing = chi0_uwing

     ! Add the intraband term if required and metallic occupation scheme is used.
     add_chi0_intraband=.FALSE. !add_chi0_intraband=.TRUE.
     if (add_chi0_intraband .and. ebands_has_metal_scheme(QP_BSt)) then

       ABI_STAT_MALLOC(chi0intra,(Ep%npwe*Ep%nI,Ep%npwe*Ep%nJ,Ep%nomega), ierr)
       ABI_CHECK(ierr==0, "Out of memory in chi0intra")

       ABI_MALLOC(chi0intra_lwing,(Ep%npwe*Ep%nI,Ep%nomega,3))
       ABI_MALLOC(chi0intra_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,3))
       ABI_MALLOC(chi0intra_head,(3,3,Ep%nomega))

       call chi0q0_intraband(Wfd,Cryst,Ep,Psps,QP_BSt,Gsph_epsG0,Pawang,Pawrad,Pawtab,Paw_ij,Paw_pwff,use_tr,Dtset%usepawu,&
       ngfft_gw,chi0intra,chi0intra_head,chi0intra_lwing,chi0intra_uwing)

       call wrtout(std_out,"Head of chi0 and chi0_intra","COLL")
       do iomega=1,Ep%nomega
         write(std_out,*)Ep%omega(iomega)*Ha_eV,REAL(chi0(1,1,iomega)),REAL(chi0intra(1,1,iomega))
         write(std_out,*)Ep%omega(iomega)*Ha_eV,AIMAG(chi0(1,1,iomega)),AIMAG(chi0intra(1,1,iomega))
       end do

       chi0       = chi0       + chi0intra
       chi0_head  = chi0_head  + chi0intra_head
       chi0_lwing = chi0_lwing + chi0intra_lwing
       chi0_uwing = chi0_uwing + chi0intra_uwing

       ABI_FREE(chi0intra)
       ABI_FREE(chi0intra_lwing)
       ABI_FREE(chi0intra_uwing)
       ABI_FREE(chi0intra_head)
     end if

     if (.False.) then
       lwl_fname = strcat(dtfil%filnam_ds(4), "_LWL")
       call lwl_write(lwl_fname,cryst,vcp,ep%npwe,ep%nomega,gsph_epsg0%gvec,chi0,chi0_head,chi0_lwing,chi0_uwing,comm)
     end if

     call timab(307,2,tsec)

   else
     ! Calculate cchi0 for q/=0.
     call timab(308,1,tsec)

     call cchi0(use_tr,Dtset,Cryst,Qmesh%ibz(:,iqibz),Ep,Psps,Kmesh,QP_BSt,Gsph_epsG0,&
&     Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,nbvw,ngfft_gw,nfftgw,ngfftf,nfftf_tot,chi0,ktabr,ktabrf,&
&     Ltg_q(iqibz),chi0_sumrule,Wfd,Wfdf)

     call timab(308,2,tsec)
   end if

   ! Print chi0(q,G,Gp,omega), then calculate epsilon and epsilon^-1 for this q-point.
   ! Only master works but this part could be parallelized over frequencies.
   call timab(309,1,tsec)

   do iomega=1,MIN(Ep%nomega,NOMEGA_PRINTED)
     write(msg,'(1x,a,i4,a,2f9.4,a)')' chi0(G,G'') at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     write(msg,'(1x,a,i3,a,i4,a)')' chi0(q =',iqibz, ', omega =',iomega,', G,G'')'
     if (Ep%nqcalc/=Ep%nqibz) write(msg,'(a,i3,a,i4,a)')'  chi0(q=',iqcalc,', omega=',iomega,', G,G'')'
     call wrtout(std_out,msg,'COLL')
!    arr99 is needed to avoid the update of all the tests. Now chi0 is divided by ucvol inside (cchi0|cchi0q0).
!    TODO should be removed but GW tests have to be updated.
     ii = MIN(9,Ep%npwe)
     ABI_MALLOC(arr_99,(ii,ii))
     arr_99 = chi0(1:ii,1:ii,iomega)*ucvol
     call print_arr(arr_99,max_r=2,unit=ab_out)
     call print_arr(arr_99,unit=std_out)
     ABI_FREE(arr_99)
     !call print_arr(chi0(:,:,iomega),max_r=2,unit=ab_out)
     !call print_arr(chi0(:,:,iomega),unit=std_out)
   end do

   if (Ep%nomega>NOMEGA_PRINTED) then
     write(msg,'(a,i3,a)')' No. of calculated frequencies > ',NOMEGA_PRINTED,', stop printing '
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end if

   ! Write chi0 to _SUSC file
   ! Master creates and write the header if this is the first q-point calculated.
   if (Dtset%prtsuscep>0 .and. my_rank==master) then
     title(1)='CHI0 file: chi0'
     title(2)=' '
     if (is_qeq0==1) then
       string='0'; if (Dtset%usepaw==0.and.Ep%inclvkb/=0) call int2char10(Ep%inclvkb,string)
       title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
     end if

     ! Open file and write header for polarizability files.
     if (is_first_qcalc) then
       ikxc=0; test_type=0; tordering=1
       hchi0 = hscr_new("polarizability",dtset,ep,hdr_local,ikxc,test_type,tordering,title,Ep%npwe,Gsph_epsG0%gvec)
       if (dtset%iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
         NCF_CHECK(nctk_open_create(unt_susc, nctk_ncify(dtfil%fnameabo_sus), xmpi_comm_self))
         NCF_CHECK(crystal_ncwrite(cryst, unt_susc))
         NCF_CHECK(ebands_ncwrite(QP_BSt, unt_susc))
#endif
       else
         unt_susc=Dtfil%unchi0
         if (open_file(dtfil%fnameabo_sus,msg,unit=unt_susc,status='unknown',form='unformatted') /= 0) then
           MSG_ERROR(msg)
         end if
       end if
       fform_chi0 = hchi0%fform
       call hscr_io(hchi0,fform_chi0,2,unt_susc,xmpi_comm_self,0,Dtset%iomode)
       call hscr_free(Hchi0)
     end if
     call write_screening("polarizability",unt_susc,Dtset%iomode,Ep%npwe,Ep%nomega,iqcalc,chi0)
   end if

!  Calculate the RPA functional correlation energy if the polarizability on a
!  Gauss-Legendre mesh along imaginary axis is available
   if (Ep%analytic_continuation .and. Dtset%gwrpacorr>0 ) then
     if (is_first_qcalc) then
       ABI_MALLOC(ec_rpa,(Dtset%gwrpacorr))
       ec_rpa(:)=zero
     end if
     call calc_rpa_functional(Dtset%gwrpacorr,label,iqibz,Ep,Vcp,Qmesh,Dtfil,gmet,chi0,comm,ec_rpa)
     if (label==Ep%nqcalc) then
       ABI_FREE(ec_rpa)
     end if
   end if

!  ==========================================================
!  === Calculate RPA \tilde\epsilon^{-1} overwriting chi0 ===
!  ==========================================================
   approx_type=0 ! RPA
   option_test=0 ! TESTPARTICLE
   dim_wing=0; if (is_qeq0==1) dim_wing=3

   if (dim_wing==0) then
     dim_wing=1
     if (.not.allocated(chi0_lwing))  then
       ABI_MALLOC(chi0_lwing,(Ep%npwe*Ep%nI,Ep%nomega,dim_wing))
     end if
     if (.not.allocated(chi0_uwing))  then
       ABI_MALLOC(chi0_uwing,(Ep%npwe*Ep%nJ,Ep%nomega,dim_wing))
     end if
     if (.not.allocated(chi0_head ))  then
       ABI_MALLOC(chi0_head,(dim_wing,dim_wing,Ep%nomega))
     end if
     dim_wing=0
   end if

#if 0
!  Using the random q for the optical limit is one of the reasons
!  why sigma breaks the initial energy degeneracies.
   chi0_lwing=czero
   chi0_uwing=czero
   chi0_head=czero
#endif

   ! Setup flags for the computation of em1
   ! If the vertex is being included for the spectrum, calculate the kernel now and pass it on
   if (dtset%gwgamma>0) rhor_kernel = rhor

   select case (dtset%gwgamma)
   case (0)
     approx_type=0; option_test=0; dim_kxcg=0
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))

   case (1, 2)
     ! ALDA TDDFT kernel vertex
     ABI_CHECK(Dtset%usepaw==0,"GWGamma=1 or 2 + PAW not available")
     MSG_WARNING('EXPERIMENTAL: Kernel is being added to screening, the SCR file will be non-standard!!')
     ikxc=7; approx_type=1; dim_kxcg=1
     if (Dtset%gwgamma==1) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==2) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))
     call kxc_driver(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,rhor_kernel,&
     Ep%npwe,dim_kxcg,kxcg,Gsph_epsG0%gvec,xmpi_comm_self)

   case (3, 4)
     ! ADA non-local kernel vertex
     ABI_CHECK(Wfd%usepaw==0,"ADA vertex + PAW not available")
     ABI_CHECK(Wfd%nsppol==1,"ADA vertex for GWGamma not available yet for spin-polarised cases")
     MSG_WARNING('EXPERIMENTAL: Kernel is being added to screening, the SCR file will be non-standard!!')
     ikxc=7; approx_type=2
     if (Dtset%gwgamma==3) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==4) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(fxc_ADA,(Ep%npwe,Ep%npwe,Ep%nqibz))
!    Use userrd to set kappa
     if (Dtset%userrd==zero) Dtset%userrd = 2.1_dp
!    Set correct value of kappa (should be scaled with alpha*r_s where)
!    r_s is Wigner-Seitz radius and alpha=(4/(9*Pi))^(1/3)
     rhoav = (omegaplasma*omegaplasma)/four_pi
     r_s = (three/(four_pi*rhoav))**third
     alpha = (four*ninth*piinv)**third
     Dtset%userrd = Dtset%userrd*alpha*r_s

     call kxc_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf,Wfd%nspden,rhor_kernel,Ep%npwe,Ep%nqibz,Ep%qibz,&
     fxc_ADA,Gsph_epsG0%gvec,xmpi_comm_self,kappa_init=Dtset%userrd)

     dim_kxcg=0
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))

!  @WC: bootstrap --
   case (-3, -4, -5, -6, -7, -8)
     ABI_CHECK(Dtset%usepaw==0,"GWGamma + PAW not available")
     if (Dtset%gwgamma>-5) then
       MSG_WARNING('EXPERIMENTAL: Bootstrap kernel is being added to screening')
       approx_type=4
     else if (Dtset%gwgamma>-7) then
       MSG_WARNING('EXPERIMENTAL: Bootstrap kernel (head-only) is being added to screening')
       approx_type=5
     else
       MSG_WARNING('EXPERIMENTAL: Bootstrap kernel (RPA-type, head-only) is being added to screening')
       approx_type=6
     end if
     dim_kxcg=0
     option_test=MOD(Dtset%gwgamma,2)
     ! 1 -> TESTELECTRON, vertex in chi0 *and* sigma
     ! 0 -> TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))
!--@WC

   case default
     MSG_ERROR(sjoin("Wrong gwgamma:", itoa(dtset%gwgamma)))
   end select

   if (approx_type<2) then
     call make_epsm1_driver(iqibz,dim_wing,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%omega,&
     approx_type,option_test,Vcp,nfftf_tot,ngfftf,dim_kxcg,kxcg,Gsph_epsG0%gvec,&
     chi0_head,chi0_lwing,chi0_uwing,chi0,spectra,comm)
   else if (approx_type<3) then !@WC: ADA
     call make_epsm1_driver(iqibz,dim_wing,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%omega,&
     approx_type,option_test,Vcp,nfftf_tot,ngfftf,dim_kxcg,kxcg,Gsph_epsG0%gvec,&
     chi0_head,chi0_lwing,chi0_uwing,chi0,spectra,comm,fxc_ADA=fxc_ADA(:,:,iqibz))
   else if (approx_type<7) then !@WC: bootstrap
     call make_epsm1_driver(iqibz,dim_wing,Ep%npwe,Ep%nI,Ep%nJ,Ep%nomega,Ep%omega,&
     approx_type,option_test,Vcp,nfftf_tot,ngfftf,dim_kxcg,kxcg,Gsph_epsG0%gvec,&
     chi0_head,chi0_lwing,chi0_uwing,chi0,spectra,comm)
   else
     MSG_ERROR(sjoin("Wrong approx_type:", itoa(approx_type)))
   end if

   ABI_FREE(chi0_lwing)
   ABI_FREE(chi0_uwing)
   ABI_FREE(chi0_head)

   if (my_rank==master .and. is_qeq0==1) then
     call spectra_repr(spectra,msg)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')

     if (Ep%nomegaer>2) then
       call spectra_write(spectra,W_EELF  ,Dtfil%fnameabo_eelf)
       call spectra_write(spectra,W_EM_LF ,Dtfil%fnameabo_em1_lf)
       call spectra_write(spectra,W_EM_NLF,Dtfil%fnameabo_em1_nlf)
     end if
   end if ! master and is_qeq0==1

   if (is_qeq0==1) call chi_free(chihw)

   call spectra_free(spectra)
   if (allocated(kxcg)) then
     ABI_FREE(kxcg)
   end if
   if (allocated(fxc_ADA)) then
     ABI_FREE(fxc_ADA)
   end if
   !
   ! Output the sum rule evaluation.
   ! Vcp%vc_sqrt(:,iqibz) Contains vc^{1/2}(q,G), complex-valued due to a possible cutoff
   epsm1 => chi0
   call output_chi0sumrule((is_qeq0==1),iqibz,Ep%npwe,omegaplasma,chi0_sumrule,epsm1(:,:,1),Vcp%vc_sqrt(:,iqibz))

   ! If input variable npvel is larger than 0, trigger the Random Stopping Power calculation
   ! Only the masternode is used
   if (my_rank==master .and. Dtset%npvel>0) then
     if (is_first_qcalc) then
       ABI_MALLOC(rspower,(Dtset%npvel))
       rspower(:)=zero
     end if
     call random_stopping_power(iqibz,Dtset%npvel,Dtset%pvelmax,Ep,Gsph_epsG0,Qmesh,Vcp,Cryst,Dtfil,epsm1,rspower)
     if (label==Ep%nqcalc) then
       ABI_FREE(rspower)
     end if
   end if

   ! Write heads and wings to main output file.
   if (is_qeq0==1) then
     write(msg,'(1x,2a)')' Heads and wings of the symmetrical epsilon^-1(G,G'') ',ch10
     call wrtout(ab_out,msg,'COLL')
     do iomega=1,Ep%nomega
       write(msg,'(2x,a,i4,a,2f9.4,a)')' Upper and lower wings at the ',iomega,' th omega',Ep%omega(iomega)*Ha_eV,' [eV]'
       call wrtout(ab_out,msg,'COLL')
       call print_arr(epsm1(1,:,iomega),max_r=9,unit=ab_out)
       call print_arr(epsm1(:,1,iomega),max_r=9,unit=ab_out)
       call wrtout(ab_out,ch10,'COLL')
     end do
   end if

   call timab(309,2,tsec)
   call timab(310,1,tsec) ! wrscr

   if (my_rank==master) then
     ! === Write the symmetrical epsilon^-1 on file ===
     title(1)='SCR file: epsilon^-1'
     if (is_qeq0==1) then
       string='0'; if (Dtset%usepaw==0.and.Ep%inclvkb/=0) call int2char10(Ep%inclvkb,string)
       title(1)=title(1)(1:21)//', calculated using inclvkb = '//string
     end if
     title(2)='TESTPARTICLE'
     ctype='RPA'
     title(2)(14:17)=ctype !this has to be modified

     if (is_first_qcalc) then
       ! === Open file and write the header for the SCR file ===
       ! * Here we write the RPA approximation for \tilde\epsilon^{-1}
       ikxc=0; test_type=0; tordering=1
       hem1 = hscr_new("inverse_dielectric_function",dtset,ep,hdr_local,ikxc,test_type,tordering,title,&
&       Ep%npwe,Gsph_epsG0%gvec)
       fform_em1 = hem1%fform
       if (dtset%iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
         NCF_CHECK(nctk_open_create(unt_em1, nctk_ncify(dtfil%fnameabo_scr), xmpi_comm_self))
         NCF_CHECK(crystal_ncwrite(cryst, unt_em1))
         NCF_CHECK(ebands_ncwrite(QP_BSt, unt_em1))
#endif
       else
         unt_em1=Dtfil%unscr
         if (open_file(dtfil%fnameabo_scr,msg,unit=unt_em1,status='unknown',form='unformatted') /= 0) then
           MSG_ERROR(msg)
         end if
       end if
       call hscr_io(hem1,fform_em1,2,unt_em1,xmpi_comm_self,0,Dtset%iomode)
       call hscr_free(Hem1)
     end if

     call write_screening("inverse_dielectric_function",unt_em1,Dtset%iomode,Ep%npwe,Ep%nomega,iqcalc,epsm1)
   end if ! my_rank==master

   call timab(310,2,tsec)
 end do ! Loop over q-points

 !----------------------------- END OF THE CALCULATION ------------------------

 ! Close Files.
 if (my_rank==master) then
   if (dtset%iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
     NCF_CHECK(nf90_close(unt_em1))
     if (dtset%prtsuscep>0) then
       NCF_CHECK(nf90_close(unt_susc))
     end if
#endif
   else
     close(unt_em1)
     if (dtset%prtsuscep>0) close(unt_susc)
   end if
 end if
!
!=====================
!==== Free memory ====
!=====================
 ABI_FREE(chi0_sumrule)
 ABI_FREE(chi0)
 if (allocated(rhor_kernel)) then
   ABI_FREE(rhor_kernel)
 end if

 ABI_FREE(rhor)
 ABI_FREE(rhog)
 ABI_FREE(ks_vbik)
 ABI_FREE(qp_vbik)
 ABI_FREE(ktabr)
 ABI_FREE(taur)
 ABI_FREE(taug)
 ABI_FREE(ks_vhartr)
 ABI_FREE(ks_vtrial)
 ABI_FREE(vpsp)
 ABI_FREE(ks_vxc)
 ABI_FREE(ph1d)
 ABI_FREE(ph1df)
 ABI_FREE(nhatgr)
 ABI_FREE(nhat)
 call pawfgr_destroy(Pawfgr)

 if (Dtset%usepaw==1) then ! Optional deallocation for PAW.
   call pawrhoij_free(Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   call paw_ij_free(Paw_ij)
   call paw_an_free(Paw_an)
   call pawpwff_free(Paw_pwff)
   if (Dtset%pawcross==1) then
     call paw_pwaves_lmn_free(Paw_onsite)
     Wfdf%bks_comm = xmpi_comm_null
     call wfd_free(Wfdf)
   end if
 end if

 ABI_DT_FREE(Pawfgrtab)
 ABI_DT_FREE(Paw_pwff)
 ABI_DT_FREE(Pawrhoij)
 ABI_DT_FREE(Paw_ij)
 ABI_DT_FREE(Paw_an)
 ABI_FREE(ktabrf)
 ABI_DT_FREE(Paw_onsite)

 call wfd_free(Wfd)
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call crystal_free(Cryst)
 call gsph_free(Gsph_epsG0)
 call gsph_free(Gsph_wfn)
 call vcoul_free(Vcp)
 call em1params_free(Ep)
 call hdr_free(Hdr_wfk)
 call hdr_free(Hdr_local)
 call ebands_free(KS_BSt)
 call ebands_free(QP_BSt)
 call littlegroup_free(Ltg_q)
 call destroy_mpi_enreg(MPI_enreg_seq)
 ABI_DT_FREE(Ltg_q)

 call timab(301,2,tsec)

 DBG_EXIT("COLL")

end subroutine screening
!!***

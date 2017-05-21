!{\src2tex{textfont=tt}}
!!****f* ABINIT/sigma
!! NAME
!! sigma
!!
!! FUNCTION
!! Calculate the matrix elements of the self-energy operator.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (GMR, VO, LR, RWG, MT, MG, RShaltaf)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! codvsn=code version
!! Dtfil<type(datafiles_type)>=variables related to files
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Pawang<type(pawang_type)>=paw angular mesh and related data
!! Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!! Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!! Psps<type(pseudopotential_type)>=variables related to pseudopotentials
!!   Before entering the first time in sigma, a significant part of Psps has been initialized :
!!   the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,mpsso,mgrid,ntypat,n1xccc,usepaw,useylm,
!!   and the arrays dimensioned to npsp. All the remaining components of Psps are to be initialized in
!!   the call to pspini. The next time the code enters screening, Psps might be identical to the
!!   one of the previous Dtset, in which case, no reinitialisation is scheduled in pspini.F90.
!! rprim(3,3)=dimensionless real space primitive translations
!!
!! OUTPUT
!!  Converged=.TRUE. if degw are within the user-specified tolerance.
!!  Output is written on the main abinit output file. Some results are stored in external files
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
!!      calc_sigc_me,calc_sigx_me,calc_ucrpa,calc_vhxc_me,chkpawovlp
!!      classify_bands,cohsex_me,crystal_free,denfgr,destroy_mpi_enreg
!!      ebands_copy,ebands_free,ebands_interpolate_kpath,ebands_report_gap
!!      ebands_update_occ,em1results_free,energies_init,esymm_free
!!      fftdatar_write,fourdp,get_gftt,getem1_from_ppm,getph,gsph_free,hdr_free
!!      init_distribfft_seq,initmpi_seq,kmesh_free,kxc_ada,kxc_driver
!!      littlegroup_free,littlegroup_init,melements_free,melements_print
!!      melements_zero,melflags_reset,metric,mkdump_er,mkrdim,nhatgrid
!!      paw_an_free,paw_an_init,paw_an_nullify,paw_check_symcprj,paw_dijhf
!!      paw_gencond,paw_ij_free,paw_ij_init,paw_ij_nullify,paw_ij_print
!!      paw_mkdijexc_core,paw_pwaves_lmn_free,paw_pwaves_lmn_init,paw_qpscgw
!!      pawcprj_alloc,pawcprj_free,pawdenpot,pawdij,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawfgrtab_print,pawinit,pawmknhat,pawprt
!!      pawpuxinit,pawpwff_free,pawpwff_init,pawrhoij_alloc,pawrhoij_copy
!!      pawrhoij_free,pawtab_get_lsize,pawtab_print,ppm_free,ppm_init
!!      prep_calc_ucrpa,print_ngfft,prtrhomxmn,pspini,rdgw,rdqps,read_rhor
!!      setsymrhoij,setup_ppmodel,setup_sigma,setvtr,show_qp,sigma_bksmask
!!      sigma_free,sigma_init,sigma_tables,sigparams_free,solve_dyson,symdij
!!      symdij_all,test_charge,timab,updt_m_lda_to_qp,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_free,wfd_get_cprj,wfd_init,wfd_mkrho
!!      wfd_print,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate,wfd_test_ortho
!!      write_sigma_header,write_sigma_results,wrqps,wrtout,xmpi_barrier
!!      xmpi_bcast,xmpi_sum
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine sigma(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,converged)

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_xomp
 use m_errors
 use m_profiling_abi
 use m_ab7_mixing
 use m_nctk
 use m_kxc
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use libxc_functionals
 use m_wfd

 use m_numeric_tools, only : imax_loc
 use m_fstrings,      only : strcat, sjoin, itoa
 use m_blas,          only : xdotc
 use m_io_tools,      only : open_file, file_exists, iomode_from_fname
 use m_mpinfo,        only : destroy_mpi_enreg
 use m_geometry,      only : normv
 use m_fftcore,       only : print_ngfft
 use m_fft_mesh,      only : get_gftt
 use m_ioarr,         only : fftdatar_write, read_rhor
 use m_crystal,       only : crystal_free, crystal_t
 use m_crystal_io,    only : crystal_ncwrite
 use m_ebands,        only : ebands_update_occ, ebands_copy, ebands_report_gap, get_valence_idx, get_bandenergy, &
&                            ebands_free, ebands_init, ebands_ncwrite, ebands_interpolate_kpath
 use m_energies,      only : energies_type, energies_init
 use m_bz_mesh,       only : kmesh_t, kmesh_free, littlegroup_t, littlegroup_init, littlegroup_free
 use m_gsphere,       only : gsphere_t, gsph_free
 use m_vcoul,         only : vcoul_t, vcoul_free
 use m_qparticles,    only : wrqps, rdqps, rdgw, show_QP, updt_m_lda_to_qp
 use m_screening,     only : mkdump_er, em1results_free, epsilonm1_results
 use m_ppmodel,       only : ppm_init, ppm_free, setup_ppmodel, getem1_from_PPm, ppmodel_t
 use m_sigma,         only : sigma_init, sigma_free, sigma_ncwrite, sigma_t, sigma_get_exene, &
&                            write_sigma_header, write_sigma_results
 use m_dyson_solver,  only : solve_dyson
 use m_esymm,         only : esymm_t, esymm_free, esymm_failed
 use m_melemts,       only : melflags_reset, melements_print, melements_free, melflags_t, melements_t, melements_zero
 use m_pawang,        only : pawang_type
 use m_pawrad,        only : pawrad_type
 use m_pawtab,        only : pawtab_type, pawtab_print, pawtab_get_lsize
 use m_paw_an,        only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,        only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_print
 use m_pawfgrtab,     only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free, pawfgrtab_print
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, symrhoij
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, paw_overlap
 use m_pawdij,        only : pawdij, symdij_all
 use m_pawfgr,        only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawpwij,       only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_paw_slater,    only : paw_mkdijexc_core, paw_dijhf
 use m_paw_dmft,      only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'sigma'
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
 logical,intent(out) :: converged
 character(len=6),intent(in) :: codvsn
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pawang_type),intent(inout) :: Pawang
 type(Pseudopotential_type),intent(inout) :: Psps
!arrays
 real(dp),intent(in) :: acell(3),rprim(3,3)
 type(Pawrad_type),intent(inout) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(inout) :: Pawtab(Psps%ntypat*Psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=40,tim_fourdp=5,master=0,cplex1=1
 integer :: approx_type,b1gw,b2gw,choice,cplex,cplex_dij,band
 integer :: dim_kxcg,gnt_option,has_dijU,has_dijso,iab,bmin,bmax,irr_idx1,irr_idx2
 integer :: iat,ib,ib1,ib2,ic,id_required,ider,idir,ii,ik,ierr,ount
 integer :: ik_bz,ikcalc,ik_ibz,ikxc,ipert,npw_k,omp_ncpus
 integer :: isp,is_idx,istep,itypat,itypatcor,izero,jj,first_band,last_band
 integer :: ks_iv,lcor,lmn2_size_max,mband,my_nband
 integer :: mgfftf,mod10,mod100,moved_atm_inside,moved_rhor,n3xccc !,mgfft
 integer :: nbsc,ndij,ndim,nfftf,nfftf_tot,nkcalc,gwc_nfft,gwc_nfftot,gwx_nfft,gwx_nfftot
 integer :: ngrvdw,nhatgrdim,nkxc,nkxc1,nprocs,nscf,nspden_rhoij,nzlmopt,optene
 integer :: optcut,optgr0,optgr1,optgr2,option,option_test,option_dij,optrad,optrhoij,psp_gencond
 integer :: my_rank,rhoxsp_method,comm,use_aerhor,use_umklp,usexcnhat
 integer :: ioe0j,spin,io,jb,nomega_sigc
 integer :: iomega,ppm_unt,temp_unt,ncid
 integer :: work_size,nstates_per_proc,my_nbks
 !integer :: jb_qp,ib_ks,ks_irr
 real(dp) :: compch_fft,compch_sph,r_s,rhoav,alpha,opt_ecut
 real(dp) :: drude_plsmf,my_plsmf,ecore,ecut_eff,ecutdg_eff,ehartree
 real(dp) :: ex_energy,gsqcutc_eff,gsqcutf_eff,gsqcut_shp,norm,oldefermi
 real(dp) :: ucvol,vxcavg,vxcavg_qp
 real(dp) :: gwc_gsq,gwx_gsq,gw_gsq
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem,ug_mem,ur_mem,cprj_mem
 complex(dpc) :: max_degw,cdummy
 logical :: use_paw_aeur,dbg_mode,pole_screening,call_pawinit
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname,pawden_fname,fname
 type(kmesh_t) :: Kmesh,Qmesh
 type(ebands_t) :: KS_BSt,QP_BSt
 type(vcoul_t) :: Vcp
 type(crystal_t) :: Cryst
 type(Energies_type) :: KS_energies,QP_energies
 type(Epsilonm1_results) :: Er
 type(gsphere_t) :: Gsph_Max,Gsph_x,Gsph_c
 type(hdr_type) :: Hdr_wfk,Hdr_sigma,hdr_rhor
 type(melflags_t) :: KS_mflags,QP_mflags
 type(melements_t) :: KS_me,QP_me
 type(MPI_type) :: MPI_enreg_seq
 type(paw_dmft_type) :: Paw_dmft
 type(pawfgr_type) :: Pawfgr
 type(ppmodel_t) :: PPm
 type(sigparams_t) :: Sigp
 type(sigma_t) :: Sr
 type(wfd_t),target :: Wfd,Wfdf
 type(wvl_data) :: Wvl
!arrays
 integer :: gwc_ngfft(18),ngfftc(18),ngfftf(18),gwx_ngfft(18)
 integer,allocatable :: nq_spl(:),nlmn_atm(:),my_spins(:)
 integer,allocatable :: tmp_gfft(:,:),ks_vbik(:,:),nband(:,:),l_size_atm(:),qp_vbik(:,:)
 integer,allocatable :: tmp_kstab(:,:,:),ks_irreptab(:,:,:),qp_irreptab(:,:,:),my_band_list(:)
 real(dp),parameter ::  k0(3)=zero
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: grchempottn(:,:),grewtn(:,:),grvdw(:,:),qmax(:)
 real(dp),allocatable :: ks_nhat(:,:),ks_nhatgr(:,:,:),ks_rhog(:,:)
 real(dp),allocatable :: ks_rhor(:,:),ks_vhartr(:),ks_vtrial(:,:),ks_vxc(:,:)
 real(dp),allocatable :: ks_taug(:,:),ks_taur(:,:)
 real(dp),allocatable :: kxc(:,:),qp_kxc(:,:),ph1d(:,:),ph1df(:,:)
 real(dp),allocatable :: prev_rhor(:,:),prev_taur(:,:),qp_nhat(:,:)
 real(dp),allocatable :: qp_nhatgr(:,:,:),qp_rhog(:,:),qp_rhor_paw(:,:)
 real(dp),allocatable :: qp_rhor_n_one(:,:),qp_rhor_nt_one(:,:)
 real(dp),allocatable :: qp_rhor(:,:),qp_vhartr(:),qp_vtrial(:,:),qp_vxc(:,:)
 real(dp),allocatable :: qp_taur(:,:),qp_taug(:,:),igwene(:,:,:)
 real(dp),allocatable :: vpsp(:),xccc3d(:),dijexc_core(:,:,:),dij_hf(:,:,:)
 complex(dpc),allocatable :: omega(:),em1_ppm(:,:,:)
 real(dp),allocatable :: ks_aepaw_rhor(:,:) !,ks_n_one_rhor(:,:),ks_nt_one_rhor(:,:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: ctmp(:,:),hbare(:,:,:,:)
 complex(dpc),target,allocatable :: sigcme(:,:,:,:,:)
 complex(dpc),allocatable :: hlda(:,:,:,:),htmp(:,:,:,:),uks2qp(:,:)
 complex(gwpc),allocatable :: kxcg(:,:),fxc_ADA(:,:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:)
 complex(dpc),pointer :: sigcme_p(:,:,:,:)
 complex(dpc), allocatable :: rhot1_q_m(:,:,:,:,:,:,:)
 complex(dpc), allocatable :: M1_q_m(:,:,:,:,:,:,:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:),bmask(:)
 type(esymm_t),target,allocatable :: KS_sym(:,:)
 type(esymm_t),pointer :: QP_sym(:,:)
 type(pawcprj_type),allocatable :: Cp1(:,:),Cp2(:,:)
 type(littlegroup_t),allocatable :: Ltg_k(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:),QP_paw_an(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:),QP_paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:),QP_pawrhoij(:),prev_Pawrhoij(:),tmp_pawrhoij(:)
 type(pawpwff_t),allocatable :: Paw_pwff(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)

!************************************************************************

 DBG_ENTER('COLL')

 call timab(401,1,tsec) ! sigma(Total)
 call timab(402,1,tsec) ! sigma(Init1)

 write(msg,'(7a)')&
& ' SIGMA: Calculation of the GW corrections ',ch10,ch10,&
& ' Based on a program developped by R.W. Godby, V. Olevano, G. Onida, and L. Reining.',ch10,&
& ' Incorporated in ABINIT by V. Olevano, G.-M. Rignanese, and M. Torrent.'
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 if(dtset%ucrpa>0) then
   write(msg,'(6a)')ch10,&
&   ' cRPA Calculation: Calculation of the screened Coulomb interaction (ucrpa/=0) ',ch10
   call wrtout(ab_out, msg,'COLL')
   call wrtout(std_out,msg,'COLL')
 end if

#if defined HAVE_GW_DPC
 if (gwpc/=8) then
   write(msg,'(6a)')ch10,&
&   'Number of bytes for double precision complex /=8 ',ch10,&
&   'Cannot continue due to kind mismatch in BLAS library ',ch10,&
&   'Some BLAS interfaces are not generated by abilint '
   MSG_ERROR(msg)
 end if
 write(msg,'(a,i2,a)')'.Using double precision arithmetic ; gwpc = ',gwpc,ch10
#else
 write(msg,'(a,i2,a)')'.Using single precision arithmetic ; gwpc = ',gwpc,ch10
#endif
 call wrtout(std_out,msg,'COLL')
 call wrtout(ab_out,msg,'COLL')

 ! Perform some additional checks for hybrid functional calculations
 if(Dtset%gwcalctyp>=100) then
   if (.not.libxc_functionals_check()) then
     msg='Hybrid functional calculations require the compilation with LIBXC library'
     MSG_ERROR(msg)
   end if
   if(MOD(Dtset%gwcalctyp,100)<10) then
     msg='gwcalctyp should enforce updated of the energies and/or wavefunctions when performing hybrid functional calculation'
     MSG_ERROR(msg)
   end if
   if(Dtset%usepaw==1) then
     msg='PAW version of hybrid functional calculations is not implemented'
     MSG_ERROR(msg)
   end if
   if(Dtset%gwcalctyp>=100 .AND. Dtset%gwcalctyp <200) then
     if( Dtset%rcut<tol6 ) then
       msg='For HSE calculation, the cutoff radius rcut has to be set to 9.090909 bohr!'
       MSG_ERROR(msg)
     end if
     if( Dtset%icutcoul /=5 .AND. Dtset%icutcoul /=15 ) then
       msg='For HSE calculation, abinit requires short-range only exchange (icutcoul=5)'
       MSG_ERROR(msg)
     end if
   end if
 end if

 !=== Initialize MPI variables, and parallelization level ===
 !* gwpara 0--> sequential run, 1--> parallelism over k-points, 2--> parallelism over bands.
 !* In case of gwpara==1 memory is not parallelized.
 !* If gwpara==2, bands are divided among processors but each proc has all the states where GW corrections are required.
 comm = xmpi_world; my_rank = xmpi_comm_rank(comm); nprocs = xmpi_comm_size(comm)
 converged = .FALSE.

 if (my_rank == master) then
   wfk_fname = dtfil%fnamewffk
   if (nctk_try_fort_or_ncfile(wfk_fname, msg) /= 0) then
     MSG_ERROR(msg)
   end if
 end if
 call xmpi_bcast(wfk_fname, master, comm, ierr)

 ! === Some variables need to be initialized/nullify at start ===
 call energies_init(KS_energies)
 usexcnhat=0
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,ab_out,rmet,rprimd,ucvol)
 !
 ! === Define FFT grid(s) sizes ===
 ! * Be careful! This mesh is only used for densities, potentials and the matrix elements of v_Hxc. It is NOT the
 ! (usually coarser) GW FFT mesh employed for the oscillator matrix elements that is defined in setmesh.F90.
 ! See also NOTES in the comments at the beginning of this file.
 ! NOTE: This mesh is defined in invars2m using ecutwfn, in GW Dtset%ecut is forced to be equal to Dtset%ecutwfn.

 call pawfgr_init(Pawfgr,Dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfftc,ngfftf,&
& gsqcutc_eff=gsqcutc_eff,gsqcutf_eff=gsqcutf_eff,gmet=gmet,k0=k0)

 ! Fake MPI_type for the sequential part.
 call initmpi_seq(MPI_enreg_seq)
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'c',ngfftc(2),ngfftc(3),'all')
 call init_distribfft_seq(MPI_enreg_seq%distribfft,'f',ngfftf(2),ngfftf(3),'all')

 call print_ngfft(ngfftf,header='Dense FFT mesh used for densities and potentials')
 nfftf_tot=PRODUCT(ngfftf(1:3))

 ! Open and read pseudopotential files ===
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,level,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)

 call timab(402,2,tsec) ! Init1
 !
 ! ===============================================
 ! ==== Initialize Sigp, Er and basic objects ====
 ! ===============================================
 ! * Sigp is completetly initialized here.
 ! * Er is only initialized with dimensions, (SCR|SUSC) file is read in mkdump_Er
 call timab(403,1,tsec) ! setup_sigma

 call setup_sigma(codvsn,wfk_fname,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&
& gwx_ngfft,gwc_ngfft,Hdr_wfk,Hdr_sigma,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_x,Gsph_c,Vcp,Er,Sigp,comm)

 call timab(403,2,tsec) ! setup_sigma
 call timab(402,1,tsec) ! Init1

!XG090617 Please, do not remove this write, unless you have checked
!that the code executes correctly on max+g95 (especially, Tv5#70).
!It is one more a silly write, perhaps needed because the compiler does not treat correctly non-nullified pointers.
 if ( sigma_needs_w(Sigp) .and. my_rank==master) then
   write(std_out,*)' screening after setup_sigma : Er%Hscr%headform=',Er%Hscr%headform
 end if
!END XG090617

 pole_screening = .FALSE.
 if (Er%fform==2002) then
   pole_screening = .TRUE.
   MSG_WARNING(' EXPERIMENTAL - Using a pole-fit screening!')
 end if

 call print_ngfft(gwc_ngfft,header='FFT mesh for oscillator strengths used for Sigma_c')
 call print_ngfft(gwx_ngfft,header='FFT mesh for oscillator strengths used for Sigma_x')

 mod10 =MOD(Sigp%gwcalctyp,10)
 mod100=MOD(Sigp%gwcalctyp,100)
 b1gw=Sigp%minbdgw
 b2gw=Sigp%maxbdgw

 gwc_nfftot=PRODUCT(gwc_ngfft(1:3))
 gwc_nfft  =gwc_nfftot  !no FFT //

 gwx_nfftot=PRODUCT(gwx_ngfft(1:3))
 gwx_nfft  =gwx_nfftot  !no FFT //
!
!TRYING TO RECREATE AN "ABINIT ENVIRONMENT"
 KS_energies%e_corepsp=ecore/Cryst%ucvol

!=== Calculate KS occupation numbers and ks_vbk(nkibz,nsppol) ====
!* ks_vbk gives the (valence|last Fermi band) index for each k and spin.
!* spinmagntarget is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 ABI_MALLOC(ks_vbik,(KS_BSt%nkpt,KS_BSt%nsppol))
 ABI_MALLOC(qp_vbik,(KS_BSt%nkpt,KS_BSt%nsppol))

 !call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=0)
 ks_vbik(:,:) = get_valence_idx(KS_BSt)

 ! ============================
 ! ==== PAW initialization ====
 ! ============================
 if (Dtset%usepaw==1) then
   call chkpawovlp(Cryst%natom,Cryst%ntypat,Dtset%pawovlp,Pawtab,Cryst%rmet,Cryst%typat,Cryst%xred)

   cplex_dij=Dtset%nspinor; cplex=1; ndij=1

   ABI_DT_MALLOC(KS_Pawrhoij,(Cryst%natom))
   nspden_rhoij=Dtset%nspden; if (Dtset%pawspnorb>0.and.Dtset%nspinor==2) nspden_rhoij=4
   call pawrhoij_alloc(KS_Pawrhoij,Dtset%pawcpxocc,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)

   ! Initialize values for several basic arrays
   gnt_option=1;if (dtset%pawxcdev==2.or.(dtset%pawxcdev==1.and.dtset%positron/=0)) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or. call_pawinit) then
     call timab(553,1,tsec)
     gsqcut_shp=two*abs(dtset%diecut)*dtset%dilatmx**2/pi**2
     call pawinit(gnt_option,gsqcut_shp,zero,Dtset%pawlcutd,Dtset%pawlmix,&
&     Psps%mpsang,Dtset%pawnphi,Cryst%nsym,Dtset%pawntheta,Pawang,Pawrad,&
&     Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,Dtset%xclevel,Dtset%usepotzero)
     call timab(553,2,tsec)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit)

   else
     if (Pawtab(1)%has_kij  ==1) Pawtab(1:Cryst%ntypat)%has_kij  =2
     if (Pawtab(1)%has_nabla==1) Pawtab(1:Cryst%ntypat)%has_nabla=2
   end if
   Psps%n1xccc=MAXVAL(Pawtab(1:Cryst%ntypat)%usetcore)

   ! Initialize optional flags in Pawtab to zero
   ! (Cannot be done in Pawinit since the routine is called only if some pars. are changed)
   Pawtab(:)%has_nabla = 0
   Pawtab(:)%usepawu   = 0
   Pawtab(:)%useexexch = 0
   Pawtab(:)%exchmix   =zero

   call setsymrhoij(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)

   ! Initialize and compute data for LDA+U
   Paw_dmft%use_dmft=Dtset%usedmft
   if (Dtset%usepawu>0.or.Dtset%useexexch>0) then
     call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla,&
&     Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,Cryst%ntypat,Pawang,Dtset%pawprtvol,&
&     Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu,dtset%ucrpa)
   end if

   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(Hdr_wfk%pawrhoij,KS_Pawrhoij)

   ! Re-symmetrize symrhoij.
   ! FIXME this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1; ipert=0; idir=0
!  call symrhoij(KS_Pawrhoij,KS_Pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
!  &             Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
!  &             Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

   !  Evaluate form factor of radial part of phi.phj-tphi.tphj.
   rhoxsp_method=1 ! Arnaud-Alouani
!  rhoxsp_method=2 ! Shiskin-Kresse
   if (Dtset%pawoptosc /= 0) rhoxsp_method = Dtset%pawoptosc

!  The q-grid must contain the FFT mesh used for sigma_c and the G-sphere for the exchange part.
!  We use the FFT mesh for sigma_c since COHSEX and the extrapolar method require oscillator
!  strengths on the FFT mesh.
   ABI_MALLOC(tmp_gfft,(3,gwc_nfftot))
   call get_gftt(gwc_ngfft,k0,gmet,gwc_gsq,tmp_gfft)
   ABI_FREE(tmp_gfft)

   gwx_gsq = Dtset%ecutsigx/(two*pi**2)
!  allocate(tmp_gfft(3,gwx_nfftot)); q0=zero
!  call get_gftt(gwx_ngfft,q0,gmet,gwx_gsq,tmp_gfft)
!  deallocate(tmp_gfft)
   gw_gsq = MAX(gwx_gsq,gwc_gsq)

!  * Set up q-grid, make qmax 20% larger than largest expected.
   ABI_MALLOC(nq_spl,(Psps%ntypat))
   ABI_MALLOC(qmax,(Psps%ntypat))
   qmax = SQRT(gw_gsq)*1.2d0 ! qmax=Psps%qgrid_ff(Psps%mqgrid_ff)
   nq_spl = Psps%mqgrid_ff
!  write(std_out,*)"using nq_spl",nq_spl,"qmax=",qmax
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
   cplex=1
   call pawfgrtab_init(Pawfgrtab,cplex,l_size_atm,Dtset%nspden,Dtset%typat)
   ABI_FREE(l_size_atm)
   compch_fft=greatest_real
   usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
   !  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
   !  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
   write(msg,'(a,i2)')' sigma : using usexcnhat = ',usexcnhat
   call wrtout(std_out,msg,'COLL')
   !
   ! Identify parts of the rectangular grid where the density has to be calculated ===
   optcut=0; optgr0=Dtset%pawstgylm; optgr1=0; optgr2=0; optrad=1-Dtset%pawstgylm
   if (Dtset%pawcross==1) optrad=1
   if (Dtset%xclevel==2.and.usexcnhat>0) optgr1=Dtset%pawstgylm

   call nhatgrid(Cryst%atindx1,gmet,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,Cryst%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,Pawfgrtab,Pawtab,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

   call pawfgrtab_print(Pawfgrtab,Cryst%natom,unit=std_out,prtvol=Dtset%pawprtvol)
 else
   ABI_DT_MALLOC(Paw_pwff,(0))
   ABI_DT_MALLOC(Pawfgrtab,(0))
 end if !End of PAW Initialization

 ! Consistency check and additional stuff done only for GW with PAW.
 ABI_DT_MALLOC(Paw_onsite,(0))

 if (Dtset%usepaw==1) then
   if (Dtset%ecutwfn < Dtset%ecut) then
     write(msg,"(5a)")&
&     "WARNING - ",ch10,&
&     "  It is highly recommended to use ecutwfn = ecut for GW calculations with PAW since ",ch10,&
&     "  an excessive truncation of the planewave basis set can lead to unphysical results."
     MSG_WARNING(msg)
   end if

   ABI_CHECK(Dtset%useexexch==0,"LEXX not yet implemented in GW")
   ABI_CHECK(Paw_dmft%use_dmft==0,"DMFT + GW not available")

   ! Optionally read core orbitals from file and calculate $ \<\phi_i|Sigma_x^\core|\phi_j\> $ for the HF decoupling.
   if (Sigp%use_sigxcore==1) then
     lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)
     ABI_MALLOC(dijexc_core,(cplex_dij*lmn2_size_max,ndij,Cryst%ntypat))

     call paw_mkdijexc_core(ndij,cplex_dij,lmn2_size_max,Cryst,Pawtab,Pawrad,dijexc_core,Dtset%prtvol,Psps%filpsp)
   end if ! HF decoupling

   if (Dtset%pawcross==1) then
     if (allocated(Paw_onsite) ) then
       ABI_DT_FREE(Paw_onsite)
     end if
     ABI_DT_MALLOC(Paw_onsite,(Cryst%natom))
     call paw_pwaves_lmn_init(Paw_onsite,Cryst%natom,Cryst%natom,Cryst%ntypat,&
&     Cryst%rprimd,Cryst%xcart,Pawtab,Pawrad,Pawfgrtab)
   end if
 end if

 ! Allocate these arrays anyway, since they are passed to subroutines.
 if (.not.allocated(ks_nhat)) then
   ABI_MALLOC(ks_nhat,(nfftf,0))
 end if
 if (.not.allocated(dijexc_core)) then
   ABI_MALLOC(dijexc_core,(1,1,0))
 end if

 ! ==================================================
 ! ==== Read KS band structure from the KSS file ====
 ! ==================================================
 !
 ! * Initialize Wfd, allocate wavefunctions and precalculate tables to do the FFT using the coarse gwc_ngfft.
 mband=Sigp%nbnds
 ABI_MALLOC(bks_mask,(mband,Kmesh%nibz,Sigp%nsppol))
 ABI_MALLOC(keep_ur ,(mband,Kmesh%nibz,Sigp%nsppol))
 keep_ur=.FALSE.; bks_mask=.FALSE.

 ABI_MALLOC(nband,(Kmesh%nibz,Sigp%nsppol))
 nband=mband

 ! autoparal section
 if (dtset%max_ncpus/=0) then
   ount =ab_out
    ! Temporary table needed to estimate memory
   ABI_MALLOC(nlmn_atm,(Cryst%natom))
   if (Dtset%usepaw==1) then
     do iat=1,Cryst%natom
       nlmn_atm(iat)=Pawtab(Cryst%typat(iat))%lmn_size
     end do
   end if

   write(ount,'(a)')"--- !Autoparal"
   write(ount,"(a)")'#Autoparal section for Sigma runs.'
   write(ount,"(a)")   "info:"
   write(ount,"(a,i0)")"    autoparal: ",dtset%autoparal
   write(ount,"(a,i0)")"    max_ncpus: ",dtset%max_ncpus
   write(ount,"(a,i0)")"    gwpara: ",dtset%gwpara
   write(ount,"(a,i0)")"    nkpt: ",dtset%nkpt
   write(ount,"(a,i0)")"    nsppol: ",dtset%nsppol
   write(ount,"(a,i0)")"    nspinor: ",dtset%nspinor
   write(ount,"(a,i0)")"    nbnds: ",Sigp%nbnds

   work_size = mband * Kmesh%nibz * Sigp%nsppol

   ! Non-scalable memory in Mb i.e. memory that is not distribute with MPI.
   nonscal_mem = (two*gwpc*Er%npwe**2*Er%nomega*(Er%mqmem+1)*b2Mb) * 1.1_dp

   ! List of configurations.
   ! Assuming an OpenMP implementation with perfect speedup!
   write(ount,"(a)")"configurations:"
   do ii=1,dtset%max_ncpus
     nstates_per_proc = 0
     eff = HUGE(one)
     max_wfsmem_mb = zero

     do my_rank=0,ii-1
       call sigma_bksmask(Dtset,Sigp,Kmesh,my_rank,ii,my_spins,bks_mask,keep_ur,ierr)
       ABI_FREE(my_spins)
       if (ierr /= 0) exit
       my_nbks = COUNT(bks_mask)
       nstates_per_proc = MAX(nstates_per_proc, my_nbks)
       eff = MIN(eff, (one * work_size) / (ii * nstates_per_proc))

       ! Memory needed for Fourier components ug.
       ug_mem = two*gwpc*Dtset%nspinor*Sigp%npwwfn*my_nbks*b2Mb

       ! Memory needed for real space ur (use gwc_nfft, instead of gwx_nfft)
       ur_mem = two*gwpc*Dtset%nspinor*gwc_nfft*COUNT(keep_ur)*b2Mb

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
   call sigma_bksmask(Dtset,Sigp,Kmesh,my_rank,nprocs,my_spins,bks_mask,keep_ur,ierr)
   ABI_CHECK(ierr==0, "Error in sigma_bksmask")
 end if

 ! Then each node owns the wavefunctions where GW corrections are required.
 do isp=1,SIZE(my_spins)
   spin = my_spins(isp)
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Irred k-point for GW
     ii=Sigp%minbnd(ikcalc,spin); jj=Sigp%maxbnd(ikcalc,spin)
     bks_mask(ii:jj,ik_ibz,spin) = .TRUE.
     if (MODULO(Dtset%gwmem,10)==1) keep_ur(ii:jj,ik_ibz,spin)=.TRUE.
   end do
 end do

 ABI_FREE(my_spins)
 opt_ecut=Dtset%ecutwfn
!opt_ecut=zero

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Sigp%npwwfn,mband,nband,Kmesh%nibz,Sigp%nsppol,bks_mask,&
& Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,gwc_ngfft,&
& Gsph_Max%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)

 if (Dtset%pawcross==1) then
   call wfd_init(Wfdf,Cryst,Pawtab,Psps,keep_ur,Dtset%paral_kgb,Sigp%npwwfn,mband,nband,Kmesh%nibz,Sigp%nsppol,bks_mask,&
&   Dtset%nspden,Dtset%nspinor,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,gwc_ngfft,&
&   Gsph_Max%gvec,Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm,opt_ecut=opt_ecut)
 end if

 ABI_FREE(bks_mask)
 ABI_FREE(nband)
 ABI_FREE(keep_ur)

 call timab(402,2,tsec) ! sigma(Init1)

 call timab(404,1,tsec) ! rdkss

 call wfd_read_wfk(Wfd,wfk_fname,iomode_from_fname(wfk_fname))

 if (Dtset%pawcross==1) then
   call wfd_copy(Wfd,Wfdf)
   call wfd_change_ngfft(Wfdf,Cryst,Psps,ngfftf)
 end if

 ! This test has been disabled (too expensive!)
 if (.False.) call wfd_test_ortho(Wfd,Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(404,2,tsec) ! rdkss

 call timab(405,1,tsec) ! Init2

!Debugging section.
!if (.TRUE.) then
 if (.FALSE.) then
!
   if (.FALSE..and.Wfd%usepaw==1) then
     ABI_DT_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
     ABI_DT_MALLOC(Cp2,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp2,0,Wfd%nlmn_atm)

     call wfd_change_ngfft(Wfd,Cryst,Psps,ngfftf)

     do spin=1,Wfd%nsppol
       do ik_bz=1,Kmesh%nbz
         ik_ibz=Kmesh%tab(ik_bz)
         do band=1,Wfd%nband(ik_ibz,spin)
           call paw_check_symcprj(Wfd,ik_bz,band,spin,1,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp1)
           call paw_check_symcprj(Wfd,ik_bz,band,spin,2,Cryst,Kmesh,Psps,Pawtab,Pawang,Cp2)

           do iat=1,Cryst%natom
             do isp=1,Wfd%nspinor
               write(789,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp1(iat,isp)%cp
               write(790,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp2(iat,isp)%cp
               write(791,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp1(iat,isp)%cp(1,:)**2 + Cp1(iat,isp)%cp(2,:)**2
               write(792,'(3i2,/,(f8.4))') band,ik_bz,spin,Cp2(iat,isp)%cp(1,:)**2 + Cp2(iat,isp)%cp(2,:)**2
             end do
           end do
         end do
       end do
     end do

     call pawcprj_free(Cp1)
     ABI_DT_FREE(Cp1)
     call pawcprj_free(Cp2)
     ABI_DT_FREE(Cp2)
   end if

 end if

 ! ==============================================================
 ! ==== Find little group of the k-points for GW corrections ====
 ! ==============================================================
 ! * The little group is calculated only if sys_sigma.
 ! * If use_umklp==1 then also symmetries requiring an umklapp to preserve k_gw are included.
 !
 ABI_DT_MALLOC(Ltg_k,(Sigp%nkptgw))
 use_umklp=1
 do ikcalc=1,Sigp%nkptgw
   if (Sigp%symsigma/=0) then
     call littlegroup_init(Sigp%kptgw(:,ikcalc),Qmesh,Cryst,use_umklp,Ltg_k(ikcalc),0)
   end if
 end do
!
!=== Compute structure factor phases and large sphere cut-off ===
!WARNING cannot use Dtset%mgfft, this has to be checked better
!mgfft=MAXVAL(ngfftc(:))
!allocate(ph1d(2,3*(2*mgfft+1)*Cryst%natom),ph1df(2,3*(2*mgfftf+1)*Cryst%natom))
!write(std_out,*)' CHECK ',Dtset%mgfftdg,mgfftf
!if (Dtset%mgfftdg/=mgfftf) then
!write(std_out,*)"WARNING Dtset%mgfftf /= mgfftf"
!write(std_out,*)'HACKING Dtset%mgfftf'
!Dtset%mgfftdg=mgfftf
!end if
 ABI_MALLOC(ph1d,(2,3*(2*Dtset%mgfft+1)*Cryst%natom))
 ABI_MALLOC(ph1df,(2,3*(2*mgfftf+1)*Cryst%natom))

 call getph(Cryst%atindx,Cryst%natom,ngfftc(1),ngfftc(2),ngfftc(3),ph1d,Cryst%xred)

 if (Psps%usepaw==1.and.Pawfgr%usefinegrid==1) then
   call getph(Cryst%atindx,Cryst%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,Cryst%xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

 !===================================================================================
 !==== Classify the GW wavefunctions according to the irreducible representation ====
 !===================================================================================
 !* Warning still under development.
 !* Only for SCGW.
 !bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw

 ABI_DT_MALLOC(KS_sym,(Wfd%nkibz,Wfd%nsppol))

 if (Sigp%symsigma==1.and.mod100>=20) then
   ! call check_zarot(Gsph_c%ng,Cryst,gwc_ngfft,Gsph_c%gvec,Psps,Pawang,Gsph_c%rottb,Gsph_c%rottbm1)
   use_paw_aeur=.FALSE. ! should pass ngfftf but the dense mesh is not forced to be symmetric
   do spin=1,Wfd%nsppol
     do ikcalc=1,Sigp%nkptgw
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       first_band = Sigp%minbnd(ikcalc,spin)
       last_band  = Sigp%maxbnd(ikcalc,spin)
!      call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,Cryst,KS_BSt,Pawtab,Pawrad,Pawang,Psps,&
       call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,Wfd%ngfft,Cryst,KS_BSt,Pawtab,Pawrad,Pawang,Psps,&
&       Dtset%tolsym,KS_sym(ik_ibz,spin))
     end do
   end do
   ! Recreate the Sig_ij tables taking advantage of the classification of the bands.
   call sigma_tables(Sigp,Kmesh,KS_sym)
 end if

 call timab(405,2,tsec) ! Init2

 call timab(406,1,tsec) ! make_vhxc
 !
 !===========================
 !=== COMPUTE THE DENSITY ===
 !===========================
 ! * Evaluate the planewave part (complete charge in case of NC pseudos).
 ABI_MALLOC(ks_rhor,(nfftf,Dtset%nspden))
 ABI_MALLOC(ks_taur,(nfftf,Dtset%nspden*Dtset%usekden))

 call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)

 if (Dtset%usekden==1) then
   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_taur,optcalc=1)
 end if
 !
 !========================================
 !==== Additional computation for PAW ====
 !========================================
 nhatgrdim=0
 if (Dtset%usepaw==1) then
   ! Calculate the compensation charge nhat.
   if (Dtset%xclevel==2) nhatgrdim=usexcnhat*Dtset%pawnhatxc
   cplex=1; ider=2*nhatgrdim; izero=0
   if (nhatgrdim>0)  then
     ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,3*nhatgrdim))
   end if
   if (nhatgrdim==0)  then
     ABI_MALLOC(ks_nhatgr,(0,0,0))
   end if

   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
&   Pawfgrtab,ks_nhatgr,ks_nhat,KS_Pawrhoij,KS_Pawrhoij,Pawtab,k0,Cryst%rprimd,&
&   Cryst%ucvol,dtset%usewvl,Cryst%xred)

   !if (nhatgrdim==0)  then
   !  ABI_FREE(ks_nhatgr)
   !end if

   ! === Evaluate onsite energies, potentials, densities ===
   !  * Initialize variables/arrays related to the PAW spheres.
   !  * Initialize also lmselect (index of non-zero LM-moments of densities).
   ABI_DT_MALLOC(KS_paw_ij,(Cryst%natom))
!  cplex=1;cplex_dij=Dtset%nspinor
   has_dijso=Dtset%pawspnorb; has_dijU=Dtset%usepawu
   call paw_ij_nullify(KS_paw_ij)
   call paw_ij_init(KS_paw_ij,cplex,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=1,&
&   has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   nkxc1=0
   ABI_DT_MALLOC(KS_paw_an,(Cryst%natom))
   call paw_an_nullify(KS_paw_an)
   call paw_an_init(KS_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,Dtset%nspden,cplex,Dtset%pawxcdev,&
&   Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=1)
!
!  Calculate onsite vxc with and without core charge.
   nzlmopt=-1; option=0; compch_sph=greatest_real
   call pawdenpot(compch_sph,KS_energies%e_paw,KS_energies%e_pawdc,ipert,&
&   Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
&   Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,KS_Paw_an,KS_Paw_an,KS_paw_ij,&
&   Pawang,Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,&
&   Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)

! TO BE REMOVED ASAP !!!
   _IBM6("Silly write for IBM6")
 else
   ABI_MALLOC(ks_nhatgr,(0,0,0))
   ABI_DT_MALLOC(KS_paw_ij,(0))
   ABI_DT_MALLOC(KS_paw_an,(0))
 end if !PAW

 !if (.not.allocated(ks_nhatgr))  then
 !  ABI_MALLOC(ks_nhatgr,(nfftf,Dtset%nspden,0))
 !end if

 call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,ks_rhor,Cryst%ucvol,&
& Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)
!
!For PAW, add the compensation charge on the FFT mesh, then get rho(G).
 if (Dtset%usepaw==1) ks_rhor=ks_rhor+ks_nhat

 call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_rhor,ucvol=ucvol)
 if(Dtset%usekden==1)then
   call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,ks_taur,optrhor=1,ucvol=ucvol)
 end if

 ABI_MALLOC(ks_rhog,(2,nfftf))
 ABI_MALLOC(ks_taug,(2,nfftf*Dtset%usekden))
 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 if(Dtset%usekden==1)then
   call fourdp(1,ks_taug,ks_taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
 end if

 !
 !The following steps have been gathered in the setvtr routine:
 !- get Ewald energy and Ewald forces
 !- compute local ionic pseudopotential vpsp
 !- eventually compute 3D core electron density xccc3d
 !- eventually compute vxc and vhartr
 !- set up ks_vtrial
 !
 !*******************************************************************
 !**** NOTE THAT HERE Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
 !*******************************************************************

 ngrvdw=0
 ABI_MALLOC(grvdw,(3,ngrvdw))
 ABI_MALLOC(grchempottn,(3,Cryst%natom))
 ABI_MALLOC(grewtn,(3,Cryst%natom))
 nkxc=0
 if (Dtset%nspden==1) nkxc=2
 if (Dtset%nspden>=2) nkxc=3 ! check GGA and spinor, quite a messy part!!!
!In case of MGGA, fxc and kxc are not available and we dont need them
!for the screening part (for now ...)
 if (Dtset%ixc<0.and.libxc_functionals_ismgga()) nkxc=0
 if (nkxc/=0)  then
   ABI_MALLOC(kxc,(nfftf,nkxc))
 end if

 n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
 ABI_MALLOC(xccc3d,(n3xccc))
 ABI_MALLOC(ks_vhartr,(nfftf))
 ABI_MALLOC(ks_vtrial,(nfftf,Dtset%nspden))
 ABI_MALLOC(vpsp,(nfftf))
 ABI_MALLOC(ks_vxc,(nfftf,Dtset%nspden))

 optene=4; moved_atm_inside=0; moved_rhor=0; istep=1

 call setvtr(Cryst%atindx1,Dtset,KS_energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,&
& istep,kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
& Cryst%nattyp,nfftf,ngfftf,ngrvdw,ks_nhat,ks_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
& optene,pawrad,Pawtab,ph1df,Psps,ks_rhog,ks_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
& Cryst%ucvol,usexcnhat,ks_vhartr,vpsp,ks_vtrial,ks_vxc,vxcavg,Wvl,xccc3d,Cryst%xred,taug=ks_taug,taur=ks_taur)

!============================
!==== Compute KS PAW Dij ====
!============================
 if (Dtset%usepaw==1) then
   !TO BE REMOVED
   _IBM6("Another silly write for IBM6")
   call timab(561,1,tsec)

   !  Calculate the unsymmetrized Dij.
   ipert=0; idir=0
   call pawdij(cplex,Dtset%enunit,Cryst%gprimd,ipert,&
&   Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
&   Dtset%nspden,Cryst%ntypat,KS_paw_an,KS_paw_ij,Pawang,Pawfgrtab,&
&   Dtset%pawprtvol,Pawrad,KS_Pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&   k0,Dtset%spnorbscl,Cryst%ucvol,dtset%charge,ks_vtrial,ks_vxc,Cryst%xred,&
&   nucdipmom=Dtset%nucdipmom)

    !  Symmetrize KS Dij
#if 0
   call symdij(Cryst%gprimd,Cryst%indsym,ipert,&
&   Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,0,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
#else
   call symdij_all(Cryst%gprimd,Cryst%indsym,ipert,&
&   Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,KS_paw_ij,Pawang,&
&   Dtset%pawprtvol,Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec)
#endif
!
!  Output the pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   call pawprt(Dtset,Cryst%natom,KS_paw_ij,KS_Pawrhoij,Pawtab)
   call timab(561,2,tsec)
 end if

 call timab(406,2,tsec) ! make_vhxc

 !=== Calculate Vxc(b1,b2,k,s)=<b1,k,s|v_{xc}|b2,k,s>  for all the states included in GW ===
 !  * ks_vxcvalme is calculated without NLCC, ks_vxcme contains NLCC (if any)
 !  * This part is parallelized within MPI_COMM_WORD since each node has all GW wavefunctions.
 !  * ks_vUme is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
 !  * Note that vH matrix elements are calculated using the true uncutted interaction.
 call timab(407,1,tsec) ! vHxc_me

 call melflags_reset(KS_mflags)
 KS_mflags%has_vhartree=1
 KS_mflags%has_vxc     =1
 KS_mflags%has_vxcval  =1
 if (Dtset%usepawu>0     )  KS_mflags%has_vu     =1
 if (Dtset%useexexch>0   )  KS_mflags%has_lexexch=1
 if (Sigp%use_sigxcore==1)  KS_mflags%has_sxcore =1
 if (mod100<10           )  KS_mflags%only_diago =1 ! off-diagonal elements only for SC on wavefunctions.

 if (.FALSE.) then ! quick and dirty hack to test HF contribution.
   MSG_WARNING("testing on-site HF")
   lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)
   ABI_MALLOC(dij_hf,(cplex_dij*lmn2_size_max,ndij,Cryst%natom))
   call paw_dijhf(ndij,cplex_dij,lmn2_size_max,Cryst%natom,Cryst%ntypat,Pawtab,Pawrad,Pawang,&
&   KS_Pawrhoij,dij_hf,Dtset%prtvol)

   do iat=1,Cryst%natom
     itypat = Cryst%typat(iat)
     ii = Pawtab(itypat)%lmn2_size
     KS_Paw_ij(iat)%dijxc(:,:) = dij_hf(1:cplex_dij*ii,:,iat)
     KS_Paw_ij(iat)%dijxc_hat(:,:) = zero
   end do
   ABI_FREE(dij_hf)

!  option_dij=3
!  call symdij(Cryst%gprimd,Cryst%indsym,ipert,&
!  &   Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
!  &   KS_paw_ij,Pawang,Dtset%pawprtvol,Pawtab,Cryst%rprimd,&
!  &   Cryst%symafm,Cryst%symrec)
 end if

 ABI_MALLOC(tmp_kstab,(2,Wfd%nkibz,Wfd%nsppol))
 tmp_kstab=0
 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw ! No spin dependent!
     ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
     tmp_kstab(1,ik_ibz,spin)=Sigp%minbnd(ikcalc,spin)
     tmp_kstab(2,ik_ibz,spin)=Sigp%maxbnd(ikcalc,spin)
   end do
 end do

 if (Dtset%gwfockmix < 0.0_dp .or. (Dtset%gwfockmix-1.0_dp) > tol8) then
   MSG_ERROR('gwfockmix is invalid.')
 end if

 call calc_vhxc_me(Wfd,KS_mflags,KS_me,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&
& ks_vtrial,ks_vhartr,ks_vxc,Psps,Pawtab,KS_paw_an,Pawang,Pawfgrtab,KS_paw_ij,dijexc_core,&
& ks_rhor,ks_rhog,usexcnhat,ks_nhat,ks_nhatgr,nhatgrdim,tmp_kstab,taug=ks_taug,taur=ks_taur)
 ABI_FREE(tmp_kstab)

!#ifdef DEV_HAVE_SCGW_SYM
!Set KS matrix elements connecting different irreps to zero. Do not touch unknown bands!.
 if (mod100>=20 .and. Sigp%symsigma > 0) then
   bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
   ABI_MALLOC(ks_irreptab,(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
   ks_irreptab=0
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       first_band = Sigp%minbnd(ikcalc,spin)
       last_band  = Sigp%maxbnd(ikcalc,spin)
       if (.not.esymm_failed(KS_sym(ik_ibz,spin))) then
         ks_irreptab(first_band:last_band,ik_ibz,spin) = KS_sym(ik_ibz,spin)%b2irrep(first_band:last_band)
!        ks_irreptab(bmin:bmax,ik_ibz,spin) = KS_sym(ik_ibz,spin)%b2irrep(bmin:bmax)
       end if
     end do
   end do
   call melements_zero(KS_me,ks_irreptab)
   ABI_FREE(ks_irreptab)
 end if
!#endif

 call melements_print(KS_me,header="Matrix elements in the KS basis set",prtvol=Dtset%prtvol)
 !
 ! If possible, calculate the EXX energy from the between the frozen core
 ! and the valence electrons using KS wavefunctions
 ! MG: BE careful here, since ex_energy is meaningful only is all occupied states are calculated.
 if( KS_mflags%has_sxcore ==1 ) then
   ! TODO
   !ex_energy = mels_get_exene_core(KS_me,kmesh,KS_BSt)
   ex_energy=zero
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       do ib=b1gw,b2gw
         if (Sigp%nsig_ab==1) then
           ex_energy = ex_energy + half*KS_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*KS_me%sxcore(ib,ib,ik,spin)
         else
           ex_energy = ex_energy + half*KS_BSt%occ(ib,ik,spin)*Kmesh%wt(ik)*SUM(KS_me%sxcore(ib,ib,ik,:))
         end if
       end do
     end do
   end do
   write(msg,'(a,2(es16.6,a))')' CORE Exchange energy with KS wavefunctions: ',ex_energy,' Ha ,',ex_energy*Ha_eV,' eV'
   call wrtout(std_out,msg,'COLL')
 end if

 call timab(407,2,tsec) ! vHxc_me

 call timab(408,1,tsec) ! hqp_init

 ! Do not break this coding! When gwcalctyp>10, the order of the bands can be interexchanged after
 ! the diagonalization. Therefore, we have to correctly assign the matrix elements to the corresponding
 ! bands and we cannot skip the following even though it looks unuseful.
 if (mod100>=10) then
   call wrtout(std_out,ch10//' *************** KS Energies *******************','COLL')
 end if

 !=== QP_BSt stores energies and occ. used for the calculation ===
 !  * Initialize QP_BSt with KS values.
 !  * In case of SC update QP_BSt using the QPS file.
 call ebands_copy(KS_BSt,QP_BSt)

 ABI_MALLOC(qp_rhor,(nfftf,Dtset%nspden))
 ABI_MALLOC(qp_taur,(nfftf,Dtset%nspden*Dtset%usekden))
 QP_sym => KS_sym

 if (mod100<10) then  ! one-shot GW, just do a copy of the KS density.
   qp_rhor=ks_rhor
   if(Dtset%usekden==1)qp_taur=ks_taur
   QP_sym => KS_sym
 else
   ! Self-consistent GW.
   !  * Read the unitary matrix and the QP energies of the previous step from the QPS file.
   call energies_init(QP_energies)
   QP_energies%e_corepsp=ecore/Cryst%ucvol

   !  m_lda_to_qp(ib,jb,k,s) := <\psi_{ib,k,s}^{KS}|\psi_{jb,k,s}^{QP}>
   ! Initialize the QP amplitudes with KS wavefunctions.
   ABI_MALLOC(Sr%m_lda_to_qp,(Sigp%nbnds,Sigp%nbnds,Kmesh%nibz,Sigp%nsppol))
   Sr%m_lda_to_qp=czero
   do ib=1,Sigp%nbnds
     Sr%m_lda_to_qp(ib,ib,:,:)=cone
   end do

   ! Now read m_lda_to_qp and update the energies in QP_BSt.
   ! TODO switch on the renormalization of n in sigma.
   ABI_MALLOC(prev_rhor,(nfftf,Dtset%nspden))
   ABI_MALLOC(prev_taur,(nfftf,Dtset%nspden*Dtset%usekden))
   ABI_DT_MALLOC(prev_Pawrhoij,(Cryst%natom*Psps%usepaw))

   call rdqps(QP_BSt,Dtfil%fnameabi_qps,Dtset%usepaw,Dtset%nspden,1,nscf,&
&   nfftf,ngfftf,Cryst%ucvol,Dtset%paral_kgb,Cryst,Pawtab,MPI_enreg_seq,nbsc,Sr%m_lda_to_qp,prev_rhor,prev_Pawrhoij)

!  Find the irreps associated to the QP amplitudes starting from the analogous table for the KS states.
!  bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
!  allocate(qp_irreptab(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
!  qp_irreptab=0
!  !qp_irreptab=ks_irreptab

!  do jb_qp=bmin,bmax
!  do ib_ks=bmin,bmax
!  if (ABS(Sr%m_lda_to_qp(ib_ks,jb_qp,ik_ibz,spin)) > tol12) then ! jb_qp has same the same character as ib_ks.
!  ks_irr = ks_irreptab(ib_ks,ib_ks,ik_ibz,spin)
!  qp_irreptab(jb_qp,jb_qp,ik_ibz,spin) = ks_irr
!  do ii=bmin,bmax
!  if (ks_irr == ks_irreptab(ii,ib_ks,ik_ibz,spin)) then
!  qp_irreptab(jb_qp,ii,ik_ibz,spin) = ks_irr
!  end if
!  end do
!  end if
!  end do
!  end do

   if (nscf==0) prev_rhor=ks_rhor
   if (nscf==0 .and. Dtset%usekden==1) prev_taur=ks_taur

   if (nscf>0.and.mod100>=20.and.wfd_iam_master(Wfd)) then ! Print the unitary transformation on std_out.
     call show_QP(QP_BSt,Sr%m_lda_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
   end if

   !  === Compute QP wfg as linear combination of KS states ===
   !  * Wfd%ug is modified inside calc_wf_qp
   !  * For PAW, update also the on-site projections.
   !  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
   !  TODO here we should use nbsc instead of nbnds

   call wfd_rotate(Wfd,Cryst,Sr%m_lda_to_qp)

   ! * Reinit the storage mode of Wfd as ug have been changed ===
   ! * Update also the wavefunctions for GW corrections on each processor
   call wfd_reset_ur_cprj(Wfd)

   ! Compute QP occupation numbers.
   call wrtout(std_out,'sigma: calculating QP occupation numbers:','COLL')

   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=0)
   qp_vbik(:,:) = get_valence_idx(QP_BSt)

!  #ifdef DEV_HAVE_SCGW_SYM
   ! Calculate the irreducible representations of the new QP amplitdues.
   if (Sigp%symsigma==1.and.mod100>=20) then
     ABI_DT_MALLOC(QP_sym,(Wfd%nkibz,Wfd%nsppol))
     use_paw_aeur=.FALSE. ! should pass ngfftf but the dense mesh is not forced to be symmetric
     do spin=1,Wfd%nsppol
       do ikcalc=1,Sigp%nkptgw
         ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
!        Quick fix for SCGW+symm TODO fix properly!
         first_band = Sigp%minbnd(ikcalc,spin)
         last_band  = Sigp%maxbnd(ikcalc,spin)
!        first_band = MINVAL(Sigp%minbnd(:,spin))
!        last_band  = MAXVAL(Sigp%maxbnd(:,spin))
!        call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,ngfftf,Cryst,QP_BSt,Pawtab,Pawrad,Pawang,Psps,&
         call classify_bands(Wfd,use_paw_aeur,first_band,last_band,ik_ibz,spin,Wfd%ngfft,Cryst,QP_BSt,Pawtab,Pawrad,Pawang,Psps,&
&         Dtset%tolsym,QP_sym(ik_ibz,spin))
       end do
     end do

     ! Recreate the Sig_ij tables taking advantage of the classification of the bands.
     call sigma_tables(Sigp,Kmesh,QP_sym)
   end if
!  #endif

   ! Compute QP density using the updated wfg.
   call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
   if(Dtset%usekden==1) then
     call wfd_mkrho(Wfd,Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_taur,optcalc=1)
   end if

   ! ========================================
   ! ==== QP self-consistent GW with PAW ====
   ! ========================================
   if (Dtset%usepaw==1) then

     ABI_MALLOC(qp_nhat,(nfftf,Dtset%nspden))
     nhatgrdim=0; if (Dtset%xclevel==2) nhatgrdim=usexcnhat
     ABI_MALLOC(qp_nhatgr,(nfftf,Dtset%nspden,3*nhatgrdim))

     ABI_DT_MALLOC(QP_pawrhoij,(Cryst%natom))
     ABI_DT_MALLOC(QP_paw_ij,(Cryst%natom))
     ABI_DT_MALLOC(QP_paw_an,(Cryst%natom))
!
!    Calculate new QP quantities: nhat, nhatgr, rho_ij, paw_ij, and paw_an.
     call paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&
&     Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&
&     QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,compch_sph,compch_fft)
   else
     ABI_MALLOC(qp_nhat,(0,0))
     ABI_MALLOC(qp_nhatgr,(0,0,0))
     ABI_DT_MALLOC(QP_pawrhoij,(0))
     ABI_DT_MALLOC(QP_paw_ij,(0))
     ABI_DT_MALLOC(QP_paw_an,(0))
   end if

!  here I should renormalize the density
   call test_charge(nfftf,KS_BSt%nelect,Dtset%nspden,qp_rhor,Cryst%ucvol,&
&   Dtset%usepaw,usexcnhat,Pawfgr%usefinegrid,compch_sph,compch_fft,drude_plsmf)

   if (Dtset%usepaw==1) qp_rhor(:,:)=qp_rhor(:,:)+qp_nhat(:,:) ! Add the "hat" term.

   call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,qp_rhor,ucvol=ucvol)
   if(Dtset%usekden==1) then
     call prtrhomxmn(std_out,MPI_enreg_seq,nfftf,ngfftf,Dtset%nspden,1,qp_taur,optrhor=1,ucvol=ucvol)
   end if

   ! Simple mixing of the PW density to damp oscillations in the Hartree potential.
   if (nscf>0 .and. (ABS(Dtset%rhoqpmix-one)>tol12) ) then
     write(msg,'(2a,f6.3)')ch10,' sigma: mixing QP densities using rhoqpmix= ',Dtset%rhoqpmix
     call wrtout(std_out,msg,'COLL')
     qp_rhor = prev_rhor + Dtset%rhoqpmix*(qp_rhor-prev_rhor)
     if(Dtset%usekden==1)qp_taur = prev_taur + Dtset%rhoqpmix*(qp_taur-prev_taur) ! mix taur.
   end if

   ABI_FREE(prev_rhor)
   ABI_FREE(prev_taur)
   if (Psps%usepaw==1.and.nscf>0) then
     call pawrhoij_free(prev_pawrhoij)
   end if
   ABI_DT_FREE(prev_pawrhoij)

   ABI_MALLOC(qp_rhog,(2,nfftf))
   ABI_MALLOC(qp_taug,(2,nfftf*Dtset%usekden))
   call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)
   if(Dtset%usekden==1)call fourdp(1,qp_taug,qp_taur(:,1),-1,MPI_enreg_seq,nfftf,ngfftf,Dtset%paral_kgb,tim_fourdp)

   ! ===========================================
   ! ==== Optional output of the QP density ====
   ! ===========================================
   if (Dtset%prtden/=0.and.wfd_iam_master(Wfd)) then
     call fftdatar_write("qp_rhor",dtfil%fnameabo_qp_den,dtset%iomode,hdr_sigma,&
     cryst,ngfftf,cplex1,nfftf,dtset%nspden,qp_rhor,mpi_enreg_seq,ebands=qp_bst)
   end if

   ! ===========================================
   ! === Optional output of the full QP density
   ! ===========================================
   if (Wfd%usepaw==1.and.Dtset%prtden==2) then
     ABI_MALLOC(qp_rhor_paw   ,(nfftf,Wfd%nspden))
     ABI_MALLOC(qp_rhor_n_one ,(nfftf,Wfd%nspden))
     ABI_MALLOC(qp_rhor_nt_one,(nfftf,Wfd%nspden))

     call denfgr(Cryst%atindx1,Cryst%gmet,comm,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,qp_nhat,&
&     Wfd%nspinor,Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,QP_pawrhoij,Pawtab,Dtset%prtvol,&
&     qp_rhor,qp_rhor_paw,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

     ABI_FREE(qp_rhor_n_one)
     ABI_FREE(qp_rhor_nt_one)
     if (Dtset%prtvol>9) then ! Print a normalisation check
       norm = SUM(qp_rhor_paw(:,1))*Cryst%ucvol/PRODUCT(Pawfgr%ngfft(1:3))
       write(msg,'(a,F8.4)') '  QUASIPARTICLE DENSITY CALCULATED - NORM OF DENSITY: ',norm
       call wrtout(std_out,msg,'PERS')
     end if

     ! Write the density to file
     if (my_rank==master) then
       call fftdatar_write("qp_pawrhor",dtfil%fnameabo_qp_pawden,dtset%iomode,hdr_sigma,&
       cryst,ngfftf,cplex1,nfftf,dtset%nspden,qp_rhor_paw,mpi_enreg_seq,ebands=qp_bst)
     end if
     ABI_FREE(qp_rhor_paw)
   end if

   nkxc=0
   if (Dtset%nspden==1) nkxc=2
   if (Dtset%nspden>=2) nkxc=3 !check GGA and spinor that is messy !!!
   !In case of MGGA, fxc and kxc are not available and we dont need them
   !for the screening part (for now ...)
   if (Dtset%ixc<0.and.libxc_functionals_ismgga()) nkxc=0
   if (nkxc/=0)  then
     ABI_MALLOC(qp_kxc,(nfftf,nkxc))
   end if

   ! **** NOTE THAT Vxc CONTAINS THE CORE-DENSITY CONTRIBUTION ****
   n3xccc=0; if (Psps%n1xccc/=0) n3xccc=nfftf
   ABI_MALLOC(qp_vhartr,(nfftf))
   ABI_MALLOC(qp_vtrial,(nfftf,Dtset%nspden))
   ABI_MALLOC(qp_vxc,(nfftf,Dtset%nspden))

   optene=4; moved_atm_inside=0; moved_rhor=0; istep=1

   call setvtr(Cryst%atindx1,Dtset,QP_energies,gmet,gprimd,grchempottn,grewtn,grvdw,gsqcutf_eff,&
&   istep,qp_kxc,mgfftf,moved_atm_inside,moved_rhor,MPI_enreg_seq,&
&   Cryst%nattyp,nfftf,ngfftf,ngrvdw,qp_nhat,qp_nhatgr,nhatgrdim,nkxc,Cryst%ntypat,Psps%n1xccc,n3xccc,&
&   optene,pawrad,Pawtab,ph1df,Psps,qp_rhog,qp_rhor,Cryst%rmet,Cryst%rprimd,strsxc,&
&   Cryst%ucvol,usexcnhat,qp_vhartr,vpsp,qp_vtrial,qp_vxc,vxcavg_qp,Wvl,&
&   xccc3d,Cryst%xred,taug=qp_taug,taur=qp_taur)

   if (allocated(qp_kxc))  then
     ABI_FREE(qp_kxc)
   end if

   if (Dtset%usepaw==1) then
     call timab(561,1,tsec)

     ! === Compute QP Dij ===
     ipert=0; idir=0
     call pawdij(cplex,Dtset%enunit,Cryst%gprimd,ipert,&
&     Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
&     Dtset%nspden,Cryst%ntypat,QP_paw_an,QP_paw_ij,Pawang,Pawfgrtab,&
&     Dtset%pawprtvol,Pawrad,QP_pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&     k0,Dtset%spnorbscl,Cryst%ucvol,dtset%charge,qp_vtrial,qp_vxc,Cryst%xred,&
&     nucdipmom=Dtset%nucdipmom)

     ! === Symmetrize total Dij ===
     option_dij=0

#if 0
     call symdij(Cryst%gprimd,Cryst%indsym,ipert,&
&     Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,option_dij,&
&     QP_paw_ij,Pawang,Dtset%pawprtvol,Pawtab,Cryst%rprimd,&
&     Cryst%symafm,Cryst%symrec)
#else
     call symdij_all(Cryst%gprimd,Cryst%indsym,ipert,&
&     Cryst%natom,Cryst%natom,Cryst%nsym,Cryst%ntypat,&
&     QP_paw_ij,Pawang,Dtset%pawprtvol,Pawtab,Cryst%rprimd,&
&     Cryst%symafm,Cryst%symrec)
#endif

     ! Output the QP pseudopotential strengths Dij and the augmentation occupancies Rhoij.
     call pawprt(Dtset,Cryst%natom,QP_paw_ij,QP_Pawrhoij,Pawtab)
     call timab(561,2,tsec)
   end if

   ehartree=half*SUM(qp_rhor(:,1)*qp_vhartr(:))/DBLE(nfftf)*Cryst%ucvol

   write(msg,'(a,80a)')ch10,('-',ii=1,80)
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(5a,f9.4,3a,es21.14,2a,es21.14)')ch10,&
&   ' QP results after the unitary transformation in the KS subspace: ',ch10,ch10,&
&   '  Number of electrons    = ',qp_rhog(1,1)*Cryst%ucvol,ch10,ch10,&
&   '  QP Band energy    [Ha] = ',get_bandenergy(QP_BSt),ch10,&
&   '  QP Hartree energy [Ha] = ',ehartree
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,80a)')ch10,('-',ii=1,80)
   call wrtout(ab_out,msg,'COLL')
!
!  TODO Since plasmonpole model 2-3-4 depend on the Fourier components of the density
!  in case of self-consistency we might calculate here the ppm coefficients using qp_rhor
 end if ! gwcalctyp>=10
!
!=== KS hamiltonian hlda(b1,b1,k,s)= <b1,k,s|H_s|b1,k,s> ===
 ABI_MALLOC(hlda,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
 hlda = czero

 if (Dtset%nspinor==1) then
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       do ib=b1gw,b2gw
         hlda(ib,ib,ik,spin) = KS_BSt%eig(ib,ik,spin)
       end do
     end do
   end do
 else
   ! Spinorial case
   !  * Note that here vxc contains the contribution of the core.
   !  * Scale ovlp if orthonormalization is not satisfied as npwwfn might be < npwvec.
   !  TODO add spin-orbit case and gwpara 2
   ABI_MALLOC(my_band_list, (wfd%mband))
   ABI_MALLOC(bmask, (wfd%mband))
   bmask = .False.; bmask(b1gw:b2gw) = .True.

   if (Wfd%usepaw==1) then
     ABI_DT_MALLOC(Cp1,(Wfd%natom,Wfd%nspinor))
     call pawcprj_alloc(Cp1,0,Wfd%nlmn_atm)
   end if

   do spin=1,Sigp%nsppol
     do ik_ibz=1,Kmesh%nibz
       ! Distribute bands in [b1gw, b2gw] range
       call wfd_distribute_bands(wfd, ik_ibz, spin, my_nband, my_band_list, bmask=bmask)
       if (my_nband == 0) cycle
       npw_k = Wfd%npwarr(ik_ibz)
       do ii=1,my_nband  ! ib=b1gw,b2gw in sequential
         ib = my_band_list(ii)
         ug1  => Wfd%Wave(ib,ik_ibz,spin)%ug
         cdummy = xdotc(npw_k*Wfd%nspinor,ug1,1,ug1,1)
         ovlp(1) = REAL(cdummy); ovlp(2) = AIMAG(cdummy)

         if (Psps%usepaw==1) then
           call wfd_get_cprj(Wfd,ib,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
           ovlp = ovlp + paw_overlap(Cp1,Cp1,Cryst%typat,Pawtab)
         end if
!        write(std_out,*)ovlp(1),ovlp(2)
         norm=DBLE(ovlp(1)+ovlp(2))
         ovlp(1)=DBLE(ovlp(1)/norm); ovlp(2)=DBLE(ovlp(2)/norm) ! ovlp(2)=cone-ovlp(1)
         hlda(ib,ib,ik_ibz,1) = KS_BSt%eig(ib,ik_ibz,1)*ovlp(1)-KS_me%vxc(ib,ib,ik_ibz,3)
         hlda(ib,ib,ik_ibz,2) = KS_BSt%eig(ib,ik_ibz,1)*ovlp(2)-KS_me%vxc(ib,ib,ik_ibz,4)
         hlda(ib,ib,ik_ibz,3) = KS_me%vxc(ib,ib,ik_ibz,3)
         hlda(ib,ib,ik_ibz,4) = KS_me%vxc(ib,ib,ik_ibz,4)
       end do
     end do
   end do

   call xmpi_sum(hlda, wfd%comm, ierr)

   ABI_FREE(my_band_list)
   ABI_FREE(bmask)
   if (Wfd%usepaw==1) then
     call pawcprj_free(Cp1)
     ABI_DT_FREE(Cp1)
   end if
 end if

!=== Initialize Sigma results ===
!TODO it is better if we use ragged arrays indexed by the k-point
 call sigma_init(Sigp,Kmesh%nibz,Dtset%usepawu,Sr)
 !
 !=== Setup of the bare Hamiltonian := T + v_{loc} + v_{nl} + v_H ===
 !* The representation depends wheter we are updating the wfs or not.
 !* ks_vUme is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
 !* Note that vH matrix elements are calculated using the true uncutted interaction.

 if (mod100<10) then  ! * For one-shot GW use the KS representation.
   Sr%hhartree=hlda-KS_me%vxcval
   !=== Additional goodies for PAW ===
   !  * LDA +U Hamiltonian
   !  * LEXX.
   !  * Core contribution estimated using Fock exchange.
   if (Dtset%usepaw==1) then
     if (Sigp%use_sigxcore==1) Sr%hhartree=hlda - (KS_me%vxc - KS_me%sxcore)
     if (Dtset%usepawu>0) Sr%hhartree=Sr%hhartree-KS_me%vu
     if (Dtset%useexexch>0) then
       MSG_ERROR("useexexch > 0 not implemented")
       Sr%hhartree = Sr%hhartree - KS_me%vlexx
     end if
   end if
 else
   ! Self-consistent on energies and|or wavefunctions.
   !   * For NC get the bare Hamiltonian  $H_{bare}= T+v_{loc}+ v_{nl}$ in the KS representation
   !   * For PAW, calculate the matrix elements of h0, store also the new Dij in QP_Paw_ij.
   !   * h0 is defined as T+vH[tn+nhat+tnZc] + vxc[tnc] + dij_eff and
   !     dij_eff = dij^0 + dij^hartree + dij^xc-dij^xc_val + dijhat - dijhat_val.
   !     In the above expression tn, tnhat are QP quantities.
   if (Dtset%usepaw==0) then
     ABI_MALLOC(hbare,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
     hbare=hlda-KS_me%vhartree-KS_me%vxcval
!
!    * Change basis from KS to QP, hbare is overwritten: A_{QP} = U^\dagger A_{KS} U
     ABI_MALLOC(htmp,(b1gw:b2gw,b1gw:b2gw,Kmesh%nibz,Sigp%nsppol*Sigp%nsig_ab))
     ABI_MALLOC(ctmp,(b1gw:b2gw,b1gw:b2gw))
     ABI_MALLOC(uks2qp,(b1gw:b2gw,b1gw:b2gw))
     htmp=hbare; hbare=czero

     do spin=1,Sigp%nsppol
       do ik=1,Kmesh%nibz
         uks2qp(:,:) = Sr%m_lda_to_qp(b1gw:b2gw,b1gw:b2gw,ik,spin)
         do iab=1,Sigp%nsig_ab
           is_idx=spin; if (Sigp%nsig_ab>1) is_idx=iab
           ctmp(:,:)=MATMUL(htmp(:,:,ik,is_idx),uks2qp(:,:))
           hbare(:,:,ik,is_idx)=MATMUL(TRANSPOSE(CONJG(uks2qp)),ctmp)
         end do
       end do
     end do
     ABI_FREE(htmp)
     ABI_FREE(ctmp)
     ABI_FREE(uks2qp)
   end if ! usepaw==0

   !=== Calculate the QP matrix elements===
   !  * This part is parallelized within MPI_COMM_WORD since each node has all GW wavefunctions.
   !  * For PAW, we have to construct the new bare Hamiltonian.
   call wrtout(std_out,ch10//' *************** QP Energies *******************','COLL')

   call melflags_reset(QP_mflags)
!  if (mod100<20) QP_mflags%only_diago=1 ! For e-only, no need of off-diagonal elements.
   QP_mflags%has_vhartree=1
   if (Dtset%usepaw==1)    QP_mflags%has_hbare  =1
!  QP_mflags%has_vxc     =1
!  QP_mflags%has_vxcval  =1
   if (Sigp%gwcalctyp >100) QP_mflags%has_vxcval_hybrid=1
!  if (Sigp%use_sigxcore==1) QP_mflags%has_sxcore =1
!  if (Dtset%usepawu>0)    QP_mflags%has_vu     =1
!  if (Dtset%useexexch>0)  QP_mflags%has_lexexch=1

   ABI_MALLOC(tmp_kstab,(2,Wfd%nkibz,Wfd%nsppol))
   tmp_kstab=0
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw ! No spin dependent!
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       tmp_kstab(1,ik_ibz,spin)=Sigp%minbnd(ikcalc,spin)
       tmp_kstab(2,ik_ibz,spin)=Sigp%maxbnd(ikcalc,spin)
     end do
   end do

   call calc_vhxc_me(Wfd,QP_mflags,QP_me,Cryst,Dtset,gsqcutf_eff,nfftf,ngfftf,&
&   qp_vtrial,qp_vhartr,qp_vxc,Psps,Pawtab,QP_paw_an,Pawang,Pawfgrtab,QP_paw_ij,dijexc_core,&
&   qp_rhor,qp_rhog,usexcnhat,qp_nhat,qp_nhatgr,nhatgrdim,tmp_kstab,taug=qp_taug,taur=qp_taur)
   ABI_FREE(tmp_kstab)

!  #ifdef DEV_HAVE_SCGW_SYM
   if (mod100>=20 .and. Sigp%symsigma>0) then
     bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
     ABI_MALLOC(qp_irreptab,(bmin:bmax,Kmesh%nibz,Sigp%nsppol))
     qp_irreptab=0
     do spin=1,Sigp%nsppol
       do ikcalc=1,Sigp%nkptgw
         ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
         first_band = Sigp%minbnd(ikcalc,spin)
         last_band  = Sigp%maxbnd(ikcalc,spin)
         if (.not.esymm_failed(QP_sym(ik_ibz,spin))) then
           qp_irreptab(first_band:last_band,ik_ibz,spin) = QP_sym(ik_ibz,spin)%b2irrep(first_band:last_band)
!          qp_irreptab(bmin:bmax,ik_ibz,spin) = QP_sym(ik_ibz,spin)%b2irrep(bmin:bmax)
         end if
       end do
     end do
     call melements_zero(QP_me,qp_irreptab)
     ABI_FREE(qp_irreptab)
   end if
!  #endif

   call melements_print(QP_me,header="Matrix elements in the QP basis set",prtvol=Dtset%prtvol)

   ! * Output the QP pseudopotential strengths Dij and the augmentation occupancies Rhoij.
   if (Dtset%usepaw==1) then
     call wrtout(std_out," *** After calc_vHxc_braket *** ",'COLL')
!    TODO terminate the implementation of this routine.
     call paw_ij_print(QP_Paw_ij,unit=std_out,pawprtvol=Dtset%pawprtvol,pawspnorb=Dtset%pawspnorb,mode_paral="COLL")
     call pawprt(Dtset,Cryst%natom,QP_paw_ij,QP_Pawrhoij,Pawtab)
   end if

   if (Dtset%usepaw==0) then
!    GA : We have an odd bug here. I have to unroll this loop, otherwise it
!    might cause segfault when running on several nodes.
!
!    Sr%hhartree = hbare + QP_me%vhartree
     if(QP_mflags%has_vxcval_hybrid==0) then
       do spin=1,Sigp%nsppol*Sr%nsig_ab
         do ikcalc=1,Sr%nkibz
           do ib1=b1gw,b2gw
             do ib2=b1gw,b2gw
               Sr%hhartree(ib2,ib1,ikcalc,spin) = hbare(ib2,ib1,ikcalc,spin) + QP_me%vhartree(ib2,ib1,ikcalc,spin)
             end do
           end do
         end do
       end do
     else
       do spin=1,Sigp%nsppol*Sr%nsig_ab
         do ikcalc=1,Sr%nkibz
           do ib1=b1gw,b2gw
             do ib2=b1gw,b2gw
               Sr%hhartree(ib2,ib1,ikcalc,spin) = hbare(ib2,ib1,ikcalc,spin) + &
&               QP_me%vhartree(ib2,ib1,ikcalc,spin) + &
&               QP_me%vxcval_hybrid(ib2,ib1,ikcalc,spin)
             end do
           end do
         end do
       end do
     end if
   else
     Sr%hhartree=QP_me%hbare
   end if

!  #ifdef DEV_HAVE_SCGW_SYM
   if (mod100>=20 .and. Sigp%symsigma > 0) then
!    bmin=Sigp%minbdgw; bmax=Sigp%maxbdgw
     do spin=1,Sigp%nsppol
       do ik_ibz=1,Kmesh%nibz
         if (.not.esymm_failed(QP_sym(ik_ibz,spin))) then
           bmin=Sigp%minbnd(ik_ibz,spin); bmax=Sigp%minbnd(ik_ibz,spin)
           do ib2=bmin,bmax
             irr_idx2 = QP_sym(ik_ibz,spin)%b2irrep(ib2)
             do ib1=bmin,bmax
               irr_idx1 = QP_sym(ik_ibz,spin)%b2irrep(ib1)
               if (irr_idx1/=irr_idx2 .and. ALL((/irr_idx1,irr_idx2/)/=0) ) Sr%hhartree(ib1,ib2,ik_ibz,spin) = czero
             end do
           end do
         end if
       end do
     end do
   end if
!  #endif

   ABI_FREE(qp_rhog)
   ABI_FREE(qp_vhartr)
   ABI_FREE(qp_vxc)
   ABI_FREE(qp_taug)
   ABI_FREE(qp_nhat)
   ABI_FREE(qp_nhatgr)
   call melements_free(QP_me)
 end if ! gwcalctyp<10

 ! Free some memory
 if (allocated(hbare)) then
   ABI_FREE(hbare)
 end if
 if (allocated(hlda )) then
   ABI_FREE(hlda)
 end if

 !=== Prepare the storage of QP amplitudes and energies ===
 !* Initialize with KS wavefunctions and energies.
 Sr%eigvec_qp=czero;  Sr%en_qp_diago=zero
 do ib=1,Sigp%nbnds
   Sr%en_qp_diago(ib,:,:)=KS_BSt%eig(ib,:,:)
   Sr%eigvec_qp(ib,ib,:,:)=cone
 end do
 !
 !=== Store <n,k,s|V_xc[n_val]|n,k,s> and <n,k,s|V_U|n,k,s> ===
 !  * Note that we store the matrix elements of V_xc in the KS basis set, not in the QP basis set
 !  * Matrix elements of V_U are zero unless we are using LDA+U as starting point
 do ib=b1gw,b2gw
   Sr%vxcme(ib,:,:)=KS_me%vxcval(ib,ib,:,:)
   if (Dtset%usepawu>0) Sr%vUme (ib,:,:)=KS_me%vu(ib,ib,:,:)
 end do

 !=== Initial guess for the GW energies ===
 !  * Save the energies of the previous iteration.
 do spin=1,Sigp%nsppol
   do ik=1,Kmesh%nibz
     do ib=1,Sigp%nbnds
       Sr%e0 (ib,ik,spin)=QP_BSt%eig(ib,ik,spin)
       Sr%egw(ib,ik,spin)=QP_BSt%eig(ib,ik,spin)
     end do
     Sr%e0gap(ik,spin)=zero
     ks_iv=ks_vbik(ik,spin)
     if (Sigp%nbnds>=ks_iv+1) Sr%e0gap(ik,spin)=Sr%e0(ks_iv+1,ik,spin)-Sr%e0(ks_iv,ik,spin)
   end do
 end do
!
!=== If required apply a scissor operator or update the energies ===
!TODO check if other Sr entries have to be updated
!moreover this part should be done only in case of semiconductors
!FIXME To me it makes more sense if we apply the scissor to KS_BS but I have to RECHECK csigme
 if (ABS(Sigp%mbpt_sciss)>tol6) then
   write(msg,'(6a,f10.5,2a)')ch10,&
&   ' sigma : performing a first self-consistency',ch10,&
&   '  update of the energies in G by a scissor operator',ch10, &
&   '  applying a scissor operator of ',Sigp%mbpt_sciss*Ha_eV,' [eV] ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   do spin=1,Sigp%nsppol
     do ik=1,Kmesh%nibz
       ks_iv=ks_vbik(ik,spin)
       if (Sigp%nbnds>=ks_iv+1) then
         Sr%egw    (ks_iv+1:Sigp%nbnds,ik,spin) = Sr%egw    (ks_iv+1:Sigp%nbnds,ik,spin)+Sigp%mbpt_sciss
         QP_BSt%eig(ks_iv+1:Sigp%nbnds,ik,spin) = QP_BSt%eig(ks_iv+1:Sigp%nbnds,ik,spin)+Sigp%mbpt_sciss
       end if
     end do
   end do
!  %call apply_scissor(QP_BSt,Sigp%mbpt_sciss)
 else if (.FALSE.) then
   write(msg,'(4a)')ch10,&
&   ' sigma : performing a first self-consistency',ch10,&
&   '  update of the energies in G by a previous GW calculation'
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
!  TODO Recheck this part, is not clear to me!
   ABI_MALLOC(igwene,(QP_Bst%mband,QP_Bst%nkpt,QP_Bst%nsppol))
   call rdgw(QP_BSt,'__in.gw__',igwene,extrapolate=.TRUE.)
   ABI_FREE(igwene)
   Sr%egw=QP_BSt%eig
   !
   ! * Recalculate the new fermi level.
   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=0)
 end if
!
!In case of AC refer all the energies wrt to the fermi level
!Take care because results from ppmodel cannot be used for AC
!FIXME check ks_energy or qp_energy (in case of SCGW?)

 if (mod10==SIG_GW_AC) then
   ! All these quantities will be passed to csigme
   ! if I skipped the self-consistent part then here I have to use fermi
   QP_BSt%eig = QP_BSt%eig -QP_BSt%fermie
   Sr%egw = Sr%egw-QP_BSt%fermie
   Sr%e0  = Sr%e0 -QP_BSt%fermie
   oldefermi=QP_BSt%fermie
   ! TODO Recheck fermi
   ! Clean EVERYTHING in particulare the treatment of E fermi
   QP_BSt%fermie=zero
 end if
 !
 ! === Setup frequencies around the KS\QP eigenvalues to compute Sigma derivatives (notice the spin) ===
 ! TODO it is better using an odd Sr%nomega4sd so that the KS\QP eigenvalue is in the middle
 ioe0j=Sr%nomega4sd/2+1
 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw
     ib1=Sigp%minbnd(ikcalc,spin)
     ib2=Sigp%maxbnd(ikcalc,spin)
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc))
     do jb=ib1,ib2
       do io=1,Sr%nomega4sd
         Sr%omega4sd(jb,ik_ibz,io,spin)=Sr%egw(jb,ik_ibz,spin)+Sigp%deltae*(io-ioe0j)
       end do
     end do
   end do
 end do

 call timab(408,2,tsec) ! hqp_init
 call timab(409,1,tsec) ! getW

 ! === Get epsilon^{-1} either from the _SCR or the _SUSC file and store it in Er%epsm1 ===
 !   * If Er%mqmem==0, allocate and read a single q-slice inside csigme.
 !     TODO Er%nomega should be initialized so that only the frequencies really needed are stored in memory

 ! TODO: The same piece of code is present in screening.
 if (sigma_needs_w(Sigp)) then
   select case (dtset%gwgamma)
   case (0)
     id_required=4; ikxc=0; approx_type=0; option_test=0; dim_kxcg=0
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))

   case (1, 2)
     ! ALDA TDDFT kernel vertex
     !ABI_CHECK(Dtset%usepaw==0,"GWGamma=1 or 2 + PAW not available")
     ABI_CHECK(Er%ID==0,"Er%ID should be 0")

     if (Dtset%usepaw==1) then
       ! If we have PAW, we need the full density on the fine grid
       ABI_MALLOC(ks_aepaw_rhor,(nfftf,Wfd%nspden))
       if (Dtset%getpawden==0 .and. Dtset%irdpawden==0) then
         MSG_ERROR("Must use get/irdpawden to provide a _PAWDEN file!")
       end if
       call wrtout(std_out,sjoin('Checking for existence of: ',Dtfil%filpawdensin),"COLL")
       if (.not. file_exists(dtfil%filpawdensin)) then
         MSG_ERROR(sjoin("Missing file:", dtfil%filpawdensin))
       end if

       ABI_DT_MALLOC(tmp_pawrhoij,(cryst%natom*wfd%usepaw))
       call read_rhor(Dtfil%filpawdensin, cplex1, nfftf_tot, Wfd%nspden, ngfftf, 1, MPI_enreg_seq, &
       ks_aepaw_rhor, hdr_rhor, tmp_pawrhoij, wfd%comm)

       call hdr_free(hdr_rhor)
       call pawrhoij_free(tmp_pawrhoij)
       ABI_DT_FREE(tmp_pawrhoij)
     end if ! Dtset%usepaw==1

     id_required=4; ikxc=7; approx_type=1; dim_kxcg=1
     if (Dtset%gwgamma==1) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==2) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))

     dbg_mode=.FALSE.
     if (Dtset%usepaw==1) then
       ! Use PAW all-electron density
       call kxc_driver(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,ks_aepaw_rhor,&
       Er%npwe,dim_kxcg,kxcg,Gsph_c%gvec,xmpi_comm_self,dbg_mode=dbg_mode)
     else
       ! Norm-conserving
       call kxc_driver(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,ks_rhor,&
       Er%npwe,dim_kxcg,kxcg,Gsph_c%gvec,xmpi_comm_self,dbg_mode=dbg_mode)
     end if

   case (3, 4)
     ! ADA non-local kernel vertex
     !ABI_CHECK(Dtset%usepaw==0,"GWGamma + PAW not available")
     ABI_CHECK(Er%ID==0,"Er%ID should be 0")
     ABI_CHECK(Sigp%nsppol==1,"ADA vertex for GWGamma not available yet for spin-polarised cases")

     if (Dtset%usepaw==1) then ! If we have PAW, we need the full density on the fine grid
       ABI_MALLOC(ks_aepaw_rhor,(nfftf,Sigp%nsppol))
       if (Dtset%getpawden==0.and.Dtset%irdpawden==0) then
         MSG_ERROR("Must use get/irdpawden to provide a _PAWDEN file!")
       end if
       call wrtout(std_out,sjoin('Checking for existence of: ',Dtfil%filpawdensin),"COLL")
       if (.not. file_exists(dtfil%filpawdensin)) then
         MSG_ERROR(sjoin("Missing file:", dtfil%filpawdensin))
       end if

       ABI_DT_MALLOC(tmp_pawrhoij,(cryst%natom*wfd%usepaw))

       call read_rhor(Dtfil%filpawdensin, cplex1, nfftf_tot, Wfd%nspden, ngfftf, 1, MPI_enreg_seq, &
       ks_aepaw_rhor, hdr_rhor, tmp_pawrhoij, wfd%comm)

       call hdr_free(hdr_rhor)
       call pawrhoij_free(tmp_pawrhoij)
       ABI_DT_FREE(tmp_pawrhoij)
     end if ! Dtset%usepaw==1

     id_required=4; ikxc=7; approx_type=2;
     if (Dtset%gwgamma==3) option_test=1 ! TESTELECTRON, vertex in chi0 *and* sigma
     if (Dtset%gwgamma==4) option_test=0 ! TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(fxc_ADA,(Er%npwe,Er%npwe,Er%nqibz))
!    Use userrd to set kappa
     if (Dtset%userrd==zero) Dtset%userrd = 2.1_dp
!    Set correct value of kappa (should be scaled with alpha*r_s where)
!    r_s is Wigner-Seitz radius and alpha=(4/(9*Pi))^(1/3)
     rhoav = (drude_plsmf*drude_plsmf)/four_pi
     r_s = (three/(four_pi*rhoav))**third
     alpha = (four*ninth*piinv)**third
     Dtset%userrd = Dtset%userrd/(alpha*r_s)

     dbg_mode=.TRUE.
     if (Dtset%usepaw==1) then
       ! Use PAW all-electron density
       call kxc_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,&
       ks_aepaw_rhor,Er%npwe,Er%nqibz,Er%qibz,&
       fxc_ADA,Gsph_c%gvec,xmpi_comm_self,kappa_init=Dtset%userrd,dbg_mode=dbg_mode)
     else
       ! Norm conserving
       call kxc_ADA(Dtset,Cryst,ikxc,ngfftf,nfftf_tot,Wfd%nspden,&
       ks_rhor,Er%npwe,Er%nqibz,Er%qibz,&
       fxc_ADA,Gsph_c%gvec,xmpi_comm_self,kappa_init=Dtset%userrd,dbg_mode=dbg_mode)
     end if

     dim_kxcg = 0
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))

!  @WC: bootstrap --
   case (-3, -4, -5, -6, -7, -8)
!    ABI_CHECK(Dtset%usepaw==0,"GWGamma=1 or 2 + PAW not available")
     ABI_CHECK(Er%ID==0,"Er%ID should be 0")

     if (Dtset%usepaw==1) then
       ! If we have PAW, we need the full density on the fine grid
       ABI_MALLOC(ks_aepaw_rhor,(nfftf,Wfd%nspden))
       if (Dtset%getpawden==0.and.Dtset%irdpawden==0) then
         MSG_ERROR("Must use get/irdpawden to provide a _PAWDEN file!")
       end if
       call wrtout(std_out,sjoin('Checking for existence of: ',Dtfil%filpawdensin),"COLL")
       if (.not. file_exists(dtfil%filpawdensin)) then
         MSG_ERROR(sjoin("Missing file:", dtfil%filpawdensin))
       end if

       ABI_DT_MALLOC(tmp_pawrhoij,(cryst%natom*wfd%usepaw))

       call read_rhor(Dtfil%filpawdensin, cplex1, nfftf_tot, Wfd%nspden, ngfftf, 1, MPI_enreg_seq, &
       ks_aepaw_rhor, hdr_rhor, tmp_pawrhoij, wfd%comm)

       call hdr_free(hdr_rhor)
       call pawrhoij_free(tmp_pawrhoij)
       ABI_DT_FREE(tmp_pawrhoij)
     end if ! Dtset%usepaw==1

     id_required=4; ikxc=7; dim_kxcg=0

     if (dtset%gwgamma>-5) then
       approx_type=4  ! full fxc(G,G')
     else if (dtset%gwgamma>-7) then
       approx_type=5  ! fxc(0,0) one-shot
     else
       approx_type=6  ! rpa-type bootstrap
     end if

     option_test=MOD(Dtset%gwgamma,2)
     ! 1 -> TESTELECTRON, vertex in chi0 *and* sigma
     ! 0 -> TESTPARTICLE, vertex in chi0 only
     ABI_MALLOC(kxcg,(nfftf_tot,dim_kxcg))
     ! --@WC

   case default
     MSG_ERROR(sjoin("Wrong gwgamma:", itoa(dtset%gwgamma)))
   end select

!  Set plasma frequency
   my_plsmf = drude_plsmf; if (Dtset%ppmfrq>tol6) my_plsmf = Dtset%ppmfrq
   Dtset%ppmfrq = my_plsmf

!   if(dtset%ucrpa==0) then
   if (Dtset%gwgamma<3) then
     call mkdump_Er(Er,Vcp,Er%npwe,Gsph_c%gvec,dim_kxcg,kxcg,id_required,&
&     approx_type,ikxc,option_test,Dtfil%fnameabo_scr,Dtset%iomode,&
&     nfftf_tot,ngfftf,comm)
   else
     call mkdump_Er(Er,Vcp,Er%npwe,Gsph_c%gvec,dim_kxcg,kxcg,id_required,&
&     approx_type,ikxc,option_test,Dtfil%fnameabo_scr,Dtset%iomode,&
&     nfftf_tot,ngfftf,comm,fxc_ADA)
   end if
!   end if
   if (allocated(kxcg))  then
     ABI_FREE(kxcg)
   end if
   if (allocated(fxc_ADA))  then
     ABI_FREE(fxc_ADA)
   end if
   if (allocated(ks_aepaw_rhor))  then
     ABI_FREE(ks_aepaw_rhor)
   end if
 end if
 !
 ! ================================================
 ! ==== Calculate plasmonpole model parameters ====
 ! ================================================
 ! TODO Maybe its better if we use mqmem as input variable
 use_aerhor=0
 ABI_MALLOC(ks_aepaw_rhor,(nfftf,Wfd%nspden*use_aerhor))

 if (sigma_needs_ppm(Sigp)) then
   my_plsmf=drude_plsmf; if (Dtset%ppmfrq>tol6) my_plsmf=Dtset%ppmfrq
   call ppm_init(PPm,Er%mqmem,Er%nqibz,Er%npwe,Sigp%ppmodel,my_plsmf,Dtset%gw_invalid_freq)

   ! PPm%force_plsmf= force_ppmfrq  ! this line to change the plasme frequency in HL expression.

   if (Wfd%usepaw==1 .and. Ppm%userho==1) then
     ! * For PAW and ppmodel 2-3-4 we need the AE rho(G) without compensation charge.
     ! * It would be possible to calculate rho(G) using Paw_pwff, though. It should be faster but
     !    results will depend on the expression used for the matrix elements. This approach is safer.
     use_aerhor=1
     ABI_FREE(ks_aepaw_rhor)
     ABI_MALLOC(ks_aepaw_rhor,(nfftf,Wfd%nspden))

     ! Check if input density file is available, otherwise compute
     pawden_fname = strcat(Dtfil%filnam_ds(3), '_PAWDEN')
     call wrtout(std_out,sjoin('Checking for existence of:',pawden_fname))
     if (file_exists(pawden_fname)) then
       ! Read density from file
       ABI_DT_MALLOC(tmp_pawrhoij,(cryst%natom*wfd%usepaw))

       call read_rhor(pawden_fname, cplex1, nfftf_tot, Wfd%nspden, ngfftf, 1, MPI_enreg_seq, &
       ks_aepaw_rhor, hdr_rhor, tmp_pawrhoij, wfd%comm)

       call hdr_free(hdr_rhor)
       call pawrhoij_free(tmp_pawrhoij)
       ABI_DT_FREE(tmp_pawrhoij)
     else
       ! Have to calculate PAW AW rhor from scratch
       ABI_MALLOC(qp_rhor_n_one,(pawfgr%nfft,Dtset%nspden))
       ABI_MALLOC(qp_rhor_nt_one,(pawfgr%nfft,Dtset%nspden))

!      FIXME
       MSG_WARNING(" denfgr in sigma seems to produce wrong results")
       write(std_out,*)" input tilde ks_rhor integrates: ",SUM(ks_rhor(:,1))*Cryst%ucvol/nfftf

       call denfgr(Cryst%atindx1,Cryst%gmet,comm,Cryst%natom,Cryst%natom,Cryst%nattyp,ngfftf,ks_nhat,Wfd%nspinor,&
       Wfd%nsppol,Wfd%nspden,Cryst%ntypat,Pawfgr,Pawrad,KS_Pawrhoij,Pawtab,Dtset%prtvol, &
       ks_rhor,ks_aepaw_rhor,qp_rhor_n_one,qp_rhor_nt_one,Cryst%rprimd,Cryst%typat,Cryst%ucvol,Cryst%xred)

       ABI_FREE(qp_rhor_n_one)
       ABI_FREE(qp_rhor_nt_one)
     end if

     write(msg,'(a,f8.4)')' sigma: PAW AE density used for PPmodel integrates to: ',SUM(ks_aepaw_rhor(:,1))*Cryst%ucvol/nfftf
     call wrtout(std_out,msg,'PERS')

     if (Er%mqmem/=0) then
       ! Calculate ppmodel parameters for all q-points.
       call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfftf,Gsph_c%gvec,ngfftf,ks_aepaw_rhor(:,1))
     end if

   else
     ! NC or PAW with PPmodel 1.
     if (Er%mqmem/=0) then
       ! Calculate ppmodel parameters for all q-points
       call setup_ppmodel(PPm,Cryst,Qmesh,Er%npwe,Er%nomega,Er%omega,Er%epsm1,nfftf,Gsph_c%gvec,ngfftf,ks_rhor(:,1))
     end if
   end if ! PAW or NC PPm and/or needs density

!  If nfreqre is set, then output the dielectric function for these frequencies to file
   write(msg,'(a,i0,a,f8.2,a)')&
&   ' Checking if PPm em1 is output, nfreqre is: ',Dtset%nfreqre,' freqremax is: ',Dtset%freqremax*Ha_eV,' eV'
   MSG_WARNING(msg)

   if (Dtset%nfreqre>1) then
     ABI_CHECK(Dtset%freqremax>tol6,'freqremax must be set to maximum frequency')
     ABI_MALLOC(em1_ppm,(1,1,Dtset%nfreqre))
     ABI_MALLOC(omega,(Dtset%nfreqre))
     do iomega=1,Dtset%nfreqre
       omega(iomega) =  CMPLX((Dtset%freqremax/REAL((Dtset%nfreqre-1)))*(iomega-1),zero)
     end do

     call getem1_from_PPm(PPm,1,1,Er%Hscr%zcut,Dtset%nfreqre,omega,Vcp,em1_ppm,only_ig1=1,only_ig2=1)

     if (open_file(Dtfil%fnameabo_em1_lf,msg,newunit=ppm_unt,form="formatted", status="unknown", action="write") /= 0) then
       MSG_ERROR(msg)
     end if
     write(ppm_unt,'(a)')'#'
     write(ppm_unt,'(a)')'# Macroscopic plasmon-pole Dielectric Function without local fields'
     write(ppm_unt,'(a)')'# Omega [eV]    Re epsilon_M       IM eps_M '
     write(ppm_unt,'(a)')'# ppmodel: ',PPm%model
     do iomega=1,Dtset%nfreqre
       write(ppm_unt,'((3es16.8))')REAL(omega(iomega))*Ha_eV,REAL(em1_ppm(1,1,iomega))
     end do
     close(ppm_unt)

     ABI_FREE(em1_ppm)
     ABI_FREE(omega)
     MSG_WARNING('Wrote epsilon-1 for PPm to file: '//TRIM(Dtfil%fnameabo_em1_lf))
   end if ! Check for output of eps^-1 for PPm

 end if ! sigma_needs_ppm

 call timab(409,2,tsec) ! getW

 if (wfd_iam_master(Wfd)) then ! Write info on the run on ab_out, then open files to store final results.
   call ebands_report_gap(KS_BSt,header='KS Band Gaps',unit=ab_out)
   if(dtset%ucrpa==0) then
     call write_sigma_header(Sigp,Er,Cryst,Kmesh,Qmesh)
   end if

   if (open_file(Dtfil%fnameabo_gw,msg,unit=unt_gw,status='unknown',form='formatted') /=0) then
     MSG_ERROR(msg)
   end if
   write(unt_gw,*)Sigp%nkptgw,Sigp%nsppol

   if (open_file(Dtfil%fnameabo_gwdiag,msg,unit=unt_gwdiag,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   write(unt_gwdiag,*)Sigp%nkptgw,Sigp%nsppol

   if (open_file(Dtfil%fnameabo_sig,msg,unit=unt_sig,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if
   if (open_file(Dtfil%fnameabo_sgr,msg,unit=unt_sgr,status='unknown',form='formatted') /= 0) then
     MSG_ERROR(msg)
   end if

   if (mod10==SIG_GW_AC) then ! Sigma along the imaginary axis.
     if (open_file(Dtfil%fnameabo_sgm,msg,unit=unt_sgm,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
   end if
 end if
!
!=======================================================================
!==== Calculate self-energy and output the results for each k-point ====
!=======================================================================
!* Here it would be possible to calculate the QP correction for the same k-point using a PPmodel
!in the first iteration just to improve the initial guess and CD or AC in the second step. Useful
!if the KS level is far from the expected QP results. Previously it was possible to do such
!calculation by simply listing the same k-point twice in kptgw. Now this trick is not allowed anymore.
!Everything, indeed, should be done in a clean and transparent way inside csigme.

 call wfd_print(Wfd,mode_paral='PERS')

 call wrtout(std_out,sigma_type_from_key(mod10),'COLL')

 if (mod100<10) then
   msg = " Perturbative Calculation"
   if (mod100==1) msg = " Newton Raphson method "
 else if (mod100<20) then
   msg = " Self-Consistent on Energies only"
 else
   msg = " Self-Consistent on Energies and Wavefunctions"
 end if
 if(Dtset%ucrpa==0) then
   call wrtout(std_out,msg,'COLL')
 end if
!
!=================================================
!==== Calculate the matrix elements of Sigma =====
!=================================================

 nomega_sigc=Sr%nomega_r+Sr%nomega4sd; if (mod10==SIG_GW_AC) nomega_sigc=Sr%nomega_i

 ib1=Sigp%minbdgw ! min and max band indeces for GW corrections.
 ib2=Sigp%maxbdgw

 !MG TODO: I don't like the fact that ib1 and ib2 are redefined here because this
 ! prevents me from refactoring the code. In particular I want to store the self-energy
 ! results inside the sigma_results datatypes hence one needs to know all the dimensions
 ! at the beginning of the execution (e.g. in setup_sigma) so that one can easily allocate the arrays in the type.
 if(Dtset%ucrpa>=1) then
   !Read the band
   if (open_file("forlb.ovlp",msg,newunit=temp_unt,form="formatted", status="unknown") /= 0) then
     MSG_ERROR(msg)
   end if
   rewind(temp_unt)
   read(temp_unt,*)
   read(temp_unt,*) msg, ib1, ib2
   close(temp_unt)
 end if

 ABI_MALLOC(sigcme,(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nkptgw,Sigp%nsppol*Sigp%nsig_ab))
 sigcme=czero
!
!==========================================================
!==== Exchange part using the dense gwx_ngfft FFT mesh ====
!==========================================================
!TODO : load distribution is not optimal if band parallelism is used.
!but this problem was also affecting the old implementation.
 call timab(421,1,tsec) ! calc_sigx_me

 ! This routine shows how the wavefunctions are distributed.
 !call wfd_show_bkstab(Wfd,unit=std_out)

 if(Dtset%ucrpa>=1) then
!  !Calculation of the Rho_n for calc_ucrpa
   call wrtout(std_out,sjoin("begin of Ucrpa calc for a nb of kpoints of: ",itoa(Sigp%nkptgw)),"COLL")
!   Wannier basis: rhot1_q_m will need an atom index in the near future.
   ic=0
   do  itypat=1,cryst%ntypat
     if(pawtab(itypat)%lpawu.ne.-1) then
       ndim=2*dtset%lpawu(itypat)+1
       ic=ic+1
       itypatcor=itypat
       lcor=dtset%lpawu(itypat)
     end if
   end do
   if(ic>1) then
     MSG_ERROR("number of correlated species is larger than one")
   end if
   ABI_ALLOCATE(rhot1_q_m,(cryst%nattyp(itypatcor),Wfd%nspinor,Wfd%nspinor,ndim,ndim,sigp%npwx,Qmesh%nibz))
   ABI_ALLOCATE(M1_q_m,(cryst%nattyp(itypatcor),Wfd%nspinor,Wfd%nspinor,ndim,ndim,sigp%npwx,Qmesh%nibz))

   M1_q_m=czero
   rhot1_q_m=czero

!   do ikcalc=1,Sigp%nkptgw

!   if(cryst%nsym==1) nkcalc=Kmesh%nbz
!   if(cryst%nsym>1)  nkcalc=Kmesh%nibz
   nkcalc=Kmesh%nbz
   do ikcalc=1,nkcalc ! for the oscillator strengh, spins are identical without SOC
!    if(Sigp%nkptgw/=Kmesh%nbz) then
!      write(msg,'(6a)')ch10,&
!&      ' nkptgw and nbz differs: this is not allowed to compute U in cRPA '
!      MSG_ERROR(msg)
!    endif
     write(std_out,*) "ikcalc",ikcalc,size(Sigp%kptgw2bz),nkcalc

     if(mod(ikcalc-1,nprocs)==Wfd%my_rank) then
       if(cryst%nsym==1) ik_ibz=Kmesh%tab(ikcalc) ! Index of the irred k-point for GW
       if(cryst%nsym>1) ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW

       call prep_calc_ucrpa(ik_ibz,ikcalc,itypatcor,ib1,ib2,Cryst,QP_bst,Sigp,Gsph_x,Vcp,&
&       Kmesh,Qmesh,lcor,M1_q_m,Pawtab,Pawang,Paw_pwff,&
&       Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,gwx_ngfft,&
&       ngfftf,Dtset%prtvol,Dtset%pawcross,rhot1_q_m)
     end if
   end do

   call xmpi_sum(rhot1_q_m,Wfd%comm,ierr)
   call xmpi_sum(M1_q_m,Wfd%comm,ierr)
!   if(Cryst%nsym==1) then
   M1_q_m=M1_q_m/Kmesh%nbz/Wfd%nsppol
   rhot1_q_m=rhot1_q_m/Kmesh%nbz/Wfd%nsppol
!   endif
!  Calculation of U in cRPA: need to treat with a different cutoff the
!  bare coulomb and the screened coulomb interaction.
   call calc_ucrpa(itypatcor,cryst,Kmesh,lcor,M1_q_m,Qmesh,Er%npwe,sigp%npwx,&
&   Cryst%nsym,rhot1_q_m,Sigp%nomegasr,Sigp%minomega_r,Sigp%maxomega_r,ib1,ib2,'Gsum',Cryst%ucvol,Wfd,Er%fname)

   ABI_DEALLOCATE(rhot1_q_m)
   ABI_DEALLOCATE(M1_q_m)

 else
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
     ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indeces for GW corrections (for this k-point)
     ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

     call calc_sigx_me(ik_ibz,ikcalc,ib1,ib2,Cryst,QP_bst,Sigp,Sr,Gsph_x,Vcp,Kmesh,Qmesh,Ltg_k(ikcalc),&
&     Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&
&     gwx_ngfft,ngfftf,Dtset%prtvol,Dtset%pawcross,Dtset%gwfockmix)
   end do

!  for the time being, do not remove this barrier!
   call xmpi_barrier(Wfd%comm)

   call timab(421,2,tsec) ! calc_sigx_me
   !
   ! ==========================================================
   ! ==== Correlation part using the coarse gwc_ngfft mesh ====
   ! ==========================================================
   if (mod10/=SIG_HF) then
     do ikcalc=1,Sigp%nkptgw
       ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
       ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indices for GW corrections (for this k-point)
       ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

       sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)

       if (any(mod10 == [SIG_SEX, SIG_COHSEX])) then
         ! Calculate static COHSEX or SEX using the coarse gwc_ngfft mesh.

         call cohsex_me(ik_ibz,ikcalc,nomega_sigc,ib1,ib2,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,Kmesh,Qmesh,&
&         Ltg_k(ikcalc),Pawtab,Pawang,Paw_pwff,Psps,Wfd,QP_sym,&
&         gwc_ngfft,Dtset%iomode,Dtset%prtvol,sigcme_p)
       else

         call calc_sigc_me(ik_ibz,ikcalc,nomega_sigc,ib1,ib2,Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,&
&         Ltg_k(ikcalc),PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&
&         gwc_ngfft,ngfftf,nfftf,ks_rhor,use_aerhor,ks_aepaw_rhor,sigcme_p)
       end if

     end do
   end if

   call xmpi_barrier(Wfd%comm)
   !
   !  =====================================================
   !  ==== Solve Dyson equation storing results in Sr% ====
   !  =====================================================
   !  * Use perturbative approach or AC to find QP corrections.
   !  * If qp-GW, diagonalize also H0+Sigma in the KS basis set to get the
   !  new QP amplitudes and energies (Sr%eigvec_qp and Sr%en_qp_diago.
   !  TODO AC with spinor not implemented yet.
   !  TODO Diagonalization of Sigma+hhartre with AC is wrong.
   !
   call timab(425,1,tsec) ! solve_dyson
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
     ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indeces for GW corrections (for this k-point)
     ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

     sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)

     call solve_dyson(ikcalc,ib1,ib2,nomega_sigc,Sigp,Kmesh,sigcme_p,QP_BSt%eig,Sr,Dtset%prtvol,Dtfil,Wfd%comm)
     !
     ! === Calculate direct gap for each spin and print out final results ===
     !   * We use the valence index of the KS system because we still do not know
     !     all the QP corrections. Ideally one should use the QP valence index
     do spin=1,Sigp%nsppol
       if ( Sigp%maxbnd(ikcalc,spin) >= ks_vbik(ik_ibz,spin)+1 .and. &
&       Sigp%minbnd(ikcalc,spin) <= ks_vbik(ik_ibz,spin)  ) then
         ks_iv=ks_vbik(ik_ibz,spin)
         Sr%egwgap (ik_ibz,spin)=  Sr%egw(ks_iv+1,ik_ibz,spin) -  Sr%egw(ks_iv,ik_ibz,spin)
         Sr%degwgap(ik_ibz,spin)= Sr%degw(ks_iv+1,ik_ibz,spin) - Sr%degw(ks_iv,ik_ibz,spin)
       else ! The "gap" cannot be computed
         Sr%e0gap  (ik_ibz,spin)=zero
         Sr%egwgap (ik_ibz,spin)=zero
         Sr%degwgap(ik_ibz,spin)=zero
       end if
     end do

     if (wfd_iam_master(Wfd)) then
       call write_sigma_results(ikcalc,ik_ibz,Sigp,Sr,KS_BSt)
     end if
   end do !ikcalc

   call timab(425,2,tsec) ! solve_dyson
   call timab(426,1,tsec) ! finalize
   !
   ! === Update the energies in QP_BSt ===
   !  * If QPSCGW, use diagonalized eigenvalues otherwise perturbative results.
   if (mod100>=10) then
     do ib=1,Sigp%nbnds
       QP_BSt%eig(ib,:,:)=Sr%en_qp_diago(ib,:,:)
     end do
   else
     QP_BSt%eig=Sr%egw
   end if

   ! ================================================================================
   ! ==== This part is done only if all k-points in the IBZ have been calculated ====
   ! ================================================================================
   if (Sigp%nkptgw==Kmesh%nibz) then

     ! Recalculate new occupations and Fermi level.
     call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=Dtset%prtvol)
     qp_vbik(:,:) = get_valence_idx(QP_BSt)

     write(msg,'(2a,3x,2(es16.6,a))')ch10,' New Fermi energy : ',QP_BSt%fermie,' Ha ,',QP_BSt%fermie*Ha_eV,' eV'
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')

     ! === If all k-points and all occupied bands are calculated, output EXX ===
     ! FIXME here be careful about the check on  ks_vbik in case of metals
     ! if (my_rank==master.and.Sigp%nkptgw==Kmesh%nibz.and.ALL(Sigp%minbnd(:)==1).and.ALL(Sigp%maxbnd(:)>=MAXVAL(nbv(:)))) then
     ! if (ALL(Sigp%minbnd==1).and. ALL(Sigp%maxbnd>=MAXVAL(MAXVAL(ks_vbik(:,:),DIM=1))) ) then
     if (ALL(Sigp%minbnd==1).and. ALL(Sigp%maxbnd>=ks_vbik) ) then  ! FIXME here the two arrays use a different indexing.

       ex_energy = sigma_get_exene(Sr,Kmesh,QP_BSt)
       write(msg,'(a,2(es16.6,a))')' New Exchange energy : ',ex_energy,' Ha ,',ex_energy*Ha_eV,' eV'
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
     end if
     !
     ! * Report the QP gaps (Fundamental and Optical)
     call ebands_report_gap(QP_BSt,header='QP Band Gaps',unit=ab_out)

     ! Band structure interpolation from QP energies computed on the k-mesh.
     if (nint(dtset%einterp(1)) /= 0 .and. all(sigp%minbdgw == sigp%minbnd) .and. all(sigp%maxbdgw == sigp%maxbnd)) then
       call ebands_interpolate_kpath(QP_BSt, dtset, cryst, [sigp%minbdgw, sigp%maxbdgw], dtfil%filnam_ds(4), comm)
     end if
   end if ! Sigp%nkptgw==Kmesh%nibz
   !
   ! === Write SCF data in case of self-consistent calculation ===
   !  * Save Sr%en_qp_diago, Sr%eigvec_qp and m_lda_to_qp in the _QPS file.
   !  * Note that in the first iteration qp_rhor contains KS rhor, then the mixed rhor.
   if (mod100>=10) then

     ! Calculate the new m_lda_to_qp
     call updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,Sr%m_lda_to_qp)

     if (wfd_iam_master(Wfd)) then
       call wrqps(Dtfil%fnameabo_qps,Sigp,Cryst,Kmesh,Psps,Pawtab,QP_Pawrhoij,&
&       Dtset%nspden,nscf,nfftf,ngfftf,Sr,QP_BSt,Sr%m_lda_to_qp,qp_rhor)
     end if

     ! === Report the MAX variation for each kptgw and spin ===
     call wrtout(ab_out,ch10//' Convergence of QP corrections ','COLL')
     converged=.TRUE.
     do spin=1,Sigp%nsppol
       write(msg,'(a,i2,a)')' >>>>> For spin ',spin,' <<<<< '
       call wrtout(ab_out,msg,'COLL')
       do ikcalc=1,Sigp%nkptgw
         ib1    = Sigp%minbnd(ikcalc,spin)
         ib2    = Sigp%maxbnd(ikcalc,spin)
         ik_bz  = Sigp%kptgw2bz(ikcalc)
         ik_ibz = Kmesh%tab(ik_bz)
         ii     = imax_loc( ABS(Sr%degw(ib1:ib2,ik_ibz,spin)) )
         max_degw = Sr%degw(ii,ik_ibz,spin)
         write(msg,('(a,i3,a,2f8.3,a,i3)'))'.  kptgw no:',ikcalc,'; Maximum DeltaE = (',max_degw*Ha_eV,') for band index:',ii
         call wrtout(ab_out,msg,'COLL')
         converged = (converged .and. ABS(max_degw) < Dtset%gw_toldfeig)
       end do
     end do
   end if
   ! ==============================================
   ! ==== Save the GW results in a NETCDF file ====
   ! ==============================================
   if (wfd_iam_master(Wfd)) then
     fname = TRIM(Dtfil%filnam_ds(4))//'_SIGRES.nc'
#ifdef HAVE_NETCDF
     NCF_CHECK(nctk_open_create(ncid, fname, xmpi_comm_self))
     NCF_CHECK(nctk_defnwrite_ivars(ncid, ["sigres_version"], [1]))
     NCF_CHECK(crystal_ncwrite(Cryst, ncid))
     NCF_CHECK(ebands_ncwrite(KS_Bst, ncid))
     NCF_CHECK(sigma_ncwrite(Sigp, Er, Sr, ncid))
     ! Add qp_rhor. Note that qp_rhor == ks_rhor if wavefunctions are not updated.
     !ncerr = nctk_write_datar("qp_rhor",path,ngfft,cplex,nfft,nspden,&
     !                          comm_fft,fftn3_distrib,ffti3_local,datar,action)
     NCF_CHECK(nf90_close(ncid))
#endif
   end if

 end if ! ucrpa

 !----------------------------- END OF THE CALCULATION ------------------------
 !
 !=====================
 !==== Close Files ====
 !=====================
 if (wfd_iam_master(Wfd)) then
   close(unt_gw )
   close(unt_gwdiag)
   close(unt_sig)
   close(unt_sgr)
   if (mod10==SIG_GW_AC) close(unt_sgm)
 end if
 !
 !===============================================
 !==== Free arrays and local data structures ====
 !===============================================
 ABI_FREE(ks_vbik)
 ABI_FREE(qp_vbik)
 ABI_FREE(ph1d)
 ABI_FREE(ph1df)
 ABI_FREE(qp_rhor)
 ABI_FREE(ks_rhor)
 ABI_FREE(ks_rhog)
 ABI_FREE(qp_taur)
 ABI_FREE(ks_taur)
 ABI_FREE(ks_taug)
 ABI_FREE(ks_vhartr)
 ABI_FREE(ks_vtrial)
 ABI_FREE(vpsp)
 ABI_FREE(ks_vxc)
 ABI_FREE(xccc3d)
 ABI_FREE(grchempottn)
 ABI_FREE(grewtn)
 ABI_FREE(grvdw)
 ABI_FREE(sigcme)

 if (allocated(kxc)) then
   ABI_FREE(kxc)
 end if
 if (allocated(qp_vtrial)) then
   ABI_FREE(qp_vtrial)
 end if

 ABI_FREE(ks_nhat)
 ABI_FREE(ks_nhatgr)
 ABI_FREE(dijexc_core)
 ABI_FREE(ks_aepaw_rhor)
 call pawfgr_destroy(Pawfgr)

 ! Deallocation for PAW.
 if (Dtset%usepaw==1) then
   call pawrhoij_free(KS_Pawrhoij)
   ABI_DT_FREE(KS_Pawrhoij)
   call pawfgrtab_free(Pawfgrtab)
   call paw_ij_free(KS_paw_ij)
   call paw_an_free(KS_paw_an)
   call pawpwff_free(Paw_pwff)
   if (mod100>=10) then
     call pawrhoij_free(QP_pawrhoij)
     call paw_ij_free(QP_paw_ij)
     call paw_an_free(QP_paw_an)
   end if
   if (Dtset%pawcross==1) then
     call paw_pwaves_lmn_free(Paw_onsite)
     Wfdf%bks_comm = xmpi_comm_null
     call wfd_free(Wfdf)
   end if
 end if
 ABI_DT_FREE(Paw_onsite)
 ABI_DT_FREE(Pawfgrtab)
 ABI_DT_FREE(Paw_pwff)
 ABI_DT_FREE(KS_paw_ij)
 ABI_DT_FREE(KS_paw_an)
 if (mod100>=10) then
   ABI_DT_FREE(QP_pawrhoij)
   ABI_DT_FREE(QP_paw_an)
   ABI_DT_FREE(QP_paw_ij)
 end if

 call wfd_free(Wfd)
 call destroy_mpi_enreg(MPI_enreg_seq)
 call littlegroup_free(Ltg_k)
 ABI_DT_FREE(Ltg_k)
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call gsph_free(Gsph_Max)
 call gsph_free(Gsph_x)
 call gsph_free(Gsph_c)
 call vcoul_free(Vcp)
 call crystal_free(Cryst)
 call sigma_free(Sr)
 call em1results_free(Er)
 call ppm_free(PPm)
 call hdr_free(Hdr_sigma)
 call hdr_free(Hdr_wfk)
 call ebands_free(KS_BSt)
 call ebands_free(QP_BSt)
 call melements_free(KS_me)

 call esymm_free(KS_sym)
 ABI_DT_FREE(KS_sym)

 if (Sigp%symsigma==1.and.mod100>=20) then
   call esymm_free(QP_sym)
   ABI_DT_FREE(QP_sym)
 end if

 call sigparams_free(Sigp)

 call timab(426,2,tsec) ! finalize
 call timab(401,2,tsec)

 DBG_EXIT('COLL')

end subroutine sigma
!!***

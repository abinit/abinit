!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_sigma_driver
!! NAME
!!  m_sigma_driver
!!
!! FUNCTION
!! Calculate the matrix elements of the self-energy operator.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (MG, GMR, VO, LR, RWG, MT)
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

module m_sigma_driver

 use defs_basis
 use m_gwdefs
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_xmpi
 use m_xomp
 use m_errors
 use m_abicore
 use m_ab7_mixing
 use m_nctk
 use m_kxc
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use libxc_functionals
 use m_wfd

 use m_time,          only : timab,cwtime
 use m_numeric_tools, only : imax_loc
 use m_fstrings,      only : strcat, sjoin, itoa, basename, ktoa, ltoa
 use m_hide_blas,     only : xdotc
 use m_io_tools,      only : open_file, file_exists, iomode_from_fname
 use m_mpinfo,        only : destroy_mpi_enreg, initmpi_seq
 use m_geometry,      only : normv, mkrdim, metric
 use m_fftcore,       only : print_ngfft
 use m_fft_mesh,      only : get_gftt, setmesh
 use m_fft,           only : fourdp
 use m_ioarr,         only : fftdatar_write, read_rhor
 use m_crystal,       only : crystal_t, crystal_print
 use m_ebands,        only : ebands_update_occ, ebands_copy, ebands_report_gap, get_valence_idx, get_bandenergy, &
&                            ebands_free, ebands_init, ebands_ncwrite, ebands_interpolate_kpath, get_eneocc_vect, &
                             enclose_degbands, get_gaps, gaps_t
 use m_energies,      only : energies_type, energies_init
 use m_bz_mesh,       only : kmesh_t, kmesh_free, littlegroup_t, littlegroup_init, littlegroup_free, &
                             kmesh_init, has_BZ_item, isamek, get_ng0sh, kmesh_print, &
                             get_bz_item, has_IBZ_item, find_qmesh
 use m_gsphere,       only : gsphere_t, gsph_init, gsph_free, merge_and_sort_kg, gsph_extend, setshells
 use m_kg,            only : getph
 use m_xcdata,        only : get_xclevel
 use m_vcoul,         only : vcoul_t, vcoul_init, vcoul_free
 use m_qparticles,    only : wrqps, rdqps, rdgw, show_QP, updt_m_lda_to_qp
 use m_screening,     only : mkdump_er, em1results_free, epsilonm1_results, init_er_from_file
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
 use m_pawrhoij,      only : pawrhoij_type, pawrhoij_alloc, pawrhoij_copy, pawrhoij_free, &
&                            pawrhoij_inquire_dim, pawrhoij_symrhoij, pawrhoij_unpack
 use m_pawcprj,       only : pawcprj_type, pawcprj_alloc, pawcprj_free, paw_overlap
 use m_pawdij,        only : pawdij, symdij_all
 use m_pawfgr,        only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_pwaves_lmn,only : paw_pwaves_lmn_t, paw_pwaves_lmn_init, paw_pwaves_lmn_free
 use m_pawpwij,       only : pawpwff_t, pawpwff_init, pawpwff_free
 use m_paw_slater,    only : paw_mkdijexc_core, paw_dijhf
 use m_paw_dmft,      only : paw_dmft_type
 use m_paw_sphharm,   only : setsym_ylm
 use m_paw_mkrho,     only : denfgr
 use m_paw_nhat,      only : nhatgrid,pawmknhat
 use m_paw_tools,     only : chkpawovlp,pawprt
 use m_paw_denpot,    only : pawdenpot
 use m_paw_init,      only : pawinit,paw_gencond
 use m_classify_bands,only : classify_bands
 use m_wfk,           only : wfk_read_eigenvalues
 use m_io_kss,        only : make_gvec_kss
 use m_vhxc_me,       only : calc_vhxc_me
 use m_cohsex,        only : cohsex_me
 use m_sigx,          only : calc_sigx_me
 use m_sigc,          only : calc_sigc_me
 use m_setvtr,        only : setvtr
 use m_mkrho,         only : prtrhomxmn
 use m_pspini,        only : pspini
 use m_calc_ucrpa,    only : calc_ucrpa
 use m_prep_calc_ucrpa,only : prep_calc_ucrpa
 use m_paw_correlations,only : pawpuxinit
 use m_plowannier,only : operwan_realspace_type,plowannier_type,init_plowannier,get_plowannier,&
                         &fullbz_plowannier,init_operwan_realspace,destroy_operwan_realspace,&
                         &destroy_plowannier,zero_operwan_realspace

 implicit none

 private
!!***

 public :: sigma
!!***

contains
!!***

!!****f* m_sigma_driver/sigma
!! NAME
!! sigma
!!
!! FUNCTION
!! Calculate the matrix elements of the self-energy operator.
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
!!      classify_bands,cohsex_me,denfgr,destroy_mpi_enreg
!!      ebands_copy,ebands_free,ebands_interpolate_kpath,ebands_report_gap
!!      ebands_update_occ,em1results_free,energies_init,esymm_free
!!      fftdatar_write,fourdp,get_gftt,getph,gsph_free,hdr_free
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
!!      setsym_ylm,setup_ppmodel,setup_sigma,setvtr,show_qp,sigma_bksmask
!!      sigma_free,sigma_init,sigma_tables,sigparams_free,solve_dyson,symdij
!!      symdij_all,test_charge,timab,updt_m_lda_to_qp,vcoul_free
!!      wfd_change_ngfft,wfd_copy,wfd_distribute_bands,wfd_free,wfd_get_cprj
!!      wfd_init,wfd_mkrho,wfd_print,wfd_read_wfk,wfd_reset_ur_cprj,wfd_rotate
!!      wfd_test_ortho,write_sigma_header,write_sigma_results,wrqps,wrtout
!!      xmpi_barrier,xmpi_bcast,xmpi_sum
!!
!! SOURCE

subroutine sigma(acell,codvsn,Dtfil,Dtset,Pawang,Pawrad,Pawtab,Psps,rprim,converged)

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
 integer,parameter :: tim_fourdp5=5,master=0,cplex1=1
 integer :: approx_type,b1gw,b2gw,choice,cplex,cplex_dij,cplex_rhoij !,band
 integer :: dim_kxcg,gwcalctyp,gnt_option,has_dijU,has_dijso,iab,bmin,bmax,irr_idx1,irr_idx2
 integer :: iat,ib,ib1,ib2,ic,id_required,ider,idir,ii,ik,ierr,ount
 integer :: ik_bz,ikcalc,ik_ibz,ikxc,ipert,npw_k,omp_ncpus,pwx,ibz
 integer :: isp,is_idx,istep,itypat,itypatcor,izero,jj,first_band,last_band
 integer :: ks_iv,lcor,lmn2_size_max,mband,my_nband
 integer :: mgfftf,mod10,moved_atm_inside,moved_rhor,n3xccc !,mgfft
 integer :: nbsc,ndij,ndim,nfftf,nfftf_tot,nkcalc,gwc_nfft,gwc_nfftot,gwx_nfft,gwx_nfftot
 integer :: ngrvdw,nhatgrdim,nkxc,nkxc1,nprocs,nscf,nspden_rhoij,nzlmopt,optene
 integer :: optcut,optgr0,optgr1,optgr2,option,option_test,option_dij,optrad,optrhoij,psp_gencond
 integer :: my_rank,rhoxsp_method,comm,use_aerhor,use_umklp,usexcnhat
 integer :: ioe0j,spin,io,jb,nomega_sigc
 integer :: temp_unt,ncid,dim,nnn,dummy
 integer :: work_size,nstates_per_proc,my_nbks
 integer :: iatom1,iatom2,il1,il2,pos1,pos2,im1,im2,ispinor1,ispinor2
 !integer :: jb_qp,ib_ks,ks_irr
 real(dp) :: compch_fft,compch_sph,r_s,rhoav,alpha
 real(dp) :: drude_plsmf,my_plsmf,ecore,ecut_eff,ecutdg_eff,ehartree
 real(dp) :: ex_energy,gsqcutc_eff,gsqcutf_eff,gsqcut_shp,norm,oldefermi
 real(dp) :: ucvol,vxcavg,vxcavg_qp
 real(dp) :: gwc_gsq,gwx_gsq,gw_gsq,cpu,wall,gflops
 real(dp):: eff,mempercpu_mb,max_wfsmem_mb,nonscal_mem,ug_mem,ur_mem,cprj_mem
 complex(dpc) :: max_degw,cdummy,xx
 logical :: use_paw_aeur,dbg_mode,pole_screening,call_pawinit,is_dfpt=.false.
 character(len=500) :: msg
 character(len=fnlen) :: wfk_fname,pawden_fname
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
 !real(dp),allocatable :: osoc_bks(:, :, :)
 real(dp),allocatable :: ks_aepaw_rhor(:,:) !,ks_n_one_rhor(:,:),ks_nt_one_rhor(:,:)
 complex(dpc) :: ovlp(2)
 complex(dpc),allocatable :: ctmp(:,:),hbare(:,:,:,:)
 complex(dpc),target,allocatable :: sigcme(:,:,:,:,:)
 complex(dpc),allocatable :: hlda(:,:,:,:),htmp(:,:,:,:),uks2qp(:,:)
 complex(gwpc),allocatable :: kxcg(:,:),fxc_ADA(:,:,:)
 complex(gwpc),ABI_CONTIGUOUS pointer :: ug1(:)
!complex(dpc),pointer :: sigcme_p(:,:,:,:)
 complex(dpc),allocatable :: sigcme_k(:,:,:,:)
 complex(dpc), allocatable :: rhot1_q_m(:,:,:,:,:,:,:)
 complex(dpc), allocatable :: M1_q_m(:,:,:,:,:,:,:)
 complex(dpc),allocatable :: buffer(:)
 logical,allocatable :: bks_mask(:,:,:),keep_ur(:,:,:),bmask(:)
 type(esymm_t),target,allocatable :: KS_sym(:,:)
 type(esymm_t),pointer :: QP_sym(:,:)
 type(pawcprj_type),allocatable :: Cp1(:,:) !,Cp2(:,:)
 type(littlegroup_t),allocatable :: Ltg_k(:)
 type(Paw_an_type),allocatable :: KS_paw_an(:),QP_paw_an(:)
 type(Paw_ij_type),allocatable :: KS_paw_ij(:),QP_paw_ij(:)
 type(Pawfgrtab_type),allocatable :: Pawfgrtab(:)
 type(Pawrhoij_type),allocatable :: KS_Pawrhoij(:),QP_pawrhoij(:),prev_Pawrhoij(:),tmp_pawrhoij(:)
 type(pawpwff_t),allocatable :: Paw_pwff(:)
 type(paw_pwaves_lmn_t),allocatable :: Paw_onsite(:)
 type(plowannier_type) :: wanbz,wanibz
 type(operwan_realspace_type), allocatable :: rhot1(:,:)
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

 gwcalctyp=Dtset%gwcalctyp
 mod10 =MOD(Dtset%gwcalctyp,10)

 ! Perform some additional checks for hybrid functional calculations
 if(mod10==5) then
   if (Dtset%ixc_sigma<0 .and. .not.libxc_functionals_check()) then
     msg='Hybrid functional calculations require the compilation with LIBXC library'
     MSG_ERROR(msg)
   end if
!  XG 20171116 : I do not agree with this condition, as one might like to do a one-shot hybrid functional calculation
!  on top of a LDA/GGA calculation ... give the power (and risks) to the user !
!  if(gwcalctyp<10) then
!    msg='gwcalctyp should enforce update of the energies and/or wavefunctions when performing hybrid functional calculation'
!    MSG_ERROR(msg)
!  end if
   if(Dtset%usepaw==1) then
     msg='PAW version of hybrid functional calculations is not implemented'
     MSG_ERROR(msg)
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
 call pspini(Dtset,Dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcutf_eff,Pawrad,Pawtab,Psps,rprimd,comm_mpi=comm)

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
 if (sigma_needs_w(Sigp) .and. my_rank==master) then
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
   call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&              nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
   call pawrhoij_alloc(KS_Pawrhoij,cplex_rhoij,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)

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

   call setsym_ylm(gprimd,Pawang%l_max-1,Cryst%nsym,Dtset%pawprtvol,Cryst%rprimd,Cryst%symrec,Pawang%zarot)

   ! Initialize and compute data for LDA+U
   Paw_dmft%use_dmft=Dtset%usedmft
   call pawpuxinit(Dtset%dmatpuopt,Dtset%exchmix,Dtset%f4of2_sla,Dtset%f6of2_sla,&
&     is_dfpt,Dtset%jpawu,Dtset%lexexch,Dtset%lpawu,Cryst%ntypat,Pawang,Dtset%pawprtvol,&
&     Pawrad,Pawtab,Dtset%upawu,Dtset%usedmft,Dtset%useexexch,Dtset%usepawu,dtset%ucrpa)

   if (my_rank == master) call pawtab_print(Pawtab)

   ! Get Pawrhoij from the header of the WFK file.
   call pawrhoij_copy(Hdr_wfk%pawrhoij,KS_Pawrhoij)

   ! Re-symmetrize rhoij.
   ! FIXME this call leads to a SIGFAULT, likely some pointer is not initialized correctly
   choice=1; optrhoij=1; ipert=0; idir=0
!  call pawrhoij_symrhoij(KS_Pawrhoij,KS_Pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
!  &                      Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
!  &                      Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

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
 ! Initialize Wfd, allocate wavefunctions and precalculate tables to do the FFT using the coarse gwc_ngfft.
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

 call wfd_init(Wfd,Cryst,Pawtab,Psps,keep_ur,mband,nband,Kmesh%nibz,Sigp%nsppol,bks_mask,&
   Dtset%nspden,Dtset%nspinor,Dtset%ecutwfn,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,gwc_ngfft,&
   Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm)

 if (Dtset%pawcross==1) then
   call wfd_init(Wfdf,Cryst,Pawtab,Psps,keep_ur,mband,nband,Kmesh%nibz,Sigp%nsppol,bks_mask,&
     Dtset%nspden,Dtset%nspinor,dtset%ecutwfn,Dtset%ecutsm,Dtset%dilatmx,Hdr_wfk%istwfk,Kmesh%ibz,gwc_ngfft,&
     Dtset%nloalg,Dtset%prtvol,Dtset%pawprtvol,comm)
 end if

 ABI_FREE(bks_mask)
 ABI_FREE(nband)
 ABI_FREE(keep_ur)

 call timab(402,2,tsec) ! sigma(Init1)
 call timab(404,1,tsec) ! rdkss

 call wfd%read_wfk(wfk_fname,iomode_from_fname(wfk_fname))

 if (Dtset%pawcross==1) then
   call wfd_copy(Wfd,Wfdf)
   call wfdf%change_ngfft(Cryst,Psps,ngfftf)
 end if

 ! This test has been disabled (too expensive!)
 if (.False.) call wfd%test_ortho(Cryst,Pawtab,unit=ab_out,mode_paral="COLL")

 call timab(404,2,tsec) ! rdkss
 call timab(405,1,tsec) ! Init2

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

 if (Sigp%symsigma==1.and.gwcalctyp>=20) then
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

 !===========================
 !=== COMPUTE THE DENSITY ===
 !===========================
 ! * Evaluate the planewave part (complete charge in case of NC pseudos).
 ABI_MALLOC(ks_rhor,(nfftf,Dtset%nspden))
 ABI_MALLOC(ks_taur,(nfftf,Dtset%nspden*Dtset%usekden))

 call wfd%mkrho(Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_rhor)

 if (Dtset%usekden==1) then
   call wfd%mkrho(Cryst,Psps,Kmesh,KS_BSt,ngfftf,nfftf,ks_taur,optcalc=1)
 end if

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
   has_dijso=Dtset%pawspnorb; has_dijU=merge(0,1,Dtset%usepawu==0)
   call paw_ij_nullify(KS_paw_ij)
   call paw_ij_init(KS_paw_ij,cplex,Dtset%nspinor,Dtset%nsppol,&
&   Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=1,&
&   has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

   nkxc1=0
   ABI_DT_MALLOC(KS_paw_an,(Cryst%natom))
   call paw_an_nullify(KS_paw_an)
   call paw_an_init(KS_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,0,Dtset%nspden,&
&   cplex,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=1)
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
 call fourdp(1,ks_rhog,ks_rhor(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,tim_fourdp5)
 if(Dtset%usekden==1)then
   call fourdp(1,ks_taug,ks_taur(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,tim_fourdp5)
 end if

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
 if (Dtset%usepawu/=0    )  KS_mflags%has_vu     =1
 if (Dtset%useexexch/=0  )  KS_mflags%has_lexexch=1
 if (Sigp%use_sigxcore==1)  KS_mflags%has_sxcore =1
 if (gwcalctyp<10           )  KS_mflags%only_diago =1 ! off-diagonal elements only for SC on wavefunctions.

 if (.FALSE.) then ! quick and dirty hack to test HF contribution.
   MSG_WARNING("testing on-site HF")
   lmn2_size_max=MAXVAL(Pawtab(:)%lmn2_size)
   ABI_MALLOC(dij_hf,(cplex_dij*lmn2_size_max,ndij,Cryst%natom))
   call paw_dijhf(ndij,cplex_dij,1,lmn2_size_max,Cryst%natom,Cryst%ntypat,Pawtab,Pawrad,Pawang,&
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

 call calc_vhxc_me(Wfd,KS_mflags,KS_me,Cryst,Dtset,nfftf,ngfftf,&
& ks_vtrial,ks_vhartr,ks_vxc,Psps,Pawtab,KS_paw_an,Pawang,Pawfgrtab,KS_paw_ij,dijexc_core,&
& ks_rhor,usexcnhat,ks_nhat,ks_nhatgr,nhatgrdim,tmp_kstab,taug=ks_taug,taur=ks_taur)
 ABI_FREE(tmp_kstab)

!#ifdef DEV_HAVE_SCGW_SYM
!Set KS matrix elements connecting different irreps to zero. Do not touch unknown bands!.
 if (gwcalctyp>=20 .and. Sigp%symsigma > 0) then
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
 if (gwcalctyp>=10) then
   call wrtout(std_out,ch10//' *************** KS Energies *******************','COLL')
 end if

 !=== QP_BSt stores energies and occ. used for the calculation ===
 !  * Initialize QP_BSt with KS values.
 !  * In case of SC update QP_BSt using the QPS file.
 call ebands_copy(KS_BSt,QP_BSt)

 ABI_MALLOC(qp_rhor,(nfftf,Dtset%nspden))
 ABI_MALLOC(qp_taur,(nfftf,Dtset%nspden*Dtset%usekden))
 QP_sym => KS_sym

 if (gwcalctyp<10) then
   ! one-shot GW, just do a copy of the KS density.
   qp_rhor=ks_rhor
   if(Dtset%usekden==1)qp_taur=ks_taur
   QP_sym => KS_sym
 else
   ! Self-consistent GW.
   ! Read the unitary matrix and the QP energies of the previous step from the QPS file.
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
              nfftf,ngfftf,Cryst%ucvol,Cryst,Pawtab,MPI_enreg_seq,nbsc,Sr%m_lda_to_qp,prev_rhor,prev_Pawrhoij)

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

   if (nscf>0.and.gwcalctyp>=20.and.wfd%iam_master()) then
     ! Print the unitary transformation on std_out.
     call show_QP(QP_BSt,Sr%m_lda_to_qp,fromb=Sigp%minbdgw,tob=Sigp%maxbdgw,unit=std_out,tolmat=0.001_dp)
   end if

   ! Compute QP wfg as linear combination of KS states ===
   !  * Wfd%ug is modified inside calc_wf_qp
   !  * For PAW, update also the on-site projections.
   !  * WARNING the first dimension of MPI_enreg MUST be Kmesh%nibz
   !  TODO here we should use nbsc instead of nbnds

   call wfd%rotate(Cryst,Sr%m_lda_to_qp)

   ! Reinit the storage mode of Wfd as ug have been changed ===
   ! Update also the wavefunctions for GW corrections on each processor
   call wfd%reset_ur_cprj()

   ! This test has been disabled (too expensive!)
   if (.False.) call wfd%test_ortho(Cryst, Pawtab, unit=std_out)

   ! Compute QP occupation numbers.
   call wrtout(std_out,'sigma: calculating QP occupation numbers:','COLL')
   call ebands_update_occ(QP_BSt,Dtset%spinmagntarget,prtvol=0)
   qp_vbik(:,:) = get_valence_idx(QP_BSt)

!  #ifdef DEV_HAVE_SCGW_SYM
   ! Calculate the irreducible representations of the new QP amplitdues.
   if (Sigp%symsigma==1.and.gwcalctyp>=20) then
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
   call wfd%mkrho(Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_rhor)
   if (Dtset%usekden==1) then
     call wfd%mkrho(Cryst,Psps,Kmesh,QP_BSt,ngfftf,nfftf,qp_taur,optcalc=1)
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

     ! Calculate new QP quantities: nhat, nhatgr, rho_ij, paw_ij, and paw_an.
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
   call fourdp(1,qp_rhog,qp_rhor(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,tim_fourdp5)
   if(Dtset%usekden==1)call fourdp(1,qp_taug,qp_taur(:,1),-1,MPI_enreg_seq,nfftf,1,ngfftf,tim_fourdp5)

   ! ===========================================
   ! ==== Optional output of the QP density ====
   ! ===========================================
   if (Dtset%prtden/=0.and.wfd%iam_master()) then
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

   if (allocated(qp_kxc)) then
     ABI_FREE(qp_kxc)
   end if

   if (Dtset%usepaw==1) then
     call timab(561,1,tsec)

     ! Compute QP Dij
     ipert=0; idir=0
     call pawdij(cplex,Dtset%enunit,Cryst%gprimd,ipert,&
&     Cryst%natom,Cryst%natom,nfftf,ngfftf(1)*ngfftf(2)*ngfftf(3),&
&     Dtset%nspden,Cryst%ntypat,QP_paw_an,QP_paw_ij,Pawang,Pawfgrtab,&
&     Dtset%pawprtvol,Pawrad,QP_pawrhoij,Dtset%pawspnorb,Pawtab,Dtset%pawxcdev,&
&     k0,Dtset%spnorbscl,Cryst%ucvol,dtset%charge,qp_vtrial,qp_vxc,Cryst%xred,&
&     nucdipmom=Dtset%nucdipmom)

     ! Symmetrize total Dij
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
       call wfd%distribute_bands(ik_ibz, spin, my_nband, my_band_list, bmask=bmask)
       if (my_nband == 0) cycle
       npw_k = Wfd%npwarr(ik_ibz)
       do ii=1,my_nband  ! ib=b1gw,b2gw in sequential
         ib = my_band_list(ii)
         ug1  => Wfd%Wave(ib,ik_ibz,spin)%ug
         cdummy = xdotc(npw_k*Wfd%nspinor,ug1,1,ug1,1)
         ovlp(1) = REAL(cdummy); ovlp(2) = AIMAG(cdummy)

         if (Psps%usepaw==1) then
           call wfd%get_cprj(ib,ik_ibz,spin,Cryst,Cp1,sorted=.FALSE.)
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

 ! Setup of the bare Hamiltonian := T + v_{loc} + v_{nl} + v_H.
 ! * The representation depends wheter we are updating the wfs or not.
 ! * ks_vUme is zero unless we are using LDA+U as starting point, see calc_vHxc_braket
 ! * Note that vH matrix elements are calculated using the true uncutted interaction.

 if (gwcalctyp<10) then
   ! For one-shot GW use the KS representation.
   Sr%hhartree=hlda-KS_me%vxcval
   ! Additional goodies for PAW
   !  * LDA +U Hamiltonian
   !  * LEXX.
   !  * Core contribution estimated using Fock exchange.
   if (Dtset%usepaw==1) then
     if (Sigp%use_sigxcore==1) Sr%hhartree=hlda - (KS_me%vxc - KS_me%sxcore)
     if (Dtset%usepawu/=0) Sr%hhartree=Sr%hhartree-KS_me%vu
     if (Dtset%useexexch/=0) then
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

     ! Change basis from KS to QP, hbare is overwritten: A_{QP} = U^\dagger A_{KS} U
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

   ! Calculate the QP matrix elements
   ! This part is parallelized within MPI_COMM_WORD since each node has all GW wavefunctions.
   ! For PAW, we have to construct the new bare Hamiltonian.
   call wrtout(std_out,ch10//' *************** QP Energies *******************','COLL')

   call melflags_reset(QP_mflags)
!  if (gwcalctyp<20) QP_mflags%only_diago=1 ! For e-only, no need of off-diagonal elements.
   QP_mflags%has_vhartree=1
   if (Dtset%usepaw==1)    QP_mflags%has_hbare  =1
!  QP_mflags%has_vxc     =1
!  QP_mflags%has_vxcval  =1
!  if (Sigp%gwcalctyp >100) QP_mflags%has_vxcval_hybrid=1
   if (mod10==5 .and. &
&   (Dtset%ixc_sigma==-402 .or. Dtset%ixc_sigma==-406 .or. Dtset%ixc_sigma==-427 .or. Dtset%ixc_sigma==-428 .or.&
&   Dtset%ixc_sigma==-456 .or. Dtset%ixc_sigma==41 .or. Dtset%ixc_sigma==42)) then
     QP_mflags%has_vxcval_hybrid=1
   end if
!  if (Sigp%use_sigxcore==1) QP_mflags%has_sxcore =1
!  if (Dtset%usepawu/=0)    QP_mflags%has_vu     =1
!  if (Dtset%useexexch/=0)  QP_mflags%has_lexexch=1

   ABI_MALLOC(tmp_kstab,(2,Wfd%nkibz,Wfd%nsppol))
   tmp_kstab=0
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw ! No spin dependent!
       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
       tmp_kstab(1,ik_ibz,spin)=Sigp%minbnd(ikcalc,spin)
       tmp_kstab(2,ik_ibz,spin)=Sigp%maxbnd(ikcalc,spin)
     end do
   end do

   call calc_vhxc_me(Wfd,QP_mflags,QP_me,Cryst,Dtset,nfftf,ngfftf,&
&   qp_vtrial,qp_vhartr,qp_vxc,Psps,Pawtab,QP_paw_an,Pawang,Pawfgrtab,QP_paw_ij,dijexc_core,&
&   qp_rhor,usexcnhat,qp_nhat,qp_nhatgr,nhatgrdim,tmp_kstab,taug=qp_taug,taur=qp_taur)
   ABI_FREE(tmp_kstab)

!  #ifdef DEV_HAVE_SCGW_SYM
   if (gwcalctyp>=20 .and. Sigp%symsigma>0) then
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
   if (gwcalctyp>=20 .and. Sigp%symsigma > 0) then
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

 ! Prepare the storage of QP amplitudes and energies ===
 ! Initialize with KS wavefunctions and energies.
 Sr%eigvec_qp=czero;  Sr%en_qp_diago=zero
 do ib=1,Sigp%nbnds
   Sr%en_qp_diago(ib,:,:)=KS_BSt%eig(ib,:,:)
   Sr%eigvec_qp(ib,ib,:,:)=cone
 end do

 ! Store <n,k,s|V_xc[n_val]|n,k,s> and <n,k,s|V_U|n,k,s> ===
 ! Note that we store the matrix elements of V_xc in the KS basis set, not in the QP basis set
 ! Matrix elements of V_U are zero unless we are using LDA+U as starting point
 do ib=b1gw,b2gw
   Sr%vxcme(ib,:,:)=KS_me%vxcval(ib,ib,:,:)
   if (Dtset%usepawu/=0) Sr%vUme (ib,:,:)=KS_me%vu(ib,ib,:,:)
 end do

 ! Initial guess for the GW energies
 ! Save the energies of the previous iteration.
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
 end if ! sigma_needs_ppm

 call timab(409,2,tsec) ! getW

 if (wfd%iam_master()) then
   ! Write info on the run on ab_out, then open files to store final results.
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

   if (mod10==SIG_GW_AC) then
     ! Sigma along the imaginary axis.
     if (open_file(Dtfil%fnameabo_sgm,msg,unit=unt_sgm,status='unknown',form='formatted') /= 0) then
       MSG_ERROR(msg)
     end if
   end if
 end if

!=======================================================================
!==== Calculate self-energy and output the results for each k-point ====
!=======================================================================
!* Here it would be possible to calculate the QP correction for the same k-point using a PPmodel
!in the first iteration just to improve the initial guess and CD or AC in the second step. Useful
!if the KS level is far from the expected QP results. Previously it was possible to do such
!calculation by simply listing the same k-point twice in kptgw. Now this trick is not allowed anymore.
!Everything, indeed, should be done in a clean and transparent way inside csigme.

 call wfd%print(mode_paral='PERS')

 call wrtout(std_out,sigma_type_from_key(mod10),'COLL')

 if (gwcalctyp<10) then
   msg = " Perturbative Calculation"
   if (gwcalctyp==1) msg = " Newton Raphson method "
 else if (gwcalctyp<20) then
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

 ! min and max band indeces for GW corrections.
 ib1=Sigp%minbdgw; ib2=Sigp%maxbdgw

 !MG TODO: I don't like the fact that ib1 and ib2 are redefined here because this
 ! prevents me from refactoring the code. In particular I want to store the self-energy
 ! results inside the sigma_results datatypes hence one needs to know all the dimensions
 ! at the beginning of the execution (e.g. in setup_sigma) so that one can easily allocate the arrays in the type.
 
 if(Dtset%ucrpa>=1) then
   !Read the band
   if (dtset%plowan_compute<10)then
     if (open_file("forlb.ovlp",msg,newunit=temp_unt,form="formatted", status="unknown") /= 0) then
       MSG_ERROR(msg)
     end if
     rewind(temp_unt)
     read(temp_unt,*)
     read(temp_unt,*) msg, ib1, ib2
     close(temp_unt)
   else
   ib1=dtset%plowan_bandi
   ib2=dtset%plowan_bandf
 endif
endif

 ABI_MALLOC(sigcme,(nomega_sigc,ib1:ib2,ib1:ib2,Sigp%nkptgw,Sigp%nsppol*Sigp%nsig_ab))
 sigcme=czero

! if (.False. .and. psps%usepaw == 0 .and. wfd%nspinor == 1 .and. any(dtset%so_psp /= 0)) then
!   call wrtout(std_out, "Computing SOC contribution with first-order perturbation theory")
!   ABI_MALLOC(bks_mask, (wfd%mband, wfd%nkibz, wfd%nsppol))
!   bks_mask = .False.
!   do spin=1,wfd%nsppol
!     do ikcalc=1,Sigp%nkptgw
!       ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Irred k-point for GW
!       ii=Sigp%minbnd(ikcalc, spin); jj=Sigp%maxbnd(ikcalc, spin)
!       bks_mask(ii:jj, ik_ibz, spin) = .True.
!     end do
!   end do
!
!   call wfd_get_socpert(wfd, cryst, psps, pawtab, bks_mask, osoc_bks)
!   ABI_FREE(bks_mask)
!   ABI_FREE(osoc_bks)
! end if

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
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Initialization of wan objects and getting psichies
!Allocation of rhot1
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   if (dtset%plowan_compute>=10)then
     write(msg,'(a)')" cRPA calculations using wannier weights from data.plowann"
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     call init_plowannier(dtset%plowan_bandf,dtset%plowan_bandi,dtset%plowan_compute,dtset%plowan_iatom,&
       &dtset%plowan_it,dtset%plowan_lcalc,dtset%plowan_natom,dtset%plowan_nbl,dtset%plowan_nt,&
       &dtset%plowan_projcalc,dtset%acell_orig,dtset%kpt,dtset%nimage,dtset%nkpt,dtset%nspinor,&
       &dtset%nsppol,dtset%wtk,wanibz)
     call get_plowannier(wanibz)
     call fullbz_plowannier(dtset,kmesh,cryst,wanibz,wanbz)
     ABI_DATATYPE_ALLOCATE(rhot1,(sigp%npwx,Qmesh%nibz))
     do pwx=1,sigp%npwx
       do ibz=1,Qmesh%nibz
         call init_operwan_realspace(wanbz,rhot1(pwx,ibz))
         call zero_operwan_realspace(wanbz,rhot1(pwx,ibz))
       enddo
     enddo
   endif



!   do ikcalc=1,Sigp%nkptgw

!   if(cryst%nsym==1) nkcalc=Kmesh%nbz
!   if(cryst%nsym>1)  nkcalc=Kmesh%nibz
   nkcalc=Kmesh%nbz
   !if(1==1)then!DEBUG
   !open(67,file="test.rhot1",status="REPLACE")
   do ikcalc=1,nkcalc ! for the oscillator strengh, spins are identical without SOC
!    if(Sigp%nkptgw/=Kmesh%nbz) then
!      write(msg,'(6a)')ch10,&
!&      ' nkptgw and nbz differs: this is not allowed to compute U in cRPA '
!      MSG_ERROR(msg)
!    endif
     write(std_out,*) " ikcalc",ikcalc,size(Sigp%kptgw2bz),nkcalc

     if(mod(ikcalc-1,nprocs)==Wfd%my_rank) then
       if(cryst%nsym==1) ik_ibz=Kmesh%tab(ikcalc) ! Index of the irred k-point for GW
       if(cryst%nsym>1) ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW

       call prep_calc_ucrpa(ik_ibz,ikcalc,itypatcor,ib1,ib2,Cryst,QP_bst,Sigp,Gsph_x,Vcp,&
&       Kmesh,Qmesh,lcor,M1_q_m,Pawtab,Pawang,Paw_pwff,&
&       Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,gwx_ngfft,&
&       ngfftf,Dtset%prtvol,Dtset%pawcross,Dtset%plowan_compute,rhot1_q_m,wanbz,rhot1)
     end if
   end do
   if (dtset%plowan_compute<10)then
     call xmpi_sum(rhot1_q_m,Wfd%comm,ierr)
     call xmpi_sum(M1_q_m,Wfd%comm,ierr)
     M1_q_m=M1_q_m/Kmesh%nbz/Wfd%nsppol
     rhot1_q_m=rhot1_q_m/Kmesh%nbz/Wfd%nsppol
     ! do pwx=1,sigp%npwx
     !   do ibz=1,Qmesh%nibz
     !     do ispinor1=1,dtset%nspinor
     !       do ispinor2=1,dtset%nspinor
     !         do iatom1=1,cryst%nattyp(itypatcor)
     !           do im1=1,2*lcor+1
     !             do im2=1,2*lcor+1
     !               write(67,*)ibz,im1,im2,rhot1_q_m(iatom1,ispinor1,ispinor2,im1,im2,pwx,ibz)
     !             enddo
     !           enddo
     !         enddo
     !       enddo
     !     enddo
     !   enddo
     ! enddo
   else
     !call cwtime(cpu,wall,gflops,"start") !reduction of rhot1
     dim=0
     do pwx=1,sigp%npwx
     do ibz=1,Qmesh%nibz
       do spin=1,wanbz%nsppol
       do ispinor1=1,wanbz%nspinor
       do ispinor2=1,wanbz%nspinor
         do iatom1=1,wanbz%natom_wan
         do iatom2=1,wanbz%natom_wan
           do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
           do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
             do il1=1,wanbz%nbl_atom_wan(iatom1)
             do il2=1,wanbz%nbl_atom_wan(iatom2)
               do im1=1,2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
               do im2=1,2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
     dim=dim+1
               enddo!im2
               enddo!im1
             enddo!il2
             enddo!il1
           enddo!pos2
           enddo!pos1
         enddo!iatom2
         enddo!iatom1
       enddo!ispinor2
       enddo!ispinor1
       enddo!spin
     enddo!ibz
     enddo!pwx
     ABI_ALLOCATE(buffer,(dim))
     nnn=0
     do pwx=1,sigp%npwx
     do ibz=1,Qmesh%nibz
         do spin=1,wanbz%nsppol
         do ispinor1=1,wanbz%nspinor
         do ispinor2=1,wanbz%nspinor
           do iatom1=1,wanbz%natom_wan
           do iatom2=1,wanbz%natom_wan
             do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
             do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
               do il1=1,wanbz%nbl_atom_wan(iatom1)
               do il2=1,wanbz%nbl_atom_wan(iatom2)
                 do im1=1,2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
                 do im2=1,2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
     nnn=nnn+1
     buffer(nnn)=rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,spin,ispinor1,ispinor2)
                 enddo!im2
                 enddo!im1
               enddo!il2
               enddo!il1
             enddo!pos2
             enddo!pos1
           enddo!iatom2
           enddo!iatom1
         enddo!ispinor2
         enddo!ispinor1
       enddo!spin
     enddo!ibz  
     enddo!pwx
     call xmpi_barrier(Wfd%comm)
     call xmpi_sum(buffer,Wfd%comm,ierr)
     call xmpi_barrier(Wfd%comm)
     buffer=buffer/Kmesh%nbz/Wfd%nsppol
     nnn=0
     do pwx=1,sigp%npwx
     do ibz=1,Qmesh%nibz
       do spin=1,wanbz%nsppol
       do ispinor1=1,wanbz%nspinor
       do ispinor2=1,wanbz%nspinor
         do iatom1=1,wanbz%natom_wan
         do iatom2=1,wanbz%natom_wan
           do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
           do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
             do il1=1,wanbz%nbl_atom_wan(iatom1)
             do il2=1,wanbz%nbl_atom_wan(iatom2)
               do im1=1,2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
               do im2=1,2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
      nnn=nnn+1
      rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,spin,ispinor1,ispinor2)=buffer(nnn)
      !write(67,*)ibz,im1,im2,rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,spin,ispinor1,ispinor2)
               enddo!im2
               enddo!im1
             enddo!il2
             enddo!il1
           enddo!pos2
           enddo!pos1
         enddo!iatom2
         enddo!iatom1
       enddo!ispinor2
       enddo!ispinor1
       enddo!spin
     enddo!ibz
     enddo!pwx
     ABI_DEALLOCATE(buffer)
   endif
   
   !close(67)
   ! call xmpi_barrier(Wfd%comm)
   ! call cwtime(cpu,wall,gflops,"stop")!reduction of rhot1
   ! write(6,*)cpu,wall,gflops
!   if(Cryst%nsym==1) then
!   M1_q_m=M1_q_m/Kmesh%nbz/Wfd%nsppol
!   rhot1_q_m=rhot1_q_m/Kmesh%nbz/Wfd%nsppol
!   endif
!  Calculation of U in cRPA: need to treat with a different cutoff the
!  bare coulomb and the screened coulomb interaction.
   ! else !DEBUG
   !   write(6,*)"DEBUGGING calc_ucrpa"
   !   open(67,file="test.rhot1",status="OLD")
   !   rewind(67)
   !    do pwx=1,sigp%npwx
   !     do ibz=1,Qmesh%nibz
   !       do spin=1,wanbz%nsppol
   !       do ispinor1=1,wanbz%nspinor
   !         do ispinor2=1,wanbz%nspinor
   !           do iatom1=1,wanbz%natom_wan
   !             do iatom2=1,wanbz%natom_wan
   !               do pos1=1,size(wanbz%nposition(iatom1)%pos,1)
   !                 do pos2=1,size(wanbz%nposition(iatom2)%pos,1)
   !                   do il1=1,wanbz%nbl_atom_wan(iatom1)
   !                     do il2=1,wanbz%nbl_atom_wan(iatom2)
   !                       do im1=1,2*wanbz%latom_wan(iatom1)%lcalc(il1)+1
   !                         do im2=1,2*wanbz%latom_wan(iatom2)%lcalc(il2)+1
   !    read(67,*)dummy,dummy,dummy,xx
   !    rhot1(pwx,ibz)%atom_index(iatom1,iatom2)%position(pos1,pos2)%atom(il1,il2)%matl(im1,im2,spin,ispinor1,ispinor2)=xx
   !                         enddo!im2
   !                       enddo!im1
   !                     enddo!il2
   !                   enddo!il1
   !                 enddo!pos2
   !               enddo!pos1
   !             enddo!iatom2
   !           enddo!iatom1
   !         enddo!ispinor2
   !       enddo!ispinor1
   !       enddo!spin
   !     enddo!ibz
   !   enddo!pwx
   !   close(67)
   ! endif!DEBUG
   call calc_ucrpa(itypatcor,cryst,Kmesh,lcor,M1_q_m,Qmesh,Er%npwe,sigp%npwx,&
&   Cryst%nsym,rhot1_q_m,Sigp%nomegasr,Sigp%minomega_r,Sigp%maxomega_r,ib1,ib2,'Gsum',Cryst%ucvol,Wfd,Er%fname,dtset%plowan_compute,rhot1,wanbz)

!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!Deallocation of wan and rhot1
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
   if (dtset%plowan_compute >=10) then
     do pwx=1,sigp%npwx
       do ibz=1,Qmesh%nibz
         call destroy_operwan_realspace(wanbz,rhot1(pwx,ibz))
       enddo
     enddo
     ABI_DATATYPE_DEALLOCATE(rhot1)
     call destroy_plowannier(wanbz)
   endif
   ABI_DEALLOCATE(rhot1_q_m)
   ABI_DEALLOCATE(M1_q_m)
   
 else
   do ikcalc=1,Sigp%nkptgw
     ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
     ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indices for GW corrections (for this k-point)
     ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

     call calc_sigx_me(ik_ibz,ikcalc,ib1,ib2,Cryst,QP_bst,Sigp,Sr,Gsph_x,Vcp,Kmesh,Qmesh,Ltg_k(ikcalc),&
&     Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&
&     gwx_ngfft,ngfftf,Dtset%prtvol,Dtset%pawcross)
   end do

   ! for the time being, do not remove this barrier!
   call xmpi_barrier(Wfd%comm)

   call timab(421,2,tsec) ! calc_sigx_me

   ! ==========================================================
   ! ==== Correlation part using the coarse gwc_ngfft mesh ====
   ! ==========================================================
   if (mod10/=SIG_HF) then
     do ikcalc=1,Sigp%nkptgw
       ik_ibz=Kmesh%tab(Sigp%kptgw2bz(ikcalc)) ! Index of the irred k-point for GW
       ib1=MINVAL(Sigp%minbnd(ikcalc,:)) ! min and max band indices for GW corrections (for this k-point)
       ib2=MAXVAL(Sigp%maxbnd(ikcalc,:))

!      sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:) ! Causes annoying Fortran runtime warning on abiref.
       ABI_ALLOCATE(sigcme_k,(nomega_sigc,ib2-ib1+1,ib2-ib1+1,Sigp%nsppol*Sigp%nsig_ab))
       sigcme_k=czero
       if (any(mod10 == [SIG_SEX, SIG_COHSEX])) then
         ! Calculate static COHSEX or SEX using the coarse gwc_ngfft mesh.
         call cohsex_me(ik_ibz,ikcalc,nomega_sigc,ib1,ib2,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_c,Vcp,Kmesh,Qmesh,&
&         Ltg_k(ikcalc),Pawtab,Pawang,Paw_pwff,Psps,Wfd,QP_sym,&
!&         gwc_ngfft,Dtset%iomode,Dtset%prtvol,sigcme_p)
&         gwc_ngfft,Dtset%iomode,Dtset%prtvol,sigcme_k)
       else
         ! Compute correlated part using the coarse gwc_ngfft mesh.
         call calc_sigc_me(ik_ibz,ikcalc,nomega_sigc,ib1,ib2,Dtset,Cryst,QP_BSt,Sigp,Sr,Er,Gsph_Max,Gsph_c,Vcp,Kmesh,Qmesh,&
&         Ltg_k(ikcalc),PPm,Pawtab,Pawang,Paw_pwff,Pawfgrtab,Paw_onsite,Psps,Wfd,Wfdf,QP_sym,&
!&         gwc_ngfft,ngfftf,nfftf,ks_rhor,use_aerhor,ks_aepaw_rhor,sigcme_p)
&         gwc_ngfft,ngfftf,nfftf,ks_rhor,use_aerhor,ks_aepaw_rhor,sigcme_k)
       end if
       sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)=sigcme_k
       ABI_DEALLOCATE(sigcme_k)

     end do
   end if

   call xmpi_barrier(Wfd%comm)

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

!    sigcme_p => sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)   ! Causes annoying Fortran runtime warning on abiref.
     ABI_ALLOCATE(sigcme_k,(nomega_sigc,ib2-ib1+1,ib2-ib1+1,Sigp%nsppol*Sigp%nsig_ab))
     sigcme_k=sigcme(:,ib1:ib2,ib1:ib2,ikcalc,:)

!    call solve_dyson(ikcalc,ib1,ib2,nomega_sigc,Sigp,Kmesh,sigcme_p,QP_BSt%eig,Sr,Dtset%prtvol,Dtfil,Wfd%comm)
     call solve_dyson(ikcalc,ib1,ib2,nomega_sigc,Sigp,Kmesh,sigcme_k,QP_BSt%eig,Sr,Dtset%prtvol,Dtfil,Wfd%comm)
     ABI_DEALLOCATE(sigcme_k)
     !
     ! Calculate direct gap for each spin and print out final results.
     ! We use the valence index of the KS system because we still do not know
     ! all the QP corrections. Ideally one should use the QP valence index
     do spin=1,Sigp%nsppol
       if ( Sigp%maxbnd(ikcalc,spin) >= ks_vbik(ik_ibz,spin)+1 .and. &
&       Sigp%minbnd(ikcalc,spin) <= ks_vbik(ik_ibz,spin)  ) then
         ks_iv=ks_vbik(ik_ibz,spin)
         Sr%egwgap (ik_ibz,spin)= Sr%egw(ks_iv+1,ik_ibz,spin) -  Sr%egw(ks_iv,ik_ibz,spin)
         Sr%degwgap(ik_ibz,spin)= Sr%degw(ks_iv+1,ik_ibz,spin) - Sr%degw(ks_iv,ik_ibz,spin)
       else
         ! The "gap" cannot be computed
         Sr%e0gap(ik_ibz,spin)=zero; Sr%egwgap (ik_ibz,spin)=zero; Sr%degwgap(ik_ibz,spin)=zero
       end if
     end do

     if (wfd%iam_master()) call write_sigma_results(ikcalc,ik_ibz,Sigp,Sr,KS_BSt,dtset%use_yaml)
   end do !ikcalc

   call timab(425,2,tsec) ! solve_dyson
   call timab(426,1,tsec) ! finalize

   ! Update the energies in QP_BSt
   ! If QPSCGW, use diagonalized eigenvalues otherwise perturbative results.
   if (gwcalctyp>=10) then
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

     ! Report the QP gaps (Fundamental and direct)
     call ebands_report_gap(QP_BSt,header='QP Band Gaps',unit=ab_out)

     ! Band structure interpolation from QP energies computed on the k-mesh.
     if (nint(dtset%einterp(1)) /= 0 .and. all(sigp%minbdgw == sigp%minbnd) .and. all(sigp%maxbdgw == sigp%maxbnd)) then
       call ebands_interpolate_kpath(QP_BSt, dtset, cryst, [sigp%minbdgw, sigp%maxbdgw], dtfil%filnam_ds(4), comm)
     end if
   end if ! Sigp%nkptgw==Kmesh%nibz
   !
   ! Write SCF data in case of self-consistent calculation ===
   ! Save Sr%en_qp_diago, Sr%eigvec_qp and m_lda_to_qp in the _QPS file.
   ! Note that in the first iteration qp_rhor contains KS rhor, then the mixed rhor.
   if (gwcalctyp>=10) then
     ! Calculate the new m_lda_to_qp
     call updt_m_lda_to_qp(Sigp,Kmesh,nscf,Sr,Sr%m_lda_to_qp)

     if (wfd%iam_master()) then
       call wrqps(Dtfil%fnameabo_qps,Sigp,Cryst,Kmesh,Psps,Pawtab,QP_Pawrhoij,&
&       Dtset%nspden,nscf,nfftf,ngfftf,Sr,QP_BSt,Sr%m_lda_to_qp,qp_rhor)
     end if

     ! Report the MAX variation for each kptgw and spin
     call wrtout(ab_out,ch10//' Convergence of QP corrections ','COLL')
     converged=.TRUE.
     do spin=1,Sigp%nsppol
       write(msg,'(a,i2,a)')' >>>>> For spin ',spin,' <<<<< '
       call wrtout(ab_out,msg,'COLL')
       do ikcalc=1,Sigp%nkptgw
         ib1 = Sigp%minbnd(ikcalc,spin); ib2 = Sigp%maxbnd(ikcalc,spin)
         ik_bz = Sigp%kptgw2bz(ikcalc); ik_ibz = Kmesh%tab(ik_bz)
         ii = imax_loc(ABS(Sr%degw(ib1:ib2,ik_ibz,spin)))
         max_degw = Sr%degw(ii,ik_ibz,spin)
         write(msg,('(a,i3,a,2f8.3,a,i3)'))'.  kptgw no:',ikcalc,'; Maximum DeltaE = (',max_degw*Ha_eV,') for band index:',ii
         call wrtout(ab_out,msg,'COLL')
         converged = (converged .and. ABS(max_degw) < Dtset%gw_toldfeig)
       end do
     end do
   end if

   ! ============================================
   ! ==== Save the GW results in NETCDF file ====
   ! ============================================
#ifdef HAVE_NETCDF
   if (wfd%iam_master()) then
     NCF_CHECK(nctk_open_create(ncid, strcat(dtfil%filnam_ds(4), '_SIGRES.nc'), xmpi_comm_self))
     NCF_CHECK(nctk_defnwrite_ivars(ncid, ["sigres_version"], [1]))
     NCF_CHECK(cryst%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(KS_Bst, ncid))
     NCF_CHECK(sigma_ncwrite(Sigp, Er, Sr, ncid))
     ! Add qp_rhor. Note that qp_rhor == ks_rhor if wavefunctions are not updated.
     !ncerr = nctk_write_datar("qp_rhor",path,ngfft,cplex,nfft,nspden,&
     !                          comm_fft,fftn3_distrib,ffti3_local,datar,action)
     NCF_CHECK(nf90_close(ncid))
   end if
#endif

 end if ! ucrpa

 !----------------------------- END OF THE CALCULATION ------------------------
 !
 !=====================
 !==== Close Files ====
 !=====================
 if (wfd%iam_master()) then
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
   if (gwcalctyp>=10) then
     call pawrhoij_free(QP_pawrhoij)
     call paw_ij_free(QP_paw_ij)
     call paw_an_free(QP_paw_an)
   end if
   if (Dtset%pawcross==1) then
     call paw_pwaves_lmn_free(Paw_onsite)
     Wfdf%bks_comm = xmpi_comm_null
     call wfdf%free()
   end if
 end if
 ABI_DT_FREE(Paw_onsite)
 ABI_DT_FREE(Pawfgrtab)
 ABI_DT_FREE(Paw_pwff)
 ABI_DT_FREE(KS_paw_ij)
 ABI_DT_FREE(KS_paw_an)
 if (gwcalctyp>=10) then
   ABI_DT_FREE(QP_pawrhoij)
   ABI_DT_FREE(QP_paw_an)
   ABI_DT_FREE(QP_paw_ij)
 end if

 call wfd%free()
 call destroy_mpi_enreg(MPI_enreg_seq)
 call littlegroup_free(Ltg_k)
 ABI_DT_FREE(Ltg_k)
 call kmesh_free(Kmesh)
 call kmesh_free(Qmesh)
 call gsph_free(Gsph_Max)
 call gsph_free(Gsph_x)
 call gsph_free(Gsph_c)
 call vcoul_free(Vcp)
 call cryst%free()
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

 if (Sigp%symsigma==1.and.gwcalctyp>=20) then
   call esymm_free(QP_sym)
   ABI_DT_FREE(QP_sym)
 end if

 call sigparams_free(Sigp)

 call timab(426,2,tsec) ! finalize
 call timab(401,2,tsec)

 DBG_EXIT('COLL')

end subroutine sigma
!!***

!!****f* ABINIT/setup_sigma
!! NAME
!! setup_sigma
!!
!! FUNCTION
!!  Initialize the data type containing parameters for a sigma calculation.
!!
!! INPUTS
!! acell(3)=length scales of primitive translations (bohr)
!! wfk_fname=Name of the WFK file.
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Dtfil<type(datafiles_type)>=variables related to files
!! rprim(3,3)=dimensionless real space primitive translations
!! ngfft(18)=information on the (fine) FFT grid used for the density.
!! Psps <Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the WFK file
!!
!! OUTPUT
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! Qmesh <kmesh_t>=Structure describing the q-point sampling.
!! Cryst<crystal_t>=Info on unit cell and symmetries.
!! Gsph_Max<gsphere_t>=Info on the G-sphere
!! Gsph_c<gsphere_t>=Info on the G-sphere for W and Sigma_c
!! Gsph_x<gsphere_t>=Info on the G-sphere for and Sigma_x
!! Hdr_wfk<hdr_type>=The header of the WFK file
!! Hdr_out<hdr_type>=The header to be used for the results of sigma calculations.
!! Vcp<vcoul_t>= Datatype gathering information on the coulombian interaction and the cutoff technique.
!! Er<Epsilonm1_results>=Datatype storing data used to construct the screening (partially Initialized in OUTPUT)
!! KS_BSt<ebands_t>=The KS energies and occupation factors.
!! gwc_ngfft(18), gwx_ngfft(18)= FFT meshes for the oscillator strengths used for the correlated and the
!!   exchange part of the self-energy, respectively.
!! comm=MPI communicator.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine setup_sigma(codvsn,wfk_fname,acell,rprim,ngfftf,Dtset,Dtfil,Psps,Pawtab,&
& gwx_ngfft,gwc_ngfft,Hdr_wfk,Hdr_out,Cryst,Kmesh,Qmesh,KS_BSt,Gsph_Max,Gsph_x,Gsph_c,Vcp,Er,Sigp,comm)

 !use m_gwdefs,        only : GW_Q0_DEFAULT, SIG_GW_AC, sigparams_t, sigma_is_herm, sigma_needs_w

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: comm
 character(len=6),intent(in) :: codvsn
 character(len=*),intent(in) :: wfk_fname
 type(Datafiles_type),intent(in) :: Dtfil
 type(Dataset_type),intent(inout) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Dtset%usepaw)
 type(sigparams_t),intent(out) :: Sigp
 type(Epsilonm1_results),intent(out) :: Er
 type(ebands_t),intent(out) :: KS_BSt
 type(kmesh_t),intent(out) :: Kmesh,Qmesh
 type(crystal_t),intent(out) :: Cryst
 type(gsphere_t),intent(out) :: Gsph_Max,Gsph_x,Gsph_c
 type(Hdr_type),intent(out) :: Hdr_wfk,Hdr_out
 type(vcoul_t),intent(out) :: Vcp
!arrays
 integer,intent(in) :: ngfftf(18)
 integer,intent(out) :: gwc_ngfft(18),gwx_ngfft(18)
 real(dp),intent(in) :: acell(3),rprim(3,3)

!Local variables-------------------------------
!scalars
 integer,parameter :: pertcase0=0,master=0
 integer :: bantot,enforce_sym,gwcalctyp,ib,ibtot,icutcoul_eff,ii,ikcalc,ikibz,io,isppol,itypat,jj,method
 integer :: mod10,mqmem,mband,ng_kss,nsheps,ikcalc2bz,ierr,gap_err,ng
 integer :: gwc_nfftot,gwx_nfftot,nqlwl,test_npwkss,my_rank,nprocs,ik,nk_found,ifo,timrev,usefock_ixc
 integer :: iqbz,isym,iq_ibz,itim,ic,pinv,ig1,ng_sigx,spin,gw_qprange,ivcoul_init,nvcoul_init,xclevel_ixc
 real(dp),parameter :: OMEGASIMIN=0.01d0,tol_enediff=0.001_dp*eV_Ha
 real(dp) :: domegas,domegasi,ucvol,rcut
 logical,parameter :: linear_imag_mesh=.TRUE.
 logical :: ltest,remove_inv,changed,found
 character(len=500) :: msg
 character(len=fnlen) :: fname,fcore,string
 type(wvl_internal_type) :: wvl
 type(gaps_t) :: gaps
!arrays
 integer :: ng0sh_opt(3),G0(3),q_umklp(3),kpos(6)
 integer,allocatable :: npwarr(:),val_indeces(:,:)
 integer,pointer :: gvec_kss(:,:),gsphere_sigx_p(:,:)
 integer,pointer :: test_gvec_kss(:,:)
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),sq(3),q_bz(3),gamma_point(3,1)
 real(dp),pointer :: energies_p(:,:,:)
 real(dp),allocatable :: doccde(:),eigen(:),occfact(:),qlwl(:,:)
 type(Pawrhoij_type),allocatable :: Pawrhoij(:)
 type(vcoul_t) :: Vcp_ks

! *************************************************************************

 DBG_ENTER('COLL')

 ! Check for calculations that are not implemented
 ltest=ALL(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol)==Dtset%nband(1))
 ABI_CHECK(ltest,'Dtset%nband(:) must be constant')

 my_rank = xmpi_comm_rank(comm); nprocs  = xmpi_comm_size(comm)

 ! Basic parameters
 Sigp%ppmodel    = Dtset%ppmodel
 Sigp%gwcalctyp  = Dtset%gwcalctyp
 Sigp%nbnds      = Dtset%nband(1)
 Sigp%symsigma   = Dtset%symsigma
 Sigp%zcut       = Dtset%zcut
 Sigp%mbpt_sciss = Dtset%mbpt_sciss

 timrev=  2 ! This information is not reported in the header
            ! 1 => do not use time-reversal symmetry
            ! 2 => take advantage of time-reversal symmetry
 if (any(dtset%kptopt == [3, 4])) timrev = 1

 ! === For HF, SEX or COHSEX use Hybertsen-Louie PPM (only $\omega=0$) ===
 ! * Use fake screening for HF.
 ! FIXME Why, we should not redefine Sigp%ppmodel
 gwcalctyp=Sigp%gwcalctyp
 mod10 =MOD(Sigp%gwcalctyp,10)
 if (mod10==5.or.mod10==6.or.mod10==7) Sigp%ppmodel=2
 if (mod10<5.and.MOD(Sigp%gwcalctyp,1)/=1) then ! * One shot GW (PPM or contour deformation).
   if (Dtset%nomegasrd==1) then ! avoid division by zero!
     Sigp%nomegasrd  =1
     Sigp%maxomega4sd=zero
     Sigp%deltae     =zero
   else
     Sigp%nomegasrd   = Dtset%nomegasrd
     Sigp%maxomega4sd = Dtset%omegasrdmax
     Sigp%deltae     = (2*Sigp%maxomega4sd)/(Sigp%nomegasrd-1)
   endif
 else
   ! For AC no need to evaluate derivative by finite differences.
   Sigp%nomegasrd  =1
   Sigp%maxomega4sd=zero
   Sigp%deltae     =zero
 end if

 ! For analytic continuation define the number of imaginary frequencies for Sigma
 ! Tests show than more than 12 freqs in the Pade approximant worsen the results!
 Sigp%nomegasi=0

 if (mod10==1) then
   Sigp%nomegasi  =Dtset%nomegasi
   Sigp%omegasimax=Dtset%omegasimax
   Sigp%omegasimin=OMEGASIMIN
   write(msg,'(4a,i3,2(2a,f8.3),a)')ch10,&
&    ' Parameters for analytic continuation : ',ch10,&
&    '  number of imaginary frequencies for sigma =  ',Sigp%nomegasi,ch10,&
&    '  min frequency for sigma on imag axis [eV] =  ',Sigp%omegasimin*Ha_eV,ch10,&
&    '  max frequency for sigma on imag axis [eV] =  ',Sigp%omegasimax*Ha_eV,ch10
   call wrtout(std_out,msg,'COLL')

   !TODO this should not be done here but in init_sigma_t
   ABI_MALLOC(Sigp%omegasi,(Sigp%nomegasi))

   if (linear_imag_mesh) then  ! * Linear mesh along the imaginary axis.
     domegasi=Sigp%omegasimax/(Sigp%nomegasi-1)
     do io=1,Sigp%nomegasi
       Sigp%omegasi(io)=CMPLX(zero,(io-1)*domegasi)
     end do
   else ! * Logarithmic mesh along the imaginary axis.
     MSG_ERROR("AC + log mesh not implemented")
     !domegasi=(Sigp%omegasimax/Sigp%omegasimin)**(one/(Sigp%nomegasi-1))
     !Sigp%omegasi(1)=czero; ldi=domegasi
     !do io=2,Sigp%nomegasi
     ! omega(io)=CMPLX(zero,ldi*Sigp%omegasimin)
     ! Sigp%omegasi(io)=ldi*domegasi
     !end do
   end if

   write(msg,'(4a)')ch10,&
&    ' setup_sigma : calculating Sigma(iw)',&
&    ' at imaginary frequencies [eV] (Fermi Level set to 0) ',ch10
   call wrtout(std_out,msg,'COLL')
   call wrtout(ab_out,msg,'COLL')
   do io=1,Sigp%nomegasi
     write(msg,'(2(f10.3,2x))')Sigp%omegasi(io)*Ha_eV
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
   end do

   ltest=(Sigp%omegasimax>0.1d-4.and.Sigp%nomegasi>0)
   ABI_CHECK(ltest,'Wrong value of omegasimax or nomegasi')
   if (Sigp%gwcalctyp/=1) then ! only one shot GW is allowed for AC.
     MSG_ERROR("SC-GW with analytic continuation is not coded")
   end if
 end if

 if (Sigp%symsigma/=0.and.gwcalctyp>=20) then
   MSG_WARNING("SC-GW with symmetries is still under development. Use at your own risk!")
 end if

 ! Setup parameters for Spectral function.
 if (Dtset%gw_customnfreqsp/=0) then
   Sigp%nomegasr = Dtset%gw_customnfreqsp
   MSG_WARNING('Custom grid for spectral function specified. Assuming experienced user.')
   if (Dtset%gw_customnfreqsp/=0) then
     Dtset%nfreqsp = Dtset%gw_customnfreqsp
     MSG_WARNING('nfreqsp has been set to the same number as gw_customnfreqsp')
   end if
 else
   Sigp%nomegasr  =Dtset%nfreqsp
   Sigp%minomega_r=Dtset%freqspmin
   Sigp%maxomega_r=Dtset%freqspmax
 end if

 if (Sigp%nomegasr>0) then
   if (Dtset%gw_customnfreqsp==0) then
     ! Check
     if (Sigp%minomega_r >= Sigp%maxomega_r) then
       MSG_ERROR('freqspmin must be smaller than freqspmax!')
     end if
     if(Sigp%nomegasr==1) then
      domegas=0.d0
     else
      domegas=(Sigp%maxomega_r-Sigp%minomega_r)/(Sigp%nomegasr-1)
     endif
     !TODO this should be moved to Sr% and done in init_sigma_t
     ABI_MALLOC(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Sigp%minomega_r + domegas*(io-1),zero)
     end do
     write(msg,'(4a,i8,3(2a,f8.3),a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  Frequency step [eV] = ',domegas*Ha_eV,ch10
     call wrtout(std_out,msg,'COLL')
   else
     Sigp%minomega_r = MINVAL(Dtset%gw_freqsp(:))
     Sigp%maxomega_r = MAXVAL(Dtset%gw_freqsp(:))
     !TODO this should be moved to Sr% and done in init_sigma_t
     ABI_MALLOC(Sigp%omega_r,(Sigp%nomegasr))
     do io=1,Sigp%nomegasr
       Sigp%omega_r(io) = CMPLX(Dtset%gw_freqsp(io),zero)
     end do
     write(msg,'(4a,i8,2(2a,f8.3),3a)')ch10,&
&      ' Parameters for the calculation of the spectral function : ',ch10,&
&      '  Number of points    = ',Sigp%nomegasr,ch10,&
&      '  Min frequency  [eV] = ',Sigp%minomega_r*Ha_eV,ch10,&
&      '  Max frequency  [eV] = ',Sigp%maxomega_r*Ha_eV,ch10,&
&      '  A custom set of frequencies is used! See the input file for values.',ch10
     call wrtout(std_out,msg,'COLL')
   end if
 else
   !In indefo all these quantities are set to zero
   !Sigp%nomegasr=1
   !allocate(Sigp%omega_r(Sigp%nomegasr))
   !Sigp%omega_r(1)=0
 end if

 ! Dimensional primitive translations rprimd (from input), gprimd, metrics and unit cell volume
 call mkrdim(acell,rprim,rprimd)
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 Sigp%npwwfn=Dtset%npwwfn
 Sigp%npwx  =Dtset%npwsigx

 ! Read parameters of the WFK, verifify them and retrieve all G-vectors.
 call wfk_read_eigenvalues(wfk_fname,energies_p,Hdr_wfk,comm)
 mband = MAXVAL(Hdr_wfk%nband)

 remove_inv = .FALSE.
 call hdr_vs_dtset(Hdr_wfk,Dtset)

 test_npwkss = 0
 call make_gvec_kss(Dtset%nkpt,Dtset%kptns,Hdr_wfk%ecut_eff,Dtset%symmorphi,Dtset%nsym,Dtset%symrel,Dtset%tnons,&
&  gprimd,Dtset%prtvol,test_npwkss,test_gvec_kss,ierr)
 ABI_CHECK(ierr==0,"Fatal error in make_gvec_kss")

 ABI_MALLOC(gvec_kss,(3,test_npwkss))
 gvec_kss = test_gvec_kss
 ng_kss = test_npwkss

 ng = MIN(SIZE(gvec_kss,DIM=2),SIZE(test_gvec_kss,DIM=2))
 ierr = 0
 do ig1=1,ng
   if (ANY(gvec_kss(:,ig1)/=test_gvec_kss(:,ig1))) then
     ierr=ierr+1
     write(std_out,*)" gvec_kss ",ig1,"/",ng,gvec_kss(:,ig1),test_gvec_kss(:,ig1)
   end if
 end do
 ABI_CHECK(ierr == 0, "Mismatch between gvec_kss and test_gvec_kss")

 ABI_FREE(test_gvec_kss)

 ! Get important dimensions from the WFK header
 Sigp%nsppol =Hdr_wfk%nsppol
 Sigp%nspinor=Hdr_wfk%nspinor
 Sigp%nsig_ab=Hdr_wfk%nspinor**2  ! TODO Is it useful calculating only diagonal terms?

 if (Sigp%nbnds>mband) then
   Sigp%nbnds    =mband
   Dtset%nband(:)=mband
   Dtset%mband   =MAXVAL(Dtset%nband)
   write(msg,'(3a,i4,a)')&
&    'Number of bands found less then required',ch10,&
&    'calculation will proceed with nbnds = ',mband,ch10
   MSG_WARNING(msg)
 end if

 ! Check input
 if (Sigp%ppmodel==3.or.Sigp%ppmodel==4) then
   if (gwcalctyp>=10) then
     write(msg,'(a,i3,a)')' The ppmodel chosen and gwcalctyp ',Dtset%gwcalctyp,' are not compatible. '
     MSG_ERROR(msg)
   end if
   if (Sigp%nspinor==2) then
     write(msg,'(a,i3,a)')' The ppmodel chosen and nspinor ',Sigp%nspinor,' are not compatible. '
     MSG_ERROR(msg)
   end if
 end if

 ! Create crystal_t data type
 cryst = hdr_get_crystal(Hdr_wfk, timrev, remove_inv)
 call crystal_print(Cryst)

 if (Sigp%npwwfn>ng_kss) then ! cannot use more G"s for the wfs than those stored on file
   Sigp%npwwfn  =ng_kss
   Dtset%npwwfn =ng_kss
   write(msg,'(2a,(a,i8,a))')&
&    'Number of G-vectors for WFS found in the KSS file is less than required',ch10,&
&    'calculation will proceed with npwwfn  = ',Sigp%npwwfn,ch10
   MSG_WARNING(msg)
 end if

 if (Sigp%npwx>ng_kss) then
   ! Have to recalcuate the (large) sphere for Sigma_x.
   pinv=1; if (remove_inv.and.Cryst%timrev==2) pinv=-1
   gamma_point(:,1) = (/zero,zero,zero/); nullify(gsphere_sigx_p)

   call merge_and_sort_kg(1,gamma_point,Dtset%ecutsigx,Cryst%nsym,pinv,Cryst%symrel,&
&   Cryst%gprimd,gsphere_sigx_p,Dtset%prtvol)

   ng_sigx=SIZE(gsphere_sigx_p,DIM=2)
   Sigp%npwx     = ng_sigx
   Dtset%npwsigx = ng_sigx

   write(msg,'(2a,(a,i8,a))')&
&    'Number of G-vectors for Sigma_x found in the KSS file is less than required',ch10,&
&    'calculation will proceed with npwsigx = ',Sigp%npwx,ch10
   MSG_WARNING(msg)

   ltest = (Sigp%npwx >= ng_kss)
   ABI_CHECK(ltest,"Sigp%npwx<ng_kss!")

   ! * Fill gvec_kss with larger sphere.
   ABI_FREE(gvec_kss)
   ABI_MALLOC(gvec_kss,(3,Sigp%npwx))
   gvec_kss = gsphere_sigx_p
   ABI_FREE(gsphere_sigx_p)
 end if

 ! Set up of the k-points and tables in the whole BZ ===
 ! TODO Recheck symmorphy and inversion
 call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.FALSE.)
 !call kmesh_init(Kmesh,Cryst,Hdr_wfk%nkpt,Hdr_wfk%kptns,Dtset%kptopt,wrap_1zone=.TRUE.)

 ! Some required information are not filled up inside kmesh_init
 ! So doing it here, even though it is not clean
 Kmesh%kptrlatt(:,:) =Dtset%kptrlatt(:,:)
 Kmesh%nshift        =Dtset%nshiftk
 ABI_ALLOCATE(Kmesh%shift,(3,Kmesh%nshift))
 Kmesh%shift(:,:)    =Dtset%shiftk(:,1:Dtset%nshiftk)

 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Kmesh,"K-mesh for the wavefunctions",ab_out, 0,           "COLL")

 ! === Initialize the band structure datatype ===
 ! * Copy WFK energies and occupations up to Sigp%nbnds==Dtset%nband(:)
 bantot=SUM(Dtset%nband(1:Dtset%nkpt*Dtset%nsppol))
 ABI_MALLOC(doccde,(bantot))
 ABI_MALLOC(eigen,(bantot))
 ABI_MALLOC(occfact,(bantot))
 doccde(:)=zero; eigen(:)=zero; occfact(:)=zero

 jj=0; ibtot=0
 do isppol=1,Dtset%nsppol
   do ikibz=1,Dtset%nkpt
     do ib=1,Hdr_wfk%nband(ikibz+(isppol-1)*Dtset%nkpt)
       ibtot=ibtot+1
       if (ib<=Sigp%nbnds) then
         jj=jj+1
         occfact(jj)=Hdr_wfk%occ(ibtot)
         eigen  (jj)=energies_p(ib,ikibz,isppol)
       end if
     end do
   end do
 end do
 ABI_FREE(energies_p)

 ! Make sure that Dtset%wtk==Kmesh%wt due to the dirty treatment of
 ! symmetry operations in the old GW code (symmorphy and inversion)
 ltest=(ALL(ABS(Dtset%wtk(1:Kmesh%nibz)-Kmesh%wt(1:Kmesh%nibz))<tol6))
 ABI_CHECK(ltest,'Mismatch between Dtset%wtk and Kmesh%wt')

 ABI_MALLOC(npwarr,(Dtset%nkpt))
 npwarr(:)=Sigp%npwwfn

 call ebands_init(bantot,KS_BSt,Dtset%nelect,doccde,eigen,Dtset%istwfk,Kmesh%ibz,Dtset%nband,&
&  Kmesh%nibz,npwarr,Dtset%nsppol,Dtset%nspinor,Dtset%tphysel,Dtset%tsmear,Dtset%occopt,occfact,Kmesh%wt,&
&  dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig,&
&  dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)

 ABI_FREE(doccde)
 ABI_FREE(eigen)
 ABI_FREE(npwarr)

 ! Calculate KS occupation numbers and ks_vbk(nkibz,nsppol) ====
 ! ks_vbk gives the (valence|last Fermi band) index for each k and spin.
 ! spinmagntarget is passed to fermi.F90 to fix the problem with newocc in case of magnetic metals
 call ebands_update_occ(KS_BSt,Dtset%spinmagntarget,prtvol=0)

 gap_err = get_gaps(KS_BSt, gaps)
 call gaps%print(unit=std_out)
 call ebands_report_gap(KS_BSt, unit=std_out)

 ABI_MALLOC(val_indeces,(KS_BSt%nkpt,KS_BSt%nsppol))
 val_indeces = get_valence_idx(KS_BSt)

 ! Create Sigma header
 ! TODO Fix problems with symmorphy and k-points
 call hdr_init(KS_BSt,codvsn,Dtset,Hdr_out,Pawtab,pertcase0,Psps,wvl)

 ! Get Pawrhoij from the header of the WFK file
 ABI_DT_MALLOC(Pawrhoij,(Cryst%natom*Dtset%usepaw))
 if (Dtset%usepaw==1) then
   call pawrhoij_alloc(Pawrhoij,1,Dtset%nspden,Dtset%nspinor,Dtset%nsppol,Cryst%typat,pawtab=Pawtab)
   call pawrhoij_copy(Hdr_wfk%Pawrhoij,Pawrhoij)
 end if
 call hdr_update(hdr_out,bantot,1.0d20,1.0d20,1.0d20,Cryst%rprimd,occfact,Pawrhoij,Cryst%xred,dtset%amu_orig(:,1))

 ABI_FREE(occfact)
 call pawrhoij_free(Pawrhoij)
 ABI_DT_FREE(Pawrhoij)

 ! ===========================================================
 ! ==== Setup of k-points and bands for the GW corrections ====
 ! ===========================================================
 ! * maxbdgw and minbdgw are the Max and min band index for GW corrections over k-points.
 !   They are used to dimension the wavefunctions and to calculate the matrix elements.
 !
 if (dtset%nkptgw==0) then
   ! Use qp_range to select the interesting k-points and the corresponing bands.
   !
   !    0 --> Compute the QP corrections only for the fundamental and the direct gap.
   ! +num --> Compute the QP corrections for all the k-points in the irreducible zone and include `num`
   !           bands above and below the Fermi level.
   ! -num --> Compute the QP corrections for all the k-points in the irreducible zone.
   !          Include all occupied states and `num` empty states.

   call wrtout(std_out, "nkptgw == 0 ==> Automatic selection of k-points and bands for the corrections.")
   if (gap_err /=0 .and. dtset%gw_qprange==0) then
     msg = "Problem while computing the fundamental and direct gap (likely metal). Will replace gw_qprange=0 with gw_qprange=1"
     MSG_WARNING(msg)
     dtset%gw_qprange = 1
   end if
   gw_qprange = dtset%gw_qprange

   if (dtset%ucrpa>0) then
     dtset%nkptgw=Kmesh%nbz
     Sigp%nkptgw =dtset%nkptgw
     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
     Sigp%kptgw(:,:)=Kmesh%bz(:,:)
     Sigp%minbnd=1
     Sigp%maxbnd=Sigp%nbnds

   else if (gw_qprange/=0) then
     ! Include all the k-points in the IBZ.
     ! Note that kptgw == ebands%kptns so we can use a single ik index in the loop over k-points
     ! No need to map kptgw onto ebands%kptns.
     dtset%nkptgw=Kmesh%nibz
     Sigp%nkptgw =dtset%nkptgw
     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))
     Sigp%kptgw(:,:)=Kmesh%ibz(:,:)
     Sigp%minbnd=1
     Sigp%maxbnd=Sigp%nbnds

     if (gw_qprange>0) then
       ! All k-points: Add buffer of bands above and below the Fermi level.
       do spin=1,Sigp%nsppol
         do ik=1,Sigp%nkptgw
           Sigp%minbnd(ik,spin) = MAX(val_indeces(ik,spin) - gw_qprange, 1)
           Sigp%maxbnd(ik,spin) = MIN(val_indeces(ik,spin) + gw_qprange + 1, Sigp%nbnds)
         end do
       end do

     else
       ! All k-points: include all occupied states and -gw_qprange empty states.
       Sigp%minbnd = 1
       do spin=1,Sigp%nsppol
         do ik=1,Sigp%nkptgw
           Sigp%maxbnd(ik,spin) = MIN(val_indeces(ik,spin) - gw_qprange, Sigp%nbnds)
         end do
       end do
     end if

   else
     ! gw_qprange is not specified in the input.
     ! Include the direct and the fundamental KS gap.
     ! The main problem here is that kptgw and nkptgw do not depend on the spin and therefore
     ! we have compute the union of the k-points where the fundamental and the direct gaps are located.
     !
     ! Find the list of `interesting` kpoints.
     ABI_CHECK(gap_err == 0, "gw_qprange 0 cannot be used because I cannot find the gap (gap_err !=0)")
     nk_found = 1; kpos(1) = gaps%fo_kpos(1,1)

     do spin=1,Sigp%nsppol
       do ifo=1,3
         ik = gaps%fo_kpos(ifo, spin)
         found = .FALSE.; jj = 0
         do while (.not. found .and. jj < nk_found)
           jj = jj + 1; found = (kpos(jj) == ik)
         end do
         if (.not. found) then
           nk_found = nk_found + 1; kpos(nk_found) = ik
         end if
       end do
     end do

     ! Now we can define the list of k-points and the bands range.
     dtset%nkptgw=nk_found
     Sigp%nkptgw =dtset%nkptgw

     ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
     ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
     ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))

     do ii=1,Sigp%nkptgw
       ik = kpos(ii)
       Sigp%kptgw(:,ii)=Kmesh%ibz(:,ik)
       do spin=1,Sigp%nsppol
         Sigp%minbnd(ii,spin) = val_indeces(ik,spin)
         Sigp%maxbnd(ii,spin) = val_indeces(ik,spin) + 1
       end do
     end do
   end if

 else
   ! Treat only the k-points and bands specified in the input file.
   Sigp%nkptgw=dtset%nkptgw
   ABI_MALLOC(Sigp%kptgw,(3,Sigp%nkptgw))
   ABI_MALLOC(Sigp%minbnd,(Sigp%nkptgw,Sigp%nsppol))
   ABI_MALLOC(Sigp%maxbnd,(Sigp%nkptgw,Sigp%nsppol))

   do spin=1,Sigp%nsppol
     Sigp%minbnd(:,spin)=dtset%bdgw(1,:,spin)
     Sigp%maxbnd(:,spin)=dtset%bdgw(2,:,spin)
   end do

   do ii=1,3
     do ikcalc=1,Sigp%nkptgw
       Sigp%kptgw(ii,ikcalc)=Dtset%kptgw(ii,ikcalc)
     end do
   end do

   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
       if (Dtset%bdgw(2,ikcalc,spin)>Sigp%nbnds) then
         write(msg,'(a,2i0,2(a,i0),2a,i0)')&
&          "For (k,s) ",ikcalc,spin," bdgw= ",Dtset%bdgw(2,ikcalc,spin), " > nbnds=",Sigp%nbnds,ch10,&
&          "Calculation will continue with bdgw =",Sigp%nbnds
         MSG_COMMENT(msg)
         Dtset%bdgw(2,ikcalc,spin)=Sigp%nbnds
       end if
     end do
   end do

 end if

 ! Make sure that all the degenerate states are included.
 ! * We will have to average the GW corrections over degenerate states if symsigma=1 is used.
 ! * KS states belonging to the same irreducible representation should be included in the basis set used for SCGW.
 if (Sigp%symsigma/=0 .or. gwcalctyp>=10) then
   do isppol=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw

       if (has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikibz,G0)) then
         call enclose_degbands(KS_BSt,ikibz,isppol,Sigp%minbnd(ikcalc,isppol),Sigp%maxbnd(ikcalc,isppol),changed,tol_enediff)
         if (changed) then
           write(msg,'(2(a,i0),2a,2(1x,i0))')&
&            "Not all the degenerate states at ikcalc= ",ikcalc,", spin= ",isppol,ch10,&
&            "were included in the bdgw set. bdgw has been changed to: ",Sigp%minbnd(ikcalc,isppol),Sigp%maxbnd(ikcalc,isppol)
           MSG_COMMENT(msg)
         end if
       else
         MSG_ERROR(sjoin('k-point', ktoa(Sigp%kptgw(:,ikcalc)), 'not in IBZ'))
       end if

     end do
   end do
 end if

 !if (.not. associated(Dtset%bdgw)) then
 !  ABI_MALLOC(Dtset%bdgw, (2,Sigp%nkptgw,Sigp%nsppol))
 !end if
 !do spin=1,Sigp%nsppol
 !  Dtset%bdgw(1,:,spin) = Sigp%minbnd(:,spin)
 !  Dtset%bdgw(2,:,spin) = Sigp%maxbnd(:,spin)
 !end do

 Sigp%minbdgw=MINVAL(Sigp%minbnd)
 Sigp%maxbdgw=MAXVAL(Sigp%maxbnd)

 ABI_MALLOC(Sigp%kptgw2bz,(Sigp%nkptgw))
 !
 !=== Check if the k-points are in the BZ ===
 !FB TODO Honestly the code is not able to treat k-points, which are not in the IBZ.
 !This extension should require to change the code in different places.
 !Therefore, one should by now prevent the user from calculating sigma for a k-point not in the IBZ.

 do ikcalc=1,Sigp%nkptgw
   if (has_BZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)) then
     !found = has_IBZ_item(Kmesh,Sigp%kptgw(:,ikcalc),ikcalc2bz,G0)
     Sigp%kptgw2bz(ikcalc) = ikcalc2bz
   else
     MSG_ERROR(sjoin('k-point:', ktoa(Sigp%kptgw(:,ikcalc)), 'not in the kbz set'))
   end if
 end do

 ! Check if there are duplicated k-point in Sigp%
 do ii=1,Sigp%nkptgw
   do jj=ii+1,Sigp%nkptgw
     if (isamek(Sigp%kptgw(:,ii),Sigp%kptgw(:,jj),G0)) then
       write(msg,'(5a)')&
&        'kptgw contains duplicated k-points. This is not allowed since ',ch10,&
&        'the QP corrections for this k-point will be calculated more than once. ',ch10,&
&        'Check your input file. '
       MSG_ERROR(msg)
     end if
   end do
 end do
 !
 ! Warn the user if SCGW run and not all the k-points are included.
 if (gwcalctyp>=10 .and. Sigp%nkptgw/=Hdr_wfk%nkpt) then
   write(msg,'(3a,2(a,i0),2a)')ch10,&
&    " COMMENT: In a self-consistent GW run, the QP corrections should be calculated for all the k-points of the KSS file ",ch10,&
&    " but nkptgw= ",Sigp%nkptgw," and WFK nkpt= ",Hdr_wfk%nkpt,ch10,&
&    " Assuming expert user. Execution will continue. "
   call wrtout(ab_out,msg,"COLL")
 end if

 ! Setup of the table used in the case of SCGW on wavefunctions to reduce the number
 ! of elements <i,kgw,s|\Sigma|j,kgw,s> that have to be calculated. No use of symmetries, except for Hermiticity.
 call sigma_tables(Sigp,Kmesh)

 ! === Read external file and initialize basic dimension of Er% ===
 ! TODO use mqmem as input variable instead of gwmem

 ! === If required, use a matrix for $\Sigma_c$ which is smaller than that stored on file ===
 ! * By default the entire matrix is read and used,
 ! * Define consistently npweps and ecuteps for \Sigma_c according the input
 if (Dtset%npweps>0.or.Dtset%ecuteps>0) then
   ! This should not happen: the Dtset array should not be modified after having been initialized.
   if (Dtset%npweps>0) Dtset%ecuteps=zero
   nsheps=0
   call setshells(Dtset%ecuteps,Dtset%npweps,nsheps,Dtset%nsym,gmet,gprimd,Dtset%symrel,'eps',ucvol)
 end if

 mqmem=0; if (Dtset%gwmem/10==1) mqmem=1

 if (Dtset%getscr/=0.or.Dtset%irdscr/=0) then
   fname=Dtfil%fnameabi_scr
 else if (Dtset%getsuscep/=0.or.Dtset%irdsuscep/=0) then
   fname=Dtfil%fnameabi_sus
 else
   fname=Dtfil%fnameabi_scr
   !FIXME this has to be cleaned, in tgw2_3 Dtset%get* and Dtset%ird* are  not defined
   !MSG_ERROR("getsuscep or irdsuscep are not defined")
 end if
 !
 ! === Setup of q-mesh in the whole BZ ===
 ! * Stop if a nonzero umklapp is needed to reconstruct the BZ. In this case, indeed,
 !   epsilon^-1(Sq) should be symmetrized in csigme using a different expression (G-G_o is needed)
 !
 if (sigma_needs_w(Sigp)) then
   if (.not. file_exists(fname)) then
     fname = nctk_ncify(fname)
     MSG_COMMENT(sjoin("File not found. Will try netcdf file:", fname))
   end if

   call init_Er_from_file(Er,fname,mqmem,Dtset%npweps,comm)

   Sigp%npwc=Er%npwe
   if (Sigp%npwc>Sigp%npwx) then
     Sigp%npwc=Sigp%npwx
     MSG_COMMENT("Found npw_correlation > npw_exchange, Imposing npwc=npwx")
     ! There is a good reason for doing so, see csigme.F90 and the size of the arrays
     ! rhotwgp and rhotwgp: we need to define a max size and we opt for Sigp%npwx.
   end if
   Er%npwe=Sigp%npwc
   Dtset%npweps=Er%npwe
   call kmesh_init(Qmesh,Cryst,Er%nqibz,Er%qibz,Dtset%kptopt)

 else
   Er%npwe     =1
   Sigp%npwc   =1
   Dtset%npweps=1
   call find_qmesh(Qmesh,Cryst,Kmesh)
   ABI_MALLOC(Er%gvec,(3,1))
   Er%gvec(:,1) = (/0,0,0/)
 end if

 call kmesh_print(Qmesh,"Q-mesh for screening function",std_out,Dtset%prtvol,"COLL")
 call kmesh_print(Qmesh,"Q-mesh for screening function",ab_out ,0           ,"COLL")

 do iqbz=1,Qmesh%nbz
   call get_BZ_item(Qmesh,iqbz,q_bz,iq_ibz,isym,itim,umklp=q_umklp)

   if (ANY(q_umklp/=0)) then
     sq = (3-2*itim)*MATMUL(Cryst%symrec(:,:,isym),Qmesh%ibz(:,iq_ibz))
     write(std_out,*) sq,Qmesh%bz(:,iqbz)
     write(msg,'(a,3f6.3,a,3f6.3,2a,9i3,a,i2,2a)')&
&      'qpoint ',Qmesh%bz(:,iqbz),' is the symmetric of ',Qmesh%ibz(:,iq_ibz),ch10,&
&      'through operation ',Cryst%symrec(:,:,isym),' and itim ',itim,ch10,&
&      'however a non zero umklapp G_o vector is required and this is not yet allowed'
     MSG_ERROR(msg)
   end if
 end do
 !
 ! === Find optimal value for G-sphere enlargment due to oscillator matrix elements ===
 ! * Here I have to be sure that Qmesh%bz is always inside the BZ, not always true size bz is buggy
 ! * -one is used because we loop over all the possibile differences, unlike screening

 call get_ng0sh(Sigp%nkptgw,Sigp%kptgw,Kmesh%nbz,Kmesh%bz,Qmesh%nbz,Qmesh%bz,-one,ng0sh_opt)
 call wrtout(std_out, sjoin(' Optimal value for ng0sh ', ltoa(ng0sh_opt)), "COLL")
 Sigp%mG0=ng0sh_opt

! G-sphere for W and Sigma_c is initialized from the SCR file.
 call gsph_init(Gsph_c,Cryst,Er%npwe,gvec=Er%gvec)
 call gsph_init(Gsph_x,Cryst,Sigp%npwx,gvec=gvec_kss)
 Sigp%ecuteps = Gsph_c%ecut
 Dtset%ecuteps = Sigp%ecuteps
 
! === Make biggest G-sphere of Sigp%npwvec vectors ===
 Sigp%npwvec=MAX(Sigp%npwwfn,Sigp%npwx)
 call gsph_init(Gsph_Max,Cryst,Sigp%npwvec,gvec=gvec_kss)
!BEGINDEBUG
 ! Make sure that the two G-spheres are equivalent.
 ierr=0
 if (sigma_needs_w(Sigp)) then
   ng = MIN(SIZE(Gsph_c%gvec,DIM=2),SIZE(gvec_kss,DIM=2))
   do ig1=1,ng
     if (ANY(Gsph_c%gvec(:,ig1)/=gvec_kss(:,ig1))) then
       ierr=ierr+1
       write(std_out,*)" Gsph_c, gvec_kss ",ig1,"/",ng,Gsph_c%gvec(:,ig1),gvec_kss(:,ig1)
     end if
   end do
   ABI_CHECK(ierr==0,"Mismatch between Gsph_c and gvec_kss")
 end if
 ierr=0
 ng = MIN(SIZE(Gsph_x%gvec,DIM=2),SIZE(gvec_kss,DIM=2))
 do ig1=1,ng
   if (ANY(Gsph_x%gvec(:,ig1)/=gvec_kss(:,ig1))) then
     ierr=ierr+1
     write(std_out,*)" Gsph_x, gvec_kss ",ig1,"/",ng,Gsph_x%gvec(:,ig1),gvec_kss(:,ig1)
   end if
 end do
 ABI_CHECK(ierr==0,"Mismatch between Gsph_x and gvec_kss")
!ENDDEBUG

 ABI_FREE(gvec_kss)
 !
 ! === Get Fourier components of the Coulombian for all q-points in the IBZ ===
 ! * If required, use a cutoff in the interaction
 ! * Pcv%vc_sqrt contains Vc^{-1/2}
 ! * Setup also the analytical calculation of the q->0 component
 ! FIXME recheck ngfftf since I got different charge outside the cutoff region

 if (Dtset%gw_nqlwl==0) then
   nqlwl=1
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,1)= GW_Q0_DEFAULT
 else
   nqlwl=Dtset%gw_nqlwl
   ABI_MALLOC(qlwl,(3,nqlwl))
   qlwl(:,:)=Dtset%gw_qlwl(:,1:nqlwl)
 end if

!The Coulomb interaction used here might have two terms :
!the first term generates the usual sigma self-energy, but possibly, one should subtract
!from it the Coulomb interaction already present in the Kohn-Sham basis,
!if the usefock associated to ixc is one.
!The latter excludes (in the present implementation) mod(Dtset%gwcalctyp,10)==5
 nvcoul_init=1
 call get_xclevel(Dtset%ixc,xclevel_ixc,usefock_ixc)
 if(usefock_ixc==1)then
   nvcoul_init=2
   if(mod(Dtset%gwcalctyp,10)==5)then
     write(msg,'(4a,i3,a,i3,4a,i5)')ch10,&
&     ' The starting wavefunctions were obtained from self-consistent calculations in the planewave basis set',ch10,&
&     ' with ixc = ',Dtset%ixc,' associated with usefock =',usefock_ixc,ch10,&
&     ' In this case, the present implementation does not allow that the self-energy for sigma corresponds to',ch10,&
&     '  mod(gwcalctyp,10)==5, while your gwcalctyp= ',Dtset%gwcalctyp
     MSG_ERROR(msg)
   endif
 endif

 do ivcoul_init=1,nvcoul_init
   rcut = Dtset%rcut
   icutcoul_eff=Dtset%icutcoul
   Sigp%sigma_mixing=one
   if( mod(Dtset%gwcalctyp,10)==5 .or. ivcoul_init==2)then
     if(abs(Dtset%hyb_mixing)>tol8)then
!      Warning : the absolute value is needed, because of the singular way used to define the default for this input variable
       Sigp%sigma_mixing=abs(Dtset%hyb_mixing)
     else if(abs(Dtset%hyb_mixing_sr)>tol8)then
       Sigp%sigma_mixing=abs(Dtset%hyb_mixing_sr)
       icutcoul_eff=5
     endif
     if(abs(rcut)<tol6 .and. abs(Dtset%hyb_range_fock)>tol8)rcut=one/Dtset%hyb_range_fock
   endif

!#if 1
   if(ivcoul_init==1)then

     if (Gsph_x%ng > Gsph_c%ng) then
       call vcoul_init(Vcp,Gsph_x,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
     else
       call vcoul_init(Vcp,Gsph_c,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
     end if

   else

!    Use a temporary Vcp_ks to compute the Coulomb interaction already present in the Fock part of the Kohn-Sham Hamiltonian
     if (Gsph_x%ng > Gsph_c%ng) then
       call vcoul_init(Vcp_ks,Gsph_x,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_x%ng,nqlwl,qlwl,ngfftf,comm)
     else
       call vcoul_init(Vcp_ks,Gsph_c,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,Dtset%vcutgeo,&
&        Dtset%ecutsigx,Gsph_c%ng,nqlwl,qlwl,ngfftf,comm)
     end if

!    Now compute the residual Coulomb interaction
     Vcp%vc_sqrt_resid=sqrt(Vcp%vc_sqrt**2-Sigp%sigma_mixing*Vcp_ks%vc_sqrt**2)
     Vcp%i_sz_resid=Vcp%i_sz-Sigp%sigma_mixing*Vcp_ks%i_sz
!    The mixing factor has already been accounted for, so set it back to one
     Sigp%sigma_mixing=one
     call vcoul_free(Vcp_ks)

!DEBUG
     write(std_out,'(a)')' setup_sigma : the residual Coulomb interaction has been computed'
!ENDDEBUG

   endif
!#else
!   call vcoul_init(Vcp,Gsph_Max,Cryst,Qmesh,Kmesh,rcut,icutcoul_eff,ivcoul_init,Dtset%vcutgeo,&
!&    Dtset%ecutsigx,Sigp%npwx,nqlwl,qlwl,ngfftf,comm)
!#endif

 enddo

#if 0
 ! Using the random q for the optical limit is one of the reasons
 ! why sigma breaks the initial energy degeneracies.
 Vcp%i_sz=zero
 Vcp%vc_sqrt(1,:)=czero
 Vcp%vcqlwl_sqrt(1,:)=czero
#endif

 ABI_FREE(qlwl)

 Sigp%ecuteps = Dtset%ecuteps
 Sigp%ecutwfn = Dtset%ecutwfn
 Sigp%ecutsigx = Dtset%ecutsigx

 ! === Setup of the FFT mesh for the oscilator strengths ===
 ! * gwc_ngfft(7:18)==Dtset%ngfft(7:18) which is initialized before entering screening.
 ! * Here we redefine gwc_ngfft(1:6) according to the following options :
 !
 ! method==0 --> FFT grid read from fft.in (debugging purpose)
 ! method==1 --> Normal FFT mesh
 ! method==2 --> Slightly augmented FFT grid to calculate exactly rho_tw_g (see setmesh.F90)
 ! method==3 --> Doubled FFT grid, same as the the FFT for the density,
 !
 ! enforce_sym==1 ==> Enforce a FFT mesh compatible with all the symmetry operation and FFT library
 ! enforce_sym==0 ==> Find the smallest FFT grid compatbile with the library, do not care about symmetries
 !
 gwc_ngfft(1:18)=Dtset%ngfft(1:18)
 gwx_ngfft(1:18)=Dtset%ngfft(1:18)

 method=2
 if (Dtset%fftgw==00 .or. Dtset%fftgw==01) method=0
 if (Dtset%fftgw==10 .or. Dtset%fftgw==11) method=1
 if (Dtset%fftgw==20 .or. Dtset%fftgw==21) method=2
 if (Dtset%fftgw==30 .or. Dtset%fftgw==31) method=3
 enforce_sym=MOD(Dtset%fftgw,10)

 ! FFT mesh for sigma_x.
 call setmesh(gmet,Gsph_Max%gvec,gwx_ngfft,Sigp%npwvec,Sigp%npwx,Sigp%npwwfn,&
&  gwx_nfftot,method,Sigp%mG0,Cryst,enforce_sym)

 ! FFT mesh for sigma_c.
 call setmesh(gmet,Gsph_Max%gvec,gwc_ngfft,Sigp%npwvec,Er%npwe,Sigp%npwwfn,&
&  gwc_nfftot,method,Sigp%mG0,Cryst,enforce_sym,unit=dev_null)

 !call new_setmesh(Cryst,ecut_osc,ecutwfn,nkpt,kpoints,method,Sigp%mG0,enforce_sym,gwx_ngfft,gwx_nfftot)
 !call new_setmesh(Cryst,ecut_osc,ecutwfn,nkpt,kpoints,method,Sigp%mG0,enforce_sym,gwc_ngfft,gwc_nfftot)

 ! ======================================================================
 ! ==== Check for presence of files with core orbitals, for PAW only ====
 ! ======================================================================
 Sigp%use_sigxcore=0
 if (Dtset%usepaw==1.and.Dtset%gw_sigxcore==1) then
   ii = 0
   do itypat=1,Cryst%ntypat
     string = Psps%filpsp(itypat)
     fcore = "CORE_"//TRIM(basename(string))
     ic = INDEX (TRIM(string), "/" , back=.TRUE.)
     if (ic>0 .and. ic<LEN_TRIM(string)) then
       ! string defines a path, prepend path to fcore
       fcore = Psps%filpsp(itypat)(1:ic)//TRIM(fcore)
     end if
     if (file_exists(fcore)) then
       ii = ii+1
     else
       MSG_WARNING(sjoin("HF decoupling is required but cannot find file:", fcore))
     end if
   end do

   Sigp%use_sigxcore=1
   if (ii/=Cryst%ntypat) then
     MSG_ERROR("Files with core orbitals not found")
   end if
 end if ! PAW+HF decoupling
 !
 ! ==============================
 ! ==== Extrapolar technique ====
 ! ==============================
 Sigp%gwcomp   = Dtset%gwcomp
 Sigp%gwencomp = Dtset%gwencomp

 if (Sigp%gwcomp==1) then
   write(msg,'(6a,e11.4,a)')ch10,&
&    'Using the extrapolar approximation to accelerate convergence',ch10,&
&    'with respect to the number of bands included',ch10,&
&    'with gwencomp: ',Sigp%gwencomp*Ha_eV,' [eV]'
   call wrtout(std_out,msg,'COLL')
 end if
 !
 ! ===================================
 ! ==== Final compatibility tests ====
 ! ===================================
 if (ANY(KS_BSt%istwfk/=1)) then
   MSG_WARNING('istwfk/=1 is still under development')
 end if

 ltest=(KS_BSt%mband==Sigp%nbnds.and.ALL(KS_BSt%nband==Sigp%nbnds))
 ABI_CHECK(ltest,'BUG in definition of KS_BSt%nband')

 ! FIXME
 if (Dtset%symsigma/=0 .and. Sigp%nomegasr/=0) then
   if (cryst%idx_spatial_inversion() == 0) then
     write(msg,'(5a)')' setup_sigma : BUG :',ch10,&
&      'It is not possible to use symsigma/=0 to calculate the spectral function ',ch10,&
&      'when the system does not have the spatial inversion. Please use symsigma=0 '
     MSG_WARNING(msg)
   end if
 end if

 if (mod10==SIG_GW_AC) then
   if (Sigp%gwcalctyp/=1) MSG_ERROR("Self-consistency with AC not implemented")
   if (Sigp%gwcomp==1) MSG_ERROR("AC with extrapolar technique not implemented")
 end if

 call gaps%free()

 ABI_FREE(val_indeces)

 DBG_EXIT('COLL')

end subroutine setup_sigma
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_tables
!! NAME
!! sigma_tables
!!
!! FUNCTION
!!  Build symmetry tables used to speedup self-consistent GW calculations.
!!
!! INPUTS
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! [Bnd_sym(Kmesh%nibz,Sigp%nsppol)] <type(Bands_Symmetries)>
!!
!! SiDE EFFECTS
!! Sigp<sigparams_t>=This routine initializes the tables:
!!   %Sigcij_tab
!!   %Sigxij_tab
!!  that are used to select the matrix elements of the self-energy that have to be calculated.
!!
!! PARENTS
!!      setup_sigma,sigma
!!
!! CHILDREN
!!
!! SOURCE


subroutine sigma_tables(Sigp,Kmesh,Bnd_sym)

!Arguments ------------------------------------
!scalars
 type(sigparams_t),intent(inout) :: Sigp
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 type(esymm_t),optional,intent(in) :: Bnd_sym(Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: gwcalctyp,spin,ikcalc,ik_ibz,bmin,bmax,bcol,brow
 integer :: ii,idx_x,idx_c,irr_idx1,irr_idx2
!arrays
 integer,allocatable :: sigc_bidx(:),sigx_bidx(:)
 logical :: use_sym_at(Kmesh%nibz,Sigp%nsppol)

! *************************************************************************

 gwcalctyp=Sigp%gwcalctyp

 ! Recreate the Sig_ij tables taking advantage of the classification of the bands.
 if (allocated(Sigp%Sigxij_tab)) then
   call sigijtab_free(Sigp%Sigxij_tab)
   ABI_DT_FREE(Sigp%Sigxij_tab)
 end if
 if (allocated(Sigp%Sigcij_tab)) then
   call sigijtab_free(Sigp%Sigcij_tab)
   ABI_DT_FREE(Sigp%Sigcij_tab)
 end if

 ABI_DT_MALLOC(Sigp%Sigcij_tab,(Sigp%nkptgw,Sigp%nsppol))
 ABI_DT_MALLOC(Sigp%Sigxij_tab,(Sigp%nkptgw,Sigp%nsppol))

 use_sym_at=.FALSE.
 if (PRESENT(Bnd_sym)) then
   do spin=1,Sigp%nsppol
     do ikcalc=1,Sigp%nkptgw
      ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))
      use_sym_at(ik_ibz,spin) = ( .not.esymm_failed(Bnd_sym(ik_ibz,spin)) )
     end do
   end do
 end if

 do spin=1,Sigp%nsppol
   do ikcalc=1,Sigp%nkptgw
     ik_ibz = Kmesh%tab(Sigp%kptgw2bz(ikcalc))

     if (use_sym_at(ik_ibz,spin)) then
       if (gwcalctyp<20) then
         MSG_ERROR("You should not be here!")
       end if

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_DT_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_DT_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))

       do bcol=bmin,bmax
         ABI_MALLOC(sigc_bidx,(bmax-bmin+1))
         ABI_MALLOC(sigx_bidx,(bmax-bmin+1))

         if (Bnd_sym(ik_ibz,spin)%err_status/=0) then   ! Band classification failed.
           sigc_bidx = (/(ii,ii=bmin,bmax)/)
           idx_c = bmax-bmin+1
           sigx_bidx = (/(ii,ii=bmin,bcol)/) ! Hermitian
           idx_x = bcol-bmin+1
         else
           irr_idx2 = Bnd_sym(ik_ibz,spin)%b2irrep(bcol)
           idx_c = 0
           do brow=bmin,bmax
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow)
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle for HF, SEX, or COHSEX.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_c = idx_c +1
               sigc_bidx(idx_c) = brow
             end if
           end do
           idx_x = 0
           do brow=bmin,bcol
             irr_idx1 = Bnd_sym(ik_ibz,spin)%b2irrep(brow)
             if (bcol<brow) CYCLE  ! Sig_x is always Hermitian.
             if (irr_idx1 == irr_idx2) then ! same character, add this row to the list.
               idx_x = idx_x +1
               sigx_bidx(idx_x) = brow
             end if
           end do
         end if
         !
         ! Table for Sigma_x matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_x))

         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= idx_x
         Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigx_bidx(1:idx_x)
         !write(std_out,*)" Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_x,(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(ii),ii=1,idx_x)
         !
         ! Table for Sigma_c matrix elements taking into account symmetries of the bands.
         ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))

         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c
         Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
         !write(std_out,*)" Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol
         !write(std_out,*)" size: ",idx_c,(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(ii), ii=1,idx_c)

         ABI_FREE(sigx_bidx)
         ABI_FREE(sigc_bidx)
       end do ! bcol

     else  ! Symmetries cannot be used for this (k,s).

       bmin=Sigp%minbnd(ikcalc,spin); bmax=Sigp%maxbnd(ikcalc,spin)
       ABI_DT_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col,(bmin:bmax))
       ABI_DT_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col,(bmin:bmax))

       if (gwcalctyp<20) then  ! QP wavefunctions == KS, therefore only diagonal elements are calculated.
         do bcol=bmin,bmax
           ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
           ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(1:1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= 1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(1) = bcol
         end do
       else
         ! Use QP wavefunctions, Sigma_ij matrix is sparse but we have to classify the states in sigma.
         ! The only thing we can do here is filling the entire matrix taking advantage of Hermiticity (if any).
         do bcol=bmin,bmax
           ABI_MALLOC(Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx,(bcol-bmin+1))
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%size1= bcol-bmin+1
           Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:) = (/(ii,ii=bmin,bcol)/) ! Sigma_x is Hermitian.
           !write(std_out,*)"Sigxij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigxij_tab(ikcalc,spin)%col(bcol)%bidx(:)

           ABI_MALLOC(sigc_bidx,(bmax-bmin+1))
           idx_c = 0
           do brow=bmin,bmax
             if (sigma_is_herm(Sigp).and.bcol<brow) CYCLE  ! Only the upper triangle of Sigc_ij is needed (SEX, COHSEX).
             idx_c = idx_c +1
             sigc_bidx(idx_c) = brow
           end do
           ABI_MALLOC(Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx,(idx_c))
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%size1= idx_c
           Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:) = sigc_bidx(1:idx_c)
           ABI_FREE(sigc_bidx)
           !write(std_out,*)"Sigcij_tab: ikcalc, spin, bcol ",ikcalc,spin,bcol,Sigp%Sigcij_tab(ikcalc,spin)%col(bcol)%bidx(:)
         end do
       end if
     end if

   end do !ikcalc
 end do !spin

end subroutine sigma_tables
!!***

!----------------------------------------------------------------------

!!****f* ABINIT/sigma_bksmask
!! NAME
!! sigma_bksmask
!!
!! FUNCTION
!!  Compute tables for the distribution and the storage of the wavefunctions in the SIGMA code.
!!
!! INPUTS
!! Dtset<type(dataset_type)>=all input variables for this dataset
!! Sigp<sigparams_t>=Parameters governing the self-energy calculation.
!! Kmesh <kmesh_t>=Structure describing the k-point sampling.
!! nprocs=Total number of MPI processors
!! my_rank=Rank of this this processor.
!!
!! OUTPUT
!! my_spins(:)=Pointer to NULL in input. In output: list of spins treated by this node.
!! bks_mask(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will treat this state.
!! keep_ur(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)=True if this node will store this state in real space.
!! ierr=Exit status.
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!
!! SOURCE

subroutine sigma_bksmask(Dtset,Sigp,Kmesh,my_rank,nprocs,my_spins,bks_mask,keep_ur,ierr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_rank,nprocs
 integer,intent(out) :: ierr
 type(Dataset_type),intent(in) :: Dtset
 type(sigparams_t),intent(in) :: Sigp
 type(kmesh_t),intent(in) :: Kmesh
!arrays
 integer,allocatable,intent(out) :: my_spins(:)
 logical,intent(out) :: bks_mask(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)
 logical,intent(out) :: keep_ur(Sigp%nbnds,Kmesh%nibz,Sigp%nsppol)

!Local variables-------------------------------
!scalars
 integer :: my_nspins,my_maxb,my_minb,isp,spin,band,nsppol,rank_spin
 logical :: store_ur
!arrays
 integer :: tmp_spins(Sigp%nsppol),nprocs_spin(Sigp%nsppol)

! *************************************************************************

 ierr=0; nsppol=Sigp%nsppol

 ! List of spins for each node, number of processors per each spin
 ! and the MPI rank in the "spin" communicator.
 my_nspins=nsppol
 ABI_MALLOC(my_spins, (nsppol))
 my_spins= [(isp, isp=1,nsppol)]
 nprocs_spin = nprocs; rank_spin = my_rank

 if (nsppol==2 .and. nprocs>1) then
   ! Distribute spins (optimal distribution if nprocs is even)
   nprocs_spin(1) = nprocs/2
   nprocs_spin(2) = nprocs - nprocs/2
   my_nspins=1
   my_spins(1)=1
   if (my_rank+1>nprocs/2) then
     ! I will treate spin=2, compute shifted rank.
     my_spins(1)=2
     rank_spin = my_rank - nprocs/2
   end if
 end if

 store_ur = (MODULO(Dtset%gwmem,10)==1)
 keep_ur=.FALSE.; bks_mask=.FALSE.

 select case (Dtset%gwpara)
 case (1)
   ! Parallelization over transitions **without** memory distributions (Except for the spin).
   my_minb=1; my_maxb=Sigp%nbnds
   if (dtset%ucrpa>0) then
     bks_mask(my_minb:my_maxb,:,:)=.TRUE.
     if (store_ur) keep_ur(my_minb:my_maxb,:,:)=.TRUE.
   else
     do isp=1,my_nspins
       spin = my_spins(isp)
       bks_mask(my_minb:my_maxb,:,spin)=.TRUE.
       if (store_ur) keep_ur(my_minb:my_maxb,:,spin)=.TRUE.
     end do
   end if

 case (2)
   ! Distribute bands and spin (alternating planes of bands)
   do isp=1,my_nspins
     spin = my_spins(isp)
     do band=1,Sigp%nbnds
       if (xmpi_distrib_with_replicas(band,Sigp%nbnds,rank_spin,nprocs_spin(spin))) then
       !if (MODULO(band-1,nprocs_spin(spin))==rank_spin) then
         bks_mask(band,:,spin)=.TRUE.
         if (store_ur) keep_ur(band,:,spin)=.TRUE.
       end if
     end do
   end do

#if 0
   ! Each node store the full set of occupied states to speed up Sigma_x.
   do isp=1,my_nspins
     spin = my_spins(isp)
     do ik_ibz=1,Kmesh%nibz
       ks_iv=ks_vbik(ik_ibz,spin) ! Valence index for this (k,s)
       bks_mask(1:ks_iv,:,spin)=.TRUE.
       if (store_ur) keep_ur(1:ks_iv,:,spin)=.TRUE.
     end do
   end do
#endif

 case default
   MSG_WARNING("Wrong value for gwpara")
   ierr = 1
 end select

 ! Return my_spins with correct size.
 tmp_spins(1:my_nspins) = my_spins(1:my_nspins)

 ABI_FREE(my_spins)
 ABI_MALLOC(my_spins, (my_nspins))
 my_spins = tmp_spins(1:my_nspins)

end subroutine sigma_bksmask
!!***

!!****f* ABINIT/paw_qpscgw
!! NAME
!! paw_qpscgw
!!
!! FUNCTION
!!  This routine is called during QP self-consistent GW calculations. It calculates the new QP on-site quantities
!!  using the QP amplitudes read from the QPS file.
!!
!! INPUTS
!!  Wfd<wfd_t>=Datatype gathering data on QP amplitudes.
!!  nscf=Number of QPSCGW iterations done so far (read from the QPS file).
!!  nfftf=Number of points in the fine FFT grid.
!!  ngfft(18)=information on the fine FFT grid used for densities and potentials.
!!  Dtset<dataset_type>=All input variables for this dataset.
!!  Cryst<crystal_t>=Info on unit cell and symmetries.
!!  Kmesh<kmesh_t>=Structure describing the k-point sampling.
!!  Psps<Pseudopotential_type)>=Info on pseudopotential, only for consistency check of the KSS file
!!  Pawang<pawang_type>=PAW angular mesh and related data.
!!  Pawrad(ntypat*usepaw)<type(pawrad_type)>=paw radial mesh and related data
!!  Pawtab(ntypat*usepaw)<type(pawtab_type)>=paw tabulated starting data
!!  Pawfgrtab(natom)<Pawfgrtab_type>= For PAW, various arrays giving data related to fine grid for a given atom.
!!  prev_Pawrhoij(Cryst%natom))<Pawrhoij_type>=Previous QP rhoij used for mixing if nscf>0 and rhoqpmix/=one.
!!  MPI_enreg=Information about MPI parallelization.
!!  QP_BSt<ebands_t>=QP band structure.
!!  QP_energies<Energies_type>=Simple datastructure to gather all part of total energy.
!!  nhatgrdim= 0 if pawgrnhat array is not used ; 1 otherwise
!!
!! OUTPUT
!!  QP_pawrhoij(Cryst%natom))<Pawrhoij_type>=on-site densities calculated from the QP amplitudes.
!!  qp_nhat(nfftf,Dtset%nspden)=Compensation charge density calculated from the QP amplitudes.
!!  qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)=Derivatives of the QP nhat on fine rectangular grid (and derivatives).
!!  qp_compch_sph=QP compensation charge integral inside spheres computed over spherical meshes.
!!  qp_compch_fft=QP compensation charge inside spheres computed over fine fft grid.
!!  QP_paw_ij(Cryst%natom)<Paw_ij_type>=Non-local D_ij strengths of the QP Hamiltonian.
!!  QP_paw_an(Cryst%natom)<Paw_an_type>=Various arrays related to the Hamiltonian given
!!   on ANgular mesh or ANgular moments.
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      sigma
!!
!! CHILDREN
!!      paw_an_init,paw_an_nullify,paw_ij_init,paw_ij_nullify,pawdenpot
!!      pawmknhat,pawrhoij_alloc,pawrhoij_unpack,pawrhoij_symrhoij,wfd_pawrhoij,wrtout
!!
!! SOURCE

subroutine paw_qpscgw(Wfd,nscf,nfftf,ngfftf,Dtset,Cryst,Kmesh,Psps,QP_BSt,&
&  Pawang,Pawrad,Pawtab,Pawfgrtab,prev_Pawrhoij,&
&  QP_pawrhoij,QP_paw_ij,QP_paw_an,QP_energies,qp_nhat,nhatgrdim,qp_nhatgr,qp_compch_sph,qp_compch_fft)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfftf,nscf,nhatgrdim
 real(dp),intent(out) :: qp_compch_fft,qp_compch_sph
 type(kmesh_t),intent(in) :: Kmesh
 type(crystal_t),intent(in) :: Cryst
 type(Dataset_type),intent(in) :: Dtset
 type(Pseudopotential_type),intent(in) :: Psps
 type(Pawang_type),intent(in) :: Pawang
 type(ebands_t),intent(in) :: QP_BSt
 type(Energies_type),intent(inout) :: QP_energies
 type(wfd_t),intent(inout) :: Wfd
!arrays
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(out) :: qp_nhat(nfftf,Dtset%nspden)
 real(dp),intent(out) :: qp_nhatgr(nfftf,Dtset%nspden,3*nhatgrdim)
 type(Pawfgrtab_type),intent(inout) :: Pawfgrtab(Cryst%natom)
 type(Pawrad_type),intent(in) :: Pawrad(Psps%ntypat*Psps%usepaw)
 type(Pawtab_type),intent(in) :: Pawtab(Psps%ntypat*Psps%usepaw)
 type(Pawrhoij_type),intent(inout) :: prev_Pawrhoij(Cryst%natom)
 type(Pawrhoij_type),intent(out) :: QP_pawrhoij(Cryst%natom)
 type(Paw_ij_type),intent(out) :: QP_paw_ij(Cryst%natom)
 type(Paw_an_type),intent(inout) :: QP_paw_an(Cryst%natom)

!Local variables ------------------------------
!scalars
 integer :: choice,cplex,cplex_rhoij,has_dijU,has_dijso,iat,ider,idir,ipert
 integer :: izero,nkxc1,nspden_rhoij,nzlmopt
 integer :: option,optrhoij,usexcnhat
 character(len=500) :: msg
!arrays
 real(dp) :: k0(3)

!************************************************************************

 ABI_UNUSED(Kmesh%nibz)
 !
 !  * 0 if Vloc in atomic data is Vbare    (Blochl s formulation)
 !  * 1 if Vloc in atomic data is VH(tnzc) (Kresse s formulation)
 usexcnhat=MAXVAL(Pawtab(:)%usexcnhat)
 !
 ! Calculate new rhoij_qp from updated Cprj_ibz, note use_rhoij_=1.
 call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&                  nspden=Dtset%nspden,spnorb=Dtset%pawspnorb,cpxocc=Dtset%pawcpxocc)
 call pawrhoij_alloc(QP_pawrhoij,cplex_rhoij,nspden_rhoij,Dtset%nspinor,Dtset%nsppol,Cryst%typat,&
&                    pawtab=Pawtab,use_rhoij_=1,use_rhoijres=1)

 ! FIXME kptop should be passed via Kmesh, in GW time reversal is always assumed.
 call wfd%pawrhoij(Cryst,QP_Bst,Dtset%kptopt,QP_pawrhoij,Dtset%pawprtvol)
 !
 ! * Symmetrize QP $\rho_{ij}$.
 choice=1; optrhoij=1; ipert=0
 call pawrhoij_symrhoij(QP_pawrhoij,QP_pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
&              Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
&              Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)

 ! ======================
 ! ==== Make QP nhat ====
 ! ======================
 cplex=1; ider=2*nhatgrdim; idir=0; ipert=0; izero=0; k0(:)=zero

 call pawmknhat(qp_compch_fft,cplex,ider,idir,ipert,izero,Cryst%gprimd,&
&  Cryst%natom,Cryst%natom,nfftf,ngfftf,nhatgrdim,Dtset%nspden,Cryst%ntypat,Pawang,&
&  Pawfgrtab,qp_nhatgr,qp_nhat,QP_Pawrhoij,QP_Pawrhoij,Pawtab,k0,Cryst%rprimd,&
&  Cryst%ucvol,Dtset%usewvl,Cryst%xred)

 ! Allocate quantities related to the PAW spheres for the QP Hamiltonian.
 ! TODO call paw_ij_init in scfcv and respfn, fix small issues
 has_dijso=Dtset%pawspnorb; has_dijU=merge(0,1,Dtset%usepawu==0)

 call paw_ij_nullify(QP_paw_ij)
 call paw_ij_init(QP_paw_ij,cplex,Dtset%nspinor,Dtset%nsppol,&
&  Dtset%nspden,Dtset%pawspnorb,Cryst%natom,Cryst%ntypat,Cryst%typat,Pawtab,&
&  has_dij=1,has_dijhartree=1,has_dijhat=1,has_dijxc=1,has_dijxc_hat=1,has_dijxc_val=1,&
&  has_dijso=has_dijso,has_dijU=has_dijU,has_exexch_pot=1,has_pawu_occ=1)

 call paw_an_nullify(QP_paw_an); nkxc1=0 ! No kernel
 call paw_an_init(QP_paw_an,Cryst%natom,Cryst%ntypat,nkxc1,0,Dtset%nspden,&
&     cplex,Dtset%pawxcdev,Cryst%typat,Pawang,Pawtab,has_vxc=1,has_vxcval=1)

 ! =====================================================
 ! ==== Optional mixing of the PAW onsite densities ====
 ! =====================================================
 if (nscf>0 .and. (ABS(Dtset%rhoqpmix-one)>tol12) ) then
   write(msg,'(a,f6.3)')' sigma: mixing on-site QP rho_ij densities using rhoqpmix= ',Dtset%rhoqpmix
   call wrtout(std_out,msg,'COLL')
   ! qp_rhor = prev_rhor + Dtset%rhoqpmix*(qp_rhor-prev_rhor)

   call pawrhoij_unpack(QP_Pawrhoij)   ! Unpack new QP %rhoijp
   call pawrhoij_unpack(prev_Pawrhoij) ! Unpack previous QP %rhoijp

   do iat=1,Cryst%natom
     QP_pawrhoij(iat)%rhoij_ = prev_Pawrhoij(iat)%rhoij_ &
&      + Dtset%rhoqpmix * (QP_pawrhoij(iat)%rhoij_ - prev_pawrhoij(iat)%rhoij_)

     prev_pawrhoij(iat)%use_rhoij_=0
     ABI_DEALLOCATE(prev_pawrhoij(iat)%rhoij_)
   end do
   !
   ! * Re-Symmetrize mixed QP $\rho_{ij}$.
   choice=1; optrhoij=1; ipert=0
   call pawrhoij_symrhoij(QP_pawrhoij,QP_pawrhoij,choice,Cryst%gprimd,Cryst%indsym,ipert,&
&                Cryst%natom,Cryst%nsym,Cryst%ntypat,optrhoij,Pawang,Dtset%pawprtvol,&
&                Pawtab,Cryst%rprimd,Cryst%symafm,Cryst%symrec,Cryst%typat)
 end if

 do iat=1,Cryst%natom
   QP_pawrhoij(iat)%use_rhoij_=0
   ABI_DEALLOCATE(QP_pawrhoij(iat)%rhoij_)
 end do
 !
 ! =================================================================================
 ! ==== Evaluate on-site energies, potentials, densities using (mixed) QP rhoij ====
 ! =================================================================================
 ! * Initialize also "lmselect" (index of non-zero LM-moments of densities).

 nzlmopt=-1; option=0; qp_compch_sph=greatest_real

 call pawdenpot(qp_compch_sph,QP_energies%e_paw,QP_energies%e_pawdc,&
&  ipert,Dtset%ixc,Cryst%natom,Cryst%natom,Dtset%nspden,&
&  Cryst%ntypat,Dtset%nucdipmom,nzlmopt,option,QP_paw_an,QP_paw_an,&
&  QP_paw_ij,Pawang,Dtset%pawprtvol,Pawrad,QP_pawrhoij,Dtset%pawspnorb,&
&  Pawtab,Dtset%pawxcdev,Dtset%spnorbscl,Dtset%xclevel,Dtset%xc_denpos,Cryst%ucvol,Psps%znuclpsp)

end subroutine paw_qpscgw
!!***

end module m_sigma_driver
!!***

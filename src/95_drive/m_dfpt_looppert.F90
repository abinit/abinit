!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_loopert
!! NAME
!!  m_dfpt_loopert
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (XG, DRH, MB, XW, MT, SPr, MJV)
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

module m_dfpt_loopert

 use defs_basis
 use m_dtset
 use m_dtfil
 use defs_wvltypes
 use m_efmas_defs
 use m_abicore
 use m_xmpi
 use m_errors
 use m_wfk
 use m_wffile
 use m_io_redirect
 use m_paral_pert
 use m_nctk
 use m_ddb
 use m_wfd
#ifdef HAVE_NETCDF
 use netcdf
#endif
 use m_hdr
 use m_ebands

 use defs_datatypes, only : pseudopotential_type, ebands_t
 use defs_abitypes, only : MPI_type
 use m_occ,        only : getnel
 use m_ddb_hdr,    only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_io_tools,   only : file_exists, iomode_from_fname, get_unit

 use m_time,       only : timab
 use m_fstrings,   only : strcat
 use m_geometry,   only : mkrdim, metric, littlegroup_pert
 use m_exit,       only : exit_check, disable_timelimit
 use m_atomdata,   only : atom_gauss
 use m_eig2d,      only : eigr2d_init,eigr2d_t, eigr2d_ncwrite,eigr2d_free, &
                          gkk_t, gkk_init, gkk_ncwrite,gkk_free, outbsd, eig2stern
 use m_crystal,    only : crystal_init, crystal_t
 use m_efmas,      only : efmas_main, efmas_analysis, print_efmas
 use m_fft,        only : fourdp
 use m_fftcore,    only : fftcore_set_mixprec
 use m_kg,         only : getcut, getmpw, kpgio, getph
 use m_iowf,       only : outwf
 use m_ioarr,      only : read_rhor
 use m_pawang,     only : pawang_type, pawang_init, pawang_free
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_paw_ij,     only : paw_ij_type
 use m_pawfgrtab,  only : pawfgrtab_type
 use m_pawrhoij,   only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_bcast, pawrhoij_copy, &
                          pawrhoij_nullify, pawrhoij_redistribute, pawrhoij_inquire_dim
 use m_pawcprj,    only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_copy, pawcprj_getdim
 use m_pawfgr,     only : pawfgr_type
 use m_paw_sphharm,only : setsym_ylm
 use m_rf2,        only : rf2_getidirs
 use m_iogkk,      only : outgkk
 use m_inwffil,    only : inwffil
 use m_spacepar,   only : rotate_rho, setsym
 use m_initylmg,   only : initylmg
 use m_dfpt_scfcv, only : dfpt_scfcv
 use m_dfpt_mkrho, only : dfpt_mkrho
 use m_mpinfo,     only : initmpi_band, distrb2, proc_distrb_cycle, proc_distrb_nband
 use m_atm2fft,    only : dfpt_atm2fft
 use m_berrytk,    only : smatrix
 use m_common,     only : prteigrs
 use m_fourier_interpol, only : transgrid
 use m_mkcore,     only : dfpt_mkcore
 use m_mklocl,     only : dfpt_vlocal, vlocalstr
 use m_cgprj,      only : ctocprj
 use m_symkpt,     only : symkpt

 implicit none

 private
!!***

 public :: dfpt_looppert
 public :: eigen_meandege
!!***

contains
!!***

!!****f* ABINIT/dfpt_looppert
!! NAME
!! dfpt_looppert
!!
!! FUNCTION
!! Loop over perturbations
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  codvsn=code version
!!  cpus=cpu time limit in seconds
!!  dim_eigbrd=1 if eigbrd is to be computed
!!  dim_eig2nkq=1 if eig2nkq is to be computed
!!  doccde(mband*nkpt*nsppol)=derivative of occupancies wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  dyew(2,3,natom,3,natom)=Ewald part of the dynamical matrix
!!  dyfrlo(3,3,natom)=frozen wavefunctions local part of the dynamical matrix
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=frozen wavefunctions non-local part of the dynamical matrix
!!  dyfrx1(2,3,natom,3,natom)=frozen wf nonlin. xc core corr.(2) part of the dynamical matrix
!!  dyfrx2(3,3,natom)=frozen wf nonlin. xc core corr.(2) part of the dynamical matrix
!!  dyvdw(2,3,natom,3,natom*usevdw)=vdw DFT-D part of the dynamical matrix
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  eigbrd(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eigbrd)=boradening factors for the electronic eigenvalues
!!  eig2nkq(2,mband*nsppol,nkpt,3,natom,3,natom*dim_eig2nkq)=second derivatives of the electronic eigenvalues
!!  eltcore(6,6)=core contribution to the elastic tensor
!!  elteew(6+3*natom,6)=Ewald contribution to the elastic tensor
!!  eltfrhar(6,6)=Hartree contribution to the elastic tensor
!!  eltfrkin(6,6)=kinetic contribution to the elastic tensor
!!  eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!!  eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!!  eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!!  eltvdw(6+3*natom,6*usevdw)=vdw DFT-D part of the elastic tensor
!!  fermie=fermi energy (Hartree)
!!  iexit=index of "exit" on first line of file (0 if not found)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhotoxc.f)
!!  mkmem =Number of k points treated by this node (GS data)
!!  mkqmem=Number of k+q points treated by this node (GS data)
!!  mk1mem=Number of k points treated by this node (RF data)
!!  mpert=maximum number of ipert
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)= # atoms of each type.
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  nkpt=number of k points
!!  nkxc=second dimension of the kxc array
!!  nspden=number of spin-density components
!!  nsym=number of symmetry elements in space group
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pertsy(3,mpert)=set of perturbations that form a basis for all other perturbations
!!  prtbbb=if 1, bbb decomposition, also dimension d2bbb
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rfpert(mpert)=array defining the type of perturbations that have to be computed
!!                1   ->   element has to be computed explicitely
!!               -1   ->   use symmetry operations to obtain the corresponding element
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrec(3,3,nsym)=3x3 matrices of the group symmetries (reciprocal space)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!!  vtrial(nfftf,nspden)=GS potential (Hartree)
!!  vxc(nfftf,nspden)=Exchange-Correlation GS potential (Hartree)
!!  vxcavg=average of vxc potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  ddkfil(3)=unit numbers for the three possible ddk files
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some second order derivatives
!!  d2lo(2,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  d2ovl(2,mpert,3,mpert*usepaw)=1st-order change of WF overlap contributions to the 2DTEs
!!  etotal=total energy (sum of 8 contributions) (hartree)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      appdig,atom_gauss,crystal_free,crystal_init,ctocprj,ddb_free
!!      ddb_from_file,ddb_hdr_free,ddb_hdr_init,ddb_hdr_open_write,dfpt_atm2fft
!!      dfpt_init_mag1,dfpt_mkcore,dfpt_mkrho,dfpt_prtene,dfpt_scfcv
!!      dfpt_vlocal,disable_timelimit,distrb2,dtset_copy,dtset_free,ebands_free
!!      ebands_init,efmas_main,efmas_analysis,eig2stern,eigen_meandege,eigr2d_free,eigr2d_init
!!      eigr2d_ncwrite,exit_check,fourdp,getcgqphase,getcut,getmpw,getnel,getph
!!      gkk_free,gkk_init,gkk_ncwrite,hdr_copy,hdr_free,hdr_init,hdr_update
!!      initmpi_band,initylmg,inwffil,kpgio,littlegroup_pert,localfilnam
!!      localrdfile,localredirect,localwrfile,metric,mkrdim,outbsd,outgkk,outwf
!!      pawang_free,pawang_init,pawcprj_alloc,pawcprj_copy,pawcprj_free
!!      pawcprj_getdim,pawrhoij_alloc,pawrhoij_copy,pawrhoij_free
!!      pawrhoij_nullify,prteigrs,put_eneocc_vect,read_rhor,rf2_getidirs
!!      rotate_rho,set_pert_comm,set_pert_paw,setsym,setsym_ylm,status,symkpt
!!      timab,transgrid,unset_pert_comm,unset_pert_paw,vlocalstr,wffclose
!!      wfk_open_read,wfk_read_eigenvalues,wrtout,xmpi_sum
!!
!! SOURCE

subroutine dfpt_looppert(atindx,blkflg,codvsn,cpus,dim_eigbrd,dim_eig2nkq,doccde,&
&  ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,dyvdw,&
&  dyfr_cplex,dyfr_nondiag,d2bbb,d2lo,d2nl,d2ovl,efmasdeg,efmasval,eigbrd,eig2nkq,&
&  eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
&  etotal,fermie,iexit,indsym,kxc,&
&  mkmem,mkqmem,mk1mem,mpert,mpi_enreg,my_natom,nattyp,&
&  nfftf,nhat,nkpt,nkxc,nspden,nsym,occ,&
&  paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&  pertsy,prtbbb,psps,rfpert,rf2_dirs_from_rfpert_nl,rhog,rhor,symq,symrec,timrev,&
&  usecprj,usevdw,vtrial,vxc,vxcavg,xred,clflg,occ_rbz_pert,eigen0_pert,eigenq_pert,&
&  eigen1_pert,nkpt_rbz,eigenq_fine,hdr_fine,hdr0)

!Arguments ------------------------------------
 integer, intent(in) :: dim_eigbrd,dim_eig2nkq,dyfr_cplex,dyfr_nondiag,mk1mem,mkmem,mkqmem,mpert
 integer, intent(in) :: nfftf,nkpt,nkxc,nspden,nsym,prtbbb,timrev,usecprj,usevdw
 integer, intent(out) :: iexit
 integer, intent(inout) :: my_natom
 real(dp), intent(in) :: cpus,vxcavg
 real(dp), intent(inout) :: fermie
 real(dp), intent(inout) :: etotal
 character(len=6), intent(in) :: codvsn
 type(MPI_type), intent(inout) :: mpi_enreg
 type(datafiles_type), intent(in) :: dtfil
 type(dataset_type), intent(in), target :: dtset
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type), intent(inout) :: psps
 integer, intent(in) :: atindx(dtset%natom),indsym(4,nsym,dtset%natom)
 integer, intent(in) :: nattyp(dtset%ntypat),pertsy(3,dtset%natom+6)
 integer, intent(in) :: rfpert(mpert),rf2_dirs_from_rfpert_nl(3,3),symq(4,2,nsym),symrec(3,3,nsym)
 integer, intent(out) :: ddkfil(3)
 integer, intent(inout) :: blkflg(3,mpert,3,mpert)
 integer, intent(out) :: clflg(3,mpert)
 real(dp), intent(in) :: doccde(dtset%mband*nkpt*dtset%nsppol)
 real(dp), intent(in) :: dyew(2,3,dtset%natom,3,dtset%natom)
 real(dp), intent(in) :: dyfrlo(3,3,dtset%natom)
 real(dp), intent(in) :: dyfrnl(dyfr_cplex,3,3,dtset%natom,1+(dtset%natom-1)*dyfr_nondiag)
 real(dp), intent(in) :: dyfrx1(2,3,dtset%natom,3,dtset%natom)
 real(dp), intent(in) :: dyfrx2(3,3,dtset%natom),dyvdw(2,3,dtset%natom,3,dtset%natom*usevdw)
 real(dp), intent(in) :: eltcore(6,6),elteew(6+3*dtset%natom,6),eltfrhar(6,6)
 real(dp), intent(in) :: eltfrkin(6,6),eltfrloc(6+3*dtset%natom,6)
 real(dp), intent(in) :: eltfrnl(6+3*dtset%natom,6)
 real(dp), intent(in) :: eltfrxc(6+3*dtset%natom,6),eltvdw(6+3*dtset%natom,6*usevdw)
 real(dp), intent(in) :: kxc(nfftf,nkxc),nhat(nfftf,nspden)
 real(dp), intent(in) :: occ(dtset%mband*nkpt*dtset%nsppol)
 real(dp), intent(in) :: rhog(2,nfftf),rhor(nfftf,nspden),vxc(nfftf,nspden)
 real(dp), intent(in) :: vtrial(nfftf,nspden)
 real(dp), intent(inout) :: xred(3,dtset%natom)
 real(dp), intent(inout) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)!vz_i
 real(dp), intent(inout) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert) !vz_i
 real(dp), intent(inout) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw) !vz_i
 real(dp), intent(out) :: eigbrd(2,dtset%mband*dtset%nsppol,nkpt,3,dtset%natom,3,dtset%natom*dim_eigbrd)
 real(dp), intent(out) :: eig2nkq(2,dtset%mband*dtset%nsppol,nkpt,3,dtset%natom,3,dtset%natom*dim_eig2nkq)
 type(efmasdeg_type),allocatable,intent(in) :: efmasdeg(:)
 type(efmasval_type),allocatable,intent(inout) :: efmasval(:,:)
 type(paw_an_type),allocatable,target,intent(inout) :: paw_an(:)
 type(paw_ij_type),allocatable,target,intent(inout) :: paw_ij(:)
 type(pawfgrtab_type),allocatable,target,intent(inout) :: pawfgrtab(:)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),allocatable,target,intent(inout) :: pawrhoij(:)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 real(dp),pointer :: eigen1_pert(:,:,:)
 real(dp),intent(out) :: occ_rbz_pert(:),eigen0_pert(:),eigenq_pert(:)
 real(dp),pointer :: eigenq_fine(:,:,:)
 integer, intent(out) :: nkpt_rbz
 type(hdr_type),intent(out) :: hdr0,hdr_fine

!Local variables-------------------------------
!scalars
 integer,parameter :: level=11,response=1,formeig1=1,master=0,fake_unit=-666
 integer :: ask_accurate,band_index,bantot,bantot_rbz,bdeigrf,bdtot1_index,nsppol,nspinor,band2tot_index
 integer :: bdtot_index,choice,cplex,cplex_rhoij,dim_eig2rf,formeig
 integer :: gscase,iband,iblok,icase,icase_eq,idir,idir0,idir1,idir2,idir_eq,idir_dkdk,ierr
 integer :: ii,ikpt,ikpt1,jband,initialized,iorder_cprj,ipert,ipert_cnt,ipert_eq,ipert_me,ireadwf0
 integer :: iscf_mod,iscf_mod_save,isppol,istr,isym,mcg,mcgq,mcg1,mcprj,mcprjq,mband
 integer :: mcgmq,mcg1mq,mpw1_mq !+/-q duplicates
 integer :: nband_me, wfk_unt, icg, ibd, ikg
 integer :: maxidir,me,mgfftf,mkmem_rbz,mk1mem_rbz,mkqmem_rbz,mpw,mpw1,my_nkpt_rbz
 integer :: n3xccc,nband_k,ncpgr,ndir,nkpt_eff,nkpt_max,nline_save,nmatel,npert_io,npert_me,nspden_rhoij
 integer :: nstep_save,nsym1,ntypat,nwffile,nylmgr,nylmgr1,old_comm_atom,openexit,option,optorth,optthm,pertcase
 integer :: qphase_rhoij,rdwr,rdwrpaw,spaceComm,smdelta,timrev_pert,timrev_kpt,to_compute_this_pert
 integer :: unitout,useylmgr,useylmgr1,vrsddb,dfpt_scfcv_retcode,optn2
#ifdef HAVE_NETCDF
 integer :: ncerr,ncid
#endif
 real(dp) :: boxcut,dosdeltae,eberry,ecore,ecut_eff,ecutf,edocc,eei,eeig0,eew,efrhar,efrkin,efrloc
 real(dp) :: efrnl,efrx1,efrx2,ehart,ehart01,ehart1,eii,ek,ek0,ek1,ek2,eloc0
 real(dp) :: elpsp1,enl,enl0,enl1,entropy,enxc,eovl1,epaw1,evdw,exc1,fsum,gsqcut,maxocc,nelectkq
 real(dp) :: residm,tolwfr,tolwfr_save,toldfe_save,toldff_save,tolrff_save,tolvrs_save
 real(dp) :: ucvol, eig1_r, eig1_i
 real(dp) :: residm_mq !+/-q duplicates
 logical,parameter :: paral_pert_inplace=.true.,remove_inv=.false.
 logical :: first_entry,found_eq_gkk,t_exist,paral_atom,write_1wfk,init_rhor1
 logical :: kramers_deg
 character(len=fnlen) :: dscrpt,fiden1i,fiwf1i,fiwf1o,fiwfddk,fnamewff(4),gkkfilnam,fname,filnam
 character(len=500) :: message
 type(crystal_t) :: crystal,ddb_crystal
 type(dataset_type), pointer :: dtset_tmp
 type(ebands_t) :: ebands_k,ebands_kq,gkk_ebands
 type(ebands_t) :: ebands_kmq !+/-q duplicates
 type(eigr2d_t)  :: eigr2d,eigi2d
 type(gkk_t)     :: gkk2d
 type(hdr_type) :: hdr,hdr_den,hdr_tmp
 type(ddb_hdr_type) :: ddb_hdr
 type(pawang_type) :: pawang1
 type(wffile_type) :: wff1,wffgs,wffkq,wffnow,wfftgs,wfftkq
 type(wfk_t) :: ddk_f(4)
 type(wfk_t) :: wfk0, wfkq, wfk1
 type(wvl_data) :: wvl
!arrays
 integer :: eq_symop(3,3),ngfftf(18),file_index(4),rfdir(9),rf2dir(9),rf2_dir1(3),rf2_dir2(3)
 integer,allocatable :: blkflg_save(:,:,:,:),dimcprj_srt(:),dummy(:),dyn(:),indkpt1(:),indkpt1_tmp(:)
 integer,allocatable :: indsy1(:,:,:),irrzon1(:,:,:),istwfk_rbz(:),istwfk_pert(:,:,:)
 integer,allocatable :: kg(:,:),kg1(:,:),nband_rbz(:),npwar1(:),npwarr(:),npwtot(:)
 integer,allocatable :: kg1_mq(:,:),npwar1_mq(:),npwtot1_mq(:) !+q/-q duplicates
 integer,allocatable :: npwtot1(:),npwar1_pert(:,:),npwarr_pert(:,:),npwtot_pert(:,:)
 integer,allocatable :: pert_calc(:,:),pert_tmp(:,:),bz2ibz_smap(:,:)
 integer,allocatable :: symaf1(:),symaf1_tmp(:),symrc1(:,:,:),symrl1(:,:,:),symrl1_tmp(:,:,:)
 integer,allocatable :: kpt_tmp (:,:)
 integer, pointer :: old_atmtab(:)
 logical, allocatable :: distrbflags(:,:,:)
 real(dp) :: dielt(3,3),gmet(3,3),gprimd(3,3),rmet(3,3),rprimd(3,3),tsec(2)
 real(dp),allocatable :: buffer1(:,:,:,:,:),cg(:,:),cg1(:,:),cg1_active(:,:),cg0_pert(:,:)
 real(dp),allocatable :: cg1_pert(:,:,:,:),cgq(:,:),gh0c1_pert(:,:,:,:)
 real(dp),allocatable :: doccde_rbz(:),docckqde(:)
 real(dp),allocatable :: gh1c_pert(:,:,:,:),eigen0(:),eigen0_copy(:),eigen1(:),eigen1_mean(:)
 real(dp),allocatable :: eigenq(:),gh1c_set(:,:),gh0c1_set(:,:),kpq(:,:)
 real(dp),allocatable :: kpq_rbz(:,:),kpt_rbz(:,:),occ_pert(:),occ_rbz(:),occkq(:),kpt_rbz_pert(:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),phnons1(:,:,:),resid(:),rhog1(:,:)
 real(dp),allocatable :: rhor1_save(:,:,:)
 real(dp),allocatable :: rhor1(:,:),rho1wfg(:,:),rho1wfr(:,:),tnons1(:,:),tnons1_tmp(:,:)
 real(dp),allocatable :: rhor1_pq(:,:),rhor1_mq(:,:),rhog1_pq(:,:),rhog1_mq(:,:)          !+q/-q duplicates
 real(dp),allocatable :: cg_mq(:,:),cg1_mq(:,:),resid_mq(:)                   !
 real(dp),allocatable :: cg1_active_mq(:,:),occk_mq(:)                 !
 real(dp),allocatable :: kmq(:,:),kmq_rbz(:,:),gh0c1_set_mq(:,:)        !
 real(dp),allocatable :: eigen_mq(:),gh1c_set_mq(:,:),docckde_mq(:),eigen1_mq(:)          !
 real(dp),allocatable :: vpsp1(:),work(:),wtk_folded(:),wtk_rbz(:),xccc3d1(:)
 real(dp),allocatable :: ylm(:,:),ylm1(:,:),ylmgr(:,:,:),ylmgr1(:,:,:),zeff(:,:,:)
 real(dp),allocatable :: phasecg(:,:),gauss(:,:)
 real(dp),allocatable :: gkk(:,:,:,:,:)
 type(pawcprj_type),allocatable :: cprj(:,:),cprjq(:,:)
 type(paw_ij_type),pointer :: paw_ij_pert(:)
 type(paw_an_type),pointer :: paw_an_pert(:)
 type(pawfgrtab_type),pointer :: pawfgrtab_pert(:)
 type(pawrhoij_type),allocatable :: pawrhoij1(:)
 type(pawrhoij_type),pointer :: pawrhoij_pert(:)
 type(ddb_type) :: ddb


! ***********************************************************************

 DBG_ENTER("COLL")

 _IBM6("In dfpt_looppert")

 call timab(141,1,tsec)

!Structured debugging if prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,' dfpt_looppert : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

 dfpt_scfcv_retcode = -1
 nsppol = dtset%nsppol; nspinor = dtset%nspinor

 kramers_deg=.true.
 if (dtset%tim1rev==0) then
   kramers_deg=.false.
 end if

!Obtain dimensional translations in reciprocal space gprimd,
!metrics and unit cell volume, from rprimd. Also output rprimd, gprimd and ucvol
 call mkrdim(dtset%acell_orig(1:3,1),dtset%rprim_orig(1:3,1:3,1),rprimd)
 call metric(gmet,gprimd,std_out,rmet,rprimd,ucvol)

 call crystal_init(dtset%amu_orig(:,1),crystal,dtset%spgroup,dtset%natom,dtset%npsp,&
& psps%ntypat,dtset%nsym,rprimd,dtset%typat,xred,dtset%ziontypat,dtset%znucl,1,&
& dtset%nspden==2.and.dtset%nsppol==1,remove_inv,psps%title,&
& symrel=dtset%symrel,tnons=dtset%tnons,symafm=dtset%symafm)

!Get FFT grid(s) sizes (be careful !) See NOTES in the comments at the beginning of respfn.F90
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   mgfftf=pawfgr%mgfft;ngfftf(:)=pawfgr%ngfft(:)
   ecutf=dtset%pawecutdg
 else
   mgfftf=dtset%mgfft;ngfftf(:)=dtset%ngfft(:)
   ecutf=dtset%ecut
 end if
 ecut_eff=dtset%ecut*(dtset%dilatmx)**2

!Compute large sphere cut-off gsqcut
 if (psps%usepaw==1) then
   call wrtout(std_out,ch10//' FFT (fine) grid used for densities/potentials:','COLL')
 end if
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,dtset%qptn,ngfftf)

!Various initializations/allocations
 iscf_mod=dtset%iscf
 ntypat=psps%ntypat
 nkpt_max=50;if (xmpi_paral==1) nkpt_max=-1
 paral_atom=(dtset%natom/=my_natom)
 cplex=2-timrev !cplex=2 ! DEBUG: impose cplex=2
 first_entry=.true.
 initialized=0
 ecore=zero ; ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero ; eii=zero
 clflg(:,:)=0 ! Array on calculated perturbations for eig2rf
 if (psps%usepaw==1) then
   ABI_ALLOCATE(dimcprj_srt,(dtset%natom))
   call pawcprj_getdim(dimcprj_srt,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
 end if

!Save values of SCF cycle parameters
 iscf_mod_save = iscf_mod
 nstep_save = dtset%nstep
 nline_save = dtset%nline
 tolwfr_save = dtset%tolwfr
 toldfe_save = dtset%toldfe
 toldff_save = dtset%toldff
 tolrff_save = dtset%tolrff
 tolvrs_save = dtset%tolvrs

!This dtset will be used in dfpt_scfcv to force non scf calculations for equivalent perturbations
 nullify(dtset_tmp)
 if (dtset%prepgkk/=0) then ! .and. dtset%use_nonscf_gkk==1) then !Later uncomment this - in scf case rhor1_save is used below only for testing
   ABI_DATATYPE_ALLOCATE(dtset_tmp,)
   dtset_tmp = dtset%copy()
 else
   dtset_tmp => dtset
 end if

!Generate the 1-dimensional phases
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*dtset%natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*dtset%natom))
 call getph(atindx,dtset%natom,dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),ph1d,xred)
 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(atindx,dtset%natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

!!Determine existence of pertubations and of pertubation symmetries
!!Create array with pertubations which have to be calculated
! ABI_ALLOCATE(pert_tmp,(3*mpert))
! ipert_cnt=0
! do ipert=1,mpert
!   do idir=1,3
!     if( rfpert(ipert)==1 .and. dtset%rfdir(idir) == 1 )then
!       if ((pertsy(idir,ipert)==1).or.&
!&       ((dtset%prepanl == 1).and.(ipert == dtset%natom+2)).or.&
!&       ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
!         ipert_cnt = ipert_cnt+1;
!         pert_tmp(ipert_cnt) = idir+(ipert-1)*3
!       else
!         write(message, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
!&         ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
!&         ' symmetric of a previously calculated perturbation.',ch10,&
!&         ' So, its SCF calculation is not needed.',ch10
!         call wrtout(std_out,message,'COLL')
!         call wrtout(ab_out,message,'COLL')
!       end if ! Test of existence of symmetry of perturbation
!     end if ! Test of existence of perturbation
!   end do
! end do
! ABI_ALLOCATE(pert_calc,(ipert_cnt))
! do icase=1,ipert_cnt
!   pert_calc(icase)=pert_tmp(icase)
! end do
! ABI_DEALLOCATE(pert_tmp)

!Initialize rf2dir :
 rf2_dir1(1:3)=dtset%rf2_pert1_dir(1:3)
 rf2_dir2(1:3)=dtset%rf2_pert2_dir(1:3)
 if (sum(rf2_dir1)==3.and.sum(rf2_dir2)==3.and.dtset%prepanl==1) then
!  Diagonal terms :
   rf2dir(1) = rf2_dirs_from_rfpert_nl(1,1)
   rf2dir(2) = rf2_dirs_from_rfpert_nl(2,2)
   rf2dir(3) = rf2_dirs_from_rfpert_nl(3,3)
!  Upper triangular terms :
   rf2dir(4) = rf2_dirs_from_rfpert_nl(2,3)
   rf2dir(5) = rf2_dirs_from_rfpert_nl(1,3)
   rf2dir(6) = rf2_dirs_from_rfpert_nl(1,2)
!  Lower triangular terms :
   rf2dir(7) = rf2_dirs_from_rfpert_nl(3,2)
   rf2dir(8) = rf2_dirs_from_rfpert_nl(3,1)
   rf2dir(9) = rf2_dirs_from_rfpert_nl(2,1)
 else
!  Diagonal terms :
   rf2dir(1) = rf2_dir1(1)*rf2_dir2(1)
   rf2dir(2) = rf2_dir1(2)*rf2_dir2(2)
   rf2dir(3) = rf2_dir1(3)*rf2_dir2(3)
!  Upper triangular terms :
   rf2dir(4) = rf2_dir1(2)*rf2_dir2(3)
   rf2dir(5) = rf2_dir1(1)*rf2_dir2(3)
   rf2dir(6) = rf2_dir1(1)*rf2_dir2(2)
!  Lower triangular terms :
   rf2dir(7) = rf2_dir1(3)*rf2_dir2(2)
   rf2dir(8) = rf2_dir1(3)*rf2_dir2(1)
   rf2dir(9) = rf2_dir1(2)*rf2_dir2(1)
 end if

!Determine existence of pertubations and of pertubation symmetries
!Create array with pertubations which have to be calculated
 ABI_ALLOCATE(pert_tmp,(3,3*(dtset%natom+6)+18))
 ipert_cnt=0
 do ipert=1,mpert
   if (ipert<dtset%natom+10) then
     maxidir = 3
     rfdir(1:3) = dtset%rfdir(:)
     rfdir(4:9) = 0
   else
     maxidir = 9
     rfdir(1:9) = rf2dir(:)
   end if
   do idir=1,maxidir
     to_compute_this_pert = 0
     if(ipert<dtset%natom+10 .and. rfpert(ipert)==1 .and. rfdir(idir) == 1 ) then
       if ((pertsy(idir,ipert)==1).or.&
&       ((dtset%prepanl == 1).and.(ipert == dtset%natom+2)).or.&
&       ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
         to_compute_this_pert = 1
       else
         write(message, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
&         ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
&         ' symmetric of a previously calculated perturbation.',ch10,&
&         ' So, its SCF calculation is not needed.',ch10
         call wrtout(std_out,message,'COLL')
         call wrtout(ab_out,message,'COLL')
       end if ! Test of existence of symmetry of perturbation
     else if (ipert==dtset%natom+11 .and. rfpert(ipert)==1 .and. rfdir(idir) == 1 ) then
       to_compute_this_pert = 1
     else if (ipert==dtset%natom+10 .and. rfpert(ipert)==1 .and. idir <= 6) then
       if (idir<=3) then
         if (rfdir(idir) == 1) to_compute_this_pert = 1
       else
         if (rfdir(idir) == 1 .or. rfdir(idir+3) == 1) to_compute_this_pert = 1
       end if
     else if (ipert==dtset%natom+11 .and. rfpert(ipert)==1) then
       if (idir <= 3 .and. rfdir(idir) == 1) then
         to_compute_this_pert = 1
       else if (idir>=4.and.idir<=6) then
         if (rfdir(idir) == 1 .or. rfdir(idir+3) == 1) to_compute_this_pert = 1
       else if (idir>=7.and.idir<=9) then
         if (rfdir(idir) == 1 .or. rfdir(idir-3) == 1) to_compute_this_pert = 1
       end if
     end if
     if (to_compute_this_pert /= 0) then
       ipert_cnt = ipert_cnt+1;
       pert_tmp(1,ipert_cnt) = ipert
       pert_tmp(2,ipert_cnt) = idir
!      Store "pertcase" in pert_tmp(3,ipert_cnt)
       if (ipert<dtset%natom+10) then
         pert_tmp(3,ipert_cnt) = idir + (ipert-1)*3
       else
         pert_tmp(3,ipert_cnt) = idir + (ipert-dtset%natom-10)*9 + (dtset%natom+6)*3
       end if
     end if
   end do ! idir
 end do !ipert
! do ipert=1,mpert
!   if (ipert<dtset%natom+10) then
!     maxidir = 3
!     rfdir(1:3) = dtset%rfdir(:)
!     rfdir(4:9) = 0
!   else
!     maxidir = 9
!     rfdir(1:9) = rf2dir(:)
!   end if
!   do idir=1,maxidir
!     if( rfpert(ipert)==1 .and. rfdir(idir) == 1 )then
!       to_compute_this_pert = 0
!       if (ipert>=dtset%natom+10) then
!         to_compute_this_pert = 1
!       else if ((pertsy(idir,ipert)==1).or.&
!&         ((dtset%prepanl == 1).and.(ipert == dtset%natom+2)).or.&
!&         ((dtset%prepgkk == 1).and.(ipert <= dtset%natom))  ) then
!         to_compute_this_pert = 1
!       end if
!       if (to_compute_this_pert /= 0) then
!         ipert_cnt = ipert_cnt+1;
!         pert_tmp(1,ipert_cnt) = ipert
!         pert_tmp(2,ipert_cnt) = idir
!!        Store "pertcase" in pert_tmp(3,ipert_cnt)
!         if (ipert<dtset%natom+10) then
!           pert_tmp(3,ipert_cnt) = idir + (ipert-1)*3
!         else
!           pert_tmp(3,ipert_cnt) = idir + (ipert-dtset%natom-10)*9 + (dtset%natom+6)*3
!         end if
!       else
!         write(message, '(a,a,i4,a,i4,a,a,a,a,a,a)' )ch10,&
!&         ' The perturbation idir=',idir,'  ipert=',ipert,' is',ch10,&
!&         ' symmetric of a previously calculated perturbation.',ch10,&
!&         ' So, its SCF calculation is not needed.',ch10
!         call wrtout(std_out,message,'COLL')
!         call wrtout(ab_out,message,'COLL')
!       end if ! Test of existence of symmetry of perturbation
!     end if ! Test of existence of perturbation
!   end do
! end do
 ABI_ALLOCATE(pert_calc,(3,ipert_cnt))
 do icase=1,ipert_cnt
   pert_calc(:,icase)=pert_tmp(:,icase)
 end do
 ABI_DEALLOCATE(pert_tmp)

 if (dtset%prepgkk/=0) then ! .and. dtset%use_nonscf_gkk==1) then !Later uncomment this - in scf case rhor1_save is used below only for testing
   ABI_ALLOCATE(rhor1_save,(cplex*nfftf,nspden,ipert_cnt))
   rhor1_save=zero
   ABI_ALLOCATE(blkflg_save,(3,mpert,3,mpert))
 end if

! Initialize quantities for netcdf print
 ABI_ALLOCATE(eigen0_copy,(dtset%mband*nkpt*dtset%nsppol))
 eigen0_copy(:)=zero

! SP : Retreval of the DDB information and computing of effective charge and
! dielectric tensor
 ABI_ALLOCATE(zeff,(3,3,dtset%natom))
 if (dtset%getddb .ne. 0 .or. dtset%irdddb .ne. 0 ) then
   filnam = dtfil%filddbsin  !'test_DDB'
   ABI_ALLOCATE(dummy,(dtset%natom))
   call ddb_from_file(ddb,filnam,1,dtset%natom,0,dummy,ddb_crystal,mpi_enreg%comm_world)
!  Get Dielectric Tensor and Effective Charges
!  (initialized to one_3D and zero if the derivatives are not available in the DDB file)
   iblok = ddb%get_dielt_zeff(ddb_crystal,1,0,0,dielt,zeff)
   call ddb_crystal%free()
   call ddb%free()
   ABI_DEALLOCATE(dummy)
 end if

!%%%% Parallelization over perturbations %%%%%
!*Define file output/log file names
 npert_io=ipert_cnt;if (dtset%nppert<=1) npert_io=0
 call localfilnam(mpi_enreg%comm_pert,mpi_enreg%comm_cell_pert,mpi_enreg%comm_world,dtfil%filnam_ds,'_PRT',npert_io)
!Compute the number of perturbation done by the current cpu
 if(mpi_enreg%paral_pert==1) then
   npert_me = 0 ; ipert_me = 0
   do icase=1,ipert_cnt
     if (mpi_enreg%distrb_pert(icase)==mpi_enreg%me_pert) npert_me=npert_me +1
   end do
 end if

!*Redefine communicators
 call set_pert_comm(mpi_enreg,dtset%nppert)

!*Redistribute PAW on-site data
 nullify(old_atmtab,pawfgrtab_pert,pawrhoij_pert,paw_an_pert,paw_ij_pert)
 if (paral_pert_inplace) then
   call set_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,&
&   paw_an,paw_ij,pawfgrtab,pawrhoij)
   pawfgrtab_pert=>pawfgrtab ; pawrhoij_pert=>pawrhoij
   paw_an_pert   =>paw_an    ; paw_ij_pert  =>paw_ij

 else
   call set_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,&
&   paw_an,paw_ij,pawfgrtab,pawrhoij,&
&   paw_an_out=paw_an_pert,paw_ij_out=paw_ij_pert,&
&   pawfgrtab_out=pawfgrtab_pert,pawrhoij_out=pawrhoij_pert)

 end if

 ! We can handle the time limit in dfpt_scfcv in a robust manner only if we have one perturbation.
 if (ipert_cnt > 1 .or. mpi_enreg%paral_pert == 1) call disable_timelimit()

!Loop on perturbations
!==========================================================================
 do icase=1,ipert_cnt
   _IBM6("In loop on perts")

!  %%%% Parallelization over perturbations %%%%%
!  Select the perturbations treated by curent processor
   if(mpi_enreg%paral_pert==1) then
     if (mpi_enreg%distrb_pert(icase)/=mpi_enreg%me_pert) cycle
   end if

   ! Redefine output/log files
   call localwrfile(mpi_enreg%comm_cell,icase,npert_io,mpi_enreg%paral_pert,0)

   ! Set precision for FFT libs.
   ii = fftcore_set_mixprec(dtset%mixprec)

!!  Retrieve type and direction of the perturbation
!   if (pert_calc(icase) <= dtset%natom*3) then
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!     ipert=( (pert_calc(icase)-idir) / 3 + 1)
!   else if (pert_calc(icase) <= dtset%natom*3+4) then
!     ipert = dtset%natom + ((pert_calc(icase) - 3*dtset%natom - 1) / 3) + 1
!     idir = mod(pert_calc(icase),3)
!     if (idir==0) idir=3
!   else
!     ipert = dtset%natom + ((pert_calc(icase) - 3*(dtset%natom+4) - 1) / 9) + 1
!   end if
!   pertcase=idir+(ipert-1)*3

!  Retrieve type and direction of the perturbation
   ipert = pert_calc(1,icase)
   idir = pert_calc(2,icase)
   istr=idir
   pertcase = pert_calc(3,icase)

!  Init MPI communicator
   spaceComm=mpi_enreg%comm_cell
   me=mpi_enreg%me_cell

!  ===== Describe the perturbation in output/log file
   _IBM6("IBM6 before print perturbation")

   write(message, '(a,80a,a,a,3f10.6)' ) ch10,('-',ii=1,80),ch10,&
&   ' Perturbation wavevector (in red.coord.) ',dtset%qptn(:)
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   if(ipert>=1 .and. ipert<=dtset%natom)then
     write(message, '(a,i4,a,i4)' )' Perturbation : displacement of atom',ipert,'   along direction',idir
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     if(iscf_mod == -3)then
       write(message, '(a,a,a,a,a,a,a,a)' )ch10,&
&       ' dfpt_looppert : COMMENT -',ch10,&
&       '  The first-order density is imposed to be zero (iscf=-3).',ch10,&
&       '  Although this is strange in the case of phonons,',ch10,&
&       '  you are allowed to do so.'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
     end if
   else if(ipert==dtset%natom+1)then
     write(message,'(a,i4)')' Perturbation : derivative vs k along direction',idir
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     if( iscf_mod /= -3 )then
       write(message, '(4a)' )ch10,&
&       ' dfpt_looppert : COMMENT -',ch10,&
&       '  In a d/dk calculation, iscf is set to -3 automatically.'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       iscf_mod=-3
     end if
     if( abs(dtset%dfpt_sciss) > 1.0d-8 )then
       write(message, '(a,a,a,a,f14.8,a,a)' )ch10,&
&       ' dfpt_looppert : WARNING -',ch10,&
&       '  Value of dfpt_sciss=',dtset%dfpt_sciss,ch10,&
&       '  Scissor with d/dk calculation : you are using a "naive" approach !'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
     end if
   else if(ipert==dtset%natom+2)then
     write(message, '(a,i4)' )' Perturbation : homogeneous electric field along direction',idir
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     if( iscf_mod == -3 )then
       write(message, '(a,a,a,a,a,a)' )ch10,&
&       ' dfpt_looppert : COMMENT -',ch10,&
&       '  The first-order density is imposed to be zero (iscf=-3).',ch10,&
&       '  This corresponds to a calculation without local fields.'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
     end if
   else if(ipert==dtset%natom+10.or.ipert==dtset%natom+11)then
     call rf2_getidirs(idir,idir1,idir2)
     if(ipert==dtset%natom+10)then
       write(message,'(2(a,i1))') ' Perturbation : 2nd derivative wrt k, idir1 = ',idir1,&
&       ' idir2 = ',idir2
     else
       write(message,'(2(a,i1),a)') ' Perturbation : 2nd derivative wrt k (idir1 =',idir1,&
&       ') and Efield (idir2 =',idir2,')'
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
     if( iscf_mod /= -3 )then
       write(message, '(4a)' )ch10,&
&       ' dfpt_looppert : COMMENT -',ch10,&
&       '  In this case, iscf is set to -3 automatically.'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
       iscf_mod=-3
     end if
     if( abs(dtset%dfpt_sciss) > 1.0d-8 )then
       write(message, '(a,a,a,a,f14.8,a,a)' )ch10,&
&       ' dfpt_looppert : WARNING -',ch10,&
&       '  Value of dfpt_sciss=',dtset%dfpt_sciss,ch10,&
&       '  Scissor with d/dk calculation : you are using a "naive" approach !'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
     end if
     ABI_ALLOCATE(occ_pert,(dtset%mband*nkpt*dtset%nsppol))
     occ_pert(:) = occ(:) - occ(1)
     maxocc = maxval(abs(occ_pert))
     if (maxocc>1.0d-6.and.abs(maxocc-occ(1))>1.0d-6) then ! True if non-zero occupation numbers are not equal
       write(message, '(3a)' ) ' ipert=natom+10 or 11 does not work for a metallic system.',ch10,&
       ' This perturbation will not be computed.'
       MSG_WARNING(message)
       ABI_DEALLOCATE(occ_pert)
       cycle
     end if
     ABI_DEALLOCATE(occ_pert)
   else if(ipert>dtset%natom+11 .or. ipert<=0 )then
     write(message, '(a,i0,3a)' ) &
&     'ipert= ',ipert,' is outside the [1,dtset%natom+11] interval.',ch10,&
&     'This perturbation is not (yet) allowed.'
     MSG_BUG(message)
   end if
!  Initialize the diverse parts of energy :
   eew=zero ; evdw=zero ; efrloc=zero ; efrnl=zero ; efrx1=zero ; efrx2=zero
   efrhar=zero ; efrkin=zero
   if(ipert<=dtset%natom)then
     eew=dyew(1,idir,ipert,idir,ipert)
     if (usevdw==1) evdw=dyvdw(1,idir,ipert,idir,ipert)
     efrloc=dyfrlo(idir,idir,ipert)
     if (dyfr_nondiag==0) efrnl=dyfrnl(1,idir,idir,ipert,1)
     if (dyfr_nondiag/=0) efrnl=dyfrnl(1,idir,idir,ipert,ipert)
     efrx1=dyfrx1(1,idir,ipert,idir,ipert)
     efrx2=dyfrx2(idir,idir,ipert)
   else if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
!    istr = 1,2,...,6 and indicates the cartesian strain component
     if(ipert==dtset%natom+4) istr=idir+3
     eii=eltcore(istr,istr)
     eew=elteew(istr,istr)
     if (usevdw==1) evdw=eltvdw(istr,istr)
     efrhar=eltfrhar(istr,istr)
     efrkin=eltfrkin(istr,istr)
     efrloc=eltfrloc(istr,istr)
     efrnl=eltfrnl(istr,istr)
     efrx1=eltfrxc(istr,istr)
   end if

!  Determine the subset of symmetry operations (nsym1 operations)
!  that leaves the perturbation invariant, and initialize corresponding arrays
!  symaf1, symrl1, tnons1 (and pawang1%zarot, if PAW)..
   ABI_ALLOCATE(symaf1_tmp,(nsym))
   ABI_ALLOCATE(symrl1_tmp,(3,3,nsym))
   ABI_ALLOCATE(tnons1_tmp,(3,nsym))
   if (dtset%prepanl/=1.and.&
&   dtset%berryopt/= 4.and.dtset%berryopt/= 6.and.dtset%berryopt/= 7.and.&
&   dtset%berryopt/=14.and.dtset%berryopt/=16.and.dtset%berryopt/=17) then
     call littlegroup_pert(gprimd,idir,indsym,ab_out,ipert,dtset%natom,nsym,nsym1,2,&
&     dtset%symafm,symaf1_tmp,symq,symrec,dtset%symrel,symrl1_tmp,0,dtset%tnons,tnons1_tmp)
   else
     nsym1 = 1
     symaf1_tmp(1) = 1
     symrl1_tmp(:,:,1) = dtset%symrel(:,:,1)
     tnons1_tmp(:,1) = 0_dp
   end if
   ABI_ALLOCATE(indsy1,(4,nsym1,dtset%natom))
   ABI_ALLOCATE(symrc1,(3,3,nsym1))
   ABI_ALLOCATE(symaf1,(nsym1))
   ABI_ALLOCATE(symrl1,(3,3,nsym1))
   ABI_ALLOCATE(tnons1,(3,nsym1))
   symaf1(1:nsym1)=symaf1_tmp(1:nsym1)
   symrl1(:,:,1:nsym1)=symrl1_tmp(:,:,1:nsym1)
   tnons1(:,1:nsym1)=tnons1_tmp(:,1:nsym1)
   ABI_DEALLOCATE(symaf1_tmp)
   ABI_DEALLOCATE(symrl1_tmp)
   ABI_DEALLOCATE(tnons1_tmp)

!  Set up corresponding symmetry data
   ABI_ALLOCATE(irrzon1,(dtset%nfft**(1-1/nsym1),2,(nspden/dtset%nsppol)-3*(nspden/4)))
   ABI_ALLOCATE(phnons1,(2,dtset%nfft**(1-1/nsym1),(nspden/dtset%nsppol)-3*(nspden/4)))
   call setsym(indsy1,irrzon1,1,dtset%natom,dtset%nfft,dtset%ngfft,nspden,dtset%nsppol,&
&   nsym1,phnons1,symaf1,symrc1,symrl1,tnons1,dtset%typat,xred)
   if (psps%usepaw==1) then
!    Allocate/initialize only zarot in pawang1 datastructure
     call pawang_init(pawang1,0,0,0,pawang%l_max-1,0,nsym1,0,1,0,0,0)
     call setsym_ylm(gprimd,pawang1%l_max-1,pawang1%nsym,0,rprimd,symrc1,pawang1%zarot)
   end if

!  Initialize k+q array
   ABI_ALLOCATE(kpq,(3,nkpt))
   if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
     kpq(:,1:nkpt)=dtset%kptns(:,1:nkpt) ! Do not modify, needed for gfortran
   else
     do ikpt=1,nkpt
       kpq(:,ikpt)=dtset%qptn(:)+dtset%kptns(:,ikpt)
     end do
   end if
!  In case wf1 at +q and -q are not related by time inversion symmetry, compute k-q as well for initializations
   if (.not.kramers_deg) then
     ABI_ALLOCATE(kmq,(3,nkpt))
     if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
       kmq(:,1:nkpt)=dtset%kptns(:,1:nkpt) ! Do not modify, needed for gfortran
     else
       do ikpt=1,nkpt
         kmq(:,ikpt)=-dtset%qptn(:)+dtset%kptns(:,ikpt) ! kmq <= k-q
       end do
     end if
   end if

!  Determine the subset of k-points needed in the "reduced Brillouin zone" and initialize other quantities
   ABI_ALLOCATE(indkpt1_tmp,(nkpt))
   ABI_ALLOCATE(wtk_folded,(nkpt))
   ABI_ALLOCATE(bz2ibz_smap, (6, nkpt))
   indkpt1_tmp(:)=0 ; optthm=0
   timrev_pert=timrev
   if(dtset%ieig2rf>0) then
     timrev_pert=0
     call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
&     1,symrc1,timrev_pert,dtset%wtk,wtk_folded, bz2ibz_smap, xmpi_comm_self)
   else
!    For the time being, the time reversal symmetry is not used
!    for ddk, elfd, mgfd perturbations.
     timrev_pert=timrev
     if(ipert==dtset%natom+1.or.ipert==dtset%natom+2.or.&
&      ipert==dtset%natom+10.or.ipert==dtset%natom+11.or. &
&      dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.  &
&      dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17.or.  &
&      ipert==dtset%natom+5.or.dtset%prtfull1wf==1) timrev_pert=0
     timrev_kpt = timrev_pert
!    The time reversal symmetry is not used for the BZ sampling when kptopt=3 or 4
     if (dtset%kptopt==3.or.dtset%kptopt==4) timrev_kpt = 0
     call symkpt(0,gmet,indkpt1_tmp,ab_out,dtset%kptns,nkpt,nkpt_rbz,&
     nsym1,symrc1,timrev_kpt,dtset%wtk,wtk_folded, bz2ibz_smap, xmpi_comm_self)
   end if

   ABI_ALLOCATE(doccde_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(indkpt1,(nkpt_rbz))
   ABI_ALLOCATE(istwfk_rbz,(nkpt_rbz))
   ABI_ALLOCATE(kpq_rbz,(3,nkpt_rbz))
   if (.not.kramers_deg) then
     ABI_ALLOCATE(kmq_rbz,(3,nkpt_rbz))
   end if
   ABI_ALLOCATE(kpt_rbz,(3,nkpt_rbz))
   ABI_ALLOCATE(nband_rbz,(nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(occ_rbz,(dtset%mband*nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(wtk_rbz,(nkpt_rbz))
   indkpt1(:)=indkpt1_tmp(1:nkpt_rbz)
   do ikpt=1,nkpt_rbz
     istwfk_rbz(ikpt)=dtset%istwfk(indkpt1(ikpt))
     kpq_rbz(:,ikpt)=kpq(:,indkpt1(ikpt))
     kpt_rbz(:,ikpt)=dtset%kptns(:,indkpt1(ikpt))
     wtk_rbz(ikpt)=wtk_folded(indkpt1(ikpt))
   end do
   if (.not.kramers_deg) then
     do ikpt=1,nkpt_rbz
       kmq_rbz(:,ikpt)=kmq(:,indkpt1(ikpt))
     end do
   end if
   ABI_DEALLOCATE(indkpt1_tmp)
   ABI_DEALLOCATE(wtk_folded)

!  Transfer occ to occ_rbz and doccde to doccde_rbz :
!  this is a more delicate issue
!  NOTE : this takes into account that indkpt1 is ordered
!  MG: What about using occ(band,kpt,spin) ???
   bdtot_index=0;bdtot1_index=0
   do isppol=1,dtset%nsppol
     ikpt1=1
     do ikpt=1,nkpt
       nband_k=dtset%nband(ikpt+(isppol-1)*nkpt)
!      Must test against ikpt1/=nkpt_rbz+1, before evaluate indkpt1(ikpt1)
       if(ikpt1/=nkpt_rbz+1)then
         if(ikpt==indkpt1(ikpt1))then
           nband_rbz(ikpt1+(isppol-1)*nkpt_rbz)=nband_k
           occ_rbz(1+bdtot1_index:nband_k+bdtot1_index)    = occ(1+bdtot_index:nband_k+bdtot_index)
           doccde_rbz(1+bdtot1_index:nband_k+bdtot1_index) = doccde(1+bdtot_index:nband_k+bdtot_index)
           ikpt1=ikpt1+1
           bdtot1_index=bdtot1_index+nband_k
         end if
       end if
       bdtot_index=bdtot_index+nband_k
     end do
   end do

   _IBM6("IBM6 in dfpt_looppert before getmpw")

!  Compute maximum number of planewaves at k
   call timab(142,1,tsec)
   call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpt_rbz,mpi_enreg,mpw,nkpt_rbz)
   call timab(142,2,tsec)

!  Allocate some k-dependent arrays at k
   ABI_ALLOCATE(kg,(3,mpw*nkpt_rbz))
   ABI_ALLOCATE(npwarr,(nkpt_rbz))
   ABI_ALLOCATE(npwtot,(nkpt_rbz))

!  Determine distribution of k-points/bands over MPI processes
   if (allocated(mpi_enreg%my_kpttab)) then
     ABI_DEALLOCATE(mpi_enreg%my_kpttab)
   end if
   ABI_ALLOCATE(mpi_enreg%my_kpttab,(nkpt_rbz))
   if(xmpi_paral==1) then
     ABI_ALLOCATE(mpi_enreg%proc_distrb,(nkpt_rbz,dtset%mband,dtset%nsppol))
     call distrb2(dtset%mband,nband_rbz,nkpt_rbz,mpi_enreg%nproc_cell,dtset%nsppol,mpi_enreg)
   else
     mpi_enreg%my_kpttab(:)=(/(ii,ii=1,nkpt_rbz)/)
   end if
   my_nkpt_rbz=maxval(mpi_enreg%my_kpttab)
   call initmpi_band(mpi_enreg,nband_rbz,nkpt_rbz,dtset%nsppol)
   mkmem_rbz =my_nkpt_rbz ; mkqmem_rbz=my_nkpt_rbz ; mk1mem_rbz=my_nkpt_rbz
   ABI_UNUSED((/mkmem,mk1mem,mkqmem/))

! given number of reduced kpt, store distribution of bands across procs
   ABI_ALLOCATE(distrbflags,(nkpt_rbz,dtset%mband,dtset%nsppol))
   distrbflags = (mpi_enreg%proc_distrb == mpi_enreg%me_kpt)

   _IBM6("IBM6 before kpgio")

!  Set up the basis sphere of planewaves at k
   call timab(143,1,tsec)
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg,&
&   kpt_rbz,mkmem_rbz,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw,npwarr,npwtot,dtset%nsppol)
   call timab(143,2,tsec)

!  Set up the spherical harmonics (Ylm) at k
   useylmgr=0; option=0 ; nylmgr=0
   if (psps%useylm==1.and. &
&   (ipert==dtset%natom+1.or.ipert==dtset%natom+3.or.ipert==dtset%natom+4.or. &
&   (psps%usepaw==1.and.ipert==dtset%natom+2))) then
     useylmgr=1; option=1 ; nylmgr=3
   else if (psps%useylm==1.and.(ipert==dtset%natom+10.or.ipert==dtset%natom+11)) then
     useylmgr=1; option=2 ; nylmgr=9
   end if
   ABI_ALLOCATE(ylm,(mpw*mkmem_rbz,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr,(mpw*mkmem_rbz,nylmgr,psps%mpsang*psps%mpsang*psps%useylm*useylmgr))
   if (psps%useylm==1) then
     call initylmg(gprimd,kg,kpt_rbz,mkmem_rbz,mpi_enreg,psps%mpsang,mpw,nband_rbz,nkpt_rbz,&
&     npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
   end if

   _IBM6("Before ieig2rf > 0")

!  Set up occupations for this perturbation
   if (dtset%ieig2rf>0) then
     if (.not.allocated(istwfk_pert)) then
       ABI_ALLOCATE(istwfk_pert,(nkpt,3,mpert))
       ABI_ALLOCATE(occ_pert,(dtset%mband*nkpt*dtset%nsppol))
       istwfk_pert(:,:,:)=0 ; occ_pert(:)=zero
     end if
     istwfk_pert(:,idir,ipert)=istwfk_rbz(:)
     occ_pert(:)=occ_rbz(:)
   end if
   if (dtset%efmas>0) then
     if (.not.allocated(istwfk_pert)) then
       ABI_ALLOCATE(istwfk_pert,(nkpt,3,mpert))
       istwfk_pert(:,:,:)=0
     end if
     istwfk_pert(:,idir,ipert)=istwfk_rbz(:)
   end if

!  Print a separator in output file
   write(message,'(3a)')ch10,'--------------------------------------------------------------------------------',ch10
   call wrtout(ab_out,message,'COLL')

!  Initialize band structure datatype at k
   bantot_rbz=sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(eigen0,(bantot_rbz))
   eigen0(:)=zero
   call ebands_init(bantot_rbz,ebands_k,dtset%nelect,doccde_rbz,eigen0,istwfk_rbz,kpt_rbz,&
&   nband_rbz,nkpt_rbz,npwarr,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,occ_rbz,wtk_rbz,&
&   dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
&   dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)
   ABI_DEALLOCATE(eigen0)

!  Initialize header, update it with evolving variables
   gscase=0 ! A GS WF file is read
   call hdr_init(ebands_k,codvsn,dtset,hdr0,pawtab,gscase,psps,wvl%descr,&
     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   call hdr0%update(bantot_rbz,etotal,fermie,&
     residm,rprimd,occ_rbz,pawrhoij_pert,xred,dtset%amu_orig(:,1),&
     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!  Initialize GS wavefunctions at k
   ireadwf0=1; formeig=0 ; ask_accurate=1 ; optorth=0
   mcg=mpw*dtset%nspinor*dtset%mband_mem*mkmem_rbz*dtset%nsppol
   if (one*mpw*dtset%nspinor*dtset%mband_mem*mkmem_rbz*dtset%nsppol > huge(1)) then
     write (message,'(4a, 5(a,i0), 2a)')&
&     "Default integer is not wide enough to store the size of the GS wavefunction array (WF0, mcg).",ch10,&
&     "Action: increase the number of processors. Consider also OpenMP threads.",ch10,&
&     "nspinor: ",dtset%nspinor, "mpw: ",mpw, "mband: ",dtset%mband, "mkmem_rbz: ",&
&     mkmem_rbz, "nsppol: ",dtset%nsppol,ch10,&
&     'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS/LAPACK...) compiled in int64 mode'
     MSG_ERROR(message)
   end if
   ABI_MALLOC_OR_DIE(cg,(2,mcg), ierr)

   ABI_ALLOCATE(eigen0,(dtset%mband*nkpt_rbz*dtset%nsppol))
   call timab(144,1,tsec)

! Initialize the wave function type and read GS WFK
   call wfk_read_my_kptbands(dtfil%fnamewffk, dtset, distrbflags, spacecomm, &
&            formeig, istwfk_rbz, kpt_rbz, nkpt_rbz, npwarr, &
&            cg, eigen=eigen0, occ=occ_rbz)
  
   call timab(144,2,tsec)

   ! Update energies GS energies at k
   call put_eneocc_vect(ebands_k, "eig", eigen0)

!  PAW: compute on-site projections of GS wavefunctions (cprj) (and derivatives) at k
   ncpgr=0
   ABI_DATATYPE_ALLOCATE(cprj,(0,0))
   if (psps%usepaw==1) then
     ncpgr=3 ! Valid for ipert<=natom (phonons), ipert=natom+2 (elec. field)
             ! or for ipert==natom+10,11
     if (ipert==dtset%natom+1) ncpgr=1
     if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) ncpgr=1
     if (usecprj==1) then
!TODO MJV: PAW case also needs porting to mband_mem loops
       mcprj=dtset%nspinor*dtset%mband*mkmem_rbz*dtset%nsppol
       ABI_DATATYPE_DEALLOCATE(cprj)
       ABI_DATATYPE_ALLOCATE(cprj,(dtset%natom,mcprj))
       call pawcprj_alloc(cprj,ncpgr,dimcprj_srt)
       if (ipert<=dtset%natom) then
         choice=2; iorder_cprj=0; idir0=0
       else if (ipert==dtset%natom+1) then
         choice=5; iorder_cprj=0; idir0=idir
       else if (ipert==dtset%natom+2) then
         choice=5; iorder_cprj=0; idir0=0
       else if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
         choice=3; iorder_cprj=0; idir0=istr
       else if (ipert==dtset%natom+10.or.ipert==dtset%natom+11) then
         choice=5; iorder_cprj=0; idir0=0 ! Compute all first derivatives
       else
         choice=1; iorder_cprj=0; idir0=0
       end if
       call ctocprj(atindx,cg,choice,cprj,gmet,gprimd,-1,idir0,iorder_cprj,istwfk_rbz,&
&       kg,kpt_rbz,mcg,mcprj,dtset%mgfft,mkmem_rbz,mpi_enreg,psps%mpsang,mpw,&
&       dtset%natom,nattyp,nband_rbz,dtset%natom,dtset%ngfft,nkpt_rbz,dtset%nloalg,&
&       npwarr,dtset%nspinor,dtset%nsppol,ntypat,dtset%paral_kgb,ph1d,psps,&
&       rmet,dtset%typat,ucvol,dtfil%unpaw,xred,ylm,ylmgr)
     end if
   end if

!  Compute maximum number of planewaves at k+q
!  Will be useful for both GS wfs at k+q and RF wavefunctions
   call timab(143,1,tsec)
   call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kpq_rbz,mpi_enreg,mpw1,nkpt_rbz)
   if (.not.kramers_deg) then
     call getmpw(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kmq_rbz,mpi_enreg,mpw1_mq,nkpt_rbz)
     !number of plane waves at k+q and k-q should be in principle the same to reconstruct rhor1_pq (?)
     !mpw1=max(mpw1,mpw1_tmp)
   else
     mpw1_mq=0
   end if
   call timab(143,2,tsec)

!  Allocate some arrays at k+q
   ABI_ALLOCATE(kg1,(3,mpw1*mk1mem_rbz))
   ABI_ALLOCATE(npwar1,(nkpt_rbz))
   ABI_ALLOCATE(npwtot1,(nkpt_rbz))
!  In case Kramers degeneracy is broken, do the same for k-q
   if (.not.kramers_deg) then
     ABI_ALLOCATE(kg1_mq,(3,mpw1_mq*mk1mem_rbz))
     ABI_ALLOCATE(npwar1_mq,(nkpt_rbz))
     ABI_ALLOCATE(npwtot1_mq,(nkpt_rbz))
   end if

!  Set up the basis sphere of planewaves at k+q
!  Will be useful for both GS wfs at k+q and RF wavefunctions
   call timab(142,1,tsec)
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg1,&
&   kpq_rbz,mk1mem_rbz,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw1,&
&   npwar1,npwtot1,dtset%nsppol)
   if (.not.kramers_deg) then
     call kpgio(ecut_eff,dtset%exchn2n3d,gmet,istwfk_rbz,kg1_mq,&
&     kmq_rbz,mk1mem_rbz,nband_rbz,nkpt_rbz,'PERS',mpi_enreg,mpw1_mq,&
&     npwar1_mq,npwtot1_mq,dtset%nsppol)
   end if
   call timab(142,2,tsec)

!  Set up the spherical harmonics (Ylm) at k+q
   useylmgr1=0; option=0 ; nylmgr1=0
   if (psps%useylm==1.and. &
&   (ipert==dtset%natom+1.or.ipert==dtset%natom+3.or.ipert==dtset%natom+4.or. &
&   (psps%usepaw==1.and.ipert==dtset%natom+2))) then
     useylmgr1=1; option=1; nylmgr1=3
   else if (psps%useylm==1.and.(ipert==dtset%natom+10.or.ipert==dtset%natom+11)) then
     useylmgr1=1; option=2; nylmgr1=9
   end if
   ABI_ALLOCATE(ylm1,(mpw1*mk1mem_rbz,psps%mpsang*psps%mpsang*psps%useylm))
   ABI_ALLOCATE(ylmgr1,(mpw1*mk1mem_rbz,nylmgr1,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))
   if (psps%useylm==1) then
     call initylmg(gprimd,kg1,kpq_rbz,mk1mem_rbz,mpi_enreg,psps%mpsang,mpw1,nband_rbz,nkpt_rbz,&
&     npwar1,dtset%nsppol,option,rprimd,ylm1,ylmgr1)
   end if
!  SPr: do the same for k-q if kramers_deg=.false.

!  Print a separator in output file
   write(message, '(a,a)' )'--------------------------------------------------------------------------------',ch10
   call wrtout(ab_out,message,'COLL')

!  Initialize band structure datatype at k+q
   ABI_ALLOCATE(eigenq,(bantot_rbz))
   eigenq(:)=zero
   call ebands_init(bantot_rbz,ebands_kq,dtset%nelect,doccde_rbz,eigenq,istwfk_rbz,kpq_rbz,&
&   nband_rbz,nkpt_rbz,npwar1,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,occ_rbz,wtk_rbz,&
&   dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
&   dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)
   if (.not.kramers_deg) then
     eigenq(:)=zero
     call ebands_init(bantot_rbz,ebands_kmq,dtset%nelect,doccde_rbz,eigenq,istwfk_rbz,kmq_rbz,&
&     nband_rbz,nkpt_rbz,npwar1_mq,dtset%nsppol,dtset%nspinor,dtset%tphysel,dtset%tsmear,dtset%occopt,occ_rbz,wtk_rbz,&
&     dtset%charge, dtset%kptopt, dtset%kptrlatt_orig, dtset%nshiftk_orig, dtset%shiftk_orig, &
&     dtset%kptrlatt, dtset%nshiftk, dtset%shiftk)
   end if
   ABI_DEALLOCATE(eigenq)

!  Initialize header
   call hdr_init(ebands_kq,codvsn,dtset,hdr,pawtab,pertcase,psps,wvl%descr, &
     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab )
   if (.not.kramers_deg) then
     call hdr_init(ebands_kmq,codvsn,dtset,hdr,pawtab,pertcase,psps,wvl%descr, &
       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab )
   end if

!  Initialize wavefunctions at k+q
!  MG: Here it is possible to avoid the extra reading if the same k mesh can be used.
   ireadwf0=1 ; formeig=0 ; ask_accurate=1 ; optorth=0
   mcgq=mpw1*dtset%nspinor*dtset%mband_mem*mkqmem_rbz*dtset%nsppol
   !SPr: verified until here, add mcgq for -q
   if (one*mpw1*dtset%nspinor*dtset%mband*mkqmem_rbz*dtset%nsppol > huge(1)) then
     write (message,'(4a, 5(a,i0), 2a)')&
&     "Default integer is not wide enough to store the size of the GS wavefunction array (WFKQ, mcgq).",ch10,&
&     "Action: increase the number of processors. Consider also OpenMP threads.",ch10,&
&     "nspinor: ",dtset%nspinor, "mpw1: ",mpw1, "mband: ",dtset%mband, "mkqmem_rbz: ",&
&     mkqmem_rbz, "nsppol: ",dtset%nsppol,ch10,&
&     'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS/LAPACK...) compiled in int64 mode'
     MSG_ERROR(message)
   end if
   ABI_MALLOC_OR_DIE(cgq,(2,mcgq), ierr)

   ABI_ALLOCATE(eigenq,(dtset%mband*nkpt_rbz*dtset%nsppol))
   if (.not.kramers_deg) then
     !ABI_MALLOC_OR_DIE(cg_pq,(2,mcgq), ierr)
     !ABI_ALLOCATE(eigen_pq,(dtset%mband*nkpt_rbz*dtset%nsppol))
     mcgmq=mpw1_mq*dtset%nspinor*dtset%mband_mem*mkqmem_rbz*dtset%nsppol
     ABI_MALLOC_OR_DIE(cg_mq,(2,mcgmq), ierr)

     ABI_ALLOCATE(eigen_mq,(dtset%mband*nkpt_rbz*dtset%nsppol))
   end if

   !if (sum(dtset%qptn(1:3)**2)>=1.d-14) then ! non-zero q
!TODO: for many other q this should be avoidable, in principle all if qptrlatt is a subgrid of kprtlatt
!  Or at very least make a pointer instead of a full copy!!!
   if (dtfil%fnamewffq == dtfil%fnamewffk .and. sum(dtset%qptn(1:3)**2) < 1.d-14) then
     call wrtout(std_out, " qpt is Gamma, psi_k+q initialized from psi_k in memory")
     cgq = cg
print *, ' mpw, mpw1 ', mpw, mpw1
print *, 'cgq ', cgq
     eigenq = eigen0
   else
     call timab(144,1,tsec)
     call wfk_read_my_kptbands(dtfil%fnamewffq, dtset, distrbflags, spacecomm, &
&            formeig, istwfk_rbz, kpq_rbz, nkpt_rbz, npwar1, &
&            cgq, eigen=eigenq, occ=occ_rbz)
     call timab(144,2,tsec)

     if (.not.kramers_deg) then
       !SPr: later "make" a separate WFQ file for "-q"
       call timab(144,1,tsec)
       call wfk_read_my_kptbands(dtfil%fnamewffq, dtset, distrbflags, spacecomm, &
&            formeig, istwfk_rbz, kpq_rbz, nkpt_rbz, npwar1_mq, &
&            cg_mq, eigen=eigen_mq, occ=occ_rbz)
       call timab(144,2,tsec)

     end if
   end if
   ! Update energies GS energies at k + q
   call put_eneocc_vect(ebands_kq, "eig", eigenq)
   if (.not.kramers_deg) then
     call put_eneocc_vect(ebands_kmq, "eig", eigen_mq)
   end if

!  PAW: compute on-site projections of GS wavefunctions (cprjq) (and derivatives) at k+q
   ABI_DATATYPE_ALLOCATE(cprjq,(0,0))
   if (psps%usepaw==1) then
     if (usecprj==1) then
!TODO MJV : PAW
       mcprjq=dtset%nspinor*dtset%mband*mkqmem_rbz*dtset%nsppol
       ABI_DATATYPE_DEALLOCATE(cprjq)
       ABI_DATATYPE_ALLOCATE(cprjq,(dtset%natom,mcprjq))
       call pawcprj_alloc(cprjq,0,dimcprj_srt)
       if (ipert<=dtset%natom.and.(sum(dtset%qptn(1:3)**2)>=1.d-14)) then ! phonons at non-zero q
         choice=1 ; iorder_cprj=0 ; idir0=0
         call ctocprj(atindx,cgq,choice,cprjq,gmet,gprimd,-1,idir0,0,istwfk_rbz,&
&         kg1,kpq_rbz,mcgq,mcprjq,dtset%mgfft,mkqmem_rbz,mpi_enreg,psps%mpsang,mpw1,&
&         dtset%natom,nattyp,nband_rbz,dtset%natom,dtset%ngfft,nkpt_rbz,dtset%nloalg,&
&         npwar1,dtset%nspinor,dtset%nsppol,ntypat,dtset%paral_kgb,ph1d,&
&         psps,rmet,dtset%typat,ucvol,dtfil%unpawq,xred,ylm1,ylmgr1)
       else if (mcprjq>0) then
         call pawcprj_copy(cprj,cprjq)
       end if
     end if
   end if

!  ===== Report on eigenq values
   if (dtset%ieig2rf>0.and.icase==ipert_cnt) then
     eigen0_pert(:) = eigen0(:)
     eigenq_pert(:) = eigenq(:)
     occ_rbz_pert(:) = occ_rbz(:)
   end if
   if (dtset%efmas>0.and.icase==ipert_cnt) then
     eigen0_pert(:) = eigen0(:)
   end if
   !call wrtout(std_out,ch10//' dfpt_looppert: eigenq array',"COLL")
   nkpt_eff=nkpt
   if( (dtset%prtvol==0.or.dtset%prtvol==1.or.dtset%prtvol==2) .and. nkpt>nkpt_max ) nkpt_eff=nkpt_max
   band_index=0
   do isppol=1,dtset%nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       if(ikpt<=nkpt_eff)then
         write(message, '(a,i2,a,i5)' )'  isppol=',isppol,', k point number',ikpt
         call wrtout(std_out,message,'COLL')
         do iband=1,nband_k,4
           write(message, '(a,4es16.6)')'  ',eigenq(iband+band_index:min(iband+3,nband_k)+band_index)
           call wrtout(std_out,message,'COLL')
         end do
       else if(ikpt==nkpt_eff+1)then
         write(message,'(a,a)' )'  respfn : prtvol=0, 1 or 2, stop printing eigenq.',ch10
         call wrtout(std_out,message,'COLL')
       end if
       band_index=band_index+nband_k
     end do
   end do

!  Generate occupation numbers for the reduced BZ at k+q
   ABI_ALLOCATE(docckqde,(dtset%mband*nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(occkq,(dtset%mband*nkpt_rbz*dtset%nsppol))
   if (.not.kramers_deg) then
     ABI_ALLOCATE(docckde_mq,(dtset%mband*nkpt_rbz*dtset%nsppol))
     ABI_ALLOCATE(occk_mq,(dtset%mband*nkpt_rbz*dtset%nsppol))
   end if

   if(0<=dtset%occopt .and. dtset%occopt<=2)then
!    Same occupation numbers at k and k+q (usually, insulating)
     occkq(:)=occ_rbz(:)
     docckqde(:)=zero  ! docckqde is irrelevant in this case
     if(.not.kramers_deg) then
       occk_mq(:)=occ_rbz(:)
       docckde_mq(:)=zero
     end if
   else
!    Metallic occupation numbers
     option=1
     dosdeltae=zero ! the DOS is not computed with option=1
     maxocc=two/(dtset%nspinor*dtset%nsppol)
     call getnel(docckqde,dosdeltae,eigenq,entropy,fermie,maxocc,dtset%mband,&
&     nband_rbz,nelectkq,nkpt_rbz,dtset%nsppol,occkq,dtset%occopt,option,&
&     dtset%tphysel,dtset%tsmear,fake_unit,wtk_rbz)
!    Compare nelect at k and nelelect at k+q
     write(message, '(a,a,a,es16.6,a,es16.6,a)')&
&     ' dfpt_looppert : total number of electrons, from k and k+q',ch10,&
&     '  fully or partially occupied states are',dtset%nelect,' and',nelectkq,'.'
     call wrtout(ab_out,message,'COLL')
     call wrtout(std_out,message,'COLL')
     if (.not.kramers_deg) then
       call getnel(docckde_mq,dosdeltae,eigen_mq,entropy,fermie,maxocc,dtset%mband,&
&       nband_rbz,nelectkq,nkpt_rbz,dtset%nsppol,occk_mq,dtset%occopt,option,&
&       dtset%tphysel,dtset%tsmear,fake_unit,wtk_rbz)
!      Compare nelect at k and nelelect at k-q
       write(message, '(a,a,a,es16.6,a,es16.6,a)')&
&       ' dfpt_looppert : total number of electrons, from k and k-q',ch10,&
&       '  fully or partially occupied states are',dtset%nelect,' and',nelectkq,'.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if
   end if

!  Debug message
   if(dtset%prtvol==-level) call wrtout(std_out,'dfpt_looppert: initialisation of q part done.','COLL')

!  Initialisation of first-order wavefunctions
   write(message,'(3a,i4)')' Initialisation of the first-order wave-functions :',ch10,'  ireadwf=',dtfil%ireadwf
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
   call appdig(pertcase,dtfil%fnamewff1,fiwf1i)
   call appdig(pertcase,dtfil%fnameabo_1wf,fiwf1o)

!  Allocate 1st-order PAW occupancies (rhoij1)
   if (psps%usepaw==1) then
     ABI_DATATYPE_ALLOCATE(pawrhoij1,(my_natom))
     call pawrhoij_nullify(pawrhoij1)
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                          nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij1,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,&
&                        dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,&
&                        comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     if (cplex_rhoij/=hdr%pawrhoij(1)%cplex_rhoij.or.qphase_rhoij/=hdr%pawrhoij(1)%qphase) then
!      Eventually reallocate hdr%pawrhoij
       call pawrhoij_free(hdr%pawrhoij)
       call pawrhoij_alloc(hdr%pawrhoij,cplex_rhoij,nspden_rhoij,dtset%nspinor,dtset%nsppol,&
&                          dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,&
&                          comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
     end if
   else
     ABI_DATATYPE_ALLOCATE(pawrhoij1,(0))
   end if

!  Initialize 1st-order wavefunctions
   formeig=1; ask_accurate=0; optorth=0
! NB: 4 Sept 2013: this was introducing a bug - for ieig2rf ==0 the dim was being set to 1 and passed to dfpt_vtowfk
   dim_eig2rf=0
   if ((dtset%ieig2rf > 0 .and. dtset%ieig2rf/=2) .or. dtset%efmas > 0) then
     dim_eig2rf=1
   end if
   mcg1=mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol
   if (one*mpw1*dtset%nspinor*dtset%mband*mk1mem_rbz*dtset%nsppol > huge(1)) then
     write (message,'(4a, 5(a,i0), 2a)')&
&     "Default integer is not wide enough to store the size of the GS wavefunction array (WFK1, mcg1).",ch10,&
&     "Action: increase the number of processors. Consider also OpenMP threads.",ch10,&
&     "nspinor: ",dtset%nspinor, "mpw1: ",mpw1, "mband: ",dtset%mband, "mk1mem_rbz: ",&
&     mk1mem_rbz, "nsppol: ",dtset%nsppol,ch10,&
&     'Note: Compiling with large int (int64) requires a full software stack (MPI/FFTW/BLAS/LAPACK...) compiled in int64 mode'
     MSG_ERROR(message)
   end if
   ABI_MALLOC_OR_DIE(cg1,(2,mcg1), ierr)
   if (.not.kramers_deg) then
     mcg1mq=mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol
     ABI_MALLOC_OR_DIE(cg1_mq,(2,mcg1mq), ierr)
   end if

   ABI_ALLOCATE(cg1_active,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
   ABI_ALLOCATE(gh1c_set,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
   ABI_ALLOCATE(gh0c1_set,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
   if (.not.kramers_deg) then
     ABI_ALLOCATE(cg1_active_mq,(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
     ABI_ALLOCATE(gh1c_set_mq,(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
     ABI_ALLOCATE(gh0c1_set_mq,(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
   end if
!  XG090606 This is needed in the present 5.8.2 , for portability for the pathscale machine.
!  However, it is due to a bug to be corrected by Paul Boulanger. When the bug will be corrected,
!  this line should be removed.
   if(mk1mem_rbz/=0 .and. dtset%ieig2rf/=0)then
     cg1_active=zero
     gh1c_set=zero
     gh0c1_set=zero
     if (.not.kramers_deg) then
       cg1_active_mq=zero
       gh1c_set_mq=zero
       gh0c1_set_mq=zero
     end if
   end if
   ABI_ALLOCATE(eigen1,(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol))
   ABI_ALLOCATE(resid,(dtset%mband*nkpt_rbz*dtset%nsppol))
   call timab(144,1,tsec)
   if (file_exists(nctk_ncify(fiwf1i)) .or. file_exists(fiwf1i)) then
     call wfk_read_my_kptbands(fiwf1i, dtset, distrbflags, spacecomm, &
&            formeig, istwfk_rbz, kpq_rbz, nkpt_rbz, npwar1, &
&            cg1, eigen=eigen1, occ=occ_rbz)
   else
     cg1 = zero
     eigen1 = zero
!     call inwffil(ask_accurate,cg1,dtset,dtset%ecut,ecut_eff,eigen1,dtset%exchn2n3d,&
!&     formeig,hdr,&
!&     dtfil%ireadwf,istwfk_rbz,kg1,kpq_rbz,dtset%localrdwf,&
!&     dtset%mband,mcg_tmp,mk1mem_rbz,mpi_enreg,mpw1,nband_rbz,dtset%ngfft,nkpt_rbz,npwar1,&
!&     dtset%nsppol,nsym1,occ_rbz,optorth,&
!&     symaf1,symrl1,tnons1,dtfil%unkg1,wff1,wffnow,dtfil%unwff1,&
!&     fiwf1i,wvl)
   end if
   call timab(144,2,tsec)

   if(.not.kramers_deg) then
     ABI_ALLOCATE(eigen1_mq,(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol))
     ABI_ALLOCATE(resid_mq,(dtset%mband*nkpt_rbz*dtset%nsppol))
     !initialize cg1_mq:
     call timab(144,1,tsec)
     if (file_exists(nctk_ncify(fiwf1i)) .or. file_exists(fiwf1i)) then
       call wfk_read_my_kptbands(fiwf1i, dtset, distrbflags, spacecomm, &
&            formeig, istwfk_rbz, kmq_rbz, nkpt_rbz, npwar1_mq, &
&            cg1_mq, eigen=eigen1_mq, occ=occ_rbz)
     else
       cg1_mq = zero
       eigen1_mq = zero
     end if
     call timab(144,2,tsec)
   end if

!  Eventually reytrieve 1st-order PAW occupancies from file header
   if (psps%usepaw==1.and.dtfil%ireadwf/=0) then
     call pawrhoij_copy(hdr%pawrhoij,pawrhoij1,comm_atom=mpi_enreg%comm_atom , &
&     mpi_atmtab=mpi_enreg%my_atmtab)
   end if

!  In case of electric field, or 2nd order perturbation : open the ddk (or dE) wf file(s)
   if ((ipert==dtset%natom+2.and.sum((dtset%qptn(1:3))**2)<1.0d-7.and. &
&   (dtset%berryopt/= 4.and.dtset%berryopt/= 6.and. &
&   dtset%berryopt/= 7.and.dtset%berryopt/=14.and. &
&   dtset%berryopt/=16.and.dtset%berryopt/=17)) .or. &
&   ipert==dtset%natom+10.or.ipert==dtset%natom+11) then
     if (ipert<dtset%natom+10) then
       ! 1st order or berry phase, one direction
       nwffile=1 ! one direction needed
       file_index(1)=idir+dtset%natom*3
       fnamewff(1)=dtfil%fnamewffddk
     else if (ipert==dtset%natom+10) then
       ! 2nd order (k,k)
       if (idir<=3) then ! one direction needed
         nwffile=1
         file_index(1)=idir+dtset%natom*3
         fnamewff(1)=dtfil%fnamewffddk
       else ! two directions needed
         nwffile=2
         file_index(1)=idir1+dtset%natom*3
         file_index(2)=idir2+dtset%natom*3
         fnamewff(1)=dtfil%fnamewffddk
         fnamewff(2)=dtfil%fnamewffddk
       end if
     else if (ipert==dtset%natom+11) then
       ! 2nd order (k,E)
       nwffile=3 ! dk, dE and dkdk
       idir_dkdk = idir
       if(idir_dkdk>6) idir_dkdk = idir_dkdk - 3
       file_index(1)=idir_dkdk +(dtset%natom+6)*3 ! dkdk
       file_index(2)=idir2+(dtset%natom+1)*3 ! defld file (dir2)
       file_index(3)=idir1+dtset%natom*3     ! ddk file (dir1)
       fnamewff(1)=dtfil%fnamewffdkdk
       fnamewff(2)=dtfil%fnamewffdelfd
       fnamewff(3)=dtfil%fnamewffddk
       if (idir>3) then
         nwffile=4
         file_index(4)=idir2+dtset%natom*3   ! ddk file (dir2)
         fnamewff(4)=dtfil%fnamewffddk
       end if
     end if
     do ii=1,nwffile
       call appdig(file_index(ii),fnamewff(ii),fiwfddk)
       ! Checking the existence of data file
       if (.not. file_exists(fiwfddk)) then
         ! Trick needed to run Abinit test suite in netcdf mode.
         if (file_exists(nctk_ncify(fiwfddk))) then
           write(message,"(3a)")"- File: ",trim(fiwfddk)," does not exist but found netcdf file with similar name."
           call wrtout(std_out,message,'COLL')
           fiwfddk = nctk_ncify(fiwfddk)
         end if
         if (.not. file_exists(fiwfddk)) then
           MSG_ERROR('Missing file: '//TRIM(fiwfddk))
         end if
       end if
       write(message,'(2a)')'-dfpt_looppert : read the wavefunctions from file: ',trim(fiwfddk)
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')
!      Note that the unit number for these files is 50,51,52 or 53 (dtfil%unddk=50)
       call wfk_open_read(ddk_f(ii),fiwfddk,formeig1,dtset%iomode,dtfil%unddk+(ii-1),spaceComm)
     end do
   end if

!  Get first-order local potentials and 1st-order core correction density change
!  (do NOT include xccc3d1 in vpsp1 : this will be done in dfpt_scfcv because vpsp1
!  might become spin-polarized)

   n3xccc=0;if(psps%n1xccc/=0)n3xccc=nfftf
   ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
   ABI_ALLOCATE(vpsp1,(cplex*nfftf))

!  PAW: compute Vloc(1) and core(1) together in reciprocal space
!  --------------------------------------------------------------
   if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
     ndir=1
     call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,istr,ipert,&
&     mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,ntypat,&
&     ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&     atmrhor1=xccc3d1,atmvlocr1=vpsp1,optn_in=n3xccc/nfftf,optn2_in=1,vspl=psps%vlspl)
   else

!    Norm-conserving psp: compute Vloc(1) in reciprocal sp. and core(1) in real sp.
!    ------------------------------------------------------------------------------

     if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) then
!      Section for strain perturbation
       call vlocalstr(gmet,gprimd,gsqcut,istr,mgfftf,mpi_enreg,&
&       psps%mqgrid_vl,dtset%natom,nattyp,nfftf,ngfftf,ntypat,ph1df,psps%qgrid_vl,&
&       ucvol,psps%vlspl,vpsp1)
     else
       call dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir,ipert,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
&       nattyp,nfftf,ngfftf,ntypat,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,psps%qgrid_vl,&
&       dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
       !SPr: need vpsp1 for -q as well, but for magnetic field it's zero, to be done later
     end if

     if(psps%n1xccc/=0)then
       call dfpt_mkcore(cplex,idir,ipert,dtset%natom,ntypat,ngfftf(1),psps%n1xccc,&
&       ngfftf(2),ngfftf(3),dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
       !SPr: same here, need xccc3d1 for -q as well for phonon pert.. to be done later
     end if ! psps%n1xccc/=0
   end if ! usepaw

   eigen1(:)=zero; resid(:)=zero
   if(.not.kramers_deg) then
     eigen1_mq(:)=zero
     resid_mq(:)=zero
   end if
!  Get starting charge density and Hartree + xc potential
   ABI_ALLOCATE(rhor1,(cplex*nfftf,nspden))
   ABI_ALLOCATE(rhog1,(2,nfftf))

   if(.not.kramers_deg) then
   !Case when first order spinors at both +q and -q are not related by symmetry (time and/or space inversion)
     ABI_ALLOCATE(rhor1_pq,(cplex*nfftf,nspden))
     ABI_ALLOCATE(rhog1_pq,(2,nfftf))

     ABI_ALLOCATE(rhor1_mq,(cplex*nfftf,nspden))
     ABI_ALLOCATE(rhog1_mq,(2,nfftf))
   end if

!  can we get this set of gkk matrices from previously calculated rhog1 through a non-scf calculation?
   found_eq_gkk=.false.
   if (dtset%prepgkk == 1 .and. ipert <= dtset%natom) then
!    NOTE: this does not take into account combinations e.g. x+y -> z
!    if rhor1 add linearly this could be done...
     do icase_eq = 1, icase-1
       idir_eq = mod(icase_eq,3)
       if (idir_eq==0) idir_eq=3
       ipert_eq = ( (icase_eq-idir_eq) / 3 + 1)

! find sym which links old perturbation to present one
       do isym=1, nsym
         ! check that isym preserves qpt to begin with
         if (symq(4,1,isym) /= 1 .or. &
&         symq(1,1,isym) /= 0 .or. &
&         symq(2,1,isym) /= 0 .or. &
&         symq(3,1,isym) /= 0      ) cycle

         eq_symop = dtset%symrel(:,:,isym)
         if (indsym(4,isym,ipert) == ipert_eq .and. &
&         abs(eq_symop(idir,idir_eq)) == 1 .and. &
&         sum(abs(eq_symop(:,idir_eq))) == 1) then
           found_eq_gkk = .true.
           exit
         end if
       end do ! isym
       if (found_eq_gkk) exit
     end do ! icase_eq
   end if ! check for prepgkk with symmetric pert

   if (found_eq_gkk) then
     write (message, '(a,l6,i6,2a,3i6,2a,3i6,a)')  &
&     ' found_eq_gkk,isym = ', found_eq_gkk, isym, ch10, &
&     ' idir,  ipert,  icase   =  ', idir, ipert, icase, ch10, &
&     ' idireq iperteq icaseeq =  ', idir_eq, ipert_eq, icase_eq, ch10
     call wrtout(std_out,message,'COLL')
!
!    Make density for present perturbation, which is symmetric of 1 or more previous perturbations:
!    rotate 1DEN arrays with symrel(isym) to produce rhog1_eq rhor1_eq
!
     if (dtset%use_nonscf_gkk == 1) then
       call rotate_rho(cplex, timrev_pert, mpi_enreg, nfftf, ngfftf, nspden, &
&       rhor1_save(:,:,icase_eq), rhog1, rhor1, eq_symop, dtset%tnons(:,isym))

       rhor1 = rhor1 * eq_symop(idir,idir_eq)

! TODO: rotate rhoij in PAW case

       blkflg_save = blkflg
       dtset_tmp%iscf = -2
       iscf_mod = -2
       dtset_tmp%nstep = 1
       dtset_tmp%nline = 1
       if (abs(dtset_tmp%tolwfr) < 1.e-24) dtset_tmp%tolwfr = 1.e-24
       dtset_tmp%toldfe = zero
       dtset_tmp%toldff = zero
       dtset_tmp%tolrff = zero
       dtset_tmp%tolvrs = zero
       write (message, '(a,i6,a)') ' NOTE: doing GKK calculation for icase ', icase, ' with non-SCF calculation'
       call wrtout(std_out,message,'COLL')
       !call wrtout(ab_out,message,'COLL') ! decomment and update output files

     else ! do not use non-scf shortcut, but save rotated 1DEN for comparison
! saves the rotated rho, for later comparison with the full SCF rhor1: comment lines below for iscf = -2
       call rotate_rho(cplex, timrev_pert, mpi_enreg, nfftf, ngfftf, nspden, &
&       rhor1_save(:,:,icase_eq), rhog1, rhor1_save(:,:,icase), eq_symop, &
&       dtset%tnons(:,isym))
       rhor1_save(:,:,icase) = rhor1_save(:,:,icase) * eq_symop(idir,idir_eq)

     end if ! force non scf calculation of other gkk, or not
   end if ! found an equiv perturbation for the gkk

   if ( (dtfil%ireadwf==0 .and. iscf_mod/=-4 .and. dtset%get1den==0 .and. dtset%ird1den==0) &
   .or. (iscf_mod== -3 .and. ipert/=dtset%natom+11) ) then
!    NOTE : For ipert==natom+11, we want to read the 1st order density from a previous calculation
     rhor1(:,:)=zero ; rhog1(:,:)=zero
!    PAW: rhoij have been set to zero in call to pawrhoij_alloc above

     init_rhor1 = ((ipert>=1 .and. ipert<=dtset%natom).or.ipert==dtset%natom+5)
     ! This section is needed in order to maintain the old behavior and pass the automatic tests
     if (psps%usepaw == 0) then
       init_rhor1 = init_rhor1 .and. all(psps%nctab(:ntypat)%has_tvale)
     else
       init_rhor1 = .False.
     end if

     if (init_rhor1) then

       if(ipert/=dtset%natom+5) then

         ! Initialize rhor1 and rhog1 from the derivative of atomic densities/gaussians.
         ndir = 1; optn2 = 3
         if (psps%usepaw==1) then
           ! FIXME: Here there's a bug because has_tvale == 0 or 1 instead of 2
           if (all(pawtab(:ntypat)%has_tvale/=0)) optn2=2
         else if (psps%usepaw==0) then
           if (all(psps%nctab(:ntypat)%has_tvale)) optn2=2
         end if

         if (optn2 == 3) then
           call wrtout(std_out," Initializing rhor1 from atom-centered gaussians", "COLL")
           ABI_ALLOCATE(gauss,(2,ntypat))
           call atom_gauss(ntypat, dtset%densty, psps%ziontypat, psps%znucltypat, gauss)

           call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,idir,ipert,&
           mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,ntypat,&
           ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
           atmrhor1=rhor1,optn_in=1,optn2_in=3,gauss=gauss)

           ABI_FREE(gauss)
         else
           call wrtout(std_out," Initializing rhor1 from valence densities taken from pseudopotential files", "COLL")
           call dfpt_atm2fft(atindx,cplex,gmet,gprimd,gsqcut,idir,ipert,&
           mgfftf,psps%mqgrid_vl,dtset%natom,ndir,nfftf,ngfftf,ntypat,&
           ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
           atmrhor1=rhor1,optn_in=1,optn2_in=2)
         end if

       else
         ! Magnetic field perturbation
         call wrtout(std_out," Initializing rhor1 guess based on the ground state XC magnetic field", "COLL")

         call dfpt_init_mag1(ipert,idir,rhor1,rhor,cplex,nfftf,nspden,vxc,kxc,nkxc)

         if(.not.kramers_deg) then
           rhor1_pq=rhor1
           call dfpt_init_mag1(ipert,idir,rhor1_mq,rhor,cplex,nfftf,nspden,vxc,kxc,nkxc)
         end if
       end if

       call fourdp(cplex,rhog1,rhor1,-1,mpi_enreg,nfftf,1,ngfftf,0)
       if (.not.kramers_deg) then
         !call fourdp(cplex,rhog1_pq,rhor1_pq,-1,mpi_enreg,nfftf,1,ngfftf,0)
         rhog1_pq=rhog1
         call fourdp(cplex,rhog1_mq,rhor1_mq,-1,mpi_enreg,nfftf,1,ngfftf,0)
       end if
     end if

   else
     ! rhor1 not being forced to 0.0
     if(iscf_mod>0) then
!      cplex=2 gets the complex density, =1 only real part
       if (psps%usepaw==1) then
!        Be careful: in PAW, rho does not include the 1st-order compensation density (to be added in dfpt_scfcv.F90) !
         ABI_ALLOCATE(rho1wfg,(2,dtset%nfft))
         ABI_ALLOCATE(rho1wfr,(dtset%nfft,nspden))
         call dfpt_mkrho(cg,cg1,cplex,gprimd,irrzon1,istwfk_rbz,&
           kg,kg1,dtset%mband,dtset%mgfft,mkmem_rbz,mk1mem_rbz,mpi_enreg,mpw,mpw1,nband_rbz,&
           dtset%nfft,dtset%ngfft,nkpt_rbz,npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,nsym1,&
           occ_rbz,phnons1,rho1wfg,rho1wfr,rprimd,symaf1,symrl1,ucvol,wtk_rbz)
         call transgrid(cplex,mpi_enreg,nspden,+1,1,1,dtset%paral_kgb,pawfgr,rho1wfg,rhog1,rho1wfr,rhor1)
         ABI_DEALLOCATE(rho1wfg)
         ABI_DEALLOCATE(rho1wfr)
       else
         !SPr: need to modify dfpt_mkrho to taken into account q,-q and set proper formulas when +q and -q spinors are related
         call dfpt_mkrho(cg,cg1,cplex,gprimd,irrzon1,istwfk_rbz,&
           kg,kg1,dtset%mband,dtset%mgfft,mkmem_rbz,mk1mem_rbz,mpi_enreg,mpw,mpw1,nband_rbz,&
           dtset%nfft,dtset%ngfft,nkpt_rbz,npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,nsym1,&
           occ_rbz,phnons1,rhog1,rhor1,rprimd,symaf1,symrl1,ucvol,wtk_rbz)
       end if

     else if (.not. found_eq_gkk) then
       ! negative iscf_mod and no symmetric rotation of rhor1
       ! Read rho1(r) from a disk file and broadcast data.
       rdwr=1;rdwrpaw=psps%usepaw;if(dtfil%ireadwf/=0) rdwrpaw=0
       if (ipert/=dtset%natom+11) then
         call appdig(pertcase,dtfil%fildens1in,fiden1i)
       else
         ! For ipert==natom+11, we want to read the 1st order density from a previous calculation
         call appdig(idir2+(dtset%natom+1)*3,dtfil%fildens1in,fiden1i)
       end if
!       call appdig(pertcase,dtfil%fildens1in,fiden1i)
       call read_rhor(fiden1i, cplex, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rhor1, &
       hdr_den, pawrhoij1, spaceComm, check_hdr=hdr)
       etotal = hdr_den%etot; call hdr_den%free()

!      Compute up+down rho1(G) by fft
       ABI_ALLOCATE(work,(cplex*nfftf))
       work(:)=rhor1(:,1)
       call fourdp(cplex,rhog1,work,-1,mpi_enreg,nfftf,1,ngfftf,0)
       ABI_DEALLOCATE(work)
     end if ! rhor1 generated or read in from file

   end if ! rhor1 set to 0 or read in from file

!  Check whether exiting was required by the user.
!  If found then do not start minimization steps
   openexit=1 ; if(dtset%chkexit==0) openexit=0
   call exit_check(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg%comm_cell,openexit)
!  If immediate exit, and wavefunctions were not read, must zero eigenvalues
   if (iexit/=0) eigen1(:)=zero
   if (iexit/=0.and.(.not.kramers_deg)) eigen1_mq(:)=zero

   if (iexit==0) then
     _IBM6("before dfpt_scfcv")

!    Main calculation: get 1st-order wavefunctions from Sternheimer equation (SCF cycle)
!    if ipert==natom+10 or natom+11 : get 2nd-order wavefunctions
     if (kramers_deg) then
       call dfpt_scfcv(atindx,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&
&       dielt,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset_tmp,&
&       d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&       ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&
&       enl0,enl1,eovl1,epaw1,etotal,evdw,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,&
&       indkpt1,indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&       kg,kg1,kpt_rbz,kxc,mgfftf,mkmem_rbz,mkqmem_rbz,mk1mem_rbz,&
&       mpert,mpi_enreg,mpw,mpw1,mpw1_mq,my_natom,&
&       nattyp,nband_rbz,ncpgr,nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,&
&       npwarr,npwar1,nspden,&
&       nsym1,n3xccc,occkq,occ_rbz,&
&       paw_an_pert,paw_ij_pert,pawang,pawang1,pawfgr,pawfgrtab_pert,pawrad,pawrhoij_pert,pawrhoij1,pawtab,&
&       pertcase,phnons1,ph1d,ph1df,prtbbb,psps,&
&       dtset%qptn,resid,residm,rhog,rhog1,&
&       rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&
&       usecprj,useylmgr,useylmgr1,ddk_f,vpsp1,vtrial,vxc,&
&       wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1,zeff,dfpt_scfcv_retcode,&
&       kramers_deg)
     else
       call dfpt_scfcv(atindx,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&
&       dielt,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset_tmp,&
&       d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&       ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&
&       enl0,enl1,eovl1,epaw1,etotal,evdw,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,&
&       indkpt1,indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&       kg,kg1,kpt_rbz,kxc,mgfftf,mkmem_rbz,mkqmem_rbz,mk1mem_rbz,&
&       mpert,mpi_enreg,mpw,mpw1,mpw1_mq,my_natom,&
&       nattyp,nband_rbz,ncpgr,nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,&
&       npwarr,npwar1,nspden,&
&       nsym1,n3xccc,occkq,occ_rbz,&
&       paw_an_pert,paw_ij_pert,pawang,pawang1,pawfgr,pawfgrtab_pert,pawrad,pawrhoij_pert,pawrhoij1,pawtab,&
&       pertcase,phnons1,ph1d,ph1df,prtbbb,psps,&
&       dtset%qptn,resid,residm,rhog,rhog1,&
&       rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&
&       usecprj,useylmgr,useylmgr1,ddk_f,vpsp1,vtrial,vxc,&
&       wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1,zeff,dfpt_scfcv_retcode,&
&       kramers_deg,&
&       cg_mq=cg_mq,cg1_mq=cg1_mq,cg1_active_mq=cg1_active_mq,docckde_mq=docckde_mq,eigen_mq=eigen_mq,&
&       eigen1_mq=eigen1_mq,gh0c1_set_mq=gh0c1_set_mq,gh1c_set_mq=gh1c_set_mq,&
&       kg1_mq=kg1_mq,npwar1_mq=npwar1_mq,occk_mq=occk_mq,resid_mq=resid_mq,residm_mq=residm_mq,&
&       rhog1_pq=rhog1_pq,rhog1_mq=rhog1_mq,rhor1_pq=rhor1_pq,rhor1_mq=rhor1_mq)
     end if

     _IBM6("after dfpt_scfcv")

!    2nd-order eigenvalues stuff
     if (dtset%ieig2rf>0) then
       if (first_entry) then
         nullify(eigen1_pert)
         first_entry = .false.
       end if
       if (.not.associated(eigen1_pert)) then
         ABI_ALLOCATE(eigen1_pert,(2*dtset%mband**2*nkpt*dtset%nsppol,3,mpert))
         ABI_MALLOC_OR_DIE(cg1_pert,(2,mpw1*nspinor*dtset%mband_mem*mk1mem_rbz*nsppol*dim_eig2rf,3,mpert),ierr)
         ABI_ALLOCATE(gh0c1_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf,3,mpert))
         ABI_ALLOCATE(gh1c_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf,3,mpert))
         ABI_ALLOCATE(kpt_rbz_pert,(3,nkpt_rbz))
         ABI_ALLOCATE(npwarr_pert,(nkpt_rbz,mpert))
         ABI_ALLOCATE(npwar1_pert,(nkpt_rbz,mpert))
         ABI_ALLOCATE(npwtot_pert,(nkpt_rbz,mpert))
         eigen1_pert(:,:,:) = zero
         cg1_pert(:,:,:,:) = zero
         gh0c1_pert(:,:,:,:) = zero
         gh1c_pert(:,:,:,:) = zero
         npwar1_pert (:,:) = 0
         npwarr_pert (:,:) = 0
         kpt_rbz_pert = kpt_rbz
       end if
       clflg(idir,ipert)=1
       eigen1_pert(1:2*dtset%mband**2*nkpt_rbz*dtset%nsppol,idir,ipert) = eigen1(:)
       if(dtset%ieig2rf==1.or.dtset%ieig2rf==3.or.dtset%ieig2rf==4.or.dtset%ieig2rf==5) then
         cg1_pert(:,:,idir,ipert)=cg1_active(:,:)
         gh0c1_pert(:,:,idir,ipert)=gh0c1_set(:,:)
         gh1c_pert(:,:,idir,ipert)=gh1c_set(:,:)
       end if
       npwarr_pert(:,ipert)=npwarr(:)
       npwar1_pert(:,ipert)=npwar1(:)
       npwtot_pert(:,ipert)=npwtot(:)
     end if
!    2nd-order eigenvalues stuff for EFMAS
     if (dtset%efmas>0) then
       if (first_entry) then
         nullify(eigen1_pert)
         first_entry = .false.
       end if
       if (.not.associated(eigen1_pert)) then
         ABI_ALLOCATE(eigen1_pert,(2*dtset%mband**2*nkpt*dtset%nsppol,3,mpert))
         ABI_ALLOCATE(cg1_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf,3,mpert))
         ABI_ALLOCATE(gh0c1_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf,3,mpert))
         ABI_ALLOCATE(gh1c_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf,3,mpert))
         ABI_ALLOCATE(kpt_rbz_pert,(3,nkpt_rbz))
         ABI_ALLOCATE(npwarr_pert,(nkpt_rbz,mpert))
         eigen1_pert(:,:,:) = zero
         cg1_pert(:,:,:,:) = zero
         gh0c1_pert(:,:,:,:) = zero
         gh1c_pert(:,:,:,:) = zero
         npwarr_pert (:,:) = 0
         kpt_rbz_pert = kpt_rbz
         ABI_ALLOCATE(cg0_pert,(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem_rbz*dtset%nsppol*dim_eig2rf))
         cg0_pert = cg
       end if
       eigen1_pert(1:2*dtset%mband**2*nkpt_rbz*dtset%nsppol,idir,ipert) = eigen1(:)
       cg1_pert(:,:,idir,ipert)=cg1_active(:,:)
       gh0c1_pert(:,:,idir,ipert)=gh0c1_set(:,:)
       gh1c_pert(:,:,idir,ipert)=gh1c_set(:,:)
       npwarr_pert(:,ipert)=npwarr(:)
     end if
     ABI_DEALLOCATE(gh1c_set)
     ABI_DEALLOCATE(gh0c1_set)
     ABI_DEALLOCATE(cg1_active)

!deallocate bit arrays for case without Kramers' degeneracy
     if(.not.kramers_deg) then
       ABI_DEALLOCATE(gh1c_set_mq)
       ABI_DEALLOCATE(gh0c1_set_mq)
       ABI_DEALLOCATE(cg1_active_mq)
     end if

   end if ! End of the check of hasty exit

   call timab(146,1,tsec)

!  Print out message at the end of the iterations
   write(message, '(80a,a,a,a,a)' ) ('=',ii=1,80),ch10,ch10,&
&   ' ----iterations are completed or convergence reached----',ch10
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

!  Print _gkk file for this perturbation
   if (dtset%prtgkk == 1) then
     call appdig(3*(ipert-1)+idir,dtfil%fnameabo_gkk,gkkfilnam)
     nmatel = dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol
     ABI_ALLOCATE(phasecg, (2, nmatel))
!     call getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg, nkpt_rbz, npwarr, npwar1, phasecg)
     phasecg(1,:) = one
     phasecg(2,:) = zero
! NB: phasecg not actually used in outgkk for the moment (2013/08/15)
     call outgkk(bantot_rbz, nmatel,gkkfilnam,eigen0,eigen1,hdr0,hdr,mpi_enreg,phasecg)
     ABI_DEALLOCATE(phasecg)

#ifdef HAVE_NETCDF
     ! Reshape eigen1 into gkk for netCDF output
     ABI_MALLOC_OR_DIE(gkk,(2*dtset%mband*dtset%nsppol,dtset%nkpt,1,1,dtset%mband), ierr)
     gkk(:,:,:,:,:) = zero
     mband = dtset%mband
     band_index = 0
     band2tot_index = 0
     do isppol=1,dtset%nsppol
       do ikpt =1,nkpt_rbz
         do iband=1,dtset%mband
           do jband=1,dtset%mband
             eig1_r = eigen1(2*jband-1+(iband-1)*2*mband+band2tot_index)
             eig1_i = eigen1(2*jband+(iband-1)*2*mband+band2tot_index)
             gkk(2*iband-1+2*band_index,ikpt,1,1,jband) = &
&             gkk(2*iband-1+2*band_index,ikpt,1,1,jband) + eig1_r
             gkk(2*iband+2*band_index,ikpt,1,1,jband) = &
&             gkk(2*iband+2*band_index,ikpt,1,1,jband) + eig1_i
           end do !jband
         end do !iband
         band2tot_index = band2tot_index + 2*mband**2
       end do !ikpt
       band_index = band_index + mband
     end do !isppol

     ! Initialize ggk_ebands to write in the GKK.nc file
     ! MG FIXME: Here there's a bug because eigen0 is dimensioned with nkpt_rbz i.e. IBZ(q)
     ! but the ebands_t object is constructed with dimensions taken from hdr0 i.e. the IBZ(q=0).
     bantot= dtset%mband*dtset%nkpt*dtset%nsppol
     call ebands_init(bantot,gkk_ebands,dtset%nelect,doccde,eigen0,hdr0%istwfk,hdr0%kptns,&
&     hdr0%nband, hdr0%nkpt,hdr0%npwarr,hdr0%nsppol,hdr0%nspinor,&
&     hdr0%tphysel,hdr0%tsmear,hdr0%occopt,hdr0%occ,hdr0%wtk,&
&     hdr0%charge, hdr0%kptopt, hdr0%kptrlatt_orig, hdr0%nshiftk_orig, hdr0%shiftk_orig, &
&     hdr0%kptrlatt, hdr0%nshiftk, hdr0%shiftk)

     ! Init a gkk_t object
     call gkk_init(gkk,gkk2d,dtset%mband,dtset%nsppol,nkpt_rbz,1,1)

     ! Write the netCDF file.
     fname = strcat(gkkfilnam,".nc")
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating GKK file")
     NCF_CHECK(crystal%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(gkk_ebands, ncid))
     call gkk_ncwrite(gkk2d,dtset%qptn(:),dtset%wtq, ncid)
     NCF_CHECK(nf90_close(ncid))

     ! Free memory
     ABI_DEALLOCATE(gkk)
     call gkk_free(gkk2d)
     call ebands_free(gkk_ebands)
#endif
   end if

   if (dtset%prepgkk == 1 .and. found_eq_gkk) then
     if (dtset%use_nonscf_gkk == 1) then
!      Restore old values of SCF cycle parameters
       iscf_mod = iscf_mod_save
       dtset_tmp%iscf = iscf_mod_save
       dtset_tmp%nstep = nstep_save
       dtset_tmp%nline = nline_save
       dtset_tmp%tolwfr = tolwfr_save
       dtset_tmp%toldfe = toldfe_save
       dtset_tmp%toldff = toldff_save
       dtset_tmp%tolrff = tolrff_save
       dtset_tmp%tolvrs = tolvrs_save
       blkflg = blkflg_save ! this ensures we do not use the (unconverged) 2DTE from this non scf run
!      Save density for present perturbation, for future use in symmetric perturbations
       rhor1_save(:,:,icase) = rhor1
     else
       write (message, '(a,3E20.10)') 'norm diff = ', sum(abs(rhor1_save(:,:,icase) - rhor1)), &
&       sum(abs(rhor1)), sum(abs(rhor1_save(:,:,icase) - rhor1))/sum(abs(rhor1))
       call wrtout(std_out,  message,'COLL')
     end if
   end if

   ! Write wavefunctions file only if convergence was not achieved.
   write_1wfk = .True.
   if (dtset%prtwf==-1 .and. dfpt_scfcv_retcode == 0) then
     write_1wfk = .False.
     call wrtout(ab_out," dfpt_looppert: DFPT cycle converged with prtwf=-1. Will skip output of the 1st-order WFK file.","COLL")
   end if

   if (write_1wfk) then
     ! Output 1st-order wavefunctions in file
     call outwf(cg1,dtset,psps,eigen1,fiwf1o,hdr,kg1,kpt_rbz,&
&     dtset%mband,mcg1,mk1mem_rbz,mpi_enreg,mpw1,dtset%natom,nband_rbz,&
&     nkpt_rbz,npwar1,dtset%nsppol,&
&     occ_rbz,resid,response,dtfil%unwff2,wvl%wfs,wvl%descr)
   end if

#ifdef HAVE_NETCDF
  ! Output DDK file in netcdf format.
   if (me == master .and. ipert == dtset%natom + 1) then
     fname = strcat(dtfil%filnam_ds(4), "_EVK.nc")
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EVK.nc file")
    ! Have to build hdr on k-grid with info about perturbation.
     call hdr_copy(hdr0, hdr_tmp)
     hdr_tmp%kptopt = dtset%kptopt
     hdr_tmp%pertcase = pertcase
     hdr_tmp%qptn = dtset%qptn(1:3)
     NCF_CHECK(hdr_tmp%ncwrite(ncid, 43, nc_define=.True.))
     call hdr_tmp%free()
     NCF_CHECK(crystal%ncwrite(ncid))
     NCF_CHECK(ebands_ncwrite(ebands_k, ncid))
     ncerr = nctk_def_arrays(ncid, [ &
     nctkarr_t('h1_matrix_elements', "dp", &
       "two, max_number_of_states, max_number_of_states, number_of_kpoints, number_of_spins")], defmode=.True.)
     NCF_CHECK(ncerr)
     NCF_CHECK(nctk_set_datamode(ncid))
     ncerr = nf90_put_var(ncid, nctk_idname(ncid, "h1_matrix_elements"), eigen1, &
       count=[2, dtset%mband, dtset%mband, nkpt_rbz, dtset%nsppol])
     NCF_CHECK(ncerr)
     NCF_CHECK(nf90_close(ncid))
   end if
#endif


!  If the perturbation is d/dk, evaluate the f-sum rule.
   if (ipert==dtset%natom+1 )then
!    Note : the factor of two is related to the difference
!    between Taylor expansion and perturbation expansion
!    Note : this expression should be modified for ecutsm.
!    Indeed, the present one will NOT tend to 1.0_dp.
     ek2=gmet(idir,idir)*(two_pi**2)*2.0_dp*dtset%nelect
     fsum=-ek1/ek2
     if(dtset%ecutsm<tol6)then
       write(message, '(a,es20.10,a,a,es20.10)' ) &
&       ' dfpt_looppert : ek2=',ek2,ch10,&
&       '          f-sum rule ratio=',fsum
     else
       write(message, '(a,es20.10,a,a,es20.10,a)' ) &
&       ' dfpt_looppert : ek2=',ek2,ch10,&
&       '          f-sum rule ratio=',fsum,' (note : ecutsm/=0)'
     end if
     call wrtout(std_out,message,'COLL')
     call wrtout(ab_out,message,'COLL')
!    Write the diagonal elements of the dH/dk operator, after averaging over degenerate states
     ABI_ALLOCATE(eigen1_mean,(dtset%mband*nkpt_rbz*dtset%nsppol))
     call eigen_meandege(eigen0,eigen1,eigen1_mean,dtset%mband,nband_rbz,nkpt_rbz,dtset%nsppol,1)
     option=4
     if (me == master) then
       call prteigrs(eigen1_mean,dtset%enunit,fermie,dtfil%fnametmp_1wf1_eig,ab_out,iscf_mod,kpt_rbz,dtset%kptopt,&
&       dtset%mband,nband_rbz,nkpt_rbz,dtset%nnsclo,dtset%nsppol,occ_rbz,dtset%occopt,&
&       option,dtset%prteig,dtset%prtvol,resid,tolwfr,vxcavg,wtk_rbz)
     end if
     ABI_DEALLOCATE(eigen1_mean)
   end if

!  Print the energies
   if (dtset%nline/=0 .or. dtset%nstep/=0)then
     call dfpt_prtene(dtset%berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&     ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,eovl1,epaw1,evdw,exc1,ab_out,&
&     ipert,dtset%natom,psps%usepaw,usevdw)
   end if

   if(mpi_enreg%paral_pert==1) then
     if (ipert_me < npert_me -1) then
       call hdr0%free()
     else
       eigen0_copy(1:dtset%mband*nkpt_rbz*dtset%nsppol) = eigen0
     end if
     ipert_me = ipert_me +1
   else
     if (icase == ipert_cnt) then
       eigen0_copy(1:dtset%mband*nkpt_rbz*dtset%nsppol) = eigen0
     else
       call hdr0%free()
     end if
   end if

!Release the temporary arrays (for k, k+q and 1st-order)
   ABI_DEALLOCATE(cg)
   ABI_DEALLOCATE(cgq)
   ABI_DEALLOCATE(cg1)
   ABI_DEALLOCATE(docckqde)
   if(.not.kramers_deg) then
     ABI_DEALLOCATE(cg_mq)
     ABI_DEALLOCATE(cg1_mq)
     ABI_DEALLOCATE(docckde_mq)
   end if
   ABI_DEALLOCATE(doccde_rbz)
   ABI_DEALLOCATE(eigen0)
   ABI_DEALLOCATE(eigenq)
   ABI_DEALLOCATE(eigen1)
   ABI_DEALLOCATE(kpq)
   if(.not.kramers_deg) then
     ABI_DEALLOCATE(eigen_mq)
     ABI_DEALLOCATE(eigen1_mq)
     ABI_DEALLOCATE(kmq)
   end if
   ABI_DEALLOCATE(indkpt1)
   ABI_DEALLOCATE(indsy1)
   ABI_DEALLOCATE(istwfk_rbz)
   ABI_DEALLOCATE(irrzon1)
   ABI_DEALLOCATE(kg)
   ABI_DEALLOCATE(kg1)
   ABI_DEALLOCATE(kpq_rbz)
   if(.not.kramers_deg) then
     ABI_DEALLOCATE(kg1_mq)
     ABI_DEALLOCATE(kmq_rbz)
   end if
   ABI_DEALLOCATE(kpt_rbz)
   ABI_DEALLOCATE(nband_rbz)
   ABI_DEALLOCATE(npwarr)
   ABI_DEALLOCATE(npwar1)
   ABI_DEALLOCATE(npwtot)
   ABI_DEALLOCATE(npwtot1)
   ABI_DEALLOCATE(occkq)
   ABI_DEALLOCATE(occ_rbz)
   ABI_DEALLOCATE(phnons1)
   ABI_DEALLOCATE(resid)
   ABI_DEALLOCATE(rhog1)
   ABI_DEALLOCATE(rhor1)
   ABI_DEALLOCATE(symaf1)
   ABI_DEALLOCATE(symrc1)
   ABI_DEALLOCATE(symrl1)
   ABI_DEALLOCATE(tnons1)
   ABI_DEALLOCATE(wtk_rbz)
   ABI_DEALLOCATE(xccc3d1)
   ABI_DEALLOCATE(vpsp1)
   ABI_DEALLOCATE(ylm)
   ABI_DEALLOCATE(ylm1)
   ABI_DEALLOCATE(ylmgr)
   ABI_DEALLOCATE(ylmgr1)
   if(.not.kramers_deg) then
     ABI_DEALLOCATE(npwar1_mq)
     ABI_DEALLOCATE(npwtot1_mq)
     ABI_DEALLOCATE(occk_mq)
     ABI_DEALLOCATE(resid_mq)
     ABI_DEALLOCATE(rhor1_pq)
     ABI_DEALLOCATE(rhor1_mq)
     ABI_DEALLOCATE(rhog1_pq)
     ABI_DEALLOCATE(rhog1_mq)
   end if
   if (psps%usepaw==1) then
     call pawang_free(pawang1)
     call pawrhoij_free(pawrhoij1)
     if (usecprj==1) then
       call pawcprj_free(cprj)
       call pawcprj_free(cprjq)
     end if
   end if
   ABI_DATATYPE_DEALLOCATE(pawrhoij1)
   ABI_DATATYPE_DEALLOCATE(cprjq)
   ABI_DATATYPE_DEALLOCATE(cprj)
   if(xmpi_paral==1)  then
     ABI_DEALLOCATE(mpi_enreg%proc_distrb)
     ABI_DEALLOCATE(mpi_enreg%my_kpttab)
   end if
   call hdr%free()

   ! Clean band structure datatypes (should use it more in the future !)
   call ebands_free(ebands_k)
   call ebands_free(ebands_kq)
   if(.not.kramers_deg) then
     call ebands_free(ebands_kmq)
   end if

!  %%%% Parallelization over perturbations %%%%%
!  *Redefine output/log files
   call localredirect(mpi_enreg%comm_cell,mpi_enreg%comm_world,npert_io,mpi_enreg%paral_pert,0)

   ABI_DEALLOCATE(bz2ibz_smap)
   ABI_DEALLOCATE(distrbflags)

   call timab(146,2,tsec)
   if(iexit/=0) exit
 end do ! End loop on perturbations
 ABI_DEALLOCATE(zeff)

!%%%% Parallelization over perturbations %%%%%
!*Restore default communicators
 call unset_pert_comm(mpi_enreg)
!*Gather output/log files
 ABI_ALLOCATE(dyn,(npert_io))
 if (npert_io>0) dyn=1
 call localrdfile(mpi_enreg%comm_pert,mpi_enreg%comm_world,.true.,npert_io,mpi_enreg%paral_pert,0,dyn)
 ABI_DEALLOCATE(dyn)
!*Restore PAW on-site data
 if (paral_pert_inplace) then
   call unset_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,paw_an,paw_ij,pawfgrtab,pawrhoij)
 else
   call unset_pert_paw(dtset,mpi_enreg,my_natom,old_atmtab,old_comm_atom,&
&   paw_an,paw_ij,pawfgrtab,pawrhoij,&
&   paw_an_out=paw_an_pert,paw_ij_out=paw_ij_pert,&
&   pawfgrtab_out=pawfgrtab_pert,pawrhoij_out=pawrhoij_pert)
 end if

!#################################################################################
!Calculate the second-order eigenvalues for a wavevector Q

 call timab(147,1,tsec)
 smdelta = dtset%smdelta
 bdeigrf = dtset%bdeigrf
 if(dtset%bdeigrf == -1) bdeigrf = dtset%mband

 if(dtset%ieig2rf > 0) then

   if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&   (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
     call wrtout(std_out,'Reading the dense grid WF file',"COLL")
!    We get the Abinit header of the file hdr_fine as ouput
!    We get eigenq_fine(mband,hdr_fine%nkpt,hdr_fine%nsppol) as ouput
     fname = dtfil%fnameabi_wfkfine
     if (dtset%iomode == IO_MODE_ETSF) fname = nctk_ncify(dtfil%fnameabi_wfkfine)

     call wfk_read_eigenvalues(fname,eigenq_fine,hdr_fine,mpi_enreg%comm_world)
     ABI_CHECK(SIZE(eigenq_fine,DIM=1)==Dtset%mband,"Size eigenq_fine != mband")
   end if
! DBSP ==> Has been changed to be able to make Bandstructure calculation
!   if(dtset%kptopt==3 .or. dtset%kptopt==0)then
   if(dtset%kptopt==3 .or. dtset%kptopt==0 .or. dtset%kptopt < -4 .or. dtset%nsym==1) then
!END
     if (dtset%nsym > 1) then ! .and. dtset%efmas==0) then
       MSG_ERROR("Symmetries are not implemented for temperature dependence calculations")
     end if
     write(std_out,*) 'Entering: eig2stern'
     if(smdelta>0)then
       if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&       (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
         call eig2stern(occ_pert,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0_pert,eigenq_pert,&
&         eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,gh0c1_pert,gh1c_pert,&
&         dtset%ieig2rf,istwfk_pert,dtset%mband,mk1mem_rbz,mpert,dtset%natom,mpi_enreg,mpw1,nkpt_rbz,&
&         npwar1_pert,dtset%nspinor,dtset%nsppol,smdelta,dtset,eigbrd,eigenq_fine,hdr_fine,hdr0)
       else
         call eig2stern(occ_pert,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0_pert,eigenq_pert,&
&         eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,gh0c1_pert,gh1c_pert,&
&         dtset%ieig2rf,istwfk_pert,dtset%mband,mk1mem_rbz,mpert,dtset%natom,mpi_enreg,mpw1,nkpt_rbz,&
&         npwar1_pert,dtset%nspinor,dtset%nsppol,smdelta,dtset,eigbrd)
       end if
     else
       if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&       (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
         call eig2stern(occ_pert,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0_pert,eigenq_pert,&
&         eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,gh0c1_pert,gh1c_pert,&
&         dtset%ieig2rf,istwfk_pert,dtset%mband,mk1mem_rbz,mpert,dtset%natom,mpi_enreg,mpw1,nkpt_rbz,&
&         npwar1_pert,dtset%nspinor,dtset%nsppol,smdelta,dtset,eigbrd,eigenq_fine,hdr_fine,hdr0)
       else
         call eig2stern(occ_pert,bdeigrf,clflg,cg1_pert,dim_eig2nkq,dim_eig2rf,eigen0_pert,eigenq_pert,&
&         eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,gh0c1_pert,gh1c_pert,&
&         dtset%ieig2rf,istwfk_pert,dtset%mband,mk1mem_rbz,mpert,dtset%natom,mpi_enreg,mpw1,nkpt_rbz,&
&         npwar1_pert,dtset%nspinor,dtset%nsppol,smdelta,dtset)
       end if
     end if
     call wrtout(std_out, 'Leaving: eig2stern', "COLL")

     if (dtset%ieig2rf==1.or.dtset%ieig2rf==2) then
       if (me==master) then
!        print _EIGR2D file for this perturbation
         dscrpt=' Note : temporary (transfer) database '
         unitout = dtfil%unddb
         vrsddb=100401

         call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&         1,xred=xred,occ=occ_pert)

         call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_eigr2d, dtfil%unddb)

         call outbsd(bdeigrf,dtset,eig2nkq,dtset%natom,nkpt_rbz,unitout)
!        print _EIGI2D file for this perturbation
         if(smdelta>0) then

           unitout = dtfil%unddb
           call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_eigi2d, unitout)

           call outbsd(bdeigrf,dtset,eigbrd,dtset%natom,nkpt_rbz,unitout)
         end if !smdelta

         call ddb_hdr_free(ddb_hdr)

!        Output of the EIGR2D.nc file.
         fname = strcat(dtfil%filnam_ds(4),"_EIGR2D.nc")

!        Electronic band energies.
         bantot= dtset%mband*dtset%nkpt*dtset%nsppol
         call ebands_init(bantot,gkk_ebands,dtset%nelect,doccde,eigen0_pert,hdr0%istwfk,hdr0%kptns,&
&         hdr0%nband, hdr0%nkpt,hdr0%npwarr,hdr0%nsppol,hdr0%nspinor,&
&         hdr0%tphysel,hdr0%tsmear,hdr0%occopt,hdr0%occ,hdr0%wtk,&
&         hdr0%charge, hdr0%kptopt, hdr0%kptrlatt_orig, hdr0%nshiftk_orig, hdr0%shiftk_orig, &
&         hdr0%kptrlatt, hdr0%nshiftk, hdr0%shiftk)

!        Second order derivative EIGR2D (real and Im)
         call eigr2d_init(eig2nkq,eigr2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
#ifdef HAVE_NETCDF
         NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGR2D file")
         NCF_CHECK(crystal%ncwrite(ncid))
         NCF_CHECK(ebands_ncwrite(gkk_ebands, ncid))
         call eigr2d_ncwrite(eigr2d,dtset%qptn(:),dtset%wtq,ncid)
         NCF_CHECK(nf90_close(ncid))
#endif
         if(smdelta>0) then
!          Output of the EIGI2D.nc file.
           fname = strcat(dtfil%filnam_ds(4),"_EIGI2D.nc")
!          Broadening EIGI2D (real and Im)
           call eigr2d_init(eigbrd,eigi2d,dtset%mband,hdr0%nsppol,nkpt_rbz,dtset%natom)
#ifdef HAVE_NETCDF
           NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EIGI2D file")
           NCF_CHECK(crystal%ncwrite(ncid))
           NCF_CHECK(ebands_ncwrite(gkk_ebands, ncid))
           call eigr2d_ncwrite(eigi2d,dtset%qptn(:),dtset%wtq,ncid)
           NCF_CHECK(nf90_close(ncid))
#endif
         end if
       end if
     end if !ieig2rf==1.or.ieig2rf==2
   else
     write(message,'(3a)')&
&     'K point grids must be the same for every perturbation: eig2stern not called',ch10,&
&     'Action: Put kptopt=3 '
     MSG_WARNING(message)
   end if !kptopt
   ABI_DEALLOCATE(gh1c_pert)
   ABI_DEALLOCATE(gh0c1_pert)
   ABI_DEALLOCATE(cg1_pert)
   ABI_DEALLOCATE(kpt_rbz_pert)
   ABI_DEALLOCATE(istwfk_pert)
   ABI_DEALLOCATE(npwarr_pert)
   ABI_DEALLOCATE(npwar1_pert)
   ABI_DEALLOCATE(npwtot_pert)
   ABI_DEALLOCATE(occ_pert)
 end if  !if dtset%ieig2rf

 ! Calculation of effective masses.
 if(dtset%efmas == 1) then
   call efmas_main(cg0_pert,cg1_pert,dim_eig2rf,dtset,efmasdeg,efmasval,eigen0_pert,&
&   eigen1_pert,gh0c1_pert,gh1c_pert,istwfk_pert,mpert,mpi_enreg,nkpt_rbz,npwarr_pert,rprimd)

   ABI_DEALLOCATE(gh1c_pert)
   ABI_DEALLOCATE(gh0c1_pert)
   ABI_DEALLOCATE(cg1_pert)
   ABI_DEALLOCATE(istwfk_pert)
   ABI_DEALLOCATE(npwarr_pert)
   ABI_DEALLOCATE(cg0_pert)

   if(dtset%prtefmas==1)then
     fname = strcat(dtfil%filnam_ds(4),"_EFMAS.nc")
     !write(std_out,*)' dfpt_looppert: will write ',fname
#ifdef HAVE_NETCDF
     NCF_CHECK_MSG(nctk_open_create(ncid, fname, xmpi_comm_self), "Creating EFMAS file")
     NCF_CHECK(crystal%ncwrite(ncid))
!    NCF_CHECK(ebands_ncwrite(ebands_k, ncid)) ! At this stage, ebands_k is not available
     call print_efmas(efmasdeg,efmasval,kpt_rbz_pert,ncid)
     NCF_CHECK(nf90_close(ncid))
#endif
   endif

   call efmas_analysis(dtset,efmasdeg,efmasval,kpt_rbz_pert,mpi_enreg,nkpt_rbz,rprimd)

   ABI_DEALLOCATE(kpt_rbz_pert)
 end if

!Free memory.
 if(dtset%ieig2rf /= 3 .and. dtset%ieig2rf /= 4 .and. dtset%ieig2rf /= 5) call hdr0%free()
 ABI_DEALLOCATE(eigen0_copy)
 call crystal%free()

 ! GKK stuff (deprecated)
 call ebands_free(gkk_ebands)
 call eigr2d_free(eigr2d)
 call eigr2d_free(eigi2d)

 call timab(147,2,tsec)
!######################################################################################

!Get ddk file information, for later use in dfpt_dyout
 ddkfil(:)=0
 do idir=1,3
   file_index(1)=idir+dtset%natom*3
   call appdig(file_index(1),dtfil%fnamewffddk,fiwfddk)
!  Check that ddk file exists
   t_exist = file_exists(fiwfddk)
   if (.not. t_exist) then
     ! Trick needed to run Abinit test suite in netcdf mode.
     t_exist = file_exists(nctk_ncify(fiwfddk))
     if (t_exist) then
       write(message,"(3a)")"- File: ",trim(fiwfddk)," does not exist but found netcdf file with similar name."
       call wrtout(std_out,message,'COLL')
       fiwfddk = nctk_ncify(fiwfddk)
     end if
   end if

!  If the file exists set ddkfil to a non-zero value
   if (t_exist) ddkfil(idir)=20+idir
 end do

 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(pert_calc)
 if (psps%usepaw==1) then
   ABI_DEALLOCATE(dimcprj_srt)
 end if

!destroy dtset_tmp
 if (dtset%prepgkk /= 0) then ! .and. dtset%use_nonscf_gkk == 1) then !Later uncomment this - in scf case rhor1_save is used below only for testing
   ABI_DEALLOCATE(rhor1_save)
   ABI_DEALLOCATE(blkflg_save)
   call dtset_tmp%free()
   ABI_DATATYPE_DEALLOCATE(dtset_tmp)
 end if

!In paral_pert-case some array's have to be reconstructed
 if(mpi_enreg%paral_pert==1) then
   ABI_ALLOCATE(buffer1,(2,3,mpert,3,mpert*(2+psps%usepaw)))
   buffer1(:,:,:,:,1:mpert)=d2lo(:,:,:,:,:)
   buffer1(:,:,:,:,1+mpert:2*mpert)=d2nl(:,:,:,:,:)
   if (psps%usepaw==1) then
     buffer1(:,:,:,:,1+2*mpert:3*mpert)=d2ovl(:,:,:,:,:)
   end if
   call xmpi_sum(buffer1,mpi_enreg%comm_pert,ierr)
   call xmpi_sum(blkflg,mpi_enreg%comm_pert,ierr)
   if(dtset%prtbbb==1) then
     call xmpi_sum(d2bbb,mpi_enreg%comm_pert,ierr)
   end if
   d2lo(:,:,:,:,:)=buffer1(:,:,:,:,1:mpert)
   d2nl(:,:,:,:,:)=buffer1(:,:,:,:,1+mpert:2*mpert)
   if (psps%usepaw==1) then
     d2ovl(:,:,:,:,:)=buffer1(:,:,:,:,1+2*mpert:3*mpert)
   end if
   ABI_DEALLOCATE(buffer1)
 end if

 if ( associated(old_atmtab)) then
   ABI_DEALLOCATE(old_atmtab)
   nullify(old_atmtab)
 end if

 call timab(141,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_looppert
!!***

!!****f* ABINIT/getcgqphase
!! NAME
!! getcgqphase
!!
!! FUNCTION
!! extract phases from wave functions, to cancel contributions to gkk matrix elements
!!
!! INPUTS
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  timrev = flag for use of time reversal symmetry
!!  cg = input wavefunctions
!!  mcg = dimension of cg = nspinor*mband*mpw*mkmem
!!  cgq = input wavefunctions at k+q
!!  mcgq = dimension of cgq = nspinor*mband*mpw*mkmem
!!  mpi_enreg = datastructure for mpi communication
!!  nkpt_rbz = number of k-points in reduced zone for present q point
!!  npwarr = array of numbers of plane waves for each k-point
!!  npwar1 = array of numbers of plane waves for each k+q point
!!
!! OUTPUT
!!  phasecg = phase of different wavefunction products <k,n | k+q,n'>
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      smatrix,wrtout,xmpi_barrier,xmpi_sum_master
!!
!! SOURCE

subroutine getcgqphase(dtset, timrev, cg,  mcg,  cgq, mcgq, mpi_enreg, nkpt_rbz, npwarr, npwar1, phasecg)

!Arguments -------------------------------
 ! scalars
 integer, intent(in) :: mcg, mcgq, timrev
 integer, intent(in) :: nkpt_rbz
 type(dataset_type), intent(in) :: dtset

 ! arrays
 integer, intent(in) :: npwarr(nkpt_rbz)
 integer, intent(in) :: npwar1(nkpt_rbz)
 real(dp), intent(in) :: cg(2,mcg)
 real(dp), intent(in) :: cgq(2,mcgq)
 type(MPI_type), intent(in) :: mpi_enreg
 real(dp),intent(out) :: phasecg(2, dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)

!Local variables -------------------------
 ! local vars
 integer :: icg, icgq, isppol, ikpt, ipw
 integer :: istate, iband1, iband2
 integer :: npw_k, npw_q
 integer :: me, ierr, master, spaceComm, nprocs
 integer :: usepaw
 integer :: ddkflag, itrs, job, maxbd, mcg1_k, minbd, shiftbd
 real(dp) :: normsmat

 integer, allocatable :: sflag_k(:)
 integer, allocatable :: pwind_k(:)

 real(dp) :: cg1_dummy(1,1)
 real(dp) :: smat_inv_dummy(1,1,1)
 real(dp) :: smat_k_paw_dummy(1,1,1)
 real(dp) :: dtm_k_dummy(2)
 real(dp), allocatable :: smat_k(:,:,:)
 real(dp), allocatable :: pwnsfac_k(:,:)
 logical, allocatable :: my_kpt(:,:)
 character(len=500) :: message

! *********************************************************************

 ABI_ALLOCATE(smat_k,(2,dtset%mband,dtset%mband))
 ABI_ALLOCATE(sflag_k,(dtset%mband))

!dummy use of timrev so abirules stops complaining.
 icg = timrev

!!MPI data for future use
 spaceComm=mpi_enreg%comm_cell
 nprocs=xmpi_comm_size(spaceComm)
 master=0
 me=mpi_enreg%me_kpt

 ABI_ALLOCATE(my_kpt, (nkpt_rbz, dtset%nsppol))
 my_kpt = .true.
 if (mpi_enreg%nproc_kpt > 1) then
   do isppol = 1, dtset%nsppol
     do ikpt = 1, nkpt_rbz
       my_kpt(ikpt, isppol) = .not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,&
&       dtset%nband(ikpt),isppol,me))
     end do
   end do
 end if


!make trivial association of G vectors: we just want <psi_k| psi_k+q>
!TODO: check this is correct wrt arrangement of kg vectors for k+q
!looks ok : usually made in initberry, from the translations associated
!to the symops, scalar product with the G vectors. The symop is the one
!used to go from the irreducible k to the full zone k. In present context
!we should be using only the reduced zone, and anyhow have the same k-grid
!for the gkk matrix elements and for the cg here...
 ABI_ALLOCATE(pwind_k,(dtset%mpw))
 ABI_ALLOCATE(pwnsfac_k,(4,dtset%mpw))
 do ipw = 1, dtset%mpw
   pwind_k(ipw) = ipw
   pwnsfac_k(1,ipw) = one
   pwnsfac_k(2,ipw) = zero
   pwnsfac_k(3,ipw) = one
   pwnsfac_k(4,ipw) = zero
 end do

!flags for call to smatrix
 usepaw = 0 ! for now
 ddkflag = 0
 itrs = 0
 job = 0
 maxbd = 1
 mcg1_k = 1
 minbd = 1
 shiftbd = 1

!from overlap matrix for each wavefunction, extract phase
 icg = 0
 icgq = 0
 istate = 0

 phasecg = zero
 do isppol = 1, dtset%nsppol
   do ikpt = 1, nkpt_rbz
!    each proc only has certain k
     if (.not. my_kpt(ikpt, isppol)) then
       istate = istate +  dtset%nband(ikpt)*dtset%nband(ikpt)
       cycle
     end if

     npw_k = npwarr(ikpt)
     npw_q= npwar1(ikpt)

!    TODO: question: are the k-points in the ibz correctly ordered in cg and cgq? if not the icg below have to be adapted.
     sflag_k = 0 ! make sure all elements are calculated
     smat_k = zero

     call smatrix(cg, cgq, cg1_dummy, ddkflag, dtm_k_dummy, icg, icgq,&
&     itrs, job, maxbd, mcg, mcgq, mcg1_k, minbd,dtset%mpw, dtset%mband, dtset%mband,&
&     npw_k, npw_q, dtset%nspinor, pwind_k, pwnsfac_k, sflag_k, shiftbd,&
&     smat_inv_dummy, smat_k, smat_k_paw_dummy, usepaw)

     icg  = icg  + npw_k*dtset%nspinor*dtset%nband(ikpt)
     icgq = icgq + npw_q*dtset%nspinor*dtset%nband(ikpt)

     do iband1 = 1, dtset%nband(ikpt)
       do iband2 = 1, dtset%nband(ikpt)
         istate = istate + 1
!        normalise the overlap matrix element to get just the phase difference phi_k - phi_k+q
         normsmat = sqrt(smat_k(1,iband2, iband1)**2 &
&         + smat_k(2,iband2, iband1)**2)
         if (normsmat > tol12) then
           phasecg(:,istate) = smat_k(:,iband2, iband1) / normsmat
!          NOTE: 21/9/2011 these appear to be always 1, i, or -i, to within 1.e-5 at worst!
         end if
       end do
     end do
   end do
 end do

!eventually do an mpi allreduce over the k-points for phasecg
 if (nprocs>1) then
   call xmpi_barrier(spaceComm)
   call xmpi_sum_master(phasecg,master,spaceComm,ierr)
   call xmpi_barrier(spaceComm)
   if (1==1) then
     write(message,'(a)') '  In getcgqphase - contributions to phasecg collected'
     call wrtout(std_out,message,'PERS')
   end if
 end if

 ABI_DEALLOCATE(sflag_k)
 ABI_DEALLOCATE(smat_k)
 ABI_DEALLOCATE(pwind_k)
 ABI_DEALLOCATE(pwnsfac_k)
 ABI_DEALLOCATE(my_kpt)

end subroutine getcgqphase
!!***

!!****f* ABINIT/dfpt_prtene
!!
!! NAME
!! dfpt_prtene
!!
!! FUNCTION
!! Print components of second derivative of total energy in nice format
!!
!! INPUTS
!! eberry=energy associated with Berry phase
!! edocc=correction to 2nd-order total energy coming from changes of occupation
!! eeig0=0th-order eigenenergies part of 2nd-order total energy
!! eew=Ewald part of 2nd-order total energy
!! efrhar=hartree frozen-wavefunction part of 2nd-order tot. en.
!! efrkin=kinetic frozen-wavefunction part of 2nd-order tot. en.
!! efrloc=local psp. frozen-wavefunction part of 2nd-order tot. en.
!! efrnl=nonlocal psp. frozen-wavefunction part of 2nd-order tot. en
!! efrx1=xc core corr.(1) frozen-wavefunction part of 2nd-order tot. en
!! efrx2=xc core corr.(2) frozen-wavefunction part of 2nd-order tot. en
!! ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!   for strain perturbation only (zero otherwise, and not used)
!! ehart1=1st-order Hartree part of 2nd-order total energy
!! eii=pseudopotential core part of 2nd-order total energy
!! ek0=0th-order kinetic energy part of 2nd-order total energy.
!! ek1=1st-order kinetic energy part of 2nd-order total energy.
!! eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!! elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!! enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!! enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!! eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!       PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008) [[cite:Audouze2008]]
!! epaw1=1st-order PAW on-site part of 2nd-order total energy.
!! evdw=DFT-D semi-empirical part of 2nd-order total energy
!! exc1=1st-order exchange-correlation part of 2nd-order total energy
!! iout=unit number to which output is written
!! ipert=type of the perturbation
!! natom=number of atoms in unit cell
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!!
!! OUTPUT
!!  (only writing)
!!
!! NOTES
!! all energies in Hartree
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine dfpt_prtene(berryopt,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,eovl1,epaw1,evdw,exc1,iout,ipert,natom,&
&  usepaw,usevdw)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,iout,ipert,natom,usepaw,usevdw
 real(dp),intent(in) :: eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1
 real(dp),intent(in) :: efrx2,ehart01,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1
 real(dp),intent(in) :: eovl1,epaw1,evdw,exc1

!Local variables -------------------------
!scalars
 integer :: nn
 logical :: berry_activated
 real(dp) :: enl1_effective,erelax,etotal
 character(len=10) :: numb
 character(len=10),parameter :: numbstr(20) = &
&  (/'One       ','Two       ','Three     ','Four      ','Five      ', &
&    'Six       ','Seven     ','Eight     ','Nine      ','Ten       ', &
&    'Eleven    ','Twelve    ','Thirteen  ','Fourteen  ','Fifteen   ', &
&    'Sixteen   ','Seventeen ','Eighteen  ','Nineteen  ','Twenty    '/)
 character(len=500) :: message

! *********************************************************************

!Count and print the number of components of 2nd-order energy
!MT feb 2015: this number is wrong! Should change it but
!             need to change a lot of ref. files
 berry_activated=(berryopt== 4.or.berryopt== 6.or.berryopt== 7.or. &
& berryopt==14.or.berryopt==16.or.berryopt==17)
 if (ipert==natom+1) nn=8
 if (ipert==natom+5) nn=8
 if (ipert==natom+2) nn=7
 if (ipert>=1.and.ipert<=natom) nn=13
 if (ipert==natom+3.or.ipert==natom+4) nn=17
 if (ipert==natom+2.and.berry_activated) nn=nn+1
 if (ipert==natom+10.or.ipert==natom+11) nn=1 ! means nothing,
! because we do not compute derivatives of the energy in this case
 if (usepaw==1) nn=nn+1
 if (usevdw==1) nn=nn+1
 write(message, '(4a)' ) ch10,&
& ' ',trim(numbstr(nn)),' components of 2nd-order total energy (hartree) are '
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 numb='1,2,3'
 write(message, '(3a)' )&
& ' ',trim(numb),': 0th-order hamiltonian combined with 1st-order wavefunctions'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')
 write(message, '(a,es17.8,a,es17.8,a,es17.8)' )&
& '     kin0=',ek0,   ' eigvalue=',eeig0,'  local=',eloc0
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 numb='4,5,6';if( ipert==natom+3.or.ipert==natom+4) numb='4,5,6,7'
 write(message, '(3a)' )&
& ' ',trim(numb),': 1st-order hamiltonian combined with 1st and 0th-order wfs'
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')
 if(ipert/=natom+1.and.ipert/=natom+2)then
   write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a)' ) &
&   ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&   ' note that "loc psp" includes a xc core correction that could be resolved'
 else if(ipert==natom+1) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '     kin1=',ek1,   '  Hartree=',ehart1,'     xc=',exc1
 else if(ipert==natom+2) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '    dotwf=',enl1,  '  Hartree=',ehart1,'     xc=',exc1
 end if
 if(ipert==natom+3 .or. ipert==natom+4) then
   write(message, '(a,es17.8,a,es17.8,a,es17.8,a,a,es17.8)' ) &
&   ' loc psp =',elpsp1,'  Hartree=',ehart1,'     xc=',exc1,ch10,&
&   '     kin1=',ek1
 end if
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 enl1_effective=enl1;if (ipert==natom+2) enl1_effective=zero
 numb='7,8,9';if( ipert==natom+3.or.ipert==natom+4) numb='8,9,10'
 write(message, '(5a,es17.8,a,es17.8,a,es17.8)' )&
& ' ',trim(numb),': eventually, occupation + non-local contributions',ch10,&
& '    edocc=',edocc,'     enl0=',enl0,'   enl1=',enl1_effective
 call wrtout(iout,message,'COLL')
 call wrtout(std_out,message,'COLL')

 if (usepaw==1) then
   numb='10';if( ipert==natom+3.or.ipert==natom+4) numb='11'
   write(message, '(3a,es17.8)' )&
&   ' ',trim(numb),': eventually, PAW "on-site" Hxc contribution: epaw1=',epaw1
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 if(ipert/=natom+10 .and.ipert/=natom+11) then
   erelax=0.0_dp
   if(ipert>=1.and.ipert<=natom)then
     erelax=ek0+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1
   else if(ipert==natom+1.or.ipert==natom+2)then
     erelax=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1+epaw1
   else if(ipert==natom+3.or.ipert==natom+4)then
     erelax=ek0+edocc+eeig0+eloc0+ek1+elpsp1+ehart1+exc1+enl0+enl1+epaw1
   else if(ipert==natom+5)then
     erelax=ek0+edocc+eeig0+eloc0+ek1+elpsp1+ehart1+exc1+enl0+enl1+epaw1
   end if
   enl1_effective=enl1
   if (ipert==natom+1.or.ipert==natom+2) then
     if (1.0_dp+enl1/10.0_dp==1.0_dp) enl1_effective=zero
   end if

   numb='1-9';if (usepaw==1) numb='1-10'
   if( ipert==natom+3.or.ipert==natom+4) then
     numb='1-10';if (usepaw==1) numb='1-11'
   end if
   write(message, '(5a,es17.8)' )&
&   ' ',trim(numb),' gives the relaxation energy (to be shifted if some occ is /=2.0)',&
&   ch10,'   erelax=',erelax
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

 if(ipert>=1.and.ipert<=natom)then

   numb='10,11,12';if (usepaw==1) numb='11,12,13'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ',&
&   'frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   ' fr.local=',efrloc,' fr.nonlo=',efrnl,'  Ewald=',eew
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   write(message, '(a,es16.6)' )' dfpt_prtene : non-relax=',efrloc+efrnl+eew
   call wrtout(std_out,message,'COLL')

   numb='13,14';if (usepaw==1) numb='14,15'
   write(message, '(3a)' )&
&   ' ',trim(numb),' Frozen wf xc core corrections (1) and (2)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8)' ) &
&   ' frxc 1  =',efrx1,'  frxc 2 =',efrx2
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if (usepaw==1) then
     numb='16'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   if (usevdw==1) then
     numb='15';if (usepaw==1) numb='17'
     write(message, '(3a,es17.8)' )&
&     ' ',trim(numb),' Contribution from van der Waals DFT-D: evdw =',evdw
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if

   write(message, '(a)' )' Resulting in : '
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax+eew+efrloc+efrnl+efrx1+efrx2+evdw
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,es20.10,a)' ) &
&   '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,a)' ) &
&   '    (  non-var. 2DEtotal :',&
&   0.5_dp*(elpsp1+enl1)+eovl1+eew+efrloc+efrnl+efrx1+efrx2+evdw,' Ha)',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else if(ipert==natom+1.or.ipert==natom+2)then
   if (ipert==natom+1.and.usepaw==1) then
     numb='11'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   write(message,*)' No Ewald or frozen-wf contrib.:',' the relaxation energy is the total one'
   if(berry_activated) then
     numb='10';
     write(message,'(3a,es20.10)')' ',trim(numb),' Berry phase energy :',eberry
   end if
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a)' ) &
&   '    (  non-var. 2DEtotal :',0.5_dp*(ek1+enl1_effective)+eovl1,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

 else if(ipert==natom+3 .or. ipert==natom+4) then
   numb='11,12,13';if (usepaw==1) numb='12,13,14'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ','frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '  fr.hart=',efrhar,'   fr.kin=',efrkin,' fr.loc=',efrloc
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   numb='14,15,16';if (usepaw==1) numb='15,16,17'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ','frozen-wavefunctions and Ewald'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8,a,es17.8,a,es17.8)' ) &
&   '  fr.nonl=',efrnl,'    fr.xc=',efrx1,'  Ewald=',eew
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')

   numb='17';if (usepaw==1) numb='18'
   write(message, '(4a)' )&
&   ' ',trim(numb),' Non-relaxation  contributions : ','pseudopotential core energy'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es17.8)' ) '  pspcore=',eii
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   if (usepaw==1) then
     numb='19'
     write(message, '(5a,es17.8)' )&
&     ' ',trim(numb),' Contribution from 1st-order change of wavefunctions overlap',&
&     ch10,' eovl1 =',eovl1
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if
   if (usevdw==1) then
     numb='18';if (usepaw==1) numb='20'
     write(message, '(3a,es17.8)' )&
&     ' ',trim(numb),' Contribution from van der Waals DFT-D: evdw =',evdw
     call wrtout(iout,message,'COLL')
     call wrtout(std_out,message,'COLL')
   end if

   write(message, '(a,es16.6)' )' dfpt_prtene : non-relax=',&
&   efrhar+efrkin+efrloc+efrnl+efrx1+eew+evdw
   call wrtout(std_out,message,'COLL')
   write(message, '(a)' )' Resulting in : '
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   etotal=erelax+efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii+evdw
   write(message, '(a,e20.10,a,e22.12,a)' ) &
&   ' 2DEtotal=',etotal,' Ha. Also 2DEtotal=',etotal*Ha_eV,' eV'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,es20.10,a)' ) &
&   '    (2DErelax=',erelax,' Ha. 2DEnonrelax=',etotal-erelax,' Ha)'
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
   write(message, '(a,es20.10,a,a)' ) &
&   '    (  non-var. 2DEtotal :',&
&   0.5_dp*(elpsp1+enl1+ek1+ehart01)+eovl1+&
&   efrhar+efrkin+efrloc+efrnl+efrx1+eew+eii+evdw,' Ha)',ch10
   call wrtout(iout,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

end subroutine dfpt_prtene
!!***

!!****f* ABINIT/eigen_meandege
!! NAME
!! eigen_meandege
!!
!! FUNCTION
!! This routine takes the mean values of the responses
!! for the eigenstates that are degenerate in energy.
!!
!! INPUTS
!!  eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)= input eigenresp
!!       eigenrep(2*mband**2*nkpt*nsppol) for first-order derivatives of the eigenvalues
!!       eigenrep(mband*nkpt*nsppol) for Fan or Debye-Waller second-order derivatives of the eigenvalues
!!  mband= maximum number of bands
!!  natom= number of atoms in the unit cell
!!  nkpt= number of k-points
!!  nsppol= 1 for unpolarized, 2 for spin-polarized
!!  option= 1 for eigen(1), 2 for eigen(2) - Fan or Debye-Waller
!!
!! OUTPUT
!!  eigenresp_mean(mband*nkpt*nsppol)= eigenresp, averaged over degenerate states
!!
!! PARENTS
!!      dfpt_looppert,respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine eigen_meandege(eigen0,eigenresp,eigenresp_mean,mband,nband,nkpt,nsppol,option)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,nkpt,nsppol,option
 integer,intent(in) :: nband(nkpt*nsppol)

!arrays
 real(dp),intent(in) :: eigen0(mband*nkpt*nsppol)
 real(dp),intent(in) :: eigenresp((3-option)*mband**(3-option)*nkpt*nsppol)
 real(dp),intent(out) :: eigenresp_mean(mband*nkpt*nsppol)

!Local variables-------------------------------
!scalars
 integer :: bdtot_index,bd2tot_index,iband,ii,ikpt,isppol,nband_k
 real(dp) :: eig0,mean
 character(len=500) :: message
!arrays

! *********************************************************************

 if(option/=1 .and. option/=2)then
   write(message, '(a,i0)' )' The argument option should be 1 or 2, while it is found that option=',option
   MSG_BUG(message)
 end if

 bdtot_index=0 ; bd2tot_index=0
 do isppol=1,nsppol
   do ikpt=1,nkpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     if(option==1)then
       do iband=1,nband_k
         eigenresp_mean(iband+bdtot_index)=&
&         eigenresp(2*iband-1 + (iband-1)*2*nband_k + bd2tot_index)
       end do
     else if(option==2)then
       do iband=1,nband_k
         eigenresp_mean(iband+bdtot_index)=eigenresp(iband+bdtot_index)
       end do
     end if

!    Treat the case of degeneracies : take the mean of degenerate states
     if(nband_k>1)then
       eig0=eigen0(1+bdtot_index)
       ii=1
       do iband=2,nband_k
         if(eigen0(iband+bdtot_index)-eig0<tol8)then
           ii=ii+1
         else
           mean=sum(eigenresp_mean(iband-ii+bdtot_index:iband-1+bdtot_index))/ii
           eigenresp_mean(iband-ii+bdtot_index:iband-1+bdtot_index)=mean
           ii=1
         end if
         eig0=eigen0(iband+bdtot_index)
         if(iband==nband_k)then
           mean=sum(eigenresp_mean(iband-ii+1+bdtot_index:iband+bdtot_index))/ii
           eigenresp_mean(iband-ii+1+bdtot_index:iband+bdtot_index)=mean
         end if
       end do
     end if

     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2
   end do
 end do

end subroutine eigen_meandege
!!***

!!****f* ABINIT/dfpt_init_mag1
!! NAME
!!  dfpt_init_mag1
!!
!! FUNCTION
!!  Initial guess of the first order magnetization/density for magnetic field perturbation.
!!  The first order magnetization is set so as to zero out the first order XC magnetic field, which
!!  should minimize the second order XC energy (without taking self-consistency into account).
!!
!! INPUTS
!!  ipert = perturbtation type (works only for ipert==natom+5)
!!  idir  = direction of the applied magnetic field
!!  cplex = complex or real first order density and magnetization
!!  nfft  = dimension of the fft grid
!!  nspden= number of density matrix components
!!  nkxc  = number of kxc components
!!  vxc0(nfft,nspden)  = GS XC potential
!!  kxc0(nfft,nspden)  = GS XC derivatives
!!  rhor0(nfft,nspden) = GS density matrix
!!
!! OUTPUT
!!  rhor1(cplex*nfft) = first order density magnetization guess
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_init_mag1(ipert,idir,rhor1,rhor0,cplex,nfft,nspden,vxc0,kxc0,nkxc)

!Arguments ------------------------------------
 integer , intent(in)    :: ipert,idir,cplex,nfft,nspden,nkxc
 real(dp), intent(in)    :: vxc0(nfft,nspden),rhor0(nfft,nspden)
 real(dp), intent(in)    :: kxc0(nfft,nkxc)
 real(dp), intent(out)   :: rhor1(cplex*nfft,nspden)

!Local variables-------------------------------
 integer  :: ipt
 real(dp) :: bxc0,bxc1
 real(dp) :: m1_norm,m0_norm
 real(dp) :: f_dot_m
 real(dp) :: mdir(3),fdir(3)

! *************************************************************************
 ABI_UNUSED(ipert)

 if (nspden==2) then

   if(cplex==1) then
     do ipt=1,nfft
       bxc1=half*(half*(kxc0(ipt,1)+kxc0(ipt,3))-kxc0(ipt,2)) ! d/dm Bxc
       !this overestimates the first order magnetization because of n1 not taken into account
       m1_norm=-half*(1/bxc1)
       rhor1(ipt,1)=zero             ! rho_up+rho_dwn    => charge density
       rhor1(ipt,2)=half*m1_norm     ! rho_up=1/2(rho+m) => half*m
     end do
   else
     do ipt=1,cplex*nfft
       rhor1(ipt,:)=zero
     end do
   end if

 else if(nspden==4) then

   fdir=zero
   fdir(idir)= 1.0d0
   do ipt=1,nfft
     m0_norm=sqrt(rhor0(ipt,2)**2+rhor0(ipt,3)**2+rhor0(ipt,4)**2)
     mdir(1)=rhor0(ipt,2)/m0_norm
     mdir(2)=rhor0(ipt,3)/m0_norm
     mdir(3)=rhor0(ipt,4)/m0_norm
     f_dot_m=fdir(1)*mdir(1)+fdir(2)*mdir(2)+fdir(3)*mdir(3) ! projection of the field direction on m0

     bxc1=half*(half*(kxc0(ipt,1)+kxc0(ipt,3))-kxc0(ipt,2))  ! d/dm Bxc
     m1_norm=(-half/bxc1)*f_dot_m                            ! get an estimate of the norm of m1

     bxc0=-sqrt((half*(vxc0(ipt,1)-vxc0(ipt,2)))**2+vxc0(ipt,3)**2+vxc0(ipt,4)**2)
     if(cplex==1) then
       rhor1(ipt,1)=zero       ! rho_up+rho_dwn    => charge density
       rhor1(ipt,2)=m1_norm*mdir(1)-half*m0_norm/bxc0*(fdir(1)-f_dot_m*mdir(1))   ! m1x
       rhor1(ipt,3)=m1_norm*mdir(2)-half*m0_norm/bxc0*(fdir(2)-f_dot_m*mdir(2))   ! m1y
       rhor1(ipt,4)=m1_norm*mdir(3)-half*m0_norm/bxc0*(fdir(3)-f_dot_m*mdir(3))   ! m1z
       rhor1(ipt,:)=zero
     else
       rhor1(2*ipt-1,1)=zero       ! Re rho_up+rho_dwn
       rhor1(2*ipt-1,2)=m1_norm*mdir(1)-half*m0_norm/bxc0*(fdir(1)-f_dot_m*mdir(1))   ! m1x
       rhor1(2*ipt-1,3)=m1_norm*mdir(2)-half*m0_norm/bxc0*(fdir(2)-f_dot_m*mdir(2))   ! m1x
       rhor1(2*ipt-1,4)=m1_norm*mdir(3)-half*m0_norm/bxc0*(fdir(3)-f_dot_m*mdir(3))   ! m1x
       rhor1(2*ipt  ,1)=zero       ! Im rho_up+rho_dwn
       rhor1(2*ipt  ,2)=zero
       rhor1(2*ipt  ,3)=zero
       rhor1(2*ipt  ,4)=zero

       rhor1(2*ipt-1,1)=zero; rhor1(2*ipt,1)=zero
       rhor1(2*ipt-1,2)=zero; rhor1(2*ipt,2)=zero
       rhor1(2*ipt-1,3)=zero; rhor1(2*ipt,3)=zero
       rhor1(2*ipt-1,4)=zero; rhor1(2*ipt,4)=zero

     end if
   end do
 end if

end subroutine dfpt_init_mag1
!!***

end module m_dfpt_loopert
!!***

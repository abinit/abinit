!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_respfn_driver
!! NAME
!!  m_respfn_driver
!!
!! FUNCTION
!!  Subdriver for DFPT calculations.
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (XG, DRH, MT, MKV)
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

module m_respfn_driver

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_efmas_defs
 use m_abicore
 use m_xmpi
 use m_exit
 use m_wffile
 use m_errors
 use m_ebands
 use m_results_respfn
 use m_hdr
 use m_crystal
 use m_xcdata
 use m_dtset

 use m_time,        only : timab
 use m_fstrings,    only : strcat
 use m_symtk,       only : matr3inv, littlegroup_q, symmetrize_xred
 use m_fft,         only : zerosym, fourdp
 use m_kpts,        only : symkchk
 use m_geometry,    only : irreducible_set_pert
 use m_dynmat,      only : chkph3, d2sym3, q0dy3_apply, q0dy3_calc, wings3, dfpt_phfrq, sytens, dfpt_prtph, &
                           asria_calc, asria_corr, cart29, cart39, chneu9, dfpt_sydy
 use m_ddb,         only : DDB_VERSION
 use m_ddb_hdr,     only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_ddb_interpolate, only : outddbnc
 use m_occ,         only : newocc
 use m_efmas,       only : efmasdeg_free_array, efmasval_free_array
 use m_wfk,         only : wfk_read_eigenvalues
 use m_ioarr,       only : read_rhor
 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type
 use m_pawtab,      only : pawtab_type, pawtab_get_lsize
 use m_paw_an,      only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,      only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,   only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_copy, &
                           pawrhoij_bcast, pawrhoij_nullify, pawrhoij_inquire_dim
 use m_pawdij,      only : pawdij, symdij
 use m_pawfgr,      only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_finegrid,only : pawexpiqr
 use m_pawxc,       only : pawxc_get_nkxc
 use m_paw_dmft,    only : paw_dmft_type
 use m_paw_sphharm, only : setsym_ylm
 use m_paw_nhat,    only : nhatgrid,pawmknhat
 use m_paw_tools,   only : chkpawovlp
 use m_paw_denpot,  only : pawdenpot
 use m_paw_init,    only : pawinit,paw_gencond
 use m_kg,          only : getcut, getph, kpgio
 use m_eig2d,       only : eig2tot, elph2_fanddw
 use m_inwffil,     only : inwffil
 use m_spacepar,    only : hartre, setsym
 use m_mkrho,       only : mkrho
 use m_vdw_dftd2,   only : vdw_dftd2
 use m_vdw_dftd3,   only : vdw_dftd3
 use m_initylmg,    only : initylmg
 use m_pspini,      only : pspini
 use m_atm2fft,     only : atm2fft
 use m_dfpt_loopert,only : dfpt_looppert, eigen_meandege
 use m_rhotoxc,     only : rhotoxc
 use m_drivexc,     only : check_kxc
 use m_mklocl,      only : mklocl, mklocl_recipspace
 use m_common,      only : setup1, prteigrs
 use m_fourier_interpol, only : transgrid
 use m_paral_atom,     only : get_my_atmtab, free_my_atmtab
 use m_paw_occupancies, only : initrhoij
 use m_paw_correlations,only : pawpuxinit
 use m_mkcore,     only : mkcore, dfpt_mkcore
 use m_dfpt_elt,   only : dfpt_eltfrxc, dfpt_eltfrloc, dfpt_eltfrkin, dfpt_eltfrhar, elt_ewald, dfpt_ewald
 use m_d2frnl,     only : d2frnl

#if defined HAVE_GPU_CUDA
 use m_manage_cuda
#endif

 implicit none

 private
!!***

 public :: respfn
!!***

contains
!!***

!!****f* m_respfn_driver/respfn
!! NAME
!! respfn
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of Response functions.
!!
!! INPUTS
!!  codvsn=code version
!!  cpui=initial cpu time
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum single fft dimension
!!   | mkmem=Number of k points treated by this node
!!   | mpw=maximum number of planewaves in basis sphere (large number)
!!   | natom=number of atoms in unit cell
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | nkpt=number of k points
!!   | nspden=number of spin-density components
!!   | nsppol=number of channels for spin-polarization (1 or 2)
!!   | nsym=number of symmetry elements in space group
!!  mkmems(3)=array containing the tree values of mkmem (see above) (k-GS, k+q-GS and RF)
!!  mpi_enreg=information about MPI parallelization
!!  npwtot(nkpt)=number of planewaves in basis and boundary at each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  etotal=total energy (sum of 7 or 8 contributions) (hartree)
!!
!! SIDE EFFECTS
!!  iexit=index of "exit" on first line of file (0 if not found)
!!  occ(mband*nkpt*nsppol)=occup number for each band (often 2) at each k point
!!    Occupations number may have been read from a previous dataset...
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!    Some dimensions in pawrad have been set in driver.f
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!    Some dimensions in pawtab have been set in driver.f
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!    Before entering the first time in respfn, a significant part of psps
!!    has been initialized: the integers dimekb,lmnmax,lnmax,mpssang,mpssoang,
!!    mpsso,mgrid,ntypat,n1xccc,usepaw,useylm, and the arrays dimensioned to npsp
!!    All the remaining components of psps are to be initialized in the call
!!    to pspini.  The next time the code enters respfn, psps might be identical
!!    to the one of the previous dtset, in which case, no reinitialisation
!!    is scheduled in pspini.f .
!!  results_respfn <type(results_respfn_type)>=stores some results of respfn calls
!!
!! NOTES
!! USE OF FFT GRIDS:
!! =================
!! In case of PAW:
!! ---------------
!!    Two FFT grids are used:
!!    - A "coarse" FFT grid (defined by ecut)
!!      for the application of the Hamiltonian on the plane waves basis.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      Hamiltonian, wave-functions, density related to WFs (rhor here), ...
!!      are expressed on this grid.
!!    - A "fine" FFT grid (defined) by ecutdg)
!!      for the computation of the density inside PAW spheres.
!!      It is defined by nfftf, ngfftf, mgfftf, ...
!!      Total density, potentials, ...
!!      are expressed on this grid.
!! In case of norm-conserving:
!! ---------------------------
!!    - Only the usual FFT grid (defined by ecut) is used.
!!      It is defined by nfft, ngfft, mgfft, ...
!!      For compatibility reasons, (nfftf,ngfftf,mgfftf)
!!      are set equal to (nfft,ngfft,mgfft) in that case.
!!
!! PARENTS
!!      driver
!!
!! CHILDREN
!!      alloc_hamilt_gpu,atm2fft,check_kxc,chkpawovlp,chkph3,crystal_free
!!      crystal_init,d2frnl,d2sym3,ddb_hdr_free,ddb_hdr_init,ddb_hdr_open_write
!!      dealloc_hamilt_gpu,dfpt_dyfro,dfpt_dyout,dfpt_dyxc1,dfpt_eltfrhar
!!      dfpt_eltfrkin,dfpt_eltfrloc,dfpt_eltfrxc,dfpt_ewald,dfpt_gatherdy
!!      dfpt_looppert,dfpt_phfrq,dfpt_prtph,ebands_free,efmasdeg_free_array
!!      efmasval_free_array,eig2tot,eigen_meandege,elph2_fanddw,elt_ewald
!!      exit_check,fourdp,getcut,getph,hartre,hdr_free,hdr_init,hdr_update
!!      initrhoij,initylmg,inwffil,irreducible_set_pert,kpgio,littlegroup_q
!!      matr3inv,mkcore,mklocl,mkrho,newocc,nhatgrid,outddbnc,paw_an_free
!!      paw_an_init,paw_an_nullify,paw_gencond,paw_ij_free,paw_ij_init
!!      paw_ij_nullify,pawdenpot,pawdij,pawexpiqr,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawinit,pawmknhat,pawpuxinit
!!      pawrhoij_alloc,pawrhoij_bcast,pawrhoij_copy,pawrhoij_free
!!      pawrhoij_nullify,pawtab_get_lsize,prteigrs,pspini,q0dy3_apply
!!      q0dy3_calc,read_rhor,rhotoxc,setsym,setsym_ylm,setup1,status,symdij
!!      symmetrize_xred,sytens,timab,transgrid,vdw_dftd2,vdw_dftd3,wffclose
!!      wings3,wrtloctens,wrtout,xcdata_init,xmpi_bcast
!!
!! SOURCE

subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&
&  mkmems,mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,results_respfn,xred)

!Arguments ------------------------------------
 integer,intent(inout) :: iexit
 real(dp),intent(in) :: cpui
 real(dp),intent(inout) :: etotal !vz_i
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
 integer,intent(in) :: mkmems(3)
 integer,intent(inout) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)
 type(results_respfn_type),intent(inout) :: results_respfn

!Local variables-------------------------------
 integer,parameter :: formeig=0,level=10
 integer,parameter :: response=1,syuse=0,master=0,cplex1=1
 integer :: nk3xc
 integer :: analyt,ask_accurate,band_index,bantot,bdeigrf,coredens_method,cplex,cplex_rhoij
 integer :: dim_eig2nkq,dim_eigbrd,dyfr_cplex,dyfr_nondiag,gnt_option
 integer :: gscase,has_dijnd,has_diju,has_kxc,iatom,iatom_tot,iband,idir,ider,ierr,ifft,ii,ikpt,indx
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer :: initialized,ipert,ipert2,ireadwf0,iscf,iscf_eff,ispden,isppol
 integer :: itypat,izero,mcg,me,mgfftf,mk1mem,mkqmem,mpert,mu
 integer :: my_natom,n1,natom,n3xccc,nband_k,nfftf,nfftot,nfftotf,nhatdim,nhatgrdim
 integer :: nkpt_eff,nkpt_max,nkpt_rbz,nkxc,nkxc1,nspden_rhoij,ntypat,nzlmopt,openexit
 integer :: optcut,option,optgr0,optgr1,optgr2,optorth,optrad
 integer :: optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv
 integer :: outd2,pawbec,pawpiezo,prtbbb,psp_gencond,qzero,rdwr,rdwrpaw
 integer :: rfasr,rfddk,rfelfd,rfphon,rfstrs,rfuser,rf2_dkdk,rf2_dkde,rfmagn
 integer :: spaceworld,sumg0,sz1,sz2,tim_mkrho,timrev,usecprj,usevdw
 integer :: usexcnhat,use_sym,vloc_method,zero_by_symm
 logical :: has_full_piezo,has_allddk,is_dfpt=.true.,non_magnetic_xc
 logical :: paral_atom,qeq0,use_nhat_gga,call_pawinit
 real(dp) :: boxcut,compch_fft,compch_sph,cpus,ecore,ecut_eff,ecutdg_eff,ecutf
 real(dp) :: eei,eew,ehart,eii,ek,enl,entropy,enxc
 real(dp) :: epaw,epawdc,etot,evdw,fermie,gsqcut,gsqcut_eff,gsqcutc_eff,qphnrm,residm
 real(dp) :: ucvol,vxcavg
 character(len=fnlen) :: dscrpt
 character(len=fnlen) :: filename
 character(len=500) :: message
 type(ebands_t) :: bstruct
 type(hdr_type) :: hdr,hdr_fine,hdr0,hdr_den
 type(ddb_hdr_type) :: ddb_hdr
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 type(crystal_t) :: Crystal
 type(xcdata_type) :: xcdata
 integer :: ddkfil(3),ngfft(18),ngfftf(18),rfdir(3),rf2_dirs_from_rfpert_nl(3,3)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:,:,:,:),blkflgfrx1(:,:,:,:),blkflg1(:,:,:,:)
 integer,allocatable :: blkflg2(:,:,:,:),carflg(:,:,:,:),clflg(:,:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),l_size_atm(:),nattyp(:),npwarr(:)
 integer,allocatable :: pertsy(:,:),rfpert(:),rfpert_nl(:,:,:,:,:,:),symq(:,:,:),symrec(:,:,:)
 real(dp) :: dum_gauss(0),dum_dyfrn(0),dum_dyfrv(0),dum_eltfrxc(0)
 real(dp) :: dum_grn(0),dum_grv(0),dum_rhog(0),dum_vg(0)
 real(dp) :: dummy6(6),gmet(3,3),gmet_for_kg(3,3),gprimd(3,3),gprimd_for_kg(3,3),qphon(3)
 real(dp) :: rmet(3,3),rprimd(3,3),rprimd_for_kg(3,3),strn_dummy6(6),strv_dummy6(6),strsxc(6),tsec(2)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: becfrnl(:,:,:),cg(:,:),d2bbb(:,:,:,:,:,:),d2cart(:,:,:,:,:)
 real(dp),allocatable :: d2cart_bbb(:,:,:,:,:,:),d2eig0(:,:,:,:,:)
 real(dp),allocatable :: d2k0(:,:,:,:,:),d2lo(:,:,:,:,:),d2loc0(:,:,:,:,:)
 real(dp),allocatable :: d2matr(:,:,:,:,:),d2nfr(:,:,:,:,:),d2nl(:,:,:,:,:),d2ovl(:,:,:,:,:)
 real(dp),allocatable :: d2nl0(:,:,:,:,:),d2nl1(:,:,:,:,:),d2tmp(:,:,:,:,:),d2vn(:,:,:,:,:)
 real(dp),allocatable :: displ(:),doccde(:)
 real(dp),allocatable :: dyew(:,:,:,:,:),dyewq0(:,:,:),dyfrlo(:,:,:),dyfrlo_indx(:,:,:)
 real(dp),allocatable :: dyfrnl(:,:,:,:,:),dyfrwf(:,:,:,:,:),dyfrx1(:,:,:,:,:),dyvdw(:,:,:,:,:)
 real(dp),allocatable :: dyfrx2(:,:,:),eigen0(:),eigval(:),eigvec(:)
 real(dp),allocatable :: eig2nkq(:,:,:,:,:,:,:),eigbrd(:,:,:,:,:,:,:)
 real(dp),allocatable :: eigen_fan(:),eigen_ddw(:),eigen_fanddw(:)
 real(dp),allocatable :: eigen_fan_mean(:),eigen_ddw_mean(:)
 real(dp),allocatable :: eltcore(:,:),elteew(:,:),eltfrhar(:,:),eltfrkin(:,:)
 real(dp),allocatable :: eltfrloc(:,:),eltfrnl(:,:),eltfrxc(:,:),eltvdw(:,:),grtn_indx(:,:)
 real(dp),allocatable :: grxc(:,:),kxc(:,:),nhat(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: ph1d(:,:),ph1df(:,:),phfrq(:),phnons(:,:,:),piezofrnl(:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhowfg(:,:),rhowfr(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:)
 real(dp),allocatable :: vxc(:,:),work(:),xccc3d(:),ylm(:,:),ylmgr(:,:,:)
 real(dp),pointer :: eigenq_fine(:,:,:),eigen1_pert(:,:,:)
 real(dp),allocatable :: eigen0_pert(:),eigenq_pert(:),occ_rbz_pert(:)
 type(efmasdeg_type),allocatable :: efmasdeg(:)
 type(efmasval_type),allocatable :: efmasval(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij_read(:)

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(132,1,tsec)
 call timab(133,1,tsec)

!Some data for parallelism
 nkpt_max=50;if(xmpi_paral==1)nkpt_max=-1
 my_natom=mpi_enreg%my_natom
 paral_atom=(my_natom/=dtset%natom)
!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

!Structured debugging if dtset%prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,' respfn : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!Option input variables
 iscf=dtset%iscf

!Respfn input variables
 rfasr=dtset%rfasr   ; rfdir(1:3)=dtset%rfdir(1:3)
 rfddk=dtset%rfddk   ; rfelfd=dtset%rfelfd ; rfmagn=dtset%rfmagn
 rfphon=dtset%rfphon ; rfstrs=dtset%rfstrs
 rfuser=dtset%rfuser ; rf2_dkdk=dtset%rf2_dkdk ; rf2_dkde=dtset%rf2_dkde

 pawbec=0  ; if(psps%usepaw==1.and.(rfphon==1.or.(rfelfd==1.or.rfelfd==3))) pawbec=1
 pawpiezo=0; if(psps%usepaw==1.and.(rfstrs/=0.or.(rfelfd==1.or.rfelfd==3))) pawpiezo=1
!AM 10152015 -- WARNING --- the full calculation of the piezoelectric tensor
!from electric field perturbation is only available
!if nsym/=1 (strain perturbation is not symmetrized):
 has_full_piezo=.False. ; if(pawpiezo==1.and.dtset%nsym==1)  has_full_piezo=.True.
 usevdw=0;if (dtset%vdw_xc>=5.and.dtset%vdw_xc<=7) usevdw=1
!mkmem variables (mkmem is already argument)
 mkqmem=mkmems(2) ; mk1mem=mkmems(3)

 ntypat=psps%ntypat
 natom=dtset%natom
 nfftot=product(ngfft(1:3))
 nfftotf=product(ngfftf(1:3))

!LIKELY TO BE TAKEN AWAY
 initialized=0
 ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero
 eii=zero ; eew=zero ; ecore=zero

!Set up for iterations
 call setup1(dtset%acell_orig(1:3,1),bantot,dtset,&
  ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
  ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
  response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!In some cases (e.g. getcell/=0), the plane wave vectors have
! to be generated from the original simulation cell
 rprimd_for_kg=rprimd
 if (dtset%getcell/=0.and.dtset%usewvl==0) rprimd_for_kg=dtset%rprimd_orig(:,:,1)
 call matr3inv(rprimd_for_kg,gprimd_for_kg)
 gmet_for_kg=matmul(transpose(gprimd_for_kg),gprimd_for_kg)

!Define the set of admitted perturbations
 mpert=natom+7
 if (rf2_dkdk>0.or.rf2_dkde>0) mpert=natom+11

!Initialize the list of perturbations rfpert
 ABI_ALLOCATE(rfpert,(mpert))
 rfpert(:)=0
 if(rfphon==1)rfpert(dtset%rfatpol(1):dtset%rfatpol(2))=1

 if(rfddk==1)rfpert(natom+1)=1
 if(rfddk==2)rfpert(natom+6)=1

 if(rf2_dkdk/=0)rfpert(natom+10)=1
 if(rf2_dkde/=0)rfpert(natom+11)=1

 if(rfelfd==1.or.rfelfd==2)rfpert(natom+1)=1
 if(rfelfd==1.or.rfelfd==3)rfpert(natom+2)=1

 if(rfstrs==1.or.rfstrs==3)rfpert(natom+3)=1
 if(rfstrs==2.or.rfstrs==3)rfpert(natom+4)=1

 if(rfuser==1.or.rfuser==3)rfpert(natom+6)=1
 if(rfuser==2.or.rfuser==3)rfpert(natom+7)=1

 if(rfmagn==1) rfpert(natom+5)=1

 qeq0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2<1.d-14)

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell
 me = xmpi_comm_rank(spaceworld)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet_for_kg,dtset%istwfk,kg,&
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,npwarr,npwtot,&
& dtset%nsppol)

!Set up the Ylm for each k point
 ABI_ALLOCATE(ylm,(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm))
 if (rfstrs/=0.or.pawbec==1.or.pawpiezo==1.or.dtset%efmas>0) then
   ABI_ALLOCATE(ylmgr,(dtset%mpw*dtset%mkmem,9,psps%mpsang*psps%mpsang*psps%useylm))
 else
   ABI_ALLOCATE(ylmgr,(0,0,psps%useylm))
 end if
 if (psps%useylm==1) then
   option=0
   if (rfstrs/=0.or.pawbec==1.or.pawpiezo==1.or.dtset%efmas>0) option=2
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
 end if

 call timab(133,2,tsec)
 call timab(134,1,tsec)

!Open and read pseudopotential files
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,pawrad,pawtab,&
& psps,rprimd,comm_mpi=mpi_enreg%comm_cell)

 call timab(134,2,tsec)
 call timab(135,1,tsec)

!Initialize band structure datatype
 bstruct = ebands_from_dtset(dtset, npwarr)

!Initialize PAW atomic occupancies
 if (psps%usepaw==1) then
   ABI_DATATYPE_ALLOCATE(pawrhoij,(my_natom))
   call pawrhoij_nullify(pawrhoij)
   call initrhoij(dtset%pawcpxocc,dtset%lexexch,dtset%lpawu, &
&   my_natom,natom,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%ntypat,&
&   pawrhoij,dtset%pawspnorb,pawtab,cplex1,dtset%spinat,dtset%typat,&
&   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)
 else
   ABI_DATATYPE_ALLOCATE(pawrhoij,(0))
 end if

!Initialize header
 gscase=0
 call hdr_init(bstruct,codvsn,dtset,hdr,pawtab,gscase,psps,wvl%descr, &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)

!Update header, with evolving variables, when available
!Here, rprimd, xred and occ are available
 etot=hdr%etot ; fermie=hdr%fermie ; residm=hdr%residm
!If parallelism over atom, hdr is distributed
 call hdr_update(hdr,bantot,etot,fermie,&
& residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1), &
& comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab)

!Clean band structure datatype (should use it more in the future !)
 call ebands_free(bstruct)

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 mcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_MALLOC_OR_DIE(cg,(2,mcg), ierr)

 ABI_ALLOCATE(eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen0(:)=zero ; ask_accurate=1
 optorth=0

 hdr%rprimd=rprimd_for_kg ! We need the rprimd that was used to generate de G vectors
 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%fnamewffk,wvl)
 hdr%rprimd=rprimd

!Close wffgs, if it was ever opened (in inwffil)
 if (ireadwf0==1) then
   call WffClose(wffgs,ierr)
 end if

 if (psps%usepaw==1.and.ireadwf0==1) then
!  if parallelism, pawrhoij is distributed, hdr%pawrhoij is not
   call pawrhoij_copy(hdr%pawrhoij,pawrhoij,comm_atom=mpi_enreg%comm_atom,&
&   mpi_atmtab=mpi_enreg%my_atmtab)
 end if

 call timab(135,2,tsec)
 call timab(136,1,tsec)

!Report on eigen0 values   ! Should use prteigrs.F90
 write(message, '(a,a)' )
 call wrtout(std_out,ch10//' respfn : eigen0 array','COLL')
 nkpt_eff=dtset%nkpt
 if( (dtset%prtvol==0.or.dtset%prtvol==1.or.dtset%prtvol==2) .and. dtset%nkpt>nkpt_max ) nkpt_eff=nkpt_max
 band_index=0
 do isppol=1,dtset%nsppol
   do ikpt=1,dtset%nkpt
     nband_k=dtset%nband(ikpt+(isppol-1)*dtset%nkpt)
     if(ikpt<=nkpt_eff)then
       write(message, '(a,i2,a,i5)' )'  isppol=',isppol,', k point number',ikpt
       call wrtout(std_out,message,'COLL')
       do iband=1,nband_k,4
         write(message, '(a,4es16.6)')'  ',eigen0(iband+band_index:min(iband+3,nband_k)+band_index)
         call wrtout(std_out,message,'COLL')
       end do
     else if(ikpt==nkpt_eff+1)then
       write(message,'(a,a)' )'  respfn : prtvol=0, 1 or 2, stop printing eigen0.',ch10
       call wrtout(std_out,message,'COLL')
     end if
     band_index=band_index+nband_k
   end do
 end do

!Allocation for forces and atomic positions (should be taken away, also argument ... )
 ABI_ALLOCATE(grxc,(3,natom))

!Do symmetry stuff
 ABI_ALLOCATE(irrzon,(nfftot**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons,(2,nfftot**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(indsym,(4,dtset%nsym,natom))
 ABI_ALLOCATE(symrec,(3,3,dtset%nsym))
 irrzon=0;indsym=0;symrec=0;phnons=zero
!If the density is to be computed by mkrho, need irrzon and phnons
 iscf_eff=0;if(dtset%getden==0)iscf_eff=1
 call setsym(indsym,irrzon,iscf_eff,natom,&
& nfftot,ngfft,dtset%nspden,dtset%nsppol,dtset%nsym,&
& phnons,dtset%symafm,symrec,dtset%symrel,dtset%tnons,dtset%typat,xred)

!Symmetrize atomic coordinates over space group elements:
 call symmetrize_xred(indsym,natom,dtset%nsym,dtset%symrel,dtset%tnons,xred)

!Examine the symmetries of the q wavevector
 ABI_ALLOCATE(symq,(4,2,dtset%nsym))
 timrev=1

! By default use symmetries.
 use_sym = 1
 if (dtset%prtgkk == 1)then
   use_sym = 0
   call littlegroup_q(dtset%nsym,dtset%qptn,symq,symrec,dtset%symafm,timrev,prtvol=dtset%prtvol,use_sym=use_sym)
 else
   call littlegroup_q(dtset%nsym,dtset%qptn,symq,symrec,dtset%symafm,timrev,prtvol=dtset%prtvol)
 end if

!Generate an index table of atoms, in order for them to be used
!type after type.
 ABI_ALLOCATE(atindx,(natom))
 ABI_ALLOCATE(atindx1,(natom))
 ABI_ALLOCATE(nattyp,(ntypat))
 indx=1
 do itypat=1,ntypat
   nattyp(itypat)=0
   do iatom=1,natom
     if(dtset%typat(iatom)==itypat)then
       atindx(iatom)=indx
       atindx1(indx)=iatom
       indx=indx+1
       nattyp(itypat)=nattyp(itypat)+1
     end if
   end do
 end do

!Here allocation of GPU for fft calculations
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call alloc_hamilt_gpu(atindx1,dtset,gprimd,mpi_enreg,nattyp,npwarr,0,psps,dtset%use_gpu_cuda)
 end if
#endif

!Compute structure factor phases for current atomic pos:
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*natom))
 call getph(atindx,natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)

 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(atindx,natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

!Compute occupation numbers and fermi energy, in case occupation scheme is metallic.
 ABI_ALLOCATE(doccde,(dtset%mband*dtset%nkpt*dtset%nsppol))
 if( dtset%occopt>=3.and.dtset%occopt<=8 ) then

   call newocc(doccde,eigen0,entropy,fermie,dtset%spinmagntarget,dtset%mband,dtset%nband,&
&   dtset%nelect,dtset%nkpt,dtset%nspinor,dtset%nsppol,occ,dtset%occopt,dtset%prtvol,dtset%stmbias,&
&   dtset%tphysel,dtset%tsmear,dtset%wtk)

!  Update fermie and occ
   etot=hdr%etot ; residm=hdr%residm
   call hdr_update(hdr,bantot,etot,fermie,residm,rprimd,occ,pawrhoij,xred,dtset%amu_orig(:,1))

 else
!  doccde is irrelevant in this case
   doccde(:)=zero
 end if

!Recompute first large sphere cut-off gsqcut, without taking into account dilatmx
 ecutf=dtset%ecut
 if (psps%usepaw==1) then
   ecutf=dtset%pawecutdg
   call wrtout(std_out,ch10//' FFT (fine) grid used in SCF cycle:','COLL')
 end if

 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!PAW: 1- Initialize values for several arrays depending only on atomic data
!2- Check overlap
!3- Identify FFT points in spheres and compute g_l(r).Y_lm(r) (and exp(-i.q.r) if needed)
!4- Allocate PAW specific arrays
!5- Compute perturbed local potential inside spheres
!6- Eventually open temporary storage files
 if(psps%usepaw==1) then
!  1-Initialize values for several arrays depending only on atomic data
   gnt_option=1
   if (dtset%pawxcdev==2.or.dtset%rfphon/=0.or.dtset%rfstrs/=0.or.dtset%rfelfd==1.or.&
   dtset%rfelfd==3.or.dtset%rf2_dkde==1) gnt_option=2

   ! Test if we have to call pawinit
   call paw_gencond(Dtset,gnt_option,"test",call_pawinit)

   if (psp_gencond==1.or.call_pawinit) then
!    Some gen-cond have to be added...
     call timab(553,1,tsec)
     call pawinit(gnt_option,zero,zero,dtset%pawlcutd,dtset%pawlmix,&
&     psps%mpsang,dtset%pawnphi,dtset%nsym,dtset%pawntheta,&
&     pawang,pawrad,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%xclevel,dtset%usepotzero)
     call setsym_ylm(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&     rprimd,symrec,pawang%zarot)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit)

     call timab(553,2,tsec)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:psps%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:psps%ntypat)%has_nabla=2
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsym_ylm(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,rprimd,symrec,pawang%zarot)
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
&   is_dfpt,dtset%jpawu,dtset%lexexch,dtset%lpawu,ntypat,pawang,dtset%pawprtvol,pawrad,&
&   pawtab,dtset%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu)
   compch_fft=-1.d5;compch_sph=-1.d5
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   usecprj=dtset%pawusecp
!  2-Check overlap
   call chkpawovlp(natom,psps%ntypat,dtset%pawovlp,pawtab,rmet,dtset%typat,xred)
!  3-Identify FFT points in spheres and compute g_l(r).Y_lm(r) and exp(-i.q.r)
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(my_natom))
   if (my_natom>0) then
     call pawtab_get_lsize(pawtab,l_size_atm,my_natom,dtset%typat,&
&     mpi_atmtab=mpi_enreg%my_atmtab)
     call pawfgrtab_init(pawfgrtab,1,l_size_atm,pawrhoij(1)%nspden,dtset%typat,&
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     ABI_DEALLOCATE(l_size_atm)
   end if
   use_nhat_gga=(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)
   optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
   if (use_nhat_gga) then
     optgr1=dtset%pawstgylm
     if (rfphon==1) optgr2=1
   end if
   if (rfphon==1.or.rfstrs/=0.or.rfelfd==3.or.rf2_dkde==1) then ! LB2016-11-28 : Why not rfelfd==1?
     if (optgr1==0) optgr1=dtset%pawstgylm
     if (optgr2==0) optgr2=dtset%pawstgylm
     if (optrad==0.and.(.not.qeq0.or.rfstrs/=0.or.rfelfd==3.or.rf2_dkde==1)) optrad=1  ! LB2016-11-28 : Why not rfelfd==1?
   end if
   if (rfelfd==1.or.rfelfd==3.or.rf2_dkde==1) then
     if (optgr1==0) optgr1=dtset%pawstgylm
   end if
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfftf,psps%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred,&
&   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab )
!  Compute exp(iq.r) factors around the atoms
   if (.not.qeq0) then
     do iatom=1,my_natom
       iatom_tot=iatom; if (paral_atom) iatom_tot=mpi_enreg%my_atmtab(iatom)
       if (allocated(pawfgrtab(iatom)%expiqr)) then
         ABI_DEALLOCATE(pawfgrtab(iatom)%expiqr)
       end if
       ABI_ALLOCATE(pawfgrtab(iatom)%expiqr,(2,pawfgrtab(iatom)%nfgd))
       call pawexpiqr(pawfgrtab(iatom)%expiqr,gprimd,pawfgrtab(iatom)%nfgd,dtset%qptn,&
&       pawfgrtab(iatom)%rfgd,xred(:,iatom_tot))
       pawfgrtab(iatom)%expiqr_allocated=1
     end do
   end if
!  4-Allocate PAW specific arrays
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_kxc=0;nkxc1=0;cplex=1
   has_dijnd=0; if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_diju=merge(0,1,dtset%usepawu==0)
   if (rfphon/=0.or.rfelfd==1.or.rfelfd==3.or.rfstrs/=0.or.rf2_dkde/=0) then
     has_kxc=1
     call pawxc_get_nkxc(nkxc1,dtset%nspden,dtset%xclevel)
   end if
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,nkxc1,0,dtset%nspden,&
&   cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_kxc=has_kxc,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%pawspnorb,&
&   natom,dtset%ntypat,dtset%typat,pawtab,has_dij=1,has_dijhartree=1,has_dijnd=has_dijnd,&
&   has_dijso=1,has_dijU=has_diju,has_pawu_occ=1,has_exexch_pot=1,nucdipmom=dtset%nucdipmom,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

 else ! PAW vs NCPP
   usexcnhat=0;usecprj=0
   use_nhat_gga=.false.
   ABI_DATATYPE_ALLOCATE(paw_an,(0))
   ABI_DATATYPE_ALLOCATE(paw_ij,(0))
   ABI_DATATYPE_ALLOCATE(pawfgrtab,(0))
 end if

 ABI_ALLOCATE(rhog,(2,nfftf))
 ABI_ALLOCATE(rhor,(nfftf,dtset%nspden))

!Read ground-state charge density from diskfile in case getden /= 0
!or compute it from wfs that were read previously : rhor as well as rhog

 if (dtset%getden /= 0 .or. dtset%irdden /= 0) then
   ! Read rho1(r) from a disk file and broadcast data.
   ! This part is not compatible with MPI-FFT (note single_proc=.True. below)

   rdwr=1;rdwrpaw=psps%usepaw;if(ireadwf0/=0) rdwrpaw=0
   if (rdwrpaw/=0) then
     ABI_DATATYPE_ALLOCATE(pawrhoij_read,(natom))
     call pawrhoij_nullify(pawrhoij_read)
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,nspden_rhoij=nspden_rhoij,&
&              nspden=dtset%nspden,spnorb=dtset%pawspnorb,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_read,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&                        dtset%nsppol,dtset%typat,pawtab=pawtab)
   else
     ABI_DATATYPE_ALLOCATE(pawrhoij_read,(0))
   end if

!    MT july 2013: Should we read rhoij from the density file ?
   call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rhor, &
   hdr_den, pawrhoij_read, spaceworld, check_hdr=hdr)
   etotal = hdr_den%etot; call hdr_free(hdr_den)

   if (rdwrpaw/=0) then
     call pawrhoij_bcast(pawrhoij_read,hdr%pawrhoij,0,spaceworld)
     call pawrhoij_free(pawrhoij_read)
   end if
   ABI_DATATYPE_DEALLOCATE(pawrhoij_read)

!  Compute up+down rho(G) by fft
   ABI_ALLOCATE(work,(nfftf))
   work(:)=rhor(:,1)
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,1,ngfftf,0)
   ABI_DEALLOCATE(work)

 else
   izero=0
!  Obtain the charge density from read wfs
!  Be careful: in PAW, compensation density has to be added !
   tim_mkrho=4
   paw_dmft%use_sc_dmft=0 ! respfn with dmft not implemented
   paw_dmft%use_dmft=0 ! respfn with dmft not implemented
   if (psps%usepaw==1) then
     ABI_ALLOCATE(rhowfg,(2,dtset%nfft))
     ABI_ALLOCATE(rhowfr,(dtset%nfft,dtset%nspden))
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhowfg,rhowfr,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
     call transgrid(1,mpi_enreg,dtset%nspden,+1,1,1,dtset%paral_kgb,pawfgr,rhowfg,rhog,rhowfr,rhor)
     ABI_DEALLOCATE(rhowfg)
     ABI_DEALLOCATE(rhowfr)
   else
     call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,&
&     mpi_enreg,npwarr,occ,paw_dmft,phnons,rhog,rhor,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs)
   end if
 end if ! getden

!In PAW, compensation density has eventually to be added
 nhatgrdim=0;nhatdim=0
 ABI_ALLOCATE(nhatgr,(0,0,0))
 if (psps%usepaw==1.and. ((usexcnhat==0).or.(dtset%getden==0).or.dtset%xclevel==2)) then
   nhatdim=1
   ABI_ALLOCATE(nhat,(nfftf,dtset%nspden))
   call timab(558,1,tsec)
   nhatgrdim=0;if (dtset%xclevel==2.and.dtset%pawnhatxc>0) nhatgrdim=usexcnhat
   ider=2*nhatgrdim
   if (nhatgrdim>0)  then
     ABI_DEALLOCATE(nhatgr)
     ABI_ALLOCATE(nhatgr,(nfftf,dtset%nspden,3))
   end if
   izero=0;cplex=1;ipert=0;idir=0;qphon(:)=zero
   call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,natom,&
&   nfftf,ngfftf,nhatgrdim,dtset%nspden,psps%ntypat,pawang,pawfgrtab,&
&   nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred, &
&   mpi_atmtab=mpi_enreg%my_atmtab, comm_atom=mpi_enreg%comm_atom)
   if (dtset%getden==0) then
     rhor(:,:)=rhor(:,:)+nhat(:,:)
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,1,ngfftf,0)
   end if
   call timab(558,2,tsec)
 else
   ABI_ALLOCATE(nhat,(0,0))
 end if

!The GS irrzon and phnons were only needed to symmetrize the GS density
 ABI_DEALLOCATE(irrzon)
 ABI_DEALLOCATE(phnons)

!jmb 2012 write(std_out,'(a)')' ' ! needed to make ibm6_xlf12 pass tests. No idea why this works. JWZ 5 Sept 2011
!Will compute now the total potential

!Compute local ionic pseudopotential vpsp and core electron density xccc3d:
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(vpsp,(nfftf))

!Determine by which method the local ionic potential and/or
! the pseudo core charge density have to be computed
!Local ionic potential:
! Method 1: PAW ; Method 2: Norm-conserving PP
 vloc_method=1;if (psps%usepaw==0) vloc_method=2
!Pseudo core charge density:
! Method 1: PAW, nc_xccc_gspace ; Method 2: Norm-conserving PP
 coredens_method=1;if (psps%usepaw==0) coredens_method=2
 if (psps%nc_xccc_gspace==1) coredens_method=1
 if (psps%nc_xccc_gspace==0) coredens_method=2

!Local ionic potential and/or pseudo core charge by method 1
 if (vloc_method==1.or.coredens_method==1) then
   call timab(562,1,tsec)
   optv=0;if (vloc_method==1) optv=1
   optn=0;if (coredens_method==1) optn=n3xccc/nfftf
   optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optn2=1
   call atm2fft(atindx1,xccc3d,vpsp,dum_dyfrn,dum_dyfrv,dum_eltfrxc,dum_gauss,gmet,gprimd,&
&   dum_grn,dum_grv,gsqcut,mgfftf,psps%mqgrid_vl,natom,nattyp,nfftf,ngfftf,&
&   ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1df,psps%qgrid_vl,&
&   dtset%qprtrb,dum_rhog,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,dum_vg,dum_vg,dum_vg,dtset%vprtrb,psps%vlspl)
   call timab(562,2,tsec)
 end if

!Local ionic potential by method 2
 if (vloc_method==2) then
   option=1
   ABI_ALLOCATE(dyfrlo_indx,(3,3,natom))
   ABI_ALLOCATE(grtn_indx,(3,natom))
   call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,&
&   grtn_indx,gsqcut,dummy6,mgfftf,mpi_enreg,natom,nattyp,&
&   nfftf,ngfftf,dtset%nspden,ntypat,option,pawtab,ph1df,psps,&
&   dtset%qprtrb,rhog,rhor,rprimd,ucvol,dtset%vprtrb,vpsp,wvl%descr,wvl%den,xred)
   ABI_DEALLOCATE(dyfrlo_indx)
   ABI_DEALLOCATE(grtn_indx)
 end if

!Pseudo core electron density by method 2
 if (coredens_method==2.and.psps%n1xccc/=0) then
   option=1
   ABI_ALLOCATE(dyfrx2,(3,3,natom))
   ABI_ALLOCATE(vxc,(0,0)) ! dummy
   call mkcore(dummy6,dyfrx2,grxc,mpi_enreg,natom,nfftf,dtset%nspden,ntypat,&
&   ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&   psps%xcccrc,psps%xccc1d,xccc3d,xred)
   ABI_DEALLOCATE(dyfrx2)
   ABI_DEALLOCATE(vxc) ! dummy
 end if

!Set up hartree and xc potential. Compute kxc here.
 ABI_ALLOCATE(vhartr,(nfftf))

 call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfftf,ngfftf,rhog,rprimd,vhartr)

 option=2 ; nk3xc=1
 nkxc=2*min(dtset%nspden,2)-1;if(dtset%xclevel==2)nkxc=12*min(dtset%nspden,2)-5
 call check_kxc(dtset%ixc,dtset%optdriver)
 ABI_ALLOCATE(kxc,(nfftf,nkxc))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))

 _IBM6("Before rhotoxc")

 call xcdata_init(xcdata,dtset=dtset)
 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
 call rhotoxc(enxc,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,non_magnetic_xc,n3xccc,option,rhor,&
& rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,vhartr=vhartr)

!Compute local + Hxc potential, and subtract mean potential.
 ABI_ALLOCATE(vtrial,(nfftf,dtset%nspden))
 do ispden=1,min(dtset%nspden,2)
   do ifft=1,nfftf
     vtrial(ifft,ispden)=vhartr(ifft)+vxc(ifft,ispden)+vpsp(ifft)
   end do
 end do
 if (dtset%nspden==4) then
   do ispden=3,4
     do ifft=1,nfftf
       vtrial(ifft,ispden)=vxc(ifft,ispden)
     end do
   end do
 end if
 ABI_DEALLOCATE(vhartr)


 if(dtset%prtvol==-level)then
   call wrtout(std_out,' respfn: ground-state density and potential set up.','COLL')
 end if

!PAW: compute Dij quantities (psp strengths)
 if (psps%usepaw==1)then
   cplex=1;ipert=0;option=1
   nzlmopt=0;if (dtset%pawnzlm>0) nzlmopt=-1

   call pawdenpot(compch_sph,epaw,epawdc,ipert,dtset%ixc,my_natom,natom,dtset%nspden,&
&   ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,&
&   pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&   dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp, &
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

   call timab(561,1,tsec)
   call pawdij(cplex,dtset%enunit,gprimd,ipert,my_natom,natom,nfftf,nfftotf,&
&   dtset%nspden,ntypat,paw_an,paw_ij,pawang,pawfgrtab,dtset%pawprtvol,&
&   pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,k0,&
&   dtset%spnorbscl,ucvol,dtset%charge,vtrial,vxc,xred,nucdipmom=dtset%nucdipmom,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call symdij(gprimd,indsym,ipert,my_natom,natom,dtset%nsym,ntypat,0,&
&   paw_ij,pawang,dtset%pawprtvol,pawtab,rprimd,dtset%symafm,symrec,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call timab(561,2,tsec)

 end if

!-----2. Frozen-wavefunctions and Ewald(q=0) parts of 2DTE

 dyfr_nondiag=0;if (psps%usepaw==1.and.rfphon==1) dyfr_nondiag=1
 dyfr_cplex=1;if (psps%usepaw==1.and.rfphon==1.and.(.not.qeq0)) dyfr_cplex=2
 ABI_ALLOCATE(dyew,(2,3,natom,3,natom))
 ABI_ALLOCATE(dyewq0,(3,3,natom))
 ABI_ALLOCATE(dyfrlo,(3,3,natom))
 ABI_ALLOCATE(dyfrx2,(3,3,natom))
 ABI_ALLOCATE(dyfrnl,(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag))
 ABI_ALLOCATE(dyfrwf,(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag))
 ABI_ALLOCATE(becfrnl,(3,natom,3*pawbec))
 ABI_ALLOCATE(piezofrnl,(6,3*pawpiezo))
 ABI_ALLOCATE(dyvdw,(2,3,natom,3,natom*usevdw))
 dyew(:,:,:,:,:)=zero
 dyewq0(:,:,:)=zero
 dyfrnl(:,:,:,:,:)=zero
 dyfrwf(:,:,:,:,:)=zero
 dyfrlo(:,:,:)=zero
 dyfrx2(:,:,:)=zero
 if (usevdw==1) dyvdw(:,:,:,:,:)=zero
 if (pawbec==1) becfrnl(:,:,:)=zero
 if (pawpiezo==1) piezofrnl(:,:)=zero

 ABI_ALLOCATE(eltcore,(6,6))
 ABI_ALLOCATE(elteew,(6+3*natom,6))
 ABI_ALLOCATE(eltfrhar,(6,6))
 ABI_ALLOCATE(eltfrnl,(6+3*natom,6))
 ABI_ALLOCATE(eltfrloc,(6+3*natom,6))
 ABI_ALLOCATE(eltfrkin,(6,6))
 ABI_ALLOCATE(eltfrxc,(6+3*natom,6))
 ABI_ALLOCATE(eltvdw,(6+3*natom,6*usevdw))
 eltcore(:,:)=zero
 elteew(:,:)=zero
 eltfrnl(:,:)=zero
 eltfrloc(:,:)=zero
 eltfrkin(:,:)=zero
 eltfrhar(:,:)=zero
 eltfrxc(:,:)=zero
 if (usevdw==1) eltvdw(:,:)=zero

!Section common to all perturbations
!Compute the nonlocal part of the elastic tensor and/or dynamical matrix
 if (rfstrs/=0.or.rfphon==1.or.dtset%efmas>0.or.pawbec==1.or.pawpiezo==1)then
   call d2frnl(becfrnl,cg,dtfil,dtset,dyfrnl,dyfr_cplex,&
&   dyfr_nondiag,efmasdeg,efmasval,eigen0,eltfrnl,gsqcut,has_allddk,indsym,kg,mgfftf,&
&   mpi_enreg,psps%mpsang,my_natom,natom,nfftf,ngfft,ngfftf,&
&   npwarr,occ,paw_ij,pawang,pawbec,pawfgrtab,pawpiezo,pawrad,&
&   pawrhoij,pawtab,ph1d,ph1df,piezofrnl,psps,rprimd,rfphon,&
&   rfstrs,symrec,vtrial,vxc,xred,ylm,ylmgr)
 end if

!No more need of these local derivatives
 if (rfphon==1.and.psps%usepaw==1.and.(.not.use_nhat_gga)) then
   do iatom=1,my_natom
     if (allocated(pawfgrtab(iatom)%gylmgr2)) then
       ABI_DEALLOCATE(pawfgrtab(iatom)%gylmgr2)
     end if
     pawfgrtab(iatom)%gylmgr2_allocated=0
   end do
 end if

!Section for the atomic displacement/electric field perturbations
 if (rfphon==1) then

!  Compute the local of the dynamical matrix
!  dyfrnl has not yet been symmetrized, but will be in the next routine
   call dfpt_dyfro(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrx2,dyfr_cplex,dyfr_nondiag,&
&   gmet,gprimd,gsqcut,indsym,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&   natom,nattyp, nfftf,ngfftf,dtset%nspden,dtset%nsym,ntypat,&
&   psps%n1xccc,n3xccc,psps,pawtab,ph1df,psps%qgrid_vl,&
&   dtset%qptn,rhog,rprimd,symq,symrec,dtset%typat,ucvol,&
&   psps%usepaw,psps%vlspl,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

   _IBM6("Before dfpt_ewald")

!  Compute Ewald (q=0) contribution
   sumg0=0;qphon(:)=zero
   call dfpt_ewald(dyew,gmet,my_natom,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   option=1
   call q0dy3_calc(natom,dyewq0,dyew,option)
!  Calculate the DFT-D2 vdw part of the dynamical matrix
   if (usevdw==1.and.dtset%vdw_xc==5) then
     call vdw_dftd2(evdw,dtset%ixc,natom,ntypat,1,dtset%typat,rprimd,dtset%vdw_tol,&
&     xred,psps%znuclpsp,dyn_vdw_dftd2=dyvdw,qphon=dtset%qptn)
   end if
! Calculate the DFT-D3/D3(BJ) vdw part of the dynamical matrix
   if (usevdw==1.and.(dtset%vdw_xc==6.or.dtset%vdw_xc==7)) then
     call vdw_dftd3(evdw,dtset%ixc,natom,ntypat,1,dtset%typat,rprimd,&
&     dtset%vdw_xc,dtset%vdw_tol,dtset%vdw_tol_3bt,xred,psps%znuclpsp,&
&     dyn_vdw_dftd3=dyvdw,qphon=dtset%qptn)
   end if
!The frozen-wavefunction part of the dynamical matrix is now:
!  d2frnl:  non-local contribution
!  dyfrlo:  local contribution
!  dyfrx2:  2nd order xc core correction contribution
!  dyew  : Ewald contribution
!  dyvdw : vdw DFT-D contribution
!  dyfrwf:  all contributions
!  In case of PAW, it misses a term coming from the perturbed overlap operator
 end if

!Section for the strain perturbation
 if(rfstrs/=0) then

!  Verify that k-point set has full space-group symmetry; otherwise exit
   timrev=1
   if (symkchk(dtset%kptns,dtset%nkpt,dtset%nsym,symrec,timrev,message) /= 0) then
     MSG_ERROR(message)
   end if

!  Calculate the kinetic part of the elastic tensor
   call dfpt_eltfrkin(cg,eltfrkin,dtset%ecut,dtset%ecutsm,dtset%effmass_free,&
&   dtset%istwfk,kg,dtset%kptns,dtset%mband,dtset%mgfft,dtset%mkmem,mpi_enreg,&
&   dtset%mpw,dtset%nband,dtset%nkpt,ngfft,npwarr,&
&   dtset%nspinor,dtset%nsppol,occ,rprimd,dtset%wtk)

!  Calculate the hartree part of the elastic tensor
   call dfpt_eltfrhar(eltfrhar,rprimd,gsqcut,mpi_enreg,nfftf,ngfftf,rhog)

!  Calculate the xc part of the elastic tensor
   call dfpt_eltfrxc(atindx,dtset,eltfrxc,enxc,gsqcut,kxc,mpi_enreg,mgfftf,&
&   nattyp,nfftf,ngfftf,ngfftf,nhat,nkxc,n3xccc,pawtab,ph1df,psps,rhor,rprimd,&
&   usexcnhat,vxc,xccc3d,xred)

!  Calculate the local potential part of the elastic tensor
   call dfpt_eltfrloc(atindx,eltfrloc,gmet,gprimd,gsqcut,mgfftf,mpi_enreg,psps%mqgrid_vl,&
&   natom,nattyp,nfftf,ngfftf,ntypat,ph1df,psps%qgrid_vl,rhog,psps%vlspl)

!  Calculate the Ewald part of the elastic tensor
   call elt_ewald(elteew,gmet,gprimd,my_natom,natom,ntypat,rmet,rprimd,&
&   dtset%typat,ucvol,xred,psps%ziontypat,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

!  Calculate the DFT-D2 vdw part of the elastic tensor
   if (usevdw==1.and.dtset%vdw_xc==5) then
     option=1; if (rfphon==1) option=0
     call vdw_dftd2(evdw,dtset%ixc,natom,ntypat,option,dtset%typat,rprimd,dtset%vdw_tol,&
&     xred,psps%znuclpsp,elt_vdw_dftd2=eltvdw)
   end if
!  Calculate the DFT-D3/D3(BJ) vdw part of the elastic tensor
   if (usevdw==1.and.(dtset%vdw_xc==6.or.dtset%vdw_xc==7)) then
     option=1; if (rfphon==1) option=0
     call vdw_dftd3(evdw,dtset%ixc,natom,ntypat,option,dtset%typat,rprimd,&
&     dtset%vdw_xc,dtset%vdw_tol,dtset%vdw_tol_3bt,xred,psps%znuclpsp,&
&     elt_vdw_dftd3=eltvdw)
   end if
!  Calculate the psp core energy part of elastic tensor (trivial)
   eltcore(1:3,1:3)=ecore/ucvol

!The frozen-wavefunction part of the elastic tensor is now:
!  eltfrnl:  non-local contribution
!  eltfrloc: local contribution
!  eltfrkin: kinetic contribution
!  eltfrhar: Hartree contribution
!  eltfrx:   XC contribution
!  eltcore:  psps core contribution
!  elteew:   Ewald contribution
!  eltvdw:   vdw DFT-D contribution
!  In case of PAW, it misses a term coming from the perturbed overlap operator
 end if

 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(xccc3d)

 if(dtset%prtvol==-level)then
   call wrtout(std_out,' respfn: frozen wavef. and Ewald(q=0) part of 2DTE done.','COLL')
 end if

 call timab(136,2,tsec)

!-----3. Initialisation of 1st response, taking into account the q vector.

 call timab(137,1,tsec)

 write(message,'(3a)')ch10,' ==>  initialize data related to q vector <== ',ch10
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 qphon(:)=dtset%qptn(:)
 sumg0=1

!Treat the case of q=0 or q too close to 0
 qzero=0
 if(qeq0)then
   qphon(:)=zero
   write(message,'(3a)')&
&   ' respfn : the norm of the phonon wavelength (as input) was small (<1.d-7).',ch10,&
&   '  q has been set exactly to (0 0 0)'
   call wrtout(std_out,message,'COLL')
   sumg0=0
   qzero=1
 else
   if(rfelfd/=0 .or. rfstrs/=0 .or. rfddk /= 0  .or. rf2_dkdk /= 0 .or. rf2_dkde /= 0) then
!    Temporarily, ...
     write(message, '(a,a,a,3es16.6,a,a,5(a,i2),a,a,a)' )ch10,&
&     'The treatment of non-zero wavevector q is restricted to phonons.',&
&     'However, the input normalized qpt is',qphon(:),',',ch10,&
&     'while rfelfd=',rfelfd,', rfddk=',rfddk,', rf2_dkdk=',rf2_dkdk,', rf2_dkde=',rf2_dkde,&
&     ' and rfstrs=',rfstrs,'.',ch10,&
&     'Action: change qpt, or rfelfd, or rfstrs in the input file.'
     MSG_ERROR(message)
   else if(rfasr.eq.2)then
     write(message,'(2a)')ch10,' rfasr=2 not allowed with q/=0 => rfasr was reset to 0.'
     MSG_WARNING(message)
     rfasr=0
   end if
 end if

 _IBM6("Before irreducible_set_pert")

!Determine the symmetrical perturbations
 ABI_ALLOCATE(pertsy,(3,natom+6))
 call irreducible_set_pert(indsym,natom+6,natom,dtset%nsym,pertsy,rfdir,rfpert,symq,symrec,dtset%symrel)
 write(message,'(a)') ' The list of irreducible perturbations for this q vector is:'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 ii=1
 do ipert=1,natom+6
   do idir=1,3
     if(rfpert(ipert)==1.and.rfdir(idir)==1)then
       if( pertsy(idir,ipert)==1 )then
         write(message, '(i5,a,i2,a,i4)' )ii,')    idir=',idir,'    ipert=',ipert
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')
         ii=ii+1
       end if
     end if
   end do
 end do

!test if the user left default rfdir 0 0 0
 if (ii==1 .and. rf2_dkdk==0 .and. rf2_dkde==0) then
   write(message,'(5a)')ch10,&
&   ' WARNING: no perturbations to be done at this q-point.',ch10,&
&   ' You may have forgotten to set the rfdir or rfatpol variables. Continuing normally.',ch10
   call wrtout(ab_out,message,'COLL')
   MSG_WARNING(message)
 end if

 if (dtset%prepanl==1.and.(rf2_dkdk/=0 .or. rf2_dkde/=0)) then
   ABI_ALLOCATE(rfpert_nl,(3,natom+2,3,natom+2,3,natom+2))
   rfpert_nl = 0
   rfpert_nl(:,natom+2,:,natom+2,:,natom+2) = 1
   rfpert_nl(:,1:natom,:,natom+2,:,natom+2) = 1
   rfpert_nl(:,natom+2,:,1:natom,:,natom+2) = 1
   rfpert_nl(:,natom+2,:,natom+2,:,1:natom) = 1
   call sytens(indsym,natom+2,natom,dtset%nsym,rfpert_nl,symrec,dtset%symrel)
   write(message, '(a,a,a,a,a)' ) ch10, &
&   ' The list of irreducible elements of the Raman and non-linear',&
&   ch10,' optical susceptibility tensors is:',ch10
   call wrtout(std_out,message,'COLL')

   write(message,'(12x,a)')&
&   'i1pert  i1dir   i2pert  i2dir   i3pert  i3dir'
   call wrtout(std_out,message,'COLL')
   n1 = 0
   rf2_dirs_from_rfpert_nl(:,:) = 0
   do i1pert = 1, natom + 2
     do i1dir = 1, 3
       do i2pert = 1, natom + 2
         do i2dir = 1, 3
           do i3pert = 1, natom + 2
             do i3dir = 1,3
               if (rfpert_nl(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
                 n1 = n1 + 1
                 write(message,'(2x,i4,a,6(5x,i3))') n1,')', &
&                 i1pert,i1dir,i2pert,i2dir,i3pert,i3dir
                 call wrtout(std_out,message,'COLL')
                 if (i2pert==natom+2) then
                   if (i3pert==natom+2) then
                     rf2_dirs_from_rfpert_nl(i3dir,i2dir) = 1
                   else if (i1pert==natom+2) then
                     rf2_dirs_from_rfpert_nl(i1dir,i2dir) = 1
                   end if
                 end if
               end if
             end do
           end do
         end do
       end do
     end do
   end do
   write(message,'(a,a)') ch10,ch10
   call wrtout(std_out,message,'COLL')

   write(message,'(a)') 'rf2_dirs_from_rfpert_nl :'
   call wrtout(std_out,message,'COLL')
   do i1dir = 1, 3
     do i2dir = 1, 3
       write(message,'(3(a,i1))') ' ',i1dir,' ',i2dir,' : ',rf2_dirs_from_rfpert_nl(i1dir,i2dir)
       call wrtout(std_out,message,'COLL')
     end do
   end do
 end if

!Contribution to the dynamical matrix from ion-ion energy
 if(rfphon==1)then
   call dfpt_ewald(dyew,gmet,my_natom,natom,qphon,rmet,sumg0,dtset%typat,ucvol,xred,psps%ziontypat, &
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call q0dy3_apply(natom,dyewq0,dyew)
 end if

!1-order contribution of the xc core correction to the dynamical matrix
 ABI_ALLOCATE(dyfrx1,(2,3,natom,3,natom))
 dyfrx1(:,:,:,:,:)=zero
 if(rfphon==1.and.psps%n1xccc/=0)then
   ABI_ALLOCATE(blkflgfrx1,(3,natom,3,natom))
!FR non-collinear magnetism
   if (dtset%nspden==4) then
     call dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,dtset%ixc,kxc,mgfftf,mpert,mpi_enreg,&
&     psps%mqgrid_vl,natom,nfftf,ngfftf,nkxc,non_magnetic_xc,dtset%nspden,&
&     ntypat,psps%n1xccc,psps,pawtab,ph1df,psps%qgrid_vl,qphon,&
&     rfdir,rfpert,rprimd,timrev,dtset%typat,ucvol,psps%usepaw,psps%xcccrc,psps%xccc1d,xred,rhor=rhor,vxc=vxc)
   else
     call dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,dtset%ixc,kxc,mgfftf,mpert,mpi_enreg,&
&     psps%mqgrid_vl,natom,nfftf,ngfftf,nkxc,non_magnetic_xc,dtset%nspden,&
&     ntypat,psps%n1xccc,psps,pawtab,ph1df,psps%qgrid_vl,qphon,&
&     rfdir,rfpert,rprimd,timrev,dtset%typat,ucvol,psps%usepaw,psps%xcccrc,psps%xccc1d,xred)
   end if
 end if

!Deallocate the arrays that were needed only for the frozen wavefunction part
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(npwarr)
 if(xmpi_paral==1) then
   ABI_DEALLOCATE(mpi_enreg%proc_distrb)
 end if

 ABI_ALLOCATE(blkflg,(3,mpert,3,mpert))
 ABI_ALLOCATE(d2eig0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2k0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2lo,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2loc0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nfr,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl0,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2nl1,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2vn,(2,3,mpert,3,mpert))
 ABI_ALLOCATE(d2ovl,(2,3,mpert,3,mpert*psps%usepaw))
 blkflg(:,:,:,:)=0
 d2eig0(:,:,:,:,:)=zero ; d2k0(:,:,:,:,:)=zero
 d2lo(:,:,:,:,:)=zero   ; d2loc0(:,:,:,:,:)=zero
 d2nfr(:,:,:,:,:)=zero  ; d2nl(:,:,:,:,:)=zero
 d2nl0(:,:,:,:,:)=zero  ; d2nl1(:,:,:,:,:)=zero
 d2vn(:,:,:,:,:)=zero
 if (psps%usepaw==1) d2ovl(:,:,:,:,:)=zero

 prtbbb=dtset%prtbbb
 ABI_ALLOCATE(d2bbb,(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 ABI_ALLOCATE(d2cart_bbb,(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb))
 if(prtbbb==1)then
   d2cart_bbb(:,:,:,:,:,:)=zero
   d2bbb(:,:,:,:,:,:)=zero
 end if

 dim_eig2nkq = 0
 if(dtset%ieig2rf /= 0) dim_eig2nkq = 1
 ABI_ALLOCATE(eig2nkq,(2,dtset%mband*dtset%nsppol,dtset%nkpt,3,natom,3,natom*dim_eig2nkq))
 dim_eigbrd=0
 if(dtset%ieig2rf /= 0 .and. dtset%smdelta>0 ) dim_eigbrd = 1
 ABI_ALLOCATE(eigbrd,(2,dtset%mband*dtset%nsppol,dtset%nkpt,3,natom,3,natom*dim_eigbrd))

 call timab(137,2,tsec)


!Check whether exiting was required by the user.
!If found then do not start minimization steps
!At this first call to exit_check, initialize cpus
 cpus=dtset%cpus
 if(abs(cpus)>1.0d-5)cpus=cpus+cpui
 openexit=1; if(dtset%chkexit==0) openexit=0
 call exit_check(cpus,dtfil%filnam_ds(1),iexit,ab_out,mpi_enreg%comm_cell,openexit)

!TEMPORARY: for testing purpose only
! if (rfstrs/=0.and.dtset%usepaw==1) iexit=1

 _IBM6("Before dfpt_looppert")

 if (iexit==0) then
!  #######################################################################
   write(message,'(a,80a)')ch10,('=',mu=1,80)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')

   ddkfil(:)=0

   ! MGNAG: WHY THIS?  all these variables are declared as optional pointers in dfpt_looppert!
   ! but they are allocated here so why using pointers! Moreover OPTIONAL arguments MUST
   ! be passed by keyword for better clarity and robustness!
   ! People should learn how to program Fort90 before being allowed to change the code!
   ! v5[26] crashes in dfpt_looppert
   ! The best solution would be using a datatype to gather the results!
   ABI_ALLOCATE(clflg,(3,mpert))
   if(dtset%ieig2rf > 0) then
     ABI_ALLOCATE(eigen0_pert,(dtset%mband*dtset%nkpt*dtset%nsppol))
     ABI_ALLOCATE(eigenq_pert,(dtset%mband*dtset%nkpt*dtset%nsppol))
     ABI_ALLOCATE(occ_rbz_pert,(dtset%mband*dtset%nkpt*dtset%nsppol))
   end if
   if(dtset%efmas > 0) then
     ABI_ALLOCATE(eigen0_pert,(dtset%mband*dtset%nkpt*dtset%nsppol))
   end if
!  Note that kg, cg, eigen0, mpw and npwarr are NOT passed to dfpt_looppert :
!  they will be reinitialized for each perturbation, with an eventually
!  reduced set of k point, thanks to the use of symmetry operations.
   call dfpt_looppert(atindx,blkflg,codvsn,cpus,dim_eigbrd,dim_eig2nkq,doccde,&
&   ddkfil,dtfil,dtset,dyew,dyfrlo,dyfrnl,dyfrx1,dyfrx2,dyvdw,dyfr_cplex,dyfr_nondiag,&
&   d2bbb,d2lo,d2nl,d2ovl,efmasdeg,efmasval,eigbrd,eig2nkq,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
&   etotal,fermie,iexit,indsym,kxc,&
&   dtset%mkmem,mkqmem,mk1mem,mpert,mpi_enreg,my_natom,nattyp,&
&   nfftf,nhat,dtset%nkpt,nkxc,dtset%nspden,dtset%nsym,occ,&
&   paw_an,paw_ij,pawang,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&   pertsy,prtbbb,psps,rfpert,rf2_dirs_from_rfpert_nl,rhog,rhor,symq,symrec,timrev,&
&   usecprj,usevdw,vtrial,vxc,vxcavg,xred,clflg,occ_rbz_pert,eigen0_pert,eigenq_pert,&
&   eigen1_pert,nkpt_rbz,eigenq_fine,hdr_fine,hdr0)

!  #####################################################################
 end if

 call timab(138,1,tsec)

 write(message, '(80a,a,a,a,a)' ) ('=',mu=1,80),ch10,ch10,&
& ' ---- first-order wavefunction calculations are completed ----',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 ABI_DEALLOCATE(vxc)

 if (dtset%prepanl==1.and.(rf2_dkdk/=0 .or. rf2_dkde/=0)) then
   ABI_DEALLOCATE(rfpert_nl)
 end if

!Output of the localization tensor
 if ( rfpert(natom+1) /= 0 .and. (me == 0) .and. dtset%occopt<=2) then
   call wrtloctens(blkflg,d2bbb,d2nl,dtset%mband,mpert,natom,dtset%prtbbb,rprimd,psps%usepaw)
 end if

!The perturbation  natom+1 was only an auxiliary perturbation,
!needed to construct the electric field response, so its flag is now set to 0.
!rfpert(natom+1)=0

!Were 2DTE computed ?
 if(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0 .or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and. rfuser==0 .and. rfmagn==0)then

   write(message,'(a,a)' )ch10,' respfn : d/dk was computed, but no 2DTE, so no DDB output.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  If 2DTE were computed, only one processor must output them and compute
!  frequencies.
 else if(me==0)then

   write(message,'(a,a)' )ch10,' ==> Compute Derivative Database <== '
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  In the RESPFN code, dfpt_nstdy and stady3 were called here
   d2nfr(:,:,:,:,:)=d2lo(:,:,:,:,:)+d2nl(:,:,:,:,:)
   if (psps%usepaw==1) d2nfr(:,:,:,:,:)=d2nfr(:,:,:,:,:)+d2ovl(:,:,:,:,:)

   zero_by_symm=1
   if(dtset%rfmeth<0)zero_by_symm=0

!  In case of bbb decomposition
   if(prtbbb==1)then
     ABI_ALLOCATE(blkflg1,(3,mpert,3,mpert))
     ABI_ALLOCATE(blkflg2,(3,mpert,3,mpert))
     blkflg2(:,:,:,:) = blkflg(:,:,:,:)
     do ipert = 1, mpert
       do ipert2 = 1, mpert
         if ((ipert /= natom + 2).and.(ipert>natom).and.(ipert2/=natom+2)) then
           blkflg2(:,ipert2,:,ipert) = 0
         end if
       end do
     end do
     ABI_ALLOCATE(d2tmp,(2,3,mpert,3,mpert))
     do iband = 1,dtset%mband
       d2tmp(:,:,:,:,:)=zero
       blkflg1(:,:,:,:) = blkflg2(:,:,:,:)
       d2tmp(:,:,natom+2,:,:) = d2bbb(:,:,:,:,iband,iband)
       call d2sym3(blkflg1,d2tmp,indsym,mpert,natom,dtset%nsym,qphon,symq,&
&       symrec,dtset%symrel,timrev,zero_by_symm)
       d2bbb(:,:,:,:,iband,iband) = d2tmp(:,:,natom+2,:,:)
     end do
     ABI_DEALLOCATE(blkflg1)
     ABI_DEALLOCATE(blkflg2)
     ABI_DEALLOCATE(d2tmp)
   end if

!  Complete the d2nfr matrix by symmetrization of the existing elements
   !write(std_out,*)"blkflg before d2sym3: ", blkflg
   call d2sym3(blkflg,d2nfr,indsym,mpert,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev,zero_by_symm)
   !write(std_out,*)"blkflg after d2sym3: ", blkflg

   if(rfphon==1.and.psps%n1xccc/=0)then
!    Complete the dyfrx1 matrix by symmetrization of the existing elements
     call d2sym3(blkflgfrx1,dyfrx1,indsym,natom,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev,zero_by_symm)
   end if

!  Note that d2sym3 usually complete the 2nd-order matrix
!  with elements that are zero by symmetry, automatically,
!  unless it has been explicitly asked not to do so.
!  blkflg is then set to 1 for these matrix elements, even if there has be no calculation.

!  Add the frozen-wf (dyfrwf) part to the ewald part (dyew),
!  the part 1 of the frozen wf part of the xc core correction
!  (dyfrx1) and the non-frozen part (dynfr) to get the second-order
!  derivative matrix (d2matr), then
!  take account of the non-cartesian coordinates (d2cart).
   ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
   ABI_ALLOCATE(carflg,(3,mpert,3,mpert))
   ABI_ALLOCATE(d2matr,(2,3,mpert,3,mpert))
   outd2=1
   call dfpt_gatherdy(becfrnl,dtset%berryopt,blkflg,carflg,&
&   dyew,dyfrwf,dyfrx1,dyfr_cplex,dyfr_nondiag,dyvdw,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
&   gprimd,dtset%mband,mpert,natom,ntypat,outd2,pawbec,pawpiezo,piezofrnl,dtset%prtbbb,&
&   rfasr,rfpert,rprimd,dtset%typat,ucvol,usevdw,psps%ziontypat)

   dscrpt=' Note : temporary (transfer) database '

!  Initialize the header of the DDB file
   call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,&
&   1,xred=xred,occ=occ,ngfft=ngfft)

!  Open the formatted derivative database file, and write the header
   call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_ddb, dtfil%unddb)

   call ddb_hdr_free(ddb_hdr)

!  Output of the dynamical matrix (master only)
   call dfpt_dyout(becfrnl,dtset%berryopt,blkflg,carflg,dtfil%unddb,ddkfil,dyew,dyfrlo,&
&   dyfrnl,dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,dyvdw,d2cart,d2cart_bbb,d2eig0,&
&   d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
&   has_full_piezo,has_allddk,ab_out,dtset%mband,mpert,natom,ntypat,&
&   outd2,pawbec,pawpiezo,piezofrnl,dtset%prtbbb,dtset%prtvol,qphon,qzero,&
&   dtset%typat,rfdir,rfpert,rfphon,rfstrs,psps%usepaw,usevdw,psps%ziontypat)

   close(dtfil%unddb)

#ifdef HAVE_NETCDF
   ! Output dynamical matrix in NetCDF format.
   call crystal_init(dtset%amu_orig(:,1), Crystal, &
&   dtset%spgroup, dtset%natom, dtset%npsp, psps%ntypat, &
&   dtset%nsym, rprimd, dtset%typat, xred, dtset%ziontypat, dtset%znucl, 1, &
&   dtset%nspden==2.and.dtset%nsppol==1, .false., hdr%title, &
&   dtset%symrel, dtset%tnons, dtset%symafm)

   filename = strcat(dtfil%filnam_ds(4),"_DDB.nc")

   call outddbnc(filename, mpert, d2matr, blkflg, dtset%qptn, Crystal)
   call crystal%free()
#endif


!  In case of phonons, diagonalize the dynamical matrix
   if(rfphon==1)then

!    First, suppress the 'wings' elements,
!    for which the diagonal element is not known
     call wings3(carflg,d2cart,mpert)

!    Check the analyticity of the dynamical matrix
     analyt=0
     if (rfpert(natom+2)==0 .or. rfpert(natom+2)==2 .or. sumg0==1 ) analyt=1

!    Diagonalize the analytic part
     ABI_ALLOCATE(displ,(2*3*natom*3*natom))
     ABI_ALLOCATE(eigval,(3*natom))
     ABI_ALLOCATE(eigvec,(2*3*natom*3*natom))
     ABI_ALLOCATE(phfrq,(3*natom))
     qphnrm=one
     call dfpt_phfrq(dtset%amu_orig(:,1),displ,d2cart,eigval,eigvec,indsym,mpert,&
&     dtset%nsym,natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,&
&     dtset%rprimd_orig(1:3,1:3,1),0,dtset%symrel,dtset%symafm,dtset%typat,ucvol)

!    Print the phonon frequencies
     call dfpt_prtph(displ,0,dtset%enunit,ab_out,natom,phfrq,qphnrm,qphon)

!    Check the completeness of the dynamical matrix and eventually send a warning
     call chkph3(carflg,0,mpert,natom)
   end if ! end case of phonons
 end if !end me == 0

!Compute the other terms for AHC dynamic and AHC full
 if (.not.(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0.or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and. rfuser==0 &
& .and. rfmagn==0)) then
   if(rfphon==1) then ! AHC can only be computed in case of phonons

!    Stuff for parallelism
     if(master /= me) then
       ABI_ALLOCATE(phfrq,(3*natom))
       ABI_ALLOCATE(displ,(2*3*natom*3*natom))
     end if
     call xmpi_bcast (phfrq,master,mpi_enreg%comm_cell,ierr) !Broadcast phfrq and displ
     call xmpi_bcast (displ,master,mpi_enreg%comm_cell,ierr) !to all processes

     if(dtset%ieig2rf == 3 .or. dtset%ieig2rf == 4 .or. dtset%ieig2rf == 5 ) then
       bdeigrf = dtset%bdeigrf
       if(dtset%bdeigrf == -1) bdeigrf = dtset%mband
!      if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
!      &         (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
!      write(std_out,*)'Reading the dense grid WF file'
!      call wfk_read_eigenvalues(dtfil%fnameabi_wfkfine,eigenq_fine,hdr_fine,mpi_enreg%comm_world)
!      ABI_CHECK(SIZE(eigenq_fine,DIM=1)==Dtset%mband,"Size eigenq_fine != mband")
!      endif
       if(dtset%kptopt==3 .or. dtset%kptopt==0 .or. dtset%nsym==1)then
         write(std_out,*) 'Entering: eig2tot'
         if(dtset%smdelta>0)then
           if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&           (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
             call eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0_pert,&
&             eigenq_pert,eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,dtset%ieig2rf,dtset%mband,mpert,&
&             dtset%natom,mpi_enreg,doccde,nkpt_rbz,dtset%nsppol,dtset%smdelta,rprimd,dtset,&
&             occ_rbz_pert,hdr0,eigbrd,eigenq_fine,hdr_fine)
           else
             call eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0_pert,&
&             eigenq_pert,eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,dtset%ieig2rf,dtset%mband,mpert,&
&             dtset%natom,mpi_enreg,doccde,nkpt_rbz,dtset%nsppol,dtset%smdelta,rprimd,dtset,&
&             occ_rbz_pert,hdr0,eigbrd)
           end if
         else
           if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&           (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
             call eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0_pert,&
&             eigenq_pert,eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,dtset%ieig2rf,dtset%mband,mpert,&
&             dtset%natom,mpi_enreg,doccde,nkpt_rbz,dtset%nsppol,dtset%smdelta,rprimd,dtset,&
&             occ_rbz_pert,hdr0,eigbrd,eigenq_fine,hdr_fine)
           else
             call eig2tot(dtfil,xred,psps,pawtab,natom,bdeigrf,clflg,dim_eig2nkq,eigen0_pert,&
&             eigenq_pert,eigen1_pert,eig2nkq,dtset%elph2_imagden,dtset%esmear,dtset%ieig2rf,dtset%mband,mpert,&
&             dtset%natom,mpi_enreg,doccde,nkpt_rbz,dtset%nsppol,dtset%smdelta,rprimd,dtset,&
&             occ_rbz_pert,hdr0)
           end if
         end if
         write(std_out,*) 'Leaving: eig2tot'
       end if
     end if
     if (dtset%ieig2rf > 0) then
       ABI_DEALLOCATE(eigen0_pert)
       ABI_DEALLOCATE(eigenq_pert)
       ABI_DEALLOCATE(occ_rbz_pert)
       ABI_DEALLOCATE(eigen1_pert)
       call hdr_free(hdr0)
       if ((dtset%getwfkfine /= 0 .and. dtset%irdwfkfine ==0) .or.&
&       (dtset%getwfkfine == 0 .and. dtset%irdwfkfine /=0) )  then
!         call hdr_free(hdr0)
         call hdr_free(hdr_fine)
         ABI_DEALLOCATE(eigenq_fine)
       end if
     end if ! ieig2rf == 3  or %ieig2rf == 4 or %ieig2rf == 5
   end if ! rfphon==1
 end if
 if(dtset%efmas>0) then
   ABI_DEALLOCATE(eigen0_pert)
   ABI_DEALLOCATE(eigen1_pert)
 end if
 ABI_DEALLOCATE(doccde)


 if(me==0)then
   if (.not.(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0 .or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and.rfuser==0 &
&   .and. rfmagn==0) )then
     if(rfphon==1)then
!      Compute and print the T=0 Fan, and possibly DDW contributions to the eigenenergies.
       if(dtset%ieig2rf > 0) then
         write(message, '(80a,9a)' ) ('=',mu=1,80),ch10,ch10,&
&         ' ---- T=0 shift of eigenenergies due to electron-phonon interation at q ---- ',ch10,&
&         ' Warning : the total shift must be computed through anaddb,                  ',ch10,&
&         ' here, only the contribution of one q point is printed.                      ',ch10,&
&         ' Print first the electronic eigenvalues, then the q-dependent Fan shift of eigenvalues.'
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,  message,'COLL')

         if(qeq0)then
           write(message, '(a)' )&
&           ' Phonons at gamma, also compute the Diagonal Debye-Waller shift of eigenvalues.'
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')
         end if

         write(message, '(a)' ) ' '
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')

         call prteigrs(eigen0,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,-1,dtset%kptns,dtset%kptopt,&
&         dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,3,0,dtset%prtvol,&
&         eigen0,zero,zero,dtset%wtk)

         write(message, '(a)' ) ch10
         call wrtout(ab_out,message,'COLL')
         call wrtout(std_out,message,'COLL')

!        Compute and print Fan contribution
         ABI_ALLOCATE(eigen_fan,(dtset%mband*dtset%nkpt*dtset%nsppol))
         ABI_ALLOCATE(eigen_fan_mean,(dtset%mband*dtset%nkpt*dtset%nsppol))
         call elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_fan,gprimd,&
&         dtset%mband,natom,dtset%nkpt,dtset%nsppol,1,phfrq,dtset%prtvol)
         call eigen_meandege(eigen0,eigen_fan,eigen_fan_mean,dtset%mband,dtset%nband,dtset%nkpt,dtset%nsppol,2)
         call prteigrs(eigen_fan_mean,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,-1,dtset%kptns,dtset%kptopt,&
&         dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,5,0,dtset%prtvol,&
&         eigen0,zero,zero,dtset%wtk)

         if(qeq0 .or. dtset%getgam_eig2nkq>0)then

           write(message, '(a)' ) ch10
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')

!          Compute and print Diagonal Debye-Waller contribution
           ABI_ALLOCATE(eigen_ddw,(dtset%mband*dtset%nkpt*dtset%nsppol))
           ABI_ALLOCATE(eigen_ddw_mean,(dtset%mband*dtset%nkpt*dtset%nsppol))
           if(qeq0)then
             call elph2_fanddw(dim_eig2nkq,displ,eig2nkq,eigen_ddw,gprimd,&
&             dtset%mband,natom,dtset%nkpt,dtset%nsppol,2,phfrq,dtset%prtvol)
             if(results_respfn%gam_jdtset == -dtset%jdtset)then
               sz1=dtset%mband*dtset%nsppol
               sz2=natom*dim_eig2nkq
               ABI_ALLOCATE(results_respfn%gam_eig2nkq,(2,sz1,dtset%nkpt,3,natom,3,sz2))
               results_respfn%gam_eig2nkq(:,:,:,:,:,:,:)=eig2nkq(:,:,:,:,:,:,:)
               results_respfn%gam_jdtset=dtset%jdtset
             end if
           else if(dtset%getgam_eig2nkq>0)then
             if(results_respfn%gam_jdtset==dtset%getgam_eig2nkq)then
               call elph2_fanddw(dim_eig2nkq,displ,results_respfn%gam_eig2nkq,eigen_ddw,&
&               gprimd,dtset%mband,natom,dtset%nkpt,dtset%nsppol,2,phfrq,dtset%prtvol)
             else
               write(message,'(a,i0,2a,i0,2a)')&
&               'results_respfn%gam_jdtset=',results_respfn%gam_jdtset,ch10,&
&               'dtset%getgam_eig2nkq=',dtset%getgam_eig2nkq,ch10,&
&               'So, it seems eig2nkq at gamma has not yet been computed, while it is needed now.'
               MSG_BUG(message)
             end if
           end if
           call eigen_meandege(eigen0,eigen_ddw,eigen_ddw_mean,dtset%mband,dtset%nband,dtset%nkpt,dtset%nsppol,2)
           call prteigrs(eigen_ddw_mean,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,-1,dtset%kptns,dtset%kptopt,&
&           dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,6,0,dtset%prtvol,&
&           eigen0,zero,zero,dtset%wtk)

           write(message, '(a)' ) ch10
           call wrtout(ab_out,message,'COLL')
           call wrtout(std_out,message,'COLL')

!          Print sum of mean Fan and DDW
           ABI_ALLOCATE(eigen_fanddw,(dtset%mband*dtset%nkpt*dtset%nsppol))
           eigen_fanddw=eigen_fan_mean+eigen_ddw_mean
           call prteigrs(eigen_fanddw,dtset%enunit,fermie,dtfil%fnameabo_eig,ab_out,-1,dtset%kptns,dtset%kptopt,&
&           dtset%mband,dtset%nband,dtset%nkpt,1,dtset%nsppol,occ,dtset%occopt,7,0,dtset%prtvol,&
&           eigen0,zero,zero,dtset%wtk)

           ABI_DEALLOCATE(eigen_ddw)
           ABI_DEALLOCATE(eigen_ddw_mean)
           ABI_DEALLOCATE(eigen_fanddw)

         end if

         ABI_DEALLOCATE(eigen_fan)
         ABI_DEALLOCATE(eigen_fan_mean)
       end if

!      In case of a non-analytical part,
!      get the phonon frequencies for three different directions (in cartesian coordinates)
       if(analyt==0)then
         qphnrm=zero
         do idir=1,3
!          Need to know the corresponding dielectric constant
           if(carflg(idir,natom+2,idir,natom+2)==1)then
             qphon(:)=zero ; qphon(idir)=one
!            Get the phonon frequencies
             call dfpt_phfrq(dtset%amu_orig(:,1),displ,d2cart,eigval,eigvec,indsym,mpert,&
&             dtset%nsym,natom,dtset%nsym,ntypat,phfrq,qphnrm,qphon,&
&             dtset%rprimd_orig(1:3,1:3,1),0,dtset%symrel,dtset%symafm,dtset%typat,ucvol)
!            Print the phonon frequencies
             call dfpt_prtph(displ,0,dtset%enunit,ab_out,natom,phfrq,qphnrm,qphon)
!            Check the completeness of the dynamical matrix
!            and eventually send a warning
             call chkph3(carflg,idir,mpert,natom)
           end if
         end do
         if (idir < 4) then
           qphon(idir)=zero
         end if
       end if

       ABI_DEALLOCATE(displ)
       ABI_DEALLOCATE(eigval)
       ABI_DEALLOCATE(eigvec)
       ABI_DEALLOCATE(phfrq)
     end if ! rfphon == 1
     ABI_DEALLOCATE(carflg)
     ABI_DEALLOCATE(d2cart)
     ABI_DEALLOCATE(d2matr)
   end if ! End condition on if.not.
 end if ! master node

!Deallocate arrays
 if (allocated(displ)) then
   ABI_DEALLOCATE(displ)
 end if
 if (allocated(eigval)) then
   ABI_DEALLOCATE(eigval)
 end if
 if (allocated(eigvec)) then
   ABI_DEALLOCATE(eigvec)
 end if
 if (allocated(phfrq)) then
   ABI_DEALLOCATE(phfrq)
 end if

 ABI_DEALLOCATE(clflg)
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(blkflg)
 ABI_DEALLOCATE(dyew)
 ABI_DEALLOCATE(dyewq0)
 ABI_DEALLOCATE(dyfrlo)
 ABI_DEALLOCATE(dyfrnl)
 ABI_DEALLOCATE(dyfrwf)
 ABI_DEALLOCATE(dyfrx1)
 ABI_DEALLOCATE(dyfrx2)
 ABI_DEALLOCATE(dyvdw)
 ABI_DEALLOCATE(d2bbb)
 ABI_DEALLOCATE(d2cart_bbb)
 ABI_DEALLOCATE(d2eig0)
 ABI_DEALLOCATE(d2k0)
 ABI_DEALLOCATE(d2lo)
 ABI_DEALLOCATE(d2loc0)
 ABI_DEALLOCATE(d2nfr)
 ABI_DEALLOCATE(d2nl)
 ABI_DEALLOCATE(d2nl0)
 ABI_DEALLOCATE(d2nl1)
 ABI_DEALLOCATE(d2ovl)
 ABI_DEALLOCATE(d2vn)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(eig2nkq)
 ABI_DEALLOCATE(eigbrd)
 ABI_DEALLOCATE(eltcore)
 ABI_DEALLOCATE(elteew)
 ABI_DEALLOCATE(eltfrhar)
 ABI_DEALLOCATE(eltfrnl)
 ABI_DEALLOCATE(eltfrloc)
 ABI_DEALLOCATE(eltfrkin)
 ABI_DEALLOCATE(eltfrxc)
 ABI_DEALLOCATE(eltvdw)
 ABI_DEALLOCATE(becfrnl)
 ABI_DEALLOCATE(piezofrnl)
 call efmasdeg_free_array(efmasdeg)
 call efmasval_free_array(efmasval)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(pertsy)
 ABI_DEALLOCATE(rfpert)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(ylm)
 ABI_DEALLOCATE(ylmgr)
 call pawfgr_destroy(pawfgr)
 if (psps%usepaw==1) then
   call pawrhoij_free(pawrhoij)
   call paw_an_free(paw_an)
   call paw_ij_free(paw_ij)
   call pawfgrtab_free(pawfgrtab)
 end if
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(nhatgr)
 ABI_DATATYPE_DEALLOCATE(pawrhoij)
 ABI_DATATYPE_DEALLOCATE(paw_an)
 ABI_DATATYPE_DEALLOCATE(paw_ij)
 ABI_DATATYPE_DEALLOCATE(pawfgrtab)
 if(rfphon==1.and.psps%n1xccc/=0)then
   ABI_DEALLOCATE(blkflgfrx1)
 end if

!Clean the header
 call hdr_free(hdr)

!Clean GPU data
#if defined HAVE_GPU_CUDA
 if (dtset%use_gpu_cuda==1) then
   call dealloc_hamilt_gpu(0,dtset%use_gpu_cuda)
 end if
#endif

 call timab(138,2,tsec)
 call timab(132,2,tsec)

 DBG_EXIT("COLL")

end subroutine respfn
!!***

!!****f* m_respfn_driver/wrtloctens
!! NAME
!! wrtloctens
!!
!! FUNCTION
!! Output of the localisation tensor
!!
!! INPUTS
!! blkflg = flags for each element of the 2DTE (=1 if computed)
!! d2bbb  = band by band decomposition of second order derivatives
!! d2nl   = non-local contributions to the 2DTEs
!! mband  = maximum number of bands
!! mpert  = maximum number of ipert
!! natom  = number of atoms in unit cell
!! prtbbb = if = 1, write the band by band decomposition of the localization tensor
!! rprimd = dimensional primitive translations for real space (bohr)
!! usepaw = flag for PAW
!!
!! OUTPUT
!!  (only writing)
!!
!! TODO
!!  The localization tensor cannot be defined in the metallic case. It should not be computed.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      wrtout
!!
!! SOURCE

subroutine wrtloctens(blkflg,d2bbb,d2nl,mband,mpert,natom,prtbbb,rprimd,usepaw)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mpert,natom,prtbbb,usepaw
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)

!Local variables ------------------------------
!scalars
 integer :: flag,iband,idir,idir2,jband
 character(len=500) :: message
!arrays
 real(dp) :: loctenscart(2,3,3)
 real(dp),allocatable :: loctenscart_bbb(:,:,:,:,:)

! *********************************************************************

!This feature is disabled in the PAW case
 if (usepaw==1) return

 if(prtbbb==1)then
   ABI_ALLOCATE(loctenscart_bbb,(2,3,3,mband,mband*prtbbb))
   loctenscart_bbb(:,:,:,:,:)=zero
 end if

!complete missing elements
 flag = 0
 do idir2 = 1,3
   do idir = 1,3

     if (blkflg(idir2,natom+1,idir,natom+1)==0) then
       if (blkflg(idir,natom+1,idir2,natom+1)==0) then
         flag = 1
       else
         d2nl(1,idir2,natom+1,idir,natom+1) = d2nl(1,idir,natom+1,idir2,natom+1)
         d2nl(2,idir2,natom+1,idir,natom+1) =-d2nl(2,idir,natom+1,idir2,natom+1)
         if(prtbbb==1)then
           d2bbb(1,idir2,idir,natom+1,:,:) = d2bbb(1,idir,idir2,natom+1,:,:)
           d2bbb(2,idir2,idir,natom+1,:,:) =-d2bbb(2,idir,idir2,natom+1,:,:)
         end if

       end if
     end if

   end do ! idir=1,3
 end do ! idir2=1,3

!Transform the tensor to cartesian coordinates

 loctenscart(1,:,:) = matmul(rprimd,d2nl(1,:,natom+1,:,natom+1))
 loctenscart(2,:,:) = matmul(rprimd,d2nl(2,:,natom+1,:,natom+1))

 loctenscart(1,:,:) = matmul(loctenscart(1,:,:),transpose(rprimd))
 loctenscart(2,:,:) = matmul(loctenscart(2,:,:),transpose(rprimd))

 loctenscart(:,:,:) = loctenscart(:,:,:)/(two_pi**2)

 if (prtbbb == 1) then

   write(message,'(a,a)')ch10, &
&   ' Band by band decomposition of the localisation tensor (bohr^2)'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

   do iband = 1,mband
     do jband = 1,mband

       loctenscart_bbb(1,:,:,iband,jband) = matmul(rprimd,d2bbb(1,:,:,natom+1,iband,jband))
       loctenscart_bbb(2,:,:,iband,jband) = matmul(rprimd,d2bbb(2,:,:,natom+1,iband,jband))

       loctenscart_bbb(1,:,:,iband,jband) = matmul(loctenscart_bbb(1,:,:,iband,jband),transpose(rprimd))
       loctenscart_bbb(2,:,:,iband,jband) = matmul(loctenscart_bbb(2,:,:,iband,jband),transpose(rprimd))

       loctenscart_bbb(:,:,:,iband,jband) = loctenscart_bbb(:,:,:,iband,jband)/(two_pi**2)


       write(message,'(a,a,i5,a,i5,a)')ch10, &
&       ' Localisation tensor (bohr^2) for band ',iband,',',jband, &
&       ' in cartesian coordinates'
       call wrtout(std_out,message,'COLL')
       call wrtout(ab_out,message,'COLL')

       write(ab_out,*)'     direction              matrix element'
       write(ab_out,*)'  alpha     beta       real part   imaginary part'
       do idir2 = 1,3
         do idir = 1,3
           write(ab_out,'(5x,i1,8x,i1,3x,2f16.10)')idir2,idir,&
&           loctenscart_bbb(1,idir2,idir,iband,jband),&
&           loctenscart_bbb(2,idir2,idir,iband,jband)
         end do
       end do

     end do !jband
   end do !iband

 end if  !prtbbb

 if (usepaw==0) then
   write(message,'(a,a,a,a)')ch10, &
&   ' Total localisation tensor (bohr^2) in cartesian coordinates',ch10,&
&   '  WARNING : still subject to testing - especially symmetries.'
 else
   write(message,'(a,a,a,a,a,a)')ch10, &
&   ' Total localisation tensor (bohr^2) in cartesian coordinates',ch10,&
&   '  WARNING : probably wrong for PAW (printing for testing purpose)',ch10,&
&   '  WARNING : still subject to testing - especially symmetries.'
 end if
 call wrtout(std_out,message,'COLL')
 call wrtout(ab_out,message,'COLL')

 write(ab_out,*)'     direction              matrix element'
 write(ab_out,*)'  alpha     beta       real part   imaginary part'
 do idir2 = 1,3
   do idir = 1,3
     write(ab_out,'(5x,i1,8x,i1,3x,2f16.10)')idir2,idir,&
&     loctenscart(1,idir2,idir),&
&     loctenscart(2,idir2,idir)
   end do
 end do

 if (flag == 1) then
   write(message,'(6a)')ch10,&
&   ' WARNING : Localization tensor calculation (this does not apply to other properties).',ch10,&
&   '  Not all d/dk perturbations were computed. So the localization tensor in reciprocal space is incomplete,',ch10,&
&   '  and transformation to cartesian coordinates may be wrong.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')
 end if

 if (prtbbb == 1) then
   ABI_DEALLOCATE(loctenscart_bbb)
 end if

end subroutine wrtloctens
!!***

!!****f* ABINIT/dfpt_dyout
!! NAME
!! dfpt_dyout
!!
!!
!! FUNCTION
!! Output of all quantities related to the 2nd-order matrix:
!! Ewald part, local and non-local frozen wf part,
!! core contributions,
!! local and non-local variational part, 2nd-order
!! matrix itself, and, for the phonon part,
!! eigenfrequencies, in Hartree, meV and cm-1.
!! Also output unformatted 2nd-order matrix for later
!! use in the Brillouin-zone interpolation
!!
!! INPUTS
!!  becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!!  blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!  matrix has been calculated ; 0 otherwise )
!!  carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!!  ddboun=unit number for the derivative database output
!!  ddkfil(3)=components are 1 if corresponding d/dk file exists, otherwise 0
!!  (in what follows, DYMX means dynamical matrix, and D2MX means 2nd-order matrix)
!!  dyew(2,3,natom,3,natom)=Ewald part of the DYMX
!!  dyfrlo(3,3,natom)=frozen wf local part of the DYMX
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=frozen wf nonloc part of the DYMX
!!  dyfrx1(2,3,natom,3,natom)=frozen wf nonlin. xc core corr.(1)
!!    part of the DYMX
!!  dyfrx2(3,3,natom)=frozen wf nonlin. xc core corr.(2) part of the DYMX
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl is non diagonal with respect to atoms; 0 otherwise
!!  dyvdw(2,3,natom,3,natom*usevdw)=vdw DFT-D part of the dynamical matrix
!!  d2cart(2,3,mpert,3,mpert)=D2MX in cartesian coordinates
!!  d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)=
!!    band by band decomposition of Born effective charges
!!    (calculated from phonon-type perturbation) in cartesian coordinates
!!  d2eig0(2,3,mpert,3,mpert)=0-order eigen. station. part of the D2MX
!!  d2k0(2,3,mpert,3,mpert)=0-order kinet. station. part of the D2MX
!!  d2lo(2,3,mpert,3,mpert)=nonstation. local part of the D2MX
!!  d2loc0(2,3,mpert,3,mpert)=0-order loc station. part of the D2MX
!!  d2matr(2,3,mpert,3,mpert)=D2MX in non-cartesian coordinates
!!  d2nl(2,3,mpert,3,mpert)=nonstation. nonloc part of the D2MX
!!  d2nl0(2,3,mpert,3,mpert)=0-order nonloc station. part of the D2MX
!!  d2nl1(2,3,mpert,3,mpert)=1-order nonloc station. part of the D2MX
!!  d2ovl(2,mpert,3,mpert*usepaw)=1st-order change of WF overlap contributions to the 2DTEs (PAW)
!!  d2vn(2,3,mpert,3,mpert)=potential*dens station. part of the D2MX and without masses included)
!!  eltcore(6,6)=core contribution to the elastic tensor
!!  elteew(6+3*natom,6)=Ewald contribution to the elastic tsenor
!!  eltfrhar(6,6)=hartree contribution to the elastic tensor
!!  eltfrkin(6,6)=kinetic contribution to the elastic tensor
!!  eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!!  eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!!  eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!!  eltvdw(6+3*natom,6*usevdw)=vdw DFT-D part of the elastic tensor
!!  has_full_piezo=the full calculation of the piezoelectric tensor from electric field perturbation
!!                 is only available if nsym=1 (strain perturbation is not symmetrized)
!!  has_allddk= True if all ddk file are present on disk, so the effective charge or piezzo
!!              electric tensor are correctly computed (PAW ONLY)
!!  iout=unit number for long write-up
!!  mband=maximum number of bands
!!  mpert =maximum number of ipert
!!  natom=number of atoms
!!  ntypat=number of atom types
!!  outd2=option for the output of the 2nd-order matrix :
!!   if outd2=1, non-stationary part
!!   if outd2=2, stationary part.
!!  pawbec= flag for the computation of frozen part of Born Effective Charges (PAW only)
!!  pawpiezo= flag for the computation of frozen part of Piezoelectric tensor (PAW only)
!!  prtbbb=if 1, print the band-by-band decomposition
!!  prtvol=print volume
!!  qphon(3)=phonon wavelength, in reduced coordinates
!!  qzero=1 if zero phonon wavevector
!!  rfdir(3)=defines the directions for the perturbations
!!  rfpert(mpert)=defines the perturbations
!!  rfphon=if 1, there are phonon perturbations
!!  rfstrs=if 1,2,3 there are strain perturbations
!!  typat(natom)=integer label of each type of atom (1,2,...)
!!  usepaw=1 if PAW, 0 otherwise
!!  usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!!  zion(ntypat)=charge corresponding to the atom type
!!
!! SIDE EFFECTS
!!  d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
!!
!! NOTES
!! This routine is called only by the processor me==0 .
!! In consequence, no use of message and wrtout routine.
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_dyout(becfrnl,berryopt,blkflg,carflg,ddboun,ddkfil,dyew,dyfrlo,dyfrnl,&
& dyfrx1,dyfrx2,dyfr_cplex,dyfr_nondiag,dyvdw,d2cart,d2cart_bbb,&
& d2eig0,d2k0,d2lo,d2loc0,d2matr,d2nl,d2nl0,d2nl1,d2ovl,d2vn,&
& eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
& has_full_piezo,has_allddk,iout,mband,mpert,natom,ntypat,&
& outd2,pawbec,pawpiezo,piezofrnl,prtbbb,prtvol,qphon,qzero,typat,rfdir,&
& rfpert,rfphon,rfstrs,usepaw,usevdw,zion)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,ddboun,dyfr_cplex,dyfr_nondiag,iout,mband,mpert
 integer,intent(in) :: natom,ntypat,outd2,pawbec,pawpiezo,prtbbb,prtvol,qzero
 integer, intent(in) :: rfphon,rfstrs,usepaw,usevdw
!arrays
 integer,intent(in) :: blkflg(3,mpert,3,mpert),carflg(3,mpert,3,mpert)
 integer,intent(in) :: ddkfil(3),rfdir(3),rfpert(mpert),typat(natom)
 real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(in) :: d2cart(2,3,mpert,3,mpert),d2eig0(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2k0(2,3,mpert,3,mpert),d2lo(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2loc0(2,3,mpert,3,mpert),d2matr(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2nl(2,3,mpert,3,mpert),d2nl0(2,3,mpert,3,mpert)
 real(dp),intent(in) :: d2nl1(2,3,mpert,3,mpert),d2ovl(2,3,mpert,3,mpert*usepaw)
 real(dp),intent(in) :: d2vn(2,3,mpert,3,mpert)
 real(dp),intent(in) :: dyew(2,3,natom,3,natom),dyfrlo(3,3,natom)
 real(dp),intent(in) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom),dyfrx2(3,3,natom)
 real(dp),intent(in) :: dyvdw(2,3,natom,3,natom*usevdw)
 real(dp),intent(in) :: eltcore(6,6),elteew(6+3*natom,6)
 real(dp),intent(in) :: eltfrhar(6,6),eltfrkin(6,6),eltfrloc(6+3*natom,6)
 real(dp),intent(in) :: eltfrnl(6+3*natom,6),eltfrxc(6+3*natom,6)
 real(dp),intent(in) :: eltvdw(6+3*natom,6*usevdw),piezofrnl(6,3*pawpiezo)
 real(dp),intent(in) :: qphon(3),zion(ntypat)
 real(dp),intent(inout) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
 logical,intent(in) :: has_allddk,has_full_piezo

!Local variables -------------------------
!scalars
 integer :: iband,idir1,idir2,ii,ipert1,ipert2,jj,nelmts,nline
 real(dp) :: qptnrm,zi,zr
!arrays
 real(dp) :: delta(3,3)

! *********************************************************************

!DEBUG
!write(std_out,*)' dfpt_dyout : enter '
!write(std_out,*)ddkfil
!ENDDEBUG

!Long print : includes detail of every part of the 2nd-order energy
 if(prtvol>=10)then

!  In case of phonon
   if (rfphon==1)then

!    write the Ewald part of the dynamical matrix
     write(iout,*)' '
     write(iout,*)' Ewald part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 dyew(1,idir1,ipert1,idir2,ipert2),&
&                 dyew(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the local frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyfrlo(idir1,idir2,ipert2),zero
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlo frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf non-local part of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2.or.dyfr_nondiag==1)then
                   if (dyfr_cplex==1) then
                     write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                     dyfrnl(1,idir1,idir2,ipert1,1+(ipert2-1)*dyfr_nondiag),zero
                   else
                     write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                     dyfrnl(:,idir1,idir2,ipert1,1+(ipert2-1)*dyfr_nondiag)
                   end if
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlinear xc core correction(1) part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc core (1) part',&
&     ' of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 dyfrx1(1,idir1,ipert1,idir2,ipert2),&
&                 dyfrx1(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlinear xc core correction(2) part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc core (2) part',&
&     ' of the dynamical matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           write(iout,*)' '
           do ipert2=1,natom
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 if(ipert1==ipert2)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyfrx2(idir1,idir2,ipert2),zero
                 else
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   zero,zero
                 end if
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the DFT-D vdw part of the dynamical matrix
     if (usevdw==1) then
       write(iout,*)' '
       write(iout,*)' DFT-D van der Waals part of the dynamical matrix'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=1,natom
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             write(iout,*)' '
             do ipert2=1,natom
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   dyvdw(1,idir1,ipert1,idir2,ipert2),&
&                   dyvdw(2,idir1,ipert1,idir2,ipert2)
                 end if
               end do
             end do
           end if
         end do
       end do
     end if

!    End of the phonon condition
   end if

!  In case of atom. strain/electric field perturbation (piezoelectric tensor)
   if (pawpiezo==1.and.(rfpert(natom+2)==1.or.rfstrs/=0).and.outd2==1)then
     write(iout,*)' '
     write(iout,*)' Frozen wf part of the piezoelectric tensor'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     ipert1=natom+2
     do idir1=1,3
       write(iout,*)' '
       ii=1
       do ipert2=natom+3,natom+4
         do idir2=1,3
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           piezofrnl(ii,idir1),zero
           ii=ii+1
         end do
       end do
     end do
   end if

!  In case of atom. displ/electric field perturbation (Born Effective Charges)
   if (pawbec==1.and.(rfpert(natom+2)==1.or.rfphon==1).and.outd2==1)then
     write(iout,*)' '
     write(iout,*)' Frozen wf part of the Born Effective Charges'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     ipert1 = natom+2
     do ipert2=1,natom
       do idir1=1,3
         write(iout,*)' '
         do idir2=1,3
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           becfrnl(idir2,ipert2,idir1),zero
         end do
       end do
     end do
   end if

!  In case of strain
   if (rfstrs/=0)then

!    Write the Ewald part of the elastic tensor
     write(iout,*)' '
     write(iout,*)' Ewald part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 elteew(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Write the Ewald part of the internal strain coupling parameters
     write(iout,*)' '
     write(iout,*)' Ewald part of the internal strain coupling parameters'
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 elteew(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the local frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrloc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf local part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrloc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the nonlo frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf nonlocal part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrnl(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf nonlocal part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrnl(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the xc part
     write(iout,*)' '
     write(iout,*)' Frozen wf xc part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrxc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

     write(iout,*)' '
     write(iout,*)' Frozen wf xc part of the internal strain coupling parameters '
     write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,natom
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+6+3*(ipert1-1)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrxc(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the kinetic frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf kinetic part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrkin(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the hartree frozen wf part
     write(iout,*)' '
     write(iout,*)' Frozen wf hartree part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltfrhar(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the psp core part
     write(iout,*)' '
     write(iout,*)' Psp core part of the elastic tensor in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=natom+3,natom+4
       do idir1=1,3
         if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&         .or.   outd2==1                           )then
           ii=idir1+3*(ipert1-natom-3)
           write(iout,*)' '
           do ipert2=natom+3,natom+4
             do idir2=1,3
               if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 jj=idir2+3*(ipert2-natom-3)
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 eltcore(ii,jj),zero
               end if
             end do
           end do
         end if
       end do
     end do

!    Now the DFT-D vdw part
     if (usevdw==1) then
       write(iout,*)' '
       write(iout,*)' DFT-D van der Waals part of the elastic tensor in cartesian coordinates'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=natom+3,natom+4
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             ii=idir1+3*(ipert1-natom-3)
             write(iout,*)' '
             do ipert2=natom+3,natom+4
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   jj=idir2+3*(ipert2-natom-3)
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   eltvdw(ii,jj),zero
                 end if
               end do
             end do
           end if
         end do
       end do

       write(iout,*)' '
       write(iout,*)' DFT-D2 van der Waals part of the internal strain coupling parameters'
       write(iout,*)'  (cartesian strain, reduced atomic coordinates)'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part   imaginary part'
       do ipert1=1,natom
         do idir1=1,3
           if ( (rfpert(ipert1)==1.and.rfdir(idir1)==1)&
&           .or.   outd2==1                           )then
             ii=idir1+6+3*(ipert1-1)
             write(iout,*)' '
             do ipert2=natom+3,natom+4
               do idir2=1,3
                 if (rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                   jj=idir2+3*(ipert2-natom-3)
                   write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                   eltvdw(ii,jj),zero
                 end if
               end do
             end do
           end if
         end do
       end do
     end if ! usevdw

!    End of the strain condition
   end if

!  Now the local nonstationary nonfrozenwf part
   if (outd2==1)then
     write(iout,*)' '
     write(iout,*)' Non-stationary local part of the 2-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2lo(1,idir1,ipert1,idir2,ipert2),&
&                 d2lo(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the nonlocal nonstationary nonfrozenwf part
   if (outd2==1)then
     write(iout,*)' '
     write(iout,*)' Non-stationary non-local part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the overlap change nonstationnary nonfrozenwf part (PAW only)
   if (outd2==1.and.usepaw==1)then
     write(iout,*)' '
     write(iout,*)' PAW: Non-stationary WF-overlap part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if ((ipert1<=natom .or.&
&         (ipert1==natom+2.and.qzero==1.and.ddkfil(idir1)/=0)).or.&
&         ((ipert1==natom+3.or.ipert1==natom+4).and.&
&         (rfpert(natom+3)==1.or.rfpert(natom+4)==1)))then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2ovl(1,idir1,ipert1,idir2,ipert2),&
&                 d2ovl(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the 0-order local stationary nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order local part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2loc0(1,idir1,ipert1,idir2,ipert2),&
&                 d2loc0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order kinetic nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order kinetic part of the 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2k0(1,idir1,ipert1,idir2,ipert2),&
&                 d2k0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order eigenvalue nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order eigenvalue part of the'&
&     ,' 2nd-order matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2eig0(1,idir1,ipert1,idir2,ipert2),&
&                 d2eig0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary potential-density nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Station. potential-density part of the ',&
&     ' 2nd-order matrix'
     write(iout,*)'  (Note : include some xc core-correction) '
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2vn(1,idir1,ipert1,idir2,ipert2),&
&                 d2vn(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 0-order nonloc nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 0-order nonlocal part of the 2-order'&
&     ,' matrix'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl0(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl0(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  Now the stationary 1-order nonloc nonfrozenwf part
   if (outd2==2)then
     write(iout,*)' '
     write(iout,*)' Stationary 1-order nonlocal part of the'&
&     ,' 2nd-order matrix'
     write(iout,*)' (or the ddk wf part of it, in case of',&
&     ' an electric field perturbation )'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part   imaginary part'
     do ipert1=1,mpert
       do idir1=1,3
         if(rfpert(ipert1)==1.and.rfdir(idir1)==1)then
           write(iout,*)' '
           do ipert2=1,mpert
             do idir2=1,3
               if(rfpert(ipert2)==1.and.rfdir(idir2)==1)then
                 write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&                 d2nl1(1,idir1,ipert1,idir2,ipert2),&
&                 d2nl1(2,idir1,ipert1,idir2,ipert2)
               end if
             end do
           end do
         end if
       end do
     end do
   end if

!  End of the long print out condition
 end if


!Derivative database initialisation

!Calculation of the number of elements
 nelmts=0
 do ipert1=1,mpert
   do idir1=1,3
     do ipert2=1,mpert
       do idir2=1,3
         nelmts=nelmts+blkflg(idir1,ipert1,idir2,ipert2)
       end do
     end do
   end do
 end do

 if(outd2==2)then
   write(ddboun, '(/,a,i8)' ) ' 2nd derivatives (stationary) - # elements :',nelmts
 else if(outd2==1)then
   write(ddboun, '(/,a,i8)' ) ' 2nd derivatives (non-stat.)  - # elements :',nelmts
 end if

!Phonon wavevector
 qptnrm=1.0_dp

!Note : if qptnrm should assume another value, it should
!be checked if the f6.1 format is OK.
 write(ddboun, '(a,3es16.8,f6.1)' ) ' qpt',(qphon(ii),ii=1,3),qptnrm

!Now the whole 2nd-order matrix, but not in cartesian coordinates,
!and masses not included
 write(iout,*)' '
 write(iout,*)' 2nd-order matrix (non-cartesian coordinates,',' masses not included,'
 write(iout,*)'  asr not included )'
 if(rfstrs/=0) then
   write(iout,*)' cartesian coordinates for strain terms (1/ucvol factor '
   write(iout,*)'  for elastic tensor components not included) '
 end if
 write(iout,*)'    j1       j2             matrix element'
 write(iout,*)' dir pert dir pert     real part     imaginary part'
 nline=1
 do ipert1=1,mpert
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=1,mpert
       do idir2=1,3
         if(blkflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2matr(1,idir1,ipert1,idir2,ipert2),&
&           d2matr(2,idir1,ipert1,idir2,ipert2)
           write(ddboun, '(4i4,2d22.14)' )idir1,ipert1,idir2,ipert2,&
&           d2matr(1,idir1,ipert1,idir2,ipert2),&
&           d2matr(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end do

!Now the dynamical matrix
 if(rfphon==1)then
   write(iout,*)' '
   write(iout,*)' Dynamical matrix, in cartesian coordinates,'
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=1,natom
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=1,natom
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the dielectric tensor ! normal case
 if(rfpert(natom+2)==1)then

   write(iout,*)' '
   write(iout,*)' Dielectric tensor, in cartesian coordinates,'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert1=natom+2
   ipert2=natom+2
   nline=1
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do idir2=1,3
       if(carflg(idir1,ipert1,idir2,ipert2)==1)then
         nline=nline+1
         write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&         d2cart(1,idir1,ipert1,idir2,ipert2),&
&         d2cart(2,idir1,ipert1,idir2,ipert2)
       end if
     end do
   end do

   if (prtbbb == 1) then

     delta(:,:) = zero
     delta(1,1) = one ; delta(2,2) = one ; delta(3,3) = one

     write(iout,*)
     write(iout,*)'Band by band decomposition of the dielectric tensor'
     write(iout,*)' '

     write(iout,*)' Vacuum polarization'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     nline=1
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do idir2=1,3
         nline=nline+1
         write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&         delta(idir2,idir1),zero
       end do
     end do

     do iband = 1,mband
       write(iout,*)' '
       write(iout,*)' Dielectric tensor, in cartesian coordinates, for band',iband
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part    imaginary part'
       ipert1 = natom + 2
       ipert2 = natom + 2
       nline=1
       do idir1=1,3
         if(nline/=0)write(iout,*)' '
         nline=0
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
!            substract vacuum polarization
             if (idir1 == idir2) then
               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband) = &
&               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband) - 1
             end if
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart_bbb(1,idir1,idir2,ipert2,iband,iband),&
&             d2cart_bbb(2,idir1,idir2,ipert2,iband,iband)
           end if
         end do
       end do
     end do !iband

   end if !prtbbb

 end if ! end natom+2 dielectric output

!Now the effective charges
!In case of the stationary calculation
 if(outd2==2 .and. rfpert(natom+2)==1 .and.rfphon==1)then
   write(iout,*)' '
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert1=natom+2
   nline=1
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=1,natom
       do idir2=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

!Now in case of the non-stationary calculation
 if(outd2==1 .and. rfpert(natom+2)==1)then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Born effectives charges are not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)' (from electric field response) '
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   ipert2=natom+2
   nline=1
   do idir2=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert1=1,natom
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

 if(outd2==1 .and. rfphon==1 .and. qzero==1&
& .and. ( (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0) .or.   &
& berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. berryopt==14 .or. berryopt==16 .or. berryopt==17 ) )then  !!HONG  need to test for fixed E and D
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Born effectives charges are not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Effective charges, in cartesian coordinates,'
   write(iout,*)' (from phonon response) '
   write(iout,*)'  if specified in the inputs, asr has been imposed'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert2=1,natom
     do idir2=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       ipert1=natom+2
       do idir1=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
   write(iout,*)' '
   write(iout,*)' '
   write(iout,*)' '

   if (prtbbb == 1) then

     write(iout,*)'Band by band decomposition of the Born effective charges'
     write(iout,*)' '
     write(iout,*)'Ionic charges in cartesian coordinates'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     zr = zero
     zi = zero
     do ipert2=1,natom
       do idir2=1,3
         if(nline/=0)write(iout,*)' '
         nline=0
         ipert1=natom+2
         do idir1=1,3
           zr = zero
           if (idir1 == idir2) then
             zr = zion(typat(ipert2))
           end if
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           zr,zi
         end do
       end do
     end do

     do iband = 1,mband
       write(iout,*)' '
       write(iout,*)' Effective charges, in cartesian coordinates, for band',iband
       write(iout,*)' (from phonon response) '
       write(iout,*)'  if specified in the inputs, asr has been imposed'
       write(iout,*)'    j1       j2             matrix element'
       write(iout,*)' dir pert dir pert     real part    imaginary part'
       nline=1
       do ipert2=1,natom
         do idir2=1,3
           if(nline/=0)write(iout,*)' '
           nline=0
           ipert1=natom+2
           do idir1=1,3
             if(carflg(idir1,ipert1,idir2,ipert2)==1)then
               nline=nline+1
               write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&               d2cart_bbb(1,idir1,idir2,ipert2,iband,iband),&
&               d2cart_bbb(2,idir1,idir2,ipert2,iband,iband)
             end if
           end do
         end do
       end do
     end do !iband
   end if !prtbbb
 end if ! end of print effective charges

!Now the elastic tensor
 if(rfstrs/=0) then
   write(iout,*)' '
   write(iout,*)' Rigid-atom elastic tensor , in cartesian coordinates,'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=natom+3,natom+4
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the internal strain coupling parameters
 if(rfstrs/=0) then
   write(iout,*)' '
   write(iout,*)' Internal strain coupling parameters, in cartesian coordinates,'
   write(iout,*)'  zero average net force deriv. has been imposed  '
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   do ipert1=1,natom
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&             d2cart(1,idir1,ipert1,idir2,ipert2),&
&             d2cart(2,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
     end do
   end do
 end if

!Now the piezoelectric tensor
 if(rfstrs/=0 .and. (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0))then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Rigid-atom proper piezoelectric tensor is not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   write(iout,*)' Rigid-atom proper piezoelectric tensor, in cartesian coordinates,'
   write(iout,*)' (from strain response)'
   write(iout,*)'    j1       j2             matrix element'
   write(iout,*)' dir pert dir pert     real part    imaginary part'
   nline=1
   ipert1=natom+2
   do idir1=1,3
     if(nline/=0)write(iout,*)' '
     nline=0
     do ipert2=natom+3,natom+4
       do idir2=1,3
         if(carflg(idir1,ipert1,idir2,ipert2)==1)then
           nline=nline+1
           write(iout,'(2(i4,i5),2(1x,f20.10))')idir1,ipert1,idir2,ipert2,&
&           d2cart(1,idir1,ipert1,idir2,ipert2),&
&           d2cart(2,idir1,ipert1,idir2,ipert2)
         end if
       end do
     end do
   end do
 end if

!Now the piezoelectric tensor
 if(outd2==1 .and. (pawpiezo==1.and.rfpert(natom+2)==1)&
& .and. (ddkfil(1)/=0.or.ddkfil(2)/=0.or.ddkfil(3)/=0)) then
   write(iout,*)' '
   if(usepaw==1.and..not.(has_allddk))then
     write(iout,*)' Warning: Rigid-atom proper piezoelectric tensor is not correctly computed'
     write(iout,*)' you need all ddk perturbations!'
   end if
   if(usepaw==1.and..not.has_full_piezo)then
     write(iout,*)' Warning: The rigid-atom proper piezoelectric tensor'
     write(iout,*)' from  electric field response requires nsym=1'
   end if
   if (has_full_piezo) then
     write(iout,*)' Rigid-atom proper piezoelectric tensor, in cartesian coordinates,'
     write(iout,*)' (from electric field response)'
     write(iout,*)'    j1       j2             matrix element'
     write(iout,*)' dir pert dir pert     real part    imaginary part'
     nline=1
     ipert1=natom+2
     do idir1=1,3
       if(nline/=0)write(iout,*)' '
       nline=0
       do ipert2=natom+3,natom+4
         do idir2=1,3
           if(carflg(idir2,ipert2,idir1,ipert1)==1)then
             nline=nline+1
             write(iout,'(2(i4,i5),2(1x,f20.10))')idir2,ipert2,idir1,ipert1,&
&             d2cart(1,idir2,ipert2,idir1,ipert1),&
&             d2cart(2,idir2,ipert2,idir1,ipert1)
           end if
         end do
       end do
     end do
   end if
 end if

end subroutine dfpt_dyout
!!***

!!****f* ABINIT/dfpt_gatherdy
!!
!! NAME
!! dfpt_gatherdy
!!
!! FUNCTION
!! Sum (gather) the different parts of the 2nd-order matrix,
!! to get the matrix of second-order derivatives (d2matr)
!! Then, generates the dynamical matrix, not including the masses,
!! but the correct non-cartesian coordinates ( => d2cart)
!!
!! INPUTS
!! becfrnl(3,natom,3*pawbec)=NL frozen contribution to Born Effective Charges (PAW only)
!! berryopt=option for berry phase treatment
!! blkflg(3,mpert,3,mpert)= ( 1 if the element of the dynamical
!!  matrix has been calculated ; 0 otherwise )
!! dyew(2,3,natom,3,natom)=Ewald part of the dyn.matrix
!! dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=frozen wf part of the dyn.matrix (except xc1)
!! dyfrx1(2,3,natom,3,natom)=xc core correction (1) part of the frozen-wf
!!  part of the dynamical matrix.
!! dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!! dyfr_nondiag=1 if dyfrwf is non diagonal with respect to atoms; 0 otherwise
!! dyvdw(2,3,natom,3,natom*usevdw)=vdw DFT-D part of the dynamical matrix
!! d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!! d2nfr(2,3,mpert,3,mpert)=non-frozen wf part of the 2nd-order matr
!! eltcore(6,6)=core contribution to the elastic tensor
!! elteew(6+3*natom,6)=Ewald contribution to the elastic tsenor
!! eltfrhar(6,6)=hartree contribution to the elastic tensor
!! eltfrkin(6,6)=kinetic contribution to the elastic tensor
!! eltfrloc(6+3*natom,6)=local psp contribution to the elastic tensor
!! eltfrnl(6+3*natom,6)=non-local psp contribution to the elastic tensor
!! eltfrxc(6+3*natom,6)=exchange-correlation contribution to the elastic tensor
!! eltvdw(6+3*natom,6*usevdw)=vdw DFT-D part of the elastic tensor
!! gprimd(3,3)=basis vector in the reciprocal space
!! mband=maximum number of bands
!! mpert =maximum number of ipert
!! natom=number of atoms in unit cell
!! ntypat=number of atom types
!! outd2=option for the output of the 2nd-order matrix :
!!  if outd2=1, non-stationary part
!!  if outd2=2, stationary part.
!! pawbec= flag for the computation of frozen part of Born Effective Charges (PAW only)
!! prtbbb=if 1, print the band-by-band decomposition, otherwise, prtbbb=0
!! rfasr= (0=> no acoustic sum rule [asr] imposed), (1=> asr is imposed,
!!  in the democratic way for the effective charges),
!! (2=> asr is imposed, in the aristocratic way for the effective
!!  charges)
!! rfpert(mpert)=define the perturbations
!! rprimd(3,3)=dimensional primitive translations (bohr)
!! typat(natom)=integer label of each type of atom (1,2,...)
!! ucvol=unit cell volume
!! usevdw= flag set to 1 if vdw DFT-D semi-empirical potential is in use
!! zion(ntypat)=charge corresponding to the atom type
!!
!! OUTPUT
!! carflg(3,mpert,3,mpert)= ( 1 if the element of the cartesian
!!  2DTE matrix has been calculated correctly ; 0 otherwise )
!! d2cart(2,3,mpert,3,mpert)=
!!  dynamical matrix, effective charges, dielectric tensor,....
!!  all in cartesian coordinates
!! d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)=
!!  band by band decomposition of Born effective charges
!!  (calculated from phonon-type perturbation) in cartesian coordinates
!! d2matr(2,3,mpert,3,mpert)=2nd-order matrix (masses non included,
!!  no cartesian coordinates : simply second derivatives)
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      asria_calc,asria_corr,cart29,cart39,chneu9
!!
!! SOURCE

subroutine dfpt_gatherdy(becfrnl,berryopt,blkflg,carflg,dyew,dyfrwf,dyfrx1,&
& dyfr_cplex,dyfr_nondiag,dyvdw,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
& eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
& gprimd,mband,mpert,natom,ntypat,outd2,pawbec,pawpiezo,piezofrnl,prtbbb,&
& rfasr,rfpert,rprimd,typat,ucvol,usevdw,zion)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: berryopt,dyfr_cplex,dyfr_nondiag,mband,mpert,natom,ntypat,outd2
 integer,intent(in) :: pawbec,pawpiezo,prtbbb,rfasr,usevdw
 real(dp),intent(in) :: ucvol
!arrays
 integer,intent(in) :: rfpert(mpert),typat(natom)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 integer,intent(out) :: carflg(3,mpert,3,mpert)
 real(dp),intent(in) :: becfrnl(3,natom,3*pawbec)
 real(dp),intent(in) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(in) :: d2nfr(2,3,mpert,3,mpert),dyew(2,3,natom,3,natom)
 real(dp),intent(in) :: dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(in) :: dyfrx1(2,3,natom,3,natom),dyvdw(2,3,natom,3,natom*usevdw)
 real(dp),intent(in) :: eltcore(6,6),elteew(6+3*natom,6)
 real(dp),intent(in) :: eltfrhar(6,6),eltfrkin(6,6),eltfrloc(6+3*natom,6)
 real(dp),intent(in) :: eltfrnl(6+3*natom,6),eltfrxc(6+3*natom,6)
 real(dp),intent(in) :: eltvdw(6+3*natom,6*usevdw),gprimd(3,3)
 real(dp),intent(in) :: piezofrnl(6,3*pawpiezo),rprimd(3,3),zion(ntypat)
 real(dp),intent(out) :: d2cart(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2cart_bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(out) :: d2matr(2,3,mpert,3,mpert)

!Local variables -------------------------
!scalars
 integer :: chneut,iband,iblok,idir,idir1,idir2,ii,ipert,ipert1,ipert2
 integer :: jj,nblok,selectz
 character(len=500) :: message
!arrays
 integer :: flg1(3),flg2(3)
 real(dp) :: vec1(3),vec2(3)
! real(dp) :: ter(3,3) ! this variable appears commented out below
 real(dp),allocatable :: d2tmp(:,:,:,:,:),d2work(:,:,:,:,:),elfrtot(:,:)

! *********************************************************************


!DEBUG
!write(std_out,*)' dfpt_gatherdy : enter '
!write(std_out,*)' outd2,mpert =',outd2,mpert
!write(std_out,*)' blkflg(:,natom+2,:,natom+2)=',blkflg(:,natom+2,:,natom+2)
!ENDDEBUG

 if(outd2/=3)then

!  Initialise the 2nd-derivative matrix
   d2matr(:,:,:,:,:)=0.0_dp

!  Add the non-frozen-part, the
!  Ewald part and the xc1 part of the frozen-wf part
!  Add the vdw part (if any)
   do ipert2=1,mpert
     do idir2=1,3
       do ipert1=1,mpert
         do idir1=1,3
           if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
             do ii=1,2
               d2matr(ii,idir1,ipert1,idir2,ipert2)=&
&               d2nfr(ii,idir1,ipert1,idir2,ipert2)
               if(ipert1<=natom .and. ipert2<=natom) then
                 d2matr(ii,idir1,ipert1,idir2,ipert2)=&
&                 d2matr(ii,idir1,ipert1,idir2,ipert2)+&
&                 dyew(ii,idir1,ipert1,idir2,ipert2)  +&
&                 dyfrx1(ii,idir1,ipert1,idir2,ipert2)
                 if (usevdw==1) then
                   d2matr(ii,idir1,ipert1,idir2,ipert2)=&
&                   d2matr(ii,idir1,ipert1,idir2,ipert2)+&
&                   dyvdw(ii,idir1,ipert1,idir2,ipert2)
                 end if
               end if
             end do
           end if
         end do
       end do
     end do
   end do

!  Add the frozen-wavefunction part
   if (dyfr_nondiag==0) then
     do ipert2=1,natom
       do idir2=1,3
         do idir1=1,3
           if( blkflg(idir1,ipert2,idir2,ipert2)==1 ) then
             d2matr(1:dyfr_cplex,idir1,ipert2,idir2,ipert2)=&
&             d2matr(1:dyfr_cplex,idir1,ipert2,idir2,ipert2)&
&             +dyfrwf(1:dyfr_cplex,idir1,idir2,ipert2,1)
           end if
         end do
       end do
     end do
   else
     do ipert2=1,natom
       do ipert1=1,natom
         do idir2=1,3
           do idir1=1,3
             if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
               d2matr(1:dyfr_cplex,idir1,ipert1,idir2,ipert2)=&
&               d2matr(1:dyfr_cplex,idir1,ipert1,idir2,ipert2)&
&               +dyfrwf(1:dyfr_cplex,idir1,idir2,ipert1,ipert2)
             end if
           end do
         end do
       end do
     end do
   end if

!  Add the frozen-wavefunction part of Born Effective Charges
   if (pawbec==1) then
     ipert2=natom+2
     do idir2=1,3            ! Direction of electric field
       do ipert1=1,natom     ! Atom
         do idir1=1,3        ! Direction of atom
           if(blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
             d2matr(1,idir1,ipert1,idir2,ipert2)=&
&             d2matr(1,idir1,ipert1,idir2,ipert2)+becfrnl(idir1,ipert1,idir2)
           end if
           if(blkflg(idir2,ipert2,idir1,ipert1)==1 ) then
             d2matr(1,idir2,ipert2,idir1,ipert1)=&
&             d2matr(1,idir2,ipert2,idir1,ipert1)+becfrnl(idir1,ipert1,idir2)
           end if
         end do
       end do
     end do
   end if

!  Section for piezoelectric tensor (from electric field response only for PAW)
   if(rfpert(natom+2)==1.and.pawpiezo==1) then
     ipert2=natom+2
     do idir2=1,3            ! Direction of electric field
       do ipert1=natom+3,natom+4     ! Strain
         do idir1=1,3
           ii=idir1+3*(ipert1-natom-3)
           if(blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
             d2matr(1,idir1,ipert1,idir2,ipert2)=&
&             d2nfr(1,idir1,ipert1,idir2,ipert2)+piezofrnl(ii,idir2)
           end if
         end do
       end do
     end do
   end if

!  Section for strain perturbation
   if(rfpert(natom+3)==1 .or. rfpert(natom+4)==1) then
!    Make sure relevant columns of output are nulled
     d2matr(:,:,:,:,natom+3:natom+4)=0.0_dp
!    Accumulate all frozen parts of the elastic tensor
     ABI_ALLOCATE(elfrtot,(6+3*natom,6))
     elfrtot(:,:)=elteew(:,:)+eltfrloc(:,:)+eltfrnl(:,:)+eltfrxc(:,:)
     elfrtot(1:6,1:6)=elfrtot(1:6,1:6)+eltcore(:,:)+eltfrhar(:,:)+eltfrkin(:,:)
     if (usevdw==1) elfrtot(:,:)=elfrtot(:,:)+eltvdw(:,:)

     do ipert2=natom+3,natom+4
       do idir2=1,3
!        Internal strain components first
         do ipert1=1,natom
           do idir1=1,3
             if( blkflg(1,ipert1,idir2,ipert2)==1 ) then
               ii=idir1+6+3*(ipert1-1)
               jj=idir2+3*(ipert2-natom-3)
               d2matr(1,idir1,ipert1,idir2,ipert2)=&
&               d2nfr(1,idir1,ipert1,idir2,ipert2)+elfrtot(ii,jj)
             end if
           end do
         end do
!        Now, electric field - strain mixed derivative (piezoelectric tensor)
         ipert1=natom+2
         do idir1=1,3
           if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
             d2matr(1,idir1,ipert1,idir2,ipert2)=&
&             d2nfr(1,idir1,ipert1,idir2,ipert2)
             if (pawpiezo==1) then
               ii=idir2+3*(ipert2-natom-3)
               d2matr(1,idir1,ipert1,idir2,ipert2)=&
&               d2matr(1,idir1,ipert1,idir2,ipert2)+piezofrnl(ii,idir1)
             end if
           end if
         end do
!        Now, strain-strain 2nd derivatives
         do ipert1=natom+3,natom+4
           do idir1=1,3
             if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
               ii=idir1+3*(ipert1-natom-3)
               jj=idir2+3*(ipert2-natom-3)
               d2matr(1,idir1,ipert1,idir2,ipert2)=&
&               d2nfr(1,idir1,ipert1,idir2,ipert2)+elfrtot(ii,jj)
             end if
           end do
         end do
       end do
     end do
     ABI_DEALLOCATE(elfrtot)
   end if
!  End section for strain perturbation

!  The second-order matrix has been computed.

!  Filter now components smaller in absolute value than 1.0d-20,
!  for automatic testing reasons
   do ipert2=1,mpert
     do idir2=1,3
       do ipert1=1,mpert
         do idir1=1,3
           if( blkflg(idir1,ipert1,idir2,ipert2)==1 ) then
             do ii=1,2
               if( d2matr(ii,idir1,ipert1,idir2,ipert2)**2 < 1.0d-40)then
                 d2matr(ii,idir1,ipert1,idir2,ipert2)=zero
               end if
             end do
           end if
         end do
       end do
     end do
   end do

!  DEBUG
!  write(std_out,*)' d2matr '
!  ipert2=natom+2
!  do idir2=1,3
!  ipert1=natom+2
!  do idir1=1,3
!  write(std_out,'(4i4,2es16.6)' )idir1,ipert1,idir2,ipert2,&
!  &       d2matr(1,idir1,ipert1,idir2,ipert2),&
!  &       d2matr(2,idir1,ipert1,idir2,ipert2)
!  end do
!  end do
!  ENDDEBUG

!  Cartesian coordinates transformation
   iblok=1 ; nblok=1

!  In the case of finite electric field, the convention for the
!  direction of the electric field perturbation was NOT the usual convention ...
!  So, there is a transformation to the usual conventions
!  to be done first ...
   if((berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. berryopt==14 .or. berryopt==16 .or. berryopt==17 )  &
&   .and. minval(abs(blkflg(:,natom+2,:,natom+2)))/=0)then   !!HONG  need to check for fixed D and E calculation
     if(minval(abs(blkflg(:,natom+2,:,natom+2)-1))/=0)then
       write(message,'(5a)')&
&       '  In case of finite electric field, and electric field perturbation,',ch10,&
&       '  the three directions for the perturbations must be treated.',ch10,&
&       '  Action : set idir to 1 1 1, or forget about finite electric field.'
       MSG_ERROR(message)
     end if
     do ipert=1,mpert
       do idir=1,3
         do ii=1,2
           vec1(:)=d2matr(ii,idir,ipert,:,natom+2)
           flg1(:)=blkflg(idir,ipert,:,natom+2)
           call cart39(flg1,flg2,gprimd,1,1,rprimd,vec1,vec2)
           d2matr(ii,idir,ipert,:,natom+2)=vec2(:)*two_pi
           blkflg(idir,ipert,:,natom+2)=flg2(:)
         end do
       end do
     end do
     do ipert=1,mpert
       do idir=1,3
         do ii=1,2
           vec1(:)=d2matr(ii,:,natom+2,idir,ipert)
           flg1(:)=blkflg(:,natom+2,idir,ipert)
           call cart39(flg1,flg2,gprimd,1,1,rprimd,vec1,vec2)
           d2matr(ii,:,natom+2,idir,ipert)=vec2(:)*two_pi
           blkflg(:,natom+2,idir,ipert)=flg2(:)
         end do
       end do
     end do
!    Also to be done, a change of sign, that I do not understand (XG071110)
!    Perhaps due to d/dk replacing id/dk ? !
     d2matr(:,:,natom+2,:,natom+2)=-d2matr(:,:,natom+2,:,natom+2)
   end if

   call cart29(blkflg,d2matr,carflg,d2cart,&
&   gprimd,iblok,mpert,natom,nblok,ntypat,rprimd,typat,ucvol,zion)

!  Band by band decomposition of the Born effective charges
   if(prtbbb==1)then
     ABI_ALLOCATE(d2work,(2,3,mpert,3,mpert))
     ABI_ALLOCATE(d2tmp,(2,3,mpert,3,mpert))
     do iband=1,mband
       d2work(:,:,:,:,:)=0.0_dp
       d2tmp(:,:,:,:,:)=0.0_dp
       d2work(:,:,natom+2,:,:) = d2bbb(:,:,:,:,iband,iband)
       call cart29(blkflg,d2work,carflg,d2tmp,&
&       gprimd,iblok,mpert,natom,nblok,ntypat,rprimd,typat,ucvol,zion)

!      Remove the ionic part
       do ipert1=1,natom
         do idir1=1,3
           d2tmp(1,idir1,natom+2,idir1,ipert1) = &
&           d2tmp(1,idir1,natom+2,idir1,ipert1) - zion(typat(ipert1))
         end do
       end do

       d2cart_bbb(:,:,:,:,iband,iband)=d2tmp(:,:,natom+2,:,:)

     end do
     ABI_DEALLOCATE(d2tmp)
     ABI_DEALLOCATE(d2work)
   end if ! prtbbb==1

!
!  Now, the cartesian elements are ready for output
!  carflg give the information on their correctness
 end if

!  Imposition of the ASR on the analytical part of the DynMat
!  Assume that if rfasr/=0, the whole cartesian matrix is correct
 if(rfasr/=0)then

   ABI_ALLOCATE(d2work,(2,3,mpert,3,mpert))
   call asria_calc(rfasr,d2work,d2cart,mpert,natom)
!  The following line imposes ASR:
   call asria_corr(rfasr,d2work,d2cart,mpert,natom)

   ABI_DEALLOCATE(d2work)

!  Imposition of the ASR on the effective charges.
   if(rfpert(natom+2)==1)then
     chneut=rfasr
     selectz=0
     call chneu9(chneut,d2cart,mpert,natom,ntypat,selectz,typat,zion)
   end if

 end if

!Additional operations on cartesian strain derivatives
 if(rfpert(natom+3)==1 .or. rfpert(natom+4)==1) then
!  Impose zero-net-force condition on internal strain tensor
   do ipert2=natom+3,natom+4
     do idir2=1,3
       vec1(:)=0.0_dp
       do ipert1=1,natom
         do idir1=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1) then
             vec1(idir1)=vec1(idir1)+d2cart(1,idir1,ipert1,idir2,ipert2)
           end if
         end do
       end do
       vec1(:)=vec1(:)/dble(natom)
       do ipert1=1,natom
         do idir1=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1) then
!            Note minus sign to convert gradients to forces
             d2cart(1,idir1,ipert1,idir2,ipert2)=&
&             -(d2cart(1,idir1,ipert1,idir2,ipert2)-vec1(idir1))
           end if
         end do
       end do
     end do
   end do
!  Divide strain 2nd deriviative by ucvol to give elastic tensor
   do ipert2=natom+3,natom+4
     do idir2=1,3
       do ipert1=natom+3,natom+4
         do idir1=1,3
           if(carflg(idir1,ipert1,idir2,ipert2)==1) then
             d2cart(1,idir1,ipert1,idir2,ipert2)=&
&             d2cart(1,idir1,ipert1,idir2,ipert2)/ucvol
           end if
         end do
       end do
     end do
   end do
 end if

!calculate Born effective charges from electric field perturbation
!do ipert1=1,natom
!ter(:,:)=zero
!do idir1=1,3
!do ii=1,3
!do jj=1,3
!if(abs(gprimd(idir1,ii))>1.0d-10)then
!ter(idir1,ii)=ter(idir1,ii)+ d2nfr(1,idir1,natom+2,jj,ipert1)*gprimd(jj,ii)
!endif
!enddo
!enddo
!add zion to bec
!ter(idir1,idir1)=ter(idir1,idir1)+zion(typat(ipert1))
!enddo
!d2cart(1,:,ipert1,:,natom+2)=ter(:,:)
!enddo
!carflg(:,1:natom,:,natom+2)=1

!Born effective charges from phonon perturbation
!do ipert1=1,natom
!ter(:,:)=zero
!do idir1=1,3
!do ii=1,3
!do jj=1,3
!if(abs(gprimd(idir1,ii))>1.0d-10)then
!ter(idir1,ii)=ter(idir1,ii)+ d2nfr(1,jj,ipert1,idir1,natom+2)*gprimd(jj,ii)
!endif
!enddo
!enddo
!! add zion to bec
!ter(idir1,idir1)=ter(idir1,idir1)+zion(typat(ipert1))
!enddo
!d2cart(1,:,natom+2,:,ipert1)=ter(:,:)
!enddo
!carflg(:,natom+2,:,1:natom)=1


!!Dielectric constant
!do ii=1,3
!do jj=1,3
!ter(ii,jj)=d2nfr(1,ii,natom+2,jj,natom+2)
!end do
!end do
!ter(:,:)=pi*four*ter(:,:)/ucvol
!
!do ii=1,3
!ter(ii,ii)=ter(ii,ii)+one
!end do
!d2cart(1,:,natom+2,:,natom+2)=ter(:,:)
!carflg(:,natom+2,:,1:natom+2)=1

!DEBUG
!Debugging, but needed to make the stuff work on the IBM Dirac ? !
!write(std_out,*)' d2cart '
!ipert2=natom+2
!do idir2=1,3
!ipert1=natom+2
!do idir1=1,3
!write(std_out,'(5i4,2d20.10)' )idir1,ipert1,idir2,ipert2,&
!&      carflg(idir1,ipert1,idir2,ipert2),&
!&      d2cart(1,idir1,ipert1,idir2,ipert2),&
!&      d2cart(2,idir1,ipert1,idir2,ipert2)
!end do
!end do
!ENDDEBUG

!DEBUG
!write(std_out,*)' dfpt_gatherdy : exit '
!ENDDEBUG

end subroutine dfpt_gatherdy
!!***

!!****f* ABINIT/dfpt_dyfro
!! NAME
!! dfpt_dyfro
!!
!! FUNCTION
!! Compute the different parts of the frozen-wavefunction part of
!! the dynamical matrix, except the non-local one, computed previously.
!! Also (when installed) symmetrize the different part and their sum.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dyfr_cplex=1 if dyfrnl is real, 2 if it is complex
!!  dyfr_nondiag=1 if dyfrnl and dyfrwf are non diagonal with respect to atoms; 0 otherwise
!!  gmet(3,3)=reciprocal space metric (bohr^-2)
!!  gprimd(3,3)=dimensional primitive translations for reciprocal space (bohr**-1)
!!  gsqcut=cutoff on G^2 based on ecut
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!   under symmetry operations (computed in symatm)
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  natom=number of atoms in unit cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!     see ~abinit/doc/variables/vargs.htm#ngfft
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  qphon(3)=wavevector of the phonon
!!  rhog(2,nfft)=electron density in G space
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  symq(4,2,nsym)=1 if symmetry preserves present qpoint. From littlegroup_q
!!  symrec(3,3,nsym)=symmetries in reciprocal space
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vlspl(mqgrid,2,ntypat)=q^2 v(q) spline for each type of atom.
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real
!!   space--only used when n1xccc/=0
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced coordinates for atoms in unit cell
!!
!! OUTPUT
!!  dyfrlo(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (local only)
!!  dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!                    frozen wavefunctions part of the dynamical matrix
!!                    (local + non-local)
!!                    If NCPP, it depends on one atom
!!                    If PAW,  it depends on two atoms
!!  dyfrxc(3,3,natom)=frozen wavefunctions part of the dynamical matrix
!!                    (non-linear xc core correction)
!!
!! SIDE EFFECTS
!! Input/Output
!!  dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)=
!!                    frozen wavefunctions part of the dynamical matrix
!!                    (non-local only)
!!                    If NCPP, it depends on one atom
!!                    If PAW,  it depends on two atoms
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      atm2fft,dfpt_sydy,fourdp,mkcore,mklocl_recipspace,timab,zerosym
!!
!! SOURCE

subroutine dfpt_dyfro(atindx1,dyfrnl,dyfrlo,dyfrwf,dyfrxc,dyfr_cplex,dyfr_nondiag,&
&  gmet,gprimd,gsqcut,indsym,mgfft,mpi_enreg,mqgrid,natom,nattyp,&
&  nfft,ngfft,nspden,nsym,ntypat,n1xccc,n3xccc,psps,pawtab,ph1d,qgrid,&
&  qphon,rhog,rprimd,symq,symrec,typat,ucvol,usepaw,vlspl,vxc,&
&  xcccrc,xccc1d,xccc3d,xred)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: dyfr_cplex,dyfr_nondiag,mgfft,mqgrid,n1xccc,n3xccc,natom,nfft,nspden
 integer,intent(in) :: nsym,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx1(natom),indsym(4,nsym,natom),nattyp(ntypat)
 integer,intent(in) :: ngfft(18),symq(4,2,nsym),symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: qgrid(mqgrid),qphon(3),rhog(2,nfft),rprimd(3,3)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat),vxc(nfft,nspden)
 real(dp),intent(in) :: xccc1d(n1xccc,6,ntypat),xcccrc(ntypat),xred(3,natom)
 real(dp),intent(inout) :: dyfrnl(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag),xccc3d(n3xccc)
 real(dp),intent(out) :: dyfrlo(3,3,natom),dyfrwf(dyfr_cplex,3,3,natom,1+(natom-1)*dyfr_nondiag)
 real(dp),intent(inout) :: dyfrxc(3,3,natom) !vz_i
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 logical, parameter :: do_final_sym=.true.
 integer :: iatom,jatom,n1,n2,n3,optatm,optdyfr,opteltfr,optgr,option
 integer :: optn,optn2,optstr,optv
 real(dp) :: eei
!arrays
 integer  :: qprtrb(3)
 real(dp) :: dummy6(6),dum_strn(6),dum_strv(6)
 real(dp) :: tsec(2),vprtrb(2)
 real(dp) :: dum_atmrho(0),dum_atmvloc(0),dum_gauss(0),dum_grn(0),dum_grv(0),dum_eltfrxc(0)
 real(dp),allocatable :: dyfrlo_tmp1(:,:,:),dyfrlo_tmp2(:,:,:,:,:),dyfrsym_tmp(:,:,:,:,:)
 real(dp),allocatable :: gr_dum(:,:),v_dum(:),vxctotg(:,:)

! *************************************************************************

 if(nspden==4)then
   MSG_WARNING('dfpt_dyfro : DFPT with nspden=4 works at the moment just for insulators and norm-conserving psp!')
 end if

 n1=ngfft(1); n2=ngfft(2); n3=ngfft(3)

 if (usepaw==1 .or. psps%nc_xccc_gspace==1) then

!  PAW or NC with nc_xccc_gspace: compute local psp and core charge contribs together
!  in reciprocal space
!  -----------------------------------------------------------------------
   call timab(563,1,tsec)
   if (n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,1,ngfft,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_DEALLOCATE(v_dum)
   else
     ABI_ALLOCATE(vxctotg,(0,0))
   end if
   optatm=0;optdyfr=1;optgr=0;optstr=0;optv=1;optn=n3xccc/nfft;optn2=1;opteltfr=0
   call atm2fft(atindx1,dum_atmrho,dum_atmvloc,dyfrxc,dyfrlo,dum_eltfrxc,&
&   dum_gauss,gmet,gprimd,dum_grn,dum_grv,gsqcut,mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,&
&   optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,qgrid,qprtrb,rhog,&
&   dum_strn,dum_strv,ucvol,usepaw,vxctotg,vxctotg,vxctotg,vprtrb,vlspl)
   ABI_DEALLOCATE(vxctotg)
   if (n3xccc==0) dyfrxc=zero
 else

!  Norm-conserving: compute local psp contribution in reciprocal space
!  and core charge contribution in real space
!  -----------------------------------------------------------------------
   option=4
   ABI_ALLOCATE(dyfrlo_tmp1,(3,3,natom))
   ABI_ALLOCATE(gr_dum,(3,natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl_recipspace(dyfrlo_tmp1,eei,gmet,gprimd,&
&   gr_dum,gsqcut,dummy6,mgfft,mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,&
&   ntypat,option,ph1d,qgrid,qprtrb,rhog,ucvol,vlspl,vprtrb,v_dum)
   do iatom=1,natom
!    Reestablish correct order of atoms
     dyfrlo(1:3,1:3,atindx1(iatom))=dyfrlo_tmp1(1:3,1:3,iatom)
   end do
   ABI_DEALLOCATE(dyfrlo_tmp1)
   ABI_DEALLOCATE(v_dum)
   if(n1xccc/=0)then
     call mkcore(dummy6,dyfrxc,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,&
&     n1,n1xccc,n2,n3,option,rprimd,typat,ucvol,vxc,xcccrc,xccc1d,xccc3d,xred)
   end if
   ABI_DEALLOCATE(gr_dum)
 end if

!Symmetrize dynamical matrix explicitly for given space group:

!Symmetrize local part of the dynamical matrix dyfrlo:
 ABI_ALLOCATE(dyfrsym_tmp,(1,3,3,natom,1))
 ABI_ALLOCATE(dyfrlo_tmp2,(1,3,3,natom,1))
 dyfrsym_tmp(1,:,:,:,1)=dyfrlo(:,:,:)
 call dfpt_sydy(1,dyfrsym_tmp,indsym,natom,0,nsym,qphon,dyfrlo_tmp2,symq,symrec)
 dyfrlo(:,:,:)=dyfrlo_tmp2(1,:,:,:,1)
 if (do_final_sym) then
   dyfrsym_tmp(1,:,:,:,1)=dyfrxc(:,:,:)
   call dfpt_sydy(1,dyfrsym_tmp,indsym,natom,0,nsym,qphon,dyfrlo_tmp2,symq,symrec)
   dyfrxc(:,:,:)=dyfrlo_tmp2(1,:,:,:,1)
 end if
 ABI_DEALLOCATE(dyfrsym_tmp)
 ABI_DEALLOCATE(dyfrlo_tmp2)

!Symmetrize nonlocal part of the dynamical matrix dyfrnl:
!atindx1 is used to reestablish the correct order of atoms
 if (dyfr_nondiag==0) then
   ABI_ALLOCATE(dyfrsym_tmp,(dyfr_cplex,3,3,natom,1))
   do iatom=1,natom
     dyfrsym_tmp(:,:,:,atindx1(iatom),1)=dyfrnl(:,:,:,iatom,1)
   end do
 else
   ABI_ALLOCATE(dyfrsym_tmp,(dyfr_cplex,3,3,natom,natom))
   do jatom=1,natom
     do iatom=1,natom
       dyfrsym_tmp(:,:,:,atindx1(iatom),atindx1(jatom))=dyfrnl(:,:,:,iatom,jatom)
     end do
   end do
 end if
 call dfpt_sydy(dyfr_cplex,dyfrsym_tmp,indsym,natom,dyfr_nondiag,nsym,qphon,dyfrnl,symq,symrec)
 ABI_DEALLOCATE(dyfrsym_tmp)

!Collect local, nl xc core, and non-local part
!of the frozen wf dynamical matrix.
 dyfrwf(:,:,:,:,:)=dyfrnl(:,:,:,:,:)
 if (dyfr_nondiag==0) then
   dyfrwf(1,:,:,:,1)=dyfrwf(1,:,:,:,1)+dyfrlo(:,:,:)+dyfrxc(:,:,:)
 else
   do iatom=1,natom
     dyfrwf(1,:,:,iatom,iatom)=dyfrwf(1,:,:,iatom,iatom)+dyfrlo(:,:,iatom)+dyfrxc(:,:,iatom)
   end do
 end if

end subroutine dfpt_dyfro
!!***

!!****f* ABINIT/dfpt_dyxc1
!! NAME
!! dfpt_dyxc1
!!
!! FUNCTION
!! Compute 2nd-order non-linear xc core-correction (part1) to the dynamical matrix.
!! In case of derivative with respect to k or electric field perturbation,
!! the 1st-order local potential vanishes.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms
!!  gmet(3,3)=metrix tensor in G space in Bohr**-2.
!!  gsqcut=cutoff value on G**2 for sphere inside fft box.
!!  ixc= choice of exchange-correlation scheme
!!  kxc(nfft,nkxc)=first-order derivative of the xc potential
!!    if (nkxc=1) LDA kxc(:,1)= d2Exc/drho2
!!    if (nkxc=2) LDA kxc(:,1)= d2Exc/drho_up drho_up
!!                    kxc(:,2)= d2Exc/drho_up drho_dn
!!                    kxc(:,3)= d2Exc/drho_dn drho_dn
!!    if (nkxc=7) GGA kxc(:,1)= d2Exc/drho2
!!                    kxc(:,2)= 1/|grad(rho)| dExc/d|grad(rho)|
!!                    kxc(:,3)= 1/|grad(rho)| d2Exc/d|grad(rho)| drho
!!                    kxc(:,4)= 1/|grad(rho)| * d/d|grad(rho)| ( 1/|grad(rho)| dExc/d|grad(rho)| )
!!                    kxc(:,5)= gradx(rho)
!!                    kxc(:,6)= grady(rho)
!!                    kxc(:,7)= gradz(rho)
!!    if (nkxc=19) spin-polarized GGA case (same as nkxc=7 with up and down components)
!!  mgfft=maximum size of 1D FFTs
!!  mpert=maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=number of grid pts in q array for f(q) spline.
!!  natom=number of atoms in cell.
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(3)=fft grid dimensions.
!!  nkxc=second dimension of the kxc array
!!   (=1 for non-spin-polarized case, =3 for spin-polarized case)
!!  nmxc= if true, handle density/potential as non-magnetic (even if it is)
!!  nspden=number of spin-density components
!!  ntypat=number of types of atoms in cell.
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phase information
!!  qgrid(mqgrid)=q grid for spline from 0 to qmax.
!!  qphon(3)=wavevector of the phonon
!!  rfdir(3)=array that define the directions of perturbations
!!  rfpert(mpert)=array defining the type of perturbations that have to be computed
!!  rprimd(3,3)=dimensional primitive translation vectors (bohr)
!!  timrev=1 if time-reversal preserves the q wavevector; 0 otherwise.
!!  typat(natom)=integer type for each atom in cell
!!  ucvol=unit cell volume (bohr**3).
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xccc1d(n1xccc,6,ntypat)=1D core charge function and five derivatives,
!!   for each type of atom, from psp
!!  xred(3,natom)=fractional coordinates for atoms in unit cell
!!
!! OUTPUT!  non_magnetic_xc= if true, handle density/potential as non-magnetic (even if it is)
!!  blkflgfrx1(3,natom,3,natom)=flag to indicate whether an element has been computed or not
!!  dyfrx1(2,3,natom,3,natom)=2nd-order non-linear xc
!!    core-correction (part1) part of the dynamical matrix
!!
!! PARENTS
!!      respfn
!!
!! CHILDREN
!!      dfpt_atm2fft,dfpt_mkcore,dfpt_mkvxc,dfpt_mkvxc_noncoll,dotprod_vn,timab
!!      xmpi_sum
!!
!! SOURCE

subroutine dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,ixc,kxc,mgfft,mpert,mpi_enreg,mqgrid,&
&          natom,nfft,ngfft,nkxc,nmxc,nspden,ntypat,n1xccc,psps,pawtab,&
&          ph1d,qgrid,qphon,rfdir,rfpert,rprimd,timrev,typat,ucvol,usepaw,xcccrc,xccc1d,xred,rhor,vxc)

 use m_cgtools,       only : dotprod_vn
 use m_atm2fft,       only : dfpt_atm2fft
 use m_dfpt_mkvxc,    only : dfpt_mkvxc, dfpt_mkvxc_noncoll

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ixc,mgfft,mpert,mqgrid,n1xccc,natom,nfft,nkxc,nspden,ntypat
 integer,intent(in) :: timrev,usepaw
 logical,intent(in) :: nmxc
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: atindx(natom),ngfft(18),rfdir(3),rfpert(mpert),typat(natom)
 real(dp),intent(in) :: gmet(3,3),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid),qphon(3)
 real(dp),intent(in) :: rprimd(3,3),xccc1d(n1xccc,6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom)
 integer,intent(out) :: blkflgfrx1(3,natom,3,natom)
 real(dp),intent(out) :: dyfrx1(2,3,natom,3,natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)
!optional
 real(dp),optional,intent(in) :: rhor(nfft,nspden)
 real(dp),optional,intent(in) :: vxc(nfft,nspden)

!Local variables-------------------------------
!scalars
 integer :: cplex,iat1,iatom1,iatom2,idir1,idir2,ierr,ifft,my_natom,comm_atom
 integer :: n1,n2,n3,n3xccc,nfftot,option,upperdir,optnc
 logical :: paral_atom
 real(dp) :: valuei,valuer
!arrays
 integer,pointer :: my_atmtab(:)
 real(dp) :: tsec(2),gprimd_dummy(3,3)
 real(dp) :: dum_nhat(0)
 real(dp),allocatable :: rhor1(:,:),vxc10(:,:),xcccwk1(:),xcccwk2(:)
! *********************************************************************

 call timab(182,1,tsec)

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 nfftot=n1*n2*n3

!Set up parallelism over atoms
 my_natom=mpi_enreg%my_natom
 my_atmtab=>mpi_enreg%my_atmtab
 comm_atom=mpi_enreg%comm_atom
 paral_atom=(my_natom/=natom)

!Zero out the output arrays :
 blkflgfrx1(:,:,:,:)=0
 dyfrx1(:,:,:,:,:)=zero

 cplex=2-timrev ; n3xccc=nfft
 ABI_ALLOCATE(vxc10,(cplex*nfft,nspden))


!Loop on the perturbation j1
 do iat1=1,my_natom
   iatom1=iat1; if(paral_atom)iatom1=my_atmtab(iat1)
   do idir1=1,3

!    Compute the derivative of the core charge with respect to j1
     ABI_ALLOCATE(xcccwk1,(cplex*n3xccc))

!    PAW or NC with nc_xccc_gspace: 1st-order core charge in reciprocal space
     if (usepaw==1 .or. psps%nc_xccc_gspace==1) then
       call dfpt_atm2fft(atindx,cplex,gmet,gprimd_dummy,gsqcut,idir1,iatom1,&
&       mgfft,mqgrid,natom,1,nfft,ngfft,ntypat,ph1d,qgrid,&
&       qphon,typat,ucvol,usepaw,xred,psps,pawtab,atmrhor1=xcccwk1,optn2_in=1)

!      Norm-conserving psp: 1st-order core charge in real space
     else
       call dfpt_mkcore(cplex,idir1,iatom1,natom,ntypat,n1,n1xccc,&
&       n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk1,xred)
     end if

!    Compute the corresponding potential
     option=0
     ABI_ALLOCATE(rhor1,(cplex*nfft,nspden))
     rhor1=zero
!FR SPr EB Non-collinear magnetism
     if (nspden==4.and.present(rhor).and.present(vxc)) then
       optnc=1
       call dfpt_mkvxc_noncoll(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,dum_nhat,0,dum_nhat,0,dum_nhat,0,nkxc,&
&       nmxc,nspden,n3xccc,optnc,option,qphon,rhor,rhor1,rprimd,0,vxc,vxc10,xcccwk1)
     else
       call dfpt_mkvxc(cplex,ixc,kxc,mpi_enreg,nfft,ngfft,dum_nhat,0,dum_nhat,0,nkxc,&
&       nmxc,nspden,n3xccc,option,qphon,rhor1,rprimd,0,vxc10,xcccwk1)
     end if
     ABI_DEALLOCATE(rhor1)
     ABI_DEALLOCATE(xcccwk1)

!    vxc10 will couple with xcccwk2, that behaves like
!    a total density (ispden=1). Only the spin-up + spin-down
!    average of vxc10 is needed.
     if (nspden/=1)then
       do ifft=1,cplex*nfft
         vxc10(ifft,1)=(vxc10(ifft,1)+vxc10(ifft,2))*half
       end do
     end if

!    Loop on the perturbation j2
     do iatom2=1,iatom1
       upperdir=3
       if(iatom1==iatom2)upperdir=idir1
       do idir2=1,upperdir
         if( (rfpert(iatom1)==1 .and. rfdir(idir1) == 1) .or. &
&         (rfpert(iatom2)==1 .and. rfdir(idir2) == 1)    )then

!          Compute the derivative of the core charge with respect to j2
           ABI_ALLOCATE(xcccwk2,(cplex*n3xccc))

!          PAW or NC with nc_xccc_gspace: 1st-order core charge in reciprocal space
           if (usepaw==1 .or. psps%nc_xccc_gspace==1) then
             call dfpt_atm2fft(atindx,cplex,gmet,gprimd_dummy,gsqcut,idir2,iatom2,&
&             mgfft,mqgrid,natom,1,nfft,ngfft,ntypat,ph1d,qgrid,&
&             qphon,typat,ucvol,usepaw,xred,psps,pawtab,atmrhor1=xcccwk2,optn2_in=1)

!            Norm-conserving psp: 1st-order core charge in real space
           else
             call dfpt_mkcore(cplex,idir2,iatom2,natom,ntypat,n1,n1xccc,&
&             n2,n3,qphon,rprimd,typat,ucvol,xcccrc,xccc1d,xcccwk2,xred)
           end if

!          Get the matrix element j1,j2

           call dotprod_vn(cplex,xcccwk2,valuer,valuei,nfft,nfftot,1,2,vxc10,ucvol)

           ABI_DEALLOCATE(xcccwk2)

           dyfrx1(1,idir1,iatom1,idir2,iatom2)= valuer
           dyfrx1(2,idir1,iatom1,idir2,iatom2)= valuei
           dyfrx1(1,idir2,iatom2,idir1,iatom1)= valuer
           dyfrx1(2,idir2,iatom2,idir1,iatom1)=-valuei
           blkflgfrx1(idir1,iatom1,idir2,iatom2)=1
           blkflgfrx1(idir2,iatom2,idir1,iatom1)=1
         end if
       end do
     end do
   end do
 end do

 if (paral_atom) then
   call timab(48,1,tsec)
   call xmpi_sum(dyfrx1,comm_atom,ierr)
   call xmpi_sum(blkflgfrx1,comm_atom,ierr)
   call timab(48,2,tsec)
 end if

 ABI_DEALLOCATE(vxc10)

 call timab(182,2,tsec)

end subroutine dfpt_dyxc1
!!***

end module m_respfn_driver
!!***

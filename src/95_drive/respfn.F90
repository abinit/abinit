!{\src2tex{textfont=tt}}
!!****f* ABINIT/respfn
!! NAME
!! respfn
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of Response functions.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, DRH, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!      alloc_hamilt_gpu,atm2fft,check_kxc,chkpawovlp,chkph3,d2frnl,d2sym3
!!      ddb_io_out,dealloc_hamilt_gpu,dfpt_dyfro,dfpt_dyout,dfpt_dyxc1
!!      dfpt_eltfrhar,dfpt_eltfrkin,dfpt_eltfrloc,dfpt_eltfrxc,dfpt_ewald
!!      dfpt_gatherdy,dfpt_looppert,dfpt_phfrq,dfpt_prtph,ebands_free
!!      efmasdeg_free_array,efmasfr_free_array,eig2tot,eigen_meandege
!!      elph2_fanddw,elt_ewald,exit_check,fourdp,getcut,getph,hdr_free,hdr_init
!!      hdr_update,initrhoij,initylmg,inwffil,irreducible_set_pert,kpgio
!!      littlegroup_q,mkcore,mklocl,mkrho,newocc,nhatgrid,outddbnc,paw_an_free
!!      paw_an_init,paw_an_nullify,paw_gencond,paw_ij_free,paw_ij_init
!!      paw_ij_nullify,pawdenpot,pawdij,pawexpiqr,pawfgr_destroy,pawfgr_init
!!      pawfgrtab_free,pawfgrtab_init,pawinit,pawmknhat,pawpuxinit
!!      pawrhoij_alloc,pawrhoij_bcast,pawrhoij_copy,pawrhoij_free
!!      pawrhoij_nullify,pawtab_get_lsize,prteigrs,psddb8,pspini,q0dy3_apply
!!      q0dy3_calc,read_rhor,rhohxc,setsym,setsymrhoij,setup1,status,symdij
!!      symmetrize_xred,sytens,timab,transgrid,vdw_dftd2,vdw_dftd3,wffclose
!!      wings3,wrtloctens,wrtout,xmpi_bcast
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine respfn(codvsn,cpui,dtfil,dtset,etotal,iexit,&
&  mkmems,mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,results_respfn,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_efmas_defs
 use m_profiling_abi
 use m_xmpi
 use m_exit
 use m_wffile
 use m_errors
 use m_ebands
 use m_results_respfn
 use m_hdr

 use m_dynmat,      only : chkph3, d2sym3, q0dy3_apply, q0dy3_calc, wings3, dfpt_phfrq
 use m_ddb,         only : psddb8, DDB_VERSION
 use m_efmas,       only : efmasdeg_free_array, efmasfr_free_array
 use m_wfk,         only : wfk_read_eigenvalues
 use m_ioarr,       only : read_rhor
 use m_pawang,      only : pawang_type
 use m_pawrad,      only : pawrad_type
 use m_pawtab,      only : pawtab_type,pawtab_get_lsize
 use m_paw_an,      only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify
 use m_paw_ij,      only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify
 use m_pawfgrtab,   only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_copy, &
&                          pawrhoij_bcast, pawrhoij_nullify
 use m_pawdij,      only : pawdij, symdij
 use m_pawfgr,      only : pawfgr_type, pawfgr_init, pawfgr_destroy
 use m_paw_finegrid,only : pawexpiqr
 use m_paw_dmft,    only : paw_dmft_type

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'respfn'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_41_xc_lowlevel
#if defined HAVE_GPU_CUDA
 use interfaces_52_manage_cuda
#endif
 use interfaces_53_ffts
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_61_occeig
 use interfaces_64_psp
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_72_response
 use interfaces_79_seqpar_mpi
 use interfaces_95_drive, except_this_one => respfn
!End of the abilint section

 implicit none

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
 integer :: analyt,ask_accurate,band_index,bantot,bdeigrf,choice,coredens_method,cplex
 integer :: dim_eig2nkq,dim_eigbrd,dyfr_cplex,dyfr_nondiag,fullinit,gnt_option
 integer :: gscase,has_dijnd,has_kxc,iatom,iatom_tot,iband,idir,ider,ierr,ifft,ii,ikpt,indx
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer :: initialized,ipert,ipert2,ireadwf0,iscf,iscf_eff,ispden,isppol
 integer :: itypat,izero,mcg,me,mgfftf,mk1mem,mkqmem,mpert,mu
 integer :: my_natom,n1,natom,n3xccc,nband_k,nblok,nfftf,nfftot,nfftotf,nhatdim,nhatgrdim
 integer :: nkpt_eff,nkpt_max,nkpt_rbz,nkxc,nkxc1,nspden_rhoij,ntypat,nzlmopt,openexit
 integer :: optcut,option,optgr0,optgr1,optgr2,optorth,optrad
 integer :: optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv
 integer :: outd2,pawbec,pawpiezo,prtbbb,psp_gencond,qzero,rdwr,rdwrpaw
 integer :: req_cplex_dij,rfasr,rfddk,rfelfd,rfphon,rfstrs,rfuser,rf2_dkdk,rf2_dkde
 integer :: spaceworld,sumg0,sz1,sz2,tim_mkrho,timrev,usecprj,usevdw
 integer :: usexcnhat,use_sym,vloc_method
 logical :: has_full_piezo,has_allddk,paral_atom,qeq0,use_nhat_gga,call_pawinit
 real(dp) :: boxcut,compch_fft,compch_sph,cpus,ecore,ecut_eff,ecutdg_eff,ecutf
 real(dp) :: eei,eew,ehart,eii,ek,enl,entropy,enxc
 real(dp) :: epaw,epawdc,etot,evdw,fermie,gsqcut,gsqcut_eff,gsqcutc_eff,qphnrm,residm
 real(dp) :: tolwfr
 real(dp) :: ucvol,vxcavg
 character(len=fnlen) :: dscrpt
 character(len=500) :: message
 type(ebands_t) :: bstruct
 type(hdr_type) :: hdr,hdr_fine,hdr0,hdr_den
 type(paw_dmft_type) :: paw_dmft
 type(pawfgr_type) :: pawfgr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 integer :: ddkfil(3),ngfft(18),ngfftf(18),rfdir(3),rf2_dirs_from_rfpert_nl(3,3)
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:,:,:,:),blkflgfrx1(:,:,:,:),blkflg1(:,:,:,:)
 integer,allocatable :: blkflg2(:,:,:,:),carflg(:,:,:,:),clflg(:,:),indsym(:,:,:)
 integer,allocatable :: irrzon(:,:,:),kg(:,:),l_size_atm(:),nattyp(:),npwarr(:)
 integer,allocatable :: pertsy(:,:),rfpert(:),rfpert_nl(:,:,:,:,:,:),symq(:,:,:),symrec(:,:,:)
 real(dp) :: dum_gauss(0),dum_dyfrn(0),dum_dyfrv(0),dum_eltfrxc(0)
 real(dp) :: dum_grn(0),dum_grv(0),dum_rhog(0),dum_vg(0)
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qphon(3)
 real(dp) :: rmet(3,3),rprimd(3,3),strn_dummy6(6),strv_dummy6(6),strsxc(6),tsec(2)
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp),allocatable :: amass(:),becfrnl(:,:,:),cg(:,:),d2bbb(:,:,:,:,:,:),d2cart(:,:,:,:,:)
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
 type(efmasfr_type),allocatable :: efmasfr(:,:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij_read(:)
! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(132,1,tsec)
 call timab(133,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'init          ')

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
 rfddk=dtset%rfddk   ; rfelfd=dtset%rfelfd
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

 call status(0,dtfil%filstat,iexit,level,'call setup1   ')

!LIKELY TO BE TAKEN AWAY
 initialized=0
 ek=zero ; ehart=zero ; enxc=zero ; eei=zero ; enl=zero
 eii=zero ; eew=zero ; ecore=zero

!Set up for iterations
 ABI_ALLOCATE(amass,(natom))
 call setup1(dtset%acell_orig(1:3,1),amass,dtset%amu_orig(:,1),bantot,dtset,&
& ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
& natom,ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
& response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!Define the set of admitted perturbations
 mpert=natom+6
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

 if(rfuser==1.or.rfuser==3)rfpert(natom+5)=1
 if(rfuser==2.or.rfuser==3)rfpert(natom+6)=1

 qeq0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2<1.d-14)

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell
 me = xmpi_comm_rank(spaceworld)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 call status(0,dtfil%filstat,iexit,level,'call kpgio(1) ')
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,&
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
   call status(0,dtfil%filstat,iexit,level,'call initylmg ')
   option=0
   if (rfstrs/=0.or.pawbec==1.or.pawpiezo==1.or.dtset%efmas>0) option=2
   call initylmg(gprimd,kg,dtset%kptns,dtset%mkmem,mpi_enreg,psps%mpsang,dtset%mpw,dtset%nband,dtset%nkpt,&
&   npwarr,dtset%nsppol,option,rprimd,ylm,ylmgr)
 end if

 call timab(133,2,tsec)
 call timab(134,1,tsec)

!Open and read pseudopotential files
 call status(0,dtfil%filstat,iexit,level,'call pspini(1)')
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,level,pawrad,pawtab,&
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
&   pawrhoij,dtset%pawspnorb,pawtab,dtset%spinat,dtset%typat,&
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

 call status(0,dtfil%filstat,iexit,level,'call inwffil(1')

!Initialize wavefunction files and wavefunctions.
 ireadwf0=1

 mcg=dtset%mpw*dtset%nspinor*dtset%mband*dtset%mkmem*dtset%nsppol
 ABI_STAT_ALLOCATE(cg,(2,mcg), ierr)
 ABI_CHECK(ierr==0, "out-of-memory in cg")

 ABI_ALLOCATE(eigen0,(dtset%mband*dtset%nkpt*dtset%nsppol))
 eigen0(:)=zero ; ask_accurate=1
 optorth=0

 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,gmet,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,rprimd,dtset%symafm,dtset%symrel,dtset%tnons,&
& dtfil%unkg,wffgs,wfftgs,dtfil%unwffgs,dtfil%fnamewffk,wvl)

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
     call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,&
&     rprimd,symrec,pawang%zarot)

     ! Update internal values
     call paw_gencond(Dtset,gnt_option,"save",call_pawinit)

     call timab(553,2,tsec)
   else
     if (pawtab(1)%has_kij  ==1) pawtab(1:psps%ntypat)%has_kij  =2
     if (pawtab(1)%has_nabla==1) pawtab(1:psps%ntypat)%has_nabla=2
   end if
   psps%n1xccc=maxval(pawtab(1:psps%ntypat)%usetcore)
   call setsymrhoij(gprimd,pawang%l_max-1,dtset%nsym,dtset%pawprtvol,rprimd,symrec,pawang%zarot)
   pawtab(:)%usepawu=0
   pawtab(:)%useexexch=0
   pawtab(:)%exchmix=zero
!  if (dtset%usepawu>0.or.dtset%useexexch>0) then
   call pawpuxinit(dtset%dmatpuopt,dtset%exchmix,dtset%f4of2_sla,dtset%f6of2_sla,&
&   dtset%jpawu,dtset%lexexch,dtset%lpawu,ntypat,pawang,dtset%pawprtvol,pawrad,&
&   pawtab,dtset%upawu,dtset%usedmft,dtset%useexexch,dtset%usepawu)
!  end if
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
   call status(0,dtfil%filstat,iexit,level,'call nhatgrid ')
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
   has_dijnd=0; req_cplex_dij=1
   if(any(abs(dtset%nucdipmom)>tol8)) then
     has_dijnd=1; req_cplex_dij=2
   end if
   if (rfphon/=0.or.rfelfd==1.or.rfelfd==3.or.rfstrs/=0.or.rf2_dkde/=0) then
     has_kxc=1;nkxc1=2*dtset%nspden-1 ! LDA only
     if(dtset%xclevel==2.and.dtset%pawxcdev==0) nkxc1=23
   end if
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,nkxc1,dtset%nspden,cplex,dtset%pawxcdev,&
&   dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_kxc=has_kxc,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%pawspnorb,&
&   natom,dtset%ntypat,dtset%typat,pawtab,has_dij=1,has_dijhartree=1,has_dijnd=has_dijnd,&
&   has_dijso=1,has_pawu_occ=1,has_exexch_pot=1,req_cplex_dij=req_cplex_dij,&
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
     nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
     call pawrhoij_alloc(pawrhoij_read,dtset%pawcpxocc,nspden_rhoij,dtset%nspinor,&
&     dtset%nsppol,dtset%typat,pawtab=pawtab)
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
   call fourdp(1,rhog,work,-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
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
     call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
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
&     ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&     psps%xcccrc,psps%xccc1d,xccc3d,xred)
   ABI_DEALLOCATE(dyfrx2)
   ABI_DEALLOCATE(vxc) ! dummy
 end if

!Set up hartree and xc potential. Compute kxc here.
 option=2 ; nk3xc=1
 nkxc=2*min(dtset%nspden,2)-1;if(dtset%xclevel==2)nkxc=23
 call check_kxc(dtset%ixc,dtset%optdriver)
 ABI_ALLOCATE(kxc,(nfftf,nkxc))
 ABI_ALLOCATE(vhartr,(nfftf))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))

 _IBM6("Before rhohxc")

 call rhohxc(dtset,enxc,gsqcut,psps%usepaw,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,dtset%nspden,n3xccc,option,rhog,rhor,&
& rprimd,strsxc,usexcnhat,vhartr,vxc,vxcavg,xccc3d)

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
&   dyfr_nondiag,efmasdeg,efmasfr,eigen0,eltfrnl,gsqcut,has_allddk,indsym,kg,mgfftf,&
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
&   psps%n1xccc,n3xccc,dtset%paral_kgb,psps,pawtab,ph1df,psps%qgrid_vl,&
&   dtset%qptn,rhog,rprimd,symq,symrec,dtset%typat,ucvol,&
&   psps%usepaw,psps%vlspl,vxc,psps%xcccrc,psps%xccc1d,xccc3d,xred)

   _IBM6("Before dfpt_ewald")

!  Compute Ewald (q=0) contribution
   sumg0=0;qphon(:)=zero
   call status(0,dtfil%filstat,iexit,level,'call dfpt_ewald(1)')
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
   call status(0,dtfil%filstat,iexit,level,'call symkchk ')
   timrev=1
   if (symkchk(dtset%kptns,dtset%nkpt,dtset%nsym,symrec,timrev,message) /= 0) then
     MSG_ERROR(message)
   end if

!  Calculate the kinetic part of the elastic tensor
   call dfpt_eltfrkin(cg,eltfrkin,dtset%ecut,dtset%ecutsm,dtset%effmass,&
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
     write(message, '(a,a,a,3es16.6,a,a,6(a,i2),a,a,a)' )ch10,&
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
   call status(0,dtfil%filstat,iexit,level,'call dfpt_ewald(2)')
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
&     psps%mqgrid_vl,natom,nfftf,ngfftf,nkxc,dtset%nspden,&
&     ntypat,psps%n1xccc,dtset%paral_kgb,psps,pawtab,ph1df,psps%qgrid_vl,qphon,&
&     rfdir,rfpert,rprimd,timrev,dtset%typat,ucvol,psps%usepaw,psps%xcccrc,psps%xccc1d,xred,rhor=rhor)
   else
     call dfpt_dyxc1(atindx,blkflgfrx1,dyfrx1,gmet,gsqcut,dtset%ixc,kxc,mgfftf,mpert,mpi_enreg,&
&     psps%mqgrid_vl,natom,nfftf,ngfftf,nkxc,dtset%nspden,&
&     ntypat,psps%n1xccc,dtset%paral_kgb,psps,pawtab,ph1df,psps%qgrid_vl,qphon,&
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
&   d2bbb,d2lo,d2nl,d2ovl,efmasdeg,efmasfr,eigbrd,eig2nkq,&
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
 if(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0 .or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and. rfuser==0)then

   write(message,'(a,a)' )ch10,' respfn : d/dk was computed, but no 2DTE, so no DDB output.'
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  If 2DTE were computed, only one processor must output them and compute
!  frequencies.
 else if(me==0)then

   write(message,'(a,a)' )ch10,' ==> Compute Derivative Database <== '
   call wrtout(std_out,message,'COLL')
   call wrtout(ab_out,message,'COLL')

!  Open the formatted derivative database file, and write the
!  preliminary information
   call status(0,dtfil%filstat,iexit,level,'call ddb_io_out')
   dscrpt=' Note : temporary (transfer) database '
!  tolwfr must be initialized here, but it is a dummy value
   tolwfr=1.0_dp
   call ddb_io_out (dscrpt,dtfil%fnameabo_ddb,natom,dtset%mband,&
&   dtset%nkpt,dtset%nsym,ntypat,dtfil%unddb,DDB_VERSION,&
&   dtset%acell_orig(1:3,1),dtset%amu_orig(:,1),dtset%dilatmx,dtset%ecut,dtset%ecutsm,&
&   dtset%intxc,iscf,dtset%ixc,dtset%kpt,dtset%kptnrm,&
&   natom,dtset%nband,ngfft,dtset%nkpt,dtset%nspden,dtset%nspinor,&
&   dtset%nsppol,dtset%nsym,ntypat,occ,dtset%occopt,dtset%pawecutdg,&
&   dtset%rprim_orig(1:3,1:3,1),dtset%dfpt_sciss,dtset%spinat,dtset%symafm,dtset%symrel,&
&   dtset%tnons,tolwfr,dtset%tphysel,dtset%tsmear,&
&   dtset%typat,dtset%usepaw,dtset%wtk,xred,psps%ziontypat,dtset%znucl)

   nblok=1 ; fullinit=1 ; choice=2
   call psddb8 (choice,psps%dimekb,psps%ekb,fullinit,psps%indlmn,&
&   psps%lmnmax,nblok,ntypat,dtfil%unddb,pawtab,&
&   psps%pspso,psps%usepaw,psps%useylm,DDB_VERSION)

!  In the RESPFN code, dfpt_nstdy and stady3 were called here
   d2nfr(:,:,:,:,:)=d2lo(:,:,:,:,:)+d2nl(:,:,:,:,:)
   if (psps%usepaw==1) d2nfr(:,:,:,:,:)=d2nfr(:,:,:,:,:)+d2ovl(:,:,:,:,:)

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
&       symrec,dtset%symrel,timrev)
       d2bbb(:,:,:,:,iband,iband) = d2tmp(:,:,natom+2,:,:)
     end do
     ABI_DEALLOCATE(blkflg1)
     ABI_DEALLOCATE(blkflg2)
     ABI_DEALLOCATE(d2tmp)
   end if

!  Complete the d2nfr matrix by symmetrization of the existing elements
   call d2sym3(blkflg,d2nfr,indsym,mpert,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev)

   if(rfphon==1.and.psps%n1xccc/=0)then
!    Complete the dyfrx1 matrix by symmetrization of the existing elements
     call d2sym3(blkflgfrx1,dyfrx1,indsym,natom,natom,dtset%nsym,qphon,symq,symrec,dtset%symrel,timrev)
   end if

!  Note that there is a bug in d2sym3 which will set some elements of
!  blkflg to 1 even when no corresponding symmetry-related element
!  has been computed.  This has the effect of producing spurious extra
!  output lines in the 2nd-order matrix listing in the .out file
!  and in the DDB file. The suprious matrix elements are all zero,
!  so this is primarily an annoyance.(DRH)


!  Add the frozen-wf (dyfrwf) part to the ewald part (dyew),
!  the part 1 of the frozen wf part of the xc core correction
!  (dyfrx1) and the non-frozen part (dynfr) to get the second-order
!  derivative matrix (d2matr), then
!  take account of the non-cartesian coordinates (d2cart).
   ABI_ALLOCATE(d2cart,(2,3,mpert,3,mpert))
   ABI_ALLOCATE(carflg,(3,mpert,3,mpert))
   ABI_ALLOCATE(d2matr,(2,3,mpert,3,mpert))
   outd2=1
   call status(0,dtfil%filstat,iexit,level,'call dfpt_gatherdy    ')
   call dfpt_gatherdy(becfrnl,dtset%berryopt,blkflg,carflg,&
&   dyew,dyfrwf,dyfrx1,dyfr_cplex,dyfr_nondiag,dyvdw,d2bbb,d2cart,d2cart_bbb,d2matr,d2nfr,&
&   eltcore,elteew,eltfrhar,eltfrkin,eltfrloc,eltfrnl,eltfrxc,eltvdw,&
&   gprimd,dtset%mband,mpert,natom,ntypat,outd2,pawbec,pawpiezo,piezofrnl,dtset%prtbbb,&
&   rfasr,rfpert,rprimd,dtset%typat,ucvol,usevdw,psps%ziontypat)
!  Output of the dynamical matrix
!  (Note : remember, previously, the processor me=0 has been selected)
   call status(0,dtfil%filstat,iexit,level,'call dfpt_dyout   ')
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
   call outddbnc(dtfil, dtset, hdr, psps, natom, mpert, rprimd, xred, dtset%qptn, d2matr, blkflg)
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
 if (.not.(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0.or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and. rfuser==0)) then
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
   if (.not.(rfphon==0 .and. (rf2_dkdk/=0 .or. rf2_dkde/=0 .or. rfddk/=0 .or. rfelfd==2) .and. rfstrs==0 .and.rfuser==0) )then
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
 ABI_DEALLOCATE(amass)
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
 call efmasfr_free_array(efmasfr)
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

 call status(0,dtfil%filstat,iexit,level,'exit')

 call timab(138,2,tsec)
 call timab(132,2,tsec)

 DBG_EXIT("COLL")

end subroutine respfn
!!***

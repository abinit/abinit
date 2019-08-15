!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_nonlinear
!! NAME
!!  m_nonlinear
!!
!! FUNCTION
!! DFT calculations of non linear response functions.
!!
!! COPYRIGHT
!!  Copyright (C) 2002-2019 ABINIT group (MVeithen,MB,LB)
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

module m_nonlinear

 use defs_basis
 use defs_datatypes
 use defs_wvltypes
 use m_wffile
 use m_errors
 use m_abicore
 use m_xmpi
 use m_hdr
 use m_ebands
 use m_xcdata
 use m_dtset
 use m_dtfil

 use defs_abitypes, only : MPI_type
 use m_fstrings, only : sjoin, itoa
 use m_time,     only : timab
 use m_symtk,    only : symmetrize_xred, littlegroup_q
 use m_dynmat,   only : d3sym, sytens
 use m_ddb,      only : nlopt, DDB_VERSION, dfptnl_doutput
 use m_ddb_hdr,  only : ddb_hdr_type, ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
 use m_ioarr,    only : read_rhor
 use m_kg,       only : getcut, kpgio, getph
 use m_fft,      only : fourdp
 use m_kpts,     only : getkgrid
 use m_inwffil,  only : inwffil
 use m_spacepar, only : hartre, setsym
 use m_pawfgr,      only : pawfgr_type,pawfgr_init, pawfgr_destroy
 use m_pawang,      only : pawang_type, pawang_init, pawang_free
 use m_pawrad,      only : pawrad_type
 use m_pawtab,      only : pawtab_type,pawtab_get_lsize
 use m_paw_an,      only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify, paw_an_print
 use m_paw_ij,      only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_print
 use m_pawfgrtab,   only : pawfgrtab_type, pawfgrtab_init, pawfgrtab_free
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_copy, &
&                          pawrhoij_bcast, pawrhoij_nullify, pawrhoij_inquire_dim
 use m_pawdij,      only : pawdij, symdij
 use m_paw_finegrid,only : pawexpiqr
 use m_pawxc,       only : pawxc_get_nkxc
 use m_paw_dmft,    only : paw_dmft_type
 use m_paw_sphharm, only : setsym_ylm
 use m_paw_nhat,    only : nhatgrid,pawmknhat
 use m_paw_denpot,  only : pawdenpot
 use m_paw_init,    only : pawinit,paw_gencond
 use m_paw_tools,   only : chkpawovlp
 use m_mkrho,       only : mkrho
 use m_getshell,    only : getshell
 use m_pspini,      only : pspini
 use m_atm2fft,     only : atm2fft
 use m_rhotoxc,     only : rhotoxc
 use m_mpinfo,      only : proc_distrb_cycle
 use m_mklocl,      only : mklocl
 use m_common,      only : setup1
 use m_fourier_interpol, only : transgrid
 use m_paw_occupancies,  only : initrhoij
 use m_paw_correlations, only : pawpuxinit
 use m_mkcore,           only : mkcore
 use m_pead_nl_loop,     only : pead_nl_loop
 use m_dfptnl_loop,      only : dfptnl_loop

 implicit none

 private
!!***

 public :: nonlinear
!!***

contains
!!***

!!****f* ABINIT/nonlinear
!! NAME
!! nonlinear
!!
!! FUNCTION
!! Primary routine for conducting DFT calculations of non linear response functions.
!!
!! INPUTS
!!  codvsn = code version
!!  dtfil <type(datafiles_type)> = variables related to files
!!  dtset <type(dataset_type)> = all input variables for this dataset
!!  etotal = new total energy (no meaning at output)
!!  mpi_enreg=informations about MPI pnarallelization
!!  occ(mband*nkpt*nsppol) = occupation number for each band and k
!!  xred(3,natom) = reduced atomic coordinates
!!
!! OUTPUT
!!
!!  npwtot(nkpt) = total number of plane waves at each k point
!!
!! SIDE EFFECTS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)> = variables related to pseudopotentials
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
!!      d3sym,dfptnl_doutput,dfptnl_loop,ebands_free,fourdp,getcut
!!      getkgrid,getshell,hdr_free,hdr_init,hdr_update,initmv,inwffil,kpgio
!!      mkcore,nlopt,pspini,read_rhor,rhotoxc,setsym,setup1,status
!!      ddb_hdr_init, ddb_hdr_free, ddb_hdr_open_write
!!      symmetrize_xred,sytens,timab,wffclose,wrtout
!!
!! SOURCE

subroutine nonlinear(codvsn,dtfil,dtset,etotal,mpi_enreg,npwtot,occ,pawang,pawrad,pawtab,psps,xred)

!Arguments ------------------------------------
!scalars
 real(dp),intent(inout) :: etotal
 character(len=6),intent(in) :: codvsn
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(pawang_type),intent(inout) :: pawang
 type(pseudopotential_type),intent(inout) :: psps
!arrays
 integer,intent(out) :: npwtot(dtset%nkpt)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),xred(3,dtset%natom)
 type(pawrad_type),intent(inout) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawtab_type),intent(inout) :: pawtab(psps%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 logical :: paral_atom,call_pawinit,qeq0
 integer,parameter :: level=50,formeig=0,response=1,cplex1=1
 integer :: ask_accurate,band_index,bantot,cplex,cplex_rhoij,dum_nshiftk,flag,gnt_option,gscase
 integer :: has_dijnd,has_diju,has_kxc,has_k3xc
 integer :: i1dir,i1pert,i2dir,i2pert,i3dir,i3pert
 integer :: iatom,indx,iband,ider,idir,ierr,ifft,ikpt,ipert,isppol
 integer :: ireadwf0,iscf_eff,ispden,itypat,izero,mcg,me,mgfftf,mkmem_max,mpert,my_natom
 integer :: n1,n3xccc,natom,nband_k,nfftf,nfftot,nfftotf,nhatdim,nhatgrdim
 integer :: nkpt_eff,nkpt_max,nkpt3,nkxc,nkxc1,nk3xc,nk3xc1,nneigh,ntypat,nsym1,nspden_rhoij,nzlmopt
 integer :: optcut,optgr0,optgr1,optgr2,optrad,option,optorth
 integer :: optatm,optdyfr,opteltfr,optgr,optstr,optv,optn,optn2
 integer :: psp_gencond,pead,qphase_rhoij,rdwr,rdwrpaw,spaceworld,tim_mkrho,timrev
 integer :: use_sym,usecprj,usexcnhat
 logical :: is_dfpt=.true.,nmxc
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: boxcut,compch_fft,compch_sph,ecore,ecut_eff,ecutdg_eff,ecutf
 real(dp) :: eei,epaw,epawdc,enxc,etot,fermie
 real(dp) :: gsqcut,gsqcut_eff,gsqcutc_eff
 real(dp) :: rdum,residm,ucvol,vxcavg
 character(len=500) :: message
 character(len=30) :: small_msg
 character(len=fnlen) :: dscrpt
 type(pawang_type) :: pawang1
 type(ebands_t) :: bstruct
 type(hdr_type) :: hdr,hdr_den
 type(ddb_hdr_type) :: ddb_hdr
 type(wffile_type) :: wffgs,wfftgs
 type(wvl_data) :: wvl
 type(xcdata_type) :: xcdata
!arrays
 integer :: dum_kptrlatt(3,3),dum_vacuum(3),ngfft(18),ngfftf(18),perm(6),ii,theunit
 integer,allocatable :: atindx(:),atindx1(:),blkflg(:,:,:,:,:,:),carflg(:,:,:,:,:,:),cgindex(:,:)
 integer,allocatable :: blkflg_tmp(:,:,:,:,:,:),blkflg_sav(:,:,:,:,:,:)
 integer,allocatable :: carflg_tmp(:,:,:,:,:,:),carflg_sav(:,:,:,:,:,:)
 integer,allocatable :: d3e_pert1(:),d3e_pert2(:),d3e_pert3(:)
 integer,allocatable :: indsym(:,:,:),indsy1(:,:,:),irrzon(:,:,:),irrzon1(:,:,:)
 integer,allocatable :: kg(:,:),kneigh(:,:),kg_neigh(:,:,:)
 integer,allocatable :: kptindex(:,:),l_size_atm(:)
 integer,allocatable :: npwarr(:),nattyp(:),pwind(:,:,:),rfpert(:,:,:,:,:,:)
 integer,allocatable :: symq(:,:,:),symrec(:,:,:),symaf1(:),symrc1(:,:,:),symrl1(:,:,:)
 real(dp) :: dum_gauss(0),dum_dyfrn(0),dum_dyfrv(0),dum_eltfrxc(0)
 real(dp) :: dum_grn(0),dum_grv(0),dum_rhog(0),dum_vg(0)
 real(dp) :: dum_shiftk(3,MAX_NSHIFTK),dummy6(6),other_dummy6(6),gmet(3,3),gprimd(3,3)
 real(dp) :: qphon(3),rmet(3,3),rprimd(3,3),strsxc(6),tsec(2)
 real(dp),allocatable :: cg(:,:),d3cart(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3etot(:,:,:,:,:,:,:),dum_kptns(:,:)
! We need all these arrays instead of one because in Fortran the maximum number of dimensions is 7...
 real(dp),allocatable :: d3e_1(:,:,:,:,:,:,:),d3cart_1(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_2(:,:,:,:,:,:,:),d3cart_2(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_3(:,:,:,:,:,:,:),d3cart_3(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_4(:,:,:,:,:,:,:),d3cart_4(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_5(:,:,:,:,:,:,:),d3cart_5(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_6(:,:,:,:,:,:,:),d3cart_6(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_7(:,:,:,:,:,:,:),d3cart_7(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_8(:,:,:,:,:,:,:),d3cart_8(:,:,:,:,:,:,:)
 real(dp),allocatable :: d3e_9(:,:,:,:,:,:,:),d3cart_9(:,:,:,:,:,:,:)
 real(dp),allocatable :: dum_wtk(:),dyfrlo_indx(:,:,:),dyfrx2(:,:,:),eigen0(:)
 real(dp),allocatable :: grtn_indx(:,:),grxc(:,:),k3xc(:,:),kpt3(:,:),kxc(:,:)
 real(dp),allocatable :: mvwtk(:,:),nhat(:,:),nhatgr(:,:,:),ph1d(:,:),ph1df(:,:),phnons(:,:,:),phnons1(:,:,:)
 real(dp),allocatable :: rhog(:,:),rhor(:,:),rhowfg(:,:),rhowfr(:,:),tnons1(:,:)
 real(dp),allocatable :: vhartr(:),vpsp(:),vtrial(:,:),vxc(:,:),work(:),xccc3d(:)
 type(pawfgr_type) :: pawfgr
 type(pawrhoij_type),allocatable :: pawrhoij(:),pawrhoij_read(:)
 type(pawfgrtab_type),allocatable,save :: pawfgrtab(:)
 type(paw_an_type),allocatable :: paw_an(:)
 type(paw_ij_type),allocatable :: paw_ij(:)
 type(paw_dmft_type) :: paw_dmft

! ***********************************************************************

 DBG_ENTER("COLL")

 call timab(501,1,tsec)

!Structured debugging if dtset%prtvol==-level
 if(dtset%prtvol==-level)then
   write(message,'(80a,a,a)')  ('=',ii=1,80),ch10,' nonlinear : enter , debug mode '
   call wrtout(std_out,message,'COLL')
 end if

!Check if the perturbations asked in the input file can be computed

 if (((dtset%d3e_pert1_phon == 1).and.(dtset%d3e_pert2_phon == 1)).or. &
& ((dtset%d3e_pert1_phon == 1).and.(dtset%d3e_pert3_phon == 1)).or. &
& ((dtset%d3e_pert2_phon == 1).and.(dtset%d3e_pert3_phon == 1))) then
   write(message,'(7a)')&
&   'You have asked for a third-order derivative with respect to',ch10,&
&   '2 or more atomic displacements.',ch10,&
&   'This is not allowed yet.',ch10,&
&   'Action : change d3e_pert1_phon, d3e_pert2_phon or d3e_pert3_phon in your input file.'
   MSG_ERROR(message)
 end if

!Computation of third order derivatives from PEAD (pead=1) or full DPFT formalism (pead=0):
 pead = dtset%usepead
 if (pead==0) then
   write(message, '(2a)' ) ch10,'NONLINEAR : PEAD=0, full DFPT computation of third order derivatives'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,message,'COLL')
 end if

!Some data for parallelism
 nkpt_max=50;if(xmpi_paral==1)nkpt_max=-1
 my_natom=mpi_enreg%my_natom
 paral_atom=(my_natom/=dtset%natom)
 if (paral_atom) then
   MSG_BUG(" Nonlinear routine is not available yet with parallelization over atoms...")
 end if

!Init spaceworld
 spaceworld=mpi_enreg%comm_cell
 me = xmpi_comm_rank(spaceworld)

!Define FFT grid(s) sizes (be careful !)
!See NOTES in the comments at the beginning of this file.
 call pawfgr_init(pawfgr,dtset,mgfftf,nfftf,ecut_eff,ecutdg_eff,ngfft,ngfftf)

 ntypat=psps%ntypat
 natom=dtset%natom
 nfftot=product(ngfft(1:3))
 nfftotf=product(ngfftf(1:3))

!Define the set of admitted perturbations taking into account
!the possible permutations
 mpert=natom+6
 ABI_ALLOCATE(blkflg,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(carflg,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(rfpert,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_pert1,(mpert))
 ABI_ALLOCATE(d3e_pert2,(mpert))
 ABI_ALLOCATE(d3e_pert3,(mpert))
 ABI_ALLOCATE(d3etot,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(blkflg_tmp,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(blkflg_sav,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(carflg_tmp,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(carflg_sav,(3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_1,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_2,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_3,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_4,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_5,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_6,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_7,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_8,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3e_9,(2,3,mpert,3,mpert,3,mpert))
 d3e_1(:,:,:,:,:,:,:) = 0_dp
 d3e_2(:,:,:,:,:,:,:) = 0_dp
 d3e_3(:,:,:,:,:,:,:) = 0_dp
 d3e_4(:,:,:,:,:,:,:) = 0_dp
 d3e_5(:,:,:,:,:,:,:) = 0_dp
 d3e_6(:,:,:,:,:,:,:) = 0_dp
 d3e_7(:,:,:,:,:,:,:) = 0_dp
 d3e_8(:,:,:,:,:,:,:) = 0_dp
 d3e_9(:,:,:,:,:,:,:) = 0_dp
 ABI_ALLOCATE(d3cart_1,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_2,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_3,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_4,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_5,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_6,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_7,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_8,(2,3,mpert,3,mpert,3,mpert))
 ABI_ALLOCATE(d3cart_9,(2,3,mpert,3,mpert,3,mpert))
 blkflg(:,:,:,:,:,:) = 0
 d3etot(:,:,:,:,:,:,:) = 0_dp
 rfpert(:,:,:,:,:,:) = 0
 d3e_pert1(:) = 0 ; d3e_pert2(:) = 0 ; d3e_pert3(:) = 0

 if (dtset%d3e_pert1_phon==1) d3e_pert1(dtset%d3e_pert1_atpol(1):dtset%d3e_pert1_atpol(2))=1
 if (dtset%d3e_pert2_phon==1) d3e_pert2(dtset%d3e_pert2_atpol(1):dtset%d3e_pert2_atpol(2))=1
 if (dtset%d3e_pert3_phon==1) d3e_pert3(dtset%d3e_pert3_atpol(1):dtset%d3e_pert3_atpol(2))=1
 if (dtset%d3e_pert1_elfd/=0) d3e_pert1(natom+2)=1
 if (dtset%d3e_pert2_elfd/=0) d3e_pert2(natom+2)=1
 if (dtset%d3e_pert3_elfd/=0) d3e_pert3(natom+2)=1

 do i1pert = 1, mpert
   do i1dir = 1, 3
     do i2pert = 1, mpert
       do i2dir = 1, 3
         do i3pert = 1, mpert
           do i3dir = 1, 3
             perm(1) = &
&              d3e_pert1(i1pert)*dtset%d3e_pert1_dir(i1dir) &
&             *d3e_pert2(i2pert)*dtset%d3e_pert2_dir(i2dir) &
&             *d3e_pert3(i3pert)*dtset%d3e_pert3_dir(i3dir)
             perm(2) = &
&              d3e_pert1(i1pert)*dtset%d3e_pert1_dir(i1dir) &
&             *d3e_pert2(i3pert)*dtset%d3e_pert2_dir(i3dir) &
&             *d3e_pert3(i2pert)*dtset%d3e_pert3_dir(i2dir)
             perm(3) = &
&              d3e_pert1(i2pert)*dtset%d3e_pert1_dir(i2dir) &
&             *d3e_pert2(i1pert)*dtset%d3e_pert2_dir(i1dir) &
&             *d3e_pert3(i3pert)*dtset%d3e_pert3_dir(i3dir)
             perm(4) = &
&              d3e_pert1(i2pert)*dtset%d3e_pert1_dir(i2dir) &
&             *d3e_pert2(i3pert)*dtset%d3e_pert2_dir(i3dir) &
&             *d3e_pert3(i1pert)*dtset%d3e_pert3_dir(i1dir)
             perm(5) = &
&              d3e_pert1(i3pert)*dtset%d3e_pert1_dir(i3dir) &
&             *d3e_pert2(i2pert)*dtset%d3e_pert2_dir(i2dir) &
&             *d3e_pert3(i1pert)*dtset%d3e_pert3_dir(i1dir)
             perm(6) = &
&              d3e_pert1(i3pert)*dtset%d3e_pert1_dir(i3dir) &
&             *d3e_pert2(i1pert)*dtset%d3e_pert2_dir(i1dir) &
&             *d3e_pert3(i2pert)*dtset%d3e_pert3_dir(i2dir)
             if (sum(perm(:)) > 0) rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
           end do
         end do
       end do
     end do
   end do
 end do

! call timab(134,2,tsec)
! call timab(135,1,tsec)

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

 call sytens(indsym,mpert,natom,dtset%nsym,rfpert,symrec,dtset%symrel)

 write(message, '(a,a,a,a,a)' ) ch10, &
& ' The list of irreducible elements of the Raman and non-linear',&
& ch10,' optical susceptibility tensors is:',ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

 write(message,'(12x,a)')&
& 'i1pert  i1dir   i2pert  i2dir   i3pert  i3dir'
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')
 n1 = 0
 do i1pert = 1, natom + 2
   do i1dir = 1, 3
     do i2pert = 1, natom + 2
       do i2dir = 1,3
         do i3pert = 1, natom + 2
           do i3dir = 1, 3
             if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
               n1 = n1 + 1
               write(message,'(2x,i4,a,6(5x,i3))') n1,')', &
&               i1pert,i1dir,i2pert,i2dir,i3pert,i3dir
               call wrtout(ab_out,message,'COLL')
               call wrtout(std_out,message,'COLL')
             else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-2) then
               blkflg(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert) = 1
               if (dtset%nonlinear_info>0) then
!                 n1 = n1 + 1
                 write(message,'(2x,i4,a,6(5x,i3),a)') n1,')', &
  &               i1pert,i1dir,i2pert,i2dir,i3pert,i3dir,' => must be zero, not computed'
                 call wrtout(ab_out,message,'COLL')
                 call wrtout(std_out,message,'COLL')
               end if
             else if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==-1) then
               if (dtset%nonlinear_info>0) then
!                 n1 = n1 + 1
                 write(message,'(2x,i4,a,6(5x,i3),a)') n1,')', &
  &               i1pert,i1dir,i2pert,i2dir,i3pert,i3dir,' => symmetric of an other element, not computed'
                 call wrtout(ab_out,message,'COLL')
                 call wrtout(std_out,message,'COLL')
               end if
             end if
           end do
         end do
       end do
     end do
   end do
 end do
 write(message,'(a,a)') ch10,ch10
 call wrtout(ab_out,message,'COLL')
 call wrtout(std_out,message,'COLL')

! For abipy :
 if (dtset%paral_rf == -1) then
   write(std_out,'(a)')"--- !IrredPerts"
   write(std_out,'(a)')'# List of irreducible perturbations for nonlinear'
   write(std_out,'(a)')'irred_perts:'

   n1 = 0
   do i1pert = 1, natom + 2
     do i1dir = 1, 3
       do i2pert = 1, natom + 2
         do i2dir = 1, 3
           do i3pert = 1, natom + 2
             do i3dir = 1,3
               if (rfpert(i1dir,i1pert,i2dir,i2pert,i3dir,i3pert)==1) then
                 n1 = n1 + 1
                 write(std_out,'(a,i0)')"   - i1pert: ",i1pert
                 write(std_out,'(a,i0)')"     i1dir: ",i1dir
                 write(std_out,'(a,i0)')"     i2pert: ",i2pert
                 write(std_out,'(a,i0)')"     i2dir: ",i2dir
                 write(std_out,'(a,i0)')"     i3pert: ",i3pert
                 write(std_out,'(a,i0)')"     i3dir: ",i3dir
               end if
             end do
           end do
         end do
       end do
     end do
   end do
   write(std_out,'(a)')"..."
   MSG_ERROR_NODUMP("aborting now")
 end if

!Set up for iterations
 call setup1(dtset%acell_orig(1:3,1),bantot,dtset,&
  ecutdg_eff,ecut_eff,gmet,gprimd,gsqcut_eff,gsqcutc_eff,&
   ngfftf,ngfft,dtset%nkpt,dtset%nsppol,&
   response,rmet,dtset%rprim_orig(1:3,1:3,1),rprimd,ucvol,psps%usepaw)

!Set up the basis sphere of planewaves
 ABI_ALLOCATE(kg,(3,dtset%mpw*dtset%mkmem))
 ABI_ALLOCATE(npwarr,(dtset%nkpt))
 call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg,&
& dtset%kptns,dtset%mkmem,dtset%nband,dtset%nkpt,'PERS',mpi_enreg,dtset%mpw,npwarr,npwtot,&
& dtset%nsppol)

!Recompute first large sphere cut-off gsqcut, without taking into account dilatmx
 ecutf=dtset%ecut
 if (psps%usepaw==1) then
   ecutf=dtset%pawecutdg
   call wrtout(std_out,ch10//' FFT (fine) grid used in SCF cycle:','COLL')
 end if
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,k0,ngfftf)

!Open and read pseudopotential files
 ecore = 0_dp
 call pspini(dtset,dtfil,ecore,psp_gencond,gsqcutc_eff,gsqcut_eff,pawrad,pawtab,&
& psps,rprimd,comm_mpi=mpi_enreg%comm_cell)

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

 call inwffil(ask_accurate,cg,dtset,dtset%ecut,ecut_eff,eigen0,dtset%exchn2n3d,&
& formeig,hdr,ireadwf0,dtset%istwfk,kg,dtset%kptns,&
& dtset%localrdwf,dtset%mband,mcg,dtset%mkmem,mpi_enreg,dtset%mpw,&
& dtset%nband,ngfft,dtset%nkpt,npwarr,dtset%nsppol,dtset%nsym,&
& occ,optorth,dtset%symafm,dtset%symrel,dtset%tnons,&
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

! call timab(135,2,tsec)
! call timab(136,1,tsec)

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

!Compute structure factor phases for current atomic pos:
 ABI_ALLOCATE(ph1d,(2,3*(2*dtset%mgfft+1)*natom))
 ABI_ALLOCATE(ph1df,(2,3*(2*mgfftf+1)*natom))
 call getph(atindx,natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)

 if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
   call getph(atindx,natom,ngfftf(1),ngfftf(2),ngfftf(3),ph1df,xred)
 else
   ph1df(:,:)=ph1d(:,:)
 end if

qeq0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2<1.d-14)
if (.not.qeq0) then
  MSG_BUG('NONLINEAR with dtset%qptn!=0 is not implemented yet')
end if

!PAW: 1- Initialize values for several arrays depending only on atomic data
!2- Check overlap
!3- Identify FFT points in spheres and compute g_l(r).Y_lm(r) (and exp(-i.q.r) if needed)
!4- Allocate PAW specific arrays
!5- Compute perturbed local potential inside spheres
!6- Eventually open temporary storage files
 if(psps%usepaw==1) then
!  1-Initialize values for several arrays depending only on atomic data

   gnt_option=2

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

!  Note: many derivatives of cprj are needed and used a few times only, so for simplicity the
!  computation of all needed derivatives will be done on-the-fly.
   usecprj=0

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
   optcut=0;optgr0=dtset%pawstgylm;optgr1=0;optgr2=0;optrad=1-dtset%pawstgylm
   optgr1=dtset%pawstgylm
   optgr2=dtset%pawstgylm
   call nhatgrid(atindx1,gmet,my_natom,natom,nattyp,ngfftf,psps%ntypat,&
&   optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred,&
&   comm_atom=mpi_enreg%comm_atom, mpi_atmtab=mpi_enreg%my_atmtab )
   ABI_DATATYPE_ALLOCATE(paw_an,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_ij,(my_natom))
   call paw_an_nullify(paw_an)
   call paw_ij_nullify(paw_ij)
   has_kxc=0;nkxc1=0;cplex=1
   has_dijnd=0;if(any(abs(dtset%nucdipmom)>tol8)) has_dijnd=1
   has_diju=merge(0,1,dtset%usepawu==0)
   has_kxc=1;nkxc1=2*dtset%nspden-1 ! LDA only
   call pawxc_get_nkxc(nkxc1,dtset%nspden,dtset%xclevel)
   has_k3xc=1; nk3xc1=3*min(dtset%nspden,2)-2 ! LDA only
   call paw_an_init(paw_an,dtset%natom,dtset%ntypat,nkxc1,nk3xc1,dtset%nspden,cplex,dtset%pawxcdev,&
&   dtset%typat,pawang,pawtab,has_vxc=1,has_vxc_ex=1,has_kxc=has_kxc,has_k3xc=has_k3xc,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   call paw_ij_init(paw_ij,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,dtset%pawspnorb,&
&   natom,dtset%ntypat,dtset%typat,pawtab,has_dij=1,has_dijhartree=1,has_dijnd=has_dijnd,&
&   has_dijso=1,has_dijU=has_diju,has_pawu_occ=1,has_exexch_pot=1,nucdipmom=dtset%nucdipmom,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
 else ! PAW vs NCPP
   usexcnhat=0;usecprj=0
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
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                          nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoij_read,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&                        dtset%nsppol,dtset%typat,qphase=qphase_rhoij,pawtab=pawtab)
   else
     ABI_DATATYPE_ALLOCATE(pawrhoij_read,(0))
   end if

!    MT july 2013: Should we read rhoij from the density file ?
   call read_rhor(dtfil%fildensin, cplex1, dtset%nspden, nfftf, ngfftf, rdwrpaw, mpi_enreg, rhor, &
   hdr_den, pawrhoij_read, spaceworld, check_hdr=hdr)
   call hdr_free(hdr_den)

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

!!jmb 2012 write(std_out,'(a)')' ' ! needed to make ibm6_xlf12 pass tests. No idea why this works. JWZ 5 Sept 2011
!!Will compute now the total potential

!Compute local ionic pseudopotential vpsp and core electron density xccc3d:
 n3xccc=0;if (psps%n1xccc/=0) n3xccc=nfftf
 ABI_ALLOCATE(xccc3d,(n3xccc))
 ABI_ALLOCATE(vpsp,(nfftf))

 eei = zero
 if (psps%usepaw==1 .or. psps%nc_xccc_gspace==1) then
!  PAW or NC with nc_xccc_gspace: compute Vloc and core charge together in reciprocal space
   call timab(562,1,tsec)
   optatm=1;optdyfr=0;opteltfr=0;optgr=0;optstr=0;optv=1;optn=n3xccc/nfftf;optn2=1
   call atm2fft(atindx1,xccc3d,vpsp,dum_dyfrn,dum_dyfrv,dum_eltfrxc,dum_gauss,gmet,gprimd,&
&   dum_grn,dum_grv,gsqcut,mgfftf,psps%mqgrid_vl,natom,nattyp,nfftf,ngfftf,&
&   ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1df,psps%qgrid_vl,&
&   dtset%qprtrb,dum_rhog,dummy6,other_dummy6,ucvol,psps%usepaw,dum_vg,dum_vg,dum_vg,dtset%vprtrb,psps%vlspl)
   call timab(562,2,tsec)
 else
!  Norm-cons.: compute Vloc in reciprocal space and core charge in real space
   option=1
   ABI_ALLOCATE(dyfrlo_indx,(3,3,natom))
   ABI_ALLOCATE(grtn_indx,(3,natom))
   call mklocl(dtset,dyfrlo_indx,eei,gmet,gprimd,&
&   grtn_indx,gsqcut,dummy6,mgfftf,mpi_enreg,natom,nattyp,&
&   nfftf,ngfftf,dtset%nspden,ntypat,option,pawtab,ph1df,psps,&
&   dtset%qprtrb,rhog,rhor,rprimd,ucvol,dtset%vprtrb,vpsp,wvl%descr,wvl%den,xred)
   ABI_DEALLOCATE(dyfrlo_indx)
   ABI_DEALLOCATE(grtn_indx)
   if (psps%n1xccc/=0) then
     ABI_ALLOCATE(dyfrx2,(3,3,natom))
     ABI_ALLOCATE(vxc,(0,0)) ! dummy
     call mkcore(dummy6,dyfrx2,grxc,mpi_enreg,natom,nfftf,dtset%nspden,ntypat,&
&     ngfftf(1),psps%n1xccc,ngfftf(2),ngfftf(3),option,rprimd,dtset%typat,ucvol,vxc,&
&     psps%xcccrc,psps%xccc1d,xccc3d,xred)
     ABI_DEALLOCATE(dyfrx2)
     ABI_DEALLOCATE(vxc) ! dummy
   end if
 end if

!Set up hartree and xc potential. Compute kxc here.
 ABI_ALLOCATE(vhartr,(nfftf))
 call hartre(1,gsqcut,psps%usepaw,mpi_enreg,nfftf,ngfftf,rhog,rprimd,vhartr)

 option=3
 nkxc=2*dtset%nspden-1 ! LDA
 if(dtset%xclevel==2.and.dtset%nspden==1) nkxc=7  ! non-polarized GGA
 if(dtset%xclevel==2.and.dtset%nspden==2) nkxc=19 ! polarized GGA
 nk3xc=3*dtset%nspden-2
 ABI_ALLOCATE(kxc,(nfftf,nkxc))
 ABI_ALLOCATE(k3xc,(nfftf,nk3xc))
 ABI_ALLOCATE(vxc,(nfftf,dtset%nspden))

 _IBM6("Before rhotoxc")

 call xcdata_init(xcdata,dtset=dtset)
 nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
 call rhotoxc(enxc,kxc,mpi_enreg,nfftf,ngfftf,&
& nhat,nhatdim,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc,option,rhor,&
& rprimd,strsxc,usexcnhat,vxc,vxcavg,xccc3d,xcdata,k3xc=k3xc,vhartr=vhartr)

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
 ABI_DEALLOCATE(vpsp)
 ABI_DEALLOCATE(vhartr)

 if(dtset%prtvol==-level)then
   call wrtout(std_out,' nonlinear : ground-state density and potential set up.','COLL')
 end if

 epaw = zero ; epawdc = zero
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

 ABI_DEALLOCATE(xccc3d)

!  Determine the subset of symmetry operations (nsym1 operations)
!  that leaves the perturbation invariant, and initialize corresponding arrays
!  symaf1, symrl1, tnons1 (and pawang1%zarot, if PAW)..
 nsym1 = 1
! symaf1_tmp(1) = 1
! symrl1_tmp(:,:,1) = dtset%symrel(:,:,1)
! tnons1_tmp(:,1) = 0_dp
   ABI_ALLOCATE(indsy1,(4,nsym1,dtset%natom))
   ABI_ALLOCATE(symrc1,(3,3,nsym1))
   ABI_ALLOCATE(symaf1,(nsym1))
   ABI_ALLOCATE(symrl1,(3,3,nsym1))
   ABI_ALLOCATE(tnons1,(3,nsym1))
   symaf1(1)= 1 !symaf1_tmp(1:nsym1)
   symrl1(:,:,1)= dtset%symrel(:,:,1) !symrl1_tmp(:,:,1:nsym1)
   tnons1(:,1)= 0_dp !tnons1_tmp(:,1:nsym1)
!   ABI_DEALLOCATE(symaf1_tmp)
!   ABI_DEALLOCATE(symrl1_tmp)
!   ABI_DEALLOCATE(tnons1_tmp)

!  Set up corresponding symmetry data
 ABI_ALLOCATE(irrzon1,(dtset%nfft**(1-1/nsym1),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 ABI_ALLOCATE(phnons1,(2,dtset%nfft**(1-1/nsym1),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4)))
 call setsym(indsy1,irrzon1,1,dtset%natom,dtset%nfft,dtset%ngfft,dtset%nspden,dtset%nsppol,&
& nsym1,phnons1,symaf1,symrc1,symrl1,tnons1,dtset%typat,xred)
 if (psps%usepaw==1) then
!  Allocate/initialize only zarot in pawang1 datastructure
   call pawang_init(pawang1,0,pawang%l_max-1,0,nsym1,0,1,0,0,0)
   call setsym_ylm(gprimd,pawang1%l_max-1,pawang1%nsym,0,rprimd,symrc1,pawang1%zarot)
 end if

 if (pead/=0) then
!  Initialize finite difference calculation of the ddk

   nkpt3 = 0

!  Prepare first call to getkgrid (obtain number of k points in FBZ)
   dum_kptrlatt(:,:) = dtset%kptrlatt(:,:)
   dum_nshiftk = dtset%nshiftk
   ABI_CHECK(dum_nshiftk <= MAX_NSHIFTK, sjoin("dum_nshiftk must be <= ", itoa(MAX_NSHIFTK)))
   dum_shiftk(:,:) = zero
   dum_shiftk(:,1:dtset%nshiftk) = dtset%shiftk(:,1:dtset%nshiftk)
   dum_vacuum(:) = 0

   ABI_ALLOCATE(dum_kptns,(3,0))
   ABI_ALLOCATE(dum_wtk,(0))
   call getkgrid(0,0,dtset%iscf,dum_kptns,3,dum_kptrlatt,&
&   rdum,dtset%nsym,0,nkpt3,dum_nshiftk,dtset%nsym,&
&   rprimd,dum_shiftk,dtset%symafm,dtset%symrel,&
&   dum_vacuum,dum_wtk)
   ABI_DEALLOCATE(dum_kptns)
   ABI_DEALLOCATE(dum_wtk)

!   write(std_out,*) 'nonlinear : nkpt, nkpt3 = ',dtset%nkpt,nkpt3
!call flush(6)
!jmb : malloc() problem with gcc461_openmpi under max2 : change order of allocations works ?!?
!allocate(kneigh(30,nkpt),kg_neigh(30,nkpt,3),mvwtk(30,nkpt))
   ABI_ALLOCATE(kg_neigh,(30,dtset%nkpt,3))
   ABI_ALLOCATE(mvwtk,(30,dtset%nkpt))
   ABI_ALLOCATE(kneigh,(30,dtset%nkpt))

   ABI_ALLOCATE(kptindex,(2,nkpt3))
   ABI_ALLOCATE(kpt3,(3,nkpt3))

   call getshell(gmet,kneigh,kg_neigh,kptindex,dtset%kptopt,&
&   dtset%kptrlatt,dtset%kptns,kpt3,dtset%mkmem,mkmem_max,mvwtk,&
&   dtset%nkpt,nkpt3,nneigh,dtset%nshiftk,rmet,rprimd,dtset%shiftk,dtset%wtk, mpi_enreg%comm_cell)

   ABI_ALLOCATE(pwind,(dtset%mpw,nneigh,dtset%mkmem))
   ABI_ALLOCATE(cgindex,(dtset%nkpt,dtset%nsppol))
   ABI_ALLOCATE(mpi_enreg%kpt_loc2ibz_sp,(0:mpi_enreg%nproc-1,1:mkmem_max, 1:2))
   ABI_ALLOCATE(mpi_enreg%mkmem,(0:mpi_enreg%nproc-1))

   call initmv(cgindex,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&
&   kpt3,dtset%mband,dtset%mkmem,mpi_enreg,dtset%mpw,dtset%nband,dtset%nkpt,&
&   nkpt3,nneigh,npwarr,dtset%nsppol,occ,pwind)

  call pead_nl_loop(blkflg,cg,cgindex,dtfil,dtset,d3etot,gmet,gprimd,gsqcut,&
&  hdr,kg,kneigh,kg_neigh,kptindex,kpt3,kxc,k3xc,dtset%mband,dtset%mgfft,&
&  dtset%mkmem,mkmem_max,dtset%mk1mem,mpert,mpi_enreg,dtset%mpw,mvwtk,natom,nfftf,&
&  dtset%nkpt,nkpt3,nkxc,nk3xc,nneigh,dtset%nspinor,dtset%nsppol,npwarr,occ,psps,pwind,&
&  rfpert,rprimd,ucvol,xred)

 else ! pead=0 in this case

   call dfptnl_loop(atindx,blkflg,cg,dtfil,dtset,d3etot,eigen0,gmet,gprimd,gsqcut,&
&   hdr,kg,kxc,k3xc,dtset%mband,dtset%mgfft,mgfftf,&
&   dtset%mkmem,dtset%mk1mem,mpert,mpi_enreg,dtset%mpw,natom,nattyp,ngfftf,nfftf,nhat,&
&   dtset%nkpt,nkxc,nk3xc,dtset%nspinor,dtset%nsppol,npwarr,occ,&
&   paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawtab,&
&   ph1d,ph1df,psps,rfpert,rhog,rhor,rprimd,ucvol,usecprj,vtrial,vxc,xred,&
&   nsym1,indsy1,symaf1,symrc1,&
&   d3e_1,d3e_2,d3e_3,d3e_4,d3e_5,d3e_6,d3e_7,d3e_8,d3e_9)

 end if ! end pead/=0

 write(message,'(a,a,a)')ch10,&
& ' --- Third order energy calculation completed --- ',ch10
 call wrtout(ab_out,message,'COLL')

!Complete missing elements using symmetry operations
 blkflg_sav = blkflg
 call d3sym(blkflg,d3etot,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)

 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_1,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_2,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_3,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_4,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_5,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_6,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_7,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_8,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)
 blkflg_tmp = blkflg_sav
 call d3sym(blkflg_tmp,d3e_9,indsym,mpert,natom,dtset%nsym,symrec,dtset%symrel)

!Open the formatted derivative database file, and write the
!preliminary information
 if (mpi_enreg%me == 0) then
   dscrpt=' Note : temporary (transfer) database '

   call ddb_hdr_init(ddb_hdr,dtset,psps,pawtab,DDB_VERSION,dscrpt,1,xred=xred,occ=occ)

   call ddb_hdr_open_write(ddb_hdr, dtfil%fnameabo_ddb, dtfil%unddb)

   call ddb_hdr_free(ddb_hdr)

!  Call main output routine
   call dfptnl_doutput(blkflg,d3etot,dtset%mband,mpert,dtset%nkpt,dtset%natom,dtset%ntypat,dtfil%unddb)

!  Close DDB
   close(dtfil%unddb)

!  Compute tensors related to third-order derivatives
   call nlopt(blkflg,carflg,d3etot,d3cart,gprimd,mpert,natom,rprimd,ucvol)
!  Note that the imaginary part is not transformed into cartesian coordinates
   call nlopt(blkflg,carflg_tmp,d3e_1,d3cart_1,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_2,d3cart_2,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_3,d3cart_3,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_4,d3cart_4,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_5,d3cart_5,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_6,d3cart_6,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_7,d3cart_7,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_8,d3cart_8,gprimd,mpert,natom,rprimd,ucvol)
   call nlopt(blkflg,carflg_tmp,d3e_9,d3cart_9,gprimd,mpert,natom,rprimd,ucvol)

   if ((d3e_pert1(natom+2)==1).and.(d3e_pert2(natom+2)==1).and. &
&   (d3e_pert3(natom+2)==1)) then

     flag = 1
     i1pert = natom+2

     d3cart(:,:,i1pert,:,i1pert,:,i1pert) = &
&     d3cart(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

     write(ab_out,*)ch10
     write(ab_out,*)' Non-linear optical susceptibility tensor d (pm/V)'
     write(ab_out,*)' in cartesian coordinates'
     write(ab_out,*)'  i1dir  i2dir  i3dir             d'

     do i1dir = 1, 3
       do i2dir = 1, 3
         do i3dir = 1, 3
           write(ab_out,'(3(5x,i2),5x,f16.9)') i1dir,i2dir,i3dir,&
&           d3cart(1,i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)
           if ((blkflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1).or.&
&           (carflg(i1dir,i1pert,i2dir,i1pert,i3dir,i1pert)/=1)) flag = 0
         end do
       end do
     end do

     if (pead==0.and.(dtset%nonlinear_info>0)) then

       d3cart_1(:,:,i1pert,:,i1pert,:,i1pert) = &
&       d3cart_1(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb
       d3cart_2(:,:,i1pert,:,i1pert,:,i1pert) = &
&       d3cart_2(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb
       d3cart_8(:,:,i1pert,:,i1pert,:,i1pert) = &
&       d3cart_8(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb
       d3cart_9(:,:,i1pert,:,i1pert,:,i1pert) = &
&       d3cart_9(:,:,i1pert,:,i1pert,:,i1pert)*16*(pi**2)*(Bohr_Ang**2)*1.0d-8*eps0/e_Cb

       theunit = ab_out

       write(small_msg,'(a)') ' ** Total :'
       call print_chi2(d3cart,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_psi1H1psi1 :'
       call print_chi2(d3cart_1,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_lambda1psi1psi1 :'
       call print_chi2(d3cart_2,small_msg,theunit)

       write(small_msg,'(a)') ' ** exc3 :'
       call print_chi2(d3cart_8,small_msg,theunit)

       write(small_msg,'(a)') ' ** exc3_paw :'
       call print_chi2(d3cart_9,small_msg,theunit)

     end if ! nonlinear_info > 0

     if (flag == 0) then
       write(message,'(a,a,a,a,a,a)')ch10,&
&       ' dfptnl_doutput: WARNING -',ch10,&
&       '  matrix of third-order energies incomplete,',ch10,&
&       '  non-linear optical coefficients may be wrong.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

   end if  ! d3e_pert1,d3e_pert2,d3e_pert3

   if (((maxval(d3e_pert1(1:natom))/=0).and.(d3e_pert2(natom+2)/=0).and. &
&   (d3e_pert3(natom+2)/=0)).or.&
   ((maxval(d3e_pert2(1:natom))/=0).and.(d3e_pert1(natom+2)/=0).and. &
&   (d3e_pert3(natom+2)/=0)).or.&
   ((maxval(d3e_pert3(1:natom))/=0).and.(d3e_pert2(natom+2)/=0).and. &
&   (d3e_pert1(natom+2)/=0))) then
!    Perform a check if all relevant elements are available

     flag = 1
     do i1pert = 1, natom
       do i1dir = 1, 3
         do i2dir = 1, 3
           do i3dir = 1, 3
             if ((blkflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
             (blkflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
             (blkflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
             if ((carflg(i1dir,i1pert,i2dir,natom+2,i3dir,natom+2) /= 1).or.&
             (carflg(i1dir,natom+2,i2dir,i1pert,i3dir,natom+2) /= 1).or.&
             (carflg(i1dir,natom+2,i2dir,natom+2,i3dir,i1pert) /= 1)) flag = 0
           end do
         end do
       end do
     end do

     write(ab_out,*)ch10
     write(ab_out,*)' First-order change in the electronic dielectric '
     write(ab_out,*)' susceptibility tensor (Bohr^-1)'
     write(ab_out,*)' induced by an atomic displacement'
     write(ab_out,*)'  atom  displacement'

     do i1pert = 1,natom
       do i1dir = 1,3
         write(ab_out,'(1x,i4,9x,i2,3(3x,f16.9))')i1pert,i1dir,&
&         d3cart(1,i1dir,i1pert,1,natom+2,:,natom+2)
         write(ab_out,'(16x,3(3x,f16.9))')&
&         d3cart(1,i1dir,i1pert,2,natom+2,:,natom+2)
         write(ab_out,'(16x,3(3x,f16.9))')&
&         d3cart(1,i1dir,i1pert,3,natom+2,:,natom+2)
       end do
       write(ab_out,*)
     end do

     if (flag == 0) then
       write(message,'(a,a,a,a,a,a)')ch10,&
&       ' dfptnl_doutput: WARNING -',ch10,&
&       '  matrix of third-order energies incomplete,',ch10,&
&       '  changes in the dielectric susceptibility may be wrong.'
       call wrtout(ab_out,message,'COLL')
       call wrtout(std_out,message,'COLL')
     end if

     if (pead==0.and.(dtset%nonlinear_info>0)) then
       theunit = ab_out

       write(small_msg,'(a)') ' ** Total :'
       call print_dchidtau(d3cart,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_psi1H1psi1 :'
       call print_dchidtau(d3cart_1,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_lambda1psi1psi1 :'
       call print_dchidtau(d3cart_2,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_lambda1psi0S1psi1 :'
       call print_dchidtau(d3cart_3,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_psi0H2psi1a :'
       call print_dchidtau(d3cart_4,small_msg,theunit)

       write(small_msg,'(a)') ' ** sum_psi0H2psi1b :'
       call print_dchidtau(d3cart_5,small_msg,theunit)

       write(small_msg,'(a)') ' ** eHxc21_paw :'
       call print_dchidtau(d3cart_6,small_msg,theunit)

       write(small_msg,'(a)') ' ** eHxc21_nhat :'
       call print_dchidtau(d3cart_7,small_msg,theunit)

       write(small_msg,'(a)') ' ** exc3 :'
       call print_dchidtau(d3cart_8,small_msg,theunit)

       write(small_msg,'(a)') ' ** exc3_paw :'
       call print_dchidtau(d3cart_9,small_msg,theunit)

     end if ! nonlinear_info > 0

   end if  ! d3e_pert1,d3e_pert2,d3e_pert3
 end if   ! mpi_enreg%me

! TO OPTIMIZE DEALLOCATION !
 if (pead/=0) then
   ABI_DEALLOCATE(cgindex)
   ABI_DEALLOCATE(kg_neigh)
   ABI_DEALLOCATE(kneigh)
   ABI_DEALLOCATE(kptindex)
   ABI_DEALLOCATE(kpt3)
   ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   ABI_DEALLOCATE(mpi_enreg%mkmem)
   ABI_DEALLOCATE(mvwtk)
   ABI_DEALLOCATE(pwind)
 end if
 ABI_DEALLOCATE(atindx)
 ABI_DEALLOCATE(atindx1)
 ABI_DEALLOCATE(blkflg)
 ABI_DEALLOCATE(blkflg_tmp)
 ABI_DEALLOCATE(blkflg_sav)
 ABI_DEALLOCATE(carflg)
 ABI_DEALLOCATE(carflg_tmp)
 ABI_DEALLOCATE(carflg_sav)
 ABI_DEALLOCATE(cg)
 ABI_DEALLOCATE(d3cart)
 ABI_DEALLOCATE(d3etot)
 ABI_DEALLOCATE(d3cart_1)
 ABI_DEALLOCATE(d3cart_2)
 ABI_DEALLOCATE(d3cart_3)
 ABI_DEALLOCATE(d3cart_4)
 ABI_DEALLOCATE(d3cart_5)
 ABI_DEALLOCATE(d3cart_6)
 ABI_DEALLOCATE(d3cart_7)
 ABI_DEALLOCATE(d3cart_8)
 ABI_DEALLOCATE(d3cart_9)
 ABI_DEALLOCATE(d3e_1)
 ABI_DEALLOCATE(d3e_2)
 ABI_DEALLOCATE(d3e_3)
 ABI_DEALLOCATE(d3e_4)
 ABI_DEALLOCATE(d3e_5)
 ABI_DEALLOCATE(d3e_6)
 ABI_DEALLOCATE(d3e_7)
 ABI_DEALLOCATE(d3e_8)
 ABI_DEALLOCATE(d3e_9)
 ABI_DEALLOCATE(d3e_pert1)
 ABI_DEALLOCATE(d3e_pert2)
 ABI_DEALLOCATE(d3e_pert3)
 ABI_DEALLOCATE(eigen0)
 ABI_DEALLOCATE(rhog)
 ABI_DEALLOCATE(rhor)
 ABI_DEALLOCATE(nhat)
 ABI_DEALLOCATE(nhatgr)
 ABI_DEALLOCATE(rfpert)
 ABI_DEALLOCATE(grxc)
 ABI_DEALLOCATE(kg)
 ABI_DEALLOCATE(kxc)
 ABI_DEALLOCATE(k3xc)
 ABI_DEALLOCATE(indsym)
 ABI_DEALLOCATE(indsy1)
 ABI_DEALLOCATE(nattyp)
 ABI_DEALLOCATE(npwarr)
 ABI_DEALLOCATE(symrec)
 ABI_DEALLOCATE(symrc1)
 ABI_DEALLOCATE(symaf1)
 ABI_DEALLOCATE(symrl1)
 ABI_DEALLOCATE(tnons1)
 ABI_DEALLOCATE(irrzon1)
 ABI_DEALLOCATE(phnons1)
 ABI_DEALLOCATE(symq)
 ABI_DEALLOCATE(ph1d)
 ABI_DEALLOCATE(ph1df)
 ABI_DEALLOCATE(vtrial)
 ABI_DEALLOCATE(vxc)
 call pawfgr_destroy(pawfgr)
 if (psps%usepaw==1) then
   call pawang_free(pawang1)
   call pawrhoij_free(pawrhoij)
   call paw_an_free(paw_an)
   call paw_ij_free(paw_ij)
   call pawfgrtab_free(pawfgrtab)
 end if
 ABI_DATATYPE_DEALLOCATE(pawrhoij)
 ABI_DATATYPE_DEALLOCATE(paw_an)
 ABI_DATATYPE_DEALLOCATE(paw_ij)
 ABI_DATATYPE_DEALLOCATE(pawfgrtab)

!Clean the header
 call hdr_free(hdr)

!As the etotal energy has no meaning here, we set it to zero
!(to avoid meaningless side-effects when comparing ouputs...)
 etotal = zero

 call timab(501,2,tsec)

 DBG_EXIT("COLL")

 contains
!!***

!!****f* nonlinear/print_chi2
!! NAME
!! print_chi2
!!
!! FUNCTION
!! Print a third derivative tensor. Used only in nonlinear
!!
!! INPUTS
!!  d3cart0 = the tensor to print
!!  msg = a short message printed before the tensor
!!  theunit = unit where the tensor is written
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!
!! SOURCE

   subroutine print_chi2(d3cart0,msg,theunit)

     integer,intent(in) :: theunit
     character(len=30) :: msg
     real(dp) :: elem1,elem2
     real(dp),intent(in) :: d3cart0(2,3,mpert,3,mpert,3,mpert)

! *************************************************************************

     write(theunit,'(2a)') ch10,msg
     do i1dir = 1, 3
       do i2dir = 1, 3
         do i3dir = 1, 3
           elem1 = d3cart0(1,i1dir,natom+2,i2dir,natom+2,i3dir,natom+2)
           elem2 = d3cart0(2,i1dir,natom+2,i2dir,natom+2,i3dir,natom+2)
           write(theunit,'(3(5x,i2),5x,f16.9,2x,f16.9)') i1dir,i2dir,i3dir,elem1,elem2
         end do
       end do
     end do

   end subroutine print_chi2
!!***

!!****f* nonlinear/print_dchidtau
!! NAME
!! print_dchidtau
!!
!! FUNCTION
!! Print a third derivative tensor. Used only in nonlinear
!!
!! INPUTS
!!  d3cart0 = the tensor to print
!!  msg = a short message printed before the tensor
!!  theunit = unit where the tensor is written
!!
!! OUTPUT
!!
!! SIDE EFFECTS
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!
!! SOURCE

   subroutine print_dchidtau(d3cart0,msg,theunit)

     integer,intent(in) :: theunit
     character(len=30) :: msg
     real(dp),intent(in) :: d3cart0(2,3,mpert,3,mpert,3,mpert)

! *************************************************************************

     write(theunit,'(a)') msg
     do i1pert = 1,natom
       do i1dir = 1,3
         write(theunit,'(1x,i4,9x,i2,3(3x,f16.9),3(3x,f16.9))')i1pert,i1dir,&
  &      d3cart0(1,i1dir,i1pert,1,natom+2,:,natom+2),d3cart0(2,i1dir,i1pert,1,natom+2,:,natom+2)
         write(theunit,'(16x,3(3x,f16.9),3(3x,f16.9))')&
  &      d3cart0(1,i1dir,i1pert,2,natom+2,:,natom+2),d3cart0(2,i1dir,i1pert,2,natom+2,:,natom+2)
         write(theunit,'(16x,3(3x,f16.9),3(3x,f16.9))')&
  &      d3cart0(1,i1dir,i1pert,3,natom+2,:,natom+2),d3cart0(2,i1dir,i1pert,3,natom+2,:,natom+2)
       end do
     end do

  end subroutine print_dchidtau
!!***

end subroutine nonlinear
!!***

!!****f* ABINIT/initmv
!! NAME
!! initmv
!!
!! FUNCTION
!! Initialize finite difference calculation of the ddk im dfptnl_mv.f
!!
!! INPUTS
!!  dtset <type(dataset_type)> = all input variables in this dataset
!!  gmet(3,3) = reciprocal space metric tensor in bohr**-2
!!  kg(3,mpw*mkmem) = reduced (integer) coordinates of G vecs in basis sphere
!!  kneigh(30,nkpt2) = index of the neighbours of each k-point
!!  kg_neigh(30,nkpt2,3) = necessary to construct the vector joining a k-point
!!                         to its nearest neighbour in case of a single k-point,
!!                         a line of k-points or a plane of k-points.
!!                         See getshell.F90 for details
!!  kptindex(2,nkpt3)= index of the k-points in the reduced BZ
!!                     related to a k-point in the full BZ
!!  kpt3(3,nkpt3) = reduced coordinates of k-points in the full BZ
!!  mband = maximum number of bands
!!  mkmem = number of k points which can fit in memory
!!  mpi_enreg = informations about MPI parallelization
!!  mpw = maximum number of plane waves
!!  nband(nkpt*nsppol)=number of bands at each k point, for each polarization
!!  nkpt2 = number of k-points in the reduced BZ
!!  nkpt3 = number of k-points in the full BZ
!!  nneigh = total number of neighbours required to evaluate the finite
!!          difference formula
!!  npwarr(nkpt2)=number of planewaves at each k point
!!  nsppol = number of spin polarizations
!!  occ(mband*nkpt*nsppol) = occupation number for each band for each k
!!
!! OUTPUT
!! cgindex(nkpt2,nsppol) = for each k-point, cgindex tores the location of the WF in the cg array
!!        me = index of the current processor
!!        ineigh = index of a neighbour
!!        ikpt_loc = index of the iteration on ikpt on the current processor
!!        ikpt_rbz = index of a k-point in the reduced BZ
!! pwind(mpw,nneigh,mkmem) = array used to compute the overlap matrix smat between k-points
!!                           (see initberry.f for more explanations)
!!
!! PARENTS
!!      nonlinear
!!
!! CHILDREN
!!      kpgio,xmpi_sum
!!
!! SOURCE

subroutine initmv(cgindex,dtset,gmet,kg,kneigh,kg_neigh,kptindex,&
&  kpt3,mband,mkmem,mpi_enreg,mpw,nband,nkpt2,&
&  nkpt3,nneigh,npwarr,nsppol,occ,pwind)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mkmem,mpw,nkpt2,nkpt3,nneigh,nsppol
 type(MPI_type),intent(inout) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: kg(3,mpw*mkmem),kneigh(30,nkpt2),kg_neigh(30,nkpt2,3)
 integer,intent(in) :: nband(nkpt2*nsppol),npwarr(nkpt2),kptindex(2,nkpt3)
 integer,intent(out) :: cgindex(nkpt2,nsppol),pwind(mpw,nneigh,mkmem)
 real(dp),intent(in) :: gmet(3,3),kpt3(3,nkpt3),occ(mband*nkpt2*nsppol)

!Local variables-------------------------------
!scalars
 integer :: flag,iband,icg,ierr,ikg,ikg1,ikpt,ikpt2,ikpt_loc,ikpt_rbz
 integer :: index,ineigh,ipw,isppol,jpw,nband_k,mband_occ,mband_occ_k,npw_k
 integer :: npw_k1,orig,spaceComm
 real(dp) :: ecut_eff,sdeg
 character(len=500) :: message
!arrays
 integer :: dg(3)
 integer,allocatable :: kg1(:,:),kg1_k(:,:),npwar1(:),npwtot(:)
 real(dp) :: dk(3),dk_(3)
 real(dp),allocatable :: kpt1(:,:)

!************************************************************************

!DEBUG
!write(std_out,*)' initmv : enter '
!ENDDEBUG

 if (xmpi_paral== 1) then
!  BEGIN TF_CHANGES
   spaceComm=mpi_enreg%comm_cell
!  END TF_CHANGES
   mpi_enreg%kpt_loc2ibz_sp(:,:,:) = 0
   mpi_enreg%mkmem(:) = 0
 end if

 ecut_eff = dtset%ecut*(dtset%dilatmx)**2
 ABI_ALLOCATE(kg1_k,(3,mpw))
 ABI_ALLOCATE(kg1,(3,mkmem*mpw))
 ABI_ALLOCATE(kpt1,(3,nkpt2))
 ABI_ALLOCATE(npwar1,(nkpt2))
 ABI_ALLOCATE(npwtot,(nkpt2))
 kg1_k(:,:) = 0
 pwind(:,:,:) = 0
 cgindex(:,:) = 0

!Compute the number of occupied bands.
!Check that it is the same for every k-point and that
!nband(ikpt) is equal to this value

 if (nsppol == 1) then
   sdeg = two
 else if (nsppol == 2) then
   sdeg = one
 end if

!DEBUG
!write(std_out,*)' list of nband '
!do isppol = 1, nsppol
!do ikpt = 1, nkpt2
!
!nband_k = nband(ikpt + (isppol - 1)*nkpt2)
!write(std_out,*)' isppol, ikpt, nband_k=',isppol, ikpt, nband_k
!end do
!end do
!ENDDEBUG

 index = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt2

     mband_occ_k = 0
     nband_k = nband(ikpt + (isppol - 1)*nkpt2)

     do iband = 1, nband_k
       index = index + 1
       if (abs(occ(index) - sdeg) < tol8) mband_occ_k = mband_occ_k + 1
     end do

     if (nband_k /= mband_occ_k) then
       write(message,'(a,a,a)')&
&       '  In a non-linear response calculation, nband must be equal ',ch10,&
&       '  to the number of valence bands.'
       MSG_ERROR(message)
     end if

!    Note that the number of bands can be different for spin up and spin down
     if (ikpt > 1) then
       if (mband_occ /= mband_occ_k) then
         message = 'The number of valence bands is not the same for every k-point'
         MSG_ERROR(message)
       end if
     else
       mband_occ = mband_occ_k
     end if

   end do                ! close loop over ikpt
 end do                ! close loop over isppol

!Find the location of each wavefunction

 icg = 0
 do isppol = 1, nsppol
   do ikpt = 1, nkpt2


!    fab: inserted the shift due to the spin...

     nband_k = dtset%nband(ikpt+(isppol - 1)*nkpt2)
     npw_k = npwarr(ikpt)

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,mpi_enreg%me)) cycle

     cgindex(ikpt,isppol) = icg
     icg = icg + dtset%nspinor*npw_k*nband_k

   end do
 end do


!Build pwind

 do ineigh = 1, nneigh

   do ikpt = 1, nkpt2
     ikpt2  = kneigh(ineigh,ikpt)
     ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ
     kpt1(:,ikpt) = dtset%kptns(:,ikpt_rbz)
   end do

!  Set up the basis sphere of plane waves at kpt1
   kg1(:,:) = 0
   call kpgio(ecut_eff,dtset%exchn2n3d,gmet,dtset%istwfk,kg1,&
&   kpt1,mkmem,dtset%nband,nkpt2,'PERS',mpi_enreg,mpw,&
&   npwar1,npwtot,dtset%nsppol)

   ikg = 0 ; ikg1 = 0 ; ikpt_loc = 0

   if(dtset%nsppol/=1)then
     if(mpi_enreg%nproc/=1)then
       message = ' At present, non-linear response calculations for spin-polarized system cannot be done in parallel.'
       MSG_ERROR(message)
     else
       isppol=1
     end if
   else
     isppol=1
   end if

   do ikpt = 1, nkpt2

     nband_k = dtset%nband(ikpt+(isppol - 1)*nkpt2)
     ikpt2  = kneigh(ineigh,ikpt)
     ikpt_rbz = kptindex(1,ikpt2)   ! index of the k-point in the reduced BZ

     if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,-1,mpi_enreg%me))  cycle

     ikpt_loc = ikpt_loc + 1

     mpi_enreg%kpt_loc2ibz_sp(mpi_enreg%me, ikpt_loc, 1) = ikpt

     flag = 0
     npw_k = npwarr(ikpt)
     npw_k1 = npwarr(ikpt_rbz)
     dk_(:) = kpt3(:,ikpt2) - dtset%kptns(:,ikpt)
     dk(:)  = dk_(:) - nint(dk_(:)) + real(kg_neigh(ineigh,ikpt,:),dp)
     dg(:)  = nint(dk(:) - dk_(:))


     if (kptindex(2,ikpt2) == 0) then
       kg1_k(:,1:npw_k1) = kg1(:,ikg1+1:ikg1+npw_k1)
       if (dg(1)==0.and.dg(2)==0.and.dg(3)==0) flag = 1
     else
       kg1_k(:,1:npw_k1) = -1*kg1(:,ikg1+1:ikg1+npw_k1)
     end if

     orig = 1
     do ipw = 1, npw_k
       do jpw = orig, npw_k1

         if ((kg(1,ikg + ipw) == kg1_k(1,jpw) - dg(1)).and. &
&         (kg(2,ikg + ipw) == kg1_k(2,jpw) - dg(2)).and. &
&         (kg(3,ikg + ipw) == kg1_k(3,jpw) - dg(3)))  then

           pwind(ipw,ineigh,ikpt_loc) = jpw
           if (flag == 1)  orig = jpw + 1
           exit

         end if

       end do
     end do

     ikg = ikg + npw_k
     ikg1 = ikg1 + npw_k1

   end do     ! close loop over k-points

 end do    ! close loop over ineigh
 mpi_enreg%mkmem(mpi_enreg%me) = mkmem

 call xmpi_sum(mpi_enreg%kpt_loc2ibz_sp,spaceComm,ierr)
 call xmpi_sum(mpi_enreg%mkmem,spaceComm,ierr)

!----------------------------------------------------------------------------

 ABI_DEALLOCATE(kg1)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(kpt1)
 ABI_DEALLOCATE(npwar1)
 ABI_DEALLOCATE(npwtot)

end subroutine initmv
!!***

end module m_nonlinear
!!***

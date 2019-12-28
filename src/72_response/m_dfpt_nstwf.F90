!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_nstwf
!! NAME
!!  m_dfpt_nstwf
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 2008-2019 ABINIT group ()
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

module m_dfpt_nstwf

 use defs_basis
 use m_xmpi
 use m_errors
 use m_abicore
 use m_wfk
 use m_hamiltonian
 use m_cgtools
 use m_nctk
 use m_dtset
 use m_dtfil


 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_time,     only : timab
 use m_io_tools, only : file_exists
 use m_geometry, only : stresssym
 use m_dynmat,   only : dfpt_sygra
 use m_mpinfo,   only : destroy_mpi_enreg, initmpi_seq, proc_distrb_cycle, proc_distrb_band
 use m_hdr,      only : hdr_skip
 use m_occ,      only : occeig
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_paw_an,   only : paw_an_type, paw_an_reset_flags
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_copy, pawrhoij_nullify, &
&                       pawrhoij_init_unpacked, pawrhoij_mpisum_unpacked, pawrhoij_inquire_dim
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_get, pawcprj_copy
 use m_pawdij,   only : pawdijfr
 use m_pawfgr,   only : pawfgr_type
 use m_paw_mkrho,only : pawmkrho
 use m_paw_nhat, only : pawnhatfr
 use m_paw_dfpt, only : pawdfptenergy
 use m_kg,       only : mkkin, kpgstr, mkkpg
 use m_fft,      only : fftpac
 use m_spacepar, only : hartrestr, symrhg
 use m_initylmg, only : initylmg
 use m_mkffnl,   only : mkffnl
 use m_getgh1c,  only : getgh1c, getdc1
 use m_dfpt_mkrho, only : dfpt_accrho
 use m_atm2fft,    only : dfpt_atm2fft
 use m_mkcore,     only : dfpt_mkcore
 use m_dfpt_mkvxc,    only : dfpt_mkvxc, dfpt_mkvxc_noncoll
 use m_dfpt_mkvxcstr, only : dfpt_mkvxcstr
 use m_mklocl,     only : dfpt_vlocal, vlocalstr
 use m_cgprj,           only : getcprj

 implicit none

 private
!!***

 public :: dfpt_nstpaw
 public :: dfpt_nstwf
!!***

contains
!!***

!!****f* ABINIT/dfpt_nstpaw
!! NAME
!! dfpt_nstpaw
!!
!! FUNCTION
!! Initially designed for PAW approach, but works also for NCPP.
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, for a whole row of
!! mixed derivatives (including diagonal terms contributing
!! to non-stationnary 2nd-order total energy).
!! Compared with NC-pseudopotentials, PAW contributions include:
!!  - changes of the overlap between 0-order wave-functions,
!!  - on-site contributions.
!!
!! COPYRIGHT
!! Copyright (C) 2010-2019 ABINIT group (MT, AM)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex=if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k projected with non-local projectors
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q projected with non-local projectors
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=1st-order eigenvalues at k,q (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for RF symmetries
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfftf,nkxc)=exchange and correlation kernel
!!  mgfftf=maximum size of 1D FFTs for the "fine" grid (see NOTES in respfn.F90)
!!  mkmem =number of k points treated by this node.
!!  mkqmem =number of k+q points treated by this node (GS data).
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid
!!  nhat(nfft,nspden*nhatdim)= -PAW only- compensation density
!!  nhat1(cplex*nfftf,nspden*usepaw)=1st-order compensation charge density (PAW)
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with i perturbation
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used otherwise, cplex*nfftf
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band at each k+q point of the reduced BZ
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k in the reduced BZ
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_an1(natom) <type(paw_an_type)>=1st-order paw arrays given on angular mesh for the perturbation (j1)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  paw_ij1(natom) <type(paw_ij_type)>=1st-order paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastructure containing only the symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic transl. phases, for RF symmetries
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information for the "fine" grid
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor(nfft,nspden)=array for GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfftf,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space in terms
!!  ucvol=unit cell volume in bohr**3.
!!  usecprj= 1 if cprj, cprjq arrays are stored in memory
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  vhartr1(cplex*nfft)=1-order Hartree potential
!!  vpsp1(cplex*nfftf)=first-order derivative of the ionic potential
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vtrial1(cplex*nfftf,nspden)= RF 1st-order potential (Hartree).
!!  vxc(nfftf,nspden)=XC GS potential
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  d2ovl(2,3,mpert,3,mpert*usepaw)=overlap contributions to the 2DTEs (PAW only)
!!  eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!        PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008) [[cite:Audouze2008]]
!!
!! NOTES
!!   We perform here the computation of
!!     delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!!     see PRB 78, 035105 (2008), Eq. (42) [[cite:Audouze2008]]
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      appdig,destroy_mpi_enreg
!!      dfpt_accrho,dfpt_atm2fft,dfpt_mkcore,dfpt_mkvxc,dfpt_mkvxc_noncoll
!!      dfpt_mkvxcstr,dfpt_sygra,dfpt_vlocal,dotprod_g,dotprod_vn,fftpac
!!      getcprj,getdc1,getgh1c,hartrestr,init_hamiltonian,init_rf_hamiltonian
!!      initmpi_seq,initylmg,kpgstr
!!      mkffnl,mkkin,mkkpg,occeig
!!      paw_an_reset_flags,paw_ij_free,paw_ij_init,paw_ij_nullify
!!      paw_ij_reset_flags,pawcprj_alloc,pawcprj_copy,pawcprj_free,pawcprj_get
!!      pawdfptenergy,pawdij2e1kb,pawdijfr,pawmkrho,pawnhatfr,pawrhoij_alloc
!!      pawrhoij_free,pawrhoij_init_unpacked,pawrhoij_mpisum_unpacked,projbd
!!      stresssym,symrhg,timab,vlocalstr,wfk_close,wfk_open_read,wfk_read_bks
!!      wrtout,xmpi_barrier,xmpi_sum
!!
!! SOURCE

subroutine dfpt_nstpaw(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&
&                  eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&
&                  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&                  mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,ncpgr,nfftf,ngfftf,nhat,nhat1,&
&                  nkpt_rbz,nkxc,npwarr,npwar1,nspden,nspinor,nsppol,nsym1,n3xccc,occkq,occ_rbz,&
&                  paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,&
&                  pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhog,rhor,rhor1,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,&
&                  ucvol,usecprj,usepaw,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,&
&                  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr1)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mgfftf,mkmem,mkqmem,mk1mem,mpert,mpw,mpw1
 integer,intent(in) :: ncpgr,nfftf,nkpt_rbz,nkxc,nspden,nspinor,nsppol,nsym1
 integer,intent(in) :: n3xccc,usecprj,usepaw,usexcnhat,useylmgr1
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(out) :: eovl1
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(MPI_type),intent(in) :: mpi_enreg
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: nattyp(dtset%ntypat),nband_rbz(nkpt_rbz*nsppol)
 integer,intent(in) :: indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: ngfftf(18),npwarr(nkpt_rbz),npwar1(nkpt_rbz)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*dtset%mband*mkmem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*nspinor*dtset%mband*mkqmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*dtset%mband*mk1mem*nsppol)
 real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*nsppol),eigen0(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: kxc(nfftf,nkxc),nhat(nfftf,nspden),nhat1(cplex*nfftf,nspden*usepaw)
 real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom),rhog(2,nfftf)
 real(dp),intent(in) :: rhor(cplex*nfftf,nspden),rhor1(cplex*nfftf,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: tnons1(3,nsym1),vhartr1(cplex*nfftf),vtrial1(cplex*nfftf,nspden),vxc(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),target,intent(in) :: vpsp1(cplex*nfftf),vtrial(nfftf,nspden),xccc3d1(cplex*n3xccc)
 real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert),d2ovl(2,3,mpert,3,mpert*usepaw)
 type(pawcprj_type),intent(in) :: cprj(dtset%natom,nspinor*dtset%mband*mkmem*nsppol*usecprj)
 type(pawcprj_type),intent(in) :: cprjq(dtset%natom,nspinor*dtset%mband*mkqmem*nsppol*usecprj)
 type(paw_an_type),intent(in) :: paw_an(:)
 type(paw_an_type),intent(inout) :: paw_an1(:)
 type(paw_ij_type),intent(in) :: paw_ij(:)
 type(paw_ij_type),intent(inout) :: paw_ij1(:)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(:)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(:)
 type(pawrhoij_type),intent(in) :: pawrhoij1(:)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: tim_fourwf=18,tim_getgh1c=2,tim_projbd=3,formeig1=1
 integer :: bd2tot_index,bdtot_index,berryopt,bufsz,choice,cpopt,cplex_rhoij,ddkcase
 integer :: dimffnl,dimffnl1,dimffnl1_idir1,dimylmgr1
 integer :: ia,iatom,iband,ibg,ibgq,ibg1,icg,icgq,icg1,ider,idir0,idir1,idir_cprj
 integer :: ierr,ii,ikg,ikg1,ikpt,ikpt_me,ilmn,iorder_cprj,ipert1
 integer :: ispden,isppol,istwf_k,istr,istr1,itypat,jband,jj,kdir1,kpert1,master,mcgq,mcprjq
 integer :: mdir1,me,mpert1,my_natom,my_comm_atom,my_nsppol,nband_k,nband_kocc,need_ylmgr1
 integer :: nfftot,nkpg,nkpg1,nkpt_me,npw_,npw_k,npw1_k,nspden_rhoij
 integer :: nvh1,nvxc1,nzlmopt_ipert,nzlmopt_ipert1,optlocal,optnl
 integer :: option,opt_gvnlx1,qphase_rhoij,sij_opt,spaceworld,usevnl,wfcorr,ik_ddk
 real(dp) :: arg,doti,dotr,dot1i,dot1r,dot2i,dot2r,dot3i,dot3r,elfd_fact,invocc,lambda,wtk_k
 logical :: force_recompute,has_dcwf,has_dcwf2,has_drho,has_ddk_file
 logical :: is_metal,is_metal_or_qne0,need_ddk_file,need_pawij10
 logical :: need_wfk,need_wf1,nmxc,paral_atom,qne0,t_exist
 character(len=500) :: msg
 character(len=fnlen) :: fiwfddk(3)
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
 type(MPI_type) :: mpi_enreg_seq
!arrays
 integer :: ddkfil(3),my_spintab(2),nband_tmp(1),npwar1_tmp(1)
 integer,allocatable :: jpert1(:),jdir1(:),kg1_k(:,:),kg_k(:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dum1(1,1),dum2(1,1),dum3(1,1),epawnst(2),kpoint(3),kpq(3)
 real(dp) :: sumelfd(2),symfact(3),tsec(2),ylmgr_dum(1,3,1)
 real(dp),allocatable :: buffer(:),ch1c(:,:,:,:),cs1c(:,:,:,:)
 real(dp),allocatable :: cwave0(:,:),cwavef(:,:),dcwavef(:,:)
 real(dp),allocatable :: doccde_k(:),doccde_kq(:)
 real(dp),allocatable :: dnhat1(:,:),drhoaug1(:,:,:,:)
 real(dp),allocatable :: drhor1(:,:),drho1wfg(:,:),drho1wfr(:,:,:)
 real(dp),allocatable :: d2nl_elfd(:,:),dkinpw(:)
 real(dp),allocatable :: d2nl_k(:,:),d2ovl_drho(:,:,:,:,:),d2ovl_k(:,:)
 real(dp),allocatable :: eig_k(:),eig_kq(:),eig1_k(:)
 real(dp),allocatable,target :: e1kbfr_spin(:,:,:,:,:,:),ffnlk(:,:,:,:),ffnl1(:,:,:,:)
 real(dp),allocatable :: gh1(:,:),gs1(:,:),gvnlx1(:,:),kinpw1(:),kpg_k(:,:),kpg1_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),ph3d(:,:,:),ph3d1(:,:,:),rhotmp(:,:),rocceig(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylm1_k(:,:),ylmgr1_k(:,:,:),vtmp1(:,:),vxc10(:,:)
 real(dp),allocatable,target :: work(:,:,:),e1kb_work(:,:,:,:)
 real(dp),pointer :: e1kbfr(:,:,:,:,:),e1kb_ptr(:,:,:,:)
 real(dp),pointer :: ffnl1_idir1(:,:,:,:),vhartr01(:),vpsp1_idir1(:),xccc3d1_idir1(:)
 type(pawcprj_type),allocatable :: dcwaveprj(:,:)
 type(pawcprj_type),allocatable,target :: cwaveprj0(:,:)
 type(pawcprj_type),pointer :: cwaveprj0_idir1(:,:)
 type(paw_ij_type),allocatable :: paw_ij10(:,:)
 type(pawrhoij_type),target,allocatable :: pawdrhoij1(:,:)
 type(pawrhoij_type),pointer :: pawdrhoij1_unsym(:,:)
 type(wfk_t) :: ddks(3)

! *********************************************************************

 DBG_ENTER("COLL")

!Keep track of total time spent in dfpt_nstpaw
 call timab(566,1,tsec)

!Not valid for PrintBandByBand
 if (dtset%prtbbb/=0) then
   MSG_BUG('not yet valid for prtbbb/=0!')
 end if

!NCPP restrictions
 if (usepaw==0) then
!  cprj cannot be used
   if (usecprj/=0) then
     MSG_BUG('NCPP: usecprj should be 0!')
   end if
!  d2ovl cannot be used
   if (size(d2ovl)/=0) then
     MSG_BUG('NCPP: d2ovl should not be allocated!')
   end if
 end if

!PAW restrictions
 if (usepaw==1) then
!  Test on FFT grid sizes
   if (pawfgr%nfft/=nfftf) then
     MSG_BUG('PAW: wrong values for nfft, nfftf!')
   end if
!  Test gradients of cprj
   if (ipert<=dtset%natom.and.ncpgr/=3) then
     MSG_BUG('PAW: wrong value of ncpgr for ipert<=natom!')
   end if
   if (ipert==dtset%natom+1.and.ncpgr/=1) then
     MSG_BUG('PAW: wrong value of ncpgr for ipert=natom+1!')
   end if
   if (ipert==dtset%natom+2.and.ncpgr/=3) then
     MSG_BUG('PAW: wrong value of ncpgr for ipert=natom+2!')
   end if
   if ((ipert==dtset%natom+3.or.ipert==dtset%natom+4).and.ncpgr/=1) then
     MSG_BUG('PAW: wrong value of ncpgr for ipert=natom+3 or 4!')
   end if
!  Test on availability of DijHartree and XC on-site potentials
   if (mpi_enreg%my_natom>0.and.ipert/=dtset%natom+1)  then
     if (paw_ij1(1)%has_dijhartree==0.or.paw_an1(1)%has_vxc==0) then
       msg='PAW: paw_ij1%dijhartree and paw=_an1%vxc1 should be allocated !'
     end if
   end if
 end if

!Set up parallelism
 master=0;me=mpi_enreg%me_kpt
 spaceworld=mpi_enreg%comm_cell
 paral_atom=(mpi_enreg%my_natom/=dtset%natom)
 my_comm_atom=mpi_enreg%comm_atom
 my_natom=mpi_enreg%my_natom
 my_atmtab=>mpi_enreg%my_atmtab
 my_spintab=mpi_enreg%my_isppoltab
 my_nsppol=count(my_spintab==1)

!Fake MPI data to be used in sequential calls to parallel routines
 call initmpi_seq(mpi_enreg_seq)
 mpi_enreg_seq%my_natom=dtset%natom

!Compute effective number of k-points
 nkpt_me=nkpt_rbz
 if(xmpi_paral==1)then
   nkpt_me=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       if (.not.(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me))) nkpt_me=nkpt_me+1
     end do
   end do
 end if

!Sizes for WF at k+q
 mcgq=mpw1*nspinor*dtset%mband*mkqmem*nsppol
 mcprjq=nspinor*dtset%mband*mkqmem*nsppol*usecprj

!Check ddk files (needed to compute electric field perturbations)
 ddkfil(:)=0
 do idir1=1,3
   ddkcase=idir1+dtset%natom*3
   call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk(idir1))
   t_exist = file_exists(fiwfddk(idir1))

   if (.not. t_exist) then
     ! Try netcdf file.
     t_exist = file_exists(nctk_ncify(fiwfddk(idir1)))
     if (t_exist) then
       fiwfddk(idir1) = nctk_ncify(fiwfddk(idir1))
       write(msg,"(3a)")"- File: ",trim(fiwfddk(idir1))," does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   if (t_exist) then
     ddkfil(idir1)=20+idir1 ! Note the use of unit numbers 21, 22 and 23
   end if
 end do
 has_ddk_file=(any(ddkfil(:)>0))

!Define the set of perturbations (j1)=(ipert1,idir1)
!The first perturbation must be (j1)=(j2)=(ipert,idir)
!because we need to compute <g|H^(j2)-Eps.S^(j2)|u0> first.
 if (ipert/=dtset%natom+1) then
   mpert1=0
   ABI_ALLOCATE(jpert1,(mpert))
   jpert1 = 0
   if (ipert/=dtset%natom+2.or.has_ddk_file) then
     mpert1=mpert1+1;jpert1(mpert1)=ipert
   end if
   do ipert1=1,mpert
     if (ipert1/=ipert.and.&
&     (ipert1<=dtset%natom.or.&
&     (ipert1==dtset%natom+2.and.has_ddk_file).or.&
&     ((ipert>dtset%natom.and.ipert/=dtset%natom+5).and.(ipert1==dtset%natom+3.or.ipert1==dtset%natom+4)).or. &
&     ((ipert1==dtset%natom+2).and.has_ddk_file))) then
       mpert1=mpert1+1;jpert1(mpert1)=ipert1
     end if
   end do
 else
   mpert1=1
   ABI_ALLOCATE(jpert1,(mpert1))
   jpert1(1)=dtset%natom+1
 end if
 mdir1=3
 ABI_ALLOCATE(jdir1,(mdir1))
 jdir1(1:3)= (/ (idir1,idir1=1,3) /)
 jdir1(1)=idir;jdir1(idir)=1

!Index of strain perturbation, if any
 istr=idir;if (ipert==dtset%natom+4) istr=idir+3

!Open ddk WF file(s)
 if (has_ddk_file) then
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if (ddkfil(idir1)/=0) then
       write(msg, '(a,a)') '-open ddk wf file :',trim(fiwfddk(idir1))
       call wrtout(std_out,msg,'COLL')
       call wrtout(ab_out,msg,'COLL')
       call wfk_open_read(ddks(idir1),fiwfddk(idir1),formeig1,dtset%iomode,ddkfil(idir1),spaceworld)
     end if
   end do
 end if

!Zero only portion of matrix to be computed here
 d2nl(:,:,1:dtset%natom+4,idir,ipert)=zero
 if (usepaw==1) d2ovl(:,:,1:dtset%natom+4,idir,ipert)=zero

!Update list of computed matrix elements
 do kpert1=1,mpert1
   ipert1=jpert1(kpert1)
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if ((ipert1<=dtset%natom).or.ipert1==dtset%natom+3.or.ipert1==dtset%natom+4.or.&
&     (ipert1==dtset%natom+1.and.((ddkfil(idir1)/=0).or.(dtset%rfdir(idir1)/=0.and.idir1<=idir))).or.&
&     (ipert1==dtset%natom+2.and.ddkfil(idir1)/=0).or.&
&     ((ipert1==dtset%natom+3.or.ipert1==dtset%natom+4).and.ddkfil(idir1)/=0)) then
       blkflg(idir1,ipert1,idir,ipert)=1
     end if
   end do
 end do

!Initialize most of the (1st-order) Hamiltonian
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
 call init_hamiltonian(gs_hamkq,psps,pawtab,nspinor,nsppol,nspden,dtset%natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& paw_ij=paw_ij,mpi_atmtab=my_atmtab,comm_atom=my_comm_atom,mpi_spintab=mpi_enreg%my_isppoltab,&
& usecprj=usecprj,nucdipmom=dtset%nucdipmom,use_gpu_cuda=dtset%use_gpu_cuda)

!Variables common to all perturbations
 arg=maxval(occ_rbz)-minval(occ_rbz)
 qne0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2>=tol14)
 is_metal=((dtset%occopt>=3.and.dtset%occopt<=8).or.(abs(arg)>tol8))
 is_metal_or_qne0=((is_metal).or.(qne0))
 nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
 ABI_ALLOCATE(ch1c,(2,dtset%mband,dtset%mband,nkpt_me))
 ch1c(:,:,:,:)=zero
 nzlmopt_ipert=0;nzlmopt_ipert1=0
 if (usepaw==1) then
   ABI_ALLOCATE(d2ovl_drho,(2,3,mpert,3,mpert))
   d2ovl_drho=zero
   if (dtset%pawnzlm/=0) then
     nzlmopt_ipert=1;if (dtset%nstep<2) nzlmopt_ipert=-1
     nzlmopt_ipert1=-1
   end if
   if (is_metal_or_qne0) then
     ABI_ALLOCATE(cs1c,(2,dtset%mband,dtset%mband,nkpt_me))
     cs1c(:,:,:,:)=zero
   end if
 end if

!Force the recomputation of on-site potentials and DijHartree
 if (usepaw==1) then
   call paw_an_reset_flags(paw_an1)
   call paw_ij_reset_flags(paw_ij1,dijhartree=.true.)
 end if

!LOOP OVER PERTURBATION TYPES (j1)
 do kpert1=1,mpert1
   ipert1=jpert1(kpert1)

!  Flag for use of DDK file
   need_ddk_file=(has_ddk_file.and.(ipert1==dtset%natom+1.or.ipert1==dtset%natom+2))

!  Factor to be applied for electric Field (Eff. charges and piezo. tensor are "minus" d2E)
   elfd_fact=one
   if ((ipert <=dtset%natom.or.ipert ==dtset%natom+3.or.ipert ==dtset%natom+4).and. &
&   (ipert1==dtset%natom+2)) elfd_fact=-one
   if ((ipert1<=dtset%natom.or.ipert1==dtset%natom+3.or.ipert1==dtset%natom+4).and. &
&   (ipert ==dtset%natom+2)) elfd_fact=-one

!  We want to compute delta_u^(j1))=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!  see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42)
   has_dcwf=.false.;has_dcwf2=.false.;has_drho=.false.
   if (usepaw==1) then
     has_dcwf =(ipert1/=dtset%natom+2)
     has_dcwf2=(ipert /=dtset%natom+2)
     has_drho =(has_dcwf.and.ipert1/=dtset%natom+1)
   end if

!  Select which WF are needed
   need_wfk=(.true.)
   need_wf1=(.true.)

!  Initialize data for NL 1st-order (j1) hamiltonian
   call init_rf_hamiltonian(cplex,gs_hamkq,ipert1,rf_hamkq,mpi_spintab=[0,0])

!  The following contributions are needed only for non-DDK perturbation:
!  - Frozen part of 1st-order Dij
!  - Contribution from local potential to dynamical matrix (due to Vxc^(j1)(tild_nc)+VH^(j1)(tild_nZc))
   if (ipert/=dtset%natom+1.and.ipert1/=dtset%natom+1) then

!    Allocations
     force_recompute=(usepaw==0) ! This is dangerous...
     nfftot=ngfftf(1)*ngfftf(2)*ngfftf(3)
     nvxc1=0;if (ipert1<=dtset%natom.or.ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) nvxc1=nspden
     nvh1=0;if (ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) nvh1=1
     ABI_ALLOCATE(vxc10,(cplex*nfftf,nvxc1))
     ABI_ALLOCATE(vhartr01,(cplex*nfftf*nvh1))
     need_pawij10=(usepaw==1)
     if (need_pawij10) then
       ABI_DATATYPE_ALLOCATE(paw_ij10,(my_natom,mdir1))
       ABI_ALLOCATE(e1kbfr_spin,(rf_hamkq%dime1kb1,rf_hamkq%dime1kb2,nspinor**2,cplex,mdir1,my_nsppol))
     else
       ABI_DATATYPE_ALLOCATE(paw_ij10,(0,0))
     end if

!    LOOP OVER PERTURBATION DIRECTIONS
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       istr1=idir1;if (ipert1==dtset%natom+4) istr1=idir1+3

!      Get first-order local potential and first-order pseudo core density
       if (ipert==ipert1.and.idir==idir1.and.(.not.force_recompute)) then
         vpsp1_idir1 => vpsp1
         xccc3d1_idir1 => xccc3d1
       else
         ABI_ALLOCATE(vpsp1_idir1,(cplex*nfftf))
         ABI_ALLOCATE(xccc3d1_idir1,(cplex*n3xccc))
         if (usepaw==1) then
           call dfpt_atm2fft(gs_hamkq%atindx,cplex,gmet,gprimd,gsqcut,istr1,ipert1,&
&           mgfftf,psps%mqgrid_vl,dtset%natom,1,nfftf,ngfftf,dtset%ntypat,&
&           ph1df,psps%qgrid_vl,dtset%qptn,dtset%typat,ucvol,psps%usepaw,xred,psps,pawtab,&
&           atmrhor1=xccc3d1_idir1,atmvlocr1=vpsp1_idir1,optv_in=1,optn_in=n3xccc/nfftf,optn2_in=1,&
&           vspl=psps%vlspl)
         else
           if(ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
             call vlocalstr(gmet,gprimd,gsqcut,istr1,mgfftf,mpi_enreg,psps%mqgrid_vl,dtset%natom,&
&             nattyp,nfftf,ngfftf,dtset%ntypat,ph1df,psps%qgrid_vl,ucvol,&
&             psps%vlspl,vpsp1_idir1)
           else
             call dfpt_vlocal(gs_hamkq%atindx,cplex,gmet,gsqcut,idir1,ipert1,mpi_enreg,psps%mqgrid_vl,&
&             dtset%natom,nattyp,nfftf,ngfftf,dtset%ntypat,ngfftf(1),ngfftf(2),ngfftf(3),&
&             ph1df,psps%qgrid_vl,dtset%qptn,ucvol,psps%vlspl,vpsp1_idir1,xred)
           end if
           if(psps%n1xccc/=0)then
             call dfpt_mkcore(cplex,idir1,ipert1,dtset%natom,dtset%ntypat,ngfftf(1),psps%n1xccc,&
&             ngfftf(2),ngfftf(3),dtset%qptn,rprimd,dtset%typat,ucvol,&
&             psps%xcccrc,psps%xccc1d,xccc3d1_idir1,xred)
           end if
         end if
       end if

!      Compute 1st-order non-local factors (Dij^(j1)_fr)
       if (need_pawij10) then
         call paw_ij_nullify(paw_ij10(:,idir1))
         call paw_ij_init(paw_ij10(:,idir1),cplex,nspinor,dtset%nsppol,dtset%nspden,&
&         0,dtset%natom,dtset%ntypat,dtset%typat,pawtab,has_dijfr=1,&
&         mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
         if (ipert/=ipert1.or.idir/=idir1.or.force_recompute) then
           option=0
           call pawdijfr(gprimd,idir1,ipert1,my_natom,dtset%natom,nfftf,ngfftf,&
&           nspden,nsppol,dtset%ntypat,option,paw_ij10(:,idir1),pawang,pawfgrtab,pawrad,&
&           pawtab,cplex,dtset%qptn,rprimd,ucvol,vpsp1_idir1,vtrial,vxc,xred,&
&           mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
         else
           do iatom=1,my_natom
             paw_ij10(iatom,idir1)%has_dijfr=paw_ij1(iatom)%has_dijfr
             if (paw_ij1(iatom)%has_dijfr==2) paw_ij10(iatom,idir1)%dijfr=paw_ij1(iatom)%dijfr
           end do
         end if
       end if

!      Get first-order exchange-correlation potential (core-correction contribution only)
       if (nvxc1>0) then
         if(ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
           option=0
           call dfpt_mkvxcstr(cplex,idir1,ipert1,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,&
&           nhat,nhat1,nkxc,nmxc,nspden,n3xccc,option,dtset%qptn,rhor,rhor1,rprimd,&
&           usepaw,usexcnhat,vxc10,xccc3d1_idir1)
         else
!          Non-collinear magnetism (should the second nkxc be nkxc_cur ?)
           if (nspden==4) then
             option=0
             call dfpt_mkvxc_noncoll(cplex,dtset%ixc,kxc,mpi_enreg,nfftf,ngfftf,dum1,0,dum2,0,dum3,0,nkxc,&
&             nmxc,nspden,n3xccc,1,option,dtset%qptn,dum1,dum1,rprimd,0,vxc,&
&             vxc10,xccc3d1_idir1)
           else
             call dfpt_mkvxc(cplex,dtset%ixc,kxc,mpi_enreg,nfftf,ngfftf,dum2,0,dum3,0,nkxc,&
&             nmxc,nspden,n3xccc,0,dtset%qptn,dum1,rprimd,0,vxc10,xccc3d1_idir1)
           end if
         end if
       end if

!      Get first-order Hartree potential (metric tensor contribution only)
       if (nvh1>0) then
         call hartrestr(gsqcut,idir1,ipert1,mpi_enreg,dtset%natom,&
&         nfftf,ngfftf,rhog,rprimd,vhartr01)
       end if

!      Get Hartree + xc + local contributions to dynamical matrix or elastic tensor
       if (ipert1<=dtset%natom.or.ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
         if (usepaw==0) then
!          vxc1 is integrated with the total 1st-order density (rhor1 )
!          vpsp1 is integrated with the 1st-order pseudo density (rhor1)
           call dotprod_vn(cplex,rhor1,dot1r,dot1i,nfftf,nfftot,nspden,2,vxc10,ucvol)
           call dotprod_vn(cplex,rhor1,dot2r,dot2i,nfftf,nfftot,1,2,vpsp1_idir1,ucvol)
         else if (usexcnhat/=0) then
!          vxc1 is integrated with the total 1st-order density (rhor1 including nhat1)
!          vpsp1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
           ABI_ALLOCATE(rhotmp,(cplex*nfftf,1))
           rhotmp(:,1)=rhor1(:,1)-nhat1(:,1)
           call dotprod_vn(cplex,rhor1,dot1r,dot1i,nfftf,nfftot,nspden,2,vxc10,ucvol)
           call dotprod_vn(cplex,rhotmp,dot2r,dot2i,nfftf,nfftot,1,2,vpsp1_idir1,ucvol)
           ABI_DEALLOCATE(rhotmp)
         else
!          vxc1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
!          vpsp1 is integrated with the 1st-order pseudo density (rhor1 without nhat1)
           ABI_ALLOCATE(rhotmp,(cplex*nfftf,nspden))
           rhotmp(:,:)=rhor1(:,:)-nhat1(:,:)
           call dotprod_vn(cplex,rhotmp,dot1r,dot1i,nfftf,nfftot,nspden,2,vxc10,ucvol)
           call dotprod_vn(cplex,rhotmp,dot2r,dot2i,nfftf,nfftot,1,2,vpsp1_idir1,ucvol)
           ABI_DEALLOCATE(rhotmp)
         end if
         if (nvh1>0) then
           call dotprod_vn(cplex,rhor1,dot3r,dot3i,nfftf,nfftot,1,2,vhartr01,ucvol)
         else
           dot3r=zero ; dot3i=zero
         end if
!        Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
         dotr=dot1r+dot2r+dot3r;doti=dot1i+dot2i+dot3i
!        In case ipert = natom+2, these lines compute the local part
!        of the Born effective charges from phonon and electric
!        field type perturbations, see eq. 43 of X. Gonze and C. Lee, PRB 55, 10355 (1997) [[cite:Gonze1997a]]
!        The minus sign is due to the fact that the effective charges
!        are minus the second derivatives of the energy
         if (ipert/=dtset%natom+1) then
           d2lo(1,idir1,ipert1,idir,ipert)=elfd_fact*dotr
           d2lo(2,idir1,ipert1,idir,ipert)=elfd_fact*doti
         end if
       end if ! ipert1<=natom

       if (ipert/=ipert1.or.idir/=idir1.or.force_recompute)  then
         ABI_DEALLOCATE(vpsp1_idir1)
         ABI_DEALLOCATE(xccc3d1_idir1)
       end if

!      End loop on directions
     end do

!    Free memory
     ABI_DEALLOCATE(vxc10)
     ABI_DEALLOCATE(vhartr01)

   else ! ddk perturbation
     d2lo(1:2,1:mdir1,ipert1,idir,ipert)=zero
     need_pawij10=.false.
     ABI_DATATYPE_ALLOCATE(paw_ij10,(0,0))
   end if

!  Prepare RF PAW files for reading and writing if mkmem, mkqmem or mk1mem==0
   iorder_cprj=0

!  Allocate arrays used to accumulate density change due to overlap
   if (has_drho) then
     ABI_ALLOCATE(drhoaug1,(cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),mdir1))
     ABI_ALLOCATE(drho1wfr,(cplex*dtset%nfft,dtset%nspden,mdir1))
     ABI_DATATYPE_ALLOCATE(pawdrhoij1,(my_natom,mdir1))
     if (paral_atom) then
       ABI_DATATYPE_ALLOCATE(pawdrhoij1_unsym,(dtset%natom,mdir1))
     else
       pawdrhoij1_unsym => pawdrhoij1
     end if
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       drho1wfr(:,:,idir1)=zero
       call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                            nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
       call pawrhoij_alloc(pawdrhoij1(:,idir1),cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&       dtset%nsppol,dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,use_rhoijp=1,use_rhoij_=0,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=my_atmtab)
       if (paral_atom) then
         call pawrhoij_alloc(pawdrhoij1_unsym(:,idir1),cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&         dtset%nsppol,dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,&
&         use_rhoijp=0,use_rhoij_=1)
       else
         call pawrhoij_init_unpacked(pawdrhoij1_unsym(:,idir1))
       end if
     end do
   end if

!  Initialize shifts for global arrays
   bdtot_index=0
   bd2tot_index=0
   ibg=0;icg=0
   ibg1=0;icg1=0
   ibgq=0;icgq=0

!  Has to get 1st-order non-local factors before the loop over spins
!  because this needs a communication over comm_atom (=comm_spinkpt)
   if (need_pawij10) then
     if (my_nsppol<nsppol) then
       ABI_ALLOCATE(e1kb_work,(rf_hamkq%dime1kb1,rf_hamkq%dime1kb2,nspinor**2,cplex))
     end if
     ii=0
     do isppol=1,nsppol
       if (my_spintab(isppol)==1) ii=ii+1
       if (my_spintab(isppol)/=1) e1kb_ptr => e1kb_work
       do kdir1=1,mdir1
         idir1=jdir1(kdir1)
         if (my_spintab(isppol)==1) e1kb_ptr => e1kbfr_spin(:,:,:,:,idir1,ii)
         call pawdij2e1kb(paw_ij10(:,idir1),isppol,my_comm_atom,e1kbfr=e1kb_ptr,mpi_atmtab=my_atmtab)
       end do
     end do
     if (my_nsppol<nsppol) then
       ABI_DEALLOCATE(e1kb_work)
     end if
   end if

!  LOOP OVER SPINS
   do isppol=1,nsppol

!    either this has to be inside the sppol loop or the size of ch1c and cs1c has to be nkpt_me*nsppol
     ikpt_me=0

!    Rewind (k+G) data if needed
     ikg=0;ikg1=0

!    Continue to initialize the GS/RF Hamiltonian
     call gs_hamkq%load_spin(isppol,with_nonlocal=.true.)
     if (need_pawij10) then
       ii=min(isppol,size(e1kbfr_spin,6))
       if (ii>0) e1kbfr => e1kbfr_spin(:,:,:,:,:,ii)
     end if

!    Initialize accumulation of density
     if (has_drho) drhoaug1(:,:,:,:)=zero

!    LOOP OVER K-POINTS
     do ikpt=1,nkpt_rbz

!      Load dimensions for this k-point
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       istwf_k=istwfk_rbz(ikpt)
       npw_k=npwarr(ikpt)
       npw1_k=npwar1(ikpt)

!      Skip loop if this k-point is not to be treated by this proc
       if (proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
         bdtot_index=bdtot_index+nband_k
         bd2tot_index=bd2tot_index+2*nband_k**2
         cycle ! Skip the rest of the k-point loop
       end if

!      Allocate/initialize local arrays and scalars for this k-point
       ABI_ALLOCATE(d2nl_k,(2,3))
       ABI_ALLOCATE(d2ovl_k,(2,3))
       ABI_ALLOCATE(eig_k,(nband_k))
       ABI_ALLOCATE(eig_kq,(nband_k))
       ABI_ALLOCATE(eig1_k,(2*nband_k**2))
       ABI_ALLOCATE(occ_k,(nband_k))
       d2nl_k(:,:)=zero;d2ovl_k(:,:)=zero
       eig_k (:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
       eig_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
       eig1_k(:)=eigen1(1+bd2tot_index:2*nband_k**2+bd2tot_index)
       occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
       nband_kocc=count(abs(occ_k(:))>tol8)
       kpoint(:)=kpt_rbz(:,ikpt)
       kpq(:)=kpoint(:);if (ipert1<dtset%natom+3) kpq(:)=kpq(:)+dtset%qptn(1:3)
       wtk_k=wtk_rbz(ikpt)
       need_ylmgr1=0;dimylmgr1=0
       nkpg=0;nkpg1=0
       ikpt_me=ikpt_me+1
       if (is_metal) then
!        For each pair of active bands (m,n), generates the ratios
!        rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
         ABI_ALLOCATE(doccde_k,(nband_k))
         ABI_ALLOCATE(doccde_kq,(nband_k))
         ABI_ALLOCATE(occ_kq,(nband_k))
         ABI_ALLOCATE(rocceig,(nband_k,nband_k))
         doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
         doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)
         occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
         call occeig(doccde_k,doccde_kq,eig_k,eig_kq,nband_k,dtset%occopt,occ_k,occ_kq,rocceig)
       end if

!      Take care of the npw and kg records in WF and DDK files
       if (need_ddk_file) then
         do kdir1=1,mdir1
           idir1=jdir1(kdir1)
           if (ddkfil(idir1)/=0)then
             !ik_ddk = wfk_findk(ddks(idir1), kpt_rbz(:,ikpt)
             ik_ddk = indkpt1(ikpt)
             npw_ = ddks(idir1)%hdr%npwarr(ik_ddk)
             if (npw_/=npw_k) then
               write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&               'For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&               'the number of plane waves in the ddk file is equal to', npw_,ch10,&
&               'while it should be ',npw_k
               MSG_ERROR(msg)
             end if

           end if
         end do
       end if

!      Allocate arrays used for NL form factors
       ABI_ALLOCATE(kg_k,(3,npw_k))
       ABI_ALLOCATE(kg1_k,(3,npw1_k))
       ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
       ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))
       if (psps%useylm==1.and.(need_ddk_file.or.ipert1==dtset%natom+1)) need_ylmgr1=1
       if (psps%useylm==1.and.(ipert1==dtset%natom+3.or.ipert1==dtset%natom+4)) need_ylmgr1=1
       dimylmgr1=max(useylmgr1,need_ylmgr1)
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*dimylmgr1))

!      Get plane-wave vectors and related data at k
       kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
       if (psps%useylm==1) then
         do jj=1,psps%mpsang*psps%mpsang
           ylm_k(1:npw_k,jj)=ylm(1+ikg:npw_k+ikg,jj)
         end do
       end if

!      Get plane-wave vectors and related data at k+q
       kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
       if (psps%useylm==1) then
         do jj=1,psps%mpsang*psps%mpsang
           ylm1_k(1:npw1_k,jj)=ylm1(1+ikg1:npw1_k+ikg1,jj)
         end do
         if (need_ylmgr1==1.and.useylmgr1/=0) then
           do jj=1,psps%mpsang*psps%mpsang
             do ia=1,3
               ylmgr1_k(1:npw1_k,ia,jj)=ylmgr1(1+ikg1:npw1_k+ikg1,ia,jj)
             end do
           end do
         end if
       end if

!      If Ylm gradients at k+q are needed and not in memory, compute them
       if (need_ylmgr1==1.and.useylmgr1==0) then
         option=-1;npwar1_tmp(1)=npw1_k;nband_tmp(1)=nband_k
         !Subtlety: initylmg is called in sequential mode
         call initylmg(gprimd,kg1_k,kpq,1,mpi_enreg_seq,psps%mpsang,&
&         npw1_k,nband_tmp,1,npwar1_tmp,nsppol,option,rprimd,&
&         ylm1_k,ylmgr1_k)
       end if

!      Compute (k+G) vectors
       nkpg=0;if(ipert1<=dtset%natom) nkpg=3*dtset%nloalg(3)
       ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
       if (nkpg>0) then
         call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
       end if

!      Compute (k+q+G) vectors
       nkpg1=0;if(ipert1<=dtset%natom.or.need_ylmgr1==1) nkpg1=3*dtset%nloalg(3)
       ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
       if (nkpg1>0) then
         call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)
       end if

!      Allocate kinetic contributions
       ABI_ALLOCATE(dkinpw,(npw_k))
       ABI_ALLOCATE(kinpw1,(npw1_k))
       dkinpw=zero
!      Compute (1/2) (2 Pi)**2 (k+q+G)**2:
!       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k)
       call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k,0,0)

!      Compute nonlocal form factors ffnl at (k+G), for all atoms
       ider=0;idir0=0
       dimffnl=0;if (ipert1<=dtset%natom) dimffnl=1
       ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
       if (ipert1<=dtset%natom) then
         call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,&
&         ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,&
&         psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,usepaw,&
&         psps%useylm,ylm_k,ylmgr_dum)
       end if

!      Compute nonlocal form factors ffnl1 at (k+q+G), for all atoms
       ider=1;if (ipert1<=dtset%natom) ider=0
       if(ipert1==dtset%natom+3.or.ipert1==dtset%natom+4)then
         dimffnl1=1;if (ider>=1) dimffnl1=2+5*psps%useylm
         idir0=0;if (ider>0.and.psps%useylm==1) idir0=-7
       else
         dimffnl1=1;if (ider>=1) dimffnl1=2+2*psps%useylm
         idir0=0;if (ider>0.and.psps%useylm==1) idir0=4
       end if
       ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
       call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,gs_hamkq%gprimd,&
&       ider,idir0,psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,&
&       psps%mqgrid_ff,nkpg1,npw1_k,psps%ntypat,psps%pspso,psps%qgrid_ff,rmet,usepaw,&
&       psps%useylm,ylm1_k,ylmgr1_k)

!      Extract non-local form factors for H^(j1)
       if (ipert1<=dtset%natom) then
         ffnl1_idir1 => ffnl1(:,:,:,:)
         dimffnl1_idir1=dimffnl1
       else
         dimffnl1_idir1=1+ider
         ABI_ALLOCATE(ffnl1_idir1,(npw1_k,dimffnl1_idir1,psps%lmnmax,psps%ntypat))
         ii=1;if (psps%useylm==0) ii=1+ider
         do itypat=1,psps%ntypat
           do ilmn=1,psps%lmnmax
             ffnl1_idir1(1:npw1_k,1:ii,ilmn,itypat)=ffnl1(1:npw1_k,1:ii,ilmn,itypat)
           end do
         end do
       end if

!      Load k-dependent part in the Hamiltonian datastructure
       ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
       call gs_hamkq%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
&       ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
       if (size(ffnlk)>0) then
         call gs_hamkq%load_k(ffnl_k=ffnlk)
       else
         call gs_hamkq%load_k(ffnl_k=ffnl1)
       end if

!      Load k+q-dependent part in the Hamiltonian datastructure
!          Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
       call gs_hamkq%load_kprime(kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
&       kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1_idir1,&
&       compute_gbound=.true.)
       if (qne0) then
         ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
         call gs_hamkq%load_kprime(ph3d_kp=ph3d1,compute_ph3d=.true.)
       end if

!      Load k-dependent part in the 1st-order Hamiltonian datastructure
       call rf_hamkq%load_k(npw_k=npw_k,dkinpw_k=dkinpw)

!      Allocate memory space for one band
       if (need_wfk)  then
         ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
       end if
       if (need_wf1)  then
         ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))
       end if
       ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
       nullify(cwaveprj0_idir1)
       if (usecprj==1) then
         ABI_DATATYPE_ALLOCATE(cwaveprj0,(dtset%natom,nspinor))
         call pawcprj_alloc(cwaveprj0,ncpgr,gs_hamkq%dimcprj)
       end if
       if (has_dcwf) then
         ABI_ALLOCATE(gs1,(2,npw1_k*nspinor))
       else
         ABI_ALLOCATE(gs1,(0,0))
       end if

!      LOOP OVER BANDS
       do iband=1,nband_k

!        Skip band if not to be treated by this proc
         if (xmpi_paral==1) then
           if (mpi_enreg%proc_distrb(ikpt,iband,isppol)/=me) cycle
         end if

!        Extract GS wavefunctions
         if (need_wfk) then
           cwave0(:,:)=cg(:,1+(iband-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
           if (usecprj==1) then
             call pawcprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,dtset%natom,iband,ibg,ikpt,iorder_cprj,&
&             isppol,dtset%mband,mkmem,dtset%natom,1,nband_k,nspinor,nsppol,dtfil%unpaw,&
&             mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
           end if
         end if

!        Extract 1st-order wavefunctions
         if (need_wf1) then
           cwavef(:,:)=cg1(:,1+(iband-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)
         end if

!        LOOP OVER PERTURBATION DIRECTIONS
         do kdir1=1,mdir1
           idir1=jdir1(kdir1)
           istr1=idir1;if(ipert1==dtset%natom+4) istr1=idir1+3

!          Not able to compute if ipert1=(Elect. field) and no ddk WF file
           if (ipert1==dtset%natom+2.and.ddkfil(idir1)==0) cycle

!          Extract 1st-order NL form factors derivatives for this idir1
           if (dimffnl1_idir1>=2.and.psps%useylm==1.and.ipert1>dtset%natom) then
             do itypat=1,psps%ntypat
               do ilmn=1,psps%lmnmax
                 ffnl1_idir1(1:npw1_k,2,ilmn,itypat)=ffnl1(1:npw1_k,1+istr1,ilmn,itypat)
               end do
             end do
           end if

!          Extract ground state projected WF and derivatives in idir1 direction
           if (usecprj==1) then
             cpopt=1

!            === Atomic displ. perturbation
             if (ipert<=dtset%natom) then
               ABI_DATATYPE_ALLOCATE(cwaveprj0_idir1,(dtset%natom,nspinor))
               call pawcprj_alloc(cwaveprj0_idir1,1,gs_hamkq%dimcprj)
               if (ipert1<=dtset%natom) then
                 call pawcprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=idir1)
               else
                 if (ipert1==dtset%natom+2) then
                   idir_cprj=idir1;choice=5
                 end if
                 if (ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
                   idir_cprj=istr1;choice=3
                 end if
                 call pawcprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=-1)
                 call getcprj(choice,cpopt,cwave0,cwaveprj0_idir1,&
&                 gs_hamkq%ffnl_kp,idir_cprj,gs_hamkq%indlmn,gs_hamkq%istwf_kp,&
&                 gs_hamkq%kg_kp,gs_hamkq%kpg_kp,gs_hamkq%kpt_kp,gs_hamkq%lmnmax,&
&                 gs_hamkq%mgfft,mpi_enreg,gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,&
&                 gs_hamkq%nloalg,gs_hamkq%npw_kp,gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkpxred,&
&                 gs_hamkq%ph1d,gs_hamkq%ph3d_kp,gs_hamkq%ucvol,gs_hamkq%useylm)
               end if

!            === Wave-vector perturbation
             else if (ipert==dtset%natom+1) then
               cwaveprj0_idir1 => cwaveprj0

!            == Electric field perturbation
             else if (ipert==dtset%natom+2) then
               ABI_DATATYPE_ALLOCATE(cwaveprj0_idir1,(dtset%natom,nspinor))
               call pawcprj_alloc(cwaveprj0_idir1,1,gs_hamkq%dimcprj)
               if (ipert1==dtset%natom+2) then
                 call pawcprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=idir1)
               else
                 if (ipert1<=dtset%natom) then
                   idir_cprj=idir1;choice=2
                 end if
                 if (ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
                   idir_cprj=istr1;choice=3
                 end if
                 call pawcprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=-1)
                 call getcprj(choice,cpopt,cwave0,cwaveprj0_idir1,&
&                 gs_hamkq%ffnl_kp,idir_cprj,gs_hamkq%indlmn,gs_hamkq%istwf_kp,&
&                 gs_hamkq%kg_kp,gs_hamkq%kpg_kp,gs_hamkq%kpt_kp,gs_hamkq%lmnmax,&
&                 gs_hamkq%mgfft,mpi_enreg,gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,gs_hamkq%nloalg,&
&                 gs_hamkq%npw_kp,gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkpxred,gs_hamkq%ph1d,&
&                 gs_hamkq%ph3d_kp,gs_hamkq%ucvol,gs_hamkq%useylm)
               end if

!            === Strain perturbation
             else if (ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
               if ((ipert1==dtset%natom+3.or.ipert1==dtset%natom+4).and.(istr==istr1)) then
                 cwaveprj0_idir1 => cwaveprj0
               else
                 ABI_DATATYPE_ALLOCATE(cwaveprj0_idir1,(dtset%natom,nspinor))
                 call pawcprj_alloc(cwaveprj0_idir1,ncpgr,gs_hamkq%dimcprj)
                 if (ipert1<=dtset%natom) then
                   idir_cprj=idir1;choice=2
                 end if
                 if (ipert1==dtset%natom+2) then
                   idir_cprj=idir1;choice=5
                 end if
                 if (ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
                   idir_cprj=istr1;choice=3
                 end if
                 call pawcprj_copy(cwaveprj0,cwaveprj0_idir1,icpgr=-1)
                 call getcprj(choice,cpopt,cwave0,cwaveprj0_idir1,&
&                 gs_hamkq%ffnl_kp,idir_cprj,gs_hamkq%indlmn,gs_hamkq%istwf_kp,gs_hamkq%kg_kp,&
&                 gs_hamkq%kpg_kp,gs_hamkq%kpt_kp,gs_hamkq%lmnmax,gs_hamkq%mgfft,mpi_enreg,&
&                 gs_hamkq%natom,gs_hamkq%nattyp,gs_hamkq%ngfft,gs_hamkq%nloalg,gs_hamkq%npw_kp,&
&                 gs_hamkq%nspinor,gs_hamkq%ntypat,gs_hamkq%phkpxred,gs_hamkq%ph1d,gs_hamkq%ph3d_kp,&
&                 gs_hamkq%ucvol,gs_hamkq%useylm)
               end if
             end if ! ipert

           else ! usecprj=0: cwaveprj0_idir1 is not used
             cwaveprj0_idir1 => cwaveprj0
           end if

!          Eventually compute 1st-order kinetic operator
           if (ipert1==dtset%natom+1) then
             call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg_k,dkinpw,kpoint,npw_k,idir1,0)
           else if (ipert1==dtset%natom+3.or.ipert1==dtset%natom+4) then
             call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,gprimd,istr1,kg_k,kpoint,npw_k)
           end if

!          Finalize initialization of 1st-order NL hamiltonian
           if (need_pawij10) rf_hamkq%e1kbfr => e1kbfr(:,:,:,:,idir1)

!          Read DDK wave function (if ipert1=electric field)
           if (ipert1==dtset%natom+2) then
             usevnl=1
             ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor*usevnl))
             if (need_ddk_file) then
               if (ddkfil(idir1)/=0) then
                 !ik_ddk = wfk_findk(ddks(idir1), kpt_rbz(:,ikpt)
                 ik_ddk = indkpt1(ikpt)
                 call ddks(idir1)%read_bks(iband, ik_ddk, isppol, xmpio_single, cg_bks=gvnlx1)
               else
                 gvnlx1=zero
               end if
               if (ipert1==dtset%natom+2) then
                 do ii=1,npw1_k*nspinor ! Multiply ddk by +i (to be consistent with getgh1c)
                   arg=gvnlx1(1,ii)
                   gvnlx1(1,ii)=-gvnlx1(2,ii)
                   gvnlx1(2,ii)=arg
                 end do
               end if
             else
               gvnlx1=zero
             end if
           else
             usevnl=0
             ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor*usevnl))
           end if

!          Get |H^(j2)-Eps_k_i.S^(j2)|u0_k_i> (VHxc-dependent part not taken into account) and S^(j2)|u0>
           lambda=eig_k(iband);berryopt=0;optlocal=0
           optnl=0;if (ipert1/=dtset%natom+1.or.idir==idir1) optnl=1
           opt_gvnlx1=0;if (ipert1==dtset%natom+2) opt_gvnlx1=2
           sij_opt=-1;if (has_dcwf) sij_opt=1
           if (usepaw==0) sij_opt=0
           call getgh1c(berryopt,cwave0,cwaveprj0_idir1,gh1,dum1,gs1,gs_hamkq,gvnlx1,idir1,ipert1,&
&           lambda,mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)
           if (sij_opt==1.and.optnl==1) gh1=gh1-lambda*gs1
           ABI_DEALLOCATE(gvnlx1)

!          If needed, compute here <delta_u^(j1)_k_i|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>
!          with delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!          (see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42))
!          This can be rewritten as:
!          -1/2.<u0_k_i|S^(j1)| Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
!          The sum over j can be computed with a single call to projbd routine
!          At first call (when j1=j2), ch1c=<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i> is stored
!          For the next calls, it is re-used.
           if (has_dcwf.or.(ipert==ipert1.and.idir==idir1.and.usepaw==1)) then
!            note: gvnlx1 used as temporary space
             ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor))
             if (ipert==ipert1.and.idir==idir1) then
               option=0;gvnlx1=gh1
             else
               option=1;gvnlx1=zero
             end if
!            Compute -Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
             call projbd(cgq,gvnlx1,-1,icgq,0,istwf_k,mcgq,0,nband_k,npw1_k,nspinor,&
&             dum1,ch1c(:,1:nband_k,iband,ikpt_me),option,tim_projbd,0,mpi_enreg%me_g0,mpi_enreg%comm_fft)
             if (has_dcwf) then
               if (ipert==ipert1.and.idir==idir1) gvnlx1=gvnlx1-gh1
               dotr=zero;doti=zero
               if (abs(occ_k(iband))>tol8) then
!                Compute: -<u0_k_i|S^(j1)| Sum_{j}[<u0_k+q_j|H^(j2)-Eps_k_i.S^(j2)|u0_k_i>.|u0_k+q_j>
                 call dotprod_g(dotr,doti,istwf_k,npw1_k*nspinor,2,gs1,gvnlx1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!                Add contribution to DDB
!                Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
!                (-1) factor already present
                 d2ovl_k(1,idir1)=d2ovl_k(1,idir1)+wtk_k*occ_k(iband)*elfd_fact*dotr
                 d2ovl_k(2,idir1)=d2ovl_k(2,idir1)+wtk_k*occ_k(iband)*elfd_fact*doti
               end if
             end if
             ABI_DEALLOCATE(gvnlx1)
           end if

!          If needed, compute here <delta_u^(j1)_k_i|H-Eps_k_i.S|u^(j2)_k_i>
!          This is equal to <delta_u^(j1)_k_i|H-Eps_k_i.S|delta_u^(j2)_k_i>  (I)
!                          +<delta_u^(j1)_k_i|H-Eps_k_i.S|u^paral^(j2)_k_i>  (II)
!          (u^paral^(j2)_k_i is the part of u^(j2)_k_i parallel to active space : metals)
!          (I) can be rewritten as:
!          Sum_j{ 1/4.<u0_k_i|S^(j1)|u0_k+q_j>.<u0_k+q_j|S^(j2)|u0_k_i>.(Eps_k+q_j-Eps_k_i) }
!          (II) can be rewritten as:
!          Sum_j{1/2.(occ_kq_j-occ_k_i).Eps1_k,q_ij.<u0_k_i|S^(j1)|u0_k+q_j> }
!          where Eps1_k,q_ij=<u0_k+q_j|H^(j2)-1/2(Eps_k+q_j-Eps_k_i)S^(j2)|u0_k_i>
!          At first call (when j1=j2), cs1c=<u0_k_i|S^(j1)|u0_k+q_j> is stored
!          For the next calls, it is re-used.
           if (has_dcwf.and.is_metal_or_qne0) then
             ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor))
             dotr=zero;doti=zero
             invocc=zero ; if (abs(occ_k(iband))>tol8) invocc=two/occ_k(iband)
             do jband=1,nband_k
!              Computation of cs1c=<u0_k_i|S^(j1)|u0_k+q_j>
               if ((ipert==ipert1.and.idir==idir1).or.(abs(occ_k(iband))>tol8)) then
                 gvnlx1(:,1:npw1_k*nspinor)=cgq(:,1+npw1_k*nspinor*(jband-1)+icgq:npw1_k*nspinor*jband+icgq)
                 call dotprod_g(dot1r,dot1i,istwf_k,npw1_k*nspinor,2,gs1,gvnlx1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
                 if (ipert==ipert1.and.idir==idir1.and.has_dcwf2) then
                   cs1c(1,jband,iband,ikpt_me)=dot1r
                   cs1c(2,jband,iband,ikpt_me)=dot1i
                 end if
               end if
               if (abs(occ_k(iband))>tol8) then
!                Computation of term (I)
                 if (has_dcwf2) then
                   arg=eig_kq(jband)-eig_k(iband)
                   dot2r=cs1c(1,jband,iband,ikpt_me)
                   dot2i=cs1c(2,jband,iband,ikpt_me)
                   dotr=dotr+(dot1r*dot2r+dot1i*dot2i)*arg
                   doti=doti+(dot1i*dot2r-dot1r*dot2i)*arg
                 end if
!                Computation of term (II)
                 if (is_metal) then
                   if (abs(rocceig(jband,iband))>tol8) then
                     ii=2*jband-1+(iband-1)*2*nband_k
                     arg=invocc*rocceig(jband,iband)*(eig_k(iband)-eig_kq(jband))
                     dot2r=eig1_k(ii);dot2i=eig1_k(ii+1)
                     dotr=dotr+arg*(dot1r*dot2r-dot1i*dot2i)
                     doti=doti+arg*(dot1r*dot2i+dot2r*dot1i)
                   end if
                 end if
               end if
             end do
             dotr=quarter*dotr;doti=quarter*doti
!            Note: factor 2 (from d2E/dj1dj2=2E^(j1j2))
             d2ovl_k(1,idir1)=d2ovl_k(1,idir1)+wtk_k*occ_k(iband)*two*elfd_fact*dotr
             d2ovl_k(2,idir1)=d2ovl_k(2,idir1)+wtk_k*occ_k(iband)*two*elfd_fact*doti
             ABI_DEALLOCATE(gvnlx1)
           end if

!          Build the matrix element <u0_k_i|H^(j1)-Eps_k_i.S^(j1)|u^(j2)_k,q_i>
!          and add contribution to DDB
           if (ipert1/=dtset%natom+1) then
             if (abs(occ_k(iband))>tol8) then
               call dotprod_g(dotr,doti,istwf_k,npw1_k*nspinor,2,gh1,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
!              Case ipert1=natom+2 (electric field):
!              gh1 contains H^(j1)|u0_k_i> (VHxc constant) which corresponds
!              to i.d/dk in Eq. (38) of Gonze, PRB 55, 10355 (1997) [[cite:Gonze1997a]].
!              * if ipert==natom+2, we apply directly Eq. (38)
!              * if ipert/=natom+2, Born effective charges are minus D2E
               d2nl_k(1,idir1)=d2nl_k(1,idir1)+wtk_k*occ_k(iband)*two*elfd_fact*dotr
               d2nl_k(2,idir1)=d2nl_k(2,idir1)+wtk_k*occ_k(iband)*two*elfd_fact*doti
             end if

!            Or compute localisation tensor (ddk)
!            See M. Veithen thesis Eq(2.5)
!            MT jan-2010: this is probably not correctly implemented for PAW !!!
!            missing terms due to S^(1) and S^(2)
           else
!            note: gh1 used as temporary space (to store idir ddk WF)
             if (idir==idir1) then
               gh1=cwavef
             else
               if (need_ddk_file.and.ddkfil(idir1)/=0) then
                 !ik_ddk = wfk_findk(ddks(idir1), kpt_rbz(:,ikpt)
                 ik_ddk = indkpt1(ikpt)
                 call ddks(idir1)%read_bks(iband, ik_ddk, isppol, xmpio_single, cg_bks=gh1)
               else
                 gh1=zero
               end if
             end if
             if (abs(occ_k(iband))>tol8) then
               call dotprod_g(dotr,doti,istwf_k,npw1_k*nspinor,2,gh1,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               call dotprod_g(dot1r,dot1i,istwf_k,npw1_k*nspinor,2,cwave0,gh1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               call dotprod_g(dot2r,dot2i,istwf_k,npw1_k*nspinor,2,cwavef,cwave0,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
               dotr=dotr-(dot1r*dot2r-dot1i*dot2i)
               doti=doti-(dot1r*dot2i+dot1i*dot2r)
               d2nl_k(1,idir1)=d2nl_k(1,idir1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
               d2nl_k(2,idir1)=d2nl_k(2,idir1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)
             end if
           end if

!          Accumulate here 1st-order density change due to overlap operator changes (if any)
           if (has_drho) then
             if (abs(occ_k(iband))>tol8) then
!              Compute here delta_u^(j1)=-1/2 Sum_{j}[<u0_k+q_j|S^(j1)|u0_k_i>.|u0_k+q_j>]
!              (see PRB 78, 035105 (2008) [[cite:Audouze2008]], Eq. (42))
               ABI_ALLOCATE(dcwavef,(2,npw1_k*nspinor))
               ABI_DATATYPE_ALLOCATE(dcwaveprj,(dtset%natom,nspinor))
               call pawcprj_alloc(dcwaveprj,0,gs_hamkq%dimcprj)
               call getdc1(cgq,cprjq,dcwavef,dcwaveprj,ibgq,icgq,istwf_k,mcgq,&
&               mcprjq,mpi_enreg,dtset%natom,nband_k,npw1_k,nspinor,1,gs1)

!              Accumulate 1st-order density due to delta_u^(j1)
               option=1;wfcorr=0
               call dfpt_accrho(cplex,cwave0,dcwavef,dcwavef,cwaveprj0_idir1,dcwaveprj,&
&               lambda,gs_hamkq,iband,idir1,ipert1,isppol,dtset%kptopt,&
&               mpi_enreg,dtset%natom,nband_k,1,npw_k,npw1_k,nspinor,occ_k,option,&
&               pawdrhoij1_unsym(:,idir1),drhoaug1(:,:,:,idir1),tim_fourwf,wfcorr,wtk_k)
               call pawcprj_free(dcwaveprj)
               ABI_DEALLOCATE(dcwavef)
               ABI_DATATYPE_DEALLOCATE(dcwaveprj)
             end if
           end if

           if((usecprj==1).and..not.(associated(cwaveprj0_idir1,cwaveprj0)))then
             call pawcprj_free(cwaveprj0_idir1)
             ABI_DATATYPE_DEALLOCATE(cwaveprj0_idir1)
           end if

!          End of loops
         end do   ! idir1
       end do     ! iband

!      Accumulate contribution of this k-point
       d2nl (:,:,ipert1,idir,ipert)=d2nl (:,:,ipert1,idir,ipert)+d2nl_k (:,:)
       if (usepaw==1) d2ovl(:,:,ipert1,idir,ipert)=d2ovl(:,:,ipert1,idir,ipert)+d2ovl_k(:,:)

!      Deallocations of arrays used for this k-point
       ABI_DEALLOCATE(gh1)
       ABI_DEALLOCATE(gs1)
       if (need_wfk)  then
         ABI_DEALLOCATE(cwave0)
       end if
       if (need_wf1)  then
         ABI_DEALLOCATE(cwavef)
       end if
       ABI_DEALLOCATE(kg_k)
       ABI_DEALLOCATE(kg1_k)
       ABI_DEALLOCATE(ylm_k)
       ABI_DEALLOCATE(ylm1_k)
       ABI_DEALLOCATE(ylmgr1_k)
       ABI_DEALLOCATE(kpg_k)
       ABI_DEALLOCATE(kpg1_k)
       ABI_DEALLOCATE(d2nl_k)
       ABI_DEALLOCATE(d2ovl_k)
       ABI_DEALLOCATE(eig_k)
       ABI_DEALLOCATE(eig_kq)
       ABI_DEALLOCATE(eig1_k)
       ABI_DEALLOCATE(occ_k)
       if (is_metal)  then
         ABI_DEALLOCATE(doccde_k)
         ABI_DEALLOCATE(doccde_kq)
         ABI_DEALLOCATE(occ_kq)
         ABI_DEALLOCATE(rocceig)
       end if
       ABI_DEALLOCATE(dkinpw)
       ABI_DEALLOCATE(kinpw1)
       ABI_DEALLOCATE(ph3d)
       if (allocated(ph3d1)) then
         ABI_DEALLOCATE(ph3d1)
       end if
       ABI_DEALLOCATE(ffnlk)
       ABI_DEALLOCATE(ffnl1)
       if (ipert1>dtset%natom)  then
         ABI_DEALLOCATE(ffnl1_idir1)
       end if
       nullify(ffnl1_idir1)
       if (usecprj==1) then
         call pawcprj_free(cwaveprj0)
         ABI_DATATYPE_DEALLOCATE(cwaveprj0)
       end if
       nullify(cwaveprj0_idir1)
!      Shift arrays
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
       if (mkmem/=0) then
         ibg=ibg+nspinor*nband_k
         icg=icg+npw_k*nspinor*nband_k
         ikg=ikg+npw_k
       end if
       if (mkqmem/=0) then
         ibgq=ibgq+nspinor*nband_k
         icgq=icgq+npw1_k*nspinor*nband_k
       end if
       if (mk1mem/=0) then
         ibg1=ibg1+nspinor*nband_k
         icg1=icg1+npw1_k*nspinor*nband_k
         ikg1=ikg1+npw1_k
       end if

     end do ! End loop over K-POINTS

!    Transfer 1st-order density change due to overlap; also take into account the spin.
     if(has_drho) then
       do kdir1=1,mdir1
         idir1=jdir1(kdir1)
         call fftpac(isppol,mpi_enreg,nspden,cplex*dtset%ngfft(1),dtset%ngfft(2),dtset%ngfft(3),&
&         cplex*dtset%ngfft(4),dtset%ngfft(5),dtset%ngfft(6),&
&         dtset%ngfft,drho1wfr(:,:,idir1),drhoaug1(:,:,:,idir1),1)
       end do
     end if

   end do ! End loop over SPINS

!  Free memory used for this type of perturbation
   call rf_hamkq%free()
   if (has_drho)  then
     ABI_DEALLOCATE(drhoaug1)
   end if
   if (need_pawij10) then
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       call paw_ij_free(paw_ij10(:,idir1))
     end do
     ABI_DEALLOCATE(e1kbfr_spin)
   end if
   ABI_DATATYPE_DEALLOCATE(paw_ij10)

!  In case of parallelism, sum 1st-order density and occupation matrix over processors
   if (has_drho.and.xmpi_paral==1) then

!    Accumulate 1st-order density
     call timab(48,1,tsec)
     bufsz=cplex*dtset%nfft*nspden*mdir1
     ABI_ALLOCATE(buffer,(bufsz))
     buffer(1:bufsz)=reshape(drho1wfr,(/bufsz/))
     call xmpi_sum(buffer,bufsz,spaceworld,ierr)
     drho1wfr(:,:,:)=reshape(buffer(1:bufsz),(/cplex*dtset%nfft,nspden,mdir1/))
     ABI_DEALLOCATE(buffer)
     call timab(48,2,tsec)

!    Accumulate 1st-order PAW occupancies
     if (usepaw==1) then
       call pawrhoij_mpisum_unpacked(pawdrhoij1_unsym,spaceworld)
     end if

   end if

!  Compute second part of overlap contribution (due to VHxc^(j2)(tild_n+hat_n))
   if (has_drho) then

     ABI_ALLOCATE(drhor1,(cplex*nfftf,nspden))
     ABI_ALLOCATE(dnhat1,(cplex*nfftf,nspden))

!    LOOP OVER PERTURBATION DIRECTIONS
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)

!      Build and symmetrize 1st-order density change due to change of overlap
       ABI_ALLOCATE(drho1wfg,(2,dtset%nfft))
       call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,&
&       nspden,nsppol,nsym1,phnons1,drho1wfg,drho1wfr(:,:,idir1),&
&       rprimd,symaf1,symrl1,tnons1)
       if (dtset%pawstgylm/=0) then
         option=0
         call pawnhatfr(option,idir1,ipert1,my_natom,dtset%natom,nspden,dtset%ntypat,&
&         pawang,pawfgrtab,pawrhoij,pawtab,rprimd,&
&         mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)
       end if
       call pawmkrho(1,arg,cplex,gprimd,idir1,indsy1,ipert1,&
&       mpi_enreg,my_natom,dtset%natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,&
&       pawfgr,pawfgrtab,-10001,pawdrhoij1(:,idir1),pawdrhoij1_unsym(:,idir1),pawtab,&
&       dtset%qptn,drho1wfg,drho1wfr(:,:,idir1),drhor1,rprimd,symaf1,symrc1,dtset%typat,&
&       ucvol,dtset%usewvl,xred,pawang_sym=pawang1,pawnhat=dnhat1,pawrhoij0=pawrhoij)
       ABI_DEALLOCATE(drho1wfg)

!      Compute plane-wave contribution to overlap contribution
!      This is subtle as it is a mix of Eq(79) and Eq(80) of PRB 78, 035105 (2008) [[cite:Audouze2008]]
!      Details:
!      The VH(tild_nZc)^(1) term of Eq(79) is:
!      <VH(tild_nZc)^(j2)|delta_tild_rho^(j1)>            = <vpsp1|drhor1-dnhat1>
!      The first term of Eq(80) is:
!      <VHxc^(j2)|delta_tild_rho^(j1)+delta_hat_rho^(j1)> = <vtrial1-vpsp1|drhor1>
!      The addition of these two terms gives:
!      <vtrial1|drhor1>-<vpsp1|dnhat1>
!      And this is more subtle when usexcnhat=0
       call dotprod_vn(cplex,drhor1,dot1r,dot1i,nfftf,nfftot,nspden,2,vtrial1,ucvol)
       if (usexcnhat/=0) then
         call dotprod_vn(cplex,dnhat1,dot2r,dot2i,nfftf,nfftot,1   ,2,vpsp1,ucvol)
       else
         ABI_ALLOCATE(vtmp1,(cplex*nfftf,nspden))
         do ispden=1,nspden
           vtmp1(:,ispden)=vtrial1(:,ispden)-vhartr1(:)
         end do
         call dotprod_vn(cplex,dnhat1,dot2r,dot2i,nfftf,nfftot,nspden,2,vtmp1,ucvol)
         ABI_DEALLOCATE(vtmp1)
       end if
       dotr=dot1r-dot2r;doti=dot1i-dot2i

!      Compute on-site contributions to overlap contribution
!      (two last terms of Eq(80) of PRB 78, 035105 (2008)) [[cite:Audouze2008]]
!      (note: Dij^(j2) and Vxc^(j2) are computed for ipert at first call)
       call pawdfptenergy(epawnst,ipert,ipert1,dtset%ixc,my_natom,dtset%natom,dtset%ntypat,&
&       nzlmopt_ipert,nzlmopt_ipert1,paw_an,paw_an1,paw_ij1,pawang,dtset%pawprtvol,&
&       pawrad,pawrhoij1,pawdrhoij1(:,idir1),pawtab,dtset%pawxcdev,dtset%xclevel,&
&       mpi_atmtab=my_atmtab,comm_atom=my_comm_atom)

!      Accumulate in 2nd-order matrix:
!      Note: factor 2 (from d2E/dj1dj2=2E^(j1j2)) eliminated by factor 1/2
!      has to take the complex conjugate because we want here Int[VHxc^(j1)^*.delta_rho^(j2)]
       dotr=dotr+epawnst(1);doti=-(doti+epawnst(2))
       d2ovl_drho(1,idir1,ipert1,idir,ipert)=elfd_fact*dotr
       d2ovl_drho(2,idir1,ipert1,idir,ipert)=elfd_fact*doti

     end do ! End loop over perturbation directions

!    Free no more needed memory
     ABI_DEALLOCATE(drhor1)
     ABI_DEALLOCATE(dnhat1)
     ABI_DEALLOCATE(drho1wfr)
     do iatom=1,my_natom
       if (pawfgrtab(iatom)%nhatfr_allocated>0)  then
         ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
       end if
       pawfgrtab(iatom)%nhatfr_allocated=0
     end do
     if (paral_atom) then
       do kdir1=1,mdir1
         idir1=jdir1(kdir1)
         call pawrhoij_free(pawdrhoij1_unsym(:,idir1))
       end do
       ABI_DATATYPE_DEALLOCATE(pawdrhoij1_unsym)
     end if
     do kdir1=1,mdir1
       idir1=jdir1(kdir1)
       call pawrhoij_free(pawdrhoij1(:,idir1))
     end do
     ABI_DATATYPE_DEALLOCATE(pawdrhoij1)
   end if ! has_drho

!  End loop over perturbations (j1)
 end do

!Final deallocations
 ABI_DEALLOCATE(ch1c)
 if (usepaw==1.and.is_metal_or_qne0) then
   ABI_DEALLOCATE(cs1c)
 end if
 call gs_hamkq%free()

!In case of parallelism, sum over processors
 if (xmpi_paral==1)then
   call timab(161,1,tsec)
   call xmpi_barrier(spaceworld)
   call timab(161,2,tsec)
   ABI_ALLOCATE(buffer,(6*mpert*(1+usepaw)))
   buffer(1:6*mpert)=reshape(d2nl(:,:,:,idir,ipert),(/6*mpert/))
   if (usepaw==1) buffer(6*mpert+1:6*mpert+6*mpert)=reshape(d2ovl(:,:,:,idir,ipert),(/6*mpert/))
   call timab(48,1,tsec)
   call xmpi_sum(buffer,6*mpert*(1+usepaw),spaceworld,ierr)
   call timab(48,2,tsec)
   d2nl (:,:,:,idir,ipert)=reshape(buffer(1:6*mpert),(/2,3,mpert/))
   if (usepaw==1) d2ovl(:,:,:,idir,ipert)=reshape(buffer(6*mpert+1:6*mpert+6*mpert),(/2,3,mpert/))
   ABI_DEALLOCATE(buffer)
 end if

!Build complete d2ovl matrix
 if (usepaw==1) d2ovl(:,:,:,idir,ipert)=d2ovl(:,:,:,idir,ipert)+d2ovl_drho(:,:,:,idir,ipert)

 if (usepaw==1) then
   ABI_DEALLOCATE(d2ovl_drho)
 end if

!Close the ddk WF files
 if (has_ddk_file) then
   do kdir1=1,mdir1
     idir1=jdir1(kdir1)
     if (ddkfil(idir1)/=0)then
       call ddks(idir1)%close()
     end if
   end do
 end if
 ABI_DEALLOCATE(jpert1)
 ABI_DEALLOCATE(jdir1)

!Symmetrize the phonons contributions, as was needed for the forces in a GS calculation
 ABI_ALLOCATE(work,(2,3,dtset%natom))
 do ipert1=1,dtset%natom
   do idir1=1,3
     work(:,idir1,ipert1)=d2nl(:,idir1,ipert1,idir,ipert)
   end do
 end do
 call dfpt_sygra(dtset%natom,d2nl(:,:,:,idir,ipert),work,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 if (usepaw==1) then
   do ipert1=1,dtset%natom
     do idir1=1,3
       work(:,idir1,ipert1)=d2ovl(:,idir1,ipert1,idir,ipert)
     end do
   end do
   call dfpt_sygra(dtset%natom,d2ovl(:,:,:,idir,ipert),work,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 end if
 ABI_DEALLOCATE(work)

!In the case of the strain perturbation time-reversal symmetry will always
!be true so imaginary part of d2nl will be must be set to zero here since
!the symmetry-reduced kpt set will leave a non-zero imaginary part.
 if(ipert==dtset%natom+3.or.ipert==dtset%natom+4) then
   d2nl(2,:,:,idir,ipert)=zero
   if (usepaw==1) d2ovl(2,:,:,idir,ipert)=zero
 else
   d2nl(2,:,dtset%natom+3:dtset%natom+4,idir,ipert)=zero
   if (usepaw==1) d2ovl(2,:,dtset%natom+3:dtset%natom+4,idir,ipert)=zero
 end if


!Symmetrize the strain perturbation contributions, as was needed for the stresses in a GS calculation
!if (ipert==dtset%natom+3.or.ipert==dtset%natom+4)then
 if (nsym1>1) then
   ABI_ALLOCATE(work,(6,1,1))
   ii=0
   do ipert1=dtset%natom+3,dtset%natom+4
     do idir1=1,3
       ii=ii+1
       work(ii,1,1)=d2nl(1,idir1,ipert1,idir,ipert)
     end do
   end do
   call stresssym(gprimd,nsym1,work(:,1,1),symrc1)
   ii=0
   do ipert1=dtset%natom+3,dtset%natom+4
     do idir1=1,3
       ii=ii+1
       d2nl(1,idir1,ipert1,idir,ipert)=work(ii,1,1)
     end do
   end do
   if (usepaw==1) then
     ii=0
     do ipert1=dtset%natom+3,dtset%natom+4
       do idir1=1,3
         ii=ii+1
         work(ii,1,1)=d2ovl(1,idir1,ipert1,idir,ipert)
       end do
     end do
     call stresssym(gprimd,nsym1,work(:,1,1),symrc1)
     ii=0
     do ipert1=dtset%natom+3,dtset%natom+4
       do idir1=1,3
         ii=ii+1
         d2ovl(1,idir1,ipert1,idir,ipert)=work(ii,1,1)
       end do
     end do
     ABI_DEALLOCATE(work)
   end if
 end if
!end if

!Must also symmetrize the electric field perturbation response !
!Note: d2ovl is not symetrized because it is zero for electric field perturbation
 if (has_ddk_file) then
   ABI_ALLOCATE(d2nl_elfd,(2,3))
!  There should not be any imaginary part, but stay general (for debugging)
   d2nl_elfd (:,:)=d2nl(:,:,dtset%natom+2,idir,ipert)
   do ii=1,3
     sumelfd(:)=zero
     do ia=1,nsym1
       do jj=1,3
         if(symrl1(ii,jj,ia)/=0)then
           if(ddkfil(jj)==0)then
             blkflg(ii,dtset%natom+2,idir,ipert)=0
           end if
         end if
       end do
       symfact(1)=dble(symrl1(ii,1,ia))
       symfact(2)=dble(symrl1(ii,2,ia))
       symfact(3)=dble(symrl1(ii,3,ia))
       sumelfd(:)=sumelfd(:)+symfact(1)*d2nl_elfd(:,1) &
&       +symfact(2)*d2nl_elfd(:,2)+symfact(3)*d2nl_elfd(:,3)
     end do
     d2nl(:,ii,dtset%natom+2,idir,ipert)=sumelfd(:)/dble(nsym1)
   end do
   ABI_DEALLOCATE(d2nl_elfd)
 end if

!Overlap: store the diagonal part of the matrix in the
!         2nd-order energy non-stationnary expression
 eovl1=zero;if (usepaw==1) eovl1=d2ovl(1,idir,ipert,idir,ipert)

 call destroy_mpi_enreg(mpi_enreg_seq)
 call timab(566,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_nstpaw
!!***

!!****f* ABINIT/dfpt_nstwf
!! NAME
!! dfpt_nstwf
!!
!! FUNCTION
!! This routine computes the non-local contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!! Only for norm-conserving pseudopotentials (no PAW)
!!
!! COPYRIGHT
!! Copyright (C) 1999-2019 ABINIT group (XG,AR,MB,MVer,MT, MVeithen)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  ddkfil(3)=unit numbers for the three possible ddk files for ipert1
!!       equal to 0 if no dot file is available for this direction
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=parameter that describes the storage of wfs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpt(3)=reduced coordinates of k point
!!  kpq(3)=reduced coordinates of k+q point
!!  mkmem =number of k points treated by this node
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  ddks(3)<wfk_t>=struct info for for the three possible DDK files for ipert1
!!  wtk_k=weight assigned to the k point.
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!
!! OUTPUT
!!  d2bbb_k(2,3,mband,mband*prtbbb)=band by band decomposition of the second
!!   order derivatives, for the present k point, and perturbation idir, ipert
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! TODO
!!  XG 20141103 The localization tensor cannot be defined in the metallic case. It should not be computed.
!!
!! PARENTS
!!      dfpt_nstdy
!!
!! CHILDREN
!!      dotprod_g,gaugetransfo,getgh1c
!!      init_rf_hamiltonian
!!      mkffnl,mkkpg,timab,wfk_read_bks
!!
!! SOURCE

subroutine dfpt_nstwf(cg,cg1,ddkfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&
&                 icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpt,kpq,&
&                 mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&
&                 occ_k,psps,rmet,ddks,wtk_k,ylm,ylm1)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,idir,ikpt,ipert,isppol,istwf_k
 integer,intent(in) :: mkmem,mk1mem,mpert,mpw,mpw1,nsppol
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ddkfil(3),kg1_k(3,npw1_k)
 integer,intent(in) :: kg_k(3,npw_k)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband_mem*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*nsppol)
 real(dp),intent(in) :: eig_k(dtset%mband*nsppol),kpt(3),kpq(3),occ_k(nband_k),rmet(3,3)
 real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(npw1_k,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: eig1_k(2*nsppol*dtset%mband**2)
 real(dp),intent(out) :: d2bbb_k(2,3,dtset%mband,dtset%mband*dtset%prtbbb)
 real(dp),intent(inout) :: d2nl_k(2,3,mpert)
 type(wfk_t),intent(inout) :: ddks(3)

!Local variables-------------------------------
!scalars
 integer :: berryopt,dimffnl,dimffnl1,dimph3d
 integer :: iband,ider,idir1,ipert1,ipw,jband,nband_kocc,nkpg,nkpg1
 integer :: ierr, iband_me, jband_me
 integer :: npw_disk,nsp,optlocal,optnl,opt_gvnlx1,sij_opt,tim_getgh1c,usevnl
 logical :: ddk
 real(dp) :: aa,dot1i,dot1r,dot2i,dot2r,dot_ndiagi,dot_ndiagr,doti,dotr,lambda
 character(len=500) :: msg
 type(rf_hamiltonian_type) :: rf_hamkq
!arrays
 integer :: ik_ddks(3)
 integer :: band_procs(nband_k)
 logical :: distrb_cycle(nband_k)
 real(dp) :: dum_grad_berry(1,1),dum_gvnlx1(1,1),dum_gs1(1,1),dum_ylmgr(1,3,1),tsec(2)
 real(dp),allocatable :: cg_k(:,:),cwave0(:,:),cwavef(:,:),cwavef_da(:,:)
 real(dp),allocatable :: cwavef_db(:,:),dkinpw(:),eig2_k(:),ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: gvnlx1(:,:),kinpw1(:),kpg1_k(:,:),kpg_k(:,:),ph3d(:,:,:)
 type(pawcprj_type),allocatable :: dum_cwaveprj(:,:)

! *********************************************************************

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='  This routine cannot be used for PAW (use pawnst3 instead) !'
   MSG_BUG(msg)
 end if

!Keep track of total time spent in dfpt_nstwf
 call timab(102,1,tsec)
 tim_getgh1c=2

!Miscelaneous inits
 ABI_DATATYPE_ALLOCATE(dum_cwaveprj,(0,0))
 ddk=(ipert==dtset%natom+1.or.ipert==dtset%natom+10.or.ipert==dtset%natom+11)

! filter for bands on this cpu for cg cg1 etc.
 distrb_cycle = (mpi_enreg%proc_distrb(ikpt,1:nband_k,isppol) /= mpi_enreg%me_kpt)
 call proc_distrb_band(band_procs,mpi_enreg%proc_distrb,ikpt,isppol,nband_k,&
&  mpi_enreg%me_band,mpi_enreg%me_kpt,mpi_enreg%comm_band)


!Additional allocations
 if (.not.ddk) then
   ABI_ALLOCATE(dkinpw,(npw_k))
   ABI_ALLOCATE(kinpw1,(npw1_k))
   kinpw1=zero;dkinpw=zero
 else
   ABI_ALLOCATE(dkinpw,(0))
   ABI_ALLOCATE(kinpw1,(0))
 end if
 ABI_ALLOCATE(gvnlx1,(2,npw1_k*dtset%nspinor))
 ABI_ALLOCATE(eig2_k,(2*nsppol*dtset%mband**2))
 ABI_ALLOCATE(cwave0,(2,npw_k*dtset%nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*dtset%nspinor))

!Compute (k+G) vectors
 nkpg=0;if (.not.ddk) nkpg=3*gs_hamkq%nloalg(3)
 ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
 if (nkpg>0) then
   call mkkpg(kg_k,kpg_k,kpt,nkpg,npw_k)
 end if

!Compute (k+q+G) vectors
 nkpg1=0;if (.not.ddk) nkpg1=3*gs_hamkq%nloalg(3)
 ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
 if (nkpg1>0) then
   call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)
 end if

!Compute nonlocal form factors ffnl at (k+G)
 dimffnl=0;if (.not.ddk) dimffnl=1
 ABI_ALLOCATE(ffnlk,(npw_k,dimffnl,psps%lmnmax,psps%ntypat))
 if (.not.ddk) then
   ider=0
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnlk,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg_k,kpg_k,kpt,psps%lmnmax,&
&   psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,psps%ntypat,psps%pspso,&
&   psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,dum_ylmgr)
 end if

!Compute nonlocal form factors ffnl1 at (k+q+G)
 dimffnl1=0;if (.not.ddk) dimffnl1=1
 ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,psps%ntypat))
 if (.not.ddk) then
   ider=0
   call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gs_hamkq%gmet,&
&   gs_hamkq%gprimd,ider,ider,psps%indlmn,kg1_k,kpg1_k,kpq,&
&   psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,npw1_k,psps%ntypat,&
&   psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1,dum_ylmgr)
 end if

!Load k-dependent part in the Hamiltonian datastructure
 call gs_hamkq%load_k(kpt_k=kpt,npw_k=npw_k,istwf_k=istwf_k,&
& kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnlk)

!Load k+q-dependent part in the Hamiltonian datastructure
!    Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
 dimph3d=0;if (.not.ddk) dimph3d=gs_hamkq%matblk
 ABI_ALLOCATE(ph3d,(2,npw1_k,dimph3d))
 call gs_hamkq%load_kprime(kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
& kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
& ph3d_kp=ph3d,compute_ph3d=(.not.ddk))

!Load k-dependent part in the 1st-order Hamiltonian datastructure
 call rf_hamkq%load_k(npw_k=npw_k,dkinpw_k=dkinpw)

!Take care of the npw and kg records
!NOTE : one should be able to modify the rwwf routine to take care
!of the band parallelism, which is not the case yet ...
 ik_ddks = 0
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
!    Read npw record
     nsp=dtset%nspinor
     ik_ddks(idir1) = ddks(idir1)%findk(kpt)
     ABI_CHECK(ik_ddks(idir1) /= -1, "Cannot find kpt")
     npw_disk = ddks(idir1)%hdr%npwarr(ik_ddks(idir1))
     if (npw_k /= npw_disk) then
       write(unit=msg,fmt='(a,i3,a,i5,a,i3,a,a,i5,a,a,i5)')&
&       'For isppol = ',isppol,', ikpt = ',ikpt,' and idir = ',idir,ch10,&
&       'the number of plane waves in the ddk file is equal to', npw_disk,ch10,&
&       'while it should be ',npw_k
       MSG_BUG(msg)
     end if
   end if
 end do

 if (ipert==dtset%natom+1) then
   nband_kocc = 0
   do iband = 1,nband_k
     if (abs(occ_k(iband)) > tol8) nband_kocc = nband_kocc + 1
     nband_kocc = max (nband_kocc, 1)
   end do
 end if

 if(dtset%prtbbb==1)then
   ABI_ALLOCATE(cwavef_da,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cwavef_db,(2,npw1_k*dtset%nspinor))
   ABI_ALLOCATE(cg_k,(2,npw_k*dtset%nspinor*dtset%mband_mem))
   if ((ipert == dtset%natom + 1).or.(ipert <= dtset%natom).or. &
&   (ipert == dtset%natom + 2).or.(ipert == dtset%natom + 5)) then
     cg_k(:,:) = cg(:,1+icg:icg+dtset%mband_mem*npw_k*dtset%nspinor)
   end if
   d2bbb_k(:,:,:,:) = zero
 end if

!Loop over bands
 iband_me = 0
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt,iband,isppol) == mpi_enreg%me_kpt) then
     iband_me = iband_me + 1

!  Read ground-state wavefunction for iband
     if (dtset%prtbbb==0 .or. ipert==dtset%natom+2) then
       cwave0(:,:)=cg(:,1+(iband_me-1)*npw_k*dtset%nspinor+icg:iband_me*npw_k*dtset%nspinor+icg)
     else    ! prtbbb==1 and ipert<=natom , already in cg_k
       cwave0(:,:)=cg_k(:,1+(iband_me-1)*npw_k*dtset%nspinor:iband_me*npw_k*dtset%nspinor)
     end if

!  Get first-order wavefunctions for iband
     cwavef(:,:)=cg1(:,1+(iband_me-1)*npw1_k*dtset%nspinor+icg1:iband_me*npw1_k*dtset%nspinor+icg1)

   end if 
   call xmpi_bcast(cwave0, band_procs(iband), mpi_enreg%comm_band, ierr)
   call xmpi_bcast(cwavef, band_procs(iband), mpi_enreg%comm_band, ierr)

!  In case non ddk perturbation
   if (ipert /= dtset%natom + 1) then

     do ipert1=1,mpert

       if( ipert1<=dtset%natom .or. ipert1==dtset%natom+2 )then

!        Initialize data for NL 1st-order hamiltonian
         call init_rf_hamiltonian(1,gs_hamkq,ipert1,rf_hamkq)

         if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)) &
&         .and.(ipert1 == dtset%natom+2).and. dtset%prtbbb==1) then
           call gaugetransfo(cg_k,cwavef,cwavef_db,mpi_enreg%comm_band,distrb_cycle,eig_k,eig1_k,iband,nband_k, &
&            dtset%mband,dtset%mband_mem,npw_k,npw1_k,dtset%nspinor,nsppol,mpi_enreg%nproc_band,occ_k)
           cwavef(:,:) = cwavef_db(:,:)
         end if

!        Define the direction along which to move the atom :
!        the polarisation (ipert1,idir1) is refered as j1.
         do idir1=1,3
           if (ipert1<=dtset%natom.or.(ipert1==dtset%natom+2.and.ddkfil(idir1)/=0)) then

!            Get |Vnon-locj^(1)|u0> :
!            First-order non-local, applied to zero-order wavefunction
!            This routine gives MINUS the non-local contribution

!            ==== Atomic displ. perturbation
             if( ipert1<=dtset%natom )then
               lambda=eig_k((isppol-1)*nband_k+iband)
               berryopt=1;optlocal=0;optnl=1;usevnl=0;opt_gvnlx1=0;sij_opt=0
               call getgh1c(berryopt,cwave0,dum_cwaveprj,gvnlx1,dum_grad_berry,&
&               dum_gs1,gs_hamkq,dum_gvnlx1,idir1,ipert1,lambda,mpi_enreg,optlocal,&
&               optnl,opt_gvnlx1,rf_hamkq,sij_opt,tim_getgh1c,usevnl)

!              ==== Electric field perturbation
             else if( ipert1==dtset%natom+2 )then
               ! TODO: Several tests fail here ifdef HAVE_MPI_IO_DEFAULT
               ! The problem is somehow related to the use of MPI-IO file views!.
!TODO MJV: this needs to be band parallelized as well, probably
               call ddks(idir1)%read_bks(iband, ik_ddks(idir1), isppol, xmpio_single, cg_bks=gvnlx1, &
               eig1_bks=eig2_k(1+(iband-1)*2*nband_k:))
                 !eig1_bks=eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k))
               !write(777,*)"eig2_k, gvnlx1 for band: ",iband,", ikpt: ",ikpt
               !do ii=1,2*nband_k
               !  write(777,*)eig2_k(ii+(iband-1))
               !end do
               !write(777,*)gvnlx1

!              In case of band-by-band,
!              construct the first-order wavefunctions in the diagonal gauge
               if (((ipert <= dtset%natom).or.(ipert == dtset%natom + 2)).and.(dtset%prtbbb==1)) then
                 call gaugetransfo(cg_k,gvnlx1,cwavef_da,mpi_enreg%comm_band,distrb_cycle,eig_k,eig2_k,iband,nband_k, &
&                  dtset%mband,dtset%mband_mem,npw_k,npw1_k,dtset%nspinor,nsppol,mpi_enreg%nproc_band,occ_k)
                 gvnlx1(:,:) = cwavef_da(:,:)
               end if
!              Multiplication by -i
               do ipw=1,npw1_k*dtset%nspinor
                 aa=gvnlx1(1,ipw)
                 gvnlx1(1,ipw)=gvnlx1(2,ipw)
                 gvnlx1(2,ipw)=-aa
               end do
             end if

!            MVeithen 021212 :
!            1) Case ipert1 = natom + 2 and ipert = natom + 2:
!            the second derivative of the energy with respect to an electric
!            field is computed from Eq. (38) of X. Gonze, PRB 55 ,10355 (1997) [[cite:Gonze1997a]].
!            The evaluation of this formula needs the operator $i \frac{d}{dk}.
!            2) Case ipert1 = natom + 2 and ipert < natom:
!            the computation of the Born effective charge tensor uses
!            the operator $-i \frac{d}{dk}.
             if (ipert==dtset%natom+2) gvnlx1(:,:) = -gvnlx1(:,:)

!            <G|Vnl1|Cnk> is contained in gvnlx1
!            construct the matrix element (<uj2|vj1|u0>)complex conjug and add it to the 2nd-order matrix
             call dotprod_g(dotr,doti,istwf_k,npw1_k*dtset%nspinor,2,cwavef,gvnlx1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*two*dotr
             d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*two*doti

!            Band by band decomposition of the Born effective charges
!            calculated from a phonon perturbation
!            d2bbb_k will be mpisummed below so only keep my iband indices on the diagonal
             if(dtset%prtbbb==1 .and. mpi_enreg%proc_distrb(ikpt,iband,isppol) == mpi_enreg%me_kpt)then
               d2bbb_k(1,idir1,iband,iband) =      wtk_k*occ_k(iband)*two*dotr
               d2bbb_k(2,idir1,iband,iband) = -one*wtk_k*occ_k(iband)*two*doti
             end if

           end if
         end do ! idir

         call rf_hamkq%free()
       end if   
     end do     ! ipert1
   end if     ! ipert /= natom +1

!  Compute the localization tensor

   if (ipert==dtset%natom+1) then

     ipert1=dtset%natom+1
     if(dtset%prtbbb==1)then
       call gaugetransfo(cg_k,cwavef,cwavef_db,mpi_enreg%comm_band,distrb_cycle,eig_k,eig1_k,iband,nband_k, &
&       dtset%mband,dtset%mband_mem,npw_k,npw1_k,dtset%nspinor,nsppol,mpi_enreg%nproc_band,occ_k)
       cwavef(:,:) = cwavef_db(:,:)
     end if

     do idir1 = 1,3
       eig2_k(:) = zero
       gvnlx1(:,:) = zero
       if (idir == idir1) then
         gvnlx1(:,:) = cwavef(:,:)
         eig2_k(:) = eig1_k(:)
       else
         if (ddkfil(idir1) /= 0) then
           call ddks(idir1)%read_bks(iband, ik_ddks(idir1), isppol, xmpio_single, cg_bks=gvnlx1, &
           eig1_bks=eig2_k(1+(iband-1)*2*nband_k:))
             !eig1_bks=eig2_k(1+(iband-1)*2*nband_k:2*iband*nband_k))
           !write(778,*)"eig2_k, gvnlx1 for band: ",iband,", ikpt: ",ikpt
           !do ii=1,2*nband_k
           !  write(778,*)eig2_k(ii+(iband-1))
           !end do
           !write(778,*)gvnlx1

           if(dtset%prtbbb==1)then
             call gaugetransfo(cg_k,gvnlx1,cwavef_da,mpi_enreg%comm_band,distrb_cycle,eig_k,eig2_k,iband,nband_k, &
&             dtset%mband,dtset%mband_mem,npw_k,npw1_k,dtset%nspinor,nsppol,mpi_enreg%nproc_band,occ_k)

             gvnlx1(:,:) = cwavef_da(:,:)
           end if

         end if    !ddkfil(idir1)
       end if    !idir == idir1

!      <G|du/dqa> is contained in gvnlx1 and <G|du/dqb> in cwavef
!      construct the matrix elements <du/dqa|du/dqb> -> dot
!      <u|du/dqa> -> dot1
!      <du/dqb|u> -> dot2
!      and add them to the 2nd-order matrix

       call dotprod_g(dotr,doti,istwf_k,npw1_k*dtset%nspinor,2,gvnlx1,cwavef,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
       d2nl_k(1,idir1,ipert1)=d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*dotr/(nband_kocc*two)
       d2nl_k(2,idir1,ipert1)=d2nl_k(2,idir1,ipert1)+wtk_k*occ_k(iband)*doti/(nband_kocc*two)


!      XG 020216 : Marek, could you check the next forty lines
!      In the parallel gauge, dot1 and dot2 vanishes
       if(dtset%prtbbb==1)then
         ! d2bbb_k will be mpisummed below - only save my local band indices for diagonal contribution
         if (mpi_enreg%proc_distrb(ikpt,iband,isppol) == mpi_enreg%me_kpt) then
           d2bbb_k(1,idir1,iband,iband)=d2bbb_k(1,idir1,iband,iband)+dotr
           d2bbb_k(2,idir1,iband,iband)=d2bbb_k(2,idir1,iband,iband)+doti
         end if
         dot_ndiagr=zero ; dot_ndiagi=zero
         jband_me = 0
         do jband = 1,nband_k              !compute dot1 and dot2
           if (mpi_enreg%proc_distrb(ikpt,iband,isppol) /= mpi_enreg%me_kpt) cycle
           jband_me = jband_me + 1

           if (abs(occ_k(jband)) > tol8) then
             dot1r=zero ; dot1i=zero
             dot2r=zero ; dot2i=zero
             cwave0(:,:)=cg_k(:,1+(jband_me-1)*npw_k*dtset%nspinor:jband_me*npw_k*dtset%nspinor)

             call dotprod_g(dot1r,dot1i,istwf_k,npw1_k*dtset%nspinor,2,cwave0,gvnlx1,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)
             call dotprod_g(dot2r,dot2i,istwf_k,npw1_k*dtset%nspinor,2,cwavef,cwave0,mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

             dot_ndiagr = dot_ndiagr + dot1r*dot2r - dot1i*dot2i
             dot_ndiagi = dot_ndiagi + dot1r*dot2i + dot1i*dot2r
! this should fill all of the iband but only the local cpu jband indices
             d2bbb_k(1,idir1,iband,jband) = d2bbb_k(1,idir1,iband,jband) - &
&             (dot1r*dot2r - dot1i*dot2i)
             d2bbb_k(2,idir1,iband,jband) = d2bbb_k(2,idir1,iband,jband) - &
&             (dot1r*dot2i + dot1i*dot2r)
           end if  ! occ_k
         end do !jband
         ! complete over jband index
         call xmpi_sum(d2bbb_k, mpi_enreg%comm_band, ierr)

         d2bbb_k(:,idir1,iband,:)=d2bbb_k(:,idir1,iband,:)*wtk_k*occ_k(iband)*half
         d2nl_k(1,idir1,ipert1)= &
&         d2nl_k(1,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagr/(nband_kocc*two)
         d2nl_k(2,idir1,ipert1)=&
&         d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*dot_ndiagi/(nband_kocc*two)
       end if ! prtbbb==1

     end do  ! idir1
   end if   ! Compute localization tensor, ipert=natom+1

 end do !  End loop over bands

!Final deallocations
 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)
 ABI_DEALLOCATE(eig2_k)
 ABI_DEALLOCATE(gvnlx1)
 ABI_DEALLOCATE(ffnlk)
 ABI_DEALLOCATE(ffnl1)
 ABI_DEALLOCATE(dkinpw)
 ABI_DEALLOCATE(kinpw1)
 ABI_DEALLOCATE(kpg_k)
 ABI_DEALLOCATE(kpg1_k)
 ABI_DEALLOCATE(ph3d)
 ABI_DATATYPE_DEALLOCATE(dum_cwaveprj)
 if(dtset%prtbbb==1)  then
   ABI_DEALLOCATE(cg_k)
   ABI_DEALLOCATE(cwavef_da)
   ABI_DEALLOCATE(cwavef_db)
 end if

 call timab(102,2,tsec)

 DBG_EXIT("COLL")

  contains
!!***

!!****f* ABINIT/gaugetransfo
!! NAME
!! gaugetransfo
!!
!! FUNCTION
!! This routine allows the passage from the parallel-transport gauge
!! to the diagonal gauge for the first-order wavefunctions
!!
!! INPUTS
!!  cg_k(2,mpw*nspinor*mband*nsppol)=planewave coefficients of wavefunctions
!!                                   for a particular k point.
!!  cwavef(2,npw1_k*nspinor)=first order wavefunction for a particular k point
!!                           in the parallel gauge
!!  comm=mpi communicator for bands
!!  distrb_cycle=array of logical flags to skip certain bands in parallelization scheme
!!  eig_k(mband*nsppol)=GS eigenvalues at k (hartree)
!!  eig1_k(2*nsppol*mband**2)=matrix of first-order eigenvalues (hartree)
!!  iband=band index of the 1WF for which the transformation has to be applied
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  nband_k=number of bands for this k point
!!  npw_k=maximum dimensioned size of npw or wfs at k
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k
!!
!! OUTPUT
!!  cwavef_d(2,npw1_k*nspinor)=first order wavefunction for a particular k point
!!                             in the diagonal gauge
!!
!! PARENTS
!!      dfpt_nstwf
!!
!! CHILDREN
!!
!! SOURCE

subroutine gaugetransfo(cg_k,cwavef,cwavef_d,comm,distrb_cycle,eig_k,eig1_k,iband,nband_k, &
&                      mband,mband_mem,npw_k,npw1_k,nspinor,nsppol,nproc_band,occ_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: iband,mband,mband_mem,nband_k,npw1_k,npw_k,nspinor,nsppol
 integer,intent(in) :: comm, nproc_band
!arrays
 logical, intent(in) :: distrb_cycle(nband_k)
 real(dp),intent(in) :: cg_k(2,npw_k*nspinor*mband_mem),cwavef(2,npw1_k*nspinor)
 real(dp),intent(in) :: eig1_k(2*nsppol*mband**2),eig_k(mband*nsppol)
 real(dp),intent(in) :: occ_k(nband_k)
 real(dp),intent(out) :: cwavef_d(2,npw1_k*nspinor)

!Local variables-------------------------------
!tolerance for non degenerated levels
!scalars
 integer :: jband,jband_me
 real(dp),parameter :: etol=1.0d-3
!arrays
 real(dp) :: cwave0(2,npw1_k*nspinor),eig1(2)

! *********************************************************************

   cwavef_d(:,:) = cwavef(:,:)

   jband_me = 0
   do jband = 1,nband_k !loop over bands
     if (distrb_cycle(jband)) cycle
     jband_me = jband_me + 1

     if ((abs(eig_k(iband)-eig_k(jband)) > etol).and.(abs(occ_k(jband)) > tol8 )) then

       cwave0(:,:) = cg_k(:,1+(jband_me-1)*npw_k*nspinor:jband_me*npw_k*nspinor)

       eig1(1) = eig1_k(2*jband-1+(iband-1)*2*nband_k)
       eig1(2) = eig1_k(2*jband +(iband-1)*2*nband_k)

       cwavef_d(1,:)=cwavef_d(1,:) &
&       - (eig1(1)*cwave0(1,:)-eig1(2)*cwave0(2,:))/(eig_k(jband)-eig_k(iband))
       cwavef_d(2,:)=cwavef_d(2,:) &
&       - (eig1(1)*cwave0(2,:)+eig1(2)*cwave0(1,:))/(eig_k(jband)-eig_k(iband))

     end if

   end do    !loop over bands
   call xmpi_sum(cwavef_d, comm, ierr)
   ! here we have summed the cwavef N times (N-1 too many), but the correction is completed over bands
   cwavef_d = cwavef_d - dble(nproc_band-1)*cwavef

  end subroutine gaugetransfo
!!***

end subroutine dfpt_nstwf
!!***

end module m_dfpt_nstwf
!!***

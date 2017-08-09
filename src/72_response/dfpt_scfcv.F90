!{\src2tex{textfont=tt}}
!!****f* ABINIT/dfpt_scfcv
!! NAME
!! dfpt_scfcv
!!
!! FUNCTION
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! optimum and optionally to compute mixed derivatives of energy.
!!
!! COPYRIGHT
!! Copyright (C) 1999-2017 ABINIT group (XG, DRH, MB, XW, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband*mkmem*nsppol)=pw coefficients of GS wavefunctions at k.
!!  cgq(2,mpw1*nspinor*mband*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  cpus= cpu time limit in seconds
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eew=2nd derivative of Ewald energy (hartree)
!!  efrhar=Contribution from frozen-wavefunction, hartree energy,
!!           to the second-derivative of total energy.
!!  efrkin=Contribution from frozen-wavefunction, kinetic energy,
!!           to the second-derivative of total energy.
!!  efrloc=Contribution from frozen-wavefunction, local potential,
!!           to the second-derivative of total energy.
!!  efrnl=Contribution from frozen-wavefunction, non-local potential,
!!           to the second-derivative of total energy.
!!  efrx1=Contribution from frozen-wavefunction, xc core correction(1),
!!           to the second-derivative of total energy.
!!  efrx2=Contribution from frozen-wavefunction, xc core correction(2),
!!           to the second-derivative of total energy.
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eii=2nd derivative of pseudopotential core energy (hartree)
!!  evdw=DFT-D semi-empirical part of 2nd-order total energy
!!  fermie=fermi energy (Hartree)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  idir=direction of the current perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data for RF symmetries
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates at k
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhohxc.f)
!!  mgfftf=maximum size of 1D FFTs for the "fine" grid (see NOTES in respfn.F90)
!!  mkmem =number of k points treated by this node (GS data)
!!  mkqmem =number of k+q points which can fit in memory (GS data); 0 if use disk
!!  mk1mem =number of k points which can fit in memory (RF data); 0 if use disk
!!  mpert=maximum number of ipert
!!  mpw=maximum dimensioned size of npw for wfs at k.
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point, for each polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf=(effective) number of FFT grid points (for this proc) for the "fine" grid (see NOTES in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid (see NOTES in respfn.F90)
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array.
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsym1=number of symmetry elements in space group consistent with perturbation
!!  n3xccc=dimension of xccc3d1 ; 0 if no XC core correction is used otherwise, nfftf
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k point of the reduced Brillouin zone.
!!  paw_an(natom) <type(paw_an_type)>=paw arrays given on angular mesh for the GS
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastr. containing only symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data for the GS
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  pertcase=fuill index of the perturbation
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic transl. phases, for RF symmetries
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  ph1df(2,3*(2*mgfftf+1)*natom)=one-dimensional structure factor information for the "fine" grid
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfftf)=array for Fourier transform of GS electron density
!!  rhor(nfftf,nspden)=array for GS electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symaf1(nsym1)=anti(ferromagnetic) part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=symmetry operations in real space in terms
!!   of primitive translations
!!  usecprj= 1 if cprj, cprjq arrays are stored in memory
!!  useylmgr = 1 if ylmgr  array is allocated
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  ddk<wfk_t>=ddk file
!!  vpsp1(cplex*nfftf)=first-order derivative of the ionic potential
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vxc(nfftf,nspden)=Exchange-Correlation GS potential (Hartree)
!!  wtk_rbz(nkpt_rbz)=weight for each k point in the reduced Brillouin zone
!!  xccc3d1(cplex*n3xccc)=3D change in core charge density, see n3xccc
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm*useylmgr)= gradients of real spherical harmonics at k
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm*useylmgr1)= gradients of real spherical harmonics at k+q
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  cg1_active(2,mpw1*nspinor*mband*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q. They are orthogonalized to the active.
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,mpert,3,mpert)=non-local contributions to the 2DTEs
!!  d2ovl(2,mpert,3,mpert*usepaw)=1st-order change of WF overlap contributions to the 2DTEs
!!  eberry=energy associated with Berry phase
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  ehart01=inhomogeneous 1st-order Hartree part of 2nd-order total energy
!!    for strain perturbation only (zero otherwise, and not used)
!!  ehart1=1st-order Hartree part of 2nd-order total energy
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues (hartree)
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy.
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  eovl1=1st-order change of wave-functions overlap, part of 2nd-order energy
!!        PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008)
!!  epaw1=1st-order PAW on-site part of 2nd-order total energy.
!!  etotal=total energy (sum of 7 contributions) (hartree)
!!  exc1=1st-order exchange-correlation part of 2nd-order total energy.
!!  gh1c_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(1)}|nK>
!!  gh0c1_set(2,mpw1*nspinor*mband*mk1mem*nsppol*dim_eig2rf)= set of <G|H^{(0)}|\Psi^{(1)}>
!!      The wavefunction is orthogonal to the active space (for metals). It is not
!!      coherent with cg1.
!!  resid(mband*nkpt_rbz*nsppol)=residuals for each band over all k points
!!   of the reduced Brillouin zone, and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  conv_retcode=return code, 0 if convergence was achieved.
!!
!! SIDE EFFECTS
!!  cg1(2,mpw1*nspinor*mband*mk1mem*nsppol)=updated wavefunctions (ortho. to occ. states);
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!  mpi_enreg=informations about MPI parallelization
!!  rhog1(2,nfftf)=array for Fourier transform of RF electron density
!!  rhor1(cplex*nfftf,nspden)=array for RF electron density in electrons/bohr**3.
!!  === if psps%usepaw==1
!!    pawrhoij1(natom) <type(pawrhoij_type)>= 1st-order paw rhoij occupancies and related data
!!
!! PARENTS
!!      dfpt_looppert
!!
!! CHILDREN
!!      ab7_mixing_deallocate,ab7_mixing_new,ab7_mixing_use_disk_cache,appdig
!!      calcdensph,destroy_efield,dfpt_etot,dfpt_newvtr,dfpt_nselt,dfpt_nstdy
!!      dfpt_nstpaw,dfpt_rhofermi,dfpt_rhotov,dfpt_vtorho,dfptff_bec,dfptff_die
!!      dfptff_ebp,dfptff_edie,dfptff_initberry,fftdatar_write_from_hdr,fourdp
!!      getcut,hdr_update,metric,newfermie1,paw_an_free,paw_an_init
!!      paw_an_nullify,paw_an_reset_flags,paw_ij_free,paw_ij_init
!!      paw_ij_nullify,paw_ij_reset_flags,pawcprj_alloc,pawcprj_free
!!      pawcprj_getdim,pawdenpot,pawdij,pawdijfr,pawmknhat,pawnhatfr
!!      pawrhoij_alloc,pawrhoij_free,qmatrix,rf2_getidirs,scprqt,status,symdij
!!      timab,wfk_close,wrtout,xmpi_barrier,xmpi_isum,xmpi_sum,xmpi_wait
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine dfpt_scfcv(atindx,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&
&  dielt,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset,&
&  d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&
&  enl0,enl1,eovl1,epaw1,etotal,evdw,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,indkpt1,&
&  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&  mpert,mpi_enreg,mpw,mpw1,my_natom,nattyp,nband_rbz,ncpgr,&
&  nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&  nsym1,n3xccc,occkq,occ_rbz,&
&  paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&
&  pertcase,phnons1,ph1d,ph1df,&
&  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&
&  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,&
&  usecprj,useylmgr,useylmgr1,ddk_f,vpsp1,vtrial,vxc,&
&  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1,zeff,conv_retcode)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_ab7_mixing
 use m_efield
 use m_errors
 use m_profiling_abi
 use m_wfk
 use m_xmpi
 use m_nctk
 use m_hdr
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use m_cgtools,  only : mean_fftr
 use m_fstrings, only : int2char4, sjoin
 use m_time,     only : abi_wtime, sec2str
 use m_io_tools, only : open_file
 use m_exit,     only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_mpinfo,   only : iwrite_fftdatar
 use m_ioarr,    only : ioarr, fftdatar_write_from_hdr, fort_denpot_skip
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_paw_an,   only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify, paw_an_reset_flags
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij, only : pawrhoij_type, pawrhoij_alloc, pawrhoij_free, pawrhoij_io
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_getdim
 use m_pawdij,   only : pawdij, pawdijfr, symdij
 use m_pawfgr,   only : pawfgr_type
 use m_rf2,      only : rf2_getidirs

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'dfpt_scfcv'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_32_util
 use interfaces_41_geometry
 use interfaces_53_ffts
 use interfaces_54_abiutil
 use interfaces_56_recipspace
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_72_response, except_this_one => dfpt_scfcv
!End of the abilint section

 implicit none

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 integer,intent(in) :: cplex,dim_eig2rf,idir,ipert,mgfftf,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mpert,mpw,mpw1,my_natom,n3xccc,ncpgr,nfftf
 integer,intent(in) :: nkpt,nkpt_rbz,nkxc,nspden
 integer,intent(in) :: nsym1,pertcase,prtbbb,usecprj,useylmgr,useylmgr1
 integer,intent(inout) :: initialized
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: atindx(dtset%natom)
 integer,intent(out) :: blkflg(3,mpert,3,mpert)
 integer,intent(in) :: indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem),nattyp(psps%ntypat)
 integer,intent(in) :: nband_rbz(nkpt_rbz*dtset%nsppol)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 integer,intent(out) :: conv_retcode
 real(dp),intent(in) :: cpus,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,eii
 real(dp),intent(out) :: eberry,edocc,eeig0,ehart01,ehart1,ek0,ek1,eloc0,elpsp1,enl0
 real(dp),intent(out) :: enl1,eovl1,epaw1,etotal,evdw,exc1,residm
 real(dp),intent(inout) :: fermie
 real(dp),intent(in) :: qphon(3)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband*mkmem*dtset%nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol)
 real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*prtbbb)
 real(dp),intent(out) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert)
 real(dp),intent(out) :: d2ovl(2,3,mpert,3,mpert*psps%usepaw)
 real(dp),intent(in) :: dielt(3,3)
 real(dp),intent(in) :: doccde_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: docckqde(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(out) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: eigenq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz),kxc(nfftf,nkxc)
 real(dp),intent(in) :: nhat(nfftf,dtset%nspden)
 real(dp),intent(in) :: occ_rbz(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: occkq(dtset%mband*nkpt_rbz*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom),ph1df(2,3*(2*mgfftf+1)*dtset%natom)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp),intent(out) :: resid(dtset%mband*nkpt_rbz*nspden)
 real(dp),intent(in) :: rhog(2,nfftf),rhor(nfftf,nspden),rprimd(3,3)
 real(dp),intent(inout) :: rhog1(2,nfftf),rhor1(cplex*nfftf,nspden),xred(3,dtset%natom)
 real(dp),target,intent(in) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: vpsp1(cplex*nfftf),vxc(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xccc3d1(cplex*n3xccc)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3+6*((ipert-dtset%natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(in) :: zeff(3,3,dtset%natom)
 type(pawcprj_type),intent(in) :: cprj(dtset%natom,dtset%nspinor*dtset%mband*mkmem*dtset%nsppol*usecprj)
 type(pawcprj_type),intent(in) :: cprjq(dtset%natom,dtset%nspinor*dtset%mband*mkqmem*dtset%nsppol*usecprj)
 type(datafiles_type),intent(in) :: dtfil
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(paw_an_type),intent(in) :: paw_an(my_natom*psps%usepaw)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(psps%ntypat*psps%usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij1(my_natom*psps%usepaw)
 type(pawtab_type), intent(in) :: pawtab(psps%ntypat*psps%usepaw)
 type(MPI_type),intent(inout) :: mpi_enreg
 type(wfk_t),intent(inout) :: ddk_f(4)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=12,response=1
 integer :: afford,bantot_rbz,choice,cplex_rhoij,dbl_nnsclo
 integer :: has_dijfr,iatom,ider,idir_dum,idir_paw1,ierr,iexit,errid,denpot
 integer :: iprcel,iscf10_mod,iscf_mod,ispden,ispmix
 integer :: istep,itypat,izero,lmn2_size,me,mgfftdiel,mvdum
 integer :: nfftdiel,nfftmix,nfftotf,nhat1grdim,npawmix,npwdiel,nspden_rhoij,nstep,nzlmopt
 integer :: optene,optfr,option,optres,prtfor,quit,quit_sum,qzero
 integer :: my_quit,quitsum_request,timelimit_exit,varid,ncerr,ncid
 integer ABI_ASYNC :: quitsum_async
 integer :: rdwrpaw,spaceComm,sz1,sz2,usexcnhat,Z_kappa
 integer :: mpi_comm_sphgrid
 logical :: need_fermie1,paral_atom,use_nhat_gga
 real(dp) :: wtime_step,now,prev
 real(dp) :: born,born_bar,boxcut,deltae,diffor,diel_q,dum,ecut,ecutf,elast
 real(dp) :: epawdc1_dum,evar,fe1fixed,fermie1,gsqcut,qphon_norm,maxfor,renorm,res2,res3,residm2
 real(dp) :: ucvol,vxcavg,elmag1
 character(len=500) :: msg
 character(len=fnlen) :: fi1o
 character(len=fnlen) :: fi1o_vtk
 type(ab7_mixing_object) :: mix
 type(efield_type) :: dtefield
!arrays
 integer :: ngfftmix(18)
 integer,allocatable :: dimcprj(:),pwindall(:,:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dielar(7)
 real(dp) :: favg(3),gmet(3,3),gprimd(3,3),q_cart(3),qphon2(3),qred2cart(3,3)
 real(dp) :: rmet(3,3),tollist(12),tsec(2)
 real(dp) :: zeff_red(3),zeff_bar(3,3)
 real(dp) :: intgden(dtset%nspden,dtset%natom),dentot(dtset%nspden)
!real(dp) :: zdmc_red(3),zdmc_bar(3,3),mean_rhor1(1) !dynamic magnetic charges and mean density 
 real(dp),allocatable :: dielinv(:,:,:,:,:)
 real(dp),allocatable :: fcart(:,:),nhat1(:,:),nhat1gr(:,:,:),nhatfermi(:,:),nvresid1(:,:),nvresid2(:,:)
 real(dp),allocatable :: qmat(:,:,:,:,:,:),resid2(:),rhog2(:,:),rhor2(:,:),rhorfermi(:,:)
 real(dp),allocatable :: susmat(:,:,:,:,:),vhartr1(:),vxc1(:,:)
 real(dp),allocatable :: vhartr1_tmp(:,:)
 real(dp),allocatable,target :: vtrial1(:,:),vtrial2(:,:)
 real(dp),pointer :: vtrial1_tmp(:,:)
 type(pawcprj_type),allocatable :: cprj1(:,:)
 type(paw_an_type),allocatable :: paw_an1(:)
 type(paw_ij_type),allocatable :: paw_ij1(:)
 type(pawrhoij_type),allocatable :: pawrhoijfermi(:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(120,1,tsec)
 call timab(154,1,tsec)

 call status(0,dtfil%filstat,iexit,level,'init          ')

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(ABI_FUNC) == ABI_FUNC) then
   write(std_out,*)"Enabling timelimit check in function: ",trim(ABI_FUNC)," with timelimit: ",trim(sec2str(get_timelimit()))
 end if

!Parallelism data
 spaceComm=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 paral_atom=(my_natom/=dtset%natom)
 my_atmtab=>mpi_enreg%my_atmtab

 _IBM6("XLF in dfpt_scfcv")

!Save some variables from dataset definition
 ecut=dtset%ecut
 ecutf=ecut;if (psps%usepaw==1) ecutf=dtset%pawecutdg
 iprcel=dtset%iprcel
 tollist(1)=dtset%tolmxf;tollist(2)=dtset%tolwfr
 tollist(3)=dtset%toldff;tollist(4)=dtset%toldfe
 tollist(6)=dtset%tolvrs;tollist(7)=dtset%tolrff
 nfftotf=product(ngfftf(1:3))
 nstep=dtset%nstep
 iscf_mod=dtset%iscf
 iscf10_mod=mod(iscf_mod,10)

 qzero=0; if(qphon(1)**2+qphon(2)**2+qphon(3)**2 < tol14) qzero=1

 need_fermie1=((qzero==1.and.dtset%frzfermi==0.and.nstep>0).and.&
& (dtset%occopt>=3.and.dtset%occopt<=8).and. &
& (ipert<=dtset%natom.or.ipert==dtset%natom+3.or.ipert==dtset%natom+4.or.ipert==dtset%natom+5)) !rfmagn deb

!The value of iscf must be modified if ddk perturbation, see dfpt_looppert.f
 if (ipert==dtset%natom+1.or.ipert==dtset%natom+10.or.ipert==dtset%natom+11) iscf_mod=-3

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute large sphere cut-off gsqcut
 qphon2(:)=zero;if (psps%usepaw==1) qphon2(:)=qphon(:)
 call getcut(boxcut,ecutf,gmet,gsqcut,dtset%iboxcut,std_out,qphon2,ngfftf)

!Some variables need to be initialized/nullify at start
 quit=0 ; dbl_nnsclo=0 ; elast=zero; conv_retcode = -1
 optres=merge(0,1,abs(iscf_mod)<10)
 usexcnhat=0
!This might be taken away later
 edocc=zero ; eeig0=zero ; ehart01=zero ; ehart1=zero ; ek0=zero ; ek1=zero
 eloc0=zero ; elpsp1=zero ; enl0=zero ; enl1=zero ; eovl1=zero; exc1=zero
 deltae=zero ; fermie1=zero ; epaw1=zero ; eberry=zero ; elmag1=zero

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor,res2,prtfor,fcart are here initialized to 0)
 choice=1 ; prtfor=0 ; diffor=zero ; res2=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))

 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,0,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)

!Allocations/initializations for PAW only
 if(psps%usepaw==1) then
   usexcnhat=maxval(pawtab(:)%usexcnhat)
   use_nhat_gga=(dtset%xclevel==2.and.dtset%pawnhatxc>0.and.usexcnhat>0)
!  1st-order compensation density
   ABI_ALLOCATE(nhat1,(cplex*nfftf,dtset%nspden))
   nhat1=zero
!  Projections of 1-st order WF on nl projectors
   ABI_DATATYPE_ALLOCATE(cprj1,(dtset%natom,dtset%nspinor*dtset%mband*mk1mem*dtset%nsppol*usecprj))
   if (usecprj==1.and.mk1mem/=0) then
     !cprj ordered by atom-type
     ABI_ALLOCATE(dimcprj,(dtset%natom))
     call pawcprj_getdim(dimcprj,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
     call pawcprj_alloc(cprj1,0,dimcprj)
     ABI_DEALLOCATE(dimcprj)
   end if
!  1st-order arrays/variables related to the PAW spheres
   ABI_DATATYPE_ALLOCATE(paw_an1,(my_natom))
   ABI_DATATYPE_ALLOCATE(paw_ij1,(my_natom))
   call paw_an_nullify(paw_an1)
   call paw_ij_nullify(paw_ij1)

   has_dijfr=0;if (ipert/=dtset%natom+1.and.ipert/=dtset%natom+10) has_dijfr=1
   call paw_an_init(paw_an1,dtset%natom,dtset%ntypat,0,dtset%nspden,cplex,dtset%pawxcdev,&
&   dtset%typat,pawang,pawtab,has_vxc=1,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

   call paw_ij_init(paw_ij1,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,0,dtset%natom,&
&   dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijfr=has_dijfr,&
&   mpi_atmtab=mpi_enreg%my_atmtab, comm_atom=mpi_enreg%comm_atom)
 else
   ABI_ALLOCATE(nhat1,(0,0))
   ABI_DATATYPE_ALLOCATE(cprj1,(0,0))
   ABI_DATATYPE_ALLOCATE(paw_an1,(0))
   ABI_DATATYPE_ALLOCATE(paw_ij1,(0))
 end if ! PAW

!Various allocations (potentials)
 ABI_ALLOCATE(vhartr1,(cplex*nfftf))
 ABI_ALLOCATE(vtrial1,(cplex*nfftf,nspden))
! TODO: for non collinear case this should always be nspden, in NCPP case as well!!!
 ABI_ALLOCATE(vxc1,(cplex*nfftf,nspden*(1-usexcnhat)*psps%usepaw)) ! Not always needed
 vtrial1_tmp => vtrial1   ! this is to avoid errors when vtrial1_tmp is unused

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (iscf_mod>0.or.iscf_mod==-3) then
   ABI_ALLOCATE(nvresid1,(cplex*nfftf,dtset%nspden))
   if (nstep==0) nvresid1=zero
   if ((dtset%getddb .ne. 0 .or. dtset%irdddb .ne.0) .and. qzero .ne. 1) then
     ABI_ALLOCATE(nvresid2,(cplex*nfftf,dtset%nspden))
     if (nstep==0) nvresid2=zero
   end if
 else
   ABI_ALLOCATE(nvresid1,(0,0))
 end if
 if(nstep>0 .and. iscf_mod>0) then
   dielar(1)=dtset%diecut;dielar(2)=dtset%dielng
   dielar(3)=dtset%diemac;dielar(4)=dtset%diemix
   dielar(5)=dtset%diegap;dielar(6)=dtset%dielam
   dielar(7)=dtset%diemix;if (dtset%iscf>=10) dielar(7)=dtset%diemixmag
!  Additional allocation for mixing within PAW
   npawmix=0
   if(psps%usepaw==1) then
     do iatom=1,my_natom
       itypat=pawrhoij1(iatom)%itypat
       lmn2_size=pawtab(itypat)%lmn2_size
       pawrhoij1(iatom)%use_rhoijres=1
       sz1=pawrhoij1(iatom)%cplex*lmn2_size;sz2=pawrhoij1(iatom)%nspden
       ABI_ALLOCATE(pawrhoij1(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,pawrhoij1(iatom)%nspden
         pawrhoij1(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij1(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij1(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij1(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij1(iatom)%nspden*pawtab(itypat)%lmnmix_sz*pawrhoij1(iatom)%cplex
     end do
   end if
   denpot = AB7_MIXING_POTENTIAL
   if (dtset%iscf > 10) denpot = AB7_MIXING_DENSITY
   if (psps%usepaw==1.and.dtset%pawmixdg==0) then
     ispmix=AB7_MIXING_FOURRIER_SPACE;nfftmix=dtset%nfft;ngfftmix(:)=dtset%ngfft(:)
   else
     ispmix=AB7_MIXING_REAL_SPACE;nfftmix=nfftf;ngfftmix(:)=ngfftf(:)
   end if
   if (iscf10_mod == 5 .or. iscf10_mod == 6) then
     call ab7_mixing_new(mix, iscf10_mod, denpot, cplex, &
&     nfftf, dtset%nspden, npawmix, errid, msg, dtset%npulayit)
   else
     call ab7_mixing_new(mix, iscf10_mod, denpot, max(cplex, ispmix), &
&     nfftmix, dtset%nspden, npawmix, errid, msg, dtset%npulayit)
   end if
   if (errid /= AB7_NO_ERROR) then
     MSG_ERROR(msg)
   end if
   if (dtset%mffmem == 0) then
     call ab7_mixing_use_disk_cache(mix, dtfil%fnametmp_fft)
   end if
 end if ! iscf, nstep

!Here, allocate arrays for computation of susceptibility and dielectric matrix or for TDDFT
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
!  Here, for TDDFT, artificially set iprcel . Also set a variable to reduce the memory needs.
   afford=1
   if(iscf_mod==-1) then
     iprcel=21
     afford=0
   end if
   npwdiel=1
   mgfftdiel=1
   nfftdiel=1
!  Now, performs allocation
!  CAUTION : the dimensions are still those of GS, except for phnonsdiel
   ABI_ALLOCATE(dielinv,(2,npwdiel*afford,nspden,npwdiel,nspden))
   ABI_ALLOCATE(susmat,(2,npwdiel*afford,nspden,npwdiel,nspden))
 end if

!Initialize Berry-phase related stuffs
 if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
   ABI_ALLOCATE(pwindall,(max(mpw,mpw1)*mkmem,8,3))
   call dfptff_initberry(dtefield,dtset,gmet,kg,kg1,dtset%mband,mkmem,mpi_enreg,&
&   mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,occ_rbz,pwindall,rprimd)
!  calculate inverse of the overlap matrix
   ABI_ALLOCATE(qmat,(2,dtefield%mband_occ,dtefield%mband_occ,nkpt,2,3))
   call qmatrix(cg,dtefield,qmat,mpw,mpw1,mkmem,dtset%mband,npwarr,nkpt,dtset%nspinor,dtset%nsppol,pwindall)
 else
   ABI_ALLOCATE(pwindall,(0,0,0))
   ABI_ALLOCATE(qmat,(0,0,0,0,0,0))
 end if

 call timab(154,2,tsec)

!######################################################################
!PERFORM ELECTRONIC ITERATIONS
!######################################################################

!Offer option of computing 2nd-order total energy with existing
!wavefunctions when nstep<=0, else do nstep iterations
!Note that for non-self-consistent calculations, this loop will be exited
!after the first call to dfpt_vtorho

!Pass through the first routines even when nstep==0
!write(std_out,*) 'dfpt_scfcv, nstep=', max(1,nstep)

 quitsum_request = xmpi_request_null; timelimit_exit = 0

 do istep=1,max(1,nstep)

   ! Handle time limit condition.
   if (istep == 1) prev = abi_wtime()
   if (istep  > 1) then
     now = abi_wtime()
     wtime_step = now - prev
     prev = now
     call wrtout(std_out,sjoin("dfpt_scfcv: previous iteration took ",sec2str(wtime_step)))

     if (have_timelimit_in(ABI_FUNC)) then
       if (istep > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then
           write(msg,"(3a)")"Approaching time limit ",trim(sec2str(get_timelimit())),". Will exit istep loop in dfpt_scfcv."
           MSG_COMMENT(msg)
           call wrtout(ab_out, msg, "COLL")
           timelimit_exit = 1
           exit
         end if
       end if

       my_quit = 0; if (now - get_start_time() + 2.15 * wtime_step > get_timelimit()) my_quit = 1
       call xmpi_isum(my_quit,quitsum_async,spacecomm,quitsum_request,ierr)
     end if
   end if

!  ######################################################################
!  The following steps are done once
!  ----------------------------------------------------------------------
   if (istep==1)then

!    PAW only: compute frozen part of 1st-order compensation density
!    and frozen part of psp strengths Dij
!    ----------------------------------------------------------------------
     if (psps%usepaw==1) then
       optfr=0
       idir_paw1 = idir
       if (ipert==dtset%natom+11) then
         call rf2_getidirs(idir,idir_dum,idir_paw1)
       end if
       call pawdijfr(cplex,gprimd,idir_paw1,ipert,my_natom,dtset%natom,nfftf,ngfftf,nspden,&
&       psps%ntypat,optfr,paw_ij1,pawang,pawfgrtab,pawrad,pawtab,qphon,&
&       rprimd,ucvol,vpsp1,vtrial,vxc,xred,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)

       if ((iscf_mod>=0.or.usexcnhat==0).and.(dtset%pawstgylm/=0)) then
         ider=0;if ((ipert<=dtset%natom).and.(use_nhat_gga)) ider=1
         call pawnhatfr(ider,idir_paw1,ipert,my_natom,dtset%natom,nspden,psps%ntypat,&
&         pawang,pawfgrtab,pawrhoij,pawtab,rprimd,&
&         mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
       end if
     end if
!    PAW only: we sometimes have to compute 1st-order compensation density
!    and eventually add it to density from 1st-order WFs
!    ----------------------------------------------------------------------
     nhat1grdim=0
     ABI_ALLOCATE(nhat1gr,(0,0,0))
     if (psps%usepaw==1.and.ipert/=dtset%natom+1.and.ipert/=dtset%natom+10) then
       call timab(564,1,tsec)
       nhat1grdim=0;if (dtset%xclevel==2) nhat1grdim=usexcnhat*dtset%pawnhatxc
       ider=2*nhat1grdim;izero=0
       if (nhat1grdim>0)   then
         ABI_DEALLOCATE(nhat1gr)
         ABI_ALLOCATE(nhat1gr,(cplex*nfftf,dtset%nspden,3*nhat1grdim))
       end if
       call pawmknhat(dum,cplex,ider,idir_paw1,ipert,izero,gprimd,my_natom,dtset%natom,&
&       nfftf,ngfftf,nhat1grdim,nspden,psps%ntypat,pawang,pawfgrtab,nhat1gr,nhat1,&
&       pawrhoij1,pawrhoij,pawtab,qphon,rprimd,ucvol,dtset%usewvl,xred,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
       if (dtfil%ireadwf/=0.and.dtset%get1den==0.and.dtset%ird1den==0.and.initialized==0) then
         rhor1(:,:)=rhor1(:,:)+nhat1(:,:)
         call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,ngfftf,dtset%paral_kgb,0)
       end if
       call timab(564,2,tsec)
     end if
!    Set initial guess for 1st-order potential
!    ----------------------------------------------------------------------
     call status(istep,dtfil%filstat,iexit,level,'get vtrial1   ')
     option=1;optene=0;if (iscf_mod==-2) optene=1
     call dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,elmag1,gmet,gprimd,gsqcut,idir,ipert,&
&     dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1gr,nhat1grdim,&
&     nkxc,nspden,n3xccc,optene,option,dtset%paral_kgb,dtset%qptn,&
&     rhog,rhog1,rhor,rhor1,rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,&
&     nvresid1,res2,vtrial1,vxc1,xccc3d1)

!    For Q=0 and metallic occupation, initialize quantities needed to
!    compute the first-order Fermi energy
!    ----------------------------------------------------------------------
     if (need_fermie1) then
       ABI_ALLOCATE(rhorfermi,(cplex*nfftf,nspden))
       if (psps%usepaw==1.and.usexcnhat==0) then
         ABI_ALLOCATE(nhatfermi,(cplex*nfftf,nspden))
       else
         ABI_ALLOCATE(nhatfermi,(0,0))
       end if
       ABI_DATATYPE_ALLOCATE(pawrhoijfermi,(my_natom*psps%usepaw))
       if (psps%usepaw==1) then
         cplex_rhoij=max(cplex,dtset%pawcpxocc)
         nspden_rhoij=dtset%nspden;if (dtset%pawspnorb>0.and.dtset%nspinor==2) nspden_rhoij=4
         call pawrhoij_alloc(pawrhoijfermi,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&         dtset%nsppol,dtset%typat,pawtab=pawtab,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_atom=mpi_enreg%comm_atom)
       end if
       
       call dfpt_rhofermi(cg,cgq,cplex,cprj,cprjq,&
&       doccde_rbz,docckqde,dtfil,dtset,eigenq,eigen0,eigen1,fe1fixed,gmet,gprimd,idir,&
&       indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,mkmem,mkqmem,mk1mem,mpi_enreg,&
&       mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1,&
&       nspden,dtset%nsppol,nsym1,occkq,occ_rbz,&
&       paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&
&       phnons1,ph1d,dtset%prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrc1,symrl1,&
&       ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)
     end if

   end if ! End the condition of istep==1

!  ######################################################################
!  The following steps are done at every iteration
!  ----------------------------------------------------------------------

   if (psps%usepaw==1)then
!    Computation of "on-site" 2nd-order energy, first-order potentials, first-order densities
     nzlmopt=0;if (istep==2.and.dtset%pawnzlm>0) nzlmopt=-1
     if (istep>2) nzlmopt=dtset%pawnzlm
     call paw_an_reset_flags(paw_an1) ! Force the recomputation of on-site potentials
     call paw_ij_reset_flags(paw_ij1,self_consistent=.true.) ! Force the recomputation of Dij
     option=0;if (dtset%iscf>0.and.dtset%iscf<10.and.nstep>0) option=1
     call status(istep,dtfil%filstat,iexit,level,'call pawdenpot')
     call pawdenpot(dum,epaw1,epawdc1_dum,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an1,paw_an,paw_ij1,pawang,&
&     dtset%pawprtvol,pawrad,pawrhoij1,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&     dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp, &
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!    First-order Dij computation
     call status(istep,dtfil%filstat,iexit,level,'call pawdij   ')
     call timab(561,1,tsec)
     if (has_dijfr>0) then
       !vpsp1 contribution to Dij already stored in frozen part of Dij
       ABI_ALLOCATE(vtrial1_tmp,(cplex*nfftf,nspden))
       vtrial1_tmp=vtrial1
       do ispden=1,min(dtset%nspden,2)
         vtrial1_tmp(:,ispden)=vtrial1_tmp(:,ispden)-vpsp1(:)
       end do
     else
       vtrial1_tmp => vtrial1
     end if
     call pawdij(cplex,dtset%enunit,gprimd,ipert,my_natom,dtset%natom,&
&     nfftf,nfftotf,dtset%nspden,psps%ntypat,paw_an1,paw_ij1,pawang,&
&     pawfgrtab,dtset%pawprtvol,pawrad,pawrhoij1,dtset%pawspnorb,pawtab,&
&     dtset%pawxcdev,qphon,dtset%spnorbscl,ucvol,dtset%charge,vtrial1_tmp,vxc1,xred,&
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     if (has_dijfr>0) then
       ABI_DEALLOCATE(vtrial1_tmp)
     end if

     call status(istep,dtfil%filstat,iexit,level,'call symdij   ')
     call symdij(gprimd,indsy1,ipert,my_natom,dtset%natom,nsym1,psps%ntypat,0,&
&     paw_ij1,pawang1,dtset%pawprtvol,pawtab,rprimd,symaf1,symrc1, &
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,&
&     qphon=qphon)
     call timab(561,2,tsec)
   end if ! end usepaw section

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------
   call status(istep,dtfil%filstat,iexit,level,'loop istep    ')

   if(iscf_mod>0.and.nstep>0)then
     write(msg, '(a,a,i4)' )ch10,' ITER STEP NUMBER  ',istep
     call wrtout(std_out,msg,'COLL')
   end if

!  For Q=0 and metallic occupation, calculate the first-order Fermi energy
   if (need_fermie1) then
     call newfermie1(cplex,fermie1,fe1fixed,ipert,istep,dtset%ixc,my_natom,dtset%natom,&
&     nfftf,nfftotf,nhatfermi,nspden,dtset%ntypat,dtset%occopt,paw_an,paw_an1,paw_ij1,pawang,&
&     dtset%pawnzlm,pawrad,pawrhoij1,pawrhoijfermi,pawtab,dtset%pawxcdev,&
&     dtset%prtvol,rhorfermi,ucvol,psps%usepaw,usexcnhat,vtrial1,vxc1,dtset%xclevel,&
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   end if

!  No need to continue and call dfpt_vtorho, when nstep==0
   if(nstep==0) exit

!  #######################e1magh###############################################
!  Compute the 1st-order density rho1 from the 1st-order trial potential
!  ----------------------------------------------------------------------

   call dfpt_vtorho(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,&
&   dbl_nnsclo,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,edocc,&
&   eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,fermie1,gh0c1_set,gh1c_set,&
&   gmet,gprimd,idir,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,&
&   mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,&
&   nhat1,nkpt_rbz,npwarr,npwar1,res2,nspden,dtset%nsppol,nsym1,dtset%ntypat,nvresid1,&
&   occkq,occ_rbz,optres,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,&
&   pawrhoij1,pawtab,phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid,residm,rhog1,&
&   rhor1,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,usecprj,useylmgr1,ddk_f,&
&   vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)

   if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
&   dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then

!    calculate \Omega E \cdot P term
     if (ipert<=dtset%natom) then
!      phonon perturbation
       call  dfptff_ebp(cg,cg1,dtefield,eberry,dtset%mband,mkmem,&
&       mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat)
     else if (ipert==dtset%natom+2) then
!      electric field perturbation
       call  dfptff_edie(cg,cg1,dtefield,eberry,idir,dtset%mband,mkmem,&
&       mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
     end if
   end if

   !if (ipert==dtset%natom+5) then
   !calculate 1st order magnetic potential contribution to the energy
   !  call dfpt_e1mag(e1mag,rhor1,rhog1);
   !endif

!  ######################################################################
!  Skip out of step loop if non-SCF (completed)
!  ----------------------------------------------------------------------

!  Indeed, nstep loops have been done inside dfpt_vtorho
   if (iscf_mod<=0 .and. iscf_mod/=-3) exit

!  ######################################################################
!  In case of density mixing , compute the total 2nd-order energy,
!  check the exit criterion, then mix the 1st-order density
!  ----------------------------------------------------------------------

   if (iscf_mod>=10) then
     optene = 1 ! use double counting scheme
     call dfpt_etot(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&     efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,elmag1,ipert,dtset%natom,optene)

     call timab(152,1,tsec)
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,kpt_rbz,maxfor,&
&     mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&     nstep,occ_rbz,0,prtfor,0,&
&     quit,res2,resid,residm,response,&
&     tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)
     call timab(152,2,tsec)

     if (istep==nstep) quit=1
!    If criteria in scprqt say to quit, then exit the loop over istep
     quit_sum=quit
     call xmpi_sum(quit_sum,spaceComm,ierr)

     if (quit_sum>0) exit
     call status(istep,dtfil%filstat,iexit,level,'call newrho   ')
!    INSERT HERE CALL TO NEWRHO3 : to be implemented
     if (psps%usepaw==1) then
       MSG_BUG("newrho3 not implemented: use potential mixing!")
     end if
     initialized=1
   end if

!  ######################################################################
!  Compute the new 1st-order potential from the 1st-order density
!  ----------------------------------------------------------------------

   if (ipert<dtset%natom+10) then
     optene=1
     call status(istep,dtfil%filstat,iexit,level,'call dfpt_rhotov   ')
     call dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,elmag1,gmet,gprimd,gsqcut,idir,ipert,&
&     dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
&     nspden,n3xccc,optene,optres,dtset%paral_kgb,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&     rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,vxc1,xccc3d1)
   end if

!  ######################################################################
!  In case of potential mixing , compute the total 2nd-order energy,
!  check the exit criterion, then mix the 1st-order potential
!  ----------------------------------------------------------------------

   if (iscf_mod<10) then

!    PAW: has to compute here the "on-site" 2nd-order energy
     if (psps%usepaw==1) then
       nzlmopt=0;if (istep==1.and.dtset%pawnzlm>0) nzlmopt=-1
       if (istep>1) nzlmopt=dtset%pawnzlm
       call paw_an_reset_flags(paw_an1) ! Force the recomputation of on-site potentials
       option=2
       call pawdenpot(dum,epaw1,epawdc1_dum,ipert,dtset%ixc,my_natom,dtset%natom,dtset%nspden,&
&       psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an1,paw_an,paw_ij1,pawang,dtset%pawprtvol,&
&       pawrad,pawrhoij1,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,&
&       dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     end if

     optene = 0 ! use direct scheme
     call dfpt_etot(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
&     efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,elmag1,ipert,dtset%natom,optene)

     call timab(152,1,tsec)
     choice=2
     call status(istep,dtfil%filstat,iexit,level,'print info    ')
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,kpt_rbz,maxfor,&
&     mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&     nstep,occ_rbz,0,prtfor,0,&
&     quit,res2,resid,residm,response,&
&     tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)
     call timab(152,2,tsec)

!    If criteria in scprqt say to quit, then exit the loop over istep
     quit_sum=quit
     call xmpi_sum(quit_sum,spaceComm,ierr)
     if (quit_sum>0) exit

     ! TODO
     ! Better error handling is the SCF cycle goes bananas:
     !   Write a BIG warning in the output file and save the wavefunctions.
     !   so that we can restart.
     if(iscf_mod/=-3)then
!      Note that nvresid1 and vtrial1 are called vresid and vtrial inside this routine
       call dfpt_newvtr(cplex,dbl_nnsclo,dielar,dtset,etotal,pawfgr%fintocoa,&
&       initialized,iscf_mod,ispmix,istep,mix,pawfgr%coatofin,&
&       mpi_enreg,my_natom,nfftf,nfftmix,ngfftf,ngfftmix,npawmix,pawrhoij1,&
&       qphon,rhor1,rprimd,psps%usepaw,nvresid1,vtrial1)
       initialized=1
     end if
   end if

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  Note that there are different "exit" instructions within the loop
!  ######################################################################
 end do ! istep

 ! Avoid pending requests if itime == ntime.
 call xmpi_wait(quitsum_request,ierr)
 if (timelimit_exit == 1) istep = istep - 1

!SP : Here read the _DDB file and extract the Born effective charge and
!     dielectric constant.
! The idea is to supress the divergence due to a residual Born effective charge
! by renormalizing the v_hart1. For this, the difference between the ionic
! Z_kappa and the Born effective charge divided by the dielectric constant is used.
! ---------------------------------------------------------------------------------
 if ((dtset%getddb .ne. 0 .or. dtset%irdddb .ne.0) .and. qzero .ne. 1) then
   ABI_ALLOCATE(rhor2,(cplex*nfftf,nspden))
   ABI_ALLOCATE(resid2,(dtset%mband*nkpt_rbz*nspden))
   ABI_ALLOCATE(rhog2,(2,nfftf))
   ABI_ALLOCATE(vtrial2,(cplex*nfftf,nspden))

   Z_kappa = nint(psps%ziontypat(dtset%typat(ipert))) ! Charge ionic from the psp
   qred2cart = two_pi*gprimd
   q_cart = MATMUL(qred2cart,qphon)
   q_cart = q_cart/SQRT(dot_product(q_cart,q_cart))
   diel_q = dot_product(MATMUL(dielt,q_cart),q_cart)
   zeff_bar = SUM(zeff(:,:,:),DIM=3)/dtset%natom
   zeff_red = MATMUL(zeff_bar(:,:),rprimd(:,idir))/two_pi
   qphon_norm = SQRT(dot_product(qphon,qphon))
   q_cart = MATMUL(qred2cart,qphon)
   born_bar = dot_product(q_cart,zeff_red(:))
   zeff_red = MATMUL(zeff(:,:,ipert),rprimd(:,idir))/two_pi
   born = dot_product(q_cart,zeff_red(:))

! To avoid problem of divergence (0/0) we add a small value to qphon
   qphon2 = qphon + tol6
   renorm = (1-(qphon2(idir)*Z_kappa-(born-born_bar)/diel_q)/(qphon2(idir)*Z_kappa-born/diel_q))

   vtrial2(:,1) = vtrial1(:,1) -renorm*vhartr1

   call dfpt_vtorho(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,&
&   dbl_nnsclo,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,edocc,&
&   eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,fermie1,gh0c1_set,gh1c_set,&
&   gmet,gprimd,idir,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,&
&   mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,&
&   nhat1,nkpt_rbz,npwarr,npwar1,res3,nspden,dtset%nsppol,nsym1,dtset%ntypat,nvresid2,&
&   occkq,occ_rbz,optres,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,&
&   pawrhoij1,pawtab,phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid2,residm2,rhog2,&
&   rhor2,rmet,rprimd,symaf1,symrc1,symrl1,ucvol,usecprj,useylmgr1,ddk_f,&
&   vtrial,vtrial2,wtk_rbz,xred,ylm,ylm1,ylmgr1,1)

   write(msg,'(a)') ' '//NEW_LINE('A')//'&
&   ---------------------------------'
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,a)')'  The charge sum rule is activated'//NEW_LINE('A')//'&
&   ---------------------------------'
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i4)') ' Z_ion (psp):',Z_kappa
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f12.8)') ' Residual Born effective charge: ',born
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f12.8)') ' Renormalisation: ',renorm
   call wrtout(ab_out,msg,'COLL')
   if (renorm > 0.01 ) then
     write(msg,'(a,a)')'   WARNING: The renormalisation seems large (> 0.01).'//NEW_LINE('A')//'&
&     You might consider increasing the k-point grid.'
     MSG_WARNING(msg)
     call wrtout(ab_out,msg,'COLL')
   end if
   write(msg,'(a)') ' '
   call wrtout(ab_out,msg,'COLL')

   ABI_DEALLOCATE(nvresid2)
   ABI_DEALLOCATE(rhor2)
   ABI_DEALLOCATE(resid2)
   ABI_DEALLOCATE(rhog2)
   ABI_DEALLOCATE(vtrial2)
 end if

 if (iscf_mod>0.or.iscf_mod==-3)  then
   ABI_DEALLOCATE(nvresid1)
 end if

!######################################################################
!Additional steps after SC iterations
!----------------------------------------------------------------------

 call timab(160,1,tsec)

!Compute Dynamic magnetic charges (dmc) in case of rfphon, 
!and magnetic susceptibility in case of rfmagn from first order density
!(results to be comapred to dmc from d2e)
!SPr deb
!if (ipert<=dtset%natom.and.dtset%nspden>=2) then
!
!  mpi_comm_sphgrid=mpi_enreg%comm_fft
!  call mean_fftr(rhor1(:,1),mean_rhor1,nfftf,nfftotf,1,mpi_comm_sphgrid)
!  write(*,*) '   Mean 1st order density: ', mean_rhor1
!  call mean_fftr(rhor1(:,2),mean_rhor1,nfftf,nfftotf,1,mpi_comm_sphgrid)
!  if (dtset%nspden==2) then
!    write(*,*) '        1st order m_z    : ', mean_rhor1
!  else !nspden==4
!    write(*,*) '        1st order m_x    : ', mean_rhor1
!    call mean_fftr(rhor1(:,3),mean_rhor1,nfftf,nfftotf,1,mpi_comm_sphgrid)
!    write(*,*) '        1st order m_y    : ', mean_rhor1
!    call mean_fftr(rhor1(:,4),mean_rhor1,nfftf,nfftotf,1,mpi_comm_sphgrid)
!    write(*,*) '        1st order m_z    : ', mean_rhor1
!  endif
! 
!endif


!Eventually close the dot file, before calling dfpt_nstdy
 if ((ipert==dtset%natom+2.and.sum((dtset%qptn(1:3))**2)<=1.0d-7.and.&
& (dtset%berryopt/=4 .and.dtset%berryopt/= 6.and.dtset%berryopt/= 7.and.&
& dtset%berryopt/=14.and.dtset%berryopt/=16.and.dtset%berryopt/=17)).or.&
& ipert==dtset%natom+10.or.ipert==dtset%natom+11) then
   call wfk_close(ddk_f(1))
 end if
 if ((ipert==dtset%natom+10 .and. idir>3) .or. ipert==dtset%natom+11) then
   call wfk_close(ddk_f(2))
 end if
 if (ipert==dtset%natom+11) then
   call wfk_close(ddk_f(3))
   if(idir>3) call wfk_close(ddk_f(4))
 end if

!Deallocate the no more needed arrays
 if (iscf_mod>0.and.nstep>0) then
   call ab7_mixing_deallocate(mix)
 end if
 if( (nstep>0 .and. iscf_mod>0) .or. iscf_mod==-1 ) then
   ABI_DEALLOCATE(dielinv)
   ABI_DEALLOCATE(susmat)
 end if
 if(allocated(rhorfermi))  then
   ABI_DEALLOCATE(rhorfermi)
 end if
 if(allocated(nhatfermi))  then
   ABI_DEALLOCATE(nhatfermi)
 end if
 if(allocated(pawrhoijfermi))  then
   call pawrhoij_free(pawrhoijfermi)
   ABI_DATATYPE_DEALLOCATE(pawrhoijfermi)
 end if
 if(psps%usepaw==1) then
   if (mk1mem/=0.and.usecprj==1) then
     call pawcprj_free(cprj1)
   end if
   do iatom=1,my_natom
     if (pawfgrtab(iatom)%nhatfr_allocated>0)  then
       ABI_DEALLOCATE(pawfgrtab(iatom)%nhatfr)
     end if
     pawfgrtab(iatom)%nhatfr_allocated=0
   end do
   if (nstep>0.and.iscf_mod>0) then
     do iatom=1,my_natom
       pawrhoij1(iatom)%lmnmix_sz=0
       pawrhoij1(iatom)%use_rhoijres=0
       ABI_DEALLOCATE(pawrhoij1(iatom)%kpawmix)
       ABI_DEALLOCATE(pawrhoij1(iatom)%rhoijres)
     end do
   end if
 end if ! PAW
 ABI_DATATYPE_DEALLOCATE(cprj1)
 ABI_DEALLOCATE(nhat1gr)

 call timab(160,2,tsec)
 call timab(150,1,tsec)

 if (psps%usepaw==0.and.dtset%userie/=919.and. &
& (ipert==dtset%natom+3.or.ipert==dtset%natom+4)) then
   call status(0,dtfil%filstat,iexit,level,'enter dfpt_nselt  ')
   call dfpt_nselt(blkflg,cg,cg1,cplex,&
&   d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,dtset%effmass_free,&
&   gmet,gprimd,gsqcut,idir,&
&   ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,mgfftf,&
&   mkmem,mk1mem,mpert,mpi_enreg,psps%mpsang,mpw,mpw1,&
&   dtset%natom,nband_rbz,nfftf,ngfftf,&
&   nkpt_rbz,nkxc,dtset%nloalg,&
&   npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,&
&   nsym1,dtset%ntypat,occ_rbz,&
&   dtset%paral_kgb,ph1d,dtset%prtbbb,psps,dtset%qptn,rhog,&
&   rhor,rhor1,rmet,rprimd,symrc1,dtset%typat,ucvol,&
&   wtk_rbz,xred,ylm,ylm1,ylmgr,ylmgr1)
 end if

!Use of NSTPAW3 for NCPP (instead of DFPT_NSELT/DFPT_NSTDY) can be forced with userie=919
!MT oct. 2015: this works perfectly on all automatic tests
 if(ipert<=dtset%natom+4)then
   if (psps%usepaw==1.or.dtset%userie==919) then
     call status(0,dtfil%filstat,iexit,level,'enter dfpt_nstpaw ')
     call dfpt_nstpaw(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&
&     eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&
&     kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,ncpgr,&
&     nfftf,ngfftf,nhat,nhat1,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,&
&     nsym1,n3xccc,occkq,occ_rbz,paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,&
&     pawrhoij,pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhog,rhor,rhor1,rmet,rprimd,symaf1,symrc1,&
&     symrl1,ucvol,usecprj,psps%usepaw,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,wtk_rbz,&
&     xccc3d1,xred,ylm,ylm1,ylmgr1)
   else
     call status(0,dtfil%filstat,iexit,level,'enter dfpt_nstdy  ')
     if (dtset%nspden==4) then
       call dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,gmet,&
&       gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,mpert,mpi_enreg,&
&       mpw,mpw1,nattyp,nband_rbz,nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&       dtset%nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,symrc1,ucvol,&
&       wtk_rbz,xred,ylm,ylm1,rhor=rhor)
     else
       call dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,gmet,&
&       gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,mpert,mpi_enreg,&
&       mpw,mpw1,nattyp,nband_rbz,nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&       dtset%nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,symrc1,ucvol,&
&       wtk_rbz,xred,ylm,ylm1)
     end if
   end if
 end if

 call timab(150,2,tsec)
 call timab(160,1,tsec)

!calculate Born effective charge and store it in d2lo
 if ((dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17).and.&
& ipert<=dtset%natom) then
   call dfptff_bec(cg,cg1,dtefield,dtset%natom,d2lo,idir,ipert,dtset%mband,mkmem,&
&   mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,1:dtset%natom)=1
 end if

!calculate dielectric tensor and store it in d2lo
 if ((dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17).and.&
& ipert==dtset%natom+2) then
   call dfptff_die(cg,cg1,dtefield,d2lo,idir,ipert,dtset%mband,mkmem,&
&   mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,dtset%natom+2)=1
 end if

!If SCF convergence was not reached (for nstep>0),
!print a warning to the output file (non-dummy arguments: nstep,
!residm, diffor - infos from tollist have been saved inside )
!Set also the value of conv_retcode
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,0,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)

!Update the content of the header (evolving variables)
 bantot_rbz = sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))
 call hdr_update(hdr,bantot_rbz,etotal,fermie,&
& residm,rprimd,occ_rbz,pawrhoij1,xred,dtset%amu_orig(:,1),&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab )

!Optionally provide output of charge density and/or potential in real space,
!as well as analysis of geometrical factors (bond lengths and bond angles).
!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.

 if(ipert==dtset%natom+5)then
  !debug: write out the vtk first-order density components
  !call appdig(pertcase,dtfil%fnameabo_den,fi1o_vtk)
  !call printmagvtk(mpi_enreg,nspden,nfftf,ngfftf,rhor1,rprimd,adjustl(adjustr(fi1o_vtk)//".vtk"))
  !compute the contributions to susceptibility from different attomic spheres:
   call calcdensph(gmet,mpi_enreg,dtset%natom,nfftf,ngfftf,nspden,&
&   dtset%ntypat,ab_out,dtset%ratsph,rhor1,rprimd,dtset%typat,ucvol,xred,&
&   idir+1,intgden,dentot)

 end if

 if (iwrite_fftdatar(mpi_enreg)) then
   if (dtset%prtden>0) then
     rdwrpaw=0
     call appdig(pertcase,dtfil%fnameabo_den,fi1o)
     ! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
     call fftdatar_write_from_hdr("first_order_density",fi1o,dtset%iomode,hdr,&
     ngfftf,cplex,nfftf,dtset%nspden,rhor1,mpi_enreg)
   end if

   ! first order potentials are always written because the eph code requires them
   ! the files are small (much much smaller that 1WFK, actually we should avoid writing 1WFK)
   rdwrpaw=0
   call appdig(pertcase,dtfil%fnameabo_pot,fi1o)
   ! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
   call fftdatar_write_from_hdr("first_order_potential",fi1o,dtset%iomode,hdr,&
   ngfftf,cplex,nfftf,dtset%nspden,vtrial1,mpi_enreg)

! output files for perturbed potential components: vhartr1,vpsp1,vxc
! NB: only 1 spin for these
   if (dtset%prtvha > 0) then
     rdwrpaw=0
     ABI_ALLOCATE(vhartr1_tmp, (cplex*nfftf, dtset%nspden))
     vhartr1_tmp = zero
     vhartr1_tmp(:,1) = vhartr1(:) 
     call appdig(pertcase,dtfil%fnameabo_vha,fi1o)
     ! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
     call fftdatar_write_from_hdr("first_order_vhartree",fi1o,dtset%iomode,hdr,&
     ngfftf,cplex,nfftf,dtset%nspden,vhartr1_tmp,mpi_enreg)
     ABI_DEALLOCATE(vhartr1_tmp)
   end if
   

! vpsp1 needs to be copied to a temp array - intent(inout) in fftdatar_write_from_hdr though I do not know why
!   if (dtset%prtvpsp > 0) then
!     rdwrpaw=0
!     call appdig(pertcase,dtfil%fnameabo_vpsp,fi1o)
!     ! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
!     call fftdatar_write_from_hdr("first_order_vpsp",fi1o,dtset%iomode,hdr,&
!       ngfftf,cplex,nfftf,1,vpsp1,mpi_enreg)
!   end if
   
   if (dtset%prtvxc > 0) then
     rdwrpaw=0
     call appdig(pertcase,dtfil%fnameabo_vxc,fi1o)
     ! TODO: should we write pawrhoij1 or pawrhoij. Note that ioarr writes hdr%pawrhoij
     call fftdatar_write_from_hdr("first_order_vxc",fi1o,dtset%iomode,hdr,&
     ngfftf,cplex,nfftf,dtset%nspden,vxc1,mpi_enreg)
   end if


   ! Add rhog1(G=0) to file
   if (mpi_enreg%me_g0 == 1) then
     if (dtset%iomode == IO_MODE_ETSF) then
#ifdef HAVE_NETCDF
       NCF_CHECK(nctk_open_modify(ncid, nctk_ncify(fi1o), xmpi_comm_self))
       ncerr = nctk_def_one_array(ncid, nctkarr_t('rhog1_g0', "dp", "two"), varid=varid)
       NCF_CHECK(ncerr)
       NCF_CHECK(nctk_set_datamode(ncid))
       NCF_CHECK(nf90_put_var(ncid, varid, rhog1(:,1)))
       NCF_CHECK(nf90_close(ncid))
#endif
     else
       ! Handle Fortran files.
       if (open_file(fi1o, msg, newunit=ncid, form='unformatted', status='old', action="readwrite") /= 0) then
         MSG_ERROR(msg)
       end if
       if (fort_denpot_skip(ncid, msg) /= 0) MSG_ERROR(msg)
       write(ncid) rhog1(:,1)
       close(ncid)
     end if
   end if

 end if ! iwrite_fftdatar(mpi_enreg)

!All procs waiting here...
 if(mpi_enreg%paral_kgb==1)then
   call timab(61,1,tsec)
   call xmpi_barrier(spaceComm)
   call timab(61,2,tsec)
 end if

!Deallocate arrays
 ABI_DEALLOCATE(fcart)
 ABI_DEALLOCATE(vtrial1)
 ABI_DEALLOCATE(vhartr1)
 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(pwindall)
 ABI_DEALLOCATE(qmat)
 if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then
   call destroy_efield(dtefield)
   if(allocated(mpi_enreg%kpt_loc2ibz_sp))  then
     ABI_DEALLOCATE(mpi_enreg%kpt_loc2ibz_sp)
   end if
 end if

 if(psps%usepaw==1) then
   call paw_an_free(paw_an1)
   call paw_ij_free(paw_ij1)
 end if
 ABI_DATATYPE_DEALLOCATE(paw_an1)
 ABI_DATATYPE_DEALLOCATE(paw_ij1)
 ABI_DEALLOCATE(nhat1)

 call status(0,dtfil%filstat,iexit,level,'exit')

 call timab(160,2,tsec)
 call timab(120,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_scfcv
!!***

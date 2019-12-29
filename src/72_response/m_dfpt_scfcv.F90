!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_dfpt_scfcv
!! NAME
!!  m_dfpt_scfcv
!!
!! FUNCTION
!!
!! COPYRIGHT
!!  Copyright (C) 1999-2019 ABINIT group (XG, DRH, MB, XW, MT, SPr, XW, MV, MM, AR)
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

module m_dfpt_scfcv

 use defs_basis
 use m_ab7_mixing
 use m_efield
 use m_errors
 use m_dtset
 use m_abicore
 use m_wfk
 use m_wffile
 use m_xmpi
 use m_nctk
 use m_hdr
 use m_dtfil
 use m_hamiltonian
#ifdef HAVE_NETCDF
 use netcdf
#endif

 use defs_datatypes, only : pseudopotential_type
 use defs_abitypes, only : MPI_type
 use m_cgtools,  only : mean_fftr, overlap_g, dotprod_vn, dotprod_vn, dotprod_g
 use m_fstrings, only : int2char4, sjoin
 use m_geometry, only : metric, stresssym
 use m_time,     only : abi_wtime, sec2str, timab
 use m_io_tools, only : open_file, file_exists, get_unit, iomode_from_fname
 use m_exit,     only : get_start_time, have_timelimit_in, get_timelimit, enable_timelimit_in
 use m_mpinfo
 use m_kg,       only : getcut, mkkin, kpgstr, mkkpg
 use m_fft,      only : fftpac, fourdp
 use m_symtk,     only : mati3inv
 use m_dynmat,    only : dfpt_sygra
 use m_occ,         only : occeig
 use m_paw_mkrho,   only : pawmkrho
 use m_mkffnl,      only : mkffnl
 use m_getgh1c,     only : getgh1c
 use m_dfpt_mkrho,  only : dfpt_accrho
 use m_nonlop,      only : nonlop
 use m_ioarr,    only : ioarr, fftdatar_write_from_hdr, fort_denpot_skip
 use m_pawang,   only : pawang_type
 use m_pawrad,   only : pawrad_type
 use m_pawtab,   only : pawtab_type
 use m_paw_an,   only : paw_an_type, paw_an_init, paw_an_free, paw_an_nullify, paw_an_reset_flags
 use m_paw_ij,   only : paw_ij_type, paw_ij_init, paw_ij_free, paw_ij_nullify, paw_ij_reset_flags
 use m_pawfgrtab,only : pawfgrtab_type
 use m_pawrhoij,    only : pawrhoij_type, pawrhoij_init_unpacked, pawrhoij_gather, pawrhoij_filter, &
                           pawrhoij_alloc, pawrhoij_free, pawrhoij_nullify, &
                           pawrhoij_free_unpacked, pawrhoij_mpisum_unpacked, pawrhoij_inquire_dim
 use m_pawcprj,  only : pawcprj_type, pawcprj_alloc, pawcprj_get, pawcprj_copy, pawcprj_axpby, pawcprj_free, pawcprj_getdim
 use m_pawdij,   only : pawdij, pawdijfr, symdij
 use m_pawfgr,   only : pawfgr_type
 use m_paw_denpot,  only : pawdenpot
 use m_paw_dfpt,    only : pawdfptenergy
 use m_paw_nhat,    only : pawmknhat,pawnhatfr
 use m_rf2,         only : rf2_getidirs
 use m_dens,        only : calcdenmagsph, prtdenmagsph
 use m_dfpt_fef,    only : dfptff_initberry, qmatrix, dfptff_edie, dfptff_ebp, dfptff_die, dfptff_bec
 use m_dfpt_vtorho, only : dfpt_vtorho
 use m_paral_atom,  only : get_my_atmtab, free_my_atmtab
 use m_common,      only : scprqt
 use m_prcref,      only : moddiel
 use m_dfpt_rhotov, only : dfpt_rhotov
 use m_dfpt_mkvxc,    only : dfpt_mkvxc, dfpt_mkvxc_noncoll
 use m_dfpt_mkvxcstr, only : dfpt_mkvxcstr
 use m_mklocl,     only : dfpt_vlocal, vlocalstr
 use m_dfpt_nstwf,   only : dfpt_nstpaw, dfpt_nstwf
 use m_mkcore,         only : dfpt_mkcore
 use m_spacepar,   only : hartrestr, symrhg

 implicit none

 private
!!***

 public :: dfpt_scfcv
!!***

contains
!!***

!!****f* ABINIT/dfpt_scfcv
!! NAME
!! dfpt_scfcv
!!
!! FUNCTION
!! Conducts set of passes or overall iterations of preconditioned
!! conjugate gradient algorithm to converge wavefunctions to
!! optimum and optionally to compute mixed derivatives of energy.
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
!!  kxc(nfftf,nkxc)=exchange and correlation kernel (see rhotoxc.f)
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
!!  tnons1(3,nsym1)=non-symmorphic translations
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
!!        PAW only - Eq(79) and Eq(80) of PRB 78, 035105 (2008) [[cite:Audouze2008]]
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
!!      calcdenmagsph,destroy_efield,dfpt_etot,dfpt_newvtr,dfpt_nselt,dfpt_nstdy
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

subroutine dfpt_scfcv(atindx,blkflg,cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cpus,&
&  dielt,dim_eig2rf,doccde_rbz,docckqde,dtfil,dtset,&
&  d2bbb,d2lo,d2nl,d2ovl,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,&
&  ehart01,ehart1,eigenq,eigen0,eigen1,eii,ek0,ek1,eloc0,elpsp1,&
&  enl0,enl1,eovl1,epaw1,etotal,evdw,exc1,fermie,gh0c1_set,gh1c_set,hdr,idir,indkpt1,&
&  indsy1,initialized,ipert,irrzon1,istwfk_rbz,&
&  kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,&
&  mpert,mpi_enreg,mpw,mpw1,mpw1_mq,my_natom,nattyp,nband_rbz,ncpgr,&
&  nfftf,ngfftf,nhat,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&  nsym1,n3xccc,occkq,occ_rbz,&
&  paw_an,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoij,pawrhoij1,pawtab,&
&  pertcase,phnons1,ph1d,ph1df,&
&  prtbbb,psps,qphon,resid,residm,rhog,rhog1,&
&  rhor,rhor1,rprimd,symaf1,symrc1,symrl1,tnons1,&
&  usecprj,useylmgr,useylmgr1,ddk_f,vpsp1,vtrial,vxc,&
&  wtk_rbz,xccc3d1,xred,ylm,ylm1,ylmgr,ylmgr1,zeff,conv_retcode,&
&  kramers_deg,&
&  cg_mq,cg1_mq,cg1_active_mq,docckde_mq,eigen_mq,eigen1_mq,gh0c1_set_mq,gh1c_set_mq,&
&  kg1_mq,npwar1_mq,occk_mq,resid_mq,residm_mq,rhog1_pq,rhog1_mq,rhor1_pq,rhor1_mq)

!Arguments ------------------------------------
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
 integer,intent(in) :: cplex,dim_eig2rf,idir,ipert,mgfftf,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mpert,mpw,mpw1,my_natom,n3xccc,ncpgr,nfftf
 integer,intent(in) :: mpw1_mq !-q duplicate
 integer,intent(in) :: nkpt,nkpt_rbz,nkxc,nspden
 integer,intent(in) :: nsym1,pertcase,prtbbb,usecprj,useylmgr,useylmgr1
 logical,intent(in) :: kramers_deg
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
 integer,optional,intent(in) :: npwar1_mq(nkpt_rbz)     !-q duplicate
 integer,optional,intent(in) :: kg1_mq(3,mpw1_mq*mk1mem)!
 integer,intent(in) :: symaf1(nsym1),symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 integer,intent(out) :: conv_retcode
 real(dp),intent(in) :: cpus,eew,efrhar,efrkin,efrloc,efrnl,efrx1,efrx2,eii
 real(dp),intent(out) :: eberry,edocc,eeig0,ehart01,ehart1,ek0,ek1,eloc0,elpsp1,enl0
 real(dp),intent(out) :: enl1,eovl1,epaw1,etotal,evdw,exc1,residm
 real(dp),optional,intent(out) :: residm_mq       !-q duplicate
 real(dp),intent(inout) :: fermie
 real(dp),intent(in) :: qphon(3)
! nfft**(1-1/nsym1) is 1 if nsym1==1, and nfft otherwise
 integer,intent(in) :: ngfftf(18)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband_mem*mkmem*dtset%nsppol)
 real(dp),intent(inout) :: cg1(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol)
 real(dp),intent(out) :: cg1_active(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh1c_set(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(out) :: gh0c1_set(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf)
 real(dp),intent(in)  :: cgq(2,mpw1*dtset%nspinor*dtset%mband_mem*mkqmem*dtset%nsppol)
 real(dp),optional,intent(inout) :: cg1_mq(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol)                  !start -q duplicates
 real(dp),optional,intent(out)   :: cg1_active_mq(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf)!
 real(dp),optional,intent(out)   :: gh1c_set_mq(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf)  !
 real(dp),optional,intent(out)   :: gh0c1_set_mq(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mk1mem*dtset%nsppol*dim_eig2rf) !
 real(dp),optional,intent(in)    :: cg_mq(2,mpw1_mq*dtset%nspinor*dtset%mband_mem*mkqmem*dtset%nsppol)                   !
 real(dp),optional,intent(in)    :: eigen_mq(dtset%mband*nkpt_rbz*dtset%nsppol)                                      !
 real(dp),optional,intent(in)    :: docckde_mq(dtset%mband*nkpt_rbz*dtset%nsppol)                                    !
 real(dp),optional,intent(out)   :: eigen1_mq(2*dtset%mband*dtset%mband*nkpt_rbz*dtset%nsppol)                       !
 real(dp),optional,intent(in)    :: occk_mq(dtset%mband*nkpt_rbz*dtset%nsppol)                                       !
 real(dp),optional,intent(out)   :: resid_mq(dtset%mband*nkpt_rbz*nspden)                                            !end
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
 real(dp),optional,intent(inout) :: rhog1_pq(2,nfftf),rhor1_pq(cplex*nfftf,nspden)                                 !+q/-q duplicates
 real(dp),optional,intent(inout) :: rhog1_mq(2,nfftf),rhor1_mq(cplex*nfftf,nspden)                                 !
 real(dp),intent(in) :: tnons1(3,nsym1)
 real(dp),target,intent(in) :: vtrial(nfftf,nspden)
 real(dp),intent(in) :: vpsp1(cplex*nfftf),vxc(nfftf,nspden)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xccc3d1(cplex*n3xccc)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3+6*((ipert-dtset%natom)/10),psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(in) :: zeff(3,3,dtset%natom)
!TODO MJV : PAW
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
 integer :: has_dijfr,has_diju,iatom,ider,idir_dum,idir_paw1,ierr,errid,denpot
 integer :: iprcel,iscf10_mod,iscf_mod,ispden,ispmix
 integer :: istep,istep_fock_outer,istep_mix,itypat,izero,me,mgfftdiel,mvdum !lmn2_size,
 integer :: nfftdiel,nfftmix,nfftotf,nhat1grdim,npawmix,npwdiel,nspden_rhoij,nstep,nzlmopt
 integer :: optene,optfr,option,optres,prtfor,qphase_rhoij,quit,quit_sum,qzero
 integer :: my_quit,quitsum_request,timelimit_exit,varid,ncerr,ncid
 integer ABI_ASYNC :: quitsum_async
 integer :: rdwrpaw,spaceComm,sz1,sz2,usexcnhat,Z_kappa
 integer :: dbl_nnsclo_mq,ifft !-q duplicate for dbl_nnsclo
!integer :: pqmq ! pqmq = indicator for potential mixing
 logical :: need_fermie1,nmxc,paral_atom,use_nhat_gga
 real(dp) :: wtime_step,now,prev
 real(dp) :: born,born_bar,boxcut,deltae,diffor,diel_q,dum,ecut,ecutf,elast
 real(dp) :: epawdc1_dum,evar,fe1fixed,fermie1,gsqcut,qphon_norm,maxfor,renorm,res2,res3,residm2
 real(dp) :: ucvol,vxcavg,elmag1
 real(dp) :: res2_mq,fe1fixed_mq,elast_mq
 real(dp) :: eberry_mq,edocc_mq,eeig0_mq,ehart01_mq,ehart1_mq,ek0_mq,ek1_mq,eloc0_mq,elpsp1_mq,enl0_mq
 real(dp) :: enl1_mq,eovl1_mq,epaw1_mq,exc1_mq,fermie1_mq,deltae_mq,elmag1_mq
 character(len=500) :: msg
 character(len=500),parameter :: MY_NAME="dfpt_scfcv"
 character(len=fnlen) :: fi1o
!character(len=fnlen) :: fi1o_vtk
 integer  :: prtopt
 type(ab7_mixing_object) :: mix
 type(efield_type) :: dtefield
!arrays
 integer :: ngfftmix(18)
 integer,allocatable :: dimcprj(:),pwindall(:,:,:)
 integer,pointer :: my_atmtab(:)
 real(dp) :: dielar(7)
 real(dp) :: favg(3),gmet(3,3),gprimd(3,3),q_cart(3),qphon2(3),qred2cart(3,3)
 real(dp) :: rhomag(2,nspden),rmet(3,3),tollist(12),tsec(2)
 real(dp) :: zeff_red(3),zeff_bar(3,3)
 real(dp) :: intgden(dtset%nspden,dtset%natom),dentot(dtset%nspden)
!real(dp) :: zdmc_red(3),zdmc_bar(3,3),mean_rhor1(1) !dynamic magnetic charges and mean density
 real(dp),allocatable :: dielinv(:,:,:,:,:)
 real(dp),allocatable :: fcart(:,:),nhat1(:,:),nhat1gr(:,:,:),nhatfermi(:,:),nvresid1(:,:),nvresid2(:,:)
 real(dp),allocatable :: qmat(:,:,:,:,:,:),resid2(:),rhog2(:,:),rhor2(:,:),rhorfermi(:,:)
 real(dp),allocatable :: susmat(:,:,:,:,:),vhartr1(:),vxc1(:,:)
 real(dp),allocatable :: vhartr1_tmp(:,:)
 real(dp),allocatable,target :: vtrial1(:,:),vtrial2(:,:)
 real(dp),allocatable :: vtrial1_pq(:,:),vtrial1_mq(:,:),rhorfermi_mq(:,:)
 real(dp),allocatable :: nvresid1_mq(:,:),vxc1_mq(:,:),vhartr1_mq(:)
 real(dp),pointer :: vtrial1_tmp(:,:)
 type(pawcprj_type),allocatable :: cprj1(:,:)
 type(paw_an_type),allocatable :: paw_an1(:)
 type(paw_ij_type),allocatable :: paw_ij1(:)
 type(pawrhoij_type),allocatable :: pawrhoijfermi(:)

! *********************************************************************

 DBG_ENTER("COLL")

 call timab(120,1,tsec)
 call timab(154,1,tsec)

 ! intel 18 really needs this to be initialized
 maxfor = zero

 ! enable time limit handler if not done in callers.
 if (enable_timelimit_in(MY_NAME) == MY_NAME) then
   write(std_out,*)"Enabling timelimit check in function: ",trim(MY_NAME)," with timelimit: ",trim(sec2str(get_timelimit()))
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
& (ipert<=dtset%natom.or.ipert==dtset%natom+3.or.ipert==dtset%natom+4.or.ipert==dtset%natom+5))

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
 nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
 usexcnhat=0
!This might be taken away later
 edocc=zero ; eeig0=zero ; ehart01=zero ; ehart1=zero ; ek0=zero ; ek1=zero
 eloc0=zero ; elpsp1=zero ; enl0=zero ; enl1=zero ; eovl1=zero; exc1=zero
 deltae=zero ; fermie1=zero ; epaw1=zero ; eberry=zero ; elmag1=zero
 elast_mq=zero
 dbl_nnsclo_mq=0
!This might be taken away later
 edocc_mq=zero ; eeig0_mq=zero ; ehart01_mq=zero ; ehart1_mq=zero ; ek0_mq=zero ; ek1_mq=zero
 eloc0_mq=zero ; elpsp1_mq=zero ; enl0_mq=zero ; enl1_mq=zero ; eovl1_mq=zero; exc1_mq=zero
 deltae_mq=zero ; fermie1_mq=zero ; epaw1_mq=zero ; eberry_mq=zero ; elmag1_mq=zero
 res2_mq=zero

!Examine tolerance criteria, and eventually  print a line to the output
!file (with choice=1, the only non-dummy arguments of scprqt are
!nstep, tollist and iscf - still, diffor,res2,prtfor,fcart are here initialized to 0)
 choice=1 ; prtfor=0 ; diffor=zero ; res2=zero
 ABI_ALLOCATE(fcart,(3,dtset%natom))

!At present, no double loop
 istep_mix=1 ; istep_fock_outer=1

 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,istep_fock_outer,istep_mix,kpt_rbz,maxfor,&
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
!TODO MJV : PAW
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
   has_diju=merge(0,1,dtset%usepawu==0)
   call paw_an_init(paw_an1,dtset%natom,dtset%ntypat,0,0,dtset%nspden,&
&   cplex,dtset%pawxcdev,dtset%typat,pawang,pawtab,has_vxc=1,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   call paw_ij_init(paw_ij1,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,0,dtset%natom,&
&   dtset%ntypat,dtset%typat,pawtab,&
&   has_dij=1,has_dijhartree=1,has_dijfr=has_dijfr,has_dijU=has_diju,&
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
 if(.not.kramers_deg) then
   ABI_ALLOCATE(vhartr1_mq,(cplex*nfftf))
   ABI_ALLOCATE(vtrial1_pq,(cplex*nfftf,nspden))
   ABI_ALLOCATE(vtrial1_mq,(cplex*nfftf,nspden))
 end if
! TODO: for non collinear case this should always be nspden, in NCPP case as well!!!
 ABI_ALLOCATE(vxc1,(cplex*nfftf,nspden*(1-usexcnhat))) ! Not always needed
 vtrial1_tmp => vtrial1   ! this is to avoid errors when vtrial1_tmp is unused

 if (.not.kramers_deg) then
   ABI_ALLOCATE(vxc1_mq,(cplex*nfftf,nspden*(1-usexcnhat)))
 end if

!Several parameters and arrays for the SCF mixing:
!These arrays are needed only in the self-consistent case
 if (iscf_mod>0.or.iscf_mod==-3) then
   ABI_ALLOCATE(nvresid1,(cplex*nfftf,dtset%nspden))
   if (nstep==0) nvresid1=zero
   if ((dtset%getddb .ne. 0 .or. dtset%irdddb .ne.0) .and. qzero .ne. 1) then
     ABI_ALLOCATE(nvresid2,(cplex*nfftf,dtset%nspden))
     if (nstep==0) nvresid2=zero
   end if
   if (.not.kramers_deg) then
     ABI_ALLOCATE(nvresid1_mq,(cplex*nfftf,dtset%nspden))
     if (nstep==0) nvresid1_mq=zero
   end if
 else
   ABI_ALLOCATE(nvresid1,(0,0))
   if(.not.kramers_deg) then
     ABI_ALLOCATE(nvresid1_mq,(0,0))
   end if
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
       pawrhoij1(iatom)%use_rhoijres=1
       sz1=pawrhoij1(iatom)%cplex_rhoij*pawrhoij1(iatom)%qphase*pawrhoij1(iatom)%lmn2_size
       sz2=pawrhoij1(iatom)%nspden
       ABI_ALLOCATE(pawrhoij1(iatom)%rhoijres,(sz1,sz2))
       do ispden=1,pawrhoij1(iatom)%nspden
         pawrhoij1(iatom)%rhoijres(:,ispden)=zero
       end do
       ABI_ALLOCATE(pawrhoij1(iatom)%kpawmix,(pawtab(itypat)%lmnmix_sz))
       pawrhoij1(iatom)%lmnmix_sz=pawtab(itypat)%lmnmix_sz
       pawrhoij1(iatom)%kpawmix=pawtab(itypat)%kmix
       npawmix=npawmix+pawrhoij1(iatom)%nspden*pawtab(itypat)%lmnmix_sz &
&                     *pawrhoij1(iatom)%cplex_rhoij*pawrhoij1(iatom)%qphase
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
   call qmatrix(cg,dtefield,qmat,mpi_enreg,mpw,mpw1,mkmem,dtset%mband,dtset%mband_mem,&
&    npwarr,nkpt,dtset%nspinor,dtset%nsppol,pwindall)
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
     call wrtout(std_out,sjoin(" dfpt_scfcv: previous iteration took ",sec2str(wtime_step)))

     if (have_timelimit_in(MY_NAME)) then
       if (istep > 2) then
         call xmpi_wait(quitsum_request,ierr)
         if (quitsum_async > 0) then
           write(msg,"(3a)")" Approaching time limit ",trim(sec2str(get_timelimit())),". Will exit istep loop in dfpt_scfcv."
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
       call pawdijfr(gprimd,idir_paw1,ipert,my_natom,dtset%natom,nfftf,ngfftf,nspden,dtset%nsppol,&
&       psps%ntypat,optfr,paw_ij1,pawang,pawfgrtab,pawrad,pawtab,cplex,qphon,&
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
         call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,1, ngfftf,0)
       end if
       call timab(564,2,tsec)
     end if
!    Set initial guess for 1st-order potential
!    ----------------------------------------------------------------------
     option=1;optene=0;if (iscf_mod==-2) optene=1
     call dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,elmag1,gsqcut,idir,ipert,&
&     dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1gr,nhat1grdim,&
&     nkxc,nspden,n3xccc,nmxc,optene,option,dtset%qptn,&
&     rhog,rhog1,rhor,rhor1,rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,&
&     nvresid1,res2,vtrial1,vxc,vxc1,xccc3d1,dtset%ixcrot)

     if(.not.kramers_deg) then
       vtrial1_pq=vtrial1 !save trial potential at +q
       !rhor1_mq=rhor1
       !rhog1_mq=rhog1
       !get initial guess for vtrial1 at -q
       do ifft=1,nfftf
         vtrial1_mq(2*ifft-1,1)=+vtrial1(2*ifft-1,1)
         vtrial1_mq(2*ifft-1,2)=+vtrial1(2*ifft-1,2)
         vtrial1_mq(2*ifft  ,1)=-vtrial1(2*ifft  ,1)
         vtrial1_mq(2*ifft  ,2)=-vtrial1(2*ifft  ,2)
         vtrial1_mq(2*ifft-1,3)= vtrial1(2*ifft  ,4) !Re[V^12]
         vtrial1_mq(2*ifft  ,3)= vtrial1(2*ifft-1,4) !Im[V^12],see definition of v(:,4) cplex=2 case
         vtrial1_mq(2*ifft  ,4)= vtrial1(2*ifft-1,3) !Re[V^21]=Re[V^12]
         vtrial1_mq(2*ifft-1,4)= vtrial1(2*ifft  ,3) !Re[V^21]=Re[V^12]
       end do
     end if

!    For Q=0 and metallic occupation, initialize quantities needed to
!    compute the first-order Fermi energy
!    ----------------------------------------------------------------------
     if (need_fermie1) then
       ABI_ALLOCATE(rhorfermi,(cplex*nfftf,nspden))
       if(.not.kramers_deg) then
         ABI_ALLOCATE(rhorfermi_mq,(cplex*nfftf,nspden))
       end if
       if (psps%usepaw==1.and.usexcnhat==0) then
         ABI_ALLOCATE(nhatfermi,(cplex*nfftf,nspden))
       else
         ABI_ALLOCATE(nhatfermi,(0,0))
       end if
       ABI_DATATYPE_ALLOCATE(pawrhoijfermi,(my_natom*psps%usepaw))
       if (psps%usepaw==1) then
         !Q phase should be 1 because Q=0
         call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                              nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
         call pawrhoij_alloc(pawrhoijfermi,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&         dtset%nsppol,dtset%typat,pawtab=pawtab,mpi_atmtab=mpi_enreg%my_atmtab,&
&         comm_atom=mpi_enreg%comm_atom)
       end if

       call dfpt_rhofermi(cg,cgq,cplex,cprj,cprjq,&
&       doccde_rbz,docckqde,dtfil,dtset,eigenq,eigen0,eigen1,fe1fixed,gmet,gprimd,idir,&
&       indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,dtset%mband_mem,mkmem,mkqmem,mk1mem,mpi_enreg,&
&       mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1,&
&       nspden,dtset%nsppol,nsym1,occkq,occ_rbz,&
&       paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&
&       phnons1,ph1d,dtset%prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,&
&       ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)
       if (.not.kramers_deg) then
         call dfpt_rhofermi(cg,cg_mq,cplex,cprj,cprjq,&
&         doccde_rbz,docckde_mq,dtfil,dtset,eigen_mq,eigen0,eigen1_mq,fe1fixed_mq,gmet,gprimd,idir,&
&         indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1_mq,kpt_rbz,dtset%mband,dtset%mband_mem,mkmem,mkqmem,mk1mem,mpi_enreg,&
&         mpw,mpw1_mq,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1_mq,&
&         nspden,dtset%nsppol,nsym1,occk_mq,occ_rbz,&
&         paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&
&         phnons1,ph1d,dtset%prtvol,psps,rhorfermi_mq,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,&
&         ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)
       end if

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
     call pawdenpot(dum,epaw1,epawdc1_dum,ipert,dtset%ixc,my_natom,dtset%natom,&
&     dtset%nspden,psps%ntypat,dtset%nucdipmom,nzlmopt,option,paw_an1,paw_an,paw_ij1,pawang,&
&     dtset%pawprtvol,pawrad,pawrhoij1,dtset%pawspnorb,pawtab,dtset%pawxcdev,&
&     dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp, &
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)

!    First-order Dij computation
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

     call symdij(gprimd,indsy1,ipert,my_natom,dtset%natom,nsym1,psps%ntypat,0,&
&     paw_ij1,pawang1,dtset%pawprtvol,pawtab,rprimd,symaf1,symrc1, &
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,&
&     qphon=qphon)
     call timab(561,2,tsec)
   end if ! end usepaw section

!  ######################################################################
!  The following steps are done only when nstep>0
!  ----------------------------------------------------------------------

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
     if (.not.kramers_deg) then
       !fermie1_mq is updated as well at "-q"
       call newfermie1(cplex,fermie1_mq,fe1fixed_mq,ipert,istep,dtset%ixc,my_natom,dtset%natom,&
&       nfftf,nfftotf,nhatfermi,nspden,dtset%ntypat,dtset%occopt,paw_an,paw_an1,paw_ij1,pawang,&
&       dtset%pawnzlm,pawrad,pawrhoij1,pawrhoijfermi,pawtab,dtset%pawxcdev,&
&       dtset%prtvol,rhorfermi_mq,ucvol,psps%usepaw,usexcnhat,vtrial1_mq,vxc1_mq,dtset%xclevel,&
&       mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
     end if
   end if

!  No need to continue and call dfpt_vtorho, when nstep==0
   if(nstep==0) exit

!  #######################e1magh###############################################
!  Compute the 1st-order density rho1 from the 1st-order trial potential
!  ----------------------------------------------------------------------

   call dfpt_vtorho(cg,cgq,cg1,cg1_active,cplex,cprj,cprjq,cprj1,&
&   dbl_nnsclo,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,dtset%qptn,edocc,&
&   eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,fermie1,gh0c1_set,gh1c_set,&
&   gmet,gprimd,idir,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,dtset%mband_mem,&
&   mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,&
&   nhat1,nkpt_rbz,npwarr,npwar1,res2,nspden,dtset%nsppol,nsym1,dtset%ntypat,nvresid1,&
&   occkq,occ_rbz,optres,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,&
&   pawrhoij1,pawtab,phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid,residm,rhog1,&
&   rhor1,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,ucvol,usecprj,useylmgr1,ddk_f,&
&   vtrial,vtrial1,wtk_rbz,xred,ylm,ylm1,ylmgr1)
   if (.not.kramers_deg) then
     rhor1_pq=rhor1 !at this stage rhor1_pq contains only one term of the 1st order density at +q
     rhog1_pq=rhog1 !same for rhog1_pq
     !get the second term related to 1st order wf at -q
     call dfpt_vtorho(cg,cg_mq,cg1_mq,cg1_active_mq,cplex,cprj,cprjq,cprj1,&
&     dbl_nnsclo_mq,dim_eig2rf,doccde_rbz,docckde_mq,dtefield,dtfil,dtset,-dtset%qptn,edocc_mq,&
&     eeig0_mq,eigen_mq,eigen0,eigen1_mq,ek0_mq,ek1_mq,eloc0_mq,enl0_mq,enl1_mq,fermie1_mq,gh0c1_set_mq,gh1c_set_mq,&
&     gmet,gprimd,idir,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1_mq,kpt_rbz,dtset%mband,dtset%mband_mem,&
&     mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1_mq,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,&
&     nhat1,nkpt_rbz,npwarr,npwar1_mq,res2_mq,nspden,dtset%nsppol,nsym1,dtset%ntypat,nvresid1_mq,&
&     occk_mq,occ_rbz,optres,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,&
&     pawrhoij1,pawtab,phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid_mq,residm_mq,rhog1_mq,&
&     rhor1_mq,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,ucvol,usecprj,useylmgr1,ddk_f,&
&     vtrial,vtrial1_mq,wtk_rbz,xred,ylm,ylm1,ylmgr1)
     !reconstruct the +q and -q densities, this might bug if fft parallelization is used, todo...
     do ifft=1,nfftf
       rhor1_pq(2*ifft-1,:) = half*(rhor1(2*ifft-1,:)+rhor1_mq(2*ifft-1,:))
       rhor1_pq(2*ifft  ,:) = half*(rhor1(2*ifft  ,:)-rhor1_mq(2*ifft  ,:))
       rhor1_mq(2*ifft-1,:) = rhor1_pq(2*ifft-1,:)
       rhor1_mq(2*ifft  ,:) =-rhor1_pq(2*ifft  ,:)
     end do
     rhor1=rhor1_pq
     call fourdp(cplex,rhog1,rhor1(:,1),-1,mpi_enreg,nfftf,1, ngfftf, 0)
     call fourdp(cplex,rhog1_mq,rhor1_mq(:,1),-1,mpi_enreg,nfftf,1, ngfftf, 0)
   end if

   if (dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
&   dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17) then

!    calculate \Omega E \cdot P term
     if (ipert<=dtset%natom) then
!      phonon perturbation
       call  dfptff_ebp(cg,cg1,dtefield,eberry,dtset%mband,dtset%mband_mem,mkmem,&
&       mpi_enreg,mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat)
     else if (ipert==dtset%natom+2) then
!      electric field perturbation
       call  dfptff_edie(cg,cg1,dtefield,eberry,idir,dtset%mband,dtset%mband_mem,mkmem,&
&       mpi_enreg,mpw,mpw1,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
     end if
   end if

!  SPr: don't remove the following comments for debugging
!  call calcdenmagsph(gmet,mpi_enreg,dtset%natom,nfftf,ngfftf,nspden,&
!&   dtset%ntypat,dtset%ratsm,dtset%ratsph,rhor1,rprimd,dtset%typat,ucvol,xred,&
!&   idir+1,cplex,intgden=intgden,rhomag=rhomag)
!  call  prtdenmagsph(cplex,intgden,dtset%natom,nspden,dtset%ntypat,ab_out,idir+1,dtset%ratsm,dtset%ratsph,rhomag,dtset%typat)

!     write(*,*) ' n ( 1,2)',intgden(1,1),' ',intgden(1,2)
!     write(*,*) ' mx( 1,2)',intgden(2,1),' ',intgden(2,2)
!     write(*,*) ' my( 1,2)',intgden(3,1),' ',intgden(3,2)
!     write(*,*) ' mz( 1,2)',intgden(4,1),' ',intgden(4,2)
!  call dfpt_etot(dtset%berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,&
!&     efrloc,efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
!&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,elmag1,ipert,dtset%natom,optene)
!     write(*,*) 'SPr: ek1=',ek1,'  exc1=',exc1,' elmag1=',elmag
!  if (ipert==dtset%natom+5) then
!  !calculate 1st order magnetic potential contribution to the energy
!    call dfpt_e1mag(e1mag,rhor1,rhog1);
!  endif

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
&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,ipert,dtset%natom,optene)
     call timab(152,1,tsec)
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,istep_fock_outer,istep_mix,kpt_rbz,maxfor,&
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
!TODO MJV: check for paralbd
     call dfpt_rhotov(cplex,ehart01,ehart1,elpsp1,exc1,elmag1,gsqcut,idir,ipert,&
&     dtset%ixc,kxc,mpi_enreg,dtset%natom,nfftf,ngfftf,nhat,nhat1,nhat1gr,nhat1grdim,nkxc,&
&     nspden,n3xccc,nmxc,optene,optres,dtset%qptn,rhog,rhog1,rhor,rhor1,&
&     rprimd,ucvol,psps%usepaw,usexcnhat,vhartr1,vpsp1,nvresid1,res2,vtrial1,vxc,vxc1,xccc3d1,dtset%ixcrot)
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
!&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,elmag1,ipert,dtset%natom,optene)
!    !debug: compute the d2E/d-qd+q energy, should be equal to the one from previous line
!    if(.not.kramers_deg) then
!      call dfpt_etot(dtset%berryopt,deltae_mq,eberry_mq,edocc_mq,eeig0_mq,eew,efrhar,efrkin,&
!&       efrloc,efrnl,efrx1,efrx2,ehart1_mq,ek0_mq,ek1_mq,eii,elast_mq,eloc0_mq,elpsp1_mq,&
!&       enl0_mq,enl1_mq,epaw1_mq,etotal_mq,evar_mq,evdw,exc1_mq,elmag1_mq,ipert,dtset%natom,optene)
!     end if
&     enl0,enl1,epaw1,etotal,evar,evdw,exc1,ipert,dtset%natom,optene)

     call timab(152,1,tsec)
     choice=2
     call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
&     etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
&     1,iscf_mod,istep,istep_fock_outer,istep_mix,kpt_rbz,maxfor,&
&     mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
&     nstep,occ_rbz,0,prtfor,0,&
&     quit,res2,resid,residm,response,&
&     tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)
!     !debug: print the information about residuals at "-q"
!     if(.not.kramers_deg) then
!       call scprqt(choice,cpus,deltae_mq,diffor,dtset,eigen0,&
!&       etotal_mq,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
!&       1,iscf_mod,istep,istep_fock_outer,istep_mix,kpt_rbz,maxfor,&
!&       mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
!&       nstep,occ_rbz,0,prtfor,0,&
!&       quit,res2_mq,resid_mq,residm_mq,response,&
!&       tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)
!     endif
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
       if (.not.kramers_deg) then
       !same problem as with density reconstruction, TODO proper fft parallelization...
         do ifft=1,nfftf
           vtrial1_mq(2*ifft-1,1)=+vtrial1(2*ifft-1,1)
           vtrial1_mq(2*ifft-1,2)=+vtrial1(2*ifft-1,2)
           vtrial1_mq(2*ifft  ,1)=-vtrial1(2*ifft  ,1)
           vtrial1_mq(2*ifft  ,2)=-vtrial1(2*ifft  ,2)
           vtrial1_mq(2*ifft-1,3)= vtrial1(2*ifft  ,4) !Re[V^12]
           vtrial1_mq(2*ifft  ,3)= vtrial1(2*ifft-1,4) !Im[V^12],see definition of v(:,4) cplex=2 case
           vtrial1_mq(2*ifft  ,4)= vtrial1(2*ifft-1,3) !Re[V^21]=Re[V^12]
           vtrial1_mq(2*ifft-1,4)= vtrial1(2*ifft  ,3) !Re[V^21]=Re[V^12]
         end do
       end if
       initialized=1
     end if
   end if

!  ######################################################################
!  END MINIMIZATION ITERATIONS
!  Note that there are different "exit" instructions within the loop
!  ######################################################################
 end do ! istep
print *, 'istep loop finished'

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
&   dbl_nnsclo,dim_eig2rf,doccde_rbz,docckqde,dtefield,dtfil,dtset,dtset%qptn,edocc,&
&   eeig0,eigenq,eigen0,eigen1,ek0,ek1,eloc0,enl0,enl1,fermie1,gh0c1_set,gh1c_set,&
&   gmet,gprimd,idir,indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,dtset%mband,dtset%mband_mem,&
&   mkmem,mkqmem,mk1mem,mpi_enreg,mpw,mpw1,my_natom,dtset%natom,nband_rbz,ncpgr,nfftf,&
&   nhat1,nkpt_rbz,npwarr,npwar1,res3,nspden,dtset%nsppol,nsym1,dtset%ntypat,nvresid2,&
&   occkq,occ_rbz,optres,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrhoij,&
&   pawrhoij1,pawtab,phnons1,ph1d,dtset%prtvol,psps,pwindall,qmat,resid2,residm2,rhog2,&
&   rhor2,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,ucvol,usecprj,useylmgr1,ddk_f,&
&   vtrial,vtrial2,wtk_rbz,xred,ylm,ylm1,ylmgr1,1)

   write(msg,'(a)') ' '//char(10)//&
'   ---------------------------------'
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,a)')'  The charge sum rule is activated'//char(10)//&
'   ---------------------------------'
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,i4)') ' Z_ion (psp):',Z_kappa
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f12.8)') ' Residual Born effective charge: ',born
   call wrtout(ab_out,msg,'COLL')
   write(msg,'(a,f12.8)') ' Renormalisation: ',renorm
   call wrtout(ab_out,msg,'COLL')
   if (renorm > 0.01 ) then
     write(msg,'(a,a)')'   WARNING: The renormalisation seems large (> 0.01).'//char(10)//&
'     You might consider increasing the k-point grid.'
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
   if (.not.kramers_deg) then
     ABI_DEALLOCATE(nvresid1_mq)
   end if
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
   call ddk_f(1)%close()
 end if
 if ((ipert==dtset%natom+10 .and. idir>3) .or. ipert==dtset%natom+11) then
   call ddk_f(2)%close()
 end if
 if (ipert==dtset%natom+11) then
   call ddk_f(3)%close()
   if(idir>3) call ddk_f(4)%close()
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
 if(allocated(rhorfermi_mq)) then
   ABI_DEALLOCATE(rhorfermi_mq)
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
print *, 'calling nselt ', ipert
   call dfpt_nselt(blkflg,cg,cg1,cplex,&
&   d2bbb,d2lo,d2nl,ecut,dtset%ecutsm,dtset%effmass_free,&
&   gmet,gprimd,gsqcut,idir,&
&   ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,dtset%mband,dtset%mband_mem,mgfftf,&
&   mkmem,mk1mem,mpert,mpi_enreg,psps%mpsang,mpw,mpw1,&
&   dtset%natom,nband_rbz,nfftf,ngfftf,&
&   nkpt_rbz,nkxc,dtset%nloalg,&
&   npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,&
&   nsym1,dtset%ntypat,occ_rbz,&
&   ph1d,dtset%prtbbb,psps,dtset%qptn,rhog,&
&   rhor,rhor1,rmet,rprimd,symrc1,dtset%typat,ucvol,&
&   wtk_rbz,xred,ylm,ylm1,ylmgr,ylmgr1)
 end if

!Use of NSTPAW3 for NCPP (instead of DFPT_NSELT/DFPT_NSTDY) can be forced with userie=919
!MT oct. 2015: this works perfectly on all automatic tests
 if(ipert<=dtset%natom+4)then
   if (psps%usepaw==1.or.dtset%userie==919) then
!TODO MJV : PAW
     call dfpt_nstpaw(blkflg,cg,cgq,cg1,cplex,cprj,cprjq,docckqde,doccde_rbz,dtfil,dtset,d2lo,d2nl,d2ovl,&
&     eigenq,eigen0,eigen1,eovl1,gmet,gprimd,gsqcut,idir,indkpt1,indsy1,ipert,irrzon1,istwfk_rbz,&
&     kg,kg1,kpt_rbz,kxc,mgfftf,mkmem,mkqmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,ncpgr,&
&     nfftf,ngfftf,nhat,nhat1,nkpt_rbz,nkxc,npwarr,npwar1,nspden,dtset%nspinor,dtset%nsppol,&
&     nsym1,n3xccc,occkq,occ_rbz,paw_an,paw_an1,paw_ij,paw_ij1,pawang,pawang1,pawfgr,pawfgrtab,pawrad,&
&     pawrhoij,pawrhoij1,pawtab,phnons1,ph1d,ph1df,psps,rhog,rhor,rhor1,rmet,rprimd,symaf1,symrc1,&
&     symrl1,tnons1,ucvol,usecprj,psps%usepaw,usexcnhat,useylmgr1,vhartr1,vpsp1,vtrial,vtrial1,vxc,wtk_rbz,&
&     xccc3d1,xred,ylm,ylm1,ylmgr1)
   else
     if (dtset%nspden==4) then
       call dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,gmet,&
&       gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,mpert,mpi_enreg,&
&       mpw,mpw1,nattyp,nband_rbz,nfftf,ngfftf,nkpt,nkpt_rbz,nkxc,npwarr,npwar1,nspden,&
&       dtset%nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,symrc1,ucvol,&
&       wtk_rbz,xred,ylm,ylm1,rhor=rhor,vxc=vxc)
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
   call dfptff_bec(cg,cg1,dtefield,dtset%natom,d2lo,idir,ipert,dtset%mband,dtset%mband_mem,mkmem,&
&   mpi_enreg,mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,1:dtset%natom)=1
 end if

!calculate dielectric tensor and store it in d2lo
 if ((dtset%berryopt== 4.or.dtset%berryopt== 6.or.dtset%berryopt== 7.or.&
& dtset%berryopt==14.or.dtset%berryopt==16.or.dtset%berryopt==17).and.&
& ipert==dtset%natom+2) then
   call dfptff_die(cg,cg1,dtefield,d2lo,idir,ipert,dtset%mband,dtset%mband_mem,mkmem,&
&   mpi_enreg,mpw,mpw1,mpert,nkpt,npwarr,npwar1,dtset%nsppol,dtset%nspinor,pwindall,qmat,rprimd)
   blkflg(:,dtset%natom+2,:,dtset%natom+2)=1
 end if

!If SCF convergence was not reached (for nstep>0),
!print a warning to the output file (non-dummy arguments: nstep,
!residm, diffor - infos from tollist have been saved inside )
!Set also the value of conv_retcode
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,eigen0,&
& etotal,favg,fcart,fermie,dtfil%fnametmp_eig,dtfil%filnam_ds(1),&
& 1,iscf_mod,istep,istep_fock_outer,istep_mix,kpt_rbz,maxfor,&
& mvdum,mpi_enreg,nband_rbz,nkpt_rbz,&
& nstep,occ_rbz,0,prtfor,0,&
& quit,res2,resid,residm,response,&
& tollist,psps%usepaw,vxcavg,wtk_rbz,xred,conv_retcode)

!Update the content of the header (evolving variables)
 bantot_rbz = sum(nband_rbz(1:nkpt_rbz*dtset%nsppol))
 call hdr%update(bantot_rbz,etotal,fermie,&
& residm,rprimd,occ_rbz,pawrhoij1,xred,dtset%amu_orig(:,1),&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab )

!Optionally provide output of charge density and/or potential in real space,
!as well as analysis of geometrical factors (bond lengths and bond angles).
!Warnings :
!- core charge is excluded from the charge density;
!- the potential is the INPUT vtrial.

 if(ipert==dtset%natom+5.or.ipert<=dtset%natom)then
   prtopt=1
   if(ipert==dtset%natom+5) then
     prtopt=idir+1;
     call calcdenmagsph(gmet,mpi_enreg,dtset%natom,nfftf,ngfftf,nspden,&
&     dtset%ntypat,dtset%ratsm,dtset%ratsph,rhor1,rprimd,dtset%typat,ucvol,xred,&
&     prtopt,cplex,intgden=intgden,dentot=dentot,rhomag=rhomag)
     call  prtdenmagsph(cplex,intgden,dtset%natom,nspden,dtset%ntypat,ab_out,prtopt,dtset%ratsm,dtset%ratsph,rhomag,dtset%typat)
     !debug: write out the vtk first-order density components
!    call appdig(pertcase,dtfil%fnameabo_den,fi1o_vtk)
!    call printmagvtk(mpi_enreg,cplex,nspden,nfftf,ngfftf,rhor1,rprimd,adjustl(adjustr(fi1o_vtk)//"_PQ"))
!    call printmagvtk(mpi_enreg,cplex,nspden,nfftf,ngfftf,rhor1,rprimd,adjustl(adjustr(fi1o_vtk)//"_MQ"))
     !SPr: add calculation of the contributions to susceptibility from all atomic spheres
   end if
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
 if (.not.kramers_deg) then
   ABI_DEALLOCATE(vhartr1_mq)
   ABI_DEALLOCATE(vxc1_mq)
   ABI_DEALLOCATE(vtrial1_pq)
   ABI_DEALLOCATE(vtrial1_mq)
 end if
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

 call timab(160,2,tsec)
 call timab(120,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_scfcv
!!***

!!****f* ABINIT/dfpt_etot
!! NAME
!! dfpt_etot
!!
!! FUNCTION
!! Assemble different contributions to the variational part of the
!! 2nd derivative of total energy
!!
!! INPUTS
!!  berryopt= 4/14: electric field is on; berryopt = 6/7/16/17: electric displacement field is on;
!!  eberry=energy associated with Berry phase
!!  edocc=correction to 2nd-order total energy coming from changes of occupation
!!  ehart1=1st-order Hartree part of 2nd-order total energy
!!  eeig0=0th-order eigenenergies part of 2nd-order total energy
!!  eew=2nd derivative of Ewald energy (hartree)
!!  efrhar=contrib. from frozen-wavefunction, hartree energy, to the 2nd-derivative of total energy
!!  efrkin=contrib. from frozen-wavefunction, kinetic energy, to the 2nd-derivative of total energy
!!  efrloc=contrib. from frozen-wavefunction, local potential, to the 2nd-derivative of total energy
!!  efrnl=contribution from frozen-wavefunction, non-local potential, to the 2nd-derivative of total energy
!!  efrx1=contrib. from frozen-wavefunction, xc core correction(1), to the 2nd-derivative of total energy
!!  efrx2=contribution from frozen-wavefunction, xc core correction(2),
!!           to the second-derivative of total energy.
!!  ek0=0th-order kinetic energy part of 2nd-order total energy.
!!  ek1=1st-order kinetic energy part of 2nd-order total energy.
!!  eii=2nd derivative of pseudopotential core energy (hartree)
!!  eloc0=0th-order local (psp+vxc+Hart) part of 2nd-order total energy
!!  elpsp1=1st-order local pseudopot. part of 2nd-order total energy.
!!  enl0=0th-order nonlocal pseudopot. part of 2nd-order total energy.
!!  enl1=1st-order nonlocal pseudopot. part of 2nd-order total energy.
!!  epaw1=1st-order PAW on-sitew part of 2nd-order total energy.
!!  evdw=DFT-D semi-empirical part of 2nd-order total energy
!!  exc1=1st-order exchange-correlation part of 2nd-order total energy
!!  ipert=type of the perturbation
!!  natom=number of atoms
!!  optene=option for the computation of 2nd-order total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!
!! OUTPUT
!!  deltae=change in energy between the previous and present SCF cycle
!!         and previous SCF cycle.
!!  etotal=2nd-order total energy
!!  evar=variational part of the 2nd-order total energy
!!
!! SIDE EFFECTS
!! input/output
!! elast=previous value of the 2nd-order total energy, needed to compute deltae,
!!      then updated (cannot simply be saved, because set to zero
!!      at each new call of dfpt_scfcv).
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!
!! SOURCE

subroutine dfpt_etot(berryopt,deltae,eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,&
&                efrnl,efrx1,efrx2,ehart1,ek0,ek1,eii,elast,eloc0,elpsp1,&
&                enl0,enl1,epaw1,etotal,evar,evdw,exc1,ipert,natom,optene)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,ipert,natom,optene
 real(dp),intent(in) :: eberry,edocc,eeig0,eew,efrhar,efrkin,efrloc,efrnl,efrx1
 real(dp),intent(in) :: efrx2,ehart1,eii,ek0,ek1,eloc0,elpsp1,enl0,enl1,epaw1
 real(dp),intent(in) :: evdw,exc1
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,etotal,evar

!Local variables-------------------------------
!scalars
! character(len=500) :: message

! *********************************************************************

 if (optene==1) then
   MSG_BUG('Double-counting scheme not yet allowed!')
 end if

 if (optene>-1) then

!  Compute 2nd-order variational energy by direct scheme
   if (optene==0) then

!    Atomic displ. perturbation
     if ( ipert>=1 .and. ipert<=natom  ) then
       evar=ek0+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1

     else if (ipert==natom+1) then
       evar=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1

     else if (ipert==natom+10 .or. ipert==natom+11) then
       evar=ek0+edocc+eeig0+eloc0+enl0+ek1 ! here ek1 contains a lot of contributions

!      For ipert==natom+2, some contributions vanishes, noticeably ek1
     else if (ipert==natom+2) then
       evar=ek0+edocc+eeig0+eloc0+ek1+ehart1+exc1+enl0+enl1+epaw1

!      All terms enter for strain perturbation
     else if ( ipert==natom+3 .or. ipert==natom+4 ) then
       evar=ek0+ek1+edocc+eeig0+eloc0+elpsp1+ehart1+exc1+enl0+enl1+epaw1

!    terms for Zeeman perturbation, SPr 2deb
     else if ( ipert==natom+5 ) then
       evar=ek0+edocc+eeig0+eloc0+ehart1+exc1+enl0+epaw1
     end if
   end if

!  Compute energy residual
   deltae=evar-elast
   elast=evar

!  Compute 2nd-order total energy by direct scheme
   if (optene==0) then
     if (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. berryopt==14 .or. berryopt==16 .or. berryopt==17) then
       if (ipert<=natom) then
         etotal=evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2+two*eberry
       else if (ipert==natom+2) then
         etotal=half*evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2+two*eberry
       end if
     else
       if (ipert/=natom+10 .and. ipert/=natom+11) then
         etotal=evar+eew+evdw+eii+efrhar+efrkin+efrloc+efrnl+efrx1+efrx2
       else
         etotal=evar ! For 2nd order sternheimer equations, the total (4th order) energy is not used (yet)
       end if
     end if
   end if

 end if

end subroutine dfpt_etot
!!***

!!****f* ABINIT/newfermie1
!! NAME
!! newfermie1
!!
!! FUNCTION
!! This routine computes the derivative of the fermi energy wrt
!! the active perturbation for use in evaluating the edocc term
!! and active subspace contribution to the first-order wavefunctions
!! in the case of metals. This is presently used only for the
!! strain and magnetic field perturbations, and only for Q = 0.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!  ipert=index of perturbation
!!  istep=index of the number of steps in the routine scfcv
!!  ixc= choice of exchange-correlation scheme
!!  mpi_atmtab(:)=--optional-- indexes of the atoms treated by current proc
!!  comm_atom=--optional-- MPI communicator over atoms
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftot= total number of FFT grid points
!!  nhatfermi(nfft,nspden)=fermi-level compensation charge density (PAW only)
!!  nspden=number of spin-density components
!!  ntypat=number of atom types
!!  occopt=option for occupancies
!!  paw_an(natom) <type(paw_an_type)>=paw arrays for 0th-order quantities given on angular mesh
!!  paw_an1(natom) <type(paw_an_type)>=paw arrays for 1st-order quantities given on angular mesh
!!  paw_ij1(natom) <type(paw_ij_type)>=(1st-order) paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawnzlm=-- PAW only -- option for the computation of non-zero
!!          lm moments of the on-sites densities
!!  pawrad(ntypat) <type(pawrad_type)>=paw radial mesh and related data
!!  pawrhoij1(natom) <type(pawrhoij_type)>= paw rhoij 1st-order occupancies
!!  pawrhoijfermi(natom) <type(pawrhoij_type)>=paw rhoij occupancies at Fermi level
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pawxcdev=Choice of XC development (0=no dev. (use of angular mesh) ; 1 or 2=dev. on moments)
!!  xclevel= XC functional level
!!  prtvol=control print volume and debugging output
!!  rhorfermi(nfft,nspden)=fermi-level electronic density
!!  ucvol=unit cell volume in bohr**3
!!  usepaw=1 if PAW is activated
!!  usexcnhat= -PAW only- flag controling use of compensation density in Vxc
!!  vtrial1(cplex*nfft,nspden)=1-st order potential
!!  vxc1(cplex*nfft,nspden)=1-st order XC potential
!!
!! OUTPUT
!!  (see side effects)
!!
!! SIDE EFFECTS
!!  fermie1=derivative of fermi energy wrt perturbation
!!   at input  : old value
!!   at output : updated value
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dotprod_vn,free_my_atmtab,get_my_atmtab,pawdfptenergy,wrtout
!!
!! SOURCE

subroutine newfermie1(cplex,fermie1,fe1fixed,ipert,istep,ixc,my_natom,natom,nfft,nfftot,&
&                     nhatfermi,nspden,ntypat,occopt,paw_an,paw_an1,paw_ij1,pawang,pawnzlm,pawrad,&
&                     pawrhoij1,pawrhoijfermi,pawtab,pawxcdev,prtvol,rhorfermi,&
&                     ucvol,usepaw,usexcnhat,vtrial1,vxc1,xclevel,&
&                     mpi_atmtab,comm_atom) ! optional arguments (parallelism)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,ipert,istep,ixc,my_natom,natom,nfft,nfftot,nspden,ntypat
 integer,intent(in) :: occopt,pawnzlm,pawxcdev,prtvol,usepaw,usexcnhat,xclevel
 integer,optional,intent(in) :: comm_atom
 real(dp),intent(in) :: fe1fixed,ucvol
 real(dp),intent(inout) :: fermie1
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,optional,target,intent(in) :: mpi_atmtab(:)
 real(dp),intent(in) :: rhorfermi(nfft,nspden),vtrial1(cplex*nfft,nspden)
 real(dp),intent(in) :: nhatfermi(:,:),vxc1(:,:)
 type(paw_an_type),intent(in) :: paw_an(my_natom*usepaw)
 type(paw_an_type),intent(inout) :: paw_an1(my_natom*usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij1(my_natom*usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij1(my_natom*usepaw),pawrhoijfermi(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: ipert0,my_comm_atom,nzlmopt,nzlmopt_fermi,option,pawprtvol
 logical :: my_atmtab_allocated,paral_atom
 real(dp) :: doti,fe1_scf,fe1_tmp,fermie1_new,fermie1rs
 character(len=500) :: msg
!arrays
 integer, pointer :: my_atmtab(:)
 real(dp) :: fe1_paw(2)
 real(dp), allocatable :: rhor_nonhat(:,:),vtrial1_novxc(:,:)

! *********************************************************************

!Tests
 if (cplex==2) then
   msg='Not compatible with cplex=2!'
   MSG_BUG(msg)
 end if
 if (usepaw==1.and.usexcnhat==0.and.(size(nhatfermi)<=0.or.size(vxc1)<=0)) then
   msg='Should have nhatfermi and vxc1 allocated with usexcnhat=0!'
   MSG_BUG(msg)
 end if

!Set up parallelism over atoms
 paral_atom=(present(comm_atom).and.(my_natom/=natom))
 nullify(my_atmtab);if (present(mpi_atmtab)) my_atmtab => mpi_atmtab
 my_comm_atom=xmpi_comm_self;if (present(comm_atom)) my_comm_atom=comm_atom
 call get_my_atmtab(my_comm_atom,my_atmtab,my_atmtab_allocated,paral_atom,natom,my_natom_ref=my_natom)

 if(occopt>=3 .and. occopt <=8) then

!  The product of the current trial potential and the so-called Fermi level
!  density is integrated to give the local potential contributions to the
!  first-order Fermi level.
   option=1
   if (usepaw==1.and.usexcnhat==0) then
     ABI_ALLOCATE(rhor_nonhat,(nfft,nspden))
     ABI_ALLOCATE(vtrial1_novxc,(nfft,nspden))
     rhor_nonhat(1:nfft,1:nspden)=rhorfermi(1:nfft,1:nspden)-nhatfermi(1:nfft,1:nspden)
     vtrial1_novxc(1:nfft,1:nspden)=vtrial1(1:nfft,1:nspden)-vxc1(1:nfft,1:nspden)
     call dotprod_vn(cplex,rhor_nonhat,fe1_scf,doti,nfft,nfftot,&
&     nspden,option,vtrial1,ucvol)
     call dotprod_vn(cplex,nhatfermi,fe1_tmp,doti,nfft,nfftot,&
&     nspden,option,vtrial1_novxc,ucvol)
     fe1_scf=fe1_scf+fe1_tmp
     ABI_DEALLOCATE(rhor_nonhat)
     ABI_DEALLOCATE(vtrial1_novxc)
   else
     call dotprod_vn(cplex,rhorfermi,fe1_scf,doti,nfft,nfftot,&
&     nspden,option,vtrial1,ucvol)
   end if

   fe1_paw(:)=zero
!  PAW on-site contribution (use Fermi level occupation matrix)
   if (usepaw==1) then
     ipert0=0;pawprtvol=0
     nzlmopt=0;if (istep>1) nzlmopt=pawnzlm
     if (istep==1.and.pawnzlm>0) nzlmopt=-1
     nzlmopt_fermi=0;if (pawnzlm>0) nzlmopt_fermi=-1
     call pawdfptenergy(fe1_paw,ipert,ipert0,ixc,my_natom,natom,ntypat,nzlmopt,&
&     nzlmopt_fermi,paw_an,paw_an1,paw_ij1,pawang,pawprtvol,pawrad,&
&     pawrhoij1,pawrhoijfermi,pawtab,pawxcdev,xclevel,&
&     mpi_atmtab=my_atmtab, comm_atom=my_comm_atom)
   end if

!  The fixed contributions consisting of non-local potential and kinetic terms
!  are added
   fermie1_new=fe1fixed+fe1_scf+fe1_paw(1)
   fermie1rs=(fermie1-fermie1_new)**2
   fermie1=fermie1_new

   if(prtvol>=10)then
     write(msg, '(a,i5,2es18.8)' ) ' fermie1, residual squared',istep,fermie1,fermie1rs
     call wrtout(std_out,msg,'COLL')
   end if

 else
   fermie1=zero
 end if

!Destroy atom table used for parallelism
 call free_my_atmtab(my_atmtab,my_atmtab_allocated)

end subroutine newfermie1
!!***

!!****f* ABINIT/dfpt_newvtr
!! NAME
!! dfpt_newvtr
!!
!! FUNCTION
!! Compute new first-order trial potential by mixing new and old values.
!! First, compute preconditioned residual first-order potential.
!! Then, call one of the self-consistency drivers, and  update vtrial.
!!
!! INPUTS
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  dielar(7)=input parameters for dielectric matrix:
!!                diecut,dielng,diemac,diemix,diegap,dielam,diemixmag.
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | isecur=level of security of the computation
!!   | mffmem=governs the number of FFT arrays which are fit in core memory
!!   |          it is either 1, in which case the array f_fftgr is used,
!!   |          or 0, in which case the array f_fftgr_disk is used
!!   | natom=number of atoms
!!   | nspden=number of spin-density components
!!   | paral_kgb=option for (kpt,g vectors,bands) parallelism
!!   | pawoptmix= - PAW only - 1 if the computed residuals include the PAW (rhoij) part
!!  etotal=the total energy obtained from the input vtrial
!!  ffttomix(nfft*(1-nfftmix/nfft))=Index of the points of the FFT (fine) grid on the grid used for mixing (coarse)
!!  initialized= if 0 the initialization of the RF run is not yet finished
!!   iscf=( <= 0 =>non-SCF), >0 => SCF)
!!    iscf =1 => determination of the largest eigenvalue of the SCF cycle
!!    iscf =2 => SCF cycle, simple mixing
!!    iscf =3 => SCF cycle, Anderson mixing
!!    iscf =4 => SCF cycle, Anderson mixing (order 2)
!!    iscf =5 => SCF cycle, CG based on the minimization of the energy
!!    iscf =7 => SCF cycle, Pulay mixing
!!  ispmix=1 if mixing is done in real space, 2 if mixing is done in reciprocal space
!!  istep= number of the step in the SCF cycle
!!  mixtofft(nfftmix*(1-nfftmix/nfft))=Index of the points of the FFT grid used for mixing (coarse) on the FFT (fine) grid
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  nfftmix=dimension of FFT grid used to mix the densities (used in PAW only)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngfftmix(18)=contain all needed information about 3D FFT, for the grid corresponding to nfftmix
!!  npawmix=-PAW only- number of spherical part elements to be mixed
!!  qphon(3)=reduced coordinates for the phonon wavelength (needed if cplex==2).
!!  pawrhoij(natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!                                         Use here rhoij residuals (and gradients)
!!  rhor(cplex*nfft,nspden)=array for 1st-order electron density
!!    in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vresid(cplex*nfft,nspden)=array for the residual of the potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  dbl_nnsclo=1 if nnsclo has to be doubled to secure the convergence.
!!
!! SIDE EFFECTS
!!  vtrial(cplex*nfft,nspden)= at input, it is the "in" trial potential that gave vresid=(v_out-v_in)
!!       at output, it is an updated "mixed" trial potential
!!  ==== if usepaw==1
!!    pawrhoij(natom)%nrhoijsel,rhoijselect,rhoijp= several arrays
!!                containing new values of rhoij (augmentation occupancies)
!!
!! NOTES
!!  In case of PAW calculations:
!!    Computations are done either on the fine FFT grid or the coarse grid (depending on dtset%pawmixdg)
!!    All variables (nfft,ngfft,mgfft) refer to the fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!    Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!!  Subtility in PAW and non-collinear magnetism:
!!    Potentials are stored in (up-up,dn-dn,Re[up-dn],Im[up-dn]) format
!!    On-site occupancies (rhoij) are stored in (n,mx,my,mz)
!!    This is compatible provided that the mixing factors for n and m are identical
!!    and that the residual is not a combination of V_res and rhoij_res (pawoptmix=0).
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      ab7_mixing_copy_current_step,ab7_mixing_eval,ab7_mixing_eval_allocate
!!      ab7_mixing_eval_deallocate,fourdp,metric,moddiel,timab
!!
!! SOURCE

subroutine dfpt_newvtr(cplex,dbl_nnsclo,dielar,dtset,etotal,ffttomix,&
&          initialized,iscf,ispmix,istep,mix,mixtofft,&
&          mpi_enreg,my_natom,nfft,nfftmix,ngfft,ngfftmix,npawmix,pawrhoij,&
&          qphon,rhor,rprimd,usepaw,vresid,vtrial)

!Arguments-------------------------------
!scalars
 integer,intent(in) :: cplex,initialized,iscf,ispmix,istep,my_natom,nfft
 integer,intent(in) :: nfftmix,npawmix,usepaw
 integer,intent(inout) :: dbl_nnsclo !vz_i
 real(dp),intent(in) :: etotal
 type(MPI_type),intent(in) :: mpi_enreg
 type(ab7_mixing_object), intent(inout) :: mix
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ffttomix(nfft*(1-nfftmix/nfft))
 integer,intent(in) :: mixtofft(nfftmix*(1-nfftmix/nfft)),ngfft(18)
 integer,intent(in) :: ngfftmix(18)
 real(dp),intent(in) :: dielar(7),qphon(3)
 real(dp), intent(in), target :: rhor(cplex*nfft,dtset%nspden)
 real(dp),intent(in) :: rprimd(3,3)
 real(dp),intent(inout) :: vresid(cplex*nfft,dtset%nspden)
 real(dp),intent(inout) :: vtrial(cplex*nfft,dtset%nspden)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*usepaw)

!Local variables-------------------------------
!scalars
 integer :: cplex_mix,cplex_rhoij,dplex,i_vresid1,i_vrespc1,iatom,ifft,indx,iq,iq0
 integer :: irhoij,ispden,jfft,jrhoij,klmn,kklmn,kmix,moved_atm_inside,nfftot,qphase
 integer :: mpicomm,errid
 logical :: mpi_summarize,reset
 real(dp) :: fact,mixfac,mixfac_eff,mixfacmag,ucvol
 character(len=500) :: msg
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 real(dp),allocatable :: rhoijrespc(:),rhoijtmp(:,:)
 real(dp),allocatable :: vresid0(:,:),vrespc(:,:),vreswk(:,:)
 real(dp), pointer :: vtrial0(:,:),vpaw(:)
 real(dp),allocatable :: vtrialg(:,:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(158,1,tsec)

!Compatibility tests
 if(usepaw==1) then
   if(dtset%nspden==4.and.dtset%pawoptmix==1) then
     MSG_ERROR('pawoptmix=1 is not compatible with nspden=4 !')
   end if
   if (my_natom>0) then
     if (pawrhoij(1)%qphase<cplex) then
       MSG_ERROR('pawrhoij()%qphase must be >=cplex !')
     end if
   end if
 end if

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 cplex_mix=max(cplex,ispmix)
 if (usepaw==1.and.my_natom>0) then
   cplex_rhoij=pawrhoij(1)%cplex_rhoij
   qphase=pawrhoij(1)%qphase
 end if

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)
 moved_atm_inside=0

!Select components of potential to be mixed
 ABI_ALLOCATE(vtrial0,(cplex_mix*nfftmix,dtset%nspden))
 ABI_ALLOCATE(vresid0,(cplex_mix*nfftmix,dtset%nspden))
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial0=vtrial;vresid0=vresid
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(cplex,vtrial0(:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,1, ngfft, 0)
     call fourdp(cplex,vresid0(:,ispden),vresid(:,ispden),-1,mpi_enreg,nfft,1, ngfft, 0)
   end do
 else
   ABI_ALLOCATE(vtrialg,(2,nfft,dtset%nspden))
   ABI_ALLOCATE(vreswk,(2,nfft))
   do ispden=1,dtset%nspden
     fact=dielar(4);if (ispden>1) fact=dielar(7)
     call fourdp(cplex,vtrialg(:,:,ispden),vtrial(:,ispden),-1,mpi_enreg,nfft,1, ngfft, 0)
     call fourdp(cplex,vreswk,vresid(:,ispden),-1,mpi_enreg,nfft,1, ngfft, 0)
     do ifft=1,nfft
       if (ffttomix(ifft)>0) then
         jfft=2*ffttomix(ifft)
         vtrial0(jfft-1,ispden)=vtrialg(1,ifft,ispden)
         vtrial0(jfft  ,ispden)=vtrialg(2,ifft,ispden)
         vresid0(jfft-1,ispden)=vreswk(1,ifft)
         vresid0(jfft  ,ispden)=vreswk(2,ifft)
       else
         vtrialg(:,ifft,ispden)=vtrialg(:,ifft,ispden)+fact*vreswk(:,ifft)
       end if
     end do
   end do
   ABI_DEALLOCATE(vreswk)
 end if

!Precondition the potential residual:
!Use a model dielectric function preconditioning, or simple mixing
 ABI_ALLOCATE(vrespc,(cplex_mix*nfftmix,dtset%nspden))
 call moddiel(cplex_mix,dielar,mpi_enreg,nfftmix,ngfftmix,dtset%nspden,ispmix,0,qphon,rprimd,vresid0,vrespc)

!PAW only : precondition the rhoij quantities (augmentation occupancies) residuals.
!Use a simple preconditionning with the same mixing factor
!as the model dielectric function.
 if (usepaw==1.and.my_natom>0) then
   ABI_ALLOCATE(rhoijrespc,(npawmix))
   mixfac=dielar(4);mixfacmag=abs(dielar(7))
   if (cplex_rhoij==1) then
     indx=0
     do iatom=1,my_natom
       do iq=1,qphase
         iq0=merge(0,cplex_rhoij*pawrhoij(iatom)%lmn2_size,iq==1)
         do ispden=1,pawrhoij(iatom)%nspden
           mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             indx=indx+1;klmn=iq0+pawrhoij(iatom)%kpawmix(kmix)
             rhoijrespc(indx)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn,ispden)
           end do
         end do
       end do
     end do
   else
     indx=-1
     do iatom=1,my_natom
       do iq=1,qphase
         iq0=merge(0,cplex_rhoij*pawrhoij(iatom)%lmn2_size,iq==1)
         do ispden=1,pawrhoij(iatom)%nspden
           mixfac_eff=mixfac;if (ispden>1) mixfac_eff=mixfacmag
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             indx=indx+2;klmn=iq0+2*pawrhoij(iatom)%kpawmix(kmix)-1
             rhoijrespc(indx:indx+1)=mixfac_eff*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
           end do
         end do
       end do
     end do
   end if
 end if

!------Compute new vtrial

 i_vresid1=mix%i_vresid(1)
 i_vrespc1=mix%i_vrespc(1)

!Initialise working arrays for the mixing object.
 call ab7_mixing_eval_allocate(mix, istep)

!Copy current step arrays.
 call ab7_mixing_copy_current_step(mix, vresid0, errid, msg, arr_respc = vrespc)

 if (errid /= AB7_NO_ERROR) then
   MSG_ERROR(msg)
 end if

 ABI_DEALLOCATE(vrespc)
 ABI_DEALLOCATE(vresid0)

!PAW: either use the array f_paw or the array f_paw_disk
 ABI_ALLOCATE(vpaw,(npawmix*usepaw))
 if (usepaw==1.and.my_natom>0) then
   dplex=cplex_rhoij-1 ; indx=-dplex
   do iatom=1,my_natom
     ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*pawrhoij(iatom)%lmn2_size,1))
     do iq=1,qphase
       iq0=merge(0,cplex_rhoij*pawrhoij(iatom)%lmn2_size,iq==1)
       do ispden=1,pawrhoij(iatom)%nspden
         rhoijtmp=zero ; jrhoij=iq0+1
         do irhoij=1,pawrhoij(iatom)%nrhoijsel
           klmn=cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
           rhoijtmp(klmn:klmn+dplex,1)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
           jrhoij=jrhoij+cplex_rhoij
         end do
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           indx=indx+cplex_rhoij;klmn=cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex ; kklmn=iq0+klmn
           vpaw(indx:indx+dplex)=rhoijtmp(klmn:klmn+dplex,1)-pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vresid1)=pawrhoij(iatom)%rhoijres(kklmn:kklmn+dplex,ispden)
           mix%f_paw(indx:indx+dplex,i_vrespc1)=rhoijrespc(indx:indx+dplex)
         end do
       end do
     end do
     ABI_DEALLOCATE(rhoijtmp)
   end do
 end if

!Unlike for GS, no need to modify the mean of vtrial

 mpicomm=0;mpi_summarize=.false.
 reset=.false.;if (initialized==0) reset=.true.
 call ab7_mixing_eval(mix, vtrial0, istep, nfftot, ucvol, &
& mpicomm, mpi_summarize, errid, msg, &
& reset = reset, isecur = dtset%isecur, &
& pawopt = dtset%pawoptmix, response = 1, pawarr = vpaw, &
& etotal = etotal, potden = rhor, comm_atom=mpi_enreg%comm_atom)

 if (errid == AB7_ERROR_MIXING_INC_NNSLOOP) then
   dbl_nnsclo = 1
 else if (errid /= AB7_NO_ERROR) then
   ! MG FIXME, Why this?
   ! One should propagate the error so that we can handle it
   ! in the caller!
   MSG_ERROR(msg)
 end if

!Do here the mixing of the potential
 if(iscf==2 .or. iscf==3 .or. iscf==7)then
!  PAW: restore rhoij from compact storage
   if (usepaw==1.and.my_natom>0) then
     dplex=cplex_rhoij-1 ; indx=-dplex
     do iatom=1,my_natom
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*qphase*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       do iq=1,qphase
         iq0=merge(0,cplex_rhoij*pawrhoij(iatom)%lmn2_size,iq==1)
         if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
           do ispden=1,pawrhoij(iatom)%nspden
             jrhoij=iq0+1
             do irhoij=1,pawrhoij(iatom)%nrhoijsel
               klmn=iq0+cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
               rhoijtmp(klmn:klmn+dplex,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
               jrhoij=jrhoij+cplex_rhoij
             end do
           end do
         end if
         do ispden=1,pawrhoij(iatom)%nspden
           do kmix=1,pawrhoij(iatom)%lmnmix_sz
             indx=indx+cplex_rhoij;klmn=iq0+cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=vpaw(indx:indx+dplex)
           end do
         end do
       end do
       call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&           pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&           pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if

 else if(iscf==5 .or. iscf==6)then
   if(ispmix/=1) then
     MSG_ERROR('Mixing on reciprocal space not allowed with iscf=5 or 6.')
   end if
!  PAW: apply a simple mixing to rhoij (this is temporary)
   if (usepaw==1.and.my_natom>0) then
     indx=1-cplex_rhoij
     do iatom=1,my_natom
       ABI_ALLOCATE(rhoijtmp,(cplex_rhoij*qphase*pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden))
       rhoijtmp=zero
       do iq=1,qphase
         iq0=merge(0,cplex_rhoij*pawrhoij(iatom)%lmn2_size,iq==1)
         if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
           do ispden=1,pawrhoij(iatom)%nspden
             do kmix=1,pawrhoij(iatom)%lmnmix_sz
               indx=indx+cplex_rhoij;klmn=iq0+cplex_rhoij*pawrhoij(iatom)%kpawmix(kmix)-dplex
               rhoijtmp(klmn:klmn+dplex,ispden)=rhoijrespc(indx:indx+dplex) &
&               -pawrhoij(iatom)%rhoijres(klmn:klmn+dplex,ispden)
             end do
           end do
         end if
         do ispden=1,pawrhoij(iatom)%nspden
           jrhoij=iq0+1
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=iq0+cplex_rhoij*pawrhoij(iatom)%rhoijselect(irhoij)-dplex
             rhoijtmp(klmn:klmn+dplex,ispden)=rhoijtmp(klmn:klmn+dplex,ispden) &
&             +pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+dplex,ispden)
             jrhoij=jrhoij+cplex_rhoij
           end do
         end do
       end do
       call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,&
&           pawrhoij(iatom)%nrhoijsel,pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,&
&           pawrhoij(iatom)%lmn2_size,pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
       ABI_DEALLOCATE(rhoijtmp)
     end do
   end if
 end if

 ABI_DEALLOCATE(vpaw)
 if (usepaw==1.and.my_natom>0)  then
   ABI_DEALLOCATE(rhoijrespc)
 end if

!Eventually write the data on disk and deallocate f_fftgr_disk
 call ab7_mixing_eval_deallocate(mix)

!Restore potential
 if (ispmix==1.and.nfft==nfftmix) then
   vtrial=vtrial0
 else if (nfft==nfftmix) then
   do ispden=1,dtset%nspden
     call fourdp(cplex,vtrial0(:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,1, ngfft,0)
   end do
 else
   do ispden=1,dtset%nspden
     do ifft=1,nfftmix
       jfft=mixtofft(ifft)
       vtrialg(1,jfft,ispden)=vtrial0(2*ifft-1,ispden)
       vtrialg(2,jfft,ispden)=vtrial0(2*ifft  ,ispden)
     end do
     call fourdp(cplex,vtrialg(:,:,ispden),vtrial(:,ispden),+1,mpi_enreg,nfft,1,ngfft,0)
   end do
   ABI_DEALLOCATE(vtrialg)
 end if
 ABI_DEALLOCATE(vtrial0)

 call timab(158,2,tsec)

 DBG_ENTER("COLL")

end subroutine dfpt_newvtr
!!***

!!****f* ABINIT/dfpt_nselt
!! NAME
!! dfpt_nselt
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, wrt strain for a whole row of
!! mixed strain derivatives.
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL,
!!    if 2, COMPLEX
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  ipert=type of the perturbation
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the
!!     storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  mgfft=maximum size of 1D FFTs
!!  mkmem =number of k points treated by this node.
!!  mk1mem =number of k points treated by this node  (RF data).
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT,
!!    see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0,
!!   the exchange-correlation kernel must be computed.
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  ntypat=number of types of atoms in unit cell.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtbbb=if 1, band-by-band decomposition (also dim of d2bbb)
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  qphon(3)=reduced coordinates for the phonon wavelength
!!  rhog(2,nfft)=array for Fourier transform of GS electron density
!!  rhor(nfft,nspden)=GS electron density in electrons/bohr**3.
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  typat(natom)=type integer for each atom in cell
!!  ucvol=unit cell volume in bohr**3.
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang)= real spherical harmonics for each G and k+q point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!  ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k+g point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!       second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      dfpt_mkcore,dfpt_mkvxcstr,dfpt_nsteltwf,dotprod_vn
!!      hartrestr,init_hamiltonian,stresssym,vlocalstr,wrtout,xmpi_barrier
!!
!! SOURCE

subroutine dfpt_nselt(blkflg,cg,cg1,cplex,&
& d2bbb,d2lo,d2nl,ecut,ecutsm,effmass_free,&
& gmet,gprimd,gsqcut,idir,&
& ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mband,mband_mem,mgfft,&
& mkmem,mk1mem,mpert,mpi_enreg,mpsang,mpw,mpw1,&
& natom,nband_rbz,nfft,ngfft,&
& nkpt_rbz,nkxc,nloalg,npwarr,npwar1,nspden,nspinor,nsppol,&
& nsym1,ntypat,occ_rbz,&
& ph1d,prtbbb,psps,qphon,rhog,&
& rhor,rhor1,rmet,rprimd,symrc1,typat,ucvol,&
& wtk_rbz,&
& xred,ylm,ylm1,ylmgr,ylmgr1)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mgfft,mk1mem
 integer,intent(in) :: mband_mem
 integer,intent(in) :: mkmem,mpert,mpsang,mpw,mpw1,natom,nfft,nkpt_rbz
 integer,intent(in) :: nkxc,nspden,nspinor,nsppol,nsym1,ntypat
 integer,intent(in) :: prtbbb
 real(dp),intent(in) :: ecut,ecutsm,effmass_free,gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk_rbz(nkpt_rbz)
 integer,intent(in) :: kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: nloalg(3),npwar1(nkpt_rbz),npwarr(nkpt_rbz)
 integer,intent(in) :: symrc1(3,3,nsym1),typat(natom)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)
 real(dp),intent(in) :: gmet(3,3)
 real(dp),intent(in) :: gprimd(3,3),kpt_rbz(3,nkpt_rbz),kxc(nfft,nkxc)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qphon(3),rhog(2,nfft)
 real(dp),intent(in) :: rhor(nfft,nspden)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: d2bbb(2,3,3,mpert,mband,mband*prtbbb)
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert)
 real(dp),intent(inout) :: d2nl(2,3,mpert,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: ban2tot,bantot,bd2tot_index,bdtot_index
 integer :: icg,icg1,idir1,ifft,ii,ikg,ikg1,ikpt,comm
 integer :: ilm,ipert1,ispden,isppol,istr1,istwf_k
 integer :: mbd2kpsp,mbdkpsp,me,n1,n2,n3,n3xccc,n4,n5,n6
 integer :: nband_k,nfftot,npw1_k,npw_k,option
 logical :: nmxc=.false.
 real(dp) :: doti,dotr
 real(dp) :: wtk_k
 character(len=500) :: message
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 integer :: ikpt_fbz(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),restr(6),dummy(0,0)
 real(dp),allocatable :: d2bbb_k(:,:,:,:),d2nl_k(:,:,:)
 real(dp),allocatable :: occ_k(:)
 real(dp),allocatable :: vhartr01(:),vpsp1(:),vxc1(:,:),xccc3d1(:),ylm1_k(:,:)
 real(dp),allocatable :: ylm_k(:,:),ylmgr1_k(:,:,:),ylmgr_k(:,:,:)
 type(pawtab_type) :: pawtab_dum(0)


! *********************************************************************

!Init me
 comm = mpi_enreg%comm_cell
 me   = mpi_enreg%me_kpt

!Unit numbers

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,natom+3:natom+4,idir,ipert)=zero
 bdtot_index=0
 bd2tot_index=0
 icg=0
 icg1=0
 mbdkpsp=mband*nkpt_rbz*nsppol
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol

!Update list of computed matrix elements
 if((ipert==natom+3) .or. (ipert==natom+4)) then
!  Eventually expand when strain coupling to other perturbations is implemented
   do ipert1=natom+3,natom+4
     do idir1=1,3
       blkflg(idir1,ipert1,idir,ipert)=1
     end do
   end do
 end if

 ABI_ALLOCATE(d2bbb_k,(2,3,mband,mband*prtbbb))
 ABI_ALLOCATE(d2nl_k,(2,3,mpert))

 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))

 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 n4=ngfft(4) ; n5=ngfft(5) ; n6=ngfft(6)
 nfftot=n1*n2*n3

!Initialize Hamiltonian (k-independent terms) - NCPP only
 call init_hamiltonian(gs_hamk,psps,pawtab_dum,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,ph1d=ph1d)

 bantot = 0
 ban2tot = 0

!LOOP OVER SPINS
 do isppol=1,nsppol

   if (nsppol/=1) then
     write(message,*)' ****  In dfpt_nselt for isppol=',isppol
     call wrtout(std_out,message,'COLL')
   end if

   ikg=0
   ikg1=0

   ikpt_fbz(1:3)=0

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)
     kpoint(:)=kpt_rbz(:,ikpt)

     bantot = bantot + nband_k
     ban2tot = ban2tot + 2*nband_k**2

! asserts at least 1 band of the current k and spin is on present processor
print *, ' ik, cycle, nband ', ikpt, proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me),&
&        proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(occ_k,(nband_k))

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang))
     ABI_ALLOCATE(ylm1_k,(npw1_k,mpsang*mpsang))
     if (ipert==natom+3.or.ipert==natom+4) then
       ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang))
       ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,mpsang*mpsang))
     end if

!    enl1_k(:)=zero
     d2nl_k(:,:,:)=zero
     if(prtbbb==1)d2bbb_k(:,:,:,:)=zero
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)

     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
       if (ipert==natom+3.or.ipert==natom+4) then
         do ilm=1,mpsang*mpsang
           do ii=1,3
             ylmgr_k(1:npw_k,ii,ilm)=ylmgr(1+ikg:npw_k+ikg,ii,ilm)
           end do
         end do
       end if
     end if

     wtk_k=wtk_rbz(ikpt)

     kg1_k(:,:) = 0

     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,mpsang*mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (ipert==natom+3.or.ipert==natom+4) then
         do ilm=1,mpsang*mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Compute the eigenvalues, wavefunction,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor1 to this k-point and this spin polarization.

!    Note that dfpt_nsteltwf is called with kpoint, while kpt is used inside dfpt_vtowfk
     call dfpt_nsteltwf(cg,cg1,d2nl_k,ecut,ecutsm,effmass_free,gs_hamk,icg,icg1,ikpt,isppol,&
&     istwf_k,kg_k,kg1_k,kpoint,mband,mband_mem,mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,natom,nband_k,&
&     npw_k,npw1_k,nspinor,nsppol,ntypat,occ_k,psps,rmet,wtk_k,ylm_k,ylmgr_k)
     d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
     if(prtbbb==1)then
       d2bbb(:,:,idir,ipert,:,:) = d2bbb(:,:,idir,ipert,:,:) + &
&       d2bbb_k(:,:,:,:)
     end if

     ABI_DEALLOCATE(occ_k)

!    Keep track of total number of bands (all k points so far, even for
!    k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       icg=icg+npw_k*nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
       ikg1=ikg1+npw1_k
     end if
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     if (ipert==natom+3.or.ipert==natom+4)  then
       ABI_DEALLOCATE(ylmgr_k)
       ABI_DEALLOCATE(ylmgr1_k)
     end if

   end do ! End big k point loop
 end do ! End loop over spins

 if(xmpi_paral==1)then
   call xmpi_barrier(comm)
   call wrtout(std_out,' dfpt_nselt: loop on k-points and spins done in parallel','COLL')
 end if

!Treat now varying occupation numbers
!if(occopt>=3 .and. occopt <=8) then
!SUPPRESSED metallic coding of vtorho

!Treat fixed occupation numbers
!else

!Accumulation over parallel processed now carried out for all terms
!in dfpt_nstdy.f

!End of test on varying or fixed occupation numbers
!end if

!The imaginary part of d2nl will be must be set to zero here since
!time-reversal symmetry will always be true for the strain peturbation.
!The symmetry-reduced kpt set will leave a non-zero imaginary part.

 d2nl(2,:,natom+3:natom+4,idir,ipert)=zero

!Symmetrize the non-local contributions,
!as was needed for the stresses in a ground-state calculation

 if (nsym1>1) then
!  Pack like symmetric-storage cartesian stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       restr(ii)=d2nl(1,idir1,ipert1,idir,ipert)
     end do
   end do
!  Do the symmetrization using the ground state routine
   call stresssym(gprimd,nsym1,restr,symrc1)
!  Unpack symmetrized stress tensor
   ii=0
   do ipert1=natom+3,natom+4
     do idir1=1,3
       ii=ii+1
       d2nl(1,idir1,ipert1,idir,ipert)=restr(ii)
     end do
   end do
 end if !nsym>1

!----------------------------------------------------------------------------
!Now, treat the local contribution

 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 n3xccc=0
 if(psps%n1xccc/=0)n3xccc=nfft
 ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
 ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))
 ABI_ALLOCATE(vhartr01,(nfft))
 xccc3d1(:)=zero

!Double loop over strain perturbations
 do ipert1=natom+3,natom+4
   do idir1=1,3
     if(ipert1==natom+3) then
       istr1=idir1
     else
       istr1=idir1+3
     end if

!    Get first-order local potential.
     call vlocalstr(gmet,gprimd,gsqcut,istr1,mgfft,mpi_enreg,&
&     psps%mqgrid_vl,natom,gs_hamk%nattyp,nfft,ngfft,ntypat,ph1d,psps%qgrid_vl,&
&     ucvol,psps%vlspl,vpsp1)

!    Get first-order hartree potential.
     call hartrestr(gsqcut,idir1,ipert1,mpi_enreg,natom,nfft,ngfft,&
&     rhog,rprimd,vhartr01)

!    Get first-order exchange-correlation potential
     if(psps%n1xccc/=0)then
       call dfpt_mkcore(cplex,idir1,ipert1,natom,ntypat,n1,psps%n1xccc,&
&       n2,n3,qphon,rprimd,typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
     end if ! psps%n1xccc/=0

     option=0
     call dfpt_mkvxcstr(cplex,idir1,ipert1,kxc,mpi_enreg,natom,nfft,ngfft,&
&     dummy,dummy,nkxc,nmxc,nspden,n3xccc,option,qphon,rhor,rhor1,rprimd,&
&     0,0,vxc1,xccc3d1)

!    Combines density j2 with local potential j1
     do ispden=1,min(nspden,2)
       do ifft=1,cplex*nfft
         vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)+vhartr01(ifft)
       end do
     end do
     call dotprod_vn(cplex,rhor1,dotr,doti,nfft,nfftot,nspden,2,vxc1,ucvol)
     write(std_out,*)
     d2lo(1,idir1,ipert1,idir,ipert)=dotr
     d2lo(2,idir1,ipert1,idir,ipert)=doti
   end do ! istr1
 end do ! ipert1

 call gs_hamk%free()

 ABI_DEALLOCATE(vxc1)
 ABI_DEALLOCATE(xccc3d1)
 ABI_DEALLOCATE(vhartr01)

 ABI_DEALLOCATE(d2bbb_k)
 ABI_DEALLOCATE(d2nl_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(vpsp1)

end subroutine dfpt_nselt
!!***

!!****f* ABINIT/dfpt_nsteltwf
!! NAME
!! dfpt_nsteltwf
!!
!! FUNCTION
!! This routine computes the non-local and kinetic contribution to the
!! 2DTE matrix elements, in the non-stationary formulation
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)=pw coefficients of RF
!!    wavefunctions at k,q.
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)  (NOT NEEDED !)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  gs_hamk <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k
!!  icg=shift to be applied on the location of data in the array cg
!!  icg1=shift to be applied on the location of data in the array cg1
!!  ikpt=number of the k-point
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  istwf_k=flag controlling the storage of WFs
!!  kg_k(3,npw_k)=reduced planewave coordinates.
!!  kg1_k(3,npw1_k)=reduced planewave coordinates at k+q, with RF k points
!!  kpoint(3)=k-point in reduced coordinates
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  mkmem =number of k points treated by this node.
!!  mk1mem =number of k points treated by this node (RF data).
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  natom=number of atoms in cell.
!!  nband_k=number of bands at this k point for that spin polarization
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  ntypat=number of types of atoms in unit cell.
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  wtk_k=weight assigned to the k point.
!!  ylm(npw_k,mpsang*mpsang)= real spherical harmonics for each G and k point
!!  ylmgr(npw_k,3,mpsang*mpsang*useylm)= gradients of real spherical for each G and k point
!!
!! OUTPUT
!!  d2nl_k(2,3,mpert)=non-local contributions to
!!   non-stationary 2DTE, for the present k point, and perturbation idir, ipert
!!
!! PARENTS
!!      dfpt_nselt
!!
!! CHILDREN
!!      dotprod_g,kpgstr,mkffnl,mkkin,nonlop
!!
!! SOURCE

subroutine dfpt_nsteltwf(cg,cg1,d2nl_k,ecut,ecutsm,effmass_free,gs_hamk,icg,icg1,ikpt,isppol,&
&  istwf_k,kg_k,kg1_k,kpoint,mband,mband_mem,mkmem,mk1mem,mpert,mpi_enreg,mpw,mpw1,natom,nband_k,&
&  npw_k,npw1_k,nspinor,nsppol,ntypat,occ_k,psps,rmet,wtk_k,ylm,ylmgr)



!Arguments ------------------------------------
!scalars
 integer,intent(in) :: icg,icg1,ikpt,isppol,istwf_k,mband,mk1mem,mkmem,mpert,mpw,mpw1,natom
 integer,intent(in) :: mband_mem
 integer,intent(in) :: nspinor,nsppol,ntypat
 integer,intent(inout) :: nband_k,npw1_k,npw_k
 real(dp),intent(in) :: ecut,ecutsm,effmass_free,wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: kg1_k(3,npw1_k),kg_k(3,npw_k)
 real(dp),intent(in) :: kpoint(3)
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)
 real(dp),intent(in) :: occ_k(nband_k),rmet(3,3)
 real(dp),intent(in) :: ylm(npw_k,psps%mpsang*psps%mpsang)
 real(dp),intent(in) :: ylmgr(npw_k,3,psps%mpsang*psps%mpsang)
 real(dp),intent(inout) :: d2nl_k(2,3,mpert)

!Local variables-------------------------------
!scalars
 integer :: choice,cpopt,dimffnl,dimffnl2,iband
 integer :: iband_me
 integer :: ider,idir0,idir1,ilmn,ipert1,ipw,ipws,ispinor,istr1,itypat
 integer :: nkpg,nnlout,paw_opt,signs,tim_nonlop
 real(dp) :: doti,dotr
 type(gs_hamiltonian_type) :: gs_hamk
!arrays
 real(dp) :: enlout(6),dum_svectout(1,1),dum(1),kpg_dum(0,0)
 real(dp),allocatable :: cwave0(:,:),cwavef(:,:),dkinpw(:),eig2_k(:)
 real(dp),allocatable :: ffnl(:,:,:,:),ffnl_ylm(:,:,:,:),ghc(:,:)
 real(dp),allocatable :: gvnlx1(:,:),gvnlxc(:,:),kinpw1(:),ph3d(:,:,:)
 type(pawcprj_type) :: cprj_dum(0,0)

! *********************************************************************

!Init me
 ABI_ALLOCATE(ghc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlxc,(2,npw1_k*nspinor))
 ABI_ALLOCATE(gvnlx1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(eig2_k,(2*nsppol*mband**2))
 ABI_ALLOCATE(kinpw1,(npw1_k))
 ABI_ALLOCATE(dkinpw,(npw_k))
 nkpg=0

!Compute nonlocal form factors ffnl at (k+G), for all atoms
 dimffnl=2
 ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
 if (psps%useylm==0) then
   ider=1;idir0=0
   call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,ider,idir0,&
&   psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,&
&   npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
 else
   ider=1;idir0=-7;dimffnl2=7
   ABI_ALLOCATE(ffnl_ylm,(npw_k,dimffnl2,psps%lmnmax,ntypat))
   call mkffnl(psps%dimekb,dimffnl2,psps%ekb,ffnl_ylm,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&   ider,idir0,psps%indlmn,kg_k,kpg_dum,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&   nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm,ylmgr)
   do itypat=1,ntypat
     do ilmn=1,psps%lmnmax
       ffnl(:,1,ilmn,itypat)=ffnl_ylm(:,1,ilmn,itypat)
     end do
   end do
 end if

!Compute kinetic contributions (1/2) (2 Pi)**2 (k+G)**2:
! call mkkin(ecut,ecutsm,effmass_free,gs_hamk%gmet,kg1_k,kinpw1,kpoint,npw1_k)
 call mkkin(ecut,ecutsm,effmass_free,gs_hamk%gmet,kg1_k,kinpw1,kpoint,npw1_k,0,0)

!Load k/k+q-dependent part in the Hamiltonian datastructure
 ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
 call gs_hamk%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,ffnl_k=ffnl,&
& ph3d_k=ph3d,compute_ph3d=.true.)

 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwavef,(2,npw1_k*nspinor))

!Loop over bands
 iband_me = 0
 do iband=1,nband_k

   if(mpi_enreg%proc_distrb(ikpt, iband, isppol) /= mpi_enreg%me_kpt) then
!    Skip the eigenvalue and the gvnl records of this band
     cycle
   end if
   iband_me = iband_me + 1

!  Get ground-state and first-order wavefunctions
   cwave0(:,:)=cg(:,1+(iband_me-1)*npw_k*nspinor+icg:iband*npw_k*nspinor+icg)
   cwavef(:,:)=cg1(:,1+(iband_me-1)*npw1_k*nspinor+icg1:iband*npw1_k*nspinor+icg1)

!  Double loop over strain perturbations
   do ipert1=natom+3,natom+4
     do idir1=1,3
       if (ipert1==natom+3) istr1=idir1
       if (ipert1==natom+4) istr1=idir1+3

!      Compute the derivative of the kinetic operator vs strain in dkinpw
       call kpgstr(dkinpw,ecut,ecutsm,effmass_free,gs_hamk%gmet,gs_hamk%gprimd,istr1,&
&       kg1_k,kpoint,npw1_k)

!      Get |vnon-locj1|u0> :
!      first-order non-local, applied to zero-order wavefunction
!      (??) this routine gives MINUS the non-local contribution

!      When using Ylms, load the correct ffnl derivative
       if (psps%useylm==1) then
         do itypat=1,ntypat
           do ilmn=1,psps%lmnmax
             ffnl(:,2,ilmn,itypat)=ffnl_ylm(:,1+istr1,ilmn,itypat)
           end do
         end do
       end if

       signs=2 ; choice=3 ; nnlout=6 ; paw_opt=0 ; cpopt=-1 ; tim_nonlop=5
       call nonlop(choice,cpopt,cprj_dum,enlout,gs_hamk,istr1,dum,mpi_enreg,1,nnlout,paw_opt,&
&       signs,dum_svectout,tim_nonlop,cwave0,gvnlx1)
!      <G|Vnl1|Cnk> is contained in gvnlx1

!      Kinetic contribution
       do ispinor=1,nspinor
         do ipw=1,npw1_k
           ipws=ipw+npw1_k*(ispinor-1)
           if(kinpw1(ipw)<huge(0.0_dp)*1.d-11)then
             gvnlx1(1,ipws)=gvnlx1(1,ipws)+dkinpw(ipw)*cwave0(1,ipws)
             gvnlx1(2,ipws)=gvnlx1(2,ipws)+dkinpw(ipw)*cwave0(2,ipws)
           else
             gvnlx1(1,ipws)=0.0_dp
             gvnlx1(2,ipws)=0.0_dp
           end if
         end do
       end do

!      construct the matrix element (<uj2|vj1|u0>)complex conjug.
!      and add it to the 2nd-order matrix
!      imaginary term should be zero for strain-strain 2nd derivatives,
!      but keep it as a test for now
       call dotprod_g(dotr,doti,gs_hamk%istwf_k,npw1_k*nspinor,2,cwavef,gvnlx1,&
&         mpi_enreg%me_g0,mpi_enreg%comm_spinorfft)

print *, 'nselt d2nl_k, wtk_k, occ_k, dotr, doti ', d2nl_k, ' %%% ', wtk_k, ' %%% ', occ_k, ' %%% ', dotr, ' %%% ', doti
       d2nl_k(1,idir1,ipert1)= d2nl_k(1,idir1,ipert1)+wtk_k*occ_k(iband)*2.0_dp*dotr
       d2nl_k(2,idir1,ipert1)= d2nl_k(2,idir1,ipert1)-wtk_k*occ_k(iband)*2.0_dp*doti

     end do !idir1
   end do !ipert1

!  UNTIL NOW, DO NOT TAKE INTO ACCOUNT istwf_k

!  End loop over bands
 end do

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwavef)

!###################################################################

 ABI_DEALLOCATE(eig2_k)
 ABI_DEALLOCATE(ghc)
 ABI_DEALLOCATE(gvnlxc)
 ABI_DEALLOCATE(gvnlx1)
 ABI_DEALLOCATE(kinpw1)
 ABI_DEALLOCATE(dkinpw)
 ABI_DEALLOCATE(ffnl)
 ABI_DEALLOCATE(ph3d)
 if (psps%useylm==1)  then
   ABI_DEALLOCATE(ffnl_ylm)
 end if

end subroutine dfpt_nsteltwf
!!***


!!****f* ABINIT/dfpt_nstdy
!! NAME
!! dfpt_nstdy
!!
!! FUNCTION
!! This routine compute the non-stationary expression for the
!! second derivative of the total energy, for a whole row of
!! mixed derivatives.
!! Only for norm-conserving pseudopotentials (no PAW)
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions at k
!!  cg1(2,mpw1*nspinor*mband_mem*mk1mem*nsppol)=pw coefficients of RF wavefunctions at k,q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL, if 2, COMPLEX
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  idir=direction of the perturbation
!!  indkpt1(nkpt_rbz)=non-symmetrized indices of the k-points
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points in the reduced BZ
!!  kxc(nfft,nkxc)=exchange and correlation kernel
!!  mkmem =number of k points treated by this node (GS data)
!!  mk1mem =number of k points treated by this node (RF data)
!!  mpert =maximum number of ipert
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  nattyp(ntypat)= # atoms of each type.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  nfft=(effective) number of FFT grid points (for this proc)
!!  ngfft(1:18)=integer array with FFT box dimensions and other
!!  nkpt=number of k points in the full BZ
!!  nkpt_rbz=number of k points in the reduced BZ for this perturbation
!!  nkxc=second dimension of the kxc array. If /=0, the XC kernel must be computed.
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with i perturbation
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band
!!   and k in the reduced Brillouin zone (usually =2)
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhor1(cplex*nfft,nspden)=RF electron density in electrons/bohr**3.
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume in bohr**3.
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point in the reduced BZ
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k+q point
!!
!! OUTPUT
!!  blkflg(3,mpert,3,mpert)=flags for each element of the 2DTE (=1 if computed)
!!  d2bbb(2,3,3,mpert,mband,mband*prtbbb)=band by band decomposition of some
!!                                        second order derivatives
!!  d2lo(2,3,mpert,3,mpert)=local contributions to the 2DTEs
!!  d2nl(2,3,mpert,3,mpert)=non-local contributions to the 2DTEs
!!
!! NOTES
!! Note that the ddk perturbation should not be treated here.
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      appdig,dfpt_mkcore,dfpt_mkvxc,dfpt_mkvxc_noncoll
!!      dfpt_nstwf,dfpt_sygra,dfpt_vlocal,dotprod_vn,init_hamiltonian
!!      mati3inv,timab,wfk_close,wfk_open_read,wrtout
!!      xmpi_sum
!!
!! SOURCE

subroutine dfpt_nstdy(atindx,blkflg,cg,cg1,cplex,dtfil,dtset,d2bbb,d2lo,d2nl,eigen0,eigen1,&
&          gmet,gsqcut,idir,indkpt1,indsy1,ipert,istwfk_rbz,kg,kg1,kpt_rbz,kxc,mkmem,mk1mem,&
&          mpert,mpi_enreg,mpw,mpw1,nattyp,nband_rbz,nfft,ngfft,nkpt,nkpt_rbz,nkxc,&
&          npwarr,npwar1,nspden,nsppol,nsym1,occ_rbz,ph1d,psps,rhor1,rmet,rprimd,&
&          symrc1,ucvol,wtk_rbz,xred,ylm,ylm1,rhor,vxc)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mk1mem,mkmem,mpert,mpw,mpw1,nfft,nkpt,nkpt_rbz,nkxc,nspden,nsppol,nsym1
 real(dp),intent(in) :: gsqcut,ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: atindx(dtset%natom),indkpt1(nkpt_rbz),indsy1(4,nsym1,dtset%natom)
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nattyp(dtset%ntypat),nband_rbz(nkpt_rbz*nsppol),ngfft(18)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz),symrc1(3,3,nsym1)
 integer,intent(inout) :: blkflg(3,mpert,3,mpert) !vz_i
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*dtset%mband_mem*mkmem*nsppol)
 real(dp),intent(in) :: cg1(2,mpw1*dtset%nspinor*dtset%mband_mem*mk1mem*nsppol)
 real(dp),intent(in) :: eigen0(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen1(2*dtset%mband*dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: gmet(3,3),kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: kxc(nfft,nkxc),occ_rbz(dtset%mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhor1(cplex*nfft,nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: wtk_rbz(nkpt_rbz),xred(3,dtset%natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: d2bbb(2,3,3,mpert,dtset%mband,dtset%mband*dtset%prtbbb)!vz_i
 real(dp),intent(inout) :: d2lo(2,3,mpert,3,mpert),d2nl(2,3,mpert,3,mpert) !vz_i
! optional
 real(dp),optional,intent(in) :: rhor(nfft,nspden)
 real(dp),optional,intent(in) :: vxc(cplex*nfft,nspden)

!Local variables-------------------------------
!scalars
 integer,parameter :: formeig1=1
 integer :: ban2tot,bantot,bdtot_index,ddkcase,iband,icg,icg1,idir1
 integer :: ierr,ifft,ii,ikg,ikg1,ikpt,ilm,ipert1,ispden,isppol
 integer :: istwf_k,isym,jj,master,me,n1,n2,n3,n3xccc,n4,n5,n6
 integer :: nband_k,nfftot,npw1_k,npw_k,nspinor_,option,spaceworld,optnc
 real(dp) :: doti,dotr,wtk_k
 logical :: nmxc=.false.,t_exist
 character(len=500) :: msg
 character(len=fnlen) :: fiwfddk
 type(gs_hamiltonian_type) :: gs_hamkq
!arrays
 integer :: ddkfil(3)
 integer,allocatable :: kg1_k(:,:),kg_k(:,:),symrl1(:,:,:)
 real(dp) :: d2nl_elfd(2,3),d2nl_mgfd(2,3),kpoint(3),kpq(3),sumelfd(2),summgfd(2),tsec(2)
 real(dp),allocatable :: buffer1(:),buffer2(:),d2bbb_k(:,:,:,:),d2nl_k(:,:,:)
 real(dp),allocatable :: eig1_k(:),eig_k(:),occ_k(:)
 real(dp) :: rhodummy(0,0)
 real(dp),allocatable :: vpsp1(:),vxc1(:,:),work1(:,:,:),xccc3d1(:),ylm1_k(:,:),ylm_k(:,:)
 type(pawtab_type) :: pawtab(dtset%ntypat*psps%usepaw)
 type(wfk_t) :: ddks(3)


! *********************************************************************

 ABI_UNUSED(nkpt)

 DBG_ENTER("COLL")

!Not valid for PAW
 if (psps%usepaw==1) then
   msg='This routine cannot be used for PAW (use dfpt_nstpaw instead) !'
   MSG_BUG(msg)
 end if


!Keep track of total time spent in dfpt_nstdy
 call timab(101,1,tsec)

!Init parallelism
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt

 master =0

!Zero only portion of nonlocal matrix to be computed here
 d2nl(:,:,1:dtset%natom+2,idir,ipert)=zero

 ABI_ALLOCATE(d2bbb_k,(2,3,dtset%mband,dtset%mband*dtset%prtbbb))
 ABI_ALLOCATE(d2nl_k,(2,3,mpert))
 ABI_ALLOCATE(eig_k,(nsppol*dtset%mband))
 ABI_ALLOCATE(eig1_k,(2*nsppol*dtset%mband**2))
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))

!Do not try to open electric field file
 ddkfil(:)=0
!The treatment of homogeneous electric field potential need the existence of d/dk files.
 do idir1=1,3
   ddkcase=idir1+dtset%natom*3
   call appdig(ddkcase,dtfil%fnamewffddk,fiwfddk)

!  Check that ddk file exists
   t_exist = file_exists(fiwfddk)
   if (.not. t_exist) then
     ! Try netcdf file.
     t_exist = file_exists(nctk_ncify(fiwfddk))
     if (t_exist) then
       fiwfddk = nctk_ncify(fiwfddk)
       write(msg,"(3a)")"- File: ",trim(fiwfddk)," does not exist but found netcdf file with similar name."
       call wrtout(std_out,msg,'COLL')
     end if
   end if

   if (t_exist) then
!    Note the use of unit numbers 21, 22 and 23
     ddkfil(idir1)=20+idir1
     write(msg, '(a,a)') '-open ddk wf file :',trim(fiwfddk)
     call wrtout(std_out,msg,'COLL')
     call wrtout(ab_out,msg,'COLL')
     call wfk_open_read(ddks(idir1),fiwfddk,formeig1,dtset%iomode,ddkfil(idir1),spaceworld)
   end if
 end do

!Update list of computed matrix elements
 if (ipert /= dtset%natom + 1) then
   do ipert1=1,mpert
     do idir1=1,3
       if(ipert1 <= dtset%natom .or. ipert1==dtset%natom+2 .and. ddkfil(idir1)/=0) then
         blkflg(idir1,ipert1,idir,ipert)=1
       end if
     end do
   end do
 else
   ipert1 = dtset%natom + 1
   do idir1=1,3
!    If was already computed in another run or dataset, or if is to be computed in the present one
     if ((ddkfil(idir1) /= 0).or. (dtset%rfdir(idir1)/=0.and. idir1<=idir) ) then
!      if ((ddkfil(idir1) /= 0).or. (idir1==idir) ) then
       blkflg(idir1,ipert1,idir,ipert)=1
     end if
   end do
 end if

 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 nspinor_=dtset%nspinor

 bantot = 0
 ban2tot = 0

!==== Initialize most of the Hamiltonian ====
!1) Allocate all arrays and initialize quantities that do not depend on k and spin.
!2) Perform the setup needed for the non-local factors:
!3) Constant kleimann-Bylander energies are copied from psps to gs_hamk.
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nsppol,nspden,dtset%natom,&
& dtset%typat,xred,nfft,dtset%mgfft,ngfft,rprimd,dtset%nloalg,ph1d=ph1d,&
& use_gpu_cuda=dtset%use_gpu_cuda)

!LOOP OVER SPINS
 bdtot_index=0
 icg=0;icg1=0
 do isppol=1,nsppol

   ikg=0;ikg1=0

!  Continue to initialize the Hamiltonian
   call gs_hamkq%load_spin(isppol,with_nonlocal=.true.)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz

     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)

     eig_k(1:nband_k) = eigen0(1+bantot:nband_k+bantot)
     eig1_k(1:2*nband_k**2) = eigen1(1+ban2tot:2*nband_k**2+ban2tot)
     bantot = bantot + nband_k
     ban2tot = ban2tot + 2*nband_k**2

print *, 'nstdy cycle, ik,isppol, nband ', proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me),&
&        ikpt, isppol, proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       bdtot_index=bdtot_index+nband_k
!      The wavefunction blocks for ddk file is skipped elsewhere in the loop
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))

!    In case of electric field pert1, read ddk wfs file
!    Note that the symmetries are not used for ddk, so read each k point
!    Also take into account implicitely the parallelism over k points

     do idir1=1,3
       if (ddkfil(idir1)/=0) then
         ii = ddks(idir1)%findk(kpt_rbz(:, ikpt))
         ABI_CHECK(ii == indkpt1(ikpt),  "ii !=  indkpt1")
       end if
     end do

     ABI_ALLOCATE(occ_k,(nband_k))
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:)+dtset%qptn(:)
     wtk_k=wtk_rbz(ikpt)
     d2nl_k(:,:,:)=zero
     if(dtset%prtbbb==1)d2bbb_k(:,:,:,:)=zero

!    Get plane-wave vectors and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get plane-wave vectors and related data at k+q
     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
     end if

!    Compute the eigenvalues, wavefunction,
!    contributions to kinetic energy, nonlocal energy, forces,
!    and update of rhor1 to this k-point and this spin polarization.
!    Note that dfpt_nstwf is called with kpoint, while kpt is used inside dfpt_vtowfk
     call dfpt_nstwf(cg,cg1,ddkfil,dtset,d2bbb_k,d2nl_k,eig_k,eig1_k,gs_hamkq,&
&     icg,icg1,idir,ikpt,ipert,isppol,istwf_k,kg_k,kg1_k,kpoint,kpq,mkmem,mk1mem,mpert,&
&     mpi_enreg,mpw,mpw1,nband_k,npw_k,npw1_k,nsppol,&
&     occ_k,psps,rmet,ddks,wtk_k,ylm_k,ylm1_k)

     d2nl(:,:,:,idir,ipert)=d2nl(:,:,:,idir,ipert)+d2nl_k(:,:,:)
     if(dtset%prtbbb==1)d2bbb(:,:,idir,ipert,:,:)=d2bbb(:,:,idir,ipert,:,:)+d2bbb_k(:,:,:,:)

!    Keep track of total number of bands
     bdtot_index=bdtot_index+nband_k

!    Shift arrays memory
     if (mkmem/=0) then
       icg=icg+npw_k*dtset%nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
       ikg=ikg+npw_k
     end if
     if (mk1mem/=0) then
       icg1=icg1+npw1_k*dtset%nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
       ikg1=ikg1+npw1_k
     end if

     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)

!    End big k point loop
   end do

!  End loop over spins
 end do

 call gs_hamkq%free()

!Treat fixed occupation numbers (as in vtorho)
 if(xmpi_paral==1)then
   ABI_ALLOCATE(buffer1,(2*3*mpert))
   ABI_ALLOCATE(buffer2,(2*3*mpert))
!  Pack d2nl
   buffer1(1:2*3*mpert)=reshape(d2nl(:,:,:,idir,ipert),(/2*3*mpert/))
!  Build sum of everything
   call timab(48,1,tsec)
   call xmpi_sum(buffer1,buffer2,2*3*mpert,spaceworld,ierr)
   call timab(48,2,tsec)
!  Unpack the final result
   d2nl(:,:,:,idir,ipert)=reshape(buffer2(:),(/2,3,mpert/))
   ABI_DEALLOCATE(buffer1)
   ABI_DEALLOCATE(buffer2)
   if(dtset%prtbbb==1)then
     ABI_ALLOCATE(buffer1,(2*3*dtset%mband*dtset%mband))
     ABI_ALLOCATE(buffer2,(2*3*dtset%mband*dtset%mband))
!    Pack d2bbb
     buffer1(1:2*3*dtset%mband*dtset%mband)=reshape(d2bbb(:,:,idir,ipert,:,:),(/2*3*dtset%mband*dtset%mband/))
!    Build sum of everything
     call timab(48,1,tsec)
     call xmpi_sum(buffer1,buffer2,2*3*dtset%mband*dtset%mband,spaceworld,ierr)
     call timab(48,2,tsec)
!    Unpack the final result
     d2bbb(:,:,idir,ipert,:,:)=reshape(buffer2(:),(/2,3,dtset%mband,dtset%mband/))
     ABI_DEALLOCATE(buffer1)
     ABI_DEALLOCATE(buffer2)
   end if
 end if ! xmpi_paral==1

!In the case of the strain perturbation time-reversal symmetry will always
!be true so imaginary part of d2nl will be must be set to zero here since
!the symmetry-reduced kpt set will leave a non-zero imaginary part.
 if(ipert==dtset%natom+3 .or. ipert==dtset%natom+4) d2nl(2,:,:,idir,ipert)=zero

!In case of electric field ipert1, close the ddk wf files
 do idir1=1,3
   if (ddkfil(idir1)/=0)then
     call ddks(idir1)%close()
   end if
 end do

!Symmetrize the non-local contributions,
!as was needed for the forces in a ground-state calculation
!However, here the quantity is complex, and there are phases !

!Do the transform
 ABI_ALLOCATE(work1,(2,3,dtset%natom))
 do ipert1=1,dtset%natom
   do idir1=1,3
     work1(1,idir1,ipert1)=d2nl(1,idir1,ipert1,idir,ipert)
     work1(2,idir1,ipert1)=d2nl(2,idir1,ipert1,idir,ipert)
   end do
 end do
 call dfpt_sygra(dtset%natom,d2nl(:,:,:,idir,ipert),work1,indsy1,ipert,nsym1,dtset%qptn,symrc1)
 ABI_DEALLOCATE(work1)

!Must also symmetrize the electric/magnetic field perturbation response !
!(XG 000803 This was not implemented until now)
 if(sum(ddkfil(:))/=0)then
!  Get the symmetry matrices in terms of real space basis
   ABI_ALLOCATE(symrl1,(3,3,nsym1))
   do isym=1,nsym1
     call mati3inv(symrc1(:,:,isym),symrl1(:,:,isym))
   end do
!  There should not be any imaginary part, but stay general (for debugging)
   d2nl_elfd(:,:)=d2nl(:,:,dtset%natom+2,idir,ipert)
   do ii=1,3
     sumelfd(:)=zero
     summgfd(:)=zero
     do isym=1,nsym1
       do jj=1,3
         if(symrl1(ii,jj,isym)/=0)then
           if(ddkfil(jj)==0)then
             blkflg(ii,dtset%natom+2,idir,ipert)=0
           end if
         end if
       end do
       sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&       dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&       dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
       summgfd(:)=summgfd(:)+dble(symrl1(ii,1,isym))*d2nl_mgfd(:,1)+&
&       dble(symrl1(ii,2,isym))*d2nl_mgfd(:,2)+&
&       dble(symrl1(ii,3,isym))*d2nl_mgfd(:,3)
     end do
     d2nl(:,ii,dtset%natom+2,idir,ipert)=sumelfd(:)/dble(nsym1)
   end do

   if ((dtset%prtbbb==1).and.(ipert<=dtset%natom)) then
     do iband = 1,dtset%mband
       d2nl_elfd(:,:)=d2bbb(:,:,idir,ipert,iband,iband)
       do ii=1,3
         sumelfd(:)=zero
         do isym=1,nsym1
           sumelfd(:)=sumelfd(:)+dble(symrl1(ii,1,isym))*d2nl_elfd(:,1)+&
&           dble(symrl1(ii,2,isym))*d2nl_elfd(:,2)+&
&           dble(symrl1(ii,3,isym))*d2nl_elfd(:,3)
         end do
         d2bbb(:,ii,idir,ipert,iband,iband)=sumelfd(:)/dble(nsym1)
       end do
     end do  !iband
   end if

   ABI_DEALLOCATE(symrl1)
 end if

!----------------------------------------------------------------------------
!Now, treat the local contribution

 nfftot=ngfft(1)*ngfft(2)*ngfft(3)
 ABI_ALLOCATE(vpsp1,(cplex*nfft))
 if (ipert /= dtset%natom + 1) then
   n3xccc=0;if(psps%n1xccc/=0) n3xccc=nfft
   ABI_ALLOCATE(xccc3d1,(cplex*n3xccc))
   ABI_ALLOCATE(vxc1,(cplex*nfft,nspden))

   do ipert1=1,mpert
     do idir1=1,3
       if(ipert1 <= dtset%natom)then

!        Get first-order local potential and first-order pseudo core density
         call dfpt_vlocal(atindx,cplex,gmet,gsqcut,idir1,ipert1,mpi_enreg,psps%mqgrid_ff,dtset%natom,&
&         nattyp,nfft,ngfft,dtset%ntypat,n1,n2,n3,ph1d,psps%qgrid_ff,&
&         dtset%qptn,ucvol,psps%vlspl,vpsp1,xred)
         if(psps%n1xccc/=0)then
           call dfpt_mkcore(cplex,idir1,ipert1,dtset%natom,dtset%ntypat,n1,psps%n1xccc,&
&           n2,n3,dtset%qptn,rprimd,dtset%typat,ucvol,psps%xcccrc,psps%xccc1d,xccc3d1,xred)
         end if

!        Get first-order exchange-correlation potential (core-correction contribution only !)
         if(psps%n1xccc/=0)then
           option=0
!FR SPr EB non-collinear magnetism
           if (nspden==4.and.present(rhor).and.present(vxc)) then
             optnc=1
             call dfpt_mkvxc_noncoll(cplex,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,rhodummy,0,rhodummy,0,rhodummy,0,&
&             nkxc,nmxc,nspden,n3xccc,optnc,option,dtset%qptn,rhor,rhor1,&
&             rprimd,0,vxc,vxc1,xccc3d1)
           else
             call dfpt_mkvxc(cplex,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,rhodummy,0,rhodummy,0,&
&             nkxc,nmxc,nspden,n3xccc,option,dtset%qptn,rhodummy,&
&             rprimd,0,vxc1,xccc3d1)
           end if
         else
           vxc1(:,:)=zero
         end if

!        Norm-conserving pseudpopotential case:
!        Combines density j2 with local potential j1 (vpsp1 and vxc1)
!        XG030514 : this is a first possible coding, however, each dotprod contains
!        a parallel section (reduction), so it is better to use only one dotprod ...
!        call dotprod_vn(cplex,rhor1,dr_psp1,di_psp1,mpi_enreg,nfft,nfftot,1,2,vpsp1,ucvol)
!        call dotprod_vn(cplex,rhor1,dr_xc1,di_xc1,mpi_enreg,nfft,nfftot,nspden,2,vxc1,ucvol)
!        dotr=dr_psp1+dr_xc1;doti=di_psp1+di_xc1... but then, one needs to overload vxc1
         do ispden=1,min(nspden,2)
           do ifft=1,cplex*nfft
             vxc1(ifft,ispden)=vxc1(ifft,ispden)+vpsp1(ifft)
           end do
         end do
         call dotprod_vn(cplex,rhor1,dotr,doti,nfft,nfftot,nspden,2,vxc1,ucvol)

!        MVeithen 021212 : in case ipert = 2, these lines compute the local part
!        of the Born effective charges from phonon and electric
!        field type perturbations, see eq. 43 of
!        X. Gonze and C. Lee, PRB 55, 10355 (1997) [[cite:Gonze1997a]]
!        The minus sign is due to the fact that the effective charges
!        are minus the second derivatives of the energy
         if (ipert == dtset%natom+2) then
           d2lo(1,idir1,ipert1,idir,ipert)=-dotr
           d2lo(2,idir1,ipert1,idir,ipert)=-doti
         else
           d2lo(1,idir1,ipert1,idir,ipert)=dotr
           d2lo(2,idir1,ipert1,idir,ipert)=doti
         end if
!        Endif ipert1<=natom
       end if
     end do
   end do

   ABI_DEALLOCATE(vxc1)
   ABI_DEALLOCATE(xccc3d1)

 end if ! ipert /= natom +1

 ABI_DEALLOCATE(d2bbb_k)
 ABI_DEALLOCATE(d2nl_k)
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)
 ABI_DEALLOCATE(vpsp1)
 ABI_DEALLOCATE(eig_k)
 ABI_DEALLOCATE(eig1_k)

 call timab(101,2,tsec)

 DBG_EXIT("COLL")

end subroutine dfpt_nstdy
!!***

!!****f* ABINIT/dfpt_rhofermi
!! NAME
!! dfpt_rhofermi
!!
!! FUNCTION
!! This routine computes the fixed contribution to the first-order
!! Fermi energy for metallic occupation and Q=0, as well as the
!! Fermi level charge density needed to compute the remainder of the
!! first-order Fermi energy from the self-consistent local potential
!! at each step in the iteration process.
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions.
!!  cgq(2,mpw1*nspinor*mband_mem*mkqmem*nsppol)=pw coefficients of GS wavefunctions at k+q.
!!  cplex: if 1, real space 1-order functions on FFT grid are REAL; if 2, COMPLEX
!!  cprj(natom,nspinor*mband_mem*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband_mem*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  doccde_rbz(mband*nkpt_rbz*nsppol)=derivative of occ_rbz wrt the energy
!!  docckqde(mband*nkpt_rbz*nsppol)=derivative of occkq wrt the energy
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables for this dataset
!!  eigenq(mband*nkpt_rbz*nsppol)=GS eigenvalues at k+q (hartree)
!!  eigen0(mband*nkpt_rbz*nsppol)=GS eigenvalues at k (hartree)
!!  gmet(3,3)=reciprocal space metric tensor in bohr**-2.
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  idir=direction of the perturbation
!!  indsy1(4,nsym1,natom)=indirect indexing array for atom labels
!!  ipert=type of the perturbation
!!  irrzon1(nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istwfk_rbz(nkpt_rbz)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced planewave coordinates.
!!  kg1(3,mpw1*mk1mem)=reduced planewave coordinates at k+q, with RF k points
!!  kpt_rbz(3,nkpt_rbz)=reduced coordinates of k points.
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  mkmem =number of k points treated by this node (GS data)
!!  mkqmem =number of k+q points treatede by this node (GS data)
!!  mk1mem =number of k points treated by this node.
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  mpw1=maximum dimensioned size of npw for wfs at k+q (also for 1-order wfs).
!!  maximum dimension for q points in grids for nonlocal form factors
!!  natom=number of atoms in cell.
!!  nband_rbz(nkpt_rbz*nsppol)=number of bands at each RF k point for each spin
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs - see comment in respfn.F90)
!!  ngfftf(1:18)=integer array with FFT box dimensions and other for the "fine" grid
!!  nhatfermi(nfft,nspden)=array for fermi-level compensation charge density (PAW only)
!!  nkpt_rbz=number of k points in the IBZ for this perturbation
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  npwarr(nkpt_rbz)=number of planewaves in basis at this GS k point
!!  npwar1(nkpt_rbz)=number of planewaves in basis at this RF k+q point
!!  nspden=number of spin-density components
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym1=number of symmetry elements in space group consistent with
!!    perturbation
!!  occkq(mband*nkpt_rbz*nsppol)=occupation number for each band (often 2)
!!   at each k+q point of the reduced Brillouin zone.
!!  occ_rbz(mband*nkpt_rbz*nsppol)=occupation number for each band and k (usually 2)
!!  paw_ij(natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels for the GS
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawang1 <type(pawang_type)>=pawang datastr. containing only symmetries preserving the perturbation
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid for the GS
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  phnons1(2,nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  prtvol=control print volume and debugging output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rmet(3,3)=real space metric (bohr**2)
!!  rprimd(3,3)=dimensional real space primitive translations
!!  symaf1(nsym1)=(anti)ferromagnetic part of symmetry operations
!!  symrc1(3,3,nsym1)=symmetry operations in reciprocal space
!!  symrl1(3,3,nsym1)=3x3 matrices of the group symmetries
!!  tnons1(3,nsym1)=non-symmorphic translations
!!  ucvol=volume of the unit cell
!!  usecprj= 1 if cprj, cprjq, cprj1 arrays are stored in memory
!!  useylmgr1= 1 if ylmgr1 array is allocated
!!  vtrial(nfftf,nspden)=GS potential (Hartree).
!!  vxc(nfftf,nspden)=XC potential (Hartree).
!!  wtk_rbz(nkpt_rbz)=weight assigned to each k point.
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylm1(mpw1*mk1mem,mpsang*mpsang*useylm)= spherical harmonics for each G and k+g point
!!  ylmgr1(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics at k+q
!!
!!
!! OUTPUT
!!  eigen1(2*mband*mband*nkpt_rbz*nsppol)=array for holding eigenvalues
!!   (hartree) - only digonal elements computed here
!!  fe1fixed=fixed contribution to the first-order Fermi energy
!!   (nonlocal and kinetic in the case of strain)
!!  nhatfermi(cplex*nfftf,nspden)=fermi-level compensation charge density (PAW only)
!!  rhorfermi(cplex*nfftf,nspden)=fermi-level electronic density
!!
!! NOTES
!!  This routine will NOT work with nspden==4:
!!    at least the use of fftpac should be modified.
!!
!! PARENTS
!!      dfpt_scfcv
!!
!! CHILDREN
!!      destroy_rf_hamiltonian,dfpt_wfkfermi,fftpac
!!      init_hamiltonian,init_rf_hamiltonian,kpgstr
!!      mkffnl,mkkin,mkkpg,occeig,paw_ij_free
!!      paw_ij_init,paw_ij_nullify,pawdijfr,pawmkrho,pawrhoij_alloc
!!      pawrhoij_free,pawrhoij_free_unpacked,pawrhoij_init_unpacked
!!      pawrhoij_mpisum_unpacked,status,symrhg,timab,xmpi_sum
!!
!! SOURCE

subroutine dfpt_rhofermi(cg,cgq,cplex,cprj,cprjq,&
& doccde_rbz,docckqde,dtfil,dtset,eigenq,eigen0,eigen1,fe1fixed,gmet,gprimd,idir,&
& indsy1,ipert,irrzon1,istwfk_rbz,kg,kg1,kpt_rbz,mband,mband_mem,mkmem,mkqmem,mk1mem,mpi_enreg,&
& mpw,mpw1,my_natom,natom,nband_rbz,ncpgr,nfftf,ngfftf,nhatfermi,nkpt_rbz,npwarr,npwar1,nspden,&
& nsppol,nsym1,occkq,occ_rbz,paw_ij,pawang,pawang1,pawfgr,pawfgrtab,pawrad,pawrhoijfermi,pawtab,&
& phnons1,ph1d,prtvol,psps,rhorfermi,rmet,rprimd,symaf1,symrc1,symrl1,tnons1,&
& ucvol,usecprj,useylmgr1,vtrial,vxc,wtk_rbz,xred,ylm,ylm1,ylmgr1)

!Arguments -------------------------------
!scalars
 integer,intent(in) :: cplex,idir,ipert,mband,mk1mem,mkmem,mkqmem
 integer,intent(in) :: mband_mem
 integer,intent(in) :: mpw,mpw1,my_natom,natom,ncpgr,nfftf,nkpt_rbz,nspden,nsppol,nsym1
 integer,intent(in) :: prtvol,usecprj,useylmgr1
 real(dp),intent(in) :: ucvol
 real(dp),intent(out) :: fe1fixed
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang,pawang1
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: indsy1(4,nsym1,natom)
 integer,intent(in) :: irrzon1(dtset%nfft**(1-1/nsym1),2,(nspden/nsppol)-3*(nspden/4))
 integer,intent(in) :: istwfk_rbz(nkpt_rbz),kg(3,mpw*mkmem),kg1(3,mpw1*mk1mem)
 integer,intent(in) :: nband_rbz(nkpt_rbz*nsppol),ngfftf(18)
 integer,intent(in) :: npwar1(nkpt_rbz),npwarr(nkpt_rbz),symaf1(nsym1)
 integer,intent(in) :: symrc1(3,3,nsym1),symrl1(3,3,nsym1)
 real(dp),intent(in) :: cg(2,mpw*dtset%nspinor*mband_mem*mkmem*nsppol)
 real(dp),intent(in) :: cgq(2,mpw1*dtset%nspinor*mband_mem*mkqmem*nsppol)
 real(dp),intent(in) :: doccde_rbz(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: docckqde(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigen0(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: eigenq(mband*nkpt_rbz*nsppol),gmet(3,3),gprimd(3,3)
 real(dp),intent(in) :: kpt_rbz(3,nkpt_rbz)
 real(dp),intent(in) :: occ_rbz(mband*nkpt_rbz*nsppol),occkq(mband*nkpt_rbz*nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*natom)
 real(dp),intent(in) :: phnons1(2,dtset%nfft**(1-1/nsym1),(nspden/nsppol)-3*(nspden/4))
 real(dp),intent(in) :: rmet(3,3),rprimd(3,3),tnons1(3,nsym1)
 real(dp),intent(in) :: vtrial(nfftf,nspden),vxc(nfftf,nspden),wtk_rbz(nkpt_rbz)
 real(dp),intent(in) :: xred(3,natom),ylm(mpw*mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylm1(mpw1*mk1mem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr1(mpw1*mk1mem,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1)
 real(dp),intent(out) :: eigen1(2*mband*mband*nkpt_rbz*nsppol)
 real(dp),intent(out) :: nhatfermi(:,:)
 real(dp),intent(out) :: rhorfermi(cplex*nfftf,nspden)
!TODO MJV : PAW
 type(pawcprj_type),intent(in) :: cprj (natom,dtset%nspinor*mband*mkmem *nsppol*usecprj)
 type(pawcprj_type),intent(in) :: cprjq(natom,dtset%nspinor*mband*mkqmem*nsppol*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(dtset%ntypat*psps%usepaw)
 type(pawrhoij_type),target,intent(inout)::pawrhoijfermi(my_natom*psps%usepaw)!vz_i
 type(pawtab_type), intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=17
 integer :: bd2tot_index,bdtot_index,buffer_size,cplex_rhoij
 integer :: dimffnl1,dimffnlk,iatom,iband,ibg,ibgq
 integer :: icg,icgq,ider,idir0,ierr,ii,ikg,ikg1,ikpt,ilm,ilmn,indx
 integer :: ispden,isppol,istr,istwf_k
 integer :: mbd2kpsp,mcgq,mcgq_disk,mcprjq,mcprjq_disk
 integer :: me,n1,n2,n3,n4,n5,n6,nband_k,nkpg,nkpg1,npw1_k,npw_k,nspden_rhoij
 integer :: optfr,qphase_rhoij,spaceworld
 logical :: paral_atom,qne0
 real(dp) :: arg,fe1norm,invfe1norm,wtk_k
 type(gs_hamiltonian_type) :: gs_hamkq
 type(rf_hamiltonian_type) :: rf_hamkq
!arrays
 integer,allocatable :: kg1_k(:,:),kg_k(:,:)
 real(dp) :: kpoint(3),kpq(3),tsec(2)
 real(dp) :: ylmgr_dum(1,1,1)
 real(dp),allocatable :: buffer1(:),dkinpw(:),doccde_k(:)
 real(dp),allocatable :: doccde_kq(:),eig0_k(:),eig0_kq(:),eig1_k(:)
 real(dp),allocatable :: fe1fixed_k(:),fe1norm_k(:)
 real(dp),allocatable :: ffnl1(:,:,:,:),ffnlk(:,:,:,:)
 real(dp),allocatable :: kinpw1(:),kpg1_k(:,:),kpg_k(:,:)
 real(dp),allocatable :: occ_k(:),occ_kq(:),ph3d(:,:,:),ph3d1(:,:,:)
 real(dp),allocatable :: rhoaug(:,:,:),rhogfermi(:,:),rhowfr(:,:)
 real(dp),allocatable :: rhoaug4(:,:,:,:)
 real(dp),allocatable :: rocceig(:,:),ylm1_k(:,:),ylm_k(:,:),ylmgr1_k(:,:,:)
 type(paw_ij_type),allocatable :: paw_ij1fr(:)
 type(pawrhoij_type),pointer :: pawrhoijfermi_unsym(:)
! real(dp),allocatable :: vlocal1(:,:,:,:),vlocal_tmp(:,:,:,:)
! real(dp),allocatable :: v1zeeman(:,:),vtrial_tmp(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Check arguments validity
 if (ipert>natom.and.ipert/=natom+3.and.ipert/=natom+4.and.ipert/=natom+5) then
   MSG_BUG('wrong ipert argument!')
 end if
 if (cplex/=1) then
   MSG_BUG('wrong cplex/=1 argument !')
 end if

!Keep track of total time spent in this routine
 call timab(121,1,tsec)
 call timab(124,1,tsec)

!Retrieve parallelism data
 spaceworld=mpi_enreg%comm_cell
 me=mpi_enreg%me_kpt
 paral_atom=(my_natom/=dtset%natom)

!Initialize output variables
 fe1fixed=zero
 if (psps%usepaw==0) rhorfermi(:,:)=zero

!Initialisations/allocation of temporary variables
 n1=dtset%ngfft(1) ; n2=dtset%ngfft(2) ; n3=dtset%ngfft(3)
 n4=dtset%ngfft(4) ; n5=dtset%ngfft(5) ; n6=dtset%ngfft(6)
 bdtot_index=0 ; bd2tot_index=0 ; ibg=0 ; ibgq=0 ; icg=0 ; icgq=0
 qne0=(dtset%qptn(1)**2+dtset%qptn(2)**2+dtset%qptn(3)**2>=tol14)
 mbd2kpsp=2*mband**2*nkpt_rbz*nsppol
 fe1norm=zero
 if (nspden/=4) then
   ABI_ALLOCATE(rhoaug,(cplex*n4,n5,n6))
 else
   ABI_ALLOCATE(rhoaug4,(cplex*n4,n5,n6,nspden))
 end if
 ABI_ALLOCATE(kg_k,(3,mpw))
 ABI_ALLOCATE(kg1_k,(3,mpw1))
 if (psps%usepaw==1) then
   ABI_ALLOCATE(rhowfr,(cplex*dtset%nfft,dtset%nspden))
   rhowfr(:,:)=zero
 end if

 mcgq=mpw1*dtset%nspinor*mband_mem*mkqmem*nsppol;mcgq_disk=0

!Prepare RF PAW files for reading and writing if mkmem, mkqmem or mk1mem==0
 if (psps%usepaw==1) then
!TODO MJV : PAW
   mcprjq=dtset%nspinor*mband*mkqmem*nsppol*usecprj;mcprjq_disk=0
 else
   mcprjq=0;mcprjq_disk=0
 end if

!PAW:has to compute frozen part of Dij^(1) (without Vpsp(1) contribution)
 if (psps%usepaw==1) then
   ABI_DATATYPE_ALLOCATE(paw_ij1fr,(my_natom))
   call paw_ij_nullify(paw_ij1fr)
   call paw_ij_init(paw_ij1fr,cplex,dtset%nspinor,dtset%nsppol,dtset%nspden,0,&
&   dtset%natom,dtset%ntypat,dtset%typat,pawtab,has_dijfr=1,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom )
   optfr=1
   ABI_ALLOCATE(buffer1,(0))
   call pawdijfr(gprimd,idir,ipert,my_natom,natom,nfftf,ngfftf,dtset%nspden,dtset%nsppol,&
&   dtset%ntypat,optfr,paw_ij1fr,pawang,pawfgrtab,pawrad,pawtab,&
&   cplex,dtset%qptn,rprimd,ucvol,buffer1,vtrial,vxc,xred,&
&   mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom)
   ABI_DEALLOCATE(buffer1)
 end if

!PAW:allocate memory for non-symetrized occupancies matrix at EFermi (pawrhoijfermi)
 pawrhoijfermi_unsym => pawrhoijfermi
 if (psps%usepaw==1) then
   if (paral_atom) then
     ABI_DATATYPE_ALLOCATE(pawrhoijfermi_unsym,(natom))
     !Q phase should be 1 because Q=0
     call pawrhoij_inquire_dim(cplex_rhoij=cplex_rhoij,qphase_rhoij=qphase_rhoij,nspden_rhoij=nspden_rhoij,&
&                              nspden=dtset%nspden,spnorb=dtset%pawspnorb,cplex=cplex,cpxocc=dtset%pawcpxocc)
     call pawrhoij_alloc(pawrhoijfermi_unsym,cplex_rhoij,nspden_rhoij,dtset%nspinor,&
&     dtset%nsppol,dtset%typat,qphase=qphase_rhoij,pawtab=pawtab,use_rhoijp=0,use_rhoij_=1)
   else
     call pawrhoij_init_unpacked(pawrhoijfermi_unsym)
   end if
 end if

!Initialize most of the Hamiltonian (arrays and quantities that do not depend on k + nl form factors)
 call init_hamiltonian(gs_hamkq,psps,pawtab,dtset%nspinor,nsppol,nspden,natom,&
& dtset%typat,xred,dtset%nfft,dtset%mgfft,dtset%ngfft,rprimd,dtset%nloalg,&
& paw_ij=paw_ij,usecprj=usecprj,ph1d=ph1d,use_gpu_cuda=dtset%use_gpu_cuda,&
& mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,mpi_spintab=mpi_enreg%my_isppoltab)
 call init_rf_hamiltonian(cplex,gs_hamkq,ipert,rf_hamkq,paw_ij1=paw_ij1fr,&
& mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,mpi_spintab=mpi_enreg%my_isppoltab)



!LOOP OVER SPINS
 do isppol=1,nsppol
   ikg=0;ikg1=0
!  Continue to initialize the Hamiltonian at k+q
   call gs_hamkq%load_spin(isppol,with_nonlocal=.true.)

   call rf_hamkq%load_spin(isppol,with_nonlocal=.true.)


!  Nullify contribution to density at EFermi from this k-point
   if (nspden/=4) then
     rhoaug(:,:,:)=zero
   else
     rhoaug4(:,:,:,:)=zero
   end if
   call timab(125,1,tsec)

!  BIG FAT k POINT LOOP
   do ikpt=1,nkpt_rbz
     nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
     istwf_k=istwfk_rbz(ikpt)
     npw_k=npwarr(ikpt)
     npw1_k=npwar1(ikpt)
     wtk_k=wtk_rbz(ikpt)

print *, ' cycle, nband ', proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me),&
&        proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me)) then
       eigen1(1+bd2tot_index : 2*nband_k**2+bd2tot_index) = zero
       bdtot_index=bdtot_index+nband_k
       bd2tot_index=bd2tot_index+2*nband_k**2
!      Skip the rest of the k-point loop
       cycle
     end if

     ABI_ALLOCATE(ylm_k,(npw_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylm1_k,(npw1_k,psps%mpsang*psps%mpsang*psps%useylm))
     ABI_ALLOCATE(ylmgr1_k,(npw1_k,3,psps%mpsang*psps%mpsang*psps%useylm*useylmgr1))

!    Continue to initialize the Hamiltonian at k+q
     kpoint(:)=kpt_rbz(:,ikpt)
     kpq(:)=kpoint(:)+dtset%qptn(1:3)

     ABI_ALLOCATE(doccde_k,(nband_k))
     ABI_ALLOCATE(doccde_kq,(nband_k))
     ABI_ALLOCATE(eig0_k,(nband_k))
     ABI_ALLOCATE(eig0_kq,(nband_k))
     ABI_ALLOCATE(eig1_k,(2*nband_k**2))
     ABI_ALLOCATE(fe1fixed_k,(nband_k))
     ABI_ALLOCATE(fe1norm_k,(nband_k))
     ABI_ALLOCATE(occ_k,(nband_k))
     ABI_ALLOCATE(occ_kq,(nband_k))
     ABI_ALLOCATE(rocceig,(nband_k,nband_k))

     eig1_k(:)=zero
     eig0_k(:)=eigen0(1+bdtot_index:nband_k+bdtot_index)
     eig0_kq(:)=eigenq(1+bdtot_index:nband_k+bdtot_index)
     occ_k(:)=occ_rbz(1+bdtot_index:nband_k+bdtot_index)
     occ_kq(:)=occkq(1+bdtot_index:nband_k+bdtot_index)
     doccde_k(:)=doccde_rbz(1+bdtot_index:nband_k+bdtot_index)
     doccde_kq(:)=docckqde(1+bdtot_index:nband_k+bdtot_index)

!    For each pair of active bands (m,n), generates the ratios
!    rocceig(m,n)=(occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n))
!    and decide to which band to attribute it.
     call occeig(doccde_k,doccde_kq,eig0_k,eig0_kq,nband_k,&
&     dtset%occopt,occ_k,occ_kq,rocceig)

!    Get plane-wave coeffs and related data at k
     kg_k(:,1:npw_k)=kg(:,1+ikg:npw_k+ikg)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm_k(1:npw_k,ilm)=ylm(1+ikg:npw_k+ikg,ilm)
       end do
     end if

!    Get plane-wave coeffs and related data at k+q
     kg1_k(:,1:npw1_k)=kg1(:,1+ikg1:npw1_k+ikg1)
     if (psps%useylm==1) then
       do ilm=1,psps%mpsang*psps%mpsang
         ylm1_k(1:npw1_k,ilm)=ylm1(1+ikg1:npw1_k+ikg1,ilm)
       end do
       if (useylmgr1==1) then
         do ilm=1,psps%mpsang*psps%mpsang
           do ii=1,3
             ylmgr1_k(1:npw1_k,ii,ilm)=ylmgr1(1+ikg1:npw1_k+ikg1,ii,ilm)
           end do
         end do
       end if
     end if

!    Set up the ground-state Hamiltonian, and some parts of the 1st-order Hamiltonian

!    Compute (k+G) vectors
     nkpg=0;if(ipert>=1.and.ipert<=natom) nkpg=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute (k+q+G) vectors
     nkpg1=0;if(ipert>=1.and.ipert<=natom) nkpg1=3*dtset%nloalg(3)
     ABI_ALLOCATE(kpg1_k,(npw1_k,nkpg1))
     if (nkpg1>0) then
       call mkkpg(kg1_k,kpg1_k,kpq,nkpg1,npw1_k)
     end if

!    ===== Preparation of non-local contributions

     dimffnlk=0;if (ipert<=natom) dimffnlk=1
     ABI_ALLOCATE(ffnlk,(npw_k,dimffnlk,psps%lmnmax,dtset%ntypat))

!    Compute nonlocal form factors ffnlk at (k+G)
     if (ipert<=natom ) then
       ider=0;idir0=0
       call mkffnl(psps%dimekb,dimffnlk,psps%ekb,ffnlk,psps%ffspl,&
&       gmet,gprimd,ider,idir0,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,&
&       psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg,npw_k,dtset%ntypat,&
&       psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_dum)
     end if

!    Compute nonlocal form factors ffnl1 at (k+q+G)
     !-- Atomic displacement perturbation
     if (ipert<=natom) then
       ider=0;idir0=0
     !-- Strain perturbation
     else if (ipert==natom+3.or.ipert==natom+4) then
       if (ipert==natom+3) istr=idir
       if (ipert==natom+4) istr=idir+3
       ider=1;idir0=-istr
     else if (ipert==natom+5) then !SPr deb rfmagn
       ider=0;idir0=0
     end if
     dimffnl1=1+ider;if (ider==1.and.idir0==0) dimffnl1=dimffnl1+2*psps%useylm
     ABI_ALLOCATE(ffnl1,(npw1_k,dimffnl1,psps%lmnmax,dtset%ntypat))
     call mkffnl(psps%dimekb,dimffnl1,psps%ekb,ffnl1,psps%ffspl,gmet,gprimd,ider,idir0,&
&     psps%indlmn,kg1_k,kpg1_k,kpq,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,nkpg1,&
&     npw1_k,dtset%ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm1_k,ylmgr1_k)

!    ===== Preparation of kinetic contributions

     ABI_ALLOCATE(dkinpw,(npw_k))
     ABI_ALLOCATE(kinpw1,(npw1_k))

!    Compute the derivative of the kinetic operator vs strain in dkinpw
     if (ipert==natom+3.or.ipert==natom+4) then
       if (ipert==natom+3) istr=idir
       if (ipert==natom+4) istr=idir+3
       call kpgstr(dkinpw,dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,gprimd,istr,&
&       kg_k,kpoint,npw_k)
     end if

!    Compute (1/2) (2 Pi)**2 (k+q+G)**2:
!     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k)
     call mkkin(dtset%ecut,dtset%ecutsm,dtset%effmass_free,gmet,kg1_k,kinpw1,kpq,npw1_k,0,0)

!    ===== Load the k/k+q dependent parts of the Hamiltonian

!    Load k-dependent part in the Hamiltonian datastructure
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamkq%matblk))
     call gs_hamkq%load_k(kpt_k=kpoint,npw_k=npw_k,istwf_k=istwf_k,kg_k=kg_k,kpg_k=kpg_k,&
&     ph3d_k=ph3d,compute_ph3d=.true.,compute_gbound=.true.)
     if (size(ffnlk)>0) then
       call gs_hamkq%load_k(ffnl_k=ffnlk)
     else
       call gs_hamkq%load_k(ffnl_k=ffnl1)
     end if

!    Load k+q-dependent part in the Hamiltonian datastructure
!        Note: istwf_k is imposed to 1 for RF calculations (should use istwf_kq instead)
     call gs_hamkq%load_kprime(kpt_kp=kpq,npw_kp=npw1_k,istwf_kp=istwf_k,&
&     kinpw_kp=kinpw1,kg_kp=kg1_k,kpg_kp=kpg1_k,ffnl_kp=ffnl1,&
&     compute_gbound=.true.)
     if (qne0) then
       ABI_ALLOCATE(ph3d1,(2,npw1_k,gs_hamkq%matblk))
       call gs_hamkq%load_kprime(ph3d_kp=ph3d1,compute_ph3d=.true.)
     end if

!    Load k-dependent part in the 1st-order Hamiltonian datastructure
     call rf_hamkq%load_k(npw_k=npw_k,dkinpw_k=dkinpw)

!    Compute fixed contributions to 1st-order Fermi energy
!    and Fermi level charge density
     fe1fixed_k(:)=zero ; fe1norm_k(:)=zero

!    Note that dfpt_wfkfermi is called with kpoint, while kpt is used inside dfpt_wfkfermi
     if (nspden/=4) then
       call dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,dtfil,eig0_k,eig1_k,fe1fixed_k,&
&       fe1norm_k,gs_hamkq,ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,dtset%kptopt,mband,mband_mem,&
&       mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k,&
&       pawrhoijfermi_unsym,prtvol,rf_hamkq,rhoaug,rocceig,wtk_k)
     else
       call dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,dtfil,eig0_k,eig1_k,fe1fixed_k,&
&       fe1norm_k,gs_hamkq,ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,dtset%kptopt,mband,mband_mem,&
&       mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,npw_k,npw1_k,dtset%nspinor,nsppol,occ_k,&
&       pawrhoijfermi_unsym,prtvol,rf_hamkq,rhoaug4,rocceig,wtk_k)
     end if
!    Free temporary storage
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(kpg1_k)
     ABI_DEALLOCATE(dkinpw)
     ABI_DEALLOCATE(ffnlk)
     ABI_DEALLOCATE(ffnl1)
     ABI_DEALLOCATE(kinpw1)
     ABI_DEALLOCATE(doccde_k)
     ABI_DEALLOCATE(doccde_kq)
     ABI_DEALLOCATE(eig0_k)
     ABI_DEALLOCATE(eig0_kq)
     ABI_DEALLOCATE(occ_kq)
     ABI_DEALLOCATE(rocceig)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylm1_k)
     ABI_DEALLOCATE(ylmgr1_k)
     ABI_DEALLOCATE(ph3d)
     if (allocated(ph3d1)) then
       ABI_DEALLOCATE(ph3d1)
     end if

!    Save eigenvalues (hartree)
     eigen1 (1+bd2tot_index : 2*nband_k**2+bd2tot_index) = eig1_k(:)

!    Accumulate sum over k points for 1st-order Fermi energy components
     do iband=1,nband_k
       fe1fixed=fe1fixed+wtk_k*occ_k(iband)*fe1fixed_k(iband)
       fe1norm=fe1norm+wtk_k*occ_k(iband)*fe1norm_k(iband)
     end do

     ABI_DEALLOCATE(eig1_k)
     ABI_DEALLOCATE(occ_k)
     ABI_DEALLOCATE(fe1fixed_k)
     ABI_DEALLOCATE(fe1norm_k)

!    Keep track of total number of bands
!    (all k points so far, even for k points not treated by me)
     bdtot_index=bdtot_index+nband_k
     bd2tot_index=bd2tot_index+2*nband_k**2

!    Shift array memory
     if (mkmem/=0) then
       ibg=ibg+nband_k
       icg=icg+npw_k*dtset%nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
       ikg=ikg+npw_k
     end if
     if (mkqmem/=0) then
       ibgq=ibgq+dtset%nspinor*nband_k
       icgq=icgq+npw1_k*dtset%nspinor*proc_distrb_nband(mpi_enreg%proc_distrb,ikpt,isppol,me)
     end if
     if (mk1mem/=0) then
       ikg1=ikg1+npw1_k
     end if

!    End big k point loop
   end do

   call timab(125,2,tsec)

!  Transfer density on augmented fft grid to normal fft grid in real space
!  Also take into account the spin.
   if (nspden/=4) then
     if (psps%usepaw==0) then
       call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhorfermi,rhoaug,1)
     else
       call fftpac(isppol,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhowfr   ,rhoaug,1)
     end if
   else
     if (psps%usepaw==0) then
       do ispden=1,4
         call fftpac(ispden,mpi_enreg,nspden,cplex*n1,n2,n3,cplex*n4,n5,n6,dtset%ngfft,rhorfermi,rhoaug4(:,:,:,ispden),1)
       end do
     end if
   end if

 end do ! End loop over spins

!More memory cleaning
 call gs_hamkq%free()
 call rf_hamkq%free()
 if(psps%usepaw==1) then
   call paw_ij_free(paw_ij1fr)
   ABI_DATATYPE_DEALLOCATE(paw_ij1fr)
 end if
 if (nspden/=4) then
   ABI_DEALLOCATE(rhoaug)
 else
   ABI_DEALLOCATE(rhoaug4)
 end if
 ABI_DEALLOCATE(kg_k)
 ABI_DEALLOCATE(kg1_k)

 call timab(124,2,tsec)


!=== MPI communications ==================
 if(xmpi_paral==1)then
   call timab(129,1,tsec)

!  Identify MPI buffer size
   buffer_size=cplex*dtset%nfft*nspden+2+mbd2kpsp
   ABI_ALLOCATE(buffer1,(buffer_size))

!  Pack rhorfermi, fe1fixed, fe1norm
   indx=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
     buffer1(1:indx)=reshape(rhorfermi,(/indx/))
   else
     buffer1(1:indx)=reshape(rhowfr,(/indx/))
   end if
   buffer1(indx+1)=fe1fixed ; buffer1(indx+2)=fe1norm
   indx=indx+2 ; bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       buffer1(indx+1:indx+2*nband_k**2)=eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       indx=indx+2*nband_k**2
     end do
   end do
   if(indx<buffer_size)buffer1(indx+1:buffer_size)=zero

!  Build sum of everything
   call timab(48,1,tsec)
   call xmpi_sum(buffer1,buffer_size,spaceworld,ierr)
   call timab(48,2,tsec)

!  Unpack the final result
   indx=cplex*dtset%nfft*nspden
   if (psps%usepaw==0) then
     rhorfermi(:,:)=reshape(buffer1(1:indx),(/cplex*dtset%nfft,nspden/))
   else
     rhowfr(:,:)=reshape(buffer1(1:indx),(/cplex*dtset%nfft,nspden/))
   end if
   fe1fixed=buffer1(indx+1) ; fe1norm =buffer1(indx+2)
   indx=indx+2 ; bd2tot_index=0
   do isppol=1,nsppol
     do ikpt=1,nkpt_rbz
       nband_k=nband_rbz(ikpt+(isppol-1)*nkpt_rbz)
       eigen1(bd2tot_index+1:bd2tot_index+2*nband_k**2)=buffer1(indx+1:indx+2*nband_k**2)
       bd2tot_index=bd2tot_index+2*nband_k**2
       indx=indx+2*nband_k**2
     end do
   end do
   ABI_DEALLOCATE(buffer1)

!  Accumulate PAW occupancies
   if (psps%usepaw==1) then
     call pawrhoij_mpisum_unpacked(pawrhoijfermi_unsym,spaceworld)
   end if

   call timab(129,2,tsec)
 end if ! if kpt parallel
!=== End communications ==================

 call timab(127,1,tsec)

!Normalize the fixed part of fermie1
 invfe1norm = zero ; if (abs(fe1norm) > tol10) invfe1norm=one/fe1norm
 fe1fixed=fe1fixed*invfe1norm


 if(nspden==4) then
! FR SPr symrhg will manage correctly this rearrangement
   rhorfermi(:,2)=rhorfermi(:,2)+(rhorfermi(:,1)+rhorfermi(:,4))    !(n+mx)
   rhorfermi(:,3)=rhorfermi(:,3)+(rhorfermi(:,1)+rhorfermi(:,4))    !(n+my)
   call timab(17,2,tsec)
 end if

!Symmetrize the density
!In order to have the symrhg working in parallel on FFT coefficients, the size
!of irzzon1 and phnons1 should be set to nfftot. Therefore, nsym\=1 does not work.
!We also have the spin-up density, symmetrized, in rhorfermi(:,2).
 ABI_ALLOCATE(rhogfermi,(2,dtset%nfft))
 if (psps%usepaw==0) then
   call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,nspden,&
&   nsppol,nsym1,phnons1,rhogfermi,rhorfermi,rprimd,symaf1,symrl1,tnons1)
 else
   call symrhg(cplex,gprimd,irrzon1,mpi_enreg,dtset%nfft,dtset%nfft,dtset%ngfft,nspden,&
&   nsppol,nsym1,phnons1,rhogfermi,rhowfr,rprimd,symaf1,symrl1,tnons1)
 end if

!PAW: Build new rhoij quantities then symetrize them
!Compute and add the compensation density to rhowfr to get the total density
 if (psps%usepaw == 1) then
   if (size(nhatfermi)>0) then
     call pawmkrho(1,arg,cplex,gprimd,0,indsy1,0,mpi_enreg,&
&     my_natom,natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,&
&     pawfgrtab,-10001,pawrhoijfermi,pawrhoijfermi_unsym,pawtab,dtset%qptn,&
&     rhogfermi,rhowfr,rhorfermi,rprimd,symaf1,symrc1,dtset%typat,ucvol,&
&     dtset%usewvl,xred,pawang_sym=pawang1,pawnhat=nhatfermi)
   else
     call pawmkrho(1,arg,cplex,gprimd,0,indsy1,0,mpi_enreg,&
&     my_natom,natom,nspden,nsym1,dtset%ntypat,dtset%paral_kgb,pawang,pawfgr,&
&     pawfgrtab,-10001,pawrhoijfermi,pawrhoijfermi_unsym,pawtab,dtset%qptn,&
&     rhogfermi,rhowfr,rhorfermi,rprimd,symaf1,symrc1,dtset%typat,ucvol,&
&     dtset%usewvl,xred,pawang_sym=pawang1)
   end if
   ABI_DEALLOCATE(rhowfr)
   call pawrhoij_free_unpacked(pawrhoijfermi_unsym)
   if (paral_atom) then
     call pawrhoij_free(pawrhoijfermi_unsym)
     ABI_DATATYPE_DEALLOCATE(pawrhoijfermi_unsym)
   end if
 end if
 ABI_DEALLOCATE(rhogfermi)

!Normalize the Fermi level charge density (and associated PAW occupancies)
 rhorfermi(:,:)=invfe1norm*rhorfermi(:,:)
 if (psps%usepaw==1) then
   if (size(nhatfermi)>0) nhatfermi(:,:)=invfe1norm*nhatfermi(:,:)
   do iatom=1,my_natom
     do ispden=1,nspden
       do ilmn=1,pawrhoijfermi(iatom)%nrhoijsel
         pawrhoijfermi(iatom)%rhoijp(ilmn,ispden)=&
&         pawrhoijfermi(iatom)%rhoijp(ilmn,ispden)*invfe1norm
       end do
     end do
   end do
 end if

 call timab(127,2,tsec)
 call timab(121,2,tsec)

 DBG_EXIT('COLL')

end subroutine dfpt_rhofermi
!!***

!!****f* ABINIT/dfpt_wfkfermi
!! NAME
!! dfpt_wfkfermi
!!
!! FUNCTION
!! This routine computes the partial Fermi-level density at a given k-point,
!! and the fixed contribution to the 1st-order Fermi energy (nonlocal and kinetic)
!!
!! INPUTS
!!  cg(2,mpw*nspinor*mband_mem*mkmem*nsppol)=planewave coefficients of wavefunctions
!!  cgq(2,mcgq)=array for planewave coefficients of wavefunctions.
!!  cplex=1 if rhoaug is real, 2 if rhoaug is complex
!TODO MJV : PAW
!!  cprj(natom,nspinor*mband*mkmem*nsppol*usecprj)= wave functions at k
!!              projected with non-local projectors: cprj=<p_i|Cnk>
!!  cprjq(natom,nspinor*mband*mkqmem*nsppol*usecprj)= wave functions at k+q
!!              projected with non-local projectors: cprjq=<p_i|Cnk+q>
!!  dtfil <type(datafiles_type)>=variables related to files
!!  eig0_k(nband_k)=GS eigenvalues at k (hartree)
!!  gs_hamkq <type(gs_hamiltonian_type)>=all data for the Hamiltonian at k+q
!!  ibg=shift to be applied on the location of data in the array cprj
!!  ibgq=shift to be applied on the location of data in the array cprjq
!!  icg=shift to be applied on the location of data in the array cg
!!  icgq=shift to be applied on the location of data in the array cgq
!!  idir=direction of the current perturbation
!!  ikpt=number of the k-point
!!  ipert=type of the perturbation
!!  isppol=1 for unpolarized, 2 for spin-polarized
!!  kptopt=option for the generation of k points
!!  mband=maximum number of bands
!!  mband_mem=maximum number of bands on this cpu
!!  mcgq=second dimension of the cgq array
!!  mcprjq=second dimension of the cprjq array
!!  mkmem =number of k points treated by this node.
!!  mpi_enreg=information about MPI parallelization
!!  mpw=maximum dimensioned size of npw or wfs at k
!!  nband_k=number of bands at this k point for that spin polarization
!!  ncpgr=number of gradients stored in cprj array (cprj=<p_i|Cnk>)
!!  npw_k=number of plane waves at this k point
!!  npw1_k=number of plane waves at this k+q point
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  occ_k(nband_k)=occupation number for each band (usually 2) for each k.
!!  prtvol=control print volume and debugging output
!!  rf_hamkq <type(gs_hamiltonian_type)>=all data for the 1st-order Hamiltonian at k,q
!!  rhoaug(cplex*n4,n5,n6)= density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output)
!!  rocceig(nband_k,nband_k)= (occ_kq(m)-occ_k(n))/(eig0_kq(m)-eig0_k(n)),
!!   if this ratio has been attributed to the band n (second argument), zero otherwise
!!  wtk_k=weight assigned to the k point.
!!
!! OUTPUT
!!  eig1_k(2*nband_k**2)=first-order eigenvalues (hartree)
!!  fe1fixed_k(nband_k)=contribution to 1st-order Fermi energy
!!      from changes of occupation from all bands at this k point.
!!  fe1norm_k(nband_k)=contribution to normalization for above
!!  rhoaug(cplex*n4,n5,n6)= Fermi-level density in electrons/bohr**3,
!!   on the augmented fft grid. (cumulative, so input as well as output).
!!  ==== if (gs_hamkq%usepaw==1) ====
!!    pawrhoijfermi(natom) <type(pawrhoij_type)>= paw rhoij occupancies
!!       at Fermi level (cumulative, so input as well as output)
!!
!! PARENTS
!!      dfpt_rhofermi
!!
!! CHILDREN
!!      dfpt_accrho,dotprod_g,getgh1c,pawcprj_alloc,pawcprj_axpby,pawcprj_copy
!!      pawcprj_free,pawcprj_get,status,timab,wrtout
!!
!! SOURCE

subroutine dfpt_wfkfermi(cg,cgq,cplex,cprj,cprjq,&
&          dtfil,eig0_k,eig1_k,fe1fixed_k,fe1norm_k,gs_hamkq,&
&          ibg,ibgq,icg,icgq,idir,ikpt,ipert,isppol,&
&          kptopt,mband,mband_mem,mcgq,mcprjq,mkmem,mpi_enreg,mpw,nband_k,ncpgr,&
&          npw_k,npw1_k,nspinor,nsppol,occ_k,pawrhoijfermi,prtvol,&
&          rf_hamkq,rhoaug,rocceig,wtk_k)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: cplex,ibg,ibgq,icg,icgq,idir,ikpt
 integer,intent(in) :: ipert,isppol,kptopt,mband,mcgq,mcprjq,mkmem,mpw,ncpgr
 integer,intent(in) :: mband_mem
 integer,intent(in) :: npw1_k,nspinor,nsppol,prtvol
 integer,intent(inout) :: nband_k,npw_k
 real(dp),intent(in) :: wtk_k
 type(MPI_type),intent(in) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(gs_hamiltonian_type),intent(inout) :: gs_hamkq
 type(rf_hamiltonian_type),intent(inout) :: rf_hamkq
!arrays
 real(dp),intent(in) :: cg(2,mpw*nspinor*mband_mem*mkmem*nsppol),cgq(2,mcgq)
 real(dp),intent(in) :: eig0_k(nband_k),occ_k(nband_k),rocceig(nband_k,nband_k)
 real(dp),intent(inout) :: rhoaug(cplex*gs_hamkq%n4,gs_hamkq%n5,gs_hamkq%n6,gs_hamkq%nvloc)
 real(dp),intent(inout) :: eig1_k(2*nband_k**2)
 real(dp),intent(out) :: fe1fixed_k(nband_k)
 real(dp),intent(out) :: fe1norm_k(nband_k)
!TODO MJV : PAW
 type(pawcprj_type),intent(in) :: cprj(gs_hamkq%natom,nspinor*mband_mem*mkmem*nsppol*gs_hamkq%usecprj)
 type(pawcprj_type),intent(in) :: cprjq(gs_hamkq%natom,mcprjq)
 type(pawrhoij_type),intent(inout) :: pawrhoijfermi(gs_hamkq%natom*gs_hamkq%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: level=18
 integer :: berryopt,iband,ii,indx,iorder_cprj
 integer :: iband_me
 integer :: ipw,me,nkpt_max,optlocal,optnl,opt_accrho,opt_corr
 integer :: opt_gvnlx1,sij_opt,tim_fourwf,tim_getgh1c,usevnl
 real(dp) :: dotr,lambda,wtband
 character(len=500) :: msg
!arrays
 real(dp) :: dum_grad_berry(1,1),dum_gvnlx1(1,1),dum_gs1(1,1),tsec(2)
 real(dp),allocatable :: cwave0(:,:),cwaveq(:,:),gh1(:,:)
 type(pawcprj_type),allocatable :: cwaveprj0(:,:),cwaveprjq(:,:),cwaveprj_tmp(:,:)

! *********************************************************************

 DBG_ENTER('COLL')

!Check arguments validity
 if (ipert>gs_hamkq%natom.and.ipert/=gs_hamkq%natom+3.and.ipert/=gs_hamkq%natom+4.and.ipert/=gs_hamkq%natom+5) then !SPr rfmagn deb
   MSG_BUG('wrong ipert argument !')
 end if
 if (cplex/=1) then
   MSG_BUG('wrong cplex/=1 argument !')
 end if

!Debugging statements
 if(prtvol==-level)then
   write(msg,'(80a,a,a)') ('=',ii=1,80),ch10,'dfpt_wfkfermi : enter'
   call wrtout(std_out,msg,'PERS')
 end if
 nkpt_max=50;if(xmpi_paral==1)nkpt_max=-1

 if(prtvol>2 .or. ikpt<=nkpt_max)then
   write(msg, '(a,a,i5,2x,a,3f9.5,2x,a)' ) ch10,&
&   ' Non-SCF iterations; k pt #',ikpt,'k=',gs_hamkq%kpt_k(:),' band residuals:'
   call wrtout(std_out,msg,'PERS')
 end if

!Retrieve parallelism data
 me=mpi_enreg%me_kpt
!Initializations and allocations

 ABI_ALLOCATE(gh1,(2,npw1_k*nspinor))
 ABI_ALLOCATE(cwave0,(2,npw_k*nspinor))
 ABI_ALLOCATE(cwaveq,(2,npw1_k*nspinor))
 iorder_cprj=0 ; eig1_k(:)=zero
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(gs_hamkq%natom,nspinor))
   ABI_DATATYPE_ALLOCATE(cwaveprjq,(gs_hamkq%natom,nspinor))
   call pawcprj_alloc(cwaveprj0,1,gs_hamkq%dimcprj)
   call pawcprj_alloc(cwaveprjq,0,gs_hamkq%dimcprj)
 else
   ABI_DATATYPE_ALLOCATE(cwaveprj0,(0,0))
   ABI_DATATYPE_ALLOCATE(cwaveprjq,(0,0))
 end if
!Arguments of getgh1c routine (want only (NL+kin) frozen H(1))
 berryopt=0;usevnl=0;sij_opt=-gs_hamkq%usepaw;tim_getgh1c=3
 optlocal=0;optnl=1;opt_gvnlx1=0
 if(ipert==gs_hamkq%natom+5) optnl=0;    ! no 1st order NL in H(1), also no kin, but this will be taken into account later
!if(ipert==gs_hamkq%natom+5) optlocal=0; ! 1st order LOCAL potential present

!Arguments of the dfpt_accrho routine
 tim_fourwf=5 ; opt_accrho=1 ; opt_corr=0
!Null potentially unassigned output variables
 fe1fixed_k(:)=zero; fe1norm_k(:)=zero

!Read the npw and kg records of wf files
 call timab(139,1,tsec)

!Loop over bands
 iband_me = 0
 do iband=1,nband_k

!  Skip bands not treated by current proc
   if(mpi_enreg%proc_distrb(ikpt, iband,isppol)/=me) cycle
   iband_me = iband_me + 1

!  Select occupied bands
   if(abs(occ_k(iband))>tol8.and.abs(rocceig(iband,iband))>tol8)then

     wtband=rocceig(iband,iband)/occ_k(iband)
!    Get ground-state wavefunctions at k
     do ipw=1,npw_k*nspinor
       cwave0(1,ipw)=cg(1,ipw+(iband_me-1)*npw_k*nspinor+icg)
       cwave0(2,ipw)=cg(2,ipw+(iband_me-1)*npw_k*nspinor+icg)
     end do

     if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!      Read PAW ground state projected WF (cprj)
!TODO MJV: PAW
       call pawcprj_get(gs_hamkq%atindx1,cwaveprj0,cprj,gs_hamkq%natom,iband,ibg,ikpt,iorder_cprj,&
&       isppol,mband,mkmem,gs_hamkq%natom,1,nband_k,nspinor,nsppol,dtfil%unpaw,&
&       mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb,&
&       icpgr=idir,ncpgr=ncpgr)
     end if

!    Read ground-state wavefunctions at k+q
     indx=npw1_k*nspinor*(iband_me-1)+icgq
     cwaveq(:,1:npw_k*nspinor)=wtband*cgq(:,1+indx:npw_k*nspinor+indx)
     if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
!      Read PAW ground state projected WF (cprj)
!TODO MJV: PAW
       indx=nspinor*(iband-1)+ibgq
       call pawcprj_copy(cprjq(:,1+indx:nspinor+indx),cwaveprjq)
       call pawcprj_axpby(zero,wtband,cwaveprj_tmp,cwaveprjq)
     end if

!    Apply H^(1)-Esp.S^(1) to Psi^(0) (H(^1)=only (NL+kin) frozen part)
     lambda=eig0_k(iband)
     call getgh1c(berryopt,cwave0,cwaveprj0,gh1,dum_grad_berry,dum_gs1,gs_hamkq,dum_gvnlx1,&
&     idir,ipert,lambda,mpi_enreg,optlocal,optnl,opt_gvnlx1,rf_hamkq,sij_opt,&
&     tim_getgh1c,usevnl)
!    Compute Eig1=<Psi^(0)|H^(1)-Eps.S^(1)|Psi(0)>
     call dotprod_g(dotr,lambda,gs_hamkq%istwf_k,npw_k*nspinor,1,cwave0,gh1,mpi_enreg%me_g0, &
&     mpi_enreg%comm_spinorfft)
     indx=2*iband-1+(iband-1)*2*nband_k
     eig1_k(indx)=dotr
!    Compute the fixed contribution to the 1st-order Fermi energy
     fe1fixed_k(iband)=two*wtband*eig1_k(indx)
     fe1norm_k(iband) =two*wtband

!    Accumulate contribution to density and PAW occupation matrix

     call dfpt_accrho(cplex,cwave0,cwaveq,cwaveq,cwaveprj0,cwaveprjq,dotr,&
       gs_hamkq,iband,0,0,isppol,kptopt,mpi_enreg,gs_hamkq%natom,nband_k,ncpgr,&
       npw_k,npw1_k,nspinor,occ_k,opt_accrho,pawrhoijfermi,rhoaug,tim_fourwf,&
       opt_corr,wtk_k)
   end if ! End of non-zero occupation and rocceig

 end do ! End loop over bands

 call timab(139,2,tsec)
 call timab(130,1,tsec)

 ABI_DEALLOCATE(cwave0)
 ABI_DEALLOCATE(cwaveq)
 ABI_DEALLOCATE(gh1)
 if (gs_hamkq%usepaw==1.and.gs_hamkq%usecprj==1) then
   call pawcprj_free(cwaveprj0)
   call pawcprj_free(cwaveprjq)
 end if
 ABI_DATATYPE_DEALLOCATE(cwaveprj0)
 ABI_DATATYPE_DEALLOCATE(cwaveprjq)

!Structured debugging : if prtvol=-level, stop here.
 if(prtvol==-level)then
   write(msg,'(a,a1,a,i2,a)')' fermie3 : exit prtvol=-',level,', debugging mode => stop '
   MSG_ERROR(msg)
 end if

 call timab(130,2,tsec)

 DBG_EXIT('COLL')

end subroutine dfpt_wfkfermi
!!***

end module m_dfpt_scfcv
!!***

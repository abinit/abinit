!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_forstr
!! NAME
!!  m_forstr
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (DCA, XG, GMR, AF, AR, MB, MT)
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

module m_forstr

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_efield
 use m_errors
 use m_xmpi
 use m_fock
 use m_cgtools
 use m_xcdata
 use m_dtset
 use m_hightemp

 use defs_datatypes,     only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : xred2xcart, metric, stresssym
 use m_energies,         only : energies_type
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawcprj,          only : pawcprj_type, pawcprj_alloc, pawcprj_free, pawcprj_get, pawcprj_reorder, pawcprj_getdim
 use m_paw_dfpt,         only : pawgrnl
 use libxc_functionals,  only : libxc_functionals_is_hybrid
 use m_stress,           only : stress
 use m_forces,           only : forces
 use m_initylmg,         only : initylmg
 use m_xchybrid,         only : xchybrid_ncpp_cc
 use m_kg,               only : mkkpg
 use m_hamiltonian,      only : init_hamiltonian, gs_hamiltonian_type !,K_H_KPRIME
 use m_electronpositron, only : electronpositron_type, electronpositron_calctype
 use m_bandfft_kpt,      only : bandfft_kpt, bandfft_kpt_type, prep_bandfft_tabs, &
&                               bandfft_kpt_savetabs, bandfft_kpt_restoretabs
 use m_spacepar,         only : meanvalue_g, hartre
 use m_mkffnl,           only : mkffnl
 use m_mpinfo,           only : proc_distrb_cycle
 use m_nonlop,           only : nonlop
 use m_fock_getghc,      only : fock_getghc
 use m_prep_kgb,         only : prep_nonlop
 use m_paw_nhat,         only : pawmknhat
 use m_rhotoxc,          only : rhotoxc
 use m_dfpt_mkvxc,       only : dfpt_mkvxc, dfpt_mkvxc_noncoll
 use m_cgprj,            only : ctocprj
 use m_psolver,          only : psolver_hartree
 use m_wvl_psi,          only : wvl_nl_gradient
 use m_fft,              only : fourdp

 implicit none

 private
!!***

 public :: forstr     ! Drives the computation of forces and/or stress tensor
 public :: nres2vres  ! Convert a density residual into a potential residual
!!***

contains
!!***

!!****f* ABINIT/forstr
!! NAME
!! forstr
!!
!! FUNCTION
!! Drives the computation of forces and/or stress tensor
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  cg(2,mcg)=wavefunctions (may be read from disk instead of input)
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   |  from Etot(npw) data (at fixed geometry), used for making
!!   |  Pulay correction to stress tensor (hartree).  Should be <=0.
!!   | ecut=cut-off energy for plane wave basis sphere (Ha)
!!   | ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!   | effmass_free=effective mass for electrons (1. in common case)
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | ionmov=governs the movement of atoms (see help file)
!!   | densfor_pred=governs the mixed electronic-atomic part of the preconditioner
!!   | istwfk(nkpt)=input option parameter that describes the storage of wfs
!!   | kptns(3,nkpt)=reduced coordinates of k points in Brillouin zone
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs
!!   | mkmem=maximum number of k points in core memory
!!   | mpw = maximum number of plane waves
!!   | natom=number of atoms in cell
!!   | nband(nkpt*nsppol)=number of bands to be included in summation at each k point
!!   | nfft=(effective) number of FFT grid points (for this processor)
!!   | ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!   | nkpt=number of k points in Brillouin zone
!!   | nloalg(3)=governs the choice of the algorithm for non-local operator.
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | pawprtvol=control print volume and debugging output for PAW
!!   | prtvol=integer controlling volume of printed output
!!   | symafm(nsym)=(anti)ferromagnetic part of symmetry operations
!!   | tfkinfunc=1 if use of Thomas-Fermi kinetic functional
!!   |          =2 if use of recursion method
!!   | typat(natom)=type integer for each atom in cell
!!   | wtk(nkpt)=weights associated with various k points
!!   | nsym=number of symmetries in space group
!!  energies <type(energies_type)>=all part of total energy.
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grchempottn(3,natom)=d(E_chemical potential)/d(xred) (hartree)
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftf= -PAW ONLY- maximum size of 1D FFTs for the fine grid
!!         (mgfftf=mgfft for norm-conserving potential runs)
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc=dimension of the xccc3d array (0 or nfftf).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfftf= -PAW ONLY- number of FFT grid points for the fine grid
!!         (nfftf=nfft for norm-conserving potential runs)
!!  ngfftf(18)= -PAW ONLY- contain all needed information about 3D FFT for the fine grid
!!              (ngfftf=ngfft for norm-conserving potential runs)
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  nhat(nfftf,nspden*psps%usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  ntypat=number of types of atoms
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  occ(mband*nkpt*nsppol)=occupancies of bands at various k points
!!  optfor=1 if computation of forces is required
!!  optres=0 if the potential residual has to be used for forces corrections
!!        =1 if the density residual has to be used for forces corrections
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases
!!  ph1df(2,3*(2*mgfftf+1)*natom)=-PAW only- 1-dim structure factor phases for the fine grid
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  rhog(2,nfftf)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfftf,nspden)=array for electron density in electrons/bohr**3.
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  strsxc(6)=xc correction to stress
!!  stress_needed=1 if computation of stress tensor is required
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  ucvol=unit cell volume in bohr**3
!!  usecprj=1 if cprj datastructure is stored in memory
!!  vhartr(nfftf)=array for holding Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  ==== if (optfor==1) ====
!!   diffor=maximal absolute value of changes in the components of
!!          force between the input and the output.
!!   favg(3)=mean of the forces before correction for translational symmetry
!!   fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!     at input, previous value of forces,
!!     at output, new value.
!!     Note : unlike fred, this array has been corrected by enforcing
!!     the translational symmetry, namely that the sum of force
!!     on all atoms is zero.
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!   maxfor=maximal absolute value of the output array force.
!!   synlgr(3,natom)=symmetrized gradients of energy due to nonlocal contributions
!!  ==== if (stress_needed==1) ====
!!   strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  ===== if psps%usepaw==1
!!  pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!!  Be careful to the meaning of nfft (size of FFT grids):
!!   - In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!   - In case of PAW calculations:
!!     Two FFT grids are used; one with nfft points (coarse grid) for
!!     the computation of wave functions ; one with nfftf points
!!     (fine grid) for the computation of total density.
!!
!! PARENTS
!!      afterscfloop,setup_positron
!!
!! CHILDREN
!!      ctocprj,forces,forstrnps,initylmg,metric,nres2vres,pawcprj_alloc
!!      pawcprj_free,pawcprj_getdim,pawgrnl,stress,timab,wvl_nl_gradient
!!      xchybrid_ncpp_cc,xred2xcart
!!
!! SOURCE

subroutine forstr(atindx1,cg,cprj,diffor,dtefield,dtset,eigen,electronpositron,energies,favg,fcart,fock,&
&                 forold,fred,grchempottn,gresid,grewtn,grhf,grvdw,grxc,gsqcut,hightemp,indsym,&
&                 kg,kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,&
&                 nfftf,ngfftf,ngrvdw,nhat,nkxc,npwarr,&
&                 ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&
&                 pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,rhog,rhor,rprimd,stress_needed,&
&                 strsxc,strten,symrec,synlgr,ucvol,usecprj,vhartr,vpsp,&
&                 vxc,wvl,xccc3d,xred,ylm,ylmgr,qvpotzero)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mcg,mcprj,mgfftf,my_natom,n3xccc,nfftf,ngrvdw,nkxc,ntypat,optfor,optres
 integer,intent(in) :: stress_needed,usecprj
 real(dp),intent(in) :: gsqcut,qvpotzero,ucvol
 real(dp),intent(inout) :: diffor,maxfor
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(in) :: energies
 type(hightemp_type),pointer,intent(inout) :: hightemp
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_data),intent(inout) :: wvl
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(ntypat),ngfftf(18)
 integer,intent(in) :: npwarr(dtset%nkpt),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: grchempottn(3,dtset%natom),grewtn(3,dtset%natom),grvdw(3,ngrvdw),kxc(dtset%nfft,nkxc)
 real(dp),intent(in) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(in) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom)
 real(dp),intent(in) :: rhog(2,nfftf),rprimd(3,3),strsxc(6),vhartr(nfftf)
 real(dp),intent(in) :: vpsp(nfftf),vxc(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: forold(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw),rhor(nfftf,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(inout),target :: nvresid(nfftf,dtset%nspden)
 real(dp),intent(out) :: favg(3)
 real(dp),intent(inout) :: fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(inout) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(inout) :: grxc(3,dtset%natom),strten(6),synlgr(3,dtset%natom)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars

 integer :: comm_grid,ifft,ispden,ncpgr,occopt_,optgr,optgr2,option,optnc,optstr,optstr2,iorder_cprj,ctocprj_choice
 integer :: idir,iatom,unpaw,mcgbz
 integer,allocatable :: dimcprj(:)
 real(dp) ::dum,dum1,ucvol_
 logical :: apply_residual
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: kinstr(6),nlstr(6),tsec(2),strdum(6),gmet(3,3),gprimd(3,3),rmet(3,3)
 real(dp) :: dummy(0)
 real(dp),allocatable :: grnl(:),vlocal(:,:),vxc_hf(:,:),xcart(:,:),ylmbz(:,:),ylmgrbz(:,:,:)
 real(dp), ABI_CONTIGUOUS pointer :: resid(:,:)


! *************************************************************************

 call timab(910,1,tsec)
 call timab(911,1,tsec)

!Do nothing if nothing is required
 if (optfor==0.and.stress_needed==0) return

!Test size of FFT grids (1 grid in norm-conserving, 2 grids in PAW)
 if (dtset%usewvl==0) then
   if ((psps%usepaw==1.and.pawfgr%nfft/=nfftf).or.(psps%usepaw==0.and.dtset%nfft/=nfftf)) then
     MSG_BUG(' wrong values for nfft, nfftf !')
   end if
   if ((psps%usepaw==1.and.pawfgr%mgfft/=mgfftf).or.(psps%usepaw==0.and.dtset%mgfft/=mgfftf)) then
     MSG_BUG('wrong values for mgfft, mgfftf!')
   end if
 end if

!==========================================================================
!Here compute terms common to forces and stresses
!==========================================================================

 !output only if (optfor==1) but we have to allocate it
 ABI_ALLOCATE(grnl,(3*dtset%natom*optfor))
 grnl(:)=zero

!Compute nonlocal psp + potential Fock ACE parts of forces and stress tensor
!-involves summations over wavefunctions at all k points
 if (dtset%tfkinfunc>0.and.stress_needed==1) then
   kinstr(1:3)=-two/three*energies%e_kinetic/ucvol ; kinstr(4:6)=zero
   nlstr(1:6)=zero
 else if (dtset%usewvl==0) then
   occopt_=0 ! This means that occ are now fixed
   if(dtset%usefock==1 .and. associated(fock)) then
!     if((dtset%optstress/=0).and.(psps%usepaw==1)) then
     if((psps%usepaw==1).and.((dtset%optstress/=0).or.(dtset%optforces==2))) then
       if(dtset%optstress==0) then
         ctocprj_choice=2
         ncpgr=3
       end if
       if(dtset%optstress/=0) then
         ctocprj_choice=20*optfor+3*dtset%optstress
         ncpgr=6*dtset%optstress+3*optfor
       end if
       if (allocated(fock%fock_BZ%cwaveocc_prj)) then
         call pawcprj_free(fock%fock_BZ%cwaveocc_prj)
         ABI_DATATYPE_DEALLOCATE(fock%fock_BZ%cwaveocc_prj)
         ABI_DATATYPE_ALLOCATE(fock%fock_BZ%cwaveocc_prj,(dtset%natom,fock%fock_BZ%mcprj))
         ABI_ALLOCATE(dimcprj,(dtset%natom))
         call pawcprj_getdim(dimcprj,dtset%natom,nattyp,dtset%ntypat,dtset%typat,pawtab,'O')
         call pawcprj_alloc(fock%fock_BZ%cwaveocc_prj,ncpgr,dimcprj)
         ABI_DEALLOCATE(dimcprj)
       end if
       iatom=-1;idir=0;iorder_cprj=0;unpaw=26
       call metric(gmet,gprimd,-1,rmet,rprimd,dum)
       if (fock%fock_BZ%mkpt/=dtset%mkmem.or.(fock%fock_BZ%mpi_enreg%paral_hf ==1)) then
         ABI_ALLOCATE(ylmbz,(dtset%mpw*fock%fock_BZ%mkpt,psps%mpsang*psps%mpsang*psps%useylm))
         ABI_ALLOCATE(ylmgrbz,(dtset%mpw*fock%fock_BZ%mkpt,3,psps%mpsang*psps%mpsang*psps%useylm))
         option=1; mcgbz=dtset%mpw*fock%fock_BZ%mkptband*fock%fock_common%my_nsppol
         call initylmg(gprimd,fock%fock_BZ%kg_bz,fock%fock_BZ%kptns_bz,fock%fock_BZ%mkpt,fock%fock_BZ%mpi_enreg,&
&         psps%mpsang,dtset%mpw,fock%fock_BZ%nbandocc_bz,fock%fock_BZ%mkpt,&
&         fock%fock_BZ%npwarr,dtset%nsppol,option,rprimd,ylmbz,ylmgrbz)
         call ctocprj(fock%fock_common%atindx,fock%fock_BZ%cgocc,ctocprj_choice,fock%fock_BZ%cwaveocc_prj,gmet,gprimd,iatom,idir,&
&         iorder_cprj,fock%fock_BZ%istwfk_bz,fock%fock_BZ%kg_bz,fock%fock_BZ%kptns_bz,mcgbz,&
&         fock%fock_BZ%mcprj,dtset%mgfft,fock%fock_BZ%mkpt,fock%fock_BZ%mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,nattyp,fock%fock_BZ%nbandocc_bz,dtset%natom,dtset%ngfft,fock%fock_BZ%mkpt,&
&         dtset%nloalg,fock%fock_BZ%npwarr,dtset%nspinor,&
&         dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,unpaw,&
&         xred,ylmbz,ylmgrbz)
         ABI_DEALLOCATE(ylmbz)
         ABI_DEALLOCATE(ylmgrbz)
       else
         call ctocprj(fock%fock_common%atindx,fock%fock_BZ%cgocc,ctocprj_choice,fock%fock_BZ%cwaveocc_prj,gmet,gprimd,iatom,idir,&
&         iorder_cprj,fock%fock_BZ%istwfk_bz,fock%fock_BZ%kg_bz,fock%fock_BZ%kptns_bz,mcg,&
&         fock%fock_BZ%mcprj,dtset%mgfft,fock%fock_BZ%mkpt,mpi_enreg,psps%mpsang,&
&         dtset%mpw,dtset%natom,nattyp,fock%fock_BZ%nbandocc_bz,dtset%natom,dtset%ngfft,fock%fock_BZ%mkpt,&
&         dtset%nloalg,fock%fock_BZ%npwarr,dtset%nspinor,&
&         dtset%nsppol,dtset%ntypat,dtset%paral_kgb,ph1d,psps,rmet,dtset%typat,ucvol,unpaw,&
&         xred,ylm,ylmgr)
       end if
     end if
   end if
   call forstrnps(cg,cprj,dtset%ecut,dtset%ecutsm,dtset%effmass_free,eigen,electronpositron,fock,grnl,&
&   dtset%istwfk,kg,kinstr,nlstr,dtset%kptns,dtset%mband,mcg,mcprj,dtset%mgfft,dtset%mkmem,&
&   mpi_enreg,psps%mpsang,dtset%mpw,my_natom,dtset%natom,dtset%nband,dtset%nfft,dtset%ngfft,&
&   dtset%nkpt,dtset%nloalg,npwarr,dtset%nspden,dtset%nspinor,dtset%nsppol,dtset%nsym,ntypat,&
&   dtset%nucdipmom,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,stress_needed,symrec,dtset%typat,&
&   usecprj,dtset%usefock,dtset%use_gpu_cuda,dtset%wtk,xred,ylm,ylmgr)
 else if (optfor>0) then !WVL
   ABI_ALLOCATE(xcart,(3, dtset%natom))
   call xred2xcart(dtset%natom, rprimd, xcart, xred)
   call wvl_nl_gradient(grnl, mpi_enreg, dtset%natom, rprimd, wvl, xcart)
   ABI_DEALLOCATE(xcart)
 end if

 call timab(911,2,tsec)
 call timab(912,1,tsec)

!PAW: add gradients due to Dij derivatives to non-local term
 if (psps%usepaw==1) then
   ABI_ALLOCATE(vlocal,(nfftf,dtset%nspden))

!$OMP PARALLEL DO COLLAPSE(2)
   do ispden=1,min(dtset%nspden,2)
     do ifft=1,nfftf
       vlocal(ifft,ispden)=vhartr(ifft)+vxc(ifft,ispden)+vpsp(ifft)
     end do
   end do
   if (dtset%nspden==4) then
!$OMP PARALLEL DO COLLAPSE(2)
     do ispden=3,4
       do ifft=1,nfftf
         vlocal(ifft,ispden)=vxc(ifft,ispden)
       end do
     end do
   end if
   ucvol_=ucvol
#if defined HAVE_BIGDFT
   if (dtset%usewvl==1) ucvol_=product(wvl%den%denspot%dpbox%hgrids)*real(product(wvl%den%denspot%dpbox%ndims),dp)
#endif
   optgr=optfor;optgr2=0;optstr=stress_needed;optstr2=0
   comm_grid=mpi_enreg%comm_fft;if(dtset%usewvl==1) comm_grid=mpi_enreg%comm_wvl
   call pawgrnl(atindx1,dtset%nspden,dummy,1,dummy,grnl,gsqcut,mgfftf,my_natom,dtset%natom,&
&   nattyp,nfftf,ngfftf,nhat,nlstr,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,optstr2,&
&   pawang,pawfgrtab,pawrhoij,pawtab,ph1df,psps,k0,rprimd,symrec,dtset%typat,ucvol_,vlocal,vxc,xred,&
&   mpi_atmtab=mpi_enreg%my_atmtab, comm_atom=mpi_enreg%comm_atom,mpi_comm_grid=comm_grid)
   ABI_DEALLOCATE(vlocal)

 end if
 call timab(912,2,tsec)
 call timab(913,1,tsec)

!==========================================================================
!Here compute forces (if required)
!==========================================================================
 if (optfor==1) then
   apply_residual=(optres==1 .and. dtset%usewvl==0.and.abs(dtset%densfor_pred)>=1 .and. &
&   abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5)

!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residual before calling forces routine
   if (apply_residual) then
     ABI_ALLOCATE(resid,(nfftf,dtset%nspden))
     option=0; if (dtset%densfor_pred<0) option=1
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) optnc=2
     call nres2vres(dtset,gsqcut,psps%usepaw,kxc,mpi_enreg,my_natom,nfftf,ngfftf,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,psps%usepaw,resid,xccc3d,xred,vxc)
   else
     resid => nvresid
   end if

   call forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,forold,fred,grchempottn,gresid,grewtn,&
&   grhf,grnl,grvdw,grxc,gsqcut,indsym,maxfor,mgfftf,&
&   mpi_enreg,psps%n1xccc,n3xccc,nattyp,&
&   nfftf,ngfftf,ngrvdw,ntypat,pawrad,pawtab,ph1df,psps,rhog,&
&   rhor,rprimd,symrec,synlgr,dtset%usefock,resid,vxc,wvl%descr,wvl%den,xred,&
&   electronpositron=electronpositron)

   if (apply_residual) then
     ABI_DEALLOCATE(resid)
   end if
 end if

 call timab(913,2,tsec)
 call timab(914,1,tsec)

!==========================================================================
!Here compute stress tensor (if required)
!==========================================================================

 if (stress_needed==1.and.dtset%usewvl==0) then
!   if (dtset%usefock==1 .and. associated(fock).and.fock%fock_common%optstr.and.psps%usepaw==0) then
   if (dtset%usefock==1 .and. associated(fock)) then
     if (fock%fock_common%optstr) then
       fock%fock_common%stress(1:3)=fock%fock_common%stress(1:3)-(two*energies%e_fock-energies%e_fock0)/ucvol
       if (n3xccc>0.and.psps%usepaw==0 .and. &
&        (dtset%ixc==41.or.dtset%ixc==42.or.libxc_functionals_is_hybrid())) then
         ABI_ALLOCATE(vxc_hf,(nfftf,dtset%nspden))
!compute Vxc^GGA(rho_val)
         call xchybrid_ncpp_cc(dtset,dum,mpi_enreg,nfftf,ngfftf,n3xccc,rhor,rprimd,strdum,dum1,xccc3d,vxc=vxc_hf,optstr=1)
       end if
     end if
   end if
   call stress(atindx1,dtset%berryopt,dtefield,energies%e_localpsp,dtset%efield,&
&   energies%e_hartree,energies%e_corepsp,fock,gsqcut,hightemp,dtset%ixc,kinstr,mgfftf,&
&   mpi_enreg,psps%mqgrid_vl,psps%n1xccc,n3xccc,dtset%natom,nattyp,&
&   nfftf,ngfftf,nlstr,dtset%nspden,dtset%nsym,ntypat,psps,pawrad,pawtab,ph1df,&
&   dtset%prtvol,psps%qgrid_vl,dtset%red_efieldbar,rhog,rprimd,strten,strsxc,symrec,&
&   dtset%typat,dtset%usefock,psps%usepaw,&
&   dtset%vdw_tol,dtset%vdw_tol_3bt,dtset%vdw_xc,psps%vlspl,vxc,vxc_hf,psps%xccc1d,xccc3d,psps%xcccrc,xred,&
&   psps%ziontypat,psps%znucltypat,qvpotzero,&
&   electronpositron=electronpositron)
 end if

!Memory deallocation
 ABI_DEALLOCATE(grnl)
 if (allocated(vxc_hf)) then
   ABI_DEALLOCATE(vxc_hf)
 end if


 call timab(914,2,tsec)
 call timab(910,2,tsec)

end subroutine forstr
!!***

!!****f* ABINIT/forstrnps
!! NAME
!! forstrnps
!!
!! FUNCTION
!! Compute nonlocal pseudopotential energy contribution to forces and/or stress tensor
!! as well as kinetic energy contribution to stress tensor.
!!
!! INPUTS
!!  cg(2,mcg)=wavefunctions (may be read from disk file)
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
!!  ecut=cut-off energy for plane wave basis sphere (Ha)
!!  ecutsm=smearing energy for plane wave kinetic energy (Ha)
!!  effmass_free=effective mass for electrons (1. in common case)
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  istwfk(nkpt)=input option parameter that describes the storage of wfs
!!  kg(3,mpw*mkmem)=reduced coordinates (integers) of G vecs in basis
!!  kpt(3,nkpt)=k points in reduced coordinates
!!  mband=maximum number of bands
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfft=maximum size of 1D FFTs
!!  mkmem=number of k points treated by this node.
!!  mpi_enreg=informations about MPI parallelization
!!  mpsang= 1+maximum angular momentum for nonlocal pseudopotentials
!!  mpw= maximum number of plane waves
!!  my_natom=number of atoms treated by current processor
!!  natom=number of atoms in cell.
!!  nband(nkpt)=number of bands at each k point
!!  nfft=number of FFT grid points
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkpt=number of k points in Brillouin zone
!!  nloalg(3)=governs the choice of the algorithm for non-local operator.
!!  npwarr(nkpt)=number of planewaves in basis and boundary at each k
!!  nspden=Number of spin Density components
!!  nspinor=number of spinorial components of the wavefunctions
!!  nsppol=1 for unpolarized, 2 for spin-polarized
!!  nsym=number of elements in symmetry group
!!  ntypat=number of types of atoms
!!  nucdipmom(3,my_natom)= nuclear dipole moments
!!  occ(mband*nkpt*nsppol)=occupation numbers for each band over all k points
!!  optfor=1 if computation of forces is required
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=one-dimensional structure factor information
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=1 if computation of stress tensor is required
!!  symrec(3,3,nsym)=symmetries in reciprocal space (dimensionless)
!!  typat(natom)=type of each atom
!!  usecprj=1 if cprj datastructure has been allocated
!!  use_gpu_cuda= 0 or 1 to know if we use cuda for nonlop call
!!  wtk(nkpt)=weight associated with each k point
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  if (optfor==1)
!!   grnl(3*natom*optfor)=stores grads of nonlocal energy wrt atomic coordinates
!!  if (stress_needed==1)
!!   kinstr(6)=kinetic energy part of stress tensor (hartree/bohr^3)
!!   Store 6 unique components of symmetric 3x3 tensor in the order
!!   11, 22, 33, 32, 31, 21
!!   npsstr(6)=nonlocal pseudopotential energy part of stress tensor
!!    (hartree/bohr^3)
!!
!! PARENTS
!!      forstr
!!
!! CHILDREN
!!      bandfft_kpt_restoretabs,bandfft_kpt_savetabs,destroy_hamiltonian
!!      fock_getghc,init_hamiltonian,load_k_hamiltonian,load_spin_hamiltonian
!!      meanvalue_g,mkffnl,mkkpg,nonlop,pawcprj_alloc,pawcprj_free,pawcprj_get
!!      pawcprj_reorder,prep_bandfft_tabs,prep_nonlop,stresssym,timab,xmpi_sum
!!
!! SOURCE

subroutine forstrnps(cg,cprj,ecut,ecutsm,effmass_free,eigen,electronpositron,fock,&
&  grnl,istwfk,kg,kinstr,npsstr,kpt,mband,mcg,mcprj,mgfft,mkmem,mpi_enreg,mpsang,&
&  mpw,my_natom,natom,nband,nfft,ngfft,nkpt,nloalg,npwarr,nspden,nspinor,nsppol,nsym,&
&  ntypat,nucdipmom,occ,optfor,paw_ij,pawtab,ph1d,psps,rprimd,&
&  stress_needed,symrec,typat,usecprj,usefock,use_gpu_cuda,wtk,xred,ylm,ylmgr)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mband,mcg,mcprj,mgfft,mkmem,mpsang,mpw,my_natom,natom,nfft,nkpt
 integer,intent(in) :: nspden,nsppol,nspinor,nsym,ntypat,optfor,stress_needed
 integer,intent(in) :: usecprj,usefock,use_gpu_cuda
 real(dp),intent(in) :: ecut,ecutsm,effmass_free
 type(electronpositron_type),pointer :: electronpositron
 type(MPI_type),intent(inout) :: mpi_enreg
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: istwfk(nkpt),kg(3,mpw*mkmem),nband(nkpt*nsppol)
 integer,intent(in) :: ngfft(18),nloalg(3),npwarr(nkpt)
 integer,intent(in) :: symrec(3,3,nsym),typat(natom)
 real(dp),intent(in) :: cg(2,mcg)
 real(dp),intent(in) :: eigen(mband*nkpt*nsppol),kpt(3,nkpt),nucdipmom(3,my_natom)
 real(dp),intent(in) :: occ(mband*nkpt*nsppol),ph1d(2,3*(2*mgfft+1)*natom)
 real(dp),intent(in) :: rprimd(3,3),wtk(nkpt),xred(3,natom)
 real(dp),intent(in) :: ylm(mpw*mkmem,mpsang*mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(mpw*mkmem,3,mpsang*mpsang*psps%useylm)
 real(dp),intent(out) :: grnl(3*natom*optfor),kinstr(6),npsstr(6)
 type(pawcprj_type),intent(inout) :: cprj(natom,mcprj*usecprj)
 type(paw_ij_type),intent(in) :: paw_ij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)
 type(fock_type),pointer, intent(inout) :: fock
!Local variables-------------------------------
!scalars
 integer,parameter :: tim_rwwf=7
 integer :: bandpp,bdtot_index,choice,cpopt,dimffnl,dimffnl_str,iband,iband_cprj,iband_last,ibg,icg,ider,ider_str
 integer :: idir,idir_str,ierr,ii,ikg,ikpt,ilm,ipositron,ipw,ishift,isppol,istwf_k
 integer :: mband_cprj,me_distrb,my_ikpt,my_nspinor,nband_k,nband_cprj_k,ndat,nkpg
 integer :: nnlout,npw_k,paw_opt,signs,spaceComm
 integer :: tim_nonlop,tim_nonlop_prep,usecprj_local,use_ACE_old
 integer :: blocksize,iblock,iblocksize,ibs,nblockbd
 real(dp) :: ar,renorm_factor,dfsm,ecutsm_inv,fact_kin,fsm,htpisq,kgc1
 real(dp) :: kgc2,kgc3,kin,xx
 type(gs_hamiltonian_type) :: gs_hamk
 logical :: compute_gbound,usefock_loc
 character(len=500) :: msg
 type(fock_common_type),pointer :: fockcommon
!arrays
 integer,allocatable :: kg_k(:,:)
 real(dp) :: kpoint(3),nonlop_dum(1,1),rmet(3,3),tsec(2)
 real(dp),allocatable :: cwavef(:,:),enlout(:),ffnl_sav(:,:,:,:),ffnl_str(:,:,:,:)
 real(dp),allocatable :: ghc_dum(:,:),gprimd(:,:),kpg_k(:,:),kpg_k_sav(:,:)
 real(dp),allocatable :: kstr1(:),kstr2(:),kstr3(:),kstr4(:),kstr5(:),kstr6(:)
 real(dp),allocatable :: lambda(:),occblock(:),ph3d(:,:,:),ph3d_sav(:,:,:)
 real(dp),allocatable :: weight(:),ylm_k(:,:),ylmgr_k(:,:,:)
 real(dp),allocatable,target :: ffnl(:,:,:,:)
 type(bandfft_kpt_type),pointer :: my_bandfft_kpt => null()
 type(pawcprj_type),target,allocatable :: cwaveprj(:,:)
 type(pawcprj_type),pointer :: cwaveprj_idat(:,:)


!*************************************************************************

 call timab(920,1,tsec)
 call timab(921,1,tsec)

!Init mpicomm and me
 if(mpi_enreg%paral_kgb==1)then
   spaceComm=mpi_enreg%comm_kpt
   me_distrb=mpi_enreg%me_kpt
 else
!* In case of HF calculation
   if (mpi_enreg%paral_hf==1) then
     spaceComm=mpi_enreg%comm_kpt
     me_distrb=mpi_enreg%me_kpt
   else
     spaceComm=mpi_enreg%comm_cell
     me_distrb=mpi_enreg%me_cell
   end if
 end if

!Some constants
 ipositron=abs(electronpositron_calctype(electronpositron))
 my_nspinor=max(1,nspinor/mpi_enreg%nproc_spinor)
!Smearing of plane wave kinetic energy
 ecutsm_inv=zero;if( ecutsm>1.0d-20) ecutsm_inv=1/ecutsm
!htpisq is (1/2) (2 Pi) **2:
 htpisq=0.5_dp*(two_pi)**2

!Check that fock is present if want to use fock option
 compute_gbound=.false.
 usefock_loc = (usefock==1 .and. associated(fock))
!Arrays initializations
 grnl(:)=zero
 if (usefock_loc) then
   fockcommon =>fock%fock_common
   fockcommon%optfor=.false.
   fockcommon%optstr=.false.
   use_ACE_old=fockcommon%use_ACE
   fockcommon%use_ACE=0
   if (fockcommon%optfor) compute_gbound=.true.
 end if
 if (stress_needed==1) then
   kinstr(:)=zero;npsstr(:)=zero
   if (usefock_loc) then
     fockcommon%optstr=.TRUE.
     fockcommon%stress=zero
     compute_gbound=.true.
   end if
 end if

 usecprj_local=usecprj

 if ((usefock_loc).and.(psps%usepaw==1)) then
   usecprj_local=1
   if(optfor==1)then
     fockcommon%optfor=.true.
     if (.not.allocated(fockcommon%forces_ikpt)) then
       ABI_ALLOCATE(fockcommon%forces_ikpt,(3,natom,mband))
     end if
     if (.not.allocated(fockcommon%forces)) then
       ABI_ALLOCATE(fockcommon%forces,(3,natom))
     end if
     fockcommon%forces=zero
     compute_gbound=.true.
   end if
 end if

!Initialize Hamiltonian (k-independent terms)

 call init_hamiltonian(gs_hamk,psps,pawtab,nspinor,nsppol,nspden,natom,&
& typat,xred,nfft,mgfft,ngfft,rprimd,nloalg,usecprj=usecprj_local,&
& comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,mpi_spintab=mpi_enreg%my_isppoltab,&
& paw_ij=paw_ij,ph1d=ph1d,electronpositron=electronpositron,fock=fock,&
& nucdipmom=nucdipmom,use_gpu_cuda=use_gpu_cuda)
 rmet = MATMUL(TRANSPOSE(rprimd),rprimd)
 call timab(921,2,tsec)

!need to reorder cprj=<p_lmn|Cnk> (from unsorted to atom-sorted)
 if (psps%usepaw==1.and.usecprj_local==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx)
 end if

!Common data for "nonlop" routine
 signs=1 ; idir=0  ; ishift=0 ; tim_nonlop=4 ; tim_nonlop_prep=12
 choice=2*optfor;if (stress_needed==1) choice=10*choice+3
 if (optfor==1.and.stress_needed==1)  ishift=6
 nnlout=max(1,6*stress_needed+3*natom*optfor)
 if (psps%usepaw==0) then
   paw_opt=0 ; cpopt=-1
 else
   paw_opt=2 ; cpopt=-1+3*usecprj_local
 end if

 call timab(921,2,tsec)

!LOOP OVER SPINS
 bdtot_index=0;ibg=0;icg=0
 do isppol=1,nsppol

!  Continue to initialize the Hamiltonian (PAW DIJ coefficients)
   call gs_hamk%load_spin(isppol,with_nonlocal=.true.)
   if (usefock_loc) fockcommon%isppol=isppol

!  Loop over k points
   ikg=0
   do ikpt=1,nkpt
     if (usefock_loc) fockcommon%ikpt=ikpt
     nband_k=nband(ikpt+(isppol-1)*nkpt)
     istwf_k=istwfk(ikpt)
     npw_k=npwarr(ikpt)
     kpoint(:)=kpt(:,ikpt)

     if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,1,nband_k,isppol,me_distrb)) then
       bdtot_index=bdtot_index+nband_k
       cycle
     end if

     call timab(922,1,tsec)

!    Parallelism over FFT and/or bands: define sizes and tabs
     if (mpi_enreg%paral_kgb==1) then
       my_ikpt=mpi_enreg%my_kpttab(ikpt)
       nblockbd=nband_k/(mpi_enreg%nproc_band*mpi_enreg%bandpp)
       bandpp=mpi_enreg%bandpp
       my_bandfft_kpt => bandfft_kpt(my_ikpt)
     else
       my_ikpt=ikpt
       bandpp=mpi_enreg%bandpp
       nblockbd=nband_k/bandpp
     end if
     blocksize=nband_k/nblockbd
     mband_cprj=mband/mpi_enreg%nproc_band
     nband_cprj_k=nband_k/mpi_enreg%nproc_band

     ABI_ALLOCATE(cwavef,(2,npw_k*my_nspinor*blocksize))
     if (psps%usepaw==1.and.usecprj_local==1) then
       ABI_DATATYPE_ALLOCATE(cwaveprj,(natom,my_nspinor*bandpp))
       call pawcprj_alloc(cwaveprj,0,gs_hamk%dimcprj)
     else
       ABI_DATATYPE_ALLOCATE(cwaveprj,(0,0))
     end if

     if (stress_needed==1) then
       ABI_ALLOCATE(kstr1,(npw_k))
       ABI_ALLOCATE(kstr2,(npw_k))
       ABI_ALLOCATE(kstr3,(npw_k))
       ABI_ALLOCATE(kstr4,(npw_k))
       ABI_ALLOCATE(kstr5,(npw_k))
       ABI_ALLOCATE(kstr6,(npw_k))
     end if

     ABI_ALLOCATE(kg_k,(3,mpw))
!$OMP PARALLEL DO
     do ipw=1,npw_k
       kg_k(:,ipw)=kg(:,ipw+ikg)
     end do

     ABI_ALLOCATE(ylm_k,(npw_k,mpsang*mpsang*psps%useylm))
     if (stress_needed==1) then
       ABI_ALLOCATE(ylmgr_k,(npw_k,3,mpsang*mpsang*psps%useylm))
     else
       ABI_ALLOCATE(ylmgr_k,(0,0,0))
     end if
     if (psps%useylm==1) then
!$OMP PARALLEL DO COLLAPSE(2)
       do ilm=1,mpsang*mpsang
         do ipw=1,npw_k
           ylm_k(ipw,ilm)=ylm(ipw+ikg,ilm)
         end do
       end do
       if (stress_needed==1) then
!$OMP PARALLEL DO COLLAPSE(2)
         do ilm=1,mpsang*mpsang
           do ii=1,3
             do ipw=1,npw_k
               ylmgr_k(ipw,ii,ilm)=ylmgr(ipw+ikg,ii,ilm)
             end do
           end do
         end do
       end if
     end if

!    Prepare kinetic contribution to stress tensor (Warning : the symmetry
!    has not been broken, like in mkkin.f or kpg3.f . It should be, in order to be coherent).
     if (stress_needed==1) then
       ABI_ALLOCATE(gprimd,(3,3))
       gprimd=gs_hamk%gprimd
!$OMP PARALLEL DO PRIVATE(fact_kin,ipw,kgc1,kgc2,kgc3,kin,xx,fsm,dfsm) &
!$OMP&SHARED(ecut,ecutsm,ecutsm_inv,gs_hamk,htpisq,kg_k,kpoint,kstr1,kstr2,kstr3,kstr4,kstr5,kstr6,npw_k)
       do ipw=1,npw_k
!        Compute Cartesian coordinates of (k+G)
         kgc1=gprimd(1,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(1,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(1,3)*(kpoint(3)+kg_k(3,ipw))
         kgc2=gprimd(2,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(2,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(2,3)*(kpoint(3)+kg_k(3,ipw))
         kgc3=gprimd(3,1)*(kpoint(1)+kg_k(1,ipw))+&
&         gprimd(3,2)*(kpoint(2)+kg_k(2,ipw))+&
&         gprimd(3,3)*(kpoint(3)+kg_k(3,ipw))
         kin=htpisq* ( kgc1**2 + kgc2**2 + kgc3**2 )
         fact_kin=1.0_dp
         if(kin>ecut-ecutsm)then
           if(kin>ecut)then
             fact_kin=0.0_dp
           else
!            See the routine mkkin.f, for the smearing procedure
             xx=(ecut-kin)*ecutsm_inv
!            This kinetic cutoff smoothing function and its xx derivatives
!            were produced with Mathematica and the fortran code has been
!            numerically checked against Mathematica.
             fsm=1.0_dp/(xx**2*(3+xx*(1+xx*(-6+3*xx))))
             dfsm=-3.0_dp*(-1+xx)**2*xx*(2+5*xx)*fsm**2
!            d2fsm=6.0_dp*xx**2*(9+xx*(8+xx*(-52+xx*(-3+xx*(137+xx*&
!            &                         (-144+45*xx))))))*fsm**3
             fact_kin=fsm+kin*(-ecutsm_inv)*dfsm
           end if
         end if
         kstr1(ipw)=fact_kin*kgc1*kgc1
         kstr2(ipw)=fact_kin*kgc2*kgc2
         kstr3(ipw)=fact_kin*kgc3*kgc3
         kstr4(ipw)=fact_kin*kgc3*kgc2
         kstr5(ipw)=fact_kin*kgc3*kgc1
         kstr6(ipw)=fact_kin*kgc2*kgc1
       end do ! ipw
       ABI_DEALLOCATE(gprimd)
     end if

!    Compute (k+G) vectors
     nkpg=3*nloalg(3)
     ABI_ALLOCATE(kpg_k,(npw_k,nkpg))
     if (nkpg>0) then
       call mkkpg(kg_k,kpg_k,kpoint,nkpg,npw_k)
     end if

!    Compute nonlocal form factors ffnl at all (k+G)
     ider=0;idir=0;dimffnl=1
     if (stress_needed==1) then
       ider=1;dimffnl=2+2*psps%useylm
     end if
     ABI_ALLOCATE(ffnl,(npw_k,dimffnl,psps%lmnmax,ntypat))
     call mkffnl(psps%dimekb,dimffnl,psps%ekb,ffnl,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&     ider,idir,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&     nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
     if ((stress_needed==1).and.(usefock_loc).and.(psps%usepaw==1))then
       ider_str=1; dimffnl_str=7;idir_str=-7
       ABI_ALLOCATE(ffnl_str,(npw_k,dimffnl_str,psps%lmnmax,ntypat))
       call mkffnl(psps%dimekb,dimffnl_str,psps%ekb,ffnl_str,psps%ffspl,gs_hamk%gmet,gs_hamk%gprimd,&
&       ider_str,idir_str,psps%indlmn,kg_k,kpg_k,kpoint,psps%lmnmax,psps%lnmax,psps%mpsang,psps%mqgrid_ff,&
&       nkpg,npw_k,ntypat,psps%pspso,psps%qgrid_ff,rmet,psps%usepaw,psps%useylm,ylm_k,ylmgr_k)
     end if

!    Load k-dependent part in the Hamiltonian datastructure
!     - Compute 3D phase factors
!     - Prepare various tabs in case of band-FFT parallelism
!     - Load k-dependent quantities in the Hamiltonian
     ABI_ALLOCATE(ph3d,(2,npw_k,gs_hamk%matblk))
     call gs_hamk%load_k(kpt_k=kpoint,istwf_k=istwf_k,npw_k=npw_k,&
&     kg_k=kg_k,kpg_k=kpg_k,ffnl_k=ffnl,ph3d_k=ph3d,compute_gbound=compute_gbound,compute_ph3d=.true.)

!    Load band-FFT tabs (transposed k-dependent arrays)
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_savetabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kpg=kpg_k_sav)
       call prep_bandfft_tabs(gs_hamk,ikpt,mkmem,mpi_enreg)
       call gs_hamk%load_k(npw_fft_k=my_bandfft_kpt%ndatarecv, &
&       kg_k     =my_bandfft_kpt%kg_k_gather, &
&       kpg_k    =my_bandfft_kpt%kpg_k_gather, &
       ffnl_k   =my_bandfft_kpt%ffnl_gather, &
       ph3d_k   =my_bandfft_kpt%ph3d_gather,compute_gbound=compute_gbound)
     end if

     call timab(922,2,tsec)

!    Loop over (blocks of) bands; accumulate forces and/or stresses
!    The following is now wrong. In sequential, nblockbd=nband_k/bandpp
!    blocksize= bandpp (JB 2016/04/16)
!    Note that in sequential mode iblock=iband, nblockbd=nband_k and blocksize=1
     ABI_ALLOCATE(lambda,(blocksize))
     ABI_ALLOCATE(occblock,(blocksize))
     ABI_ALLOCATE(weight,(blocksize))
     ABI_ALLOCATE(enlout,(nnlout*blocksize))
     occblock=zero;weight=zero;enlout(:)=zero
     if (usefock_loc) then
       if (fockcommon%optstr) then
         ABI_ALLOCATE(fockcommon%stress_ikpt,(6,nband_k))
         fockcommon%stress_ikpt=zero
       end if
     end if
     if ((usefock_loc).and.(psps%usepaw==1)) then
       if (fockcommon%optfor) then
         fockcommon%forces_ikpt=zero
       end if
     end if

     do iblock=1,nblockbd

       iband=(iblock-1)*blocksize+1;iband_last=min(iband+blocksize-1,nband_k)
       iband_cprj=(iblock-1)*bandpp+1
       if(proc_distrb_cycle(mpi_enreg%proc_distrb,ikpt,iband,iband_last,isppol,me_distrb)) cycle

!      Select occupied bandsddk
       occblock(:)=occ(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
       if( abs(maxval(occblock))>=tol8 ) then
         call timab(923,1,tsec)
         weight(:)=wtk(ikpt)*occblock(:)

!        gs_hamk%ffnl_k is changed in fock_getghc, so that it is necessary to restore it when stresses are to be calculated.
         if ((stress_needed==1).and.(usefock_loc).and.(psps%usepaw==1))then
           call gs_hamk%load_k(ffnl_k=ffnl)
         end if

!        Load contribution from n,k
         cwavef(:,1:npw_k*my_nspinor*blocksize)=&
&         cg(:,1+(iblock-1)*npw_k*my_nspinor*blocksize+icg:iblock*npw_k*my_nspinor*blocksize+icg)
         if (psps%usepaw==1.and.usecprj_local==1) then
           call pawcprj_get(gs_hamk%atindx1,cwaveprj,cprj,natom,iband_cprj,ibg,ikpt,0,isppol,&
&           mband_cprj,mkmem,natom,bandpp,nband_cprj_k,my_nspinor,nsppol,0,&
&           mpicomm=mpi_enreg%comm_kpt,proc_distrb=mpi_enreg%proc_distrb)
         end if

         call timab(923,2,tsec)
         call timab(926,1,tsec)

         lambda(1:blocksize)= eigen(1+(iblock-1)*blocksize+bdtot_index:iblock*blocksize+bdtot_index)
         if (mpi_enreg%paral_kgb/=1) then
           call nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,lambda,mpi_enreg,blocksize,nnlout,&
&           paw_opt,signs,nonlop_dum,tim_nonlop,cwavef,cwavef)
         else
           call prep_nonlop(choice,cpopt,cwaveprj,enlout,gs_hamk,idir,lambda,blocksize,&
&           mpi_enreg,nnlout,paw_opt,signs,nonlop_dum,tim_nonlop_prep,cwavef,cwavef)
         end if
         if ((stress_needed==1).and.(usefock_loc).and.(psps%usepaw==1))then
           call gs_hamk%load_k(ffnl_k=ffnl_str)
         end if
         call timab(926,2,tsec)

!        Accumulate non-local contributions from n,k
         if (optfor==1) then
           do iblocksize=1,blocksize
             ibs=nnlout*(iblocksize-1)
             grnl(1:3*natom)=grnl(1:3*natom)+weight(iblocksize)*enlout(ibs+1+ishift:ibs+3*natom+ishift)
           end do
         end if
         if (stress_needed==1) then
           do iblocksize=1,blocksize
             ibs=nnlout*(iblocksize-1)
             npsstr(1:6)=npsstr(1:6) + weight(iblocksize)*enlout(ibs+1:ibs+6)
           end do
         end if

!        Accumulate stress tensor kinetic contributions
         if (stress_needed==1) then
           call timab(924,1,tsec)
           do iblocksize=1,blocksize
             call meanvalue_g(ar,kstr1,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(1)=kinstr(1)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr2,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(2)=kinstr(2)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr3,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(3)=kinstr(3)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr4,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(4)=kinstr(4)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr5,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(5)=kinstr(5)+weight(iblocksize)*ar
             call meanvalue_g(ar,kstr6,0,istwf_k,mpi_enreg,npw_k,my_nspinor,&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),&
&             cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),0)
             kinstr(6)=kinstr(6)+weight(iblocksize)*ar
           end do
           call timab(924,2,tsec)
         end if

!        Accumulate stress tensor and forces for the Fock part
         if (usefock_loc) then
           if(fockcommon%optstr.or.fockcommon%optfor) then
             if (mpi_enreg%paral_kgb==1) then
               msg='forsrtnps: Paral_kgb is not yet implemented for fock stresses'
               MSG_BUG(msg)
             end if
             ndat=mpi_enreg%bandpp
             if (gs_hamk%usepaw==0) cwaveprj_idat => cwaveprj
             ABI_ALLOCATE(ghc_dum,(0,0))
             do iblocksize=1,blocksize
               fockcommon%ieigen=(iblock-1)*blocksize+iblocksize
               fockcommon%iband=(iblock-1)*blocksize+iblocksize
               if (gs_hamk%usepaw==1) then
                 cwaveprj_idat => cwaveprj(:,(iblocksize-1)*my_nspinor+1:iblocksize*my_nspinor)
               end if
               call fock_getghc(cwavef(:,1+(iblocksize-1)*npw_k*my_nspinor:iblocksize*npw_k*my_nspinor),cwaveprj_idat,&
&               ghc_dum,gs_hamk,mpi_enreg)
               if (fockcommon%optstr) then
                 fockcommon%stress(:)=fockcommon%stress(:)+weight(iblocksize)*fockcommon%stress_ikpt(:,fockcommon%ieigen)
               end if
               if (fockcommon%optfor) then
                 fockcommon%forces(:,:)=fockcommon%forces(:,:)+weight(iblocksize)*fockcommon%forces_ikpt(:,:,fockcommon%ieigen)
               end if
             end do
             ABI_DEALLOCATE(ghc_dum)
           end if
         end if
       end if

     end do ! End of loop on block of bands

!    Restore the bandfft tabs
     if (mpi_enreg%paral_kgb==1) then
       call bandfft_kpt_restoretabs(my_bandfft_kpt,ffnl=ffnl_sav,ph3d=ph3d_sav,kpg=kpg_k_sav)
     end if

!    Increment indexes
     bdtot_index=bdtot_index+nband_k
     if (mkmem/=0) then
       ibg=ibg+my_nspinor*nband_cprj_k
       icg=icg+npw_k*my_nspinor*nband_k
       ikg=ikg+npw_k
     end if

     if (usefock_loc) then
       if (fockcommon%optstr) then
         ABI_DEALLOCATE(fockcommon%stress_ikpt)
       end if
     end if

     if (psps%usepaw==1) then
       call pawcprj_free(cwaveprj)
     end if
     ABI_DATATYPE_DEALLOCATE(cwaveprj)
     ABI_DEALLOCATE(cwavef)

     ABI_DEALLOCATE(lambda)
     ABI_DEALLOCATE(occblock)
     ABI_DEALLOCATE(weight)
     ABI_DEALLOCATE(enlout)
     ABI_DEALLOCATE(ffnl)
     ABI_DEALLOCATE(kg_k)
     ABI_DEALLOCATE(kpg_k)
     ABI_DEALLOCATE(ylm_k)
     ABI_DEALLOCATE(ylmgr_k)
     ABI_DEALLOCATE(ph3d)
     if (stress_needed==1) then
       ABI_DEALLOCATE(kstr1)
       ABI_DEALLOCATE(kstr2)
       ABI_DEALLOCATE(kstr3)
       ABI_DEALLOCATE(kstr4)
       ABI_DEALLOCATE(kstr5)
       ABI_DEALLOCATE(kstr6)
     end if
     if ((stress_needed==1).and.(usefock_loc).and.(psps%usepaw==1))then
       ABI_DEALLOCATE(ffnl_str)
     end if

   end do ! End k point loop
 end do ! End loop over spins

!Stress is equal to dE/d_strain * (1/ucvol)
 npsstr(:)=npsstr(:)/gs_hamk%ucvol

!Parallel case: accumulate (n,k) contributions
 if (xmpi_paral==1) then
!  Forces
   if (optfor==1) then
     call timab(65,1,tsec)
     call xmpi_sum(grnl,spaceComm,ierr)
     call timab(65,2,tsec)
     if ((usefock_loc).and.(psps%usepaw==1)) then
       call xmpi_sum(fockcommon%forces,spaceComm,ierr)
     end if
   end if
!  Stresses
   if (stress_needed==1) then
     call timab(65,1,tsec)
     call xmpi_sum(kinstr,spaceComm,ierr)
     call xmpi_sum(npsstr,spaceComm,ierr)
     if ((usefock_loc).and.(fockcommon%optstr)) then
       call xmpi_sum(fockcommon%stress,spaceComm,ierr)
     end if
     call timab(65,2,tsec)
   end if
 end if

 call timab(925,1,tsec)

!Do final normalizations and symmetrizations of stress tensor contributions
 if (stress_needed==1) then
   renorm_factor=-(two_pi**2)/effmass_free/gs_hamk%ucvol
   kinstr(:)=kinstr(:)*renorm_factor
   if (nsym>1) then
     call stresssym(gs_hamk%gprimd,nsym,kinstr,symrec)
     call stresssym(gs_hamk%gprimd,nsym,npsstr,symrec)
     if ((usefock_loc).and.(fockcommon%optstr)) then
       call stresssym(gs_hamk%gprimd,nsym,fockcommon%stress,symrec)
     end if
   end if
 end if

!Need to reorder cprj=<p_lmn|Cnk> (from atom-sorted to unsorted)
 if (psps%usepaw==1.and.usecprj_local==1) then
   call pawcprj_reorder(cprj,gs_hamk%atindx1)
 end if

!Deallocate temporary space
 call gs_hamk%free()
 if (usefock_loc) then
   fockcommon%use_ACE=use_ACE_old
 end if
 call timab(925,2,tsec)
 call timab(920,2,tsec)

end subroutine forstrnps
!!***

!!****f* ABINIT/nres2vres
!!
!! NAME
!! nres2vres
!!
!! FUNCTION
!! Convert a density residual into a potential residual
!! using a first order formula:
!!     V^res(r)=dV/dn.n^res(r)
!!             =V_hartree(n^res)(r) + Kxc.n^res(r)
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom= number of atoms in cell
!!  | nspden=number of spin-density components
!!  | ntypat=number of atom types
!!  | typat(natom)=type (integer) for each atom
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! izero=if 1, unbalanced components of Vhartree(g) are set to zero
!! kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!! mpi_enreg=information about MPI parallelization
!! my_natom=number of atoms treated by current processor
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT
!! nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!! nkxc=second dimension of the array kxc, see rhotoxc.F90 for a description
!! nresid(nfft,nspden)= the input density residual
!! n3xccc=dimension of the xccc3d array (0 or nfft).
!! optnc=option for non-collinear magnetism (nspden=4):
!!       1: the whole 2x2 Vres matrix is computed
!!       2: only Vres^{11} and Vres^{22} are computed
!! optxc=0 if LDA part of XC kernel has only to be taken into account (even for GGA)
!!       1 if XC kernel has to be fully taken into
!!      -1 if XC kernel does not have to be taken into account
!! pawang <type(pawang_type)>=paw angular mesh and related data
!! pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!! pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhor(nfft,nspden)=electron density in real space
!!                   (used only if Kxc was not computed before)
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! xccc3d(n3xccc)=3D core electron density for XC core correction (bohr^-3)
!! xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! === optional inputs ===
!! vxc(cplex*nfft,nspden)=XC GS potential
!!
!! OUTPUT
!! vresid(nfft,nspden)= the output potential residual
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      dfpt_mkvxc,dfpt_mkvxc_noncoll,fourdp,hartre,metric,pawmknhat
!!      psolver_hartree,rhotoxc,xcdata_init
!!
!! SOURCE

subroutine nres2vres(dtset,gsqcut,izero,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&
&                 nkxc,nresid,n3xccc,optnc,optxc,pawang,pawfgrtab,pawrhoij,pawtab,&
&                 rhor,rprimd,usepaw,vresid,xccc3d,xred,vxc)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: izero,my_natom,n3xccc,nfft,nkxc,optnc,optxc,usepaw
 real(dp),intent(in) :: gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(pawang_type),intent(in) :: pawang
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: kxc(nfft,nkxc),nresid(nfft,dtset%nspden)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3),xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(out) :: vresid(nfft,dtset%nspden)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*usepaw)
 type(pawrhoij_type),intent(in) :: pawrhoij(my_natom*usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*usepaw)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden) !FR TODO:cplex?

!Local variables-------------------------------
!scalars
 integer :: cplex,ider,idir,ipert,ispden,nhatgrdim,nkxc_cur,option,me,nproc,comm,usexcnhat
 logical :: has_nkxc_gga,non_magnetic_xc
 real(dp) :: dum,energy,m_norm_min,ucvol,vxcavg
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 integer :: nk3xc
 real(dp) :: dummy6(6),gmet(3,3),gprimd(3,3),qq(3),rmet(3,3)
 real(dp),allocatable :: dummy(:),kxc_cur(:,:),nhatgr(:,:,:)
 real(dp),allocatable :: nresg(:,:),rhor0(:,:),vhres(:)

! *************************************************************************

!Compatibility tests:
 has_nkxc_gga=(nkxc==7.or.nkxc==19)

 if(optxc<-1.or.optxc>1)then
   write(message,'(a,i0)')' Wrong value for optxc ',optxc
   MSG_BUG(message)
 end if

 if((optnc/=1.and.optnc/=2).or.(dtset%nspden/=4.and.optnc/=1))then
   write(message,'(a,i0)')' Wrong value for optnc ',optnc
   MSG_BUG(message)
 end if

 if(dtset%icoulomb==1.and.optxc/=-1)then
   write(message,'(a)')' This routine is not compatible with icoulomb==1 and optxc/=-1 !'
   MSG_BUG(message)
 end if

 if(dtset%nspden==4.and.dtset%xclevel==2.and.optxc==1.and.(.not.has_nkxc_gga))then
   MSG_ERROR(' Wrong values for optxc and nkxc !')
 end if

 qq=zero
 nkxc_cur=0
 m_norm_min=EPSILON(0.0_dp)**2
 usexcnhat=0;if (usepaw==1) usexcnhat=maxval(pawtab(1:dtset%ntypat)%usexcnhat)
 non_magnetic_xc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)
 if (dtset%xclevel==1.or.optxc==0) nkxc_cur= 2*min(dtset%nspden,2)-1 ! LDA: nkxc=1,3
 if (dtset%xclevel==2.and.optxc==1)nkxc_cur=12*min(dtset%nspden,2)-5 ! GGA: nkxc=7,19
 ABI_ALLOCATE(vhres,(nfft))

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Compute density residual in reciprocal space
 if (dtset%icoulomb==0) then
   ABI_ALLOCATE(nresg,(2,nfft))
   ABI_ALLOCATE(dummy,(nfft))
   dummy(:)=nresid(:,1)
   call fourdp(1,nresg,dummy,-1,mpi_enreg,nfft,1,ngfft,0)
   ABI_DEALLOCATE(dummy)
 end if

!For GGA, has to recompute gradients of nhat
 nhatgrdim=0
 if ((nkxc==nkxc_cur.and.has_nkxc_gga).or.(optxc==-1.and.has_nkxc_gga).or.&
& (optxc/=-1.and.nkxc/=nkxc_cur)) then
   if (usepaw==1.and.dtset%xclevel==2.and.usexcnhat>0.and.dtset%pawnhatxc>0) then
     nhatgrdim=1
     ABI_ALLOCATE(nhatgr,(nfft,dtset%nspden,3))
     ider=1;cplex=1;ipert=0;idir=0
     call pawmknhat(dum,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&     nfft,ngfft,nhatgrdim,dtset%nspden,dtset%ntypat,pawang,pawfgrtab,&
&     nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qq,rprimd,ucvol,dtset%usewvl,xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&     distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
   else
     ABI_ALLOCATE(nhatgr,(0,0,0))
   end if
 else
   ABI_ALLOCATE(nhatgr,(0,0,0))
 end if

 ABI_ALLOCATE(dummy,(0))
!First case: Kxc has already been computed
!-----------------------------------------
 if (nkxc==nkxc_cur.or.optxc==-1) then

!  Compute VH(n^res)(r)
   if (dtset%icoulomb == 0) then
     call hartre(1,gsqcut,izero,mpi_enreg,nfft,ngfft,nresg,rprimd,vhres)
   else
     comm=mpi_enreg%comm_cell
     nproc=xmpi_comm_size(comm)
     me=xmpi_comm_rank(comm)
     call psolver_hartree(energy, (/ rprimd(1,1) / dtset%ngfft(1), &
&     rprimd(2,2) / dtset%ngfft(2), rprimd(3,3) / dtset%ngfft(3) /), dtset%icoulomb, &
&     me, comm, dtset%nfft, dtset%ngfft(1:3), nproc, dtset%nscforder, dtset%nspden, &
&     nresid(:,1), vhres, dtset%usewvl)
   end if

!  Compute Kxc(r).n^res(r)
   if (optxc/=-1) then

!    Collinear magnetism or non-polarized
     if (dtset%nspden/=4) then
       !Note: imposing usexcnhat=1 avoid nhat to be substracted
       call dfpt_mkvxc(1,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&       nkxc,non_magnetic_xc,dtset%nspden,0,2,qq,nresid,rprimd,1,vresid,dummy)
     else
!FR    call routine for Non-collinear magnetism
       ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
       rhor0(:,:)=rhor(:,:)-nresid(:,:)
       !Note: imposing usexcnhat=1 avoid nhat to be substracted
       call dfpt_mkvxc_noncoll(1,dtset%ixc,kxc,mpi_enreg,nfft,ngfft,nhat,usepaw,nhat,usepaw,nhatgr,nhatgrdim,&
&       nkxc,non_magnetic_xc,dtset%nspden,0,2,2,qq,rhor0,nresid,rprimd,1,vxc,vresid,xccc3d)
       ABI_DEALLOCATE(rhor0)
     end if

   else
     vresid=zero
   end if

 end if

!2nd case: Kxc has to be computed
!--------------------------------
 if (nkxc/=nkxc_cur.and.optxc/=-1) then

!  Has to use the "initial" density to compute Kxc
   ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
   rhor0(:,:)=rhor(:,:)-nresid(:,:)

!  Compute VH(n^res) and XC kernel (Kxc) together
   ABI_ALLOCATE(kxc_cur,(nfft,nkxc_cur))

   option=2;if (dtset%xclevel==2.and.optxc==0) option=12

   call hartre(1,gsqcut,izero,mpi_enreg,nfft,ngfft,nresg,rprimd,vhres)
   call xcdata_init(xcdata,dtset=dtset)

!  To be adjusted for the call to rhotoxc
   nk3xc=1
   call rhotoxc(energy,kxc_cur,mpi_enreg,nfft,ngfft,&
&   nhat,usepaw,nhatgr,nhatgrdim,nkxc_cur,nk3xc,non_magnetic_xc,n3xccc,option,&
&   rhor0,rprimd,dummy6,usexcnhat,vresid,vxcavg,xccc3d,xcdata,vhartr=vhres)  !vresid=work space
   if (dtset%nspden/=4)  then
     ABI_DEALLOCATE(rhor0)
   end if

!  Compute Kxc(r).n^res(r)

   if (dtset%nspden/=4) then
!    Collinear magnetism or non-polarized
     call dfpt_mkvxc(1,dtset%ixc,kxc_cur,mpi_enreg,nfft,ngfft,nhat,usepaw,nhatgr,nhatgrdim,&
&     nkxc_cur,non_magnetic_xc,dtset%nspden,0,2,qq,nresid,rprimd,1,vresid,dummy)
   else
!    Non-collinear magnetism
     ABI_ALLOCATE(rhor0,(nfft,dtset%nspden))
     rhor0(:,:)=rhor(:,:)-nresid(:,:)
     call dfpt_mkvxc_noncoll(1,dtset%ixc,kxc_cur,mpi_enreg,nfft,ngfft,nhat,usepaw,nhat,usepaw,nhatgr,nhatgrdim,&
&     nkxc,non_magnetic_xc,dtset%nspden,0,2,2,qq,rhor0,nresid,rprimd,1,vxc,vresid,xccc3d)
     ABI_DEALLOCATE(rhor0)
   end if

   ABI_DEALLOCATE(kxc_cur)
 end if

 !if (nhatgrdim>0)  then
 ABI_DEALLOCATE(nhatgr)
 !end if

!Assemble potential residual: V^res(r)=VH(n^res)(r) + Kxc(r).n^res(r)
!--------------------------------------------------------------------
 do ispden=1,dtset%nspden/optnc
   vresid(:,ispden)=vresid(:,ispden)+vhres(:)
 end do

 if (dtset%icoulomb==0)  then
   ABI_DEALLOCATE(nresg)
 end if
 ABI_DEALLOCATE(vhres)
 ABI_DEALLOCATE(dummy)

end subroutine nres2vres
!!***

end module m_forstr
!!***

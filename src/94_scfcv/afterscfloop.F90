!{\src2tex{textfont=tt}}
!!****f* ABINIT/afterscfloop
!! NAME
!! afterscfloop
!!
!! FUNCTION
!! Perform all calculations needed after the SCF loop, independent of the
!! call to scfcv (with or without atomic displacements), and exclusive
!! of print or write purposes, or deallocations.
!!
!! COPYRIGHT
!! Copyright (C) 2005-2017 ABINIT group (XG)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! INPUTS
!!  atindx(natom)=index table for atoms (see gstate.f)
!!  atindx1(natom)=index table for atoms, inverse of atindx (see gstate.f)
!!  cg(2,mcg)=wavefunctions (may be read from disk instead of input)
!!  cprj(natom,mcprj*usecprj)=<p_lmn|Cnk> coefficients for each WF |Cnk> and each NL proj |p_lmn>
!!  cpus= cpu time limit in seconds
!!  deltae=change in energy between the previous and present SCF cycle
!!  dtfil <type(datafiles_type)>=variables related to files
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | mband=maximum number of bands
!!   | mgfft=maximum size of 1D FFTs (see NOTES at beginning of scfcv)
!!   | mkmem=maximum number of k points in core memory
!!   | mpw = maximum number of plane waves
!!   | natom=number of atoms in cell
!!   | nfft=(effective) number of FFT grid points (for this processor) (see NOTES at beginning of scfcv)
!!   | nkpt=number of k points in Brillouin zone
!!   | nspden=number of spin-density components
!!   | nsppol=1 for unpolarized, 2 for spin-polarized
!!   | nsym=number of symmetries in space group
!!  eigen(mband*nkpt*nsppol)=array for holding eigenvalues (hartree)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D2 dispersion (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*dtset%ecut/(2._dp*(Pi**2)
!!  hdr <type(hdr_type)>=the header of wf, den and pot files
!!  indsym(4,nsym,natom)=index showing transformation of atom labels
!!                       under symmetry operations (computed in symatm)
!!  irrzon(nfft**(1-1/nsym),2,(nspden/nsppol)-3*(nspden/4))=irreducible zone data
!!  istep=number of the SCF iteration
!!  kg(3,mpw*mkmem)=reduced (integer) coordinates of G vecs in basis sphere
!!  kxc(nfftf,nkxc)=XC kernel
!!  mcg=size of wave-functions array (cg) =mpw*nspinor*mband*mkmem*nsppol
!!  mcprj=size of projected wave-functions array (cprj) =nspinor*mband*mkmem*nsppol
!!  mgfftf= - PAW only - maximum size of 1D FFTs for the "fine" grid (see NOTES at beginning of scfcv)
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  moved_atm_inside: if==1, the atoms are allowed to move.
!!  mpi_enreg=informations about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(dtset%ntypat)=number of atoms of each type
!!  nfftf= - PAW only - number of FFT grid points for the "fine" grid (see NOTES at beginning of scfcv)
!!  ngfft(18)=contain all needed information about 3D FFT (see NOTES at beginning of scfcv)
!!  ngfftf(18)= - PAW only - contain all needed information about 3D FFT  for the "fine" grid
!!  nhat(nfftf,nspden*psps%usepaw)= -PAW only- compensation density
!!  nkxc=dimension of kxc
!!  npwarr(nkpt)=number of planewaves in basis and on boundary for each k
!!  nvresid(nfftf,nspden)=array for the residual of the density/potential
!!  occ(mband*nkpt*nsppol)=occupancies of bands at various k points
!!  optres=0: the potential residual has been computed in scfcv
!!         1: the density residual has been computed in scfcv
!!  paw_an(my_natom*usepaw) <type(paw_an_type)>=paw arrays given on angular mesh
!!  paw_ij(my_natom*usepaw) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgr <type(pawfgr_type)>=fine grid parameters and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrhoij(my_natom*usepaw) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!  pawtab(dtset%ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  pel(3)=reduced coordinates of the electronic polarization (a. u.)
!!  pel_cg(3) = reduced coordinates of the electronic polarization (a. u.)
!!             computed in the SCF loop
!!  ph1df(2,3*(2*mgfftf+1)*natom)= - PAW only - 1-dim structure factor phases for the "fine" grid
!!      (see NOTES at beginning of scfcv)
!!  phnons(2,nfft**(1-1/nsym),(nspden/nsppol)-3*(nspden/4))=nonsymmorphic translation phases
!!  pion(3)=reduced coordinates of the ionic polarization (a. u.)
!!  prtfor=1 only if forces have to be printed (0 otherwise)
!!  prtxml=1 if XML file has to be output
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!   | mpsang= 1+maximum angular momentum
!!  pwind(pwind_alloc,2,3) = array used to compute
!!           the overlap matrix smat between k-points (see initberry.f)
!!  pwind_alloc = first dimension of pwind
!!  pwnsfac(2,pwind_alloc) = phase factors for non-symmorphic translations
!!  res2=density/potential residual (squared)
!!  resid(mband*nkpt*nsppol)=residuals for each band over all k points and spins
!!  residm=maximum value from resid array (except for nbdbuf highest bands)
!!  rhog(2,nfftf)=Fourier transform of total electron density (including compensation density in PAW)
!!  rhor(nfftf,nspden)=total electron density (including compensation density in PAW)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  stress_needed=1 if stresses are needed, 0 otherwise
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  tollist(12)=list of tolerances
!!  usecprj=1 if cprj datastructure has been allocated
!!  vhartr(nfftf)=Hartree potential
!!  vpsp(nfftf)=array for holding local psp
!!  vxc(nfftf,nspden)=exchange-correlation potential (hartree) in real space
!!  vxcavg=vxc average
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  ylm(mpw*mkmem,mpsang*mpsang*useylm)= real spherical harmonics for each G and k point
!!  ylmgr(mpw*mkmem,3,mpsang*mpsang*useylm)= gradients of real spherical harmonics
!!
!! OUTPUT
!!  conv_retcode= Non-zero if convergence is not achieved.
!!  elfr(nfftf,nspden)=electron localization function
!!  grhor(nfft,nspden,3)= gradient of electron density in electrons/bohr**4, real space
!!  lrhor(nfft,nspden)= Laplacian of electron density in electrons/bohr**5, real space
!!  results_gs <type(results_gs_type)>=results (energy and its components,
!!   forces and its components, the stress tensor) of a ground-state computation
!!   (should be made a pure output quantity)
!!  taug(2,nfftf)=Fourier transform of total kinetic energy density
!!  taur(nfftf,nspden)=total kinetic energy density in real space
!!  ==== if forces are required ====
!!   diffor=maximal absolute value of changes in the components of
!!          force between the input and the output.
!!   favg(3)=mean of the forces before correction for translational symmetry
!!   fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!     at input, previous value of forces,
!!     at output, new value.
!!     Note : unlike fred, this array has been corrected by enforcing
!!     the translational symmetry, namely that the sum of force
!!     on all atoms is zero.
!!   fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!   gresid(3,natom)=forces due to the residual of the potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!   maxfor=maximal absolute value of the output array force.
!!   synlgr(3,natom)=symmetrized gradients of energy due to nonlocal contributions
!!  ==== if stress tensor is required ====
!!   strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!
!! SIDE EFFECTS
!! computed_forces=1 if forces have been computed, 0 otherwise
!! dtefield <type(efield_type)> = variables related to Berry phase
!!       and electric field calculations (see initberry.f).
!!       In case dtset%berryopt = 4/6/7/14/16/17, the overlap matrices computed
!!       in this routine are stored in dtefield%smat in order
!!       to be used in the electric field calculation.
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation
!!  energies <type(energies_type)>=all part of total energy.
!!   | entropy(IN)=entropy due to the occupation number smearing (if metal)
!!   | e_localpsp(IN)=local psp energy (hartree)
!!   | e_eigenvalues(IN)=Sum of the eigenvalues - Band energy (Hartree)
!!   | e_chempot(IN)=energy from spatially varying chemical potential (hartree)
!!   | e_ewald(IN)=Ewald energy (hartree)
!!   | e_vdw_dftd(IN)=VdW DFT-D energy
!!   | e_hartree(IN)=Hartree part of total energy (hartree units)
!!   | e_corepsp(IN)=psp core-core energy
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(IN)=nonlocal pseudopotential part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely    !!HONG
!!   |                  on the electric field:
!!   |                 enefield =  -ebar_i p_i - Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j for fixed E/ebar
!!   |                          =  Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  for fixed D/d
!! etotal=total energy, might be correct by improved polarization computation
!! forold(3,natom)=old forces
!! xred(3,natom)=reduced dimensionless atomic coordinates
!! ===== if dtset%densfor_pred==3 .and. moved_atm_inside==1 =====
!!   ph1d(2,3*(2*mgfft+1)*natom)=1-dim structure factor phases (coarse grid)
!!   ph1df(2,3*(2*mgfftf+1)*natom)=1-dim structure factor phases (fine PAW grid)
!!  wvl <type(wvl_data)>=all wavelets data.
!!
!! NOTES
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      denspot_free_history,eigensystem_info,elpolariz,energies_copy
!!      exchange_electronpositron,forces,forstr,getph,hdr_update
!!      kswfn_free_scf_data,last_orthon,metric,mkrho,nhatgrid,prtposcar
!!      prtrhomxmn,scprqt,setnoccmmp,spin_current,timab,total_energies
!!      write_energies,wrtout,wvl_eigen_abi2big,wvl_mkrho,wvl_nl_gradient
!!      wvl_occ_abi2big,wvl_psitohpsi,wvl_tail_corrections,wvl_vtrial_abi2big
!!      xcden,xmpi_sum,xred2xcart
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine afterscfloop(atindx,atindx1,cg,computed_forces,cprj,cpus,&
& deltae,diffor,dtefield,dtfil,dtset,eigen,electronpositron,elfr,&
& energies,etotal,favg,fcart,fock,forold,fred,gresid,grewtn,grhf,grhor,grvdw,&
& grxc,gsqcut,hdr,indsym,irrzon,istep,kg,kxc,lrhor,maxfor,mcg,mcprj,mgfftf,&
& moved_atm_inside,mpi_enreg,my_natom,n3xccc,nattyp,nfftf,ngfft,ngfftf,ngrvdw,nhat,&
& nkxc,npwarr,nvresid,occ,optres,paw_an,paw_ij,pawang,pawfgr,&
& pawfgrtab,pawrhoij,pawtab,pel,pel_cg,ph1d,ph1df,phnons,pion,prtfor,prtxml,&
& psps,pwind,pwind_alloc,pwnsfac,res2,resid,residm,results_gs,&
& rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,taug,&
& taur,tollist,usecprj,vhartr,vpsp,vtrial,vxc,vxcavg,wvl,&
& xccc3d,xred,ylm,ylmgr,qvpotzero,conv_retcode)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_parameters
 use defs_wvltypes
 use m_energies
 use m_errors
 use m_profiling_abi
 use m_efield
 use m_ab7_mixing
 use m_hdr

 use m_xmpi,             only : xmpi_sum, xmpi_comm_rank,xmpi_comm_size
 use m_results_gs ,      only : results_gs_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype,exchange_electronpositron
 use m_dtset,            only : dtset_copy, dtset_free
 use m_abi2big,          only : wvl_eigen_abi2big,wvl_occ_abi2big,wvl_vtrial_abi2big
 use m_paw_dmft,         only : paw_dmft_type
 use m_pawang,           only : pawang_type
 use m_pawtab,           only : pawtab_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_paw_an,           only : paw_an_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawcprj,          only : pawcprj_type
 use m_pawfgr,           only : pawfgr_type
 use m_fock,             only : fock_type

#ifdef HAVE_BIGDFT
 use BigDFT_API, only : last_orthon, &
      & kswfn_free_scf_data, denspot_free_history,&
      & write_energies, total_energies,&
      & eigensystem_info
#endif

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'afterscfloop'
 use interfaces_14_hidewrite
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_56_recipspace
 use interfaces_56_xc
 use interfaces_62_wvl_wfs
 use interfaces_65_paw
 use interfaces_67_common
 use interfaces_94_scfcv, except_this_one => afterscfloop
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: istep,mcg,mcprj,mgfftf,moved_atm_inside,my_natom,n3xccc,nfftf,ngrvdw,nkxc
 integer,intent(in) :: optres,prtfor,prtxml,pwind_alloc,stress_needed,usecprj
 integer,intent(inout) :: computed_forces
 real(dp),intent(in) :: cpus,deltae,gsqcut,res2,residm
 real(dp),intent(in) :: qvpotzero
 real(dp),intent(inout) :: diffor,etotal,maxfor,vxcavg
 type(MPI_type),intent(inout) :: mpi_enreg
 type(datafiles_type),intent(in) :: dtfil
 type(dataset_type),intent(inout) :: dtset
 type(efield_type),intent(inout) :: dtefield
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(hdr_type),intent(inout) :: hdr
 type(pawang_type),intent(in) :: pawang
 type(pawfgr_type),intent(in) :: pawfgr
 type(pseudopotential_type),intent(in) :: psps
 type(results_gs_type),intent(inout) :: results_gs
 type(wvl_data),intent(inout) :: wvl
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx(dtset%natom),atindx1(dtset%natom)
 integer,intent(in) :: indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: irrzon(dtset%nfft**(1-1/dtset%nsym),2,(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 integer,intent(in) :: kg(3,dtset%mpw*dtset%mkmem),nattyp(dtset%ntypat)
 integer,intent(in) :: ngfft(18),ngfftf(18),npwarr(dtset%nkpt)
 integer,intent(in) :: pwind(pwind_alloc,2,3),symrec(3,3,dtset%nsym)
 integer,intent(out) :: conv_retcode
 real(dp),intent(in) :: grewtn(3,dtset%natom),grvdw(3,ngrvdw)
 real(dp),intent(in) :: phnons(2,dtset%nfft**(1-1/dtset%nsym),(dtset%nspden/dtset%nsppol)-3*(dtset%nspden/4))
 real(dp),intent(in) :: pwnsfac(2,pwind_alloc)
 real(dp),intent(in) :: resid(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(in) :: rprimd(3,3),tollist(12),vpsp(nfftf)
 real(dp),intent(inout) :: vtrial(nfftf,dtset%nspden)
 real(dp),intent(in) :: ylm(dtset%mpw*dtset%mkmem,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(in) :: ylmgr(dtset%mpw*dtset%mkmem,3,psps%mpsang*psps%mpsang*psps%useylm)
 real(dp),intent(inout) :: cg(2,mcg)
 real(dp),intent(inout) :: eigen(dtset%mband*dtset%nkpt*dtset%nsppol)
 real(dp),intent(inout) :: forold(3,dtset%natom)
 real(dp),intent(inout) :: nhat(nfftf,dtset%nspden*psps%usepaw)
 real(dp),intent(inout) :: nvresid(nfftf,dtset%nspden),pel(3)
 real(dp),intent(inout) :: occ(dtset%mband*dtset%nkpt*dtset%nsppol),pel_cg(3)
 real(dp),intent(inout) :: ph1d(2,3*(2*dtset%mgfft+1)*dtset%natom)
 real(dp),intent(inout) :: ph1df(2,3*(2*mgfftf+1)*dtset%natom),pion(3)
 real(dp),intent(inout) :: rhog(2,nfftf),rhor(nfftf,dtset%nspden),strsxc(6)
 real(dp),intent(inout) :: vhartr(nfftf),vxc(nfftf,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc),xred(3,dtset%natom)
 real(dp),intent(inout) :: favg(3),fcart(3,dtset%natom),fred(3,dtset%natom)
 real(dp),intent(inout) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(inout) :: grxc(3,dtset%natom),kxc(nfftf,nkxc),strten(6)
 real(dp),intent(inout) :: synlgr(3,dtset%natom)
 real(dp),pointer :: elfr(:,:),grhor(:,:,:),lrhor(:,:),taug(:,:),taur(:,:)
 type(pawcprj_type),intent(inout) :: cprj(dtset%natom,mcprj*usecprj)
 type(paw_an_type),intent(inout) :: paw_an(my_natom*psps%usepaw)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom*psps%usepaw)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(dtset%ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer,parameter :: response=0
 integer :: bantot,bufsz,choice,cplex,ierr,ifft,igrad,ishift,ispden,ngrad
 integer :: optcut,optfor,optgr0,optgr1,optgr2,optrad,quit
 integer :: spaceComm_fft,tim_mkrho,me,nproc,comm
 logical :: test_gylmgr,test_nfgd,test_rfgd
 logical :: wvlbigdft=.false.
 real(dp) :: dtaur,dtaurzero,c_fermi,ucvol
 character(len=500) :: message
 type(paw_dmft_type) :: paw_dmft
#if defined HAVE_BIGDFT
 logical :: do_last_ortho
 real(dp) :: eproj,raux
#endif
!arrays
 real(dp) :: gmet(3,3),gprimd(3,3),pelev(3),rmet(3,3),tsec(2)
 real(dp) :: dmatdum(0,0,0,0)
 real(dp),allocatable :: grnl(:,:),mpibuf(:,:)
 real(dp),allocatable :: qphon(:),rhonow(:,:,:),sqnormgrhor(:,:)
 real(dp),allocatable :: xcart(:,:)

! *************************************************************************

 DBG_ENTER("COLL")

 call timab(250,1,tsec)
 call timab(251,1,tsec)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!MPI FFT communicator
 spaceComm_fft=mpi_enreg%comm_fft

!Recompute structure factor phases if atomic positions have changed
 if (moved_atm_inside==1) then
   if (psps%usepaw==1.and.pawfgr%usefinegrid==1) then
     call getph(atindx,dtset%natom,ngfft(1),ngfft(2),ngfft(3),ph1d,xred)
   else
     ph1d(:,:)=ph1df(:,:)
   end if
 end if

!----------------------------------------------------------------------
!Wavelet case: transform psi to KS orbitals
!----------------------------------------------------------------------
 if (dtset%usewvl == 1) then

!  wvlbigdft indicates that the BigDFT workflow will be followed
   if(dtset%wvl_bigdft_comp==1) wvlbigdft=.true.

#if defined HAVE_BIGDFT
!  Transform to KS orbitals
   comm=mpi_enreg%comm_cell
   me=xmpi_comm_rank(comm)
   nproc=xmpi_comm_size(comm)
!  do_last_ortho in case of direct minimization, since
!  we never diagonalized the hamiltonian and the eigenvalues are unknown.
   if(wvlbigdft) then
     do_last_ortho=dtset%iscf==0 !direct minimization scheme in case of wvlbigdft.
   else
     do_last_ortho=dtset%nnsclo==0 !direct minimization scheme.
   end if
   if ( do_last_ortho ) then
     call total_energies(wvl%e%energs, istep, me)
     call write_energies(istep,0,wvl%e%energs,zero,zero,"FINAL")
!    XG 121101 Suppressed, because wvl%energs%eKS was NaN on shiva distchk, for nearly all automatic tests...
!    etotal = wvl%energs%eKS
     if(.not. wvlbigdft) then
!      Store xcart for each atom
       ABI_ALLOCATE(xcart,(3, dtset%natom))
       call xred2xcart(dtset%natom, rprimd, xcart, xred)
!      If ISCF>10, we exit scfcv at a place where bigdft objects
!      do not contain the KS potential. Hence, we copy vtrial to wvl%den
       if(dtset%iscf>=10) then
         call wvl_vtrial_abi2big(1,vtrial,wvl%den)
       end if
!      hpsi is lost in hpsitopsi,
!      so we recalculate hpsi, last_orthon will be called afterwards.
       call wvl_psitohpsi(dtset%diemix,energies%e_exactX,energies%e_xc,energies%e_hartree,&
&       energies%e_kinetic, energies%e_localpsp, eproj, energies%e_sicdc,&
&       istep,1, dtset%iscf, mpi_enreg%me_wvl,dtset%natom,&
&       nfftf,mpi_enreg%nproc_wvl,dtset%nspden,&
&       raux,.false.,energies%e_xcdc,wvl,wvlbigdft,xcart,strsxc)
       ABI_DEALLOCATE(xcart)
     end if

     call last_orthon(me, nproc, istep, wvl%wfs%ks, wvl%e%energs%evsum, .true.)
     if (nproc == 1) nullify(wvl%wfs%ks%psit)

     call eigensystem_info(me,nproc,0.d0,&
     wvl%wfs%ks%Lzd%Glr%wfd%nvctr_c+7*wvl%wfs%ks%Lzd%Glr%wfd%nvctr_f,&
     wvl%wfs%ks%orbs,wvl%wfs%ks%psi)

!    Copy eigenvalues from BigDFT object to "eigen"
     call wvl_eigen_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,eigen,2,wvl%wfs)
!    Copy occupations from BigDFT objects to ABINIT
     call wvl_occ_abi2big(dtset%mband,dtset%nkpt,dtset%nsppol,occ,2,wvl%wfs)

   end if

!  Clean KSwfn parts only needed in the SCF loop.
   call kswfn_free_scf_data(wvl%wfs%ks, (nproc > 1))
!  Clean denspot parts only needed in the SCF loop.
   call denspot_free_history(wvl%den%denspot)
#else
   BIGDFT_NOTENABLED_ERROR()
#endif

!  Tail corrections, pending for wvlbigdft==.false.
!  TODO put it at the end of gstate.
!  WVL - maybe compute the tail corrections to energy
   if (dtset%tl_radius > tol12 .and. wvlbigdft) then
!    Store xcart for each atom
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     comm=mpi_enreg%comm_cell
     me=xmpi_comm_rank(comm)
     nproc=xmpi_comm_size(comm)
!    Use the tails to improve energy precision.
     call wvl_tail_corrections(dtset, energies, etotal, mpi_enreg, psps, wvl, xcart)
     ABI_DEALLOCATE(xcart)
   end if

   if(wvlbigdft) then
!    Change the density according to the KS projection.
     call wvl_mkrho(dtset, irrzon, mpi_enreg, phnons, rhor, wvl%wfs, wvl%den)
   end if
 end if

 call timab(251,2,tsec)
 call timab(252,1,tsec)

!----------------------------------------------------------------------
!Polarization Calculation
!----------------------------------------------------------------------

 if(dtset%berryopt/=0)then
   call elpolariz(atindx1,cg,cprj,dtefield,dtfil,dtset,etotal,energies%e_elecfield,gprimd,hdr,&
&   kg,dtset%mband,mcg,mcprj,dtset%mkmem,mpi_enreg,dtset%mpw,my_natom,dtset%natom,nattyp,dtset%nkpt,&
&   npwarr,dtset%nsppol,psps%ntypat,pawrhoij,pawtab,pel,pel_cg,pelev,pion,&
&   psps,pwind,pwind_alloc,pwnsfac,rprimd,ucvol,usecprj,xred)
 end if


 call timab(252,2,tsec)
 call timab(253,1,tsec)

!----------------------------------------------------------------------
!Gradient and Laplacian of the Density Calculation
!----------------------------------------------------------------------

!We use routine xcden which get gradient of rhor (grhor), and eventually laplacian of rhor (lrhor).
 if(dtset%prtgden/=0 .or. dtset%prtlden/=0)then

!  Compute gradient of the electron density
   ngrad=2
   cplex=1
   ishift=0
   ABI_ALLOCATE(rhonow,(nfftf,dtset%nspden,ngrad*ngrad))
   if(dtset%prtlden/=0)then
     nullify(lrhor)
     ABI_ALLOCATE(lrhor,(nfftf,dtset%nspden))
   end if
   write(message,'(a,a)') ch10, " Compute gradient of the electron density"
   call wrtout(ab_out,message,'COLL')
   if(dtset%prtlden/=0) then
     write(message,'(a)') " and also Compute Laplacian of the electron density"
     call wrtout(ab_out,message,'COLL')
   end if
   write(message,'(a)') "--------------------------------------------------------------------------------"
   call wrtout(ab_out,message,'COLL')

   ABI_ALLOCATE(qphon,(3))
   qphon(:)=zero
   if(dtset%prtlden/=0) then
     call xcden (cplex,gprimd,ishift,mpi_enreg,nfftf,ngfftf,ngrad,dtset%nspden,dtset%paral_kgb,qphon,rhor,rhonow,lrhonow=lrhor)
   else
     call xcden (cplex,gprimd,ishift,mpi_enreg,nfftf,ngfftf,ngrad,dtset%nspden,dtset%paral_kgb,qphon,rhor,rhonow)
   end if
   ABI_DEALLOCATE(qphon)

!  Copy gradient which has been output in rhonow to grhor (and free rhonow)
   nullify(grhor)
   ABI_ALLOCATE(grhor,(nfftf,dtset%nspden,3))
   do ispden=1,dtset%nspden
     do ifft=1,nfftf
       grhor(ifft,ispden,1:3) = rhonow(ifft,ispden,2:4)
     end do
   end do
   ABI_DEALLOCATE(rhonow)

   if(dtset%prtgden/=0) then
!    Print result for grhor
     write(message,'(a,a)') ch10, " Result for gradient of the electron density for each direction (1,2,3):"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,a,a)') ch10," 1rst direction:",ch10,&
&     "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,grhor(:,:,1),optrhor=2,ucvol=ucvol)
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,a,a)') ch10," 2nd direction:",ch10,&
&     "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,grhor(:,:,2),optrhor=2,ucvol=ucvol)
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a,a,a)') ch10," 3rd direction:",ch10,&
&     "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,grhor(:,:,3),optrhor=2,ucvol=ucvol)
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
   end if

   if(dtset%prtlden/=0) then
!    Print result for lrhor
     write(message,'(a,a)') ch10, " Result for Laplacian of the electron density :"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a)') ch10, "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,lrhor,optrhor=3,ucvol=ucvol)
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
   end if

   write(message,'(a)') "--------------------------------------------------------------------------------"
   call wrtout(ab_out,message,'COLL')
 end if

!----------------------------------------------------------------------
!Kinetic Energy Density Calculation
!----------------------------------------------------------------------

 call timab(253,2,tsec)
 call timab(254,1,tsec)

!We use routine mkrho with option=1 to compute kinetic energy density taur (and taug)
 if(dtset%prtkden/=0 .or. dtset%prtelf/=0)then
   nullify(taug,taur)
!  tauX are reused in outscfcv for output
!  should be deallocated there
   ABI_ALLOCATE(taug,(2,nfftf))
   ABI_ALLOCATE(taur,(nfftf,dtset%nspden))
   tim_mkrho=5
   if(dtset%prtelf/=0) then
     write(message,'(a,a)') ch10, " Compute ELF"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
   end if
   write(message,'(a,a)') ch10, " Compute kinetic energy density"
   call wrtout(ab_out,message,'COLL')
   paw_dmft%use_sc_dmft=0 ! dmft not used here
   paw_dmft%use_dmft=0 ! dmft not used here
   call mkrho(cg,dtset,gprimd,irrzon,kg,mcg,mpi_enreg,&
&   npwarr,occ,paw_dmft,phnons,taug,taur,rprimd,tim_mkrho,ucvol,wvl%den,wvl%wfs,option=1)
!  Print result
   if(dtset%prtkden/=0) then
     write(message,'(a,a)') ch10, "Result for kinetic energy density :"
     call wrtout(ab_out,message,'COLL')
     write(message,'(a,a)') ch10, "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,taur,optrhor=1,ucvol=ucvol)
     write(message,'(a)') "--------------------------------------------------------------------------------"
     call wrtout(ab_out,message,'COLL')
   end if
   ABI_DEALLOCATE(taug)
 end if

!----------------------------------------------------------------------
!Electron Localization Function (ELF) Calculation
!----------------------------------------------------------------------

 call timab(254,2,tsec)
 call timab(255,1,tsec)

!We use routine xcden to compute gradient of electron density (grhor),
!NOTE: If GGA is used, gradient of electron density is already computed
!and it is stored in exchange correlation kernel kxc(:,18:23).
!In order to save memory and do not have the same quantity twice
!in memory we should use kxc.
!But unfortunately only spin up ans spin down component are stored in kxc.
!So far we thus use grhor which contains all component (nspden=4)
!just like rhonow in xcden instead of kxc.

 if((dtset%prtelf/=0))then
   if(dtset%nspden<=2) then

     ngrad=2
     cplex=1
     if((cplex*dtset%nfft)/=nfftf)then
       write(message, '(a,a,a,a)' ) ch10,&
&       ' afterscfloop: ERROR -', ch10, &
&       '   The density is complex, ELF analysis cannot be performed.'
       call wrtout(std_out,message,'COLL')
!      MSG_ERROR(message)
     end if

     if((dtset%prtgden==0) .and. (dtset%prtlden==0)) then
!      Compute gradient of the electron density
       ishift=0
       ABI_ALLOCATE(rhonow,(nfftf,dtset%nspden,ngrad*ngrad))
       write(message,'(a,a)') ch10, " Compute gradient of the electron density"
       call wrtout(ab_out,message,'COLL')
       ABI_ALLOCATE(qphon,(3))
       qphon(:)=zero
       call xcden (cplex,gprimd,ishift,mpi_enreg,nfftf,ngfftf,ngrad,dtset%nspden,dtset%paral_kgb,qphon,rhor,rhonow)
       ABI_DEALLOCATE(qphon)
!      Copy gradient which has been output in rhonow to grhor (and free rhonow)
       ABI_ALLOCATE(grhor,(nfftf,dtset%nspden,3))
       do ispden=1,dtset%nspden
         do ifft=1,nfftf
           grhor(ifft,ispden,1:3) = rhonow(ifft,ispden,2:4)
         end do
       end do
       ABI_DEALLOCATE(rhonow)
     end if
!    Compute square norm of gradient of the electron density (|grhor|**2)
     if(dtset%nspden==1)then
       ABI_ALLOCATE(sqnormgrhor,(nfftf,dtset%nspden))
       do ifft=1,nfftf
         sqnormgrhor(ifft,1) = zero
       end do
     elseif(dtset%nspden==2)then
       ABI_ALLOCATE(sqnormgrhor,(nfftf,dtset%nspden+1))
!      because we not only want (total and up quantities, but also down)
!      Indeed after having token the square norm we can not recover the
!      down quantity by substracting total and up quantities (as we do for usual densities)
       do ispden=1,dtset%nspden+1
         do ifft=1,nfftf
           sqnormgrhor(ifft,ispden) = zero
         end do
       end do
     end if

     do igrad=1,3
       do ispden=1,dtset%nspden !total (and up)
         do ifft=1,nfftf
           sqnormgrhor(ifft,ispden) = sqnormgrhor(ifft,ispden) + grhor(ifft,ispden,igrad)**2
         end do
       end do
       if(dtset%nspden==2)then
         do ifft=1,nfftf !down
           sqnormgrhor(ifft,3) = sqnormgrhor(ifft,3) + (grhor(ifft,1,igrad)-grhor(ifft,2,igrad))**2
         end do
       end if
     end do

!    Compute electron localization function (ELF) (here it is elfr)

     nullify(elfr)
     if(dtset%nspden==1)then
       ABI_ALLOCATE(elfr,(nfftf,dtset%nspden))
     elseif(dtset%nspden==2)then
       ABI_ALLOCATE(elfr,(nfftf,dtset%nspden+1))
!      1rst is total elf, 2nd is spin-up elf, and 3rd is spin-down elf. (elf_tot /= elf_up + elf_down)
     end if
     c_fermi = 3.0d0/10.0d0*((3.0d0*pi**2)**(2.0d0/3.0d0))

!    First compute total elf
     ispden=1
     do ifft=1,nfftf
       dtaurzero = c_fermi*rhor(ifft,ispden)**(5.0d0/3.0d0)
       dtaur = taur(ifft,ispden)
       dtaur = dtaur - (1.0d0/8.0d0)*(sqnormgrhor(ifft,ispden)/rhor(ifft,ispden))
!      Ensure that dtaur is always positive or zero, as it should be.
       if(dtaur<0.0d0)dtaur=0.0d0
!      To avoid NaN values we check that dtaurzero is not to small compare to dtaur
       if(dtaurzero<(1.0d-20*dtaur)) then
         elfr(ifft,ispden) = 0.0d0
       else
         elfr(ifft,ispden) = 1.0d0/(1.0d0 + (dtaur/dtaurzero)**2)
!        For atomic shell studies we could also add the condition that when (dtaur/dtaurzero)
!        is very close to zero we actually set it exactly to zero (or --> elfr = 1.0d0
!        which is the upper limit of elfr.)
       end if
     end do

!    If spin-dependent densities are avalaible, compute spin-dependent elf
!    and offer the possibility to compute total elf in an alternative approach
!    (see doc/theory/ELF)
     if(dtset%nspden==2)then

!      alternative approach to the total elf
       if(dtset%prtelf==2)then
         ispden=1
         do ifft=1,nfftf
           dtaurzero = 2.0d0**(2.0d0/3.0d0)*c_fermi*( rhor(ifft,ispden+1)**(5.0d0/3.0d0) + &
&           (rhor(ifft,ispden) - rhor(ifft,ispden+1))**(5.0d0/3.0d0) )
           dtaur = taur(ifft,ispden)
           dtaur = dtaur - (1.0d0/8.0d0)*(sqnormgrhor(ifft,ispden+1)/rhor(ifft,ispden+1)) &
&           - (1.0d0/8.0d0)*(sqnormgrhor(ifft,ispden+2)/(rhor(ifft,ispden)-rhor(ifft,ispden+1)))
           if(dtaur<0.0d0)dtaur=0.0d0
!          To avoid NaN values we check that dtaurzero is not to small compare to dtaur
           if(abs(dtaurzero)<abs(1.0d-20*dtaur)) then
             elfr(ifft,ispden) = 0.0d0
           else
             elfr(ifft,ispden) = 1.0d0/(1.0d0 + (dtaur/dtaurzero)**2)
           end if
         end do
       end if

!      elf_up
       ispden=2
       do ifft=1,nfftf
         dtaurzero = 2.0d0**(2.0d0/3.0d0)*c_fermi*rhor(ifft,ispden)**(5.0d0/3.0d0)
         dtaur = taur(ifft,ispden)
         dtaur = dtaur - (1.0d0/8.0d0)*(sqnormgrhor(ifft,ispden)/rhor(ifft,ispden))
         if(dtaur<0.0d0)dtaur=0.0d0
!        To avoid NaN values we check that dtaurzero is not to small compare to dtaur
         if(abs(dtaurzero)<abs(1.0d-20*dtaur)) then
           elfr(ifft,ispden) = 0.0d0
         else
           elfr(ifft,ispden) = 1.0d0/(1.0d0 + (dtaur/dtaurzero)**2)
         end if
       end do

!      elf_down
       ispden=1
       do ifft=1,nfftf
         dtaurzero = 2.0d0**(2.0d0/3.0d0)*c_fermi*(rhor(ifft,ispden)-rhor(ifft,ispden+1))**(5.0d0/3.0d0)
         dtaur = taur(ifft,ispden)-taur(ifft,ispden+1)
         dtaur = dtaur - (1.0d0/8.0d0)*(sqnormgrhor(ifft,ispden+2)/(rhor(ifft,ispden)-rhor(ifft,ispden+1)))
         if(dtaur<0.0d0)dtaur=0.0d0
!        To avoid NaN values we check that dtaurzero is not to small compare to dtaur
         if(abs(dtaurzero)<abs(1.0d-20*dtaur)) then
           elfr(ifft,ispden+2) = 0.0d0
         else
           elfr(ifft,ispden+2) = 1.0d0/(1.0d0 + (dtaur/dtaurzero)**2)
         end if
       end do

     end if !endif dtset%nspden==2

!    Print result for elfr
     call prtrhomxmn(ab_out,mpi_enreg,nfftf,ngfftf,dtset%nspden,1,elfr,optrhor=4,ucvol=ucvol)

     ABI_DEALLOCATE(grhor)
     ABI_DEALLOCATE(sqnormgrhor)

   else
     message ='ELF is not yet implemented for non collinear spin cases.'
     MSG_WARNING(message)

     ABI_ALLOCATE(elfr,(nfftf,dtset%nspden))
     do ispden=1,dtset%nspden
       do ifft=1,nfftf
         elfr(ifft,ispden) = -2.0d0
       end do
     end do
!    even if elf is not computed we want to finish the abinit run.
!    To ensure that users won't use the _ELF file which will be produced
!    we set elf to -2.0 (a meaningless value)

   end if ! endif dtset%nspden<=2

   write(message,'(a,a)') ch10, "--------------------------------------------------------------------------------"
   call wrtout(ab_out,message,'COLL')
   write(message,'(a)') " End of ELF section"
   call wrtout(ab_out,message,'COLL')

 end if !endif prtelf/=0

!######################################################################
!Compute forces (if they were not computed during the elec. iterations)
!and stresses (if requested by user)
!----------------------------------------------------------------------

 call timab(255,2,tsec)
 call timab(256,1,tsec)

 optfor=0

 if (computed_forces==0.and.dtset%optforces>0.and.dtset%iscf>=0) then
   if (dtset%nstep>0.or.dtfil%ireadwf==1) optfor=1
 end if

 if (optfor>0.or.stress_needed>0) then

!  PAW: eventually, compute g_l(r).Y_lm(r) gradients (if not already done)
   if (psps%usepaw==1) then
     test_nfgd  =any(pawfgrtab(:)%nfgd==0)
     test_rfgd  =any(pawfgrtab(:)%rfgd_allocated==0)
     test_gylmgr=any(pawfgrtab(:)%gylmgr_allocated==0)
     if (test_nfgd.or.&
&     (test_gylmgr.and.dtset%pawstgylm==1).or.&
&     (test_rfgd.and.stress_needed==1.and.dtset%pawstgylm==1).or.&
     (test_rfgd.and.dtset%pawstgylm==0)) then
       optcut=0;optgr0=0;optgr1=dtset%pawstgylm;optgr2=0
       optrad=1-dtset%pawstgylm;if (stress_needed==1) optrad=1
       call nhatgrid(atindx1,gmet,my_natom,dtset%natom,nattyp,ngfftf,dtset%ntypat,&
&       optcut,optgr0,optgr1,optgr2,optrad,pawfgrtab,pawtab,rprimd,dtset%typat,ucvol,xred,&
&       comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&       comm_fft=spaceComm_fft,distribfft=mpi_enreg%distribfft)
     end if
   end if

   if (dtset%usewvl == 0 .or. dtset%berryopt == 4 .or. dtset%berryopt == 6 .or. dtset%berryopt == 7 .or.  &
&   dtset%berryopt == 14 .or. dtset%berryopt == 16 .or. dtset%berryopt == 17)  then

     call forstr(atindx1,cg,cprj,diffor,dtefield,dtset,&
&     eigen,electronpositron,energies,favg,fcart,fock,forold,fred,gresid,grewtn,&
&     grhf,grvdw,grxc,gsqcut,indsym,&
&     kg,kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,&
&     n3xccc,nattyp,nfftf,ngfftf,ngrvdw,nhat,nkxc,&
&     npwarr,dtset%ntypat,nvresid,occ,optfor,optres,&
&     paw_ij,pawang,pawfgr,pawfgrtab,pawrhoij,pawtab,ph1d,ph1df,&
&     psps,rhog,rhor,rprimd,stress_needed,strsxc,strten,symrec,synlgr,&
&     ucvol,usecprj,vhartr,vpsp,vxc,wvl%descr,wvl%den,&
&     xccc3d,xred,ylm,ylmgr,qvpotzero)
   else
     ABI_ALLOCATE(xcart,(3, dtset%natom))
     call xred2xcart(dtset%natom, rprimd, xcart, xred)
     ABI_ALLOCATE(grnl,(3, dtset%natom))
     call wvl_nl_gradient(grnl, mpi_enreg, dtset%natom, rprimd, wvl, xcart)
     ABI_DEALLOCATE(xcart)
     call forces(atindx1, diffor, dtefield, dtset, favg, fcart,fock, forold, fred, gresid, grewtn,&
&     grhf, grnl, grvdw, grxc, gsqcut, indsym, maxfor, mgfftf, mpi_enreg, &
&     0, n3xccc, nattyp, dtset%nfft, dtset%ngfft, ngrvdw, &
&     dtset%ntypat, pawtab, ph1d, psps, rhog, rhor, rprimd, symrec, &
&     synlgr, dtset%usefock, nvresid, vxc, wvl%descr,wvl%den, xred)
     ABI_DEALLOCATE(grnl)
   end if
 end if
 if (optfor==1) computed_forces=1
 if (optfor==1) diffor=9.9999999999D99
 if (stress_needed==0) strten(:)=9.9999999999D99
 if (computed_forces==0) fcart(:,:)=9.9999999999D99
 if (dtset%prtstm/=0) strten(:)=zero

 call timab(256,2,tsec)
 call timab(257,1,tsec)

!If SCF convergence was not reached (for dtset%nstep>0),
!print a warning to the output file (non-dummy arguments: dtset%nstep,
!residm, diffor - infos from tollist have been saved inside )
 choice=3
 call scprqt(choice,cpus,deltae,diffor,dtset,&
& eigen,etotal,favg,fcart,energies%e_fermie,dtfil%fnameabo_app_eig,dtfil%filnam_ds(1),&
& 1,dtset%iscf,istep,dtset%kptns,maxfor,&
& moved_atm_inside,mpi_enreg,dtset%nband,dtset%nkpt,&
& dtset%nstep,occ,optres,prtfor,prtxml,quit,&
& res2,resid,residm,response,tollist,psps%usepaw,vxcavg,dtset%wtk,xred,conv_retcode,&
& electronpositron=electronpositron)

!output POSCAR and FORCES files, VASP style, for PHON code and friends.
 if (dtset%prtposcar == 1) then
   call prtposcar(fcart, dtfil%filnam_ds(4), dtset%natom, dtset%ntypat, rprimd, dtset%typat, ucvol, xred, dtset%znucl)
 end if ! prtposcar

 if(allocated(qphon))   then
   ABI_DEALLOCATE(qphon)
 end if

!get current operator on wavefunctions
 if (dtset%prtspcur == 1) then
   call spin_current(cg,dtfil,dtset,gprimd,hdr,kg,mcg,mpi_enreg,psps)
 end if

!Electron-positron stuff: if last calculation was a positron minimization,
!exchange electron and positron data in order to
!get electronic quantities in global variables
 if (dtset%positron/=0) then
   electronpositron%scf_converged=.false.
   if (dtset%positron<0.and.electronpositron_calctype(electronpositron)==1) then
     call exchange_electronpositron(cg,cprj,dtset,eigen,electronpositron,energies,fred,mcg,mcprj,&
&     mpi_enreg,my_natom,nfftf,ngfftf,nhat,npwarr,occ,paw_an,pawrhoij,rhog,rhor,strten,usecprj,vhartr)
   end if
 end if

!If PAW+U and density mixing, has to update nocc_mmp
 if (psps%usepaw==1.and.dtset%usepawu>0.and.(dtset%iscf>0.or.dtset%iscf==-3)) then
   call setnoccmmp(1,0,dmatdum,0,0,indsym,my_natom,dtset%natom,dtset%natpawu,&
&   dtset%nspinor,dtset%nsppol,dtset%nsym,dtset%ntypat,paw_ij,pawang,dtset%pawprtvol,&
&   pawrhoij,pawtab,dtset%spinat,dtset%symafm,dtset%typat,0,dtset%usepawu,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

!Update the content of the header (evolving variables)
 bantot=hdr%bantot
 if (dtset%positron==0) then
   call hdr_update(hdr,bantot,etotal,energies%e_fermie,residm,rprimd,occ,&
&   pawrhoij,xred,dtset%amu_orig(:,1),&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 else
   call hdr_update(hdr,bantot,electronpositron%e0,energies%e_fermie,residm,rprimd,occ,&
&   pawrhoij,xred,dtset%amu_orig(:,1),&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
 end if

#ifdef HAVE_LOTF
 if(dtset%ionmov==23 .and. mpi_enreg%nproc_band>1) then
   bufsz=2+2*dtset%natom;if (moved_atm_inside==1) bufsz=bufsz+dtset%natom
   ABI_ALLOCATE(mpibuf,(3,bufsz))
   mpibuf(:,1:dtset%natom)=fred(:,1:dtset%natom)
   mpibuf(:,dtset%natom+1:2*dtset%natom)=fcart(:,1:dtset%natom)
   if (moved_atm_inside==1) mpibuf(:,2*dtset%natom+1:3*dtset%natom)=xred(:,1:dtset%natom)
   mpibuf(1:3,bufsz-1:bufsz)=reshape(strten(1:6),(/3,2/))
   call xmpi_sum(mpibuf,mpi_enreg%comm_band,ierr)
   fred(:,1:dtset%natom)=mpibuf(:,1:dtset%natom)/mpi_enreg%nproc_band
   fcart(:,1:dtset%natom)=mpibuf(:,dtset%natom+1:2*dtset%natom)/mpi_enreg%nproc_band
   if (moved_atm_inside==1) xred(:,1:dtset%natom)=mpibuf(:,2*dtset%natom+1:3*dtset%natom)/mpi_enreg%nproc_band
   strten(1:6)=reshape(mpibuf(1:3,bufsz-1:bufsz),(/6/))/mpi_enreg%nproc_band
   ABI_DEALLOCATE(mpibuf)
 end if
#endif

!In case of FFT parallelisation, has to synchronize positions and forces
!to avoid numerical noise
 if (mpi_enreg%nproc_fft>1) then
   bufsz=2+2*dtset%natom;if (moved_atm_inside==1) bufsz=bufsz+dtset%natom
   ABI_ALLOCATE(mpibuf,(3,bufsz))
   mpibuf(:,1:dtset%natom)=fred(:,1:dtset%natom)
   mpibuf(:,dtset%natom+1:2*dtset%natom)=fcart(:,1:dtset%natom)
   if (moved_atm_inside==1) mpibuf(:,2*dtset%natom+1:3*dtset%natom)=xred(:,1:dtset%natom)
   mpibuf(1:3,bufsz-1:bufsz)=reshape(strten(1:6),(/3,2/))
   call xmpi_sum(mpibuf,mpi_enreg%comm_fft,ierr)
   fred(:,1:dtset%natom)=mpibuf(:,1:dtset%natom)/mpi_enreg%nproc_fft
   fcart(:,1:dtset%natom)=mpibuf(:,dtset%natom+1:2*dtset%natom)/mpi_enreg%nproc_fft
   if (moved_atm_inside==1) xred(:,1:dtset%natom)=mpibuf(:,2*dtset%natom+1:3*dtset%natom)/mpi_enreg%nproc_fft
   strten(1:6)=reshape(mpibuf(1:3,bufsz-1:bufsz),(/6/))/mpi_enreg%nproc_fft
   ABI_DEALLOCATE(mpibuf)
 end if

!results_gs%energies   = energies
 call energies_copy(energies,results_gs%energies)
 results_gs%etotal     =etotal
 results_gs%deltae     =deltae
 results_gs%diffor     =diffor
 results_gs%residm     =residm
 results_gs%res2       =res2
 results_gs%fcart(:,:) =fcart(:,:)
 results_gs%fred(:,:)  =fred(:,:)
 results_gs%gresid(:,:)=gresid(:,:)
 results_gs%grewtn(:,:)=grewtn(:,:)
 results_gs%grxc(:,:)  =grxc(:,:)
 results_gs%pel(1:3)   =pel(1:3)
 results_gs%strten(1:6)=strten(1:6)
 results_gs%synlgr(:,:)=synlgr(:,:)
 results_gs%vxcavg     =vxcavg
 if (ngrvdw>0) results_gs%grvdw(1:3,1:ngrvdw)=grvdw(1:3,1:ngrvdw)
 if (dtset%nstep == 0 .and. dtset%occopt>=3.and.dtset%occopt<=8) then
   results_gs%etotal = results_gs%etotal - dtset%tsmear * results_gs%entropy
 end if

 DBG_EXIT("COLL")

 call timab(257,2,tsec)
 call timab(250,2,tsec)

#if !defined HAVE_BIGDFT
 if (.false.) write(std_out,*) vtrial(1,1)
#endif
end subroutine afterscfloop
!!***

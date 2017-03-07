!{\src2tex{textfont=tt}}
!!****f* ABINIT/forstr
!! NAME
!! forstr
!!
!! FUNCTION
!! Drives the computation of forces and/or stress tensor
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (DCA, XG, GMR, MB, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
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
!!   | effmass=effective mass for electrons (1. in common case)
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
!!   | ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/input_variables/vargs.htm#ngfft
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
!!  nkxc=second dimension of the array kxc, see rhohxc.f for a description
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
!!      fock_updatecwaveocc,forces,forstrnps,nres2vres,pawgrnl,stress,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine forstr(atindx1,cg,cprj,diffor,dtefield,dtset,eigen,electronpositron,energies,favg,fcart,fock,&
&                 forold,fred,grchempottn,gresid,grewtn,grhf,grvdw,grxc,gsqcut,indsym,&
&                 kg,kxc,maxfor,mcg,mcprj,mgfftf,mpi_enreg,my_natom,n3xccc,nattyp,&
&                 nfftf,ngfftf,ngrvdw,nhat,nkxc,npwarr,&
&                 ntypat,nvresid,occ,optfor,optres,paw_ij,pawang,pawfgr,&
&                 pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,ph1df,psps,rhog,rhor,rprimd,stress_needed,&
&                 strsxc,strten,symrec,synlgr,ucvol,usecprj,vhartr,vpsp,&
&                 vxc,wvl,xccc3d,xred,ylm,ylmgr,qvpotzero)


 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_profiling_abi
 use m_efield
 use m_errors

 use m_electronpositron, only : electronpositron_type
 use m_energies,         only : energies_type
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_paw_ij,           only : paw_ij_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_pawfgr,           only : pawfgr_type
 use m_pawcprj,          only : pawcprj_type
 use m_fock,             only : fock_type,fock_updatecwaveocc

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'forstr'
 use interfaces_18_timing
 use interfaces_41_geometry
 use interfaces_62_wvl_wfs
 use interfaces_65_paw
 use interfaces_67_common, except_this_one => forstr
!End of the abilint section

 implicit none

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
 integer :: comm_grid,ifft,ispden,occopt_,optgr,optgr2,option,optnc,optstr,optstr2
 real(dp) :: dum,ucvol_
 logical :: apply_residual
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: kinstr(6),nlstr(6),tsec(2)
 real(dp) :: dummy(0)
 real(dp),allocatable :: grnl(:),vlocal(:,:),xcart(:,:)
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

!Compute nonlocal pseudopotential parts of forces and stress tensor
!-involves summations over wavefunctions at all k points
 if (dtset%tfkinfunc>0.and.stress_needed==1) then
   kinstr(1:3)=-two/three*energies%e_kinetic/ucvol ; kinstr(4:6)=zero
   nlstr(1:6)=zero
 else if (dtset%usewvl==0) then
   occopt_=0 ! This means that occ are now fixed
   if(dtset%usefock==1 .and. associated(fock)) then
     if((dtset%optforces/=0).or.(dtset%optstress/=0)) then
       call fock_updatecwaveocc(cg,cprj,dtset,fock,dum,indsym,fock%nnsclo_hf+1,mcg,mcprj,&
&       mpi_enreg,nattyp,npwarr,occ,ucvol)
     end if
   end if
   call forstrnps(cg,cprj,dtset%ecut,dtset%ecutsm,dtset%effmass,eigen,electronpositron,fock,grnl,&
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
#  if defined HAVE_BIGDFT
   if (dtset%usewvl==1) ucvol_=product(wvl%den%denspot%dpbox%hgrids)*real(product(wvl%den%denspot%dpbox%ndims),dp)
#  endif
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
&                  abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5)

!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residual before calling forces routine
   if (apply_residual) then
     ABI_ALLOCATE(resid,(nfftf,dtset%nspden))
     option=0; if (dtset%densfor_pred<0) option=1
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) optnc=2
     call nres2vres(dtset,gsqcut,psps%usepaw,kxc,mpi_enreg,my_natom,nfftf,ngfftf,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,psps%usepaw,resid,xccc3d,xred)
   else
     resid => nvresid
   end if

   call forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,forold,fred,gresid,grewtn,&
&     grhf,grnl,grvdw,grxc,gsqcut,indsym,maxfor,mgfftf,&
&     mpi_enreg,psps%n1xccc,n3xccc,nattyp,&
&     nfftf,ngfftf,ngrvdw,ntypat,pawrad,pawtab,ph1df,psps,rhog,&
&     rhor,rprimd,symrec,synlgr,dtset%usefock,resid,vxc,wvl%descr,wvl%den,xred,&
&     electronpositron=electronpositron)

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
   if (dtset%usefock==1 .and. associated(fock).and.fock%optstr) then
!write(80,*)fock%stress
     fock%stress(1:3)=fock%stress(1:3)-energies%e_fock/ucvol
!write(80,*) "forstr",-energies%e_fock/ucvol, energies%e_fock
   end if
   call stress(atindx1,dtset%berryopt,dtefield,energies%e_localpsp,dtset%efield,&
&   energies%e_hartree,energies%e_corepsp,fock,gsqcut,dtset%ixc,kinstr,mgfftf,&
&   mpi_enreg,psps%mqgrid_vl,psps%n1xccc,n3xccc,dtset%natom,nattyp,&
&   nfftf,ngfftf,nlstr,dtset%nspden,dtset%nsym,ntypat,dtset%paral_kgb,psps,pawrad,pawtab,ph1df,&
&   dtset%prtvol,psps%qgrid_vl,dtset%red_efieldbar,rhog,rprimd,strten,strsxc,symrec,&
&   dtset%typat,dtset%usefock,psps%usepaw,dtset%vdw_tol,dtset%vdw_tol_3bt,dtset%vdw_xc,psps%vlspl,&
&   vxc,psps%xccc1d,xccc3d,psps%xcccrc,xred,psps%ziontypat,psps%znucltypat,qvpotzero,&
&   electronpositron=electronpositron)
 end if

!Memory deallocation
 ABI_DEALLOCATE(grnl)

 call timab(914,2,tsec)
 call timab(910,2,tsec)

end subroutine forstr
!!***

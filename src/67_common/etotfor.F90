!{\src2tex{textfont=tt}}
!!****f* ABINIT/etotfor
!! NAME
!! etotfor
!!
!! FUNCTION
!! This routine is called to compute the total energy and various parts of it.
!! The routine computes -if requested- the forces.
!!
!! COPYRIGHT
!! Copyright (C) 1998-2017 ABINIT group (XG, GMR, MT)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt.
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtset <type(dataset_type)>=all input variables in this dataset
!!   | berryopt  = 4: electric field is on -> add the contribution of the
!!   |                - \Omega E.P term to the total energy
!!   |          /= 4: electric field is off
!!   | bfield = cartesian coordinates of magnetic field in atomic units
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | densfor_pred=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  gmet(3,3)=metric tensor for G vecs (in bohr**-2)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grchempottn(3,natom)=grads of spatially-varying chemical potential energy (hartree)
!!  grewtn(3,natom)=grads of Ewald energy (hartree)
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  nhat(nfft,nspden*usepaw)= -PAW only- compensation density
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optene=option for the computation of total energy
!!         (-1=no computation; 0=direct scheme; 1=double-counting scheme)
!!  optforces=option for the computation of forces
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom*usepaw) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) information.
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetry operations in reciprocal space
!!  ucvol=unit cell volume
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!  ===== if optforces==1
!!   diffor=maximum absolute change in component of forces between present and previous SCF cycle.
!!   favg(3)=mean of fcart before correction for translational symmetry
!!   fcart(3,natom)=cartesian forces from fred (hartree/bohr)
!!   fred(3,natom)=symmetrized form of grtn (grads of Etot) (hartree)
!!   gresid(3,natom)=forces due to the residual of the density/potential
!!   grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!   grxc(3,natom)=d(Exc)/d(xred) derivatives (0 without core charges)
!!   maxfor=maximum absolute value of force
!!   synlgr(3,natom)=symmetrized form of grads of Enl (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
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
!!   | e_hybcomp_E0(IN)=energy compensation energy for the hybrid functionals at frozen density
!!   | e_hybcomp_v0(IN)=potential compensation energy for the hybrid functionals at frozen density
!!   | e_hybcomp_v (IN)=potential compensation energy for the hybrid functionals at self-consistent density
!!   | e_kinetic(IN)=kinetic energy part of total energy.
!!   | e_nonlocalpsp(IN)=nonlocal pseudopotential part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!!   | e_magfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the magnetic field:  e_magfield = -ucvol*E*P
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  ===== if optforces==1
!!   forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!   grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!                 Input for norm-conserving psps, output for PAW
!!  ===== if psps%usepaw==1
!!   pawrhoij(my_natom) <type(pawrhoij_type)>= paw rhoij occupancies and related data
!!    (gradients of rhoij for each atom with respect to atomic positions are computed here)
!!
!! NOTES
!!  In case of PAW calculations:
!!    All computations are done on the fine FFT grid.
!!    All variables (nfft,ngfft,mgfft) refer to this fine FFT grid.
!!    All arrays (densities/potentials...) are computed on this fine FFT grid.
!!  ! Developpers have to be careful when introducing others arrays:
!!      they have to be stored on the fine FFT grid.
!!  In case of norm-conserving calculations the FFT grid is the usual FFT grid.
!!
!! PARENTS
!!      scfcv
!!
!! CHILDREN
!!      forces,nres2vres,pawgrnl,timab
!!
!! SOURCE

#if defined HAVE_CONFIG_H
#include "config.h"
#endif

#include "abi_common.h"

subroutine etotfor(atindx1,deltae,diffor,dtefield,dtset,&
&  elast,electronpositron,energies,&
&  etotal,favg,fcart,fock,forold,fred,gmet,grchempottn,gresid,grewtn,grhf,grnl,grvdw,&
&  grxc,gsqcut,indsym,kxc,maxfor,mgfft,mpi_enreg,my_natom,nattyp,&
&  nfft,ngfft,ngrvdw,nhat,nkxc,ntypat,nvresid,n1xccc,n3xccc,optene,optforces,optres,&
&  pawang,pawfgrtab,pawrad,pawrhoij,pawtab,ph1d,red_ptot,psps,rhog,rhor,rmet,rprimd,&
&  symrec,synlgr,ucvol,usepaw,vhartr,vpsp,vxc,wvl,wvl_den,xccc3d,xred)

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use defs_wvltypes
 use m_efield
 use m_profiling_abi
 use m_fock,             only : fock_type
 use m_pawang,           only : pawang_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_pawfgrtab,        only : pawfgrtab_type
 use m_pawrhoij,         only : pawrhoij_type
 use m_energies,         only : energies_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype

!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
#undef ABI_FUNC
#define ABI_FUNC 'etotfor'
 use interfaces_18_timing
 use interfaces_65_paw
 use interfaces_67_common, except_this_one => etotfor
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,mgfft,n1xccc,n3xccc,nfft,ngrvdw,nkxc,ntypat,optene,optforces
 integer,intent(in) :: optres,usepaw
 real(dp),intent(in) :: gsqcut
 real(dp),intent(inout) :: elast,ucvol
 real(dp),intent(out) :: deltae,diffor,etotal,maxfor
 type(MPI_type),intent(in) :: mpi_enreg
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer :: electronpositron
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: gmet(3,3),grchempottn(3,dtset%natom),grewtn(3,dtset%natom),grvdw(3,ngrvdw),kxc(nfft,nkxc)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom),red_ptot(3)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden),rmet(3,3),rprimd(3,3)
 real(dp),intent(in) :: vhartr(nfft),vpsp(nfft),vxc(nfft,dtset%nspden)
 real(dp),intent(in) :: xccc3d(n3xccc)
 real(dp),intent(inout) :: forold(3,dtset%natom),grnl(3*dtset%natom)
 real(dp),intent(inout) :: nhat(nfft,dtset%nspden*psps%usepaw)
 real(dp),intent(inout),target :: nvresid(nfft,dtset%nspden)
 real(dp),intent(inout) :: xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom)
 real(dp),intent(inout) :: fcart(3,dtset%natom)
 real(dp),intent(out) :: gresid(3,dtset%natom),grhf(3,dtset%natom)
 real(dp),intent(inout) :: grxc(3,dtset%natom)
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom*psps%usepaw)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom*psps%usepaw)
 type(pawrad_type),intent(in) :: pawrad(ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: comm_grid,dimnhat,ifft,ipositron,ispden,optgr,optgr2,option,optnc,optstr,optstr2,iir,jjr,kkr
 logical :: apply_residual
 real(dp) :: eenth,ucvol_
!arrays
 real(dp),parameter :: k0(3)=(/zero,zero,zero/)
 real(dp) :: tsec(2),A(3,3),A1(3,3),A_new(3,3),efield_new(3)
 real(dp) :: dummy(0),nhat_dum(0,0)
 real(dp),allocatable :: vlocal(:,:)
 real(dp), ABI_CONTIGUOUS pointer :: resid(:,:)

! *********************************************************************

 call timab(80,1,tsec)

 ipositron=electronpositron_calctype(electronpositron)

 if (optene>-1) then

!  When the finite-temperature VG broadening scheme is used,
!  the total entropy contribution "tsmear*entropy" has a meaning,
!  and gather the two last terms of Eq.8 of VG paper
!  Warning : might have to be changed for fixed moment calculations
   if(dtset%occopt>=3 .and. dtset%occopt<=8) then
     if (abs(dtset%tphysel) < tol10) then
       energies%e_entropy = - dtset%tsmear * energies%entropy
     else
       energies%e_entropy = - dtset%tphysel * energies%entropy
     end if
   else
     energies%e_entropy = zero
   end if

!  Turn it into an electric enthalpy, refer to Eq.(33) of Suppl. of Nat. Phys. paper (5,304,2009),
!    the missing volume is added here
   energies%e_elecfield=zero
   if ((dtset%berryopt==4.or.dtset%berryopt==14).and.ipositron/=1) then
     energies%e_elecfield=-dot_product(dtset%red_efieldbar,red_ptot)  !!ebar_i p_i
     eenth=zero
     do iir=1,3
       do jjr=1,3
         eenth=eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)  !!g^{-1})_ij ebar_i ebar_j
       end do
     end do
     energies%e_elecfield=energies%e_elecfield-eenth*ucvol/(8._dp*pi)
   end if

!  Turn it into an internal energy, refer to Eq.(36) of Suppl. of Nat. Phys. paper (5,304,2009),
!    but a little different: U=E_ks + (vol/8*pi) *  g^{-1})_ij ebar_i ebar_j
   if ((dtset%berryopt==6.or.dtset%berryopt==16).and.ipositron/=1) then
     energies%e_elecfield=zero
     eenth=zero
     do iir=1,3
       do jjr=1,3
         eenth=eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)  !! g^{-1})_ij ebar_i ebar_j
       end do
     end do
     energies%e_elecfield=energies%e_elecfield+eenth*ucvol/(8._dp*pi)
   end if

!  Calculate internal energy and electric enthalpy for mixed BC case.
   if (dtset%berryopt==17.and.ipositron/=1) then
     energies%e_elecfield=zero
     A(:,:)=(four_pi/ucvol)*rmet(:,:)
     A1(:,:)=A(:,:) ; A_new(:,:)=A(:,:)
     efield_new(:)=dtset%red_efield(:)
     eenth=zero
     do kkr=1,3
       if (dtset%jfielddir(kkr)==1) then    ! fixed ebar direction
!        step 1 add -ebar*p
         eenth=eenth-dtset%red_efieldbar(kkr)*red_ptot(kkr)
!        step 2  chang to e_new (change e to ebar)
         efield_new(kkr)=dtset%red_efieldbar(kkr)
!        step 3  chang matrix A to A1
         do iir=1,3
           do jjr=1,3
             if (iir==kkr .and. jjr==kkr) A1(iir,jjr)=-1.0/A(kkr,kkr)
             if ((iir==kkr .and. jjr/=kkr) .or.  (iir/=kkr .and.  jjr==kkr)) &
&             A1(iir,jjr)=-1.0*A(iir,jjr)/A(kkr,kkr)
             if (iir/=kkr .and. jjr/=kkr) A1(iir,jjr)=A(iir,jjr)-A(iir,kkr)*A(kkr,jjr)/A(kkr,kkr)
           end do
         end do
         A(:,:)=A1(:,:) ; A_new(:,:)=A1(:,:)
       end if
     end do  ! end fo kkr
     do iir=1,3
       do jjr=1,3
         eenth= eenth+half*A_new(iir,jjr)*efield_new(iir)*efield_new(jjr)
       end do
     end do
     energies%e_elecfield=energies%e_elecfield+eenth
   end if   ! berryopt==17

!  Turn it into a magnetic enthalpy, by adding orbital electronic contribution
   energies%e_magfield = zero
!  if (dtset%berryopt == 5 .and. ipositron/=1) then
!  emag = dot_product(mag_cart,dtset%bfield)
!  energies%e_magfield = emag
!  end if

!  Compute total (free)- energy by direct scheme
   if (optene==0) then
     etotal = energies%e_kinetic + energies%e_hartree + energies%e_xc + &
&     energies%e_localpsp + energies%e_corepsp + energies%e_fock+&
&     energies%e_entropy + energies%e_elecfield + energies%e_magfield+&
&     energies%e_hybcomp_E0 - energies%e_hybcomp_v0 + energies%e_hybcomp_v
     etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
     if (usepaw==0) etotal = etotal + energies%e_nonlocalpsp
     if (usepaw/=0) etotal = etotal + energies%e_paw
   end if

!  Compute total (free) energy by double-counting scheme
   if (optene==1) then
     etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc &
&     - energies%e_xcdc + energies%e_corepsp - energies%e_corepspdc+ energies%e_fock- energies%e_fockdc &
&     + energies%e_entropy + energies%e_elecfield + energies%e_magfield &
&     + energies%e_hybcomp_E0 - energies%e_hybcomp_v0
     etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
     if (usepaw/=0) etotal = etotal + energies%e_pawdc
   end if

!  Additional stuff for electron-positron
   if (dtset%positron/=0) then
     if (ipositron==0) then
       energies%e_electronpositron  =zero
       energies%edc_electronpositron=zero
     else
       energies%e_electronpositron  =electronpositron%e_hartree+electronpositron%e_xc
       energies%edc_electronpositron=electronpositron%e_hartree+electronpositron%e_xcdc
       if (usepaw==1) then
         energies%e_electronpositron  =energies%e_electronpositron  +electronpositron%e_paw
         energies%edc_electronpositron=energies%edc_electronpositron+electronpositron%e_pawdc
       end if
     end if
     if (optene==0) electronpositron%e0=etotal
     if (optene==1) electronpositron%e0=etotal-energies%edc_electronpositron
     etotal=electronpositron%e0+energies%e0_electronpositron+energies%e_electronpositron
   end if

!  Compute energy residual
   deltae=etotal-elast
   elast=etotal
 end if !optene/=-1

 call timab(80,2,tsec)

!------Compute forces-----------------------------------------------------

 if (optforces==1) then

!  PAW: add gradients due to Dij derivatives to non-local term
   if (usepaw==1) then
     ABI_ALLOCATE(vlocal,(nfft,dtset%nspden))
     do ispden=1,min(dtset%nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,vhartr,vlocal,vpsp,vxc)
       do ifft=1,nfft
         vlocal(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
       end do
     end do

     if(dtset%nspden==4)then
       do ispden=3,4
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,vlocal,vxc)
         do ifft=1,nfft
           vlocal(ifft,ispden)=vxc(ifft,ispden)
         end do
       end do
     end if
     ucvol_=ucvol
#if defined HAVE_BIGDFT
     if (dtset%usewvl==1) ucvol_=product(wvl_den%denspot%dpbox%hgrids)*real(product(wvl_den%denspot%dpbox%ndims),dp)
#endif
     dimnhat=0;optgr=1;optgr2=0;optstr=0;optstr2=0
     comm_grid=mpi_enreg%comm_fft;if(dtset%usewvl==1) comm_grid=mpi_enreg%comm_wvl
     call pawgrnl(atindx1,dimnhat,dummy,1,dummy,grnl,gsqcut,mgfft,my_natom, &
&     dtset%natom, nattyp,nfft,ngfft,nhat_dum,dummy,dtset%nspden,dtset%nsym,ntypat,optgr,optgr2,optstr,optstr2,&
&     pawang,pawfgrtab,pawrhoij,pawtab,ph1d,psps,k0,rprimd,symrec,dtset%typat,ucvol_,vlocal,vxc,xred, &
&     mpi_atmtab=mpi_enreg%my_atmtab,comm_atom=mpi_enreg%comm_atom,mpi_comm_grid=mpi_enreg%comm_fft,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,paral_kgb=mpi_enreg%paral_kgb)
     ABI_DEALLOCATE(vlocal)
   end if

   apply_residual=(optres==1 .and. dtset%usewvl==0.and.abs(dtset%densfor_pred)>=1 .and. &
&   abs(dtset%densfor_pred)<=6.and.abs(dtset%densfor_pred)/=5)

!  If residual is a density residual (and forces from residual asked),
!  has to convert it into a potential residual before calling forces routine
   if (apply_residual) then
     ABI_ALLOCATE(resid,(nfft,dtset%nspden))
     option=0; if (dtset%densfor_pred<0) option=1
     optnc=1;if (dtset%nspden==4.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) optnc=2
     call nres2vres(dtset,gsqcut,usepaw,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&
&     nkxc,nvresid,n3xccc,optnc,option,pawang,pawfgrtab,pawrhoij,pawtab,&
&     rhor,rprimd,usepaw,resid,xccc3d,xred,vxc)
   else
     resid => nvresid
   end if
   call forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,forold,fred,grchempottn,gresid,grewtn,&
&   grhf,grnl,grvdw,grxc,gsqcut,indsym,maxfor,mgfft,mpi_enreg,&
&   n1xccc,n3xccc,nattyp,nfft,ngfft,ngrvdw,ntypat,&
&   pawrad,pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,dtset%usefock,resid,vxc,wvl,wvl_den,xred,&
&   electronpositron=electronpositron)
   if (apply_residual) then
     ABI_DEALLOCATE(resid)
   end if

!  Returned fred are full symmetrized gradients of Etotal
!  wrt reduced coordinates xred, d(Etotal)/d(xred)
!  Forces are contained in array fcart

 else   ! if optforces==0
   fcart=zero
   fred=zero
   favg=zero
   diffor=zero
   gresid=zero
   grhf=zero
   maxfor=zero
   synlgr=zero
 end if

 call timab(80,2,tsec)

end subroutine etotfor
!!***

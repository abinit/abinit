!{\src2tex{textfont=tt}}
!!****m* ABINIT/m_odamix
!! NAME
!!  m_odamix
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2019 ABINIT group (FJ, MT)
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

module m_odamix

 use defs_basis
 use defs_datatypes
 use defs_abitypes
 use m_abicore
 use m_errors
 use m_xmpi
 use m_xcdata

 use m_time,       only : timab
 use m_geometry,   only : metric
 use m_cgtools,    only : dotprod_vn
 use m_pawang,     only : pawang_type
 use m_pawrad,     only : pawrad_type
 use m_pawtab,     only : pawtab_type
 use m_paw_an,     only : paw_an_type
 use m_paw_ij,     only : paw_ij_type
 use m_pawfgrtab,  only : pawfgrtab_type
 use m_pawrhoij,   only : pawrhoij_type,pawrhoij_filter
 use m_paw_nhat,   only : pawmknhat
 use m_paw_denpot, only : pawdenpot
 use m_energies,   only : energies_type
 use m_spacepar,   only : hartre
 use m_rhotoxc,    only : rhotoxc
 use m_fft,        only : fourdp

 implicit none

 private
!!***

 public :: odamix
!!***

contains
!!***

!!****f* ABINIT/odamix
!! NAME
!! odamix
!!
!! FUNCTION
!! This routine is called to compute the total energy and various parts of it.
!! The routine computes -if requested- the forces.
!!
!! INPUTS
!!  [add_tfw]=flag controling the addition of Weiszacker gradient correction to Thomas-Fermi kin energy
!!  dtset <type(dataset_type)>=all input variables in this dataset
!! berryopt  =  4/14: electric field is on -> add the contribution of the
!!                      -ebar_i p_i - Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!   = 6/16, or 7/17: electric displacement field is on  -> add the contribution of the
!!                      Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!   | berryopt  = 5: magnetic field is on -> add the contribution of the
!!   |                - \Omega B.M term to the total energy
!!   |          /= 5: magnetic field is off
!!   | bfield     = cartesian coordinates of the magnetic field in atomic units
!!   | dfield     = cartesian coordinates of the electric displacement field in atomic units (berryopt==6/7)
!!   | efield     = cartesian coordinates of the electric field in atomic units  (berryopt==4)
!!   | red_dfield = reduced the electric displacement field  (berryopt==16/17)
!!   | red_efieldbar = reduced the electric field (ebar)  (berryopt==14)
!!   | iatfix(3,natom)=1 for frozen atom along some direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | natom=number of atoms in cell.
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetry elements in space group
!!   | occopt=option for occupancies
!!   | prtvol=integer controlling volume of printed output
!!   | tsmear=smearing energy or temperature (if metal)
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!   | xclevel= XC functional level
!!  gprimd(3,3)=dimensional reciprocal space primitive translations
!!  gsqcut=cutoff on (k+G)^2 (bohr^-2)
!!  mpi_enreg=information about MPI parallelization
!!  my_natom=number of atoms treated by current processor
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nkxc=second dimension of the array kxc, see rhotoxc.f for a description
!!  ntypat=number of types of atoms in unit cell.
!!  nvresid(nfft,nspden)=potential or density residual
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  optres=0 if residual array (nvresid) contains the potential residual
!!        =1 if residual array (nvresid) contains the density residual
!!  paw_ij(my_natom) <type(paw_ij_type)>=paw arrays given on (i,j) channels
!!  paw_an(my_natom) <type(paw_an_type)>=paw arrays given on angular mesh
!!  pawang <type(pawang_type)>=paw angular mesh and related data
!!  pawfgrtab(my_natom) <type(pawfgrtab_type)>=atomic data given on fine rectangular grid
!!  pawrad
!!  pawtab(ntypat) <type(pawtab_type)>=paw tabulated starting data
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=array for Fourier transform of electron density
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  [taur(nfftf,nspden*dtset%usekden)]=array for kinetic energy density
!!  ucvol = unit cell volume (Bohr**3)
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vhartr(nfft)=array for holding Hartree potential
!!  vpsp(nfft)=array for holding local psp
!!  vxc(nfft,nspden)=array for holding XC potential
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  deltae=change in total energy
!!         between the previous and present SCF cycle
!!  etotal=total energy (hartree)
!!
!! SIDE EFFECTS
!! Input/Output:
!!  elast=previous value of the energy,
!!        needed to compute deltae, then updated.
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
!!   | e_nlpsp_vfock(IN)=nonlocal psp + potential Fock ACE part of total energy.
!!   | e_xc(IN)=exchange-correlation energy (hartree)
!!   | e_xcdc(IN)=exchange-correlation double-counting energy (hartree)
!!   | e_paw(IN)=PAW spherical part energy
!!   | e_pawdc(IN)=PAW spherical part double-counting energy
!!   | e_elecfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the electric field:  enefield = -ucvol*E*P
!!   | e_magfield(OUT)=the term of the energy functional that depends explicitely
!!   |                  on the magnetic field:  enmagfield = -ucvol*B*M
!!   | e_entropy(OUT)=entropy energy due to the occupation number smearing (if metal)
!!   |                this value is %entropy * dtset%tsmear (hartree).
!!  kxc(nfft,nkxc)=exchange-correlation kernel, needed only if nkxc>0
!!  [vxctau(nfftf,dtset%nspden,4*dtset%usekden)]=derivative of XC energy density with respect to
!!      kinetic energy density (metaGGA cases) (optional output)
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
!!      dotprod_vn,fourdp,hartre,metric,pawdenpot,pawmknhat,rhotoxc,timab
!!      xcdata_init,xmpi_sum
!!
!! SOURCE

subroutine odamix(deltae,dtset,elast,energies,etotal,&
&          gprimd,gsqcut,kxc,mpi_enreg,my_natom,nfft,ngfft,nhat,&
&          nkxc,ntypat,nvresid,n3xccc,optres,paw_ij,&
&          paw_an,pawang,pawfgrtab,pawrad,pawrhoij,pawtab,&
&          red_ptot,psps,rhog,rhor,rprimd,strsxc,ucvol,usepaw,&
&          usexcnhat,vhartr,vpsp,vtrial,vxc,vxcavg,xccc3d,xred,&
&          taur,vxctau,add_tfw) ! optional arguments

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: my_natom,n3xccc,nfft,nkxc,ntypat,optres
 integer,intent(in) :: usepaw,usexcnhat
 logical,intent(in),optional :: add_tfw
 real(dp),intent(in) :: gsqcut,ucvol
 real(dp),intent(inout) :: elast
 real(dp),intent(out) :: deltae,etotal,vxcavg
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
 type(energies_type),intent(inout) :: energies
 type(pawang_type),intent(in) :: pawang
 type(pseudopotential_type),intent(in) :: psps
!arrays
 integer,intent(in) :: ngfft(18)
 logical :: add_tfw_
 real(dp),intent(in) :: gprimd(3,3)
 real(dp),intent(in) :: red_ptot(3),rprimd(3,3),vpsp(nfft),xred(3,dtset%natom)
 real(dp),intent(in),optional :: taur(nfft,dtset%nspden*dtset%usekden)
 real(dp),intent(inout) :: kxc(nfft,nkxc),nhat(nfft,dtset%nspden*usepaw)
 real(dp),intent(inout) :: nvresid(nfft,dtset%nspden),rhog(2,nfft)
 real(dp),intent(inout) :: rhor(nfft,dtset%nspden),vhartr(nfft)
 real(dp),intent(inout) :: vtrial(nfft,dtset%nspden),vxc(nfft,dtset%nspden)
 real(dp),intent(inout) :: xccc3d(n3xccc)
 real(dp),intent(out) :: strsxc(6)
 real(dp),intent(inout),optional :: vxctau(nfft,dtset%nspden,4*dtset%usekden)
 type(paw_an_type),intent(inout) :: paw_an(my_natom)
 type(paw_ij_type),intent(inout) :: paw_ij(my_natom)
 type(pawfgrtab_type),intent(inout) :: pawfgrtab(my_natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat)
 type(pawrhoij_type),intent(inout) :: pawrhoij(my_natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat)

!Local variables-------------------------------
!scalars
 integer :: cplex,iatom,ider,idir,ierr,ifft,ipert,irhoij,ispden,itypat,izero,iir,jjr,kkr
 integer :: jrhoij,klmn,klmn1,kmix,nfftot,nhatgrdim,nzlmopt,nk3xc,option,optxc
 logical :: nmxc,with_vxctau
 real(dp) :: alphaopt,compch_fft,compch_sph,doti,e1t10,e_ksnm1,e_xcdc_vxctau
 real(dp) :: eenth,fp0,gammp1,ro_dlt,ucvol_local
 character(len=500) :: message
 type(xcdata_type) :: xcdata
!arrays
 real(dp) :: A(3,3),A1(3,3),A_new(3,3),efield_new(3)
 real(dp) :: gmet(3,3),gprimdlc(3,3),qpt(3),rmet(3,3),tsec(2)
 real(dp),allocatable :: nhatgr(:,:,:),rhoijtmp(:,:)

! *********************************************************************

!DEBUG
!write(std_out,*)' odamix : enter'
!ENDDEBUG

 call timab(80,1,tsec)

!Check that usekden is not 0 if want to use vxctau
 with_vxctau = (present(vxctau).and.present(taur).and.(dtset%usekden/=0))

!To be adjusted for the call to rhotoxc
 add_tfw_=.false.;if (present(add_tfw)) add_tfw_=add_tfw
 nk3xc=1;nmxc=(dtset%usepaw==1.and.mod(abs(dtset%usepawu),10)==4)

!faire un test sur optres=1, usewvl=0, nspden=1,nhatgrdim
 if(optres/=1)then
   write(message,'(a,i0,a)')' optres=',optres,', not allowed in oda => stop '
   MSG_ERROR(message)
 end if

 if(dtset%usewvl/=0)then
   write(message,'(a,i0,a)')' usewvl=',dtset%usewvl,', not allowed in oda => stop '
   MSG_ERROR(message)
 end if

 if(dtset%nspden/=1)then
   write(message,'(a,i0,a)')'  nspden=',dtset%nspden,', not allowed in oda => stop '
   MSG_ERROR(message)
 end if

 if (my_natom>0) then
   if(paw_ij(1)%has_dijhat==0)then
     message = ' dijhat variable must be allocated in odamix ! '
     MSG_ERROR(message)
   end if
   if(paw_ij(1)%cplex_dij==2.or.paw_ij(1)%qphase==2)then
     message = ' complex dij not allowed in odamix! '
     MSG_ERROR(message)
   end if
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! calculation of f'(0)= Eband_new-EH_old-E_xcdc_old-Ek_old-E_loc_old-E_nonloc_old
!!!!!!!!!! save previous energy E(rho_tild_n)

 fp0=energies%e_eigenvalues-energies%h0-two*energies%e_hartree-energies%e_xcdc
 if (usepaw==1) then
   do iatom=1,my_natom
     ABI_CHECK(pawrhoij(iatom)%qphase==1,'ODA mixing not allowed with a Q phase in PAW objects!')
     itypat=pawrhoij(iatom)%itypat
     do ispden=1,pawrhoij(iatom)%nspden
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,ispden)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,ispden)-paw_ij(iatom)%dijhat(klmn,ispden))
         jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
       end do
       klmn1=1
       do klmn=1,pawrhoij(iatom)%lmn2_size
         ro_dlt=-pawrhoij(iatom)%rhoijres(klmn1,ispden)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,ispden)-paw_ij(iatom)%dijhat(klmn,ispden))
         klmn1=klmn1+pawrhoij(iatom)%cplex_rhoij
       end do
     end do
     if (paw_ij(iatom)%ndij>=2.and.pawrhoij(iatom)%nspden==1) then
       jrhoij=1
       do irhoij=1,pawrhoij(iatom)%nrhoijsel
         klmn=pawrhoij(iatom)%rhoijselect(irhoij)
         ro_dlt=pawrhoij(iatom)%rhoijp(jrhoij,1)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,2)-paw_ij(iatom)%dijhat(klmn,2))
         jrhoij=jrhoij+pawrhoij(iatom)%cplex_rhoij
       end do
       klmn1=1
       do klmn=1,pawrhoij(iatom)%lmn2_size
         ro_dlt=-pawrhoij(iatom)%rhoijres(klmn1,1)*pawtab(itypat)%dltij(klmn)
         e1t10=e1t10+ro_dlt*(paw_ij(iatom)%dij(klmn,2)-paw_ij(iatom)%dijhat(klmn,2))
         klmn1=klmn1+pawrhoij(iatom)%cplex_rhoij
       end do
       e1t10=half*e1t10
     end if
   end do
   if (mpi_enreg%nproc_atom>1) then
     call xmpi_sum(e1t10,mpi_enreg%comm_atom,ierr)
   end if
   fp0=fp0-e1t10
 end if
 e_ksnm1=etotal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calculation of quantities that do not depend on rho_n+1

!PAW: eventually recompute compensation density (and gradients)
 nhatgrdim=0
 if (usepaw==1) then
   ider=-1;if (dtset%xclevel==2.or.usexcnhat==0) ider=0
   if (dtset%xclevel==2.and.usexcnhat==1) ider=ider+2
   if (ider>0) then
     nhatgrdim=1
     ABI_ALLOCATE(nhatgr,(nfft,dtset%nspden,3))
   end if
   if (ider>=0) then
     ider=0;izero=0;cplex=1;ipert=0;idir=0;qpt(:)=zero
     call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
     nfft,ngfft,nhatgrdim,dtset%nspden,ntypat,pawang,pawfgrtab,&
&     nhatgr,nhat,pawrhoij,pawrhoij,pawtab,qpt,rprimd,ucvol,dtset%usewvl,xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&     distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
   end if
 end if

!------Compute Hartree and xc potentials----------------------------------
 nfftot=PRODUCT(ngfft(1:3))

 call hartre(1,gsqcut,usepaw,mpi_enreg,nfft,ngfft,rhog,rprimd,vhartr)

 call xcdata_init(xcdata,dtset=dtset)

!Compute xc potential (separate up and down if spin-polarized)
 optxc=1
 call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
& nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc,optxc,rhor,rprimd,strsxc,&
& usexcnhat,vxc,vxcavg,xccc3d,xcdata,taur=taur,vhartr=vhartr,vxctau=vxctau,add_tfw=add_tfw_)

!------Compute parts of total energy depending on potentials--------

 ucvol_local=ucvol

!Compute Hartree energy energies%e_hartree
 call dotprod_vn(1,rhor,energies%e_hartree,doti,nfft,nfftot,1,1,vhartr,ucvol_local,&
& mpi_comm_sphgrid=mpi_enreg%comm_fft)
 energies%e_hartree=half*energies%e_hartree


!Compute local psp energy energies%e_localpsp
 call dotprod_vn(1,rhor,energies%e_localpsp,doti,nfft,nfftot,1,1,vpsp,ucvol_local,&
& mpi_comm_sphgrid=mpi_enreg%comm_fft)

!Compute double-counting XC energy energies%e_xcdc
 call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local,&
& mpi_comm_sphgrid=mpi_enreg%comm_fft)
 if (with_vxctau) then
   call dotprod_vn(1,taur,e_xcdc_vxctau,doti,nfft,nfftot,dtset%nspden,1,&
&   vxctau(:,:,1),ucvol_local,mpi_comm_sphgrid=mpi_enreg%comm_fft)
   energies%e_xcdc=energies%e_xcdc+e_xcdc_vxctau
 end if

 if (usepaw/=0) then
   nzlmopt=dtset%pawnzlm; option=2
   do iatom=1,my_natom
     itypat=paw_ij(iatom)%itypat
     ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(itypat)%lmn2_size))
     paw_ij(iatom)%has_dijhartree=1
   end do
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,0,dtset%ixc,my_natom,dtset%natom,dtset%nspden,ntypat,&
&   dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang,dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,&
&   pawtab,dtset%pawxcdev,dtset%spnorbscl,dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   do iatom=1,my_natom
     ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
   end do
 end if

!When the finite-temperature VG broadening scheme is used,
!the total entropy contribution "tsmear*entropy" has a meaning,
!and gather the two last terms of Eq.8 of the VG paper
!Warning : might have to be changed for fixed moment calculations
 if(dtset%occopt>=3 .and. dtset%occopt<=8) then
   energies%e_entropy = - dtset%tsmear * energies%entropy
 else
   energies%e_entropy = zero
 end if
!Turn it into an electric enthalpy,refer to Eq.(33) of Suppl. of Nat. Phys. paper (5,304,2009) [[cite:Stengel2009]]
! the missing volume is added here
 energies%e_elecfield = zero
 if (dtset%berryopt == 4 .or. dtset%berryopt == 14 ) then             !!HONG

   energies%e_elecfield = -dot_product(dtset%red_efieldbar,red_ptot)

   call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
   eenth = zero
   do iir=1,3
     do jjr=1,3
       eenth= eenth+gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)         !! HONG g^{-1})_ij ebar_i ebar_j
     end do
   end do
   eenth=-1_dp*(ucvol_local/(8.0d0*pi))*eenth
   energies%e_elecfield = energies%e_elecfield + eenth

 end if

 energies%e_magfield = zero
!if (dtset%berryopt == 5) then
!emag = dot_product(mag_cart,dtset%bfield)
!energies%e_magfield = emag
!end if

!HONG  Turn it into an internal enthalpy, refer to Eq.(36) of Suppl. of Nat. Phys. paper (5,304,2009) [[cite:Stengel2009]], 
!but a little different: U=E_ks + (vol/8*pi) *  g^{-1})_ij ebar_i ebar_j
 if (dtset%berryopt == 6 .or. dtset%berryopt == 16 )  then
   energies%e_elecfield=zero
   call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
   do iir=1,3
     do jjr=1,3
       energies%e_elecfield = energies%e_elecfield + gmet(iir,jjr)*dtset%red_efieldbar(iir)*dtset%red_efieldbar(jjr)     !! HONG g^{-1})_ij ebar_i ebar_j
     end do
   end do
   energies%e_elecfield = ucvol_local/(8.0d0*pi)*energies%e_elecfield
 end if

!HONG  calculate internal energy and electric enthalpy for mixed BC case.
 if ( dtset%berryopt == 17 ) then
   energies%e_elecfield = zero
   call metric(gmet,gprimdlc,-1,rmet,rprimd,ucvol_local)
   A(:,:)=(4*pi/ucvol_local)*rmet(:,:)
   A1(:,:)=A(:,:)
   A_new(:,:)=A(:,:)
   efield_new(:)=dtset%red_efield(:)
   eenth = zero

   do kkr=1,3
     if (dtset%jfielddir(kkr)==1) then    ! fixed ebar direction

!      step 1 add -ebar*p
       eenth=eenth - dtset%red_efieldbar(kkr)*red_ptot(kkr)

!      step 2  chang to e_new (change e to ebar)
       efield_new(kkr)=dtset%red_efieldbar(kkr)

!      step 3  chang matrix A to  A1

       do iir=1,3
         do jjr=1,3
           if (iir==kkr .and. jjr==kkr) A1(iir,jjr)=-1.0/A(kkr,kkr)
           if ((iir==kkr .and. jjr/=kkr) .or.  (iir/=kkr .and.  jjr==kkr)) &
&           A1(iir,jjr)=-1.0*A(iir,jjr)/A(kkr,kkr)
           if (iir/=kkr .and. jjr/=kkr) A1(iir,jjr)=A(iir,jjr)-A(iir,kkr)*A(kkr,jjr)/A(kkr,kkr)
         end do
       end do

       A(:,:)=A1(:,:)
       A_new(:,:)=A1(:,:)
     end if

   end do  ! end fo kkr


   do iir=1,3
     do jjr=1,3
       eenth= eenth+(1/2.0)*A_new(iir,jjr)*efield_new(iir)*efield_new(jjr)
     end do
   end do

   energies%e_elecfield=energies%e_elecfield+eenth

 end if   ! berryopt==17

 etotal = energies%e_kinetic+ energies%e_hartree + energies%e_xc + &
& energies%e_localpsp + energies%e_nlpsp_vfock - energies%e_fock0 + energies%e_corepsp + &
& energies%e_entropy + energies%e_elecfield + energies%e_magfield
!etotal = energies%e_eigenvalues - energies%e_hartree + energies%e_xc - &
!& energies%e_xcdc + energies%e_corepsp + &
!& energies%e_entropy + energies%e_elecfield
 etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
 if (usepaw==1) then
   etotal = etotal + energies%e_paw
 end if

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! now, compute mixed densities

 gammp1=etotal-e_ksnm1-fp0
 if (fp0>0.d0) then
   write(std_out,*) "fp0 est positif"
!  stop
 end if
 write(std_out,*) "fp0 ",fp0
 alphaopt=-fp0/two/gammp1

 if (alphaopt>one.or.alphaopt<0.d0) alphaopt=one
 if (abs(energies%h0)<=tol10) alphaopt=one
 write(std_out,*) " alphaopt",alphaopt

 energies%h0=(one-alphaopt)*energies%h0 + alphaopt*(energies%e_kinetic+energies%e_localpsp)
 energies%h0=energies%h0 + alphaopt*energies%e_nlpsp_vfock

 rhor= rhor+(alphaopt-one)*nvresid
 call fourdp(1,rhog,rhor(:,1),-1,mpi_enreg,nfft,1,ngfft,0)

 if (usepaw==1) then
   if (my_natom>0) then
     ABI_ALLOCATE(rhoijtmp,(pawrhoij(1)%cplex_rhoij*pawrhoij(1)%lmn2_size,pawrhoij(1)%nspden))
   end if
   do iatom=1,my_natom
     rhoijtmp=zero
     if (pawrhoij(iatom)%cplex_rhoij==1) then
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         do ispden=1,pawrhoij(iatom)%nspden
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=pawrhoij(iatom)%rhoijselect(irhoij)
             rhoijtmp(klmn,ispden)=pawrhoij(iatom)%rhoijp(irhoij,ispden)
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           klmn=pawrhoij(iatom)%kpawmix(kmix)
           rhoijtmp(klmn,ispden)=rhoijtmp(klmn,ispden)+(alphaopt-one)*pawrhoij(iatom)%rhoijres(klmn,ispden)
         end do
       end do
     else
       if (pawrhoij(iatom)%lmnmix_sz<pawrhoij(iatom)%lmn2_size) then
         jrhoij=1
         do ispden=1,pawrhoij(iatom)%nspden
           do irhoij=1,pawrhoij(iatom)%nrhoijsel
             klmn=2*pawrhoij(iatom)%rhoijselect(irhoij)-1
             rhoijtmp(klmn:klmn+1,ispden)=pawrhoij(iatom)%rhoijp(jrhoij:jrhoij+1,ispden)
             jrhoij=jrhoij+2
           end do
         end do
       end if
       do ispden=1,pawrhoij(iatom)%nspden
         do kmix=1,pawrhoij(iatom)%lmnmix_sz
           klmn=2*pawrhoij(iatom)%kpawmix(kmix)-1
           rhoijtmp(klmn:klmn+1,ispden)=rhoijtmp(klmn:klmn+1,ispden) &
&           +(alphaopt-one)*pawrhoij(iatom)%rhoijres(klmn:klmn+1,ispden)
         end do
       end do
     end if
     call pawrhoij_filter(pawrhoij(iatom)%rhoijp,pawrhoij(iatom)%rhoijselect,pawrhoij(iatom)%nrhoijsel,&
&                         pawrhoij(iatom)%cplex_rhoij,pawrhoij(iatom)%qphase,pawrhoij(iatom)%lmn2_size,&
&                         pawrhoij(iatom)%nspden,rhoij_input=rhoijtmp)
   end do ! iatom
   if (allocated(rhoijtmp)) then
     ABI_DEALLOCATE(rhoijtmp)
   end if
 end if ! usepaw

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!! Calcul des quantites qui dependent de rho_tilde_n+1 (rho apres mixing)

 if (usepaw==1) then
   if (ider>=0) then
     izero=0;cplex=1;ipert=0;idir=0;qpt(:)=zero
     call pawmknhat(compch_fft,cplex,ider,idir,ipert,izero,gprimd,my_natom,dtset%natom,&
&     nfft,ngfft,nhatgrdim,dtset%nspden,ntypat,pawang,pawfgrtab,nhatgr,&
&     nhat,pawrhoij,pawrhoij,pawtab,qpt,rprimd,ucvol,dtset%usewvl,xred,&
&     comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab,&
&     comm_fft=mpi_enreg%comm_fft,paral_kgb=dtset%paral_kgb,me_g0=mpi_enreg%me_g0,&
&     distribfft=mpi_enreg%distribfft,mpi_comm_wvl=mpi_enreg%comm_wvl)
   end if
 end if

!------Compute Hartree and xc potentials----------------------------------

 call hartre(1,gsqcut,usepaw,mpi_enreg,nfft,ngfft,rhog,rprimd,vhartr)

!Compute xc potential (separate up and down if spin-polarized)
 optxc=1;if (nkxc>0) optxc=2
 call rhotoxc(energies%e_xc,kxc,mpi_enreg,nfft,ngfft,&
& nhat,usepaw,nhatgr,nhatgrdim,nkxc,nk3xc,nmxc,n3xccc,optxc,rhor,rprimd,strsxc,&
& usexcnhat,vxc,vxcavg,xccc3d,xcdata,taur=taur,vhartr=vhartr,vxctau=vxctau,add_tfw=add_tfw_)

 if (nhatgrdim>0)  then
   ABI_DEALLOCATE(nhatgr)
 end if

!------Compute parts of total energy depending on potentials--------

 ucvol_local = ucvol

!Compute Hartree energy energies%e_hartree
 call dotprod_vn(1,rhor,energies%e_hartree,doti,nfft,nfftot,1,1,vhartr,ucvol_local,&
& mpi_comm_sphgrid=mpi_enreg%comm_fft)
 energies%e_hartree=half*energies%e_hartree

!Compute double-counting XC energy energies%e_xcdc
 call dotprod_vn(1,rhor,energies%e_xcdc,doti,nfft,nfftot,dtset%nspden,1,vxc,ucvol_local,&
& mpi_comm_sphgrid=mpi_enreg%comm_fft)

 etotal=energies%h0+energies%e_hartree+energies%e_xc+energies%e_corepsp + &
& energies%e_entropy + energies%e_elecfield + energies%e_magfield
 etotal = etotal + energies%e_ewald + energies%e_chempot + energies%e_vdw_dftd
 if (usepaw==1) then
   do iatom=1,my_natom
     itypat=paw_ij(iatom)%itypat
     ABI_ALLOCATE(paw_ij(iatom)%dijhartree,(pawtab(itypat)%lmn2_size))
     paw_ij(iatom)%has_dijhartree=1
   end do
   call pawdenpot(compch_sph,energies%e_paw,energies%e_pawdc,0,dtset%ixc,my_natom,dtset%natom, &
&   dtset%nspden,ntypat,dtset%nucdipmom,nzlmopt,option,paw_an,paw_an,paw_ij,pawang, &
&   dtset%pawprtvol,pawrad,pawrhoij,dtset%pawspnorb,pawtab,dtset%pawxcdev,dtset%spnorbscl,&
&   dtset%xclevel,dtset%xc_denpos,ucvol,psps%znuclpsp,&
&   comm_atom=mpi_enreg%comm_atom,mpi_atmtab=mpi_enreg%my_atmtab)
   do iatom=1,my_natom
     ABI_DEALLOCATE(paw_ij(iatom)%dijhartree)
     paw_ij(iatom)%has_dijhartree=0
   end do
   etotal=etotal+energies%e_paw
 end if
!Compute energy residual
 deltae=etotal-elast
 elast=etotal

 do ispden=1,min(dtset%nspden,2)
!$OMP PARALLEL DO PRIVATE(ifft) SHARED(ispden,nfft,vhartr,vpsp,vxc)
   do ifft=1,nfft
     vtrial(ifft,ispden)=vhartr(ifft)+vpsp(ifft)+vxc(ifft,ispden)
   end do
 end do
 if(dtset%nspden==4) vtrial(:,3:4)=vxc(:,3:4)

 call timab(80,2,tsec)

!DEBUG
!write(std_out,*) 'eeig-ehart+enxc-enxcdc+eew+eii+eent+enefield+epawdc',eeig,ehart,enxc,enxcdc,eew,eii,eent,enefield,epawdc
!ENDEBUG

end subroutine odamix
!!***

end module m_odamix
!!***

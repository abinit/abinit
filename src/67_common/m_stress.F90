!!****m* ABINIT/m_stress
!! NAME
!!  m_stress
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!!  Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, FJ, MT)
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

module m_stress

 use defs_basis
 use m_efield
 use m_abicore
 use m_errors
 use m_xmpi

 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : metric, stresssym
 use m_fock,             only : fock_type
 use m_ewald,            only : ewald2
 use defs_datatypes,     only : pseudopotential_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use m_fft,              only : zerosym, fourdp
 use m_mpinfo,           only : ptabs_fourdp
 use m_vdw_dftd2,        only : vdw_dftd2
 use m_vdw_dftd3,        only : vdw_dftd3
 use m_atm2fft,          only : atm2fft
 use m_mklocl,           only : mklocl_recipspace
 use m_mkcore,           only : mkcore, mkcore_alt

 implicit none

 private
!!***

 public :: stress
!!***

contains
!!***

!!****f* ABINIT/stress
!!
!! NAME
!! stress
!!
!! FUNCTION
!! Compute the stress tensor
!! strten(i,j) = (1/ucvol)*d(Etot)/(d(eps(i,j)))
!! where Etot is energy per unit cell, ucvol is the unstrained unit cell
!! volume, r(i,iat) is the ith position of atom iat,
!! and eps(i,j) is an infinitesimal strain which maps each
!! point r to r(i) -> r(i) + Sum(j) [eps(i,j)*r(j)].
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!! berryopt    =  4/14: electric field is on -> add the contribution of the
!!                      -ebar_i p_i - Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!     = 6/16, or 7/17: electric displacement field is on  -> add the contribution of the
!!                      Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!   from Etot(npw) data (at fixed geometry), used for making
!!   Pulay correction to stress tensor (hartree).  Should be <=0.
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  eei=local pseudopotential part of Etot (hartree)
!!  efield = cartesian coordinates of the electric field in atomic units
!!  ehart=Hartree energy (hartree)
!!  eii=pseudoion core correction energy part of Etot (hartree)
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  ixc = choice of exchange-correlation functional
!!  kinstr(6)=kinetic energy part of stress tensor
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  mqgrid=dimensioned number of q grid points for local psp spline
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  natom=number of atoms in cell
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  nlstr(6)=nonlocal part of stress tensor
!!  nspden=number of spin-density components
!!  nsym=number of symmetries in space group
!!  ntypat=number of types of atoms
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  prtvol=integer controlling volume of printed output
!!  qgrid(mqgrid)=q point array for local psp spline fits
!!  red_efieldbar(3) = efield in reduced units relative to reciprocal lattice
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  strsxc(6)=xc correction to stress
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  typat(natom)=type integer for each atom in cell
!!  usefock=1 if fock operator is used; 0 otherwise.
!!  usekden=1 is kinetic energy density has to be taken into account, 0 otherwise
!!  usepaw= 0 for non paw calculation; =1 for paw calculation
!!  vdw_tol= Van der Waals tolerance
!!  vdw_tol_3bt= Van der Waals tolerance on the 3-body term (only effective
!!               vdw_xc=6)
!!  vdw_xc= Van der Waals correction flag
!!  vlspl(mqgrid,2,ntypat)=local psp spline
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  vxctau(nfft,nspden,4*usekden)=(only for meta-GGA) derivative of XC energy density
!!                                wrt kinetic energy density (depsxcdtau)
!!  vxc_hf(nfft,nspden)=exchange-correlation potential (hartree) in real space for Hartree-Fock corrections
!!  xccc1d(n1xccc*(1-usepaw),6,ntypat)=1D core charge function and five derivatives,
!!                          for each type of atom, from psp (used in Norm-conserving only)
!!  xccc3d(n3xccc)=3D core electron density for XC core correction, bohr^-3
!!  xcctau3d(n3xccc*usekden)=(only for meta-GGA): 3D core electron kinetic energy density for XC core correction
!!  xcccrc(ntypat)=XC core correction cutoff radius (bohr) for each atom type
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  zion(ntypat)=valence charge of each type of atom
!!  znucl(ntypat)=atomic number of atom type
!!
!! OUTPUT
!!  strten(6)=components of the stress tensor (hartree/bohr^3) for the
!!    6 unique components of this symmetric 3x3 tensor:
!!    Given in order (1,1), (2,2), (3,3), (3,2), (3,1), (2,1).
!!    The diagonal components of the returned stress tensor are
!!    CORRECTED for the Pulay stress.
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!
!! NOTES
!! * Concerning the stress tensor:
!!   See O. H. Nielsen and R. M. Martin, PRB 32, 3792 (1985) [[cite:Nielsen1985a]].
!!   Note that first term in equation (2) should have minus sign
!!   (for kinetic energy contribution to stress tensor).
!!   Normalizations in this code differ somewhat from those employed
!!   by Nielsen and Martin.
!!   For the stress tensor contribution from the nonlocal Kleinman-Bylander
!!   separable pseudopotential, see D. M. Bylander, L. Kleinman, and
!!   S. Lee, PRB 42, 1394 (1990) [[cite:Bylander1990]].
!!   Again normalization conventions differ somewhat.
!!   See Doug Allan s notes starting page 795 (13 Jan 1992).
!! * This subroutine calls different subroutines to compute the stress
!!   tensor contributions from the following parts of the total energy:
!!   (1) kinetic energy, (2) exchange-correlation energy,
!!   (3) Hartree energy, (4) local pseudopotential energy,
!!   (5) pseudoion core correction energy, (6) nonlocal pseudopotential energy,
!!   (7) Ewald energy.
!!
!! PARENTS
!!      m_forstr
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

 subroutine stress(atindx1,berryopt,dtefield,eei,efield,ehart,eii,fock,gsqcut,ixc,kinstr,&
&                  mgfft,mpi_enreg,mqgrid,n1xccc,n3xccc,natom,nattyp,&
&                  nfft,ngfft,nlstr,nspden,nsym,ntypat,psps,pawrad,pawtab,ph1d,&
&                  prtvol,qgrid,red_efieldbar,rhog,rprimd,strten,strsxc,symrec,&
&                  typat,usefock,usekden,usepaw,vdw_tol,vdw_tol_3bt,vdw_xc,&
&                  vlspl,vxc,vxctau,vxc_hf,xccc1d,xccc3d,xcctau3d,xcccrc,xred,zion,znucl,qvpotzero,&
&                  electronpositron) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: berryopt,ixc,mgfft,mqgrid,n1xccc,n3xccc,natom,nfft,nspden
 integer,intent(in) :: nsym,ntypat,prtvol,usefock,usekden,usepaw,vdw_xc
 real(dp),intent(in) :: eei,ehart,eii,gsqcut,vdw_tol,vdw_tol_3bt,qvpotzero
 type(efield_type),intent(in) :: dtefield
 type(pseudopotential_type),intent(in) :: psps
 type(electronpositron_type),pointer,optional :: electronpositron
 type(MPI_type),intent(in) :: mpi_enreg
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(natom),nattyp(ntypat),ngfft(18),symrec(3,3,nsym)
 integer,intent(in) :: typat(natom)
 real(dp),intent(in) :: efield(3),kinstr(6),nlstr(6)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*natom),qgrid(mqgrid)
 real(dp),intent(in) :: red_efieldbar(3),rhog(2,nfft),strsxc(6)
 real(dp),intent(in) :: vlspl(mqgrid,2,ntypat),vxc(nfft,nspden),vxctau(nfft,nspden,4*usekden)
 real(dp),allocatable,intent(in) :: vxc_hf(:,:)
 real(dp),intent(in) :: xccc1d(n1xccc*(1-usepaw),6,ntypat),xcccrc(ntypat)
 real(dp),intent(in) :: xred(3,natom),zion(ntypat),znucl(ntypat)
 real(dp),intent(inout) :: xccc3d(n3xccc),xcctau3d(n3xccc*usekden),rprimd(3,3)
 real(dp),intent(out) :: strten(6)
 type(pawrad_type),intent(in) :: pawrad(ntypat*usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: coredens_method,coretau_method,iatom,icoulomb,idir,ii,ipositron,mu,nkpt=1
 integer :: optatm,optdyfr,opteltfr,opt_hybr,optgr,option,optn,optn2,optstr,optv,sdir,vloc_method
 real(dp),parameter :: tol=1.0d-15
 real(dp) :: e_dum,dum_rcut=zero,strsii,ucvol,vol_element
 character(len=500) :: message
 logical :: calc_epaw3_stress, efield_flag
!arrays
 integer :: qprtrb_dum(3),icutcoul=3
 real(dp) :: corstr(6),dumstr(6),ep3(3),epaws3red(6),ewestr(6),gmet(3,3),vcutgeo(3)
 real(dp) :: gprimd(3,3),harstr(6),lpsstr(6),rmet(3,3),taustr(6),tsec(2),uncorr(3)
 real(dp) :: vdwstr(6),vprtrb_dum(2)
 real(dp) :: Maxstr(6),ModE !Maxwell-stress constribution, and magnitude of efield
 real(dp) :: dummy_in(0)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0),dummy_out7(0)
 real(dp),allocatable :: dummy(:),dyfr_dum(:,:,:),gr_dum(:,:),rhog_ep(:,:),v_dum(:)
 real(dp),allocatable :: vxctotg(:,:)
 character(len=10) :: EPName(1:2)=(/"Electronic","Positronic"/)

! *************************************************************************

 call timab(37,1,tsec)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 opt_hybr=0;if (allocated(vxc_hf)) opt_hybr=1
 icoulomb=0 ! not yet compatible with icoulomb

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

!Determine by which method the local ionic potential and/or the pseudo core charge density
! contributions have to be computed
!Local ionic potential:
! Method 1: PAW
! Method 2: Norm-conserving PP, icoulomb>0, wavelets
 vloc_method=1;if (usepaw==0) vloc_method=2
 if (psps%usewvl==1) vloc_method=2
!Pseudo core charge density:
! Method 1: PAW, nc_xccc_gspace
! Method 2: Norm-conserving PP, wavelets
 coredens_method=1;if (usepaw==0) coredens_method=2
 if (psps%nc_xccc_gspace==1) coredens_method=1
 if (psps%nc_xccc_gspace==0) coredens_method=2
 if (psps%usewvl==1) coredens_method=2
 coretau_method=0
 if (usekden==1.and.psps%usepaw==1) then
   coretau_method=1;if (psps%nc_xccc_gspace==0) coretau_method=2
 end if

!Local ionic potential and/or pseudo core charge by method 1
 if (vloc_method==1.or.coredens_method==1.or.coretau_method==1) then
   call timab(551,1,tsec)
!  Compute Vxc in reciprocal space
   if (coredens_method==1.and.n3xccc>0) then
     ABI_MALLOC(v_dum,(nfft))
     ABI_MALLOC(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,1,ngfft,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_FREE(v_dum)
   else
     ABI_MALLOC(vxctotg,(0,0))
   end if
!  Compute contribution to stresses from Vloc and/or pseudo core density
   optv=0;if (vloc_method==1) optv=1
   optn=0;if (coredens_method==1) optn=n3xccc/nfft
   optatm=0;optdyfr=0;opteltfr=0;optgr=0;optstr=1;optn2=1
   if (vloc_method==1.or.coredens_method==1) then
     call atm2fft(atindx1,dummy_out1,dummy_out2,dummy_out3,dummy_out4,&
&     dummy_out5,dummy_in,gmet,gprimd,dummy_out6,dummy_out7,gsqcut,&
&     mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
&     psps,pawtab,ph1d,qgrid,qprtrb_dum,dum_rcut,rhog,rprimd,corstr,lpsstr,ucvol,usepaw,vxctotg,vxctotg,vxctotg,vprtrb_dum,vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
   end if
   if (n3xccc==0.and.coredens_method==1) corstr=zero
   ABI_FREE(vxctotg)
   if (usekden==1.and.coretau_method==1..and.n3xccc>0) then
!    Compute contribution to stresses from pseudo kinetic energy core density
     optv=0;optn=1;optn2=4
     ABI_MALLOC(v_dum,(nfft))
     ABI_MALLOC(vxctotg,(2,nfft))
     v_dum(:)=vxctau(:,1,1);if (nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxctau(:,2,1))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,1,ngfft,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_FREE(v_dum)
     call atm2fft(atindx1,dummy_out1,dummy_out2,dummy_out3,dummy_out4,&
&     dummy_out5,dummy_in,gmet,gprimd,dummy_out6,dummy_out7,gsqcut,&
&     mgfft,mqgrid,natom,nattyp,nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
&     psps,pawtab,ph1d,qgrid,qprtrb_dum,dum_rcut,rhog,rprimd,taustr,dumstr,ucvol,usepaw,vxctotg,vxctotg,vxctotg,vprtrb_dum,vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
     corstr(1:6)=corstr(1:6)+taustr(1:6)
   end if
   call timab(551,2,tsec)
 end if

!Local ionic potential by method 2
 if (vloc_method==2) then
   option=3
   ABI_MALLOC(dyfr_dum,(3,3,natom))
   ABI_MALLOC(gr_dum,(3,natom))
   ABI_MALLOC(v_dum,(nfft))
   call mklocl_recipspace(dyfr_dum,eei,gmet,gprimd,gr_dum,gsqcut,icutcoul,lpsstr,mgfft,&
&   mpi_enreg,mqgrid,natom,nattyp,nfft,ngfft,nkpt,ntypat,option,ph1d,qgrid,&
&   qprtrb_dum,dum_rcut,rhog,rprimd,ucvol,vcutgeo,vlspl,vprtrb_dum,v_dum)
   ABI_FREE(dyfr_dum)
   ABI_FREE(gr_dum)
   ABI_FREE(v_dum)
 end if

!Pseudo core electron density by method 2
 if (coredens_method==2.or.coretau_method==2) then
   if (n1xccc/=0) then
     call timab(55,1,tsec)
     option=3
     ABI_MALLOC(dyfr_dum,(3,3,natom))
     ABI_MALLOC(gr_dum,(3,natom))
     ABI_MALLOC(v_dum,(nfft))
     if (coredens_method==2) then
       if (psps%usewvl==0.and.usepaw==0.and.icoulomb==0) then
         if(opt_hybr==0) then
           call mkcore(corstr,dyfr_dum,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,ngfft(1),&
&           n1xccc,ngfft(2),ngfft(3),option,rprimd,typat,ucvol,vxc,&
&           xcccrc,xccc1d,xccc3d,xred)
         else
           call mkcore(corstr,dyfr_dum,gr_dum,mpi_enreg,natom,nfft,nspden,ntypat,ngfft(1),&
&           n1xccc,ngfft(2),ngfft(3),option,rprimd,typat,ucvol,vxc_hf,&
&           xcccrc,xccc1d,xccc3d,xred)
         end if
       else if (psps%usewvl==0.and.(usepaw==1.or.icoulomb==1)) then
         call mkcore_alt(atindx1,corstr,dyfr_dum,gr_dum,icoulomb,mpi_enreg,natom,nfft,&
&         nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&         ucvol,vxc,xcccrc,xccc1d,xccc3d,xred,pawrad,pawtab,usepaw)
       end if
     end if
     if (usekden==1.and.coretau_method==2) then
       call mkcore_alt(atindx1,taustr,dyfr_dum,gr_dum,icoulomb,mpi_enreg,natom,nfft,&
&       nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&       ucvol,vxctau(:,:,1),xcccrc,xccc1d,xcctau3d,xred,pawrad,pawtab,usepaw,&
&       usekden=.true.)

     end if
     ABI_FREE(dyfr_dum)
     ABI_FREE(gr_dum)
     ABI_FREE(v_dum)
     call timab(55,2,tsec)
   else
     corstr(:)=zero
   end if
 end if

!=======================================================================
!======================= Hartree energy contribution ===================
!=======================================================================

 call strhar(ehart,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,rprimd)

!=======================================================================
!======================= Ewald contribution ============================
!=======================================================================

 call timab(38,1,tsec)
 call ewald2(gmet,natom,ntypat,rmet,rprimd,ewestr,typat,ucvol,xred,zion)

!=======================================================================
!================== VdW DFT-D contribution ============================
!=======================================================================

 if (vdw_xc==5) then
   call vdw_dftd2(e_dum,ixc,natom,ntypat,0,typat,rprimd,vdw_tol,&
&   xred,znucl,str_vdw_dftd2=vdwstr)
 elseif (vdw_xc==6.or.vdw_xc==7) then
   call vdw_dftd3(e_dum,ixc,natom,ntypat,0,typat,rprimd,vdw_xc,&
&   vdw_tol,vdw_tol_3bt,xred,znucl,str_vdw_dftd3=vdwstr)
 end if

 call timab(38,2,tsec)

!HONG  no Berry phase contribution if using reduced ebar or d according to
!HONG  PRL 89, 117602 (2002) [[cite:Souza2002]]
!HONG  Nature Physics: M. Stengel et.al. (2009)) [[cite:Stengel1999]]
!=======================================================================
!=================== Berry phase contribution ==========================
!=======================================================================

!if (berryopt==4) then
!berrystr_tmp(:,:) = zero
!Diagonal:
!do mu = 1, 3
!do ii = 1, 3
!berrystr_tmp(mu,mu) = berrystr_tmp(mu,mu) - &
!&       efield(mu)*rprimd(mu,ii)*(pel(ii) + pion(ii))/ucvol
!end do
!end do
!Off-diagonal (symmetrized before adding it to strten):
!do ii = 1, 3
!berrystr_tmp(3,2) = berrystr_tmp(3,2) &
!&     - efield(3)*rprimd(2,ii)*(pel(ii) + pion(ii))/ucvol
!berrystr_tmp(2,3) = berrystr_tmp(2,3) &
!&     - efield(2)*rprimd(3,ii)*(pel(ii) + pion(ii))/ucvol
!berrystr_tmp(3,1) = berrystr_tmp(3,1) &
!&     - efield(3)*rprimd(1,ii)*(pel(ii) + pion(ii))/ucvol
!berrystr_tmp(1,3) = berrystr_tmp(1,3) &
!&     - efield(1)*rprimd(3,ii)*(pel(ii) + pion(ii))/ucvol
!berrystr_tmp(2,1) = berrystr_tmp(2,1) &
!&     - efield(2)*rprimd(1,ii)*(pel(ii) + pion(ii))/ucvol
!berrystr_tmp(1,2) = berrystr_tmp(1,2) &
!&     - efield(1)*rprimd(2,ii)*(pel(ii) + pion(ii))/ucvol
!end do
!berrystr(1) = berrystr_tmp(1,1)
!berrystr(2) = berrystr_tmp(2,2)
!berrystr(3) = berrystr_tmp(3,3)
!berrystr(4) = (berrystr_tmp(3,2) + berrystr_tmp(2,3))/two
!berrystr(5) = (berrystr_tmp(3,1) + berrystr_tmp(1,3))/two
!berrystr(6) = (berrystr_tmp(2,1) + berrystr_tmp(1,2))/two
!end if

!=======================================================================
!================= Other (trivial) contributions =======================
!=======================================================================

!Nonlocal part of stress has already been computed
!(in forstrnps(norm-conserving) or pawgrnl(PAW))

!Kinetic part of stress has already been computed
!(in forstrnps)

!XC part of stress tensor has already been computed in "strsxc"

!ii part of stress (diagonal) is trivial!
 strsii=-eii/ucvol
!qvpotzero is non zero, only when usepotzero=1
 strsii=strsii+qvpotzero/ucvol

!======================================================================
!HONG  Maxwell stress when electric/displacement field is non-zero=====
!======================================================================
 efield_flag = (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or. &
& berryopt==14 .or. berryopt==16 .or. berryopt==17)
 calc_epaw3_stress = (efield_flag .and. usepaw == 1)
 if ( efield_flag ) then
   ModE=dot_product(efield,efield)
   do ii=1,3
     Maxstr(ii)=two*efield(ii)*efield(ii)-ModE
   end do
   Maxstr(4)=two*efield(3)*efield(2)
   Maxstr(5)=two*efield(3)*efield(1)
   Maxstr(6)=two*efield(2)*efield(1)
!  Converting to units of Ha/Bohr^3
!  Maxstr(:)=Maxstr(:)*e_Cb*Bohr_Ang*1.0d-10/(Ha_J*8.0d0*pi)

   Maxstr(:)=Maxstr(:)*eps0*Ha_J*Bohr_Ang*1.0d-10/(8.0d0*pi*e_Cb**2)

   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of Maxwell stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   ' Maxstr(1 1)=',Maxstr(1),' Maxstr(3 2)=',Maxstr(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   ' Maxstr(2 2)=',Maxstr(2),' Maxstr(3 1)=',Maxstr(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   ' Maxstr(3 3)=',Maxstr(3),' Maxstr(2 1)=',Maxstr(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')

 end if

! compute additional F3-type stress due to projectors for electric field with PAW
 if ( efield_flag .and. calc_epaw3_stress ) then
   do sdir = 1, 6
     ep3(:) = zero
     do idir = 1, 3
       vol_element=one/(ucvol*dtefield%nstr(idir)*dtefield%nkstr(idir))
       do iatom = 1, natom
         ep3(idir) = ep3(idir) + vol_element*dtefield%epaws3(iatom,idir,sdir)
       end do ! end loop over atoms
     end do ! end loop over idir (components of P)
! note no appearance of ucvol here unlike in forces, stress definition includes
! division by ucvol, which cancels the factor in -ucvol e . p
     epaws3red(sdir) = -dot_product(red_efieldbar(1:3),ep3(1:3))
   end do

!   write(message, '(a,a)' )ch10,&
!&   ' Cartesian components of PAW sigma_3 stress tensor (hartree/bohr^3)'
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,  message,'COLL')
!   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
!&   ' epaws3red(1 1)=',epaws3red(1),' epaws3red(3 2)=',epaws3red(4)
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,  message,'COLL')
!   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
!&   ' epaws3red(2 2)=',epaws3red(2),' epaws3red(3 1)=',epaws3red(5)
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,  message,'COLL')
!   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
!&   ' epaws3red(3 3)=',epaws3red(3),' epaws3red(2 1)=',epaws3red(6)
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,  message,'COLL')
!   write(message, '(a)' ) ' '
!   call wrtout(ab_out,message,'COLL')
!   call wrtout(std_out,  message,'COLL')

 end if

!=======================================================================
!===== Assemble the various contributions to the stress tensor =========
!=======================================================================
!In cartesian coordinates (symmetric storage)

 strten(:)=kinstr(:)+ewestr(:)+corstr(:)+strsxc(:)+harstr(:)+lpsstr(:)+nlstr(:)

 if (usefock==1 .and. associated(fock)) then
   if (fock%fock_common%optstr) then
     strten(:)=strten(:)+fock%fock_common%stress(:)
   end if
 end if

!Add contributions for constant E or D calculation.
 if ( efield_flag ) then
   strten(:)=strten(:)+Maxstr(:)
   if ( calc_epaw3_stress ) strten(:) = strten(:) + epaws3red(:)
 end if
 if (vdw_xc>=5.and.vdw_xc<=7) strten(:)=strten(:)+vdwstr(:)

!Additional stuff for electron-positron
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     if (allocated(electronpositron%stress_ep)) ipositron=electronpositron_calctype(electronpositron)
   end if
 end if
 if (abs(ipositron)==1) then
   strten(:)=strten(:)-harstr(:)-ewestr(:)-corstr(:)-lpsstr(:)
   harstr(:)=zero;ewestr(:)=zero;corstr(:)=zero;strsii=zero
   lpsstr(:)=-lpsstr(:);lpsstr(1:3)=lpsstr(1:3)-two*eei/ucvol
   strten(:)=strten(:)+lpsstr(:)
   if (vdw_xc>=5.and.vdw_xc<=7) strten(:)=strten(:)-vdwstr(:)
   if (vdw_xc>=5.and.vdw_xc<=7) vdwstr(:)=zero
 end if
 if (abs(ipositron)==2) then
   ABI_MALLOC(rhog_ep,(2,nfft))
   ABI_MALLOC(dummy,(6))
   call fourdp(1,rhog_ep,electronpositron%rhor_ep,-1,mpi_enreg,nfft,1,ngfft,0)
   rhog_ep=-rhog_ep
   call strhar(electronpositron%e_hartree,gsqcut,dummy,mpi_enreg,nfft,ngfft,rhog_ep,rprimd)
   strten(:)=strten(:)+dummy(:);harstr(:)=harstr(:)+dummy(:)
   ABI_FREE(rhog_ep)
   ABI_FREE(dummy)
 end if
 if (ipositron>0) strten(:)=strten(:)+electronpositron%stress_ep(:)

!Symmetrize resulting tensor if nsym>1
 if (nsym>1) then
   call stresssym(gprimd,nsym,strten,symrec)
 end if

!Set to zero very small values of stress
 do mu=1,6
   if (abs(strten(mu))<tol) strten(mu)=zero
 end do

!Include diagonal terms, save uncorrected stress for output
 do mu=1,3
   uncorr(mu)=strten(mu)+strsii
   strten(mu)=uncorr(mu)
 end do

!=======================================================================
!================ Print out info about stress tensor ===================
!=======================================================================
 if (prtvol>=10.and.ipositron>=0) then
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of hartree stress is',harstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of loc psp stress is',lpsstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,&
&     ' of kinetic stress is',kinstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of nonlocal ps stress is',nlstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,' of     core xc stress is',corstr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' )&
&     ' stress: component',mu,&
&     ' of Ewald energ stress is',ewestr(mu)
     call wrtout(std_out,message,'COLL')
   end do
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   do mu=1,6
     write(message, '(a,i5,a,1p,e22.12)' ) &
&     ' stress: component',mu,' of xc stress is',strsxc(mu)
     call wrtout(std_out,message,'COLL')
   end do
   if (vdw_xc>=5.and.vdw_xc<=7) then
     write(message, '(a)' ) ' '
     call wrtout(std_out,message,'COLL')
     do mu=1,6
       write(message, '(a,i5,a,1p,e22.12)' )&
&       ' stress: component',mu,&
&       ' of VdW DFT-D stress is',vdwstr(mu)
       call wrtout(std_out,message,'COLL')
     end do
   end if
   write(message, '(a)' ) ' '
   call wrtout(std_out,message,'COLL')
   write(message, '(a,1p,e22.12)' ) &
&   ' stress: ii (diagonal) part is',strsii
   call wrtout(std_out,message,'COLL')
   if (berryopt==4 .or. berryopt==6 .or. berryopt==7 .or.  &
&   berryopt==14 .or. berryopt==16 .or. berryopt==17) then  !!HONG
     write(message, '(a)' ) ' '
     call wrtout(std_out,message,'COLL')
     do mu = 1, 6
       write(message, '(a,i2,a,1p,e22.12)' )&
&       ' stress: component',mu,' of Maxwell stress is',&
&       Maxstr(mu)
       call wrtout(std_out,message,'COLL')
     end do
   end if
   if (ipositron/=0) then
     write(message, '(a)' ) ' '
     call wrtout(std_out,message,'COLL')
     do mu=1,6
       write(message, '(a,i5,3a,1p,e22.12)' ) &
&       ' stress: component',mu,' of ',EPName(abs(ipositron)), &
&       ' stress is',electronpositron%stress_ep(mu)
       call wrtout(std_out,message,'COLL')
     end do
   end if

 end if ! prtvol
 if (ipositron>=0) then
   write(message, '(a,a)' )ch10,&
&   ' Cartesian components of stress tensor (hartree/bohr^3)'
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(1 1)=',strten(1),'  sigma(3 2)=',strten(4)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(2 2)=',strten(2),'  sigma(3 1)=',strten(5)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a,1p,e16.8,a,1p,e16.8)' ) &
&   '  sigma(3 3)=',strten(3),'  sigma(2 1)=',strten(6)
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
   write(message, '(a)' ) ' '
   call wrtout(ab_out,message,'COLL')
   call wrtout(std_out,  message,'COLL')
 end if
 call timab(37,2,tsec)

end subroutine stress
!!***

!!****f* ABINIT/strhar
!!
!! NAME
!! strhar
!!
!! FUNCTION
!! Compute Hartree energy contribution to stress tensor (Cartesian coordinates).
!!
!! INPUTS
!!  ehart=Hartree energy (hartree)
!!  gsqcut=cutoff value on $G^2$ for (large) sphere inside fft box.
!!  $gsqcut=(boxcut^2)*ecut/(2._dp*(\pi^2))$
!!  mpi_enreg=information about MPI parallelization
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhog(2,nfft)= optional argument: Fourier transform of a second charge density (bohr^-3)
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!
!! OUTPUT
!!  harstr(6)=components of Hartree part of stress tensor
!!   (Cartesian coordinates, symmetric tensor) in hartree/bohr^3
!!   Definition of symmetric tensor storage: store 6 unique components
!!   in the order 11, 22, 33, 32, 31, 21 (suggested by Xavier Gonze).
!!
!! PARENTS
!!      m_stress
!!
!! CHILDREN
!!      metric,ptabs_fourdp,timab,xmpi_sum
!!
!! SOURCE

subroutine strhar(ehart,gsqcut,harstr,mpi_enreg,nfft,ngfft,rhog,rprimd,&
&                 rhog2) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft
 real(dp),intent(in) :: ehart,gsqcut
 type(MPI_type),intent(in) :: mpi_enreg
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rprimd(3,3),rhog(2,nfft)
 real(dp),intent(in),optional :: rhog2(2,nfft)
 real(dp),intent(out) :: harstr(6)

!Local variables-------------------------------
!scalars
 integer,parameter :: im=2,re=1
 integer :: i1,i2,i3,id1,id2,id3,ierr,ig1,ig2,ig3,ii,irho2,me_fft,n1,n2,n3,nproc_fft
 real(dp) :: cutoff,gsquar,rhogsq,tolfix=1.000000001_dp,ucvol
!arrays
 real(dp) :: gcart(3),gmet(3,3),gprimd(3,3),rmet(3,3),tsec(2)
 integer, ABI_CONTIGUOUS pointer :: fftn2_distrib(:),ffti2_local(:)
 integer, ABI_CONTIGUOUS pointer :: fftn3_distrib(:),ffti3_local(:)

! *************************************************************************

 call timab(568,1,tsec)

 harstr(:)=zero
!ehtest=0.0_dp (used for testing)

 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

 irho2=0;if (present(rhog2)) irho2=1

!Conduct looping over all fft grid points to find G vecs inside gsqcut
!Include G**2 on surface of cutoff sphere as well as inside:
 cutoff=gsqcut*tolfix
 n1=ngfft(1) ; n2=ngfft(2) ; n3=ngfft(3)
 me_fft=ngfft(11)
 nproc_fft=ngfft(10)
 id1=n1/2+2
 id2=n2/2+2
 id3=n3/2+2
 ii=0

 ! Get the distrib associated with this fft_grid
 call ptabs_fourdp(mpi_enreg,n2,n3,fftn2_distrib,ffti2_local,fftn3_distrib,ffti3_local)

 do i3=1,n3
   ig3=i3-(i3/id3)*n3-1
   do i2=1,n2
     ig2=i2-(i2/id2)*n2-1
     if (fftn2_distrib(i2)==me_fft) then
       do i1=1,n1
         ig1=i1-(i1/id1)*n1-1
!        ii=ii+1
         ii=i1+n1*(ffti2_local(i2)-1+(n2/nproc_fft)*(i3-1))
!        **     GET RID OF THIS IF STATEMENT LATER for speed if needed
!        Avoid G=0:
!        if (ii>1) then
         if (ig1==0 .and. ig2==0 .and. ig3==0) cycle
!        Compute cartesian components of G
         gcart(1)=gprimd(1,1)*dble(ig1)+gprimd(1,2)*dble(ig2)+gprimd(1,3)*dble(ig3)
         gcart(2)=gprimd(2,1)*dble(ig1)+gprimd(2,2)*dble(ig2)+gprimd(2,3)*dble(ig3)
         gcart(3)=gprimd(3,1)*dble(ig1)+gprimd(3,2)*dble(ig2)+gprimd(3,3)*dble(ig3)
!        Compute |G|^2
         gsquar=gcart(1)**2+gcart(2)**2+gcart(3)**2

!        Keep only G**2 inside larger cutoff (not sure this is needed):
         if (gsquar<=cutoff) then
!          take |rho(G)|^2 for complex rhog
           if (irho2==0) then
             rhogsq=rhog(re,ii)**2+rhog(im,ii)**2
           else
             rhogsq=rhog(re,ii)*rhog2(re,ii)+rhog(im,ii)*rhog2(im,ii)
           end if
           harstr(1)=harstr(1)+(rhogsq/gsquar**2)*gcart(1)*gcart(1)
           harstr(2)=harstr(2)+(rhogsq/gsquar**2)*gcart(2)*gcart(2)
           harstr(3)=harstr(3)+(rhogsq/gsquar**2)*gcart(3)*gcart(3)
           harstr(4)=harstr(4)+(rhogsq/gsquar**2)*gcart(3)*gcart(2)
           harstr(5)=harstr(5)+(rhogsq/gsquar**2)*gcart(3)*gcart(1)
           harstr(6)=harstr(6)+(rhogsq/gsquar**2)*gcart(2)*gcart(1)
         end if
!        end if
       end do
     end if
   end do
 end do

!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
#ifdef FC_IBM
 write(std_out,*)' strhar : before mpi_comm, harstr=',harstr
#endif

!Init mpi_comm
 if(mpi_enreg%nproc_fft>1)then
   call timab(48,1,tsec)
   call xmpi_sum(harstr,mpi_enreg%comm_fft ,ierr)
   call timab(48,2,tsec)
 end if

#ifdef FC_IBM
!DO not remove : seems needed to avoid problem with pathscale compiler, in parallel
 write(std_out,*)' strhar : after mpi_comm, harstr=',harstr
 write(std_out,*)' strhar : ehart,ucvol=',ehart,ucvol
#endif

!Normalize and add term -ehart/ucvol on diagonal
 harstr(1)=harstr(1)/pi-ehart/ucvol
 harstr(2)=harstr(2)/pi-ehart/ucvol
 harstr(3)=harstr(3)/pi-ehart/ucvol
 harstr(4)=harstr(4)/pi
 harstr(5)=harstr(5)/pi
 harstr(6)=harstr(6)/pi

 call timab(568,2,tsec)

end subroutine strhar
!!***

end module m_stress
!!***

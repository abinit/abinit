!!****m* ABINIT/m_forces
!! NAME
!!  m_forces
!!
!! FUNCTION
!!
!!
!! COPYRIGHT
!! Copyright (C) 1998-2020 ABINIT group (DCA, XG, GMR, FJ, MM, MT, SCE)
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

module m_forces

 use defs_basis
 use defs_wvltypes
 use m_abicore
 use m_efield
 use m_errors
 use m_atomdata
 use m_dtset

 use defs_datatypes,     only : pseudopotential_type
 use defs_abitypes,      only : MPI_type
 use m_time,             only : timab
 use m_geometry,         only : fred2fcart, metric, xred2xcart
 use m_fock,             only : fock_type
 use m_pawrad,           only : pawrad_type
 use m_pawtab,           only : pawtab_type
 use m_electronpositron, only : electronpositron_type,electronpositron_calctype
 use libxc_functionals,  only : libxc_functionals_is_hybrid
 use m_fft,              only : zerosym, fourdp
 use m_cgtools,          only : mean_fftr
 use m_mpinfo,           only : pre_gather, pre_scatter
 use m_atm2fft,          only : atm2fft
 use m_mklocl,           only : mklocl
 use m_predtk,           only : prtxvf
 use m_xchybrid,         only : xchybrid_ncpp_cc
 use m_mkcore,           only : mkcore, mkcore_alt
 use m_mkcore_wvl,       only : mkcore_wvl

 implicit none

 private
!!***

 public :: forces
 public :: fresid
!!***

contains
!!***

!!****f* ABINIT/forces
!! NAME
!! forces
!!
!! FUNCTION
!! Assemble gradients of various total energy terms with respect
!! to reduced coordinates, including possible symmetrization,
!! in order to produce forces.
!!
!!     fcart(i,iat) = d(Etot)/(d(r(i,iat)))
!!
!! INPUTS
!!  atindx1(natom)=index table for atoms, inverse of atindx
!!  dtefield <type(efield_type)> = variables related to Berry phase
!!  dtset <type(dataset_type)>=all input variables in this dataset
!! berryopt    =  4/14: electric field is on -> add the contribution of the
!!                      -ebar_i p_i - Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!     = 6/16, or 7/17: electric displacement field is on  -> add the contribution of the
!!                      Omega/(8*pi) (g^{-1})_ij ebar_i ebar_j  terms to the total energy
!!   | efield = cartesian coordinates of the electric field in atomic units
!!   | dfield = cartesian coordinates of the electric displacement field in atomic units
!!   | iatfix(3,natom)=1 for frozen atom along specified direction, 0 for unfrozen
!!   | ionmov=governs the movement of atoms (see help file)
!!   | densfor_pred=governs the mixed electronic-atomic part of the preconditioner
!!   | natom=number of atoms in cell
!!   | nconeq=number of atomic constraint equations
!!   | nspden=number of spin-density components
!!   | nsym=number of symmetries in space group
!!   | prtvol=integer controlling volume of printed output
!!   | typat(natom)=type integer for each atom in cell
!!   | wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  fock <type(fock_type)>= quantities to calculate Fock exact exchange
!!  grchempottn(3,natom)=d(E_chemical potential)/d(xred) (hartree)
!!  grcondft(3,natom)=d(E_constrainedDFT)/d(xred) (hartree)
!!  grewtn(3,natom)=d(Ewald)/d(xred) (hartree)
!!  grnl(3*natom)=gradients of Etot due to nonlocal contributions
!!  grvdw(3,ngrvdw)=gradients of energy due to Van der Waals DFT-D dispersion (hartree)
!!  gsqcut=cutoff value on G**2 for (large) sphere inside FFT box.
!!                       gsqcut=(boxcut**2)*ecut/(2._dp*(Pi**2)
!!  indsym(4,nsym,natom)=indirect indexing array for atom labels
!!  mgfft=maximum size of 1D FFTs
!!  mpi_enreg=information about MPI parallelization
!!  n1xccc=dimension of xccc1d ; 0 if no XC core correction is used
!!  n3xccc=dimension of the xccc3d array (0 or nfft).
!!  nattyp(ntypat)=number of atoms of each type
!!  nfft=(effective) number of FFT grid points (for this processor)
!!  ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!!  ngrvdw=size of grvdw(:,:); can be 0 or natom according to dtset%vdw_xc
!!  ntypat=number of types of atoms
!!  pawrad(ntypat*usepaw) <type(pawrad_type)>=paw radial mesh and related data
!!  pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!!  ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase (structure factor) array
!!  psps <type(pseudopotential_type)>=variables related to pseudopotentials
!!  rhog(2,nfft)=Fourier transform of charge density (bohr^-3)
!!  rhor(nfft,nspden)=array for electron density in electrons/bohr**3
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  symrec(3,3,nsym)=symmetries in reciprocal space, reduced coordinates
!!  usefock=1 if fock operator is used; 0 otherwise.
!!  usekden= 1 is kinetic energy density has to be computed, 0 otherwise
!!  vresid(nfft,nspden)=potential residual (if non-collinear magn., only trace of it)
!!  vxc(nfft,nspden)=exchange-correlation potential (hartree) in real space
!!  vxctau(nfft,nspden,4*usekden)=(only for meta-GGA): derivative of XC energy density
!!                                wrt kinetic energy density (depsxcdtau)
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!  xred_old(3,natom)=previous reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximal absolute value of changes in the components of
!!         force between the input and the output.
!!  favg(3)=mean of the forces before correction for translational symmetry
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!  fred(3,natom)=symmetrized grtn = d(etotal)/d(xred)
!!  gresid(3,natom)=forces due to the residual of the density/potential
!!  grhf(3,natom)=Hellman-Feynman derivatives of the total energy
!!  grxc(9+3*natom)=d(Exc)/d(xred) if core charges are used
!!  maxfor=maximal absolute value of the output array force.
!!  synlgr(3,natom)=symmetrized d(enl)/d(xred)
!!
!! SIDE EFFECTS
!!  electronpositron <type(electronpositron_type)>=quantities for the electron-positron annihilation (optional argument)
!!  fcart(3,natom)=forces in cartesian coordinates (Ha/Bohr)
!!    Note : unlike fred, this array has been corrected by enforcing
!!    the translational symmetry, namely that the sum of force
!!    on all atoms is zero.
!!
!! NOTES
!! * Symmetrization of gradients with respect to reduced
!!   coordinates xred is conducted according to the expression
!!   [d(e)/d(t(n,a))]_symmetrized = (1/Nsym) Sum(S) symrec(n,m,S)*
!!                [d(e)/d(t(m,b))]_unsymmetrized
!!   where t(m,b)= (symrel^-1)(m,n)*(t(n,a)-tnons(n)) and tnons
!!   is a possible nonsymmorphic translation.  The label "b" here
!!   refers to the atom which gets rotated into "a" under symmetry "S".
!!   symrel is the symmetry matrix in real space, which is the inverse
!!   transpose of symrec.  symrec is the symmetry matrix in reciprocal
!!   space.  sym_cartesian = R * symrel * R^-1 = G * symrec * G^-1
!!   where the columns of R and G are the dimensional primitive translations
!!   in real and reciprocal space respectively.
!! * Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      etotfor,forstr
!!
!! CHILDREN
!!      atm2fft,constrf,dgemv,fourdp,fred2fcart,fresid,fresidrsp,metric,mkcore
!!      mkcore_alt,mkcore_wvl,mklocl,sygrad,timab,xchybrid_ncpp_cc,zerosym
!!
!! SOURCE

subroutine forces(atindx1,diffor,dtefield,dtset,favg,fcart,fock,&
&                  forold,fred,grchempottn,grcondft,gresid,grewtn,&
&                  grhf,grnl,grvdw,grxc,gsqcut,indsym,&
&                  maxfor,mgfft,mpi_enreg,n1xccc,n3xccc,&
&                  nattyp,nfft,ngfft,ngrvdw,ntypat,&
&                  pawrad,pawtab,ph1d,psps,rhog,rhor,rprimd,symrec,synlgr,usefock,&
&                  vresid,vxc,vxctau,wvl,wvl_den,xred,&
&                  electronpositron) ! optional argument

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,n1xccc,n3xccc,nfft,ngrvdw,ntypat,usefock
 real(dp),intent(in) :: gsqcut
 real(dp),intent(out) :: diffor,maxfor
 type(MPI_type),intent(in) :: mpi_enreg
 type(efield_type),intent(in) :: dtefield
 type(dataset_type),intent(in) :: dtset
 type(electronpositron_type),pointer,optional :: electronpositron
 type(pseudopotential_type),intent(in) :: psps
 type(wvl_internal_type), intent(in) :: wvl
 type(wvl_denspot_type), intent(inout) :: wvl_den
 type(fock_type),pointer, intent(inout) :: fock
!arrays
 integer,intent(in) :: atindx1(dtset%natom),indsym(4,dtset%nsym,dtset%natom)
 integer,intent(in) :: nattyp(ntypat),ngfft(18),symrec(3,3,dtset%nsym)
 real(dp),intent(in) :: grchempottn(3,dtset%natom),grcondft(3,dtset%natom),grewtn(3,dtset%natom)
 real(dp),intent(in) :: grvdw(3,ngrvdw),grnl(3*dtset%natom)
 real(dp),intent(in) :: ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: rhog(2,nfft),rhor(nfft,dtset%nspden)
 real(dp),intent(in) :: vxc(nfft,dtset%nspden),vxctau(nfft,dtset%nspden,4*dtset%usekden)
 real(dp),intent(inout) :: fcart(3,dtset%natom),forold(3,dtset%natom)
 real(dp),intent(inout) :: vresid(nfft,dtset%nspden),xred(3,dtset%natom)
 real(dp),intent(out) :: favg(3),fred(3,dtset%natom),gresid(3,dtset%natom)
 real(dp),intent(out) :: grhf(3,dtset%natom),rprimd(3,3)
 real(dp),intent(inout) :: grxc(3,dtset%natom)
 real(dp),intent(out) :: synlgr(3,dtset%natom)
 type(pawrad_type),intent(in) :: pawrad(ntypat*psps%usepaw)
 type(pawtab_type),intent(in) :: pawtab(ntypat*psps%usepaw)

!Local variables-------------------------------
!scalars
 integer :: coredens_method,coretau_method,fdir,iatom,idir,indx,ipositron,itypat,mu
 integer :: optatm,optdyfr,opteltfr,optgr,option,optn,optn2,optstr,optv,vloc_method
 real(dp) :: eei_dum1,eei_dum2,ucvol,ucvol_local,vol_element
 logical :: calc_epaw3_forces, efield_flag
 logical :: is_hybrid_ncpp
!arrays
 integer :: qprtrb_dum(3)
 real(dp) :: dummy6(6),ep3(3),fioncart(3),gmet(3,3),gprimd(3,3)
 real(dp) :: rmet(3,3),strn_dummy6(6),strv_dummy6(6),tsec(2),vprtrb_dum(2)
 real(dp),allocatable :: atmrho_dum(:),atmvloc_dum(:),dyfrlo_dum(:,:,:)
 real(dp),allocatable :: dyfrn_dum(:,:,:),dyfrv_dum(:,:,:)
 real(dp),allocatable :: dyfrx2_dum(:,:,:),eltfrn_dum(:,:),gauss_dum(:,:)
 real(dp),allocatable :: epawf3red(:,:),fin(:,:),fionred(:,:),grl(:,:),grl_dum(:,:)
 real(dp),allocatable :: grnl_tmp(:,:),grtn(:,:),grtn_indx(:,:),grxctau(:,:),v_dum(:),vxctotg(:,:)
 real(dp),allocatable :: xccc3d_dum(:)

! *************************************************************************

 call timab(69,1,tsec)

!Save input value of forces
 ABI_ALLOCATE(fin,(3,dtset%natom))
 fin(:,:)=fcart(:,:)

!Compute different geometric tensor, as well as ucvol, from rprimd
 call metric(gmet,gprimd,-1,rmet,rprimd,ucvol)

!Check if we're in hybrid norm conserving pseudopotential
 is_hybrid_ncpp=(psps%usepaw==0 .and. &
& (dtset%ixc==41.or.dtset%ixc==42.or.libxc_functionals_is_hybrid()))

!=======================================================================
!========= Local pseudopotential and core charge contributions =========
!=======================================================================

 ABI_ALLOCATE(grl,(3,dtset%natom))

!Determine by which method the local ionic potential and/or the pseudo core
!  charge density contributions have to be computed
!Local ionic potential:
! Method 1: PAW
! Method 2: Norm-conserving PP, icoulomb>0, wavelets
 vloc_method=1;if (psps%usepaw==0) vloc_method=2
 if (dtset%icoulomb>0) vloc_method=2
 if (psps%usewvl==1) vloc_method=2
!Pseudo core charge density:
! Method 1: PAW, nc_xccc_gspace
! Method 2: Norm-conserving PP, wavelets
 coredens_method=1;if (psps%usepaw==0) coredens_method=2
 if (psps%nc_xccc_gspace==1) coredens_method=1
 if (psps%nc_xccc_gspace==0) coredens_method=2
 if (psps%usewvl==1) coredens_method=2
 coretau_method=0
 if (dtset%usekden==1.and.psps%usepaw==1) then
   coretau_method=1;if (psps%nc_xccc_gspace==0) coretau_method=2
 end if

!Local ionic potential and/or pseudo core charge by method 1
 if (vloc_method==1.or.coredens_method==1.or.coretau_method==1) then
   if (psps%nc_xccc_gspace==1.and.psps%usepaw==0.and.is_hybrid_ncpp) then
     MSG_BUG(' Not yet implemented !')
   end if
   call timab(550,1,tsec)
!  Allocate (unused) dummy variables, otherwise some compilers complain
   ABI_ALLOCATE(gauss_dum,(0,0))
   ABI_ALLOCATE(atmrho_dum,(0))
   ABI_ALLOCATE(atmvloc_dum,(0))
   ABI_ALLOCATE(dyfrn_dum,(0,0,0))
   ABI_ALLOCATE(dyfrv_dum,(0,0,0))
   ABI_ALLOCATE(eltfrn_dum,(0,0))
!  Compute Vxc in reciprocal space
   if (coredens_method==1.and.n3xccc>0) then
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxc(:,1);if (dtset%nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxc(:,2))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,1,ngfft,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_DEALLOCATE(v_dum)
   else
     ABI_ALLOCATE(vxctotg,(0,0))
   end if
!  Compute contribution to forces from Vloc and/or pseudo core density
   optv=0;if (vloc_method==1) optv=1
   optn=0;if (coredens_method==1) optn=n3xccc/nfft
   optatm=0;optdyfr=0;optgr=1;optstr=0;optn2=1;opteltfr=0
   if (vloc_method==1.or.coredens_method==1) then
     call atm2fft(atindx1,atmrho_dum,atmvloc_dum,dyfrn_dum,dyfrv_dum,&
&     eltfrn_dum,gauss_dum,gmet,gprimd,&
&     grxc,grl,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&     optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,qprtrb_dum,&
&     dtset%rcut,rhog,rprimd,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,vxctotg,vxctotg,vxctotg,vprtrb_dum,psps%vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
   end if
   if (n3xccc==0.and.coredens_method==1) grxc=zero
   ABI_DEALLOCATE(vxctotg)
   if (dtset%usekden==1.and.coretau_method==1..and.n3xccc>0) then
!    Compute contribution to forces from pseudo kinetic energy core density
     optv=0;optn=1;optn2=4
     ABI_ALLOCATE(grxctau,(3,dtset%natom))
     ABI_ALLOCATE(grl_dum,(0,0))
     ABI_ALLOCATE(v_dum,(nfft))
     ABI_ALLOCATE(vxctotg,(2,nfft))
     v_dum(:)=vxctau(:,1,1);if (dtset%nspden>=2) v_dum(:)=0.5_dp*(v_dum(:)+vxctau(:,2,1))
     call fourdp(1,vxctotg,v_dum,-1,mpi_enreg,nfft,1,ngfft,0)
     call zerosym(vxctotg,2,ngfft(1),ngfft(2),ngfft(3),&
&     comm_fft=mpi_enreg%comm_fft,distribfft=mpi_enreg%distribfft)
     ABI_DEALLOCATE(v_dum)
     call atm2fft(atindx1,atmrho_dum,atmvloc_dum,dyfrn_dum,dyfrv_dum,&
&     eltfrn_dum,gauss_dum,gmet,gprimd,&
&     grxctau,grl_dum,gsqcut,mgfft,psps%mqgrid_vl,dtset%natom,nattyp,nfft,ngfft,ntypat,&
&     optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,psps,pawtab,ph1d,psps%qgrid_vl,qprtrb_dum,&
&     dtset%rcut,rhog,rprimd,strn_dummy6,strv_dummy6,ucvol,psps%usepaw,vxctotg,vxctotg,vxctotg,vprtrb_dum,psps%vlspl,&
&     comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
&     paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)
     grxc(:,:)=grxc(:,:)+grxctau(:,:)
     ABI_DEALLOCATE(grl_dum)
     ABI_DEALLOCATE(grxctau)
     ABI_DEALLOCATE(vxctotg)
   end if
!  Deallocate temporary arrays
   ABI_DEALLOCATE(gauss_dum)
   ABI_DEALLOCATE(atmrho_dum)
   ABI_DEALLOCATE(atmvloc_dum)
   ABI_DEALLOCATE(dyfrn_dum)
   ABI_DEALLOCATE(dyfrv_dum)
   ABI_DEALLOCATE(eltfrn_dum)
   call timab(550,2,tsec)
 end if

!Local ionic potential by method 2
 if (vloc_method==2) then
   option=2
   ABI_ALLOCATE(dyfrlo_dum,(3,3,dtset%natom))
   ABI_ALLOCATE(grtn_indx,(3,dtset%natom))
   ABI_ALLOCATE(v_dum,(nfft))
   call mklocl(dtset,dyfrlo_dum,eei_dum1,gmet,gprimd,grtn_indx,gsqcut,dummy6,mgfft,&
&   mpi_enreg,dtset%natom,nattyp,nfft,ngfft,dtset%nspden,ntypat,option,pawtab,ph1d,psps,&
&   qprtrb_dum,rhog,rhor,rprimd,ucvol,vprtrb_dum,v_dum,wvl,wvl_den,xred)
   do iatom=1,dtset%natom
!    Has to use the indexing array atindx1
     grl(1:3,atindx1(iatom))=grtn_indx(1:3,iatom)
   end do
   ABI_DEALLOCATE(dyfrlo_dum)
   ABI_DEALLOCATE(grtn_indx)
   ABI_DEALLOCATE(v_dum)
!  If gradients are computed in real space, we need to symmetrize the system before summing.
!  Rshaltaf: I changed the following line to include surfaces BC
   if (dtset%icoulomb == 1 .or. dtset%icoulomb == 2) then
     ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
     call sygrad(grnl_tmp,dtset%natom,grl,dtset%nsym,symrec,indsym)
     grl(:, :) = grnl_tmp(:, :)
     ABI_DEALLOCATE(grnl_tmp)
   end if
 end if

!Pseudo core electron density by method 2
 if (coredens_method==2.or.coretau_method==2) then
   if (n1xccc/=0) then
     call timab(53,1,tsec)
     option=2
     ABI_ALLOCATE(dyfrx2_dum,(3,3,dtset%natom))
     ABI_ALLOCATE(xccc3d_dum,(n3xccc))
     if (coredens_method==2) then
       if (is_hybrid_ncpp) then
         call xchybrid_ncpp_cc(dtset,eei_dum1,mpi_enreg,nfft,ngfft,n3xccc,rhor,rprimd,&
&         dummy6,eei_dum2,xccc3d_dum,grxc=grxc,xcccrc=psps%xcccrc,xccc1d=psps%xccc1d,xred=xred,n1xccc=n1xccc)
       else
         if (psps%usewvl==0.and.psps%usepaw==0.and.dtset%icoulomb==0) then
           call mkcore(dummy6,dyfrx2_dum,grxc,mpi_enreg,dtset%natom,nfft,dtset%nspden,ntypat,&
&           ngfft(1),n1xccc, ngfft(2),ngfft(3),option,rprimd,dtset%typat,ucvol,vxc,&
&           psps%xcccrc,psps%xccc1d,xccc3d_dum,xred)
         else if (psps%usewvl==0.and.(psps%usepaw==1.or.dtset%icoulomb==1)) then
           call mkcore_alt(atindx1,dummy6,dyfrx2_dum,grxc,dtset%icoulomb,mpi_enreg,dtset%natom,nfft,&
&           dtset%nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&           ucvol,vxc,psps%xcccrc,psps%xccc1d,xccc3d_dum,xred,pawrad,pawtab,psps%usepaw)
         else if (psps%usewvl==1.and.psps%usepaw==1) then
           ucvol_local=ucvol
#if defined HAVE_BIGDFT
!          ucvol_local=product(wvl_den%denspot%dpbox%hgrids)*real(product(wvl_den%denspot%dpbox%ndims),dp)
!          call mkcore_wvl_old(atindx1,dummy6,dyfrx2_dum,wvl%atoms%astruct%geocode,grxc,wvl%h,dtset%natom,&
! &           nattyp,nfft,wvl_den%denspot%dpbox%nscatterarr(mpi_enreg%me_wvl,:),dtset%nspden,ntypat,&
! &           wvl%Glr%d%n1,wvl%Glr%d%n1i,wvl%Glr%d%n2,wvl%Glr%d%n2i,wvl%Glr%d%n3,wvl_den%denspot%dpbox%n3pi,&
! &           n3xccc,option,pawrad,pawtab,psps%gth_params%psppar,rprimd,ucvol_local,vxc,xccc3d_dum,xred,&
! &           mpi_comm_wvl=mpi_enreg%comm_wvl)
           call mkcore_wvl(atindx1,dummy6,grxc,dtset%natom,nattyp,nfft,dtset%nspden,ntypat,&
&           n1xccc,n3xccc,option,pawrad,pawtab,rprimd,vxc,psps%xccc1d,xccc3d_dum,&
&           psps%xcccrc,xred,wvl_den,wvl,mpi_comm_wvl=mpi_enreg%comm_wvl)
#endif
         end if
       end if
     end if
     if (dtset%usekden==1.and.coretau_method==2) then
       ABI_ALLOCATE(grxctau,(3,dtset%natom))
       call mkcore_alt(atindx1,dummy6,dyfrx2_dum,grxctau,dtset%icoulomb,mpi_enreg,dtset%natom,nfft,&
&       dtset%nspden,nattyp,ntypat,ngfft(1),n1xccc,ngfft(2),ngfft(3),option,rprimd,&
&       ucvol,vxctau(:,:,1),psps%xcccrc,psps%xccc1d,xccc3d_dum,xred,pawrad,pawtab,psps%usepaw,&
&       usekden=.true.)
       grxc(:,:)=grxc(:,:)+grxctau(:,:)
       ABI_DEALLOCATE(grxctau)
     end if
     ABI_DEALLOCATE(xccc3d_dum)
     ABI_DEALLOCATE(dyfrx2_dum)
     call timab(53,2,tsec)
   else
     grxc(:,:)=zero
   end if
 end if

!=======================================================================
!===================== Nonlocal contributions ==========================
!=======================================================================

!Only has to apply symmetries
 ABI_ALLOCATE(grnl_tmp,(3,dtset%natom))
 do iatom=1,dtset%natom
   indx=3*(iatom-1);grnl_tmp(1:3,atindx1(iatom))=grnl(indx+1:indx+3)
 end do
 if (dtset%usewvl == 0) then
   call sygrad(synlgr,dtset%natom,grnl_tmp,dtset%nsym,symrec,indsym)
 else
   synlgr = grnl_tmp
 end if
 ABI_DEALLOCATE(grnl_tmp)

!=======================================================================
!============ Density/potential residual contributions =================
!=======================================================================

 if (dtset%usewvl==0.and.abs(dtset%densfor_pred)>=1.and.abs(dtset%densfor_pred)<=3) then
   call fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,1,&
&   pawtab,rhor,rprimd,ucvol,vresid,xred,xred,psps%znuclpsp)
 else if (dtset%usewvl==0.and.(abs(dtset%densfor_pred)==4.or.abs(dtset%densfor_pred)==6)) then
   call fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,&
&   mpi_enreg,psps%mqgrid_vl,nattyp,nfft,ngfft,ntypat,psps,pawtab,ph1d,&
&   psps%qgrid_vl,rprimd,ucvol,psps%usepaw,vresid,psps%zionpsp,psps%znuclpsp)
 else
   gresid(:,:)=zero
 end if

!=======================================================================
!======================= Other contributions ===========================
!=======================================================================

!Ewald energy contribution to forces as already been computed in "ewald"

!Potential residual contribution to forces as already been computed (forstr)

!Add Berry phase contributions (berryopt == 4/6/7/14/16/17)
!(compute the electric field force on the ion cores)
 efield_flag = (dtset%berryopt==4 .or. dtset%berryopt==6 .or. dtset%berryopt==7 .or. &
& dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17)
 calc_epaw3_forces = (efield_flag .and. dtset%optforces /= 0 .and. psps%usepaw == 1)
 if ( efield_flag ) then
   ABI_ALLOCATE(fionred,(3,dtset%natom))
   fionred(:,:)=zero
   do iatom=1,dtset%natom
     itypat=dtset%typat(iatom)
! force on ion due to electric field, cartesian representation
     fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
! form fionred = rprimd^T * fioncart, note that forces transform
! oppositely to coordinates, because they are derivative with respect to
! coordinates
     call dgemv('T',3,3,one,rprimd,3,fioncart,1,zero,fionred(1:3,iatom),1)
!     do mu=1,3
!       fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
!&       +rprimd(2,mu)*fioncart(2) &
!&       +rprimd(3,mu)*fioncart(3)
!     end do
   end do
 end if

!(compute additional F3-type force due to projectors for electric field with PAW)
 if ( efield_flag .and. calc_epaw3_forces ) then
   ABI_ALLOCATE(epawf3red,(3,dtset%natom))
! dtefield%epawf3(iatom,idir,fdir) contains
   epawf3red(:,:)=zero
   do iatom=1,dtset%natom
     do fdir = 1, 3
       do idir = 1, 3
! vol_element is volume/pt for integration of epawf3. volume is BZ volume
! so 1/ucvol, and number of kpts is nstr(idir)*nkstr(idir)
         vol_element=one/(ucvol*dtefield%nstr(idir)*dtefield%nkstr(idir))
         ep3(idir) = vol_element*dtefield%epawf3(iatom,idir,fdir)
       end do
       epawf3red(fdir,iatom) = -ucvol*dot_product(dtset%red_efieldbar(1:3),ep3(1:3))
     end do
   end do ! end loop over iatom
 end if

!This was incorrect coding. Bug found by Jiawang Hong
!if (dtset%berryopt==4) then
!allocate(fionred(3,dtset%natom));fionred(:,:)=zero
!iatom = 0
!do itypat=1,ntypat
!do iattyp=1,nattyp(itypat)
!iatom=iatom+1
!fioncart(:)=psps%ziontypat(itypat)*dtset%efield(:)
!do mu=1,3
!fionred(mu,iatom)=rprimd(1,mu)*fioncart(1) &
!&         +rprimd(2,mu)*fioncart(2) &
!&         +rprimd(3,mu)*fioncart(3)
!end do
!end do
!end do
!end if

!=======================================================================
!======= Assemble the various contributions to the forces ==============
!=======================================================================

!Collect grads of etot wrt reduced coordinates
!This gives non-symmetrized Hellman-Feynman reduced gradients
 ABI_ALLOCATE(grtn,(3,dtset%natom))
 grtn(:,:)=grl(:,:)+grchempottn(:,:)+grcondft(:,:)+grewtn(:,:)+synlgr(:,:)+grxc(:,:)

 if (usefock==1 .and. associated(fock)) then
   if (fock%fock_common%optfor) then
     grtn(:,:)=grtn(:,:)+fock%fock_common%forces(:,:)
   end if
 end if

 if (ngrvdw==dtset%natom) grtn(:,:)=grtn(:,:)+grvdw(:,:)
! note that fionred is subtracted, because it really is a force and we need to
! turn it back into a gradient. The fred2fcart routine below includes the minus
! sign to convert gradients back to forces
 if ( efield_flag ) grtn(:,:)=grtn(:,:)-fionred(:,:)
! epawf3red is added, because it actually is a gradient, not a force
 if ( efield_flag .and. calc_epaw3_forces ) grtn(:,:) = grtn(:,:) + epawf3red(:,:)

!Symmetrize explicitly for given space group and store in grhf :
 call sygrad(grhf,dtset%natom,grtn,dtset%nsym,symrec,indsym)

!If residual forces are too large, there must be a problem: cancel them !
 if (dtset%usewvl==0.and.abs(dtset%densfor_pred)>0.and.abs(dtset%densfor_pred)/=5) then
   do iatom=1,dtset%natom
     do mu=1,3
       if (abs(gresid(mu,iatom))>10000._dp*abs(grtn(mu,iatom))) gresid(mu,iatom)=zero
     end do
   end do
 end if

!Add residual potential correction
 grtn(:,:)=grtn(:,:)+gresid(:,:)

!Additional stuff for electron-positron
 ipositron=0
 if (present(electronpositron)) then
   if (associated(electronpositron)) then
     if (allocated(electronpositron%fred_ep)) ipositron=electronpositron_calctype(electronpositron)
   end if
 end if
 if (abs(ipositron)==1) then
   grtn(:,:)=grtn(:,:)-grxc(:,:)-grchempottn(:,:)-grcondft(:,:)-grewtn(:,:)-gresid(:,:)-two*grl(:,:)
!  grtn(:,:)=grtn(:,:)-grxc(:,:)-grewtn(:,:)-gresid(:,:)-two*grl(:,:)
   grl(:,:)=-grl(:,:);grxc(:,:)=zero;gresid(:,:)=zero
   if (ngrvdw==dtset%natom) grtn(:,:)=grtn(:,:)-grvdw(:,:)
   if ( dtset%berryopt== 4 .or. dtset%berryopt== 6 .or. dtset%berryopt== 7 .or. &
&   dtset%berryopt==14 .or. dtset%berryopt==16 .or. dtset%berryopt==17)  then
     grtn(:,:)=grtn(:,:)+fionred(:,:)
     fionred(:,:)=zero
   end if
 end if
 if (ipositron>0) grtn(:,:)=grtn(:,:)+electronpositron%fred_ep(:,:)

!Symmetrize all grads explicitly for given space group:
 if (dtset%usewvl == 0) then
   call sygrad(fred,dtset%natom,grtn,dtset%nsym,symrec,indsym)
 else
   fred = grtn
 end if

!Conversion to cartesian coordinates (bohr) AND
!Subtract off average force from each force component
!to avoid spurious drifting of atoms across cell.
! notice that fred2fcart multiplies fred by -1 to convert it
! from a gradient (input) to a force (output)

 call fred2fcart(favg,(dtset%jellslab==0 .and. dtset%nzchempot==0),fcart,fred,gprimd,dtset%natom)

!Compute maximal force and maximal difference
 maxfor=zero;diffor=zero
 do iatom=1,dtset%natom
   do mu=1,3
     if (dtset%iatfix(mu,iatom) /= 1) then
       maxfor=max(maxfor,abs(fcart(mu,iatom)))
       diffor=max(diffor,abs(fcart(mu,iatom)-fin(mu,iatom)))
     else if (dtset%ionmov==4 .or. dtset%ionmov==5) then
!      Make the force vanish on fixed atoms when ionmov=4 or 5
!      This is because fixing of atom cannot be imposed at the
!      level of a routine similar to brdmin or moldyn for these options.
       fcart(mu,iatom)=zero
     end if
   end do
 end do

!Apply any generalized constraints to the forces
 if (dtset%nconeq>0) call constrf(diffor,fcart,forold,fred,dtset%iatfix,dtset%ionmov,maxfor,&
& dtset%natom,dtset%nconeq,dtset%prtvol,rprimd,dtset%wtatcon,xred)

!=======================================================================
!Memory deallocations
 ABI_DEALLOCATE(grl)
 ABI_DEALLOCATE(grtn)
 ABI_DEALLOCATE(fin)
 if ( efield_flag )  then
   ABI_DEALLOCATE(fionred)
   if ( calc_epaw3_forces ) then
     ABI_DEALLOCATE(epawf3red)
   end if
 end if

 call timab(69,2,tsec)

end subroutine forces
!!***

!!****f* ABINIT/sygrad
!!
!! NAME
!! sygrad
!!
!! FUNCTION
!! Symmetrize derivatives of energy with respect to coordinates.
!! Unsymmetrized gradients are input as dedt; symmetrized grads are then placed in fred.
!! If nsym=1 simply copy dedt into fred (only symmetry is identity).
!!
!! INPUTS
!!  natom=number of atoms in cell
!!  dedt(3,natom)=unsymmetrized gradients wrt dimensionless tn (hartree)
!!  nsym=number of symmetry operators in group
!!  symrec(3,3,nsym)=symmetries of group in terms of operations on
!!    reciprocal space primitive translations--see comments below
!!  indsym(4,nsym,natom)=label given by subroutine symatm, indicating atom
!!   label which gets rotated into given atom by given symmetry
!!   (first three elements are related primitive translation--
!!   see symatm where this is computed)
!!
!! OUTPUT
!! fred(3,3,natom)=symmetrized gradients wrt reduced coordinates (hartree)
!!
!! NOTES
!! symmetrization of gradients with respect to reduced
!! coordinates tn is conducted according to the expression
!! $[d(e)/d(t(n,a))]_{symmetrized} = (1/Nsym)*Sum(S)*symrec(n,m,S)*
!!              [d(e)/d(t(m,b))]_{unsymmetrized}$
!! where $t(m,b)= (symrel^{-1})(m,n)*(t(n,a)-tnons(n))$ and tnons
!! is a possible nonsymmorphic translation.  The label "b" here
!! refers to the atom which gets rotated into "a" under symmetry "S".
!! symrel is the symmetry matrix in real space, which is the inverse
!! transpose of symrec.  symrec is the symmetry matrix in reciprocal
!! space.  $sym_{cartesian} = R * symrel * R^{-1} = G * symrec * G^{-1}$
!! where the columns of R and G are the dimensional primitive translations
!! in real and reciprocal space respectively.
!! Note the use of "symrec" in the symmetrization expression above.
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!
!! SOURCE

subroutine sygrad(fred,natom,dedt,nsym,symrec,indsym)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: natom,nsym
!arrays
 integer,intent(in) :: indsym(4,nsym,natom),symrec(3,3,nsym)
 real(dp),intent(in) :: dedt(3,natom)
 real(dp),intent(out) :: fred(3,natom)

!Local variables-------------------------------
!scalars
 integer :: ia,ind,isym,mu
 real(dp),parameter :: tol=1.0d-30
 real(dp) :: summ

! *************************************************************************
!
 if (nsym==1) then
!  only symmetry is identity so simply copy
   fred(:,:)=dedt(:,:)
 else
!  actually conduct symmetrization
   do ia=1,natom
     do mu=1,3
       summ=0._dp
       do isym=1,nsym
         ind=indsym(4,isym,ia)
         summ=summ+dble(symrec(mu,1,isym))*dedt(1,ind)+&
&         dble(symrec(mu,2,isym))*dedt(2,ind)+&
&         dble(symrec(mu,3,isym))*dedt(3,ind)
       end do
       fred(mu,ia)=summ/dble(nsym)
       if(abs(fred(mu,ia))<tol)fred(mu,ia)=0.0_dp
     end do
   end do
 end if

end subroutine sygrad
!!***

!!****f* ABINIT/fresidrsp
!!
!! NAME
!! fresidrsp
!!
!! FUNCTION
!! Compute the forces due to the residual of the potential (or density)
!! in RECIPROCAL SPACE, using
!!  - the atomic density read in psp file (PAW or NC with nctval_spl e.g. psp8 format)
!!  - a gaussian atomic density (norm-conserving psps if nctval_spl is not available)
!!
!! INPUTS
!! atindx1(natom)=index table for atoms, inverse of atindx
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | densty(ntypat,4)=parameters for initialisation of the density of each atom type
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | ixc= choice of exchange-correlation scheme
!!  | natom=number of atoms in cell.
!!  | nspden=number of spin-density components
!!  | typat(natom)=integer type for each atom in cell
!! gmet(3,3)=reciprocal space metric
!! gprimd(3,3)=reciprocal space dimensional primitive translations
!! gsqcut=cutoff value on G**2 for sphere inside fft box
!! mgfft=maximum size of 1D FFTs
!! mpi_enreg=information about MPI parallelization
!! mqgrid=number of grid pts in q array for atomic density spline n^AT(q)
!! nattyp(ntypat)=number of atoms of each type in cell
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! psps <type(pseudopotential_type)>=variables related to pseudopotentials
!! pawtab(ntypat*usepaw) <type(pawtab_type)>=paw tabulated starting data
!! ph1d(2,3*(2*mgfft+1)*natom)=1-dim phase information for given atom coordinates.
!! qgrid(mqgrid)=q grid for spline atomic valence density n^AT(q) from 0 to qmax
!! ucvol=unit cell volume (bohr**3).
!! usepaw= 0 for non paw calculation; =1 for paw calculation
!! vresid(nfft,nspden)=potential residual - (non-collinear magn. : only V11 and V22 are used)
!! zion(ntypat)=charge on each type of atom (real number)
!! znucl(ntypat)=atomic number, for each type of atom
!!
!! OUTPUT
!! gresid(3,natom)=forces due to the residual of the potential
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!      atm2fft,fourdp,wrtout
!!
!! SOURCE

subroutine fresidrsp(atindx1,dtset,gmet,gprimd,gresid,gsqcut,mgfft,mpi_enreg,mqgrid,nattyp,nfft,&
&          ngfft,ntypat,psps,pawtab,ph1d,qgrid,rprimd,ucvol,usepaw,vresid,zion,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: mgfft,mqgrid,nfft,ntypat,usepaw
 real(dp),intent(in) :: gsqcut,ucvol
 type(pseudopotential_type),intent(in) :: psps
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: atindx1(dtset%natom),nattyp(ntypat),ngfft(18)
 real(dp),intent(in) :: gmet(3,3),gprimd(3,3),ph1d(2,3*(2*mgfft+1)*dtset%natom)
 real(dp),intent(in) :: qgrid(mqgrid),vresid(nfft,dtset%nspden),zion(ntypat)
 real(dp),intent(inout) :: rprimd(3,3)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(out) :: gresid(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*usepaw)

!Local variables-------------------------------
!scalars
 integer :: itypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv
 logical :: usegauss
!arrays
 integer :: dummy3(3)
 real(dp) :: dummy2(2)
 real(dp) :: dummy_in1(0),dummy_in2(0)
 real(dp) :: dummy_out1(0),dummy_out2(0),dummy_out3(0),dummy_out4(0),dummy_out5(0),dummy_out6(0)
 real(dp) :: strn_dummy6(6),strv_dummy6(6)
 real(dp),allocatable :: gauss(:,:),vresg(:,:),work(:)

! *************************************************************************

!Inits
 optatm=0;optdyfr=0;opteltfr=0;optgr=1;optstr=0;optv=0;optn=1
 ABI_ALLOCATE(vresg,(2,nfft))

!Transfer potential residual to reciprocal space
!Use only Vres=Vres11+Vres22=Vres_up+Vres_dn
 ABI_ALLOCATE(work,(nfft))
 work(:)=vresid(:,1)
 if (dtset%nspden>=2) work(:)=work(:)+vresid(:,2)
 call fourdp(1,vresg,work,-1,mpi_enreg,nfft,1,ngfft,0)
 ABI_DEALLOCATE(work)

!Determine wether a gaussan atomic density has to be used or not
 usegauss=.true.
 if (usepaw==0) usegauss = any(.not.psps%nctab(1:ntypat)%has_tvale)
 if (usepaw==1) usegauss=(minval(pawtab(1:ntypat)%has_tvale)==0)
 if (usegauss) then
   optn2=3
   ABI_ALLOCATE(gauss,(2,ntypat))
   do itypat=1,ntypat
     gauss(1,itypat)=zion(itypat)
     gauss(2,itypat) = atom_length(dtset%densty(itypat,1),zion(itypat),znucl(itypat))
   end do
   call wrtout(std_out," Computing residual forces using gaussian functions as atomic densities", "COLL")
 else
   optn2=2
   ABI_ALLOCATE(gauss,(2,0))
   call wrtout(std_out," Computing residual forces using atomic densities taken from pseudos", "COLL")
 end if

!Compute forces due to residual
 call atm2fft(atindx1,dummy_out1,dummy_out2,dummy_out3,dummy_out4,&
& dummy_out5,gauss,gmet,gprimd,gresid,dummy_out6,gsqcut,mgfft,&
& mqgrid,dtset%natom,nattyp,nfft,ngfft,ntypat,optatm,optdyfr,opteltfr,optgr,optn,optn2,optstr,optv,&
& psps,pawtab,ph1d,qgrid,dummy3,dtset%rcut,dummy_in1,rprimd,strn_dummy6,strv_dummy6,ucvol,usepaw,&
& vresg,vresg,vresg,dummy2,dummy_in2,comm_fft=mpi_enreg%comm_fft,me_g0=mpi_enreg%me_g0,&
& paral_kgb=mpi_enreg%paral_kgb,distribfft=mpi_enreg%distribfft)

!In case of nspden>=2, has to apply 1/2 factor
 if (dtset%nspden>=2) gresid=gresid*half

 ABI_DEALLOCATE(gauss)
 ABI_DEALLOCATE(vresg)

end subroutine fresidrsp
!!***

!!****f* ABINIT/fresid
!! NAME
!! fresid
!!
!! FUNCTION
!! If option=1, compute the forces due to the residual of the potential
!! If option=2, generate approximate new density from old one,
!!              old atomic positions, and new atomic positions
!!
!! INPUTS
!! dtset <type(dataset_type)>=all input variables in this dataset
!!  | icoulomb=0 periodic treatment of Hartree potential, 1 use of Poisson solver
!!  | natom=number of atoms in cell.
!!  | nspden=number of spin-density components
!!  | typat(natom)=integer type for each atom in cell
!!  | usepaw= 0 for non paw calculation; =1 for paw calculation
!!  | xclevel= level of the XC functional
!! mpi_enreg=information about MPI parallelization
!! nfft=(effective) number of FFT grid points (for this processor)
!! ngfft(18)=contain all needed information about 3D FFT, see ~abinit/doc/variables/vargs.htm#ngfft
!! ntypat=number of types of atoms in cell.
!! option=see below
!! pawtab(ntypat*dtset%usepaw) <type(pawtab_type)>=paw tabulated starting data
!! rhor(nfft,nspden)=electron density in electrons/bohr**3 (slices of it if FTT parallelism).
!! rprimd(3,3)=dimensional primitive translation vectors (bohr)
!! ucvol=unit cell volume (bohr**3).
!! xred_new(3,natom)=new reduced coordinates for atoms in unit cell
!! xred_old(3,natom)=old reduced coordinates for atoms in unit cell
!! znucl(ntypat)=real(dp), atomic number of atom type
!!
!! OUTPUT
!! gresid(3,natom)=forces due to the residual of the potential
!!
!! SIDE EFFECTS
!! work(nfft,nspden)=functions on the fft grid (slices of it if FTT parallelism):
!!  if option==1, the POTENTIAL residual is input
!!  if option==2, the interpolated density is output
!!
!! NOTES
!! FFT parallelism:
!! At the beginning of this routine, the plane-waves are ditributed over all the processors.
!! In the main part, all the processors perform the same calculations over the whole FFT grid.
!! At the end, each processor gets its part of the whole FFT grid.
!! These modifications are not efficient when large FFT grids are used.
!! So they have to be considered as a first step before a comprehensive parallelization of this routine.
!!
!! PARENTS
!!      forces,prcref,prcref_PMA,scfcv
!!
!! CHILDREN
!!      atomdata_from_znucl,mean_fftr,pre_gather,pre_scatter
!!
!! SOURCE

subroutine fresid(dtset,gresid,mpi_enreg,nfft,ngfft,ntypat,option,&
&                 pawtab,rhor,rprimd,ucvol,work,xred_new,xred_old,znucl)

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: nfft,ntypat,option
 real(dp),intent(in) :: ucvol
 type(MPI_type),intent(in) :: mpi_enreg
 type(dataset_type),intent(in) :: dtset
!arrays
 integer,intent(in) :: ngfft(18)
 real(dp),intent(in) :: rhor(nfft,dtset%nspden),rprimd(3,3)
 real(dp),intent(in) :: xred_new(3,dtset%natom),xred_old(3,dtset%natom)
 real(dp),intent(in) :: znucl(ntypat)
 real(dp),intent(inout) :: work(nfft,dtset%nspden)
 real(dp),intent(out) :: gresid(3,dtset%natom)
 type(pawtab_type),intent(in) :: pawtab(ntypat*dtset%usepaw)

!Local variables-------------------------------
!real(dp), parameter :: app_remain=0.001_dp
!scalars
 integer,parameter :: natnum=110
 integer :: atmove,i1,i1_new,i1m,i1p,i2,i2_new,i2m,i2p,i3,i3_new,i3m,i3p
 integer :: iatom,ifft,ifft_new,iloop,ind2m,ind2m3m,ind2p,ind2p3p,ind3m,ind3p
 integer :: index,index_new,ishift,ishift1,ishift2,ishift3,ispden,ixp,mshift,mu
 integer :: n1,n2,n3,n4,nfft_tmp,nfftot,nu,quit
 real(dp),parameter :: app_remain=0.01_dp
 real(dp) :: diff_rem1,diff_rem2,diff_rem3,difmag,difmag2
 real(dp) :: difmag2_fact,difmag2_part,drho1,drho11,drho12,drho13,drho14
 real(dp) :: drho1dn,drho1mx,drho1my,drho1mz,drho1tot,drho1up,drho2,drho21
 real(dp) :: drho22,drho23,drho24,drho2dn,drho2mx,drho2my,drho2mz,drho2tot
 real(dp) :: drho2up,drho3,drho31,drho32,drho33,drho34,drho3dn,drho3mx,drho3my
 real(dp) :: drho3mz,drho3tot,drho3up,drhox00,drhox01,drhox10,drhox11,drhoxy0
 real(dp) :: drhoxy1,drhoxyz,fact,range,range2,rcov,rcov2,rcovm1,rdiff1
 real(dp) :: rdiff2,rdiff3,vresid1,vresid2,vresid3,vresid4,xx
 type(atomdata_t) :: atom
!arrays
 integer :: diff_igrid(3),igrid(3),irange(3)
 integer,allocatable :: ii(:,:)
 real(dp) :: diff_grid(3),diff_rem(3),diff_tau(3),diff_xred(3),lencp(3)
 real(dp) :: rho_tot(4),rhosum(4),rmet(3,3),scale(3),tau(3)
 real(dp),allocatable :: approp(:),atmrho(:,:),rhor_tot(:,:),rrdiff(:,:)
 real(dp),allocatable :: work_tot(:,:)
 logical,allocatable :: my_sphere(:)

! *************************************************************************

!Compute lengths of cross products for pairs of primitive
!translation vectors (used in setting index search range below)
 lencp(1)=cross_fr(rprimd(1,2),rprimd(2,2),rprimd(3,2),&
& rprimd(1,3),rprimd(2,3),rprimd(3,3))
 lencp(2)=cross_fr(rprimd(1,3),rprimd(2,3),rprimd(3,3),&
& rprimd(1,1),rprimd(2,1),rprimd(3,1))
 lencp(3)=cross_fr(rprimd(1,1),rprimd(2,1),rprimd(3,1),&
& rprimd(1,2),rprimd(2,2),rprimd(3,2))

!Compute factor R1.(R2xR3)/|R2xR3| etc for 1, 2, 3
!(recall ucvol=R1.(R2xR3))
 scale(:)=ucvol/lencp(:)

!initialize diff_igrid, otherwise valgrind complains
 diff_igrid=0

!Compute metric tensor in real space rmet
 do nu=1,3
   rmet(:,nu)=rprimd(1,:)*rprimd(1,nu)+rprimd(2,:)*rprimd(2,nu)+rprimd(3,:)*rprimd(3,nu)
 end do

!FFT parallelization: Starting from now, calculations are performed on the whole FFT grid
!and no more on slices. The nfft variable becomes nfft_tmp until the end
 n1=ngfft(1);n2=ngfft(2);n3=ngfft(3)
 n4=n3/mpi_enreg%nproc_fft
 nfftot=PRODUCT(ngfft(1:3));nfft_tmp=nfftot
 if(mpi_enreg%paral_kgb==1) then
   ABI_ALLOCATE(rhor_tot,(nfftot,dtset%nspden))
   ABI_ALLOCATE(work_tot,(nfftot,dtset%nspden))
   do ispden=1,dtset%nspden
     call pre_gather(rhor(:,ispden),rhor_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
     call pre_gather(work(:,ispden),work_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
   end do
 end if

 gresid(1:3,1:dtset%natom)=0.0_dp
 quit=0

!Initialize appropriation function
 ABI_ALLOCATE(approp,(nfft_tmp))
 ABI_ALLOCATE(atmrho,(nfft_tmp,dtset%nspden))
 ABI_ALLOCATE(my_sphere,(nfft_tmp))

 approp(:)=app_remain
!First loop over atoms in unit cell : build appropriation function
!Second loop : compute forces
 do iloop=1,2

!  Take into account the remaining density
   if(option==2 .and. iloop==2)then
     if(mpi_enreg%paral_kgb==1) then
!      FFT parallelization: All the processors perform the same calculation.
!      We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
       do ispden=1,dtset%nspden
         do ifft=1,nfft_tmp
           work_tot(ifft,ispden)=rhor_tot(ifft,ispden)*approp(ifft)*app_remain
         end do
       end do
       call mean_fftr(work_tot,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
       rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
     else
       do ispden=1,dtset%nspden
         do ifft=1,nfft_tmp
           work(ifft,ispden)=rhor(ifft,ispden)*approp(ifft)*app_remain
         end do
       end do
       call mean_fftr(work,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
     end if

!    This will be used to restore proper normalization of density
     rho_tot(1:dtset%nspden)=rhosum(1:dtset%nspden)*nfftot
   end if

   do iatom=1,dtset%natom

!    Get the covalent radius
     call atomdata_from_znucl(atom,znucl(dtset%typat(iatom)))
     rcov = atom%rcov
!    PAW choose PAW radius instead...
     if (dtset%usepaw==1) rcov=max(rcov,pawtab(dtset%typat(iatom))%rpaw)

!    Set search range
     rcov2=rcov**2
     range=2._dp*rcov
     range2=range**2
     rcovm1=1.0_dp/rcov

!    Use range to compute an index range along R(1:3)
!    (add 1 to make sure it covers full range)
     irange(1)=1+nint((range/scale(1))*dble(n1))
     irange(2)=1+nint((range/scale(2))*dble(n2))
     irange(3)=1+nint((range/scale(3))*dble(n3))

!    Allocate ii and rrdiff
     mshift=2*maxval(irange(1:3))+1
     ABI_ALLOCATE(ii,(mshift,3))
     ABI_ALLOCATE(rrdiff,(mshift,3))

!    Consider each component in turn
     do mu=1,3

!      Convert reduced coord of given atom to [0,1)
       tau(mu)=mod(xred_old(mu,iatom)+1._dp-aint(xred_old(mu,iatom)),1._dp)

!      Use tau to find nearest grid point along R(mu)
!      (igrid=0 is the origin; shift by 1 to agree with usual index)
       igrid(mu)=nint(tau(mu)*dble(ngfft(mu)))

!      Set up a counter that explore the relevant range
!      of points around the atom
       ishift=0
       do ixp=igrid(mu)-irange(mu),igrid(mu)+irange(mu)
         ishift=ishift+1
         ii(ishift,mu)=1+mod(ngfft(mu)+mod(ixp,ngfft(mu)),ngfft(mu))
         rrdiff(ishift,mu)=dble(ixp)/dble(ngfft(mu))-tau(mu)
       end do

!      If option 2, set up quantities related with the change of atomic coordinates
       if(option==2 .and. iloop==2)then
         diff_xred(mu)=xred_new(mu,iatom)-xred_old(mu,iatom)
!        Convert to [0,1)
         diff_tau(mu)=mod(diff_xred(mu)+1._dp-aint(diff_xred(mu)),1._dp)
!        Convert to [0,ngfft)
         diff_grid(mu)=diff_tau(mu)*dble(ngfft(mu))
!        Integer part
         diff_igrid(mu)=int(diff_grid(mu))
!        Compute remainder
         diff_rem(mu)=diff_grid(mu)-diff_igrid(mu)

!        DEBUG
!        write(std_out,*)' mu,diff',mu,diff_igrid(mu),diff_rem(mu)
!        ENDDEBUG

       end if

!      End loop on mu
     end do

!    May be the atom is fixed
     atmove=1
     if(option==2 .and. iloop==2)then
       if(diff_xred(1)**2+diff_xred(2)**2+diff_xred(3)**2 < 1.0d-24)then
         atmove=0
       else
         diff_rem1=diff_rem(1)
         diff_rem2=diff_rem(2)
         diff_rem3=diff_rem(3)
       end if
     end if

!    If second loop, initialize atomic density, and the variable
!    that says whether a fft point belongs to the sphere of the atom
     if(iloop==2) then
       atmrho(:,:)=0.0_dp
       my_sphere(:)=.false.
     end if

!    Conduct triple loop over restricted range of grid points for iatom

     do ishift3=1,1+2*irange(3)
!      map back to [1,ngfft(3)] for usual fortran index in unit cell
       i3=ii(ishift3,3)
       i3m=i3-1 ; if(i3==1)i3m=n3
       i3p=i3+1 ; if(i3==n3)i3p=1

!      find vector from atom location to grid point (reduced)
       rdiff3=rrdiff(ishift3,3)

       do ishift2=1,1+2*irange(2)
         i2=ii(ishift2,2)
         i2m=i2-1 ; if(i2==1)i2m=n2
         i2p=i2+1 ; if(i2==n2)i2p=1
         index=n1*(i2-1+n2*(i3-1))
         ind3m=n1*(i2-1+n2*(i3m-1))
         ind3p=n1*(i2-1+n2*(i3p-1))
         ind2m=n1*(i2m-1+n2*(i3-1))
         ind2p=n1*(i2p-1+n2*(i3-1))
         ind2p3p=n1*(i2p-1+n2*(i3p-1))

         rdiff2=rrdiff(ishift2,2)
!        Prepare the computation of difmag2
         difmag2_part=rmet(3,3)*rdiff3**2+rmet(2,2)*rdiff2**2&
&         +2.0_dp*rmet(3,2)*rdiff3*rdiff2
         difmag2_fact=2.0_dp*(rmet(3,1)*rdiff3+rmet(2,1)*rdiff2)

         do ishift1=1,1+2*irange(1)
           rdiff1=rrdiff(ishift1,1)

!          Compute (rgrid-tau-Rprim)**2
           difmag2= difmag2_part+rdiff1*(difmag2_fact+rmet(1,1)*rdiff1)

!          Only accept contribution inside defined range
!          This condition means that x, calculated below, cannot exceed 2.0_dp
           if (difmag2<range2) then

!            Will compute contribution to appropriation function based on
!            rcov2, range2 and difmag2
             i1=ii(ishift1,1)
             ifft=i1+index

             if(iloop==1)then

!              Build appropriation function
               if (difmag2<rcov2)then
                 approp(ifft)=approp(ifft)+1.0_dp
               else
                 difmag=sqrt(difmag2)
                 xx=difmag*rcovm1
!                The following function is 1. at xx=1, 0. at xx=2, with vanishing
!                derivatives at these points.
                 approp(ifft)=approp(ifft)+((2.0_dp*xx-9.0_dp)*xx+12.0_dp)*xx-4.0_dp
               end if

             else

               if (difmag2<rcov2) then
                 fact=one
               else
                 difmag=sqrt(difmag2)
                 xx=difmag*rcovm1
                 fact=((2.0_dp*xx-9.0_dp)*xx+12.0_dp)*xx-4.0_dp
               end if

!              Build atomic density
               if(mpi_enreg%paral_kgb==1) then
                 atmrho(ifft,1:dtset%nspden)=atmrho(ifft,1:dtset%nspden) &
&                 +rhor_tot(ifft,1:dtset%nspden)*fact*approp(ifft)
               else
                 atmrho(ifft,1:dtset%nspden)=atmrho(ifft,1:dtset%nspden) &
&                 +rhor(ifft,1:dtset%nspden)*fact*approp(ifft)
               end if

!              Compute the sphere of the atom : it is different for
!              option 1 and for option 2
               i1p=i1+1 ; if(i1==n1)i1p=1
               if(option==1)then
                 i1m=i1-1 ; if(i1==1)i1m=n1
                 my_sphere(ifft)=.true.
                 my_sphere(i1p+index)=.true. ; my_sphere(i1m+index)=.true.
                 my_sphere(i1+ind2p)=.true. ; my_sphere(i1+ind2m)=.true.
                 my_sphere(i1+ind3p)=.true. ; my_sphere(i1+ind3m)=.true.
               else
                 my_sphere(ifft)=.true. ; my_sphere(i1p+index)=.true.
                 my_sphere(i1+ind2p)=.true. ; my_sphere(i1p+ind2p)=.true.
                 my_sphere(i1+ind3p)=.true. ; my_sphere(i1p+ind3p)=.true.
                 my_sphere(i1+ind2p3p)=.true. ; my_sphere(i1p+ind2p3p)=.true.
               end if

             end if

!            End of condition on the range
           end if

!          End loop on ishift1
         end do

!        End loop on ishift2
       end do

!      End loop on ishift3
     end do
!    At the end of the second loop for each atom, compute the force
!    from the atomic densities, or translate density.
!    In the first case, use a two-point finite-difference approximation,
!    since this calculation serves only to decrease the error,
!    and should not be very accurate, but fast.
!    In the second case, using a crude trilinear interpolation scheme
!    for the same reason.
!
!    The section is skipped if option==2 and the atom is fixed
     if(iloop==2 .and. (option==1 .or. atmove==1) )then

       do i3=1,n3
         i3m=i3-1 ; if(i3==1)i3m=n3
         i3p=i3+1 ; if(i3==n3)i3p=1
!        note: diff_igrid is only set  if(option==2 .and. iloop==2)
         i3_new=i3+diff_igrid(3) ; if(i3_new > n3)i3_new=i3_new-n3
         do i2=1,n2
           i2m=i2-1 ; if(i2==1)i2m=n2
           i2p=i2+1 ; if(i2==n2)i2p=1
           i2_new=i2+diff_igrid(2) ; if(i2_new > n2)i2_new=i2_new-n2
           index=n1*(i2-1+n2*(i3-1))
           index_new=n1*(i2_new-1+n2*(i3_new-1))
           ind3m=n1*(i2-1+n2*(i3m-1))
           ind3p=n1*(i2-1+n2*(i3p-1))
           ind2m=n1*(i2m-1+n2*(i3-1))
           ind2p=n1*(i2p-1+n2*(i3-1))
           ind2m3m=n1*(i2m-1+n2*(i3m-1))
           do i1=1,n1
             ifft=i1+index
             if(my_sphere(ifft))then

               i1m=i1-1 ; if(i1==1)i1m=n1

               if(option==1)then
!                Treat option 1 : computation of residual forces
                 i1p=i1+1 ; if(i1==n1)i1p=1
!                Distinguish spin-unpolarized and spin-polarized
                 if(dtset%nspden==1)then ! Non magnetic
!                  Note that the factor needed to obtain a true finite difference
!                  estimation of the derivative will be applied afterwards, for speed
                   drho1=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho2=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho3=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                   else
                     vresid1=work(ifft,1)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho1*vresid1
                   gresid(2,iatom)=gresid(2,iatom)+drho2*vresid1
                   gresid(3,iatom)=gresid(3,iatom)+drho3*vresid1
                 else if(dtset%nspden==2) then ! Collinear magnetism
                   drho1tot=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho2tot=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho3tot=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   drho1up=atmrho(i1p+index,2)-atmrho(i1m+index,2)
                   drho2up=atmrho(i1+ind2p,2) -atmrho(i1+ind2m,2)
                   drho3up=atmrho(i1+ind3p,2) -atmrho(i1+ind3m,2)
                   drho1dn=drho1tot-drho1up
                   drho2dn=drho2tot-drho2up
                   drho3dn=drho3tot-drho3up
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                     vresid2=work_tot(ifft,2)
                   else
                     vresid1=work(ifft,1)
                     vresid2=work(ifft,2)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho1up*vresid1+drho1dn*vresid2
                   gresid(2,iatom)=gresid(2,iatom)+drho2up*vresid1+drho2dn*vresid2
                   gresid(3,iatom)=gresid(3,iatom)+drho3up*vresid1+drho3dn*vresid2
                 else ! Non-collinear magnetism
                   drho1tot=atmrho(i1p+index,1)-atmrho(i1m+index,1)
                   drho1mx =atmrho(i1p+index,2)-atmrho(i1m+index,2)
                   drho1my =atmrho(i1p+index,3)-atmrho(i1m+index,3)
                   drho1mz =atmrho(i1p+index,4)-atmrho(i1m+index,4)
                   drho2tot=atmrho(i1+ind2p,1) -atmrho(i1+ind2m,1)
                   drho2mx =atmrho(i1+ind2p,2) -atmrho(i1+ind2m,2)
                   drho2my =atmrho(i1+ind2p,3) -atmrho(i1+ind2m,3)
                   drho2mz =atmrho(i1+ind2p,4) -atmrho(i1+ind2m,4)
                   drho3tot=atmrho(i1+ind3p,1) -atmrho(i1+ind3m,1)
                   drho3mx =atmrho(i1+ind3p,2) -atmrho(i1+ind3m,2)
                   drho3my =atmrho(i1+ind3p,3) -atmrho(i1+ind3m,3)
                   drho3mz =atmrho(i1+ind3p,4) -atmrho(i1+ind3m,4)
                   drho11=half*(drho1tot+drho1mz)
                   drho12=half*(drho1tot-drho1mz)
                   drho13= half*drho1mx
                   drho14=-half*drho1my
                   drho21=half*(drho2tot+drho2mz)
                   drho22=half*(drho2tot-drho2mz)
                   drho23= half*drho2mx
                   drho24=-half*drho2my
                   drho31=half*(drho3tot+drho3mz)
                   drho32=half*(drho3tot-drho3mz)
                   drho33= half*drho3mx
                   drho34=-half*drho3my
                   if(mpi_enreg%paral_kgb==1) then
                     vresid1=work_tot(ifft,1)
                     vresid2=work_tot(ifft,2)
                     vresid3=work_tot(ifft,3)
                     vresid4=work_tot(ifft,4)
                   else
                     vresid1=work(ifft,1)
                     vresid2=work(ifft,2)
                     vresid3=work(ifft,3)
                     vresid4=work(ifft,4)
                   end if
                   gresid(1,iatom)=gresid(1,iatom)+drho11*vresid1+drho12*vresid2+two*(drho13*vresid3+drho14*vresid4)
                   gresid(2,iatom)=gresid(2,iatom)+drho21*vresid1+drho22*vresid2+two*(drho23*vresid3+drho24*vresid4)
                   gresid(3,iatom)=gresid(3,iatom)+drho31*vresid1+drho32*vresid2+two*(drho33*vresid3+drho34*vresid4)
                 end if
!                Treat the case option==2 now : trilinear interpolation of the density
               else
                 i1_new=i1+diff_igrid(1) ; if(i1_new > n1)i1_new=i1_new-n1
                 ifft_new=i1_new+index_new
                 do ispden=1,dtset%nspden
                   drhox00=(atmrho(i1m+index,ispden)-atmrho(i1+index,ispden))*diff_rem1 &
&                   +atmrho(i1+index,ispden)
                   drhox10=(atmrho(i1m+ind2m,ispden)-atmrho(i1+ind2m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind2m,ispden)
                   drhox01=(atmrho(i1m+ind3m,ispden)-atmrho(i1+ind3m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind3m,ispden)
                   drhox11=(atmrho(i1m+ind2m3m,ispden)-atmrho(i1+ind2m3m,ispden))*diff_rem1 &
&                   +atmrho(i1+ind2m3m,ispden)
                   drhoxy0=(drhox10-drhox00)*diff_rem2+drhox00
                   drhoxy1=(drhox11-drhox01)*diff_rem2+drhox01
                   drhoxyz=(drhoxy1-drhoxy0)*diff_rem3+drhoxy0
                   if(mpi_enreg%paral_kgb==1) then
                     work_tot(ifft_new,ispden)=work_tot(ifft_new,ispden)+drhoxyz
                   else
                     work(ifft_new,ispden)=work(ifft_new,ispden)+drhoxyz
                   end if
                   rho_tot(ispden)=rho_tot(ispden)+drhoxyz
                 end do
               end if

!              End condition of belonging to the sphere of influence of the atom
             end if
           end do
         end do
       end do
!      The finite-difference factor applied here also take
!      into account diverse factors
       fact=-ucvol/dble(nfftot)
       gresid(1,iatom)=gresid(1,iatom)*dble(n1)*.5_dp*fact
       gresid(2,iatom)=gresid(2,iatom)*dble(n2)*.5_dp*fact
       gresid(3,iatom)=gresid(3,iatom)*dble(n3)*.5_dp*fact
     end if

!    Update work if the atom is fixed.
     if(iloop==2 .and. option==2 .and. atmove==0)then
       if(mpi_enreg%paral_kgb==1) then
!        FFT parallelization: All the processors perform the same calculation.
!        We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
         do ispden=1,dtset%nspden
           do ifft=1,nfft_tmp
             work_tot(ifft,ispden)=work_tot(ifft,ispden)+atmrho(ifft,ispden)
           end do
         end do
         call mean_fftr(atmrho,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
         rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
       else
         do ispden=1,dtset%nspden
           do ifft=1,nfft_tmp
             work(ifft,ispden)=work(ifft,ispden)+atmrho(ifft,ispden)
           end do
         end do
         call mean_fftr(atmrho,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
       end if

       rho_tot(1:dtset%nspden)=rho_tot(1:dtset%nspden)+rhosum(1:dtset%nspden)*nfftot
     end if

     ABI_DEALLOCATE(ii)
     ABI_DEALLOCATE(rrdiff)

!    End loop on atoms
   end do

!  DEBUG
!  if(option==2)then
!  if(iloop==1)then
!  write(std_out,*)' fresid : rhor, approp'
!  do ifft=1,n1
!  write(std_out,*)ifft,rhor(ifft,1),approp(ifft)
!  end do
!  end if
!  if(iloop==2)then
!  write(std_out,*)' fresid : rhor, approp, work(:,:)'
!  do ifft=1,n1
!  write(std_out,'(i4,3es18.8)' )ifft,rhor(ifft,1),approp(ifft),work(ifft,1)
!  end do
!  do ifft=1,nfft_tmp
!  if(work(ifft,1)<0.0_dp)then
!  write(std_out,*)' f_fft negative value',work(ifft,1),' for ifft=',ifft
!  end if
!  if(rhor(ifft,1)<0.0_dp)then
!  write(std_out,*)' rhor  negative value',rhor(ifft,1),' for ifft=',ifft
!  end if
!  end do
!  end if
!  end if
!  ENDDEBUG

   if(quit==1)exit

!  At the end of the first loop, where the appropriation function is generated,
!  invert it, to save cpu time later.
   if(iloop==1)approp(:)=1.0_dp/approp(:)

!  End first or second pass through atoms
 end do

!Restore proper normalisation of density
!(Non-collinear magnetism: n, mx,my,mz integral conservation)
 if(option==2)then
   if(mpi_enreg%paral_kgb==1) then
!    FFT parallelization: All the processors perform the same calculation.
!    We divided by nproc_fft in order to "remove" the xmpi_sum made in mean_fftr
!    Trangel: mpicomm now is optional in mean_fftr, no need to divide over nproc_fft
     call mean_fftr(rhor_tot,rhosum,nfft_tmp,nfftot,dtset%nspden)
!    rhosum(1:dtset%nspden)=rhosum(1:dtset%nspden)/mpi_enreg%nproc_fft
   else
     call mean_fftr(rhor,rhosum,nfft_tmp,nfftot,dtset%nspden,mpi_comm_sphgrid=mpi_enreg%comm_fft)
   end if
!  "!OCL NOPREEX" to avoid zero division after optimization (-Of) by MM
!  (Even if nspden=1, "1.0/rho_tot" will appear on vpp fujitsu
!  OCL NOPREEX
   if(mpi_enreg%paral_kgb==1) then
     do ispden=1,dtset%nspden
       fact=rhosum(ispden)*dble(nfftot)/rho_tot(ispden)
       work_tot(:,ispden)=fact*work_tot(:,ispden)
       call pre_scatter(work(:,ispden),work_tot(:,ispden),n1,n2,n3,n4,mpi_enreg)
     end do
   else
     do ispden=1,dtset%nspden
       fact=rhosum(ispden)*dble(nfftot)/rho_tot(ispden)
       work(:,ispden)=fact*work(:,ispden)
     end do
   end if
!  DEBUG
!  Here, zero all the hard work, for checking purposes !
!  work(:,:)=rhor(:,:)
!  ENDDEBUG
 end if

 ABI_DEALLOCATE(approp)
 ABI_DEALLOCATE(atmrho)
 ABI_DEALLOCATE(my_sphere)
 if(mpi_enreg%paral_kgb==1) then
   ABI_DEALLOCATE(rhor_tot)
   ABI_DEALLOCATE(work_tot)
 end if

!DEBUG
!write(std_out,*)' fresid : exit '
!do iatom=1,dtset%natom
!write(std_out,*)iatom,gresid(1:3,iatom)
!enddo
!ENDDEBUG

 contains

   function cross_fr(xx,yy,zz,aa,bb,cc)
!Define magnitude of cross product of two vectors
   real(dp) :: cross_fr
   real(dp),intent(in) :: xx,yy,zz,aa,bb,cc
   cross_fr=sqrt((yy*cc-zz*bb)**2+(zz*aa-xx*cc)**2+(xx*bb-yy*aa)**2)
 end function cross_fr

end subroutine fresid
!!***

!!****f* ABINIT/constrf
!! NAME
!! constrf
!!
!! FUNCTION
!! Computes projected forces, fpcart, which satisfy a set of
!! constraint equations of the form
!!  Sum[mu,iatom]: wtatcon(mu,iatom,iconeq)*fpcart(mu,iatom) = 0 (iconeq=1,nconeq).
!! These projected forces are returned in fcart and thus replace
!! the original forces.
!!
!! INPUTS
!!  iatfix(3,natom)=1 for frozen atom along each direction, 0 for unfrozen
!!  natom=number of atoms in cell
!!  nconeq=number of atomic constraint equations
!!  prtvol=control print volume and debugging
!!  rprimd(3,3)=dimensional primitive translations in real space (bohr)
!!  wtatcon(3,natom,nconeq)=weights for atomic constraints
!!  xred(3,natom)=reduced dimensionless atomic coordinates
!!
!! OUTPUT
!!  diffor=maximum absolute change in component of projected fcart between present
!!          and previous SCF cycle
!!  fred(3,natom)=grads of Etot wrt reduced coordinates (hartree)
!!  maxfor=maximum absolute value of fcart
!!
!! SIDE EFFECTS
!!  fcart(3,natom)=cartesian forces (hartree/bohr) on input, projected forces on output
!!  forold(3,natom)=cartesian forces of previous SCF cycle (hartree/bohr)
!!
!! TODO
!!
!! PARENTS
!!      forces
!!
!! CHILDREN
!!      dposv,prtxvf,wrtout,xred2xcart
!!
!! SOURCE

subroutine constrf(diffor,fcart,forold,fred,iatfix,ionmov,maxfor,natom,&
& nconeq,prtvol,rprimd,wtatcon,xred)

 use m_linalg_interfaces

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: ionmov,natom,nconeq,prtvol
 real(dp),intent(out) :: diffor,maxfor
!arrays
 integer,intent(in) :: iatfix(3,natom)
 real(dp),intent(in) :: rprimd(3,3),wtatcon(3,natom,nconeq)
 real(dp),intent(inout) :: fcart(3,natom),forold(3,natom),xred(3,natom)
 real(dp),intent(inout) :: fred(3,natom) !vz_i

!Local variables -------------------------
!scalars
 integer :: iatom,iconeq,iconeq1,iconeq2,index,info,mu,prtvel
 character(len=500) :: message
!arrays
 real(dp),allocatable :: fcartvec(:),fpcart(:,:),fvector(:),vel_dummy(:,:)
 real(dp),allocatable :: wmatrix(:,:),wtatconvec(:,:),wtcoeffs(:),xcart(:,:)

!************************************************************************

!Allocate temporary variables
 ABI_ALLOCATE(fpcart,(3,natom))
 ABI_ALLOCATE(fcartvec,(3*natom))
 ABI_ALLOCATE(fvector,(nconeq))
 ABI_ALLOCATE(vel_dummy,(3,natom))
 ABI_ALLOCATE(wmatrix,(nconeq,nconeq))
 ABI_ALLOCATE(wtatconvec,(3*natom,nconeq))
 ABI_ALLOCATE(wtcoeffs,(nconeq))
 ABI_ALLOCATE(xcart,(3,natom))

!If prtvol>10, output coordinates and forces prior to projecting
 if(prtvol>=10)then
   write(message,'(a)')' constrf - coordinates and forces prior to constraint projections:'
   call wrtout(std_out,message,'COLL')
   call xred2xcart(natom,rprimd,xcart,xred)
   prtvel=0
   call prtxvf(fcart,fred,iatfix,06,natom,prtvel,vel_dummy,xcart,xred)
 end if

!Transfer fcart and wtatcon to flat vectors
 index=0
 do iatom=1,natom
   do mu=1,3
     index=index+1
     fcartvec(index)=fcart(mu,iatom)
     wtatconvec(index,:)=wtatcon(mu,iatom,:)
   end do
 end do

!Compute a matrix (wmatrix) and vector (fvector) such that solving
!the linear equations wmatrix*wcoeffs=fvector gives the coefficients
!of wtatcon (wcoeffs) needed to compute the projected forces
 do iconeq2=1,nconeq
   fvector(iconeq2)=ddot(3*natom,fcartvec,1,wtatconvec(1,iconeq2),1)
   do iconeq1=1,nconeq
     wmatrix(iconeq1,iconeq2)=ddot(3*natom,wtatconvec(1,iconeq1),1,wtatconvec(1,iconeq2),1)
   end do
 end do

!Solve the system of linear equations, wmatrix*wcoeffs=fvector
 call dposv('U',nconeq,1,wmatrix,nconeq,fvector,nconeq,info)

 if (info/=0) then
   write(message, '(a,a,a,a,a)' )&
&   'Constraint matrix is not positive definite,',ch10,&
&   'probably because constraints are linearly dependent.',ch10,&
&   'Action: Check for linear dependence of constraints.'
   MSG_ERROR(message)
 end if

!The solution vector is returned in fvector, so copy it to a more sensible location
 wtcoeffs(:)=fvector(:)

!Compute the projected forces, which now satisfy all the constraint equations
 fpcart(:,:)=fcart(:,:)
 do iconeq=1,nconeq
   fpcart(:,:)=fpcart(:,:)-wtcoeffs(iconeq)*wtatcon(:,:,iconeq)
 end do

!Reconvert constrained forces back from fpcart to fred
 do iatom=1,natom
   do mu=1,3
     fred(mu,iatom)= - (rprimd(1,mu)*fpcart(1,iatom)+&
&     rprimd(2,mu)*fpcart(2,iatom)+&
&     rprimd(3,mu)*fpcart(3,iatom))
   end do
 end do

!If prtvol>=10, output coordinates and forces after projecting
 if(prtvol>=10)then
   write(message,'(a)')' constrf - coordinates and forces after constraint projections:'
   call wrtout(std_out,message,'COLL')
   prtvel=0
   call prtxvf(fpcart,fred,iatfix,06,natom,prtvel,vel_dummy,xcart,xred)
 end if

!Copy the constrained forces, fpcart, back to fcart
 fcart(:,:)=fpcart(:,:)

!Compute maximal force and maximal difference of the projected forces,
!overriding the values already computed in forces
 maxfor=0.0_dp
 diffor=0.0_dp
 do iatom=1,natom
   do mu=1,3
     if (iatfix(mu,iatom) /= 1) then
       maxfor=max(maxfor,abs(fcart(mu,iatom)))
       diffor=max(diffor,abs(fcart(mu,iatom)-forold(mu,iatom)))
     else if (ionmov==4 .or. ionmov==5) then
!      Make the force vanish on fixed atoms when ionmov=4 or 5
!      This is because fixing of atom cannot be imposed at the
!      level of a routine similar to brdmin or moldyn for these options.
       fcart(mu,iatom)=0.0_dp
     end if
   end do
 end do
 forold(:,:)=fcart(:,:)

!Dellocate temporary variables
 ABI_DEALLOCATE(fpcart)
 ABI_DEALLOCATE(fcartvec)
 ABI_DEALLOCATE(fvector)
 ABI_DEALLOCATE(vel_dummy)
 ABI_DEALLOCATE(wmatrix)
 ABI_DEALLOCATE(wtatconvec)
 ABI_DEALLOCATE(wtcoeffs)
 ABI_DEALLOCATE(xcart)

end subroutine constrf
!!***

end module m_forces
!!***
